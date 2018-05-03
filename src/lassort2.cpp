/*
    lastools
    Copyright (C) 2015 German Tischler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <libmaus2/dazzler/align/LasSort2.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/dazzler/align/AlignmentFile.hpp>
#include <libmaus2/dazzler/align/SimpleOverlapParser.hpp>
#include <libmaus2/dazzler/align/OverlapDataInterface.hpp>
#include <libmaus2/aio/FileRemoval.hpp>
#include <libmaus2/aio/InputOutputStream.hpp>
#include <libmaus2/aio/InputOutputStreamFactoryContainer.hpp>
#include <libmaus2/aio/OutputStreamFactoryContainer.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [-M<memory> -l -t<numthreads> -T<tmpprefix> -f<mergefanin> -s<sortorder> --index] <out.las> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -t      : number of threads (defaults to number of cores on machine)\n";
	ostr << " -T      : prefix for temporary files (default: create files in current working directory)\n";
	ostr << " -f      : merge fan in (default: 64)\n";
	ostr << " -s      : sort order\n";
	ostr << " -l      : list mode (input contains list of file names instead of actual LAS files)\n";
	ostr << " -M      : memory block size\n";
	ostr << " --index : construct index for output file\n";

	return ostr.str();
}


std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}


template<typename comparator_type>
int lassort2Template(libmaus2::util::ArgParser const & arg)
{
	if ( arg.size() < 1 )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] usage: " << arg.progname << " <out.las> <in.las> ..." << std::endl;
		lme.finish();
		throw lme;
	}

	std::string const outlas = arg[0];
	std::string const tmpbase = getTmpFileBase(arg);
	uint64_t const blocksize = arg.uniqueArgPresent("M") ? arg.getUnsignedNumericArg<uint64_t>("M") : (1024ull * 1024ull * 1024ull);
	bool const listmode = arg.uniqueArgPresent("l");
	uint64_t const fanin = 16;
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();

	std::vector < std::string > Vin;

	if ( listmode )
	{
		for ( uint64_t z = 1; z < arg.size(); ++z )
		{
			libmaus2::aio::InputStreamInstance ISI(arg[z]);

			while ( ISI )
			{
				std::string s;
				ISI >> s;

				if ( s.size() )
					Vin.push_back(s);
			}
		}

		for ( uint64_t i = 0; i < Vin.size(); ++i )
			std::cerr << "[L]\t" << Vin[i] << std::endl;
	}
	else
	{
		for ( uint64_t z = 1; z < arg.size(); ++z )
			Vin.push_back(arg[z]);
	}

	uint64_t tmpid = 0;
	bool deletein = false;
	while ( (numthreads > 1) && (Vin.size() > fanin) )
	{
		uint64_t const tsplit = (Vin.size() + numthreads -1) / numthreads;
		uint64_t const two = 2;
		uint64_t const split = std::max(two,tsplit);
		uint64_t const parts = (Vin.size() + split - 1)/split;

		std::vector < std::string > Vtmp(parts);
		for ( uint64_t i = 0; i < parts; ++i )
		{
			std::ostringstream ostr;
			ostr << tmpbase << "_threadmerge_" << (tmpid++);
			Vtmp[i] = ostr.str();
		}
		std::vector < std::string > Nin(parts);

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < parts; ++i )
		{
			uint64_t const low = i * split;
			uint64_t const high = std::min(low+split,static_cast<uint64_t>(Vin.size()));
			std::string const outfn = Vtmp[i] + "_out";
			Nin[i] = outfn;

			std::ostringstream ostr;
			libmaus2::dazzler::align::LasSort2<comparator_type>::lassort2(
				outfn,
				std::vector<std::string>(
					Vin.begin()+low,
					Vin.begin()+high
				), blocksize, fanin, Vtmp[i], false, &ostr);

			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << ostr.str();
		}

		if ( deletein )
			for ( uint64_t i = 0; i < Vin.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile(Vin[i]);

		deletein = true;
		Vin = Nin;
	}

	bool const index = arg.uniqueArgPresent("index");

	return libmaus2::dazzler::align::LasSort2<comparator_type>::lassort2(outlas, Vin, blocksize, fanin, tmpbase, index, &std::cerr);

	if ( deletein )
		for ( uint64_t i = 0; i < Vin.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(Vin[i]);
}

int lassort2(libmaus2::util::ArgParser const & arg)
{
	std::string const sortorder = arg.uniqueArgPresent("s") ? arg["s"] : std::string("canonical");

	if ( sortorder == "canonical" )
		return lassort2Template<libmaus2::dazzler::align::OverlapDataInterfaceFullComparator>(arg);
	else if ( sortorder == "ba" )
		return lassort2Template<libmaus2::dazzler::align::OverlapDataInterfaceBAComparator>(arg);
	else
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] unknown sort order " << sortorder << std::endl;
		lme.finish();
		throw lme;
	}
}

/**
 * sort a set of LAS/OVL files and merge the sorted files to a single output file
 **/

int main(int argc, char *argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.argPresent("h") || arg.argPresent("help") )
		{
			std::cerr << getUsage(arg);
			return EXIT_SUCCESS;
		}
		else if ( arg.argPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.size() < 2 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return lassort2(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
