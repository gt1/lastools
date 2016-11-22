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
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/GetFileSize.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [-t<threads> -T<tmpprefix> -f<mergefanin> -s<sortorder>] <out.las> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -T : prefix for temporary files (default: create files in current working directory)\n";
	ostr << " -f : merge fan in (default: 64)\n";
	ostr << " -s : sort order (canonical or ba, default: canonical)\n";

	return ostr.str();
}

template<typename comparator_type>
int lasmergeTemplate(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	comparator_type comp;
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : arginfo.getDefaultTmpFileName();
	uint64_t const mergefanin = arg.uniqueArgPresent("f") ? arg.getUnsignedNumericArg<uint64_t>("f") : 64;
	bool const iflag = arg.uniqueArgPresent("i");
	uint64_t nextid = 0;

	std::string const outfilename = arg[0];
	std::vector<std::string> infilenames;

	if ( iflag )
	{
		for ( uint64_t i = 1; i < arg.size(); ++i )
		{
			libmaus2::aio::InputStreamInstance ISI(arg[i]);
			std::string line;
			while ( ISI )
			{
				std::getline(ISI,line);
				if ( line.size() )
					infilenames.push_back(line);
			}
		}
	}
	else
	{
		for ( uint64_t i = 1; i < arg.size(); ++i )
			infilenames.push_back(arg[i]);
	}

	assert ( infilenames.size() );

	if ( infilenames.size() > 1 )
	{
		assert ( mergefanin > 1 );
		bool deleteinput = false;

		while ( infilenames.size() > mergefanin )
		{
			uint64_t const packs = ( infilenames.size() + mergefanin - 1 ) / mergefanin;
			uint64_t const tmergefanin = ( infilenames.size() + packs - 1 ) / packs;

			uint64_t s = 0;
			std::vector<std::string> IV(packs);
			for ( uint64_t z = 0; z < packs; ++z )
			{
				uint64_t const low  = z * tmergefanin;
				uint64_t const high = std::min(low+tmergefanin,infilenames.size());
				assert ( high > low );
				s += high - low;

				std::vector<std::string> I(infilenames.begin()+low,infilenames.begin()+high);
				std::ostringstream fnostr;
				fnostr << tmpfilebase << "_" << std::setw(6) << std::setfill('0') << (nextid++);
				std::string const O = fnostr.str();
				libmaus2::util::TempFileRemovalContainer::addTempFile(O);

				libmaus2::dazzler::align::SortingOverlapOutputBuffer<comparator_type>::mergeFiles(I,O,comp);

				if ( deleteinput )
					for ( uint64_t i = 0; i < I.size(); ++i )
						libmaus2::dazzler::align::SortingOverlapOutputBuffer<>::removeFileAndIndex(I[i]);

				IV[z] = O;
			}
			assert ( s == infilenames.size() );

			infilenames = IV;
			deleteinput = true;
		}

		assert ( infilenames.size() <= mergefanin );

		libmaus2::dazzler::align::SortingOverlapOutputBuffer<comparator_type>::mergeFiles(infilenames,outfilename,comp);

		if ( deleteinput )
			for ( uint64_t i = 0; i < infilenames.size(); ++i )
				libmaus2::aio::FileRemoval::removeFile(infilenames[i]);
	}
	else
	{
		libmaus2::util::GetFileSize::copy(infilenames[0],outfilename);
	}

	return EXIT_SUCCESS;
}

int lasmerge(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	std::string const sortorder = arg.uniqueArgPresent("s") ? arg["s"] : std::string("canonical");

	if ( sortorder == "ba" )
		return lasmergeTemplate<libmaus2::dazzler::align::OverlapComparatorBReadARead>(arg,arginfo);
	else if ( sortorder == "aidbid" )
		return lasmergeTemplate<libmaus2::dazzler::align::OverlapComparatorAIdBId>(arg,arginfo);
	else if ( sortorder == "bidaid" )
		return lasmergeTemplate<libmaus2::dazzler::align::OverlapComparatorBIdAId>(arg,arginfo);
	else
		return lasmergeTemplate<libmaus2::dazzler::align::OverlapComparator>(arg,arginfo);
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

		return lasmerge(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
