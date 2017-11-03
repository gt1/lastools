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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/lcs/NP.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [-t<numthreads>] [--tspace<value>] <out.las> <in_a.db> <in_b.db> <in.las>" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int64_t getDefaultTSpace()
{
	return libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
}

int laschangetspace(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
{
	std::string const outlas = arg[0];
	std::string const db0name = arg[1];
	std::string const db1name = arg[2];

	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	int64_t const outtspace = arg.uniqueArgPresent("tspace") ? arg.getUnsignedNumericArg<uint64_t>("tspace") : getDefaultTSpace();

	libmaus2::dazzler::db::DatabaseFile DB0(db0name);
	DB0.computeTrimVector();
	std::vector<uint64_t> RL0;
	DB0.getAllReadLengths(RL0);

	libmaus2::dazzler::db::DatabaseFile DB1(db1name);
	DB1.computeTrimVector();
	std::vector<uint64_t> RL1;
	DB1.getAllReadLengths(RL1);

	libmaus2::autoarray::AutoArray<libmaus2::lcs::AlignmentTraceContainer::unique_ptr_type> AATC(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::lcs::AlignmentTraceContainer::unique_ptr_type tptr(new libmaus2::lcs::AlignmentTraceContainer);
		AATC[i] = UNIQUE_PTR_MOVE(tptr);
	}
	libmaus2::autoarray::AutoArray<libmaus2::lcs::NP::unique_ptr_type> ANP(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::lcs::NP::unique_ptr_type tptr(new libmaus2::lcs::NP);
		ANP[i] = UNIQUE_PTR_MOVE(tptr);
	}

	libmaus2::dazzler::db::DatabaseFile::ReadDataRange::unique_ptr_type adata(
		DB0.decodeReadIntervalParallel(0,RL0.size(),numthreads,false /* rc */,0 /* term */));
	libmaus2::dazzler::db::DatabaseFile::ReadDataRange::unique_ptr_type bdata_forward(
		DB1.decodeReadIntervalParallel(0,RL1.size(),numthreads,false /* rc */,0 /* term */));
	libmaus2::dazzler::db::DatabaseFile::ReadDataRange::unique_ptr_type bdata_reverse(
		DB1.decodeReadIntervalParallel(0,RL1.size(),numthreads,true /* rc */,0 /* term */));

	libmaus2::dazzler::align::AlignmentWriter AW(outlas,outtspace);

	libmaus2::dazzler::align::Overlap OVL;
	for ( uint64_t i = 3; i < arg.size(); ++i )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(arg[i]));
		int64_t const intspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(arg[i]);
		int64_t const novl = libmaus2::dazzler::align::AlignmentFile::getNovl(arg[i]);

		uint64_t maxovl = 256*1024;
		uint64_t o = 0;
		std::vector < libmaus2::dazzler::align::Overlap > VOVL(maxovl);

		uint64_t proc = 0;

		do
		{
			o = 0;

			while ( o < maxovl && PIN->getNextOverlap(OVL) )
				VOVL[o++] = OVL;

			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
			#endif
			for ( uint64_t j = 0; j < o; ++j )
			{
				#if defined(_OPENMP)
				uint64_t const tid = omp_get_thread_num();
				#else
				uint64_t const tid = 0;
				#endif

				libmaus2::lcs::AlignmentTraceContainer & ATC = *(AATC[tid]);
				libmaus2::lcs::NP & NP = *(ANP[tid]);

				uint8_t const * ua = adata->get(VOVL[j].aread).first;
				uint8_t const * ub = VOVL[j].isInverse() ?
					bdata_reverse->get(VOVL[j].bread).first
					:
					bdata_forward->get(VOVL[j].bread).first;

				VOVL[j].computeTrace(ua,ub,intspace,ATC,NP);

				VOVL[j] = libmaus2::dazzler::align::Overlap::computeOverlap(
					VOVL[j].flags,
					VOVL[j].aread,VOVL[j].bread,
					VOVL[j].path.abpos,VOVL[j].path.aepos,
					VOVL[j].path.bbpos,VOVL[j].path.bepos,
					outtspace,
					ATC
				);
			}

			for ( uint64_t j = 0; j < o; ++j )
				AW.put(VOVL[j]);

			proc += o;

			std::cerr << "[V] " << proc << " / " << novl << std::endl;
		} while ( o );
	}

	return EXIT_SUCCESS;
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

		return laschangetspace(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
