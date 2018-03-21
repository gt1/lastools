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

	ostr << "usage: " << arg.progname << " [-t<numthreads>] [--tspace<value>] <in_a.db> <in_b.db> <in.las>" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int laschangetspace(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
{
	std::string const db0name = arg[0];
	std::string const db1name = arg[1];

	uint64_t const maxovl = arg.uniqueArgPresent("maxovl") ? arg.getUnsignedNumericArg<uint64_t>("maxovl") : 64*1024;
	bool const verbose = arg.uniqueArgPresent("verbose");

	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();

	std::cerr << "[V] loading data for " << db0name << " to memory...";
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type Pdb0data(
		libmaus2::dazzler::db::DatabaseFile::copyToArrays(db0name)
	);
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet const * db0data = Pdb0data.get();
	std::cerr << "done." << std::endl;

	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type Pdb1data;
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet const * db1data = 0;
	if ( db1name != db0name )
	{
		std::cerr << "[V] loading data for " << db1name << " to memory...";
		libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type Tdb1data(
			libmaus2::dazzler::db::DatabaseFile::copyToArrays(db1name)
		);
		Pdb1data = UNIQUE_PTR_MOVE(Tdb1data);
		db1data = Pdb1data.get();
		std::cerr << "done." << std::endl;
	}
	else
	{
		db1data = Pdb0data.get();
	}

	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB0(
		new libmaus2::dazzler::db::DatabaseFile(db0data->getDBURL())
	);
	libmaus2::dazzler::db::DatabaseFile & DB0 = *PDB0;
	DB0.computeTrimVector();
	std::vector<uint64_t> RL0;
	DB0.getAllReadLengths(RL0);

	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB1;
	libmaus2::dazzler::db::DatabaseFile * pDB1 = 0;
	if ( db1name != db0name )
	{
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type TDB1(
			new libmaus2::dazzler::db::DatabaseFile(db1data->getDBURL())
		);
		PDB1 = UNIQUE_PTR_MOVE(TDB1);
		PDB1->computeTrimVector();
		pDB1 = PDB1.get();
	}
	else
	{
		pDB1 = PDB0.get();
	}
	libmaus2::dazzler::db::DatabaseFile & DB1 = *pDB1;

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

	libmaus2::dazzler::align::Overlap OVL;
	double volatile maxerr = 0;
	for ( uint64_t i = 2; i < arg.size(); ++i )
	{
		std::cerr << "[V] processing " << arg[i] << std::endl;

		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(arg[i]));
		int64_t const intspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(arg[i]);
		int64_t const novl = libmaus2::dazzler::align::AlignmentFile::getNovl(arg[i]);

		uint64_t o = 0;
		std::vector < libmaus2::dazzler::align::Overlap > VOVL(maxovl);
		std::vector < uint64_t > Vaid(maxovl);
		std::vector < uint64_t > Vbid(maxovl);

		uint64_t proc = 0;

		do
		{
			o = 0;

			// std::cerr << "[V] loading alignments...";
			while ( o < maxovl && PIN->getNextOverlap(OVL) )
			{
				Vaid[o] = OVL.aread;
				Vbid[o] = OVL.bread;
				VOVL[o] = OVL;
				++o;
			}
			// std::cerr << "got " << o << std::endl;

			std::sort(Vaid.begin(),Vaid.begin()+o);
			uint64_t const numa = std::unique(Vaid.begin(),Vaid.begin()+o) - Vaid.begin();
			std::sort(Vbid.begin(),Vbid.begin()+o);
			uint64_t const numb = std::unique(Vbid.begin(),Vbid.begin()+o) - Vbid.begin();

			// std::cerr << "[V] numa=" << numa << " numb=" << numb << std::endl;

			libmaus2::dazzler::db::DatabaseFile::ReadDataRange::unique_ptr_type adata(
				DB0.decodeReadDataByArrayParallelDecode(Vaid.begin(),numa,numthreads,false /* rc */,0 /* term */));
			libmaus2::dazzler::db::DatabaseFile::ReadDataRange::unique_ptr_type bdata_forward(
				DB1.decodeReadDataByArrayParallelDecode(Vbid.begin(),numb,numthreads,false /* rc */,0 /* term */));
			libmaus2::dazzler::db::DatabaseFile::ReadDataRange::unique_ptr_type bdata_reverse(
				DB1.decodeReadDataByArrayParallelDecode(Vbid.begin(),numb,numthreads,true /* rc */,0 /* term */));

			// std::cerr << "[V] decoded read data" << std::endl;

			int volatile gfailed = 0;
			libmaus2::parallel::PosixSpinLock gfailedlock;

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

				std::vector<uint64_t>::const_iterator it_a = ::std::lower_bound(Vaid.begin(),Vaid.begin()+numa,VOVL[j].aread);
				assert ( it_a != Vaid.begin()+numa && static_cast<int64_t>(*it_a) == VOVL[j].aread );
				std::vector<uint64_t>::const_iterator it_b = ::std::lower_bound(Vbid.begin(),Vbid.begin()+numb,VOVL[j].bread);
				assert ( it_b != Vbid.begin()+numb && static_cast<int64_t>(*it_b) == VOVL[j].bread );

				uint64_t const arank = it_a - Vaid.begin();
				uint64_t const brank = it_b - Vbid.begin();

				uint8_t const * ua = adata->get(arank).first;
				uint8_t const * ub = VOVL[j].isInverse() ?
					bdata_reverse->get(brank).first
					:
					bdata_forward->get(brank).first;

				int64_t const l_a = adata->get(arank).second - adata->get(arank).first;
				int64_t const l_b = bdata_forward->get(brank).second - bdata_forward->get(brank).first;
				bool const contained_A = VOVL[j].path.abpos == 0 && VOVL[j].path.aepos == l_a;
				bool const contained_B = VOVL[j].path.bbpos == 0 && VOVL[j].path.bepos == l_b;

				try
				{
					VOVL[j].computeTrace(ua,ub,intspace,ATC,NP);

					if ( verbose )
					{
						{
							libmaus2::lcs::AlignmentStatistics const AS = ATC.getAlignmentStatistics();
							double const lerr = AS.getErrorRate();

							libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
							std::cerr << "[V] " << arg[i] << " " << AS;

							if ( lerr >= 0.45 )
								std::cerr << " HIGH";
							if ( contained_A )
								std::cerr << " C_A";
							if ( contained_B )
								std::cerr << " C_B";

							std::cerr << std::endl;
							if ( lerr > maxerr )
								maxerr = lerr;
						}
					}
				}
				catch(std::exception const & ex)
				{
					{
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						std::cerr << ex.what() << std::endl;
					}

					gfailedlock.lock();
					gfailed = 1;
					gfailedlock.unlock();
				}
			}

			if ( gfailed )
			{
				std::cerr << "[V] failed" << std::endl;
				return EXIT_FAILURE;
			}

			proc += o;

			// std::cerr << "[V] " << proc << " / " << novl << std::endl;
		} while ( o );
	}

	std::cout << "[V] ok" << std::endl;
	if ( verbose )
		std::cerr << "[V] maxerr=" << maxerr << std::endl;

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
