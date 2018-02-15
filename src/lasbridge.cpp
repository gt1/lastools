/*
    lastools
    Copyright (C) 2016 German Tischler

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

#include <config.h>

#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/align/SimpleOverlapParser.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/bambam/DecoderBase.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

int lasbridge(libmaus2::util::ArgParser const & arg)
{
	std::string const outfn = arg[0];
	std::string const dbfn = arg[1];
	std::string const infn = arg[2];
	// border around tandems considered in base pairs
	int64_t const border = arg.uniqueArgPresent("border") ? arg.getUnsignedNumericArg<uint64_t>("border") : 500;
	// maximum indel error accepted between A and B read in [abpos,aepos],[bbpos,bepos]
	double const maxdiferr = arg.uniqueArgPresent("maxdiferr") ? arg.getParsedArg<double>("maxdiferr") : 0.3;
	// output tspace value
	int64_t const outtspace = arg.uniqueArgPresent("tspace") ? arg.getUnsignedNumericArg<uint64_t>("tspace") : libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
	// maximum error allowed for 500 bp windows on output alignments
	double const maxwindowerror = arg.uniqueArgPresent("maxwindowerror") ? arg.getParsedArg<double>("maxdiferr") : 0.4;
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	bool const dbmem = arg.uniqueArgPresent("dbmem");

	// libmaus2::dazzler::align::AlignmentWriter AW(outfn,outtspace);

	libmaus2::dazzler::align::AlignmentWriterArray AWA(outfn + "_AWA", numthreads, outtspace);

	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	DB.computeTrimVector();

	libmaus2::dazzler::db::Track::unique_ptr_type ptrack(DB.readTrack("vtan"));

	if ( ptrack->trackannosize )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[V] track " << "vtan" << " does not look like a mask track" << std::endl;
		lme.finish();
		throw lme;
	}

	libmaus2::dazzler::db::TrackAnnoInterface const & anno = ptrack->getAnno();

	assert ( anno.size() == DB.size()+1 );

	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type dbptr;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB;
	libmaus2::dazzler::db::DatabaseFile * DBP = &DB;

	if ( dbmem )
	{
		std::cerr << "[V] copying " << dbfn << " to memory...";
		libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type tdbptr(
			libmaus2::dazzler::db::DatabaseFile::copyToArrays(dbfn)
		);
		dbptr = UNIQUE_PTR_MOVE(tdbptr);

		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type tPDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->getDBURL()));
		PDB = UNIQUE_PTR_MOVE(tPDB);
		PDB->computeTrimVector();
		
		DBP = PDB.get();

		std::cerr << "done." << std::endl;
	}


	libmaus2::dazzler::align::SimpleOverlapParser SOP(
		infn,1024*1024 /* buf size */,
		libmaus2::dazzler::align::OverlapParser::overlapparser_do_not_split_ab
	);

	while ( SOP.parseNextBlock() )
	{
		libmaus2::dazzler::align::OverlapData & data = SOP.getData();

		std::vector < std::pair<uint64_t,uint64_t > > VI;

		uint64_t ilow = 0;
		while ( ilow < data.size() )
		{
			int32_t const aread = libmaus2::dazzler::align::OverlapData::getARead(data.getData(ilow).first);
			int32_t const bread = libmaus2::dazzler::align::OverlapData::getBRead(data.getData(ilow).first);
			bool const isinv = libmaus2::dazzler::align::OverlapData::getInverseFlag(data.getData(ilow).first);

			uint64_t ihigh = ilow+1;

			while (
				ihigh < data.size()
				&&
				libmaus2::dazzler::align::OverlapData::getARead(data.getData(ihigh).first) == aread
				&&
				libmaus2::dazzler::align::OverlapData::getBRead(data.getData(ihigh).first) == bread
				&&
				libmaus2::dazzler::align::OverlapData::getInverseFlag(data.getData(ihigh).first) == isinv
			)
				++ihigh;

			VI.push_back(std::pair<uint64_t,uint64_t >(ilow,ihigh));

			ilow = ihigh;
		}

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t z = 0; z < VI.size(); ++z )
		{
			std::ostringstream logstr;

			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			libmaus2::dazzler::align::AlignmentWriter & AW = AWA[tid];

			uint64_t const ilow = VI[z].first;
			uint64_t const ihigh = VI[z].second;

			int32_t const aread = libmaus2::dazzler::align::OverlapData::getARead(data.getData(ilow).first);
			int32_t const bread = libmaus2::dazzler::align::OverlapData::getBRead(data.getData(ilow).first);
			bool const isinv = libmaus2::dazzler::align::OverlapData::getInverseFlag(data.getData(ilow).first);

			uint64_t const alow  = anno[aread];
			uint64_t const ahigh = anno[aread+1];
			uint64_t const asize = ahigh-alow;

			unsigned char const * p = ptrack->Adata->begin() + alow;
			assert ( asize % (2*sizeof(int32_t)) == 0 );
			uint64_t const anumintv = asize / (2*sizeof(int32_t));

			for ( uint64_t j = 0; j < anumintv; ++j )
			{
				int64_t const rfrom = libmaus2::bambam::DecoderBase::getLEInteger(p + (j * 2 + 0) * sizeof(uint32_t), sizeof(uint32_t));
				int64_t const rto = libmaus2::bambam::DecoderBase::getLEInteger(p + (j * 2 + 1) * sizeof(uint32_t), sizeof(uint32_t));

				logstr << aread << " " << bread << " " << ihigh-ilow << " " << j << "/" << anumintv << "=[" << rfrom << "," << rto << ")" << std::endl;

				bool spanner = false;

				// check whether there is a spanning read
				for ( uint64_t k = ilow; (!spanner) && k < ihigh; ++k )
				{
					int32_t const abpos = libmaus2::dazzler::align::OverlapData::getABPos(data.getData(k).first);
					int32_t const aepos = libmaus2::dazzler::align::OverlapData::getAEPos(data.getData(k).first);

					if ( abpos <= rfrom-border && aepos >= rto+border )
					{
						spanner = true;
					}
				}

				// if not then check whether we find something anchored on the left and something on the right
				if ( !spanner )
				{
					std::vector < uint64_t > Vanchorleft;
					std::vector < uint64_t > Vanchorright;

					for ( uint64_t k = ilow; k < ihigh; ++k )
					{
						int32_t const abpos = libmaus2::dazzler::align::OverlapData::getABPos(data.getData(k).first);
						int32_t const aepos = libmaus2::dazzler::align::OverlapData::getAEPos(data.getData(k).first);

						if ( abpos <= rfrom-border && aepos >= rfrom+border )
						{
							Vanchorleft.push_back(k);
						}
						else if ( abpos <= rto-border && aepos >= rto+border )
						{
							Vanchorright.push_back(k);
						}
					}

					for ( uint64_t i = 0; i < Vanchorleft.size(); ++i )
						for ( uint64_t j = 0; j < Vanchorright.size(); ++j )
						{
							uint64_t const ileft = Vanchorleft[i];
							uint64_t const iright = Vanchorright[j];

							int32_t const abpos = libmaus2::dazzler::align::OverlapData::getABPos(data.getData(ileft).first);
							int32_t const aepos = libmaus2::dazzler::align::OverlapData::getAEPos(data.getData(iright).first);
							int32_t const bbpos = libmaus2::dazzler::align::OverlapData::getBBPos(data.getData(ileft).first);
							int32_t const bepos = libmaus2::dazzler::align::OverlapData::getBEPos(data.getData(iright).first);

							int64_t const adif = aepos-abpos;
							int64_t const bdif = bepos-bbpos;
							int64_t const maxdif = std::max(adif,bdif);
							int64_t const mindif = std::min(adif,bdif);
							int64_t const ddif = maxdif-mindif;
							double const e = static_cast<double>(ddif) / mindif;


							if ( e <= maxdiferr )
							{
								std::string const sa = DBP->decodeRead(aread,false);
								std::string const sb = DBP->decodeRead(bread,isinv);

								libmaus2::lcs::NP np;

								np.np(
									sa.begin() + abpos,
									sa.begin() + aepos,
									sb.begin() + bbpos,
									sb.begin() + bepos
								);

								libmaus2::lcs::AlignmentTraceContainer::WindowErrorLargeResult const WELR =
									libmaus2::lcs::AlignmentTraceContainer::windowErrorLargeDetail(np.ta,np.te,500);
								libmaus2::lcs::AlignmentStatistics const ASw =
									libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(WELR.t0,WELR.t1);
								double const werr = ASw.getErrorRate();

								logstr << "\t"
									<< np.getAlignmentStatistics()
									<< " "
									<< ASw
									<< std::endl;

								#if 0
								libmaus2::lcs::AlignmentPrint::printAlignmentLines(
									logstr,
									sa.begin() + abpos,
									aepos-abpos,
									sb.begin() + bbpos,
									bepos-bbpos,
									80,
									np.ta,
									np.te
								);
								#endif

								if ( werr <= maxwindowerror )
								{
									libmaus2::dazzler::align::Overlap const OVL = libmaus2::dazzler::align::Overlap::computeOverlap(
										libmaus2::dazzler::align::OverlapData::getFlags(data.getData(ileft).first),
										aread,
										bread,
										abpos,
										aepos,
										bbpos,
										bepos,
										outtspace,
										np.getTraceContainer());

									AW.put(OVL);
								}
							}
						}
				}
			}

			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << logstr.str() << std::flush;
		}
	}

	AWA.merge(outfn,outfn + ".tmp");

	return EXIT_SUCCESS;
}

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " out.las in.db in.las" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	
	ostr << "\t--border<int>: border around tandems considered in base pairs (default 500)" << std::endl;
	ostr << "\t--maxdifferr<double>: maximum indel error accepted between A and B read in [abpos,aepos],[bbpos,bepos] (default 0.3)" << std::endl;
	ostr << "\t--outtspace<int>: output tspace (default " << libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace() << ")" << std::endl;
	ostr << "\t--maxwindowerror<double>: maximum error allowed for 500 bp windows on output alignments (default 0.4)" << std::endl;
	ostr << "\t-t<int>: number of threads used (defaults to number of CPU cores detected)" << std::endl;
	ostr << "\t--dbmem: copy dazzler database to memory for processing" << std::endl;

	return ostr.str();
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 3 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << getUsage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return lasbridge(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
