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

int lasbridge(libmaus2::util::ArgParser const & arg)
{
	std::string const outfn = arg[0];
	std::string const dbfn = arg[1];
	std::string const infn = arg[2];
	int64_t const border = arg.uniqueArgPresent("border") ? arg.getUnsignedNumericArg<uint64_t>("border") : 500;
	double const maxdiferr = arg.uniqueArgPresent("maxdiferr") ? arg.getParsedArg<double>("maxdiferr") : 0.3;

	int64_t const outtspace = libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();

	libmaus2::dazzler::align::AlignmentWriter AW(outfn,outtspace);

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

	libmaus2::dazzler::align::SimpleOverlapParser SOP(
		infn,1024*1024 /* buf size */,
		libmaus2::dazzler::align::OverlapParser::overlapparser_do_not_split_ab
	);

	while ( SOP.parseNextBlock() )
	{
		libmaus2::dazzler::align::OverlapData & data = SOP.getData();

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

				std::cerr << aread << " " << bread << " " << ihigh-ilow << " " << j << "/" << anumintv << "=[" << rfrom << "," << rto << ")" << std::endl;

				bool spanner = false;

				for ( uint64_t k = ilow; (!spanner) && k < ihigh; ++k )
				{
					int32_t const abpos = libmaus2::dazzler::align::OverlapData::getABPos(data.getData(k).first);
					int32_t const aepos = libmaus2::dazzler::align::OverlapData::getAEPos(data.getData(k).first);

					if ( abpos <= rfrom-border && aepos >= rto+border )
					{
						// std::cerr << "spanner " << abpos << "," << aepos << std::endl;
						spanner = true;
					}
				}

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
								std::string const sa = DB.decodeRead(aread,false);
								std::string const sb = DB.decodeRead(bread,isinv);

								libmaus2::lcs::NP np;

								np.np(
									sa.begin() + abpos,
									sa.begin() + aepos,
									sb.begin() + bbpos,
									sb.begin() + bepos
								);

								libmaus2::lcs::AlignmentTraceContainer::WindowErrorLargeResult const WELR =
									libmaus2::lcs::AlignmentTraceContainer::windowErrorLargeDetail(np.ta,np.te,500);

								std::cerr << "\t" << np.getAlignmentStatistics() << " " << WELR.maxerr << std::endl;

								#if 0
								libmaus2::lcs::AlignmentPrint::printAlignmentLines(
									std::cerr,
									sa.begin() + abpos,
									aepos-abpos,
									sb.begin() + bbpos,
									bepos-bbpos,
									80,
									np.ta,
									np.te
								);
								#endif

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

			ilow = ihigh;
		}
	}
	#if 0
		for ( uint64_t i = 0; i < DB.size(); ++i )
		{
			uint64_t const s = DB[i].size();

			unsigned char const * p = ptrack->Adata->begin() + low;

			assert ( size % (2*sizeof(int32_t)) == 0 );

			uint64_t const numintv = size / (2*sizeof(int32_t));

			for ( uint64_t j = 0; j < numintv; ++j )
			{
	#endif


	return EXIT_SUCCESS;
}

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [-d<reads.db> -b<reads.bam>] <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 1 )
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
