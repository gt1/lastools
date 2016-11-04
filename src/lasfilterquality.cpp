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
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <db|dam> <out.las> <in.las>" << std::endl;
	ostr << "\n";

	return ostr.str();
}

int lasfilterquality(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & /* arginfo */)
{
	libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pfile(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(arg[2]));
	int64_t const tspace = Pfile->Palgn->tspace;
	libmaus2::dazzler::align::Overlap OVL;
	libmaus2::dazzler::align::AlignmentWriter AW(arg[1],tspace,true);
	libmaus2::dazzler::db::DatabaseFile DB(arg[0]);
	DB.computeTrimVector();
	libmaus2::dazzler::db::Track::unique_ptr_type Pqual(DB.readTrack("qual"));
	libmaus2::dazzler::db::TrackAnnoInterface const & qualanno = Pqual->getAnno();
	libmaus2::autoarray::AutoArray<unsigned char> const & qualdata = *(Pqual->Adata);
	unsigned int const defq = 25;
	unsigned int const q = arg.argPresent("q") ? arg.getUnsignedNumericArg<uint64_t>("q") : defq;

	libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > I;
	libmaus2::autoarray::AutoArray < std::pair< int32_t,int32_t > > O;
	uint64_t Io = 0;
	int64_t prevaread = -1;

	while ( Pfile->getNextOverlap(OVL) )
	{
		int64_t const aread = OVL.aread;

		if ( aread != prevaread )
		{
			Io = 0;

			int64_t const alow = qualanno[aread];
			int64_t const ahigh = qualanno[aread+1];
			unsigned char const * ua = qualdata.begin() + alow;
			unsigned char const * ue = qualdata.begin() + ahigh;

			unsigned char const * up = ua;
			while ( up != ue )
			{
				while ( up != ue && *up > q )
					++up;
				assert ( up == ue || *up <= q );

				unsigned char const * ul = up;
				while ( up != ue && *up <= q )
					++up;

				if ( up != ul )
					I.push(Io,std::pair<uint64_t,uint64_t>((ul-ua)*tspace,(up-ua)*tspace));
			}

			#if 0
			std::cerr << aread << "\t";
			for ( uint64_t i = 0; i < Io; ++i )
				std::cerr << "(" << I[i].first << "," << I[i].second << ")";
			std::cerr << std::endl;
			#endif

			prevaread = aread;
		}

		OVL.alignToTracePoints(tspace);
		uint64_t const o = OVL.filterIntervals(I.begin(),I.begin()+Io,O);

		for ( uint64_t i = 0; i < o; ++i )
		{
			libmaus2::dazzler::align::Overlap fOVL = OVL.filter(O[i],tspace);
			assert ( ! fOVL.isEmpty() );
			AW.put(fOVL);
		}
	}

	return EXIT_SUCCESS;
}

/**
 *
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
		else if ( arg.size() < 3 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return lasfilterquality(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
