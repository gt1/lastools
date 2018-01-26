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
#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/random/Random.hpp>
#include <libmaus2/bitio/BitVector.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <out.las> <out.db> <in.db> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int lassubset(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
{
	std::string const outfilename = arg[0];
	std::string const outfasta = arg[1];
	std::string const dbname = arg[2];

	std::vector<std::string> Vin;
	for ( uint64_t i = 3; i < arg.size(); ++i )
		Vin.push_back(arg[i]);

	libmaus2::dazzler::db::DatabaseFile DB(dbname);
	DB.computeTrimVector();
	uint64_t const n = DB.size();

	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vin);

	libmaus2::bitio::BitVector BV(n);
	for ( uint64_t i = 0; i < n; ++i )
		BV.erase(i);

	for ( uint64_t i = 0; i < Vin.size(); ++i )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(Vin[i]));
		libmaus2::dazzler::align::Overlap OVL;

		while ( PIN->getNextOverlap(OVL) )
		{
			BV.set(OVL.aread);
			BV.set(OVL.bread);
		}
	}

	libmaus2::dazzler::align::AlignmentWriter::unique_ptr_type AW(
		new libmaus2::dazzler::align::AlignmentWriter(outfilename,tspace,false /* index */, 0 /* expt */)
	);

	libmaus2::rank::ImpCacheLineRank ICLR(n);
	libmaus2::rank::ImpCacheLineRank::WriteContext context = ICLR.getWriteContext();
	for ( uint64_t i = 0; i < n; ++i )
		context.writeBit(BV.get(i));
	context.flush();

	for ( uint64_t i = 0; i < Vin.size(); ++i )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(Vin[i]));
		libmaus2::dazzler::align::Overlap OVL;

		while ( PIN->getNextOverlap(OVL) )
		{
			OVL.aread = ICLR.rank1(OVL.aread)-1;
			OVL.bread = ICLR.rank1(OVL.bread)-1;
			AW->put(OVL);
		}
	}

	libmaus2::aio::OutputStreamInstance OSI(outfasta);
	for ( uint64_t i = 0; i < n; ++i )
		if ( BV.get(i) )
		{
			std::string const r = DB[i];
			OSI << ">L0/" << i << "/" << 0 << "_" << r.size() << "\n";

			char const * p = r.c_str();
			char const * pe = p + r.size();

			while ( p != pe )
			{
				uint64_t const rest = pe-p;
				uint64_t const cols = 80;
				uint64_t const toprint = std::min(rest,cols);

				OSI.write(p,toprint);
				OSI.put('\n');

				p += toprint;
			}
		}
	OSI.flush();

	AW.reset();

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
		else if ( arg.size() < 3 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return lassubset(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
