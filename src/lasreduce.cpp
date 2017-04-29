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
#include <libmaus2/rank/ImpCacheLineRank.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " -p<prob=1> -s<seed> <out.las> <in.db> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int lasreduce(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
{
	std::string const outfilename = arg[0];
	std::string const dbname = arg[1];

	double const p = arg.uniqueArgPresent("p") ? arg.getParsedArg<double>("p") : 1.0;
	assert ( p <= 1.0 );

	if ( arg.uniqueArgPresent("s") )
	{
		uint64_t const seed = arg.getParsedArg<uint64_t>("s");
		std::cerr << "[V] using seed " << seed << std::endl;
		::libmaus2::random::Random::setup(seed);
	}
	else
	{
		::libmaus2::random::Random::setup();
	}

	libmaus2::dazzler::db::DatabaseFile DB(dbname);
	DB.computeTrimVector();

	std::vector<uint64_t> RL;
	DB.getAllReadLengths(RL);
	uint64_t const accrl = std::accumulate(RL.begin(),RL.end(),0ull);
	uint64_t const paccrl = static_cast<uint64_t>(std::floor(p * accrl + 0.5));
	assert ( paccrl < accrl );

	std::pair<int64_t,int64_t> P = DB.getTrimmedBlockInterval(0);
	uint64_t const num = P.second-P.first;
	assert ( P.first == 0 );

	uint64_t keep = static_cast<uint64_t>(p*num+0.5);
	if ( keep > num )
		keep = num;

	std::vector<bool> B(num,false);
	uint64_t nkeep = 0;
	uint64_t accsum = 0;

	while ( accsum < paccrl )
	{
		uint64_t const i = ::libmaus2::random::Random::rand64() % num;

		if ( ! B[i] )
		{
			B[i] = true;
			accsum += RL[i];
			nkeep++;
		}
	}
	
	libmaus2::rank::ImpCacheLineRank ICLR(num);
	libmaus2::rank::ImpCacheLineRank::WriteContext WC = ICLR.getWriteContext();
	for ( uint64_t i = 0; i < num; ++i )
		WC.writeBit(B[i]);
	WC.flush();
	
	for ( uint64_t i = 0; i < num; ++i )
		if ( B[i] )
		{
			std::string const data = DB[i];
			
			std::cout << ">L0/" << i << "/0_" << data.size() << " RQ=0.851\n";

			uint64_t p = 0;
			while ( p < data.size() )
			{
				uint64_t const rest = data.size()-p;
				uint64_t const cols = 80;
				uint64_t const print = std::min(cols,rest);
				
				std::cout << data.substr(p,print) << "\n";
				
				p += print;
			}
		}

	std::cerr << "[V] keeping " << static_cast<double>(nkeep) / num << " bases " << static_cast<double>(accsum)/accrl << std::endl;

	std::vector<std::string> VI;
	for ( uint64_t i = 2 ; i < arg.size(); ++i )
		VI.push_back(arg[i]);
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(VI);

	libmaus2::dazzler::align::AlignmentWriter::unique_ptr_type AW(
		new libmaus2::dazzler::align::AlignmentWriter(outfilename,tspace,false /* index */, 0 /* expt */)
	);

	libmaus2::dazzler::align::Overlap OVL;
	for ( uint64_t i = 0; i < VI.size(); ++i )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(VI[i]));

		while ( PIN->getNextOverlap(OVL) )
			if ( B[OVL.aread] && B[OVL.bread] )
			{
				uint64_t const r0 = ICLR.rankm1(OVL.aread);
				uint64_t const r1 = ICLR.rankm1(OVL.bread);
				OVL.aread = r0;
				OVL.bread = r1;
				AW->put(OVL);
			}
	}

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
		else if ( arg.size() < 2 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return lasreduce(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
