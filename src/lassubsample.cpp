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

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " -p<prob=1> -s<seed> <out.las> <in.db> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int lassort(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
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

	std::pair<int64_t,int64_t> P = DB.getTrimmedBlockInterval(0);
	uint64_t const num = P.second-P.first;
	assert ( P.first == 0 );

	uint64_t keep = static_cast<uint64_t>(p*num+0.5);
	if ( keep > num )
		keep = num;

	std::vector<bool> B(num,false);
	uint64_t nkeep = 0;

	while ( nkeep < keep )
	{
		uint64_t const i = ::libmaus2::random::Random::rand64() % num;

		if ( ! B[i] )
		{
			B[i] = true;
			nkeep++;
		}
	}

	std::cerr << "[V] keeping " << static_cast<double>(nkeep) / num << std::endl;

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
				AW->put(OVL);
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

		return lassort(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
