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

	ostr << "usage: " << arg.progname << " <in_a.db> <in_b.db> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int lascheck(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
{
	std::string const db0name = arg[0];
	std::string const db1name = arg[1];

	libmaus2::dazzler::db::DatabaseFile DB0(db0name);
	DB0.computeTrimVector();
	std::vector<uint64_t> RL0;
	DB0.getAllReadLengths(RL0);

	libmaus2::dazzler::db::DatabaseFile DB1(db1name);
	DB1.computeTrimVector();
	std::vector<uint64_t> RL1;
	DB1.getAllReadLengths(RL1);

	libmaus2::dazzler::align::Overlap OVL;
	for ( uint64_t i = 2; i < arg.size(); ++i )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(arg[i]));

		bool ok = true;

		while ( PIN->getNextOverlap(OVL) )
		{
			std::ostringstream ostr;
			if ( OVL.path.abpos < 0 )
				ostr << "[E] broken abpos " << OVL.path.abpos << std::endl;
			if ( OVL.path.aepos > static_cast<int64_t>(RL0[OVL.aread]) )
				ostr << "[E] broken aepos " << OVL.path.aepos << " > " << RL0[OVL.aread] << std::endl;
			if ( OVL.path.bbpos < 0 )
				ostr << "[E] broken bbpos " << OVL.path.bbpos << std::endl;
			if ( OVL.path.bepos > static_cast<int64_t>(RL1[OVL.bread]) )
				ostr << "[E] broken bepos " << OVL.path.bepos << " > " << RL1[OVL.bread] << std::endl;
			if ( !OVL.checkBSpan() )
				ostr << "[E] bspan broken " << OVL.path.getBSpan() << " != " << OVL.path.bepos-OVL.path.bbpos << std::endl;

			if ( ostr.str().size() )
			{
				std::cerr << OVL << std::endl;
				std::cerr << ostr.str();
				ok = false;
			}
		}

		std::cout << arg[i] << " " << (ok?"ok":"broken") << std::endl;
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

		return lascheck(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
