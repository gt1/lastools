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

	ostr << "usage: " << arg.progname << " [-t<numthreads>] [--tspace<value>] <out.las> <in.las>" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int64_t getDefaultTSpace()
{
	return libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
}

int lasproject(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
{
	std::string const outlas = arg[0];
	std::string const inlas = arg[1];

	if ( ! arg.uniqueArgPresent("a") )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] required argument -a not present" << std::endl;
		lme.finish();
		throw lme;
	}

	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(inlas);
	int64_t const aid = arg.uniqueArgPresent("a") ? arg.getUnsignedNumericArg<int64_t>("a") : -1;
	int64_t const bid = arg.uniqueArgPresent("b") ? arg.getUnsignedNumericArg<int64_t>("b") : -1;
	libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(inlas));
	libmaus2::dazzler::align::AlignmentWriter AW(outlas,tspace);

	libmaus2::dazzler::align::Overlap OVL;

	while ( PIN->getNextOverlap(OVL) )
		if ( OVL.aread == aid && (bid < 0 || OVL.bread == bid) )
			AW.put(OVL);

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

		return lasproject(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
