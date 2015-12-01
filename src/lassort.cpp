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
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [-t<threads> -T<tmpprefix> -f<mergefanin> -s<sortorder>] <out.las> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -t : number of threads (defaults to number of cores on machine)\n";
	ostr << " -T : prefix for temporary files (default: create files in current working directory)\n";
	ostr << " -f : merge fan in (default: 64)\n";
	ostr << " -s : sort order (canonical or ba, default: canonical)\n";

	return ostr.str();
}

int lassort(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : arginfo.getDefaultTmpFileName();
	uint64_t const mergefanin = arg.uniqueArgPresent("f") ? arg.getUnsignedNumericArg<uint64_t>("f") : 64;
	std::string const sortorder = arg.uniqueArgPresent("s") ? arg["s"] : std::string("canonical");
	assert ( mergefanin );

	std::string const outfilename = arg[0];
	std::vector<std::string> infilenames;
	for ( uint64_t i = 1; i < arg.size(); ++i )
		infilenames.push_back(arg[i]);

	if ( sortorder == "ba" )
		libmaus2::dazzler::align::SortingOverlapOutputBuffer<libmaus2::dazzler::align::OverlapComparatorBReadARead>::sortAndMerge(infilenames,outfilename,tmpfilebase,mergefanin,1 /* numthreads */);
	else
		libmaus2::dazzler::align::SortingOverlapOutputBuffer<>::sortAndMerge(infilenames,outfilename,tmpfilebase,mergefanin,numthreads);

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
