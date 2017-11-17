/*
    lastools
    Copyright (C) 2017 German Tischler

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
#include <libmaus2/dazzler/align/AlignmentFileSymmetryCheck.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/aio/SerialisedPeeker.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <in.db> <in1.las> ..." << std::endl;

	return ostr.str();
}


int laschecksymmetry(libmaus2::util::ArgParser & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	std::string const dbfn = arg[0];
	std::vector<std::string> Vin;
	for ( uint64_t i = 1; i < arg.size(); ++i )
		Vin.push_back(arg[i]);

	int const r = libmaus2::dazzler::align::AlignmentFileSymmetryCheck::checkSymmetry(dbfn,Vin,tmpfilebase,&std::cerr);

	if ( r == EXIT_SUCCESS )
		std::cout << "ok" << std::endl;
	else
		std::cout << "not ok" << std::endl;

	return r;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);

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

		return laschecksymmetry(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
