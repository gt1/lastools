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

#include <config.h>

#include <runProgram.hpp>
#include <which.hpp>

#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/TempFileNameGenerator.hpp>
#include <libmaus2/aio/OutputStreamInstance.hpp>
#include <libmaus2/aio/OutputStreamFactoryContainer.hpp>
#include <libmaus2/util/WriteableString.hpp>
#include <libmaus2/network/Socket.hpp>
#include <libmaus2/util/GetFileSize.hpp>
#include <libmaus2/util/ContainerDescriptionList.hpp>
#include <libmaus2/util/CommandContainer.hpp>

struct CommandContainerView
{
	std::string const cdl;
	libmaus2::util::ContainerDescriptionList CDL;
	std::vector < libmaus2::util::ContainerDescription > & CDLV;
	std::vector < libmaus2::util::CommandContainer > VCC;

	static libmaus2::util::ContainerDescriptionList loadCDL(std::string const & cdl)
	{
		libmaus2::util::ContainerDescriptionList CDL;
		libmaus2::aio::InputStreamInstance ISI(cdl);
		CDL.deserialise(ISI);
		return CDL;
	}

	static std::vector < libmaus2::util::CommandContainer > loadVCC(std::vector < libmaus2::util::ContainerDescription > & CDLV)
	{
		std::vector < libmaus2::util::CommandContainer > VCC(CDLV.size());
		for ( uint64_t i = 0; i < CDLV.size(); ++i )
		{
			libmaus2::aio::InputStreamInstance ISI(CDLV[i].fn);
			VCC[i].deserialise(ISI);
			CDLV[i].missingdep = 0;
		}

		return VCC;
	}

	CommandContainerView(std::string const & rcdl)
	: 
	  cdl(rcdl),
	  CDL(loadCDL(cdl)),
	  CDLV(CDL.V),
	  VCC(loadVCC(CDLV))
	{
	}
};

std::ostream & operator<<(std::ostream & out, CommandContainerView const & C)
{
	for ( uint64_t i = 0; i < C.VCC.size(); ++i )
	{
		out << "CommandContainer[" << i << "]=" << C.VCC[i] << std::endl;
	}
	
	return out;
}

int showcdl(libmaus2::util::ArgParser const & arg)
{
	std::string const cdl = arg[0];

	CommandContainerView CCV(cdl);
	
	std::cout << CCV;

	return EXIT_SUCCESS;
}

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;
	ostr << "usage: " << arg.progname << " <cdl>";
	return ostr.str();
}

int main(int argc, char * argv[])
{
	try
	{
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
		else if ( arg.size() < 1 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		int const r = showcdl(arg);

		return r;
	}
	catch(std::exception const & ex)
	{
		std::cerr << "[E] exception in main: " << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
