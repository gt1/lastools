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

#include <libmaus2/util/CommandContainer.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/TempFileNameGenerator.hpp>
#include <libmaus2/aio/OutputStreamInstance.hpp>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <in.sh" << std::endl;

	return ostr.str();
}

static uint64_t getDefaultLinesPerPack()
{
	return 1000;
}

static std::string getDefaultD(libmaus2::util::ArgParser const & arg)
{
	return libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
}

std::string which(std::string const prog)
{
	if ( prog.size() && prog.find('/') != std::string::npos )
		return prog;

	std::string const path = getenv("PATH");
	std::deque<std::string> const Vpath = libmaus2::util::stringFunctions::tokenize(path,std::string(":"));

	for ( uint64_t i = 0; i < Vpath.size(); ++i )
	{
		std::string const p = Vpath[i] + '/' + prog;
		struct stat sb;
		int const r = stat(p.c_str(),&sb);

		if ( r == 0 )
		{
			if (
				(S_ISREG(sb.st_mode) || S_ISLNK(sb.st_mode))
				&&
				(access(p.c_str(),X_OK) == 0)
			)
			{
				return p;
			}
		}
	}

	libmaus2::exception::LibMausException lme;
	lme.getStream() << "[E] " << prog << " not found" << std::endl;
	lme.finish();
	throw lme;
}

void handle(libmaus2::util::TempFileNameGenerator & tgen, std::vector<std::string> & lines, uint64_t const id, uint64_t & subid)
{
	libmaus2::util::CommandContainer CN;

	for ( uint64_t i = 0; i < lines.size(); ++i )
	{
		std::string line = lines[i];

		while ( line.size() && isspace(line[0]) )
			line = line.substr(1);

		if ( ! line.size() )
			continue;

		uint64_t h = 0;
		while ( h < line.size() && !isspace(line[h]) )
			++h;

		std::string const fcom = which(line.substr(0,h));

		while ( h < line.size() && isspace(line[h]) )
			++h;

		line = fcom + " " + line.substr(h);

		if ( line.size() )
		{
			std::string const in = "/dev/null";

			std::ostringstream ostr;
			ostr << tgen.getFileName() << "_" << id << "_" << subid;

			std::string const filebase = ostr.str();
			std::string const out = filebase + ".out";
			std::string const err = filebase + ".err";
			std::string const script = filebase + ".script";
			std::string const code = filebase + ".returncode";
			std::string const command = filebase + ".com";

			{
				libmaus2::aio::OutputStreamInstance OSI(script);
				OSI << "#! /bin/bash\n";
				OSI << line << "\n";
				OSI << "RT=$?\n";
				OSI << "echo ${RT} >" << code << "\n";
				OSI << "exit ${RT}\n";
				OSI.flush();
			}

			{
				libmaus2::aio::OutputStreamInstance OSI(command);
				OSI << line << "\n";
				OSI.flush();
			}

			std::vector<std::string> tokens;

			tokens.push_back("/bin/bash");
			tokens.push_back(script);

			libmaus2::util::Command const C(in,out,err,tokens);
			CN.V.push_back(C);
		}
	}

	std::ostringstream ostr;
	ostr << tgen.getFileName() << "_" << id << "_" << subid;
	std::string const filebase = ostr.str();
	std::string const container = filebase + ".container";

	{
		libmaus2::aio::OutputStreamInstance OSI(container);
		CN.serialise(OSI);
		OSI.flush();
	}

	std::cout << container << "\t" << CN.V.size() << std::endl;

	++subid;
	lines.resize(0);
}

int commandpack(libmaus2::util::ArgParser const & arg)
{
	uint64_t linesperpack = arg.uniqueArgPresent("l") ? arg.getUnsignedNumericArg<uint64_t>("l") : getDefaultLinesPerPack();
	std::string const dn = arg.uniqueArgPresent("d") ? arg["d"] : getDefaultD(arg);
	std::vector < std::string > lines;
	uint64_t id = 0;
	uint64_t subid = 0;
	libmaus2::util::TempFileNameGenerator tgen(dn,4);

	while ( std::cin )
	{
		std::string line;
		std::getline(std::cin,line);

		if ( line.size() )
		{
			// force split at #
			if ( line.size() && line[0] == '#' )
			{
				if ( lines.size() )
					handle(tgen,lines,id,subid);

				id++;
			}
			else
			{
				lines.push_back(line);

				if ( lines.size() >= linesperpack )
				{
					handle(tgen,lines,id,subid);
				}
			}
		}
	}

	if ( lines.size() )
		handle(tgen,lines,id,subid);

	return EXIT_SUCCESS;
}

/**
 **/
int main(int argc, char *argv[])
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
		#if 0
		else if ( arg.size() < 0 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}
		#endif
		else
		{
			return commandpack(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
