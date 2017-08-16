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
#include <libmaus2/util/ContainerDescriptionList.hpp>
#include <libmaus2/aio/PosixFdInputOutputStream.hpp>

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

struct ContainerInfo
{
	libmaus2::util::CommandContainer CN;
	std::string fn;

	ContainerInfo()
	{

	}
};


ContainerInfo handle(libmaus2::util::TempFileNameGenerator & tgen, std::vector<std::string> & lines, uint64_t const id, uint64_t & subid, uint64_t const cnid, std::vector<uint64_t> const & depid)
{
	libmaus2::util::CommandContainer CN;
	CN.id = cnid;
	CN.depid = depid;

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

			libmaus2::util::Command const C(in,out,err,code,tokens);
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

	// std::cout << container << "\t" << CN.V.size() << std::endl;

	++subid;
	lines.resize(0);

	ContainerInfo CI;
	CI.CN = CN;
	CI.fn = container;

	return CI;
}

std::string getcontextdir()
{
	char const * home = getenv("HOME");

	if ( ! home )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] unable to get ${HOME} directory" << std::endl;
		lme.finish();
		throw lme;
	}

	return std::string(home) + "/.commandpack";
}

void makecontextdir()
{
	std::string const command = std::string("mkdir -p ") + getcontextdir();

	int const r = system(command.c_str());

	if ( r != 0 )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] failed to run " << command << std::endl;
		lme.finish();
		throw lme;
	}
}

bool canlock(std::string const & fn)
{
	try
	{
		libmaus2::aio::PosixFdInputOutputStream PIS(fn,std::ios::in|std::ios::out);
		libmaus2::aio::PosixFdInputOutputStreamBuffer::LockObject const LO = PIS.lock();
		PIS.unlock(LO);
		return true;
	}
	catch(...)
	{
		return false;
	}
}

int commandpack(libmaus2::util::ArgParser const & arg)
{
	uint64_t linesperpack = arg.uniqueArgPresent("l") ? arg.getUnsignedNumericArg<uint64_t>("l") : getDefaultLinesPerPack();
	std::string const dn = arg.uniqueArgPresent("d") ? arg["d"] : getDefaultD(arg);
	std::vector < std::string > lines;
	uint64_t id = 0;
	uint64_t subid = 0;
	uint64_t cnid = 0;
	libmaus2::util::TempFileNameGenerator tgen(dn,4);

	std::vector < uint64_t > depid;
	std::vector < uint64_t > ndepid;
	std::vector < ContainerInfo > containers;

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
				{
					ndepid.push_back(cnid);
					containers.push_back(handle(tgen,lines,id,subid,cnid++,depid));
				}

				depid = ndepid;
				ndepid.resize(0);
				id++;
			}
			else
			{
				lines.push_back(line);

				if ( lines.size() >= linesperpack )
				{
					ndepid.push_back(cnid);
					containers.push_back(handle(tgen,lines,id,subid,cnid++,depid));
				}
			}
		}
	}

	if ( lines.size() )
	{
		ndepid.push_back(cnid);
		containers.push_back(handle(tgen,lines,id,subid,cnid++,depid));
	}
	
	for ( uint64_t i = 0; i < containers.size(); ++i )
	{
		libmaus2::util::CommandContainer & CC = containers[i].CN;
		std::vector<uint64_t> const & depid = CC.depid;
		
		for ( uint64_t j = 0; j < depid.size(); ++j )
			containers [ depid[j] ] . CN . rdepid.push_back(i);
	}

	for ( uint64_t i = 0; i < containers.size(); ++i )
	{
		std::cerr << containers[i].CN << std::endl;
		libmaus2::aio::OutputStreamInstance OSI(containers[i].fn);
		containers[i].CN.serialise(OSI);
		OSI.flush();
	}
	
	libmaus2::util::ContainerDescriptionList CDL;
	for ( uint64_t i = 0; i < containers.size(); ++i )
	{
		ContainerInfo const & CI = containers[i];

		std::string const & fn = CI.fn;
		std::vector<uint64_t> const & dep = CI.CN.depid;

		libmaus2::util::ContainerDescription const CD(fn, false /* started */, dep.size());
		CDL.V.push_back(CD);
	}

	std::string outfn;

	{
		std::ostringstream ostrfn;
		ostrfn << tgen.getFileName() << "_CDL";

		libmaus2::aio::OutputStreamInstance OSI(ostrfn.str());
		CDL.serialise(OSI);

		if ( canlock(ostrfn.str()) )
			outfn = ostrfn.str();
		else
		{
			std::cerr << "[W] cannot lock " << ostrfn.str() << std::endl;
		}
	}

	if ( ! outfn.size() )
	{
		makecontextdir();
		std::string const contextdir = getcontextdir();
		std::string const fn = contextdir + "/" + libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname) + "_CDL";

		libmaus2::aio::OutputStreamInstance OSI(fn);
		CDL.serialise(OSI);

		if ( canlock(fn) )
			outfn = fn;
		else
		{
			std::cerr << "[W] cannot lock " << fn << std::endl;
		}

	}

	std::cout << outfn << std::endl;

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
