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
#include <libmaus2/parallel/NumCpus.hpp>

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


ContainerInfo handle(libmaus2::util::TempFileNameGenerator & tgen, std::vector<std::string> & lines, uint64_t const id, uint64_t const subid, uint64_t const cnid, std::vector<uint64_t> const & depid, uint64_t const numthreads)
{
	libmaus2::util::CommandContainer CN;
	CN.id = cnid;
	CN.depid = depid;

	uint64_t o = 0;
	for ( uint64_t i = 0; i < lines.size(); ++i )
	{
		std::string line = lines[i];

		while ( line.size() && isspace(line[0]) )
			line = line.substr(1);

		if ( line.size() )
			lines[o++] = line;
	}
	lines.resize(o);
	
	CN.V.resize(lines.size());
	
	int volatile tfailed = 0;
	libmaus2::parallel::PosixSpinLock tfailedlock;
	libmaus2::parallel::PosixSpinLock tgenlock;

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
	#endif
	for ( uint64_t i = 0; i < lines.size(); ++i )
	{
		try
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
				
				{
					libmaus2::parallel::ScopePosixSpinLock slock(tgenlock);
					ostr << tgen.getFileName() << "_" << id << "_" << subid;
				}

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
				CN.V[i] = C;
			}
		}
		catch(std::exception const & ex)
		{
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << ex.what() << std::endl;
			}
		
			tfailedlock.lock();
			tfailed = 1;
			tfailedlock.unlock();
		}
	}
	
	if ( tfailed )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] handle failed" << std::endl;
		lme.finish();
		throw lme;
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

static std::vector < std::pair< std::string, std::vector<std::string> > > parseBatches(std::istream & ISI)
{
	std::vector < std::pair< std::string, std::vector<std::string> > > Vbatch;
	
	while ( ISI )
	{
		std::string line;
		std::getline(ISI,line);
		if ( ISI )
		{
			if ( line.size() && line.at(0) == '#' )
			{
				line = line.substr(1);
				while ( line.size() && isspace(line[0]) )
					line = line.substr(1);
				Vbatch.push_back(
					std::pair< std::string, std::vector<std::string> >(
						line,
						std::vector<std::string>()
					)
				);
			}
			else
			{
				if ( Vbatch.size() == 0 )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] first line outside batch: " << line << std::endl;
					lme.finish();
					throw lme;
				}

				Vbatch.back().second.push_back(line);
			}
		}
	}
	
	return Vbatch;
}


int commandpack(libmaus2::util::ArgParser const & arg)
{
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	uint64_t linesperpack = arg.uniqueArgPresent("l") ? arg.getUnsignedNumericArg<uint64_t>("l") : getDefaultLinesPerPack();
	std::string const dn = arg.uniqueArgPresent("d") ? arg["d"] : getDefaultD(arg);
	std::vector < std::string > lines;
	uint64_t cnid = 0;
	libmaus2::util::TempFileNameGenerator tgen(dn,4,16 /* dirmod */, 16 /* filemod */);

	std::vector < uint64_t > depid;
	std::vector < ContainerInfo > containers;

	std::vector < std::pair< std::string, std::vector<std::string> > > Vbatch = parseBatches(std::cin);

	for ( uint64_t id = 0; id < Vbatch.size(); ++id )
	{
		std::vector < uint64_t > ndepid;
		std::vector<std::string> & V = Vbatch[id].second;
		
		if ( ! V.size() )
		{
			std::string const com = "echo OK";
			V.push_back(com);
		}
		
		assert ( V.size() );
		
		uint64_t const packs = (V.size() + linesperpack - 1)/linesperpack;
		
		std::cerr << "[V] processing id=" << id << " " << Vbatch[id].first << " of size " << V.size() << " with " << packs << " packages" << std::endl;
		
		for ( uint64_t subid = 0; subid < packs; ++subid )
		{
			uint64_t const low = subid * linesperpack;
			uint64_t const high = std::min(low + linesperpack, static_cast<uint64_t>(V.size()));
			
			std::cerr << "\t[V] handling [" << low << "," << high << ") for cnid=" << cnid << std::endl;
			
			std::vector < std::string > lines(V.begin() + low, V.begin()+high);
			ndepid.push_back(cnid);
			containers.push_back(handle(tgen,lines,id,subid,cnid++,depid,numthreads));
		}
		
		depid = ndepid;
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
