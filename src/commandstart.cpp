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
#include <libmaus2/util/WriteableString.hpp>
#include <libmaus2/network/Socket.hpp>
#include <sys/wait.h>

#if defined(__APPLE__)
#include <crt_externs.h>
#endif

char ** getEnviron()
{
	#if defined(__APPLE__)
	return *_NSGetEnviron();
	#else
	return environ;
	#endif
}


std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " container" << std::endl;

	return ostr.str();
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

std::string runProgram(std::vector<std::string> const & args, libmaus2::util::ArgParser const & arg)
{
	std::string const tmpprefix = std::string("/tmp/") + libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	libmaus2::autoarray::AutoArray < libmaus2::util::WriteableString::unique_ptr_type > AW(args.size());
	for ( uint64_t i = 0; i < args.size(); ++i )
	{
		libmaus2::util::WriteableString::unique_ptr_type tptr(new libmaus2::util::WriteableString(args[i]));
		AW[i] = UNIQUE_PTR_MOVE(tptr);
	}
	libmaus2::autoarray::AutoArray < char * > AA(args.size()+1);
	for ( uint64_t i = 0; i < args.size(); ++i )
		AA[i] = AW[i]->A.begin();
	AA[args.size()] = 0;

	std::string const prog = which(args[0]);

	std::string const errfile = tmpprefix + "_runprog.err";
	libmaus2::util::TempFileRemovalContainer::addTempFile(errfile);
	std::string const outfile = tmpprefix + "_runprog.out";
	libmaus2::util::TempFileRemovalContainer::addTempFile(outfile);

	pid_t const pid = fork();

	if ( pid == static_cast<pid_t>(-1) )
	{
		int const error = errno;
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] fork failed: " << strerror(error) << std::endl;
		lme.finish();
		throw lme;
	}
	if ( pid == 0 )
	{
		int const fdout = ::open(outfile.c_str(),O_CREAT|O_TRUNC|O_RDWR,0600);
		if ( fdout == -1 )
		{
			int const error = errno;
			std::cerr << "[E] failed to open " << outfile << ":" << strerror(error) << std::endl;
			_exit(EXIT_FAILURE);
		}
		int const fderr = ::open(errfile.c_str(),O_CREAT|O_TRUNC|O_RDWR,0600);
		if ( fderr == -1 )
		{
			int const error = errno;
			std::cerr << "[E] failed to open " << errfile << ":" << strerror(error) << std::endl;
			_exit(EXIT_FAILURE);
		}

		if ( ::close(STDOUT_FILENO) == -1 )
			_exit(EXIT_FAILURE);
		if ( ::close(STDERR_FILENO) == -1 )
			_exit(EXIT_FAILURE);

		if ( dup2 ( fdout, STDOUT_FILENO ) == -1 )
			_exit(EXIT_FAILURE);
		if ( dup2 ( fderr, STDERR_FILENO ) == -1 )
			_exit(EXIT_FAILURE);

		execve(prog.c_str(),AA.begin(),getEnviron());
		_exit(EXIT_FAILURE);
	}

	int status = 0;
	int const r = waitpid(pid,&status,0/*options */);
	if ( r == -1 )
	{
		int const error = errno;
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] waitpid failed: " << strerror(error) << std::endl;
		lme.finish();
		throw lme;
	}

	if ( WIFEXITED(status) && (WEXITSTATUS(status) == 0) )
	{
		std::cerr << "[V] " << prog << " finished ok" << "\n";
	}
	else
	{
		libmaus2::exception::LibMausException lme;
		std::ostream & errstr = lme.getStream();

		errstr << "[E] " << args[0] << " failed" << std::endl;

		if ( WIFEXITED(status) )
		{
			errstr << "[E] exit code " << WEXITSTATUS(status) << std::endl;
		}

		if ( libmaus2::util::GetFileSize::fileExists(errfile) )
		{
			libmaus2::aio::InputStreamInstance ISI(errfile);
			while ( ISI )
			{
				std::string line;
				std::getline(ISI,line);
				if ( ISI )
					errstr << line << std::endl;
			}
		}

		lme.finish();

		throw lme;
	}

	if ( !libmaus2::util::GetFileSize::fileExists(outfile) )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] output file " << outfile << "  does not exist." << std::endl;
		lme.finish();
		throw lme;
	}

	libmaus2::aio::InputStreamInstance ISI(outfile);
	std::istreambuf_iterator<char> ia(ISI);
	std::istreambuf_iterator<char> ib;
	std::string const res(ia,ib);

	libmaus2::aio::FileRemoval::removeFile(errfile);
	libmaus2::aio::FileRemoval::removeFile(outfile);

	return res;
}


int commandstart(libmaus2::util::ArgParser const & arg)
{
	std::string const hostname = arg[0];
	uint64_t const port = arg.getParsedRestArg<uint64_t>(1);

	std::string cocommand = which("cocommand");
	std::string commandfollowup = which("commandfollowup");
	std::string commandstart = which("commandstart");

	libmaus2::network::ClientSocket CS(port,hostname.c_str());

	std::string data = CS.readString();
	std::istringstream PIS(data);

	libmaus2::util::ContainerDescriptionList CDL(PIS);

	std::string const tmpprefix = std::string("/tmp/") + libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	std::string const jobfn = tmpprefix + "_job.sbatch";
	libmaus2::util::TempFileRemovalContainer::addTempFile(jobfn);

	for ( uint64_t i = 0; i < CDL.V.size(); ++i )
		if ( !CDL.V[i].started && !CDL.V[i].missingdep )
		{
			std::string const fn = CDL.V[i].fn;
			libmaus2::aio::InputStreamInstance CCISI(fn);
			libmaus2::util::CommandContainer CC(CCISI);
			std::cerr << "could start " << i << " " << fn << std::endl;

			int64_t njobid = -1;

			{
				std::ostringstream ostr;
				ostr << "#!/bin/bash\n";
				ostr << "#\n";
				ostr << "#SBATCH --job-name=command_" << CC.id << "\n";
				ostr << "#SBATCH --output=/dev/null\n";
				ostr << "#SBATCH --error=/dev/null\n";
				ostr << "#SBATCH --array=0-" << CC.V.size()-1 << "\n";
				ostr << "#\n";
				ostr << "#SBATCH --ntasks=1\n";
				ostr << "#SBATCH --time=1440\n";
				ostr << "#SBATCH --mem=240000\n";
				ostr << "#SBATCH --cpus-per-task=24\n";

				#if 0
				if ( depids.size() )
				{
					ostr << "#SBATCH --dependency=afterok";
					for ( uint64_t i = 0; i < depids.size(); ++i )
						ostr << ":" << depids[i];
					ostr << "\n";
				}
				#endif

				ostr << "\n";
				ostr << "srun " << cocommand << " " << fn << "\n";

				std::cerr << ostr.str();

				{
					libmaus2::aio::OutputStreamInstance OSI(jobfn);
					OSI << ostr.str();
				}

				std::vector<std::string> sbatchargs;
				sbatchargs.push_back("sbatch");
				sbatchargs.push_back(jobfn);
				std::string out;
				bool outok = false;

				for ( uint64_t qq = 0; (!outok) && (qq < 100); ++qq )
				{
					try
					{
						std::cerr << "calling " << sbatchargs[0] << " " << sbatchargs[1] << " try " << qq << std::endl;
						out = runProgram(sbatchargs,arg);
						outok = true;
					}
					catch(std::exception const & ex)
					{
						std::cerr << ex.what() << std::endl;
					}
				}

				if ( outok )
				{
					std::istringstream sistr(out);
					while ( sistr )
					{
						std::string line;
						std::getline(sistr,line);

						if (
							sistr
							&&
							line.find("Submitted") == 0
						)
						{
							std::deque<std::string> Vtokens = libmaus2::util::stringFunctions::tokenize(line,std::string(" "));
							if ( Vtokens.size() )
							{
								std::istringstream uistr(Vtokens.back());
								uint64_t jobid;
								uistr >> jobid;
								njobid = jobid;

								std::cerr << "[V] submitted job " << jobid << " for " << sbatchargs[1] << "\n";
							}
						}
					}

					if ( njobid < 0 )
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] failed to run " << sbatchargs[1] << " (job id not found)" << "\n";
						lme.finish();
						throw lme;
					}
				}
				else
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] failed to run " << sbatchargs[1] << " (job id not found)" << "\n";
					lme.finish();
					throw lme;
				}

				libmaus2::aio::FileRemoval::removeFile(jobfn);
			}

			int64_t nnjobid = -1;
			{
				std::ostringstream ostr;
				ostr << "#!/bin/bash\n";
				ostr << "#\n";
				ostr << "#SBATCH --job-name=followup_" << CC.id << "\n";
				// ostr << "#SBATCH --output=/dev/null\n";
				ostr << "#SBATCH --output=followup_" << CC.id << ".out\n";
				//ostr << "#SBATCH --error=/dev/null\n";
				ostr << "#\n";
				ostr << "#SBATCH --ntasks=1\n";
				ostr << "#SBATCH --time=1440\n";
				ostr << "#SBATCH --mem=4000\n";
				ostr << "#SBATCH --cpus-per-task=1\n";
				ostr << "#SBATCH --dependency=afterany:" << njobid << "\n";

				ostr << "\n";
				ostr << "srun " << commandfollowup << " " << i << " " << hostname << " " << port << "\n";

				std::cerr << ostr.str();

				{
					libmaus2::aio::OutputStreamInstance OSI(jobfn);
					OSI << ostr.str();
				}

				std::vector<std::string> sbatchargs;
				sbatchargs.push_back("sbatch");
				sbatchargs.push_back(jobfn);
				std::string out;
				bool outok = false;

				for ( uint64_t qq = 0; (!outok) && (qq < 100); ++qq )
				{
					try
					{
						std::cerr << "calling " << sbatchargs[0] << " " << sbatchargs[1] << " try " << qq << std::endl;
						out = runProgram(sbatchargs,arg);
						outok = true;
					}
					catch(std::exception const & ex)
					{
						std::cerr << ex.what() << std::endl;
					}
				}

				if ( outok )
				{
					std::istringstream sistr(out);
					while ( sistr )
					{
						std::string line;
						std::getline(sistr,line);

						if (
							sistr
							&&
							line.find("Submitted") == 0
						)
						{
							std::deque<std::string> Vtokens = libmaus2::util::stringFunctions::tokenize(line,std::string(" "));
							if ( Vtokens.size() )
							{
								std::istringstream uistr(Vtokens.back());
								uint64_t jobid;
								uistr >> jobid;
								nnjobid = jobid;

								std::cerr << "[V] submitted follow up job " << jobid << " for " << sbatchargs[1] << "\n";
							}
						}
					}

					if ( nnjobid < 0 )
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] failed to run " << sbatchargs[1] << " (job id not found)" << "\n";
						lme.finish();
						throw lme;
					}
				}
				else
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] failed to run " << sbatchargs[1] << " (job id not found)" << "\n";
					lme.finish();
					throw lme;
				}

				libmaus2::aio::FileRemoval::removeFile(jobfn);
			}

			CDL.V[i].started = true;
		}

	std::ostringstream OPIS;
	CDL.serialise(OPIS);

	CS.writeString(OPIS.str());
	CS.writeSingle<uint64_t>(0); // keep running

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
		else if ( arg.size() < 1 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}
		else
		{
			return commandstart(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
