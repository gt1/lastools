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
#include <libmaus2/util/Command.hpp>
#include <sys/wait.h>

#include <sys/types.h>
#include <pwd.h>

#include <slurm/slurm.h>

pid_t startCommand(libmaus2::util::Command const & C)
{
	pid_t const pid = fork();
	
	if ( pid == 0 )
	{
		try
		{
			int const r = C.dispatch();
			_exit(r);
		}
		catch(...)
		{
			_exit(EXIT_FAILURE);
		}
	}
	else if ( pid == static_cast<pid_t>(-1) )
	{
		int const error = errno;
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[V] fork failed: " << strerror(error) << std::endl;
		lme.finish();
		throw lme;
	}
	else
	{
		return pid;
	}
}

int slurmworker(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	
	std::string const hostname = arg[0];
	uint64_t const port = arg.getParsedRestArg<uint64_t>(1);
	
	char const * jobid_s = getenv("SLURM_JOB_ID");
	
	if ( ! jobid_s )
	{
		std::cerr << "[E] job id not found in SLURM_JOB_ID" << std::endl;
		return EXIT_FAILURE;
	}
	
	std::istringstream jobidistr(jobid_s);
	uint64_t jobid;
	jobidistr >> jobid;
	if ( ! jobidistr || jobidistr.peek() != std::istream::traits_type::eof() )
	{
		std::cerr << "[E] job id " << jobid_s << " found in SLURM_JOB_ID not parseable" << std::endl;
		return EXIT_FAILURE;
	}
	
	libmaus2::network::ClientSocket sockA(port,hostname.c_str());
	sockA.writeSingle<uint64_t>(jobid);
	
	uint64_t const workerid = sockA.readSingle<uint64_t>();
	
	bool running = true;
	
	enum state_type
	{
		state_idle,
		state_running
	};
	
	state_type state = state_idle;
	pid_t workpid = static_cast<pid_t>(-1);
	
	while ( running )
	{
		sockA.writeSingle<uint64_t>(2);
		sockA.readSingle<uint64_t>();
		
		switch ( state )
		{
			case state_idle:
			{
				// tell control we are idle
				sockA.writeSingle<uint64_t>(0);
				
				uint64_t const rep = sockA.readSingle<uint64_t>();
				
				// execute command
				if ( rep == 0 )
				{
					std::string const jobdesc = sockA.readString();
					std::istringstream jobdescistr(jobdesc);
					libmaus2::util::Command const com(jobdescistr);
					
					std::cerr << "[V] starting command " << com << std::endl;
					
					workpid = startCommand(com);
					state = state_running;
				}
				// stay idle
				else if ( rep == 1 )
				{
					sleep(1);
				}
				// terminate
				else if ( rep == 2 )
				{
					running = false;
				}
				
				break;
			}
			case state_running:
			{
				int status = 0;
				pid_t const wpid = waitpid(workpid, &status, WNOHANG);
				
				if ( wpid == 0 )
				{
					// still running
					::sleep(1);
				}
				else if ( wpid < 0 )
				{
					int const error = errno;
					
					switch ( error )
					{
						case EINTR:
						case EAGAIN:
							break;
						default:
						{
							std::cerr << "[E] waitpid failed: " << strerror(errno);
							return EXIT_FAILURE;
						}
					}
				}
				else
				{
					assert ( wpid == workpid );
					
					workpid = static_cast<pid_t>(-1);

					// tell control we finished a job
					sockA.writeSingle<uint64_t>(1);
					sockA.writeSingle<uint64_t>(status);
					
					std::cerr << "[V] finished with status " << status << std::endl;
					
					state = state_idle;
				}
			
				break;
			}
		}
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		
		int const r = slurmworker(arg);
		
		return r;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
