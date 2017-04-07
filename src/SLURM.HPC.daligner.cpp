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
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/GetFileSize.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/stringFunctions.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/util/WriteableString.hpp>
#include <libmaus2/util/PosixExecute.hpp>
#include <libmaus2/util/TempFileNameGenerator.hpp>
#include <cstring>

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

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
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

std::string fullpath(std::string prog)
{
	uint64_t i = 0;
	while ( i < prog.size() && !isspace(prog[i]) )
		++i;

	std::string command = prog.substr(0,i);
	std::string args = i < prog.size() ? prog.substr(i+1) : std::string();

	command = which(command);

	return command + " " + args;
}

#if 0
struct RunInfo
{
	RunInfo()
	{

	}

	RunInfo(std::istream & in)
	{

	}
};
#endif

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		std::string const tmpname = libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
		std::vector<std::string> Vsbatch;

		bool const havecoverage = arg.uniqueArgPresent("c");
		int64_t const coverage = havecoverage ? arg.getParsedArg<unsigned int>("c") : -1;

		std::string DASedit;
		try
		{
			DASedit = which("DASedit");
		}
		catch(std::exception const &)
		{

		}

		if ( DASedit.size() )
			std::cerr << "[V] found DASedit at " << DASedit << std::endl;

		std::string DASpatch;
		try
		{
			DASpatch = which("DASpatch");
		}
		catch(std::exception const &)
		{

		}

		if ( DASpatch.size() )
			std::cerr << "[V] found DASpatch at " << DASpatch << std::endl;

		std::string DAStrim;
		try
		{
			DAStrim = which("DAStrim");
		}
		catch(std::exception const &)
		{

		}

		if ( DAStrim.size() )
			std::cerr << "[V] found DAStrim at " << DAStrim << std::endl;

		std::string dasqv;
		try
		{
			dasqv = which("DASqv");
		}
		catch(std::exception const &)
		{

		}

		if ( dasqv.size() )
			std::cerr << "[V] found DASqv at " << dasqv << std::endl;

		std::string cattrack;
		try
		{
			cattrack = which("Catrack");
		}
		catch(std::exception const &)
		{

		}

		if ( cattrack.size() )
			std::cerr << "[V] found Catrack at " << cattrack << std::endl;

		int64_t sepindex = argc;
		for ( int i = 1; i < argc; ++i )
			if ( std::string(argv[i]) == "--" )
			{
				sepindex = i;
				break;
			}
			else
			{
				Vsbatch.push_back(argv[i]);
			}

		if ( ! (sepindex + 1 < argc) )
		{
			std::cerr << "[E] usage: " << argv[0] << " <opts> -- <prog> ..." << std::endl;
			return EXIT_FAILURE;
		}

		argc -= (sepindex+1);
		argv += (sepindex+1);

		std::vector<std::string> args;
		for ( int i = 0; i < argc; ++i )
			args.push_back(argv[i]);

		libmaus2::util::ArgParser const argdal(argc,argv);
		uint64_t numthreads;
		if ( argdal.uniqueArgPresent("T") )
		{
			numthreads = argdal.getParsedArg<unsigned int>("T");
		}
		else
		{
			numthreads = 4;
			std::ostringstream ostr;
			ostr << "-T" << numthreads;
			args.push_back(ostr.str());
		}

		uint64_t mem;
		if ( argdal.uniqueArgPresent("M") )
		{
			mem = argdal.getParsedArg<unsigned int>("M");
		}
		else
		{
			mem = 30;
			std::ostringstream ostr;
			ostr << "-M" << mem;
			args.push_back(ostr.str());
		}

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

		std::string const HPC_daligner = which(argv[0]);
		std::string const tmpfileprefix = getTmpFileBase(arg);
		libmaus2::util::TempFileNameGenerator tmpgen(tmpfileprefix + "_subdir",5);

		std::string const errfile = tmpfileprefix + "_HPC.err";
		libmaus2::util::TempFileRemovalContainer::addTempFile(errfile);
		std::string const outfile = tmpfileprefix + "_HPC.out";
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
				return EXIT_FAILURE;
			}
			int const fderr = ::open(errfile.c_str(),O_CREAT|O_TRUNC|O_RDWR,0600);
			if ( fderr == -1 )
			{
				int const error = errno;
				std::cerr << "[E] failed to open " << errfile << ":" << strerror(error) << std::endl;
				return EXIT_FAILURE;
			}

			if ( ::close(STDOUT_FILENO) == -1 )
				return EXIT_FAILURE;
			if ( ::close(STDERR_FILENO) == -1 )
				return EXIT_FAILURE;

			if ( dup2 ( fdout, STDOUT_FILENO ) == -1 )
				return EXIT_FAILURE;
			if ( dup2 ( fderr, STDERR_FILENO ) == -1 )
				return EXIT_FAILURE;

			execve(HPC_daligner.c_str(),AA.begin(),getEnviron());
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
			std::cerr << "[V] " << HPC_daligner << " finished ok" << std::endl;

			if ( !libmaus2::util::GetFileSize::fileExists(outfile) )
			{
				std::cerr << "[E] output file does not exist." << std::endl;
				return EXIT_FAILURE;
			}

			libmaus2::aio::InputStreamInstance::unique_ptr_type ISI(new libmaus2::aio::InputStreamInstance(outfile));
			std::vector < std::pair< std::string, std::vector<std::string> > > Vbatch;
			while ( *ISI )
			{
				std::string line;
				std::getline(*ISI,line);
				if ( *ISI )
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
							std::cerr << "[E] first line outside batch: " << line << std::endl;
							return EXIT_FAILURE;
						}

						Vbatch.back().second.push_back(line);
					}
				}
			}
			ISI.reset();

			std::ostringstream patchfnstr;
			patchfnstr << tmpgen.getFileName() << "_ndb";
			std::string const patchfn = patchfnstr.str();

			int64_t lastmerge = -1;
			for ( uint64_t i = 0; i < Vbatch.size(); ++i )
				if ( Vbatch[i].second.size() && Vbatch[i].second[0].find("LAmerge") == 0 )
					lastmerge = i;
			std::vector<std::string> Voutputfiles;
			if ( lastmerge != -1 )
			{
				for ( uint64_t i = 0; i < Vbatch[lastmerge].second.size(); ++i )
				{
					std::string s = Vbatch[lastmerge].second[i];

					if ( s.find("LAmerge ") == 0 )
					{
						s = s.substr(strlen("LAmerge "));

						if ( s.find(" ") )
						{
							s = s.substr(0,s.find(" "));
							std::cerr << "Last merge " << lastmerge << " " << s << std::endl;
							Voutputfiles.push_back(s);
						}
					}
				}
			}
			if ( dasqv.size() && havecoverage )
			{
				Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("dasqv", std::vector<std::string>()));

				for ( uint64_t i = 0; i < Voutputfiles.size(); ++i )
				{
					std::ostringstream ostr;
					ostr << dasqv << " -c" << coverage << " " << argdal[0] << " " << Voutputfiles[i];
					Vbatch.back().second.push_back(ostr.str());
				}

				if ( cattrack.size() )
				{
					Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("Catrack", std::vector<std::string>()));
					std::ostringstream ostr;
					ostr << cattrack << " " << argdal[0] << " " << "qual";
					Vbatch.back().second.push_back(ostr.str());

					if ( DAStrim.size() )
					{
						Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DAStrim", std::vector<std::string>()));

						for ( uint64_t i = 0; i < Voutputfiles.size(); ++i )
						{
							std::ostringstream ostr;
							ostr << DAStrim << " -v -g20 -b26" << " " << argdal[0] << " " << Voutputfiles[i];
							Vbatch.back().second.push_back(ostr.str());
						}

						if ( cattrack.size() )
						{
							Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("Catrack", std::vector<std::string>()));
							std::ostringstream ostr;
							ostr << cattrack << " " << argdal[0] << " " << "trim";
							Vbatch.back().second.push_back(ostr.str());
						}
					}

					if ( DASpatch.size() && cattrack.size() )
					{
						Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DASpatch", std::vector<std::string>()));

						for ( uint64_t i = 0; i < Voutputfiles.size(); ++i )
						{
							std::ostringstream ostr;
							ostr << DASpatch << " -v " << " " << argdal[0] << " " << Voutputfiles[i];
							Vbatch.back().second.push_back(ostr.str());
						}

						if ( cattrack.size() )
						{
							Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("Catrack", std::vector<std::string>()));
							std::ostringstream ostr;
							ostr << cattrack << " " << argdal[0] << " " << "patch";
							Vbatch.back().second.push_back(ostr.str());
						}
					}

					if ( DASedit.size() )
					{
						Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DASedit", std::vector<std::string>()));

						std::ostringstream ostr;
						ostr << DASedit << " -v -x2000 " << " " << argdal[0] << " " << patchfn;
						Vbatch.back().second.push_back(ostr.str());
					}
				}
			}

			std::vector< std::vector<std::string> > Vlevellogs;
			std::vector< std::vector<std::string> > Vlevelsublogs;

			std::vector<std::string> Vlevellogcat;
			std::vector<std::string> Vlevelloglistfiles;
			std::vector<std::string> Vlevelsublogcat;
			std::vector<std::string> Vlevelsubloglistfiles;

			std::ostringstream commandlistfnstr;
			commandlistfnstr << tmpgen.getFileName() << ".commandlist";
			std::string const commandlistfn = commandlistfnstr.str();

			libmaus2::aio::OutputStreamInstance::unique_ptr_type comOSI(new libmaus2::aio::OutputStreamInstance(commandlistfn));

			std::vector<uint64_t> depids;
			std::vector<uint64_t> allids;
			for ( uint64_t i = 0; i < Vbatch.size(); ++i )
			{
				Vlevellogs.push_back(std::vector<std::string>(0));
				Vlevelsublogs.push_back(std::vector<std::string>(0));

				std::ostringstream catlogfnstr;
				catlogfnstr << tmpgen.getFileName() << "_" << i << ".log";
				std::ostringstream listlogfnstr;
				listlogfnstr << tmpgen.getFileName() << "_" << i << ".list";

				Vlevellogcat.push_back(catlogfnstr.str());
				Vlevelloglistfiles.push_back(listlogfnstr.str());

				std::ostringstream catsublogfnstr;
				catsublogfnstr << tmpgen.getFileName() << "_" << i << ".sublog";
				std::ostringstream listsublogfnstr;
				listsublogfnstr << tmpgen.getFileName() << "_" << i << ".sublist";

				Vlevelsublogcat.push_back(catsublogfnstr.str());
				Vlevelsubloglistfiles.push_back(listsublogfnstr.str());

				std::cerr << Vbatch[i].first << std::endl;
				std::vector<uint64_t> ndepids;
				for ( uint64_t j = 0; j < Vbatch[i].second.size(); ++j )
				{
					std::string command = Vbatch[i].second[j];

					std::vector<std::string> Largs;
					std::string Targ;
					uint64_t l = 0;
					while ( l < command.size() )
					{
						while ( l < command.size() && isspace(command[l]) )
							++l;

						uint64_t h = l;
						while ( h < command.size() && !isspace(command[h]) )
							++h;

						if ( h != l )
						{
							std::string const arg = command.substr(l,h-l);
							Largs.push_back(arg);
							l = h;

							if ( arg.size() >= 2 && arg.substr(0,2) == "-T" )
								Targ = arg.substr(2);
						}
					}

					uint64_t lnumthreads = 1;
					if ( Targ.size() )
					{
						std::istringstream istr(Targ);
						istr >> lnumthreads;
					}
					else
					{
						if ( Largs.at(0) == "daligner" )
						{
							lnumthreads = 4;

							std::ostringstream snumthreads;
							snumthreads << lnumthreads;

							command.replace(0,strlen("daligner"),std::string("daligner -T")+snumthreads.str());
						}
						else
							lnumthreads = 1;
					}

					std::ostringstream commandfnstr;
					commandfnstr << tmpgen.getFileName() << "_" << i << "_" << j << ".command";
					std::string const commandfn = commandfnstr.str();

					std::ostringstream sublogfnstr;
					sublogfnstr << tmpgen.getFileName() << "_" << i << "_" << j << ".sublog";
					std::string const sublogfn = sublogfnstr.str();

					std::ostringstream logfnstr;
					logfnstr << tmpgen.getFileName() << "_" << i << "_" << j << ".log";
					std::ostringstream jobfnstr;
					jobfnstr << tmpgen.getFileName() << "_" << i << "_" << j << ".job";
					std::string const logfn = logfnstr.str();
					std::ostringstream jobnamestr;
					jobnamestr << tmpname << "_" << i << "_" << j;
					std::string const jobname = jobnamestr.str();

					command = fullpath(command);

					{
						libmaus2::aio::OutputStreamInstance OSI(commandfn);
						OSI << "#!/bin/bash" << std::endl;
						OSI << command << std::endl;
						// OSI << "rm -f " << commandfn < std::endl;
						OSI.flush();
					}

					Vlevellogs.back().push_back(logfn);
					Vlevelsublogs.back().push_back(sublogfn);

					std::ostringstream ostr;
					ostr << "#!/bin/bash\n";
					ostr << "#\n";
					ostr << "#SBATCH --job-name=" << jobname << "\n";
					ostr << "#SBATCH --output=" << logfn << "\n";
					ostr << "#\n";
					ostr << "#SBATCH --ntasks=1\n";
					ostr << "#SBATCH --time=24:00:00\n";
					ostr << "#SBATCH --mem=" << (mem*1024)+512 << "\n";
					ostr << "#SBATCH --cpus-per-task=" << lnumthreads << "\n";

					if ( depids.size() )
					{
						ostr << "#SBATCH --dependency=afterok";
						for ( uint64_t i = 0; i < depids.size(); ++i )
							ostr << ":" << depids[i];
						ostr << "\n";
					}

					ostr << "\n";
					ostr << "srun bash " << commandfn << " >" << sublogfn << " 2>&1" << "\n";

					(*comOSI) << commandfn << std::endl;

					libmaus2::aio::OutputStreamInstance::unique_ptr_type Pjob(new libmaus2::aio::OutputStreamInstance(jobfnstr.str()));
					(*Pjob) << ostr.str();
					Pjob.reset();

					std::cerr << ostr.str() << std::endl;

					std::ostringstream comstr;
					comstr << "sbatch " << jobfnstr.str();
					std::string const com = comstr.str();

					std::string out;
					std::string err;
					std::cerr << "[V] calling " << com << std::endl;
					int const r = libmaus2::util::PosixExecute::execute(com,out,err,false /* do not throw */);
					std::cerr << "[V] called " << com << std::endl;

					if ( r != EXIT_SUCCESS )
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] failed to run " << com << std::endl;
						lme.finish();
						throw lme;
					}

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
								ndepids.push_back(jobid);
								allids.push_back(jobid);

								std::cerr << "[V] submitted job " << jobid << " for " << com << std::endl;
							}
						}
					}

					libmaus2::aio::FileRemoval::removeFile(jobfnstr.str());
				}

				{
				libmaus2::aio::OutputStreamInstance liststr(Vlevelloglistfiles.back());
				for ( uint64_t j = 0; j < Vlevellogs.back().size(); ++j )
					liststr << Vlevellogs.back()[j] << "\n";
				}
				{
				libmaus2::aio::OutputStreamInstance liststr(Vlevelsubloglistfiles.back());
				for ( uint64_t j = 0; j < Vlevelsublogs.back().size(); ++j )
					liststr << Vlevelsublogs.back()[j] << "\n";
				}

				depids = ndepids;
			}

			comOSI.reset();

			uint64_t failjobid = 0;
			if ( allids.size() )
			{
				std::ostringstream logfnstr;
				logfnstr << tmpgen.getFileName() << "_fail_job.log";
				std::string const logfn = logfnstr.str();

				std::ostringstream jobfnstr;
				jobfnstr << tmpgen.getFileName() << "_fail_job.job";

				std::ostringstream jobnamestr;
				jobnamestr << tmpname << "_fail_job";
				std::string const jobname = jobnamestr.str();

				std::ostringstream ostr;
				ostr << "#!/bin/bash\n";
				ostr << "#\n";
				ostr << "#SBATCH --job-name=" << jobname << "\n";
				ostr << "#SBATCH --output=" << logfn << "\n";
				ostr << "#\n";
				ostr << "#SBATCH --ntasks=1\n";
				ostr << "#SBATCH --time=1:00:00\n";
				ostr << "#SBATCH --mem=" << 1000 << "\n";
				ostr << "#SBATCH --cpus-per-task=1\n";

				ostr << "#SBATCH --dependency=";
				for ( uint64_t i = 0; i < allids.size(); ++i )
					ostr << ((i==0)?"":"?") << "afternotok:" << allids[i];
				ostr << "\n";

				ostr << "\n";

				ostr << "srun bash -c \"scancel";

				for ( uint64_t i = 0; i < allids.size(); ++i )
					ostr << " " << allids[i];

				ostr << "\"";

				ostr << "\n";

				libmaus2::aio::OutputStreamInstance::unique_ptr_type Pjob(new libmaus2::aio::OutputStreamInstance(jobfnstr.str()));
				(*Pjob) << ostr.str();
				Pjob.reset();

				std::cerr << ostr.str() << std::endl;

				std::ostringstream comstr;
				comstr << "sbatch " << jobfnstr.str();
				std::string const com = comstr.str();

				std::string out;
				std::string err;
				std::cerr << "[V] calling " << com << std::endl;
				int const r = libmaus2::util::PosixExecute::execute(com,out,err,false /* do not throw */);
				std::cerr << "[V] called " << com << std::endl;

				if ( r != EXIT_SUCCESS )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] failed to run " << com << "\n" << out << err << std::endl;
					lme.finish();
					throw lme;
				}

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
							uistr >> failjobid;

							std::cerr << "[V] submitted job " << failjobid << " for " << com << std::endl;
						}
					}
				}

				libmaus2::aio::FileRemoval::removeFile(jobfnstr.str());
			}

			if ( allids.size() )
			{
				std::ostringstream logfnstr;
				logfnstr << tmpgen.getFileName() << "_success_job.log";
				std::string const logfn = logfnstr.str();

				std::ostringstream jobfnstr;
				jobfnstr << tmpgen.getFileName() << "_success_job.job";

				std::ostringstream jobnamestr;
				jobnamestr << tmpname << "_success_job";
				std::string const jobname = jobnamestr.str();

				std::ostringstream ostr;
				ostr << "#!/bin/bash\n";
				ostr << "#\n";
				ostr << "#SBATCH --job-name=" << jobname << "\n";
				ostr << "#SBATCH --output=" << logfn << "\n";
				ostr << "#\n";
				ostr << "#SBATCH --ntasks=1\n";
				ostr << "#SBATCH --time=1:00:00\n";
				ostr << "#SBATCH --mem=" << 1000 << "\n";
				ostr << "#SBATCH --cpus-per-task=1\n";

				ostr << "#SBATCH --dependency=afterok";
				for ( uint64_t i = 0; i < allids.size(); ++i )
					ostr << ":" << allids[i];
				ostr << "\n";

				ostr << "\n";

				ostr << "srun bash -c \"scancel " << failjobid;

				for ( uint64_t j = 0; j < Vlevelloglistfiles.size(); ++j )
				{
					ostr << "; xargs <" << Vlevelloglistfiles[j] << " cat > " << Vlevellogcat[j];
					ostr << "; xargs <" << Vlevelloglistfiles[j] << " rm";
					ostr << "; rm " << Vlevelloglistfiles[j];
				}

				for ( uint64_t j = 0; j < Vlevelsubloglistfiles.size(); ++j )
				{
					ostr << "; xargs <" << Vlevelsubloglistfiles[j] << " cat > " << Vlevelsublogcat[j];
					ostr << "; xargs <" << Vlevelsubloglistfiles[j] << " rm";
					ostr << "; rm " << Vlevelsubloglistfiles[j];
				}

				ostr << "; xargs <" << commandlistfn << " rm";
				ostr << "; rm -f " << commandlistfn;

				ostr << "\"\n";

				libmaus2::aio::OutputStreamInstance::unique_ptr_type Pjob(new libmaus2::aio::OutputStreamInstance(jobfnstr.str()));
				(*Pjob) << ostr.str();
				Pjob.reset();

				std::cerr << ostr.str() << std::endl;

				std::ostringstream comstr;
				comstr << "sbatch " << jobfnstr.str();
				std::string const com = comstr.str();

				std::string out;
				std::string err;
				std::cerr << "[V] calling " << com << std::endl;
				int const r = libmaus2::util::PosixExecute::execute(com,out,err,false /* do not throw */);
				std::cerr << "[V] called " << com << std::endl;

				if ( r != EXIT_SUCCESS )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] failed to run " << com << "\n" << out << err << std::endl;
					lme.finish();
					throw lme;
				}

				libmaus2::aio::FileRemoval::removeFile(jobfnstr.str());
			}
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			std::ostream & errstr = lme.getStream();

			errstr << "[E] " << argv[0] << " failed" << std::endl;

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
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
