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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <cstring>
#include <libmaus2/util/CommandContainer.hpp>
#include <slurm/slurm.h>

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

std::string whiche(std::string const & prog)
{
	try
	{
		std::string const s = which(prog);
		std::cerr << "[V] found " << prog << " at " << s << "\n";
		return s;
	}
	catch(...)
	{
		return std::string();
	}
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

std::string runProgram(std::vector<std::string> const & args, libmaus2::util::TempFileNameGenerator & tmpgen)
{
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

	std::string const HPC_daligner = which(args[0]);

	std::string const errfile = tmpgen.getFileName() + "_HPC.err";
	libmaus2::util::TempFileRemovalContainer::addTempFile(errfile);
	std::string const outfile = tmpgen.getFileName() + "_HPC.out";
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
		std::cerr << "[V] " << HPC_daligner << " finished ok" << "\n";
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

void parseBatches(std::string const & HPCdalout, std::vector < std::pair< std::string, std::vector<std::string> > > & Vbatch)
{
	std::istringstream ISI(HPCdalout);
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
}

std::string basename(std::string const & s)
{
	if ( s.find('/') == std::string::npos )
		return s;
	else
		return s.substr(s.find_last_of('/')+1);
}

std::string dirname(std::string const & s)
{
	if ( s.find('/') == std::string::npos )
		return ".";
	else
		return s.substr(0,s.find_last_of('/'));
}

std::string clipdbdam(std::string const & s)
{
	if ( s.rfind(".db") == s.size() - 3 )
		return s.substr(0,s.size()-3);
	else if ( s.rfind(".dam") == s.size() - 4 )
		return s.substr(0,s.size()-4);
	else
		return s;
}

std::string moveDB(std::string const & from, std::string const & to)
{
	std::string const fromdir = dirname(from);
	std::string const fromdbbase = basename(from);
	std::string const frombase = clipdbdam(fromdbbase);

	std::string const todir = dirname(to);
	std::string const todbbase = basename(to);
	std::string const tobase = clipdbdam(todbbase);

	std::ostringstream ostr;
	ostr << "mv " << from << ".db" << " " << to;
	ostr << " ; mv " << fromdir << "/." << frombase << ".idx " << todir << "/." << tobase << ".idx";
	ostr << " ; mv " << fromdir << "/." << frombase << ".bps " << todir << "/." << tobase << ".bps";
	ostr << " ; for i in " << fromdir << "/." << frombase << ".*.data ; do I=${i#" << fromdir << "/} ; I=${I#." << frombase << ".} ; I=${I%.data}";
	ostr << " ; mv " << fromdir << "/." << frombase << ".${I}.data " << todir << "/." << tobase << ".${I}.data";
	ostr << " ; mv " << fromdir << "/." << frombase << ".${I}.anno " << todir << "/." << tobase << ".${I}.anno";
	ostr << " ; done ";

	return ostr.str();
}

uint64_t submitJob(
	std::string command,
	libmaus2::util::TempFileNameGenerator & tmpgen,
	std::string const & tmpname,
	uint64_t const i,
	uint64_t const j,
	std::vector< std::vector<std::string> > & Vlevellogs,
	std::vector< std::vector<std::string> > & Vlevelsublogs,
	uint64_t const mem,
	uint64_t const lnumthreads,
	std::vector<uint64_t> const & depids,
	std::ostream & comOSI,
	bool const deppos
)
{
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
		OSI << "#!/bin/bash" << "\n";
		OSI << command << "\n";
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
		if ( deppos )
		{
			ostr << "#SBATCH --dependency=afterok";
			for ( uint64_t i = 0; i < depids.size(); ++i )
				ostr << ":" << depids[i];
			ostr << "\n";
		}
		else
		{
			ostr << "#SBATCH --dependency=";
			for ( uint64_t i = 0; i < depids.size(); ++i )
				ostr << ((i==0)?"":"?") << "afternotok:" << depids[i];
			ostr << "\n";
		}
	}

	ostr << "\n";
	ostr << "srun bash " << commandfn << " >" << sublogfn << " 2>&1" << "\n";

	comOSI << commandfn << "\n";

	libmaus2::aio::OutputStreamInstance::unique_ptr_type Pjob(new libmaus2::aio::OutputStreamInstance(jobfnstr.str()));
	(*Pjob) << ostr.str();
	Pjob.reset();

	std::cerr << ostr.str() << "\n";

	std::ostringstream comstr;
	comstr << "sbatch " << jobfnstr.str();
	std::string const com = comstr.str();

	std::string out;
	std::string err;
	std::cerr << "[V] calling " << com << "\n";
	int const r = libmaus2::util::PosixExecute::execute(com,out,err,false /* do not throw */);
	std::cerr << "[V] called " << com << "\n";

	if ( r != EXIT_SUCCESS )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] failed to run " << com << "\n";
		lme.finish();
		throw lme;
	}

	std::istringstream sistr(out);
	int64_t njobid = -1;
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

				std::cerr << "[V] submitted job " << jobid << " for " << com << "\n";
			}
		}
	}

	if ( njobid < 0 )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] failed to run " << com << " (job id not found)" << "\n";
		lme.finish();
		throw lme;
	}

	libmaus2::aio::FileRemoval::removeFile(jobfnstr.str());

	return njobid;
}

struct SlurmControlConfig
{
	slurm_ctl_conf_t * conf;

	SlurmControlConfig()
	: conf(0)
	{
		int const r = slurm_load_ctl_conf(static_cast<time_t>(0),&conf);

		if ( r != 0 )
		{
			int const error = slurm_get_errno();
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] slurm_load_ctl_conf failed: " << slurm_strerror(error) << std::endl;
			lme.finish();
			throw lme;
		}
	}

	~SlurmControlConfig()
	{
		slurm_free_ctl_conf(conf);
	}

	uint64_t getMaxArraySize() const
	{
		return conf->max_array_sz;
	}
};

struct SlurmPartitions
{
	partition_info_msg_t * partitions;

	SlurmPartitions()
	{
		int const r = slurm_load_partitions(static_cast<time_t>(0),&partitions,0/*flags*/);

		if ( r != 0 )
		{
			int const error = slurm_get_errno();
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] slurm_load_partitions failed: " << slurm_strerror(error) << std::endl;
			lme.finish();
			throw lme;
		}
	}

	~SlurmPartitions()
	{
		slurm_free_partition_info_msg(partitions);
	}

	uint64_t size() const
	{
		return partitions->record_count;
	}

	std::string getName(uint64_t const i) const
	{
		assert ( i < size() );
		return partitions->partition_array[i].name;
	}

	std::string getNodes(uint64_t const i) const
	{
		assert ( i < size() );
		return partitions->partition_array[i].nodes;
	}

	uint64_t getTotalCpus(uint64_t const i) const
	{
		assert ( i < size() );
		return partitions->partition_array[i].total_cpus;
	}

	uint64_t getTotalNodes(uint64_t const i) const
	{
		assert ( i < size() );
		return partitions->partition_array[i].total_nodes;
	}

	uint64_t getMaxTime(uint64_t const i) const
	{
		assert ( i < size() );
		return partitions->partition_array[i].max_time;
	}
};

static std::string updateCommand(std::string command, uint64_t & lnumthreads)
{

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

	command = fullpath(command);

	return command;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		long const slurmapi = slurm_api_version();

		std::cerr << "[V] slurm API version is " << slurmapi << std::endl;

		SlurmControlConfig slurmconf;

		std::cerr << "[V] maximum job array size is " << slurmconf.getMaxArraySize() << std::endl;

		SlurmPartitions slurmpart;
		for ( uint64_t i = 0; i < slurmpart.size(); ++i )
		{
			std::cerr << "partition " << i << " name " << slurmpart.getName(i)
				<< " nodes " << slurmpart.getNodes(i)
				<< " numnodes=" << slurmpart.getTotalNodes(i)
				<< " numcpus=" << slurmpart.getTotalCpus(i)
				<< " maxtime=" << slurmpart.getMaxTime(i)
				<< std::endl;
		}

		std::string const tmpname = libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
		std::vector<std::string> Vsbatch;

		bool const havegenomesize = arg.uniqueArgPresent("g");
		
		std::string const outdbname = arg.uniqueArgPresent("O") ? arg["O"] : std::string("out.db");
		// bool const havecoverage = arg.uniqueArgPresent("c");
		// int64_t const coverage = havecoverage ? arg.getParsedArg<unsigned int>("c") : -1;

		std::string const DASedit = whiche("DASedit");
		std::string const DASpatch = whiche("DASpatch");
		std::string const DAStrim = whiche("DAStrim");
		std::string dasqv = whiche("DASqv");
		std::string cattrack = which("Catrack");
		std::string DBsplit = which("DBsplit");
		std::string HPC_TANmask = which("HPC.TANmask");
		std::string HPC_REPmask = which("HPC.REPmask");
		std::string cocommand = which("cocommand");

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
			std::cerr << "[E] usage: " << argv[0] << " <opts> -- <prog> ..." << "\n";
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

		std::string const dbname = argdal[0];
		std::string const dbbase = clipdbdam(dbname);

		std::cerr << "[V] db=" << dbname << " dbbase=" << dbbase << "\n";

		uint64_t blocksize = 0;
		uint64_t cutoff = 0;
		uint64_t basesum = 0;
		{
			libmaus2::dazzler::db::DatabaseFile DB(dbname);
			DB.computeTrimVector();

			blocksize = DB.blocksize;
			cutoff = DB.cutoff;
			basesum = DB.computeReadLengthSum();
		}

		std::cerr << "[V] blocksize=" << blocksize << " cutoff=" << cutoff << " basesum=" << basesum << "\n";

		int64_t coverage = havegenomesize ?
			( std::floor((static_cast<double>(basesum) / arg.getUnsignedNumericArg<uint64_t>(std::string("g"))) + 0.5) ) : -1;
		bool const havecoverage = coverage > 0;

		std::cerr << "[V] basesum=" << basesum << " coverage=" << coverage << "\n";

		uint64_t const bdiv = 1000ull * 1000000ull;
		uint64_t const bs = (basesum + bdiv - 1)/bdiv;
		uint64_t const numblocks = (basesum + (bs*1000000) - 1)/(bs*1000000);

		std::cerr << "[V] expected blocks " << numblocks << "\n";

		std::string const tmpfileprefix = getTmpFileBase(arg);
		libmaus2::util::TempFileNameGenerator tmpgen(tmpfileprefix + "_subdir",5);

		std::vector<std::string> tanargs;
		tanargs.push_back(HPC_TANmask);
		{
			std::ostringstream ostr;
			ostr << "-T" << numthreads;
			tanargs.push_back(ostr.str());
		}
		tanargs.push_back(dbname);

		std::vector < std::pair< std::string, std::vector<std::string> > > Vbatch;

		std::string const HPCtanout = runProgram(tanargs,tmpgen);
		parseBatches(HPCtanout,Vbatch);

		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("Catrack", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << cattrack << " " << argdal[0] << " " << "tan"; ostr << ";";
			ostr << "rm -f ." << dbbase << ".*.tan.*"; ostr << ";";
			Vbatch.back().second.push_back(ostr.str());
		}

		// 1:   (1*2*c)/numblocks   * 333
		// 10:  (10*2*c)/numblocks  * 25
		// 100: (100*2*c)/numblocks * 1.6

		/*
		 * REP1
		 */
		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DBsplit", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << bs << " " << dbname;
			Vbatch.back().second.push_back(ostr.str());
		}

		std::vector<std::string> repargs1;
		repargs1.push_back(HPC_REPmask);
		repargs1.push_back("-d");
		repargs1.push_back("-mtan");
		{
			std::ostringstream ostr;
			ostr << "-T" << numthreads;
			repargs1.push_back(ostr.str());
		}
		repargs1.push_back("-g1");
		{
			std::ostringstream ostr;
			ostr << "-c" << (1  *2*coverage*1200)/numblocks;
			repargs1.push_back(ostr.str());
		}
		repargs1.push_back(dbname);

		{
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << bs << " " << dbname;
			int const r = system(ostr.str().c_str());
			if ( r != 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] failed to run DBsplit " << ostr.str() << "\n";
				lme.finish();
				throw lme;
			}
		}
		std::string const HPCrep1out = runProgram(repargs1,tmpgen);
		{
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << (blocksize/1000000) << " " << dbname;
			int const r = system(ostr.str().c_str());
			if ( r != 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] failed to run DBsplit" << "\n";
				lme.finish();
				throw lme;
			}
		}

		parseBatches(HPCrep1out,Vbatch);

		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("work clean", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << "rm -fR work* temp*\n";
			Vbatch.back().second.push_back(ostr.str());
		}

		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("Catrack", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << cattrack << " " << argdal[0] << " " << "rep1"; ostr << ";";
			ostr << "rm -f ." << dbbase << ".*.rep1.*"; ostr << ";";
			Vbatch.back().second.push_back(ostr.str());
		}

		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DBsplit", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << (blocksize/1000000) << " " << dbname;
			Vbatch.back().second.push_back(ostr.str());
		}


		/*
		 * REP10
		 */
		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DBsplit", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << bs << " " << dbname;
			Vbatch.back().second.push_back(ostr.str());
		}

		std::vector<std::string> repargs10;
		repargs10.push_back(HPC_REPmask);
		repargs10.push_back("-d");
		repargs10.push_back("-mtan");
		repargs10.push_back("-mrep1");
		{
			std::ostringstream ostr;
			ostr << "-T" << numthreads;
			repargs10.push_back(ostr.str());
		}
		repargs10.push_back("-g10");
		{
			std::ostringstream ostr;
			ostr << "-c" << (10 *2*coverage*90)/numblocks;
			repargs10.push_back(ostr.str());
		}
		repargs10.push_back(dbname);

		{
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << bs << " " << dbname;
			int const r = system(ostr.str().c_str());
			if ( r != 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] failed to run DBsplit " << ostr.str() << "\n";
				lme.finish();
				throw lme;
			}
			// system((std::string("cat ") + dbname).c_str());
			std::cerr << "numblocks=" << numblocks << " " << 72 * static_cast<double>(2*coverage) / numblocks << "\n";
		}
		std::string const HPCrep10out = runProgram(repargs10,tmpgen);
		{
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << (blocksize/1000000) << " " << dbname;
			int const r = system(ostr.str().c_str());
			if ( r != 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] failed to run DBsplit" << "\n";
				lme.finish();
				throw lme;
			}
		}
		parseBatches(HPCrep10out,Vbatch);

		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("work clean", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << "rm -fR work* temp*\n";
			Vbatch.back().second.push_back(ostr.str());
		}

		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("Catrack", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << cattrack << " " << argdal[0] << " " << "rep10"; ostr << ";";
			ostr << "rm -f ." << dbbase << ".*.rep10.*"; ostr << ";";
			Vbatch.back().second.push_back(ostr.str());
		}


		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DBsplit", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << (blocksize/1000000) << " " << dbname;
			Vbatch.back().second.push_back(ostr.str());
		}

		/*
		 * REP100
		 */
		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DBsplit", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << bs << " " << dbname;
			Vbatch.back().second.push_back(ostr.str());
		}

		std::vector<std::string> repargs100;
		repargs100.push_back(HPC_REPmask);
		repargs100.push_back("-d");
		repargs100.push_back("-mtan");
		repargs100.push_back("-mrep1");
		repargs100.push_back("-mrep10");
		{
			std::ostringstream ostr;
			ostr << "-T" << numthreads;
			repargs10.push_back(ostr.str());
			repargs100.push_back(ostr.str());
		}
		repargs100.push_back("-g100");
		{
			std::ostringstream ostr;
			ostr << "-c" << (100*2*coverage*6)/numblocks;
			repargs100.push_back(ostr.str());
		}
		repargs100.push_back(dbname);

		{
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << bs << " " << dbname;
			int const r = system(ostr.str().c_str());
			if ( r != 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] failed to run DBsplit " << ostr.str() << "\n";
				lme.finish();
				throw lme;
			}
		}
		std::string const HPCrep100out = runProgram(repargs100,tmpgen);
		{
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << (blocksize/1000000) << " " << dbname;
			int const r = system(ostr.str().c_str());
			if ( r != 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] failed to run DBsplit" << "\n";
				lme.finish();
				throw lme;
			}
		}
		parseBatches(HPCrep100out,Vbatch);

		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("work clean", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << "rm -fR work* temp*\n";
			Vbatch.back().second.push_back(ostr.str());
		}

		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("Catrack", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << cattrack << " " << argdal[0] << " " << "rep100"; ostr << ";";
			ostr << "rm -f ." << dbbase << ".*.rep100.*"; ostr << ";";
			Vbatch.back().second.push_back(ostr.str());
		}


		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DBsplit", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << DBsplit << " -f -x" << cutoff << " -s" << (blocksize/1000000) << " " << dbname;
			Vbatch.back().second.push_back(ostr.str());
		}

		#if 0
		std::cerr << "1:   " <<  (1  *2*coverage*1200)/numblocks << "\n";
		std::cerr << "10:  " <<  (10 *2*coverage*90)/numblocks << "\n";
		std::cerr << "100: " <<  (100*2*coverage*6)/numblocks << "\n";
		#endif

		args.push_back("-mtan");
		args.push_back("-mrep1");
		args.push_back("-mrep10");
		args.push_back("-mrep100");

		std::string const HPCdalout = runProgram(args,tmpgen);
		parseBatches(HPCdalout,Vbatch);
		{
			Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("work clean", std::vector<std::string>()));
			std::ostringstream ostr;
			ostr << "rm -fR work* temp*\n";
			Vbatch.back().second.push_back(ostr.str());
		}


		for ( uint64_t i = 0; i < Vbatch.size(); ++i )
			if ( Vbatch[i].second.size() == 0 )
				Vbatch[i].second.push_back(std::string("echo OK"));

		std::ostringstream patchfnstr;
		patchfnstr << tmpgen.getFileName() << "_ndb";
		std::string const patchfn = patchfnstr.str();

		// find last merging batch
		int64_t lastmerge = -1;
		for ( uint64_t i = 0; i < Vbatch.size(); ++i )
			if ( Vbatch[i].second.size() && Vbatch[i].second[0].find("LAmerge") == 0 )
				lastmerge = i;
		// extract output files
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
						std::cerr << "Last merge " << lastmerge << " " << s << "\n";
						Voutputfiles.push_back(s);
					}
				}
			}
		}
		// add dasQV etc. call

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
				ostr << ";";
				ostr << "rm -f ." << dbbase << ".*.qual.*";
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
						ostr << cattrack << " " << argdal[0] << " " << "trim"; ostr << ";"; ostr << "rm -f ." << dbbase << ".*.trim.*"; ostr << ";";
						ostr << cattrack << " " << argdal[0] << " " << "hq"; ostr << ";"; ostr << "rm -f ." << dbbase << ".*.hq.*"; ostr << ";";
						ostr << cattrack << " " << argdal[0] << " " << "hole"; ostr << ";"; ostr << "rm -f ." << dbbase << ".*.hole.*"; ostr << ";";
						ostr << cattrack << " " << argdal[0] << " " << "span"; ostr << ";"; ostr << "rm -f ." << dbbase << ".*.span.*"; ostr << ";";
						ostr << cattrack << " " << argdal[0] << " " << "split"; ostr << ";"; ostr << "rm -f ." << dbbase << ".*.split.*"; ostr << ";";
						ostr << cattrack << " " << argdal[0] << " " << "adapt"; ostr << ";"; ostr << "rm -f ." << dbbase << ".*.adapt.*"; ostr << ";";
						ostr << cattrack << " " << argdal[0] << " " << "keep"; ostr << ";"; ostr << "rm -f ." << dbbase << ".*.keep.*"; ostr << ";";
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
						ostr << ";";
						ostr << "rm -f ." << dbbase << ".*.patch.*";
						Vbatch.back().second.push_back(ostr.str());
					}
				}

				if ( DASedit.size() )
				{
					Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("DASedit", std::vector<std::string>()));

					std::ostringstream ostr;
					ostr << DASedit << " -v -x2000 " << " " << argdal[0] << " " << patchfn;
					Vbatch.back().second.push_back(ostr.str());

					{
						std::string const movedbfn = tmpgen.getFileName();

						{
							libmaus2::aio::OutputStreamInstance OSI(movedbfn);
							OSI << moveDB(patchfn,outdbname);
						}

						Vbatch.push_back(std::pair< std::string, std::vector<std::string> >("movedb", std::vector<std::string>()));
						Vbatch.back().second.push_back(
							std::string("bash ") + movedbfn
						);
					}
				}
			}
		}

		{
			libmaus2::aio::OutputStreamInstance OSI("batches.debug");
			for ( uint64_t i = 0; i < Vbatch.size(); ++i )
			{
				OSI << "#" << "[" << i << "]: " << Vbatch[i].first << "\n";
				for ( uint64_t j = 0; j < Vbatch[i].second.size(); ++j )
					OSI << Vbatch[i].second[j] << "\n";
			}
			OSI << std::flush;
		}

		std::vector<uint64_t> depids;
		// iterate over batches
		for ( uint64_t i = 0; i < Vbatch.size(); ++i )
		{
			// size of batch in commands
			uint64_t const batchsize = Vbatch[i].second.size();
			uint64_t const maxjobsperbatch = 512;
			uint64_t const linesperjob = (batchsize + maxjobsperbatch - 1)/maxjobsperbatch;

			std::vector< std::vector<std::string> > Vpacks;
			uint64_t lnumthreads = 1;
			uint64_t const numjobs = (batchsize + linesperjob - 1)/linesperjob;
			for ( uint64_t j = 0; j < numjobs; ++j )
			{
				uint64_t const low = j * linesperjob;
				uint64_t const high = std::min(low+linesperjob,static_cast<uint64_t>(Vbatch[i].second.size()));

				Vpacks.push_back(std::vector<std::string>(0));
				for ( uint64_t k = low; k < high; ++k )
					Vpacks.back().push_back(updateCommand(Vbatch[i].second[k],lnumthreads));
			}

			assert ( Vpacks.size() == numjobs );

			// uint64_t const numjobs = Vbatch[i].second.size();
			uint64_t const jobsperarray = slurmconf.getMaxArraySize();
			uint64_t const numarrays = (numjobs + jobsperarray - 1)/jobsperarray;

			std::vector<uint64_t> ndepids;

			for ( uint64_t zz = 0; zz < numarrays; ++zz )
			{
				uint64_t const jlow = zz * jobsperarray;
				uint64_t const jhigh = std::min(jlow + jobsperarray, numjobs);

				std::ostringstream tmppreostr;
				tmppreostr << tmpgen.getFileName() << "_" << zz;
				std::string const tmppre = tmppreostr.str();

				std::string const contfn = tmppre + ".cont";
				std::string const jobfn = tmppre + ".job";
				// libmaus2::aio::OutputStreamInstance::unique_ptr_type Pcont(new libmaus2::aio::OutputStreamInstance(contfn));

				libmaus2::util::CommandContainer CC;

				for ( uint64_t j = jlow; j < jhigh; ++j )
				{
					std::string const tmppre = tmpgen.getFileName();

					std::string const in = "/dev/null";

					std::ostringstream outstr;
					outstr << tmppre << "_" << i << "_" << j << ".out";
					std::ostringstream errstr;
					errstr << tmppre << "_" << i << "_" << j << ".err";
					std::ostringstream comstr;
					comstr << tmppre << "_" << i << "_" << j << ".command";

					std::vector<std::string> const & pack = Vpacks[j];

					{
						libmaus2::aio::OutputStreamInstance OSI(comstr.str());

						for ( uint64_t i = 0; i < pack.size(); ++i )
						{
							OSI << pack[i] << "\n";
							OSI << "RET=$?\n";
							OSI << "echo \"return code ${RET}\" >/dev/stderr\n";
							OSI << "if [ $RET -ne 0 ] ; then\n";
							OSI << "\texit ${RET}\n";
							OSI << "fi\n";
						}

						OSI << "exit 0\n";

						OSI.flush();
					}

					std::vector<std::string> args;

					args.push_back(which("bash"));
					args.push_back(comstr.str());

					libmaus2::util::Command C(in,outstr.str(),errstr.str(),args);
					CC.V.push_back(C);

					// Largs
				}

				{
					libmaus2::aio::OutputStreamInstance cOSI(contfn);
					CC.serialise(cOSI);
					cOSI.flush();
				}


				std::ostringstream ostr;
				ostr << "#!/bin/bash\n";
				ostr << "#\n";
				ostr << "#SBATCH --job-name=SLURM_HPC_daligner_" << i << "_" << zz << "\n";
				ostr << "#SBATCH --output=/dev/null\n";
				ostr << "#SBATCH --error=/dev/null\n";
				//ostr << "#SBATCH --output=out_%A_%a.out\n";
				//ostr << "#SBATCH --error=err_%A_%a.err\n";
				ostr << "#SBATCH --array=0-" << jhigh-jlow-1 << "\n";
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
				ostr << "srun " << cocommand << " " << contfn << "\n";

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
						out = runProgram(sbatchargs,tmpgen);
						outok = true;
					}
					catch(std::exception const & ex)
					{
						std::cerr << ex.what() << std::endl;
					}
				}

				std::istringstream sistr(out);
				int64_t njobid = -1;
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

				ndepids.push_back(njobid);

				libmaus2::aio::FileRemoval::removeFile(jobfn);
			}

			depids = ndepids;
		}

		#if 0
		// exit(0);

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

			std::cerr << Vbatch[i].first << "\n";
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

				uint64_t const jobid = submitJob(command,tmpgen,tmpname,i,j,Vlevellogs,Vlevelsublogs,mem,lnumthreads,depids,*comOSI,true /* deppos */);
				ndepids.push_back(jobid);
				allids.push_back(jobid);

				#if 0
				// XYZ
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
					OSI << "#!/bin/bash" << "\n";
					OSI << command << "\n";
					// OSI << "rm -f " << commandfn < "\n";
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

				(*comOSI) << commandfn << "\n";

				libmaus2::aio::OutputStreamInstance::unique_ptr_type Pjob(new libmaus2::aio::OutputStreamInstance(jobfnstr.str()));
				(*Pjob) << ostr.str();
				Pjob.reset();

				std::cerr << ostr.str() << "\n";

				std::ostringstream comstr;
				comstr << "sbatch " << jobfnstr.str();
				std::string const com = comstr.str();

				std::string out;
				std::string err;
				std::cerr << "[V] calling " << com << "\n";
				int const r = libmaus2::util::PosixExecute::execute(com,out,err,false /* do not throw */);
				std::cerr << "[V] called " << com << "\n";

				if ( r != EXIT_SUCCESS )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] failed to run " << com << "\n";
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

							std::cerr << "[V] submitted job " << jobid << " for " << com << "\n";
						}
					}
				}

				libmaus2::aio::FileRemoval::removeFile(jobfnstr.str());
				#endif
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

			assert ( ndepids.size() );

			uint64_t zz = 0;
			while ( ndepids.size() > 1 )
			{
				std::vector<uint64_t> nndepids;

				uint64_t const quant = 64;
				for ( uint64_t q = 0; q < ndepids.size(); q += quant )
				{
					uint64_t const low = q;
					uint64_t const high = std::min(q+quant,static_cast<uint64_t>(ndepids.size()));

					uint64_t const jj = Vbatch[i].second.size() + (zz++);
					std::vector<uint64_t> ldepids(ndepids.begin()+low,ndepids.begin()+high);
					uint64_t const jobid = submitJob(std::string("echo OK"),tmpgen,tmpname,i,jj,Vlevellogs,Vlevelsublogs,mem,1/*num threads */,ldepids,*comOSI,true /* deppos */);
					nndepids.push_back(jobid);
					allids.push_back(jobid);

					std::cerr << "[V] reducing [" << low << "," << high << ")" << " to " << jobid << "\n";
				}

				if ( ndepids.size() % 2 ==  1)
					nndepids.push_back(ndepids.back());

				ndepids = nndepids;
			}

			depids = ndepids;
		}

		comOSI.reset();

		std::vector<uint64_t> failids = allids;
		{

			uint64_t zz = 0;
			while ( failids.size() > 1 )
			{
				std::vector<uint64_t> nfailids;

				uint64_t const quant = 64;
				for ( uint64_t q = 0; q < failids.size(); q += quant )
				{
					uint64_t const low = q;
					uint64_t const high = std::min(q+quant,static_cast<uint64_t>(failids.size()));

					std::vector<uint64_t> ldepids(failids.begin()+low,failids.begin()+high);
					uint64_t const jobid = submitJob(std::string("echo OK && exit 1"),tmpgen,tmpname,Vbatch.size(),zz++,Vlevellogs,Vlevelsublogs,mem,1/*num threads */,ldepids,*comOSI,false /* deppos */);
					nfailids.push_back(jobid);
					failids.push_back(jobid);

					std::cerr << "[V] reducing neg [" << low << "," << high << ")" << " to " << jobid << "\n";
				}

				failids = nfailids;
			}
		}

		uint64_t failjobid = 0;
		if ( failids.size() )
		{
			assert ( failids.size() == 1 );

			std::ostringstream logfnstr;
			logfnstr << tmpgen.getFileName() << "_fail_job.log";
			std::string const logfn = logfnstr.str();

			std::ostringstream commandfnstr;
			commandfnstr << tmpgen.getFileName() << "_fail_job.command";
			std::string const commandfn = commandfnstr.str();

			std::ostringstream allpidfnstr;
			allpidfnstr << tmpgen.getFileName() << "_fail_job.allpid";
			std::string const allpidfn = allpidfnstr.str();

			std::ostringstream jobfnstr;
			jobfnstr << tmpgen.getFileName() << "_fail_job.job";

			std::ostringstream jobnamestr;
			jobnamestr << tmpname << "_fail_job";
			std::string const jobname = jobnamestr.str();

			{
				libmaus2::aio::OutputStreamInstance OSI(allpidfn);
				for ( uint64_t i = 0; i < allids.size(); ++i )
					OSI << allids[i] << "\n";
			}

			{
				libmaus2::aio::OutputStreamInstance OSI(commandfn);
				OSI << "xargs scancel <" << allpidfn << "\n";
				OSI << DBsplit << " -f -x" << cutoff << " -s" << (blocksize/1000000) << " " << dbname << "\n";
				OSI << "rm -f " << commandfn << "\n";
			}

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

			ostr << "#SBATCH --dependency=afternotok:" << failids.at(0);
			ostr << "\n";

			ostr << "srun bash " << commandfn;

			libmaus2::aio::OutputStreamInstance::unique_ptr_type Pjob(new libmaus2::aio::OutputStreamInstance(jobfnstr.str()));
			(*Pjob) << ostr.str();
			Pjob.reset();

			std::cerr << ostr.str() << "\n";

			std::ostringstream comstr;
			comstr << "sbatch " << jobfnstr.str();
			std::string const com = comstr.str();

			std::string out;
			std::string err;
			std::cerr << "[V] calling " << com << "\n";
			int const r = libmaus2::util::PosixExecute::execute(com,out,err,false /* do not throw */);
			std::cerr << "[V] called " << com << "\n";

			if ( r != EXIT_SUCCESS )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] failed to run " << com << "\n" << out << err << "\n";
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

						std::cerr << "[V] submitted job " << failjobid << " for " << com << "\n";
					}
				}
			}

			libmaus2::aio::FileRemoval::removeFile(jobfnstr.str());
		}

		if ( depids.size() )
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
			for ( uint64_t i = 0; i < depids.size(); ++i )
				ostr << ":" << depids[i];
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

			std::string const movedbfn = tmpgen.getFileName();

			{
				libmaus2::aio::OutputStreamInstance OSI(movedbfn);
				OSI << moveDB(patchfn,"out.db");
			}

			ostr << "; bash " << movedbfn;

			ostr << ";" << DBsplit << " -f -x" << cutoff << " -s" << (blocksize/1000000) << " " << dbname;

			ostr << "\"\n";

			libmaus2::aio::OutputStreamInstance::unique_ptr_type Pjob(new libmaus2::aio::OutputStreamInstance(jobfnstr.str()));
			(*Pjob) << ostr.str();
			Pjob.reset();

			std::cerr << ostr.str() << "\n";

			std::ostringstream comstr;
			comstr << "sbatch " << jobfnstr.str();
			std::string const com = comstr.str();

			std::string out;
			std::string err;
			std::cerr << "[V] calling " << com << "\n";
			int const r = libmaus2::util::PosixExecute::execute(com,out,err,false /* do not throw */);
			std::cerr << "[V] called " << com << "\n";

			if ( r != EXIT_SUCCESS )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] failed to run " << com << "\n" << out << err << "\n";
				lme.finish();
				throw lme;
			}

			libmaus2::aio::FileRemoval::removeFile(jobfnstr.str());
		}
		#endif
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << "\n";
		return EXIT_FAILURE;
	}
}
