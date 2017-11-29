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
#include <libmaus2/util/ContainerDescriptionList.hpp>
#include <libmaus2/util/CommandContainer.hpp>
#include <sys/wait.h>

#include <sys/types.h>
#include <pwd.h>

#include <slurm/slurm.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " container" << std::endl;

	return ostr.str();
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
	
	uint64_t getIdForName(std::string const & name) const
	{
		for ( uint64_t i = 0; i < size(); ++i )
			if ( getName(i) == name )
				return i;
		
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] partition named " << name << " not found" << std::endl;
		lme.finish();
		throw lme;
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

	uint64_t getMaxMemPerCpu(uint64_t const i) const
	{
		assert ( i < size() );
		return partitions->partition_array[i].max_mem_per_cpu;
	}

	uint64_t getMaxCpusPerNode(uint64_t const i) const
	{
		assert ( i < size() );
		return partitions->partition_array[i].max_cpus_per_node;
	}
};

struct SlurmJobs
{
	job_info_msg_t * jobs;

	SlurmJobs()
	: jobs(0)
	{
		int const r = slurm_load_jobs(0,&jobs,0 /* showflags */);
		
		if ( r != 0 )
		{
			int const error = slurm_get_errno();
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] slurm_load_jobs failed: " << slurm_strerror(error) << std::endl;
			lme.finish();
			throw lme;
		}
	}
	
	~SlurmJobs()
	{
		slurm_free_job_info_msg(jobs);
	}
	
	uint64_t size() const
	{
		return jobs->record_count;
	}
	
	uint64_t getUserId(uint64_t const i) const
	{
		return jobs->job_array[i].user_id;
	} 

	std::string getUserName(uint64_t const i) const
	{
		uint64_t const uid = getUserId(i);
		struct passwd * pw = getpwuid(uid);

		if ( pw )
		{
			return std::string(pw->pw_name);
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] no user name found for uid " << uid << std::endl;
			lme.finish();
			throw lme;
		}
	} 

	char const * getName(uint64_t const i) const
	{
		return jobs->job_array[i].name;
	} 
	
	uint64_t getJobId(uint64_t const i) const
	{
		return jobs->job_array[i].job_id;
	}
};

static void writeJobDescription(
	std::string const & fn,
	std::string const & jobname,
	std::string const & outfn,
	uint64_t const utime,
	uint64_t const umem,
	uint64_t const threads,
	std::string const partition,
	std::string const command
)
{
	libmaus2::aio::OutputStreamInstance OSI(fn);
	
	OSI << "#!/bin/bash\n";
	OSI << "#SBATCH --job-name=" << jobname << "\n";
	OSI << "#SBATCH --output=" << outfn << "\n";
	OSI << "#SBATCH --ntasks=1" << "\n";
	OSI << "#SBATCH --time=" << utime << "\n";
	OSI << "#SBATCH --mem=" << umem << "\n";
	OSI << "#SBATCH --cpus-per-task=" << threads << "\n";
	OSI << "#SBATCH --cpus-per-task=" << threads << "\n";
	OSI << "#SBATCH --partition=" << partition << "\n";
	OSI << "srun bash -c \"" << command << "\"\n";
}

struct WorkerInfo
{
	int64_t id;
	libmaus2::network::SocketBase::unique_ptr_type Asocket;
	bool active;
	std::pair<int64_t,int64_t> packageid;
	uint64_t workerid;
	
	WorkerInfo() { reset(); }
	
	void reset()
	{
		id = -1;
		Asocket.reset();
		active = false;
		packageid = std::pair<int64_t,int64_t>(-1,-1);
		workerid = std::numeric_limits<uint64_t>::max();
	}
};

void startWorker(
	uint64_t & nextworkerid, std::string const & tmpfilebase,
	std::string const & hostname, uint64_t const serverport,
	uint64_t const workertime, uint64_t const workermem,
	uint64_t const workerthreads,
	std::string const & partition,
	libmaus2::util::ArgParser const & arg,
	WorkerInfo * AW,
	uint64_t const i,
	std::map<uint64_t,uint64_t> & idToSlot,
	uint64_t const workers
)
{
	uint64_t const workerid = nextworkerid++;
	std::ostringstream workernamestr;
	workernamestr << "worker_" << workerid;
	std::string const workername = workernamestr.str();

	std::ostringstream outfnstr;
	outfnstr << tmpfilebase << "_" << workerid << ".out";
	std::string const outfn = outfnstr.str();
	
	std::string const descname = tmpfilebase + "_worker.sbatch";
	
	std::ostringstream commandstr;
	commandstr << "slurmworker " << hostname << " " << serverport;
	std::string command = commandstr.str();
	
	writeJobDescription(
		descname,
		workername,
		outfn,
		workertime,
		workermem,
		workerthreads,
		partition,
		command
	);
	
	std::vector<std::string> Varg;
	Varg.push_back("sbatch");
	Varg.push_back(descname);
	
	std::string const jobid_s = runProgram(Varg,arg);
	
	std::deque<std::string> Vtoken = libmaus2::util::stringFunctions::tokenize(jobid_s,std::string(" "));
	
	if ( Vtoken.size() >= 4 )
	{
		std::istringstream istr(Vtoken[3]);
		uint64_t id;
		istr >> id;
		
		if ( istr && istr.peek() == '\n' )
		{
			// std::cerr << "got job id " << id << std::endl;
			
			AW [ i ].id = id;
			AW [ i ].workerid = workerid;
			idToSlot[id] = i;
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to find job id in " << jobid_s << std::endl;
			lme.finish();
			throw lme;
		}
	}
	else
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] unable to find job id in " << jobid_s << std::endl;
		lme.finish();
		throw lme;		
	}
					
	libmaus2::aio::FileRemoval::removeFile(descname);

	std::cerr << "[V] started job " << (i+1) << " out of " << workers << " with id " << AW[i].id << std::endl;
}

int slurmcontrol(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	std::string const hostname = libmaus2::network::GetHostName::getHostName();
	std::string const cdl = arg[0];
	
	unsigned short serverport = 50000;
	uint64_t const backlog = 1024;
	uint64_t const tries = 1000;
	uint64_t nextworkerid = 0;

	uint64_t const workertime = arg.uniqueArgPresent("workertime") ? arg.getParsedArg<uint64_t>("workertime") : 1440;
	uint64_t const workermem = arg.uniqueArgPresent("workermem") ? arg.getParsedArg<uint64_t>("workermem") : 40000;
	std::string const partition = arg.uniqueArgPresent("partition") ? arg["partition"] : "haswell";
	uint64_t const workers = arg.uniqueArgPresent("workers") ? arg.getParsedArg<uint64_t>("workers") : 16;

	libmaus2::autoarray::AutoArray<WorkerInfo> AW(workers);
	std::map<uint64_t,uint64_t> idToSlot;
	// std::vector < int64_t > Vworkerid(workers,-1);
	
	libmaus2::util::ContainerDescriptionList CDL;
	{
		libmaus2::aio::InputStreamInstance ISI(cdl);
		CDL.deserialise(ISI);
	}

	
	std::vector < libmaus2::util::ContainerDescription > & CDLV = CDL.V;
	
	std::vector < libmaus2::util::CommandContainer > VCC(CDLV.size());
	std::set < std::pair<uint64_t,uint64_t> > Sunfinished;
	std::map < uint64_t, uint64_t > Munfinished;
	for ( uint64_t i = 0; i < CDLV.size(); ++i )
	{
		libmaus2::aio::InputStreamInstance ISI(CDLV[i].fn);
		VCC[i].deserialise(ISI);
		
		if ( CDLV[i].missingdep == 0 )
		{
			for ( uint64_t j = 0; j < VCC[i].V.size(); ++j )
			{
				Sunfinished.insert(std::pair<uint64_t,uint64_t>(i,j));
				Munfinished[i]++;
			}
		}
	}
	
	uint64_t maxthreads = 1;
	for ( uint64_t i = 0; i < VCC.size(); ++i )
		maxthreads = std::max(maxthreads,VCC[i].threads);

	uint64_t const workerthreads = maxthreads;
	
	libmaus2::network::ServerSocket::unique_ptr_type Pservsock(
		libmaus2::network::ServerSocket::allocateServerSocket(
			serverport,
			backlog,
			hostname,
			tries
		)
	);

	std::cerr << "[V] hostname=" << hostname << " serverport=" << serverport << " number of container " << CDLV.size() << std::endl;

	SlurmControlConfig SCC;
	SlurmPartitions SP;
	
	uint64_t const partid = SP.getIdForName(partition);
	uint64_t const maxcores = SP.getMaxCpusPerNode(partid);

	
	for ( uint64_t i = 0; i < workers; ++i )
		startWorker(
			nextworkerid,tmpfilebase,hostname,serverport,
			workertime,workermem,workerthreads,partition,arg,AW.begin(),i,
			idToSlot,workers
		);
	
	std::set<uint64_t> restartSet;
	uint64_t pending = 0;
	
	while ( Sunfinished.size() || pending )
	{
		for ( std::set<uint64_t>::const_iterator it = restartSet.begin(); it != restartSet.end(); ++it )
		{
			uint64_t const i = *it;
			startWorker(
				nextworkerid,tmpfilebase,hostname,serverport,
				workertime,workermem,workerthreads,partition,arg,AW.begin(),i,
				idToSlot,workers
			);
		}
		restartSet.clear();
	
		std::cerr << "[V] Sunfinished.size()=" << Sunfinished.size() << " pending=" << pending << std::endl;
	
		std::vector<int> Vfd;
		Vfd.push_back(Pservsock->getFD());

		int pfd = -1;
		int const r = libmaus2::network::SocketBase::multiPoll(Vfd,pfd);

		if ( r > 0 )
		{
			if ( pfd == Pservsock->getFD() )
			{
				libmaus2::network::SocketBase::unique_ptr_type nptr = Pservsock->accept();
				uint64_t const jobid = nptr->readSingle<uint64_t>();
				
				std::cerr << "[V] accepted connection for jobid=" << jobid << std::endl;
				
				if ( idToSlot.find(jobid) != idToSlot.end() )
				{
					uint64_t const slot = idToSlot.find(jobid)->second;
					nptr->writeSingle<uint64_t>(AW[slot].workerid);
					
					if ( ! AW[slot].Asocket )
					{
						AW[slot].Asocket = UNIQUE_PTR_MOVE(nptr);
						AW[slot].active = true;
					}
					else
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] erratic worker trying to open third connection" << std::endl;
						lme.finish();
						throw lme;	
					}
				}
				else
				{
					std::cerr << "[V] job id unknown" << std::endl;
				}
			}
		}
		else if ( r < 0 )
		{
			int const error = errno;
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] poll failed: " << strerror(error) << std::endl;
			lme.finish();
			throw lme;		
		}
		else
		{
			for ( uint64_t i = 0; i < AW.size(); ++i )
			{
				try
				{
					if ( AW[i].active )
					{
						// worker asking for work
						uint64_t const rd = AW[i].Asocket->readSingle<uint64_t>();
						
						if ( rd == 0 )
						{
							if ( Sunfinished.size() )
							{
								std::pair<uint64_t,uint64_t> const currentid = *(Sunfinished.begin());
								Sunfinished.erase(currentid);
								libmaus2::util::Command const com = VCC[currentid.first].V[currentid.second];
								std::ostringstream ostr;
								com.serialise(ostr);
								// process command
								AW[i].packageid = currentid;
								AW[i].Asocket->writeSingle<uint64_t>(0);
								AW[i].Asocket->writeString(ostr.str());
								pending += 1;
								
								std::cerr << "[V] started " << com << " for " << currentid.first << "," << currentid.second << " on slot " << i << std::endl;
							}
							else
							{
								// stay idle
								AW[i].Asocket->writeSingle<uint64_t>(1);								
							}
						}
						else if ( rd == 1 )
						{
							pending -= 1;
						
							uint64_t const status = AW[i].Asocket->readSingle<uint64_t>();
							int const istatus = static_cast<int>(status);
							
							std::cerr << "[V] slot " << i << " reports job ended with istatus=" << istatus << std::endl;
							
							if ( WIFEXITED(istatus) && (WEXITSTATUS(istatus) == 0) )
							{
								std::cerr << "[V] updating command container " << AW[i].packageid.first << std::endl;
								
								Munfinished [ AW[i].packageid.first ] --;
								
								if ( ! Munfinished [AW[i].packageid.first] )
								{
									std::cerr << "[V] finished command container " << AW[i].packageid.first << std::endl;
									
									libmaus2::util::CommandContainer & CC = VCC[
										AW[i].packageid.first
									];
																	
									for ( uint64_t j = 0; j < CC.rdepid.size(); ++j )
									{
										uint64_t const k = CC.rdepid[j];
										
										assert ( CDLV[k].missingdep );
										
										CDLV[k].missingdep -= 1;
										
										if ( ! CDLV[k].missingdep )
										{
											std::cerr << "[V] activating container " << k << std::endl;
											for ( uint64_t j = 0; j < VCC[k].V.size(); ++j )
											{
												Sunfinished.insert(std::pair<uint64_t,uint64_t>(k,j));
												Munfinished[k]++;
											}
										}
									}
								}
								
								AW[i].packageid.first = -1;
								AW[i].packageid.second = -1;
							}
							else
							{
								// todo: check how often job failed
								Sunfinished.insert(AW[i].packageid);
							}
						}
						else if ( rd == 2 )
						{
							// std::cerr << "[V] got ping signal from " << i << std::endl;
							AW[i].Asocket->writeSingle<uint64_t>(0);
						}
						else
						{
							uint64_t const id = AW[i].id;
							AW[i].reset();
							idToSlot.erase(id);
							restartSet.insert(i);
							
							std::cerr << "[V] process for slot " << i << " jobid " << id << " is erratic" << std::endl;
						}
					}
				}
				catch(...)
				{
					uint64_t const id = AW[i].id;
					AW[i].reset();
					idToSlot.erase(id);
					restartSet.insert(i);
					
					std::cerr << "[V] exception for slot " << i << " jobid " << id << std::endl;				
				}
			}
			
			sleep(1);
		}
	}
	
	std::set<uint64_t> Sactive;
	for ( uint64_t i = 0; i < AW.size(); ++i )
		if ( AW[i].active )
			Sactive.insert(i);

	while ( Sactive.size() )
	{
		std::vector < uint64_t > Vterm;
		
		for ( std::set<uint64_t>::const_iterator it = Sactive.begin(); it != Sactive.end(); ++it )
		{
			uint64_t const i = *it;
		
			try
			{
				// worker asking for work
				uint64_t const rd = AW[i].Asocket->readSingle<uint64_t>();
				
				if ( rd == 0 )
				{
					// request terminate
					AW[i].Asocket->writeSingle<uint64_t>(2);
					AW[i].reset();
					Vterm.push_back(i);
				}
				else if ( rd == 1 )
				{
					std::cerr << "[V] slot " << i << " reports finished job with no jobs active" << std::endl;
				}
				else if ( rd == 2 )
				{
					// std::cerr << "[V] got ping signal from " << i << std::endl;
					AW[i].Asocket->writeSingle<uint64_t>(0);
				}
				else
				{
					uint64_t const id = AW[i].id;
					AW[i].reset();
					idToSlot.erase(id);
					restartSet.insert(i);
					
					std::cerr << "[V] process for slot " << i << " jobid " << id << " is erratic" << std::endl;
				}
			}
			catch(...)
			{
				uint64_t const id = AW[i].id;
				AW[i].reset();
				idToSlot.erase(id);
				restartSet.insert(i);
				
				std::cerr << "[V] exception for slot " << i << " jobid " << id << std::endl;				
			}
		}

		for ( uint64_t i = 0; i < Vterm.size(); ++i )
			Sactive.erase(Vterm[i]);
	}

	
	#if 0
	for ( uint64_t i = 0; i < SP.size(); ++i )
	{
		std::cerr << "[V]\t" << i << "\t" << SP.getName(i) << "\t" << SP.getNodes(i) << "\t" << SP.getTotalCpus(i) << "\t"<< SP.getMaxTime(i) << "\t" << SP.getMaxMemPerCpu(i) << "\t" << SP.getMaxCpusPerNode(i) << std::endl;
	}
	#endif
	
	std::cerr << "[V] partition " << partition << " max cores per node " << maxcores << std::endl;

	#if 0	
	SlurmJobs jobs;
	
	std::cerr << "number of jobs: " << jobs.size() << std::endl;
	
	for ( uint64_t i = 0; i < jobs.size(); ++i )
	{
		std::cerr << "[V]\t" << i << "\t" << jobs.getUserName(i) << "\t" << jobs.getName(i) << "\t" << jobs.getJobId(i) << std::endl;
	}
	#endif

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		
		int const r = slurmcontrol(arg);
		
		return r;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
