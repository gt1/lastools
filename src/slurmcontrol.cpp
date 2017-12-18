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
#include <sys/wait.h>
#include <sys/epoll.h>

#include <sys/types.h>
#include <pwd.h>

#if 0
#include <slurm/slurm.h>

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
#endif

#include <FDIO.hpp>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [<parameters>] <container.cdl>" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -t          : number of threads (defaults to number of cores on machine)\n";
	ostr << " -T          : prefix for temporary files (default: create files in current working directory)\n";
	ostr << " -workertime : time for workers in (default: 1440)\n";
	ostr << " -workermem  : memory for workers (default: 40000)\n";
	ostr << " -partition  : cluster partition (default: haswell)\n";
	ostr << " -workers    : number of workers (default: 16)\n";

	return ostr.str();
}


struct SlurmControl
{
	struct JobDescription
	{
		int64_t containerid;
		int64_t subid;

		JobDescription() : containerid(-1), subid(-1)
		{
		}
		JobDescription(uint64_t const rcontainerid, uint64_t const rsubid)
		: containerid(rcontainerid), subid(rsubid)
		{

		}

		bool operator<(JobDescription const & O) const
		{
			if ( containerid != O.containerid )
				return containerid < O.containerid;
			else if ( subid != O.subid )
				return subid < O.subid;
			else
				return false;
		}

		void reset()
		{
			containerid = -1;
			subid = -1;
		}
	};

	struct WorkerInfo
	{
		int64_t id;
		libmaus2::network::SocketBase::unique_ptr_type Asocket;
		bool active;
		JobDescription packageid;
		uint64_t workerid;

		WorkerInfo() { reset(); }

		void reset()
		{
			id = -1;
			Asocket.reset();
			active = false;
			workerid = std::numeric_limits<uint64_t>::max();
			resetPackageId();
		}

		void resetPackageId()
		{
			packageid.reset();
		}
	};

	struct StartWorkerRequest
	{
		uint64_t * nextworkerid;
		std::string tmpfilebase;
		std::string hostname;
		uint64_t serverport;
		uint64_t workertime;
		uint64_t workermem;
		uint64_t workerthreads;
		std::string partition;
		libmaus2::util::ArgParser const * arg;
		WorkerInfo * AW;
		uint64_t i;
		std::map<uint64_t,uint64_t> * idToSlot;
		uint64_t workers;

		StartWorkerRequest() {}
		StartWorkerRequest(
			uint64_t & rnextworkerid,
			std::string rtmpfilebase,
			std::string rhostname,
			uint64_t rserverport,
			uint64_t rworkertime,
			uint64_t rworkermem,
			uint64_t rworkerthreads,
			std::string rpartition,
			libmaus2::util::ArgParser const & rarg,
			WorkerInfo * rAW,
			uint64_t ri,
			std::map<uint64_t,uint64_t> & ridToSlot,
			uint64_t rworkers
		) :
			nextworkerid(&rnextworkerid),
			tmpfilebase(rtmpfilebase),
			hostname(rhostname),
			serverport(rserverport),
			workertime(rworkertime),
			workermem(rworkermem),
			workerthreads(rworkerthreads),
			partition(rpartition),
			arg(&rarg),
			AW(rAW),
			i(ri),
			idToSlot(&ridToSlot),
			workers(rworkers)
		{

		}

		void dispatch()
		{
			uint64_t const workerid = (*nextworkerid)++;
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

			std::string const jobid_s = runProgram(Varg,*arg);

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
					(*idToSlot)[id] = i;
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
	};


	struct EPoll
	{
		int fd;

		EPoll(int const
		#if defined(HAVE_EPOLL_CREATE) && !defined(HAVE_EPOLL_CREATE1)
			size
		#endif
		) : fd(-1)
		{
			#if defined(HAVE_EPOLL_CREATE1)
			fd = epoll_create1(0);

			if ( fd < 0 )
			{
				int const error = errno;

				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] epoll_create1() failed: " << strerror(error) << std::endl;
				lme.finish();
				throw lme;
			}
			#elif defined(HAVE_EPOLL_CREATE)
			fd = epoll_create(size);

			if ( fd < 0 )
			{
				int const error = errno;

				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] epoll_create1() failed: " << strerror(error) << std::endl;
				lme.finish();
				throw lme;
			}

			#else
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] EPoll: epoll interface not supported " << strerror(error) << std::endl;
			lme.finish();
			throw lme;
			#endif
		}

		~EPoll()
		{
			if ( fd >= 0 )
			{
				::close(fd);
				fd = -1;
			}
		}

		#if defined(HAVE_EPOLL_CREATE) || defined(HAVE_EPOLL_CREATE1)
		void add(int const addfd)
		{
			struct epoll_event ev;
			ev.events = EPOLLIN
				#if defined(EPOLLRDHUP)
				| EPOLLRDHUP
				#endif
				;
			ev.data.fd = addfd;

			while ( true )
			{
				int const r = epoll_ctl(
					fd,
					EPOLL_CTL_ADD,
					addfd,
					&ev
				);

				if ( r == 0 )
					break;

				int const error = errno;

				switch ( error )
				{
					case EINTR:
					case EAGAIN:
						break;
					default:
					{

						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] EPoll:add: epoll_ctl() failed: " << strerror(error) << std::endl;
						lme.finish();
						throw lme;
					}
				}
			}
		}

		void remove(int const remfd)
		{
			while ( true )
			{
				int const r = epoll_ctl(
					fd,
					EPOLL_CTL_DEL,
					remfd,
					NULL
				);

				if ( r == 0 )
					break;
				else
				{
					int const error = errno;

					switch ( error )
					{
						case EINTR:
						case EAGAIN:
							break;
						default:
						{
							libmaus2::exception::LibMausException lme;
							lme.getStream() << "[E] EPoll:remove: epoll_ctl() failed: " << strerror(error) << std::endl;
							lme.finish();
							throw lme;
						}
					}
				}
			}
		}

		bool wait(int & rfd, int const timeout = 1000 /* milli seconds */)
		{
			rfd = -1;

			struct epoll_event events[1];

			while ( true )
			{
				int const nfds = epoll_wait(fd, &events[0], sizeof(events)/sizeof(events[0]), timeout);

				if ( nfds < 0 )
				{
					int const error = errno;

					switch ( error )
					{
						case EINTR:
						case EAGAIN:
							break;
						default:
						{
							libmaus2::exception::LibMausException lme;
							lme.getStream() << "[E] EPoll:wait: epoll_wait() failed: " << strerror(error) << std::endl;
							lme.finish();
							throw lme;
						}
					}
				}

				if ( nfds == 0 )
				{
					return false;
				}
				else
				{
					rfd = events[0].data.fd;
					return true;
				}
			}
		}
		#else
		void add(int const)
		{
		}
		void remove(int const)
		{
		}
		bool wait(int &, int const = 1000 /* milli seconds */)
		{
		}
		#endif
	};

	struct ProgState
	{
		uint64_t numunfinished;
		uint64_t numpending;

		ProgState() : numunfinished(0), numpending(0) {}
		ProgState(
			uint64_t const rnumunfinished,
			uint64_t const rnumpending
		) : numunfinished(rnumunfinished), numpending(rnumpending)
		{

		}

		bool operator==(ProgState const & o) const
		{
			return numunfinished == o.numunfinished && numpending == o.numpending;
		}

		bool operator!=(ProgState const & o) const
		{
			return !operator==(o);
		}
	};

	std::string const curdir;
	unsigned short serverport;
	uint64_t const backlog;
	uint64_t const tries;
	uint64_t nextworkerid;
	std::string const hostname;
	std::string const tmpfilebase;
	uint64_t const workertime;
	uint64_t const workermem;
	std::string const partition;
	uint64_t const workers;
	std::string const cdl;

	libmaus2::autoarray::AutoArray<WorkerInfo> AW;
	std::map<uint64_t,uint64_t> idToSlot;
	std::map<int,uint64_t> fdToSlot;
	std::map < JobDescription, uint64_t > Mfail;

	libmaus2::util::ContainerDescriptionList CDL;
	std::vector < libmaus2::util::ContainerDescription > & CDLV;
	std::vector < libmaus2::util::CommandContainer > VCC;


	std::set < JobDescription > Sunfinished;
	std::set < JobDescription > Srunning;
	std::set<uint64_t> Sresubmit;
	uint64_t ndeepsleep;
	std::map < uint64_t, uint64_t > Munfinished;

	uint64_t const maxthreads;
	uint64_t const workerthreads;

	EPoll EP;

	libmaus2::network::ServerSocket::unique_ptr_type Pservsock;

	std::set<uint64_t> restartSet;
	std::set<uint64_t> wakeupSet;
	// uint64_t pending;

	ProgState pstate;
	bool failed;

	std::vector < StartWorkerRequest > Vreq;

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


	void processWakeupSet()
	{
		for ( std::set < uint64_t >::const_iterator it = wakeupSet.begin();
			it != wakeupSet.end(); ++it )
		{
			uint64_t const i = *it;
			std::cerr << "[V] sending wakeup to slot " << i << std::endl;

			FDIO fdio(AW[i].Asocket->getFD());
			fdio.writeNumber(1);
		}
		wakeupSet.clear();
	}

	void processResubmitSet()
	{
		for ( std::set < uint64_t >::const_iterator it = Sresubmit.begin();
			it != Sresubmit.end(); ++it )
		{
			uint64_t const i = *it;
			std::cerr << "[V] resubmitting slot " << i << " after deep sleep" << std::endl;
			Vreq[i].dispatch();

		}
		Sresubmit.clear();
	}

	uint64_t computeMaxThreads()
	{
		uint64_t maxthreads = 1;
		for ( uint64_t i = 0; i < VCC.size(); ++i )
			maxthreads = std::max(maxthreads,VCC[i].threads);
		return maxthreads;
	}

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

	void countUnfinished()
	{
		// count number of unfinished jobs per command container
		for ( uint64_t i = 0; i < CDLV.size(); ++i )
		{
			libmaus2::util::CommandContainer & CC = VCC[i];
			uint64_t numunfinished = 0;

			for ( uint64_t j = 0; j < CC.V.size(); ++j )
				if ( CC.V[j].completed )
				{

				}
				else
				{
					Munfinished[i]++;
					numunfinished += 1;
				}

			std::cerr << "[V] container " << i << " has " << numunfinished << " unfinished jobs" << std::endl;

			if ( numunfinished )
			{
				for ( uint64_t j = 0; j < CC.rdepid.size(); ++j )
				{
					uint64_t const k = CC.rdepid[j];

					std::cerr << "[V] container " << k << " has missing dependency " << i << std::endl;

					CDLV[k].missingdep += 1;
				}
			}
		}
	}

	void addUnfinished(JobDescription const J)
	{
		Sunfinished.insert(J);
	}

	JobDescription getUnfinished()
	{
		std::set< JobDescription >::const_iterator const it = Sunfinished.begin();
		JobDescription const P = *it;
		Sunfinished.erase(it);
		return P;
	}

	void enqueUnfinished()
	{
		for ( uint64_t i = 0; i < CDLV.size(); ++i )
		{
			if ( CDLV[i].missingdep == 0 )
			{
				std::cerr << "[V] container " << i << " has no missing dependencies, enqueuing jobs" << std::endl;

				for ( uint64_t j = 0; j < VCC[i].V.size(); ++j )
				{
					addUnfinished(JobDescription(i,j));
				}
			}
		}
	}

	std::vector < StartWorkerRequest > computeStartRequests(libmaus2::util::ArgParser const & arg)
	{
		std::vector < StartWorkerRequest > Vreq(workers);
		for ( uint64_t i = 0; i < workers; ++i )
			Vreq[i] = StartWorkerRequest(
				nextworkerid,tmpfilebase,hostname,serverport,
				workertime,workermem,workerthreads,partition,arg,AW.begin(),i,
				idToSlot,workers
			);
		return Vreq;
	}

	void checkRequeue(uint64_t const slotid)
	{
		Mfail [ AW[slotid].packageid ] += 1;

		libmaus2::util::CommandContainer & CC = VCC[AW[slotid].packageid.containerid];
		libmaus2::util::Command & CO = CC.V[AW[slotid].packageid.subid];

		// mark pipeline as failed
		if ( Mfail [ AW[slotid].packageid ] >= CC.maxattempt )
		{
			std::cerr << "[V] too many failures on " << AW[slotid].packageid.containerid << "," << AW[slotid].packageid.subid << ", marking pipeline as failed" << std::endl;

			if ( !CO.ignorefail )
				failed = true;
		}
		// requeue
		else
		{
			std::cerr << "[V] requeuing " << AW[slotid].packageid.containerid << "," << AW[slotid].packageid.subid << std::endl;

			addUnfinished(AW[slotid].packageid);
			processWakeupSet();
			processResubmitSet();
		}
	}


	void resetSlot(uint64_t const slotid)
	{
		uint64_t const id = AW[slotid].id;
		EP.remove(AW[slotid].Asocket->getFD());
		fdToSlot.erase(AW[slotid].Asocket->getFD());
		AW[slotid].reset();
		idToSlot.erase(id);
		wakeupSet.erase(id);
		restartSet.insert(slotid);
	}

	void writeContainer(uint64_t const i)
	{
		std::string const fn = CDLV.at(i).fn;
		std::string const tmpfn = fn + ".tmp";

		libmaus2::aio::OutputStreamInstance::unique_ptr_type pOSI(
			new libmaus2::aio::OutputStreamInstance(tmpfn)
		);

		VCC . at ( i ) . serialise ( *pOSI );

		pOSI->flush();
		pOSI.reset();

		libmaus2::aio::OutputStreamFactoryContainer::rename(tmpfn,fn);
	}

	void handleSuccessfulCommand(
		uint64_t const slotid,
		bool const verbose = false
	)
	{
		if ( verbose )
			std::cerr << "[V] getting package id for slot " << slotid << std::endl;
		JobDescription const packageid = AW[slotid].packageid;
		if ( verbose )
			std::cerr << "[V] found package id " << packageid.containerid << "," << packageid.subid << std::endl;

		if ( verbose )
			std::cerr << "[V] getting command container" << std::endl;
		libmaus2::util::CommandContainer & CC = VCC[packageid.containerid];
		if ( verbose )
			std::cerr << "[V] got command container" << std::endl;

		libmaus2::util::Command & CO = CC.V[packageid.subid];

		if ( verbose )
			std::cerr << "[V] updating numattempts,completed" << std::endl;
		CO.numattempts += 1;
		CO.completed = true;
		if ( CO.deepsleep )
		{
			assert ( ndeepsleep > 0 );
			ndeepsleep -= 1;
		}
		Srunning.erase(packageid);
		if ( verbose )
			std::cerr << "[V] updated numattempts,completed to " << CO.numattempts << "," << CO.completed << std::endl;

		writeContainer(packageid.containerid);

		uint64_t const numunfin = --Munfinished [ packageid.containerid ];

		if ( !numunfin )
		{
			std::cerr << "[V] finished command container " << AW[slotid].packageid.containerid << std::endl;

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
						addUnfinished(JobDescription(k,j));
					}
					if ( Sunfinished.size() )
					{
						processWakeupSet();
						processResubmitSet();
					}
				}
			}
		}

		AW[slotid].resetPackageId();
	}

	void handleFailedCommand(uint64_t const slotid)
	{
		std::cerr << "[V] getting package id for slot " << slotid << std::endl;
		JobDescription const packageid = AW[slotid].packageid;
		std::cerr << "[V] found package id " << packageid.containerid << "," << packageid.subid << std::endl;

		std::cerr << "[V] getting reference to command container" << std::endl;
		libmaus2::util::CommandContainer & CC = VCC.at(packageid.containerid);
		std::cerr << "[V] got reference to command container" << std::endl;

		std::cerr << "[V] gettting reference to command" << std::endl;
		libmaus2::util::Command & CO = CC.V[packageid.subid];
		std::cerr << "[V] got reference to command" << std::endl;

		std::cerr << "[V] incrementing numattempts" << std::endl;
		CO.numattempts += 1;
		if ( CO.deepsleep )
		{
			assert ( ndeepsleep > 0 );
			ndeepsleep -= 1;
		}
		Srunning.erase(packageid);
		std::cerr << "[V] incremented numattempts to " << CO.numattempts << std::endl;

		if ( CO.numattempts >= CC.maxattempt && CO.ignorefail )
		{
			std::cerr << "[V] number of attempts reached max " << CC.maxattempt << " but container has ignorefail flag set" << std::endl;

			std::cerr << "[V] decreasing numattempts" << std::endl;
			CO.numattempts -= 1;
			if ( CO.deepsleep )
				ndeepsleep += 1;
			Srunning.insert(packageid);
			std::cerr << "[V] decreased numattempts to " << CO.numattempts << std::endl;

			std::cerr << "[V] calling handleSuccesfulCommand" << std::endl;
			handleSuccessfulCommand(slotid,true);
			std::cerr << "[V] returned from handleSuccesfulCommand" << std::endl;
		}
		else
		{
			writeContainer(packageid.containerid);
			checkRequeue(slotid);
			AW[slotid].resetPackageId();
		}
	}


	SlurmControl(
		std::string const & rtmpfilebase,
		uint64_t const rworkertime,
		uint64_t const rworkermem,
		std::string const rpartition,
		uint64_t const rworkers,
		std::string const & rcdl,
		int64_t const rworkerthreads,
		libmaus2::util::ArgParser const & rarg
	)
	: curdir(libmaus2::util::ArgInfo::getCurDir()),
	  serverport(50000), backlog(1024), tries(1000), nextworkerid(0), hostname(libmaus2::network::GetHostName::getHostName()),
	  tmpfilebase(rtmpfilebase),
	  workertime(rworkertime),
	  workermem(rworkermem),
	  partition(rpartition),
	  workers(rworkers),
	  cdl(rcdl),
	  AW(workers),
	  idToSlot(),
	  fdToSlot(),
	  Mfail(),
	  CDL(loadCDL(cdl)),
	  CDLV(CDL.V),
	  VCC(loadVCC(CDLV)),
	  Sunfinished(),
	  Srunning(),
	  ndeepsleep(0),
	  Munfinished(),
	  maxthreads(computeMaxThreads()),
	  workerthreads(rworkerthreads > 0 ? rworkerthreads : maxthreads),
	  EP(workers+1),
	  Pservsock(
		libmaus2::network::ServerSocket::allocateServerSocket(
			serverport,
			backlog,
			hostname,
			tries
		)
	  ),
	  restartSet(),
	  wakeupSet(),
	  // pending(0),
	  pstate(),
	  failed(false),
	  Vreq(computeStartRequests(rarg))
	{
		std::cerr << "[V] hostname=" << hostname << " serverport=" << serverport << " number of containers " << CDLV.size() << std::endl;

		EP.add(Pservsock->getFD());
		fdToSlot[Pservsock->getFD()] = std::numeric_limits<uint64_t>::max();

		countUnfinished();
		enqueUnfinished();
	}

	int process()
	{
		for ( uint64_t i = 0; i < workers; ++i )
			Vreq[i].dispatch();

		while ( Sunfinished.size() || Srunning.size() )
		{
			std::set<uint64_t> nrestartSet;
			for ( std::set<uint64_t>::const_iterator it = restartSet.begin(); it != restartSet.end(); ++it )
			{
				uint64_t const i = *it;

				try
				{
					Vreq[i].dispatch();
				}
				catch(std::exception const & ex)
				{
					std::cerr << "[E] job start failed:\n" << ex.what() << std::endl;
					nrestartSet.insert(i);
					AW[i].reset();
				}
			}
			restartSet.clear();
			restartSet = nrestartSet;

			ProgState npstate(
				Sunfinished.size(),
				Srunning.size()
			);

			if ( npstate != pstate )
			{
				pstate = npstate;
				std::cerr << "[V] Sunfinished.size()=" << pstate.numunfinished << " pending=" << pstate.numpending << std::endl;
			}

			int rfd = -1;
			if ( EP.wait(rfd) )
			{
				std::map<int,uint64_t>::const_iterator itslot = fdToSlot.find(rfd);

				if ( itslot == fdToSlot.end() )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] EPoll::wait returned unknown file descriptor" << std::endl;
					lme.finish();
					throw lme;
				}
				else if ( itslot->second == std::numeric_limits<uint64_t>::max() )
				{
					assert ( rfd == Pservsock->getFD() );
					int64_t slot = -1;

					try
					{
						libmaus2::network::SocketBase::unique_ptr_type nptr = Pservsock->accept();

						FDIO fdio(nptr->getFD());
						uint64_t const jobid = fdio.readNumber();

						std::cerr << "[V] accepted connection for jobid=" << jobid << std::endl;

						if ( idToSlot.find(jobid) != idToSlot.end() )
						{
							slot = idToSlot.find(jobid)->second;
							fdio.writeNumber(AW[slot].workerid);
							fdio.writeString(curdir);
							bool const curdirok = fdio.readNumber();

							if ( curdirok )
							{
								std::ostringstream tmpostr;
								tmpostr << tmpfilebase << "_container_" << AW[slot].workerid;
								fdio.writeString(tmpostr.str());

								if ( ! AW[slot].Asocket )
								{
									AW[slot].Asocket = UNIQUE_PTR_MOVE(nptr);
									EP.add(AW[slot].Asocket->getFD());
									fdToSlot[AW[slot].Asocket->getFD()] = slot;
									AW[slot].active = true;

									std::cerr << "[V] marked slot " << slot << " active for jobid " << AW[slot].id << std::endl;
								}
								else
								{
									libmaus2::exception::LibMausException lme;
									lme.getStream() << "[E] erratic worker trying to open second connection" << std::endl;
									lme.finish();
									throw lme;
								}
							}
						}
						else
						{
							std::cerr << "[V] job id unknown" << std::endl;
						}
					}
					catch(std::exception const & ex)
					{
						std::cerr << "[E] error while accepting new connection:\n" << ex.what() << std::endl;
						if ( slot >= 0 )
							AW[slot].reset();
					}
				}
				else
				{
					uint64_t const i = itslot->second;

					std::cerr << "[V] epoll returned slot " << i << " ready for reading" << std::endl;

					if ( ! AW[i].active )
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] epoll returned file descriptor for inactive slot" << std::endl;
						lme.finish();
						throw lme;
					}

					try
					{
						// read worker state
						FDIO fdio(AW[i].Asocket->getFD());
						uint64_t const rd = fdio.readNumber();

						// worker is idle
						if ( rd == 0 )
						{
							if ( Sunfinished.size() )
							{
								// get next package
								JobDescription const currentid = getUnfinished();
								// get command
								libmaus2::util::Command const com = VCC[currentid.containerid].V[currentid.subid];
								// serialise command to string
								std::ostringstream ostr;
								com.serialise(ostr);
								// process command
								AW[i].packageid = currentid;
								fdio.writeNumber(0);
								fdio.writeString(ostr.str());
								fdio.writeNumber(currentid.containerid);
								fdio.writeNumber(currentid.subid);

								Srunning.insert(currentid);
								if ( com.deepsleep )
									ndeepsleep += 1;

								std::cerr << "[V] started " << com << " for " << currentid.containerid << "," << currentid.subid << " on slot " << i << std::endl;
							}
							else
							{
								if ( ndeepsleep == Srunning.size() )
								{
									// put slot to deep sleep
									std::cerr << "[V] putting slot " << i << " to deep sleep" << std::endl;

									// request termination
									fdio.writeNumber(2);
									EP.remove(AW[i].Asocket->getFD());
									fdToSlot.erase(AW[i].Asocket->getFD());
									AW[i].reset();
									Sresubmit.insert(i);
								}
								else
								{
									std::cerr << "[V] putting slot " << i << " in wakeupSet" << std::endl;

									wakeupSet.insert(i);
								}
							}
						}
						// worker has finished a job (may or may not be succesful)
						else if ( rd == 1 )
						{
							uint64_t const status = fdio.readNumber();
							int const istatus = static_cast<int>(status);
							// acknowledge
							fdio.writeNumber(0);

							std::cerr << "[V] slot " << i << " reports job ended with istatus=" << istatus << std::endl;

							if ( WIFEXITED(istatus) && (WEXITSTATUS(istatus) == 0) )
							{
								handleSuccessfulCommand(i);
							}
							else
							{
								std::cerr << "[V] slot " << i << " failed, checking requeue " << AW[i].packageid.containerid << "," << AW[i].packageid.subid << std::endl;
								handleFailedCommand(i);
							}
						}
						// worker is still running a job
						else if ( rd == 2 )
						{
							// acknowledge
							fdio.writeNumber(0);
						}
						else
						{
							std::cerr << "[V] process for slot " << i << " jobid " << AW[i].id << " is erratic" << std::endl;

							if ( AW[i].packageid.containerid >= 0 )
							{
								handleFailedCommand(i);
							}

							resetSlot(i /* slotid */);
						}
					}
					catch(std::exception const & ex)
					{
						std::cerr << "[V] exception for slot " << i << " jobid " << AW[i].id << std::endl;
						std::cerr << ex.what() << std::endl;

						if ( AW[i].packageid.containerid >= 0 )
						{
							try
							{
								handleFailedCommand(i);
							}
							catch(std::exception const & ex)
							{
								std::cerr << "[E] exception in handleFailedCommand: " << std::endl;
								std::cerr << ex.what() << std::endl;
								throw;
							}
						}

						try
						{
							resetSlot(i /* slotid */);
						}
						catch(std::exception const & ex)
						{
							std::cerr << "[E] exception in resetSlot: " << std::endl;
							std::cerr << ex.what() << std::endl;
							throw;
						}
					}
				}
			}
		}

		processWakeupSet();
		processResubmitSet();

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
					// read worker state
					FDIO fdio(AW[i].Asocket->getFD());
					uint64_t const rd = fdio.readNumber();

					// worker is idle
					if ( rd == 0 )
					{
						// request terminate
						fdio.writeNumber(2);
						EP.remove(AW[i].Asocket->getFD());
						fdToSlot.erase(AW[i].Asocket->getFD());
						AW[i].reset();
						Vterm.push_back(i);
					}
					else if ( rd == 1 )
					{
						std::cerr << "[V] slot " << i << " reports finished job with no jobs active" << std::endl;
					}
					// worker is still running a job
					else if ( rd == 2 )
					{
						std::cerr << "[V] slot " << i << " reports job running, but we know of no such job" << std::endl;
						fdio.writeNumber(0);
					}
					else
					{
						std::cerr << "[V] process for slot " << i << " jobid " << AW[i].id << " is erratic" << std::endl;

						resetSlot(i /* slotid */);
					}
				}
				catch(...)
				{
					std::cerr << "[V] exception for slot " << i << " jobid " << AW[i].id << std::endl;

					resetSlot(i /* slotid */);
				}
			}

			for ( uint64_t i = 0; i < Vterm.size(); ++i )
				Sactive.erase(Vterm[i]);
		}

		if ( failed )
		{
			std::cerr << "[E] pipeline failed" << std::endl;
			return EXIT_FAILURE;
		}
		else
		{
			std::cerr << "[V] pipeline finished ok" << std::endl;
			return EXIT_SUCCESS;
		}
	}
};

int slurmcontrol(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	uint64_t const workertime = arg.uniqueArgPresent("workertime") ? arg.getParsedArg<uint64_t>("workertime") : 1440;
	uint64_t const workermem = arg.uniqueArgPresent("workermem") ? arg.getParsedArg<uint64_t>("workermem") : 40000;
	std::string const partition = arg.uniqueArgPresent("p") ? arg["p"] : "haswell";
	uint64_t const workers = arg.uniqueArgPresent("workers") ? arg.getParsedArg<uint64_t>("workers") : 16;

	std::string const cdl = arg[0];


	SlurmControl SC(
		tmpfilebase,workertime,workermem,partition,workers,cdl,
		arg.uniqueArgPresent("workerthreads") ? arg.getParsedArg<uint64_t>("workerthreads") : -1,
		arg
	);

	int const r = SC.process();

	return r;
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

		int const r = slurmcontrol(arg);

		return r;
	}
	catch(std::exception const & ex)
	{
		std::cerr << "[E] exception in main: " << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
