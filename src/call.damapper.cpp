/*
    lastools
    Copyright (C) 2016 German Tischler

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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/GetFileSize.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/util/WriteableString.hpp>
#include <libmaus2/util/stringFunctions.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/aio/OutputStreamInstance.hpp>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <config.h>

#include <libgen.h>

static std::string sdirname(std::string const & path)
{
	libmaus2::util::WriteableString W(path);
	return std::string(::dirname(W.A.begin()));
}

static std::string sbasename(std::string const & path)
{
	libmaus2::util::WriteableString W(path);
	return std::string(::basename(W.A.begin()));
}

static std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <ref.fasta> <reads.fasta" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n\n";
	ostr << arg.progname << " passes the arguments for the parameters k,t,M,e,s,n,B,T,b,v and p through to damapper" << std::endl;
	ostr << "\nFurther parameters\n";
	ostr << "--refblocksize: reference block size used for DBsplit (--refblocksize250 by default)" << std::endl;
	ostr << "--readblocksize: read block size used for DBsplit (--refblocksize250 by default)" << std::endl;
	ostr << "-I<dir>: directory to store the dazzler reference index (directory containing reference FastA by default)" << std::endl;
	ostr << "-W<dir>: work directory (directory for storing temporary files, by default the current directory)" << std::endl;

	return ostr.str();
}

static bool which(std::string const prog, std::string & res)
{
	char const * path = getenv("PATH");

	if ( path )
	{
		std::string const spath(path);
		std::deque<std::string> pathtokens = libmaus2::util::stringFunctions::tokenize(spath,std::string(":"));

		for ( uint64_t i = 0; i < pathtokens.size(); ++i )
		{
			std::string path = pathtokens[i];

			if ( path.size() && path[path.size()-1] == '/' )
			{
				uint64_t j = path.size();
				while ( j && path[j-1] == '/' )
					--j;

				if ( ! j )
					++j;

				path = path.substr(0,j);
			}

			res = path + "/" + prog;
			if ( libmaus2::util::GetFileSize::fileExists(res) )
				return true;
		}

		return false;
	}
	else
	{
		res = "./" + prog;
		return libmaus2::util::GetFileSize::fileExists(res);
	}
}

static std::string getCurrentDir()
{
	uint64_t const l = 2 * PATH_MAX;
	libmaus2::autoarray::AutoArray<char> A(l,false);
	char * p = getcwd(A.begin(),A.size());

	if ( p == A.begin() )
		return std::string(A.begin());
	else
	{
		int const error = errno;
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: unable to get current directory " << strerror(error) << std::endl;
		lme.finish();
		throw lme;
	}
}

static std::string makeAbsolute(std::string const & s)
{
	if ( ! s.size() )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: makeAbsolute: string is empty" << std::endl;
		lme.finish();
		throw lme;
	}

	if ( s[0] == '/' )
		return s;
	else
		return getCurrentDir() + "/" + s;
}

#include <stack>
#include <sys/stat.h>
#include <sys/types.h>

static void mkdirP(std::string dirname)
{
	std::stack<std::string> Vfn;

	Vfn.push(dirname);

	do
	{
		dirname = sdirname(dirname);
		Vfn.push(dirname);
	} while ( dirname != sdirname(dirname) );

	while ( !Vfn.empty() )
	{
		std::string const dir = Vfn.top();
		Vfn.pop();
		if ( mkdir(dir.c_str(),0700) )
		{
			int const error = errno;

			switch ( error )
			{
				case EEXIST:
					break;
				default:
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "call.damapper: mkdir(" << dir << "): " << strerror(error) << std::endl;
					lme.finish();
					throw lme;
					break;
			}
		}
	}
}

struct ProgramDescriptor
{
	std::string program;
	std::string path;
	std::string dir;

	static std::string dirname(std::string const & s)
	{
		libmaus2::util::WriteableString W(s);
		char * d = ::dirname(W.A.begin());
		return std::string(d);
	}

	ProgramDescriptor() {}
	ProgramDescriptor(std::string const & rprogram, std::string const & rcall)
	: program(rprogram)
	{
		assert ( program.size() );
		assert ( program[0] != '/' );

		std::string call = makeAbsolute(rcall);

		if ( ! which(program,path) )
		{
			path = sdirname(call) + "/" + program;

			if ( ! libmaus2::util::GetFileSize::fileExists(path) )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "call.damapper: unable to find program " << program << std::endl;
				lme.finish();
				throw lme;
			}
		}

		path = makeAbsolute(path);
		dir = dirname(path);
	}
};

static void callrename(std::string const & from, std::string const & to)
{
	if ( ::rename(from.c_str(),to.c_str()) != 0 )
	{
		int const error = errno;
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: rename(" << from << "," << to << ") failed: " << strerror(error) << std::endl;
		lme.finish();
		throw lme;
	}
}

static bool startsWith(std::string const & s, std::string const & t)
{
	return
		s.size() >= t.size() &&
		std::equal(
			s.begin(),
			s.begin()+t.size(),
			t.begin()
		);
}


int callsystem(std::string const & s)
{
	pid_t const pid = fork();

	// fork error
	if ( pid == static_cast<pid_t>(-1) )
	{
		int const error = errno;
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: fork() failed: " << strerror(error) << std::endl;
		lme.finish();
		throw lme;
	}
	// child
	else if ( pid == 0 )
	{
		// close standard output
		::close(STDOUT_FILENO);
		// redirect standard output to standard error
		dup2(STDERR_FILENO,STDOUT_FILENO);
		int const r = system(s.c_str());
		_exit(r);
	}
	// parent
	else
	{
		int status = -1;
		waitpid(pid,&status,0);
		return status;
	}
}

static int call_damapper(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	unsigned int const damapper_arg_k = arg.uniqueArgPresent("k") ? arg.getUnsignedNumericArg<uint64_t>("k") : 20;
	         int const damapper_arg_t = arg.uniqueArgPresent("t") ? arg.getParsedArg<int>("t") : -1;
	         int const damapper_arg_M = arg.uniqueArgPresent("M") ? arg.getParsedArg<int>("M") : -1;
	double const damapper_arg_e = arg.uniqueArgPresent("e") ? arg.getParsedArg<double>("e") : std::numeric_limits<double>::min();
	         int const damapper_arg_s = arg.uniqueArgPresent("s") ? arg.getParsedArg<int>("s") : -1;
	double const damapper_arg_n = arg.uniqueArgPresent("n") ? arg.getParsedArg<double>("n") : std::numeric_limits<double>::min();
	         int const damapper_arg_B = arg.uniqueArgPresent("B") ? arg.getParsedArg<int>("B") : -1;
	         int const damapper_arg_T = arg.uniqueArgPresent("T") ? arg.getParsedArg<int>("T") : -1;
	         int const damapper_arg_b = arg.uniqueArgPresent("b") ? 1 :  0; // comp bias
	         int const damapper_arg_v = arg.uniqueArgPresent("v") ? 1 :  0; // verbose
	         int const damapper_arg_p = arg.uniqueArgPresent("p") ? 1 :  0; // repeat profile (not used)

	unsigned int const numthreads = damapper_arg_T < 0 ? 4 : damapper_arg_T;

	unsigned int const refblocksize = arg.uniqueArgPresent("refblocksize") ? arg.getUnsignedNumericArg<uint64_t>("refblocksize") : 250;
	unsigned int const readblocksize = arg.uniqueArgPresent("readblocksize") ? arg.getUnsignedNumericArg<uint64_t>("readblocksize") : 250;

	std::string const argindexdir = arg.uniqueArgPresent("I") ? arg["I"] :  std::string();
	std::string const argworkdir = arg.uniqueArgPresent("W") ? arg["W"] :  std::string();
	std::string const workdir = makeAbsolute(argworkdir.size() ? argworkdir : ".");

	mkdirP(workdir);

	unsigned int const fastacolwidth = 80;
	bool const reformatref = true;
	bool const reformatreads = true;

	std::string progname = arg.progname;
	ProgramDescriptor PD_fasta2DAM("fasta2DAM",progname);
	ProgramDescriptor PD_DBsplit("DBsplit",progname);
	ProgramDescriptor PD_HPC_damapper("HPC.damapper",progname);
	ProgramDescriptor PD_damapper("damapper",progname);
	ProgramDescriptor PD_lascat("lascat",progname);
	ProgramDescriptor PD_lastobam("lastobam",progname);
	ProgramDescriptor PD_LAsort("LAsort",progname);
	ProgramDescriptor PD_LAcat("LAcat",progname);
	ProgramDescriptor PD_LAmerge("LAmerge",progname);

	std::vector<std::string> Vpathadd;
	Vpathadd.push_back(PD_LAsort.dir);
	Vpathadd.push_back(PD_LAcat.dir);
	Vpathadd.push_back(PD_LAmerge.dir);
	std::sort(Vpathadd.begin(),Vpathadd.end());
	Vpathadd.resize ( std::unique(Vpathadd.begin(),Vpathadd.end()) - Vpathadd.begin() );

	std::ostringstream npathstr;
	npathstr << getenv("PATH");
	for ( uint64_t i = 0; i < Vpathadd.size(); ++i )
		npathstr << ":" << Vpathadd[i];
	std::string const npath = npathstr.str();
	libmaus2::util::WriteableString Wnpath(npath);
	setenv("PATH",Wnpath.A.begin(),1 /* overwrite */);

	std::string const reffn = arg[0];

	if ( ! libmaus2::util::GetFileSize::fileExists(reffn) )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: reference file " << reffn << " does not exist" << std::endl;
		lme.finish();
		throw lme;
	}

	static char const * faend[] = { ".fasta",".fa",0 };

	std::string const refdir = sdirname(reffn);
	std::string const indexdir = makeAbsolute(argindexdir.size() ? argindexdir : refdir);
	std::string const reffile = sbasename(reffn);
	std::string const reffileclipped = libmaus2::util::OutputFileNameTools::endClip(reffile,&faend[0]);

	mkdirP(indexdir);

	std::string const refdam = indexdir + "/"  + reffileclipped + ".dam";
	std::string const refidx = indexdir + "/." + reffileclipped + ".idx";
	std::string const refhdr = indexdir + "/." + reffileclipped + ".hdr";
	std::string const refbps = indexdir + "/." + reffileclipped + ".bps";
	std::string const indexfasta = indexdir + "/"  + reffileclipped + ".dam.fasta";

	if (
		! libmaus2::util::GetFileSize::fileExists(refdam) || libmaus2::util::GetFileSize::isOlder(refdam,reffn) ||
		! libmaus2::util::GetFileSize::fileExists(refidx) || libmaus2::util::GetFileSize::isOlder(refidx,reffn) ||
		! libmaus2::util::GetFileSize::fileExists(refhdr) || libmaus2::util::GetFileSize::isOlder(refhdr,reffn) ||
		! libmaus2::util::GetFileSize::fileExists(refbps) || libmaus2::util::GetFileSize::isOlder(refbps,reffn) ||
		! libmaus2::util::GetFileSize::fileExists(indexfasta) || libmaus2::util::GetFileSize::isOlder(indexfasta,reffn)
	)
	{
		std::string const tmpbase = arginfo.getDefaultTmpFileName();
		std::string const reffastatmp  = indexdir + "/"  + tmpbase + ".fasta";
		std::string const refdamtmp    = indexdir + "/"  + tmpbase + ".dam";
		std::string const refdamidxtmp = indexdir + "/." + tmpbase + ".idx";
		std::string const refdamhdrtmp = indexdir + "/." + tmpbase + ".hdr";
		std::string const refdambpstmp = indexdir + "/." + tmpbase + ".bps";

		libmaus2::util::TempFileRemovalContainer::addTempFile(refdamtmp);
		libmaus2::util::TempFileRemovalContainer::addTempFile(refdamidxtmp);
		libmaus2::util::TempFileRemovalContainer::addTempFile(refdamhdrtmp);
		libmaus2::util::TempFileRemovalContainer::addTempFile(refdambpstmp);
		libmaus2::util::TempFileRemovalContainer::addTempFile(reffastatmp);

		if ( reformatref )
		{
			libmaus2::fastx::FastAReader FAreader(reffn);
			libmaus2::fastx::FastAReader::pattern_type pat;
			libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(reffastatmp));
			while ( FAreader.getNextPatternUnlocked(pat) )
				pat.printMultiLine(*OSI,fastacolwidth);
			OSI->flush();
			OSI.reset();
		}
		else
		{
			if ( symlink(reffn.c_str(),reffastatmp.c_str()) )
			{
				int const error = errno;

				libmaus2::exception::LibMausException lme;
				lme.getStream() << "call.damapper: symlink(" << reffn << "," << reffastatmp << "): " << strerror(error) << std::endl;
				lme.finish();
				throw lme;
			}
		}

		std::string const fasta2DAMcom = PD_fasta2DAM.path + " " + refdamtmp + " " + reffastatmp;

		int const fasta2DAMcomres = callsystem(fasta2DAMcom.c_str());

		if ( damapper_arg_v )
			std::cerr << "[V] " << fasta2DAMcom << std::endl;

		if ( fasta2DAMcomres != 0 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "call.damapper: " << fasta2DAMcom << " failed" << std::endl;
			lme.finish();
			throw lme;
		}

		std::ostringstream splitstr;
		splitstr << PD_DBsplit.path << " -s" << refblocksize << " -x" << damapper_arg_k << " " << refdamtmp;
		std::string const splitcom = splitstr.str();

		int const splitcomres = callsystem(splitcom.c_str());

		if ( damapper_arg_v )
			std::cerr << "[V] " << splitcom << std::endl;

		if ( splitcomres != 0 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "call.damapper: " << splitcomres << " failed" << std::endl;
			lme.finish();
			throw lme;
		}

		callrename(refdamtmp,refdam);
		callrename(refdamidxtmp,refidx);
		callrename(refdamhdrtmp,refhdr);
		callrename(refdambpstmp,refbps);
		callrename(reffastatmp,indexfasta);
	}

	std::string const tmpbase = arginfo.getDefaultTmpFileName();
	std::string const readsfastatmp = workdir + "/" + tmpbase + "_reads" + ".fasta";
	std::string const readsdamtmp = workdir + "/" + tmpbase + "_reads" + ".dam";
	std::string const readsidxtmp = workdir + "/." + tmpbase + "_reads" + ".idx";
	std::string const readsbpstmp = workdir + "/." + tmpbase + "_reads" + ".bps";
	std::string const readshdrtmp = workdir + "/." + tmpbase + "_reads" + ".hdr";
	std::string const readsdamtmpscript = workdir + "/" + tmpbase + "_reads" + ".dam.script";
	std::string const readslastmp = workdir + "/" + tmpbase + "_reads" + ".las";

	libmaus2::util::TempFileRemovalContainer::addTempFile(readsfastatmp);
	libmaus2::util::TempFileRemovalContainer::addTempFile(readsdamtmp);
	libmaus2::util::TempFileRemovalContainer::addTempFile(readsidxtmp);
	libmaus2::util::TempFileRemovalContainer::addTempFile(readsbpstmp);
	libmaus2::util::TempFileRemovalContainer::addTempFile(readshdrtmp);
	libmaus2::util::TempFileRemovalContainer::addTempFile(readsdamtmpscript);
	libmaus2::util::TempFileRemovalContainer::addTempFile(readslastmp);

	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(readsfastatmp));
	if ( reformatreads )
	{
		libmaus2::fastx::StreamFastAReaderWrapper SFAR(std::cin);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pat;
		while ( SFAR.getNextPatternUnlocked(pat) )
			pat.printMultiLine(*OSI,fastacolwidth);

	}
	else
	{
		libmaus2::autoarray::AutoArray<char> A(4096,false);
		while ( std::cin )
		{
			std::cin.read(A.begin(),A.size());
			uint64_t const r = std::cin.gcount();
			OSI->write(A.begin(),r);
		}
	}
	OSI->flush();
	OSI.reset();

	std::string const fasta2DAMcom = PD_fasta2DAM.path + " " + readsdamtmp + " " + readsfastatmp;

	int const fasta2DAMcomres = callsystem(fasta2DAMcom.c_str());

	if ( damapper_arg_v )
		std::cerr << "[V] " << fasta2DAMcom << std::endl;

	if ( fasta2DAMcomres != 0 )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: " << fasta2DAMcom << " failed" << std::endl;
		lme.finish();
		throw lme;
	}

	std::ostringstream splitstr;
	splitstr << PD_DBsplit.path << " -s" << readblocksize << " -x" << damapper_arg_k << " " << readsdamtmp;
	std::string const splitcom = splitstr.str();

	int const splitcomres = callsystem(splitcom.c_str());

	if ( damapper_arg_v )
		std::cerr << "[V] " << splitcom << std::endl;

	if ( splitcomres != 0 )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: " << splitcomres << " failed" << std::endl;
		lme.finish();
		throw lme;
	}

	std::string const calldir = getCurrentDir();
	if ( chdir(workdir.c_str()) )
	{
		int const error = errno;
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: chdir(" << workdir << "): " << strerror(error) << std::endl;
		lme.finish();
		throw lme;
	}

	std::ostringstream HPC1str;
	HPC1str << PD_HPC_damapper.path;

	HPC1str << " -k" << damapper_arg_k;
	if ( damapper_arg_t != -1 )
		HPC1str << " -t" << damapper_arg_t;
	if ( damapper_arg_M != -1 )
		HPC1str << " -M" << damapper_arg_M;
	if ( damapper_arg_e != std::numeric_limits<double>::min() )
		HPC1str << " -e" << damapper_arg_e;
	if ( damapper_arg_s != -1 )
		HPC1str << " -s" << damapper_arg_s;
	if ( damapper_arg_n != std::numeric_limits<double>::min() )
		HPC1str << " -n" << damapper_arg_n;
	if ( damapper_arg_B != -1 )
		HPC1str << " -B" << damapper_arg_B;
	HPC1str << " -T" << numthreads;
	if ( damapper_arg_b )
		HPC1str << " -b";
	if ( damapper_arg_v )
		HPC1str << " -v";
	if ( damapper_arg_p )
		HPC1str << " -p";

	HPC1str << " -C -N " << refdam << " " << readsdamtmp << " >" << readsdamtmpscript;
	std::string const HPC1 = HPC1str.str();

	int const HPC1res = system(HPC1.c_str());

	if ( damapper_arg_v )
		std::cerr << "[V] " << HPC1 << std::endl;

	if ( HPC1res != 0 )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: " << HPC1 << " failed" << std::endl;
		lme.finish();
		throw lme;
	}

	std::vector<std::string> Voutput;
	std::vector<std::string> Vdamapper;
	libmaus2::aio::InputStreamInstance::unique_ptr_type HPCisi(new libmaus2::aio::InputStreamInstance(readsdamtmpscript));
	std::string const exp = "# Catenate jobs (2)";
	std::string const expdamapper = "# Damapper jobs (1)";
	std::string const expsub = "LAsort -a ";
	std::string const expcheck = "# Check all .las files (optional but recommended)";
	std::string const expchecksub = "LAcheck -vS ";

	bool activedamapper = false;
	bool activeoutput = false;
	bool activationfound = false;
	bool activationdamapperfound = false;

	while ( *HPCisi )
	{
		std::string line;
		std::getline(*HPCisi,line);

		// std::cerr << line << std::endl;

		if ( line.size() )
		{
			if ( line[0] == '#' )
			{
				activeoutput = false;
				activedamapper = false;

				if ( line == expcheck )
				{
					activeoutput = true;
					activationfound = true;
				}
				else if ( line == expdamapper )
				{
					activedamapper = true;
					activationdamapperfound = true;
				}
			}
			else if ( activeoutput )
			{
				if (
					line.size() >= expchecksub.size()
					&&
					line.substr(0,expchecksub.size()) == expchecksub
				)
				{
					line = line.substr(expchecksub.size());

					std::vector<std::string> tokens;
					uint64_t l = 0;
					while ( l < line.size() )
					{
						while ( l < line.size() && ::std::isspace(static_cast<unsigned char>(line[l])) )
							++l;

						uint64_t h = l;
						while ( h < line.size() && ! ::std::isspace(static_cast<unsigned char>(line[h])) )
							++h;

						if ( h > l )
						{
							tokens.push_back(line.substr(l,h-l));
							// std::cerr << "token " << tokens.back() << std::endl;
						}

						l = h;
					}

					if ( tokens.size() == 3 )
					{
						Voutput.push_back(tokens[2] + ".las");
					}
					else
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "call.damapper: cannot understand line " << line << std::endl;
						lme.finish();
						throw lme;
					}
				}
				else
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "call.damapper: cannot understand line " << line << std::endl;
					lme.finish();
					throw lme;
				}
			}
			else if ( activedamapper )
			{
				if ( startsWith(line,"damapper ") )
				{
					line = line.substr(std::string("damapper ").size());
					line = PD_damapper.path + " " + line;
					Vdamapper.push_back(line);
				}
				else
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "call.damapper: cannot understand line " << line << std::endl;
					lme.finish();
					throw lme;
				}
			}
		}
	}
	HPCisi.reset();

	if ( ! activationfound )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: activation line for cat jobs not found" << std::endl;
		lme.finish();
		throw lme;
	}
	if ( ! activationdamapperfound )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: activation line for damapper not found" << std::endl;
		lme.finish();
		throw lme;
	}

	for ( uint64_t i = 0; i < Voutput.size(); ++i )
		libmaus2::util::TempFileRemovalContainer::addTempFile(Voutput[i]);

	/* call damapper */
	for ( uint64_t i = 0; i < Vdamapper.size(); ++i )
	{
		int const damapperret = callsystem(Vdamapper[i].c_str());

		if ( damapper_arg_v )
			std::cerr << "[V] " << Vdamapper[i] << std::endl;

		if ( damapperret != 0 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "call.damapper: damapper call failed" << std::endl;
			lme.finish();
			throw lme;
		}
	}

	/* concatenate all output .las files */
	std::ostringstream catostr;
	catostr << PD_lascat.path << " " << readslastmp;
	for ( uint64_t i = 0; i < Voutput.size(); ++i )
	{
		catostr << " " << Voutput[i];
	}
	std::string const catcom = catostr.str();

	if ( damapper_arg_v )
		std::cerr << "[V] " << catcom << std::endl;

	if ( system(catcom.c_str()) != 0 )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: lascat call " << catcom << " failed" << std::endl;
		lme.finish();
		throw lme;
	}

	for ( uint64_t i = 0; i < Voutput.size(); ++i )
		::remove(Voutput[i].c_str());

	/* convert .las file to BAM */
	std::ostringstream lastobamstr;
	lastobamstr << PD_lastobam.path << " -t" << numthreads << " -snone " << refdam << " " << indexfasta << " " << readsdamtmp << " " << readsfastatmp << " " << readslastmp;
	std::string const lastobamcom = lastobamstr.str();

	if ( damapper_arg_v )
		std::cerr << "[V] " << lastobamcom << std::endl;

	if ( system(lastobamcom.c_str()) != 0 )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: lastobam call " << lastobamcom << " failed" << std::endl;
		lme.finish();
		throw lme;
	}

	if ( chdir(calldir.c_str()) )
	{
		int const error = errno;
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "call.damapper: chdir(" << calldir << "): " << strerror(error) << std::endl;
		lme.finish();
		throw lme;
	}


	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
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

		return call_damapper(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
