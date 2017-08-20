/*
    libmaus2
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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/GetFileSize.hpp>
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/util/stringFunctions.hpp>

#include <unistd.h>

std::string getSelf()
{
	std::ostringstream ostr;

	ostr << "/proc/self/exe";

	if ( libmaus2::util::GetFileSize::fileExists(ostr.str()) )
	{
		libmaus2::autoarray::AutoArray<char> A(PATH_MAX+1);
		std::fill(A.begin(),A.end(),0);

		int const r = readlink(ostr.str().c_str(),A.begin(),A.size()-1);

		if ( r < 0 )
		{
			int const error = errno;
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] " << strerror(error) << std::endl;
			lme.finish();
			throw lme;
		}

		return std::string(A.begin(),A.begin()+r);
	}

	return std::string();
}

std::string relConfigPath()
{
	std::string self = getSelf();

	if ( self.find('/') != std::string::npos )
		self = self.substr(0,self.find_last_of('/'));

	if ( self.size() )
		return self + "/../share/derive.rc";
	else
		return std::string();
}

libmaus2::util::ArgParser parseArgs(std::string const & line)
{
	std::vector<std::string> V;
	uint64_t cl = 0;
	while ( cl < line.size() )
	{
		while ( cl < line.size() && isspace(line[cl]) )
			++cl;
		uint64_t ch = cl;
		while ( ch < line.size() && !isspace(line[ch]) )
			++ch;
		if ( ch != cl )
		{
			V.push_back(line.substr(cl,ch-cl));
		}

		cl = ch;
	}

	return libmaus2::util::ArgParser(V);
}

struct ConfigEntry
{
	enum config_entry_type
	{
		config_entry_type_thread,
		config_entry_type_mem
	};

	bool mandatory;
	std::string arg;
	config_entry_type type;
	int64_t defaultvalue;
	uint64_t multiplier;
	uint64_t add;

	ConfigEntry() {}
	ConfigEntry(bool const rmandatory, std::string const & rarg, config_entry_type rtype, int64_t rdefaultvalue, uint64_t rmultiplier, uint64_t radd)
	: mandatory(rmandatory), arg(rarg), type(rtype), defaultvalue(rdefaultvalue), multiplier(rmultiplier), add(radd)
	{

	}
};

std::map < std::string, std::vector<ConfigEntry> > loadConfig(libmaus2::util::ArgParser const & arg)
{
	std::string configfile;
	if (
		arg.uniqueArgPresent("c")
	)
	{
		if ( libmaus2::util::GetFileSize::fileExists(arg["c"]) )
		{
			configfile = arg["c"];
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] file " << arg["c"] << " does not exist" << std::endl;
			lme.finish();
			throw lme;
		}
	}
	else if ( libmaus2::util::GetFileSize::fileExists(std::string(getenv("HOME"))+"/.lascommand/derive.rc") )
	{
		configfile = std::string(getenv("HOME"))+"/.lascommand/derive.rc";
	}
	else if ( libmaus2::util::GetFileSize::fileExists(relConfigPath()) )
	{
		configfile = relConfigPath();
	}
	else if ( libmaus2::util::GetFileSize::fileExists(std::string(DATA_PATH) + "/" + "derive.rc") )
	{
		configfile = std::string(DATA_PATH) + "/" + "derive.rc";
	}

	if ( ! configfile.size() )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] no config file found" << std::endl;
		lme.finish();
		throw lme;
	}

	std::cerr << "[V] using config file " << configfile << std::endl;

	std::map < std::string, std::vector<ConfigEntry> > CEM;

	libmaus2::aio::InputStreamInstance ISI(configfile);

	while ( ISI )
	{
		std::string line;
		std::getline(ISI,line);

		if ( line.size() )
		{
			if ( line[0] == '#' )
			{
			}
			else
			{
				std::deque<std::string> tokens = libmaus2::util::stringFunctions::tokenize(line,std::string("\t"));

				// prog type flag
				if ( tokens.size() >= 6 )
				{
					std::string const prog = tokens[0];

					ConfigEntry::config_entry_type type;
					if ( tokens[1] == "thread" )
						type = ConfigEntry::config_entry_type_thread;
					else if ( tokens[1] == "mem" )
						type = ConfigEntry::config_entry_type_mem;
					else
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] unknown resource type " << tokens[1] << std::endl;
						lme.finish();
						throw lme;
					}

					bool mandatory;
					if ( tokens[2] == "mandatory" )
						mandatory = true;
					else if ( tokens[2] == "optional" )
						mandatory = false;
					else
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] unknown flag type " << tokens[2] << " (use either mandatory or optional)" << std::endl;
						lme.finish();
						throw lme;
					}

					std::istringstream multistr(tokens[4]);
					uint64_t multiplier;
					multistr >> multiplier;
					if ( ! ( multistr && multistr.peek() == std::istream::traits_type::eof() ) )
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] unable to parse multiplier " << tokens[4] << std::endl;
						lme.finish();
						throw lme;
					}

					std::istringstream addistr(tokens[5]);
					uint64_t add;
					addistr >> add;
					if ( ! ( addistr && addistr.peek() == std::istream::traits_type::eof() ) )
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] unable to parse add " << tokens[5] << std::endl;
						lme.finish();
						throw lme;
					}

					int64_t defaultvalue = -1;
					if ( 6 < tokens.size() )
					{
						defaultvalue = atol(tokens[4].c_str());
					}

					ConfigEntry CE(mandatory,tokens[3],type,defaultvalue,multiplier,add);

					CEM [ prog ] . push_back(CE);
				}
				else
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] insufficient token count in config line " << line << std::endl;
					lme.finish();
					throw lme;
				}
			}
		}
	}

	return CEM;
}

int commandderive(libmaus2::util::ArgParser const & arg)
{
	std::map < std::string, std::vector < ConfigEntry > > const CEM = loadConfig(arg);

	while ( std::cin )
	{
		std::string line;

		std::getline(std::cin,line);

		if ( line.size() )
		{
			if ( line[0] == '#' )
				std::cout << line << std::endl;
			else
			{
				uint64_t ch = 0;
				while ( ch < line.size() && !isspace(line[ch]) )
					++ch;

				std::string const command = line.substr(0,ch);

				libmaus2::util::ArgParser const subarg = parseArgs(line);

				uint64_t mem = 1000;
				uint64_t threads = 1;

				if ( CEM.find(command) != CEM.end() )
				{
					std::vector<ConfigEntry> const & CEV = CEM.find(command)->second;

					for ( uint64_t i = 0; i < CEV.size(); ++i )
					{
						ConfigEntry const & CE = CEV[i];

						if ( CE.mandatory && !subarg.uniqueArgPresent(CE.arg) )
						{
							libmaus2::exception::LibMausException lme;
							lme.getStream() << "[E] need " << CE.arg << " arg for " << command <<  " calls" << std::endl;
							lme.finish();
							throw lme;
						}

						if ( subarg.uniqueArgPresent(CE.arg) )
						{
							uint64_t const v = subarg.getUnsignedNumericArg<uint64_t>(CE.arg);

							if ( CE.type == ConfigEntry::config_entry_type_mem )
								mem = v * CE.multiplier + CE.add;
							else if ( CE.type == ConfigEntry::config_entry_type_thread )
								threads = v * CE.multiplier + CE.add;
						}
						else
						{
							if ( CE.defaultvalue >= 0 )
							{
								if ( CE.type == ConfigEntry::config_entry_type_mem )
									mem = CE.defaultvalue * CE.multiplier + CE.add;
								else if ( CE.type == ConfigEntry::config_entry_type_thread )
									threads = CE.defaultvalue * CE.multiplier + CE.add;
							}
						}
					}
				}

				std::cout << threads << "\t" << mem << "\t" << line << std::endl;
			}
		}
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);

		return commandderive(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
