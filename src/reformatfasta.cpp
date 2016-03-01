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
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/fastx/StreamFastQReader.hpp>
#include <libmaus2/aio/PosixFdInputStream.hpp>
#include <libmaus2/fastx/SpaceTable.hpp>
#include <libmaus2/random/Random.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;
	ostr << "usage: " << arg.progname << " [<parameters>] <in.fasta\n";
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -p : read name prolog (defaults to L)\n";
	ostr << " -i : read id base (defaults to 0)\n";
	ostr << " -c : column width (defaults to 80)\n";
	ostr << " -r : replace non A, C, G, T bases by random bases (0 or 1, defaults to 1)\n";
	return ostr.str();
}

template<typename reader_type>
void fastareformat(libmaus2::util::ArgParser const & arg, reader_type & reader)
{
	typename reader_type::pattern_type pattern;
	std::string const prolog = arg.uniqueArgPresent("p") ? arg["p"] : "L";
	uint64_t readid = arg.uniqueArgPresent("i") ? arg.getUnsignedNumericArg<uint64_t>("i") : 0;
	uint64_t const cols = arg.uniqueArgPresent("c") ? arg.getUnsignedNumericArg<uint64_t>("c") : 80;
	bool const randomN = arg.uniqueArgPresent("r") ? arg.getParsedArg<int>("r") : 1;

	libmaus2::fastx::SpaceTable const ST;
	libmaus2::random::Random::setup();
	bool replaced = false;
	bool replaceprinted = false;

	while ( reader.getNextPatternUnlocked(pattern) )
	{
		std::string s = pattern.spattern;

		if ( randomN )
			for ( uint64_t i = 0; i < s.size(); ++i )
				if ( libmaus2::fastx::mapChar(s[i]) >= 4 )
				{
					s[i] = libmaus2::fastx::remapChar(libmaus2::random::Random::rand8() & 3);
					replaced = true;
				}
		char const * c = s.c_str();

		if ( replaced && ! replaceprinted )
		{
			std::cerr << "[V] warning, replacing non ACGT symbols with random bases" << std::endl;
			replaceprinted = true;
		}

		std::cout << '>' << prolog << '/' << (readid++) << '/' << 0 << '_' << s.size() << " RQ=0.851 " << pattern.sid << "\n";
		uint64_t low = 0;

		while ( low < s.size() )
		{
			uint64_t const high = std::min(low+cols,static_cast<uint64_t>(s.size()));

			std::cout.write(c+low,high-low);
			std::cout.put('\n');

			low = high;
		}
	}
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

		libmaus2::aio::PosixFdInputStream PFIS(STDIN_FILENO,64*1024,1024);

		int const c = PFIS.peek();

		if ( c == '>' || c == std::istream::traits_type::eof() )
		{
			libmaus2::fastx::StreamFastAReaderWrapper in(PFIS);
			fastareformat(arg,in);
		}
		else if ( c == '@' )
		{
			libmaus2::fastx::StreamFastQReaderWrapper in(PFIS);
			fastareformat(arg,in);
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "Unknown input file format, first character " << static_cast<char>(c) << std::endl;
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
