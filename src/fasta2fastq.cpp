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

/**
 * convert FastA input to FastQ adding dummy quality value 'H' for each base
 **/
int main()
{
	try
	{
		libmaus2::fastx::StreamFastAReaderWrapper SFARW(std::cin);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
		libmaus2::autoarray::AutoArray<char> H(0,false);
		while ( SFARW.getNextPatternUnlocked(pattern) )
		{
			if ( pattern.spattern.size() > H.size() )
			{
				H.resize(pattern.spattern.size());
				std::fill(H.begin(),H.end(),'H');
			}

			for ( std::string::size_type i = 0; i < pattern.spattern.size(); ++i )
				pattern.spattern[i] = ::std::toupper(pattern.spattern[i]);

			std::cout << "@" << pattern.sid << "\n" << pattern.spattern << "\n" << "+\n";
			std::cout.write(H.begin(),pattern.spattern.size());
			std::cout.put('\n');
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
