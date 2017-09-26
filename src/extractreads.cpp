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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/util/ArgParser.hpp>

int extractreads(libmaus2::util::ArgParser const & arg)
{
	bool const trim = !arg.uniqueArgPresent("notrim");
	bool const rc = arg.uniqueArgPresent("rc");

	std::string const dbname = arg[0];
	libmaus2::dazzler::db::DatabaseFile DB(dbname);
	if ( trim )
		DB.computeTrimVector();
	libmaus2::aio::InputStream::unique_ptr_type Pbasestr(DB.openBaseStream());
	libmaus2::autoarray::AutoArray<char> A;

	std::vector < std::string > args;
	for ( uint64_t i = 1; i < arg.size(); ++i )
		args.push_back(arg[i]);
	if ( ! args.size() )
	{
		std::ostringstream ostr;
		ostr << 0 << "-" << DB.size();
		args.push_back(ostr.str());
	}

	for ( uint64_t a = 0; a < args.size(); ++a )
	{
		std::string const arg = args[a];
		std::pair<uint64_t,uint64_t> P;

		if ( arg.find('-') != std::string::npos )
		{
			std::string const front = arg.substr(0,arg.find('-'));
			std::string const back  = arg.substr(arg.find('-')+1);
			std::istringstream frontistr(front);
			std::istringstream backistr(back);

			frontistr >> P.first;
			backistr >> P.second;

			if ( frontistr.bad() || frontistr.fail() || frontistr.peek() != std::istream::traits_type::eof() )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "usage: cannot parse " << front << std::endl;
				lme.finish();
				throw lme;
			}
			if ( backistr.bad() || backistr.fail() || backistr.peek() != std::istream::traits_type::eof() )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "usage: cannot parse " << back << std::endl;
				lme.finish();
				throw lme;
			}
		}
		else
		{
			std::string const front = arg.substr(0,arg.find('-'));
			std::istringstream frontistr(front);

			frontistr >> P.first;

			if ( frontistr.bad() || frontistr.fail() || frontistr.peek() != std::istream::traits_type::eof() )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "usage: cannot parse " << front << std::endl;
				lme.finish();
				throw lme;
			}

			P.second = P.first+1;
		}

		std::vector<libmaus2::dazzler::db::Read> VR;
		std::cerr << "[V] extracting read interval [" << P.first << "," << P.second << ")" << std::endl;
		DB.getReadInterval(P.first,P.second,VR);

		for ( uint64_t i = P.first; i < P.second; ++i )
		{
			uint64_t const untrimmedid = trim ? DB.trimmedToUntrimmed(i) : i;
			libmaus2::dazzler::db::Read const & R = VR[i-P.first];
			Pbasestr->seekg(R.boff,std::ios::beg);
			DB.decodeRead(*Pbasestr,A,R.rlen);

			if ( rc )
			{
				std::reverse(A.begin(),A.begin()+R.rlen);
				for ( int64_t i = 0; i < R.rlen; ++i )
					A[i] = libmaus2::fastx::invertUnmapped(A[i]);
			}

			std::string const name = DB.getReadName(untrimmedid,R);

			std::cout << ">" << name << std::endl;

			char const * p = A.begin();
			uint64_t t = R.rlen;

			while ( t )
			{
				uint64_t const maxlen = 80;
				uint64_t const towrite = std::min(maxlen,t);
				std::cout.write(p,towrite);
				p += towrite;
				t -= towrite;
				std::cout.put('\n');
			}
		}
	}

	return EXIT_SUCCESS;
}

/*
 * extract a set of reads from a database and writes them as FastA files
 *
 * use trim=1 to extract from the trimmed database
 */

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		if ( arg.size() < 1 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "usage: extractreads <db> <readid> ..." << std::endl;
			lme.finish();
			throw lme;
		}

		return extractreads(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
