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
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/dazzler/align/LasIntervals.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <out.las> <in.db> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int64_t getDefaultTSpace()
{
	return libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
}

int lasprojectsorted(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
{
	std::string const outlas = arg[0];
	
	std::string const dbfn(arg[1]);
	std::vector<std::string> Vin;
	for ( uint64_t i = 2; i < arg.size(); ++i )
		Vin.push_back(arg[i]);
	
	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	DB.computeTrimVector();
	uint64_t const nreads = DB.size();
	
	libmaus2::dazzler::align::LasIntervals LAI(Vin,nreads,std::cerr);
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vin);
	libmaus2::dazzler::align::AlignmentWriter AW(outlas,tspace);
	
	std::vector<std::string> A = arg("a");

	for ( uint64_t i = 0; i < A.size(); ++i )
	{
		int64_t aid = -1;
		int64_t bid = -1;
		
		std::string a = A[i];
		
		if ( a.find(",") != std::string::npos )
		{
			std::string b = a.substr(a.find(",")+1);
			a = a.substr(0,a.find(","));
			
			std::istringstream istr(b);
			istr >> bid;
			
			if ( !istr || istr.peek() != std::istream::traits_type::eof() )
			{
				std::cerr << "[E] cannot parse " << b << ", ignoring " << A[i] << std::endl;
				continue;
			}
		}
		
		{		
			std::istringstream istr(a);
			istr >> aid;
			
			if ( !istr || istr.peek() != std::istream::traits_type::eof() )
			{
				std::cerr << "[E] cannot parse " << a << ", ignoring " << A[i] << std::endl;
				continue;
			}
		}

		libmaus2::dazzler::align::AlignmentFileCat::unique_ptr_type pdec(LAI.openRange(aid,aid+1));
	
		libmaus2::dazzler::align::Overlap OVL;

		while ( pdec->getNextOverlap(OVL) )
			if ( OVL.aread == aid && (bid < 0 || OVL.bread == bid) )
				AW.put(OVL);
	}
	
	return EXIT_SUCCESS;
}

/**
 * sort a set of LAS/OVL files and merge the sorted files to a single output file
 **/

int main(int argc, char *argv[])
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
		else if ( arg.size() < 3 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return lasprojectsorted(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
