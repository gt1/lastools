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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [-a] <in_a.db> <in_b.db> <in.las>" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

int64_t getDefaultTSpace()
{
	return libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
}

int lasshow(libmaus2::util::ArgParser const & arg)
{
	bool const showalignment = arg.uniqueArgPresent("a");
	int64_t const vaid = arg.uniqueArgPresent("aid") ? arg.getUnsignedNumericArg<uint64_t>("aid") : -1;
	int64_t const vbid = arg.uniqueArgPresent("bid") ? arg.getUnsignedNumericArg<uint64_t>("bid") : -1;

	std::string const db0name = arg[0];
	std::string const db1name = arg[1];

	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type uDB0;
	libmaus2::dazzler::db::DatabaseFile * pDB0 = 0;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type uDB1;
	libmaus2::dazzler::db::DatabaseFile * pDB1 = 0;

	{
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type tDB0(new libmaus2::dazzler::db::DatabaseFile(db0name));
		uDB0 = UNIQUE_PTR_MOVE(tDB0);
		pDB0 = uDB0.get();
	}

	if ( db1name == db0name )
	{
		pDB1 = uDB0.get();
	}
	else
	{
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type tDB1(new libmaus2::dazzler::db::DatabaseFile(db1name));
		uDB1 = UNIQUE_PTR_MOVE(tDB1);
		pDB1 = uDB1.get();
	}

	libmaus2::dazzler::db::DatabaseFile & DB0 = *pDB0;
	DB0.computeTrimVector();
	std::vector<uint64_t> RL0;
	DB0.getAllReadLengths(RL0);

	libmaus2::dazzler::db::DatabaseFile & DB1 = *pDB1;
	DB1.computeTrimVector();
	std::vector<uint64_t> RL1;
	DB1.getAllReadLengths(RL1);

	libmaus2::lcs::NP NP;
	libmaus2::lcs::AlignmentTraceContainer ATC;

	libmaus2::dazzler::align::Overlap OVL;
	for ( uint64_t i = 2; i < arg.size(); ++i )
	{
		int64_t const intspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(arg[i]);
		// int64_t const novl = libmaus2::dazzler::align::AlignmentFile::getNovl(arg[i]);
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(arg[i]));

		std::string a;
		int64_t aid = -1;
		std::string b;
		std::string br;
		int64_t bid = -1;

		while ( PIN->getNextOverlap(OVL) )
		{
			if ( OVL.aread != aid )
			{
				aid = OVL.aread;
				a = DB0[aid];

				//std::cerr << "[V] aid=" << aid << " length=" << a.size() << std::endl;
			}
			if ( OVL.bread != bid )
			{
				bid = OVL.bread;
				b = DB1[bid];
				br = libmaus2::fastx::reverseComplementUnmapped(b);
			}

			if (
				(vaid < 0 || vaid == OVL.aread)
				&&
				(vbid < 0 || vbid == OVL.bread)
			)
			{
				std::cout
					<< aid << "[" << OVL.path.abpos << "," << OVL.path.aepos << ")/" << a.size()
					<< " "
					<< bid
					<< (OVL.isInverse()?'c':'n')
					<< "[" << OVL.path.bbpos << "," << OVL.path.bepos << ")/" << b.size() << " " << OVL.path.diffs << " diffs " << OVL.getErrorSum() << " esum"
					<< " e " << (static_cast<double>(OVL.getErrorSum()) / (OVL.path.aepos-OVL.path.abpos))
					;

				if ( OVL.isTrue() )
				{
					std::cout << " " << "true";
				}
				if ( OVL.isHaploFlag() )
				{
					std::cout << " " << "haplo";
				}

				if ( showalignment )
				{

					uint8_t const * aptr = reinterpret_cast<uint8_t const *>(a.c_str());
					uint8_t const * bptr =
						OVL.isInverse()
						?
						reinterpret_cast<uint8_t const *>(br.c_str())
						:
						reinterpret_cast<uint8_t const *>( b.c_str());

					OVL.computeTrace(aptr,bptr,intspace,ATC,NP);

					std::cout << " " << ATC.getAlignmentStatistics() << "\n";
					std::cout.put('\n');

					libmaus2::lcs::AlignmentPrint::printAlignmentLines(
						std::cout,
						aptr + OVL.path.abpos, OVL.path.aepos-OVL.path.abpos,
						bptr + OVL.path.bbpos, OVL.path.bepos-OVL.path.bbpos,
						80,
						ATC.ta,
						ATC.te
					);
				}
				else
				{
					std::cout.put('\n');
				}
			}
		}
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

		return lasshow(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
