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

#include <config.h>

#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus2/dazzler/align/OverlapInfoIndexer.hpp>

void handle(
	std::vector < libmaus2::dazzler::align::Overlap > & VOVL,
	libmaus2::dazzler::align::AlignmentWriter & AW
)
{
	for ( uint64_t i = 1; i < VOVL.size(); ++i )
		assert (
			VOVL[i].aread == VOVL[0].aread
			&&
			VOVL[i].bread == VOVL[0].bread
			&&
			VOVL[i].isInverse() == VOVL[0].isInverse()
		);

	std::sort(VOVL.begin(),VOVL.end(),libmaus2::dazzler::align::OverlapFullComparator());

	if ( VOVL.size() )
	{
		uint64_t o = 1;
		libmaus2::dazzler::align::OverlapFullComparator comp;

		for ( uint64_t i = 1; i < VOVL.size(); ++i )
			if ( comp(VOVL[i-1],VOVL[i]) )
			{
				VOVL[o++] = VOVL[i];
			}
			else
			{
				#if 0
				std::cerr << "[V] removing exact copy" << std::endl;
				std::cerr << VOVL[i-1].getHeader().getInfo() << std::endl;
				std::cerr << VOVL[i-0].getHeader().getInfo() << std::endl;
				#endif
			}

		VOVL.resize(o);
	}

	for ( uint64_t i = 0; i < VOVL.size(); ++i )
		AW.put(VOVL[i]);

	VOVL.resize(0);
}

int lasdedupexact(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPV;
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : arginfo.getDefaultTmpFileName();
	std::string const symkillfn = tmpfilebase + "_symkill.tmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(symkillfn);
	libmaus2::aio::OutputStreamInstance::unique_ptr_type Psymkill(new libmaus2::aio::OutputStreamInstance(symkillfn));

	std::string const outfn = arg[0];
	std::string const dbfn = arg[1];

	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	DB.computeTrimVector();
	std::vector<uint64_t> RL;
	DB.getAllReadLengths(RL);


	std::vector<std::string> Vin;
	for ( uint64_t i = 2; i < arg.size(); ++i )
	{
		std::string const infn = arg[i];
		Vin.push_back(infn);
	}
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vin);

	libmaus2::dazzler::align::AlignmentWriter AW(outfn,tspace);

	for ( uint64_t i = 0; i < Vin.size(); ++i )
	{
		std::string const infn = Vin[i];

		std::cerr << "[V] sorting " << infn << "...";
		libmaus2::dazzler::align::SortingOverlapOutputBuffer<libmaus2::dazzler::align::OverlapFullComparator>::sort(infn,tmpfilebase+".tmp");
		std::cerr << "done." << std::endl;

		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));
		libmaus2::dazzler::align::Overlap OVL;
		int64_t prevaread = std::numeric_limits<int64_t>::max();
		int64_t prevbread = std::numeric_limits<int64_t>::max();
		bool previnverse = false;
		std::vector < libmaus2::dazzler::align::Overlap > VOVL;

		while ( Plas->getNextOverlap(OVL) )
		{
			if ( OVL.aread != prevaread || OVL.bread != prevbread || OVL.isInverse() != previnverse )
				handle(VOVL,AW);

			prevaread = OVL.aread;
			prevbread = OVL.bread;
			previnverse = OVL.isInverse();
			VOVL.push_back(OVL);
		}

		handle(VOVL,AW);
	}

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 1 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] out.las in.db in.las\n\n";
			std::cerr << "options:" << std::endl;
			std::cerr << " -T : prefix for temporary files (default: create files in current working directory)\n";
			return EXIT_SUCCESS;
		}
		else
		{
			libmaus2::util::ArgInfo const arginfo(argc,argv);
			return lasdedupexact(arg,arginfo);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
