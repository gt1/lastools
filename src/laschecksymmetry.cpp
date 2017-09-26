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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus2/aio/SerialisedPeeker.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <in.db> <in1.las> ..." << std::endl;

	return ostr.str();
}


int laschecksymmetry(libmaus2::util::ArgParser & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	std::string const symforfn = tmpfilebase + "ovl.for.sym";
	std::string const syminvfn = tmpfilebase + "ovl.inv.sym";

	libmaus2::util::TempFileRemovalContainer::addTempFile(symforfn);
	libmaus2::util::TempFileRemovalContainer::addTempFile(syminvfn);

	libmaus2::aio::OutputStreamInstance::unique_ptr_type Psymfor(
		new libmaus2::aio::OutputStreamInstance(symforfn)
	);
	libmaus2::aio::OutputStreamInstance::unique_ptr_type Psyminv(
		new libmaus2::aio::OutputStreamInstance(syminvfn)
	);

	std::string const dbfn = arg[0];

	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	DB.computeTrimVector();

	std::vector<uint64_t> RL;
	DB.getAllReadLengths(RL);

	for ( uint64_t i = 1; i < arg.size(); ++i )
	{
		std::string const lasfn = arg[i];
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type afile(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(lasfn));

		libmaus2::dazzler::align::Overlap OVL;
		uint64_t c = 0;
		while ( afile->getNextOverlap(OVL) )
		{
			libmaus2::dazzler::align::OverlapInfo info = OVL.getHeader().getInfo().swapped();

			if ( info.aread & 1 )
				info = info.inverse(RL.begin());

			info.serialise(*Psyminv);
			OVL.getHeader().getInfo().serialise(*Psymfor);

			if ( (++c % (1024*1024)) == 0 )
				std::cerr << "[V] " << c << std::endl;
		}
	}

	Psyminv->flush();
	Psyminv.reset();
	Psymfor->flush();
	Psymfor.reset();

	libmaus2::sorting::SerialisingSortingBufferedOutputFile<libmaus2::dazzler::align::OverlapInfo>::sort(symforfn,16*1024*1024);
	libmaus2::sorting::SerialisingSortingBufferedOutputFile<libmaus2::dazzler::align::OverlapInfo>::sort(syminvfn,16*1024*1024);

	libmaus2::aio::SerialisedPeeker<libmaus2::dazzler::align::OverlapInfo> peekerfor(symforfn);
	libmaus2::aio::SerialisedPeeker<libmaus2::dazzler::align::OverlapInfo> peekerinv(syminvfn);
	bool ok = true;

	libmaus2::dazzler::align::OverlapInfo infofor;
	libmaus2::dazzler::align::OverlapInfo infoinv;
	while ( peekerinv.peekNext(infoinv) && peekerfor.peekNext(infofor) )
	{
		if ( infoinv < infofor )
		{
			std::cerr << "[E] missing " << infoinv.swapped() << std::endl;
			peekerinv.getNext(infoinv);
			ok = false;
		}
		else if ( infofor < infoinv )
		{
			std::cerr << "[E] missing " << infofor.swapped() << std::endl;
			peekerfor.getNext(infofor);
			ok = false;
		}
		else
		{
			peekerinv.getNext(infoinv);
			peekerfor.getNext(infofor);
			// std::cerr << "ok " << infofor << std::endl;
		}
	}

	if ( ok )
	{
		std::cout << "ok" << std::endl;
		return EXIT_SUCCESS;
	}
	else
	{
		std::cout << "not ok" << std::endl;
		return EXIT_FAILURE;
	}
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);

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
		else if ( arg.size() < 2 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return laschecksymmetry(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
