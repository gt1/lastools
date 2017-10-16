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

int lassimplededup(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPV;
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : arginfo.getDefaultTmpFileName();

	for ( uint64_t i = 0; i < arg.size(); ++i )
	{
		std::string const infn = arg[i];

		std::cerr << "[V] sorting " << infn << "...";
		// libmaus2::dazzler::align::SortingOverlapOutputBuffer<libmaus2::dazzler::align::OverlapFullComparator>::sort(infn,tmpfilebase+".tmp");
		std::cerr << "done." << std::endl;

		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(infn);
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));
		libmaus2::dazzler::align::AlignmentWriter AW(libmaus2::util::OutputFileNameTools::clipOff(infn,".las")+"_simplededup.las",tspace);

		libmaus2::dazzler::align::OverlapFullComparator comp;
		libmaus2::dazzler::align::Overlap OVLprev;
		libmaus2::dazzler::align::Overlap OVL;
		bool prevvalid = false;

		while ( Plas->getNextOverlap(OVL) )
		{
			if ( prevvalid )
			{
				bool eq =
					!comp(OVLprev,OVL)
					&&
					!comp(OVL,OVLprev);

				assert (
					comp(OVLprev,OVL)
					||
					eq
				);
			}

			if ( !prevvalid || comp(OVLprev,OVL) )
				AW.put(OVL);

			OVLprev = OVL;
			prevvalid = true;
		}
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
			std::cerr << "usage: " << arg.progname << " [options] in.las\n";
			return EXIT_SUCCESS;
		}
		else
		{
			libmaus2::util::ArgInfo const arginfo(argc,argv);
			return lassimplededup(arg,arginfo);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
