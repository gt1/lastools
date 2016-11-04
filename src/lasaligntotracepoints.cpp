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
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <out.las> <in.las>" << std::endl;
	ostr << "\n";

	return ostr.str();
}

int lasaligntotracepoints(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & /* arginfo */)
{
	libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pfile(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(arg[1]));
	int64_t const tspace = Pfile->Palgn->tspace;
	libmaus2::dazzler::align::Overlap OVL;
	libmaus2::dazzler::align::AlignmentWriter AW(arg[0],tspace,true);

	while ( Pfile->getNextOverlap(OVL) )
	{
		OVL.alignToTracePoints(tspace);
		if ( ! OVL.isEmpty() )
			AW.put(OVL);
	}

	return EXIT_SUCCESS;
}

/**
 *
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
		else if ( arg.size() < 2 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return lasaligntotracepoints(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
