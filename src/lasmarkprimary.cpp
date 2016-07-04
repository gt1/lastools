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
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/geometry/RangeSet.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " <out.las> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
}

struct ReadInterval
{
	uint64_t from;
	uint64_t to;
	uint64_t id;

	ReadInterval() {}
	ReadInterval(uint64_t const rfrom, uint64_t const rto, uint64_t const rid) : from(rfrom), to(rto), id(rid) {}

	uint64_t getFrom() const
	{
		return from;
	}

	uint64_t getTo() const
	{
		return to;
	}
};

void handleVector(std::vector < libmaus2::dazzler::align::Overlap > & VOVL, libmaus2::dazzler::align::AlignmentWriter & AW)
{
	uint64_t maxaepos = 0;
	for ( uint64_t i = 0; i < VOVL.size(); ++i )
		maxaepos = std::max(maxaepos,static_cast<uint64_t>(VOVL[i].path.aepos));

	libmaus2::util::SimpleQueue< libmaus2::geometry::RangeSet<ReadInterval>::search_q_element > SQ;

	libmaus2::geometry::RangeSet<ReadInterval> RS(maxaepos);

	for ( uint64_t i = 0; i < VOVL.size(); ++i )
	{
		uint64_t const abpos = VOVL[i].path.abpos;
		uint64_t const aepos = VOVL[i].path.aepos;
		uint64_t const nf = RS.search(ReadInterval(abpos,aepos,i),SQ);
		bool const primary = (! nf);

		if ( primary )
		{
			VOVL[i].setPrimary();
			std::cerr << "primary " << VOVL[i] << std::endl;
			RS.insert(ReadInterval(abpos,aepos,i));
		}
		else
		{
			std::cerr << "secondary " << VOVL[i] << std::endl;
		}
	}


	for ( uint64_t i = 0; i < VOVL.size(); ++i )
		AW.put(VOVL[i]);
	VOVL.resize(0);
}

int lassort(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const &)
{
	std::string const outfilename = arg[0];
	std::vector<std::string> VI;
	for ( uint64_t i = 1 ; i < arg.size(); ++i )
		VI.push_back(arg[i]);
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(VI);

	libmaus2::dazzler::align::AlignmentWriter::unique_ptr_type AW(
		new libmaus2::dazzler::align::AlignmentWriter(outfilename,tspace,false /* index */, 0 /* expt */)
	);

	libmaus2::dazzler::align::Overlap OVL;
	for ( uint64_t i = 0; i < VI.size(); ++i )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type PIN(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(VI[i]));

		std::vector < libmaus2::dazzler::align::Overlap > VOVL;
		while ( PIN->getNextOverlap(OVL) )
		{
			if ( VOVL.size() && VOVL[0].bread != OVL.bread )
				handleVector(VOVL,*AW);
			VOVL.push_back(OVL);
		}
		handleVector(VOVL,*AW);
	}

	AW.reset();

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
		else if ( arg.size() < 2 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return lassort(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
