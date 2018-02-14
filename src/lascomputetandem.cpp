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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/bambam/BamNumericalIndexDecoder.hpp>
#include <libmaus2/bambam/BamNumericalIndexGenerator.hpp>
#include <libmaus2/bambam/BamAccessor.hpp>

int64_t imax(int64_t const a, int64_t const b)
{
	return std::max(a,b);
}

void handle(
	std::vector < libmaus2::math::IntegerInterval<int64_t> > & VI,
	std::vector < libmaus2::math::IntegerInterval<int64_t> > & VT,
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type & Pdb,
	int64_t const prevaread,
	int const verbose,
	std::string const & bamfn
)
{
	// if not empty
	if ( VI.size() )
	{
		// merge intervals
		VI = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(VI);

		// read length
		int64_t const rlen =
			VI.size() ?
			(
				Pdb ? Pdb->getRead(prevaread).rlen : -1
			)
			:
			-1;

		for ( uint64_t i = 0; i < VI.size(); ++i )
		{
			// read id
			libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,prevaread);
			// interval start
			libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,VI[i].from);
			// interval end
			libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,VI[i].to+1);

			if ( verbose )
			{
				std::cerr << prevaread << "\t" << VI[i];
				if ( rlen >= 0 )
				{
					// read length
					std::cerr << "\t" << rlen;
					// distance from start
					std::cerr << "\t" << imax(0,VI[i].from);
					// distance from end
					std::cerr << "\t" << imax(0,rlen - VI[i].to);
				}

				if ( bamfn.size() )
				{
					// get BAM accessor
					libmaus2::bambam::BamAccessor::unique_ptr_type tptr(new libmaus2::bambam::BamAccessor(bamfn,prevaread));
					// get BAM record
					libmaus2::bambam::BamAlignment const & algn = (*tptr)[prevaread];

					// check size match
					if ( rlen >= 0 )
						assert ( algn.getLseq() == rlen );

					// get reference interval for complete read
					libmaus2::math::IntegerInterval<int64_t> intv = algn.getReferenceInterval();

					// tandem interval on read
					uint64_t ifrom = VI[i].from;
					uint64_t ito = VI[i].to+1;

					if ( algn.isReverse() )
					{
						std::swap(ifrom,ito);
						ifrom = rlen - ifrom;
						ito = rlen - ito;
					}

					// tandem interval on reference
					uint64_t const tfrom = algn.getRefPosForReadPos(ifrom);
					uint64_t const tto = algn.getRefPosForReadPos(ito-1);

					libmaus2::math::IntegerInterval<int64_t> tintv(tfrom,tto);

					std::cerr << "\t" << intv;
					std::cerr << "\t" << tintv;

					VT.push_back(tintv);
				}

				std::cerr << "\n";
			}
		}
		VI.resize(0);
	}

}

int lascomputetandem(libmaus2::util::ArgParser const & arg)
{
	bool const onlytrue = arg.uniqueArgPresent("onlytrue");
	bool const verbose = arg.uniqueArgPresent("verbose");

	std::string const db = arg.uniqueArgPresent("d") ? arg["d"] : std::string();
	std::string const bamfn = arg.uniqueArgPresent("b") ? arg["b"] : std::string();

	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type Pdb;

	if ( db.size() )
	{
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type Tdb(new libmaus2::dazzler::db::DatabaseFile(db));
		Pdb = UNIQUE_PTR_MOVE(Tdb);
	}

	if ( bamfn.size() )
	{
		uint64_t const genindexmod = 128;
		libmaus2::bambam::BamNumericalIndexGenerator::indexFileCheck(bamfn,genindexmod,1 /* numthreads */);
	}

	std::vector < libmaus2::math::IntegerInterval<int64_t> > VT;

	for ( uint64_t i = 0; i < arg.size(); ++i )
	{
		std::string const infn = arg[i];
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));
		libmaus2::dazzler::align::Overlap OVL;
		int64_t prevaread = std::numeric_limits<int64_t>::max();

		std::vector < libmaus2::math::IntegerInterval<int64_t> > VI;

		while ( Plas->getNextOverlap(OVL) )
		{
			// ignore non self overlaps
			if ( OVL.aread != OVL.bread )
				continue;

			// new aread
			if ( OVL.aread != prevaread )
				handle(VI,VT,Pdb,prevaread,verbose,bamfn);

			// skip inverse overlaps
			if ( OVL.isInverse() )
				continue;

			// get A and B intervals
			libmaus2::math::IntegerInterval<int64_t> const IA(OVL.path.abpos,OVL.path.aepos-1);
			libmaus2::math::IntegerInterval<int64_t> const IB(OVL.path.bbpos,OVL.path.bepos-1);

			// if intersecting
			if ( ! IA.intersection(IB).isEmpty() )
			{
				// compute union
				libmaus2::math::IntegerInterval<int64_t> const IS = libmaus2::math::IntegerInterval<int64_t>::span(IA,IB);
				// push union
				VI.push_back(IS);
			}
			// otherwise push both
			else if ( ! onlytrue )
			{
				VI.push_back(IA);
				VI.push_back(IB);
			}

			prevaread = OVL.aread;
		}

		handle(VI,VT,Pdb,prevaread,verbose,bamfn);
	}

	VT = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(VT);
	for ( uint64_t i = 0; i < VT.size(); ++i )
		std::cerr << "VT[" << i << "]=" << VT[i] << std::endl;

	return EXIT_SUCCESS;
}

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [-d<reads.db> -b<reads.bam>] <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -d : Dazzler database for LAS files (default: none)\n";
	ostr << " -b : BAM file storing ground truth alignments for reads (default: none)\n";

	return ostr.str();
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
			std::cerr << getUsage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return lascomputetandem(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
