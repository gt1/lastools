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

int lascomputetandem(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & /* arginfo */)
{
	bool const onlytrue = arg.uniqueArgPresent("onlytrue");

	for ( uint64_t i = 0; i < arg.size(); ++i )
	{
		std::string const infn = arg[i];
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));
		libmaus2::dazzler::align::Overlap OVL;
		int64_t prevaread = std::numeric_limits<int64_t>::max();

		std::vector < libmaus2::math::IntegerInterval<int64_t> > VI;

		while ( Plas->getNextOverlap(OVL) )
		{
			assert ( OVL.aread == OVL.bread );

			if ( OVL.aread != prevaread )
			{
				VI = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(VI);
				for ( uint64_t i = 0; i < VI.size(); ++i )
				{
					libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,prevaread);
					libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,VI[i].from);
					libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,VI[i].to);
					// std::cerr << prevaread << "\t" << VI[i] << "\n";
				}
				VI.resize(0);
			}

			if ( OVL.isInverse() )
				continue;

			libmaus2::math::IntegerInterval<int64_t> const IA(OVL.path.abpos,OVL.path.aepos-1);
			libmaus2::math::IntegerInterval<int64_t> const IB(OVL.path.bbpos,OVL.path.bepos-1);

			if ( ! IA.intersection(IB).isEmpty() )
			{
				libmaus2::math::IntegerInterval<int64_t> const IS = libmaus2::math::IntegerInterval<int64_t>::span(IA,IB);
				VI.push_back(IS);
			}
			else if ( ! onlytrue )
			{
				VI.push_back(IA);
				VI.push_back(IB);
			}

			prevaread = OVL.aread;
		}

		VI = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(VI);
		for ( uint64_t i = 0; i < VI.size(); ++i )
		{
			libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,prevaread);
			libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,VI[i].from);
			libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,VI[i].to);
			//std::cerr << prevaread << "\t" << VI[i] << "\n";
		}
		VI.resize(0);
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
			return lascomputetandem(arg,arginfo);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
