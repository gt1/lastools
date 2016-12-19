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

void handle(
	std::vector < libmaus2::dazzler::align::Overlap > & VOVL, 
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> & TPV, 
	int64_t const tspace,
	libmaus2::dazzler::align::AlignmentWriter & AW
)
{
	std::sort(VOVL.begin(),VOVL.end(),libmaus2::dazzler::align::OverlapComparator());
		
	std::vector<bool> killvec(VOVL.size(),false);
	
	uint64_t l = 0;
	while ( l < VOVL.size() )
	{
		uint64_t h = l+1;
		while ( h < VOVL.size() && VOVL[h].isInverse() == VOVL[l].isInverse() )
			++h;
			
		uint64_t o = 0;
		for ( uint64_t i = l ; i < h; ++i )
			o = VOVL[i].getTracePoints(tspace,i,TPV,o);
			
		std::sort(TPV.begin(),TPV.begin()+o);
		
		uint64_t tl = 0;
		while ( tl != o )
		{
			uint64_t th = tl+1;
			while ( th != o && TPV[th].apos == TPV[tl].apos && TPV[th].bpos == TPV[tl].bpos )
				++th;
				
			if ( th-tl > 1 )
			{
				for ( uint64_t i = tl; i < th-1; ++i )
				{
					int64_t const aid = TPV[i].id;
					libmaus2::dazzler::align::Overlap const & A = VOVL[aid];
					libmaus2::math::IntegerInterval<int64_t> IA(A.path.abpos,A.path.aepos-1);
					
					for ( uint64_t j = i+1; j < th; ++j )
					{
						int64_t const bid = TPV[j].id;
						libmaus2::dazzler::align::Overlap const & B = VOVL[bid];
						libmaus2::math::IntegerInterval<int64_t> IB(B.path.abpos,B.path.aepos-1);
						
						if ( IA == IB )
						{
							assert ( aid < bid );
							killvec[bid] = true;
						}
						else
						{
							libmaus2::math::IntegerInterval<int64_t> IC = IA.intersection(IB);
							
							if ( IB.diameter() < IA.diameter() && IC.diameter() >= 0.95 * IB.diameter() )
								killvec[bid] = true;
							else if ( IA.diameter() < IB.diameter() && IC.diameter() >= 0.95 * IA.diameter() )
								killvec[aid] = true;
						}
					}
				}
			}
				
			tl = th;
		}
			
		l = h;
	}
	
	uint64_t o = 0;
	for ( uint64_t i = 0; i < VOVL.size(); ++i )
		if ( ! killvec[i] )
			VOVL[o++] = VOVL[i];
	VOVL.resize(o);
	for ( uint64_t i = 0; i < VOVL.size(); ++i )
		AW.put(VOVL[i]);
	
	VOVL.resize(0);
}

int lasdedup(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & /* arginfo */)
{
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPV;

	for ( uint64_t i = 0; i < arg.size(); ++i )
	{
		std::string const infn = arg[i];
		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(infn);
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));
		libmaus2::dazzler::align::Overlap OVL;
		libmaus2::dazzler::align::AlignmentWriter AW(libmaus2::util::OutputFileNameTools::clipOff(infn,".las")+"_dedup.las",tspace);
		int64_t prevaread = std::numeric_limits<int64_t>::max();
		int64_t prevbread = std::numeric_limits<int64_t>::max();
		std::vector < libmaus2::dazzler::align::Overlap > VOVL;

		while ( Plas->getNextOverlap(OVL) )
		{
			if ( OVL.aread != prevaread || OVL.bread != prevbread )
				handle(VOVL,TPV,tspace,AW);

			prevaread = OVL.aread;
			prevbread = OVL.bread;
			VOVL.push_back(OVL);
		}

		handle(VOVL,TPV,tspace,AW);
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
			return lasdedup(arg,arginfo);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
