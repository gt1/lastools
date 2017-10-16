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

struct Intersection
{
	uint64_t aid;
	uint64_t bid;
	uint64_t apos;
	uint64_t bpos;

	Intersection()
	{

	}

	Intersection(
		uint64_t raid,
		uint64_t rbid,
		uint64_t rapos,
		uint64_t rbpos
	) : aid(raid), bid(rbid), apos(rapos), bpos(rbpos)
	{

	}

	bool operator<(Intersection const & O) const
	{
		if ( aid != O.aid )
			return aid < O.aid;
		else if ( bid != O.bid )
			return bid < O.bid;
		else if ( apos < O.apos )
			return apos < O.apos;
		else
			return bpos < O.bpos;
	}
};

void handle(
	std::vector < libmaus2::dazzler::align::Overlap > & VOVL,
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> & TPV,
	int64_t const tspace,
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


	bool changed = true;

	while ( changed )
	{
		changed = false;

		std::set<uint64_t> killset;
		std::set< Intersection > mergeset;

		// collect trace points with ids
		uint64_t o = 0;
		for ( uint64_t i = 0; i < VOVL.size(); ++i )
			o = VOVL[i].getTracePoints(tspace,i,TPV,o);

		// sort trace points
		std::sort(TPV.begin(),TPV.begin()+o);

		// check for common trace points
		uint64_t tl = 0;
		while ( tl != o )
		{
			// same apos,bpos pair
			uint64_t th = tl+1;
			while ( th != o && TPV[th].apos == TPV[tl].apos && TPV[th].bpos == TPV[tl].bpos )
				++th;

			// more than one?
			if ( th-tl > 1 )
			{
				// check all pairs
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

						// same interval?
						if ( IA == IB )
						{
							assert ( aid < bid );
							killset.insert(bid);
						}
						else
						{
							libmaus2::math::IntegerInterval<int64_t> IC = IA.intersection(IB);

							if ( IC == IA )
							{
								killset.insert(aid);
							}
							else if ( IC == IB )
							{
								killset.insert(bid);
							}
							else
							{
								mergeset.insert(Intersection(aid,bid,TPV[tl].apos,TPV[tl].bpos));
							}
						}
					}
				}
			}

			tl = th;
		}

		std::set<uint64_t> changedset;

		std::vector < libmaus2::dazzler::align::Overlap > Vadd;
		for ( std::set< Intersection >::const_iterator it = mergeset.begin(); it != mergeset.end(); ++it )
			if (
				changedset.find(it->aid) == changedset.end()
				&&
				changedset.find(it->bid) == changedset.end()
			)
			{
				Intersection const P = *it;

				uint64_t const aid = P.aid;
				uint64_t const bid = P.bid;
				int64_t const apos = P.apos;
				int64_t const bpos = P.bpos;
				libmaus2::dazzler::align::Overlap const & A = VOVL[aid];
				libmaus2::dazzler::align::Overlap const & B = VOVL[bid];
				std::vector<libmaus2::dazzler::align::TraceBlock> const TA = A.getTraceBlocks(tspace);
				std::vector<libmaus2::dazzler::align::TraceBlock> const TB = B.getTraceBlocks(tspace);

				assert ( A.path.abpos <= B.path.abpos );

				uint64_t i = 0;
				while ( i < TA.size()
					&&
					(
						TA[i].A.second != apos
						||
						TA[i].B.second != bpos
					)
				)
					++i;

				assert ( i < TA.size() );
				assert ( TA[i].A.second == apos );
				assert ( TA[i].B.second == bpos );

				uint64_t j = 0;
				while ( j < TB.size()
					&&
					(
						TB[j].A.first != apos
						||
						TB[j].B.first != bpos
					)
				)
					++j;
				assert ( j < TB.size() );
				assert ( TB[j].A.first == apos );
				assert ( TB[j].B.first == bpos );


				libmaus2::dazzler::align::Overlap OVL;
				OVL.aread = A.aread;
				OVL.bread = A.bread;
				OVL.flags = A.flags;

				OVL.path.abpos = A.path.abpos;
				OVL.path.bbpos = A.path.bbpos;
				OVL.path.aepos = B.path.aepos;
				OVL.path.bepos = B.path.bepos;

				for ( uint64_t z = 0; z <= i; ++z )
				{
					OVL.path.path.push_back(
						std::pair<uint16_t,uint16_t>(TA[z].err,TA[z].B.second-TA[z].B.first)
					);
				}
				for ( uint64_t z = j; z < TB.size(); ++z )
				{
					OVL.path.path.push_back(
						std::pair<uint16_t,uint16_t>(TB[z].err,TB[z].B.second-TB[z].B.first)
					);
				}

				OVL.path.tlen = 2 * OVL.path.path.size();

				#if 0
				std::cerr << "should merge" << std::endl;
				std::cerr << A << std::endl;
				std::cerr << B << std::endl;
				std::cerr << "apos=" << apos << std::endl;
				std::cerr << "bpos=" << bpos << std::endl;
				std::cerr << "output " << OVL << std::endl;
				#endif

				killset.insert(it->aid);
				killset.insert(it->bid);
				changedset.insert(it->aid);
				changedset.insert(it->bid);
				Vadd.push_back(OVL);
			}

		uint64_t ko = 0;
		for ( uint64_t i = 0; i < VOVL.size(); ++i )
			if ( killset.find(i) == killset.end() )
				VOVL[ko++] = VOVL[i];
		VOVL.resize(ko);

		for ( uint64_t i = 0; i < Vadd.size(); ++i )
			VOVL.push_back(Vadd[i]);

		std::sort(
			VOVL.begin(),VOVL.end(),libmaus2::dazzler::align::OverlapFullComparator()
		);

		changed = (killset.size() != 0);
	}

	if ( VOVL.size() )
		std::cerr << "[V] left " << VOVL[0].aread << "," << VOVL[0].bread << "," << VOVL[0].isInverse() << std::endl;

	for ( uint64_t i = 0; i < VOVL.size(); ++i )
		AW.put(VOVL[i]);

	VOVL.resize(0);
}

int lasdedup(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPV;
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : arginfo.getDefaultTmpFileName();

	for ( uint64_t i = 0; i < arg.size(); ++i )
	{
		std::string const infn = arg[i];

		std::cerr << "[V] sorting " << infn << "...";
		//libmaus2::dazzler::align::SortingOverlapOutputBuffer<libmaus2::dazzler::align::OverlapFullComparator>::sort(infn,tmpfilebase+".tmp");
		std::cerr << "done." << std::endl;

		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(infn);
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));
		libmaus2::dazzler::align::Overlap OVL;
		libmaus2::dazzler::align::AlignmentWriter AW(libmaus2::util::OutputFileNameTools::clipOff(infn,".las")+"_dedup.las",tspace);
		int64_t prevaread = std::numeric_limits<int64_t>::max();
		int64_t prevbread = std::numeric_limits<int64_t>::max();
		bool previnverse = false;
		std::vector < libmaus2::dazzler::align::Overlap > VOVL;

		while ( Plas->getNextOverlap(OVL) )
		{
			if ( OVL.aread != prevaread || OVL.bread != prevbread || OVL.isInverse() != previnverse )
				handle(VOVL,TPV,tspace,AW);

			prevaread = OVL.aread;
			prevbread = OVL.bread;
			previnverse = OVL.isInverse();
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
