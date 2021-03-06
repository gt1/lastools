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
#include <libmaus2/lcs/NP.hpp>

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
	libmaus2::dazzler::db::DatabaseFile const & DB,
	std::vector<uint64_t> const & RL,
	libmaus2::dazzler::align::AlignmentWriter & AW,
	libmaus2::dazzler::align::AlignmentWriter & AWsym,
	std::ostream & symkill
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

	// remove exact copies
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

	bool changed = true;

	std::vector < libmaus2::dazzler::align::Overlap > VKILL;
	std::vector < libmaus2::dazzler::align::Overlap > VNEW;
	std::vector < libmaus2::dazzler::align::Overlap > VSYM;

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
		{
			if ( killset.find(i) == killset.end() )
			{
				VOVL[ko++] = VOVL[i];
			}
			else
			{
				#if 0
				std::cerr << "Kill " << VOVL[i].getHeader().getInfo() << std::endl;
				#endif
				VKILL.push_back(VOVL[i]);
			}
		}
		VOVL.resize(ko);

		for ( uint64_t i = 0; i < Vadd.size(); ++i )
		{
			VOVL.push_back(Vadd[i]);
			VNEW.push_back(Vadd[i]);
		}

		std::sort(VOVL.begin(),VOVL.end(),libmaus2::dazzler::align::OverlapFullComparator());

		changed = (killset.size() != 0);
	}

	// remove exact copies
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

	// remove exact copies
	if ( VNEW.size() )
	{
		uint64_t o = 1;
		libmaus2::dazzler::align::OverlapFullComparator comp;

		for ( uint64_t i = 1; i < VNEW.size(); ++i )
			if ( comp(VNEW[i-1],VNEW[i]) )
			{
				VNEW[o++] = VNEW[i];
			}
			else
			{
				#if 0
				std::cerr << "[V] removing exact copy" << std::endl;
				std::cerr << VNEW[i-1].getHeader().getInfo() << std::endl;
				std::cerr << VNEW[i-0].getHeader().getInfo() << std::endl;
				#endif
			}

		VNEW.resize(o);
	}

	std::sort(VNEW.begin(),VNEW.end(),libmaus2::dazzler::align::OverlapFullComparator());

	{
		uint64_t l0 = 0;
		uint64_t l1 = 0;

		while ( l0 < VNEW.size() && l1 < VOVL.size() )
		{
			libmaus2::dazzler::align::OverlapInfo const OInew = VNEW[l0].getHeader().getInfo();
			libmaus2::dazzler::align::OverlapInfo const OIovl = VOVL[l1].getHeader().getInfo();

			if ( OInew < OIovl )
			{
				++l0;
			}
			else if ( OIovl < OInew )
			{
				++l1;
			}
			else
			{
				assert ( OInew == OIovl );

				uint64_t h0 = l0+1;
				uint64_t h1 = l1+1;

				while ( h0 < VNEW.size() && VNEW[h0].getHeader().getInfo() == OInew )
					++h0;
				while ( h1 < VOVL.size() && VOVL[h1].getHeader().getInfo() == OIovl )
					++h1;

				libmaus2::dazzler::align::Overlap const & OVL = VNEW[l0];

				libmaus2::lcs::AlignmentTraceContainer ATC;
				libmaus2::lcs::NP aligner;

				std::basic_string<uint8_t> const ua = DB.getu(OVL.aread,false);
				std::basic_string<uint8_t> const ub = DB.getu(OVL.bread,OVL.isInverse());

				libmaus2::dazzler::align::Overlap const ROVL = OVL.getSwapped(
					tspace,
					ua.c_str(),
					ua.size(),
					ub.c_str(),
					ub.size(),
					ATC,
					aligner
				);

				VSYM.push_back(ROVL);

				l0 = h0;
				l1 = h1;
			}
		}
	}

	std::sort(VSYM.begin(),VSYM.end());

	#if 0
	if ( VOVL.size() )
		std::cerr << "[V] left " << VOVL[0].aread << "," << VOVL[0].bread << "," << VOVL[0].isInverse() << std::endl;
	#endif

	for ( uint64_t i = 0; i < VOVL.size(); ++i )
		AW.put(VOVL[i]);
	for ( uint64_t i = 0; i < VSYM.size(); ++i )
		AWsym.put(VOVL[i]);

	for ( uint64_t i = 0; i < VKILL.size(); ++i )
	{
		libmaus2::dazzler::align::OverlapInfo const finfo = VKILL[i].getHeader().getInfo();
		libmaus2::dazzler::align::OverlapInfo const info = finfo.swappedStraight(RL.begin());

		assert ( (finfo.aread & 1) == 0 );
		assert ( (info.aread & 1) == 0 );

		finfo.serialise(symkill);
		info.serialise(symkill);
	}

	VOVL.resize(0);
}

int lasdedup(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
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
	libmaus2::dazzler::align::AlignmentWriter AWsym(libmaus2::util::OutputFileNameTools::clipOff(outfn,".las") + ".sym.las",tspace);

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
				handle(VOVL,TPV,tspace,DB,RL,AW,AWsym,*Psymkill);

			prevaread = OVL.aread;
			prevbread = OVL.bread;
			previnverse = OVL.isInverse();
			VOVL.push_back(OVL);
		}

		handle(VOVL,TPV,tspace,DB,RL,AW,AWsym,*Psymkill);
	}

	Psymkill->flush();
	Psymkill.reset();

	libmaus2::sorting::SerialisingSortingBufferedOutputFile<libmaus2::dazzler::align::OverlapInfo>::sort(symkillfn,16*1024*1024);
	libmaus2::aio::OutputStreamFactoryContainer::rename (symkillfn, arg[0] + ".symkill");
	libmaus2::dazzler::align::OverlapInfoIndexer::createInfoIndex(arg[0] + ".symkill",DB.size());

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
			return lasdedup(arg,arginfo);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
