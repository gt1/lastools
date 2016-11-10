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

#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " [parameters] <out.las> <in.las> ..." << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -T : prefix for temporary files (default: create files in current working directory)\n";
	ostr << " -f : merge fan in (default: 64)\n";
	ostr << " -s : sort order (canonical or ba, default: canonical)\n";

	return ostr.str();
}

template<typename _comparator_type = libmaus2::dazzler::align::OverlapComparator>
struct OverlapBuffer
{
	typedef _comparator_type comparator_type;
	typedef OverlapBuffer<comparator_type> this_type;

	struct BlockInfo
	{
		uint64_t p;
		uint64_t numc;
		uint64_t numa;

		BlockInfo()
		{
		}

		BlockInfo(
			uint64_t const rp,
			uint64_t const rnumc,
			uint64_t const rnuma
		) : p(rp), numc(rnumc), numa(rnuma) {}
	};

	std::string tmpname;
	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI;
	uint64_t const tspace;
	bool const small;
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::Overlap> AOVL;
	uint64_t AOVLo;
	libmaus2::autoarray::AutoArray<uint64_t> I;

	comparator_type comparator;

	std::vector < BlockInfo > B;

	OverlapBuffer(std::string const & rtmpname, int64_t const rtspace, bool const rsmall, uint64_t const s)
	: tmpname(rtmpname), OSI(new libmaus2::aio::OutputStreamInstance(tmpname)), tspace(rtspace), small(rsmall), AOVL(s), AOVLo(0), comparator()
	{
		libmaus2::dazzler::align::AlignmentFile::serialiseHeader(*OSI, 0 /* novl */, tspace);
	}

	struct HeapComparator
	{
		comparator_type const & comparator;

		HeapComparator(comparator_type const & rcomparator) : comparator(rcomparator) {}

		bool operator()(
			std::pair<uint64_t,libmaus2::dazzler::align::Overlap> const & A,
			std::pair<uint64_t,libmaus2::dazzler::align::Overlap> const & B
		) const
		{
			if ( comparator(A.second,B.second) ) // A < B
				return false;
			else if ( comparator(B.second,A.second) ) // B < A
				return true;
			else
				return A.first > B.first;
		}
	};

	struct SubComparator
	{
		libmaus2::dazzler::align::Overlap const * V;
		comparator_type const & comparator;

		SubComparator(
			libmaus2::dazzler::align::Overlap const * rV,
			comparator_type const & rcomparator
		) : V(rV), comparator(rcomparator)
		{

		}

		bool operator()(uint64_t const i, uint64_t const j)
		{
			return comparator(V[i],V[j]);
		}
	};

	void implicitFlush()
	{
		uint64_t Io = 0;

		for ( uint64_t i = 0; i < AOVLo; ++i )
		{
			libmaus2::dazzler::align::Overlap const & OVL = AOVL[i];
			bool const isstart = OVL.isStart();
			bool const isnext = OVL.isNext();
			bool const isnonchain = (!isstart) && (!isnext);

			if ( isstart || isnonchain )
				I.push(Io,i);
		}

		SubComparator SC(AOVL.begin(),comparator);
		std::sort(I.begin(),I.begin()+Io,SC);

		for ( uint64_t i = 1 ; i < Io; ++i )
			assert ( !comparator(AOVL[I[i]],AOVL[I[i-1]]) );

		uint64_t const p_bef = OSI->tellp();
		uint64_t w = 0;
		for ( uint64_t i = 0; i < Io; ++i )
		{
			uint64_t j = I[i];
			assert ( j < AOVLo );
			assert ( !AOVL[j].isNext() );
			bool const isstart = AOVL[j].isStart();

			AOVL[j++].serialiseWithPath(*OSI,small);
			w += 1;

			if ( isstart )
			{
				while ( j < AOVLo && AOVL[j].isNext() )
				{
					AOVL[j++].serialiseWithPath(*OSI,small);
					w += 1;
				}
			}
		}

		assert ( w == AOVLo );

		if ( Io )
			B.push_back(BlockInfo(p_bef,Io,AOVLo));

		AOVLo = 0;
	}

	void flush()
	{
		if ( AOVLo )
		{
			implicitFlush();
		}
	}

	void putVector(std::vector < libmaus2::dazzler::align::Overlap > & TV)
	{
		if ( TV.size() + AOVLo > AOVL.size() )
		{
			implicitFlush();
			assert ( AOVLo == 0 );
		}
		if ( TV.size() + AOVLo > AOVL.size() )
		{
			assert ( AOVLo == 0 );
			AOVL = libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::Overlap>(0);
			AOVL = libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::Overlap>(TV.size());
		}
		assert ( TV.size() + AOVLo <= AOVL.size() );

		for ( uint64_t i = 0; i < TV.size(); ++i )
			AOVL.push(AOVLo,TV[i]);

		TV.resize(0);
	}

	std::string merge(std::string const & tmpfilebase, uint64_t tmpid, uint64_t const fanin)
	{
		flush();
		OSI->flush();

		if ( B.size() == 1 )
		{
			OSI->seekp(0);
			uint64_t offset = 0;
			libmaus2::dazzler::align::AlignmentFile::putLittleEndianInteger8(*OSI,B[0].numa,offset);
		}

		OSI.reset();

		HeapComparator HC(comparator);

		while ( B.size() > 1 )
		{
			std::string nexttmpfn = tmpfilebase + libmaus2::util::NumberSerialisation::formatNumber(tmpid++,6);
			libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(nexttmpfn));
			libmaus2::dazzler::align::AlignmentFile::serialiseHeader(*OSI, 0 /* novl */, tspace);

			uint64_t const numblocks = B.size();
			uint64_t const numblockruns = (numblocks + fanin - 1)/fanin;

			std::vector<BlockInfo> NB;

			for ( uint64_t r = 0; r < numblockruns; ++r )
			{
				uint64_t const blow = r * fanin;
				uint64_t const bhigh = std::min(static_cast<uint64_t>(blow+fanin),static_cast<uint64_t>(B.size()));
				uint64_t const brange = bhigh - blow;

				std::vector<BlockInfo> BB(B.begin()+blow,B.begin()+bhigh);

				std::priority_queue <
					std::pair<uint64_t,libmaus2::dazzler::align::Overlap>,
					std::vector< std::pair<uint64_t,libmaus2::dazzler::align::Overlap> >,
					HeapComparator > Q(HC);
				libmaus2::dazzler::align::Overlap OVL;

				libmaus2::autoarray::AutoArray<libmaus2::aio::InputStreamInstance::unique_ptr_type> AIN(brange);
				uint64_t s = 0;
				for ( uint64_t i = 0; i < BB.size(); ++i )
				{
					BlockInfo const & BI = BB[i];

					libmaus2::aio::InputStreamInstance::unique_ptr_type tptr(
						new libmaus2::aio::InputStreamInstance(tmpname)
					);

					tptr->seekg(BI.p);

					AIN[i] = UNIQUE_PTR_MOVE(tptr);

					if ( BI.numc )
					{
						libmaus2::dazzler::align::AlignmentFile::readOverlap(*AIN[i],OVL,s,small);
						Q.push(std::pair<uint64_t,libmaus2::dazzler::align::Overlap>(i,OVL));
					}
				}

				uint64_t p_beg = OSI->tellp();
				uint64_t numc = 0;
				uint64_t numa = 0;

				while ( !Q.empty() )
				{
					std::pair<uint64_t,libmaus2::dazzler::align::Overlap> const P = Q.top();
					Q.pop();

					uint64_t const id = P.first;

					P.second.serialiseWithPath(*OSI,small);
					numc += 1;
					numa += 1;
					BB[id].numa -= 1;
					BB[id].numc -= 1;

					while ( BB[id].numa )
					{
						libmaus2::dazzler::align::AlignmentFile::readOverlap(*AIN[id],OVL,s,small);

						if ( OVL.isNext() )
						{
							assert ( P.second.isStart() );
							OVL.serialiseWithPath(*OSI,small);
							BB[id].numa -= 1;
							numa += 1;
						}
						else
						{
							Q.push(std::pair<uint64_t,libmaus2::dazzler::align::Overlap>(id,OVL));
							break;
						}
					}
				}

				NB.push_back(BlockInfo(p_beg,numc,numa));
			}

			OSI->flush();

			if ( NB.size() == 1 )
			{
				OSI->seekp(0);
				uint64_t offset = 0;
				libmaus2::dazzler::align::AlignmentFile::putLittleEndianInteger8(*OSI,NB[0].numa,offset);
			}

			OSI.reset();
			B = NB;
			libmaus2::aio::FileRemoval::removeFile(tmpname);
			tmpname = nexttmpfn;
		}

		return tmpname;
	}
};

template<typename comparator_type>
int laschainsortTemplate(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	std::vector<std::string> Vin;
	for ( uint64_t i = 1; i < arg.size(); ++i )
		Vin.push_back(arg[i]);

	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vin);
	bool const small = libmaus2::dazzler::align::AlignmentFile::tspaceToSmall(tspace);
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : arginfo.getDefaultTmpFileName();
	uint64_t const mergefanin = arg.uniqueArgPresent("f") ? arg.getUnsignedNumericArg<uint64_t>("f") : 64;

	uint64_t tmpid = 0;
	std::string tmpfn = tmpfilebase + libmaus2::util::NumberSerialisation::formatNumber(tmpid++,6);
	uint64_t const ovlbuffersize = 1024*1024;
	// uint64_t const ovlbuffersize = 32;

	OverlapBuffer<comparator_type> OB(tmpfn,tspace,small,ovlbuffersize);

	for ( uint64_t i = 0; i < Vin.size(); ++i )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pin(
			libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(Vin[i]));

		std::vector < libmaus2::dazzler::align::Overlap > TV;
		libmaus2::dazzler::align::Overlap OVL;

		while ( Pin->getNextOverlap(OVL) )
		{
			bool const isstart = OVL.isStart();
			bool const isnext = OVL.isNext();
			bool const isnonchain = (!isstart) && (!isnext);

			if ( isstart || isnonchain )
				OB.putVector(TV);

			TV.push_back(OVL);
		}

		if ( TV.size() )
			OB.putVector(TV);
	}

	std::string const finfn = OB.merge(tmpfilebase, tmpid, mergefanin);

	libmaus2::aio::OutputStreamFactoryContainer::rename(finfn,arg[0]);

	{
		bool prevvalid = false;
		libmaus2::dazzler::align::Overlap OVL, OVLprev;
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pin(
			libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(arg[0]));
		comparator_type comparator;

		while ( Pin->getNextOverlap(OVL) )
		{
			bool const isstart = OVL.isStart();
			bool const isnext = OVL.isNext();
			bool const isnonchain = (!isstart) && (!isnext);

			if ( isstart || isnonchain )
			{
				if ( prevvalid )
				{
					bool const ok = ! comparator(OVL,OVLprev);
					assert ( ok );
				}

				prevvalid = true;
				OVLprev = OVL;
			}
			else
			{
				assert ( prevvalid );
				assert ( OVLprev.aread == OVL.aread );
				assert ( OVLprev.bread == OVL.bread );
			}
		}
	}

	return EXIT_SUCCESS;

}

int laschainsort(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	std::string const sortorder = arg.uniqueArgPresent("s") ? arg["s"] : std::string("canonical");

	if ( sortorder == "ba" )
		return laschainsortTemplate<libmaus2::dazzler::align::OverlapComparatorBReadARead>(arg,arginfo);
	else if ( sortorder == "aidbid" )
		return laschainsortTemplate<libmaus2::dazzler::align::OverlapComparatorAIdBId>(arg,arginfo);
	else if ( sortorder == "bidaid" )
		return laschainsortTemplate<libmaus2::dazzler::align::OverlapComparatorBIdAId>(arg,arginfo);
	else
		return laschainsortTemplate<libmaus2::dazzler::align::OverlapComparator>(arg,arginfo);
}

int main(int argc, char * argv[])
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

		return laschainsort(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

}
