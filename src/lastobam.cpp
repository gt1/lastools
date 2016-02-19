/*
    lastobam
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

#include <iostream>
#include <cassert>
#include <csignal>

#include <config.h>

#include <libmaus2/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/bambam/BamFlagBase.hpp>
#include <libmaus2/bambam/BamHeader.hpp>
#include <libmaus2/bambam/BamSeqEncodeTable.hpp>
#include <libmaus2/bambam/MdStringComputationContext.hpp>
#include <libmaus2/bambam/parallel/AddWritePendingBgzfBlockInterface.hpp>
#include <libmaus2/bambam/parallel/BgzfLinearMemCompressWorkPackageDispatcher.hpp>
#include <libmaus2/bambam/parallel/BgzfLinearMemCompressWorkPackageReturnInterface.hpp>
#include <libmaus2/bambam/parallel/FragmentAlignmentBuffer.hpp>
#include <libmaus2/bambam/parallel/FragmentAlignmentBufferAllocator.hpp>
#include <libmaus2/bambam/parallel/FragmentAlignmentBufferTypeInfo.hpp>
#include <libmaus2/bambam/parallel/GetBgzfDeflateZStreamBaseInterface.hpp>
#include <libmaus2/bambam/parallel/PutBgzfDeflateZStreamBaseInterface.hpp>
#include <libmaus2/bambam/parallel/SmallLinearBlockCompressionPendingObjectFinishedInterface.hpp>
#include <libmaus2/bambam/parallel/SmallLinearBlockCompressionPendingObjectHeapComparator.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/SimpleOverlapParser.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/fastx/CharBuffer.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/lcs/AlignerFactory.hpp>
#include <libmaus2/lcs/DalignerLocalAlignment.hpp>
#include <libmaus2/lcs/EditDistanceTraceContainer.hpp>
#include <libmaus2/lz/BgzfDeflate.hpp>
#include <libmaus2/lz/BgzfDeflateOutputCallbackMD5.hpp>
#include <libmaus2/lz/BgzfDeflateZStreamBase.hpp>
#include <libmaus2/lz/BgzfDeflateZStreamBaseAllocator.hpp>
#include <libmaus2/lz/BgzfDeflateZStreamBaseTypeInfo.hpp>
#include <libmaus2/lz/BgzfDeflateOutputBufferBase.hpp>
#include <libmaus2/lz/BgzfDeflateOutputBufferBaseAllocator.hpp>
#include <libmaus2/lz/BgzfDeflateOutputBufferBaseTypeInfo.hpp>
#include <libmaus2/parallel/LockedFreeList.hpp>
#include <libmaus2/parallel/LockedGrowingFreeList.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/parallel/SimpleThreadWorkPackageDispatcher.hpp>
#include <libmaus2/parallel/SimpleThreadPoolWorkPackageFreeList.hpp>
#include <libmaus2/parallel/SimpleThreadPool.hpp>
#include <libmaus2/parallel/SynchronousQueue.hpp>
#include <libmaus2/parallel/TerminatableSynchronousQueue.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>

int getDefaultLevel()
{
	return -1;
}

int getDefaultVerbose()
{
	return 0;
}

int getDefaultMD5()
{
	return 0;
}

int getDefaultCalMdNm()
{
	return 0;
}

std::string getDefaultSupStoreStrat()
{
	return "none";
}

struct RgInfo
{
	std::string ID;
	std::string CN;
	std::string DS;
	std::string DT;
	std::string FO;
	std::string KS;
	std::string LB;
	std::string PG; // = fastqtobam
	std::string PI;
	std::string PL;
	std::string PM;
	std::string PU;
	std::string SM;

	RgInfo() {}
	RgInfo(libmaus2::util::ArgParser const & arg)
	:
		ID(arg.uniqueArgPresent("RGID") ? arg.getEqArg("RGID") : std::string()),
		CN(arg.uniqueArgPresent("RGCN") ? arg.getEqArg("RGCN") : std::string()),
		DS(arg.uniqueArgPresent("RGDS") ? arg.getEqArg("RGDS") : std::string()),
		DT(arg.uniqueArgPresent("RGDT") ? arg.getEqArg("RGDT") : std::string()),
		FO(arg.uniqueArgPresent("RGFO") ? arg.getEqArg("RGFO") : std::string()),
		KS(arg.uniqueArgPresent("RGKS") ? arg.getEqArg("RGKS") : std::string()),
		LB(arg.uniqueArgPresent("RGLB") ? arg.getEqArg("RGLB") : std::string()),
		PG(arg.uniqueArgPresent("RGPG") ? arg.getEqArg("RGPG") : std::string("lastobam")),
		PI(arg.uniqueArgPresent("RGPI") ? arg.getEqArg("RGPI") : std::string()),
		PL(arg.uniqueArgPresent("RGPL") ? arg.getEqArg("RGPL") : std::string()),
		PM(arg.uniqueArgPresent("RGPM") ? arg.getEqArg("RGPM") : std::string()),
		PU(arg.uniqueArgPresent("RGPU") ? arg.getEqArg("RGPU") : std::string()),
		SM(arg.uniqueArgPresent("RGSM") ? arg.getEqArg("RGSM") : std::string())
	{

	}

	std::string toString() const
	{
		std::ostringstream ostr;

		if ( ID.size() )
		{
			ostr << "@RG\tID:" << ID;

			if ( CN.size() ) ostr << "\tCN:" << CN;
			if ( DS.size() ) ostr << "\tDS:" << DS;
			if ( DT.size() ) ostr << "\tDT:" << DT;
			if ( FO.size() ) ostr << "\tFO:" << FO;
			if ( KS.size() ) ostr << "\tKS:" << KS;
			if ( LB.size() ) ostr << "\tLB:" << LB;
			if ( PG.size() ) ostr << "\tPG:" << PG;
			if ( PI.size() ) ostr << "\tPI:" << PI;
			if ( PL.size() ) ostr << "\tPL:" << PL;
			if ( PM.size() ) ostr << "\tPM:" << PM;
			if ( PU.size() ) ostr << "\tPU:" << PU;
			if ( SM.size() ) ostr << "\tSM:" << SM;

			ostr << "\n";
		}

		return ostr.str();
	}
};

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;
	ostr << "usage: " << arg.progname << " [<parameters>] <reference.db> <reference.fasta> <reads.db> <reads.fasta> <alignments.las>\n";
	ostr << "\n";
	ostr << "parameters:\n";
	ostr << " -t : number of threads (defaults to number of cores on machine)\n";
	ostr << " -m : maximum number of alignments stored per read (default: all)\n";
	ostr << " -s : base storage strategy for secondary reads (none, soft, hard)\n";
	ostr << " -c : calculate MD, NM and AS fields (0 or 1, defaults to 1)\n";
	ostr << " -M : maximum amount of memory used for loading reads (defaults to 2GB)\n";
	ostr << " -l : zlib compression level for output BAM file (defaults to zlib default)\n";
	return ostr.str();
}

struct RElement
{
	typedef RElement this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::dazzler::align::OverlapData::shared_ptr_type odata;
	uint64_t id;
	bool final;

	std::vector < std::pair<uint64_t,uint64_t> > split;

	RElement() {}
	RElement(
		libmaus2::dazzler::align::OverlapData::shared_ptr_type rodata,
		uint64_t rid,
		bool rfinal
	) : odata(rodata), id(rid), final(rfinal) {}

	bool operator<(RElement const & rhs) const
	{
		return id < rhs.id;
	}

	void computeSplit(uint64_t const tparts)
	{
		uint64_t const tpartsize = (odata->size() + tparts - 1) / tparts;
		split.resize(0);

		uint64_t l = 0;
		while ( l < odata->size() )
		{
			// upper limit, not included
			uint64_t h = std::min(odata->size(),l + tpartsize);
			assert ( h > l );

			while (
				h < odata->size() &&
				libmaus2::dazzler::align::OverlapData::getBRead(odata->getData(h-1).first) == libmaus2::dazzler::align::OverlapData::getBRead(odata->getData(h).first)
			)
				++h;

			split.push_back(std::pair<uint64_t,uint64_t>(l,h));

			l = h;
		}

		while ( split.size() < tparts )
			split.push_back(std::pair<uint64_t,uint64_t>(odata->size(),odata->size()));

		for ( uint64_t i = 1; i < split.size(); ++i )
			if ( split[i].second != split[i].first )
			{
				assert ( split[i-1].second != split[i-1].first );

				uint8_t const * prev = odata->getData(split[i-1].second-1).first;
				uint8_t const * cur = odata->getData(split[i].first).first;

				assert ( libmaus2::dazzler::align::OverlapData::getBRead(prev) != libmaus2::dazzler::align::OverlapData::getBRead(cur) );
			}
	}
};

struct RElementTypeInfo
{
	typedef RElementTypeInfo this_type;

	typedef RElement::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct RElementAllocator
{
	RElementAllocator() {}

	RElement::shared_ptr_type operator()() const
	{
		RElement::shared_ptr_type ptr(new RElement);
		return ptr;
	}
};


struct RecodePackage
{
	typedef RecodePackage this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	int64_t p_low;
	int64_t p_high;
	std::string lasfn;

	RecodePackage() : p_low(-1), p_high(-1), lasfn()
	{}

	RecodePackage(int64_t const rp_low, int64_t const rp_high, std::string const & rlasfn)
	: p_low(rp_low), p_high(rp_high), lasfn(rlasfn)
	{

	}
};

struct RecodePackageTypeInfo
{
	typedef RecodePackageTypeInfo this_type;

	typedef RecodePackage::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct RecodePackageAllocator
{
	RecodePackageAllocator() {}

	RecodePackage::shared_ptr_type operator()() const
	{
		RecodePackage::shared_ptr_type ptr(new RecodePackage);
		return ptr;
	}
};

struct OverlapDataTypeInfo
{
	typedef OverlapDataTypeInfo this_type;

	typedef libmaus2::dazzler::align::OverlapData::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct OverlapDataAllocator
{
	OverlapDataAllocator() {}

	libmaus2::dazzler::align::OverlapData::shared_ptr_type operator()() const
	{
		libmaus2::dazzler::align::OverlapData::shared_ptr_type ptr(new libmaus2::dazzler::align::OverlapData);
		return ptr;
	}
};


struct OverlapReadWorkPackage : public libmaus2::parallel::SimpleThreadWorkPackage
{
	typedef OverlapReadWorkPackage this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	OverlapReadWorkPackage() : libmaus2::parallel::SimpleThreadWorkPackage() {}
	OverlapReadWorkPackage(uint64_t const rpriority, uint64_t const rdispatcherid, uint64_t const rpackageid = 0)
	: libmaus2::parallel::SimpleThreadWorkPackage(rpriority,rdispatcherid,rpackageid)
	{

	}
	virtual ~OverlapReadWorkPackage() {}

	virtual char const * getPackageName() const
	{
		return "OverlapReadWorkPackage";
	}
};

struct OverlapReadWorkPackageFinishedInterface
{
	virtual ~OverlapReadWorkPackageFinishedInterface() {}
	virtual void overlapReadWorkPackageFinished(OverlapReadWorkPackage * package) = 0;
};


struct OverlapDataReadInterface
{
	virtual ~OverlapDataReadInterface() {}
	virtual void overlapDataRead(RElement::shared_ptr_type R) = 0;
};

struct OverlapReadWorkPackageDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef OverlapReadWorkPackageDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	uint64_t const dispatcherid;
	libmaus2::dazzler::align::SimpleOverlapParser * PSOP;
	libmaus2::parallel::LockedFreeList<libmaus2::dazzler::align::OverlapData,OverlapDataAllocator,OverlapDataTypeInfo> & ODFL;
	libmaus2::parallel::PosixSpinLock lock;
	OverlapReadWorkPackageFinishedInterface & finishedInterface;
	OverlapDataReadInterface & enqueInterface;
	uint64_t volatile id;

	libmaus2::parallel::PosixSpinLock finishedLock;
	int volatile finished;

	libmaus2::parallel::LockedGrowingFreeList<RElement,RElementAllocator,RElementTypeInfo> & REFL;

	bool getFinished()
	{
		int lfinished;
		finishedLock.lock();
		lfinished = finished;
		finishedLock.unlock();
		return lfinished;
	}

	void setFinished()
	{
		finishedLock.lock();
		finished = true;
		finishedLock.unlock();
	}

	OverlapReadWorkPackageDispatcher(
		uint64_t const rdispatcherid,
		libmaus2::dazzler::align::SimpleOverlapParser * rPSOP,
		libmaus2::parallel::LockedFreeList<libmaus2::dazzler::align::OverlapData,OverlapDataAllocator,OverlapDataTypeInfo> & rODFL,
		libmaus2::parallel::LockedGrowingFreeList<RElement,RElementAllocator,RElementTypeInfo> & rREFL,
		OverlapReadWorkPackageFinishedInterface & rfinishedInterface,
		OverlapDataReadInterface & renqueInterface
	) : dispatcherid(rdispatcherid), PSOP(rPSOP), ODFL(rODFL), finishedInterface(rfinishedInterface), enqueInterface(renqueInterface), id(1), finished(0),
	    REFL(rREFL)
	{

	}

	virtual ~OverlapReadWorkPackageDispatcher() {}
	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		OverlapReadWorkPackage * BP = dynamic_cast<OverlapReadWorkPackage *>(P);
		assert ( BP );

		std::vector < RElement::shared_ptr_type > R;

		if ( lock.trylock() )
		{
			libmaus2::parallel::ScopePosixSpinLock slock(lock,true /* pre locked */);

			libmaus2::dazzler::align::OverlapData::shared_ptr_type POD;
			while ( (!getFinished()) && (POD=ODFL.getIf()) )
			{
				if ( PSOP->parseNextBlock() )
				{
					PSOP->getData().swap(*POD);
					RElement::shared_ptr_type PR = REFL.get();
					*PR = RElement(POD,id++,false /* final */);
					R.push_back(PR);
				}
				else
				{
					POD->overlapsInBuffer = 0;
					RElement::shared_ptr_type PR = REFL.get();
					*PR = RElement(POD,id++,true /* final */);
					R.push_back(PR);
					setFinished();
				}
			}
		}

		for ( uint64_t i = 0; i < R.size(); ++i )
			enqueInterface.overlapDataRead(R[i]);

		finishedInterface.overlapReadWorkPackageFinished(BP);
	}
};


struct AlignerTypeInfo
{
	typedef AlignerTypeInfo this_type;

	typedef libmaus2::lcs::Aligner::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct AlignerAllocator
{
	AlignerAllocator() {}

	libmaus2::lcs::Aligner::shared_ptr_type operator()() const
	{
		std::set<libmaus2::lcs::AlignerFactory::aligner_type> const S = libmaus2::lcs::AlignerFactory::getSupportedAligners();

		if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_y256_8) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_y256_8));
			libmaus2::lcs::Aligner::shared_ptr_type S(T.release());
			return S;
		}
		else if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8));
			libmaus2::lcs::Aligner::shared_ptr_type S(T.release());
			return S;
		}
		else if ( S.size() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(*(S.begin())));
			libmaus2::lcs::Aligner::shared_ptr_type S(T.release());
			return S;
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "AlignerAllocator: no aligners found" << std::endl;
			lme.finish();
			throw lme;
		}

	}
};

struct LASToBAMConverter
{
	typedef LASToBAMConverter this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::autoarray::AutoArray < std::pair<uint16_t,uint16_t> > Apath;

	std::vector<uint64_t> const & refoff;
	libmaus2::autoarray::AutoArray<char> const & refdata;
	uint64_t const ref_low;
	uint64_t const ref_high;

	std::vector<uint64_t> const & readsoff;
	libmaus2::autoarray::AutoArray<char> const & readsdata;
	uint64_t const reads_low;
	uint64_t const reads_high;

	libmaus2::autoarray::AutoArray<char const *> const & Preadnames;

	#define USE_DALIGNER

	#if defined(LIBMAUS2_HAVE_DALIGNER) && defined(USE_DALIGNER)
	libmaus2::lcs::DalignerLocalAlignment DLA;
	#endif

	libmaus2::lcs::EditDistanceTraceContainer EATC;
	libmaus2::lcs::EditDistanceTraceContainer & ATC;
	libmaus2::lcs::Aligner::unique_ptr_type Pal;
	libmaus2::lcs::Aligner & ND;
	::libmaus2::fastx::UCharBuffer ubuffer;
	::libmaus2::fastx::UCharBuffer sbuffer;
	::libmaus2::fastx::UCharBuffer tbuffer;
	::libmaus2::bambam::MdStringComputationContext context;
	::libmaus2::bambam::BamSeqEncodeTable const seqenc;

	int64_t const tspace;
	bool const small;

	bool const calmdnm;

	libmaus2::autoarray::AutoArray<std::pair<libmaus2::lcs::AlignmentTraceContainer::step_type,uint64_t> > Aopblocks;
	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Acigop;
	libmaus2::autoarray::AutoArray<uint8_t> ASQ;

	enum supplementary_seq_strategy_t {
		supplementary_seq_strategy_soft,
		supplementary_seq_strategy_hard,
		supplementary_seq_strategy_none
	};

	supplementary_seq_strategy_t const supplementary_seq_strategy;

	std::string const rgid;


	static libmaus2::lcs::Aligner::unique_ptr_type constructAligner()
	{
		std::set<libmaus2::lcs::AlignerFactory::aligner_type> const S = libmaus2::lcs::AlignerFactory::getSupportedAligners();

		if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_Daligner_NP) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_Daligner_NP));
			return UNIQUE_PTR_MOVE(T);
		}
		else if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_y256_8) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_y256_8));
			return UNIQUE_PTR_MOVE(T);
		}
		else if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8));
			return UNIQUE_PTR_MOVE(T);
		}
		else if ( S.find(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_NP) != S.end() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_NP));
			return UNIQUE_PTR_MOVE(T);
		}
		else if ( S.size() )
		{
			libmaus2::lcs::Aligner::unique_ptr_type T(libmaus2::lcs::AlignerFactory::construct(*(S.begin())));
			return UNIQUE_PTR_MOVE(T);
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "LASToBAMConverter::constructAligner: no aligners found" << std::endl;
			lme.finish();
			throw lme;
		}
	}

	LASToBAMConverter(
		std::vector<uint64_t> const & rrefoff,
		libmaus2::autoarray::AutoArray<char> const & rrefdata,
		uint64_t const rref_low,
		uint64_t const rref_high,
		std::vector<uint64_t> const & rreadsoff,
		libmaus2::autoarray::AutoArray<char> const & rreadsdata,
		uint64_t const rreads_low,
		uint64_t const rreads_high,
		libmaus2::autoarray::AutoArray<char const *> const & rPreadnames,
		int64_t const rtspace,
		bool const rcalmdnm,
		supplementary_seq_strategy_t const rsupplementaryStrategy,
		std::string const rrgid
	)
	:
	  refoff(rrefoff),
	  refdata(rrefdata),
	  ref_low(rref_low),
	  ref_high(rref_high),
	  readsoff(rreadsoff),
	  readsdata(rreadsdata),
	  reads_low(rreads_low),
	  reads_high(rreads_high),
	  Preadnames(rPreadnames),
	  #if defined(LIBMAUS2_HAVE_DALIGNER) && defined(USE_DALIGNER)
	  DLA(),
	  EATC(),
	  ATC(DLA),
	  #else
	  EATC(),
	  ATC(EATC),
	  #endif
	  Pal(constructAligner()),
	  ND(*Pal),
	  tspace(rtspace),
	  small(libmaus2::dazzler::align::AlignmentFile::tspaceToSmall(tspace)),
	  calmdnm(rcalmdnm),
	  supplementary_seq_strategy(rsupplementaryStrategy),
	  rgid(rrgid)
	{
	}

	void operator()(
		// the overlap
		uint8_t const * OVL,
		// buffer for storing bam record,
		libmaus2::bambam::parallel::FragmentAlignmentBufferFragment & FABR,
		// is this the primary alignment for this read?
		bool const primary,
		// header
		libmaus2::bambam::BamHeader const & bamheader
	)
	{
		int64_t const aread = libmaus2::dazzler::align::OverlapData::getARead(OVL);
		int64_t const bread = libmaus2::dazzler::align::OverlapData::getBRead(OVL);

		if (
			!
			(
				aread >= static_cast<int64_t>(ref_low) &&
				aread < static_cast<int64_t>(ref_high) &&
				bread >= static_cast<int64_t>(reads_low) &&
				bread < static_cast<int64_t>(reads_high)
			)
		)
		{
			return;
		}

		bool const bIsInverse = libmaus2::dazzler::align::OverlapData::getFlags(OVL) & 1;
		int64_t const abpos = libmaus2::dazzler::align::OverlapData::getABPos(OVL);
		int64_t const aepos = libmaus2::dazzler::align::OverlapData::getAEPos(OVL);
		int64_t const bbpos = libmaus2::dazzler::align::OverlapData::getBBPos(OVL);
		int64_t const bepos = libmaus2::dazzler::align::OverlapData::getBEPos(OVL);
		uint64_t const tracelen = libmaus2::dazzler::align::OverlapData::decodeTraceVector(OVL,Apath,small);

		uint64_t const refbaseoff = refoff [ aread - ref_low ];
		// length of padded ref seq
		uint64_t const reflenp = refoff [ (aread + 1) - ref_low ] - refoff [ aread - ref_low ];
		// length of ref seq
		uint64_t const reflen = reflenp - 2;
		char const * aptr = refdata.begin() + refbaseoff + 1;

		uint64_t const readbaseoff = readsoff [ bread - reads_low ];
		// length of padded read
		uint64_t const readlenp = (readsoff [ bread - reads_low + 1 ] - readsoff [ bread - reads_low ]) >> 1;
		// length of read
		uint64_t const readlen = readlenp - 2;
		// pointer to b read
		char const * bptr = readsdata.begin() + readbaseoff + (bIsInverse ? readlenp : 0) + 1;

		char const * readname = Preadnames.at(bread - reads_low);

		int32_t const blen   = (bepos - bbpos);
		int32_t const bclipleft = bbpos;
		int32_t const bclipright = readlen - blen - bclipleft;

		#if defined(LIBMAUS2_HAVE_DALIGNER) && defined(USE_DALIGNER)
		/* LocalEditDistanceResult const LEDR = */ DLA.computeDenseTracePreMapped(
			reinterpret_cast<uint8_t const *>(aptr),reflen,
			reinterpret_cast<uint8_t const *>(bptr),readlen,
			tspace,
			Apath.begin(),tracelen,
			libmaus2::dazzler::align::OverlapData::getDiffs(OVL),
			abpos,
			bbpos,
			aepos,
			bepos);
		#else
		// compute dense alignment trace
		libmaus2::dazzler::align::Overlap::computeTracePreMapped(Apath.begin(),tracelen,abpos,aepos,bbpos,bepos,reinterpret_cast<uint8_t const *>(aptr),reinterpret_cast<uint8_t const *>(bptr),tspace,ATC,ND);
		#endif

		// compute cigar operation blocks based on alignment trace
		size_t const nopblocks = ATC.getOpBlocks(Aopblocks);
		// compute final number of cigar operations based on clipping information
		size_t const cigopblocks = nopblocks + (bclipleft?1:0) + (bclipright?1:0);

		// reallocate cigar operation vector if necessary
		if ( cigopblocks > Acigop.size() )
			Acigop.resize(cigopblocks);
		uint64_t cigp = 0;

		int64_t as = 0;
		if ( bclipleft )
		{
			if ( primary )
				Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,bclipleft);
			else
			{
				switch ( supplementary_seq_strategy )
				{
					case supplementary_seq_strategy_soft:
						Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,bclipleft);
						break;
					case supplementary_seq_strategy_hard:
					case supplementary_seq_strategy_none:
						Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP,bclipleft);
						break;
				}
			}
		}
		for ( uint64_t i = 0; i < nopblocks; ++i )
		{
			switch ( Aopblocks[i].first )
			{
				case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
					Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL,Aopblocks[i].second);
					as += Aopblocks[i].second;
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
					Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF,Aopblocks[i].second);
					as -= static_cast<int64_t>(Aopblocks[i].second);
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
					Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS,Aopblocks[i].second);
					as -= static_cast<int64_t>(Aopblocks[i].second);
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
					Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL,Aopblocks[i].second);
					as -= static_cast<int64_t>(Aopblocks[i].second);
					break;
				default:
					break;
			}
		}
		if ( bclipright )
		{
			if ( primary )
				Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,bclipright);
			else
			{
				switch ( supplementary_seq_strategy )
				{
					case supplementary_seq_strategy_soft:
						Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP,bclipright);
						break;
					case supplementary_seq_strategy_hard:
					case supplementary_seq_strategy_none:
						Acigop[cigp++] = libmaus2::bambam::cigar_operation(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP,bclipright);
						break;
				}
			}
		}

		// reallocate quality value buffer if necessary
		if ( readlen > ASQ.size() )
		{
			ASQ.resize(readlen);
			std::fill(ASQ.begin(),ASQ.end(),255);
		}

		// length of stored sequence
		uint64_t storeseqlen;

		if ( primary )
		{
			storeseqlen = readlen;
		}
		else
		{
			switch ( supplementary_seq_strategy )
			{
				// hard clip
				case supplementary_seq_strategy_hard:
					storeseqlen = readlen - (bclipleft+bclipright);
					break;
				// do not store
				case supplementary_seq_strategy_none:
					storeseqlen = 0;
					break;
				// soft clip
				case supplementary_seq_strategy_soft:
				default:
					storeseqlen = readlen;
					break;
			}
		}

		libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignmentPreMapped(
			tbuffer,
			// seqenc,
			readname,
			strlen(readname),
			aread, // ref id
			abpos, // pos
			255, // mapq (none given)
			(bIsInverse ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE : 0)
			|
			(primary ? 0 : libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSECONDARY)
			,
			Acigop.begin(),
			cigp,
			-1, // next ref id
			-1, // next pos
			0, // template length (not given)
			primary ? bptr : (bptr + ((supplementary_seq_strategy==supplementary_seq_strategy_hard) ? bclipleft : 0)),
			storeseqlen,
			ASQ.begin(), /* quality (not given) */
			0, /* quality offset */
			true /* reset buffer */
		);

		libmaus2::bambam::libmaus2_bambam_alignment_validity val = libmaus2::bambam::BamAlignmentDecoderBase::valid(tbuffer.begin(),tbuffer.end() - tbuffer.begin());

		if ( val != libmaus2::bambam::libmaus2_bambam_alignment_validity_ok )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "LASToBAMConverter: constructed alignment is invalid: " << val << std::endl;
			::libmaus2::bambam::BamFormatAuxiliary auxdata;
			libmaus2::bambam::BamAlignmentDecoderBase::formatAlignment(lme.getStream(),tbuffer.begin(),tbuffer.end() - tbuffer.begin(),bamheader,auxdata);
			lme.getStream() << std::endl;
			lme.finish();
			throw lme;
		}

		if ( calmdnm )
		{
			try
			{
				if ( primary || (supplementary_seq_strategy != supplementary_seq_strategy_none) )
				{
					libmaus2::bambam::BamAlignmentDecoderBase::calculateMdMapped(
						tbuffer.begin(),tbuffer.end()-(tbuffer.begin()),
						context,aptr + abpos,false
					);
				}
				else
				{
					sbuffer.reset();
					libmaus2::bambam::BamAlignmentEncoderBase::encodeSeqPreMapped(sbuffer,bptr+bclipleft,readlen-(bclipleft+bclipright));

					libmaus2::bambam::BamAlignmentDecoderBase::calculateMdMapped(
						tbuffer.begin(),tbuffer.end()-(tbuffer.begin()),
						context,aptr + abpos,
						sbuffer.begin(),
						readlen-(bclipleft+bclipright),
						false
					);
				}
				libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(tbuffer,"MD",context.md.get());
				libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(tbuffer,"NM",'i',context.nm);
				libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(tbuffer,"AS",'i',as);
			}
			catch(std::exception const & ex)
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "LASToBAMConverter: calmdnm exception for read " << readname << "\n" << ex.what() << std::endl;
				::libmaus2::bambam::BamFormatAuxiliary auxdata;
				libmaus2::bambam::BamAlignmentDecoderBase::formatAlignment(lme.getStream(),tbuffer.begin(),tbuffer.end() - tbuffer.begin(),bamheader,auxdata);
				lme.finish();
				throw lme;
			}
		}

		if ( rgid.size() )
			libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(tbuffer,"RG",rgid.c_str());

		libmaus2::math::IntegerInterval<int64_t> const covered = libmaus2::bambam::BamAlignmentDecoderBase::getCoveredReadInterval(tbuffer.begin());

		bool const c_ok = (covered.from == bbpos) && (covered.to+1 == bepos);

		if ( ! c_ok )
		{
			std::ostringstream ostr;
			ostr << "covered " << covered << " expected " << bbpos << "," << bepos << std::endl;
			ostr << "front clipping " << libmaus2::bambam::BamAlignmentDecoderBase::getFrontClipping(tbuffer.begin()) << std::endl;
			ostr << "back clipping " << libmaus2::bambam::BamAlignmentDecoderBase::getBackClipping(tbuffer.begin()) << std::endl;
			ostr << "lseq " << libmaus2::bambam::BamAlignmentDecoderBase::getLseq(tbuffer.begin()) << std::endl;
			ostr << "from " << covered.from << std::endl;
			ostr << "to " << covered.to << std::endl;

			::libmaus2::bambam::BamFormatAuxiliary auxdata;
			libmaus2::bambam::BamAlignmentDecoderBase::formatAlignment(ostr,tbuffer.begin(),tbuffer.end() - tbuffer.begin(),bamheader,auxdata);
			ostr.put('\n');

			std::string const s = ostr.str();
			char const * c = s.c_str();
			std::cerr.write(c,s.size());

			assert ( c_ok );
		}

		FABR.pushAlignmentBlock(tbuffer.begin(),tbuffer.end() - tbuffer.begin());

		#if defined(LASTOBAM_ALIGNMENT_PRINT_DEBUG)
		// print alignment if requested
		if ( printAlignments )
		{
			libmaus2::lcs::AlignmentStatistics const AS = ATC.getAlignmentStatistics();

			if ( AS.getErrorRate() <= eratelimit )
			{
				std::cout << OVL << std::endl;
				std::cout << AS << std::endl;
				ATC.printAlignmentLines(std::cout,aptr + abpos,aepos-abpos,bptr + bbpos,bepos-bbpos,80);
			}
		}
		#endif

		#if 0
		::libmaus2::bambam::BamFormatAuxiliary auxdata;
		libmaus2::bambam::BamAlignmentDecoderBase::formatAlignment(std::cerr,tbuffer.begin(),tbuffer.end() - tbuffer.begin(),bamheader,auxdata);
		#endif
	}
};

struct LASToBAMConverterTypeInfo
{
	typedef LASToBAMConverterTypeInfo this_type;

	typedef LASToBAMConverter::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct LASToBAMConverterAllocator
{
	std::vector<uint64_t> const & refoff;
	libmaus2::autoarray::AutoArray<char> const & refdata;
	uint64_t const ref_low;
	uint64_t const ref_high;
	std::vector<uint64_t> const & readsoff;
	libmaus2::autoarray::AutoArray<char> const & readsdata;
	uint64_t const reads_low;
	uint64_t const reads_high;
	libmaus2::autoarray::AutoArray<char const *> const & Preadnames;
	int64_t const tspace;
	bool const calmdnm;
	LASToBAMConverter::supplementary_seq_strategy_t const supplementaryStrategy;
	std::string const rgid;

	LASToBAMConverterAllocator(
		std::vector<uint64_t> const & rrefoff,
		libmaus2::autoarray::AutoArray<char> const & rrefdata,
		uint64_t const rref_low,
		uint64_t const rref_high,
		std::vector<uint64_t> const & rreadsoff,
		libmaus2::autoarray::AutoArray<char> const & rreadsdata,
		uint64_t const rreads_low,
		uint64_t const rreads_high,
		libmaus2::autoarray::AutoArray<char const *> const & rPreadnames,
		int64_t const rtspace, /* algn.tspace */
		bool const rcalmdnm,
		LASToBAMConverter::supplementary_seq_strategy_t const rsupplementaryStrategy,
		std::string const rrgid
	)
	:
	  refoff(rrefoff),
	  refdata(rrefdata),
	  ref_low(rref_low),
	  ref_high(rref_high),
	  readsoff(rreadsoff),
	  readsdata(rreadsdata),
	  reads_low(rreads_low),
	  reads_high(rreads_high),
	  Preadnames(rPreadnames),
	  tspace(rtspace),
	  calmdnm(rcalmdnm),
	  supplementaryStrategy(rsupplementaryStrategy),
	  rgid(rrgid)
	{}

	LASToBAMConverter::shared_ptr_type operator()() const
	{
		LASToBAMConverter::shared_ptr_type Tptr(
			new LASToBAMConverter(
				refoff,refdata,ref_low,ref_high,readsoff,readsdata,reads_low,reads_high,Preadnames,tspace,calmdnm,supplementaryStrategy,rgid
			)
		);
		return Tptr;
	}
};

struct LasToBamConversionRequest
{
	typedef LasToBamConversionRequest this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type FAB;
	uint64_t volatile finished;
	libmaus2::parallel::PosixSpinLock lock;

	LasToBamConversionRequest() : FAB(), finished(false), lock() {}
	LasToBamConversionRequest(libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type rFAB) : FAB(rFAB), finished(false), lock()
	{

	}

	bool incrementFinished()
	{
		bool r = false;
		lock.lock();
		finished += 1;
		if ( finished == FAB->size() )
			r = true;
		lock.unlock();
		return r;
	}
};

struct LasToBamConversionRequestAllocator
{
	LasToBamConversionRequestAllocator()
	{}

	LasToBamConversionRequest::shared_ptr_type operator()() const
	{
		LasToBamConversionRequest::shared_ptr_type Tptr(new LasToBamConversionRequest);
		return Tptr;
	}
};

struct LasToBamConversionRequestTypeInfo
{
	typedef LasToBamConversionRequestTypeInfo this_type;

	typedef LasToBamConversionRequest::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct OVLDataInvAbPosComparator
{
	libmaus2::dazzler::align::OverlapData const * data;

	OVLDataInvAbPosComparator(libmaus2::dazzler::align::OverlapData const * rdata)
	: data(rdata)
	{

	}

	bool operator()(uint64_t const i, uint64_t const j) const
	{
		uint8_t const * idata = data->getData(i).first;
		uint8_t const * jdata = data->getData(j).first;

		int const inv_i = libmaus2::dazzler::align::OverlapData::getFlags(idata) & 1;
		int const inv_j = libmaus2::dazzler::align::OverlapData::getFlags(jdata) & 1;

		if ( inv_i != inv_j )
			return inv_i < inv_j;

		return
			libmaus2::dazzler::align::OverlapData::getABPos(idata) < libmaus2::dazzler::align::OverlapData::getABPos(jdata);
	}
};

struct OVLDataInvBePosComparator
{
	libmaus2::dazzler::align::OverlapData const * data;

	OVLDataInvBePosComparator(libmaus2::dazzler::align::OverlapData const * rdata)
	: data(rdata)
	{

	}

	bool operator()(uint64_t const i, uint64_t const j) const
	{
		uint8_t const * idata = data->getData(i).first;
		uint8_t const * jdata = data->getData(j).first;

		return
			libmaus2::dazzler::align::OverlapData::getBEPos(idata) < libmaus2::dazzler::align::OverlapData::getBEPos(jdata);
	}
};

struct OVLDataScoreComparator
{
	libmaus2::dazzler::align::OverlapData const * data;

	OVLDataScoreComparator(libmaus2::dazzler::align::OverlapData const * rdata)
	: data(rdata)
	{

	}

	static int64_t getScore(uint8_t const * edata)
	{
		int64_t const qlen = libmaus2::dazzler::align::OverlapData::getBEPos(edata)-libmaus2::dazzler::align::OverlapData::getBBPos(edata);
		int64_t const diffs = libmaus2::dazzler::align::OverlapData::getDiffs(edata);
		return qlen - diffs;
	}

	bool operator()(uint64_t const i, uint64_t const j) const
	{
		uint8_t const * idata = data->getData(i).first;
		uint8_t const * jdata = data->getData(j).first;
		return getScore(idata) > getScore(jdata);
	}
};

struct LasToBamConversionRequestPart
{
	typedef LasToBamConversionRequestPart this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	LasToBamConversionRequest::shared_ptr_type request;
	libmaus2::bambam::parallel::FragmentAlignmentBufferFragment * FABF;
	std::pair<uint64_t,uint64_t> range;
	RElement::shared_ptr_type relement;
	libmaus2::bambam::BamHeader const * bamheader;
	uint64_t maxconvert;

	LasToBamConversionRequestPart() : request(), FABF(0), range(), relement(), bamheader(0), maxconvert(0) {}

	bool dispatch(LASToBAMConverter::shared_ptr_type converter)
	{
		libmaus2::dazzler::align::OverlapData const * rdata = relement->odata.get();

		uint64_t low = range.first;

		while ( low < range.second )
		{
			int64_t const ref_bread = libmaus2::dazzler::align::OverlapData::getBRead(rdata->getData(low).first);

			uint64_t high = low+1;

			while (
				high < range.second &&
				libmaus2::dazzler::align::OverlapData::getBRead(rdata->getData(high).first) == ref_bread
			)
				++high;

			std::vector<int64_t> I(high-low);
			for ( uint64_t i = low; i < high; ++i )
				I[i-low] = i;
			std::sort(I.begin(),I.end(),OVLDataScoreComparator(rdata));

			#if 0
			std::cerr << std::string(80,'-') << std::endl;
			for ( uint64_t i = 0; i < I.size(); ++i )
			{
				libmaus2::dazzler::align::OverlapData::toString(std::cerr,rdata->getData(I[i]).first);
				std::cerr << " " << OVLDataScoreComparator::getScore(rdata->getData(I[i]).first);
				std::cerr << std::endl;
			}
			#endif

			if ( maxconvert )
			{
				(*converter)(relement->odata->getData(I[0]).first,*FABF,true /* primary */,*bamheader);
				for ( uint64_t i = 1; (i < maxconvert) && (i < I.size()); ++i )
					(*converter)(relement->odata->getData(I[i]).first,*FABF,false /* primary */,*bamheader);
			}

			low = high;
		}

		return request->incrementFinished();
	}
};

struct LasToBamConversionRequestPartAllocator
{
	LasToBamConversionRequestPartAllocator()
	{}

	LasToBamConversionRequestPart::shared_ptr_type operator()() const
	{
		LasToBamConversionRequestPart::shared_ptr_type Tptr(new LasToBamConversionRequestPart);
		return Tptr;
	}
};

struct LasToBamConversionRequestPartTypeInfo
{
	typedef LasToBamConversionRequestPartTypeInfo this_type;

	typedef LasToBamConversionRequestPart::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct ReturnLasToBamRequestPartInterface
{
	virtual ~ReturnLasToBamRequestPartInterface() {}
	virtual void returnLasToBamRequestPart(LasToBamConversionRequestPart::shared_ptr_type reqpart, bool const alldone) = 0;
};


struct LasToBamConversionRequestPartWorkPackage : public libmaus2::parallel::SimpleThreadWorkPackage
{
	typedef LasToBamConversionRequestPartWorkPackage this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	LasToBamConversionRequestPart::shared_ptr_type reqpart;

	LasToBamConversionRequestPartWorkPackage() : libmaus2::parallel::SimpleThreadWorkPackage() {}
	LasToBamConversionRequestPartWorkPackage(
		uint64_t const rpriority, uint64_t const rdispatcherid, uint64_t const rsubid,
		LasToBamConversionRequestPart::shared_ptr_type rreqpart
	)
	: libmaus2::parallel::SimpleThreadWorkPackage(rpriority,rdispatcherid,0 /* package id */,rsubid), reqpart(rreqpart)
	{

	}
	virtual ~LasToBamConversionRequestPartWorkPackage() {}

	virtual char const * getPackageName() const
	{
		return "LasToBamConversionRequestPartWorkPackage";
	}

};

struct LasToBamConversionRequestPartWorkPackageFinishedInterface
{
	virtual ~LasToBamConversionRequestPartWorkPackageFinishedInterface() {}
	virtual void lasToBamConversionRequestPartWorkPackageFinished(LasToBamConversionRequestPartWorkPackage * package) = 0;
};

struct LasToBamConversionRequestPartWorkPackageDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	typedef OverlapReadWorkPackageDispatcher this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	uint64_t const dispatcherid;
	libmaus2::parallel::LockedGrowingFreeList<LASToBAMConverter,LASToBAMConverterAllocator,LASToBAMConverterTypeInfo> & lastobamfreelist;
	ReturnLasToBamRequestPartInterface & returnLasToBamRequestPartInterface;
	LasToBamConversionRequestPartWorkPackageFinishedInterface & lasToBamConversionRequestPartWorkPackageFinishedInterface;

	LasToBamConversionRequestPartWorkPackageDispatcher(
		uint64_t const rdispatcherid,
		libmaus2::parallel::LockedGrowingFreeList<LASToBAMConverter,LASToBAMConverterAllocator,LASToBAMConverterTypeInfo> & rlastobamfreelist,
		ReturnLasToBamRequestPartInterface & rreturnLasToBamRequestPartInterface,
		LasToBamConversionRequestPartWorkPackageFinishedInterface & rlasToBamConversionRequestPartWorkPackageFinishedInterface
	) : dispatcherid(rdispatcherid), lastobamfreelist(rlastobamfreelist), returnLasToBamRequestPartInterface(rreturnLasToBamRequestPartInterface), lasToBamConversionRequestPartWorkPackageFinishedInterface(rlasToBamConversionRequestPartWorkPackageFinishedInterface)
	{

	}

	virtual ~LasToBamConversionRequestPartWorkPackageDispatcher() {}
	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		LasToBamConversionRequestPartWorkPackage * BP = dynamic_cast<LasToBamConversionRequestPartWorkPackage *>(P);
		assert ( BP );

		LasToBamConversionRequestPart::shared_ptr_type reqpart = BP->reqpart;

		LASToBAMConverter::shared_ptr_type converter = lastobamfreelist.get();
		bool const alldone = reqpart->dispatch(converter);
		lastobamfreelist.put(converter);

		returnLasToBamRequestPartInterface.returnLasToBamRequestPart(reqpart,alldone);
		lasToBamConversionRequestPartWorkPackageFinishedInterface.lasToBamConversionRequestPartWorkPackageFinished(BP);
	}
};

struct FragmentAlignmentBufferIdComparator
{
	bool operator()(libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type const & A, libmaus2::bambam::parallel::FragmentAlignmentBuffer::FragmentAlignmentBuffer::shared_ptr_type const & B) const
	{
		return A->id < B->id;
	}
};

struct RElementComp
{
	bool operator()(RElement::shared_ptr_type const & A, RElement::shared_ptr_type const & B) const
	{
		return *A < *B;
	}
};

struct RecodeControl :
	public OverlapReadWorkPackageFinishedInterface,
	public OverlapDataReadInterface,
	public ReturnLasToBamRequestPartInterface,
	public LasToBamConversionRequestPartWorkPackageFinishedInterface,
	public libmaus2::bambam::parallel::GetBgzfDeflateZStreamBaseInterface,
	public libmaus2::bambam::parallel::PutBgzfDeflateZStreamBaseInterface,
	public libmaus2::bambam::parallel::BgzfLinearMemCompressWorkPackageReturnInterface,
	public libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObjectFinishedInterface,
	public libmaus2::bambam::parallel::AddWritePendingBgzfBlockInterface
{
	std::ostream & out;

	libmaus2::bambam::BamHeader const & bamheader;

	std::vector<uint64_t> const & refoff;
	libmaus2::autoarray::AutoArray<char> const & refdata;
	uint64_t const ref_low;
	uint64_t const ref_high;

	std::vector<uint64_t> const & readsoff;
	libmaus2::autoarray::AutoArray<char> const & readsdata;
	uint64_t const reads_low;
	uint64_t const reads_high;

	libmaus2::parallel::SimpleThreadPool & STP;
	libmaus2::aio::InputStreamInstance::unique_ptr_type PSOPISI;
	libmaus2::dazzler::align::SimpleOverlapParser::unique_ptr_type PSOP;
	typedef libmaus2::parallel::LockedFreeList<libmaus2::dazzler::align::OverlapData,OverlapDataAllocator,OverlapDataTypeInfo> ODFL_type;
	ODFL_type ODFL;
	libmaus2::parallel::LockedGrowingFreeList<RElement,RElementAllocator,RElementTypeInfo> REFL;

	libmaus2::parallel::SimpleThreadPoolWorkPackageFreeList<OverlapReadWorkPackage> ORWPFL;
	libmaus2::parallel::SimpleThreadPoolWorkPackageFreeList<LasToBamConversionRequestPartWorkPackage> LTBCRPWPFL;
	libmaus2::parallel::SimpleThreadPoolWorkPackageFreeList<libmaus2::bambam::parallel::BgzfLinearMemCompressWorkPackage> BLMCWPFL;

	uint64_t const ORWPDid;
	OverlapReadWorkPackageDispatcher ORWPD;
	uint64_t const LTBCRPWPDid;
	LasToBamConversionRequestPartWorkPackageDispatcher LTBCRPWPD;
	uint64_t const BLMCWPDid;
	libmaus2::bambam::parallel::BgzfLinearMemCompressWorkPackageDispatcher BLMCWPD;

	std::set<RElement::shared_ptr_type,RElementComp> Rtodo;
	uint64_t volatile rnext;
	libmaus2::parallel::PosixSpinLock rlock;
	int volatile rfinalseen;

	std::deque<RElement::shared_ptr_type> rreadyforconversion;
	libmaus2::parallel::PosixSpinLock rreadyforconversionlock;

	libmaus2::autoarray::AutoArray<char const *> const & Preadnames;

	libmaus2::parallel::LockedGrowingFreeList<LASToBAMConverter,LASToBAMConverterAllocator,LASToBAMConverterTypeInfo> lastobamfreelist;
	libmaus2::parallel::LockedFreeList<
		libmaus2::bambam::parallel::FragmentAlignmentBuffer,
		libmaus2::bambam::parallel::FragmentAlignmentBufferAllocator,
		libmaus2::bambam::parallel::FragmentAlignmentBufferTypeInfo
	> fabfreelist;

	libmaus2::parallel::LockedGrowingFreeList<LasToBamConversionRequest,LasToBamConversionRequestAllocator,LasToBamConversionRequestTypeInfo> lastobamconversionrequestfreelist;
	libmaus2::parallel::LockedGrowingFreeList<LasToBamConversionRequestPart,LasToBamConversionRequestPartAllocator,LasToBamConversionRequestPartTypeInfo> lastobamconversionrequestpartfreelist;

	std::set < libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type, FragmentAlignmentBufferIdComparator > FABpending;
	libmaus2::parallel::PosixSpinLock FABpendinglock;
	uint64_t volatile FABpendingnext;

	int volatile recodingdone;
	libmaus2::parallel::PosixSpinLock recodingdonelock;

	int const zlevel;

	libmaus2::parallel::LockedGrowingFreeList<libmaus2::lz::BgzfDeflateZStreamBase,
		libmaus2::lz::BgzfDeflateZStreamBaseAllocator,
		libmaus2::lz::BgzfDeflateZStreamBaseTypeInfo> lzfreelist;

	libmaus2::parallel::LockedFreeList<
		libmaus2::lz::BgzfDeflateOutputBufferBase,
		libmaus2::lz::BgzfDeflateOutputBufferBaseAllocator,
		libmaus2::lz::BgzfDeflateOutputBufferBaseTypeInfo
	> bgzfbufferfreelist;

	std::priority_queue<
		libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObject,
		std::vector<libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObject>,
		libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObjectHeapComparator
	> smallbgzfpending;
	libmaus2::parallel::PosixSpinLock smallbgzfpendinglock;

	std::map<int64_t,uint64_t> smallbgzfunissued;
	libmaus2::parallel::PosixSpinLock smallbgzfunissuedlock;

	int64_t volatile smallbgzfpendingnextblock;
	uint64_t volatile smallbgzfpendingnextsubblock;
	libmaus2::parallel::PosixSpinLock smallbgzfpendingnextlock;

	std::map<int64_t,uint64_t> smallbgzfunfinished;
	libmaus2::parallel::PosixSpinLock smallbgzfunfinishedlock;

	std::map<int64_t,libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type> fabunfinished;
	libmaus2::parallel::PosixSpinLock fabunfinishedlock;

	std::map<int64_t,uint64_t> smallbgzfunwritten;
	libmaus2::parallel::PosixSpinLock smallbgzfunwrittenlock;

	int64_t volatile finalblock;
	libmaus2::parallel::PosixSpinLock finalblocklock;

	int64_t volatile writenextblock;
	uint64_t volatile writenextsubblock;
	libmaus2::parallel::PosixSpinLock writenextlock;

	int volatile writingdone;
	libmaus2::parallel::PosixSpinLock writingdonelock;

	std::map <
		std::pair<int64_t,uint64_t>,
		std::pair<libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type,libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo>
	> writepending;
	libmaus2::parallel::PosixSpinLock writependinglock;

	int64_t const maxconvert;

	struct SmallObjectComp
	{
		bool operator()(libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObject const & A,
			libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObject const & B)
		const
		{
			if ( A.blockid != B.blockid )
				return A.blockid < B.blockid;
			else
				return A.subid < B.subid;
		}
	};

	void printState()
	{
		std::cerr << "rnext=" << rnext << std::endl;
		std::cerr << "rfinalseen=" << rfinalseen << std::endl;
		std::cerr << "rreadyforconversion.size()=" << rreadyforconversion.size() << std::endl;
		std::cerr << "fabfreelist.empty()=" << fabfreelist.empty() << std::endl;
		std::cerr << "FABpending.size()=" << FABpending.size() << std::endl;
		std::cerr << "FABpendingnext=" << FABpendingnext << std::endl;
		std::cerr << "recodingdone=" << recodingdone << std::endl;
		std::cerr << "bgzfbufferfreelist.empty()=" << bgzfbufferfreelist.empty() << std::endl;
		std::cerr << "smallbgzfpending.size()=" << smallbgzfpending.size() << std::endl;
		if ( smallbgzfpending.size() )
			std::cerr << "smallbgzfpending.top()=" << smallbgzfpending.top().blockid << "," << smallbgzfpending.top().subid << std::endl;
		std::cerr << "smallbgzfpendingnextblock=" << smallbgzfpendingnextblock << std::endl;
		std::cerr << "smallbgzfpendingnextsubblock=" << smallbgzfpendingnextsubblock << std::endl;

		for ( std::map<int64_t,uint64_t>::const_iterator ita = smallbgzfunissued.begin(); ita != smallbgzfunissued.end(); ++ita )
			std::cerr << "smallbgzfunissued[" << ita->first << "]=" << ita->second << std::endl;
		for ( std::map<int64_t,uint64_t>::const_iterator ita = smallbgzfunfinished.begin(); ita != smallbgzfunfinished.end(); ++ita )
			std::cerr << "smallbgzfunfinished[" << ita->first << "]=" << ita->second << std::endl;
		std::cerr << "fabunfinished.size()=" << fabunfinished.size() << std::endl;
		for ( std::map<int64_t,uint64_t>::const_iterator ita = smallbgzfunwritten.begin(); ita != smallbgzfunwritten.end(); ++ita )
			std::cerr << "smallbgzfunwritten[" << ita->first << "]=" << ita->second << std::endl;
		std::cerr << "finalblock=" << finalblock << std::endl;

		std::cerr << "writenextblock=" << writenextblock << std::endl;
		std::cerr << "writenextsubblock=" << writenextsubblock << std::endl;

		std::cerr << "writingdone=" << writingdone << std::endl;

		std::cerr << "writepending.size()=" << writepending.size();
		for (
			std::map <
				std::pair<int64_t,uint64_t>,
				std::pair<libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type,libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo>
			>::const_iterator ita = writepending.begin(); ita != writepending.end(); ++ita )
		{
			std::cerr << "writepending contains " << ita->first.first << "," << ita->first.second << std::endl;
		}
	}

	void returnBgzfLinearMemCompressWorkPackage(libmaus2::bambam::parallel::BgzfLinearMemCompressWorkPackage * package)
	{
		BLMCWPFL.returnPackage(package);
	}

	void smallLinearBlockCompressionPendingObjectFinished(libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObject const & smallobj)
	{
		bool fabdone = false;
		{
			libmaus2::parallel::ScopePosixSpinLock slock(smallbgzfunfinishedlock);
			assert ( smallbgzfunfinished.find(smallobj.blockid) != smallbgzfunfinished.end() );
			if ( ! --smallbgzfunfinished[smallobj.blockid] )
			{
				smallbgzfunfinished.erase(smallbgzfunfinished.find(smallobj.blockid));
				fabdone = true;
			}
		}

		libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type FABreturn;

		if ( fabdone )
		{
			libmaus2::parallel::ScopePosixSpinLock slock(fabunfinishedlock);
			FABreturn = fabunfinished.find(smallobj.blockid)->second;
			fabunfinished.erase(fabunfinished.find(smallobj.blockid));
		}

		if ( FABreturn )
		{
			// std::cerr << "[V] compression of block " << FABreturn->id << " done." << std::endl;
			returnFragmentBuffer(FABreturn);
		}
	}

	bool getWritePendingNext(std::pair<int64_t,uint64_t> & key, std::pair<libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type,libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo> & value)
	{
		libmaus2::parallel::ScopePosixSpinLock swlock(writependinglock);
		libmaus2::parallel::ScopePosixSpinLock snlock(writenextlock);

		key = std::pair<int64_t,uint64_t>(writenextblock,writenextsubblock);

		if ( writepending.find(key) != writepending.end() )
		{
			value = writepending.find(key)->second;
			writepending.erase(writepending.find(key));
			return true;
		}
		else
		{
			return false;
		}
	}

	void checkWritePending()
	{
		int lwritingdone = 0;
		std::pair<int64_t,uint64_t> K;
		std::pair<libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type,libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo> P;

		while ( getWritePendingNext(K,P) )
		{
			libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type obuf = P.first;
			libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo const flushinfo = P.second;

			bool blockdone = false;
			{
				libmaus2::parallel::ScopePosixSpinLock slock(smallbgzfunwrittenlock);
				assert ( smallbgzfunwritten.find(K.first) != smallbgzfunwritten.end() );
				if ( ! --smallbgzfunwritten[K.first] )
				{
					blockdone = true;
					smallbgzfunwritten.erase(smallbgzfunwritten.find(K.first));
				}
			}

			char const * outp = reinterpret_cast<char const *>(obuf->outbuf.begin());
			uint64_t n = 0;

			assert ( flushinfo.blocks == 1 || flushinfo.blocks == 2 );

			if ( flushinfo.blocks == 1 )
			{
				/* write data to stream, one block */
				n = flushinfo.block_a_c;
			}
			else
			{
				assert ( flushinfo.blocks == 2 );
				/* write data to stream, two blocks */
				n = flushinfo.block_a_c + flushinfo.block_b_c;
			}

			out.write(outp, n);

			if ( ! out )
			{
				// check for output errors
			}

			// return block
			bgzfbufferfreelist.put(obuf);
			checkSmallPendingQueue();

			if ( blockdone )
			{
				// std::cerr << "[V] writing block " << K.first << " done" << std::endl;

				int64_t lfinalblock;
				{
					finalblocklock.lock();
					lfinalblock = finalblock;
					finalblocklock.unlock();
				}

				if ( K.first == lfinalblock )
				{
					lwritingdone = 1;
				}

				libmaus2::parallel::ScopePosixSpinLock snlock(writenextlock);
				writenextblock += 1;
				writenextsubblock = 0;
			}
			else
			{
				libmaus2::parallel::ScopePosixSpinLock snlock(writenextlock);
				writenextsubblock += 1;
			}
		}

		if ( lwritingdone )
		{
			writingdonelock.lock();
			writingdone = 1;
			writingdonelock.unlock();
		}
	}

	bool getWritingDone()
	{
		int lwritingdone;
		writingdonelock.lock();
		lwritingdone = writingdone;
		writingdonelock.unlock();
		return lwritingdone;
	}

	void addWritePendingBgzfBlock(
		int64_t const blockid,
		int64_t const subid,
		libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type obuf,
		libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo const & info
	)
	{
		{
			libmaus2::parallel::ScopePosixSpinLock slock(writependinglock);
			writepending[std::pair<int64_t,uint64_t>(blockid,subid)] = std::pair<libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type,libmaus2::lz::BgzfDeflateZStreamBaseFlushInfo>(obuf,info);
		}

		checkWritePending();
	}

	bool getFinalSeen()
	{
		int lfinalseen;
		rlock.lock();
		lfinalseen = rfinalseen;
		rlock.unlock();
		return lfinalseen;
	}

	static libmaus2::aio::InputStreamInstance::unique_ptr_type openFile(std::string const & fn, uint64_t const fna)
	{
		libmaus2::aio::InputStreamInstance::unique_ptr_type Tptr(new libmaus2::aio::InputStreamInstance(fn));
		Tptr->seekg(fna);
		return UNIQUE_PTR_MOVE(Tptr);
	}


	RecodeControl(
		std::ostream & rout,
		libmaus2::bambam::BamHeader const & rbamheader,
		std::vector<uint64_t> const & rrefoff,
		libmaus2::autoarray::AutoArray<char> const & rrefdata,
		uint64_t const rref_low, uint64_t const rref_high,
		std::vector<uint64_t> const & rreadsoff,
		libmaus2::autoarray::AutoArray<char> const & rreadsdata,
		uint64_t const rreads_low, uint64_t const rreads_high,
		std::string const & fn,
		int64_t const tspace,
		uint64_t const fna,
		uint64_t const fnb,
		libmaus2::parallel::SimpleThreadPool & rSTP,
		libmaus2::autoarray::AutoArray<char const *> const & rPreadnames,
		bool const rcalmdnm,
		LASToBAMConverter::supplementary_seq_strategy_t const rsupplementaryStrategy,
		std::string const & rgid,
		int64_t const rmaxconvert,
		int const rzlevel
	)
	:
	  out(rout),
	  bamheader(rbamheader),
	  refoff(rrefoff),
	  refdata(rrefdata),
	  ref_low(rref_low),
	  ref_high(rref_high),
	  readsoff(rreadsoff),
	  readsdata(rreadsdata),
	  reads_low(rreads_low),
	  reads_high(rreads_high),
	  STP(rSTP),
	  PSOPISI(openFile(fn,fna)),
	  PSOP(new libmaus2::dazzler::align::SimpleOverlapParser(*PSOPISI,tspace,1024*1024,libmaus2::dazzler::align::OverlapParser::overlapparser_do_not_split_b /* dont split */,fnb-fna)),
	  ODFL(2*STP.getNumThreads()),
	  ORWPDid(STP.getNextDispatcherId()),
	  ORWPD(ORWPDid,PSOP.get(),ODFL,REFL,*this,*this),
	  LTBCRPWPDid(STP.getNextDispatcherId()),
	  LTBCRPWPD(LTBCRPWPDid,lastobamfreelist,*this,*this),
	  BLMCWPDid(STP.getNextDispatcherId()),
	  BLMCWPD(*this,*this,*this,*this,*this),
	  Rtodo(), rnext(1), rlock(), rfinalseen(0),
	  Preadnames(rPreadnames),
	  lastobamfreelist(
		LASToBAMConverterAllocator(refoff,refdata,ref_low,ref_high,readsoff,readsdata,reads_low,reads_high,Preadnames,tspace,rcalmdnm,rsupplementaryStrategy,rgid)
	  ),
	  fabfreelist(
	  	4,libmaus2::bambam::parallel::FragmentAlignmentBufferAllocator(STP.getNumThreads(),1 /* pointer mult */)
	  ),
	  lastobamconversionrequestfreelist(),
	  FABpending(),
	  FABpendinglock(),
	  FABpendingnext(0),
	  recodingdone(0),
	  recodingdonelock(),
	  zlevel(rzlevel),
	  lzfreelist(libmaus2::lz::BgzfDeflateZStreamBaseAllocator(zlevel)),
	  bgzfbufferfreelist(4*STP.getNumThreads(),libmaus2::lz::BgzfDeflateOutputBufferBaseAllocator(zlevel)),
	  smallbgzfpendingnextblock(0),
	  smallbgzfpendingnextsubblock(0),
	  finalblock(-1),
	  writenextblock(0),
	  writenextsubblock(0),
	  writingdone(0),
	  writingdonelock(),
	  maxconvert(rmaxconvert)
	{
		STP.registerDispatcher(ORWPDid,&ORWPD);
		STP.registerDispatcher(LTBCRPWPDid,&LTBCRPWPD);
		STP.registerDispatcher(BLMCWPDid,&BLMCWPD);
	}

	libmaus2::lz::BgzfDeflateZStreamBase::shared_ptr_type getBgzfDeflateZStreamBase()
	{
		return lzfreelist.get();
	}

	void putBgzfDeflateZStreamBase(libmaus2::lz::BgzfDeflateZStreamBase::shared_ptr_type & ptr)
	{
		lzfreelist.put(ptr);
	}

	bool getRecodingDone()
	{
		bool ldone = false;
		recodingdonelock.lock();
		ldone = (recodingdone != 0);
		recodingdonelock.unlock();
		return ldone;
	}

	void wait()
	{
		while ( ! getWritingDone() && ! STP.isInPanicMode() )
		{
			sleep(1);
		}
	}

	void start(bool const writeheader)
	{
		libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type FAB = fabfreelist.get();
		FAB->reset();
		FAB->id = 0;

		if ( writeheader )
		{
			libmaus2::bambam::parallel::FragmentAlignmentBufferFragment * fragment = (*FAB)[0];
			std::ostringstream ostr;
			bamheader.serialise(ostr);
			std::string const sheader = ostr.str();
			char const * cheader = sheader.c_str();
			uint8_t const * uheader = reinterpret_cast<uint8_t const *>(cheader);
			fragment->push(uheader,sheader.size());
		}

		FABpending.insert(FAB);

		enqueReadPackage();
	}

	void enqueReadPackage()
	{
		OverlapReadWorkPackage * package = ORWPFL.getPackage();
		*package = OverlapReadWorkPackage(1 /* prio */, ORWPD.dispatcherid);
		STP.enque(package);
	}

	void overlapReadWorkPackageFinished(OverlapReadWorkPackage * package)
	{
		ORWPFL.returnPackage(package);
	}

	void lasToBamConversionRequestPartWorkPackageFinished(LasToBamConversionRequestPartWorkPackage * package)
	{
		LTBCRPWPFL.returnPackage(package);
	}

	void returnOverlapDataObject(libmaus2::dazzler::align::OverlapData::shared_ptr_type odata)
	{
		ODFL.put(odata);
		if ( ! getFinalSeen() )
			enqueReadPackage();
	}

	void returnFragmentBuffer(libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type FAB)
	{
		// return buffer
		fabfreelist.put(FAB);
		checkRReadyForConversion();
	}

	bool smallGetBuffer(
		libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObject & smallobj,
		libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type & obuffer
	)
	{
		libmaus2::parallel::ScopePosixSpinLock slock(smallbgzfpendinglock);
		libmaus2::parallel::ScopePosixSpinLock snlock(smallbgzfpendingnextlock);

		if (
			smallbgzfpending.size() &&
			smallbgzfpending.top().blockid == smallbgzfpendingnextblock &&
			smallbgzfpending.top().subid == smallbgzfpendingnextsubblock &&
			(obuffer = bgzfbufferfreelist.getIf())
		)
		{
			smallobj = smallbgzfpending.top();
			smallbgzfpending.pop();

			libmaus2::parallel::ScopePosixSpinLock sulock(smallbgzfunissuedlock);
			if ( ! --smallbgzfunissued[smallobj.blockid] )
			{
				smallbgzfunissued.erase(smallbgzfunissued.find(smallobj.blockid));
				smallbgzfpendingnextblock += 1;
				// std::cerr << "[V] switching to issue block " << smallbgzfpendingnextblock << " last previous " << smallbgzfpendingnextsubblock << std::endl;
				smallbgzfpendingnextsubblock = 0;
			}
			else
			{
				smallbgzfpendingnextsubblock += 1;
			}

			return true;
		}
		else
			return false;
	}

	libmaus2::parallel::PosixSpinLock cerrlock;

	void checkSmallPendingQueue()
	{
		libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObject smallobj;
		libmaus2::lz::BgzfDeflateOutputBufferBase::shared_ptr_type obuffer;

		while ( smallGetBuffer(smallobj,obuffer) )
		{
			libmaus2::bambam::parallel::BgzfLinearMemCompressWorkPackage * package = BLMCWPFL.getPackage();

			*package = libmaus2::bambam::parallel::BgzfLinearMemCompressWorkPackage(
				0 /* prio */,
				smallobj,
				obuffer,
				BLMCWPDid
			);

			package->subid = (smallobj.blockid << 24) | smallobj.subid;

			STP.enque(package);
		}
	}

	void checkFabPendingList()
	{
		libmaus2::parallel::ScopePosixSpinLock slock(FABpendinglock);

		bool done = false;

		while (
			FABpending.size() &&
			(*(FABpending.begin()))->id == FABpendingnext
		)
		{
			// std::cerr << "[V] recode buffer " << FABpendingnext << " finished" << std::endl;

			libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type FAB = *(FABpending.begin());
			FABpending.erase(FABpending.begin());

			if ( FAB->final )
			{
				done = true;

				finalblocklock.lock();
				finalblock = FAB->id;
				finalblocklock.unlock();
			}

			std::vector<std::pair<uint8_t *,uint8_t *> > V;
			FAB->getLinearOutputFragments(libmaus2::lz::BgzfConstants::getBgzfMaxBlockSize(), V);

			// std::cerr << "[V] using " << V.size() << " blocks for block " << FAB->id << std::endl;

			{
				libmaus2::parallel::ScopePosixSpinLock slock(fabunfinishedlock);
				fabunfinished[FAB->id] = FAB;
			}

			{
				libmaus2::parallel::ScopePosixSpinLock slock(smallbgzfunissuedlock);
				smallbgzfunissued[FAB->id] = V.size();
			}

			{
				libmaus2::parallel::ScopePosixSpinLock slock(smallbgzfunfinishedlock);
				smallbgzfunfinished[FAB->id] = V.size();
			}

			{
				libmaus2::parallel::ScopePosixSpinLock slock(smallbgzfunwrittenlock);
				smallbgzfunwritten[FAB->id] = V.size();
			}

			for ( uint64_t i = 0; i < V.size(); ++i )
			{
				libmaus2::bambam::parallel::SmallLinearBlockCompressionPendingObject const smallobj(FAB->id,i,V[i].first,V[i].second);
				libmaus2::parallel::ScopePosixSpinLock slock(smallbgzfpendinglock);
				smallbgzfpending.push(smallobj);
			}

			FABpendingnext += 1;

			// returnFragmentBuffer(FAB);
		}

		if ( done )
		{
			recodingdonelock.lock();
			recodingdone = 1;
			recodingdonelock.unlock();
		}

		checkSmallPendingQueue();
	}

	void returnLasToBamRequestPart(LasToBamConversionRequestPart::shared_ptr_type reqpart, bool const alldone)
	{
		if ( alldone )
		{
			{
			libmaus2::parallel::ScopePosixSpinLock slock(FABpendinglock);
			FABpending.insert(reqpart->request->FAB);
			}

			checkFabPendingList();

			// return data object
			returnOverlapDataObject(reqpart->relement->odata);
			// erase pointer
			reqpart->relement->odata = libmaus2::dazzler::align::OverlapData::shared_ptr_type();
			// return RElement
			REFL.put(reqpart->relement);
			// return request
			lastobamconversionrequestfreelist.put(reqpart->request);
		}

		lastobamconversionrequestpartfreelist.put(reqpart);

	}

	libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type FAB = fabfreelist.get();

	bool getRReady(RElement::shared_ptr_type & R, libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type & FAB)
	{
		libmaus2::parallel::ScopePosixSpinLock srreadyforconversionlock(rreadyforconversionlock);

		if ( rreadyforconversion.size() && (FAB=fabfreelist.getIf()) )
		{
			R = rreadyforconversion.front();
			rreadyforconversion.pop_front();
			return true;
		}
		else
		{
			return false;
		}
	}

	void checkRReadyForConversion()
	{
		libmaus2::bambam::parallel::FragmentAlignmentBuffer::shared_ptr_type FAB;
		RElement::shared_ptr_type R;

		while ( getRReady(R,FAB) )
		{
			R->computeSplit(FAB->size());
			assert ( R->split.size() == FAB->size() );

			// std::cerr << "[V] queuing " << R->id << " split " << R->split.size() << std::endl;

			FAB->reset();
			LasToBamConversionRequest::shared_ptr_type req = lastobamconversionrequestfreelist.get();
			req->FAB = FAB;
			req->FAB->id = R->id;
			req->FAB->final = R->final;
			req->finished = 0;

			for ( uint64_t z = 0; z < R->split.size(); ++z )
			{
				LasToBamConversionRequestPart::shared_ptr_type reqpart = lastobamconversionrequestpartfreelist.get();
				reqpart->request = req;
				reqpart->FABF = (*(req->FAB))[z];
				reqpart->range = R->split[z];
				reqpart->relement = R;
				reqpart->bamheader = &bamheader;
				reqpart->maxconvert = (maxconvert >= 0) ? static_cast<uint64_t>(maxconvert) : std::numeric_limits<uint64_t>::max();

				LasToBamConversionRequestPartWorkPackage * package = LTBCRPWPFL.getPackage();
				*package = LasToBamConversionRequestPartWorkPackage(
					0 /* prio */, LTBCRPWPDid, R->id /* sub id */, reqpart
				);

				STP.enque(package);
			}
		}
	}

	void overlapDataRead(RElement::shared_ptr_type R)
	{
		{
			libmaus2::parallel::ScopePosixSpinLock slock(rlock);
			Rtodo.insert(R);

			#if 0
			libmaus2::dazzler::align::OverlapData::shared_ptr_type odata = R.odata;
			assert ( odata );
			if ( odata->size() )
			{
				std::cerr
					<< libmaus2::dazzler::align::OverlapData::getARead(odata->getData(0).first)
					<< "\t"
					<< libmaus2::dazzler::align::OverlapData::getARead(odata->getData(odata->size()-1).first)
					<< std::endl;
			}
			#endif
		}

		{
			libmaus2::parallel::ScopePosixSpinLock slock(rlock);

			while ( Rtodo.size() && (*(Rtodo.begin()))->id == rnext )
			{
				RElement::shared_ptr_type R = *(Rtodo.begin());
				Rtodo.erase(Rtodo.begin());
				rnext += 1;

				{
				libmaus2::parallel::ScopePosixSpinLock srreadyforconversionlock(rreadyforconversionlock);
				rreadyforconversion.push_back(R);
				}

				if ( R->final )
					rfinalseen = 1;
			}
		}

		checkRReadyForConversion();
	}
};

#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>

int lastobam(libmaus2::util::ArgParser const & arg)
{

	uint64_t const threads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	int64_t const maxconvert = arg.uniqueArgPresent("m") ? static_cast<int64_t>(arg.getUnsignedNumericArg<uint64_t>("m")) : -1;

	try
	{
		std::ostream & out = std::cout;

		std::string const db_ref_fn = arg[0];
		std::string const reference = arg[1];
		std::string const db_reads_fn = arg[2];
		std::string const reads_fa_fn = arg[3];

		libmaus2::dazzler::db::DatabaseFile DB_ref(db_ref_fn);
		DB_ref.computeTrimVector();
		libmaus2::dazzler::db::DatabaseFile DB_reads(db_reads_fn);
		DB_reads.computeTrimVector();

		/* get the original read names */
		std::ostringstream readnameostr;
		std::vector<uint64_t> readnameoff;
		libmaus2::autoarray::AutoArray<char> Creadnames;
		libmaus2::autoarray::AutoArray<char const *> CPreadnames;

		{
			libmaus2::aio::InputStreamInstance ISI(reads_fa_fn);
			libmaus2::fastx::StreamFastAReaderWrapper SFAR(ISI);
			libmaus2::fastx::FastAReader::pattern_type pattern;
			uint64_t off = 0;
			uint64_t inid = 0;
			while ( SFAR.getNextPatternUnlocked(pattern) )
			{
				if ( DB_reads.isInTrimmed(inid++) )
				{
					std::string const id = pattern.getShortStringId();
					uint64_t const l = id.size();
					readnameostr.write(id.c_str(),l+1);
					readnameoff.push_back(off);
					off += (l+1);
				}
			}

			Creadnames.resize(off);
			std::string const readnames = readnameostr.str();
			char const * in = readnames.c_str();
			std::copy(in,in+off,Creadnames.begin());

			CPreadnames.resize(readnameoff.size());
			for ( uint64_t i = 0; i < readnameoff.size(); ++i )
				CPreadnames[i] = Creadnames.begin() + readnameoff[i];
		}

		std::string const supstorestrat_s = arg.uniqueArgPresent("s") ? arg["s"] : getDefaultSupStoreStrat();
		bool const calmdnm = arg.uniqueArgPresent("c") ? arg.getParsedArg<int>("c") : 1;
		int const zlevel = arg.uniqueArgPresent("l") ? arg.getParsedArg<int>("l") : Z_DEFAULT_COMPRESSION;

		LASToBAMConverter::supplementary_seq_strategy_t supstorestrat;

		if ( supstorestrat_s == "soft" )
			supstorestrat = LASToBAMConverter::supplementary_seq_strategy_soft;
		else if ( supstorestrat_s == "hard" )
			supstorestrat = LASToBAMConverter::supplementary_seq_strategy_hard;
		else if ( supstorestrat_s == "none" )
			supstorestrat = LASToBAMConverter::supplementary_seq_strategy_none;
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "unknown strategy for storing supplementary alignments: " << supstorestrat_s << " (available options are soft, hard and none)" << std::endl;
			lme.finish();
			throw lme;
		}

		libmaus2::fastx::FastAIndex::unique_ptr_type Prefindex(libmaus2::fastx::FastAIndex::load(reference+".fai"));
		libmaus2::fastx::FastAIndex const & refindex = *Prefindex;

		std::string const curdir = libmaus2::util::ArgInfo::getCurDir();
		std::string const absreference =  (reference.size() && reference[0] == '/') ? reference : (curdir+'/'+reference);
		std::ostringstream sqstream;

		for ( uint64_t i = 0; i < refindex.size(); ++i )
			if ( DB_ref.isInTrimmed(i) )
			{
				libmaus2::fastx::FastAIndexEntry const & entry = refindex[i];
				sqstream << "@SQ\t" << "SN:" << entry.name << "\tLN:" << entry.length << "\tUR:file://" << absreference << std::endl;

				if ( static_cast<int64_t>(entry.length) != static_cast<int64_t>(DB_ref.getRead(i).rlen) )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "reference length " << entry.length << " for id " << i << " is not in sync with database entry length " << DB_ref.getRead(i).rlen << std::endl;
					lme.finish();
					throw lme;
				}
			}

		// construct new header
		RgInfo const rginfo(arg);
		std::string const rgid = rginfo.ID;

		std::ostringstream headerostr;
		headerostr << "@HD\tVN:1.5\tSO:unknown\n";
		headerostr
			<< "@PG"<< "\t"
			<< "ID:" << "lastobam" << "\t"
			<< "PN:" << "lastobam" << "\t"
			<< "CL:" << arg.commandline << "\t"
			<< "VN:" << std::string(PACKAGE_VERSION)
			<< std::endl;
		headerostr << rginfo.toString();
		headerostr << sqstream.str();
		::libmaus2::bambam::BamHeader bamheader(headerostr.str());

		// std::cerr << headerostr.str();

		uint64_t maxreadmem = arg.uniqueArgPresent("M") ? static_cast<int64_t>(arg.getUnsignedNumericArg<uint64_t>("M")) : (2*1024ull*1024ull*1024ull);

		libmaus2::dazzler::db::DatabaseFile::SplitResult ref_split = DB_ref.splitDb(std::numeric_limits<uint64_t>::max());
		libmaus2::dazzler::db::DatabaseFile::SplitResult reads_split = DB_reads.splitDb(maxreadmem);

		bool writeheader = true;

		for ( uint64_t z_ref_split = 0; z_ref_split < ref_split.size(); ++z_ref_split )
		{
			libmaus2::dazzler::db::DatabaseFile::SplitResultElement const ref_interval = ref_split[z_ref_split];

			std::vector<uint64_t> refoff;
			libmaus2::autoarray::AutoArray<char> refdata;
			DB_ref.decodeReadsMappedTerm(ref_interval.low,ref_interval.high,refdata,refoff);

			for ( uint64_t z_reads_split = 0; z_reads_split < reads_split.size(); ++z_reads_split )
			{
				libmaus2::dazzler::db::DatabaseFile::SplitResultElement const reads_interval = reads_split[z_reads_split];
				std::vector<uint64_t> readsoff;
				libmaus2::autoarray::AutoArray<char> readsdata;
				DB_reads.decodeReadsAndReverseComplementMappedTerm(reads_interval.low,reads_interval.high,readsdata,readsoff);

				#if 0
				libmaus2::autoarray::AutoArray<char> Areadnames;
				libmaus2::autoarray::AutoArray<char const *> Preadnames;
                                DB_reads.getReadNameInterval(reads_interval.low,reads_interval.high,Areadnames,Preadnames);

                                for ( uint64_t i = 0; i < Preadnames.size(); ++i )
                                {
                                	assert ( Preadnames[i] );
                                	assert ( strlen(Preadnames[i]) );
				}
				#else
				libmaus2::autoarray::AutoArray<char const *> Preadnames(reads_interval.high-reads_interval.low);
				std::copy(
					CPreadnames.begin() + reads_interval.low,
					CPreadnames.begin() + reads_interval.high,
					Preadnames.begin()
				);
				#endif

				for ( uint64_t f = 4; f < arg.size(); ++f )
				{
					std::string const lasfn = arg[f];
					int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(lasfn);
					uint64_t const startpos = 12ull;
					uint64_t const endpos = libmaus2::util::GetFileSize::getFileSize(lasfn);

					libmaus2::parallel::SimpleThreadPool STP(threads);
					try
					{
						RecodeControl RC(
							out,
							bamheader,refoff,refdata,ref_interval.low,ref_interval.high,
							readsoff,readsdata,reads_interval.low,reads_interval.high,
							lasfn,tspace,startpos,endpos,
							STP,Preadnames,calmdnm,supstorestrat,rgid,maxconvert,zlevel);
						try
						{
							RC.start(writeheader);
							writeheader = false;
							RC.wait();
						}
						catch(std::exception const & ex)
						{
							std::cerr << ex.what() << std::endl;
							RC.wait();
							STP.terminate();
							throw;
						}
						STP.terminate();
					}
					catch(std::exception const & ex)
					{
						std::cerr << ex.what() << std::endl;
						STP.terminate();
						throw;
					}
					catch(...)
					{
						STP.terminate();
						throw;
					}
				}
			}
		}

		out << libmaus2::lz::BgzfDeflate<std::ostream>::getEOFBlock();
		out.flush();

		return EXIT_SUCCESS;
	}
	catch(...)
	{
		throw;
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
		else if ( arg.size() < 4 )
		{
			std::cerr << getUsage(arg);
			return EXIT_FAILURE;
		}

		return lastobam(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
