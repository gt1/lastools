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

#include <iostream>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/util/iterator.hpp>
#include <queue>

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

struct Kmer
{
	uint64_t v;
	uint64_t p;

	Kmer() {}
	Kmer(uint64_t const rv, uint64_t const rp) : v(rv), p(rp) {}

	bool operator<(Kmer const & K) const
	{
		if ( v != K.v )
			return v < K.v;
		else
			return p < K.p;
	}
};

std::string kmerToString(uint64_t const v, uint64_t const k)
{
	std::ostringstream ostr;

	for ( uint64_t i = 0; i < k; ++i )
		ostr.put(libmaus2::fastx::remapChar((v >> (2*(k-i-1)))&3));

	return ostr.str();
}

uint64_t getKmers(std::string const & s, uint64_t const k, libmaus2::autoarray::AutoArray<Kmer> & K)
{
	uint64_t const n = s.size();
	uint64_t const numk = (n >= k) ? (n-k+1) : 0;
	uint64_t const kmask = libmaus2::math::lowbits(2*(k-1));
	assert ( k );

	uint64_t v = 0;
	for ( uint64_t i = 0; i < k-1; ++i )
	{
		v <<= 2;
		v |= libmaus2::fastx::mapChar(s[i]);
	}

	uint64_t o = 0;
	for ( uint64_t i = 0; i < numk ; ++i )
	{
		v &= kmask;
		v <<= 2;
		v |= libmaus2::fastx::mapChar(s[i+k-1]);
		K.push(o,Kmer(v,i));

		// std::cerr << kmerToString(v,k) << std::endl;
	}

	std::sort(K.begin(),K.begin()+o);

	return o;
}

struct HeapTodo
{
	Kmer const * ka;
	Kmer const * kb_c;
	Kmer const * kb_e;

	HeapTodo()
	{

	}

	HeapTodo(
		Kmer const * rka,
		Kmer const * rkb_c,
		Kmer const * rkb_e)
	: ka(rka), kb_c(rkb_c), kb_e(rkb_e)
	{

	}

	bool hasNext() const
	{
		return kb_c+1 != kb_e;
	}

	HeapTodo getNext()
	{
		HeapTodo T = *this;
		++T.kb_c;
		return T;
	}

	int64_t getDiag() const
	{
		return static_cast<int64_t>(ka->p) - static_cast<int64_t>(kb_c->p);
	}

	bool operator<(HeapTodo const & O) const
	{
		return getDiag() < O.getDiag();
	}

	bool operator>(HeapTodo const & O) const
	{
		return getDiag() > O.getDiag();
	}
};

struct Match
{
	int64_t off;
	int64_t l;

	Match()
	{

	}
	Match(int64_t const roff, int64_t const rl)
	: off(roff), l(rl) {}

	bool operator<(Match const & M) const
	{
		return off < M.off;
	}

	int64_t getFrom() const
	{
		return off;
	}

	int64_t getTo() const
	{
		return off + l;
	}

	libmaus2::math::IntegerInterval<int64_t> getInterval() const
	{
		return libmaus2::math::IntegerInterval<int64_t>(getFrom(),getTo()-1);
	}
};

void process(
	std::string const & sa, std::string const & sb,
	uint64_t const r, uint64_t const k, int64_t const tspace,
	libmaus2::dazzler::align::AlignmentWriter & AW, int64_t const aid, int64_t const bid,
	uint64_t const minlen
)
{
	std::string const ra = sa;
	std::string const rb = r ? libmaus2::fastx::reverseComplementUnmapped(sb) : sb;

	libmaus2::autoarray::AutoArray<Kmer> KA;
	libmaus2::autoarray::AutoArray<Kmer> KB;

	// compute k-mers in read a and b
	uint64_t const oa = getKmers(ra,k,KA);
	uint64_t const ob = getKmers(rb,k,KB);

	uint64_t ia = 0;
	uint64_t ib = 0;

	std::priority_queue< HeapTodo,std::vector<HeapTodo>,std::greater<HeapTodo> > Q;

	// look for matching k-mers
	while ( ia < oa && ib < ob )
	{
		if ( KA[ia].v < KB[ib].v )
			++ia;
		else if ( KB[ib].v < KA[ia].v )
			++ib;
		else
		{
			assert ( KA[ia].v == KB[ib].v );

			uint64_t ja = ia+1;
			while ( ja < oa && KA[ja].v == KA[ia].v )
				++ja;
			uint64_t jb = ib+1;
			while ( jb < ob && KB[jb].v == KB[ib].v )
				++jb;

			// reverse region in B so positions on B are in decreasing order
			std::reverse(KB.begin()+ib,KB.begin()+jb);

			for ( uint64_t z = ia; z < ja; ++z )
				Q.push(HeapTodo(KA.begin()+z,KB.begin()+ib,KB.begin()+jb));

			ia = ja;
			ib = jb;
		}
	}

	libmaus2::lcs::NNP nnp;
	libmaus2::lcs::NNPTraceContainer nnptracecontainer;
	libmaus2::lcs::AlignmentTraceContainer ATC;
	libmaus2::autoarray::AutoArray<Match const *> QR;
	libmaus2::util::SimpleQueue<libmaus2::geometry::RangeSet<Match>::search_q_element> Rtodo;

	std::map< int64_t, libmaus2::geometry::RangeSet<Match>::shared_ptr_type > D;

	#if defined(D_DEBUG)
	std::map< int64_t, std::vector<libmaus2::math::IntegerInterval<int64_t> > > TT;
	#endif

	libmaus2::lcs::NNPAlignResult bestres(0,1,0,1,1);

	int64_t prevdiag = std::numeric_limits<int64_t>::min();

	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPV;
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPVold;
	std::set< libmaus2::dazzler::align::TracePoint > TPVseen;
	std::map< uint64_t, libmaus2::dazzler::align::Overlap > OVLseen;

	// process matching k-mers in order of increasing diagonal index
	for ( uint64_t traceid = 0; ! Q.empty(); ++traceid )
	{
		HeapTodo H = Q.top();
		Q.pop();

		if ( H.hasNext() )
			Q.push(H.getNext());

		int64_t const nd = H.getDiag();

		assert ( nd >= prevdiag );
		prevdiag = nd;

		// remove diagonals which we no longer need
		while ( D.begin() != D.end() && D.begin()->first < nd )
		{
			// std::cerr << "removing " << D.begin()->first << std::endl;
			D.erase(D.begin());
		}
		// std::cerr << nd << " " << H.ka->p << " " << H.kb_c->p << std::endl;

		// check whether seed was already processed (was used on a previous alignment)
		if ( D.find(nd) != D.end() && D.find(nd)->second->search(Match(std::min(H.ka->p,H.kb_c->p),k),Rtodo) )
			continue;

		// compute alignment
		libmaus2::lcs::NNPAlignResult res;
		if ( aid == bid )
		{
			res = nnp.align(
				ra.begin(),ra.end(),H.ka->p,
				ra.begin(),ra.end(),H.kb_c->p,
				nnptracecontainer,
				aid==bid
			);
		}
		else
		{
			res = nnp.align(
				ra.begin(),ra.end(),H.ka->p,
				rb.begin(),rb.end(),H.kb_c->p,
				nnptracecontainer,
				aid==bid
			);
		}

		// compute dense trace
		nnptracecontainer.computeTrace(ATC);

		libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = ATC.ta;
		libmaus2::lcs::AlignmentTraceContainer::step_type const * const te = ATC.te;

		// register matches so we can avoid starting from a seed which is already covered by an alignment
		int64_t apos = res.abpos;
		int64_t bpos = res.bbpos;

		uint64_t matchcount = 0;
		for ( ; ta != te; ++ta )
		{
			switch ( *ta )
			{
				case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
				case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
				case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
					if ( matchcount )
					{
						int64_t const d = apos - bpos;
						int64_t const off = std::min(apos,bpos)-matchcount;

						// std::cerr << "M " << matchcount << " " << d << " " << off << std::endl;

						if ( D.find(d) == D.end() )
						{
							libmaus2::geometry::RangeSet<Match>::shared_ptr_type R(
								new libmaus2::geometry::RangeSet<Match>(
									std::min(ra.size()+k,rb.size()+k)
								)
							);
							D[d] = R;
						}

						D.find(d)->second->insert(Match(off,matchcount));

						#if defined(D_DEBUG)
						TT[d].push_back(
							libmaus2::math::IntegerInterval<int64_t>(off,off+matchcount-1)
						);
						#endif

						matchcount = 0;
					}
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
					++matchcount;
					break;
				default:
					break;
			}
			switch ( *ta )
			{
				case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
					bpos++;
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
					apos++;
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
					apos++;
					bpos++;
					break;
				case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
					apos++;
					bpos++;
					break;
				default:
					break;
			}
		}

		if ( matchcount )
		{
			int64_t const d = apos - bpos;
			int64_t const off = std::min(apos,bpos)-matchcount;

			// std::cerr << "M " << matchcount << " " << d << " " << off << std::endl;

			if ( D.find(d) == D.end() )
			{
				libmaus2::geometry::RangeSet<Match>::shared_ptr_type R(
					new libmaus2::geometry::RangeSet<Match>(
						std::min(ra.size()+k,rb.size()+k)
					)
				);
				D[d] = R;
			}

			D.find(d)->second->insert(Match(off,matchcount));

			#if defined(D_DEBUG)
			TT[d].push_back(
				libmaus2::math::IntegerInterval<int64_t>(off,off+matchcount-1)
			);
			#endif

			matchcount = 0;
		}

		assert ( apos == static_cast<int64_t>(res.aepos) );
		assert ( bpos == static_cast<int64_t>(res.bepos) );

		// if match is sufficiently long
		if ( res.aepos - res.abpos >= minlen )
		{
			//std::cerr << res << std::endl;

			if ( res.getErrorRate() < bestres.getErrorRate() )
				bestres = res;

			// compute dazzler style overlap data structure
			libmaus2::dazzler::align::Overlap const OVL = libmaus2::dazzler::align::Overlap::computeOverlap(
				r ? libmaus2::dazzler::align::Overlap::getInverseFlag() : 0,
				aid,
				bid,
				res.abpos,
				res.aepos,
				res.bbpos,
				res.bepos,
				tspace,
				ATC
			);

			// get dazzler trace points
			uint64_t const tpvo = OVL.getTracePoints(tspace,traceid,TPV,0);
			libmaus2::math::IntegerInterval<int64_t> IA(res.abpos,res.aepos-1);

			bool dup = false;

			// check whether this or a previous alignment is a duplicate
			for ( uint64_t i = 0; (!dup) && i < tpvo; ++i )
			{
				// look for trace point
				std::set<libmaus2::dazzler::align::TracePoint>::const_iterator it = TPVseen.lower_bound(
					libmaus2::dazzler::align::TracePoint(TPV[i].apos,TPV[i].bpos,0)
				);
				std::vector<uint64_t> killlist;

				for (
					;
					(!dup)
					&&
					it != TPVseen.end()
					&&
					it->apos == TPV[i].apos
					&&
					it->bpos == TPV[i].bpos
					;
					++it
				)
				{
					uint64_t const oldtraceid = it->id;

					assert ( OVLseen.find(oldtraceid) != OVLseen.end() );

					libmaus2::dazzler::align::Overlap const & OVLold = OVLseen.find(oldtraceid)->second;

					libmaus2::math::IntegerInterval<int64_t> IO(OVLold.path.abpos,OVLold.path.aepos-1);
					libmaus2::math::IntegerInterval<int64_t> IC = IA.intersection(IO);

					if ( IA.diameter() <= IO.diameter() && IC.diameter() >= 0.95 * IA.diameter() )
					{
						dup = true;

						#if 0
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						std::cerr << "dup ?\n";
						std::cerr << OVL << std::endl;
						std::cerr << OVLold << std::endl;
						#endif
					}
					else if ( IO.diameter() <= IA.diameter() && IC.diameter() >= 0.95 * IO.diameter() )
					{
						killlist.push_back(oldtraceid);

						#if 0
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						std::cerr << "dup ?\n";
						std::cerr << OVLold << std::endl;
						std::cerr << OVL << std::endl;
						#endif
					}
				}

				for ( uint64_t i = 0; i < killlist.size(); ++i )
				{
					uint64_t const oldtraceid = killlist[i];
					libmaus2::dazzler::align::Overlap const & OVLold = OVLseen.find(oldtraceid)->second;

					uint64_t const tpvoo = OVLold.getTracePoints(tspace,oldtraceid,TPVold,0);
					for ( uint64_t j = 0; j < tpvoo; ++j )
					{
						assert ( TPVseen.find(TPVold[j]) != TPVseen.end() );
						TPVseen.erase(TPVold[j]);
					}

					OVLseen.erase(oldtraceid);
				}

				#if 0
				if ( TPVseen.find(std::pair<int64_t,int64_t>(TPV[i].apos,TPV[i].bpos)) != TPVseen.end() )
				{
					dup = true;
					break;
				}
				#endif
			}

			if ( ! dup )
			{
				for ( uint64_t i = 0; i < tpvo; ++i )
					TPVseen.insert(TPV[i]);
				OVLseen[traceid] = OVL;
			}

			// std::cerr << OVL << std::endl;

			#if 0
			libmaus2::lcs::AlignmentPrint::printAlignmentLines(
				std::cerr,
				ra.begin()+res.abpos,res.aepos-res.abpos,
				rb.begin()+res.bbpos,res.bepos-res.bbpos,
				80,
				ATC.ta,ATC.te
			);
			#endif
		}
	}

	// write the alignments we produced
	for ( std::map<uint64_t,libmaus2::dazzler::align::Overlap>::const_iterator ita = OVLseen.begin();
		ita != OVLseen.end(); ++ita )
		AW.put(ita->second);

	//std::cerr << "best " << bestres << std::endl;
}

// pair file accessor class
struct PairFileReader
{
	typedef PairFileReader this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	std::string const fn;
	mutable libmaus2::aio::InputStreamInstance ISI;

	typedef libmaus2::util::ConstIterator<PairFileReader,std::pair<uint64_t,uint64_t> > const_iterator;

	PairFileReader(std::string const & rfn)
	: fn(rfn), ISI(fn) {}

	uint64_t size() const
	{
		ISI.clear();
		ISI.seekg(0,std::ios::end);
		uint64_t const p = ISI.tellg();
		assert ( p % 2*sizeof(uint64_t) == 0 );
		return p / (2*sizeof(uint64_t));
	}

	std::pair<uint64_t,uint64_t>  operator[](uint64_t const i) const
	{
		ISI.clear();
		ISI.seekg(i*2*sizeof(uint64_t),std::ios::beg);
		std::pair<uint64_t,uint64_t> M;
		M.first  = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		M.second = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		return M;
	}

	std::pair<uint64_t,uint64_t>  get(uint64_t const i) const
	{
		return operator[](i);
	}

	const_iterator begin() const
	{
		return const_iterator(this,0);
	}

	const_iterator end() const
	{
		return const_iterator(this,size());
	}
};

int tandemaligner(libmaus2::util::ArgParser const & arg)
{
	uint64_t const k = arg.uniqueArgPresent("k") ? arg.getUnsignedNumericArg<uint64_t>("k") : 14;
	uint64_t const minlen = arg.uniqueArgPresent("l") ? arg.getUnsignedNumericArg<uint64_t>("l") : 1000;
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	int64_t const tspace = 100;

	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	std::string const outfn = arg[0];
	std::string const dbname = arg[1];
	std::string const pfn = arg[2];

	std::cerr << "[V] copying " << dbname << " to memory...";
	libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type dbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(dbname,"mem:dbprefix"));
	std::cerr << "done." << std::endl;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->fn));
	PDB->computeTrimVector();

	uint64_t start = 0;
	uint64_t stop = PairFileReader(pfn).size();

	if ( arg.uniqueArgPresent("I") )
	{
		std::string const Is = arg["I"];
		std::istringstream istr(Is);
		uint64_t cnt;
		istr >> cnt;

		if ( ! istr )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		int const c = istr.get();

		if ( ! istr || c == std::istream::traits_type::eof() || c != ',' )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		uint64_t div;
		istr >> div;

		if ( ! istr || istr.peek() != std::istream::traits_type::eof() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		if ( ! (cnt < div) )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] " << cnt << " >= " << div << std::endl;
			lme.finish();
			throw lme;
		}

		uint64_t const partsize = (stop + div - 1)/div;

		start = cnt * partsize;
		stop = std::min(start+partsize,stop);
	}

	libmaus2::autoarray::AutoArray < PairFileReader::unique_ptr_type > APFR(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		PairFileReader::unique_ptr_type tptr(new PairFileReader(pfn));
		APFR[i] = UNIQUE_PTR_MOVE(tptr);
	}

	libmaus2::dazzler::align::AlignmentWriterArray AWA(tmpfilebase + "_out.las.array",numthreads,tspace);

	std::cerr << "[V] handling [" << start << "," << stop << ") out of [" << 0 << "," << PairFileReader(pfn).size() << ")" << std::endl;

	uint64_t volatile done = 0;
	libmaus2::parallel::PosixSpinLock donelock;

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
	#endif
	for ( uint64_t i = start; i < stop; ++i )
	{
		#if defined(_OPENMP)
		uint64_t const tid = omp_get_thread_num();
		#else
		uint64_t const tid = 0;
		#endif

		PairFileReader const & PFR = *(APFR[tid]);
		std::pair<uint64_t,uint64_t> const P = PFR[i];

		std::string const ra = (*PDB)[P.first];
		std::string const rb = (*PDB)[P.second];

		process(ra,rb,0,k,tspace,AWA[tid],P.first,P.second,minlen);
		process(ra,rb,1,k,tspace,AWA[tid],P.first,P.second,minlen);

		uint64_t ldone;

		donelock.lock();
		ldone = ++done;
		donelock.unlock();

		if ( ldone % 1024 == 0 )
		{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << "[V] " << ldone << "/" << stop-start << " " << static_cast<double>(ldone) / (stop-start) << std::endl;
		}
	}

	AWA.merge(outfn,tmpfilebase + "_out.las.array.mergetmp");

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 3 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " -t<threads=numcores> -k<kmersize=14> -I<from_read,to_read> -l<minlen=1000> [options] out.las in.db in.intv\n";
			return EXIT_SUCCESS;
		}
		else
		{
			return tandemaligner(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
