/*
    lastools
    Copyright (C) 2017 German Tischler

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
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/LineBuffer.hpp>
#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/aio/ArrayFile.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSort.hpp>
#include <libmaus2/lcs/NNPCor.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/dazzler/align/AlignmentFile.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <libmaus2/bambam/CigarStringParser.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " out.las query.fasta [target.fasta] <in.paf" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";
	#if 0
	ostr << " -t : number of threads (defaults to number of cores on machine)\n";
	ostr << " -T : prefix for temporary files (default: create files in current working directory)\n";
	ostr << " -f : merge fan in (default: 64)\n";
	ostr << " -s : sort order (canonical or ba, default: canonical)\n";
	#endif

	return ostr.str();
}

static uint64_t parse(std::vector<char> const & V)
{
	std::string const s(V.begin(),V.end());
	std::istringstream istr(s);
	uint64_t u;
	istr >> u;

	if ( !(istr && istr.peek() == std::istream::traits_type::eof()) )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] unable to parse " << s << " as unsigned integer" << std::endl;
		lme.finish();
		throw lme;
	}
	return u;
}

static uint64_t parse(char const * a, char const * e)
{
	std::vector<char> V(a,e);
	return parse(V);
}

struct Meta
{
	uint64_t size;
	uint64_t od;
	uint64_t id;


	Meta()
	{

	}

	Meta(libmaus2::lf::ImpCompactHuffmanWaveletLF const & LF, uint64_t sp)
	{
		std::vector < char > VS;
		while ( LF[sp] != '\n' )
		{
			VS.push_back(LF[sp]);
			sp = LF(sp);
		}
		std::reverse(VS.begin(),VS.end());
		assert ( LF[sp] == '\n' );
		sp = LF(sp);

		std::vector < char > VOD;
		while ( LF[sp] != '\n' )
		{
			VOD.push_back(LF[sp]);
			sp = LF(sp);
		}
		std::reverse(VOD.begin(),VOD.end());
		assert ( LF[sp] == '\n' );
		sp = LF(sp);

		std::vector < char > VID;
		while ( LF[sp] != 0 )
		{
			VID.push_back(LF[sp]);
			sp = LF(sp);
		}
		std::reverse(VID.begin(),VID.end());
		assert ( LF[sp] == 0 );
		sp = LF(sp);


		id = parse(VID);
		od = parse(VOD);
		size = parse(VS);

		#if 0
		std::cerr
			<< std::string(VS.begin(),VS.end()) << "\t"
			<< std::string(VOD.begin(),VOD.end()) << "\t"
			<< std::string(VID.begin(),VID.end()) << "\n";
		std::cerr
			<< size << "\t"
			<< od << "\t"
			<< id << "\n";
		#endif
	}
};

static std::pair<char const *, char const *> arc(char const * c_ab, char const * c_ae, libmaus2::autoarray::AutoArray<char> & R)
{
	uint64_t const n = c_ae - c_ab;

	R.ensureSize(n);

	char * r = R.begin() + n;

	for ( char const * c_c = c_ab; c_c != c_ae; ++c_c )
		*(--r) = libmaus2::fastx::invertUnmapped(*c_c);

	assert ( r == R.begin() );

	return std::pair<char const *,char const *>(R.begin(), R.begin() + n);
}

struct Index
{
	typedef Index this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	libmaus2::autoarray::AutoArray<char> D;
	uint64_t numreads;
	libmaus2::lf::ImpCompactHuffmanWaveletLF::unique_ptr_type PLF;
	libmaus2::fm::FM<libmaus2::lf::ImpCompactHuffmanWaveletLF>::unique_ptr_type PFM;

	Index(std::string const & infasta) : D(), numreads(0)
	{
		std::ostringstream namestr;
		uint64_t od = 0;

		libmaus2::fastx::FastAReader FA(infasta);
		libmaus2::fastx::FastAReader::pattern_type pattern;

		for ( ; FA.getNextPatternUnlocked(pattern); ++numreads )
		{
			namestr << numreads << '\n';
			namestr << od << '\n';
			namestr << pattern.spattern.size();
			namestr.put(0);
			namestr << pattern.getShortStringId();
			namestr.put(0);

			pattern.spattern = libmaus2::fastx::remapString(libmaus2::fastx::mapString(pattern.spattern));

			D.ensureSize(od + pattern.spattern.size());

			std::copy(
				pattern.spattern.begin(),
				pattern.spattern.end(),
				D.begin() + od
			);

			od += pattern.spattern.size();
		}

		std::string const namedata = namestr.str();
		libmaus2::aio::ArrayFile<std::string::const_iterator> AF(namedata.begin(),namedata.end());

		uint64_t const numbwtthreads = 1; //libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultNumThreads();
		// uint64_t const numthreads = 1; // libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultNumThreads();
		libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions options(
			AF.getURL(),
			1*1024ull*1024ull*1024ull, /* memory 16GB */
			numbwtthreads,
			"bytestream",
			false /* bwtonly */,
			std::string("mem:tmp_"),
			std::string(), // sparse
			std::string("mem:file.bwt"),
			32 /* isa */,
			32 /* sa */
		);
		options.verbose = 0;

		// construct BWT, SA and ISA
		std::ostringstream bwtostr;
		libmaus2::suffixsort::bwtb3m::BwtMergeSortResult res = libmaus2::suffixsort::bwtb3m::BwtMergeSort::computeBwt(options,&bwtostr);

		// load LF object
		libmaus2::lf::ImpCompactHuffmanWaveletLF::unique_ptr_type TLF = res.loadLF("mem://tmp_",numbwtthreads);
		PLF = UNIQUE_PTR_MOVE(TLF);
		// uint64_t const n = PLF->W->size();

		// construct FM object
		libmaus2::fm::FM<libmaus2::lf::ImpCompactHuffmanWaveletLF>::unique_ptr_type TFM(res.loadFM("mem://tmp",numbwtthreads));
		PFM = UNIQUE_PTR_MOVE(TFM);
	}
};

int paftolas(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	libmaus2::util::LineBuffer LB(std::cin);

	//int64_t const tspace = 100;
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();
	uint64_t const tspace = arg.uniqueArgPresent("tspace") ? arg.getUnsignedNumericArg<uint64_t>("tspace") : libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
	std::cerr << "[V] using tspace " << tspace << std::endl;
	bool const small = libmaus2::dazzler::align::AlignmentFile::tspaceToSmall(tspace);
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : arginfo.getDefaultTmpFileName();

	std::string const outfn = arg[0];
	std::string const infasta = arg[1];

	Index::unique_ptr_type Pqindex(new Index(infasta));
	Index::unique_ptr_type Ptindex;

	Index const * p_qindex = Pqindex.get();
	Index const * p_tindex = Pqindex.get();

	if ( 2 < arg.size() )
	{
		std::string const infasta = arg[2];
		Index::unique_ptr_type Tindex(new Index(infasta));
		Ptindex = UNIQUE_PTR_MOVE(Tindex);
		p_tindex = Ptindex.get();
	}

	Index const & qindex = *p_qindex;
	Index const & tindex = *p_tindex;

	static char const * cigprefix = "cg:Z:";
	int64_t const cigprefixlen = strlen(cigprefix);

	libmaus2::dazzler::align::AlignmentWriterArray AWA(tmpfilebase + "_AWA",numthreads,tspace);
	#if 0
	std::string prevqsid;
	#endif

	libmaus2::autoarray::AutoArray < char > A;
	uint64_t const oAmax = 16*1024*1024;
	libmaus2::autoarray::AutoArray < uint64_t > O;

	bool running = true;

	uint64_t volatile finished = 0;
	libmaus2::parallel::PosixSpinLock finishedlock;

	libmaus2::timing::RealTimeClock rtc;
	rtc.start();

	libmaus2::autoarray::AutoArray < libmaus2::autoarray::AutoArray < std::pair<char const *, char const *> > > A_AP(numthreads);
	libmaus2::autoarray::AutoArray < libmaus2::autoarray::AutoArray < char > > A_C(numthreads);
	libmaus2::autoarray::AutoArray < libmaus2::lcs::NNPCor::unique_ptr_type > A_nnpcor(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::lcs::NNPCor::unique_ptr_type tptr(new libmaus2::lcs::NNPCor);
		A_nnpcor[i] = UNIQUE_PTR_MOVE(tptr);
	}
	libmaus2::autoarray::AutoArray < libmaus2::lcs::NNPTraceContainer > A_nnptrace(numthreads);

	libmaus2::autoarray::AutoArray < libmaus2::autoarray::AutoArray < char > > A_RC(numthreads);
	libmaus2::autoarray::AutoArray < libmaus2::autoarray::AutoArray < char > > B_RC(numthreads);

	while ( running )
	{
		uint64_t oA = 0;
		uint64_t oO = 0;

		while ( oA < oAmax )
		{
			char const * a = 0;
			char const * e = 0;
			bool const ok = LB.getline(&a,&e);

			if ( ok )
			{
				O.push(oO,oA);
				A.ensureSize(oA + e-a);
				std::copy(a,e,A.begin() + oA);
				oA += (e-a);
			}
			else
			{
				running = false;
				break;
			}
		}

		if ( oO )
		{
			uint64_t const entries = oO;
			O.push(oO,oA);

			std::cerr << "[V] " << entries << " PAF lines loaded" << std::endl;

			uint64_t volatile tfailed = 0;
			libmaus2::parallel::PosixSpinLock tfailedlock;

			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
			#endif
			for ( uint64_t zzz = 0; zzz < entries; ++zzz )
			{
				try
				{
					#if defined(_OPENMP)
					uint64_t const threadid = omp_get_thread_num();
					#else
					uint64_t const threadid = 0;
					#endif

					libmaus2::autoarray::AutoArray < std::pair<char const *, char const *> > & AP = A_AP[threadid];
					libmaus2::autoarray::AutoArray < char > & C = A_C[threadid];

					libmaus2::dazzler::align::AlignmentWriter & AW = AWA[threadid];

					char const * a = A.begin()+O[zzz];
					char const * e = A.begin()+O[zzz+1];
					char const * c = a;

					// std::cerr << std::string(a,e) << std::endl;

					uint64_t oap = 0;
					while ( c != e )
					{
						char const * ce = c;
						while ( ce != e && *ce != '\t' )
							++ce;


						AP.push(oap,std::pair<char const *,char const *>(c,ce));

						if ( ce != e )
						{
							assert ( *ce == '\t' );
							++ce;
						}

						c = ce;
					}

					if ( oap >= 12 )
					{
						// query id
						std::string const qsid(AP[0].first,AP[0].second);
						// target id
						std::string const tsid(AP[5].first,AP[5].second);

						#if 0
						if ( qsid != prevqsid )
						{
							std::cerr << qsid << std::endl;
							prevqsid = qsid;
						}
						#endif

						std::ostringstream qidstr;
						qidstr.put(0);
						qidstr << qsid;
						qidstr.put(0);
						std::string const qid = qidstr.str();

						std::ostringstream tidstr;
						tidstr.put(0);
						tidstr << tsid;
						tidstr.put(0);
						std::string const tid = tidstr.str();

						uint64_t asp, aep;
						qindex.PFM->search(qid.begin(), qid.size(), asp, aep);

						uint64_t bsp, bep;
						tindex.PFM->search(tid.begin(), tid.size(), bsp, bep);

						if ( aep-asp == 1 && bep-bsp == 1 )
						{
							Meta MQ(*(qindex.PLF),asp);
							Meta MT(*(tindex.PLF),bsp);

							int64_t cigid = -1;
							for ( uint64_t j = 12; j < oap; ++j )
								if ( AP[j].second-AP[j].first >= cigprefixlen && memcmp(cigprefix,AP[j].first,cigprefixlen) == 0 )
									cigid = j;

							uint64_t const qlen = parse(AP[1].first,AP[1].second);
							uint64_t const qstart = parse(AP[2].first,AP[2].second);
							uint64_t const qend = parse(AP[3].first,AP[3].second);

							uint64_t const tlen = parse(AP[6].first,AP[6].second);
							uint64_t const tstart = parse(AP[7].first,AP[7].second);
							uint64_t const tend = parse(AP[8].first,AP[8].second);

							if ( qlen != MQ.size || tlen != MT.size )
							{
								libmaus2::exception::LibMausException lme;
								lme.getStream() << "[E] read length is not consistent between PAF and FASTA file" << std::endl;
								lme.getStream() << "[E] qlen=" << qlen << " MQ.size=" << MQ.size << std::endl;
								lme.getStream() << "[E] tlen=" << tlen << " MT.size=" << MT.size << std::endl;
								lme.finish();
								throw lme;
							}

							bool const strandok = ( AP[4].second-AP[4].first == 1 ) && ( AP[4].first[0] == '+' || AP[4].first[0] == '-' );
							if ( !strandok )
							{
								libmaus2::exception::LibMausException lme;
								lme.getStream() << "[E] strand designator is invalid" << std::endl;
								lme.finish();
								throw lme;
							}

							bool const rc = (AP[4].first[0] == '-');

							char const * c_ab = qindex.D.begin() + MQ.od;
							char const * c_ae = c_ab + MQ.size;

							char const * c_bb = tindex.D.begin() + MT.od;
							char const * c_be = c_bb + MT.size;

							char const * c_n_ab;
							char const * c_n_ae;

							if ( rc )
							{
								std::pair<char const *, char const *> const AP = arc(c_ab,c_ae,A_RC[threadid]);
								c_n_ab = AP.first;
								c_n_ae = AP.second;
							}
							else
							{
								c_n_ab = c_ab;
								c_n_ae = c_ae;
							}

							libmaus2::lcs::AlignmentTraceContainer ATC;
							int64_t abpos;
							int64_t aepos;
							int64_t bbpos;
							int64_t bepos;

							if ( cigid < 0 )
							{
								libmaus2::lcs::NNPCor & nnpcor = *(A_nnpcor[threadid]);
								libmaus2::lcs::NNPTraceContainer & nnptrace = A_nnptrace[threadid];

								libmaus2::lcs::NNPAlignResult res = nnpcor.align(
									c_bb,c_be,tstart,
									c_n_ab,c_n_ae,rc ? (MQ.size - qend) : qstart,
									nnptrace
								);

								abpos = res.abpos;
								aepos = res.aepos;
								bbpos = res.bbpos;
								bepos = res.bepos;

								nnptrace.computeTrace(ATC);
							}
							else
							{
								abpos = tstart;
								aepos = tend;
								bbpos = rc ? (MQ.size - qend  ) : qstart;
								bepos = rc ? (MQ.size - qstart) : qend;

								char const * ita = c_bb +   abpos;
								char const * itb = c_n_ab + bbpos;

								// cigar string
								char const * cigptr = AP[cigid].first + cigprefixlen;
								uint64_t const ciglen = AP[cigid].second - cigptr;

								C.ensureSize(ciglen+1);
								std::copy(cigptr,cigptr+ciglen,C.begin());
								C[ciglen] = 0;

								libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Aop;
								size_t const numcig = libmaus2::bambam::CigarStringParser::parseCigarString(C.begin(), Aop);

								uint64_t gop = 0;
								for ( uint64_t i = 0; i < numcig; ++i )
									gop += Aop[i].second;

								if ( ATC.capacity() < gop )
									ATC.resize(gop);
								ATC.reset();

								// std::cerr << C.begin() << std::endl;

								for ( uint64_t i = 0; i < numcig; ++i )
								{
									libmaus2::bambam::BamFlagBase::bam_cigar_ops const op = static_cast<libmaus2::bambam::BamFlagBase::bam_cigar_ops>(Aop[i].first);
									uint64_t const len = Aop[i].second;

									// std::cerr << op << "," << len << std::endl;

									for ( uint64_t j = 0; j < len; ++j )
									{
										switch ( op )
										{
											case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CMATCH:
												if ( *(ita++) == *(itb++) )
												{
													*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_MATCH;
												}
												else
												{
													*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_MISMATCH;
												}
												break;
											case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
												ita++;
												*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_DEL;
												break;
											case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
												itb++;
												*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_INS;
												break;
											case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
												assert ( *ita == *itb );
												++ita;
												++itb;
												*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_MATCH;
												break;
											case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
												assert ( *ita != *itb );
												++ita;
												++itb;
												*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_MISMATCH;
												break;
											default:
											{
												libmaus2::exception::LibMausException lme;
												lme.getStream() << "[W] unsupported cigar operation " << op << " encountered" << std::endl;
												lme.finish();
												throw lme;
												break;
											}
										}
									}
								}

								assert ( ATC.ta + gop == ATC.te );

								assert ( ita == c_bb   + aepos );
								assert ( itb == c_n_ab + bepos );

								ATC.reverse();

								// std::cerr << "numcig = " << numcig << std::endl;
							}

							#if 0
							libmaus2::lcs::AlignmentPrint::printAlignmentLines(
								std::cerr,
								c_bb   + abpos, aepos-abpos,
								c_n_ab + bbpos, bepos-bbpos,
								80,
								ATC.ta,
								ATC.te,
								static_cast<uint64_t>(abpos),
								static_cast<uint64_t>(bbpos)
							);
							#endif

							libmaus2::dazzler::align::Overlap const OVL = libmaus2::dazzler::align::Overlap::computeOverlap(
								rc ? libmaus2::dazzler::align::Overlap::getInverseFlag() : 0,
								MT.id,
								MQ.id,
								abpos,
								aepos,
								bbpos,
								bepos,
								tspace,
								ATC);

							bool ok = true;
							if ( small )
							{
								for ( uint64_t i = 0; i < OVL.path.path.size(); ++i )
								{
									bool const okfirst = OVL.path.path[i].first <= std::numeric_limits<int8_t>::max();
									bool const oksecond = OVL.path.path[i].second <= std::numeric_limits<int8_t>::max();
									ok = ok && okfirst && oksecond;

									#if 0
									if ( ! okfirst )
									{
										std::cerr << OVL.path.path[i].first << " > " << std::numeric_limits<int8_t>::max() << std::endl;
									}
									if ( ! oksecond )
									{
										std::cerr << OVL.path.path[i].second << " > " << std::numeric_limits<int8_t>::max() << std::endl;
									}
									#endif
								}
							}

							// std::cerr << OVL << std::endl;

							if ( ok )
							{
								AW.put(OVL);
							}
							else
							{
								std::cerr << "[W] cannot write alignment due to excessive block values" << std::endl;
							}

						}
						else
						{
							if ( aep-asp != 1 )
								std::cerr << "[W] ignoring line " << std::string(a,e) << " because " << qsid << " is unknown" << std::endl;
							if ( bep-bsp != 1 )
								std::cerr << "[W] ignoring line " << std::string(a,e) << " because " << tsid << " is unknown" << std::endl;
						}
						// std::cerr << qsid << "\t" << tsid << std::endl;
					}
					else
					{
						std::cerr << "[W] ignoring line (too few tokens) " << std::string(a,e) << std::endl;
					}

					finishedlock.lock();
					uint64_t const lfinished = ++finished;
					finishedlock.unlock();

					if ( lfinished % 10000 == 0 )
					{
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						std::cerr << "[V] " << lfinished << " " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
					}
				}
				catch(std::exception const & ex)
				{
					tfailedlock.lock();
					tfailed = 1;
					tfailedlock.unlock();

					libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
					std::cerr << ex.what() << std::endl;
				}
			}

			if ( tfailed )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] conversion loop failed" << std::endl;
				lme.finish();
				throw lme;
			}
		}
	}
	std::cerr << "[V] finished " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

	std::cerr << "[V] merging thread las files" << std::endl;
	AWA.merge(outfn,outfn+std::string("_tmp"));
	std::cerr << "[V] merged thread las files" << std::endl;

	return EXIT_SUCCESS;
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

		return paftolas(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
