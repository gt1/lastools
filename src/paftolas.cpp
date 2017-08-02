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
#include <libmaus2/bambam/CigarStringParser.hpp>

#include <config.h>

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " out.las in.fasta <in.paf" << std::endl;
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

int paftolas(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & /*arginfo */)
{
	libmaus2::util::LineBuffer LB(std::cin);

	//int64_t const tspace = 100;
	uint64_t const tspace = arg.uniqueArgPresent("tspace") ? arg.getUnsignedNumericArg<uint64_t>("tspace") : libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
	bool const small = libmaus2::dazzler::align::AlignmentFile::tspaceToSmall(tspace);

	std::string const outfn = arg[0];
	std::string const infasta = arg[1];

	std::ostringstream namestr;
	libmaus2::autoarray::AutoArray<char> D;
	uint64_t od = 0;
	uint64_t numreads = 0;

	std::cerr << "[V] loading reads...";
	{
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

			// std::cerr << pattern.getShortStringId() << std::endl;
		}
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
	libmaus2::lf::ImpCompactHuffmanWaveletLF::unique_ptr_type PLF = res.loadLF("mem://tmp_",numbwtthreads);
	// uint64_t const n = PLF->W->size();

	// construct FM object
	libmaus2::fm::FM<libmaus2::lf::ImpCompactHuffmanWaveletLF>::unique_ptr_type PFM(res.loadFM("mem://tmp",numbwtthreads));

	std::cerr << "done, number " << numreads << " bases " << od << std::endl;

	char const * a = 0;
	char const * e = 0;
	libmaus2::autoarray::AutoArray < std::pair<char const *, char const *> > AP;
	libmaus2::autoarray::AutoArray < char > C;
	static char const * cigprefix = "cg:Z:";
	int64_t const cigprefixlen = strlen(cigprefix);

	libmaus2::dazzler::align::AlignmentWriter AW(outfn,tspace);
	std::string prevasid;

	while ( LB.getline(&a,&e) )
	{
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
			std::string const asid(AP[0].first,AP[0].second);
			std::string const bsid(AP[5].first,AP[5].second);

			if ( asid != prevasid )
			{
				std::cerr << asid << std::endl;
				prevasid = asid;
			}

			std::ostringstream aidstr;
			aidstr.put(0);
			aidstr << asid;
			aidstr.put(0);
			std::string const aid = aidstr.str();

			std::ostringstream bidstr;
			bidstr.put(0);
			bidstr << bsid;
			bidstr.put(0);
			std::string const bid = bidstr.str();

			uint64_t asp, aep;
			PFM->search(aid.begin(), aid.size(), asp, aep);

			uint64_t bsp, bep;
			PFM->search(bid.begin(), bid.size(), bsp, bep);

			if ( aep-asp == 1 && bep-bsp == 1 )
			{
				Meta MA(*PLF,asp);
				Meta MB(*PLF,bsp);

				int64_t cigid = -1;
				for ( uint64_t j = 12; j < oap; ++j )
					if ( AP[j].second-AP[j].first >= cigprefixlen && memcmp(cigprefix,AP[j].first,cigprefixlen) == 0 )
						cigid = j;

				uint64_t const alen = parse(AP[1].first,AP[1].second);
				uint64_t const astart = parse(AP[2].first,AP[2].second);
				uint64_t const aend = parse(AP[3].first,AP[3].second);

				uint64_t const blen = parse(AP[6].first,AP[6].second);
				uint64_t const bstart = parse(AP[7].first,AP[7].second);
				uint64_t const bend = parse(AP[8].first,AP[8].second);

				if ( alen != MA.size || blen != MB.size )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] read length is not consistent between PAF and FASTA file" << std::endl;
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

				std::string const a(D.begin() + MA.od, D.begin() + MA.od + MA.size);
				std::string const b(D.begin() + MB.od, D.begin() + MB.od + MB.size);
				std::string const an = rc ? libmaus2::fastx::reverseComplementUnmapped(a) : a;
				std::string const bn = rc ? libmaus2::fastx::reverseComplementUnmapped(b) : b;

				libmaus2::lcs::AlignmentTraceContainer ATC;
				int64_t abpos;
				int64_t aepos;
				int64_t bbpos;
				int64_t bepos;

				if ( cigid < 0 )
				{
					libmaus2::lcs::NNPCor nnpcor;
					libmaus2::lcs::NNPTraceContainer nnptrace;

					libmaus2::lcs::NNPAlignResult res = nnpcor.align(
						an.begin(),an.end(),rc ? (MA.size - aend) : astart,
						b.begin(),b.end(),bstart,
						nnptrace
					);

					// std::cerr << res << std::endl;

					abpos = res.abpos;
					aepos = res.aepos;
					bbpos = res.bbpos;
					bepos = res.bepos;

					nnptrace.computeTrace(ATC);
				}
				else
				{
					abpos = rc ? (MA.size - aend  ) : astart;
					aepos = rc ? (MA.size - astart) : aend;
					bbpos = bstart;
					bepos = bend;

					std::string::const_iterator ita = an.begin() + abpos;
					std::string::const_iterator itb = b.begin()  + bbpos;

					char const * cigptr = AP[cigid].first + cigprefixlen;
					uint64_t const ciglen = AP[cigid].second - cigptr;

					C.ensureSize(ciglen+1);

					std::copy(
						cigptr,cigptr+ciglen,C.begin()
					);
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

					uint64_t calen = 0;
					uint64_t cblen = 0;

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
									calen++;
									cblen++;
									break;
								case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
									cblen++;
									itb++;
									*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_INS;
									break;
								case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
									calen++;
									ita++;
									*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_DEL;
									break;
								case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
									assert ( *ita == *itb );
									calen++;
									cblen++;
									++ita;
									++itb;
									*(--ATC.ta) = libmaus2::lcs::BaseConstants::STEP_MATCH;
									break;
								case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
									assert ( *ita != *itb );
									calen++;
									cblen++;
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

					#if 0
					std::cerr << calen << " " << cblen << std::endl;
					std::cerr << (ita - an.begin()+abpos) << " " << aepos-abpos << std::endl;
					std::cerr << (itb -  b.begin()+bbpos) << " " << bepos-bbpos << std::endl;
					#endif

					assert ( ita == an.begin() + aepos );
					assert ( itb == b.begin()  + bepos );

					ATC.reverse();

					// std::cerr << "numcig = " << numcig << std::endl;
				}

				if ( rc )
				{
					std::swap(abpos,aepos);
					abpos = MA.size - abpos;
					aepos = MA.size - aepos;

					std::swap(bbpos,bepos);
					bbpos = MB.size - bbpos;
					bepos = MB.size - bepos;

					std::reverse(ATC.ta,ATC.te);
				}

				#if 0
				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					std::cerr,
					a.begin() + abpos, aepos-abpos,
					bn.begin() + bbpos, bepos-bbpos,
					80,
					ATC.ta,
					ATC.te,
					static_cast<uint64_t>(abpos),
					static_cast<uint64_t>(bbpos)
				);
				#endif

				libmaus2::dazzler::align::Overlap const OVL = libmaus2::dazzler::align::Overlap::computeOverlap(
					rc ? libmaus2::dazzler::align::Overlap::getInverseFlag() : 0,
					MA.id,
					MB.id,
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
					std::cerr << "[W] ignoring line " << std::string(a,e) << " because " << asid << " is unknown" << std::endl;
				if ( bep-bsp != 1 )
					std::cerr << "[W] ignoring line " << std::string(a,e) << " because " << bsid << " is unknown" << std::endl;
			}
			// std::cerr << asid << "\t" << bsid << std::endl;
		}
		else
		{
			std::cerr << "[W] ignoring line (too few tokens) " << std::string(a,e) << std::endl;
		}
	}

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
