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
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/util/ArgParser.hpp>

void loadNames(std::string const & fn, std::map<std::string,uint64_t> & M, std::map<std::string,uint64_t> & L)
{
	libmaus2::fastx::FastAReader reader(fn);
	libmaus2::fastx::FastAReader::pattern_type pattern;
	uint64_t id = 0;
	while ( reader.getNextPatternUnlocked(pattern) )
	{
		M [ pattern.getShortStringId() ] = id++;
		L [ pattern.getShortStringId() ] = pattern.getPatternLength();
	}
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( ! arg.size() )
		{
			std::cerr << "usage: " << argv[0] << " out.las in.fasta <in.bam" << std::endl;
			return EXIT_FAILURE;
		}

		libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(arginfo));
		libmaus2::bambam::BamAlignmentDecoder & dec = decwrapper->getDecoder();
		libmaus2::bambam::BamAlignment const & algn = dec.getAlignment();
		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Acigop;
		libmaus2::lcs::AlignmentTraceContainer ATC;
		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
		bool const small = libmaus2::dazzler::align::AlignmentFile::tspaceToSmall(tspace);

		assert ( !small );

		std::string const outfn = arginfo.getUnparsedRestArg(0);
		std::string const infasta = arginfo.getUnparsedRestArg(1);

		std::map<std::string,uint64_t> M;
		std::map<std::string,uint64_t> L;
		loadNames(infasta,M,L);

		libmaus2::dazzler::align::AlignmentWriter AW(outfn,tspace,false /* no index */);
		std::set<int> cigset;

		while ( dec.readAlignment() )
		{
			if ( algn.isMapped() )
			{
				assert ( M.find(algn.getName()) != M.end() );
				assert ( L.find(algn.getName()) != L.end() );
				uint64_t const readid = M.find(algn.getName())->second;
				int64_t const lrl = L.find(algn.getName())->second;

				// get cigar operations vector
				uint32_t const numcig = algn.getCigarOperations(Acigop);
				// convert to trace
				libmaus2::bambam::CigarStringParser::cigarToTrace(Acigop.begin(),Acigop.begin()+numcig,ATC,true /* ignore unknown */);
				// record operation set
				for ( uint64_t i = 0; i < numcig; ++i )
				{
					cigset.insert(Acigop[i].first);

					switch ( Acigop[i].first )
					{
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP:
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP:
							break;
						default:
						{
							libmaus2::exception::LibMausException lme;
							lme.getStream() << "[E] file contains unhadled CIGAR operation " << libmaus2::bambam::BamAlignmentDecoderBase::cigarOpToChar(Acigop[i].first) << std::endl;
							lme.finish();
							throw lme;
						}
					}
				}

				int64_t const abpos = algn.getPos();
				int64_t const aepos = abpos + algn.getReferenceLength();
				int64_t const leftclip = libmaus2::bambam::BamAlignmentDecoderBase::getFrontClipping(algn.D.begin());
				int64_t const rightclip = libmaus2::bambam::BamAlignmentDecoderBase::getBackClipping(algn.D.begin());
				int64_t const rl = libmaus2::bambam::BamAlignmentDecoderBase::getReadLengthByCigar(algn.D.begin());

				assert ( rl == lrl );

				int64_t const bbpos = leftclip;
				int64_t const bepos = rl - rightclip;

				std::pair<uint64_t,uint64_t> const SLU = ATC.getStringLengthUsed();
				assert ( static_cast<int64_t>(SLU.first) == aepos-abpos );
				assert ( static_cast<int64_t>(SLU.second) == bepos-bbpos );

				libmaus2::dazzler::align::Overlap OVL =libmaus2::dazzler::align::Overlap::computeOverlap(
					algn.isReverse() ? 1 : 0,
					algn.getRefID(),
					readid,
					abpos,
					aepos,
					bbpos,
					bepos,
					tspace,
					ATC);

				if (! OVL.bConsistent() )
				{
					std::cerr << "[E] broken alignment " << OVL << std::endl;
				}

				if ( OVL.path.pathValidSmall(small) )
					AW.put(OVL);
				else
				{
					std::cerr << "[E] not writing " << OVL << " for " << algn.getName() << std::endl;
					for ( uint64_t i = 0; i < OVL.path.path.size(); ++i )
						if ( OVL.path.path[i].first >= 256 || OVL.path.path[i].second >= 256 )
							std::cerr << "(" << OVL.path.path[i].first << "," << OVL.path.path[i].second << std::endl;
				}
			}
		}

		std::cerr << "[V] cigar ops encountered ";
		for ( std::set<int>::const_iterator ita = cigset.begin(); ita != cigset.end(); ++ita )
		{
			std::cerr << libmaus2::bambam::BamAlignmentDecoderBase::cigarOpToChar(static_cast<libmaus2::bambam::BamFlagBase::bam_cigar_ops>(*ita));
		}
		std::cerr << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
