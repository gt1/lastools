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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/lcs/DalignerLocalAlignment.hpp>
#include <libmaus2/lcs/LocalAlignmentPrint.hpp>
#include <libmaus2/lcs/DalignerNP.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>

uint8_t remap(uint8_t const & a)
{
	return libmaus2::fastx::remapChar(a);
}

std::pair<uint8_t const *,uint64_t> getRead(libmaus2::autoarray::AutoArray<char> const & data, std::vector<uint64_t> const & off, uint64_t const id, bool const reverse)
{
	uint64_t const dif = off.at(id+1)-off.at(id);
	assert ( dif % 2 == 0 );
	uint64_t const padlen = dif/2;
	uint64_t const rl = padlen - 2;
	uint8_t const * base = reinterpret_cast<uint8_t const *>(data.begin() + off[id]);

	if ( reverse )
		return std::pair<uint8_t const *,uint64_t>(base + padlen + 1,rl);
	else
		return std::pair<uint8_t const *,uint64_t>(base + 1,rl);
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);

		std::string const outfn = arg[0];
		std::string const outfntmp = outfn + ".tmp";
		libmaus2::util::TempFileRemovalContainer::addTempFile(outfntmp);
		std::string const db0 = arg[1];
		std::string const db1 = arg[2];

		libmaus2::dazzler::db::DatabaseFile DB0(db0);
		DB0.computeTrimVector();
		uint64_t const n0 = DB0.size();
		libmaus2::dazzler::db::DatabaseFile DB1(db1);
		DB1.computeTrimVector();
		uint64_t const n1 = DB1.size();

		libmaus2::autoarray::AutoArray<char> Adata;
		std::vector<uint64_t> Aoff;
		DB0.decodeReadsAndReverseComplementMappedTerm(0,n0,Adata,Aoff);

		libmaus2::autoarray::AutoArray<char> Bdata;
		std::vector<uint64_t> Boff;
		DB1.decodeReadsAndReverseComplementMappedTerm(0,n1,Bdata,Boff);

		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(arg[3]);
		libmaus2::dazzler::align::AlignmentWriter::unique_ptr_type AW(new libmaus2::dazzler::align::AlignmentWriter(outfntmp,tspace,false /* no index */));

		libmaus2::lcs::DalignerNP DNP;
		libmaus2::lcs::AlignmentTraceContainer ATC;

		for ( uint64_t i = 3; i < arg.size(); ++i )
		{
			std::string const lasfn = arg[i];

			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type afile(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFile(lasfn));

			libmaus2::dazzler::align::Overlap OVL;
			while ( afile->getNextOverlap(OVL) )
			{
				int64_t const aread = OVL.aread;
				int64_t const bread = OVL.bread;

				std::pair<uint8_t const *,uint64_t> DA(getRead(Adata,Aoff,aread,false));
				std::pair<uint8_t const *,uint64_t> DB(getRead(Bdata,Boff,bread,OVL.isInverse()));

				libmaus2::dazzler::align::Overlap const OVLswapped = OVL.getSwappedPreMapped(afile->Palgn->tspace,DA.first,DA.second,DB.first,DB.second,ATC,DNP);

				AW->put(OVLswapped);

				#if 0
				#if 0
				libmaus2::lcs::DalignerLocalAlignment DLA;
				libmaus2::autoarray::AutoArray<std::pair<uint16_t,uint16_t> > Atrace(OVL.path.tlen/2);
				std::copy(
					OVL.path.path.begin(),
					OVL.path.path.begin() + OVL.path.tlen/2,
					Atrace.begin()
				);

				std::cerr << OVL << std::endl;
				libmaus2::lcs::LocalEditDistanceResult const DLALEDR = DLA.computeDenseTracePreMapped(
					reinterpret_cast<uint8_t const *>(adata),alen,
					reinterpret_cast<uint8_t const *>(bdata),blen,
					afile->Palgn->tspace,
					Atrace.begin(),
					Atrace.size(),
					OVL.path.diffs,
					OVL.path.abpos,
					OVL.path.bbpos,
					OVL.path.aepos,
					OVL.path.bepos
				);
				std::cerr << "[VDLA]" << std::endl;
				libmaus2::lcs::LocalAlignmentPrint::printAlignmentLines(
					std::cerr,
					reinterpret_cast<uint8_t const *>(adata),reinterpret_cast<uint8_t const *>(adata)+alen,
					reinterpret_cast<uint8_t const *>(bdata),reinterpret_cast<uint8_t const *>(bdata)+blen,
					80,
					DLA.ta,
					DLA.te,
					DLALEDR,
					remap,
					remap
				);
				#endif

				#if 1
				libmaus2::lcs::DalignerLocalAlignment DLA;
				libmaus2::autoarray::AutoArray<std::pair<uint16_t,uint16_t> > Atrace(OVLswapped.path.tlen/2);
				std::copy(
					OVLswapped.path.path.begin(),
					OVLswapped.path.path.begin() + OVLswapped.path.tlen/2,
					Atrace.begin()
				);

				std::cerr << OVLswapped << std::endl;

				std::pair<uint8_t const *,uint64_t> RA(getRead(Bdata,Boff,bread,false));
				std::pair<uint8_t const *,uint64_t> RB(getRead(Adata,Aoff,aread,OVLswapped.isInverse()));

				// if ( ! OVLswapped.isInverse() )
				{
					libmaus2::lcs::LocalEditDistanceResult const DLALEDR = DLA.computeDenseTracePreMapped(
						RA.first,RA.second,
						RB.first,RB.second,
						afile->Palgn->tspace,
						Atrace.begin(),
						Atrace.size(),
						OVLswapped.path.diffs,
						OVLswapped.path.abpos,
						OVLswapped.path.bbpos,
						OVLswapped.path.aepos,
						OVLswapped.path.bepos
					);
					std::cerr << "[VDLA]" << std::endl;
					libmaus2::lcs::LocalAlignmentPrint::printAlignmentLines(
						std::cerr,
						RA.first,RA.first+RA.second,
						RB.first,RB.first+RB.second,
						80,
						DLA.ta,
						DLA.te,
						DLALEDR,
						remap,
						remap
					);
				}
				#endif
				#endif
			}

			#if 0
			int64_t amin = std::numeric_limits<int64_t>::max();
			int64_t amax = std::numeric_limits<int64_t>::min();
			int64_t bmin = std::numeric_limits<int64_t>::max();
			int64_t bmax = std::numeric_limits<int64_t>::min();
			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type afile(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFile(lasfn));

			libmaus2::dazzler::align::Overlap OVL;
			while ( afile->getNextOverlap(OVL) )
			{
				amin = std::min(amin,static_cast<int64_t>(OVL.aread));
				amax = std::max(amax,static_cast<int64_t>(OVL.aread));
				bmin = std::min(bmin,static_cast<int64_t>(OVL.bread));
				bmax = std::max(bmax,static_cast<int64_t>(OVL.bread));
			}
			#endif
		}

		//AW->flush();
		AW.reset();

		libmaus2::dazzler::align::SortingOverlapOutputBuffer<libmaus2::dazzler::align::OverlapFullComparator>::sortFile(outfntmp,outfn);
		libmaus2::aio::FileRemoval::removeFile(outfntmp);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
