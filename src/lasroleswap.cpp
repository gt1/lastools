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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/lcs/AlignerFactory.hpp>
#include <libmaus2/util/ArgParser.hpp>

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

		if ( ! arg.size() )
		{
			std::cerr << "usage: " << argv[0] << " out.las in_a.db in_b.db in.las" << std::endl;
			return EXIT_FAILURE;
		}

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

		libmaus2::lcs::Aligner::unique_ptr_type Palgn(constructAligner());
		libmaus2::lcs::AlignmentTraceContainer ATC;

		for ( uint64_t i = 3; i < arg.size(); ++i )
		{
			std::string const lasfn = arg[i];

			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type afile(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(lasfn));

			libmaus2::dazzler::align::Overlap OVL;
			while ( afile->getNextOverlap(OVL) )
			{
				std::pair<uint8_t const *,uint64_t> DA(getRead(Adata,Aoff,OVL.aread,false));
				std::pair<uint8_t const *,uint64_t> DB(getRead(Bdata,Boff,OVL.bread,OVL.isInverse()));
				libmaus2::dazzler::align::Overlap const OVLswapped = OVL.getSwappedPreMapped(afile->Palgn->tspace,DA.first,DA.second,DB.first,DB.second,ATC,*Palgn);
				AW->put(OVLswapped);
			}
		}

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
