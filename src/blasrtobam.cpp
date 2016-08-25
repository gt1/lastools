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
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/lcs/NNPLocalAligner.hpp>
#include <libmaus2/fastx/FastAIndexGenerator.hpp>
#include "config.h"

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

std::vector < std::string > tokenize(std::string const & line)
{
	uint64_t l = 0;
	std::vector<std::string> V;
	while ( l < line.size() )
	{
		while ( l < line.size() && isspace(line[l]) )
			++l;
		uint64_t h = l;
		while ( h < line.size() && !isspace(line[h]) )
			++h;

		if ( h > l )
			V.push_back(line.substr(l,h-l));
		l = h;
	}

	return V;
}

bool parseInteger(std::string const & s, int64_t & i)
{
	std::istringstream istr(s);
	istr >> i;
	if ( ! istr )
		return false;
	if ( istr.peek() != std::istream::traits_type::eof() )
		return false;
	return true;
}

bool parseDouble(std::string const & s, double & i)
{
	std::istringstream istr(s);
	istr >> i;
	if ( ! istr )
		return false;
	if ( istr.peek() != std::istream::traits_type::eof() )
		return false;
	return true;
}

int bamtolas(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo & arginfo)
{
	std::string const reffn = arg[0];
	std::string const readsfn = arg[1];

	arginfo.insertKey("level","0");

	libmaus2::fastx::FastAIndexGenerator::generate(reffn,reffn+".fai",true /* verbose */);
	libmaus2::fastx::FastAIndex::unique_ptr_type PFAI(libmaus2::fastx::FastAIndex::load(reffn+".fai"));

	std::map<std::string,uint64_t> nameToId;

	for ( uint64_t i = 0; i < PFAI->size(); ++i )
	{
		nameToId [ (*PFAI)[i].name ] = i;
		// std::cerr << (*PFAI)[i].name << std::endl;
	}

	std::map < uint64_t, std::string> refdb;

	libmaus2::fastx::FastAReader Rin(reffn);
	libmaus2::fastx::FastAReader::pattern_type Rpattern;
	for ( uint64_t id = 0; Rin.getNextPatternUnlocked(Rpattern); ++id )
	{
		for ( uint64_t i = 0; i < Rpattern.spattern.size(); ++i )
			Rpattern.spattern[i] = toupper(Rpattern.spattern[i]);
		refdb [ id ] = Rpattern.spattern;
	}

	std::ostringstream nameostr;
	std::vector < uint64_t > Vnameoffsets;
	std::ostringstream dataostr;
	std::vector < uint64_t > Vdataoffsets;

	libmaus2::fastx::FastAReader Fin(readsfn);
	libmaus2::fastx::FastAReader::pattern_type pattern;
	uint64_t offset = 0;
	uint64_t dataoffset = 0;

	std::map < std::string, uint64_t > readNameToId;

	std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;

	std::ostringstream headerostr;
	headerostr << "@HD\tVN:1.5\n";
	for ( uint64_t i = 0; i < PFAI->size(); ++i )
		headerostr << "@SQ\tLN:" << ((*PFAI)[i]).length << "\tSN:" << (*PFAI)[i].name << "\n";
	libmaus2::bambam::BamHeader header(headerostr.str());

	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(
		libmaus2::bambam::BamBlockWriterBaseFactory::construct(header,arginfo,Pcbs));

	for ( uint64_t id = 0; Fin.getNextPatternUnlocked(pattern); ++id )
	{
		std::string const sid = pattern.getShortStringId();
		uint64_t const l = sid.size();
		nameostr.write(sid.c_str(),l+1);
		Vnameoffsets.push_back(offset);
		offset += (l+1);

		readNameToId [ sid ] = id;

		std::string rpattern = libmaus2::fastx::remapString(libmaus2::fastx::mapString(pattern.spattern));
		uint64_t const lp = rpattern.size();
		dataostr.write(rpattern.c_str(),lp+1);
		Vdataoffsets.push_back(dataoffset);
		dataoffset += (lp+1);
	}

	std::string const sreaddata = dataostr.str();
	char const * creaddata = sreaddata.c_str();
	libmaus2::lcs::NP np;

	libmaus2::autoarray::AutoArray< std::pair<libmaus2::lcs::AlignmentTraceContainer::step_type,uint64_t> > Aopblocks;
	libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Aop;
	::libmaus2::fastx::EntityBuffer<uint8_t,libmaus2::bambam::BamAlignment::D_array_alloc_type> buffer;
	libmaus2::bambam::BamSeqEncodeTable const seqenc;
	uint64_t c = 0;

	bool prevnamevalid = false;
	std::string prevname;

	while ( std::cin )
	{
		std::string line;
		std::getline(std::cin,line);

		std::vector<std::string> V = tokenize(line);

		if ( V.size() == 13 )
		{
			std::string readname = V[0];

			if ( readNameToId.find(readname) == readNameToId.end() )
			{
				if ( readname.find_last_of('/') != std::string::npos )
					readname = readname.substr(0,readname.find_last_of('/'));
			}

			if ( readNameToId.find(readname) == readNameToId.end() )
			{
				std::cerr << "[E] cannot find read id " << readname << std::endl;
				continue;
			}

			bool const secondary = prevnamevalid && (readname == prevname);
			uint64_t const secondaryflag = secondary ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSECONDARY : 0;

			uint64_t const readid = readNameToId.find(readname)->second;
			uint64_t const dataoffset = Vdataoffsets[readid];
			char const * lcreaddata = creaddata + dataoffset;

			std::string refname = V[1];

			if ( nameToId.find(refname) == nameToId.end() )
			{
				std::cerr << "[E] cannot find ref id " << refname << std::endl;
				continue;
			}

			uint64_t const refid = nameToId.find(refname)->second;

			int64_t refrc;

			if ( ! parseInteger(V[2],refrc) )
			{
				std::cerr << "[E] cannot parse " << V[2] << std::endl;
				continue;
			}

			int64_t readrc;

			if ( ! parseInteger(V[3],readrc) )
			{
				std::cerr << "[E] cannot parse " << V[3] << std::endl;
				continue;
			}

			int64_t score;

			if ( ! parseInteger(V[4],score) )
			{
				std::cerr << "[E] cannot parse " << V[4] << std::endl;
				continue;
			}

			double persim;

			if ( ! parseDouble(V[5],persim) )
			{
				std::cerr << "[E] cannot parse " << V[5] << std::endl;
				continue;
			}

			int64_t rtstart;

			if ( ! parseInteger(V[6],rtstart) )
			{
				std::cerr << "[E] cannot parse " << V[6] << std::endl;
				continue;
			}

			int64_t rtend;

			if ( ! parseInteger(V[7],rtend) )
			{
				std::cerr << "[E] cannot parse " << V[7] << std::endl;
				continue;
			}

			int64_t tlength;

			if ( ! parseInteger(V[8],tlength) )
			{
				std::cerr << "[E] cannot parse " << V[8] << std::endl;
				continue;
			}

			int64_t rqstart;

			if ( ! parseInteger(V[9],rqstart) )
			{
				std::cerr << "[E] cannot parse " << V[9] << std::endl;
				continue;
			}

			int64_t rqend;

			if ( ! parseInteger(V[10],rqend) )
			{
				std::cerr << "[E] cannot parse " << V[10] << std::endl;
				continue;
			}

			int64_t qlength;

			if ( ! parseInteger(V[11],qlength) )
			{
				std::cerr << "[E] cannot parse " << V[11] << std::endl;
				continue;
			}

			int64_t ncells;

			if ( ! parseInteger(V[12],ncells) )
			{
				std::cerr << "[E] cannot parse " << V[12] << std::endl;
				continue;
			}

			uint64_t const qstart = readrc ? (qlength-rqend) : rqstart;
			uint64_t const qend =   readrc ? (qlength-rqstart) : rqend;
			uint64_t const tstart = readrc ? (tlength-rtend) : rtstart;
			uint64_t const tend   = readrc ? (tlength-rtstart) : rtend;

			#if 0
			std::cerr << line << std::endl;

			//for ( uint64_t i = 0; i < V.size(); ++i )
			//	std::cerr << V[i] << std::endl;
			std::cerr << "readname\t" << readname << std::endl;
			std::cerr << "refname\t" << refname << std::endl;
			std::cerr << "refrc\t" << refrc << std::endl;
			std::cerr << "readrc\t" << readrc << std::endl;
			std::cerr << "score\t" << score << std::endl;
			std::cerr << "persim\t" << persim << std::endl;
			std::cerr << "rtstart\t" << rtstart << std::endl;
			std::cerr << "rtend\t" << rtend << std::endl;
			std::cerr << "tstart\t" << tstart << std::endl;
			std::cerr << "tend\t" << tend << std::endl;
			std::cerr << "tlen\t" << tlength << std::endl;
			std::cerr << "rqstart\t" << rqstart << std::endl;
			std::cerr << "rqend\t" << rqend << std::endl;
			std::cerr << "qstart\t" << qstart << std::endl;
			std::cerr << "qend\t" << qend << std::endl;
			std::cerr << "qlen\t" << qlength << std::endl;
			std::cerr << "ncells\t" << ncells << std::endl;
			#endif

			assert ( strlen(lcreaddata) == static_cast<uint64_t>(qlength) );

			assert ( rtstart >= 0 );
			assert ( rtend >= 0 );
			assert ( rtstart <= rtend );
			assert ( rqstart >= 0 );
			assert ( rqend >= 0 );
			assert ( rqstart <= rqend );

			std::string const & refuse = refdb.find(refid)->second;
			assert ( rtstart <= static_cast<int64_t>(refuse.size()) );
			assert ( rtend <= static_cast<int64_t>(refuse.size()) );
			std::string const refused = refuse.substr(tstart,tend-tstart);

			std::string const prereaduse = std::string(lcreaddata);
			assert ( prereaduse.size() == strlen(lcreaddata) );
			assert ( rqstart <= static_cast<int64_t>(prereaduse.size()) );
			assert ( rqend <= static_cast<int64_t>(prereaduse.size()) );

			#if 0
			if ( readrc )
			{
				libmaus2::lcs::NNPLocalAligner LA(6,14,1024*1024,30,50);

				std::string const readuse = libmaus2::fastx::reverseComplementUnmapped(prereaduse);

				// std::string const A = libmaus2::fastx::mapString(refused);
				std::string const A = libmaus2::fastx::mapString(refuse);
				std::string const B = libmaus2::fastx::mapString(readuse);

				std::vector< std::pair<libmaus2::lcs::NNPAlignResult,libmaus2::lcs::NNPTraceContainer::shared_ptr_type> > const V = LA.align(A.begin(),A.end(),B.begin(),B.end());

				for ( uint64_t i = 0; i < V.size(); ++i )
					std::cerr << V[i].first << std::endl;
			}
			#endif

			std::string const readuse = readrc ? libmaus2::fastx::reverseComplementUnmapped(prereaduse) : prereaduse;
			std::string const readused = readuse.substr(qstart,qend-qstart);

			np.np(refuse.begin() + tstart, refuse.begin() + tend, readused.begin(),readused.end());

			#if 0
			libmaus2::lcs::AlignmentPrint::printAlignmentLines(
				std::cerr,
				//refused.begin(),
				//refused.size(),
				refuse.begin() + tstart, tend-tstart,
				readused.begin(), qend-qstart,
				80,
				np.getTraceContainer().ta,
				np.getTraceContainer().te
			);
			#endif

			uint64_t const ncigar = libmaus2::bambam::CigarStringParser::traceToCigar(
				np.getTraceContainer(),Aopblocks,Aop,
				qstart,
				0,
				0,
				qlength-qend
			);

			libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignment(
				buffer,
				seqenc,
				readname.c_str(),
				readname.size(),
				refid,
				tstart,
				0 /* mapq */,
				(readrc ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE : 0) | secondaryflag,
				Aop.begin(),
				ncigar,
				-1 /* next ref id */,
				-1 /* next pos */,
				0 /* tlen */,
				readused.c_str() /* seq */,
				0 /* seq len */,
				readused.c_str() /* qual */
			);

			uint64_t const frontdel = libmaus2::bambam::BamAlignmentDecoderBase::getFrontDel(buffer.buffer);
			libmaus2::bambam::BamAlignmentEncoderBase::putPos(buffer.buffer, tstart + frontdel);

			Pwriter->writeBamBlock(buffer.buffer,buffer.length);

			if ( ++c % 1024 == 0 )
			{
				std::cerr << "[V] " << c << "\t" << readname << std::endl;
			}

			prevname = readname;
			prevnamevalid = true;
		}
	}

	Pwriter.reset();

	return EXIT_SUCCESS;
}

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;
	ostr << "usage: " << arg.progname << " ref.fasta reads.fasta <in.blasr\n";
	return ostr.str();
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo arginfo(argc,argv);
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

		return bamtolas(arg,arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
