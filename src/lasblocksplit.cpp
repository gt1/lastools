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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/align/SimpleOverlapParser.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/bambam/DecoderBase.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>

int lasblocksplit(libmaus2::util::ArgParser const & arg)
{
	std::string const outprefix = arg[0];
	std::string const dbfn = arg[1];
	bool const ilist = arg.uniqueArgPresent("IL");

	std::vector<std::string> Vinfn;

	if ( ilist )
	{
		for ( uint64_t i = 2; i < arg.size(); ++i )
		{
			libmaus2::aio::InputStreamInstance ISI(arg[i]);
			while ( ISI )
			{
				std::string line;
				std::getline(ISI,line);
				if ( line.size() )
				{
					Vinfn.push_back(line);
				}
			}
		}
	}
	else
	{
		for ( uint64_t i = 2; i < arg.size(); ++i )
			Vinfn.push_back(arg[i]);
	}
	bool const splitb = arg.uniqueArgPresent("b");

	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	if ( DB.cutoff < 0 )
	{
		std::cerr << "[E] database " << dbfn << " is not split, please run DBsplit" << std::endl;
		return EXIT_FAILURE;
	}
	DB.computeTrimVector();

	std::vector < std::pair<uint64_t,uint64_t> > Vblock;
	for ( uint64_t i = 1; i <= DB.numblocks; ++i )
		Vblock.push_back(DB.getTrimmedBlockInterval(i));

	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vinfn);
	int64_t prevaread = std::numeric_limits<int64_t>::min();
	int64_t blockid = 0;
	int64_t fblockid = std::numeric_limits<int64_t>::min();
	uint64_t nal = 0;

	libmaus2::aio::OutputStreamInstance::unique_ptr_type PAW;

	for ( uint64_t z = 0; z < Vinfn.size(); ++z )
	{
		std::string const infn = Vinfn[z];

		libmaus2::dazzler::align::SimpleOverlapParser SOP(
			infn,1024*1024 /* buf size */,
			libmaus2::dazzler::align::OverlapParser::overlapparser_do_split
		);

		while ( SOP.parseNextBlock() )
		{
			libmaus2::dazzler::align::OverlapData & data = SOP.getData();

			for ( uint64_t z = 0; z < data.size(); ++z )
			{
				std::pair<uint8_t const *, uint8_t const *> const P = data.getData(z);
				int32_t const aread =
					splitb ?
						libmaus2::dazzler::align::OverlapData::getBRead(P.first)
						:
						libmaus2::dazzler::align::OverlapData::getARead(P.first)
						;

				if ( aread < prevaread )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] input file is not sorted by a-read" << std::endl;
					lme.finish();
					throw lme;
				}

				prevaread = aread;

				while ( blockid < static_cast<int64_t>(Vblock.size()) && aread >= static_cast<int64_t>(Vblock[blockid].second) )
					++blockid;

				if ( blockid == static_cast<int64_t>(Vblock.size()) )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] a-read" << aread << " is not in any database block" << std::endl;
					lme.finish();
					throw lme;
				}

				if ( blockid != fblockid )
				{
					if ( PAW )
					{
						PAW->seekp(0);

						uint64_t offset = 0;
						libmaus2::dazzler::db::OutputBase::putLittleEndianInteger8(*PAW,nal,offset);

						nal = 0;
						PAW.reset();
					}

					fblockid = blockid;
					std::ostringstream fnostr;
					fnostr << outprefix << "." << (fblockid+1) << ".las";

					libmaus2::aio::OutputStreamInstance::unique_ptr_type TAW(
						new libmaus2::aio::OutputStreamInstance(fnostr.str())
					);
					PAW = UNIQUE_PTR_MOVE(TAW);

					libmaus2::dazzler::align::AlignmentFile::serialiseHeader(*PAW,0 /* novl */,tspace);
				}

				assert ( PAW );

				PAW->write(
					reinterpret_cast<char const *>(P.first),
					P.second-P.first
				);
				nal += 1;

			}
		}
	}

	if ( PAW )
	{
		PAW->seekp(0);

		uint64_t offset = 0;
		libmaus2::dazzler::db::OutputBase::putLittleEndianInteger8(*PAW,nal,offset);

		PAW.reset();
	}

	return EXIT_SUCCESS;
}

std::string getUsage(libmaus2::util::ArgParser const & arg)
{
	std::ostringstream ostr;

	ostr << "usage: " << arg.progname << " out_prefix in.db in.las" << std::endl;
	ostr << "\n";
	ostr << "parameters:\n";

	return ostr.str();
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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 3 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << getUsage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return lasblocksplit(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
