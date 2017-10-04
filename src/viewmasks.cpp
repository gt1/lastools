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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/bambam/DecoderBase.hpp>

int viewmasks(libmaus2::util::ArgParser const & arg)
{
	bool const trim = true;
	std::string const dbname = arg[0];
	libmaus2::dazzler::db::DatabaseFile DB(dbname);
	if ( trim )
		DB.computeTrimVector();
	libmaus2::aio::InputStream::unique_ptr_type Pbasestr(DB.openBaseStream());
	libmaus2::autoarray::AutoArray<char> A;

	std::vector < std::string > args;
	for ( uint64_t i = 1; i < arg.size(); ++i )
		args.push_back(arg[i]);

	for ( uint64_t a = 0; a < args.size(); ++a )
	{
		std::string const arg = args[a];
		libmaus2::dazzler::db::Track::unique_ptr_type ptrack(DB.readTrack(arg));
		
		if ( ptrack->trackannosize )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[V] track " << arg << " does not look like a mask track" << std::endl;
			lme.finish();
			throw lme;
		}
		
		libmaus2::dazzler::db::TrackAnnoInterface const & anno = ptrack->getAnno();
		
		assert ( anno.size() == DB.size()+1 );
		
		for ( uint64_t i = 0; i < DB.size(); ++i )
		{
			uint64_t const low = anno[i];
			uint64_t const high = anno[i+1];
			uint64_t const size = high-low;
			uint64_t const s = DB[i].size();
			
			unsigned char const * p = ptrack->Adata->begin() + low;
			
			assert ( size % (2*sizeof(int32_t)) == 0 );
			
			uint64_t const numintv = size / (2*sizeof(int32_t));
			
			for ( uint64_t j = 0; j < numintv; ++j )
			{
				uint64_t const from = libmaus2::bambam::DecoderBase::getLEInteger(
					p + (j * 2 + 0) * sizeof(uint32_t), sizeof(uint32_t)
				);
				uint64_t const to = libmaus2::bambam::DecoderBase::getLEInteger(
					p + (j * 2 + 1) * sizeof(uint32_t), sizeof(uint32_t)
				);

				std::cout << i << "\t" << j << "\t" << from << "\t" << to << "\t" << s << "\n";
			}
			
			// std::cerr << i << "\t" << anno[i] << std::endl;
		}
	}

	return EXIT_SUCCESS;
}

/*
 * extract a set of reads from a database and writes them as FastA files
 *
 * use trim=1 to extract from the trimmed database
 */

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		if ( arg.size() < 1 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "usage: viewmasks <db> <maskname> ..." << std::endl;
			lme.finish();
			throw lme;
		}

		return viewmasks(arg);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
