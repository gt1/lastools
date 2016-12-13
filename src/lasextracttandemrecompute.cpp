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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/iterator.hpp>

struct MarkedInterval
{
	int64_t readid;
	libmaus2::math::IntegerInterval<int64_t> I;

	MarkedInterval() {}
	MarkedInterval(int64_t const rreadid) : readid(rreadid) {}
	MarkedInterval(int64_t const rreadid, libmaus2::math::IntegerInterval<int64_t> const & rI) : readid(rreadid), I(rI) {}

	bool operator<(MarkedInterval const & O) const
	{
		if ( readid != O.readid )
			return readid < O.readid;
		else if ( I.from != O.I.from )
			return I.from < O.I.from;
		else
			return I.to < O.I.to;
	}
};

std::ostream & operator<<(std::ostream & out, MarkedInterval const & M)
{
	out << "MarkedInterval(" << M.readid << "," << M.I << ")";
	return out;
}

struct MarkedIntervalReadIdComparator
{
	bool operator()(MarkedInterval const & A, MarkedInterval const & B) const
	{
		return A.readid < B.readid;
	}
};

bool hasOverlap(int64_t const readid, int64_t const bpos, int64_t const epos, MarkedInterval const * M, uint64_t const o)
{
	std::pair<MarkedInterval const *, MarkedInterval const *> P =
		::std::equal_range(M,M+o,MarkedInterval(readid),MarkedIntervalReadIdComparator());

	libmaus2::math::IntegerInterval<int64_t> IR(bpos,epos-1);
	for ( MarkedInterval const * p = P.first; p != P.second; ++p )
		if ( !IR.intersection(p->I).isEmpty() )
			return true;

	return false;
}


struct IntervalFileReader
{
	std::string const fn;
	mutable libmaus2::aio::InputStreamInstance ISI;
	// std::istream & istr;

	typedef libmaus2::util::ConstIterator<IntervalFileReader,MarkedInterval> const_iterator;

	IntervalFileReader(std::string const & rfn)
	: fn(rfn), ISI(fn) {}

	uint64_t size() const
	{
		ISI.clear();
		ISI.seekg(0,std::ios::end);
		uint64_t const p = ISI.tellg();
		assert ( p % 3*sizeof(uint64_t) == 0 );
		return p / (3*sizeof(uint64_t));
	}

	MarkedInterval operator[](uint64_t const i) const
	{
		ISI.clear();
		ISI.seekg(i*3*sizeof(uint64_t),std::ios::beg);
		MarkedInterval M;
		M.readid = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		M.I.from = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		M.I.to = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
		return M;
	}

	MarkedInterval get(uint64_t const i) const
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

int lasextracttandemrecompute(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & /* arginfo */)
{
	// load intervals
	IntervalFileReader IFR(arg[0]);
	IntervalFileReader::const_iterator ite = IFR.end();
	libmaus2::autoarray::AutoArray<MarkedInterval> AI;
	uint64_t AIo = 0;
	for ( IntervalFileReader::const_iterator it = IFR.begin(); it != ite; ++it )
		AI.push(AIo,*it);
	std::sort(AI.begin(),AI.begin()+AIo);

	std::cerr << "[V] number of intervals " << AIo << std::endl;

	for ( uint64_t i = 1; i < arg.size(); ++i )
	{
		std::string const infn = arg[i];

		int64_t const minaread = libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(infn);
		int64_t const maxaread = libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(infn);
		int64_t const toparead = maxaread+1;

		std::cerr << "[V] " << infn << " [" << minaread << "," << toparead << ")" << std::endl;

		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));
		libmaus2::dazzler::align::Overlap OVL;

		int64_t prevaread = std::numeric_limits<int64_t>::min();
		std::set<uint64_t> S;

		while ( Plas->getNextOverlap(OVL) )
		{
			if ( OVL.aread != prevaread )
			{
				for ( std::set<uint64_t>::const_iterator ita = S.begin(); ita != S.end(); ++ita )
				{
					libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,prevaread);
					libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,*ita);
				}
				S.clear();
				if ( prevaread >= 0 )
					std::cerr << "[V] " << prevaread << std::endl;
			}

			if (
				hasOverlap(OVL.aread,OVL.path.abpos,OVL.path.aepos,AI.begin(),AIo)
				||
				hasOverlap(OVL.bread,OVL.path.bbpos,OVL.path.bepos,AI.begin(),AIo)
			)
			{
				S.insert(OVL.bread);
				// std::cerr << "yes " << OVL.aread << " " << OVL.bread << std::endl;
			}
			else
			{
				// std::cerr << "no " << OVL.aread << " " << OVL.bread << std::endl;
			}

			prevaread = OVL.aread;
		}

		for ( std::set<uint64_t>::const_iterator ita = S.begin(); ita != S.end(); ++ita )
		{
			//std::cout << prevaread << "\t" << *ita << "\n";
			libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,prevaread);
			libmaus2::util::NumberSerialisation::serialiseNumber(std::cout,*ita);
		}
		S.clear();
		if ( prevaread >= 0 )
			std::cerr << "[V] " << prevaread << std::endl;
	}

	return EXIT_SUCCESS;
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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] intv in.las\n";
			return EXIT_SUCCESS;
		}
		else
		{
			libmaus2::util::ArgInfo const arginfo(argc,argv);
			return lasextracttandemrecompute(arg,arginfo);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
