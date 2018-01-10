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
#if ! defined(RUNINFO_HPP)
#define RUNINFO_HPP

#include <libmaus2/util/NumberSerialisation.hpp>

struct RunInfo
{
	uint64_t containerid;
	uint64_t subid;
	uint64_t outstart;
	uint64_t outend;
	uint64_t errstart;
	uint64_t errend;

	RunInfo()
	{

	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,containerid);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,subid);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,outstart);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,outend);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,errstart);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,errend);
		return out;
	}
};
#endif
