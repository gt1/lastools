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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/CommandContainer.hpp>

int64_t getTaskId()
{
	char const * s = getenv("SLURM_ARRAY_TASK_ID");

	if ( ! s )
	{
	        libmaus2::exception::LibMausException lme;
                lme.getStream() << "SLURM_ARRAY_TASK_ID not set" << std::endl;
                lme.finish();
                throw lme;
	}

	std::istringstream istr(s);
	uint64_t id;
	istr >> id;

	if ( istr && istr.peek() == std::istream::traits_type::eof() )
	        return id;
        else
        {
	        libmaus2::exception::LibMausException lme;
                lme.getStream() << "SLURM_ARRAY_TASK_ID: cannot parse " << s << std::endl;
                lme.finish();
                throw lme;
        }
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		libmaus2::aio::InputStreamInstance ISI(arg[0]);
		int64_t const taskid = getTaskId();
		return libmaus2::util::CommandContainer::dispatch(ISI,taskid);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << "\n";
		return EXIT_FAILURE;
	}
}
