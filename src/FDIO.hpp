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
#if ! defined(LASTOOLS_FDIO_HPP)
#define LASTOOLS_FDIO_HPP

#include <unistd.h>
#include <cstring>
#include <libmaus2/exception/LibMausException.hpp>

struct FDIO
{
	int const fd;

	FDIO(int const rfd) : fd(rfd)
	{
	
	}
	
	void readArray(unsigned char * A, ::ssize_t l)
	{
		while ( l )
		{
			::ssize_t const r = ::read(fd,A,l);
			
			if ( r < 0 )
			{
				int const error = errno;
				
				switch ( error )
				{
					case EINTR:
					case EAGAIN:
						break;
					default:
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] read failed with error " << strerror(error) << std::endl;
						lme.finish();
						throw lme;			
					}
				}
			}
			else if ( r == 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] read failed with EOF" << std::endl;
				lme.finish();
				throw lme;			
			}
			else
			{
				assert ( r <= l );
				
				A += r;
				l -= r;
			}
		}
	}

	void writeArray(unsigned char const * A, ::ssize_t l)
	{
		while ( l )
		{
			::ssize_t const r = ::write(fd,A,l);
			
			if ( r < 0 )
			{
				int const error = errno;
				
				switch ( error )
				{
					case EINTR:
					case EAGAIN:
						break;
					default:
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] write failed with error " << strerror(error) << std::endl;
						lme.finish();
						throw lme;			
					}
				}
			}
			else if ( r == 0 )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] write failed with EOF" << std::endl;
				lme.finish();
				throw lme;			
			}
			else
			{
				assert ( r <= l );
				
				A += r;
				l -= r;
			}
		}
	}

	void writeNumber(uint64_t v)
	{
		static unsigned int const n = 8;
		unsigned char A[n];
		for ( unsigned int i = 0; i < n; ++i )
		{
			unsigned int const shift = n - i - 1;
			unsigned int const shift8 = shift << 3;
			A[i] = (v >> shift8) & 0xFFull;
		}
		writeArray(&A[0],n);
	}
	
	uint64_t readNumber()
	{
		static unsigned int const n = 8;
		unsigned char A[n];
		readArray(&A[0],n);
		uint64_t v = 0;
		for ( unsigned int i = 0; i < n; ++i )
		{
			v <<= 8;
			v |= static_cast<uint64_t>(A[i]);
		}
		return v;
	}
	
	void writeString(std::string const & s)
	{
		writeNumber(s.size());
		writeArray(reinterpret_cast<unsigned char const *>(s.c_str()),s.size());
	}
	
	std::string readString()
	{
		uint64_t const n = readNumber();
		libmaus2::autoarray::AutoArray<unsigned char> A(n,false);
		readArray(A.begin(),n);
		return std::string(A.begin(),A.end());
	}
};
#endif
