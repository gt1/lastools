# lastools
Tools for processing DALIGNER .las files

Source
------

The lastools source code is hosted on github:

	git@github.com:gt1/lastools.git

Compilation of lastools
---------------------

lastools needs libmaus2 [https://github.com/gt1/libmaus2] . When libmaus2
is installed in ${LIBMAUSPREFIX} then lastools can be compiled and
installed in ${HOME}/lastools using

	- autoreconf -i -f
	- ./configure --with-libmaus2=${LIBMAUSPREFIX} \
		--prefix=${HOME}/lastools
	- make install
