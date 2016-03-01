# lastools
Tools for processing DALIGNER .las files

Source
------

The lastools source code is hosted on github:

	git@github.com:gt1/lastools.git

Compilation of lastools
-----------------------

lastools needs libmaus2 [https://github.com/gt1/libmaus2] . When libmaus2
is installed in ${LIBMAUSPREFIX} then lastools can be compiled and
installed in ${HOME}/lastools using

	- autoreconf -i -f
	- ./configure --with-libmaus2=${LIBMAUSPREFIX} \
		--prefix=${HOME}/lastools
	- make install

Contained programs
------------------

The command line parameters of all programs can be obtained by calling the respective program with the --help switch as only parameter.
The contained programs are:

 * bamtolas: Convert an input BAM file to the LAS format. This only converts the data which can be stored in LAS files, i.e. no auxiliary tags etc.
 * fasta2fastq: Convert an input FastA file to FastQ by adding dummy quality values ('H').
 * lasindex: Produce an index for one or more LAS files.
 * lasroleswap: Swap the roles of the A and B read in an LAS file.
 * lassort: Sort a set of LAS files and merge into a single output file.
 * lastobam: Convert a sorted input LAS file to the BAM format. This does currently not add information about unmapped reads.
 * reformatfasta: Reformats a FastA file so it follows the naming and formatting (limited column width) of the FastA expected by fasta2DB in the DAZZ_DB suite.
