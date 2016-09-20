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

Using lastobam with damapper
----------------------------

The lastobam program requires the A reads in LAS files to refer to the reference sequence database and the B reads to the read sequence database.
This is not the default output format of damapper, so make sure to call damapper (or HPC.damapper) using the the switches -C -N.
lastobam expects to see the alignments in increasing order of the B read id. This is exactly the order in which damapper outputs the alignments, so it is sufficient to concatenate
the files as an input for lastobam. A sample pipeline is thus:

```
# k-mer size used by damapper
k=20
# create dazzler DAM file for reference
fasta2DAM ref.dam ref.fasta
# call DBsplit (block size 250MB, cut off at the kmer size)
DBsplit -s250 -x${k} ref.dam
# create dazzler DAM file for the reads
fasta2DAM reads.dam reads.fasta
# call DBsplit on the reads with block size 250MB, cut nothing off
# (-x0 is important because lastobam cannot handle trimmed read databases)
DBsplit -s250 -x0 reads.dam
# run damapper
HPC.damapper -C -N -k${k} <other switches> ref.dam reads.dam | grep "^damapper" | ${SHELL}
# get output file names produced by damapper
DAMAPPEROUT=`HPC.damapper -C -N -k${k} <other switches> ref_primary_damapper.dam ref_primary_p10_damapper.dam | grep LAsort | perl -p -e "s/\s*\&.*//" | perl -p -e "s/LAsort.*?ref_/ref_/"`
# concatenate damapper output files to a single file
LAcat ${DAMAPPEROUT} >cat.las
# convert the alignments to bam format
lastobam -snone ref.dam ref.fasta reads.dam reads.fasta cat.las >out.bam
```

Using call.damapper
-------------------

call.damapper is a simplified interface to damapper and lastobam. An example call is

```
call.damapper ref.fasta <reads.fasta >reads.bam
```

call.damapper requires the programs fasta2DAM, DBsplit (both from https://github.com/thegenemyers/DAZZ_DB),
HPC.damapper, damapper (both from https://github.com/thegenemyers/DAMAPPER), 
lascat and lastobam (both in this repository) either accessible via PATH or in the same directory as call.damapper.
