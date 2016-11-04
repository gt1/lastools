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
 * laschainsort: Chain aware sorting of LAS files. Use this for sorting LAS files produced by damapper.
 * reformatfasta: Reformats a FastA file so it follows the naming and formatting (limited column width) of the FastA expected by fasta2DB in the DAZZ_DB suite.
 * call.damapper: A wrapper program for damapper and lastobam. It takes a reference and a read FastA file and produces a BAM file containing the reads as mapped by damapper.

Using lastobam with damapper
----------------------------

The lastobam program requires the A reads in LAS files to refer to the reference sequence database and the B reads to the read sequence database.
This is not the default output format of damapper, so make sure to call damapper (or HPC.damapper) using the the switches -C -N.
lastobam expects to see the alignments in increasing order of the B read id. Use laschainsort with switch -sba to obtain this order.
A sample pipeline is thus:

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
# sort
laschainsort -sba sorted.las ref.reads.las
# convert the alignments to bam format
lastobam -snone ref.dam ref.fasta reads.dam reads.fasta cat.las >out.bam
```

Using call.damapper
-------------------

call.damapper is currently broken due to incompatible changes in the damapper interface. It will be functional again soon.

call.damapper is a simplified interface to damapper and lastobam. An example call is

```
call.damapper ref.fasta <reads.fasta >reads.bam
```

call.damapper requires the programs fasta2DAM, DBsplit (both from https://github.com/thegenemyers/DAZZ_DB),
HPC.damapper, damapper (both from https://github.com/thegenemyers/DAMAPPER), 
lascat and lastobam (both in this repository) either accessible via PATH or in the same directory as call.damapper.

call.damapper passes the parameters k,t,M,e,s,n,B,T,b,v and p through to damapper. 
Please see https://github.com/thegenemyers/DAMAPPER and https://github.com/thegenemyers/DALIGNER for their meaning.
Some of the more frequently relevant ones are

 * -T: number of threads (by default -T4)
 * -M: memory limit (by default the machine's memory size)
 * -k: k-mer size (by default -k20)
 
Additionally call.damapper has the options

 * --refblocksize: reference block size used for DBsplit. The default is --refblocksize256, which sets the block size to 256MB.
 * --readblocksize: read block size used for DBsplit. The default is --readblocksize256, which sets the block size to 256MB.
 * -I: directory used for storing the reference index. By default the index is stored in the directory containing the reference FastA file.
 * -W: working directory. By default the current directory is used for storing temporary files.

Whitespace is not allowed between key and value when providing parameters (e.g. -Idir is valid, -I dir is not).

Please see https://github.com/gt1/damapper_bwt for special aux fields produced in the BAM output by lastobam/call.damapper/damapper_bwt.
