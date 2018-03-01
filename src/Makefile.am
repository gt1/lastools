ACLOCAL_AMFLAGS=-I m4

AM_CPPFLAGS = -DDATA_PATH=\"$(datadir)\"

noinst_HEADERS = which.hpp runProgram.hpp FDIO.hpp RunInfo.hpp

MANPAGES = 

man_MANS = ${MANPAGES}

data_DATA = derive.rc

EXTRA_DIST = ${MANPAGES} derive.rc
EXTRA_PROGRAMS =

bin_PROGRAMS = lassort lasindex lastobam reformatfasta fasta2fastq lasroleswap bamtolas blasrtobam lascat lasmarkprimary call.damapper lasaligntotracepoints lasfilterquality \
	laschainsort lasmerge lascomputetandem lasextracttandemrecompute tandemaligner lasdedup lassubsample lassuffixprefix cocommand lasreduce \
	lasfilterlength paftolas commandpack commandstart commandfollowup commandprime SLURM.HPC.daligner commandderive laschecksymmetry extractreads viewmasks \
	lassimplededup lascheck laschangetspace lasproject lasshow slurmworker slurmcontrol lasdeepcheck showcdl lassort2 lasprojectsorted lassubset lasbridge lasblocksplit @SLURMTEST@ \
	lasdedupexact
dist_bin_SCRIPTS = swappedLasToBAM.sh

lasindex_SOURCES = lasindex.cpp
lasindex_LDADD = ${LIBMAUS2LIBS}
lasindex_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasindex_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lassort_SOURCES = lassort.cpp
lassort_LDADD = ${LIBMAUS2LIBS}
lassort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lassort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lastobam_SOURCES = lastobam.cpp
lastobam_LDADD = ${LIBMAUS2LIBS}
lastobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lastobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

reformatfasta_SOURCES = reformatfasta.cpp
reformatfasta_LDADD = ${LIBMAUS2LIBS}
reformatfasta_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
reformatfasta_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

fasta2fastq_SOURCES = fasta2fastq.cpp
fasta2fastq_LDADD = ${LIBMAUS2LIBS}
fasta2fastq_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
fasta2fastq_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasroleswap_SOURCES = lasroleswap.cpp
lasroleswap_LDADD = ${LIBMAUS2LIBS}
lasroleswap_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasroleswap_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamtolas_SOURCES = bamtolas.cpp
bamtolas_LDADD = ${LIBMAUS2LIBS}
bamtolas_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
bamtolas_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

blasrtobam_SOURCES = blasrtobam.cpp
blasrtobam_LDADD = ${LIBMAUS2LIBS}
blasrtobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
blasrtobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lascat_SOURCES = lascat.cpp
lascat_LDADD = ${LIBMAUS2LIBS}
lascat_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lascat_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasmarkprimary_SOURCES = lasmarkprimary.cpp
lasmarkprimary_LDADD = ${LIBMAUS2LIBS}
lasmarkprimary_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasmarkprimary_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

call_damapper_SOURCES = call.damapper.cpp
call_damapper_LDADD = ${LIBMAUS2LIBS}
call_damapper_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
call_damapper_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasaligntotracepoints_SOURCES = lasaligntotracepoints.cpp
lasaligntotracepoints_LDADD = ${LIBMAUS2LIBS}
lasaligntotracepoints_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasaligntotracepoints_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasfilterquality_SOURCES = lasfilterquality.cpp
lasfilterquality_LDADD = ${LIBMAUS2LIBS}
lasfilterquality_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasfilterquality_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

laschainsort_SOURCES = laschainsort.cpp
laschainsort_LDADD = ${LIBMAUS2LIBS}
laschainsort_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
laschainsort_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasmerge_SOURCES = lasmerge.cpp
lasmerge_LDADD = ${LIBMAUS2LIBS}
lasmerge_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasmerge_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lascomputetandem_SOURCES = lascomputetandem.cpp
lascomputetandem_LDADD = ${LIBMAUS2LIBS}
lascomputetandem_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lascomputetandem_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasextracttandemrecompute_SOURCES = lasextracttandemrecompute.cpp
lasextracttandemrecompute_LDADD = ${LIBMAUS2LIBS}
lasextracttandemrecompute_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasextracttandemrecompute_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

tandemaligner_SOURCES = tandemaligner.cpp
tandemaligner_LDADD = ${LIBMAUS2LIBS}
tandemaligner_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
tandemaligner_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasdedup_SOURCES = lasdedup.cpp
lasdedup_LDADD = ${LIBMAUS2LIBS}
lasdedup_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasdedup_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lassimplededup_SOURCES = lassimplededup.cpp
lassimplededup_LDADD = ${LIBMAUS2LIBS}
lassimplededup_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lassimplededup_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lassubsample_SOURCES = lassubsample.cpp
lassubsample_LDADD = ${LIBMAUS2LIBS}
lassubsample_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lassubsample_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

SLURM_HPC_daligner_SOURCES = SLURM.HPC.daligner.cpp
SLURM_HPC_daligner_LDADD = ${LIBMAUS2LIBS}
SLURM_HPC_daligner_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
SLURM_HPC_daligner_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lassuffixprefix_SOURCES = lassuffixprefix.cpp
lassuffixprefix_LDADD = ${LIBMAUS2LIBS}
lassuffixprefix_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lassuffixprefix_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

cocommand_SOURCES = cocommand.cpp
cocommand_LDADD = ${LIBMAUS2LIBS}
cocommand_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
cocommand_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasreduce_SOURCES = lasreduce.cpp
lasreduce_LDADD = ${LIBMAUS2LIBS}
lasreduce_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasreduce_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasfilterlength_SOURCES = lasfilterlength.cpp
lasfilterlength_LDADD = ${LIBMAUS2LIBS}
lasfilterlength_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasfilterlength_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

paftolas_SOURCES = paftolas.cpp
paftolas_LDADD = ${LIBMAUS2LIBS}
paftolas_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
paftolas_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

commandpack_SOURCES = commandpack.cpp
commandpack_LDADD = ${LIBMAUS2LIBS}
commandpack_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
commandpack_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

commandstart_SOURCES = commandstart.cpp
commandstart_LDADD = ${LIBMAUS2LIBS}
commandstart_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
commandstart_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

commandfollowup_SOURCES = commandfollowup.cpp
commandfollowup_LDADD = ${LIBMAUS2LIBS}
commandfollowup_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
commandfollowup_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

commandprime_SOURCES = commandprime.cpp
commandprime_LDADD = ${LIBMAUS2LIBS}
commandprime_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
commandprime_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

commandderive_SOURCES = commandderive.cpp
commandderive_LDADD = ${LIBMAUS2LIBS}
commandderive_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
commandderive_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

laschecksymmetry_SOURCES = laschecksymmetry.cpp
laschecksymmetry_LDADD = ${LIBMAUS2LIBS}
laschecksymmetry_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
laschecksymmetry_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

extractreads_SOURCES = extractreads.cpp
extractreads_LDADD = ${LIBMAUS2LIBS}
extractreads_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
extractreads_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

viewmasks_SOURCES = viewmasks.cpp
viewmasks_LDADD = ${LIBMAUS2LIBS}
viewmasks_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
viewmasks_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lascheck_SOURCES = lascheck.cpp
lascheck_LDADD = ${LIBMAUS2LIBS}
lascheck_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lascheck_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

laschangetspace_SOURCES = laschangetspace.cpp
laschangetspace_LDADD = ${LIBMAUS2LIBS}
laschangetspace_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
laschangetspace_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasproject_SOURCES = lasproject.cpp
lasproject_LDADD = ${LIBMAUS2LIBS}
lasproject_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasproject_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasprojectsorted_SOURCES = lasprojectsorted.cpp
lasprojectsorted_LDADD = ${LIBMAUS2LIBS}
lasprojectsorted_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasprojectsorted_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasshow_SOURCES = lasshow.cpp
lasshow_LDADD = ${LIBMAUS2LIBS}
lasshow_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasshow_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasdeepcheck_SOURCES = lasdeepcheck.cpp
lasdeepcheck_LDADD = ${LIBMAUS2LIBS}
lasdeepcheck_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasdeepcheck_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

slurmcontrol_SOURCES = slurmcontrol.cpp which.cpp runProgram.cpp
slurmcontrol_LDADD = ${LIBMAUS2LIBS}
slurmcontrol_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
slurmcontrol_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

slurmworker_SOURCES = slurmworker.cpp which.cpp runProgram.cpp
slurmworker_LDADD = ${LIBMAUS2LIBS}
slurmworker_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
slurmworker_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

showcdl_SOURCES = showcdl.cpp
showcdl_LDADD = ${LIBMAUS2LIBS}
showcdl_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
showcdl_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lassort2_SOURCES = lassort2.cpp
lassort2_LDADD = ${LIBMAUS2LIBS}
lassort2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lassort2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lassubset_SOURCES = lassubset.cpp
lassubset_LDADD = ${LIBMAUS2LIBS}
lassubset_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lassubset_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasbridge_SOURCES = lasbridge.cpp
lasbridge_LDADD = ${LIBMAUS2LIBS}
lasbridge_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasbridge_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasblocksplit_SOURCES = lasblocksplit.cpp
lasblocksplit_LDADD = ${LIBMAUS2LIBS}
lasblocksplit_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasblocksplit_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasdedupexact_SOURCES = lasdedupexact.cpp
lasdedupexact_LDADD = ${LIBMAUS2LIBS}
lasdedupexact_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS} ${AM_LDFLAGS}
lasdedupexact_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}
