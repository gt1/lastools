#!/bin/bash
if [ $# -lt 3 ] ; then
	echo "usage: $0 <ref.fasta> <reads.fasta> <in1.las> ..."
	exit 1
fi

REFFASTA="$1"
shift
READSFASTA="$1"
shift

REFDB=
READSDB=
if [ -e "${REFFASTA%.fasta}.dam" ] ; then
	REFDB="${REFFASTA%.fasta}.dam"
elif [ -e "${REFFASTA%.fasta}.db" ] ; then
	REFDB="${REFFASTA%.fasta}.db"
fi
if [ -e "${READSFASTA%.fasta}.dam" ] ; then
	READSDB="${READSFASTA%.fasta}.dam"
elif [ -e "${READSFASTA%.fasta}.db" ] ; then
	READSDB="${READSFASTA%.fasta}.db"
fi

if [ -z "${REFDB}" ] ; then
	echo "$0: unable to find reference db for ${REFFASTA}"
	exit 1
fi
if [ -z "${READSDB}" ] ; then
	echo "$0: unable to find reference db for ${READSFASTA}"
	exit 1
fi

SWAPTMP=tmp_$$_swapped.las
SWAPBATMP=tmp_$$_swapped_ba.las

SCRIPTDIR=`dirname "$(readlink -f "$0")"`

"${SCRIPTDIR}/lasroleswap" ${SWAPTMP} ${READSDB} ${REFDB} $*
rm .${SWAPTMP}.bidx

"${SCRIPTDIR}/lassort" -sba ${SWAPBATMP} ${SWAPTMP}
rm -f ${SWAPTMP}
rm -f .${SWAPBATMP}.bidx

"${SCRIPTDIR}/lastobam" -snone ${REFDB} ${REFFASTA} ${READSDB} ${READSFASTA} ${SWAPBATMP}
rm -f ${SWAPBATMP}
