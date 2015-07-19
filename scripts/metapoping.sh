#!/bin/bash

#synthetic.sh <#qtl> <#markers> <type:quanti|ntrl> <optimize: Gt|NA>
#To be used only in a machine with /scratch directory

set -e

if [ $# -ne 4 ]  #Check arguments
then
	echo "Usage: $0 <#qtl> <#markers> <type:quanti|ntrl> <optimize>" 
	exit 1
fi
   
#Set arguments
QTL=$1
MARKER=$2
TYPE=$3
OPT=$4

#Some variables
CASE=${QTL}_${MARKER}
WD=${PWD} 
SCRATCH=/scratch/${USER}/${CASE}/${TYPE}/${OPT}
LOG=${WD}/LOG_${CASE}_${TYPE}_${OPT}

MTP_PREFIX=${WD}/${CASE}/subpops/${TYPE}_genotype_g50_
MTP_SUFFIX=.dat.mtp

#set up directories
mkdir -p ${SCRATCH}  
OUTDIR=${WD}/${CASE}/synthetic/${TYPE}/${OPT}
mkdir -p ${OUTDIR}
touch $LOG

#Let's go!
echo "Run started at `date`" >> $LOG
cd $SCRATCH


replicates=$(ls ${MTP_PREFIX}*.dat.mtp | wc -l)
# Iterate over replicates to get and clean the genotypes from the base populations
for ((a=1;a<=$replicates;a++))
do
	INDEX=$(printf "%03d" $a) # Indices have three digits
	FINALPOOL=${OUTDIR}/poolFile_${INDEX}
	if [ -e ${FINALPOOL} ]
	then 
		echo "${FINALPOOL} does exist" >> $LOG
	else
		echo "${FINALPOOL} does not exist" >> $LOG 
		metapop ${MTP_PREFIX}${INDEX}.dat.mtp ${WD}/templates/configMTP_${OPT}  > out
		cp poolFile ${FINALPOOL}
		echo "Replica ${INDEX} done!" >> $LOG
	fi
done

# Free scratch space
cd ${WD}
rm -r ${SCRATCH}
