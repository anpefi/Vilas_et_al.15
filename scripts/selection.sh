#!/bin/bash

#selection.sh <#qtl> <#markers> <type:quanti|ntrl> <optimize: Gt|NA>
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
LOG=${WD}/LOG_sel_${CASE}_${TYPE}_${OPT}


#set up directories  
SYN_DIR=${WD}/${CASE}/synthetic/${TYPE}/${OPT}
mkdir -p ${WD}/${CASE}/stats
mkdir -p ${SCRATCH}
touch $LOG

#Let's go!
echo "Run started at `date`" >> $LOG
cd $SCRATCH

# Generate synthetic pops from poolfiles


replicates=$(ls ${SYN_DIR}/poolFile_* | wc -l)
# Iterate over replicates to get and clean the genotypes from the base populations
	echo "rep   gen  PM_A    Gen_Val   HS_S   HS_N   A_S   A_N" > data
for ((a=1;a<=$replicates;a++))
do
	INDEX=$(printf "%03d" $a) # Indices have three digits

	NEU_GEN=${WD}/${CASE}/subpops/ntrl_genotype_g50_${INDEX}.dat
	QTL_GEN=${WD}/${CASE}/subpops/quanti_genotype_g50_${INDEX}.dat
	POOLFILE=${SYN_DIR}/poolFile_${INDEX}
	mkdir -p ${SCRATCH}/${INDEX}

	python ${WD}/scripts/synthetic.py ${QTL_GEN} ${NEU_GEN} ${POOLFILE} ${SCRATCH}/${INDEX}
	
	sed '1,15d' ${WD}/${CASE}/allelic_file.txt | cut -f3 > ${SCRATCH}/efecto
	od -vAn -N4 -tu4 </dev/urandom > seedfile

	for ((b=1;b<=100;b++))
	do
		REP=$(printf "%03d" $b) 
		echo "${QTL} 256 3" > ${SCRATCH}/quanti_genotype
		echo "${MARKER} 256 3" > ${SCRATCH}/ntrl_genotype
		cut -d' ' -f2- ${SCRATCH}/${INDEX}/ntrl_genotype_${REP}.dat >> ${SCRATCH}/ntrl_genotype
		cut -d' ' -f2-  ${SCRATCH}/${INDEX}/quanti_genotype_${REP}.dat >> ${SCRATCH}/quanti_genotype
				
		${WD}/scripts/selection >> out < ${WD}/templates/params.txt

		#sed -i "s/REP/$b/g" data  # Esto es muy costoso en tiempo y no merece la pena
	done
	echo "Finish replicate ${INDEX}" >> ${LOG}
done

# SUMMARIZE
Rscript ${WD}/scripts/summarize.R
mv data ${WD}/${CASE}/stats/data.${TYPE}.${OPT}
mv averages ${WD}/${CASE}/stats/averages.${TYPE}.${OPT}
mv errors ${WD}/${CASE}/stats/errors.${TYPE}.${OPT}

# Free scratch space
rm -r ${SCRATCH}
echo "Run finished at `date`" >> ${LOG}
