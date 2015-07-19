#!/bin/bash

#preMetapop.sh <#qtl> <#markers>
#To be used only in a machine with /scratch directory

set -e

if [ $# -ne 2 ]  #Check arguments
then
	echo "Usage: $0 <#qtl> <#markers> "
	exit 1
fi
   
#Set arguments
QTL=$1
MARKER=$2

#Some variables
CASE=${QTL}_${MARKER}
WD=${PWD} 
SCRATCH=/scratch/${USER}/${CASE}
ALL_VAR=$(bc -l <<< `echo 0.625/${QTL}`)
LOG=${WD}/LOG_${CASE}_`hostname`
GENBASE=13000

#set up directories
mkdir -p ${SCRATCH}  
mkdir -p ${WD}/${CASE}/basepop
mkdir -p ${WD}/${CASE}/subpops
mkdir -p ${SCRATCH}/basepop
mkdir -p ${SCRATCH}/subpops
touch $LOG

#Let's go!
echo "Run started at `date`" >> $LOG
cd $SCRATCH

# Use templates to generate case-specific configuration files
sed -e "s/SETquanti_loci/${QTL}/"  -e "s/SETquanti_allelic_var/${ALL_VAR}/" -e "s/SETgens/${GENBASE}/" -e "s/SETntrl_loci/${MARKER}/" ${WD}/templates/quantiNemo.ini > ${SCRATCH}/quantiNemo.ini #Adapts quantiNemo.ini to this case
$WD/scripts/genAllelicFile.sh $QTL > ${SCRATCH}/allelic_file.txt  #creates allelic_file.txt required by QuantiNemo
$WD/scripts/genAllelicMarkerFile.sh $MARKER > ${SCRATCH}/ntrl_allelic_file.txt  #creates ntrl_allelic_file.txt required by QuantiNemo
cp ${SCRATCH}/*allelic_file.txt ${WD}/${CASE}/ #Make a copy of the allelic files into the working directory to be able to check them

# Run Quantinemo to produce basepopulations
quantiNEMO >> out.${CASE} < quantiNemo.ini
echo "Base population finished at `date`" >> $LOG

#Quick path to the folders where genotypes are being generated
QUA=${SCRATCH}/simulation/quanti_genotype
NTR=${SCRATCH}/simulation/ntrl_genotype

# Create array with subpopulations size 
SIZE_SUBPOPS=(2000 100 100 100 100 100 100 100)
# Create file labelling the subpopultion number of each individuals
touch labels
for ((i=0; i<${#SIZE_SUBPOPS[@]};i++)); do for ((j=0; j<${SIZE_SUBPOPS[$i]}; j++)); do echo $(expr ${i} + 1)>>labels; done; done
TOTAL_SIZE=0; for i in ${SIZE_SUBPOPS[@]}; do TOTAL_SIZE=$(( $TOTAL_SIZE + $i ));  done;

# USe template for the second run of quantiNEMO (subpopulations)
sed -e "s/SETquanti_loci/${QTL}/"  -e "s/SETquanti_allelic_var/${ALL_VAR}/" -e "s/SETntrl_loci/${MARKER}/" ${WD}/templates/quantiNemoSubpop.ini > ${SCRATCH}/quantiNemo.ini #Adapts quantiNemoSubpop.ini to this case

# A dirty way to get the number of replicates
replicates=$(ls ${QUA}/simulation_g${GENBASE}_r*.dat | wc -l) 

# Iterate over replicates to get and clean the genotypes from the base populations
for ((a=1;a<=$replicates;a++))
do
	INDEX=$(printf "%03d" $a) # Indices have three digits

	#Deleting the header and loci names to the final base population storing location
	sed "1,$(($QTL+1))d" ${QUA}/simulation_g${GENBASE}_r${INDEX}.dat | cut -d' ' -f2-$(($QTL+1)) > ${SCRATCH}/basepop/quanti_genotype_${INDEX}
	sed "1,$(($MARKER+1))d" ${NTR}/simulation_g${GENBASE}_r${INDEX}.dat | cut -d' ' -f2-$(($MARKER+1)) > ${SCRATCH}/basepop/ntrl_genotype_${INDEX}
done
rm -r ${SCRATCH}/simulation #No need it anymore, saves disk space.

# Iterate over replicates to run the subpopulations
for ((a=1;a<=$replicates;a++))
do
	INDEX=$(printf "%03d" $a)


	#Create initial genotypes files from those in basepop by sampling and subdividing
	#New header for the genotypes files to be load in quantiNEMO
	echo "${#SIZE_SUBPOPS[@]}	${QTL}	256	3" > ${SCRATCH}/quanti_ini_genotype.dat
	echo "${#SIZE_SUBPOPS[@]}       ${MARKER}  256     3" > ${SCRATCH}/ntrl_ini_genotype.dat
	#Create new dumb loci names (as required by be loaded by quantiNEMO)
	for((i=0; i<${QTL};i++)); do echo l${i} >> ${SCRATCH}/quanti_ini_genotype.dat; done
	for((i=0; i<${MARKER};i++)); do echo l${i} >> ${SCRATCH}/ntrl_ini_genotype.dat; done
	#Get the sample of the N first genotypes for and give them a supopulation number (stored in labels)
	paste labels <(head -${TOTAL_SIZE} ${SCRATCH}/basepop/quanti_genotype_${INDEX}) >> ${SCRATCH}/quanti_ini_genotype.dat
	paste labels <(head -${TOTAL_SIZE} ${SCRATCH}/basepop/ntrl_genotype_${INDEX}) >> ${SCRATCH}/ntrl_ini_genotype.dat

	## RUNNING SUBPOPULATIONS
	quantiNEMO >> ${SCRATCH}/out.subpop.${CASE} 

	## SAVING RELEVANT FILES for each Replicate
	cp ${SCRATCH}/ntrl_ini_genotype.dat ${SCRATCH}/subpops/ntrl_ini_genotype_${INDEX}.dat
	cp ${SCRATCH}/quanti_ini_genotype.dat ${SCRATCH}/subpops/quanti_ini_genotype_${INDEX}.dat 
	cp ${SCRATCH}/simulation2/ntrl_genotype/simulation_g50.dat ${SCRATCH}/subpops/ntrl_genotype_g50_${INDEX}.dat
	cp ${SCRATCH}/simulation2/quanti_genotype/simulation_g50.dat  ${SCRATCH}/subpops/quanti_genotype_g50_${INDEX}.dat
	
	# CONVERT TO METAPOP
	python ${WD}/scripts/fs2mp.py ${SCRATCH}/subpops/ntrl_genotype_g50_${INDEX}.dat
	python ${WD}/scripts/fs2mp.py ${SCRATCH}/subpops/quanti_genotype_g50_${INDEX}.dat
	#Metapop file is subpops/ntrl_genotype_g50_${INDEX}.dat.mtp	
done

##COPYING BACK RELEVANT FILES (Removing headers and loci names in genotype files to save space)
echo "Retrieving relevant files at `date`" >> $LOG
#basepop files already trimmed, just copy to head
cp -r ${SCRATCH}/basepop ${WD}/${CASE}/
# metapop files, copy as are
cp -r ${SCRATCH}/subpops/*.mtp ${WD}/${CASE}/subpops/
# Trimming the subpops genotype files saving them in the WD
for ((a=1;a<=$replicates;a++))
do
	INDEX=$(printf "%03d" $a)
	sed "1,$(($QTL+1))d"  ${SCRATCH}/subpops/quanti_ini_genotype_${INDEX}.dat >  ${WD}/${CASE}/subpops/quanti_ini_genotype_${INDEX}.dat
	sed "1,$(($QTL+1))d"  ${SCRATCH}/subpops/quanti_genotype_g50_${INDEX}.dat >  ${WD}/${CASE}/subpops/quanti_genotype_g50_${INDEX}.dat
	sed "1,$(($MARKER+1))d"  ${SCRATCH}/subpops/ntrl_ini_genotype_${INDEX}.dat >  ${WD}/${CASE}/subpops/ntrl_ini_genotype_${INDEX}.dat
	sed "1,$(($MARKER+1))d"  ${SCRATCH}/subpops/ntrl_genotype_g50_${INDEX}.dat >  ${WD}/${CASE}/subpops/ntrl_genotype_g50_${INDEX}.dat
done		
echo "Copy finished at `date`" >> $LOG

# Free scratch space
rm -r ${SCRATCH}
