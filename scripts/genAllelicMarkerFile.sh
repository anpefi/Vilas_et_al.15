#!/bin/bash

#genAllelicFile.sh <#qtl>
# Generate an Allelic File suitable for quantiNemo with a specified number of loci using 256 alleles with the same frequency

MRK=$1

echo '#Allelic	values'			
echo  '################################################################'				
echo ' '				
echo '[FILE_INFO]{'				
echo '	col_locus	1'		
echo '	col_allele	2'		
echo '	col_mut_freq	3'		
echo '	col_ini_freq	4'		
echo '}	'			
echo '	'			
echo '	'			
echo '	'			
echo '	'			
echo '#locus	mut_freq	ini_freq'
for ((l=1;l<=$MRK;l++))
do
	for a in {1..256}
	do
		FREQ=$(bc -l <<< '1/256')
		echo "$l	$a	$FREQ	$FREQ"
	done
done

