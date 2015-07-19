#!/bin/bash

#genAllelicFile.sh <#qtl>

QTL=$1

echo '#Allelic	values'			
echo  '################################################################'				
echo ' '				
echo '[FILE_INFO]{'				
echo '	col_locus	1'		
echo '	col_allele	2'		
echo '	col_allelic_value	3'		
echo '	col_mut_freq	4'		
echo '	col_ini_freq	5'		
echo '}	'			
echo '	'			
echo '	'			
echo '	'			
echo '	'			
echo '#locus	allele	value	mut_freq	ini_freq'
for ((l=1;l<=$QTL;l++))
do
	for a in {1..256}
	do
		VALUE=$(bc -l <<< `echo "( ( $a - 1 ) * ( 30 / 255 ) - 15 ) / $QTL"`)
		FREQ=$(bc -l <<< '1/256')
		echo "$l	$a	$VALUE	$FREQ	$FREQ"
	done
done

