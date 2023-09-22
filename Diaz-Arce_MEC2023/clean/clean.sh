#!/bin/bash

###################################
#                                 #
#    CHANGE THIS for each PROJECT #
#                                 #
###################################

# first argument name of pool
# second argument - number of samples minus 1

clean=/share/projects/GBYP/clean_90/$1

temp=(`grep -v "#" $clean/$1.txt  |wc`)
n=(`expr $temp - 1`)

read1=(`grep -w $1 ruta_aceitu.txt |cut -f 3`)
read2=(`grep -w $1 ruta_aceitu.txt |cut -f 4`)

##################################################################33

names=(`grep -v "#" $clean/$1.txt |cut -f 1`)
barcodes=(` grep -v "#" $clean/$1.txt | cut -f 2`)
grep -v "#" $clean/$1.txt | cut -f 2 > $clean/barcodes.txt

process_radtags -i gzfastq -1 $read1 -2 $read2 -o $clean/ -e sbfI -b $clean/barcodes.txt -c -q -r -t 90

echo "sample retainedReads" >  $clean/sample_reads.txt

for i in `seq 0 $n`; do


	gunzip $clean/sample_${barcodes[i]}.1.fq.gz $clean/sample_${barcodes[i]}.rem.1.fq.gz 
	cat $clean/sample_${barcodes[i]}.1.fq $clean/sample_${barcodes[i]}.rem.1.fq > $clean/${names[i]}.fq_1
	gzip $clean/${names[i]}.fq_1

        gunzip $clean/sample_${barcodes[i]}.2.fq.gz $clean/sample_${barcodes[i]}.rem.2.fq.gz
        cat $clean/sample_${barcodes[i]}.2.fq $clean/sample_${barcodes[i]}.rem.2.fq > $clean/${names[i]}.fq_2
        gzip $clean/${names[i]}.fq_2

        echo -n ${names[i]} >> $clean/sample_reads.txt
        echo -n " " >> $clean/sample_reads.txt
        c="$(grep -c ^@ $clean/sample_${barcodes[i]}.1.fq)"
	echo -n $c >> $clean/sample_reads.txt
        echo -n " " >> $clean/sample_reads.txt
        c="$(zgrep -c ^@ $clean/${names[i]}.fq_1.gz)"
        echo  $c >> $clean/sample_reads.txt

	rm -f  $clean/sample_${barcodes[i]}.*

done

