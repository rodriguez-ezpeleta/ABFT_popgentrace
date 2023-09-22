src=/home/natalia/bin/picard_test
#path to ref_genome
ref=/share/projects/BFT-poolSeq/ref-genomes/mtDNA/annotated_mtDNA

names=(`cut -f 1 /share/projects/GBYP_ND/mapping/popmap.txt`)
in=/share/projects/GBYP_ND/clean_nc
out=/share/projects/GBYP_ND/mapping_mtDNA

## FORMAT MTDNA
#cd ${ref}
#bwa index NC_004901.2.fasta
#cd ${out}

for i in `seq 0 539`; do
	bwa mem -t 16 ${ref}/NC_004901.2.fasta $in/${names[i]}.fq_1_.gz $in/${names[i]}.fq_2_.gz > $out/${names[i]}.sam
	samtools view -S -b ${out}/${names[i]}.sam > ${out}/${names[i]}_unsorted.bam
	samtools sort ${out}/${names[i]}_unsorted.bam -o ${out}/${names[i]}_0.bam
        rm -rf ${out}/${names[i]}_unsorted.bam
        samtools index ${out}/${names[i]}_0.bam
#We skip mark and remove duplicates because we already did it usign clone_filter
#remove unmapped reads
	samtools view -b -F 0x0004 ${out}/${names[i]}_0.bam > ${out}/${names[i]}_2.bam
# remove mate-unmapped reads
        samtools view -b -F 0x0008 ${out}/${names[i]}_2.bam > ${out}/${names[i]}_3.bam
# remove not primary alignment
	samtools view -b -F 0x0100 ${out}/${names[i]}_3.bam > ${out}/${names[i]}_5.bam
	mv ${out}/${names[i]}_5.bam ${out}/${names[i]}.bam
	samtools index ${out}/${names[i]}.bam
# remove Reads2 # We save these in mapping/ folder, but we keep both reads when using stacks2 ##
	samtools view -b -F 0x0080 ${out}/${names[i]}.bam > ${out}/${names[i]}_R1.bam
	samtools index ${out}/${names[i]}_R1.bam
	rm -rf ${out}/${names[i]}_0.bam 
	rm -rf ${out}/${names[i]}_2.bam
	rm -rf ${out}/${names[i]}_3.bam
done

cd ${out}
