src=/home/natalia/bin/picard_test
#path to ref_genome
ref=/share/projects/BFT-poolSeq/ref-genomes/PBFT_2019
in=/share/projects/GBYP_ND/clean_nc
out=/share/projects/GBYP_ND/mapping_ABFT_R1_2
mkdir ${out}
echo "names" > ${out}/names.txt
echo "total" > ${out}/total_0.txt
echo "noUnmap" > ${out}/nounmap_2.txt
echo "noMateUnmap" > ${out}/nomateunmap_3.txt
echo "noR2" > ${out}/noR2_4.txt
echo "noPrimaryAl" > ${out}/noPrimAlig_5.txt

cd $out 
cut -f 1 popmap.txt > inds.txt

for name in $(cat inds.txt); do
	bwa mem -t 16 ${ref}/BKCK_.fasta $in/${name}.fq_1_.gz $in/${name}.fq_2_.gz > $out/${names}.sam
	samtools view -S -b ${out}/${name}.sam > ${out}/${name}_unsorted.bam
	rm -rf ${out}/${name}.sam
	samtools sort ${out}/${name}_unsorted.bam -o ${out}/${name}_0.bam
        rm -rf ${out}/${name}_unsorted.bam
        samtools index ${out}/${name}_0.bam
	echo ${name} >> ${out}/names.txt
	samtools view ${out}/${name}_0.bam | wc -l >> ${out}/total_0.txt
#remove unmapped reads
	samtools view -b -F 0x0004 ${out}/${name}_0.bam > ${out}/${name}_2.bam
	samtools view ${out}/${name}_2.bam | wc -l >> ${out}/nounmap_2.txt
# remove mate-unmapped reads 
	samtools view -b -F 0x0008 ${out}/${name}_2.bam > ${out}/${name}_3.bam
	samtools view ${out}/${name}_3.bam | wc -l >> ${out}/nomateunmap_3.txt
# remove not primary alignment
	samtools view -b -F 0x0100 ${out}/${name}_3.bam > ${out}/${name}_5.bam
	samtools view ${out}/${name}_5.bam | wc -l >> ${out}/noPrimAlig_5.txt
	mv ${out}/${name}_5.bam ${out}/${name}.bam
	samtools index ${out}/${name}.bam
# remove Reads2 # We save these in mapping/ folder, but we keep both reads when using stacks2 ##
	samtools view -b -F 0x0080 ${out}/${name}.bam > ${out2}/${name}.bam
	samtools index ${out}/${name}.bam
	rm -rf ${out}/${names[i]}_0.bam.* 
	rm -rf ${out}/${names[i]}_1.bam
	rm -rf ${out}/${names[i]}_2.bam
	rm -rf ${out}/${names[i]}_3.bam
	rm -rf ${out}/${names[i]}_4.bam
done

paste ${out}/names.txt ${out}/total_0.txt > ${out}/paste1.txt
paste paste1.txt ${out}/nounmap_2.txt > ${out}/paste3.txt
paste paste3.txt ${out}/nomateunmap_3.txt > ${out}/paste4.txt
paste paste4.txt ${out}/noR2_4.txt > ${out}/paste5.txt
paste paste5.txt ${out}/noPrimAlig_5.txt > ${out}/reads_mapping.txt
rm -rf ${out}/paste2.txt
rm -rf ${out}/paste3.txt
rm -rf ${out}/paste4.txt
rm -rf ${out}/paste5.txt



