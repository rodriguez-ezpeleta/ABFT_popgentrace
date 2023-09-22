for i in /share/projects/GBYP_ND/mapping/*.sam; do 
	samtools view -S -b $i > $i.bam
	samtools sort $i.bam -o $i.sorted.bam
	samtools index $i.sorted.bam
	rm -rf $i.bam
done
