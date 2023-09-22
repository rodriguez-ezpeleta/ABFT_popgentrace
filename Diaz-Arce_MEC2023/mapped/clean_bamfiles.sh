for f in *.bam; do 
	samtools view -b -F 80 $f > $f.int.bam
	samtools view -b -F 100 $f.int.bam > $f.fin.bam
	rm -rf $f
	rm -rf $f.int.bam
	mv $f.fin.bam $f
	samtools index -b $f
done

