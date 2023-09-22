for f in *fq_1.gz; do

        gunzip ${f%.fq_1.gz}.fq_1.gz
        gzip  ${f%.fq_1.gz}.fq_1
        gunzip ${f%.fq_1.gz}.fq_2.gz
        gzip  ${f%.fq_1.gz}.fq_2
        trimmomatic PE -trimlog log ${f%.fq_1.gz}.fq_1.gz ${f%.fq_1.gz}.fq_2.gz ../clean_90_nc/${f%.fq_1.gz}_1.fq.gz ../clean_90_nc/${f%.fq_1.gz}_1_un.fq.gz ../clean_90_nc/${f%.fq_1.gz}_2.fq.gz ../clean_90_nc/${f%.fq_1.gz}_2_un.fq.gz CROP:90
        rm -rf ../clean_90_nc/*_un.fq.gz
        clone_filter -1 ../clean_90_nc/${f%.fq_1.gz}_1.fq.gz  -2 ../clean_90_nc/${f%.fq_1.gz}_2.fq.gz -i gzfastq  -o . > ../clean_90_nc/${f%.fq_1.gz}.log
done

for f in *log; do  echo -n $f; echo -n " "; awk '{mul+=$2*$1} {sum+=$2} END {print mul, sum}' $f ; done >remainingReadsAfterCF.txt
