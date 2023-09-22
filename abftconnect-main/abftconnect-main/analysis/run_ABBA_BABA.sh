# From vcf file to geno format
vcf=
src=/home/natalia/downloaded_programs/genomics_general
python ${src}/VCF_processing/parseVCF.py -i ${vcf}.vcf.gz --skipIndels | gzip > ${vcf}_noIndels.geno.gz

## From geno format to allele frequencies
## aquí copiar el que haya sido el popmap del populations. El archivo "pops.txt" es un textfile con col1 = nombre ind y col2= nombre pop, igual que el popmap
cp ../ABFTL_ALA_MAC_PBF/pop_map.txt ./pops.txt

python ${src}/freq.py -g ${vcf}_noIndels.geno.gz -p GOMA -p GOML -p SSL -p SSY -p MEDL -p MEDY -p MEDA -p ALA -p MAC --popsFile pops.txt --target derived -o ${vcf}_noIndels.derFreq.tsv.gz
