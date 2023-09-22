#from /share/projects/GBYP_ND/pstacks_m3/csstacks_n6/
populations -b 1 -P . -M /share/projects/GBYP_ND/samples-onlyL_ALA.txt -t 8 --plink -p 4 -r 0.80
mkdir onlyL_ALA_p4_r80
mv batch_1.hapstats.tsv batch_1.sumstats.tsv batch_1.sumstats_summary.tsv batch_1.haplotypes.tsv batch_1.plink.map batch_1.populations.log batch_1.plink.ped ./onlyL_ALA_p4_r80
cd onlyL_ALA_p4_r80

grep -v "#" batch_1.sumstats.tsv | grep Talalunga | cut -f 1-8 | grep "-" | cut -f 2,5 | sed 's/\t/_/' > SNPs_homozygous_in_ALA.txt
plink --noweb --file plink --extract SNPs_homozygous_in_ALA.txt --out plink_HOMO_ALA --recode

echo "" > ./SNPs.txt
for f in MED SLOPE GOM Talalunga; do
        grep $f plink_HOMO_ALA.ped > ${f}.ped
        cp plink_HOMO_ALA.map ${f}.map
        plink --noweb --file ${f} --hwe 0.05 --out ${f}_hwe05 --recode
        sed 's/ \+/\t/g' ${f}_hwe05.map | cut -f 2 >> ./SNPs.txt
        rm -rf ${f}.*
        rm -rf ${f}_hwe05.*
done

rm -rf temp.*
sort ./SNPs.txt | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sed '1d' | sort -n -k 1 | grep "^4" > SNPs_notfailing_in_anypop.txt
rm -rf SNPs.txt

plink --noweb --file plink_HOMO_ALA --extract SNPs_notfailing_in_anypop.txt --out plink_HOMO_ALA_hwe05 --recode

grep -v ^# batch_1.phylip.log | cut -f 2,3 | sed 's/\t/_/' > SNPs_fixed_ALA_ABFT.txt

plink --noweb --file plink_HOMO_ALA_hwe05 --exclude SNPs_fixed_ALA_ABFT.txt --out plink_HOMO_ALA_hwe05_WOFIXED --recode
cut -f 2 plink_HOMO_ALA_hwe05_WOFIXED.map > SNPs_HOMO_ALA_hwe05_WOFIXED.txt

vcftools --vcf batch_1.vcf --snps SNPs_HOMO_ALA_hwe05_WOFIXED.txt --out Homo_ALA_hwe05_WOFIXED --recode
vcftools --vcf Homo_ALA_hwe05_WOFIXED.recode.vcf --mac 2 --out Homo_ALA_hwe05_WOFIXED_mac2 --recode
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2.recode.vcf | cut -f 3 | sed 's/_/\t/g' | sort -n -k 1,1 -u | sed 's/\t/_/1' > SNPs_HOMO_ALA_hwe05_WOFIXED_mac2_oneSNPperTag.txt

vcftools --vcf Homo_ALA_hwe05_WOFIXED_mac2.recode.vcf --snps SNPs_HOMO_ALA_hwe05_WOFIXED_mac2_oneSNPperTag.txt --out Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG --recode

mkdir dadi_dataset
cd dadi_dataset
cp ../Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf ./
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 1,2 > Gene_pos.txt
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 4 > ref_allele.txt
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 5 > alt_allele.txt
### here 97 is the column for the last alalunga
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 1-5,97 | sed 's/:/\t/g' | sed 's/\//\t/' | awk '$6==0 {print $4} AND $6==1 {print$5}' > alalunga_alleles.txt
# from column 10 to column 37 is the information for MED population genotypes on the vcf file. It extracts the first three characters of each column, which are thosewith the genotype info; and then counts the number of "zeros" and "ones"
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 10-37 | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) } OFS=" "; print $0} ' | awk '{print gsub(/0/,"")}' > ref_allele_counts_MED.txt
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 10-37 | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) } OFS=" "; print $0} ' | awk '{print gsub(/1/,"")}' > alt_allele_counts_MED.txt
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 38-58 | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) } OFS=" "; print $0} ' | awk '{print gsub(/0/,"")}' > ref_allele_counts_GOM.txt
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 38-58 | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) } OFS=" "; print $0} ' | awk '{print gsub(/1/,"")}' > alt_allele_counts_GOM.txt
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 59-93 | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) } OFS=" "; print $0} ' | awk '{print gsub(/0/,"")}' > ref_allele_counts_SLOPE.txt
grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 59-93 | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) } OFS=" "; print $0} ' | awk '{print gsub(/1/,"")}' > alt_allele_counts_SLOPE.txt
#grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 94-97 | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) } OFS=" "; print $0} ' | awk '{print gsub(/0/,"")}' > ref_allele_counts_ALA.txt
#grep -v ^# Homo_ALA_hwe05_WOFIXED_mac2_oneSNPperTAG.recode.vcf | cut -f 94-97 | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) } OFS=" "; print $0} ' | awk '{print gsub(/1/,"")}' > alt_allele_counts_ALA.txt

sed 's/^/-/' ref_allele.txt | sed 's/$/-/' > ref_column.txt
sed 's/^/-/' alalunga_alleles.txt | sed 's/$/-/' > alalunga_column.txt

echo "Pacific	Alalunga	Allele1	MED	GOM	SLOPE	Allele2	MED	GOM	SLOPE	Gene	Position" > header.txt
#echo "Pacific	Alalunga	Allele1	MED	GOM	SLOPE	ALA	Allele2	MED	GOM	SLOPE	ALA	Gene	Position" > header.txt
paste ref_column.txt alalunga_column.txt ref_allele.txt ref_allele_counts_MED.txt ref_allele_counts_GOM.txt ref_allele_counts_SLOPE.txt alt_allele.txt alt_allele_counts_MED.txt alt_allele_counts_GOM.txt alt_allele_counts_SLOPE.txt Gene_pos.txt > table.txt
#paste ref_column.txt alalunga_column.txt ref_allele.txt ref_allele_counts_MED.txt ref_allele_counts_GOM.txt ref_allele_counts_SLOPE.txt ref_allele_counts_ALA.txt alt_allele.txt alt_allele_counts_MED.txt alt_allele_counts_GOM.txt alt_allele_counts_SLOPE.txt alt_allele_counts_ALA.txt Gene_pos.txt > table.txt

cat header.txt table.txt > input_file.data

rm -rf alalunga_alleles.txt alalunga_column.txt alt_allele_counts_GOM.txt alt_allele_counts_MED.txt alt_allele_counts_SLOPE.txt alt_allele.txt Gene_pos.txt header.txt ref_allele_counts_GOM.txt ref_allele_counts_MED.txt ref_allele_counts_SLOPE.txt ref_allele.txt ref_column.txt table.txt
