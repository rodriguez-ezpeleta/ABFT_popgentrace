pop=/share/projects/GBYP_ND/samples_onlyL99_ALA_PBF_MAC.txt
exp=plink
n=6

populations -b 1 -P . -M $pop -t 8 --plink -p 5 -r 0.75
mkdir ABFTL_ALA_MAC_PBF
mv -f batch_1.sumstats* batch_1.haplotypes.tsv batch_1.hapstats* *plink* batch_1.populations.log ABFTL_ALA_MAC_PBF
cd ABFTL_ALA_MAC_PBF

echo "step tags snps" > population_stats.txt

cp batch_1.plink.ped ${exp}.ped
cp batch_1.plink.map ${exp}.map
plink --noweb --file ${exp} --missing --out ${exp}

tags="$(wc -l ${exp}.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}.map| cut -d "_" -f 1 | uniq|wc -l | cut -d " " -f 1)"
echo "catalog075 "  $tags " " $snps  >> population_stats.txt

plink --noweb --file  plink --mind 0.20 --out ${exp}_mind20 --recode
plink --noweb --file  ${exp}_mind20 --missing --out  ${exp}_mind20

tags="$(wc -l ${exp}_mind20.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}_mind20.map| cut -d "_" -f 1 | uniq |wc -l| cut -d " " -f 1)"
echo "mind25 " $tags " " $snps  >> population_stats.txt

plink --noweb --file  ${exp}_mind20 --geno 0.10 --out ${exp}_mind20_geno10 --recode
plink --noweb --file  ${exp}_mind20_geno10 --missing --out  ${exp}_mind20_geno10

tags="$(wc -l ${exp}_mind20_geno10.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}_mind20_geno10.map| cut -d "_" -f 1 | uniq |wc -l| cut -d " " -f 1)"
echo "mind20_geno10 " $tags " " $snps  >> population_stats.txt

plink --noweb --file  ${exp}_mind20_geno10 --maf 0.01 --out ${exp}_mind20_geno10_maf01 --recode
plink --noweb --file  ${exp}_mind20_geno10_maf01 --missing --out  ${exp}_mind20_geno10_maf01

tags="$(wc -l ${exp}_mind20_geno10_maf01.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}_mind20_geno10_maf01.map| cut -d "_" -f 1 | uniq |wc -l| cut -d " " -f 1)"
echo "mind20_geno10_maf01 " $tags " " $snps  >> population_stats.txt

ped=./${exp}_mind20_geno10.ped
map=./${exp}_mind20_geno10.map
maf=./${exp}_mind20_geno10_maf01

grep -v Tori $ped | grep -v Tmac | grep -v Tala > temp.ped
cp $map ./temp.map

echo "" > ./SNPs.txt
for f in MEDL SSL GOML; do
        grep $f temp.ped > ${f}.ped
        cp temp.map ${f}.map
        plink --noweb --file ${f} --hwe 0.05 --out ${f}_hwe05 --recode
        sed 's/ \+/\t/g' ${f}_hwe05.map | cut -f 2 >> ./SNPs.txt
        rm -rf ${f}.*
        rm -rf ${f}_hwe05.*
done

rm -rf temp.*
sort ./SNPs.txt | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sed '1d' | sort -n -k 1 | grep "^3" | cut -f 2 > SNPs_notfailing_in_anypop.txt
rm -rf SNPs.txt

plink --noweb --file ${maf} --extract SNPs_notfailing_in_anypop.txt --out ${maf}_hwe05 --recode
plink --noweb --file ${maf}_hwe05 --missing --out ${maf}_hwe05

tags="$(wc -l ${maf}_hwe05.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${maf}_hwe05.map| cut -d "_" -f 1 | uniq |wc -l| cut -d " " -f 1)"
echo "mind20_geno10_maf01_hwe05 " $tags " " $snps  >> population_stats.txt

plink --noweb --file ${maf}_hwe05 --geno 0.05 --out ${maf}_hwe05_geno05 --recode

#cut -f 2 ${maf}_hwe05.map | sed 's/_/\t/g' > whitelist_for_treemix.txt
cut -f 2 ${maf}_hwe05_geno05.map | sed 's/_/\t/g' > whitelist_for_treemix_geno05.txt
cut -d " " -f 1 ${maf}_hwe05.ped > col_2.txt
cut -d " " -f 2 ${maf}_hwe05.ped > col_1.txt
paste col_1.txt col_2.txt > pop_map.txt
rm -rf col_1.txt
rm -rf col_2.txt

cd ..

populations -b 1 -P . -M ./ABFTL_ALA_MAC_PBF/pop_map.txt -t 8 --treemix -W ./ABFTL_ALA_MAC_PBF/whitelist_for_treemix.txt

mkdir ABFTL_ALA_MAC_PBF_treemix
mv -f batch_1.sumstats* batch_1.haplotypes.tsv batch_1.hapstats* batch_1.populations.log batch_1.treemix* ABFTL_ALA_MAC_PBF_treemix
cd ABFTL_ALA_MAC_PBF_treemix

sed -i '1d' batch_1.treemix
gzip batch_1.treemix
for f in `seq 0 10`; do 
treemix -i batch_1.treemix.gz -m $f -o stem_m${f}
done 

cd ..

populations -b 1 -P . -M ./ABFTL_ALA_MAC_PBF/pop_map.txt -t 8 --treemix -W ./ABFTL_ALA_MAC_PBF/whitelist_for_treemix_geno05.txt
mkdir ABFTL_ALA_MAC_PBF_treemix_geno05
mv -f batch_1.sumstats* batch_1.haplotypes.tsv batch_1.hapstats* batch_1.populations.log batch_1.treemix* ABFTL_ALA_MAC_PBF_treemix_geno05
cd ABFTL_ALA_MAC_PBF_treemix_geno05
sed -i '1d' batch_1.treemix
gzip batch_1.treemix
for f in `seq 0 10`; do
treemix -i batch_1.treemix.gz -m $f -o stem_m${f}
done


##### LOOK AT TREEMIX RESULTS: 
#####	A)EXPLORE LIKELIHOODS FROM *.lik FILES, WHICH ARE LIKELIHOOD ASSOCIATED TO EACH TREE + MIGRATION EVENTS ######
#####		 for f in `seq 1 10`; do cat stem_m${f}.llik; done	
#####	B)PLOT MIGRATION EVENTS IN R
#####


###### ABBA / BABA ######

cd ../

populations -b 1 -P . -M ./ABFTL_ALA_MAC_PBF/pop_map.txt -t 8 --vcf -W ./ABFTL_ALA_MAC_PBF/whitelist_for_treemix.txt
mkdir ABFTL_ALA_MAC_PBF_ABBABABA
mv -f batch_1.sumstats* batch_1.haplotypes.tsv batch_1.hapstats* batch_1.populations.log batch_1.vcf ABFTL_ALA_MAC_PBF_ABBABABA
cd ABFTL_ALA_MAC_PBF_ABBABABA
gzip batch_1.vcf
path=/home/natalia/downloaded_programs/genomics_general
python $path/VCF_processing/parseVCF.py -i batch_1.vcf.gz --skipIndels | gzip > noIndels.geno.gz
python ${path}/freq.py -g noIndels.geno.gz -p GOML -p MEDL -p SSL -p Talalunga -p Tmaccoyii --popsFile ../ABFTL_ALA_MAC_PBF/pop_map.txt --target derived -o noIndels.derFreq.tsv.gz
