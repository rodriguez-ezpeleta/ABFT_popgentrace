##Convert file from plink to BayeScan. Decide missing data, and think of filtering or not for hwe.
pedfile=../plink_mind25_geno10_maf05_hwe05
file=plink_mind25_geno10_maf05_hwe05
#grep -v nwatlA ../gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped | grep -v .medA |grep -v .medY |  grep -v nwatlY | sed 's/^.MEDL/MEDL/' > gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG_onlyL.ped
sed 's/^.MED./MED/' ${pedfile}.ped | sed 's/^GOM./GOM/' | sed 's/^SS./SS/' > ${file}.ped
cp ${pedfile}.map ./${file}.map

sed "s/map_file/${file}.map/g" /share/projects/GBYP_ND/ped2bayescan.spid > ped2bayescan.spid
java -Xmx8G -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile ${file}.ped -outputfile ${file}.bay -spid ped2bayescan.spid

## Run BayeScan. Here we could modify parameters of the chain, prior_odds, etc. Everything remains default:
bayescan ${file}.bay -threads 6

## From BayeScan outputs, we should look at *.sel using R to check for convergence. Information about outlier SNPs is in the *_fst.txt file
# Decide threshold for FDR and and posterior probability values
awk ' $2 > 0.75 ' *_fst.txt | awk ' $4 < 0.05 ' > SNPs_Prob76_qvalue05.txt

# Seek for those SNPs that are ranked into the map_file
cut -f 1 SNPs_Prob76_qvalue05.txt > SNP_numbers.txt
snps=(`wc -l ${file}.map`)
seq 1 $snps > maps_column.txt
paste maps_column.txt ${file}.map > new_map.txt

#This line searches for outlier-SNP numbers in the map file and prints their contig name
awk 'FNR==NR{a[$1];next}($1 in a){print}' SNP_numbers.txt new_map.txt | cut -f 3 > SNP_outlier_contigs.txt
paste SNP_outlier_contigs.txt SNPs_Prob76_qvalue05.txt > Renamed_SNPs_Prob76_qvalue05.txt

plink --noweb --file ${file} --extract SNP_outlier_contigs.txt --out BS_outliers_76_05 --recode

sed "s/map_file/BS_outliers_76_05.map/g" /share/projects/GBYP_ND/ped2structure.spid > ped2structure.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile BS_outliers_76_05.ped -outputfile BS_outliers_76_05.str -spid ped2structure.spid
