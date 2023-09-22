popone=/share/projects/GBYP_ND/samples_allsps-ABFTonePop.txt
pop=/share/projects/GBYP_ND/samples_allsps_398ABFT.txt
exp=gbyp
n=6

mv ./*tags.tsv.gz ./csstacks_n$n/
mv ./*models.tsv.gz ./csstacks_n$n/
cd ./csstacks_n$n/
populations -b 1 -P . -M $popone -t 8 -r 0.75 -p 5
mkdir  075$exp
mv -f batch_1.sumstats* batch_1.haplotypes.tsv  batch_1.populations.log  075$exp
grep  -v "^#" 075$exp/batch_1.sumstats.tsv | cut -f 2 | sort | uniq > tagsIn075_$exp.tsv
populations -b 1 -P . -M $pop -t 8 --plink -W tagsIn075_$exp.tsv
mv -f batch_1.sumstats* batch_1.haplotypes.tsv batch_1.hapstats* *plink* batch_1.populations.log 075$exp
cd 075$exp

echo "step tags snps" > population_stats.txt

mv batch_1.plink.ped ${exp}.ped
mv batch_1.plink.map ${exp}.map
plink --noweb --file  ${exp} --missing --out  $exp

tags="$(wc -l ${exp}.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}.map| cut -d "_" -f 1 | uniq|wc -l | cut -d " " -f 1)"
echo "catalog075 "  $tags " " $snps  >> population_stats.txt

plink --noweb --file  ${exp} --mind 0.25 --out ${exp}_mind25 --recode
plink --noweb --file  ${exp}_mind25 --missing --out  ${exp}_mind25

tags="$(wc -l ${exp}_mind25.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}_mind25.map| cut -d "_" -f 1 | uniq |wc -l| cut -d " " -f 1)"
echo "mind25 " $tags " " $snps  >> population_stats.txt

plink --noweb --file  ${exp}_mind25 --geno 0.10 --out ${exp}_mind25_geno10 --recode
plink --noweb --file  ${exp}_mind25_geno10 --missing --out  ${exp}_mind25_geno10

tags="$(wc -l ${exp}_mind25_geno10.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}_mind25_geno10.map| cut -d "_" -f 1 | uniq |wc -l| cut -d " " -f 1)"
echo "mind25_geno10 " $tags " " $snps  >> population_stats.txt

plink --noweb --file  ${exp}_mind25_geno10 --maf 0.05 --out ${exp}_mind25_geno10_maf05 --recode
plink --noweb --file  ${exp}_mind25_geno10_maf05 --missing --out  ${exp}_mind25_geno10_maf05

tags="$(wc -l ${exp}_mind25_geno10_maf05.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}_mind25_geno10_maf05.map| cut -d "_" -f 1 | uniq |wc -l| cut -d " " -f 1)"
echo "mind25_geno10_maf05 " $tags " " $snps  >> population_stats.txt

ped=./gbyp_mind25_geno10.ped
map=./gbyp_mind25_geno10.map
maf=./gbyp_mind25_geno10_maf05

# HWE considering 5 locations
grep -v cmedA $ped | grep -v emedA | grep -v wmedA | grep -v nwatlA | sed 's/^CMED./CMED/'|  sed 's/^EMED./EMED/' | sed 's/WMED./WMED/' | sed 's/SLOPE./SLOPE/' | sed 's/GOM./GOM/' > temp.ped

# HWE considering 3 locations
grep -v cmedA $ped | grep -v emedA | grep -v wmedA | grep -v nwatlA | sed 's/^CMED./MED/'|  sed 's/^EMED./MED/' | sed 's/WMED./MED/' | sed 's/SS./SLOPE/' | sed 's/GOM./GOM/' > temp.ped

#sed 's/^CMED./CMED/' $ped | sed 's/^EMED./EMED/' | sed 's/WMED./WMED/' | sed 's/SLOPE./SLOPE/' | sed 's/GOM./GOM/' > temp.ped
cp $map ./temp.map

echo "" > ./SNPs.txt
# considering 3
for f in MED SLOPE GOM; do
# considering 5
for f in CMED EMED WMED SLOPE GOM; do
        grep $f temp.ped > ${f}.ped
        cp temp.map ${f}.map
        plink --noweb --file ${f} --hwe 0.05 --out ${f}_hwe05 --recode
        sed 's/ \+/\t/g' ${f}_hwe05.map | cut -f 2 >> ./SNPs.txt
        rm -rf ${f}.*
        rm -rf ${f}_hwe05.*
done

rm -rf temp.*
#
sort ./SNPs.txt | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sed '1d' | sort -n -k 1 | grep "^5" | cut -f 2 > SNPs_notfailing_in_anypop.txt
#
sort ./SNPs.txt | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sed '1d' | sort -n -k 1 | grep "^3" | cut -f 2 > SNPs_notfailing_in_anypop_HWEMED.txt
rm -rf SNPs.txt

#
plink --noweb --file ${maf} --extract SNPs_notfailing_in_anypop.txt --out ${maf}_hwe05 --recode
plink --noweb --file ${maf}_hwe05 --missing --out ${maf}_hwe05
#
plink --noweb --file ${maf} --extract SNPs_notfailing_in_anypop_HWEMED.txt  --out ${maf}_hwe05MED --recode
plink --noweb --file ${maf}_hwe05 --missing --out ${maf}_hwe05MED

tags="$(wc -l ${exp}_mind25_geno10_maf05_hwe05.map| cut -d " " -f 1)"
snps="$(cut -f 2 ${exp}_mind25_geno10_maf05_hwe05.map| cut -d "_" -f 1 | uniq |wc -l| cut -d " " -f 1)"
echo "mind25_geno10_maf05_hwe05 " $tags " " $snps  >> population_stats.txt

# convert to bayescan
sed "s/map_file/${exp}_mind25_geno10_maf05_hwe05.map/g" ../../../ped2bayescan.spid > ped2bayescan.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile ${exp}_mind25_geno10_maf05_hwe05.ped -outputfile ${exp}_mind25_geno10_maf05_hwe05.bay -spid ped2bayescan.spid

# convert to genpop
sed "s/map_file/${exp}_mind25_geno10_maf05_hwe05.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile ${exp}_mind25_geno10_maf05_hwe05.ped -outputfile ${exp}_mind25_geno10_maf05_hwe05.genepop -spid ped2genepop.spid

###only first filtered snps per tag####
sed 's/ \+/\t/g' ${maf}_hwe05.map | sed 's/_/\t/g' | sort -n -k 2,2 -u | sed 's/\t/_/2' | cut -f 2 > oneSNPperTAG_whitelist.txt
plink --noweb --file ${maf}_hwe05 --extract oneSNPperTAG_whitelist.txt --out ${maf}_hwe05_oneSNPperTAG --recode
plink --noweb --file ${maf}_hwe05_oneSNPperTAG --missing --recode
###only first filtered snps per tag 3 locations HWE ##
sed 's/ \+/\t/g' ${maf}_hwe05MED.map | sed 's/_/\t/g' | sort -n -k 2,2 -u | sed 's/\t/_/2' | cut -f 2 > oneSNPperTAG_whitelist_HWEMED.txt
plink --noweb --file ${maf}_hwe05 --extract oneSNPperTAG_whitelist_HWEMED.txt --out ${maf}_hwe05MED_oneSNPperTAG --recode
plink --noweb --file ${maf}_hwe05MED_oneSNPperTAG --missing --recode



##convert to str##
sed "s/map_file/gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.map/g" ../../../ped2structure.spid > ped2structure.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile ${maf}_hwe05_oneSNPperTAG.ped -outputfile ${maf}_hwe05_oneSNPperTAG.str -spid ped2structure.spid

## convert to BayeScan
sed "s/map_file/gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.map/g" ../../../ped2bayescan.spid > ped2bayescan.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped -outputfile gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.bay -spid ped2bayescan.spid


###converto to str and bayeScan only 3 locations###
##convert to str##
sed "s/map_file/gbyp_mind25_geno10_maf05_hwe05MED_oneSNPperTAG.map/g" ../../../ped2structure.spid > ped2structure.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile ${maf}_hwe05MED_oneSNPperTAG.ped -outputfile ${maf}_hwe05MED_oneSNPperTAG.str -spid ped2structure.spid

## convert to BayeScan
sed "s/map_file/gbyp_mind25_geno10_maf05_hwe05MED_oneSNPperTAG.map/g" ../../../ped2bayescan.spid > ped2bayescan.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile gbyp_mind25_geno10_maf05_hwe05MED_oneSNPperTAG.ped -outputfile gbyp_mind25_geno10_maf05_hwe05MED_oneSNPperTAG.bay -spid ped2bayescan.spid



##### SOLO GOM-SS #####
grep -v med gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped > atl_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped
cp gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.map atl_mind25_geno10_maf05_hwe05_oneSNPperTAG.map


### SOLO MED ####
grep med gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped > med_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped
cp gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.map med_mind25_geno10_maf05_hwe05_oneSNPperTAG.map

##### NOW RUN ADMIXTURE with run_admixture script, PCA IN R, FST AND F3 ######

### F3 

file=gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG
mkdir F3_admixtools
cd F3_admixtools
sed 's/-9/2/1' ../${file}.ped > ./${file}.ped
sed 's/^0/1/' ../${file}.map | cut -f 1-3 > ./new_map.map
snp=(`wc -l ../${file}.map`)
snps=$(($snp*100000))
seq 1 100000 $snps > maps_column.txt
paste ./new_map.map maps_column.txt > ${file}.map
rm -rf new_map.map maps_column.txt

### first we have to convert file format from ped to eigenstrat. To do so, we can use the convertf function from ADMIXTOOLS
### we need to fill the par.PED.EIGENSTRAT file, which looks like this (withouth # symbols):
#genotypename:   gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped
#snpname:        gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.map
#indivname:      gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped
#outputformat:   EIGENSTRAT
#genotypeoutname:        gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.eigenstratgeno
#snpoutname:     gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.snp
#indivoutname:   gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.ind
#familynames:    NO

/home/natalia/bin/AdmixTools/bin/convertf -p par.PED.EIGENSTRAT

# it generates the files *.eigenstratgeno, *.ind and *.snps
# Now we need the files: 
#. list_3PopF: which includes three columns, where the first two are source populations and the third one is the tested pop. Each row represents one test
#. par_qp3PopF: 
#genotypename:   gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.eigenstratgeno
#snpname:        gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.snp
#indivname:      gbyp_429_mind25_geno10_maf05_hwe05_oneSNPperTAG.ind
#popfilename:    list_3PopF
#. pops.txt, which is like the first column of the ped file:
cut -d " " -f 1 ${file}.ped | sed 's/.MEDL/MEDL/g' | sed 's/.MEDY/MEDY/' > pops.txt
cut -c 1-29 gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.ind > temp.txt
paste -d'\0' temp.txt pops.txt > ${file}.ind
rm -rf temp.txt

/home/natalia/bin/AdmixTools/bin/qp3Pop -p par_qp3PopF > qp3Pop_results.txt

