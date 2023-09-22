### F3
file=plink_5POPs_geno20_mind20_geno10_maf05_hwe05_mind20_oneSNPperTAG

mkdir F3
cd F3
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

