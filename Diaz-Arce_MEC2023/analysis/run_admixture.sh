file=gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG
outfold=admixture
K=2
snp=(`wc -l ${file}.map`)
snps=$(($snp*100000))

# multiply number of snps per 100,000
#snps=884600000

seq 1 100000 $snps > maps_column.txt
cut -f 1-3 ${file}.map > new_map.txt
mkdir ${outfold}
paste new_map.txt maps_column.txt > ./${outfold}/${file}.map
cp ./${file}.ped ./${outfold}/
rm -rf ./new_map.txt
rm -rf maps_column.txt
cd ./${outfold}

plink --noweb --file ${file} --recode12 --out ${file}

#cross-validation
mkdir ./cross-validation
cd ./cross-validation
for K in 1 2 3 4 5 6 7 8 9 10; do admixture -j16 --cv ../${file}.ped $K | tee log${K}.out; done
cd ..

## C value is the number of the steps for admixture to run over each iteration. It can be removed and left to default. However I consider it more correct to assgin a constant number
## of steps for each iteration. To do this, fist run one iteration and explore the number of steps it needs to reach convergence and use a value for c which ensure convergence for each step.
~/bin/admixture_linux-1.3.0/admixture -j16 -B5000 -c20 ${file}.ped ${K}


cut -d " " -f 1,2 ${file}.ped > names.txt
