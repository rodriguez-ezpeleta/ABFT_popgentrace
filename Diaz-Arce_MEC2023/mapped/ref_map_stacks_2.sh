prj=/share/projects/GBYP_ND
in=$prj/mapping
mkdir $prj/ref_map

gstacks -I $in -M $in/popmap.txt --rm-pcr-duplicates -O $prj/ref_map/ -t 8

##### with the new reference genome BKCK2019 [PBFT], change manually for other reference genomes

prj=/share/projects/GBYP_ND
in=$prj/mapping_BKCK_R_1_2
mkdir $prj/BKCK_mapped_catalog

gstacks -I $in -M $in/popmap_540pops.txt --ignore-pe-reads -O $prj/ref_map/ -t 28
grep -A 544 "BEGIN effective_coverages_per_sample" gstacks.log.distribs > tags_per_sample.txt
tail -n 541 tags_per_sample.txt | head -n 540 | awk '$2<25000' > Inds_with_less_than25000TAGS.txt
tail -n 541 tags_per_sample.txt | head -n 540 | awk '$2>25000' | cut -f 1 > popmap_529inds.txt
sed -i 's/$/\tABFT/' popmap_529inds.txt
grep ABFT popmap_529inds.txt > popmap_ABFT.txt
