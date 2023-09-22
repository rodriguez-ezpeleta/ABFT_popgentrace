#!/bin/bash
prj=/share/projects/GBYP_ND
sample=(`cut -f 1 $prj/samples_csstacks.txt`)
temp=(`wc -l $prj/samples_csstacks.txt`)
c=(`expr $temp - 1`)
m=3
M=2
n=6
folder=stacks_2_m${m}_M${M}_n${n}

#
# Build loci de novo in each sample for the single-end reads only. If paired-end reads are available, 
# they will be integrated in a later stage (tsv2bam stage).
# This loop will run ustacks on each sample, e.g.
#   ustacks -f ./samples/sample_01.1.fq.gz -o ./stacks -i 1 --name sample_01 -M 4 --gapped -p 8
#

mkdir $prj/${folder}
echo "" > clean_nc_reads_lines.txt
echo "" > clean_nc_reads_names.txt

id=1
for i in `seq 0 $c`; do
	ustacks -f $prj/clean_nc/${sample[i]}.fq_1.gz -o $prj/${folder} -i $id -t gzfastq --name ${sample[i]} -m ${m} -M ${M} --max_locus_stacks 2 -d -p 8 > $prj/${folder}/${sample[i]}.log
	let "id+=1"
	zcat $prj/clean_nc/${sample[i]}.fq_1.gz | wc -l >> $prj/${folder}/clean_nc_reads_lines.txt
	echo "${sample[i]}" >> clean_nc_reads_names.txt	
done

echo "Names       UsedReads       Tags    Alleles SNPs" > $prj/${folder}/ustacks_stats.txt
echo "" > Names.txt
echo "" > Used_Reads.txt
echo "" > Tags.txt
echo "" > Alleles.tt
echo "" > SNPs.txt

for i in `seq 0 $c`; do
	echo "${sample[i]}" >> Names.txt
	zcat $prj/${folder}/${sample[i]}.tags.tsv.gz | grep -c 'primary\|secondary' >> $prj/${folder}/Used_Reads.txt
	zcat ${sample[i]}.tags.tsv.gz | grep -c consensus >> $prj/${folder}/Tags.txt
	zcat ${sample[i]}.alleles.tsv.gz | wc -l >> $prj/${folder}/Alleles.txt
	zcat ${sample[i]}.snps.tsv.gz | grep -c 'E' >> $prj/${folder}/SNPs.txt
done

# 
# Build the catalog of loci available in the metapopulation from the samples contained
# in the population map. To build the catalog from a subset of individuals, supply
# a separate population map only containing those samples.
#

cstacks -n ${n} -P $prj/${folder}/ -M $prj/popmap_527ABFT.txt -p 8 

#
# Run sstacks. Match all samples supplied in the population map against the catalog.
#
sstacks --gapped -P $src/stacks/ -M $src/popmaps/popmap -p 8

# Run tsv2bam to transpose the data so it is stored by locus, instead of by sample. We will include
# paired-end reads using tsv2bam. tsv2bam expects the paired read files to be in the samples
# directory and they should be named consistently with the single-end reads,
# e.g. sample_01.1.fq.gz and sample_01.2.fq.gz, which is how process_radtags will output them.
#
tsv2bam -P $prj/${folder} -R $prj/clean_nc/ -M $prj/mapping/popmap.txt -t 8
#
# Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided),
# align reads per sample, call variant sites in the population, genotypes in each individual.
#
gstacks -P $prj/${folder} -M $prj/mapping/popmap_thynnus.txt -t 8
#
# Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics
# export several output files.
#

# A PARTIR DEL NUMERO DE TAGS DE USTACKS, HE SELECCIONADO INDIVIDUOS CON > 25,000 TAGS ###
populations -P $prj/${folder} -M $prj/${folder}/popmap_485inds_morethan25000tags_ONEPOP.txt -r 0.75 --plink -t 16
mkdir $prj/${folder}/population_r075
mv $prj/${folder}/populations.* $prj/${folder}/populations_r075/
grep  -v "^#" $prj/${folder}/populations_r075/populations.sumstats.tsv | cut -f 1 | sort | uniq > $prj/${folder}/populations_r075/tagsIn075.tsv
populations -P $prj/${folder} -M $prj/${folder}/popmap_485inds_morethan25000tags_13POP.txt -t 8 --plink -W $prj/${folder}/populations_r075/tagsIn075.tsv
rm -rf $prj/${folder}/populations_r075/populations.*
mv $prj/${folder}/populations.* $prj/${folder}/populations_r075/
cd $prj/${folder}/populations_r075/
plink --noweb --file populations.plink --missing
sed -i 's/^ \+//' plink.imiss



###### MAPPED CATALOG ######

gstacks -I $prj/mapping_R_1_2 -M $prj/mapping/popmap.txt -O $prj/mapped_stacks_v2_3e -t 8
populations -P $prj/mapped_stacks_v2_3e -M $prj/mapped_stacks_v2_3e/popmap_528_ONE.txt -r 0.75 --plink -t 16
mkdir $prj/mapped_stacks_v2_3e/population_r075
mv $prj/mapped_stacks_v2_3e/populations.* $prj/mapped_stacks_v2_3e/population_r075
grep  -v "^#" $prj/mapped_stacks_v2_3e/population_r075populations.sumstats.tsv | cut -f 1 | sort | uniq > $prj/mapped_stacks_v2_3e/population_r075/tagsIn075.tsv
populations -P $prj/mapped_stacks_v2_3e/ -M $prj/mapped_stacks_v2_3e/popmap_528.txt --plink -t 16 -W $prj/mapped_stacks_v2_3e/population_r075/tagsIn075.tsv
