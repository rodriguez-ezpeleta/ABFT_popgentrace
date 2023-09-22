gstacks -I /share/projects/GBYP_ND/mapping_mtDNA -M /share/projects/GBYP_ND/popmap_ONE_527ABFT.txt \ 
--rm-pcr-duplicates -O /share/projects/GBYP_ND/gstacks_mtDNA -t 16

mkdir /share/projects/GBYP_ND/mapped_catalog_mtDNA/populations_r075
cd /share/projects/GBYP_ND/mapped_catalog_mtDNA/populations_r075
populations -P ./ -M ../popmap_ONE_527ABFT.txt -r 0.75 --plink --vcf -t 16
plink --noweb --file populations.plink --extract diagnostic_post.txt --out plink_diagnostic --recode

