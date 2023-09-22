sed "s/map_file/${file}.map/g" /share/projects/GBYP_ND/ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile ${file}.ped -outputfile ${file}_genepop.txt -spid ped2genepop.spid

