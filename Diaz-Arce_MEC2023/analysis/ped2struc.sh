##convert to str##
sed "s/map_file/${file}.map/g" /share/projects/GBYP_ND/ped2structure.spid > ped2structure.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile ${file}.ped -outputfile ${file}.str -spid ped2structure.spid
