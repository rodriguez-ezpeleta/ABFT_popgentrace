ped=gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.ped
map=gbyp_mind25_geno10_maf05_hwe05_oneSNPperTAG.map

mkdir permutated_Fst_oneSNPperTAG
cd permutated_Fst_oneSNPperTAG
cp ../$ped ./
cp ../$map ./


### MED - GOM ###
sed 's/^.MED./MED/' $ped | sed 's/^GOM./GOM/' | grep "MED\|GOM" | cut -d " " -f 2- > incomplete_ped_MED_GOM.ped
echo "Fst_MED_GOM" > Fst_MED_GOM.txt
for f in `seq 1 10000`; do
sed 's/^.MED./MED/' $ped | sed 's/^GOM./GOM/' | grep "MED\|GOM" | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped_MED_GOM.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_MED_GOM.txt
rm -rf permut_${f}_genepop.txt.*
done


### MED - SS ###
sed 's/^.MED./MED/' $ped | sed 's/^SS./SS/' | grep "MED\|SS" | cut -d " " -f 2- > incomplete_ped_MED_SS.ped
echo "Fst_MED_SS" > Fst_MED_SS.txt
for f in `seq 1 10000`; do
sed 's/^.MED./MED/' $ped | sed 's/^SS./SS/' | grep "MED\|SS" | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped_MED_SS.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_MED_SS.txt
rm -rf permut_${f}_genepop.txt.*
done

### GOM - SS ###
sed 's/^GOM./GOM/' $ped | sed 's/^SS./SS/' | grep "GOM\|SS" | cut -d " " -f 2- > incomplete_ped_GOM_SS.ped
echo "Fst_GOM_SS" > Fst_GOM_SS.txt
for f in `seq 1 10000`; do
sed 's/^GOM./GOM/' $ped | sed 's/^SS./SS/' | grep "GOM\|SS" | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped_GOM_SS.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_GOM_SS.txt
rm -rf permut_${f}_genepop.txt.*
done


echo "Fst_MED_GOMA" > Fst_MED_GOMA_new.txt
for f in `seq 1 10000`; do
sed 's/^.MED./MED/' $ped | grep "MED\|GOMA" | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_MED_GOMA_new.txt
rm -rf permut_${f}_genepop.txt.*
done

echo "Fst_GOML_GOMA" > Fst_GOML_GOMA_new.txt
for f in `seq 1 10000`; do
grep "GOML\|GOMA" $ped | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_GOML_GOMA_new.txt
rm -rf permut_${f}_genepop.txt.*
done

echo "Fst_MED_SLOPEL" > Fst_MED_SLOPEL_new.txt
for f in `seq 1 10000`; do
sed 's/^.MED./MED/' $ped | grep "MED\|SLOPEL" | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_MED_SLOPEL_new.txt
rm -rf permut_${f}_genepop.txt.*
done

echo "Fst_MED_SLOPEY" > Fst_MED_SLOPEY_new.txt
for f in `seq 1 10000`; do
sed 's/^.MED./MED/' $ped | grep "MED\|SLOPEY" | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_MED_SLOPEY_new.txt
rm -rf permut_${f}_genepop.txt.*
done

echo "Fst_GOML_SLOPEL" > Fst_GOML_SLOPEL_new.txt
for f in `seq 1 10000`; do
grep "GOML\|SLOPEL" $ped | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_GOML_SLOPEL_new.txt
rm -rf permut_${f}_genepop.txt.*
done

echo "Fst_GOML_SLOPEY" > Fst_GOML_SLOPEY_new.txt
for f in `seq 1 10000`; do
grep "GOML\|SLOPEY" $ped | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_GOML_SLOPEY_new.txt
rm -rf permut_${f}_genepop.txt.*
done

echo "Fst_GOMA_SLOPEL" > Fst_GOMA_SLOPEL_new.txt
for f in `seq 1 10000`; do
grep "GOMA\|SLOPEL" $ped | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_GOMA_SLOPEL_new.txt
rm -rf permut_${f}_genepop.txt.*
done

echo "Fst_GOMA_SLOPEY" > Fst_GOMA_SLOPEY_new.txt
for f in `seq 1 10000`; do
grep "GOMA\|SLOPEY" $ped | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid
        rm -rf fichier.in
        echo permut_${f}_genepop.txt > genepop.cmd
        echo   >> genepop.cmd
        echo 5 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 6 >> genepop.cmd
        echo 2 >> genepop.cmd
        echo   >> genepop.cmd
        echo 9 >> genepop.cmd
        Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_GOMA_SLOPEY_new.txt
rm -rf permut_${f}_genepop.txt.*
done

echo "Fst_SLOPEL_SLOPEY" > Fst_SLOPEL_SLOPEY_new.txt
for f in `seq 1 10000`; do
grep "SLOPEL\|SLOPEY" $ped | cut -d " " -f 1 | shuf > extracolumn${f}.txt
paste extracolumn${f}.txt incomplete_ped.ped | sort > permut_${f}.ped
cp $map permut_${f}.map
sed "s/map_file/permut_${f}.map/g" ../../../ped2genepop.spid > ped2genepop.spid
java -Xmx1024m -Xms512m -jar /home/naiara/downloaded_programs/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile permut_${f}.ped -outputfile permut_${f}_genepop.txt -spid ped2genepop.spid

	rm -rf fichier.in
	echo permut_${f}_genepop.txt > genepop.cmd
	echo   >> genepop.cmd
	echo 5 >> genepop.cmd
	echo 2 >> genepop.cmd
	echo   >> genepop.cmd
	echo 6 >> genepop.cmd
	echo 2 >> genepop.cmd
	echo   >> genepop.cmd
	echo 9 >> genepop.cmd
	Genepop < genepop.cmd

rm -rf permut_${f}_genepop.txt
rm -rf permut_${f}.ped
rm -rf permut_${f}.map
rm -rf extracolumn${f}.txt
head -n 4 permut_${f}_genepop.txt.MIG | tail -n 1 >> Fst_SLOPEL_SLOPEY_new.txt
rm -rf permut_${f}_genepop.txt.*
done
