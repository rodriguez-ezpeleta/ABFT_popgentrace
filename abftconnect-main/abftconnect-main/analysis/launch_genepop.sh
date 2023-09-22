#for d in /share/projects/ACEITUNA/csstacks*/075*/; do
#        cd $d
        for f in *_genepop.txt ; do
                rm -rf fichier.in
                echo $f > genepop.cmd
                echo   >> genepop.cmd
                echo 5 >> genepop.cmd
                echo 2 >> genepop.cmd
                echo   >> genepop.cmd
                echo 6 >> genepop.cmd
                echo 2 >> genepop.cmd
                echo   >> genepop.cmd
                echo 9 >> genepop.cmd
                Genepop < genepop.cmd
        done
#done 
