for n in 1 2 3 4 5 6 7 8 9 10; do 
	cd m${n}_s*
	for f in *optimized*; do 
	grep -v "^Model" ${f} | sort -n -k 4 | head -n 1 
	done >> ../results_models.txt
	cd ../ 
done
