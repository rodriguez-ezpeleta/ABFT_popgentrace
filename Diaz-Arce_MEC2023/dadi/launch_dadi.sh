m=$1
t=$2

for i in `seq 1 $t`; do

	echo -n "Model is: "
	echo $m
	python dadi_Run_3D_Set_Para.py $m $i > ${m}_${i}.log

done

