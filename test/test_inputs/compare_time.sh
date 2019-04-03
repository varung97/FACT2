# compare_time.sh n k
# Will run oldFACT and newFACT algorithm on same input and echo timing statistics
# echo "Gen trees:"

scen=1
num_trials=5

# Fixed k, varying n
# ks=(100)
# ns=(500 1000 1500 2000 3000 4000 5000 7500 10000)

# Varying k, fixed n
ks=(4000)
ns=(1000)

echo $num_trials
for k in "${ks[@]}"; do
	for n in "${ns[@]}"; do
		echo $n, $k
		for i in $(seq 1 $num_trials); do
			if [ ! -f input_trees_${scen}_${n}_${k}_${i}.nex ]; then
				python3 gen_input_trees.py $n $k $scen input_trees_${scen}_${n}_${k}_${i}.nex
			fi
			./newFACT freq input_trees_${scen}_${n}_${k}_${i}.nex
			./oldFACT freq input_trees_${scen}_${n}_${k}_${i}.nex 1
			./oldFACT freq input_trees_${scen}_${n}_${k}_${i}.nex 2
		done
	done
done

# python3 gen_input_trees.py $1 $2 input_trees.nex
# echo
# echo "New FACT"
# ./newFACT freq input_trees.nex
# echo
# echo "Old FACT"
# ./oldFACT freq input_trees.nex
# echo
# echo "Are results same:" $(python3 compare_trees.py $(cat newoutput.txt) $(cat oldoutput.txt))
# rm newoutput.txt
# rm oldoutput.txt
