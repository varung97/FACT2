# test_correctness.sh n k num_tests
# Will run num_tests tests checking equivalence of oldFACT and newFACT with n and k
for i in $(seq 1 $3); do
	python3 gen_input_trees.py $1 $2 input_trees.nex
	./newFACT freq input_trees.nex >> newoutput.txt
	./oldFACT freq input_trees.nex >> oldoutput.txt
	echo $(python3 compare_trees.py $(cat newoutput.txt) $(cat oldoutput.txt))
	rm newoutput.txt
	rm oldoutput.txt
done
