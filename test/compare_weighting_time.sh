# compare_weighting_time.sh n k
# Will run oldFACT and newFACT algorithm on same input and echo timing statistics for only the weighting step
echo "Gen trees:"
time python3 gen_input_trees.py $1 $2 input_weights.nex
echo
echo "New FACT"
time ./FACT++ freq input_weights.nex onlyw >> newoutput.txt
echo
echo "Old FACT"
time ./oldFACT freq input_weights.nex onlyw >> oldoutput.txt
rm newoutput.txt
rm oldoutput.txt
