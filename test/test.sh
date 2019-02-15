echo "Gen trees:"
time python3 gen_input_trees.py $1 $2 input_trees.nex
echo
echo "New FACT"
time ./FACT++ freq input_trees.nex >> newoutput.txt
echo
echo "Old FACT"
time ./oldFACT freq input_trees.nex >> oldoutput.txt
echo
echo "Are results same:" $(python3 compare_trees.py $(cat newoutput.txt) $(cat oldoutput.txt))
rm newoutput.txt
rm oldoutput.txt
