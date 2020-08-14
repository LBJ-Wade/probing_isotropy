#! /bin/bash


for i in {101..1000}
do
cp  out_cone.$i inp.txt

gfortran chi_prob_for-script.f90 -I/vol/software/software/tools/nag/fll6a23dhl/nag_interface_blocks -L/vol/software/software/tools/nag/fll6a23dhl/lib -L/vol/software/software/tools/nag/fll6a23dhl/acml -lnag_nag
./a.out
rm -f inp.txt
mv out_prob.txt out_prob.$i
rm -f a.out
echo "$i"
done

