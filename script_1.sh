#! /bin/bash

for i in {101..1000}
do
mv  rand.$i input.txt

ifort cone_analysis_healp_mu_dl-file_for-script.f90
./a.out
rm -f input.txt
mv out_cone.txt out_cone.$i
rm -f a.out
echo "$i"
done

