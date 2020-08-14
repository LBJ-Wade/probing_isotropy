#! /bin/bash


for i in {1..1000}
do
rm -f a.out
rm -f outfile.txt
cp  out_prob.$i infile.txt
ifort min_of_out_probs.f90
./a.out
awk '{print $1}' outfile.txt>>map_prob.txt
rm -f infile.txt 
done
