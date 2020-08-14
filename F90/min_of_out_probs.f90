!this code finds the minimum of the chi-squared probability in the file of the 192 coordinates and the probabilities
program mi

integer,parameter:: num=192
real(8):: prob(num),minn,dumy,dumy2,dumy3
integer sne_count(num), min_num


minn=100
open(10,file="infile.txt")
!open(10,file="out_prob.11")
do i=1,num
read(10,*) dumy,dumy2,dumy3,prob(i),sne_count(i)
if(sne_count(i).ge.3)then
if(prob(i).lt.minn)then
minn=prob(i)
min_num=sne_count(i)
endif
endif
enddo

open(100,file="outfile.txt")
write(100,*) minn, min_num
 close(100)




end program mi
