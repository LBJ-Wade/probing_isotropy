!This code compute the chi squred probability for the value of the chi_sq_best 
program test
USE nag_library
integer, parameter:: num=192
integer(8) ifail
real(kind=nag_wp)  x, df,prob(num),prob_best(num)    !,dummy1, dummy2, dummy3, dummy4, dummy5,dummy6, dummy7, dummy8
 character(1) tail


integer i,j,k ,sne_count(num)
real(8)::long(num), lat(num), omega(num), error_omega(num),chi_sq(num),chi_sq_best(num) !,error_mu(num),mu(num)


open(10,file="inp.txt")


do i=1,num
read(10,*) long(i), lat(i), sne_count(i),omega(i), error_omega(i),chi_sq(i),chi_sq_best(i)
enddo

ifail=0
tail='u'

do i=1,num
df=sne_count(i)-2
x=chi_sq(i)
prob(i)= g01ecf (tail,x,df,ifail)
x=chi_sq_best(i)
prob_best(i)= g01ecf (tail,x,df,ifail)
enddo

open(400,file="out_prob.txt")
do i=1,num
!write(400,5) long(i), lat(i), sne_count(i), omega(i), error_omega(i),mu(i), error_mu(i), chi_sq(i),chi_sq_best(i),prob(i),prob_best(i),prob_best(i)-prob(i)
write(400,*) long(i), lat(i), prob(i), prob_best(i)
!5 format (F,F,I,F,F,F,F,F,F,F,F,F)
enddo
 close(400)



end program test




!gfortran test.f90 -I/vol/software/software/tools/nag/fll6a23dhl/nag_interface_blocks -L/vol/software/software/tools/nag/fll6a23dhl/lib -L/vol/software/software/tools/nag/fll6a23dhl/acml -lnag_nag

