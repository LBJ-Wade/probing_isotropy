program per

integer,parameter:: num=1000
real(8) prob(num)
integer boz

open(10,file="map_prob.txt")
do i=1,num
read(10,*) prob(i)
if(prob(i).le.0.1170300848253082)then
boz=boz+1
endif
enddo


print*,boz



end program per
