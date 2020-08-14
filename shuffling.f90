program shuffle
implicit none

integer,parameter:: sne_num=350
integer i,j,k,done_sim(sne_num),done_data(sne_num),n_run
real(8):: x,y, z_data(sne_num),dm_data(sne_num),dm_err(sne_num),long(sne_num),lat(sne_num),long_new(sne_num),lat_new(sne_num)
 character(256) ::numstr,numstring
open(10,file="cut_on_redshift.txt")
do i=1,sne_num
read(10,*) z_data(i),dm_data(i),dm_err(i),long(i),lat(i)
enddo

call Random_Seed

do n_run=1,1000
 done_sim=0
 done_data=0
 do i=1,sne_num
 do while(done_data(i).eq.0)
 call Random_Number(x)
    k=int((x*sne_num)+1)
    if(done_sim(k).eq.0)then
     long_new(k)=long(i)
     lat_new(k)=lat(i)
     done_sim(k)=1
     done_data(i)=1
    endif
 enddo
 enddo

 if(n_run.le.9)then
 write(numstr,'(i1)') n_run
 numstring='rand.'//numstr
 open(70,file=numstring)
 endif

 if(n_run.ge.10.and.n_run.le.99)then
 write(numstr,'(i2)') n_run
 numstring='rand.'//numstr
 open(70,file=numstring)
 endif

 if(n_run.ge.100.and.n_run.lt.1000)then
 write(numstr,'(i3)') n_run
 numstring='rand.'//numstr
 open(70,file=numstring)
 endif

 if(n_run.eq.1000)then
 write(numstr,'(i4)') n_run
 numstring='rand.'//numstr
 open(70,file=numstring)
 endif

 do i=1,sne_num
 write(70,5) z_data(i),dm_data(i),dm_err(i),long_new(i),lat_new(i)
 5 format (F,F,F,F,F)
 enddo

print*, n_run
enddo


!open(100,file="rand_1.txt")
!do i=1,sne_num
!write(100,5) z_data(i),dm_data(i),dm_err(i),long_new(i),lat_new(i)
!5 format (F,F,F,F,F)
!enddo






end program shuffle
