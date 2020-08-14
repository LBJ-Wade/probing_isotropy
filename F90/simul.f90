program sim
implicit none

integer, parameter:: sne_num=350
integer i,j,k,temp1,temp2,done_sim(sne_num),done_data(sne_num),n_run
real,parameter:: z_max=1.5, omeg_bin=0.01, omeg_min=0.0, omeg_max=1.0, omeg_num=1+(omeg_max-omeg_min)/omeg_bin
real(8):: x,y,z,q,z_data(sne_num),res(sne_num),dm_err(sne_num),long(sne_num),lat(sne_num),dm_best(sne_num),z_new(sne_num),dm_new(sne_num)
 character(256) ::numstr,numstring
call Random_Seed


open(10,file="z_dmbest_res_dm-err_long_lat.txt")
do i=1,sne_num
read(10,*) z_data(i),dm_best(i),res(i),dm_err(i),long(i),lat(i)
enddo


do n_run=1,1000 
 done_sim=0
 done_data=0
 do i=1,sne_num
 do while(done_data(i).eq.0)
 call Random_Number(x)
    k=int((x*sne_num)+1)
    if(done_sim(k).eq.0)then
     z_new(k)=z_data(i)
     dm_new(k)=dm_best(i)
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
 call Random_Number(y)
 temp1=int(351*y)
 write(70,5) z_new(i),dm_new(i)+res(temp1),dm_err(i),long(i),lat(i)
 5 format (F,F,F,F,F)
 enddo


print*, n_run
enddo


end program sim
