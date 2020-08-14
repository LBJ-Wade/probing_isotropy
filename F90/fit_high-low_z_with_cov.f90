program dis_mod
USE nag_library
implicit none

integer,parameter:: sne_num=580, znum=580, lda=znum+1, n=znum, ldb=znum
integer i,j,k,m, minimum,minimum_m,k_max,k_min,m_max,m_min,error_m,error_k
integer error_m_up,error_m_dn,error_k_up,error_k_dn,IFAIL,data_index(znum)
real,parameter:: omeg_bin=0.01, omeg_min=0.0, omeg_max=1.0
real,parameter:: mu0_bin=0.001, mu0_min=42.0, mu0_max=44.0
integer,parameter:: omeg_num=(omeg_max-omeg_min)/omeg_bin, mu0_num=(mu0_max-mu0_min)/mu0_bin
real(8):: z,z_data(sne_num),dm_data(sne_num),dm_err(sne_num),chi_sq(mu0_num,omeg_num),omega_m,dummy,temp1(znum)
real(8):: chi_sq_min,long,lat, dl(omeg_num,sne_num),mu0,err_om_up,err_om_dn,err_del_up,err_del_dn,vec(znum),cov(sne_num,sne_num)
REAL (KIND=nag_wp) inv(ldb,n), zz(n),cov_z(lda,n)

 Character(10)::names(sne_num)


!********************** reading the covariance matrix ********************
open(12,file="SCPUnion2.1_covmat_sys.txt")
READ(12,*) ((cov(i,j),j=1,sne_num),i=1,sne_num)


!********************** reading the covariance matrix DONE********************


open(10,file="SCPUnion2.1_mu_vs_z_without_header.txt")
do i=1,sne_num
read(10,*) names(i),z_data(i),dm_data(i),dm_err(i),dummy
enddo

open(11,file="omegaM_z_dl_original_order.txt")
do i=1,omeg_num
do j=1,sne_num
read(11,*) omega_m,z,dl(i,j)
enddo
enddo


i=1
do j=1,sne_num
if(z_data(j).lt.2)then
data_index(i)=j
i=i+1
endif
enddo

do i=1,znum
m=data_index(i)
do j=1,znum
k=data_index(j)
cov_z(i,j)=cov(m,k)
enddo
enddo

IFAIL=0
CALL f01abf(cov_z,lda,n,inv,ldb,zz,ifail)

do i=1,znum
do j=i,znum
inv(i,j)=inv(j,i)
enddo
enddo


!************************* chi_sq
do m=1,mu0_num
mu0=(m*mu0_bin)+mu0_min
do i=1,omeg_num
k=1
do j=1,sne_num
if(z_data(j).lt.2)then

vec(k)=mu0+5*log10(dl(i,j))-dm_data(j)
k=k+1
endif
enddo
temp1=matmul(vec,inv)


 chi_sq(m,i)=chi_sq(m,i)+dot_product(vec,temp1)


enddo
enddo


 chi_sq_min=10000.0
do m=1,mu0_num
do k=1,omeg_num
if(chi_sq(m,k).lt.chi_sq_min)then 
 chi_sq_min=chi_sq(m,k)
 minimum=k
minimum_m=m
endif
enddo
enddo

print*, minimum, (minimum-1)*omeg_bin, minimum_m, (minimum_m*mu0_bin)+mu0_min, chi_sq(minimum_m,minimum)

!open(300,file="conf_int_1sigma.txt")
!open(400,file="conf_int_2sigma.txt")

k_max=-1
k_min=200
m_max=-1
m_min=10000
do m=1,mu0_num
do k=1,omeg_num
if(chi_sq(m,k)-chi_sq(minimum_m,minimum).lt.1.1.and.chi_sq(m,k)-chi_sq(minimum_m,minimum).gt.0.9)then
error_m=m
error_k=k

if(k.lt.k_min)then
k_min=k
endif
if(k.gt.k_max)then
k_max=k
endif
if(m.gt.m_max)then
m_max=m
endif
if(m.lt.m_min)then
m_min=m
endif


endif
enddo
enddo


error_m_up=abs(m_max-minimum_m)
error_m_dn=abs(minimum_m-m_min)

error_k_up=abs(k_max-minimum)
error_k_dn=abs(minimum-k_min)

err_del_up=error_m_up*mu0_bin
err_del_dn=error_m_dn*mu0_bin

err_om_up=error_k_up*omeg_bin
err_om_dn=error_k_dn*omeg_bin





!error_m=abs(minimum_m-error_m)
!error_k=abs(minimum-error_k)



!do m=1,mu0_num
!do k=1,omeg_num

!if(chi_sq(m,k)-chi_sq(minimum_m,minimum).le.2.30)then
!write(300,*) 1-(k-1)*omeg_bin, (m*mu0_bin)+mu0_min! , chi_sq(m,k)
!endif
!if(chi_sq(m,k)-chi_sq(minimum_m,minimum).le.6.17)then
!write(400,*) 1-(k-1)*omeg_bin, (m*mu0_bin)+mu0_min !, chi_sq(m,k)
!endif

!enddo
!enddo
! close(300)
! close(400)


open(800,file="fit_results_z_all_3.txt")
write(800,*) "Omega_matter= ",(minimum-1)*omeg_bin, "delta= ",(minimum_m*mu0_bin)+mu0_min, "chi_square= ", chi_sq(minimum_m,minimum)
write(800,*) "omeg_err_up= ",err_om_up,"omeg_err_dn= ",err_om_dn, "delta_err_up= ",err_del_up,"delta_err_dn= ",err_del_dn
 close(800)

!print*, (error_m*mu0_bin), error_k*omeg_bin
end program dis_mod




!gfortran-4.7 fit_with_dl_file_omega_mu0_cov.f90 -I/vol/software/software/tools/nag/fll6a24dfl/nag_interface_blocks -L/vol/software/software/tools/nag/fll6a24dfl/lib -L/vol/software/software/tools/nag/fll6a24dfl/acml -lnag_nag

