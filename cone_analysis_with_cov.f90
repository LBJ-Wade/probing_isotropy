!corrected for error estimation
program cone
USE nag_library
implicit none

integer,parameter:: sne_num=580, heal_num=192

integer i,j,k,m, p,minimum, sne_count,s,minimum_m,k_max,k_min,m_max,m_min,error_m,error_k
integer error_m_up,error_m_dn,error_k_up,error_k_dn,cone_num(heal_num),ldb,n,lda,IFAIL

real,parameter:: omeg_bin=0.01, omeg_min=0.0, omeg_max=1.0, sigma_sample=0.0
real,parameter::cone_angle=30.0,mu0_bin=0.001, mu0_min=42.0, mu0_max=44.0
integer,parameter:: omeg_num=(omeg_max-omeg_min)/omeg_bin, mu0_num=(mu0_max-mu0_min)/mu0_bin

real(8):: z,z_data(sne_num),dm_data(sne_num),dm_err(sne_num),dl(omeg_num,sne_num),chi_sq(mu0_num,omeg_num),omega_m, chi_sq_min
real(8):: long(sne_num), lat(sne_num), cone_vec_long, cone_vec_lat,cone_vec_x,cone_vec_y,cone_vec_z,radius_of_separation
real(8)::sn_vec_x,sn_vec_y,sn_vec_z,sn_separation, error,heal_lat(heal_num), heal_long(heal_num),mu0,torad,delta,error_del
real(8):: error_om,err_om_up,err_om_dn,err_del_up,err_del_dn,chimin,chibest,cov(sne_num,sne_num),dummy

real*8,allocatable:: inv(:,:), zz(:),cov_cone(:,:),vec(:),temp1(:)
integer,allocatable:: data_index(:)
 character:: names
torad=3.14159265359/180.0

radius_of_separation=sin(torad*cone_angle)/sin(torad*(180-cone_angle)/2)   !checked


!********************** reading the covariance matrix ********************
!open(12,file="SCPUnion2.1_covmat_sys.txt")
!READ(12,*) ((cov(i,j),j=1,sne_num),i=1,sne_num)


!********************** reading the covariance matrix DONE********************



open(10,file="name_z_mu_err_prob_l_b.txt")
do i=1,sne_num
read(10,*) names,z_data(i),dm_data(i),dm_err(i), dummy,long(i), lat(i)
enddo
lat=90-lat

!*************** building the diagonal covariance matrix ********
cov=0
do i=1,sne_num
cov(i,i)=dm_err(i)**2
enddo

!*************** end of building the diagonal covariance matrix ********


open(11,file="Healpix_Coordinates_nside_4.txt")
do i=1,heal_num
read(11,*) heal_lat(i), heal_long(i)
enddo

open(13,file="omegaM_z_dl_original_order.txt")
do i=1,omeg_num
do j=1,sne_num
read(13,*) omega_m,z,dl(i,j)
enddo
enddo

open(14,file="cone_num_30_zgt02.txt")
do i=1,heal_num
read(14,*) cone_num(i)
enddo


open(400,file="long_lat_sn-num_omega-&er_delta-&er_chi-sq_chi-sq-best_30_z-gt-02_cov.txt")
!open(400,file="long_lat_sn-num_omega-&er_delta-&er_chi-sq_chi-sq-best_60_z-all.txt")


do s=1,heal_num
!do s=122,122   !s=122,122  s=182,182 

lda=cone_num(s)+1
n=cone_num(s)
ldb=cone_num(s)

allocate(data_index(n))
allocate(vec(n))
allocate(temp1(n))
allocate(cov_cone(lda,n))
allocate(inv(ldb,n))
allocate(zz(n))

print*, heal_lat(s), heal_long(s)

cone_vec_long=heal_long(s)
cone_vec_lat=heal_lat(s)
! ********************** calculating r_p **********************************
cone_vec_x=sin(torad*cone_vec_lat)*cos(torad*cone_vec_long)
cone_vec_y=sin(torad*cone_vec_lat)*sin(torad*cone_vec_long)
cone_vec_z=cos(torad*cone_vec_lat)
!************************ end of calculating r_p ********************************



!********************************* chi square fitting ************************
 sne_count=0
 chi_sq=0
p=1
do j=1,sne_num
if(z_data(j).gt.0.2)then
!****************** calculating r_i ******************************************
sn_vec_x=sin(torad*lat(j))*cos(torad*long(j))
sn_vec_y=sin(torad*lat(j))*sin(torad*long(j))
sn_vec_z=cos(torad*lat(j))
!****************** end of calculating r_i ******************************************

!****************** calculating the separation of the sn from the r_p ********************
sn_separation=sqrt((cone_vec_x-sn_vec_x)**2+(cone_vec_y-sn_vec_y)**2+(cone_vec_z-sn_vec_z)**2)
!****************** end calculating the separation of the sn from the r_p ********************


if(sn_separation.le.radius_of_separation)then      !checking the cone condition
data_index(p)=j   !index of the SNe in the cone
p=p+1
sne_count=sne_count+1

endif    !if(sn_separation.le.radius_of_separation)then
endif   !if(z_data(j).gt.0.2)then
enddo   !do j=1,sne_num

if(sne_count.ge.5)then


!**************************************** BUILDING THE COVARIANCE MATRIX ************************************
do i=1,n
do j=1,n
m=data_index(i)
k=data_index(j)
cov_cone(i,j)=cov(m,k)
enddo
enddo
!**************************************** BUILDING THE COVARIANCE MATRIX   DONE   ************************************

!******************************************* INVERTING THE COV MATRIX *************************************
IFAIL=0


CALL f01abf(cov_cone,lda,n,inv,ldb,zz,ifail)

do i=1,cone_num(s)
do j=i,cone_num(s)
inv(i,j)=inv(j,i)
enddo
enddo

!******************************************* INVERTING THE COV MATRIX  DONE *************************************


!**************************** chi square fitting *****************
do m=1,mu0_num
mu0=(m*mu0_bin)+mu0_min
do i=1,omeg_num

do j=1,n
p=data_index(j)
vec(j)=mu0+5*log10(dl(i,p))-dm_data(p)
enddo
temp1=matmul(vec,inv)

 chi_sq(m,i)=chi_sq(m,i)+dot_product(vec,temp1)

enddo   !do i=1,omeg_num
enddo   !do m=1,mu0_num


!**************************** end of chi square fitting *****************




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


!error_m=abs(max(m_max-minimum_m,minimum_m-m_min))
!error_k=abs(max(k_max-minimum,minimum-k_min))

error_m_up=abs(m_max-minimum_m)
error_m_dn=abs(minimum_m-m_min)

error_k_up=abs(k_max-minimum)
error_k_dn=abs(minimum-k_min)

err_del_up=error_m_up*mu0_bin
err_del_dn=error_m_dn*mu0_bin

err_om_up=error_k_up*omeg_bin
err_om_dn=error_k_dn*omeg_bin

omega_m=(minimum-1)*omeg_bin
delta=(minimum_m*mu0_bin)+mu0_min


 chimin=chi_sq(minimum_m,minimum)
 !chibest=chi_sq(1171,31)
 chibest=chi_sq(1159,29)
write(400,*) cone_vec_long, 90-cone_vec_lat, sne_count, omega_m,err_om_up,err_om_dn, delta,err_del_up,err_del_dn, chimin,chibest
!5 format(F,F,I,F,F,F,F,F,F,F,F)

!print*, cone_vec_long, 90-cone_vec_lat, sne_count, omega_m,err_om_up,err_om_dn, delta,err_del_up,err_del_dn, chimin,chibest

else

write(400,*) cone_vec_long, 90-cone_vec_lat, sne_count, 0,0,0, 0,0,0, 0,0

endif

deallocate(data_index)
deallocate(vec)
deallocate(temp1)
deallocate(cov_cone)
deallocate(inv)
deallocate(zz)
enddo   !do s=1,heal_num

 close(400)



end program cone

!gfortran-4.7 cone_analysis_with_cov.f90 -I/vol/software/software/tools/nag/fll6a24dfl/nag_interface_blocks -L/vol/software/software/tools/nag/fll6a24dfl/lib -L/vol/software/software/tools/nag/fll6a24dfl/acml -lnag_nag
