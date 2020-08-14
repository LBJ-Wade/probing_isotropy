program cone
implicit none

integer,parameter:: sne_num=580, heal_num=192

integer i,j,k,m, minimum, sne_count,s,minimum_m,k_max,k_min,m_max,m_min,error_m,error_k

real,parameter:: omeg_bin=0.01, omeg_min=0.0, omeg_max=1.0, omeg_num=(omeg_max-omeg_min)/omeg_bin, sigma_sample=0.0, cone_angle=90.0,mu0_bin=0.001, mu0_min=42.0, mu0_max=44.0, mu0_num=(mu0_max-mu0_min)/mu0_bin

real(8):: z,z_data(sne_num),dm_data(sne_num),dm_err(sne_num),dl(omeg_num,sne_num),chi_sq(mu0_num,omeg_num),omega_m, chi_sq_min, long(sne_num), lat(sne_num), cone_vec_long, cone_vec_lat,cone_vec_x,cone_vec_y,cone_vec_z,radius_of_separation,sn_vec_x,sn_vec_y,sn_vec_z,sn_separation, error,heal_lat(heal_num), heal_long(heal_num),mu0
 

radius_of_separation=sind(cone_angle)/sind((180-cone_angle)/2)   !checked

open(10,file="input.txt")
do i=1,sne_num
read(10,*) z_data(i),dm_data(i),dm_err(i), long(i), lat(i)
enddo
lat=90-lat

open(11,file="Healpix_Coordinates_nside_4.txt")
do i=1,heal_num
read(11,*) heal_lat(i), heal_long(i)
enddo

open(12,file="omegaM_z_dl.txt")
do i=1,omeg_num
do j=1,sne_num
read(12,*) omega_m,z,dl(i,j)
enddo
enddo


open(400,file="out_cone.txt")
!open(400,file="long_lat_sn-num_omega-&er_mu0-&er_chi-sq_60_z-all.txt")

do s=1,heal_num

!print*, heal_lat(s), heal_long(s)

cone_vec_long=heal_long(s)
cone_vec_lat=heal_lat(s)
! ********************** calculating r_p **********************************
cone_vec_x=sind(cone_vec_lat)*cosd(cone_vec_long)
cone_vec_y=sind(cone_vec_lat)*sind(cone_vec_long)
cone_vec_z=cosd(cone_vec_lat)
!************************ end of calculating r_p ********************************



!********************************* chi square fitting ************************
 sne_count=0
 chi_sq=0

do j=1,sne_num
if(z_data(j).gt.0.2)then
!****************** calculating r_i ******************************************
sn_vec_x=sind(lat(j))*cosd(long(j))
sn_vec_y=sind(lat(j))*sind(long(j))
sn_vec_z=cosd(lat(j))
!****************** end of calculating r_i ******************************************

!****************** calculating the separation of the sn from the r_p ********************
sn_separation=sqrt((cone_vec_x-sn_vec_x)**2+(cone_vec_y-sn_vec_y)**2+(cone_vec_z-sn_vec_z)**2)
!****************** end calculating the separation of the sn from the r_p ********************


if(sn_separation.le.radius_of_separation)then      !checking the cone condition
sne_count=sne_count+1


do m=1,mu0_num
mu0=(m*mu0_bin)+mu0_min

do i=1,omeg_num


 chi_sq(m,i)=chi_sq(m,i)+(((mu0+5*log10(dl(i,j))-dm_data(j))**2)/(dm_err(j)**2+sigma_sample**2))

!**************************** end of chi square fitting *****************

enddo   !do i=1,omeg_num
enddo   !do m=1,mu0_num

endif    !if(sn_separation.le.radius_of_separation)then
endif   !if(z_data(j).gt.0.2)then

enddo   !do j=1,sne_num



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

!write(300,*) (k-1)*omeg_bin, (m*mu0_bin)+mu0_min, chi_sq(m,k)
!if(k.lt.k_min)then
!k_min=k
!endif
!if(k.gt.k_max)then
!k_max=k
!endif
!if(m.gt.m_max)then
!m_max=m
!endif
!if(m.lt.m_min)then
!m_min=m
!endif
endif
enddo
enddo


error_m=abs(minimum_m-error_m)
error_k=abs(minimum-error_k)

!if(sne_count.gt.5)then
!write(400,5) cone_vec_long, 90-cone_vec_lat, sne_count, (minimum-1)*omeg_bin,error_k*omeg_bin, (minimum_m*mu0_bin)+mu0_min,(error_m*mu0_bin), chi_sq(minimum_m,minimum),chi_sq(1160,29)
!5 format (F,F,I,F,F,F,F,F,F)
write(400,5) cone_vec_long, 90-cone_vec_lat, sne_count,(minimum-1)*omeg_bin,error_k*omeg_bin, chi_sq(minimum_m,minimum),chi_sq(1160,29)
5 format (F,F,I,F,F,F,F)
!endif


enddo   !do s=1,heal_num

 close(400)



end program cone
