program modif                                         

USE netcdf                                            

use gsw_mod_kinds
use gsw_mod_netcdf
use gsw_mod_toolbox

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

IMPLICIT NONE                                         

INTEGER :: fidT, dimID_depth, dimID_x, dimID_y, dimID_nbounds, mz, mx, my, mnbounds, t_sd_ID, t_an_ID, s_sd_ID, s_an_ID, fidM,  &
&          z_bnds_ID, depth_ID, lon_bnds_ID, lon_ID, lat_bnds_ID, lat_ID, fidS, k30, iWOA, jWOA, kWOA, kisf, kr, ks, k01,       &
&          mean_s_an_1_ID, mean_t_an_1_ID, mean_s_an_2_ID, mean_t_an_2_ID, mean_s_an_3_ID, mean_t_an_3_ID, fidISF, status,      &
&          dimID_Nisf, Nisf, front_max_lat_ID, z_ID, front_min_lat_ID, front_max_lon_ID, front_min_lon_ID, front_ice_dep_avg_ID,&
&          front_ice_dep_min_ID, front_bot_dep_avg_ID, front_bot_dep_max_ID, melt_isf_ID, ii, jj, kk, kdepth, mdepth,           &
&          temperature_ID, salinity_ID, temperature_Con_Shlf_ID, temperature_Offshore_ID, salinity_Con_Shlf_ID,                 &
&          salinity_Offshore_ID

CHARACTER(LEN=150) :: file_isf, file_out, file_in_T, file_in_S

INTEGER, ALLOCATABLE, DIMENSION(:) :: NN, k1, k2

REAL*8,ALLOCATABLE,DIMENSION(:) :: front_max_lat, front_min_lat, front_max_lon, front_min_lon, front_ice_dep_avg, front_ice_dep_min, &
&                                  front_bot_dep_avg, front_bot_dep_max, z, depth, melt_isf
 
REAL*8,ALLOCATABLE,DIMENSION(:,:) :: z_bnds, lon, lat, tmplon

REAL*8,ALLOCATABLE,DIMENSION(:,:) :: mean_t_an_1, mean_s_an_1, &
&                                    mean_t_an_2, mean_s_an_2, &
&                                    mean_t_an_3, mean_s_an_3
 
REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: t_an, s_an, theta_an, temperature, salinity, temperature_Con_Shlf, &
&                                      temperature_Offshore, salinity_Con_Shlf, salinity_Offshore, t_sh,  &
&                                      t_of, s_sh, s_of, theta_sh, theta_of

REAL*8 :: dif, mindif, front_inf, front_sup, &
&         lon_box_1, lon_box_2, lon_box_3, lat_box_1, lat_box_2, lat_box_3

REAL*8 :: Sabs

LOGICAL :: ll_z_positive
 
file_in_T  = '/store/njourd/DATA/WOA2018_EN4_MEOP/obs_temperature_1995-2017_8km_x_60m_masked.nc'
file_in_S  = '/store/njourd/DATA/WOA2018_EN4_MEOP/obs_salinity_1995-2017_8km_x_60m_masked.nc'
file_isf   = 'RTopo-2.0.1_30sec_ice_shelf_mask_ANT3.nc'
file_out   = 'TS_per_ice_shelf_ANT3_WOA2018_EN4_MEOP.nc'

! Boxes :
! NB: 5.0 x 1.75 is the effective resolution at 70°S for a model of 1° resolution in longitude (assuming 5 delta X and a Mercator grid)
lon_box_1 =  5.0 ; lat_box_1 = 1.75 ! Continental shelf only
lon_box_2 = 10.0 ; lat_box_2 = 3.50 ! Continental shelf only
lon_box_3 = 20.0 ; lat_box_3 = 7.00 ! Offshore only

call gsw_saar_init (.true.)
 
!---------------------------------------
! Read temperature :

write(*,*) 'Reading ', TRIM(file_in_T)

status = NF90_OPEN(TRIM(file_in_T),0,fidT); call erreur(status,.TRUE.,"read WOA temperature")

status = NF90_INQ_DIMID(fidT,"z",dimID_depth); call erreur(status,.TRUE.,"inq_dimID_depth")
status = NF90_INQ_DIMID(fidT,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidT,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidT,"nbounds",dimID_nbounds); call erreur(status,.TRUE.,"inq_dimID_nbounds")

status = NF90_INQUIRE_DIMENSION(fidT,dimID_depth,len=mz); call erreur(status,.TRUE.,"inq_dim_depth")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_x,len=mx); call erreur(status,.TRUE.,"inq_dim_lon")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_y,len=my); call erreur(status,.TRUE.,"inq_dim_lat")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_nbounds,len=mnbounds); call erreur(status,.TRUE.,"inq_dim_nbounds")

ALLOCATE(  temperature(mx,my,mz)  ) 
ALLOCATE(  temperature_Con_Shlf(mx,my,mz)  ) 
ALLOCATE(  temperature_Offshore(mx,my,mz)  ) 
ALLOCATE(  z_bnds(mnbounds,mz)  ) 
ALLOCATE(  z(mz)  ) 
ALLOCATE(  lon(mx,my)  ) 
ALLOCATE(  tmplon(mx,my)  ) 
ALLOCATE(  lat(mx,my)  ) 

status = NF90_INQ_VARID(fidT,"temperature",temperature_ID); call erreur(status,.TRUE.,"inq_temperature_ID")
status = NF90_INQ_VARID(fidT,"temperature_Con_Shlf",temperature_Con_Shlf_ID); call erreur(status,.TRUE.,"inq_temperature_Con_Shlf_ID")
status = NF90_INQ_VARID(fidT,"temperature_Offshore",temperature_Offshore_ID); call erreur(status,.TRUE.,"inq_temperature_Offshore_ID")
status = NF90_INQ_VARID(fidT,"z_bnds",z_bnds_ID); call erreur(status,.TRUE.,"inq_z_bnds_ID")
status = NF90_INQ_VARID(fidT,"z",z_ID); call erreur(status,.TRUE.,"inq_z_ID")
status = NF90_INQ_VARID(fidT,"lon",lon_ID); call erreur(status,.TRUE.,"inq_lon_ID")
status = NF90_INQ_VARID(fidT,"lat",lat_ID); call erreur(status,.TRUE.,"inq_lat_ID")
 
status = NF90_GET_VAR(fidT,temperature_ID,temperature); call erreur(status,.TRUE.,"getvar_temperature")
status = NF90_GET_VAR(fidT,temperature_Con_Shlf_ID,temperature_Con_Shlf); call erreur(status,.TRUE.,"getvar_temperature_Con_Shlf")
status = NF90_GET_VAR(fidT,temperature_Offshore_ID,temperature_Offshore); call erreur(status,.TRUE.,"getvar_temperature_Offshore")
status = NF90_GET_VAR(fidT,z_bnds_ID,z_bnds); call erreur(status,.TRUE.,"getvar_z_bnds")
status = NF90_GET_VAR(fidT,z_ID,z); call erreur(status,.TRUE.,"getvar_z")
status = NF90_GET_VAR(fidT,lon_ID,lon); call erreur(status,.TRUE.,"getvar_lon")
status = NF90_GET_VAR(fidT,lat_ID,lat); call erreur(status,.TRUE.,"getvar_lat")
 
status = NF90_CLOSE(fidT); call erreur(status,.TRUE.,"close_file")

write(*,*) '  Minimum longitude for this dataset is : ', MINVAL(lon)

if ( z(2) .lt. 0.0 ) then
  z(:) = -z(:)
endif

mdepth=77
ALLOCATE( depth(mdepth), k1(mdepth), k2(mdepth) )
ALLOCATE( t_an(mx,my,mdepth) )
ALLOCATE( t_sh(mx,my,mdepth) )
ALLOCATE( t_of(mx,my,mdepth) )
depth = (/ 0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., &
       85., 90., 95., 100., 125., 150., 175., 200., 225., 250., 275., 300., 325., 350., 375., &
          400., 425., 450., 475., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., &
          1000., 1050., 1100., 1150., 1200., 1250., 1300., 1350., 1400., 1450., 1500., 1550., &
          1600., 1650., 1700., 1750., 1800., 1850., 1900., 1950., 2000., 2100., 2200., 2300., &
          2400., 2500., 2600., 2700., 2800., 2900., 3000. /)

k1(:) = 1
k2(:) = 1
do kdepth=1,mdepth
  mindif = 1.e4
  do kk=2,mz-1
    if ( abs( depth(kdepth) - z(kk) ) .lt. mindif ) then
      mindif = abs( depth(kdepth) - z(kk) )
      k1(kdepth) = kk
    endif 
  enddo
  if ( abs( depth(kdepth) - z(MIN(mz,k1(kdepth)+1)) ) .lt. abs( depth(kdepth) - z(MAX(1,k1(kdepth)-1)) ) ) then
    k2(kdepth) = MIN(mz,k1(kdepth)+1)
  else
    k2(kdepth) = MAX(1,k1(kdepth)-1)
  endif
  if ( depth(kdepth) .gt. MAX(z(k1(kdepth)),z(k2(kdepth))) ) then
    write(*,*) '## ', depth(kdepth), z(MAX(k1(kdepth),k2(kdepth)))
    t_an(:,:,kdepth) = temperature(:,:,MAX(k1(kdepth),k2(kdepth)))
    t_sh(:,:,kdepth) = temperature_Con_Shlf(:,:,MAX(k1(kdepth),k2(kdepth)))
    t_of(:,:,kdepth) = temperature_Offshore(:,:,MAX(k1(kdepth),k2(kdepth)))
  elseif ( depth(kdepth) .lt. MIN(z(k1(kdepth)),z(k2(kdepth))) ) then
    write(*,*) '## ', depth(kdepth), z(MIN(k1(kdepth),k2(kdepth)))
    t_an(:,:,kdepth) = temperature(:,:,MIN(k1(kdepth),k2(kdepth)))
    t_sh(:,:,kdepth) = temperature_Con_Shlf(:,:,MIN(k1(kdepth),k2(kdepth)))
    t_of(:,:,kdepth) = temperature_Offshore(:,:,MIN(k1(kdepth),k2(kdepth)))
  else
    write(*,*) '## ', depth(kdepth), z(k1(kdepth)), z(k2(kdepth))
    where( abs(temperature(:,:,k1(kdepth))) .lt. 50.0 .and. abs(temperature(:,:,k2(kdepth))) .lt. 50.0 )
      t_an(:,:,kdepth) = (   temperature(:,:,k1(kdepth)) * ( z(k2(kdepth)) - depth(kdepth) )   &
      &                    + temperature(:,:,k2(kdepth)) * ( depth(kdepth) - z(k1(kdepth)) ) ) &
      &                / ( z(k2(kdepth)) - z(k1(kdepth)) )
    elsewhere
      t_an(:,:,kdepth) = NF90_FILL_DOUBLE
    endwhere
    where( abs(temperature_Con_Shlf(:,:,k1(kdepth))) .lt. 50.0 .and. abs(temperature_Con_Shlf(:,:,k2(kdepth))) .lt. 50.0 )
      t_sh(:,:,kdepth) = (   temperature_Con_Shlf(:,:,k1(kdepth)) * ( z(k2(kdepth)) - depth(kdepth) )   &
      &                    + temperature_Con_Shlf(:,:,k2(kdepth)) * ( depth(kdepth) - z(k1(kdepth)) ) ) &
      &                / ( z(k2(kdepth)) - z(k1(kdepth)) )
    elsewhere
      t_sh(:,:,kdepth) = NF90_FILL_DOUBLE
    endwhere
    where( abs(temperature_Offshore(:,:,k1(kdepth))) .lt. 50.0 .and. abs(temperature_Offshore(:,:,k2(kdepth))) .lt. 50.0 )
      t_of(:,:,kdepth) = (   temperature_Offshore(:,:,k1(kdepth)) * ( z(k2(kdepth)) - depth(kdepth) )   &
      &                    + temperature_Offshore(:,:,k2(kdepth)) * ( depth(kdepth) - z(k1(kdepth)) ) ) &
      &                / ( z(k2(kdepth)) - z(k1(kdepth)) )
    elsewhere
      t_of(:,:,kdepth) = NF90_FILL_DOUBLE
    endwhere
  endif
enddo

DEALLOCATE( temperature, temperature_Con_Shlf, temperature_Offshore )

!---------------------------------------
! Read salinity :

ALLOCATE( s_an(mx,my,mdepth) )
ALLOCATE( s_sh(mx,my,mdepth) )
ALLOCATE( s_of(mx,my,mdepth) )

write(*,*) 'Reading ', TRIM(file_in_S)

status = NF90_OPEN(TRIM(file_in_S),0,fidS); call erreur(status,.TRUE.,"read WOA temperature")

ALLOCATE(  salinity(mx,my,mz)  ) 
ALLOCATE(  salinity_Con_Shlf(mx,my,mz)  ) 
ALLOCATE(  salinity_Offshore(mx,my,mz)  ) 

status = NF90_INQ_VARID(fidS,"salinity",salinity_ID); call erreur(status,.TRUE.,"inq_salinity_ID")
status = NF90_INQ_VARID(fidS,"salinity_Con_Shlf",salinity_Con_Shlf_ID); call erreur(status,.TRUE.,"inq_salinity_Con_Shlf_ID")
status = NF90_INQ_VARID(fidS,"salinity_Offshore",salinity_Offshore_ID); call erreur(status,.TRUE.,"inq_salinity_Offshore_ID")
 
status = NF90_GET_VAR(fidS,salinity_ID,salinity); call erreur(status,.TRUE.,"getvar_salinity")
status = NF90_GET_VAR(fidS,salinity_Con_Shlf_ID,salinity_Con_Shlf); call erreur(status,.TRUE.,"getvar_salinity")
status = NF90_GET_VAR(fidS,salinity_Offshore_ID,salinity_Offshore); call erreur(status,.TRUE.,"getvar_salinity")
 
status = NF90_CLOSE(fidS); call erreur(status,.TRUE.,"close_file")

do kdepth=1,mdepth
  if ( depth(kdepth) .gt. MAX(z(k1(kdepth)),z(k2(kdepth))) ) then
    s_an(:,:,kdepth) = salinity(:,:,MAX(k1(kdepth),k2(kdepth)))
    s_sh(:,:,kdepth) = salinity_Con_Shlf(:,:,MAX(k1(kdepth),k2(kdepth)))
    s_of(:,:,kdepth) = salinity_Offshore(:,:,MAX(k1(kdepth),k2(kdepth)))
  elseif ( depth(kdepth) .lt. MIN(z(k1(kdepth)),z(k2(kdepth))) ) then
    s_an(:,:,kdepth) = salinity(:,:,MIN(k1(kdepth),k2(kdepth)))
    s_sh(:,:,kdepth) = salinity_Con_Shlf(:,:,MIN(k1(kdepth),k2(kdepth)))
    s_of(:,:,kdepth) = salinity_Offshore(:,:,MIN(k1(kdepth),k2(kdepth)))
  else
    where( abs(salinity(:,:,k1(kdepth))) .lt. 50.0 .and. abs(salinity(:,:,k2(kdepth))) .lt. 50.0 )
      s_an(:,:,kdepth) = (   salinity(:,:,k1(kdepth)) * ( z(k2(kdepth)) - depth(kdepth) )   &
      &                    + salinity(:,:,k2(kdepth)) * ( depth(kdepth) - z(k1(kdepth)) ) ) &
      &                / ( z(k2(kdepth)) - z(k1(kdepth)) )
    elsewhere
      s_an(:,:,kdepth) = NF90_FILL_DOUBLE
    endwhere
    where( abs(salinity_Con_Shlf(:,:,k1(kdepth))) .lt. 50.0 .and. abs(salinity_Con_Shlf(:,:,k2(kdepth))) .lt. 50.0 )
      s_sh(:,:,kdepth) = (   salinity_Con_Shlf(:,:,k1(kdepth)) * ( z(k2(kdepth)) - depth(kdepth) )   &
      &                    + salinity_Con_Shlf(:,:,k2(kdepth)) * ( depth(kdepth) - z(k1(kdepth)) ) ) &
      &                / ( z(k2(kdepth)) - z(k1(kdepth)) )
    elsewhere
      s_sh(:,:,kdepth) = NF90_FILL_DOUBLE
    endwhere
    where( abs(salinity_Offshore(:,:,k1(kdepth))) .lt. 50.0 .and. abs(salinity_Offshore(:,:,k2(kdepth))) .lt. 50.0 )
      s_of(:,:,kdepth) = (   salinity_Offshore(:,:,k1(kdepth)) * ( z(k2(kdepth)) - depth(kdepth) )   &
      &                    + salinity_Offshore(:,:,k2(kdepth)) * ( depth(kdepth) - z(k1(kdepth)) ) ) &
      &                / ( z(k2(kdepth)) - z(k1(kdepth)) )
    elsewhere
      s_of(:,:,kdepth) = NF90_FILL_DOUBLE
    endwhere
  endif
enddo

DEALLOCATE( salinity, salinity_Con_Shlf, salinity_Offshore )

!---------------------------------------
! Convert to potential temperature :

ALLOCATE( theta_an(mx,my,mdepth)  )
ALLOCATE( theta_sh(mx,my,mdepth)  )
ALLOCATE( theta_of(mx,my,mdepth)  )

write(*,*) 'Converting in-situ temperatures to potential temperatures...'

do ii=1,mx
do jj=1,my
do kk=1,mdepth
 Sabs = gsw_sa_from_sp( DBLE(s_an(ii,jj,kk)), DBLE(depth(kk)), DBLE(lon(ii,jj)), DBLE(lat(ii,jj)) )
 theta_an(ii,jj,kk) = gsw_pt0_from_t(Sabs, DBLE(t_an(ii,jj,kk)), DBLE(depth(kk)))
 Sabs = gsw_sa_from_sp( DBLE(s_sh(ii,jj,kk)), DBLE(depth(kk)), DBLE(lon(ii,jj)), DBLE(lat(ii,jj)) )
 theta_sh(ii,jj,kk) = gsw_pt0_from_t(Sabs, DBLE(t_sh(ii,jj,kk)), DBLE(depth(kk)))
 Sabs = gsw_sa_from_sp( DBLE(s_of(ii,jj,kk)), DBLE(depth(kk)), DBLE(lon(ii,jj)), DBLE(lat(ii,jj)) )
 theta_of(ii,jj,kk) = gsw_pt0_from_t(Sabs, DBLE(t_of(ii,jj,kk)), DBLE(depth(kk)))
enddo
enddo
enddo

write(*,*) '         >> done'

DEALLOCATE( t_an, t_sh, t_of )

!---------------------------------------                   
! Read ice shelf mask file :                                 

write(*,*) 'Reading ', TRIM(file_isf)

status = NF90_OPEN(TRIM(file_isf),0,fidISF)          
call erreur(status,.TRUE.,"read ice shelf mask file") 

status = NF90_INQ_DIMID(fidISF,"Nisf",dimID_Nisf)
call erreur(status,.TRUE.,"inq_dimID_Nisf")

status = NF90_INQUIRE_DIMENSION(fidISF,dimID_Nisf,len=Nisf)
call erreur(status,.TRUE.,"inq_dim_Nisf")

ALLOCATE(  front_max_lat(Nisf)  ) 
ALLOCATE(  front_min_lat(Nisf)  ) 
ALLOCATE(  front_max_lon(Nisf)  ) 
ALLOCATE(  front_min_lon(Nisf)  ) 
ALLOCATE(  front_ice_dep_avg(Nisf)  ) 
ALLOCATE(  front_ice_dep_min(Nisf)  ) 
ALLOCATE(  front_bot_dep_avg(Nisf)  ) 
ALLOCATE(  front_bot_dep_max(Nisf)  ) 
ALLOCATE(  melt_isf(Nisf)  )

status = NF90_INQ_VARID(fidISF,"front_max_lat",front_max_lat_ID)
call erreur(status,.TRUE.,"inq_front_max_lat_ID")
status = NF90_INQ_VARID(fidISF,"front_min_lat",front_min_lat_ID)
call erreur(status,.TRUE.,"inq_front_min_lat_ID")
status = NF90_INQ_VARID(fidISF,"front_max_lon",front_max_lon_ID)
call erreur(status,.TRUE.,"inq_front_max_lon_ID")
status = NF90_INQ_VARID(fidISF,"front_min_lon",front_min_lon_ID)
call erreur(status,.TRUE.,"inq_front_min_lon_ID")
status = NF90_INQ_VARID(fidISF,"front_ice_dep_avg",front_ice_dep_avg_ID)
call erreur(status,.TRUE.,"inq_front_ice_dep_avg_ID")
status = NF90_INQ_VARID(fidISF,"front_ice_dep_min",front_ice_dep_min_ID)
call erreur(status,.TRUE.,"inq_front_ice_dep_min_ID")
status = NF90_INQ_VARID(fidISF,"front_bot_dep_avg",front_bot_dep_avg_ID)
call erreur(status,.TRUE.,"inq_front_bot_dep_avg_ID")
status = NF90_INQ_VARID(fidISF,"front_bot_dep_max",front_bot_dep_max_ID)
call erreur(status,.TRUE.,"inq_front_bot_dep_max_ID")
status = NF90_INQ_VARID(fidISF,"melt_isf",melt_isf_ID)
call erreur(status,.TRUE.,"inq_melt_isf_ID")

status = NF90_GET_VAR(fidISF,front_max_lat_ID,front_max_lat)
call erreur(status,.TRUE.,"getvar_front_max_lat")
status = NF90_GET_VAR(fidISF,front_min_lat_ID,front_min_lat)
call erreur(status,.TRUE.,"getvar_front_min_lat")
status = NF90_GET_VAR(fidISF,front_max_lon_ID,front_max_lon)
call erreur(status,.TRUE.,"getvar_front_max_lon")
status = NF90_GET_VAR(fidISF,front_min_lon_ID,front_min_lon)
call erreur(status,.TRUE.,"getvar_front_min_lon")
status = NF90_GET_VAR(fidISF,front_ice_dep_avg_ID,front_ice_dep_avg)
call erreur(status,.TRUE.,"getvar_front_ice_dep_avg")
status = NF90_GET_VAR(fidISF,front_ice_dep_min_ID,front_ice_dep_min)
call erreur(status,.TRUE.,"getvar_front_ice_dep_min")
status = NF90_GET_VAR(fidISF,front_bot_dep_avg_ID,front_bot_dep_avg)
call erreur(status,.TRUE.,"getvar_front_bot_dep_avg")
status = NF90_GET_VAR(fidISF,front_bot_dep_max_ID,front_bot_dep_max)
call erreur(status,.TRUE.,"getvar_front_bot_dep_max")
status = NF90_GET_VAR(fidISF,melt_isf_ID,melt_isf)
call erreur(status,.TRUE.,"getvar_melt_isf")

status = NF90_CLOSE(fidISF)                      
call erreur(status,.TRUE.,"fin_lecture")     

!------------------------------------------------------
! Calculate T,S profiles for individual ice shelves :

! find level closest to 100m depth :
mindif = 1.e4
do kWOA=1,mdepth
  dif = abs( depth(kWOA) - 100.0 )
  if ( dif .lt. mindif ) then
    mindif = dif
    k01 = kWOA
  endif
enddo
! find level closest to 3000m depth :
mindif = 1.e4
do kWOA=1,mdepth
  dif = abs( depth(kWOA) - 3000.0 )
  if ( dif .lt. mindif ) then
    mindif = dif
    k30 = kWOA
  endif
enddo
write(*,*) 'Index corresponding to 3000m depth : ', k30, '      check: ', depth(k30)

ALLOCATE( NN(k30) )
ALLOCATE( mean_t_an_1(Nisf,k30), mean_s_an_1(Nisf,k30) )
ALLOCATE( mean_t_an_2(Nisf,k30), mean_s_an_2(Nisf,k30) )
ALLOCATE( mean_t_an_3(Nisf,k30), mean_s_an_3(Nisf,k30) )

mean_t_an_1(:,:) = 0.e0
mean_s_an_1(:,:) = 0.e0
mean_t_an_2(:,:) = 0.e0
mean_s_an_2(:,:) = 0.e0
mean_t_an_3(:,:) = 0.e0
mean_s_an_3(:,:) = 0.e0

DO kisf=1,Nisf

 IF (       front_min_lon(kisf) .ne. 180.0 .and. front_max_lon(kisf) .ne. -180.0 &      ! to remove non-attributed ice shelves
 &    .and. front_min_lat(kisf) .ne. 90.00 .and. front_max_lat(kisf) .ne. -90.00 ) THEN

  if ( kisf .ne. 10 ) then ! not Ross IS
    tmplon(:,:) = lon(:,:)
  else ! ROSS
    where ( lon(:,:) .lt. 0.0 )
      tmplon(:,:) = lon(:,:) + 360.0
    elsewhere
      tmplon(:,:) = lon(:,:)  
    endwhere
    front_max_lon(kisf) = front_max_lon(kisf) + 360.0
    write(*,*) 'Ross front_max_lon = ', front_max_lon(kisf)
  endif

  ! closest to the ice shelf front within 2 deg in lon and lat :
  NN(:) = 0.e0
  ! trick not to llok for T,S on the other side of the Peninsula :
  if ( front_max_lon(kisf) .ge. -68.0 .and. front_max_lon(kisf) .le. -66.0 ) then
     front_sup = -66.0
  else
     front_sup = front_max_lon(kisf) + lon_box_1
  endif
  if ( front_min_lon(kisf) .ge. -66.0 .and. front_min_lon(kisf) .le. -64.0 ) then
     front_inf = -66.0
  else
     front_inf = front_min_lon(kisf) - lon_box_1
  endif
  do iWOA=1,mx
  do jWOA=1,my
    if (       tmplon(iWOA,jWOA) .ge. front_inf                        &
    &    .and. tmplon(iWOA,jWOA) .le. front_sup                        &
         .and.    lat(iWOA,jWOA) .ge. front_min_lat(kisf) - lat_box_1  &
         .and.    lat(iWOA,jWOA) .le. front_max_lat(kisf) + lat_box_1  ) then
        do kWOA=1,k30 
          if (       theta_sh(iWOA,jWOA, 1  ) .ge.  -5.0 .and. theta_sh(iWOA,jWOA, 1  ) .lt.  33.0 &  ! to remove ice shelf cavities
          &    .and. theta_sh(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_sh(iWOA,jWOA,kWOA) .lt.  33.0 &
          &    .and.     s_sh(iWOA,jWOA,kWOA) .ge.   0.0 .and.     s_sh(iWOA,jWOA,kWOA) .lt.  60.0 ) then
            mean_t_an_1(kisf,kWOA) = mean_t_an_1(kisf,kWOA) + theta_sh(iWOA,jWOA,kWOA)
            mean_s_an_1(kisf,kWOA) = mean_s_an_1(kisf,kWOA) + s_sh(iWOA,jWOA,kWOA)
            NN(kWOA) = NN(kWOA) + 1
          endif
        enddo
    endif
  enddo
  enddo
  do kWOA=1,k30
    if ( NN(kWOA) .ne. 0 ) then
      mean_t_an_1(kisf,kWOA) = mean_t_an_1(kisf,kWOA) / NN(kWOA)
      mean_s_an_1(kisf,kWOA) = mean_s_an_1(kisf,kWOA) / NN(kWOA)
    else
      mean_t_an_1(kisf,kWOA) = NF90_FILL_DOUBLE
      mean_s_an_1(kisf,kWOA) = NF90_FILL_DOUBLE
    endif
  enddo
  !! filling gaps :
  !do ks=1,5
  !do kr=ks,1,-1
  !  do kWOA=2,k30-kr
  !    if ( NN(kWOA) .eq. 0 .and. NN(kWOA+kr) .gt. 0 .and. NN(kWOA-1) .gt. 0 ) then 
  !      mean_t_an_1(kisf,kWOA) = ( mean_t_an_1(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_an_1(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
  !      mean_s_an_1(kisf,kWOA) = ( mean_s_an_1(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_an_1(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
  !      NN(kWOA)=1
  !    endif
  !  enddo
  !enddo
  !enddo
  !- extrapolate downward only if there exists a non-zero value below 100m
  do kr=1,k30-k01 
    do kWOA=k01+1,k30
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_1(kisf,kWOA) = mean_t_an_1(kisf,kWOA-1)
        mean_s_an_1(kisf,kWOA) = mean_s_an_1(kisf,kWOA-1)
        NN(kWOA)=1
      endif
    enddo
  enddo

  ! within 5deg from the ice shelf front but only on the continental shelf (sea floor < 1500m)
  NN(:) = 0.e0
  ! trick not to llok for T,S on the other side of the Peninsula :
  if ( front_max_lon(kisf) .ge. -71.0 .and. front_max_lon(kisf) .le. -66.0 ) then
     front_sup = -66.0
  else
     front_sup = front_max_lon(kisf) + lon_box_2
  endif
  if ( front_min_lon(kisf) .ge. -66.0 .and. front_min_lon(kisf) .le. -61.0 ) then
     front_inf = -66.0
  else
     front_inf = front_min_lon(kisf) - lon_box_2
  endif
  do iWOA=1,mx
  do jWOA=1,my
    if (       tmplon(iWOA,jWOA) .ge. front_inf                        &
    &    .and. tmplon(iWOA,jWOA) .le. front_sup                        &
    &    .and.    lat(iWOA,jWOA) .ge. front_min_lat(kisf) - lat_box_2  &
    &    .and.    lat(iWOA,jWOA) .le. front_max_lat(kisf) + lat_box_2  ) then
        do kWOA=1,k30 
          if (         theta_sh(iWOA,jWOA, 1  ) .ge.  -5.0 .and. theta_sh(iWOA,jWOA, 1  ) .lt.  33.0   &  ! to remove ice shelf cavities
          &    .and.   theta_sh(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_sh(iWOA,jWOA,kWOA) .lt.  33.0   &
          &    .and.       s_sh(iWOA,jWOA,kWOA) .ge.   0.0 .and.     s_sh(iWOA,jWOA,kWOA) .lt.  60.0   ) then
            mean_t_an_2(kisf,kWOA) = mean_t_an_2(kisf,kWOA) + theta_sh(iWOA,jWOA,kWOA)
            mean_s_an_2(kisf,kWOA) = mean_s_an_2(kisf,kWOA) +     s_sh(iWOA,jWOA,kWOA)
            NN(kWOA) = NN(kWOA) + 1
          endif
        enddo
    endif
  enddo
  enddo
  do kWOA=1,k30
    if ( NN(kWOA) .ne. 0 ) then
      mean_t_an_2(kisf,kWOA) = mean_t_an_2(kisf,kWOA) / NN(kWOA)
      mean_s_an_2(kisf,kWOA) = mean_s_an_2(kisf,kWOA) / NN(kWOA)
    else
      mean_t_an_2(kisf,kWOA) = NF90_FILL_DOUBLE
      mean_s_an_2(kisf,kWOA) = NF90_FILL_DOUBLE
    endif
  enddo
  !! filling gaps :
  !do ks=1,5
  !do kr=ks,1,-1
  !  do kWOA=2,k30-kr
  !    if ( NN(kWOA) .eq. 0 .and. NN(kWOA+kr) .gt. 0 .and. NN(kWOA-1) .gt. 0 ) then 
  !      mean_t_an_2(kisf,kWOA) = ( mean_t_an_2(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_an_2(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
  !      mean_s_an_2(kisf,kWOA) = ( mean_s_an_2(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_an_2(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
  !      NN(kWOA)=1
  !    endif
  !  enddo
  !enddo
  !enddo
  !- extrapolate downward only if there exists a non-zero value below 100m
  do kr=1,k30-k01+1
    do kWOA=k01+1,k30
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_2(kisf,kWOA) = mean_t_an_2(kisf,kWOA-1)
        mean_s_an_2(kisf,kWOA) = mean_s_an_2(kisf,kWOA-1)
        NN(kWOA)=1
      endif
    enddo
  enddo

  ! within 10deg from the ice shelf front but only offshore of the continental shelf (sea floor > 1500m)
  NN(:) = 0.e0
  ! trick not to llok for T,S on the other side of the Peninsula :
  if ( front_max_lon(kisf) .ge. -76.0 .and. front_max_lon(kisf) .le. -66.0 ) then
     front_sup = -66.0
  else
     front_sup = front_max_lon(kisf) + lon_box_3
  endif
  if ( front_min_lon(kisf) .ge. -66.0 .and. front_min_lon(kisf) .le. -56.0 ) then
     front_inf = -66.0
  else
     front_inf = front_min_lon(kisf) - lon_box_3
  endif
  do iWOA=1,mx
  do jWOA=1,my
    if (       tmplon(iWOA,jWOA) .ge. front_inf                        &
    &    .and. tmplon(iWOA,jWOA) .le. front_sup                        &
    &    .and.    lat(iWOA,jWOA) .ge. front_min_lat(kisf) - lat_box_3  &
    &    .and.    lat(iWOA,jWOA) .le. front_max_lat(kisf) + lat_box_3  ) then
        do kWOA=1,k30 
          if (       theta_of(iWOA,jWOA, 1  ) .ge.  -5.0 .and. theta_of(iWOA,jWOA, 1  ) .lt.  33.0 &  ! to remove ice shelf cavities
          &    .and. theta_of(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_of(iWOA,jWOA,kWOA) .lt.  33.0 &
          &    .and.     s_of(iWOA,jWOA,kWOA) .ge.   0.0 .and.     s_of(iWOA,jWOA,kWOA) .lt.  60.0 ) then
            mean_t_an_3(kisf,kWOA) = mean_t_an_3(kisf,kWOA) + theta_of(iWOA,jWOA,kWOA)
            mean_s_an_3(kisf,kWOA) = mean_s_an_3(kisf,kWOA) +     s_of(iWOA,jWOA,kWOA)
            NN(kWOA) = NN(kWOA) + 1
          endif
        enddo
    endif
  enddo
  enddo
  do kWOA=1,k30
    if ( NN(kWOA) .ne. 0 ) then
      mean_t_an_3(kisf,kWOA) = mean_t_an_3(kisf,kWOA) / NN(kWOA)
      mean_s_an_3(kisf,kWOA) = mean_s_an_3(kisf,kWOA) / NN(kWOA)
    else
      mean_t_an_3(kisf,kWOA) = NF90_FILL_DOUBLE
      mean_s_an_3(kisf,kWOA) = NF90_FILL_DOUBLE
    endif
  enddo
  !! filling gaps :
  !do ks=1,5
  !do kr=ks,1,-1
  !  do kWOA=2,k30-kr
  !    if ( NN(kWOA) .eq. 0 .and. NN(kWOA+kr) .gt. 0 .and. NN(kWOA-1) .gt. 0 ) then 
  !      mean_t_an_3(kisf,kWOA) = ( mean_t_an_3(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_an_3(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
  !      mean_s_an_3(kisf,kWOA) = ( mean_s_an_3(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_an_3(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
  !      NN(kWOA)=1
  !    endif
  !  enddo
  !enddo
  !enddo
  !- extrapolate downward only if there exists a non-zero value below 100m
  do kr=1,k30-k01+1
    do kWOA=k01+1,k30
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_3(kisf,kWOA) = mean_t_an_3(kisf,kWOA-1)
        mean_s_an_3(kisf,kWOA) = mean_s_an_3(kisf,kWOA-1)
        NN(kWOA)=1
      endif
    enddo
  enddo

 ELSE

  mean_t_an_1(kisf,:) = NF90_FILL_DOUBLE
  mean_s_an_1(kisf,:) = NF90_FILL_DOUBLE
  mean_t_an_2(kisf,:) = NF90_FILL_DOUBLE
  mean_s_an_2(kisf,:) = NF90_FILL_DOUBLE
  mean_t_an_3(kisf,:) = NF90_FILL_DOUBLE 
  mean_s_an_3(kisf,:) = NF90_FILL_DOUBLE 

 ENDIF

ENDDO


!---------------------------------------
! Writing new netcdf file :                                   

write(*,*) 'Creating ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create output file')

status = NF90_DEF_DIM(fidM,"Nisf",Nisf,dimID_Nisf); call erreur(status,.TRUE.,"def_dimID_Nisf")
status = NF90_DEF_DIM(fidM,"depth",k30,dimID_depth); call erreur(status,.TRUE.,"def_dimID_depth")
status = NF90_DEF_DIM(fidM,"nbounds",mnbounds,dimID_nbounds); call erreur(status,.TRUE.,"def_dimID_nbounds")

!status = NF90_DEF_VAR(fidM,"depth_bnds",NF90_DOUBLE,(/dimID_nbounds,dimID_depth/),depth_bnds_ID); call erreur(status,.TRUE.,"def_var_depth_bnds_ID")
status = NF90_DEF_VAR(fidM,"depth",NF90_DOUBLE,(/dimID_depth/),depth_ID); call erreur(status,.TRUE.,"def_var_depth_ID")
status = NF90_DEF_VAR(fidM,"t_an_1",NF90_DOUBLE,(/dimID_Nisf,dimID_depth/),mean_t_an_1_ID); call erreur(status,.TRUE.,"def_var_mean_t_an_1_ID")
status = NF90_DEF_VAR(fidM,"s_an_1",NF90_DOUBLE,(/dimID_Nisf,dimID_depth/),mean_s_an_1_ID); call erreur(status,.TRUE.,"def_var_mean_s_an_1_ID")
status = NF90_DEF_VAR(fidM,"t_an_2",NF90_DOUBLE,(/dimID_Nisf,dimID_depth/),mean_t_an_2_ID); call erreur(status,.TRUE.,"def_var_mean_t_an_2_ID")
status = NF90_DEF_VAR(fidM,"s_an_2",NF90_DOUBLE,(/dimID_Nisf,dimID_depth/),mean_s_an_2_ID); call erreur(status,.TRUE.,"def_var_mean_s_an_2_ID")
status = NF90_DEF_VAR(fidM,"t_an_3",NF90_DOUBLE,(/dimID_Nisf,dimID_depth/),mean_t_an_3_ID); call erreur(status,.TRUE.,"def_var_mean_t_an_3_ID")
status = NF90_DEF_VAR(fidM,"s_an_3",NF90_DOUBLE,(/dimID_Nisf,dimID_depth/),mean_s_an_3_ID); call erreur(status,.TRUE.,"def_var_mean_s_an_3_ID")

!status = NF90_PUT_ATT(fidM,depth_bnds_ID,"comment","depth bounds"); call erreur(status,.TRUE.,"put_att_depth_bnds_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"axis","Z"); call erreur(status,.TRUE.,"put_att_depth_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"units","meters"); call erreur(status,.TRUE.,"put_att_depth_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"positive","down"); call erreur(status,.TRUE.,"put_att_depth_ID")
!status = NF90_PUT_ATT(fidM,depth_ID,"bounds","depth_bnds"); call erreur(status,.TRUE.,"put_att_depth_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"standard_name","depth"); call erreur(status,.TRUE.,"put_att_depth_ID")

status = NF90_PUT_ATT(fidM,mean_t_an_1_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_1_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_1_ID,"long_name","Potential temperature in box_1 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_1_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_1_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_1_ID,"long_name","Practical Salinity in box_1 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")

status = NF90_PUT_ATT(fidM,mean_t_an_2_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_2_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_2_ID,"long_name","Potential temperature in box_2 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_2_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_2_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_2_ID,"long_name","Practical Salinity in box_2 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")

status = NF90_PUT_ATT(fidM,mean_t_an_3_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_3_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_3_ID,"long_name","Potential temperature in box_3 seaward of the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_3_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_3_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_3_ID,"long_name","Practical Salinity in box_3 seaward of the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using calculate_TS_profiles_WOA_MEOP.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lon_box_1",lon_box_1); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lat_box_1",lat_box_1); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lon_box_2",lon_box_2); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lat_box_2",lat_box_2); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lon_box_3",lon_box_3); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lat_box_3",lat_box_3); call erreur(status,.TRUE.,"put_att_GLOBAL")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 

!status = NF90_PUT_VAR(fidM,depth_bnds_ID,depth_bnds(:,1:k30)); call erreur(status,.TRUE.,"var_depth_bnds_ID")
status = NF90_PUT_VAR(fidM,depth_ID,depth(1:k30));             call erreur(status,.TRUE.,"var_depth_ID")
status = NF90_PUT_VAR(fidM,mean_t_an_1_ID,mean_t_an_1); call erreur(status,.TRUE.,"var_mean_t_an_1_ID")
status = NF90_PUT_VAR(fidM,mean_s_an_1_ID,mean_s_an_1); call erreur(status,.TRUE.,"var_mean_s_an_1_ID")
status = NF90_PUT_VAR(fidM,mean_t_an_2_ID,mean_t_an_2); call erreur(status,.TRUE.,"var_mean_t_an_2_ID")
status = NF90_PUT_VAR(fidM,mean_s_an_2_ID,mean_s_an_2); call erreur(status,.TRUE.,"var_mean_s_an_2_ID")
status = NF90_PUT_VAR(fidM,mean_t_an_3_ID,mean_t_an_3); call erreur(status,.TRUE.,"var_mean_t_an_3_ID")
status = NF90_PUT_VAR(fidM,mean_s_an_3_ID,mean_s_an_3); call erreur(status,.TRUE.,"var_mean_s_an_3_ID")

status = NF90_CLOSE(fidM); call erreur(status,.TRUE.,"final")         

end program modif



SUBROUTINE erreur(iret, lstop, chaine)
! pour les messages d'erreur
USE netcdf
INTEGER, INTENT(in)                     :: iret
LOGICAL, INTENT(in)                     :: lstop
CHARACTER(LEN=*), INTENT(in)            :: chaine
!
CHARACTER(LEN=80)                       :: message
!
IF ( iret .NE. 0 ) THEN
WRITE(*,*) 'ROUTINE: ', TRIM(chaine)
WRITE(*,*) 'ERREUR: ', iret
message=NF90_STRERROR(iret)
WRITE(*,*) 'CA VEUT DIRE:',TRIM(message)
IF ( lstop ) STOP
ENDIF
!
END SUBROUTINE erreur
