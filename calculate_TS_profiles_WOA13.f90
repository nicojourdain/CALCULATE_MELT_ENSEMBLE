program modif                                         

USE netcdf                                            

use gsw_mod_kinds
use gsw_mod_netcdf
use gsw_mod_toolbox

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

IMPLICIT NONE                                         

INTEGER :: fidT, dimID_depth, dimID_lon, dimID_lat, dimID_nbounds, mdepth, mlon, mlat, mnbounds, t_se_ID, t_an_ID, s_se_ID, s_an_ID, fidM, &
&          depth_bnds_ID, depth_ID, lon_bnds_ID, lon_ID, lat_bnds_ID, lat_ID, fidS, k15, k30, iWOA, jWOA, kWOA, j50S, kisf, kr, ks, k01,   &
&          mean_s_an_1_ID, mean_s_er_1_ID, mean_t_an_1_ID, mean_t_er_1_ID, mean_s_an_2_ID, mean_s_er_2_ID, mean_t_an_2_ID, mean_t_er_2_ID, &
&          mean_s_an_3_ID, mean_s_er_3_ID, mean_t_an_3_ID, mean_t_er_3_ID, fidISF, status, dimID_Nisf, Nisf, front_max_lat_ID,             &
&          front_min_lat_ID, front_max_lon_ID, front_min_lon_ID, front_ice_dep_avg_ID, front_ice_dep_min_ID, front_bot_dep_avg_ID,         &
&          front_bot_dep_max_ID, melt_isf_ID, ii, jj, kk

CHARACTER(LEN=150) :: file_isf, file_out, file_in_T, file_in_S

INTEGER, ALLOCATABLE, DIMENSION(:) :: NN

REAL*4,ALLOCATABLE,DIMENSION(:) :: front_max_lat, front_min_lat, front_max_lon, front_min_lon, front_ice_dep_avg, front_ice_dep_min, &
&                                  front_bot_dep_avg, front_bot_dep_max, depth, lon, lat, tmplon, melt_isf
 
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: depth_bnds, lon_bnds, lat_bnds

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: mean_t_an_1, mean_t_er_1, mean_s_an_1, mean_s_er_1, &
&                                    mean_t_an_2, mean_t_er_2, mean_s_an_2, mean_s_er_2, &
&                                    mean_t_an_3, mean_t_er_3, mean_s_an_3, mean_s_er_3 
 
REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: t_se, t_an, s_se, s_an, theta_an, theta_se

REAL*4 :: dif, mindif, front_inf, front_sup, &
&         lon_box_1, lon_box_2, lon_box_3, lat_box_1, lat_box_2, lat_box_3

REAL*8 :: Sabs
 
file_in_T  = '/store/njourd/DATA/WOA_2013/v2/woa13_v2_025_T_annual_mean.nc'
file_in_S  = '/store/njourd/DATA/WOA_2013/v2/woa13_v2_025_S_annual_mean.nc'
file_isf   = 'RTopo-2.0.1_30sec_ice_shelf_mask_ANT3.nc'
file_out   = 'TS_per_ice_shelf_ANT3_WOA13.nc'

! NB: 5.0 x 1.75 is the effective resolution at 70°S for a model of 1° resolution in longitude (assuming 5 delta X and a Mercator grid)
lon_box_1 =  5.0 ; lat_box_1 = 1.75 ! Continental shelf only
lon_box_2 = 10.0 ; lat_box_2 = 3.50 ! Continental shelf only
lon_box_3 = 20.0 ; lat_box_3 = 7.00 ! Offshore only

call gsw_saar_init (.true.)
 
!---------------------------------------
! Read WOA2013 temperature :

write(*,*) 'Reading ', TRIM(file_in_T)

status = NF90_OPEN(TRIM(file_in_T),0,fidT); call erreur(status,.TRUE.,"read WOA temperature")

status = NF90_INQ_DIMID(fidT,"depth",dimID_depth); call erreur(status,.TRUE.,"inq_dimID_depth")
status = NF90_INQ_DIMID(fidT,"lon",dimID_lon); call erreur(status,.TRUE.,"inq_dimID_lon")
status = NF90_INQ_DIMID(fidT,"lat",dimID_lat); call erreur(status,.TRUE.,"inq_dimID_lat")
status = NF90_INQ_DIMID(fidT,"nbounds",dimID_nbounds); call erreur(status,.TRUE.,"inq_dimID_nbounds")

status = NF90_INQUIRE_DIMENSION(fidT,dimID_depth,len=mdepth); call erreur(status,.TRUE.,"inq_dim_depth")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_lon,len=mlon); call erreur(status,.TRUE.,"inq_dim_lon")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_lat,len=mlat); call erreur(status,.TRUE.,"inq_dim_lat")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_nbounds,len=mnbounds); call erreur(status,.TRUE.,"inq_dim_nbounds")

ALLOCATE(  t_se(mlon,mlat,mdepth)  ) 
ALLOCATE(  t_an(mlon,mlat,mdepth)  ) 
ALLOCATE(  depth_bnds(mnbounds,mdepth)  ) 
ALLOCATE(  depth(mdepth)  ) 
ALLOCATE(  lon_bnds(mnbounds,mlon)  ) 
ALLOCATE(  lon(mlon)  ) 
ALLOCATE(  tmplon(mlon)  ) 
ALLOCATE(  lat_bnds(mnbounds,mlat)  ) 
ALLOCATE(  lat(mlat)  ) 

status = NF90_INQ_VARID(fidT,"t_se",t_se_ID); call erreur(status,.TRUE.,"inq_t_se_ID")
status = NF90_INQ_VARID(fidT,"t_an",t_an_ID); call erreur(status,.TRUE.,"inq_t_an_ID")
status = NF90_INQ_VARID(fidT,"depth_bnds",depth_bnds_ID); call erreur(status,.TRUE.,"inq_depth_bnds_ID")
status = NF90_INQ_VARID(fidT,"depth",depth_ID); call erreur(status,.TRUE.,"inq_depth_ID")
status = NF90_INQ_VARID(fidT,"lon_bnds",lon_bnds_ID); call erreur(status,.TRUE.,"inq_lon_bnds_ID")
status = NF90_INQ_VARID(fidT,"lon",lon_ID); call erreur(status,.TRUE.,"inq_lon_ID")
status = NF90_INQ_VARID(fidT,"lat_bnds",lat_bnds_ID); call erreur(status,.TRUE.,"inq_lat_bnds_ID")
status = NF90_INQ_VARID(fidT,"lat",lat_ID); call erreur(status,.TRUE.,"inq_lat_ID")
 
status = NF90_GET_VAR(fidT,t_se_ID,t_se); call erreur(status,.TRUE.,"getvar_t_se")
status = NF90_GET_VAR(fidT,t_an_ID,t_an); call erreur(status,.TRUE.,"getvar_t_an")
status = NF90_GET_VAR(fidT,depth_bnds_ID,depth_bnds); call erreur(status,.TRUE.,"getvar_depth_bnds")
status = NF90_GET_VAR(fidT,depth_ID,depth); call erreur(status,.TRUE.,"getvar_depth")
status = NF90_GET_VAR(fidT,lon_bnds_ID,lon_bnds); call erreur(status,.TRUE.,"getvar_lon_bnds")
status = NF90_GET_VAR(fidT,lon_ID,lon); call erreur(status,.TRUE.,"getvar_lon")
status = NF90_GET_VAR(fidT,lat_bnds_ID,lat_bnds); call erreur(status,.TRUE.,"getvar_lat_bnds")
status = NF90_GET_VAR(fidT,lat_ID,lat); call erreur(status,.TRUE.,"getvar_lat")
 
status = NF90_CLOSE(fidT); call erreur(status,.TRUE.,"close_file")

write(*,*) '  Minimum longitude for WOA2013 is : ', MINVAL(lon)

!---------------------------------------
! Read WOA2013 salinity :

write(*,*) 'Reading ', TRIM(file_in_S)

status = NF90_OPEN(TRIM(file_in_S),0,fidS); call erreur(status,.TRUE.,"read WOA temperature")

ALLOCATE(  s_se(mlon,mlat,mdepth)  ) 
ALLOCATE(  s_an(mlon,mlat,mdepth)  ) 

status = NF90_INQ_VARID(fidS,"s_se",s_se_ID); call erreur(status,.TRUE.,"inq_s_se_ID")
status = NF90_INQ_VARID(fidS,"s_an",s_an_ID); call erreur(status,.TRUE.,"inq_s_an_ID")
 
status = NF90_GET_VAR(fidS,s_se_ID,s_se); call erreur(status,.TRUE.,"getvar_s_se")
status = NF90_GET_VAR(fidS,s_an_ID,s_an); call erreur(status,.TRUE.,"getvar_s_an")
 
status = NF90_CLOSE(fidS); call erreur(status,.TRUE.,"close_file")

!---------------------------------------
! Convert to potential temperature :

ALLOCATE( theta_an(mlon,mlat,mdepth)  )
ALLOCATE( theta_se(mlon,mlat,mdepth)  )

write(*,*) 'Converting in-situ temperatures to potential temperatures...'

do ii=1,mlon
do jj=1,mlat
do kk=1,mdepth
 Sabs = gsw_sa_from_sp( DBLE(s_an(ii,jj,kk)), DBLE(depth(kk)), DBLE(lon(ii)), DBLE(lat(jj)) )
 theta_an(ii,jj,kk) = gsw_pt0_from_t(Sabs, DBLE(t_an(ii,jj,kk)), DBLE(depth(kk)))
 theta_se(ii,jj,kk) = gsw_pt0_from_t(Sabs, DBLE(t_an(ii,jj,kk)+t_se(ii,jj,kk)), DBLE(depth(kk)))
 theta_se(ii,jj,kk) = theta_se(ii,jj,kk) - theta_an(ii,jj,kk)
enddo
enddo
enddo

write(*,*) '         >> done'

DEALLOCATE( t_an, t_se )

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

! find WOA13 level closest to 100m depth :
mindif = 1.e4
do kWOA=1,mdepth
  dif = abs( depth(kWOA) - 100.0 )
  if ( dif .lt. mindif ) then
    mindif = dif
    k01 = kWOA
  endif
enddo
! find WOA13 level closest to 1500m depth :
mindif = 1.e4
do kWOA=1,mdepth
  dif = abs( depth(kWOA) - 1500.0 )
  if ( dif .lt. mindif ) then
    mindif = dif
    k15 = kWOA
  endif
enddo
write(*,*) 'Index corresponding to 1500m depth : ', k15, '      check: ', depth(k15)
! find WOA13 level closest to 3000m depth :
mindif = 1.e4
do kWOA=1,mdepth
  dif = abs( depth(kWOA) - 3000.0 )
  if ( dif .lt. mindif ) then
    mindif = dif
    k30 = kWOA
  endif
enddo
write(*,*) 'Index corresponding to 3000m depth : ', k30, '      check: ', depth(k30)
! find WOA13's 50degS :
mindif = 1.e4
do jWOA=1,mlat
  dif = abs( lat(jWOA) + 50.0 )
  if ( dif .lt. mindif ) then
    mindif = dif
    j50S = jWOA
  endif
enddo
write(*,*) 'Index corresponding to 50degS : ', j50S, '      check: ', lat(j50S)

ALLOCATE( NN(k30) )
ALLOCATE( mean_t_an_1(Nisf,k30), mean_t_er_1(Nisf,k30), mean_s_an_1(Nisf,k30), mean_s_er_1(Nisf,k30) )
ALLOCATE( mean_t_an_2(Nisf,k30), mean_t_er_2(Nisf,k30), mean_s_an_2(Nisf,k30), mean_s_er_2(Nisf,k30) )
ALLOCATE( mean_t_an_3(Nisf,k30), mean_t_er_3(Nisf,k30), mean_s_an_3(Nisf,k30), mean_s_er_3(Nisf,k30) )

mean_t_an_1(:,:) = 0.e0
mean_t_er_1(:,:) = 0.e0
mean_s_an_1(:,:) = 0.e0
mean_s_er_1(:,:) = 0.e0
mean_t_an_2(:,:) = 0.e0
mean_t_er_2(:,:) = 0.e0
mean_s_an_2(:,:) = 0.e0
mean_s_er_2(:,:) = 0.e0
mean_t_an_3(:,:) = 0.e0
mean_t_er_3(:,:) = 0.e0
mean_s_an_3(:,:) = 0.e0
mean_s_er_3(:,:) = 0.e0

DO kisf=1,Nisf

 IF (       front_min_lon(kisf) .ne. 180.0 .and. front_max_lon(kisf) .ne. -180.0 &      ! to remove non-attributed ice shelves
 &    .and. front_min_lat(kisf) .ne. 90.00 .and. front_max_lat(kisf) .ne. -90.00 ) THEN

  if ( kisf .ne. 10 ) then ! not Ross IS
    tmplon(:) = lon(:)
  else ! ROSS
    where ( lon(:) .lt. 0.0 )
      tmplon(:) = lon(:) + 360.0
    elsewhere
      tmplon(:) = lon(:)  
    endwhere
    front_max_lon(kisf) = front_max_lon(kisf) + 360.0
    write(*,*) 'Ross front_max_lon = ', front_max_lon(kisf)
  endif

  ! closest to the ice shelf front within 2 deg in lon and lat :
  NN(:) = 0.e0
  do iWOA=1,mlon
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
    if (       tmplon(iWOA) .ge. front_inf  &
    &    .and. tmplon(iWOA) .le. front_sup  ) then
      do jWOA=1,j50S
         if (       lat(jWOA) .ge. front_min_lat(kisf) - lat_box_1  &
         &    .and. lat(jWOA) .le. front_max_lat(kisf) + lat_box_1  ) then
           do kWOA=1,k30 
             if (       theta_an(iWOA,jWOA, 1  ) .ge.  -5.0 .and. theta_an(iWOA,jWOA, 1  ) .lt.  33.0 &  ! to remove ice shelf cavities
             &    .and. theta_an(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_an(iWOA,jWOA,kWOA) .lt.  33.0 &
             &    .and. theta_se(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_se(iWOA,jWOA,kWOA) .lt.  33.0 &
             &    .and. s_an(iWOA,jWOA,kWOA) .ge.   0.0 .and. s_an(iWOA,jWOA,kWOA) .lt.  60.0 &
             &    .and. s_se(iWOA,jWOA,kWOA) .ge.   0.0 .and. s_se(iWOA,jWOA,kWOA) .lt.  60.0 ) then
               mean_t_an_1(kisf,kWOA) = mean_t_an_1(kisf,kWOA) + theta_an(iWOA,jWOA,kWOA)
               mean_t_er_1(kisf,kWOA) = mean_t_er_1(kisf,kWOA) + theta_se(iWOA,jWOA,kWOA)
               mean_s_an_1(kisf,kWOA) = mean_s_an_1(kisf,kWOA) + s_an(iWOA,jWOA,kWOA)
               mean_s_er_1(kisf,kWOA) = mean_s_er_1(kisf,kWOA) + s_se(iWOA,jWOA,kWOA)
               NN(kWOA) = NN(kWOA) + 1
             endif
           enddo
         endif
      enddo
    endif
  enddo
  do kWOA=1,k30
    if ( NN(kWOA) .ne. 0 ) then
      mean_t_an_1(kisf,kWOA) = mean_t_an_1(kisf,kWOA) / NN(kWOA)
      mean_t_er_1(kisf,kWOA) = mean_t_er_1(kisf,kWOA) / NN(kWOA)
      mean_s_an_1(kisf,kWOA) = mean_s_an_1(kisf,kWOA) / NN(kWOA)
      mean_s_er_1(kisf,kWOA) = mean_s_er_1(kisf,kWOA) / NN(kWOA)
    else
      mean_t_an_1(kisf,kWOA) = NF90_FILL_FLOAT
      mean_t_er_1(kisf,kWOA) = NF90_FILL_FLOAT
      mean_s_an_1(kisf,kWOA) = NF90_FILL_FLOAT
      mean_s_er_1(kisf,kWOA) = NF90_FILL_FLOAT
    endif
  enddo
  ! filling gaps :
  do ks=1,5
  do kr=ks,1,-1
    do kWOA=2,k30-kr
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA+kr) .gt. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_1(kisf,kWOA) = ( mean_t_an_1(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_an_1(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_t_er_1(kisf,kWOA) = ( mean_t_er_1(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_er_1(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_s_an_1(kisf,kWOA) = ( mean_s_an_1(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_an_1(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_s_er_1(kisf,kWOA) = ( mean_s_er_1(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_er_1(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        NN(kWOA)=1
      endif
    enddo
  enddo
  enddo
  !- extrapolate downward only if there exists a non-zero value below 100m
  do kr=1,k30-k01 
    do kWOA=k01+1,k30
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_1(kisf,kWOA) = mean_t_an_1(kisf,kWOA-1)
        mean_t_er_1(kisf,kWOA) = mean_t_er_1(kisf,kWOA-1)
        mean_s_an_1(kisf,kWOA) = mean_s_an_1(kisf,kWOA-1)
        mean_s_er_1(kisf,kWOA) = mean_s_er_1(kisf,kWOA-1)
        NN(kWOA)=1
      endif
    enddo
  enddo

  ! within 5deg from the ice shelf front but only on the continental shelf (sea floor < 1500m)
  NN(:) = 0.e0
  do iWOA=1,mlon
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
    if (       tmplon(iWOA) .ge. front_inf  &
    &    .and. tmplon(iWOA) .le. front_sup  ) then
      do jWOA=1,j50S
         if (       lat(jWOA) .ge. front_min_lat(kisf) - lat_box_2  &
         &    .and. lat(jWOA) .le. front_max_lat(kisf) + lat_box_2  ) then
           do kWOA=1,k30 
             if (         theta_an(iWOA,jWOA, 1  ) .ge.  -5.0 .and. theta_an(iWOA,jWOA, 1  ) .lt.  33.0   &  ! to remove ice shelf cavities
             &    .and. ( theta_an(iWOA,jWOA, k15) .lt.  -5.0  .or. theta_an(iWOA,jWOA, k15) .gt.  33.0 ) &  ! to remove pts offshore of the continental shelf
             &    .and.   theta_an(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_an(iWOA,jWOA,kWOA) .lt.  33.0   &
             &    .and.   theta_se(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_se(iWOA,jWOA,kWOA) .lt.  33.0   &
             &    .and.   s_an(iWOA,jWOA,kWOA) .ge.   0.0 .and. s_an(iWOA,jWOA,kWOA) .lt.  60.0   &
             &    .and.   s_se(iWOA,jWOA,kWOA) .ge.   0.0 .and. s_se(iWOA,jWOA,kWOA) .lt.  60.0   ) then
               mean_t_an_2(kisf,kWOA) = mean_t_an_2(kisf,kWOA) + theta_an(iWOA,jWOA,kWOA)
               mean_t_er_2(kisf,kWOA) = mean_t_er_2(kisf,kWOA) + theta_se(iWOA,jWOA,kWOA)
               mean_s_an_2(kisf,kWOA) = mean_s_an_2(kisf,kWOA) + s_an(iWOA,jWOA,kWOA)
               mean_s_er_2(kisf,kWOA) = mean_s_er_2(kisf,kWOA) + s_se(iWOA,jWOA,kWOA)
               NN(kWOA) = NN(kWOA) + 1
             endif
           enddo
         endif
      enddo
    endif
  enddo
  do kWOA=1,k30
    if ( NN(kWOA) .ne. 0 ) then
      mean_t_an_2(kisf,kWOA) = mean_t_an_2(kisf,kWOA) / NN(kWOA)
      mean_t_er_2(kisf,kWOA) = mean_t_er_2(kisf,kWOA) / NN(kWOA)
      mean_s_an_2(kisf,kWOA) = mean_s_an_2(kisf,kWOA) / NN(kWOA)
      mean_s_er_2(kisf,kWOA) = mean_s_er_2(kisf,kWOA) / NN(kWOA)
    else
      mean_t_an_2(kisf,kWOA) = NF90_FILL_FLOAT
      mean_t_er_2(kisf,kWOA) = NF90_FILL_FLOAT
      mean_s_an_2(kisf,kWOA) = NF90_FILL_FLOAT
      mean_s_er_2(kisf,kWOA) = NF90_FILL_FLOAT
    endif
  enddo
  ! filling gaps :
  do ks=1,5
  do kr=ks,1,-1
    do kWOA=2,k30-kr
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA+kr) .gt. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_2(kisf,kWOA) = ( mean_t_an_2(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_an_2(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_t_er_2(kisf,kWOA) = ( mean_t_er_2(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_er_2(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_s_an_2(kisf,kWOA) = ( mean_s_an_2(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_an_2(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_s_er_2(kisf,kWOA) = ( mean_s_er_2(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_er_2(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        NN(kWOA)=1
      endif
    enddo
  enddo
  enddo
  !- extrapolate downward only if there exists a non-zero value below 100m
  do kr=1,k30-k01+1
    do kWOA=k01+1,k30
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_2(kisf,kWOA) = mean_t_an_2(kisf,kWOA-1)
        mean_t_er_2(kisf,kWOA) = mean_t_er_2(kisf,kWOA-1)
        mean_s_an_2(kisf,kWOA) = mean_s_an_2(kisf,kWOA-1)
        mean_s_er_2(kisf,kWOA) = mean_s_er_2(kisf,kWOA-1)
        NN(kWOA)=1
      endif
    enddo
  enddo

  ! within 10deg from the ice shelf front but only offshore of the continental shelf (sea floor > 1500m)
  NN(:) = 0.e0
  do iWOA=1,mlon
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
    if (       tmplon(iWOA) .ge. front_inf  &
    &    .and. tmplon(iWOA) .le. front_sup  ) then
      do jWOA=1,j50S
         if (       lat(jWOA) .ge. front_min_lat(kisf) - lat_box_3  &
         &    .and. lat(jWOA) .le. front_max_lat(kisf) + lat_box_3  ) then
           do kWOA=1,k30 
             if (       theta_an(iWOA,jWOA, 1  ) .ge.  -5.0 .and. theta_an(iWOA,jWOA, 1  ) .lt.  33.0 &  ! to remove ice shelf cavities
             &    .and. theta_an(iWOA,jWOA, k15) .ge.  -5.0 .and. theta_an(iWOA,jWOA, k15) .lt.  33.0 &  ! to only keep pts offshore of the continental shelf
             &    .and. theta_an(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_an(iWOA,jWOA,kWOA) .lt.  33.0 &
             &    .and. theta_se(iWOA,jWOA,kWOA) .ge.  -5.0 .and. theta_se(iWOA,jWOA,kWOA) .lt.  33.0 &
             &    .and. s_an(iWOA,jWOA,kWOA) .ge.   0.0 .and. s_an(iWOA,jWOA,kWOA) .lt.  60.0 &
             &    .and. s_se(iWOA,jWOA,kWOA) .ge.   0.0 .and. s_se(iWOA,jWOA,kWOA) .lt.  60.0 ) then
               mean_t_an_3(kisf,kWOA) = mean_t_an_3(kisf,kWOA) + theta_an(iWOA,jWOA,kWOA)
               mean_t_er_3(kisf,kWOA) = mean_t_er_3(kisf,kWOA) + theta_se(iWOA,jWOA,kWOA)
               mean_s_an_3(kisf,kWOA) = mean_s_an_3(kisf,kWOA) + s_an(iWOA,jWOA,kWOA)
               mean_s_er_3(kisf,kWOA) = mean_s_er_3(kisf,kWOA) + s_se(iWOA,jWOA,kWOA)
               NN(kWOA) = NN(kWOA) + 1
             endif
           enddo
         endif
      enddo
    endif
  enddo
  do kWOA=1,k30
    if ( NN(kWOA) .ne. 0 ) then
      mean_t_an_3(kisf,kWOA) = mean_t_an_3(kisf,kWOA) / NN(kWOA)
      mean_t_er_3(kisf,kWOA) = mean_t_er_3(kisf,kWOA) / NN(kWOA)
      mean_s_an_3(kisf,kWOA) = mean_s_an_3(kisf,kWOA) / NN(kWOA)
      mean_s_er_3(kisf,kWOA) = mean_s_er_3(kisf,kWOA) / NN(kWOA)
    else
      mean_t_an_3(kisf,kWOA) = NF90_FILL_FLOAT
      mean_t_er_3(kisf,kWOA) = NF90_FILL_FLOAT
      mean_s_an_3(kisf,kWOA) = NF90_FILL_FLOAT
      mean_s_er_3(kisf,kWOA) = NF90_FILL_FLOAT
    endif
  enddo
  ! filling gaps :
  do ks=1,5
  do kr=ks,1,-1
    do kWOA=2,k30-kr
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA+kr) .gt. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_3(kisf,kWOA) = ( mean_t_an_3(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_an_3(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_t_er_3(kisf,kWOA) = ( mean_t_er_3(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_t_er_3(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_s_an_3(kisf,kWOA) = ( mean_s_an_3(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_an_3(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        mean_s_er_3(kisf,kWOA) = ( mean_s_er_3(kisf,kWOA+kr) * (depth(kWOA)-depth(kWOA-1)) + mean_s_er_3(kisf,kWOA-1) * (depth(kWOA+kr)-depth(kWOA)) ) / (depth(kWOA+kr)-depth(kWOA-1))
        NN(kWOA)=1
      endif
    enddo
  enddo
  enddo
  !- extrapolate downward only if there exists a non-zero value below 100m
  do kr=1,k30-k01+1
    do kWOA=k01+1,k30
      if ( NN(kWOA) .eq. 0 .and. NN(kWOA-1) .gt. 0 ) then 
        mean_t_an_3(kisf,kWOA) = mean_t_an_3(kisf,kWOA-1)
        mean_t_er_3(kisf,kWOA) = mean_t_er_3(kisf,kWOA-1)
        mean_s_an_3(kisf,kWOA) = mean_s_an_3(kisf,kWOA-1)
        mean_s_er_3(kisf,kWOA) = mean_s_er_3(kisf,kWOA-1)
        NN(kWOA)=1
      endif
    enddo
  enddo

 ELSE

  mean_t_an_1(kisf,:) = NF90_FILL_FLOAT
  mean_t_er_1(kisf,:) = NF90_FILL_FLOAT
  mean_s_an_1(kisf,:) = NF90_FILL_FLOAT
  mean_s_er_1(kisf,:) = NF90_FILL_FLOAT
  mean_t_an_2(kisf,:) = NF90_FILL_FLOAT
  mean_t_er_2(kisf,:) = NF90_FILL_FLOAT
  mean_s_an_2(kisf,:) = NF90_FILL_FLOAT
  mean_s_er_2(kisf,:) = NF90_FILL_FLOAT
  mean_t_an_3(kisf,:) = NF90_FILL_FLOAT 
  mean_t_er_3(kisf,:) = NF90_FILL_FLOAT 
  mean_s_an_3(kisf,:) = NF90_FILL_FLOAT 
  mean_s_er_3(kisf,:) = NF90_FILL_FLOAT 

 ENDIF

ENDDO


!---------------------------------------
! Writing new netcdf file :                                   

write(*,*) 'Creating ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create output file')

status = NF90_DEF_DIM(fidM,"Nisf",Nisf,dimID_Nisf); call erreur(status,.TRUE.,"def_dimID_Nisf")
status = NF90_DEF_DIM(fidM,"depth",k30,dimID_depth); call erreur(status,.TRUE.,"def_dimID_depth")
status = NF90_DEF_DIM(fidM,"nbounds",mnbounds,dimID_nbounds); call erreur(status,.TRUE.,"def_dimID_nbounds")

status = NF90_DEF_VAR(fidM,"depth_bnds",NF90_FLOAT,(/dimID_nbounds,dimID_depth/),depth_bnds_ID); call erreur(status,.TRUE.,"def_var_depth_bnds_ID")
status = NF90_DEF_VAR(fidM,"depth",NF90_FLOAT,(/dimID_depth/),depth_ID); call erreur(status,.TRUE.,"def_var_depth_ID")
status = NF90_DEF_VAR(fidM,"t_an_1",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_t_an_1_ID); call erreur(status,.TRUE.,"def_var_mean_t_an_1_ID")
status = NF90_DEF_VAR(fidM,"t_er_1",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_t_er_1_ID); call erreur(status,.TRUE.,"def_var_mean_t_er_1_ID")
status = NF90_DEF_VAR(fidM,"s_an_1",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_s_an_1_ID); call erreur(status,.TRUE.,"def_var_mean_s_an_1_ID")
status = NF90_DEF_VAR(fidM,"s_er_1",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_s_er_1_ID); call erreur(status,.TRUE.,"def_var_mean_s_er_1_ID")
status = NF90_DEF_VAR(fidM,"t_an_2",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_t_an_2_ID); call erreur(status,.TRUE.,"def_var_mean_t_an_2_ID")
status = NF90_DEF_VAR(fidM,"t_er_2",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_t_er_2_ID); call erreur(status,.TRUE.,"def_var_mean_t_er_2_ID")
status = NF90_DEF_VAR(fidM,"s_an_2",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_s_an_2_ID); call erreur(status,.TRUE.,"def_var_mean_s_an_2_ID")
status = NF90_DEF_VAR(fidM,"s_er_2",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_s_er_2_ID); call erreur(status,.TRUE.,"def_var_mean_s_er_2_ID")
status = NF90_DEF_VAR(fidM,"t_an_3",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_t_an_3_ID); call erreur(status,.TRUE.,"def_var_mean_t_an_3_ID")
status = NF90_DEF_VAR(fidM,"t_er_3",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_t_er_3_ID); call erreur(status,.TRUE.,"def_var_mean_t_er_3_ID")
status = NF90_DEF_VAR(fidM,"s_an_3",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_s_an_3_ID); call erreur(status,.TRUE.,"def_var_mean_s_an_3_ID")
status = NF90_DEF_VAR(fidM,"s_er_3",NF90_FLOAT,(/dimID_Nisf,dimID_depth/),mean_s_er_3_ID); call erreur(status,.TRUE.,"def_var_mean_s_er_3_ID")

status = NF90_PUT_ATT(fidM,depth_bnds_ID,"comment","depth bounds"); call erreur(status,.TRUE.,"put_att_depth_bnds_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"axis","Z"); call erreur(status,.TRUE.,"put_att_depth_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"units","meters"); call erreur(status,.TRUE.,"put_att_depth_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"positive","down"); call erreur(status,.TRUE.,"put_att_depth_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"bounds","depth_bnds"); call erreur(status,.TRUE.,"put_att_depth_ID")
status = NF90_PUT_ATT(fidM,depth_ID,"standard_name","depth"); call erreur(status,.TRUE.,"put_att_depth_ID")

status = NF90_PUT_ATT(fidM,mean_t_an_1_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_1_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_1_ID,"long_name","Potential temperature in box_1 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_1_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_1_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_1_ID,"long_name","mean temperature stddev in box_1"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_1_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_1_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_1_ID,"long_name","Practical Salinity in box_1 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_1_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_1_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_1_ID,"long_name","mean salinity stddev in box_1"); call erreur(status,.TRUE.,"put_att_t_an_ID")

status = NF90_PUT_ATT(fidM,mean_t_an_2_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_2_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_2_ID,"long_name","Potential temperature in box_2 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_2_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_2_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_2_ID,"long_name","mean temperature stddev in box_2"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_2_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_2_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_2_ID,"long_name","Practical Salinity in box_2 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_2_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_2_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_2_ID,"long_name","mean salinity stddev in box_2"); call erreur(status,.TRUE.,"put_att_t_an_ID")

status = NF90_PUT_ATT(fidM,mean_t_an_3_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_3_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_3_ID,"long_name","Potential temperature in box_3 seaward of the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_3_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_3_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_er_3_ID,"long_name","mean temperature stddev in box_3"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_3_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_3_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_3_ID,"long_name","Practical Salinity in box_3 seaward of the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_3_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_3_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_er_3_ID,"long_name","mean salinity stddev in box_3"); call erreur(status,.TRUE.,"put_att_t_an_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using calculate_TS_profiles_WOA13.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lon_box_1",lon_box_1); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lat_box_1",lat_box_1); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lon_box_2",lon_box_2); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lat_box_2",lat_box_2); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lon_box_3",lon_box_3); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lat_box_3",lat_box_3); call erreur(status,.TRUE.,"put_att_GLOBAL")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 

status = NF90_PUT_VAR(fidM,depth_bnds_ID,depth_bnds(:,1:k30)); call erreur(status,.TRUE.,"var_depth_bnds_ID")
status = NF90_PUT_VAR(fidM,depth_ID,depth(1:k30));             call erreur(status,.TRUE.,"var_depth_ID")
status = NF90_PUT_VAR(fidM,mean_t_an_1_ID,mean_t_an_1); call erreur(status,.TRUE.,"var_mean_t_an_1_ID")
status = NF90_PUT_VAR(fidM,mean_t_er_1_ID,mean_t_er_1); call erreur(status,.TRUE.,"var_mean_t_er_1_ID")
status = NF90_PUT_VAR(fidM,mean_s_an_1_ID,mean_s_an_1); call erreur(status,.TRUE.,"var_mean_s_an_1_ID")
status = NF90_PUT_VAR(fidM,mean_s_er_1_ID,mean_s_er_1); call erreur(status,.TRUE.,"var_mean_s_er_1_ID")
status = NF90_PUT_VAR(fidM,mean_t_an_2_ID,mean_t_an_2); call erreur(status,.TRUE.,"var_mean_t_an_2_ID")
status = NF90_PUT_VAR(fidM,mean_t_er_2_ID,mean_t_er_2); call erreur(status,.TRUE.,"var_mean_t_er_2_ID")
status = NF90_PUT_VAR(fidM,mean_s_an_2_ID,mean_s_an_2); call erreur(status,.TRUE.,"var_mean_s_an_2_ID")
status = NF90_PUT_VAR(fidM,mean_s_er_2_ID,mean_s_er_2); call erreur(status,.TRUE.,"var_mean_s_er_2_ID")
status = NF90_PUT_VAR(fidM,mean_t_an_3_ID,mean_t_an_3); call erreur(status,.TRUE.,"var_mean_t_an_3_ID")
status = NF90_PUT_VAR(fidM,mean_t_er_3_ID,mean_t_er_3); call erreur(status,.TRUE.,"var_mean_t_er_3_ID")
status = NF90_PUT_VAR(fidM,mean_s_an_3_ID,mean_s_an_3); call erreur(status,.TRUE.,"var_mean_s_an_3_ID")
status = NF90_PUT_VAR(fidM,mean_s_er_3_ID,mean_s_er_3); call erreur(status,.TRUE.,"var_mean_s_er_3_ID")

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
