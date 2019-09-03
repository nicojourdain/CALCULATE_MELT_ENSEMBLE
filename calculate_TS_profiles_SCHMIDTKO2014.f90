program modif                                         

USE netcdf                                            

IMPLICIT NONE                                         

INTEGER :: fidT, dimID_depth, dimID_lon, dimID_lat, dimID_nbounds, mdepth, mlon, mlat, mnbounds, t_se_ID, theta_an_ID, s_se_ID, s_an_ID, &
&          depth_bnds_ID, depth_ID, lon_bnds_ID, lon_ID, lat_bnds_ID, lat_ID, fidS, k15, k30, iSHMID, jSHMID, kWOA, kisf, kr, ks, k01,   &
&          mean_s_an_4_ID, mean_t_an_4_ID, mean_s_an_5_ID, mean_t_an_5_ID, fidISF, status, dimID_Nisf, front_max_lat_ID, melt_isf_ID,    &
&          front_min_lat_ID, front_max_lon_ID, front_min_lon_ID, front_ice_dep_avg_ID, front_ice_dep_min_ID, front_bot_dep_avg_ID, Nisf, &
&          front_bot_dep_max_ID, ii, jj, kk, NN, fidM

CHARACTER(LEN=150) :: file_isf, file_out, file_in_T, file_in_S

REAL*8,ALLOCATABLE,DIMENSION(:) :: front_max_lat, front_min_lat, front_max_lon, front_min_lon, front_ice_dep_avg, front_ice_dep_min, &
&                                  front_bot_dep_avg, front_bot_dep_max, depth, lon, lat, tmplon, melt_isf
 
REAL*8,ALLOCATABLE,DIMENSION(:) :: mean_t_an_4, mean_s_an_4, mean_t_an_5, mean_s_an_5
 
REAL*8,ALLOCATABLE,DIMENSION(:,:) :: s_an, theta_an

REAL*8 :: dif, mindif, front_inf, front_sup, lon_box_1, lon_box_2, lat_box_1, lat_box_2

REAL*8 :: Sabs
 
file_in_T  = '/store/njourd/DATA/SCHMIDTKO_2014/ASBW_Schmidtko.nc'
file_in_S  = '/store/njourd/DATA/SCHMIDTKO_2014/ASBW_Schmidtko.nc'
file_isf   = 'RTopo-2.0.1_30sec_ice_shelf_mask_ANT3.nc'
file_out   = 'TS_per_ice_shelf_ANT3_SCHMIDTKO2014.nc'

! Boxes :
! NB: 5.0 x 1.75 is the effective resolution at 70°S for a model of 1° resolution in longitude (assuming 5 delta X and a Mercator grid)
lon_box_1 =  5.0 ; lat_box_1 = 1.75 ! Continental shelf only
lon_box_2 = 10.0 ; lat_box_2 = 3.50 ! Continental shelf only

!---------------------------------------
! Read sea floor potential temperature :

write(*,*) 'Reading ', TRIM(file_in_T)

status = NF90_OPEN(TRIM(file_in_T),0,fidT); call erreur(status,.TRUE.,"read WOA temperature")

status = NF90_INQ_DIMID(fidT,"lon",dimID_lon); call erreur(status,.TRUE.,"inq_dimID_lon")
status = NF90_INQ_DIMID(fidT,"lat",dimID_lat); call erreur(status,.TRUE.,"inq_dimID_lat")

status = NF90_INQUIRE_DIMENSION(fidT,dimID_lon,len=mlon); call erreur(status,.TRUE.,"inq_dim_lon")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_lat,len=mlat); call erreur(status,.TRUE.,"inq_dim_lat")

ALLOCATE(  theta_an(mlon,mlat)  ) 
ALLOCATE(  lon(mlon)  ) 
ALLOCATE(  tmplon(mlon)  ) 
ALLOCATE(  lat(mlat)  ) 

status = NF90_INQ_VARID(fidT,"PT",theta_an_ID); call erreur(status,.TRUE.,"inq_t_an_ID")
status = NF90_INQ_VARID(fidT,"lon",lon_ID); call erreur(status,.TRUE.,"inq_lon_ID")
status = NF90_INQ_VARID(fidT,"lat",lat_ID); call erreur(status,.TRUE.,"inq_lat_ID")
 
status = NF90_GET_VAR(fidT,theta_an_ID,theta_an); call erreur(status,.TRUE.,"getvar_t_an")
status = NF90_GET_VAR(fidT,lon_ID,lon); call erreur(status,.TRUE.,"getvar_lon")
status = NF90_GET_VAR(fidT,lat_ID,lat); call erreur(status,.TRUE.,"getvar_lat")
 
status = NF90_CLOSE(fidT); call erreur(status,.TRUE.,"close_file")

write(*,*) '  Minimum longitude for Schmidtko et al. (2014) is : ', MINVAL(lon)

! put lon in [-180:180] :
do ii=1,mlon
  if ( lon(ii) .gt. 180.e0 ) then
    lon(ii) = lon(ii) - 360.e0
  endif
enddo

!---------------------------------------
! Read WOA2013 salinity :

write(*,*) 'Reading ', TRIM(file_in_S)

status = NF90_OPEN(TRIM(file_in_S),0,fidS); call erreur(status,.TRUE.,"read WOA temperature")

ALLOCATE(  s_an(mlon,mlat)  ) 

status = NF90_INQ_VARID(fidS,"PS",s_an_ID); call erreur(status,.TRUE.,"inq_s_an_ID")
 
status = NF90_GET_VAR(fidS,s_an_ID,s_an); call erreur(status,.TRUE.,"getvar_s_an")
 
status = NF90_CLOSE(fidS); call erreur(status,.TRUE.,"close_file")

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
ALLOCATE(  melt_isf(Nisf)  )

status = NF90_INQ_VARID(fidISF,"front_max_lat",front_max_lat_ID)
call erreur(status,.TRUE.,"inq_front_max_lat_ID")
status = NF90_INQ_VARID(fidISF,"front_min_lat",front_min_lat_ID)
call erreur(status,.TRUE.,"inq_front_min_lat_ID")
status = NF90_INQ_VARID(fidISF,"front_max_lon",front_max_lon_ID)
call erreur(status,.TRUE.,"inq_front_max_lon_ID")
status = NF90_INQ_VARID(fidISF,"front_min_lon",front_min_lon_ID)
call erreur(status,.TRUE.,"inq_front_min_lon_ID")
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
status = NF90_GET_VAR(fidISF,melt_isf_ID,melt_isf)
call erreur(status,.TRUE.,"getvar_melt_isf")

status = NF90_CLOSE(fidISF)                      
call erreur(status,.TRUE.,"fin_lecture")     

!------------------------------------------------------
! Calculate T,S profiles for individual ice shelves :

ALLOCATE( mean_t_an_4(Nisf), mean_s_an_4(Nisf) )
ALLOCATE( mean_t_an_5(Nisf), mean_s_an_5(Nisf) )

mean_t_an_4(:) = 0.e0
mean_s_an_4(:) = 0.e0
mean_t_an_5(:) = 0.e0
mean_s_an_5(:) = 0.e0

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

  !== First box ==
  NN = 0.e0
  do iSHMID=1,mlon
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
    if (       tmplon(iSHMID) .ge. front_inf  &
    &    .and. tmplon(iSHMID) .le. front_sup  ) then
      do jSHMID=1,mlat
         if (       lat(jSHMID) .ge. front_min_lat(kisf) - lat_box_1  &
         &    .and. lat(jSHMID) .le. front_max_lat(kisf) + lat_box_1  ) then
           if (   theta_an(iSHMID,jSHMID) .ge.  -5.0 .and. theta_an(iSHMID,jSHMID) .lt.  33.0 &
           &    .and. s_an(iSHMID,jSHMID) .ge.   0.0 .and. s_an(iSHMID,jSHMID) .lt.  60.0 ) then
             mean_t_an_4(kisf) = mean_t_an_4(kisf) + theta_an(iSHMID,jSHMID)
             mean_s_an_4(kisf) = mean_s_an_4(kisf) + s_an(iSHMID,jSHMID)
             NN = NN + 1
           endif
         endif
      enddo
    endif
  enddo
  if ( NN .ne. 0 ) then
    mean_t_an_4(kisf) = mean_t_an_4(kisf) / NN
    mean_s_an_4(kisf) = mean_s_an_4(kisf) / NN
  else
    mean_t_an_4(kisf) = NF90_FILL_DOUBLE
    mean_s_an_4(kisf) = NF90_FILL_DOUBLE
  endif

  !== Second box ==
  NN = 0.e0
  do iSHMID=1,mlon
    ! trick not to llok for T,S on the other side of the Peninsula :
    if ( front_max_lon(kisf) .ge. -68.0 .and. front_max_lon(kisf) .le. -66.0 ) then
       front_sup = -66.0
    else
       front_sup = front_max_lon(kisf) + lon_box_2
    endif
    if ( front_min_lon(kisf) .ge. -66.0 .and. front_min_lon(kisf) .le. -64.0 ) then
       front_inf = -66.0
    else
       front_inf = front_min_lon(kisf) - lon_box_2
    endif
    if (       tmplon(iSHMID) .ge. front_inf  &
    &    .and. tmplon(iSHMID) .le. front_sup  ) then
      do jSHMID=1,mlat
         if (       lat(jSHMID) .ge. front_min_lat(kisf) - lat_box_2  &
         &    .and. lat(jSHMID) .le. front_max_lat(kisf) + lat_box_2  ) then
           if (   theta_an(iSHMID,jSHMID) .ge.  -5.0 .and. theta_an(iSHMID,jSHMID) .lt.  33.0 &
           &    .and. s_an(iSHMID,jSHMID) .ge.   0.0 .and. s_an(iSHMID,jSHMID) .lt.  60.0 ) then
             mean_t_an_5(kisf) = mean_t_an_5(kisf) + theta_an(iSHMID,jSHMID)
             mean_s_an_5(kisf) = mean_s_an_5(kisf) + s_an(iSHMID,jSHMID)
             NN = NN + 1
           endif
         endif
      enddo
    endif
  enddo
  if ( NN .ne. 0 ) then
    mean_t_an_5(kisf) = mean_t_an_5(kisf) / NN
    mean_s_an_5(kisf) = mean_s_an_5(kisf) / NN
  else
    mean_t_an_5(kisf) = NF90_FILL_DOUBLE
    mean_s_an_5(kisf) = NF90_FILL_DOUBLE
  endif

 ELSE

  mean_t_an_4(kisf) = NF90_FILL_DOUBLE
  mean_s_an_4(kisf) = NF90_FILL_DOUBLE
  mean_t_an_5(kisf) = NF90_FILL_DOUBLE
  mean_s_an_5(kisf) = NF90_FILL_DOUBLE

 ENDIF

ENDDO


!---------------------------------------
! Writing new netcdf file :                                   

write(*,*) 'Creating ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create output file')

status = NF90_DEF_DIM(fidM,"Nisf",Nisf,dimID_Nisf); call erreur(status,.TRUE.,"def_dimID_Nisf")

status = NF90_DEF_VAR(fidM,"t_an_4",NF90_DOUBLE,(/dimID_Nisf/),mean_t_an_4_ID); call erreur(status,.TRUE.,"def_var_mean_t_an_4_ID")
status = NF90_DEF_VAR(fidM,"s_an_4",NF90_DOUBLE,(/dimID_Nisf/),mean_s_an_4_ID); call erreur(status,.TRUE.,"def_var_mean_s_an_4_ID")
status = NF90_DEF_VAR(fidM,"t_an_5",NF90_DOUBLE,(/dimID_Nisf/),mean_t_an_5_ID); call erreur(status,.TRUE.,"def_var_mean_t_an_5_ID")
status = NF90_DEF_VAR(fidM,"s_an_5",NF90_DOUBLE,(/dimID_Nisf/),mean_s_an_5_ID); call erreur(status,.TRUE.,"def_var_mean_s_an_5_ID")

status = NF90_PUT_ATT(fidM,mean_t_an_4_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_4_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_4_ID,"long_name","Potential sea floor temperature in box_1 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_4_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_4_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_4_ID,"long_name","Practical sea floor salinity in box_1 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")

status = NF90_PUT_ATT(fidM,mean_t_an_5_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_5_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_t_an_5_ID,"long_name","Potential sea floor temperature in box_2 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_5_ID,"_FillValue",NF90_FILL_DOUBLE); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_5_ID,"units","psu"); call erreur(status,.TRUE.,"put_att_t_an_ID")
status = NF90_PUT_ATT(fidM,mean_s_an_5_ID,"long_name","Practical sea floor salinity in box_2 over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_an_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using calculate_TS_profiles_SCHMIDTKO2014.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lon_box_1",lon_box_1); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lat_box_1",lat_box_1); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lon_box_2",lon_box_2); call erreur(status,.TRUE.,"put_att_GLOBAL")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lat_box_2",lat_box_2); call erreur(status,.TRUE.,"put_att_GLOBAL")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 

status = NF90_PUT_VAR(fidM,mean_t_an_4_ID,mean_t_an_4); call erreur(status,.TRUE.,"var_mean_t_an_4_ID")
status = NF90_PUT_VAR(fidM,mean_s_an_4_ID,mean_s_an_4); call erreur(status,.TRUE.,"var_mean_s_an_4_ID")
status = NF90_PUT_VAR(fidM,mean_t_an_5_ID,mean_t_an_5); call erreur(status,.TRUE.,"var_mean_t_an_5_ID")
status = NF90_PUT_VAR(fidM,mean_s_an_5_ID,mean_s_an_5); call erreur(status,.TRUE.,"var_mean_s_an_5_ID")

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
