program modif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. Jourdain, IGE-CNRS, MAY-2019
!
! Put RTOPO's mon/lat grid to a stereographic projection.
!
! We use BedMachine-Antarctica's grid as a target.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE netcdf

IMPLICIT NONE

INTEGER :: fidISF, fidST, fidM, status, dimID_lat, dimID_lon, dimID_y, dimID_x, mlat, mlon, mx_new, my_new, isfmask_ID, area_isf_ID, err_melt_isf_ID, melt_isf_ID, lat_ID, lon_ID, x_ID, y_ID, rs, rx, ry, ki, kj, ineig, jneig, kinr, kjnr, iR, jR, pm
 
INTEGER*4,ALLOCATABLE,DIMENSION(:) :: y_new, x_new

INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: NN

INTEGER*2,ALLOCATABLE,DIMENSION(:,:) :: isfmask, isfmask_new, tmp_isfmask_new

REAL*4,ALLOCATABLE,DIMENSION(:) :: lat, lon

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: lat2d, lon2d, x_RTOPO, y_RTOPO

REAL*4 :: a, e, lat_c, t_c, m_c, rho, lon_0, dx, dy, deg2rad, pi

LOGICAL :: ll_check
 
CHARACTER(LEN=150) :: file_isf, file_stgrid, file_out

!-------------------------------------------------------------------------------

file_isf  = 'RTopo-2.0.1_30sec_ice_shelf_mask_ANT.nc'  ! NB: ANT = extended grid compared to ANT3
 
file_stgrid  = '/store/njourd/DATA/DATA_BATHYMETRY/BedMachineAntarctica-2018-12-22.nc'

file_out = 'RTopo-2.0.1_30sec_ice_shelf_mask_STEREOGRAPHIC.nc'

pi = dacos(-1.d0)

deg2rad = dacos(-1.d0) / 180.d0 

!-------------------------------------------------------------------------------
! Read ice shelves mask and name :
 
write(*,*) 'Reading ', TRIM(file_isf)
 
status = NF90_OPEN(TRIM(file_isf),0,fidISF); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidISF,"lat",dimID_lat); call erreur(status,.TRUE.,"inq_dimID_lat")
status = NF90_INQ_DIMID(fidISF,"lon",dimID_lon); call erreur(status,.TRUE.,"inq_dimID_lon")
 
status = NF90_INQUIRE_DIMENSION(fidISF,dimID_lat,len=mlat); call erreur(status,.TRUE.,"inq_dim_lat")
status = NF90_INQUIRE_DIMENSION(fidISF,dimID_lon,len=mlon); call erreur(status,.TRUE.,"inq_dim_lon")
  
ALLOCATE(  isfmask(mlon,mlat)  ) 
ALLOCATE(  lat(mlat)  ) 
ALLOCATE(  lon(mlon)  ) 
ALLOCATE(  lat2d(mlon,mlat)  )
ALLOCATE(  lon2d(mlon,mlat)  )
 
status = NF90_INQ_VARID(fidISF,"isfmask",isfmask_ID); call erreur(status,.TRUE.,"inq_isfmask_ID")
status = NF90_INQ_VARID(fidISF,"lat",lat_ID); call erreur(status,.TRUE.,"inq_lat_ID")
status = NF90_INQ_VARID(fidISF,"lon",lon_ID); call erreur(status,.TRUE.,"inq_lon_ID")
 
status = NF90_GET_VAR(fidISF,isfmask_ID,isfmask); call erreur(status,.TRUE.,"getvar_isfmask")
status = NF90_GET_VAR(fidISF,lat_ID,lat); call erreur(status,.TRUE.,"getvar_lat")
status = NF90_GET_VAR(fidISF,lon_ID,lon); call erreur(status,.TRUE.,"getvar_lon")
 
status = NF90_CLOSE(fidISF); call erreur(status,.TRUE.,"close_file")

do iR=1,mlon
  lat2d(iR,:)=lat(:)
enddo
do jR=1,mlat
  lon2d(:,jR)=lon(:)
enddo

DEALLOCATE(lon,lat)

!-------------------------------------------------------------------------------
! Read stereographic grid :
 
write(*,*) 'Reading ', TRIM(file_stgrid)
 
status = NF90_OPEN(TRIM(file_stgrid),0,fidST); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidST,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidST,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
 
status = NF90_INQUIRE_DIMENSION(fidST,dimID_y,len=my_new); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidST,dimID_x,len=mx_new); call erreur(status,.TRUE.,"inq_dim_x")
  
ALLOCATE(  y_new(my_new)  ) 
ALLOCATE(  x_new(mx_new)  ) 
 
status = NF90_INQ_VARID(fidST,"y",y_ID); call erreur(status,.TRUE.,"inq_y_ID")
status = NF90_INQ_VARID(fidST,"x",x_ID); call erreur(status,.TRUE.,"inq_x_ID")
 
status = NF90_GET_VAR(fidST,y_ID,y_new); call erreur(status,.TRUE.,"getvar_y")
status = NF90_GET_VAR(fidST,x_ID,x_new); call erreur(status,.TRUE.,"getvar_x")
 
status = NF90_CLOSE(fidST); call erreur(status,.TRUE.,"close_file")

!-------------------------------------------------------------------------------
!Conversion to stereographic coordinates:

write(*,*) 'Conversion to stereographic coordinates'

ALLOCATE( x_RTOPO(mlon,mlat) )
ALLOCATE( y_RTOPO(mlon,mlat) )

!- WGS84 :
a     = 6378137.0         ! Earth radius (m)
e     =       0.08181919  ! excentricity
lat_c =     -71.0         ! true latitude (deg)
lon_0 =       0.0

!- convert to radians :
lat2d(:,:) = deg2rad * lat2d(:,:)
lat_c      = deg2rad * lat_c
lon2d(:,:) = deg2rad * lon2d(:,:)
lon_0      = deg2rad * lon_0

!- if the standard parallel is in S.Hemi., switch signs.
if ( lat_c .lt. 0.0 ) then 
    pm         = -1    ! plus or minus, north lat. or south
    lat2d(:,:) = -lat2d(:,:)
    lat_c      = -lat_c
    lon2d(:,:) = -lon2d(:,:)
    lon_0      = -lon_0
else
    pm         = 1
endif

t_c = tan(pi/4-lat_c     /2) / ( (1-e*sin(lat_c     )) / (1+e*sin(lat_c     )) )**(e/2)
m_c = cos(lat_c) / sqrt( 1 - e**2 * (sin(lat_c))**2 )

do iR=1,mlon
do jR=1,mlat
  rho = a * m_c * tan(pi/4-lat2d(iR,jR)/2) / ( (1-e*sin(lat2d(iR,jR))) / (1+e*sin(lat2d(iR,jR))) )**(e/2) / t_c
  x_RTOPO(iR,jR) =  pm * rho * sin(lon2d(iR,jR)-lon_0)
  y_RTOPO(iR,jR) = -pm * rho * cos(lon2d(iR,jR)-lon_0)
enddo
enddo

write(*,*) 'min/max x_RTOPO = ', MINVAL(x_RTOPO(:,1:mlat)), MAXVAL(x_RTOPO(:,1:mlat))
write(*,*) 'min/max y_RTOPO = ', MINVAL(y_RTOPO(:,1:mlat)), MAXVAL(y_RTOPO(:,1:mlat))

!-----

write(*,*) 'Starting interpolation from RTOPO'

ALLOCATE( isfmask_new(mx_new,my_new) )
ALLOCATE( NN(mx_new,my_new) )

isfmask_new(:,:) = 0
NN(:,:) = 0

dx = FLOAT(ABS(x_new(2)-x_new(1)))
dy = FLOAT(ABS(y_new(2)-y_new(1)))

do iR=1,mlon
do jR=1,mlat

   !- nearset neighbour :
   kinr = MINLOC( abs( x_RTOPO(iR,jR)-x_new(:) ), 1 )
   kjnr = MINLOC( abs( y_RTOPO(iR,jR)-y_new(:) ), 1 )

   !- here, overwrite if several values in the save x,y mesh (no average because isfmask is discrete)
   if (       x_RTOPO(iR,jR) .ge. x_new(kinr)-0.5*dx .and. x_RTOPO(iR,jR) .le. x_new(kinr)+0.5*dx & 
   &    .and. y_RTOPO(iR,jR) .ge. y_new(kjnr)-0.5*dy .and. y_RTOPO(iR,jR) .le. y_new(kjnr)+0.5*dy ) then
        isfmask_new(kinr,kjnr) = isfmask(iR,jR) 
        NN(kinr,kjnr) = NN(kinr,kjnr) + 1
   endif

enddo
enddo

DEALLOCATE( isfmask )

!-----

write(*,*) 'Completing interpolation from RTOPO'

ALLOCATE( tmp_isfmask_new(mx_new,my_new) )
tmp_isfmask_new(:,:) = isfmask_new(:,:)

do ki=1,mx_new
do kj=1,my_new
 !- nearest-neighbour interpolation
 if ( NN(ki,kj) .eq. 0 ) then
   ineig=0
   jneig=0
   ll_check = .false.
   do rs=1,20
     do rx=0,rs
     do ry=0,rs
       if ( NN(MIN(mx_new,ki+rx),MIN(my_new,kj+ry)) .ge. 1 ) then; ineig=MIN(mx_new,ki+rx); jneig=MIN(my_new,kj+ry); ll_check=.true.; exit; endif
       if ( NN(MAX(    1 ,ki-rx),MIN(my_new,kj+ry)) .ge. 1 ) then; ineig=MAX(    1 ,ki-rx); jneig=MIN(my_new,kj+ry); ll_check=.true.; exit; endif
       if ( NN(MIN(mx_new,ki+rx),MAX(    1 ,kj-ry)) .ge. 1 ) then; ineig=MIN(mx_new,ki+rx); jneig=MAX(    1 ,kj-ry); ll_check=.true.; exit; endif
       if ( NN(MAX(    1 ,ki-rx),MAX(    1 ,kj-ry)) .ge. 1 ) then; ineig=MAX(    1 ,ki-rx); jneig=MAX(    1 ,kj-ry); ll_check=.true.; exit; endif
     enddo
     if ( ll_check ) exit
     enddo
     if ( ll_check ) exit
   enddo
   if ( ll_check ) then
     tmp_isfmask_new(ki,kj) = isfmask_new(ineig,jneig)
   else
     ! area not covered by our extraction of RTOPO2:
     if ( ( FLOAT(x_new(ki))**2 + FLOAT(y_new(kj))**2 ) .lt. 1.e12 ) then
       tmp_isfmask_new(ki,kj) = 0 ! center of ice sheet (grounded)
     else
       tmp_isfmask_new(ki,kj) = 1 ! far offshore
     endif
   endif
 endif
 ! Manual correction to remove grouned iceberg (?) in the Weddell Sea :
 if ( FLOAT(x_new(ki)) .gt. -0.9e6 .and. FLOAT(x_new(ki)) .lt. -1.1e6 .and. &
 &    FLOAT(y_new(kj)) .gt.  1.1e6 .and. FLOAT(y_new(kj)) .lt.  1.3e6       ) then
   tmp_isfmask_new(ki,kj) = 1
   NN(ki,kj) = 1
 endif
enddo
enddo

isfmask_new(:,:) = tmp_isfmask_new(:,:)

DEALLOCATE( tmp_isfmask_new )

!-------------------------------------------------------------------------------
! Writing new netcdf file :
 
write(*,*) 'Creating ', TRIM(file_out)
 
status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"y",my_new,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx_new,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")
  
status = NF90_DEF_VAR(fidM,"isfmask",NF90_FLOAT,(/dimID_x,dimID_y/),isfmask_ID); call erreur(status,.TRUE.,"def_var_isfmask_ID")
status = NF90_DEF_VAR(fidM,"y",NF90_INT,(/dimID_y/),y_ID); call erreur(status,.TRUE.,"def_var_y_ID")
status = NF90_DEF_VAR(fidM,"x",NF90_INT,(/dimID_x/),x_ID); call erreur(status,.TRUE.,"def_var_x_ID")

status = NF90_PUT_ATT(fidM,isfmask_ID,"title","isfmask"); call erreur(status,.TRUE.,"put_att_isfmask_ID")
status = NF90_PUT_ATT(fidM,isfmask_ID,"long_name","ice shelf mask (0 for grounded, 1 for ocean, isf ID for ice shelves)"); call erreur(status,.TRUE.,"put_att_isfmask_ID")
status = NF90_PUT_ATT(fidM,isfmask_ID,"units","-"); call erreur(status,.TRUE.,"put_att_isfmask_ID")
status = NF90_PUT_ATT(fidM,y_ID,"units","meter"); call erreur(status,.TRUE.,"put_att_y_ID")
status = NF90_PUT_ATT(fidM,y_ID,"standard_name","projection_y_coordinate"); call erreur(status,.TRUE.,"put_att_y_ID")
status = NF90_PUT_ATT(fidM,y_ID,"long_name","Cartesian y-coordinate"); call erreur(status,.TRUE.,"put_att_y_ID")
status = NF90_PUT_ATT(fidM,x_ID,"units","meter"); call erreur(status,.TRUE.,"put_att_x_ID")
status = NF90_PUT_ATT(fidM,x_ID,"standard_name","projection_x_coordinate"); call erreur(status,.TRUE.,"put_att_x_ID")
status = NF90_PUT_ATT(fidM,x_ID,"long_name","Cartesian x-coordinate"); call erreur(status,.TRUE.,"put_att_x_ID")
 
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using put_isfmask_on_stereographic_grid.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
 
status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidM,isfmask_ID,isfmask_new); call erreur(status,.TRUE.,"var_isfmask_ID")
status = NF90_PUT_VAR(fidM,y_ID,y_new); call erreur(status,.TRUE.,"var_y_ID")
status = NF90_PUT_VAR(fidM,x_ID,x_new); call erreur(status,.TRUE.,"var_x_ID")

status = NF90_CLOSE(fidM); call erreur(status,.TRUE.,"final")

end program modif

!======================================================

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
