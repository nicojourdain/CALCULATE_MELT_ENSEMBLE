program modif
 
USE netcdf
 
IMPLICIT NONE
 
INTEGER :: fidA, status, dimID_mpara, dimID_mstat, dimID_Nisf, dimID_lat, dimID_lon, mmpara, mmstat, mNisf,  &
&          mlat, mlon, index_tun_ID, index_WOA_ID, index_CMIP_ID, index_para_ID, Kcoef_ID, kstat2, &
&          mean_melt_futu_ID, total_melt_futu_ID, mean_melt_pres_ID, total_melt_pres_ID, Futu_Melt_ID,       &
&          Pres_Melt_ID, lat_ID, lon_ID, fidM, nn_tuning, ID_isf, kk_tuning, kk_para, imin, imax, jmin, jmax,&
&          ii, jj, kstat
 
CHARACTER(LEN=100) :: file_in, file_out
 
REAL*4,ALLOCATABLE,DIMENSION(:) :: lat, lon
 
INTEGER*1,ALLOCATABLE,DIMENSION(:) :: tmp_index_tun, tmp_index_WOA, tmp_index_CMIP, tmp_index_para

INTEGER*1,ALLOCATABLE,DIMENSION(:,:) :: index_tun, index_WOA, index_CMIP, index_para

INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:) :: Futu_NN, Pres_NN
 
REAL*4,ALLOCATABLE,DIMENSION(:) :: tmp_Kcoef, tmp_mean_melt_futu, tmp_total_melt_futu, tmp_mean_melt_pres, tmp_total_melt_pres

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: Kcoef, mean_melt_futu, total_melt_futu, mean_melt_pres, total_melt_pres
 
REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: Futu_Melt, Pres_Melt, tmp_Futu_Melt, tmp_Pres_Melt

LOGICAL :: ll_check
 
ID_isf = 11
nn_tuning = 19

101 FORMAT('melt_ensemble_',i3.3,'_',i1,'.nc')
102 FORMAT('melt_ensemble_',i3.3,'_',i2,'.nc')
103 FORMAT('melt_ensemble_',i3.3,'.nc')

kstat2 = 0

DO kk_tuning=1,nn_tuning

  if ( kk_tuning .lt. 10 ) then
    write(file_in,101) ID_isf, kk_tuning
  else
    write(file_in,102) ID_isf, kk_tuning
  endif
  write(file_out,103) ID_isf
   
  !---------------------------------------
  ! Read netcdf input file :
   
  write(*,*) 'Reading ', TRIM(file_in)
   
  status = NF90_OPEN(TRIM(file_in),0,fidA); call erreur(status,.TRUE.,"read")
   
  status = NF90_INQ_DIMID(fidA,"mpara",dimID_mpara); call erreur(status,.TRUE.,"inq_dimID_mpara")
  status = NF90_INQ_DIMID(fidA,"mstat",dimID_mstat); call erreur(status,.TRUE.,"inq_dimID_mstat")
  status = NF90_INQ_DIMID(fidA,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
  status = NF90_INQ_DIMID(fidA,"lat",dimID_lat); call erreur(status,.TRUE.,"inq_dimID_lat")
  status = NF90_INQ_DIMID(fidA,"lon",dimID_lon); call erreur(status,.TRUE.,"inq_dimID_lon")
   
  status = NF90_INQUIRE_DIMENSION(fidA,dimID_mpara,len=mmpara); call erreur(status,.TRUE.,"inq_dim_mpara")
  status = NF90_INQUIRE_DIMENSION(fidA,dimID_mstat,len=mmstat); call erreur(status,.TRUE.,"inq_dim_mstat")
  status = NF90_INQUIRE_DIMENSION(fidA,dimID_Nisf,len=mNisf); call erreur(status,.TRUE.,"inq_dim_Nisf")
  status = NF90_INQUIRE_DIMENSION(fidA,dimID_lat,len=mlat); call erreur(status,.TRUE.,"inq_dim_lat")
  status = NF90_INQUIRE_DIMENSION(fidA,dimID_lon,len=mlon); call erreur(status,.TRUE.,"inq_dim_lon")

  ALLOCATE(  tmp_index_tun(mmstat)  ) 
  ALLOCATE(  tmp_index_WOA(mmstat)  ) 
  ALLOCATE(  tmp_index_CMIP(mmstat)  ) 
  ALLOCATE(  tmp_index_para(mmstat)  ) 
  !ALLOCATE(  tmp_Ccoef(mmstat)  ) 
  ALLOCATE(  tmp_Kcoef(mmstat)  ) 
  ALLOCATE(  tmp_mean_melt_futu(mmstat)  ) 
  ALLOCATE(  tmp_total_melt_futu(mmstat)  ) 
  ALLOCATE(  tmp_mean_melt_pres(mmstat)  ) 
  ALLOCATE(  tmp_total_melt_pres(mmstat)  ) 
  ALLOCATE(  tmp_Futu_Melt(mlon,mlat,mmpara)  ) 
  ALLOCATE(  tmp_Pres_Melt(mlon,mlat,mmpara)  ) 
   
  status = NF90_INQ_VARID(fidA,"index_tun",index_tun_ID); call erreur(status,.TRUE.,"inq_index_tun_ID")
  status = NF90_INQ_VARID(fidA,"index_WOA",index_WOA_ID); call erreur(status,.TRUE.,"inq_index_WOA_ID")
  status = NF90_INQ_VARID(fidA,"index_CMIP",index_CMIP_ID); call erreur(status,.TRUE.,"inq_index_CMIP_ID")
  status = NF90_INQ_VARID(fidA,"index_para",index_para_ID); call erreur(status,.TRUE.,"inq_index_para_ID")
  !status = NF90_INQ_VARID(fidA,"Ccoef",Ccoef_ID); call erreur(status,.TRUE.,"inq_Ccoef_ID")
  status = NF90_INQ_VARID(fidA,"Kcoef",Kcoef_ID); call erreur(status,.TRUE.,"inq_Kcoef_ID")
  status = NF90_INQ_VARID(fidA,"mean_melt_futu",mean_melt_futu_ID); call erreur(status,.TRUE.,"inq_mean_melt_futu_ID")
  status = NF90_INQ_VARID(fidA,"total_melt_futu",total_melt_futu_ID); call erreur(status,.TRUE.,"inq_total_melt_futu_ID")
  status = NF90_INQ_VARID(fidA,"mean_melt_pres",mean_melt_pres_ID); call erreur(status,.TRUE.,"inq_mean_melt_pres_ID")
  status = NF90_INQ_VARID(fidA,"total_melt_pres",total_melt_pres_ID); call erreur(status,.TRUE.,"inq_total_melt_pres_ID")
  status = NF90_INQ_VARID(fidA,"Futu_Melt",Futu_Melt_ID); call erreur(status,.TRUE.,"inq_Futu_Melt_ID")
  status = NF90_INQ_VARID(fidA,"Pres_Melt",Pres_Melt_ID); call erreur(status,.TRUE.,"inq_Pres_Melt_ID")
   
  status = NF90_GET_VAR(fidA,index_tun_ID,tmp_index_tun,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_index_tun",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,index_WOA_ID,tmp_index_WOA,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_index_WOA",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,index_CMIP_ID,tmp_index_CMIP,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_index_CMIP",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,index_para_ID,tmp_index_para,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_index_para",(/ID_isf,1/),(/1,mmpara/))
  !status = NF90_GET_VAR(fidA,Ccoef_ID,tmp_Ccoef,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_Ccoef",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,Kcoef_ID,tmp_Kcoef,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_Kcoef",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,mean_melt_futu_ID,tmp_mean_melt_futu,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_mean_melt_futu",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,total_melt_futu_ID,tmp_total_melt_futu,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_total_melt_futu",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,mean_melt_pres_ID,tmp_mean_melt_pres,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_mean_melt_pres",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,total_melt_pres_ID,tmp_total_melt_pres,(/ID_isf,1/),(/1,mmpara/)); call erreur(status,.TRUE.,"getvar_total_melt_pres",(/ID_isf,1/),(/1,mmpara/))
  status = NF90_GET_VAR(fidA,Futu_Melt_ID,tmp_Futu_Melt); call erreur(status,.TRUE.,"getvar_Futu_Melt")
  status = NF90_GET_VAR(fidA,Pres_Melt_ID,tmp_Pres_Melt); call erreur(status,.TRUE.,"getvar_Pres_Melt")

  if ( kk_tuning .eq. 1 ) then
    ALLOCATE(  index_tun(mNisf,mmstat)  ) 
    ALLOCATE(  index_WOA(mNisf,mmstat)  ) 
    ALLOCATE(  index_CMIP(mNisf,mmstat)  ) 
    ALLOCATE(  index_para(mNisf,mmstat)  ) 
    !ALLOCATE(  Ccoef(mNisf,mmstat)  ) 
    ALLOCATE(  Kcoef(mNisf,mmstat)  ) 
    ALLOCATE(  mean_melt_futu(mNisf,mmstat)  ) 
    ALLOCATE(  total_melt_futu(mNisf,mmstat)  ) 
    ALLOCATE(  mean_melt_pres(mNisf,mmstat)  ) 
    ALLOCATE(  total_melt_pres(mNisf,mmstat)  ) 
    ALLOCATE(  Futu_Melt(mlon,mlat,mmpara)  ) 
    ALLOCATE(  Pres_Melt(mlon,mlat,mmpara)  )
    ALLOCATE(  Futu_NN  (mlon,mlat,mmpara)  ) 
    ALLOCATE(  Pres_NN  (mlon,mlat,mmpara)  )
    ! 
    index_tun       (:,:) = 0
    index_WOA       (:,:) = 0
    index_CMIP      (:,:) = 0
    index_para      (:,:) = 0
    !Ccoef           (:,:) = NF90_FILL_FLOAT
    Kcoef           (:,:) = NF90_FILL_FLOAT
    mean_melt_futu  (:,:) = NF90_FILL_FLOAT
    total_melt_futu (:,:) = NF90_FILL_FLOAT
    mean_melt_pres  (:,:) = NF90_FILL_FLOAT
    total_melt_pres (:,:) = NF90_FILL_FLOAT
    Futu_Melt(:,:,:) = 0.e0
    Pres_Melt(:,:,:) = 0.e0
    Futu_NN  (:,:,:) = 0
    Pres_NN  (:,:,:) = 0
    ! 
    ALLOCATE(  lat(mlat)  ) 
    ALLOCATE(  lon(mlon)  ) 
    status = NF90_INQ_VARID(fidA,"lat",lat_ID); call erreur(status,.TRUE.,"inq_lat_ID")
    status = NF90_INQ_VARID(fidA,"lon",lon_ID); call erreur(status,.TRUE.,"inq_lon_ID")
    status = NF90_GET_VAR(fidA,lat_ID,lat); call erreur(status,.TRUE.,"getvar_lat")
    status = NF90_GET_VAR(fidA,lon_ID,lon); call erreur(status,.TRUE.,"getvar_lon")
    status = NF90_GET_ATT(fidA,NF90_GLOBAL, "imin", imin)
    status = NF90_GET_ATT(fidA,NF90_GLOBAL, "imax", imax)
    status = NF90_GET_ATT(fidA,NF90_GLOBAL, "jmin", jmin)
    status = NF90_GET_ATT(fidA,NF90_GLOBAL, "jmax", jmax)
  endif
   
  status = NF90_CLOSE(fidA); call erreur(status,.TRUE.,"close_file")

  !---------------------------------------------------------------------------------
  ! Gathering data :

  ll_check = .false. 
  do kstat=1,mmstat
    if ( tmp_index_tun(kstat) .eq. kk_tuning ) then
      kstat2 = kstat2 + 1
      ll_check = .true.
      index_tun       (ID_isf,kstat2) = tmp_index_tun       (kstat)
      index_WOA       (ID_isf,kstat2) = tmp_index_WOA       (kstat)
      index_CMIP      (ID_isf,kstat2) = tmp_index_CMIP      (kstat)
      index_para      (ID_isf,kstat2) = tmp_index_para      (kstat)
      !Ccoef           (ID_isf,kstat2) = tmp_Ccoef           (kstat)
      Kcoef           (ID_isf,kstat2) = tmp_Kcoef           (kstat)
      mean_melt_futu  (ID_isf,kstat2) = tmp_mean_melt_futu  (kstat)
      total_melt_futu (ID_isf,kstat2) = tmp_total_melt_futu (kstat)
      mean_melt_pres  (ID_isf,kstat2) = tmp_mean_melt_pres  (kstat)
      total_melt_pres (ID_isf,kstat2) = tmp_total_melt_pres (kstat)
    endif
  enddo 

  if ( ll_check ) then
     do ii=1,mlon
     do jj=1,mlat
     do kk_para=1,mmpara
        if ( tmp_Futu_Melt(ii,jj,kk_para) .ne. NF90_FILL_FLOAT ) then
          Futu_Melt(ii,jj,kk_para) = Futu_Melt(ii,jj,kk_para) + tmp_Futu_Melt(ii,jj,kk_para)
          Futu_NN  (ii,jj,kk_para) = Futu_NN  (ii,jj,kk_para) + 1
        endif
        if ( tmp_Pres_Melt(ii,jj,kk_para) .ne. NF90_FILL_FLOAT ) then
          Pres_Melt(ii,jj,kk_para) = Pres_Melt(ii,jj,kk_para) + tmp_Pres_Melt(ii,jj,kk_para)
          Pres_NN  (ii,jj,kk_para) = Pres_NN  (ii,jj,kk_para) + 1
        endif
     enddo
     enddo
     enddo
  endif

  DEALLOCATE( tmp_index_tun, tmp_index_WOA, tmp_index_CMIP, tmp_index_para )
  DEALLOCATE( tmp_Kcoef, tmp_mean_melt_futu, tmp_total_melt_futu, tmp_mean_melt_pres, tmp_total_melt_pres )
  DEALLOCATE( tmp_Futu_Melt, tmp_Pres_Melt )

ENDDO

do ii=1,mlon
do jj=1,mlat
do kk_para=1,mmpara
   if ( Futu_NN(ii,jj,kk_para) .gt. 0 ) then
     Futu_Melt(ii,jj,kk_para) = Futu_Melt(ii,jj,kk_para) / Futu_NN(ii,jj,kk_para)
   else
     Futu_Melt(ii,jj,kk_para) = NF90_FILL_FLOAT
   endif
   if ( Pres_NN(ii,jj,kk_para) .gt. 0 ) then
     Pres_Melt(ii,jj,kk_para) = Pres_Melt(ii,jj,kk_para) / Pres_NN(ii,jj,kk_para)
   else
     Pres_Melt(ii,jj,kk_para) = NF90_FILL_FLOAT
   endif
enddo
enddo
enddo

!---------------------------------------------------------------------------------
! Writing new netcdf file :
 
write(*,*) 'Creating ', TRIM(file_out)
 
status = NF90_CREATE(TRIM(file_out),or(NF90_CLOBBER,NF90_64BIT_OFFSET),fidM); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"mpara",mmpara,dimID_mpara); call erreur(status,.TRUE.,"def_dimID_mpara")
status = NF90_DEF_DIM(fidM,"mstat",mmstat,dimID_mstat); call erreur(status,.TRUE.,"def_dimID_mstat")
status = NF90_DEF_DIM(fidM,"Nisf",mNisf,dimID_Nisf); call erreur(status,.TRUE.,"def_dimID_Nisf")
status = NF90_DEF_DIM(fidM,"lat",mlat,dimID_lat); call erreur(status,.TRUE.,"def_dimID_lat")
status = NF90_DEF_DIM(fidM,"lon",mlon,dimID_lon); call erreur(status,.TRUE.,"def_dimID_lon")
  
status = NF90_DEF_VAR(fidM,"index_tun",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_tun_ID); call erreur(status,.TRUE.,"def_var_index_tun_ID")
status = NF90_DEF_VAR(fidM,"index_WOA",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_WOA_ID); call erreur(status,.TRUE.,"def_var_index_WOA_ID")
status = NF90_DEF_VAR(fidM,"index_CMIP",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_CMIP_ID); call erreur(status,.TRUE.,"def_var_index_CMIP_ID")
status = NF90_DEF_VAR(fidM,"index_para",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_para_ID); call erreur(status,.TRUE.,"def_var_index_para_ID")
!status = NF90_DEF_VAR(fidM,"Ccoef",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),Ccoef_ID); call erreur(status,.TRUE.,"def_var_Ccoef_ID")
status = NF90_DEF_VAR(fidM,"Kcoef",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),Kcoef_ID); call erreur(status,.TRUE.,"def_var_Kcoef_ID")
status = NF90_DEF_VAR(fidM,"mean_melt_futu",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),mean_melt_futu_ID); call erreur(status,.TRUE.,"def_var_mean_melt_futu_ID")
status = NF90_DEF_VAR(fidM,"total_melt_futu",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),total_melt_futu_ID); call erreur(status,.TRUE.,"def_var_total_melt_futu_ID")
status = NF90_DEF_VAR(fidM,"mean_melt_pres",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),mean_melt_pres_ID); call erreur(status,.TRUE.,"def_var_mean_melt_pres_ID")
status = NF90_DEF_VAR(fidM,"total_melt_pres",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),total_melt_pres_ID); call erreur(status,.TRUE.,"def_var_total_melt_pres_ID")
status = NF90_DEF_VAR(fidM,"Futu_Melt",NF90_FLOAT,(/dimID_lon,dimID_lat,dimID_mpara/),Futu_Melt_ID); call erreur(status,.TRUE.,"def_var_Futu_Melt_ID")
status = NF90_DEF_VAR(fidM,"Pres_Melt",NF90_FLOAT,(/dimID_lon,dimID_lat,dimID_mpara/),Pres_Melt_ID); call erreur(status,.TRUE.,"def_var_Pres_Melt_ID")
status = NF90_DEF_VAR(fidM,"lat",NF90_FLOAT,(/dimID_lat/),lat_ID); call erreur(status,.TRUE.,"def_var_lat_ID")
status = NF90_DEF_VAR(fidM,"lon",NF90_FLOAT,(/dimID_lon/),lon_ID); call erreur(status,.TRUE.,"def_var_lon_ID")
 
status = NF90_PUT_ATT(fidM,index_tun_ID,"title","index_tun"); call erreur(status,.TRUE.,"put_att_index_tun_ID")
status = NF90_PUT_ATT(fidM,index_tun_ID,"long_name","index defining the sampling of Rignot (2013)'s range of uncertainty"); call erreur(status,.TRUE.,"put_att_index_tun_ID")
status = NF90_PUT_ATT(fidM,index_tun_ID,"units","-"); call erreur(status,.TRUE.,"put_att_index_tun_ID")
status = NF90_PUT_ATT(fidM,index_WOA_ID,"title","index_WOA"); call erreur(status,.TRUE.,"put_att_index_WOA_ID")
status = NF90_PUT_ATT(fidM,index_WOA_ID,"latg_name","index defining the sampling of WOA2013 range of uncertainty"); call erreur(status,.TRUE.,"put_att_index_WOA_ID")
status = NF90_PUT_ATT(fidM,index_WOA_ID,"units","-"); call erreur(status,.TRUE.,"put_att_index_WOA_ID")
status = NF90_PUT_ATT(fidM,index_CMIP_ID,"title","index_CMIP"); call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
status = NF90_PUT_ATT(fidM,index_CMIP_ID,"long_name","index defining the CMIP5 model used for the climate anomaly"); call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
status = NF90_PUT_ATT(fidM,index_CMIP_ID,"units","-"); call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
status = NF90_PUT_ATT(fidM,index_para_ID,"title","index_para"); call erreur(status,.TRUE.,"put_att_index_para_ID")
status = NF90_PUT_ATT(fidM,index_para_ID,"long_name","index defining the melt parameterization"); call erreur(status,.TRUE.,"put_att_index_para_ID")
status = NF90_PUT_ATT(fidM,index_para_ID,"units","-"); call erreur(status,.TRUE.,"put_att_index_para_ID")
!status = NF90_PUT_ATT(fidM,Ccoef_ID,"title","C coeff"); call erreur(status,.TRUE.,"put_att_Ccoef_ID")
!status = NF90_PUT_ATT(fidM,Ccoef_ID,"long_name","C coefficient (PICO's overturning)"); call erreur(status,.TRUE.,"put_att_Ccoef_ID")
!status = NF90_PUT_ATT(fidM,Ccoef_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_Ccoef_ID")
!status = NF90_PUT_ATT(fidM,Ccoef_ID,"units","TBD"); call erreur(status,.TRUE.,"put_att_Ccoef_ID")
status = NF90_PUT_ATT(fidM,Kcoef_ID,"title","K coeff"); call erreur(status,.TRUE.,"put_att_Kcoef_ID")
status = NF90_PUT_ATT(fidM,Kcoef_ID,"long_name","K coefficient (tuning factor)"); call erreur(status,.TRUE.,"put_att_Kcoef_ID")
status = NF90_PUT_ATT(fidM,Kcoef_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_Kcoef_ID")
status = NF90_PUT_ATT(fidM,Kcoef_ID,"units","TBD"); call erreur(status,.TRUE.,"put_att_Kcoef_ID")
status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"title","Future Melt Rate"); call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"long_name","Future mean melt rate over the cavity"); call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"units","m/yr"); call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"title","Future Melt Flux"); call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"long_name","Total future cavity melt flux"); call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"units","Gt/yr"); call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"title","Present Melt Rate"); call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"long_name","Present mean melt rate over the cavity"); call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"units","m/yr"); call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"title","Present Melt Flux"); call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"long_name","Total present cavity melt flux"); call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"units","Gt/yr"); call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"title","Future melt rate"); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"long_name","Future ice shelf melt rate"); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"units","m/yr"); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"coordinates","lat lon Nisf"); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"title","Present melt rate"); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"long_name","Present ice shelf melt rate"); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"units","m/yr"); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"coordinates","lat lon Nisf"); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,lat_ID,"name","lat"); call erreur(status,.TRUE.,"put_att_lat_ID")
status = NF90_PUT_ATT(fidM,lat_ID,"long_name","latitude coordinate"); call erreur(status,.TRUE.,"put_att_lat_ID")
status = NF90_PUT_ATT(fidM,lat_ID,"units","degrees_north"); call erreur(status,.TRUE.,"put_att_lat_ID")
status = NF90_PUT_ATT(fidM,lat_ID,"standard_name","latitude"); call erreur(status,.TRUE.,"put_att_lat_ID")
status = NF90_PUT_ATT(fidM,lon_ID,"name","lon"); call erreur(status,.TRUE.,"put_att_lon_ID")
status = NF90_PUT_ATT(fidM,lon_ID,"long_name","longitude coordinate"); call erreur(status,.TRUE.,"put_att_lon_ID")
status = NF90_PUT_ATT(fidM,lon_ID,"units","degrees_east"); call erreur(status,.TRUE.,"put_att_lon_ID")
status = NF90_PUT_ATT(fidM,lon_ID,"standard_name","longitude"); call erreur(status,.TRUE.,"put_att_lon_ID")
 
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"origin","Created using calculate_melt_rate_ensemble.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using rebuild_melt_large_ice_shelves.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imin",imin); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imax",imax); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmin",jmin); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmax",jmax); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
 
status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidM,index_tun_ID,index_tun); call erreur(status,.TRUE.,"var_index_tun_ID")
status = NF90_PUT_VAR(fidM,index_WOA_ID,index_WOA); call erreur(status,.TRUE.,"var_index_WOA_ID")
status = NF90_PUT_VAR(fidM,index_CMIP_ID,index_CMIP); call erreur(status,.TRUE.,"var_index_CMIP_ID")
status = NF90_PUT_VAR(fidM,index_para_ID,index_para); call erreur(status,.TRUE.,"var_index_para_ID")
!status = NF90_PUT_VAR(fidM,Ccoef_ID,Ccoef); call erreur(status,.TRUE.,"var_Ccoef_ID")
status = NF90_PUT_VAR(fidM,Kcoef_ID,Kcoef); call erreur(status,.TRUE.,"var_Kcoef_ID")
status = NF90_PUT_VAR(fidM,mean_melt_futu_ID,mean_melt_futu); call erreur(status,.TRUE.,"var_mean_melt_futu_ID")
status = NF90_PUT_VAR(fidM,total_melt_futu_ID,total_melt_futu); call erreur(status,.TRUE.,"var_total_melt_futu_ID")
status = NF90_PUT_VAR(fidM,mean_melt_pres_ID,mean_melt_pres); call erreur(status,.TRUE.,"var_mean_melt_pres_ID")
status = NF90_PUT_VAR(fidM,total_melt_pres_ID,total_melt_pres); call erreur(status,.TRUE.,"var_total_melt_pres_ID")
status = NF90_PUT_VAR(fidM,Futu_Melt_ID,Futu_Melt); call erreur(status,.TRUE.,"var_Futu_Melt_ID")
status = NF90_PUT_VAR(fidM,Pres_Melt_ID,Pres_Melt); call erreur(status,.TRUE.,"var_Pres_Melt_ID")
status = NF90_PUT_VAR(fidM,lat_ID,lat); call erreur(status,.TRUE.,"var_lat_ID")
status = NF90_PUT_VAR(fidM,lon_ID,lon); call erreur(status,.TRUE.,"var_lon_ID")
 
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
