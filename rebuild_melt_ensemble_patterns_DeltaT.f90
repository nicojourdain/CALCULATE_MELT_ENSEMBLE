program modif
 
USE netcdf
 
IMPLICIT NONE
 
INTEGER :: fidA, fidG, status, dimID_lat, dimID_lon, mlat, mlon, lat_ID, lon_ID, dimID_mpara, dimID_mstat, &
&          dimID_Nisf, mmpara, mmstat, mNisf, index_tun_ID, index_WOA_ID, imin, imax, jmin, jmax, kisf,    &
&          index_CMIP_ID, index_para_ID, Ccoef_ID, Kcoef_ID, mean_melt_futu_ID, total_melt_futu_ID, fidM,  &
&          mean_melt_pres_ID, total_melt_pres_ID, Futu_Melt_ID, Pres_Melt_ID, mlonreg, mlatreg, mmpara2,   &
&          kk_para, ii, jj
 
CHARACTER(LEN=100) :: file_melt_isf, file_out, file_G
 
REAL*4,ALLOCATABLE,DIMENSION(:) :: lat, lon
 
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: Futu_Melt, Pres_Melt, Futu_Melt_reg, Pres_Melt_reg

LOGICAL :: ex

file_G  = 'RTopo-2.0.1_30sec_ice_shelf_mask_ANT3.nc'
 
!---------------------------------------
! Reading the global grid :
 
write(*,*) 'Reading ', TRIM(file_G)
 
status = NF90_OPEN(TRIM(file_G),0,fidG); call erreur(status,.TRUE.,"read RTOPO2 grid")
 
status = NF90_INQ_DIMID(fidG,"lat",dimID_lat); call erreur(status,.TRUE.,"inq_dimID_lat")
status = NF90_INQ_DIMID(fidG,"lon",dimID_lon); call erreur(status,.TRUE.,"inq_dimID_lon")
 
status = NF90_INQUIRE_DIMENSION(fidG,dimID_lat,len=mlat); call erreur(status,.TRUE.,"inq_dim_lat")
status = NF90_INQUIRE_DIMENSION(fidG,dimID_lon,len=mlon); call erreur(status,.TRUE.,"inq_dim_lon")
  
ALLOCATE(  lat(mlat)  ) 
ALLOCATE(  lon(mlon)  ) 
 
status = NF90_INQ_VARID(fidG,"lat",lat_ID); call erreur(status,.TRUE.,"inq_lat_ID")
status = NF90_INQ_VARID(fidG,"lon",lon_ID); call erreur(status,.TRUE.,"inq_lon_ID")
 
status = NF90_GET_VAR(fidG,lat_ID,lat); call erreur(status,.TRUE.,"getvar_lat")
status = NF90_GET_VAR(fidG,lon_ID,lon); call erreur(status,.TRUE.,"getvar_lon")
 
status = NF90_CLOSE(fidG); call erreur(status,.TRUE.,"close_file")

!------------------------------------------------------------------------------------

! Reading mmpara, mmstat, mNisf :
write(*,*) 'Reading melt_ensemble_DeltaT_066.nc'
status = NF90_OPEN('melt_ensemble_DeltaT_066.nc',0,fidA); call erreur(status,.TRUE.,"read mmpara, mmstat, mNisf")
status = NF90_INQ_DIMID(fidA,"mpara",dimID_mpara); call erreur(status,.TRUE.,"inq_dimID_mpara")
status = NF90_INQ_DIMID(fidA,"mstat",dimID_mstat); call erreur(status,.TRUE.,"inq_dimID_mstat")
status = NF90_INQ_DIMID(fidA,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_mpara,len=mmpara); call erreur(status,.TRUE.,"inq_dim_mpara")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_mstat,len=mmstat); call erreur(status,.TRUE.,"inq_dim_mstat")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_Nisf,len=mNisf); call erreur(status,.TRUE.,"inq_dim_Nisf")
status = NF90_CLOSE(fidA); call erreur(status,.TRUE.,"close_file")

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
111 FORMAT('melt_ensemble_DeltaT_',i3.3,'.nc')
112 FORMAT('melt_pattern_ANT3_DeltaT_kkpara',i2.2,'.nc')

DO kk_para=1,mmpara

   ALLOCATE( Futu_Melt(mlon,mlat), Pres_Melt(mlon,mlat) )
 
   write(file_out,112) kk_para

   DO kisf=2,mNisf
  
     write(file_melt_isf,111) kisf
     INQUIRE(file=file_melt_isf,exist=ex)  
 
     if ( .not. ex ) then

       write(*,*) '@@ WARNING @@ : ', TRIM(file_melt_isf), ' DOES NOT EXIST !!!'

     else

       !------------------------------------------
       ! Read melt rates of individual ice shelves
      
       write(*,*) 'Reading ', TRIM(file_melt_isf)
        
       status = NF90_OPEN(TRIM(file_melt_isf),0,fidA); call erreur(status,.TRUE.,"read melt")
        
       status = NF90_INQ_DIMID(fidA,"mpara",dimID_mpara); call erreur(status,.TRUE.,"inq_dimID_mpara")
       status = NF90_INQ_DIMID(fidA,"lat",dimID_lat); call erreur(status,.TRUE.,"inq_dimID_lat")
       status = NF90_INQ_DIMID(fidA,"lon",dimID_lon); call erreur(status,.TRUE.,"inq_dimID_lon")
        
       status = NF90_INQUIRE_DIMENSION(fidA,dimID_mpara,len=mmpara2); call erreur(status,.TRUE.,"inq_dim_mpara")
       status = NF90_INQUIRE_DIMENSION(fidA,dimID_lat,len=mlatreg); call erreur(status,.TRUE.,"inq_dim_lat")
       status = NF90_INQUIRE_DIMENSION(fidA,dimID_lon,len=mlonreg); call erreur(status,.TRUE.,"inq_dim_lon")
        
       if ( mmpara2 .ne. mmpara ) then
         write(*,*) '~!@#$%^&* DIMENSION MISMATCH (mpara) >>>>>>>>>>> STOP !!'
         stop
       endif
      
       ALLOCATE(  Futu_Melt_reg(mlonreg,mlatreg)  ) 
       ALLOCATE(  Pres_Melt_reg(mlonreg,mlatreg)  ) 
       
       status = NF90_INQ_VARID(fidA,"Futu_Melt",Futu_Melt_ID); call erreur(status,.TRUE.,"inq_Futu_Melt_ID")
       status = NF90_INQ_VARID(fidA,"Pres_Melt",Pres_Melt_ID); call erreur(status,.TRUE.,"inq_Pres_Melt_ID")
        
       status = NF90_GET_VAR(fidA,Futu_Melt_ID,Futu_Melt_reg,(/1, 1, kk_para/),(/mlonreg, mlatreg, 1/)); call erreur(status,.TRUE.,"getvar_Futu_Melt")
       status = NF90_GET_VAR(fidA,Pres_Melt_ID,Pres_Melt_reg,(/1, 1, kk_para/),(/mlonreg, mlatreg, 1/)); call erreur(status,.TRUE.,"getvar_Pres_Melt")
       
       status = NF90_GET_ATT(fidA,NF90_GLOBAL, "imin", imin)
       status = NF90_GET_ATT(fidA,NF90_GLOBAL, "imax", imax)
       status = NF90_GET_ATT(fidA,NF90_GLOBAL, "jmin", jmin)
       status = NF90_GET_ATT(fidA,NF90_GLOBAL, "jmax", jmax)
       write(*,212) imin, imax, jmin, jmax
       212 FORMAT('   > i = ',i5,' : ',i5,'   ;   j = ',i5,' : ',i5)
      
       status = NF90_CLOSE(fidA); call erreur(status,.TRUE.,"close_file")
        
       !---------------------------------------
       ! Gather everything in single variables:
  
       !Pres_Melt(imin:imax,jmin:jmax) = Pres_Melt_reg(:,:)
       !Futu_Melt(imin:imax,jmin:jmax) = Futu_Melt_reg(:,:)
       do ii=imin,imax
       do jj=jmin,jmax 
          if ( Pres_Melt_reg(ii-imin+1,jj-jmin+1) .eq. 0.e0 .and. Futu_Melt_reg(ii-imin+1,jj-jmin+1) .eq. 0.e0 ) then
            Pres_Melt_reg(ii-imin+1,jj-jmin+1) = NF90_FILL_FLOAT
            Futu_Melt_reg(ii-imin+1,jj-jmin+1) = NF90_FILL_FLOAT
          endif
          if ( Pres_Melt_reg(ii-imin+1,jj-jmin+1) .ne. NF90_FILL_FLOAT ) then
            Pres_Melt(ii,jj) = Pres_Melt_reg(ii-imin+1,jj-jmin+1)
          endif
          if ( Futu_Melt_reg(ii-imin+1,jj-jmin+1) .ne. NF90_FILL_FLOAT ) then
            Futu_Melt(ii,jj) = Futu_Melt_reg(ii-imin+1,jj-jmin+1)
          endif
       enddo
       enddo

       DEALLOCATE( Futu_Melt_reg, Pres_Melt_reg )

     endif

  ENDDO  !! kisf

  !---------------------------------------------
  ! Writing new netcdf file with melt patterns :
   
  write(*,*) 'Creating ', TRIM(file_out)
   
  status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')
   
  status = NF90_DEF_DIM(fidM,"lat",mlat,dimID_lat); call erreur(status,.TRUE.,"def_dimID_lat")
  status = NF90_DEF_DIM(fidM,"lon",mlon,dimID_lon); call erreur(status,.TRUE.,"def_dimID_lon")
    
  status = NF90_DEF_VAR(fidM,"Futu_Melt",NF90_FLOAT,(/dimID_lon,dimID_lat/),Futu_Melt_ID); call erreur(status,.TRUE.,"def_var_Futu_Melt_ID")
  status = NF90_DEF_VAR(fidM,"Pres_Melt",NF90_FLOAT,(/dimID_lon,dimID_lat/),Pres_Melt_ID); call erreur(status,.TRUE.,"def_var_Pres_Melt_ID")
  status = NF90_DEF_VAR(fidM,"lat",NF90_FLOAT,(/dimID_lat/),lat_ID); call erreur(status,.TRUE.,"def_var_lat_ID")
  status = NF90_DEF_VAR(fidM,"lon",NF90_FLOAT,(/dimID_lon/),lon_ID); call erreur(status,.TRUE.,"def_var_lon_ID")
   
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
   
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using rebuild_melt_ensemble_patterns_DeltaT.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"kk_para",kk_para); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
   
  status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
   
  status = NF90_PUT_VAR(fidM,Futu_Melt_ID,Futu_Melt); call erreur(status,.TRUE.,"var_Futu_Melt_ID")
  status = NF90_PUT_VAR(fidM,Pres_Melt_ID,Pres_Melt); call erreur(status,.TRUE.,"var_Pres_Melt_ID")
  status = NF90_PUT_VAR(fidM,lon_ID,lon); call erreur(status,.TRUE.,"var_lon_ID")
  status = NF90_PUT_VAR(fidM,lat_ID,lat); call erreur(status,.TRUE.,"var_lat_ID")
   
  status = NF90_CLOSE(fidM); call erreur(status,.TRUE.,"final")

  DEALLOCATE( Futu_Melt, Pres_Melt )

ENDDO !! kk_para
 

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
