program modif
 
USE netcdf
 
IMPLICIT NONE
 
INTEGER :: fidA, status, dimID_x, dimID_y, mx, my, Futu_Melt_ID, Pres_Melt_ID, x_ID, y_ID, fidM, Npara, kpara
 
INTEGER*4,ALLOCATABLE,DIMENSION(:) :: x, y
 
CHARACTER(LEN=100) :: file_in, file_out
 
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: Futu_Melt, Pres_Melt, Futu_Melt_mean, Pres_Melt_mean
 
Npara = 18

!---------------------------------------
 
file_out = 'melt_pattern_ANT3_DeltaT_ALL_STEREOGRAPHIC.nc'

!---------------------------------------
do kpara=1,Npara

  write(file_in,101) kpara
  101 FORMAT('melt_pattern_ANT3_DeltaT_kkpara',i2.2,'_STEREOGRAPHIC.nc') 
   
  write(*,*) 'Reading ', TRIM(file_in)
   
  status = NF90_OPEN(TRIM(file_in),0,fidA); call erreur(status,.TRUE.,"read")
   
  status = NF90_INQ_DIMID(fidA,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
  status = NF90_INQ_DIMID(fidA,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
   
  status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx); call erreur(status,.TRUE.,"inq_dim_x")
  status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my); call erreur(status,.TRUE.,"inq_dim_y")
    
  if ( kpara .eq. 1 ) then
    ALLOCATE(  Futu_Melt_mean(mx,my)  ) 
    ALLOCATE(  Pres_Melt_mean(mx,my)  ) 
    ALLOCATE(  x(mx)  )
    ALLOCATE(  y(my)  )
    Futu_Melt_mean = 0.e0
    Pres_Melt_mean = 0.e0
  endif
   
  ALLOCATE(  Futu_Melt(mx,my)  )
  status = NF90_INQ_VARID(fidA,"Futu_Melt",Futu_Melt_ID); call erreur(status,.TRUE.,"inq_Futu_Melt_ID")
  status = NF90_GET_VAR(fidA,Futu_Melt_ID,Futu_Melt); call erreur(status,.TRUE.,"getvar_Futu_Melt")
  Futu_Melt_mean(:,:) = Futu_Melt_mean(:,:) + Futu_Melt(:,:) * 1.e0 / Npara
  DEALLOCATE( Futu_Melt )

  ALLOCATE(  Pres_Melt(mx,my)  ) 
  status = NF90_INQ_VARID(fidA,"Pres_Melt",Pres_Melt_ID); call erreur(status,.TRUE.,"inq_Pres_Melt_ID")
  status = NF90_GET_VAR(fidA,Pres_Melt_ID,Pres_Melt); call erreur(status,.TRUE.,"getvar_Pres_Melt")
  Pres_Melt_mean(:,:) = Pres_Melt_mean(:,:) + Pres_Melt(:,:) * 1.e0 / Npara
  DEALLOCATE( Pres_Melt )

  if ( kpara .eq. 1 ) then
    status = NF90_INQ_VARID(fidA,"x",x_ID); call erreur(status,.TRUE.,"inq_x_ID")
    status = NF90_GET_VAR(fidA,x_ID,x); call erreur(status,.TRUE.,"getvar_x")
    status = NF90_INQ_VARID(fidA,"y",y_ID); call erreur(status,.TRUE.,"inq_y_ID")
    status = NF90_GET_VAR(fidA,y_ID,y); call erreur(status,.TRUE.,"getvar_y")
  endif
   
  status = NF90_CLOSE(fidA); call erreur(status,.TRUE.,"close_file")
   
enddo
 
!---------------------------------------
! Writing new netcdf file :
 
write(*,*) 'Creating ', TRIM(file_out)
 
status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"x",mx,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"y",my,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
  
status = NF90_DEF_VAR(fidM,"Futu_Melt",NF90_FLOAT,(/dimID_x,dimID_y/),Futu_Melt_ID); call erreur(status,.TRUE.,"def_var_Futu_Melt_ID")
status = NF90_DEF_VAR(fidM,"Pres_Melt",NF90_FLOAT,(/dimID_x,dimID_y/),Pres_Melt_ID); call erreur(status,.TRUE.,"def_var_Pres_Melt_ID")
status = NF90_DEF_VAR(fidM,"x",NF90_INT,(/dimID_x/),x_ID); call erreur(status,.TRUE.,"def_var_x_ID")
status = NF90_DEF_VAR(fidM,"y",NF90_INT,(/dimID_y/),y_ID); call erreur(status,.TRUE.,"def_var_y_ID")
 
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"title","Future melt rate"); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"long_name","Future ice shelf melt rate"); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Futu_Melt_ID,"units","m/yr"); call erreur(status,.TRUE.,"put_att_Futu_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"title","Present melt rate"); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"long_name","Present ice shelf melt rate"); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,Pres_Melt_ID,"units","m/yr"); call erreur(status,.TRUE.,"put_att_Pres_Melt_ID")
status = NF90_PUT_ATT(fidM,x_ID,"long_name","Cartesian x-coordinate"); call erreur(status,.TRUE.,"put_att_x_ID")
status = NF90_PUT_ATT(fidM,x_ID,"standard_name","projection_x_coordinate"); call erreur(status,.TRUE.,"put_att_x_ID")
status = NF90_PUT_ATT(fidM,x_ID,"units","meter"); call erreur(status,.TRUE.,"put_att_x_ID")
status = NF90_PUT_ATT(fidM,y_ID,"long_name","Cartesian y-coordinate"); call erreur(status,.TRUE.,"put_att_y_ID")
status = NF90_PUT_ATT(fidM,y_ID,"standard_name","projection_y_coordinate"); call erreur(status,.TRUE.,"put_att_y_ID")
status = NF90_PUT_ATT(fidM,y_ID,"units","meter"); call erreur(status,.TRUE.,"put_att_y_ID")
 
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using calculate_mean_melt_all_param_DeltaT.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
 
status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidM,Futu_Melt_ID,Futu_Melt_mean); call erreur(status,.TRUE.,"var_Futu_Melt_ID")
status = NF90_PUT_VAR(fidM,Pres_Melt_ID,Pres_Melt_mean); call erreur(status,.TRUE.,"var_Pres_Melt_ID")
status = NF90_PUT_VAR(fidM,x_ID,x); call erreur(status,.TRUE.,"var_x_ID")
status = NF90_PUT_VAR(fidM,y_ID,y); call erreur(status,.TRUE.,"var_y_ID")
 
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
