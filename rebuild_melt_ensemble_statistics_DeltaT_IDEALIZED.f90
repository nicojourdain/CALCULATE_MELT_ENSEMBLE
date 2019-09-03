program modif
 
USE netcdf
 
IMPLICIT NONE
 
INTEGER :: fidA, status, dimID_mpara, dimID_mstat, imin, imax, jmin, jmax, mean_melt_pres_ID, total_melt_pres_ID,      &
&          dimID_Nisf, dimID_lat, dimID_lon, mmpara, mmstat, mNisf, mlat, mlon, index_K_ID, index_WOA_ID, kisf,        &
&          index_CMIP_ID, index_para_ID, Kcoef_ID, mean_melt_futu_ID, total_melt_futu_ID, fidSTAT, DeltaT_ID
 
CHARACTER(LEN=100) :: file_melt_isf, file_out
 
INTEGER*1,ALLOCATABLE,DIMENSION(:,:) :: index_K, index_WOA, index_CMIP, index_para,                &
&                                       index_K_reg, index_WOA_reg, index_CMIP_reg, index_para_reg
 
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: Kcoef, mean_melt_futu, total_melt_futu, mean_melt_pres, total_melt_pres,   &
&                                    Kcoef_reg, mean_melt_futu_reg, total_melt_futu_reg, mean_melt_pres_reg,&
&                                    total_melt_pres_reg, DeltaT, DeltaT_reg

LOGICAL :: ex, ll_init

mNisf = 216
 
file_out = 'melt_ensemble_statistics_ANT3_DeltaT_IDEALIZED.nc'

ll_init = .true.

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
111 FORMAT('melt_ensemble_DeltaT_IDEALIZED_',i3.3'.nc')

DO kisf=2,mNisf

  write(file_melt_isf,111) kisf

  INQUIRE(file=file_melt_isf,exist=ex)

  if ( .not. ex ) then

    write(*,*) '~!@#$%^* WARNING : ', TRIM(file_melt_isf), ' DOES NOT EXIST !!!'

  else

    !------------------------------------------
    ! Read melt rates of individual ice shelves
   
    write(*,*) 'Reading ', TRIM(file_melt_isf)
     
    status = NF90_OPEN(TRIM(file_melt_isf),0,fidA); call erreur(status,.TRUE.,"read")
     
    status = NF90_INQ_DIMID(fidA,"mpara",dimID_mpara); call erreur(status,.TRUE.,"inq_dimID_mpara")
    status = NF90_INQ_DIMID(fidA,"mstat",dimID_mstat); call erreur(status,.TRUE.,"inq_dimID_mstat")
    status = NF90_INQ_DIMID(fidA,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
     
    status = NF90_INQUIRE_DIMENSION(fidA,dimID_mpara,len=mmpara); call erreur(status,.TRUE.,"inq_dim_mpara")
    status = NF90_INQUIRE_DIMENSION(fidA,dimID_mstat,len=mmstat); call erreur(status,.TRUE.,"inq_dim_mstat")
    status = NF90_INQUIRE_DIMENSION(fidA,dimID_Nisf,len=mNisf); call erreur(status,.TRUE.,"inq_dim_Nisf")
      
    ALLOCATE(  index_K_reg(mNisf,mmstat)  ) 
    ALLOCATE(  index_WOA_reg(mNisf,mmstat)  ) 
    ALLOCATE(  index_CMIP_reg(mNisf,mmstat)  ) 
    ALLOCATE(  index_para_reg(mNisf,mmstat)  ) 
    ALLOCATE(  DeltaT_reg(mNisf,mmstat)  ) 
    ALLOCATE(  Kcoef_reg(mNisf,mmstat)  ) 
    ALLOCATE(  mean_melt_futu_reg(mNisf,mmstat)  ) 
    ALLOCATE(  total_melt_futu_reg(mNisf,mmstat)  ) 
    ALLOCATE(  mean_melt_pres_reg(mNisf,mmstat)  ) 
    ALLOCATE(  total_melt_pres_reg(mNisf,mmstat)  ) 
    
    if ( ll_init ) then
      ALLOCATE(  index_K(mNisf,mmstat)  ) 
      ALLOCATE(  index_WOA(mNisf,mmstat)  ) 
      ALLOCATE(  index_CMIP(mNisf,mmstat)  ) 
      ALLOCATE(  index_para(mNisf,mmstat)  ) 
      ALLOCATE(  DeltaT(mNisf,mmstat)  ) 
      ALLOCATE(  Kcoef(mNisf,mmstat)  ) 
      ALLOCATE(  mean_melt_futu(mNisf,mmstat)  ) 
      ALLOCATE(  total_melt_futu(mNisf,mmstat)  ) 
      ALLOCATE(  mean_melt_pres(mNisf,mmstat)  ) 
      ALLOCATE(  total_melt_pres(mNisf,mmstat)  ) 
      index_K  (:,:) = 0
      index_WOA  (:,:) = 0
      index_CMIP (:,:) = 0
      index_para (:,:) = 0
      DeltaT     (:,:) = NF90_FILL_FLOAT
      Kcoef      (:,:) = NF90_FILL_FLOAT
      mean_melt_futu (:,:) = NF90_FILL_FLOAT
      total_melt_futu(:,:) = NF90_FILL_FLOAT
      mean_melt_pres (:,:) = NF90_FILL_FLOAT
      total_melt_pres(:,:) = NF90_FILL_FLOAT
      ll_init = .false.
    endif
   
    status = NF90_INQ_VARID(fidA,"index_K",index_K_ID); call erreur(status,.TRUE.,"inq_index_K_ID")
    status = NF90_INQ_VARID(fidA,"index_WOA",index_WOA_ID); call erreur(status,.TRUE.,"inq_index_WOA_ID")
    status = NF90_INQ_VARID(fidA,"index_CMIP",index_CMIP_ID); call erreur(status,.TRUE.,"inq_index_CMIP_ID")
    status = NF90_INQ_VARID(fidA,"index_para",index_para_ID); call erreur(status,.TRUE.,"inq_index_para_ID")
    status = NF90_INQ_VARID(fidA,"DeltaT",DeltaT_ID); call erreur(status,.TRUE.,"inq_DeltaT_ID")
    status = NF90_INQ_VARID(fidA,"Kcoef",Kcoef_ID); call erreur(status,.TRUE.,"inq_Kcoef_ID")
    status = NF90_INQ_VARID(fidA,"mean_melt_futu",mean_melt_futu_ID); call erreur(status,.TRUE.,"inq_mean_melt_futu_ID")
    status = NF90_INQ_VARID(fidA,"total_melt_futu",total_melt_futu_ID); call erreur(status,.TRUE.,"inq_total_melt_futu_ID")
    status = NF90_INQ_VARID(fidA,"mean_melt_pres",mean_melt_pres_ID); call erreur(status,.TRUE.,"inq_mean_melt_pres_ID")
    status = NF90_INQ_VARID(fidA,"total_melt_pres",total_melt_pres_ID); call erreur(status,.TRUE.,"inq_total_melt_pres_ID")
     
    status = NF90_GET_VAR(fidA,index_K_ID,index_K_reg); call erreur(status,.TRUE.,"getvar_index_K")
    status = NF90_GET_VAR(fidA,index_WOA_ID,index_WOA_reg); call erreur(status,.TRUE.,"getvar_index_WOA")
    status = NF90_GET_VAR(fidA,index_CMIP_ID,index_CMIP_reg); call erreur(status,.TRUE.,"getvar_index_CMIP")
    status = NF90_GET_VAR(fidA,index_para_ID,index_para_reg); call erreur(status,.TRUE.,"getvar_index_para")
    status = NF90_GET_VAR(fidA,DeltaT_ID,DeltaT_reg); call erreur(status,.TRUE.,"getvar_DeltaT")
    status = NF90_GET_VAR(fidA,Kcoef_ID,Kcoef_reg); call erreur(status,.TRUE.,"getvar_Kcoef")
    status = NF90_GET_VAR(fidA,mean_melt_futu_ID,mean_melt_futu_reg); call erreur(status,.TRUE.,"getvar_mean_melt_futu")
    status = NF90_GET_VAR(fidA,total_melt_futu_ID,total_melt_futu_reg); call erreur(status,.TRUE.,"getvar_total_melt_futu")
    status = NF90_GET_VAR(fidA,mean_melt_pres_ID,mean_melt_pres_reg); call erreur(status,.TRUE.,"getvar_mean_melt_pres")
    status = NF90_GET_VAR(fidA,total_melt_pres_ID,total_melt_pres_reg); call erreur(status,.TRUE.,"getvar_total_melt_pres")
   
    write(*,*) '   Kcoef_reg in ', MINVAL(Kcoef_reg(kisf,:)), ' : ', MAXVAL(Kcoef_reg(kisf,:))
    write(*,*) '   mean_melt_pres_reg in ', MINVAL(mean_melt_pres_reg(kisf,:)), ' : ', MAXVAL(mean_melt_pres_reg(kisf,:))
 
    status = NF90_CLOSE(fidA); call erreur(status,.TRUE.,"close_file")
     
    !---------------------------------------
    ! Gather everything in single variables:
  
    index_K  (kisf,:) = index_K_reg  (kisf,:)
    index_WOA  (kisf,:) = index_WOA_reg  (kisf,:)
    index_CMIP (kisf,:) = index_CMIP_reg (kisf,:)
    index_para (kisf,:) = index_para_reg (kisf,:)
    DeltaT     (kisf,:) = DeltaT_reg     (kisf,:)
    Kcoef      (kisf,:) = Kcoef_reg      (kisf,:)
    mean_melt_futu (kisf,:) = mean_melt_futu_reg (kisf,:)
    total_melt_futu(kisf,:) = total_melt_futu_reg(kisf,:)
    mean_melt_pres (kisf,:) = mean_melt_pres_reg (kisf,:)
    total_melt_pres(kisf,:) = total_melt_pres_reg(kisf,:)

    write(*,*) '   Kcoef in ', MINVAL(Kcoef(kisf,:)), ' : ', MAXVAL(Kcoef(kisf,:))
    write(*,*) '   mean_melt_pres in ', MINVAL(mean_melt_pres(kisf,:)), ' : ', MAXVAL(mean_melt_pres(kisf,:))
  
    DEALLOCATE( index_K_reg, index_WOA_reg, index_CMIP_reg, index_para_reg, DeltaT_reg, Kcoef_reg )
    DEALLOCATE( mean_melt_futu_reg, total_melt_futu_reg, mean_melt_pres_reg, total_melt_pres_reg )

  endif

ENDDO

write(*,*) 'Kcoef in ', MINVAL(Kcoef), ' : ', MAXVAL(Kcoef)
write(*,*) 'mean_melt_pres in ', MINVAL(mean_melt_pres), ' : ', MAXVAL(mean_melt_pres)
 
!---------------------------------------
! Writing new netcdf file with statistics :
 
write(*,*) 'Creating ', TRIM(file_out)
 
status = NF90_CREATE(TRIM(file_out),or(NF90_CLOBBER,NF90_64BIT_OFFSET),fidSTAT); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidSTAT,"mpara",mmpara,dimID_mpara); call erreur(status,.TRUE.,"def_dimID_mpara")
status = NF90_DEF_DIM(fidSTAT,"mstat",mmstat,dimID_mstat); call erreur(status,.TRUE.,"def_dimID_mstat")
status = NF90_DEF_DIM(fidSTAT,"Nisf",mNisf,dimID_Nisf); call erreur(status,.TRUE.,"def_dimID_Nisf")
  
status = NF90_DEF_VAR(fidSTAT,"index_K",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_K_ID); call erreur(status,.TRUE.,"def_var_index_K_ID")
status = NF90_DEF_VAR(fidSTAT,"index_WOA",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_WOA_ID); call erreur(status,.TRUE.,"def_var_index_WOA_ID")
status = NF90_DEF_VAR(fidSTAT,"index_CMIP",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_CMIP_ID); call erreur(status,.TRUE.,"def_var_index_CMIP_ID")
status = NF90_DEF_VAR(fidSTAT,"index_para",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_para_ID); call erreur(status,.TRUE.,"def_var_index_para_ID")
status = NF90_DEF_VAR(fidSTAT,"DeltaT",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),DeltaT_ID); call erreur(status,.TRUE.,"def_var_DeltaT_ID")
status = NF90_DEF_VAR(fidSTAT,"Kcoef",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),Kcoef_ID); call erreur(status,.TRUE.,"def_var_Kcoef_ID")
status = NF90_DEF_VAR(fidSTAT,"mean_melt_futu",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),mean_melt_futu_ID); call erreur(status,.TRUE.,"def_var_mean_melt_futu_ID")
status = NF90_DEF_VAR(fidSTAT,"total_melt_futu",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),total_melt_futu_ID); call erreur(status,.TRUE.,"def_var_total_melt_futu_ID")
status = NF90_DEF_VAR(fidSTAT,"mean_melt_pres",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),mean_melt_pres_ID); call erreur(status,.TRUE.,"def_var_mean_melt_pres_ID")
status = NF90_DEF_VAR(fidSTAT,"total_melt_pres",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),total_melt_pres_ID); call erreur(status,.TRUE.,"def_var_total_melt_pres_ID")
 
status = NF90_PUT_ATT(fidSTAT,index_K_ID,"title","index_K"); call erreur(status,.TRUE.,"put_att_index_K_ID")
status = NF90_PUT_ATT(fidSTAT,index_K_ID,"long_name","index defining the factor value"); call erreur(status,.TRUE.,"put_att_index_K_ID")
status = NF90_PUT_ATT(fidSTAT,index_K_ID,"units","-"); call erreur(status,.TRUE.,"put_att_index_K_ID")
status = NF90_PUT_ATT(fidSTAT,index_WOA_ID,"title","index_WOA"); call erreur(status,.TRUE.,"put_att_index_WOA_ID")
status = NF90_PUT_ATT(fidSTAT,index_WOA_ID,"latg_name","index defining the sampling of WOA2013 range of uncertainty"); call erreur(status,.TRUE.,"put_att_index_WOA_ID")
status = NF90_PUT_ATT(fidSTAT,index_WOA_ID,"units","-"); call erreur(status,.TRUE.,"put_att_index_WOA_ID")
status = NF90_PUT_ATT(fidSTAT,index_CMIP_ID,"title","index_CMIP"); call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
status = NF90_PUT_ATT(fidSTAT,index_CMIP_ID,"long_name","index defining the CMIP5 model used for the climate anomaly"); call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
status = NF90_PUT_ATT(fidSTAT,index_CMIP_ID,"units","-"); call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
status = NF90_PUT_ATT(fidSTAT,index_para_ID,"title","index_para"); call erreur(status,.TRUE.,"put_att_index_para_ID")
status = NF90_PUT_ATT(fidSTAT,index_para_ID,"long_name","index defining the melt parameterization"); call erreur(status,.TRUE.,"put_att_index_para_ID")
status = NF90_PUT_ATT(fidSTAT,index_para_ID,"units","-"); call erreur(status,.TRUE.,"put_att_index_para_ID")
status = NF90_PUT_ATT(fidSTAT,DeltaT_ID,"title","DeltaT"); call erreur(status,.TRUE.,"put_att_DeltaT_ID")
status = NF90_PUT_ATT(fidSTAT,DeltaT_ID,"long_name","Thermal Forcing correction"); call erreur(status,.TRUE.,"put_att_DeltaT_ID")
status = NF90_PUT_ATT(fidSTAT,DeltaT_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_DeltaT_ID")
status = NF90_PUT_ATT(fidSTAT,DeltaT_ID,"units","K"); call erreur(status,.TRUE.,"put_att_DeltaT_ID")
status = NF90_PUT_ATT(fidSTAT,Kcoef_ID,"title","K coeff"); call erreur(status,.TRUE.,"put_att_Kcoef_ID")
status = NF90_PUT_ATT(fidSTAT,Kcoef_ID,"long_name","K coefficient (tuning factor)"); call erreur(status,.TRUE.,"put_att_Kcoef_ID")
status = NF90_PUT_ATT(fidSTAT,Kcoef_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_Kcoef_ID")
status = NF90_PUT_ATT(fidSTAT,Kcoef_ID,"units","TBD"); call erreur(status,.TRUE.,"put_att_Kcoef_ID")
status = NF90_PUT_ATT(fidSTAT,mean_melt_futu_ID,"title","Future Melt Rate"); call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
status = NF90_PUT_ATT(fidSTAT,mean_melt_futu_ID,"long_name","Future mean melt rate over the cavity"); call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
status = NF90_PUT_ATT(fidSTAT,mean_melt_futu_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
status = NF90_PUT_ATT(fidSTAT,mean_melt_futu_ID,"units","m/yr"); call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
status = NF90_PUT_ATT(fidSTAT,total_melt_futu_ID,"title","Future Melt Flux"); call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
status = NF90_PUT_ATT(fidSTAT,total_melt_futu_ID,"long_name","Total future cavity melt flux"); call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
status = NF90_PUT_ATT(fidSTAT,total_melt_futu_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
status = NF90_PUT_ATT(fidSTAT,total_melt_futu_ID,"units","Gt/yr"); call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
status = NF90_PUT_ATT(fidSTAT,mean_melt_pres_ID,"title","Present Melt Rate"); call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
status = NF90_PUT_ATT(fidSTAT,mean_melt_pres_ID,"long_name","Present mean melt rate over the cavity"); call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
status = NF90_PUT_ATT(fidSTAT,mean_melt_pres_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
status = NF90_PUT_ATT(fidSTAT,mean_melt_pres_ID,"units","m/yr"); call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
status = NF90_PUT_ATT(fidSTAT,total_melt_pres_ID,"title","Present Melt Flux"); call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
status = NF90_PUT_ATT(fidSTAT,total_melt_pres_ID,"long_name","Total present cavity melt flux"); call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
status = NF90_PUT_ATT(fidSTAT,total_melt_pres_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
status = NF90_PUT_ATT(fidSTAT,total_melt_pres_ID,"units","Gt/yr"); call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
 
status = NF90_PUT_ATT(fidSTAT,NF90_GLOBAL,"history","Created using rebuild_melt_ensemble_statistics_DeltaT_IDEALIZED.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
 
status = NF90_ENDDEF(fidSTAT); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidSTAT,index_K_ID,index_K); call erreur(status,.TRUE.,"var_index_K_ID")
status = NF90_PUT_VAR(fidSTAT,index_WOA_ID,index_WOA); call erreur(status,.TRUE.,"var_index_WOA_ID")
status = NF90_PUT_VAR(fidSTAT,index_CMIP_ID,index_CMIP); call erreur(status,.TRUE.,"var_index_CMIP_ID")
status = NF90_PUT_VAR(fidSTAT,index_para_ID,index_para); call erreur(status,.TRUE.,"var_index_para_ID")
status = NF90_PUT_VAR(fidSTAT,DeltaT_ID,DeltaT); call erreur(status,.TRUE.,"var_DeltaT_ID")
status = NF90_PUT_VAR(fidSTAT,Kcoef_ID,Kcoef); call erreur(status,.TRUE.,"var_Kcoef_ID")
status = NF90_PUT_VAR(fidSTAT,mean_melt_futu_ID,mean_melt_futu); call erreur(status,.TRUE.,"var_mean_melt_futu_ID")
status = NF90_PUT_VAR(fidSTAT,total_melt_futu_ID,total_melt_futu); call erreur(status,.TRUE.,"var_total_melt_futu_ID")
status = NF90_PUT_VAR(fidSTAT,mean_melt_pres_ID,mean_melt_pres); call erreur(status,.TRUE.,"var_mean_melt_pres_ID")
status = NF90_PUT_VAR(fidSTAT,total_melt_pres_ID,total_melt_pres); call erreur(status,.TRUE.,"var_total_melt_pres_ID")
 
status = NF90_CLOSE(fidSTAT); call erreur(status,.TRUE.,"final")

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
