program ensemble_melt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! N. JOURDAIN, IGE-CNRS, Grenoble, August 2018
!
! Purpose : Calculate ensemble ice-shelf melt rates from multiple parameterizations and inputs. 
! ~~~~~~~~~
!
! Required input : 
! ~~~~~~~~~~~~~~~~
!                  - RTopo-2.0.1_30sec_ice_shelf_mask.nc (built by identify_ice_shelf_mask.f90)
!                  - WOA13_TS_per_ice_sehlf.nc (built by calculate_TS_profiles_WOA13.f90) 
!                  - CMIP5_anom_TS_per_ice_shelf.nc (built by calculate_TS_profiles_CMIP5_anom.f90)
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


USE netcdf

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER(KIND=2),ALLOCATABLE,DIMENSION(:,:) :: isfmask,           & ! 0 grounded; 1 ocean; >1 ice shelf number
&                                             IF_mask, GL_mask,  & ! Ice shelf Front (IF) and Grounding Line (GL) mask
&                                             index_para, index_CMIP, index_WOA, index_K

!----------------------------------------------------------------------------------------------------

INTEGER(KIND=4) :: nD, ii, jj, kk, kki, kkj, max_nb_box, kkp, dimID_mstat, dimID_mpara, dimID_StrLen,stun, nn_TS,    &
&                  Pmelt_ID, Fmelt_ID, total_melt_pres_ID, mean_melt_pres_ID, total_melt_futu_ID, mean_melt_futu_ID, &
&                  index_para_ID, index_CMIP_ID, index_WOA_ID, index_K_ID, nD1, nD2, nD3, im1, ip1, jm1, jp1, rs,    &
&                  fidBOX, fidPLUME, Abox_ID, Zbox_ID, zGL_ID, alpha_ID, dGL_ID, dIF_ID, nnp, dimID_nnp, nnp_ID,     &
&                  dimID_max_nb_box1, max_nb_box1, dimID_max_nb_box2, max_nb_box2, Kcoef_ID, status, fidM, cpttmp,   &
&                  aa, nbDir, nnp2, Ntun, interval, cptglob, imin_ID, jmin_ID, imax_ID, jmax_ID, cond12, cond23,     &
&                  cond13, fidWOA, dimID_depth, mdepth, mdepth2, mNisf, mNisf2, s_er_3_ID, s_an_3_ID, t_er_3_ID,     &
&                  t_an_3_ID, s_er_2_ID, s_an_2_ID, t_er_2_ID, t_an_2_ID, s_er_1_ID, s_an_1_ID, t_er_1_ID, t_an_1_ID,&
&                  fidCMIP, dimID_Nmod, mStrLen, mNmod, s_anom_3_ID, t_anom_3_ID, s_anom_2_ID, t_anom_2_ID, kkinf,   &
&                  s_anom_1_ID, t_anom_1_ID, depth_ID, front_max_lat_ID, front_min_lat_ID, front_max_lon_ID, kksup,  &
&                  front_min_lon_ID, front_ice_dep_avg_ID, front_ice_dep_min_ID, front_bot_dep_avg_ID, nn_para, mtot,&
&                  front_bot_dep_max_ID, kisf, kk_TS_pres, nn_TS_pres, kk_tuning, nn_tuning, kk_CMIP5_anom, kk_para, &
&                  fidBED, dimID_lon, dimID_lat, mlondim, mlatdim, lon_ID, lat_ID, bedrock_topography_ID, fidISF,    &
&                  isfmask_ID, fidDFT, ice_base_topography_ID, dimID_Nisf, err_melt_isf_ID, melt_isf_ID, fidCMIPbot, &
&                  IF_mask_ID, GL_mask_ID, s_er_4_ID, s_an_4_ID, t_er_4_ID, t_an_4_ID, s_er_5_ID, s_an_5_ID, fidSTAT,&
&                  t_er_5_ID, t_an_5_ID, fidSHMIDTKO, s_anom_5_ID, t_anom_5_ID, s_anom_4_ID, t_anom_4_ID,            &
&                  kk_K, nn_K, mstat, mpara2, DeltaT_ID, min_melt_pres_ID, max_melt_pres_ID

INTEGER(KIND=4),ALLOCATABLE,DIMENSION(:) :: kstat, imin, imax, jmin, jmax, indx2

INTEGER(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: Pnnnn, Fnnnn

!----------------------------------------------------------------------------------------------------

REAL(KIND=4) :: dlon, dlat

REAL(KIND=4),ALLOCATABLE,DIMENSION(:) :: lon, lat, err_melt_isf, melt_isf, front_max_lat, front_min_lat, S_futu,     &
&                                        front_max_lon, front_min_lon, front_ice_dep_avg, front_ice_dep_min, T_futu, &
&                                        front_bot_dep_avg, front_bot_dep_max, depth, T_pres, S_pres, s_er_4, s_an_4,&
&                                        t_er_4, t_an_4, s_er_5, s_an_5, t_er_5, t_an_5, perct, Kcoef2

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: bedrock_topography, ice_base_topography, s_er_3, s_an_3, t_er_3, t_an_3,  &
&                                          s_er_2, s_an_2, t_er_2, t_an_2, s_er_1, s_an_1, t_er_1, t_an_1, Kcoef,    &
&                                          total_melt_pres, mean_melt_pres, total_melt_futu, mean_melt_futu,         &
&                                          s_anom_5, t_anom_5, s_anom_4, t_anom_4, Kcoef_pct, DeltaT_out,            &
&                                          min_melt_pres, max_melt_pres

REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: Pmelt_out, Fmelt_out, s_anom_3, t_anom_3, s_anom_2, t_anom_2,         &
&                                              s_anom_1, t_anom_1

!----------------------------------------------------------------------------------------------------

REAL(KIND=8) ::  zzz, T0c, T0, S0, Tf, lbd1, lbd2, lbd3, meltfac, K, gT, alphap, wn, x0, E0, Cd, GamT, GefT, GTS, MM,&
&                M_hat, X_hat, ll, rhostar, CC, beta, g1, g2, rr, qqq, xerr, isfslope_x, isfslope_y, tmp1, Tstar,    &
&                xbox, TF_avg, TF_tot, H0, gg, gg_SI, det, lonP2, latP2, lonP3, latP3, radius, angle, TT, SS, incr,  &
&                yearinsec, rhoi_SI, rhosw_SI, rhofw_SI, Lf_SI, cpw_SI, grav_SI, lambda1, lambda2, lambda3, DeltaT,  &
&                rhoi, rhosw, rhofw, Lf, cpw, gravity, gammaT_SI, facPDC_SI, E0_LAZER, GT_LAZER, zz1, zz2, area,     &
&                GTS0_LAZER, gammaT_PICO_SI, C_PICO_SI, alpha_PICO, beta_PICO, rhostar_SI_PICO, snt, angle1, angle2, &
&                gammaT, facPDC, gammaT_PICO, C_PICO, rhostar_PICO, dist, lonGLmin, latGLmin, zGLmin, sn, distmp,    &
&                zGLtmp, RT, deg2rad,  aainf, aasup, target_melt, mmm_tot_p, mmm_avg_p, mmm_tot_f, mmm_avg_f, aaa,   &
&                maxDeltaT, ctau, cr1, alpha_LAZER, beta_LAZER, ttt, ttt_SI, m_closest, reldif, DeltaT_closest,      &
&                mmm_min_p, mmm_max_p

REAL(KIND=8), DIMENSION(12) :: pp

REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::  Tbox, Sbox, Mbox

REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::  Zbox, Abox, dGL, dIF, Melt

REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: zGL, alpha

!----------------------------------------------------------------------------------------------------

REAL(KIND=16), DIMENSION(:,:,:), ALLOCATABLE :: Pmelt, Fmelt

!----------------------------------------------------------------------------------------------------

CHARACTER(LEN=150) :: file_out, file_box, file_plume, file_bed, file_isf, file_dft, file_in_stat, &
&                     file_TS_WOA, file_TS_SHMIDTKO, file_CMIP5_anom, file_CMIP5_anom_bot

CHARACTER(LEN=10) ::  TypeL

!----------------------------------------------------------------------------------------------------

LOGICAL :: llslp, llnD, ex, ll_curr_undr, ll_prev_undr, ll_curr_over, ll_prev_over

!----------------------------------------------------------------------------------------------------

file_dft            = '/store/njourd/DATA/DATA_BATHYMETRY/RTopo-2.0.1_30sec_ice_base_topography_ANT3.nc'
file_bed            = '/store/njourd/DATA/DATA_BATHYMETRY/RTopo-2.0.1_30sec_bedrock_topography_ANT3.nc'
file_isf            = 'RTopo-2.0.1_30sec_ice_shelf_mask_ANT3.nc'
file_in_stat        = 'melt_ensemble_statistics_ANT3.nc'
file_TS_WOA         = 'TS_per_ice_shelf_ANT3_WOA2018_EN4_MEOP.nc'
!file_TS_SHMIDTKO    = 'TS_per_ice_shelf_ANT3_SCHMIDTKO2014.nc'
file_CMIP5_anom     = 'CMIP5_anom_TS_per_ice_shelf_ANT3.nc'
file_CMIP5_anom_bot = 'CMIP5_anom_TS_bottom_per_ice_shelf_ANT3.nc'

deg2rad = dacos(-1.d0) / 180.d0
RT = 6.371d6 ! Earth radius in m
write(*,*) 'pi = ', dacos(-1.d0)

nn_K       = 17  ! =19 for every 5th percentile, NB : keep spaces "nn_K       ="
nn_para    = 18

! Nb of iterations to find deltaT
Ntun = 50

! Bounds DeltaT : abs(DeltaT) < maxDeltaT
maxDeltaT = 5.d0

!------------------------------------------------------------------------------
! 0- Read index and properties specific to individual ice shelves
!------------------------------------------------------------------------------

write(*,*) 'Reading isf properties in ', TRIM(file_isf)

status = NF90_OPEN(TRIM(file_isf),0,fidISF)          
call erreur(status,.TRUE.,"read ice shelf mask") 

status = NF90_INQ_DIMID(fidISF,"Nisf",dimID_Nisf)
call erreur(status,.TRUE.,"inq_dimID_Nisf")

status = NF90_INQUIRE_DIMENSION(fidISF,dimID_Nisf,len=mNisf)
call erreur(status,.TRUE.,"inq_dim_Nisf")

ALLOCATE(  melt_isf(mNisf)  )
ALLOCATE(  err_melt_isf(mNisf)  )
ALLOCATE(  front_max_lat(mNisf)  ) 
ALLOCATE(  front_min_lat(mNisf)  ) 
ALLOCATE(  front_max_lon(mNisf)  ) 
ALLOCATE(  front_min_lon(mNisf)  ) 
ALLOCATE(  front_ice_dep_avg(mNisf)  ) 
ALLOCATE(  front_ice_dep_min(mNisf)  ) 
ALLOCATE(  front_bot_dep_avg(mNisf)  ) 
ALLOCATE(  front_bot_dep_max(mNisf)  ) 
ALLOCATE(  imin(mNisf), imax(mNisf), jmin(mNisf), jmax(mNisf) )

status = NF90_INQ_VARID(fidISF,"err_melt_isf",err_melt_isf_ID)
call erreur(status,.TRUE.,"inq_err_melt_isf_ID")
status = NF90_INQ_VARID(fidISF,"melt_isf",melt_isf_ID)
call erreur(status,.TRUE.,"inq_melt_isf_ID")
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
status = NF90_INQ_VARID(fidISF,"imin",imin_ID)
call erreur(status,.TRUE.,"inq_imin_ID")
status = NF90_INQ_VARID(fidISF,"imax",imax_ID)
call erreur(status,.TRUE.,"inq_imax_ID")
status = NF90_INQ_VARID(fidISF,"jmin",jmin_ID)
call erreur(status,.TRUE.,"inq_jmin_ID")
status = NF90_INQ_VARID(fidISF,"jmax",jmax_ID)
call erreur(status,.TRUE.,"inq_jmax_ID")

status = NF90_GET_VAR(fidISF,err_melt_isf_ID,err_melt_isf)
call erreur(status,.TRUE.,"getvar_err_melt_isf")
status = NF90_GET_VAR(fidISF,melt_isf_ID,melt_isf)
call erreur(status,.TRUE.,"getvar_melt_isf")
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
status = NF90_GET_VAR(fidISF,imin_ID,imin)
call erreur(status,.TRUE.,"getvar_imin")
status = NF90_GET_VAR(fidISF,imax_ID,imax)
call erreur(status,.TRUE.,"getvar_imax")
status = NF90_GET_VAR(fidISF,jmin_ID,jmin)
call erreur(status,.TRUE.,"getvar_jmin")
status = NF90_GET_VAR(fidISF,jmax_ID,jmax)
call erreur(status,.TRUE.,"getvar_jmax")

status = NF90_CLOSE(fidISF)                      
call erreur(status,.TRUE.,"fin_lecture")     

!----------------------------------------------------------
! Read statistics from first calculation (with DeltaT = 0) :
 
write(*,*) 'Reading ', TRIM(file_in_stat)
 
status = NF90_OPEN(TRIM(file_in_stat),0,fidSTAT); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidSTAT,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
status = NF90_INQ_DIMID(fidSTAT,"mstat",dimID_mstat); call erreur(status,.TRUE.,"inq_dimID_mstat")
status = NF90_INQ_DIMID(fidSTAT,"mpara",dimID_mpara); call erreur(status,.TRUE.,"inq_dimID_mpara")
 
status = NF90_INQUIRE_DIMENSION(fidSTAT,dimID_Nisf,len=mNisf2); call erreur(status,.TRUE.,"inq_dim_Nisf")
status = NF90_INQUIRE_DIMENSION(fidSTAT,dimID_mstat,len=mstat); call erreur(status,.TRUE.,"inq_dim_mstat")
status = NF90_INQUIRE_DIMENSION(fidSTAT,dimID_mpara,len=mpara2); call erreur(status,.TRUE.,"inq_dim_mpara")
  
if ( mNisf .ne. mNisf2 .or. mpara2 .ne. nn_para ) then
  write(*,*) '~!@#$%^*  PB WITH THE NB OF ICE SHELVES or PARAM. BETW. THE DIFF. FILES  >>>> STOP !'
  stop
endif

ALLOCATE(  Kcoef(mNisf,mstat)  ) 
ALLOCATE(  index_para(mNisf,mstat)  ) 

status = NF90_INQ_VARID(fidSTAT,"Kcoef",Kcoef_ID); call erreur(status,.TRUE.,"inq_Kcoef_ID")
status = NF90_INQ_VARID(fidSTAT,"index_para",index_para_ID); call erreur(status,.TRUE.,"inq_index_para_ID")

status = NF90_GET_VAR(fidSTAT,Kcoef_ID,Kcoef); call erreur(status,.TRUE.,"getvar_Kcoef")
status = NF90_GET_VAR(fidSTAT,index_para_ID,index_para); call erreur(status,.TRUE.,"getvar_index_para")

status = NF90_CLOSE(fidSTAT); call erreur(status,.TRUE.,"close_file")

!----- calculate percentiles for each parameterization :

ALLOCATE( perct(nn_K), Kcoef_pct(nn_para,nn_K) )

do kk_K=1,nn_K
  perct(kk_K) = 1.e0 * kk_K / (nn_K+1) ! value of the evenly-spread percentiles
enddo
write(*,*) 'Calculating percentiles ', perct(:)

DO kk_para = 1, nn_para

  write(*,*) '    for param. ', kk_para  

  nnp=0
  !do ii=1,mNisf
  do ii=10,75   !! NB: only take coeff of known ice shelves :
  do jj=1,mstat
    if ( index_para(ii,jj) .eq. kk_para .and. Kcoef(ii,jj) .ne. NF90_FILL_FLOAT ) nnp=nnp+1
  enddo
  enddo
  write(*,*) '          over ', nnp, ' values'
  
  ALLOCATE( Kcoef2(nnp), indx2(nnp) )
  nnp=0
  !do ii=1,mNisf
  do ii=10,75   !! NB: only take coeff of known ice shelves :
  do jj=1,mstat
    if ( index_para(ii,jj) .eq. kk_para .and. Kcoef(ii,jj) .ne. NF90_FILL_FLOAT ) then
      nnp=nnp+1
      Kcoef2(nnp) = Kcoef(ii,jj)
      indx2(nnp) = nnp
    endif
  enddo
  enddo

  write(*,*) 'MIN/MAX Kcoef = ', MINVAL(Kcoef), MAXVAL(Kcoef)
  write(*,*) 'MIN/MAX Kcoef2 = ', MINVAL(Kcoef2), MAXVAL(Kcoef2)
  write(*,*) 'MIN/MAX indx2 = ', MINVAL(indx2), MAXVAL(indx2)

  CALL hpsort_eps_epw(nnp,Kcoef2,indx2,1.d-20) ! this sorts Kcoef2
  do kk_K=1,nn_K
    Kcoef_pct(kk_para,kk_K) = Kcoef2(NINT(perct(kk_K)*nnp))
    write(*,*) 'kk_para, kk_K, NINT(perct(kk_K)*nnp), Kcoef_pct = ', kk_para, kk_K, NINT(perct(kk_K)*nnp), Kcoef_pct(kk_para,kk_K)
  enddo

  DEALLOCATE( Kcoef2, indx2 )

ENDDO

DEALLOCATE( Kcoef, index_para )

!---------------------------------------
! Read 3d TS in WOA :
 
write(*,*) 'Reading ', TRIM(file_TS_WOA)
 
status = NF90_OPEN(TRIM(file_TS_WOA),0,fidWOA); call erreur(status,.TRUE.,"read WOA T,S profiles")
 
status = NF90_INQ_DIMID(fidWOA,"depth",dimID_depth); call erreur(status,.TRUE.,"inq_dimID_depth")
status = NF90_INQ_DIMID(fidWOA,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
 
status = NF90_INQUIRE_DIMENSION(fidWOA,dimID_depth,len=mdepth); call erreur(status,.TRUE.,"inq_dim_depth")
status = NF90_INQUIRE_DIMENSION(fidWOA,dimID_Nisf,len=mNisf); call erreur(status,.TRUE.,"inq_dim_Nisf")
  
ALLOCATE(  s_an_3(mNisf,mdepth)  ) 
ALLOCATE(  t_an_3(mNisf,mdepth)  ) 
ALLOCATE(  s_an_2(mNisf,mdepth)  ) 
ALLOCATE(  t_an_2(mNisf,mdepth)  ) 
ALLOCATE(  s_an_1(mNisf,mdepth)  ) 
ALLOCATE(  t_an_1(mNisf,mdepth)  ) 
ALLOCATE(  depth(mdepth)  ) 
 
status = NF90_INQ_VARID(fidWOA,"s_an_3",s_an_3_ID); call erreur(status,.TRUE.,"inq_s_an_3_ID")
status = NF90_INQ_VARID(fidWOA,"t_an_3",t_an_3_ID); call erreur(status,.TRUE.,"inq_t_an_3_ID")
status = NF90_INQ_VARID(fidWOA,"s_an_2",s_an_2_ID); call erreur(status,.TRUE.,"inq_s_an_2_ID")
status = NF90_INQ_VARID(fidWOA,"t_an_2",t_an_2_ID); call erreur(status,.TRUE.,"inq_t_an_2_ID")
status = NF90_INQ_VARID(fidWOA,"s_an_1",s_an_1_ID); call erreur(status,.TRUE.,"inq_s_an_1_ID")
status = NF90_INQ_VARID(fidWOA,"t_an_1",t_an_1_ID); call erreur(status,.TRUE.,"inq_t_an_1_ID")
status = NF90_INQ_VARID(fidWOA,"depth",depth_ID); call erreur(status,.TRUE.,"inq_depth_ID")
 
status = NF90_GET_VAR(fidWOA,s_an_3_ID,s_an_3); call erreur(status,.TRUE.,"getvar_s_an_3")
status = NF90_GET_VAR(fidWOA,t_an_3_ID,t_an_3); call erreur(status,.TRUE.,"getvar_t_an_3")
status = NF90_GET_VAR(fidWOA,s_an_2_ID,s_an_2); call erreur(status,.TRUE.,"getvar_s_an_2")
status = NF90_GET_VAR(fidWOA,t_an_2_ID,t_an_2); call erreur(status,.TRUE.,"getvar_t_an_2")
status = NF90_GET_VAR(fidWOA,s_an_1_ID,s_an_1); call erreur(status,.TRUE.,"getvar_s_an_1")
status = NF90_GET_VAR(fidWOA,t_an_1_ID,t_an_1); call erreur(status,.TRUE.,"getvar_t_an_1")
status = NF90_GET_VAR(fidWOA,depth_ID,depth); call erreur(status,.TRUE.,"getvar_depth")
 
! Close file
status = NF90_CLOSE(fidWOA); call erreur(status,.TRUE.,"close_file")

!!--------------------------------------------------
!! Read sea-floor TS from Schmidtko et al. (2014) :
! 
!write(*,*) 'Reading ', TRIM(file_TS_SHMIDTKO)
! 
!status = NF90_OPEN(TRIM(file_TS_SHMIDTKO),0,fidSHMIDTKO); call erreur(status,.TRUE.,"read Shmidtko sea-floor T,S")
! 
!ALLOCATE(  s_an_5(mNisf)  ) 
!ALLOCATE(  t_an_5(mNisf)  ) 
!ALLOCATE(  s_an_4(mNisf)  ) 
!ALLOCATE(  t_an_4(mNisf)  ) 
! 
!status = NF90_INQ_VARID(fidSHMIDTKO,"s_an_5",s_an_5_ID); call erreur(status,.TRUE.,"inq_s_an_5_ID")
!status = NF90_INQ_VARID(fidSHMIDTKO,"t_an_5",t_an_5_ID); call erreur(status,.TRUE.,"inq_t_an_5_ID")
!status = NF90_INQ_VARID(fidSHMIDTKO,"s_an_4",s_an_4_ID); call erreur(status,.TRUE.,"inq_s_an_4_ID")
!status = NF90_INQ_VARID(fidSHMIDTKO,"t_an_4",t_an_4_ID); call erreur(status,.TRUE.,"inq_t_an_4_ID")
! 
!status = NF90_GET_VAR(fidSHMIDTKO,s_an_5_ID,s_an_5); call erreur(status,.TRUE.,"getvar_s_an_5")
!status = NF90_GET_VAR(fidSHMIDTKO,t_an_5_ID,t_an_5); call erreur(status,.TRUE.,"getvar_t_an_5")
!status = NF90_GET_VAR(fidSHMIDTKO,s_an_4_ID,s_an_4); call erreur(status,.TRUE.,"getvar_s_an_4")
!status = NF90_GET_VAR(fidSHMIDTKO,t_an_4_ID,t_an_4); call erreur(status,.TRUE.,"getvar_t_an_4")
! 
!status = NF90_CLOSE(fidSHMIDTKO); call erreur(status,.TRUE.,"close_file")
!
!---------------------------------------
! Read CMIP5 3d anomalies :
 
write(*,*) 'Reading ', TRIM(file_CMIP5_anom)
 
status = NF90_OPEN(TRIM(file_CMIP5_anom),0,fidCMIP); call erreur(status,.TRUE.,"read CMIP5 anomaly")
 
status = NF90_INQ_DIMID(fidCMIP,"StrLen",dimID_StrLen); call erreur(status,.TRUE.,"inq_dimID_StrLen")
status = NF90_INQ_DIMID(fidCMIP,"depth",dimID_depth); call erreur(status,.TRUE.,"inq_dimID_depth")
status = NF90_INQ_DIMID(fidCMIP,"Nmod",dimID_Nmod); call erreur(status,.TRUE.,"inq_dimID_Nmod")
status = NF90_INQ_DIMID(fidCMIP,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
 
status = NF90_INQUIRE_DIMENSION(fidCMIP,dimID_StrLen,len=mStrLen); call erreur(status,.TRUE.,"inq_dim_StrLen")
status = NF90_INQUIRE_DIMENSION(fidCMIP,dimID_depth,len=mdepth2); call erreur(status,.TRUE.,"inq_dim_depth")
status = NF90_INQUIRE_DIMENSION(fidCMIP,dimID_Nmod,len=mNmod); call erreur(status,.TRUE.,"inq_dim_Nmod")
status = NF90_INQUIRE_DIMENSION(fidCMIP,dimID_Nisf,len=mNisf2); call erreur(status,.TRUE.,"inq_dim_Nisf")
  
if ( mdepth .ne. mdepth2 ) then
  write(*,*) '~!@#$%^* CHECK VERTICAL DIMENSION: MIOSMATCH BETWEEN WOA AND CMIP  >>>> STOP !!'
  stop
endif

if ( mNisf .ne. mNisf2 ) then
  write(*,*) '~!@#$%^*  MISMATCH IN THE NUMBER OF ICE SHELVES BETWEEN THE DIFFERENT FILES  >>>> STOP !'
  stop
endif 

ALLOCATE(  s_anom_3(mNmod,mNisf,mdepth)  ) 
ALLOCATE(  t_anom_3(mNmod,mNisf,mdepth)  ) 
ALLOCATE(  s_anom_2(mNmod,mNisf,mdepth)  ) 
ALLOCATE(  t_anom_2(mNmod,mNisf,mdepth)  ) 
ALLOCATE(  s_anom_1(mNmod,mNisf,mdepth)  ) 
ALLOCATE(  t_anom_1(mNmod,mNisf,mdepth)  ) 
 
status = NF90_INQ_VARID(fidCMIP,"s_anom_3",s_anom_3_ID); call erreur(status,.TRUE.,"inq_s_anom_3_ID")
status = NF90_INQ_VARID(fidCMIP,"t_anom_3",t_anom_3_ID); call erreur(status,.TRUE.,"inq_t_anom_3_ID")
status = NF90_INQ_VARID(fidCMIP,"s_anom_2",s_anom_2_ID); call erreur(status,.TRUE.,"inq_s_anom_2_ID")
status = NF90_INQ_VARID(fidCMIP,"t_anom_2",t_anom_2_ID); call erreur(status,.TRUE.,"inq_t_anom_2_ID")
status = NF90_INQ_VARID(fidCMIP,"s_anom_1",s_anom_1_ID); call erreur(status,.TRUE.,"inq_s_anom_1_ID")
status = NF90_INQ_VARID(fidCMIP,"t_anom_1",t_anom_1_ID); call erreur(status,.TRUE.,"inq_t_anom_1_ID")
status = NF90_INQ_VARID(fidCMIP,"depth",depth_ID); call erreur(status,.TRUE.,"inq_depth_ID")
 
status = NF90_GET_VAR(fidCMIP,s_anom_3_ID,s_anom_3); call erreur(status,.TRUE.,"getvar_s_anom_3")
status = NF90_GET_VAR(fidCMIP,t_anom_3_ID,t_anom_3); call erreur(status,.TRUE.,"getvar_t_anom_3")
status = NF90_GET_VAR(fidCMIP,s_anom_2_ID,s_anom_2); call erreur(status,.TRUE.,"getvar_s_anom_2")
status = NF90_GET_VAR(fidCMIP,t_anom_2_ID,t_anom_2); call erreur(status,.TRUE.,"getvar_t_anom_2")
status = NF90_GET_VAR(fidCMIP,s_anom_1_ID,s_anom_1); call erreur(status,.TRUE.,"getvar_s_anom_1")
status = NF90_GET_VAR(fidCMIP,t_anom_1_ID,t_anom_1); call erreur(status,.TRUE.,"getvar_t_anom_1")
status = NF90_GET_VAR(fidCMIP,depth_ID,depth); call erreur(status,.TRUE.,"getvar_depth")
 
status = NF90_CLOSE(fidCMIP); call erreur(status,.TRUE.,"close_file")

!!---------------------------------------
!! Read CMIP5 bottom anomalies :
! 
!write(*,*) 'Reading ', TRIM(file_CMIP5_anom_bot)
! 
!status = NF90_OPEN(TRIM(file_CMIP5_anom_bot),0,fidCMIPbot); call erreur(status,.TRUE.,"read CMIP5 anomaly")
! 
!status = NF90_INQ_DIMID(fidCMIPbot,"StrLen",dimID_StrLen); call erreur(status,.TRUE.,"inq_dimID_StrLen")
!status = NF90_INQ_DIMID(fidCMIPbot,"Nmod",dimID_Nmod); call erreur(status,.TRUE.,"inq_dimID_Nmod")
!status = NF90_INQ_DIMID(fidCMIPbot,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
! 
!status = NF90_INQUIRE_DIMENSION(fidCMIPbot,dimID_StrLen,len=mStrLen); call erreur(status,.TRUE.,"inq_dim_StrLen")
!status = NF90_INQUIRE_DIMENSION(fidCMIPbot,dimID_Nmod,len=mNmod); call erreur(status,.TRUE.,"inq_dim_Nmod")
!status = NF90_INQUIRE_DIMENSION(fidCMIPbot,dimID_Nisf,len=mNisf2); call erreur(status,.TRUE.,"inq_dim_Nisf")
!  
!if ( mNisf .ne. mNisf2 ) then
!  write(*,*) '~!@#$%^*  MISMATCH IN THE NUMBER OF ICE SHELVES BETWEEN THE DIFFERENT FILES  >>>> STOP !'
!  stop
!endif 
!
!ALLOCATE(  s_anom_5(mNmod,mNisf)  ) 
!ALLOCATE(  t_anom_5(mNmod,mNisf)  ) 
!ALLOCATE(  s_anom_4(mNmod,mNisf)  ) 
!ALLOCATE(  t_anom_4(mNmod,mNisf)  ) 
! 
!status = NF90_INQ_VARID(fidCMIPbot,"s_anom_5",s_anom_5_ID); call erreur(status,.TRUE.,"inq_s_anom_5_ID")
!status = NF90_INQ_VARID(fidCMIPbot,"t_anom_5",t_anom_5_ID); call erreur(status,.TRUE.,"inq_t_anom_5_ID")
!status = NF90_INQ_VARID(fidCMIPbot,"s_anom_4",s_anom_4_ID); call erreur(status,.TRUE.,"inq_s_anom_4_ID")
!status = NF90_INQ_VARID(fidCMIPbot,"t_anom_4",t_anom_4_ID); call erreur(status,.TRUE.,"inq_t_anom_4_ID")
! 
!status = NF90_GET_VAR(fidCMIPbot,s_anom_5_ID,s_anom_5); call erreur(status,.TRUE.,"getvar_s_anom_5")
!status = NF90_GET_VAR(fidCMIPbot,t_anom_5_ID,t_anom_5); call erreur(status,.TRUE.,"getvar_t_anom_5")
!status = NF90_GET_VAR(fidCMIPbot,s_anom_4_ID,s_anom_4); call erreur(status,.TRUE.,"getvar_s_anom_4")
!status = NF90_GET_VAR(fidCMIPbot,t_anom_4_ID,t_anom_4); call erreur(status,.TRUE.,"getvar_t_anom_4")
! 
!status = NF90_CLOSE(fidCMIPbot); call erreur(status,.TRUE.,"close_file")
!
!------------------------------------------------------------------------------
! 2- Define constants and parameters of the simulation :
!------------------------------------------------------------------------------

write(*,*) 'Defining constants and parameters of the simulation...'

!- MISMIP+ values :
yearinsec = 86400.0 * 365.2422 ! Year length in second
rhoi_SI   =  918.0             ! Ice density (kg/m^3)
rhosw_SI  = 1028.0             ! Sea water density (kg/m^3)
rhofw_SI  = 1000.0             ! Fresh water density (kg/m^3)
Lf_SI     =    3.34e5          ! Fusion Latent heat of Ice (J/kg)
cpw_SI    = 3974.0             ! Specific heat of sea water (J/kg/K)
grav_SI   =    9.81            ! Gravity (m/s^2)
lambda1   =   -0.0573          ! Liquidus slope  (K/psu)
lambda2   =    0.0832          ! Liquidus intercept  (K)
lambda3   =   -7.53e-8 * grav_SI * rhosw_SI ! Liquidus pressure coefficient  (K/m)

!- conversions to Elmer units:
rhoi    = rhoi_SI  / (1.0e6*yearinsec**2)
rhosw   = rhosw_SI / (1.0e6*yearinsec**2)
rhofw   = rhofw_SI / (1.0e6*yearinsec**2)
Lf      = Lf_SI   * yearinsec**2
cpw     = cpw_SI  * yearinsec**2
gravity = grav_SI * yearinsec**2 * (-1.0)

meltfac    = rhosw_SI * cpw_SI / ( rhoi_SI * Lf_SI )

!-------------------------------------------------
!- specific to individual melt parameterizations :

! Simple parameterizations :
gammaT_SI     = 3.8 * 1.0e-5            ! exchange velocity (m/s) used for Beckman & Goose 2002; ORIGINAL=1.e-5
facPDC_SI     = 35.6 * 4.997318167e-7   ! factor "K*Kt" (m/s) in Pollard & DeConto 2012; ORIGINAL=1.0*4.997318167e-7

! Plume Emulator :
Cd            = 2.5d-3            ! Drag coefficent
E0_LAZER      = 3.6e-2            ! Entrainment coeff. (no unit) in Lazeroms et al.
GT_LAZER      = 1.1e-3            ! Thermal Stanton number (no unit) in Lazeroms et al.
GTS0_LAZER    = 6.0e-4            ! Effective thermal Stanton number (no unit) in Lazeroms et al.
pp = (/  1.371330075095435d-1, &  ! p0
      &  5.527656234709359d1,  &  ! p1
      & -8.951812433987858d2,  &  ! p2
      &  8.927093637594877d3,  &  ! p3
      & -5.563863123811898d4,  &  ! p4 
      &  2.218596970948727d5,  &  ! p5
      & -5.820015295669482d5,  &  ! p6
      &  1.015475347943186d6,  &  ! p7
      & -1.166290429178556d6,  &  ! p8
      &  8.466870335320488d5,  &  ! p9
      & -3.520598035764990d5,  &  ! p10
      &  6.387953795485420d4   /) ! p11
x0            = 0.56              ! Dimensionless transition melting/freezing
alpha_LAZER   = 3.87d-5
beta_LAZER    = 7.86d-4

gg = 9.81 * yearinsec**2


! Box models :
gammaT_PICO_SI  = 2.0d-5
C_PICO_SI       = 1.0d6                 ! Circulation parameter C (m^6/kg/s) in [0.1;9]*1.e6
alpha_PICO      = 7.5d-5                ! Thermal expansion coefficient for linear EOS (K^-1)
beta_PICO       = 7.7d-4                ! Salinity contraction coefficient for linear EOS (psu^-1)
rhostar_SI_PICO = 1027.51               ! In situ density for linear EOS (kg/m3)
gg_SI           = 9.81d0                ! gravity acceleration, in m.s^-2
ttt_SI          = 5.d10                 ! tuning relationship between C and gammT (inferred from Fig.4 of Reese et al. 2018)

gammaT = gammaT_SI * yearinsec
facPDC = facPDC_SI * yearinsec
gammaT_PICO = gammaT_PICO_SI * yearinsec
C_PICO = C_PICO_SI * (1.0d6*yearinsec**2) * yearinsec
rhostar_PICO = rhostar_SI_PICO / (1.0d6*yearinsec**2)
gg = gg_SI *yearinsec**2
ttt = ttt_SI * (1.0d6*yearinsec**2)

!--------------------
lbd1 = lambda1
lbd2 = lambda2
lbd3 = lambda3

!===================================================================================
!===================================================================================

nn_TS_pres = 3   ! No need for Schmidtko or uncertainty on T here because T is corrected.
                 ! (but still keep three locations because of different CMIP anomalies)
nn_tuning  = 17  !

mtot=nn_K*mNmod*nn_para*nn_tuning*nn_TS_pres

ALLOCATE( T_pres(mdepth), S_pres(mdepth), T_futu(mdepth), S_futu(mdepth) )
ALLOCATE( total_melt_pres(mNisf,mtot), mean_melt_pres(mNisf,mtot) )
ALLOCATE( total_melt_futu(mNisf,mtot), mean_melt_futu(mNisf,mtot) )
ALLOCATE( min_melt_pres(mNisf,mtot), max_melt_pres(mNisf,mtot) )
ALLOCATE( DeltaT_out(mNisf,mtot), Kcoef(mNisf,mtot) ) !, Ccoef(mNisf,mtot) )
ALLOCATE( kstat(mNisf) )
ALLOCATE( index_para(mNisf,mtot), index_CMIP(mNisf,mtot), index_WOA(mNisf,mtot), index_K(mNisf,mtot) )

max_nb_box = 10 ! maximum nb of boxes used in PICO
ALLOCATE( Zbox(max_nb_box,max_nb_box), Abox(max_nb_box,max_nb_box) )

total_melt_pres(:,:) = NF90_FILL_FLOAT
total_melt_futu(:,:) = NF90_FILL_FLOAT
mean_melt_pres(:,:)  = NF90_FILL_FLOAT
mean_melt_futu(:,:)  = NF90_FILL_FLOAT
min_melt_pres(:,:)   = NF90_FILL_FLOAT
max_melt_pres(:,:)   = NF90_FILL_FLOAT
DeltaT_out(:,:)  = NF90_FILL_FLOAT
Kcoef(:,:)  = NF90_FILL_FLOAT
!Ccoef(:,:)  = NF90_FILL_FLOAT
index_para(:,:) = 0
index_CMIP(:,:) = 0
index_WOA (:,:) = 0
index_K   (:,:) = 0

kstat(:)=0

write(*,*) 'MAXVAL(isfmask) = ', MAXVAL(isfmask)

DO kisf=2,mNisf

 IF (       front_min_lon(kisf) .ne. 180.0 .and. front_max_lon(kisf) .ne. -180.0 &      ! to remove non-attributed ice shelves
 &    .and. front_min_lat(kisf) .ne. 90.00 .and. front_max_lat(kisf) .ne. -90.00 ) THEN

   write(*,*) ' '
   write(*,*) '==============================================================================='
   write(*,*) 'kisf = ', kisf
   write(*,*) '                      i in  [ ',imin(kisf),' , ',imax(kisf),' ]'
   write(*,*) '                      j in  [ ',jmin(kisf),' , ',jmax(kisf),' ]'

   mlondim = imax(kisf) - imin(kisf) + 1
   mlatdim = jmax(kisf) - jmin(kisf) + 1

   ALLOCATE(  lon(mlondim)  )
   ALLOCATE(  lat(mlatdim)  )
   ALLOCATE(  bedrock_topography(mlondim,mlatdim)  )
   ALLOCATE(  ice_base_topography(mlondim,mlatdim)  )
   ALLOCATE(  isfmask(mlondim,mlatdim)  )
   ALLOCATE(  IF_mask(mlondim,mlatdim)  )
   ALLOCATE(  GL_mask(mlondim,mlatdim)  )

   ALLOCATE( Melt(mlondim,mlatdim), Pmelt(mlondim,mlatdim,nn_para), Fmelt(mlondim,mlatdim,nn_para) )
   ALLOCATE( Pmelt_out(mlondim,mlatdim,nn_para), Fmelt_out(mlondim,mlatdim,nn_para) )
   ALLOCATE( Pnnnn(mlondim,mlatdim,nn_para), Fnnnn(mlondim,mlatdim,nn_para) )
   ALLOCATE( dGL(mlondim,mlatdim), dIF(mlondim,mlatdim) )
   nnp = 3 ! nb of options for Lazeroms (see TypeL values)
   ALLOCATE( zGL(mlondim,mlatdim,nnp), alpha(mlondim,mlatdim,nnp) )

   Melt(:,:)    = 0.d0
   Pmelt(:,:,:) = 0.q0
   Fmelt(:,:,:) = 0.q0
   Pmelt_out(:,:,:) = 0.e0
   Fmelt_out(:,:,:) = 0.e0
   Pnnnn(:,:,:) = 0
   Fnnnn(:,:,:) = 0
   dGL(:,:)     = 1.d9
   dIF(:,:)     = 1.d9
   zGL(:,:,:)   = 0.d0
   alpha(:,:,:) = 0.d0

   !---------------------------------------                   
   ! Read bedrock, lon, lat
   
   write(*,*) 'Reading ', TRIM(file_bed)
   
   status = NF90_OPEN(TRIM(file_bed),0,fidBED); call erreur(status,.TRUE.,"read bedrock") 
   
   status = NF90_INQ_VARID(fidBED,"lon",lon_ID); call erreur(status,.TRUE.,"inq_lon_ID")
   status = NF90_INQ_VARID(fidBED,"lat",lat_ID); call erreur(status,.TRUE.,"inq_lat_ID")
   status = NF90_INQ_VARID(fidBED,"bedrock_topography",bedrock_topography_ID); call erreur(status,.TRUE.,"inq_bedrock_topography_ID")
   
   status = NF90_GET_VAR(fidBED,lon_ID,lon,start=(/imin(kisf)/),count=(/mlondim/)); call erreur(status,.TRUE.,"getvar_lon")
   status = NF90_GET_VAR(fidBED,lat_ID,lat,start=(/jmin(kisf)/),count=(/mlatdim/)); call erreur(status,.TRUE.,"getvar_lat")
   status = NF90_GET_VAR(fidBED,bedrock_topography_ID,bedrock_topography,start=(/imin(kisf), jmin(kisf)/),count=(/mlondim, mlatdim/)); call erreur(status,.TRUE.,"getvar_bedrock_topography")
   
   status = NF90_CLOSE(fidBED); call erreur(status,.TRUE.,"fin_lecture")     
   
   dlon=abs(lon(3)-lon(2))
   dlat=abs(lat(3)-lat(2))
   write(*,*) '              dlon, dlat = ', dlon, dlat  
 
   !---------------------------------------                   
   ! Read ice draft 
   
   write(*,*) 'Reading ', TRIM(file_dft)
   
   status = NF90_OPEN(TRIM(file_dft),0,fidDFT); call erreur(status,.TRUE.,"read ice draft") 
   
   status = NF90_INQ_VARID(fidDFT,"ice_base_topography",ice_base_topography_ID); call erreur(status,.TRUE.,"inq_ice_base_topography_ID")
   
   status = NF90_GET_VAR(fidDFT,ice_base_topography_ID,ice_base_topography,start=(/imin(kisf), jmin(kisf)/),count=(/mlondim, mlatdim/)); call erreur(status,.TRUE.,"getvar_ice_base_topography")
   
   status = NF90_CLOSE(fidDFT); call erreur(status,.TRUE.,"fin_lecture")     
   
   !---------------------------------------                   
   ! Read ice shelf mask 
   
   write(*,*) 'Reading spatial fields in ', TRIM(file_isf)
   
   status = NF90_OPEN(TRIM(file_isf),0,fidISF); call erreur(status,.TRUE.,"read ice shelf mask") 
   
   status = NF90_INQ_VARID(fidISF,"isfmask",isfmask_ID); call erreur(status,.TRUE.,"inq_isfmask_ID")
   status = NF90_INQ_VARID(fidISF,"IF_mask",IF_mask_ID); call erreur(status,.TRUE.,"inq_IF_mask_ID")
   status = NF90_INQ_VARID(fidISF,"GL_mask",GL_mask_ID); call erreur(status,.TRUE.,"inq_GL_mask_ID")
   
   status = NF90_GET_VAR(fidISF,isfmask_ID,isfmask,start=(/imin(kisf), jmin(kisf)/),count=(/mlondim, mlatdim/)); call erreur(status,.TRUE.,"getvar_isfmask")
   status = NF90_GET_VAR(fidISF,IF_mask_ID,IF_mask,start=(/imin(kisf), jmin(kisf)/),count=(/mlondim, mlatdim/)); call erreur(status,.TRUE.,"getvar_IF_mask")
   status = NF90_GET_VAR(fidISF,GL_mask_ID,GL_mask,start=(/imin(kisf), jmin(kisf)/),count=(/mlondim, mlatdim/)); call erreur(status,.TRUE.,"getvar_GL_mask")
   
   status = NF90_CLOSE(fidISF); call erreur(status,.TRUE.,"fin_lecture")     

   !========== BOXES characteristics for PICO ==========

   write(file_box,381) kisf
   381 FORMAT('BOXES/box_',i3.3,'.nc')
   write(*,*) 'Reading ', TRIM(file_box)
   
   status = NF90_OPEN(TRIM(file_box),0,fidBOX); call erreur(status,.TRUE.,"read file_box") 

   status = NF90_INQ_DIMID(fidBOX,"max_nb_box1",dimID_max_nb_box1); call erreur(status,.TRUE.,"inq_dimID_max_nb_box1")
   status = NF90_INQ_DIMID(fidBOX,"max_nb_box2",dimID_max_nb_box2); call erreur(status,.TRUE.,"inq_dimID_max_nb_box2")
   
   status = NF90_INQUIRE_DIMENSION(fidBOX,dimID_max_nb_box1,len=max_nb_box1); call erreur(status,.TRUE.,"inq_dim_max_nb_box1")
   status = NF90_INQUIRE_DIMENSION(fidBOX,dimID_max_nb_box2,len=max_nb_box2); call erreur(status,.TRUE.,"inq_dim_max_nb_box2")

   if ( max_nb_box1 .ne. max_nb_box .or. max_nb_box2 .ne. max_nb_box ) then
     write(*,*) '~!@#$%^&* max_nb_box not consistent between fortran prgm and netcdf box input file  >>>> STOP !!'
     stop
   endif

   status = NF90_INQ_VARID(fidBOX,"dGL",dGL_ID)   ; call erreur(status,.TRUE.,"inq_dGL_ID")
   status = NF90_INQ_VARID(fidBOX,"dIF",dIF_ID)   ; call erreur(status,.TRUE.,"inq_dIF_ID")
   status = NF90_INQ_VARID(fidBOX,"Abox",Abox_ID) ; call erreur(status,.TRUE.,"inq_Abox_ID")
   status = NF90_INQ_VARID(fidBOX,"Zbox",Zbox_ID) ; call erreur(status,.TRUE.,"inq_Zbox_ID")
   
   status = NF90_GET_VAR(fidBOX,dGL_ID,dGL)   ; call erreur(status,.TRUE.,"getvar_dGL")
   status = NF90_GET_VAR(fidBOX,dIF_ID,dIF)   ; call erreur(status,.TRUE.,"getvar_dIF")
   status = NF90_GET_VAR(fidBOX,Abox_ID,Abox) ; call erreur(status,.TRUE.,"getvar_Abox")
   status = NF90_GET_VAR(fidBOX,Zbox_ID,Zbox) ; call erreur(status,.TRUE.,"getvar_Zbox")
  
   status = NF90_GET_ATT(fidBOX,NF90_GLOBAL,"nD1",nD1); call erreur(status,.TRUE.,"get_att_nD1")
   status = NF90_GET_ATT(fidBOX,NF90_GLOBAL,"nD2",nD2); call erreur(status,.TRUE.,"get_att_nD2")
   status = NF90_GET_ATT(fidBOX,NF90_GLOBAL,"nD3",nD3); call erreur(status,.TRUE.,"get_att_nD3")
 
   status = NF90_CLOSE(fidBOX); call erreur(status,.TRUE.,"fin_lecture")     

   !========== Slope and GL parameters for the plume model ==========

   write(file_plume,382) kisf
   382 FORMAT('PLUME/plume_',i3.3,'.nc')
   write(*,*) 'Reading ', TRIM(file_plume)
   
   status = NF90_OPEN(TRIM(file_plume),0,fidPLUME); call erreur(status,.TRUE.,"read file_plume") 

   status = NF90_INQ_DIMID(fidPLUME,"nnp",dimID_nnp); call erreur(status,.TRUE.,"inq_dimID_nnp")
   
   status = NF90_INQUIRE_DIMENSION(fidPLUME,dimID_nnp,len=nnp2); call erreur(status,.TRUE.,"inq_dim_nnp")

   if ( nnp .ne. nnp2 ) then
     write(*,*) '~!@#$%^&* nb of option for alpha & zGL not consistent between fortran prgm and netcdf box input file  >>>> STOP !!'
     stop
   endif

   status = NF90_INQ_VARID(fidPLUME,"zGL",zGL_ID)   ; call erreur(status,.TRUE.,"inq_zGL_ID")
   status = NF90_INQ_VARID(fidPLUME,"alpha",alpha_ID)   ; call erreur(status,.TRUE.,"inq_alpha_ID")
   
   status = NF90_GET_VAR(fidPLUME,zGL_ID,zGL)   ; call erreur(status,.TRUE.,"getvar_zGL")
   status = NF90_GET_VAR(fidPLUME,alpha_ID,alpha)   ; call erreur(status,.TRUE.,"getvar_alpha")
  
   status = NF90_CLOSE(fidPLUME); call erreur(status,.TRUE.,"fin_lecture")     

   !==========

   DO kk_tuning=1,nn_tuning  !== Total cavity melt used for param calibration ==

      ! Assuming that Rignot's uncertainties represent a normal distribution :
      ! Keeping a symmetric distribution to preserve the mean.
      CALL norm_sample(kk_tuning,nn_tuning,xerr)

      target_melt = melt_isf(kisf) + xerr * err_melt_isf(kisf) ! Target cavity melt in Gt/yr
      write(*,*) '*******************************'
      write(*,*) 'target melt = ', target_melt

      DO kk_TS_pres=1,nn_TS_pres  !== Present-day T,S profiles ==

        if ( kk_TS_pres .eq. 1 ) then
          ! WOA 2013 over the close continental shelf
          T_pres(:) = t_an_1(kisf,:)
          S_pres(:) = s_an_1(kisf,:)
          nn_TS = 1
        elseif ( kk_TS_pres .eq. 2 ) then
          ! WOA 2013 over the "nearby" continental shelf
          T_pres(:) = t_an_2(kisf,:)
          S_pres(:) = s_an_2(kisf,:)
          nn_TS = 2
        elseif ( kk_TS_pres .eq. 3 ) then
          ! WOA 2013 offshore "nearby"
          T_pres(:) = t_an_3(kisf,:)
          S_pres(:) = s_an_3(kisf,:)
          nn_TS = 3
        else
          write(*,*) 'case not anticipated >>>>> stop !!'
          stop
        endif

        DO kk_K=1,nn_K  ! == percentiles of the tuned coeficients ==  (keep no spaces)

          DO kk_para=1,nn_para  !== Melt parameterization ==

            SELECT CASE (kk_para)
 
            !*************************************************************************************
            CASE(1:3)     ! kk_para=1   -> Linear local               : m = k . ( T(z) - Tf(z) )
                          ! kk_para=2   -> Quadratic local            : m = k . ( T(z) - Tf(z) ) . | T(z) - Tf(z) | 
                          ! kk_para=101 -> Linear with bottom Temp    : m = k . ( Tbot - Tf(z) )
                          ! kk_para=3   -> Quadratic with bottom Temp : m = k . ( Tbot - Tf(z) ) . | Tbot - Tf(z) | 

              ll_local = .false. ; if ( kk_para .eq. 1 .or. kk_para .eq. 2 ) ll_local = .true.  ! local / bottom
              ll_quadr = .false. ; if ( kk_para .eq. 2 .or. kk_para .eq. 3 ) ll_quadr = .true.  ! quadratic / linear
 
              gT  =  Kcoef_pct(kk_para,kk_K)

              DeltaT = 0.d0
              incr = maxDeltaT / (Ntun-2)  ! initial increment for the dichotomy on DeltaT
              ll_curr_over = .false.
              ll_curr_undr = .false.
              ll_prev_over = .false.
              ll_prev_undr = .false.

              !=== Present-day melt rates ===
              m_closest = 1.d20
              DeltaT_closest = -99999.99
              DO stun=1,Ntun ! tuning iterations
               Melt(:,:) = 0.d0
               mmm_tot_p = 0.d0
               mmm_avg_p = 0.d0
               area = 0.d0
               do ii=1,mlondim
               do jj=1,mlatdim
                 if ( ll_local ) then
                   zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                 else
                   zz1=front_bot_dep_max(kisf) ! deepest point of the entrance
                 endif
                 zz2=-ice_base_topography(ii,jj) ! ice draft depth
                 if ( isfmask(ii,jj) .eq. kisf ) then
                   CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                   T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                   S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                   Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                   ! Melt in m/yr (meters of ice per year), positive if ice ablation
                   if ( ll_quadr ) then
                     Melt(ii,jj) = - gT * meltfac * meltfac * (T0-Tf+DeltaT) * abs(T0-Tf+DeltaT)
                   else
                     Melt(ii,jj) = - gT * meltfac * (T0-Tf+DeltaT)
                   endif
                   ! total melt (positive if melting) in Gt/yr
                   aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                   mmm_tot_p = mmm_tot_p - Melt(ii,jj) * aaa * 1.d-9 * rhoi_SI / rhofw_SI
                   mmm_avg_p = mmm_avg_p +               aaa * 1.d-9
                 else
                   Melt(ii,jj) = 0.d0
                 endif
               enddo
               enddo
               !
               !write(*,939) kk_para, DeltaT, mmm_tot_p
               ! keep better DeltaT if we never find a satisfactory DeltaT value:
               reldif = ABS( ( mmm_tot_p - target_melt ) / target_melt )
               if ( reldif .lt. m_closest ) then
                 m_closest = reldif
                 DeltaT_closest = DeltaT
                 mmm_avg_p = mmm_tot_p / mmm_avg_p ! average melt rate (in m.w.e/yr)
                 mmm_min_p = MINVAL( -Melt )       ! minimum melt rate in m/yr (meters of ice per year) 
                 mmm_max_p = MAXVAL( -Melt )       ! maximum melt rate in m/yr (meters of ice per year) 
               endif
               !
               IF ( reldif .LT. 1.d-2 ) THEN ! This DeltaT value is satisfactory
                 EXIT
               ELSEIF ( ABS(DeltaT) .ge. maxDeltaT ) THEN  ! Bound DeltaT
                 DeltaT = DeltaT_closest
                 EXIT
               ELSEIF ( stun .eq. Ntun ) THEN ! no possibility of tunning DeltaT to get the correct melt rate
                 DeltaT = -99999.99
               ELSEIF ( mmm_tot_p .gt. target_melt ) THEN
                 ll_curr_over = .true.
                 if ( ll_curr_over .and. ll_prev_undr ) then
                   incr = incr * 1.d-1
                   ll_curr_over = .false.
                   ll_curr_undr = .false.
                   ll_prev_over = .false.
                   ll_prev_undr = .false.
                 endif
                 DeltaT = DeltaT - incr
               ELSEIF ( mmm_tot_p .lt. target_melt ) THEN
                 ll_curr_undr = .true.
                 if ( ll_curr_undr .and. ll_prev_over ) then
                   incr = incr * 1.d-1
                   ll_curr_over = .false.
                   ll_curr_undr = .false.
                   ll_prev_over = .false.
                   ll_prev_undr = .false.
                 endif
                 DeltaT = DeltaT + incr
               ENDIF
               ll_prev_over = ll_curr_over
               ll_prev_undr = ll_curr_undr
              ENDDO ! stun
              if ( DeltaT .gt. -1.d2 ) then
                ! store mean melt pattern for each parametrization :
                !$OMP DO
                do ii=1,mlondim
                do jj=1,mlatdim
                  if( isfmask(ii,jj) .eq. kisf ) then
                    Pmelt(ii,jj,kk_para) = Pmelt(ii,jj,kk_para) - Melt(ii,jj)
                    Pnnnn(ii,jj,kk_para) = Pnnnn(ii,jj,kk_para) + 1
                  endif
                enddo
                enddo
                !$OMP END DO
              endif

              !=== Future melt rates ===
              DO kk_CMIP5_anom = 1,mNmod  !== CMIP5 anomaly ==
       
                if     ( nn_TS .eq. 1 ) then  ! T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_1(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_1(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 2 ) then  ! T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_2(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_2(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 3 ) then  ! T,S taken offshore "nearby"
                   T_futu(:) = T_pres(:) + t_anom_3(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_3(kk_CMIP5_anom,kisf,:)
                !elseif ( nn_TS .eq. 4 ) then  ! sea-floor T,S over the close continental shelf
                !   T_futu(:) = T_pres(:) + t_anom_4(kk_CMIP5_anom,kisf)
                !   S_futu(:) = S_pres(:) + s_anom_4(kk_CMIP5_anom,kisf)
                !elseif ( nn_TS .eq. 5 ) then  ! sea-floor T,S over the "nearby" continental shelf
                !   T_futu(:) = T_pres(:) + t_anom_5(kk_CMIP5_anom,kisf)
                !   S_futu(:) = S_pres(:) + s_anom_5(kk_CMIP5_anom,kisf)
                else
                   write(*,*) '~!@#$% error with T_futu and S_futu calculation >>>>>> STOP !!'
                   stop
                endif

                if ( DeltaT .gt. -1.d2 ) then 
                  Melt(:,:) = 0.d0
                  mmm_tot_f = 0.d0
                  mmm_avg_f = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    if ( ll_local ) then
                      zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                    else
                      zz1=front_bot_dep_max(kisf) ! deepest point of the entrance
                    endif
                    zz2=-ice_base_topography(ii,jj) ! ice draft depth
                    if ( isfmask(ii,jj) .eq. kisf ) then
                      CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                      T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                      S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                      Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                      ! Melt in m/yr (meters of ice per year), positive if ice ablation
                      if ( ll_quadr ) then
                        Melt(ii,jj) = - gT * meltfac * meltfac * (T0-Tf+DeltaT) * abs(T0-Tf+DeltaT)
                      else
                        Melt(ii,jj) = - gT * meltfac * (T0-Tf+DeltaT)
                      endif
                      ! total melt (positive if melting) in Gt/yr
                      aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                      mmm_tot_f = mmm_tot_f - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg_f = mmm_avg_f +               aaa * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) for future and corresponding present state:
                  mmm_avg_f = mmm_tot_f / mmm_avg_f
                  ! saved quantities for offline statistical analysis:
                  kstat(kisf) = kstat(kisf) + 1
                  index_para     (kisf,kstat(kisf)) = kk_para
                  index_CMIP     (kisf,kstat(kisf)) = kk_CMIP5_anom
                  index_WOA      (kisf,kstat(kisf)) = kk_TS_pres
                  index_K        (kisf,kstat(kisf)) = kk_K
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot_p
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot_f
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg_p
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg_f
                  min_melt_pres  (kisf,kstat(kisf)) = mmm_min_p
                  max_melt_pres  (kisf,kstat(kisf)) = mmm_max_p
                  DeltaT_out     (kisf,kstat(kisf)) = DeltaT
                  Kcoef          (kisf,kstat(kisf)) = gT
                endif

              ENDDO  !! DO kk_CMIP5_anom = 1,mNmod
  
            !*****************************************************************************************************************
            CASE(4:6)   ! kk_para=4 -> Quadratic with mixed local/mean TF                  : m = k . ( T(f) - Tf(z) ) . |< T(z) - Tf(z) >|
                        ! kk_para=5 -> Quadratic with mixed local/mean TF and local slope  : m = k . < T(z) - Tf(z) > . |< T(z) - Tf(z) >| . sin(theta)
                        ! kk_para=6 -> Quadratic with mixed local/mean TF and cavity slope : m = k . < T(z) - Tf(z) > . |< T(z) - Tf(z) >| . sin(theta)

              if     ( kk_para .eq.  4 ) then; TypeL='none__'
              elseif ( kk_para .eq.  5 ) then; TypeL='locale'
              elseif ( kk_para .eq.  6 ) then; TypeL='cavity'; endif

              gT  =  gammaT

              if ( dIFzGLmin .gt. 1.e0 ) then
                sn = sin( atan( abs(front_ice_dep_avg-zGLmin) / dIFzGLmin ) ) ! sinus of mean cavity slope
              else
                sn = 0.d0
              endif

              DeltaT = 0.d0
              incr = maxDeltaT / (Ntun-2)  ! initial increment for the dichotomy on DeltaT
              ll_curr_over = .false.
              ll_curr_undr = .false.
              ll_prev_over = .false.
              ll_prev_undr = .false.

              m_closest = 1.d20
              DeltaT_closest = -99999.99
              DO stun=1,Ntun ! tuning iterations
                !== Mean present-day Thermal Forcing ===
                TF_tot = 0.d0
                TF_avg = 0.d0
                do ii=1,mlondim
                do jj=1,mlatdim
                  zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                  zz2=-ice_base_topography(ii,jj) ! ice draft depth
                  if ( isfmask(ii,jj) .eq. kisf ) then
                    CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                    T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                    S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                    Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                    aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                    TF_tot = TF_tot + (T0-Tf+DeltaT) * aaa
                    TF_avg = TF_avg +                  aaa
                  endif
                enddo
                enddo
                TF_avg = TF_tot / TF_avg
                !=== Present-day melt rates === 
                Melt(:,:) = 0.d0
                mmm_tot_p = 0.d0
                mmm_avg_p = 0.d0
                do ii=1,mlondim
                do jj=1,mlatdim
                   if ( isfmask(ii,jj) .eq. kisf ) then
                     zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                     zz2=-ice_base_topography(ii,jj) ! ice draft depth
                     CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                     T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                     S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                     Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                     ! Melt in m/yr (meters of ice per year), positive if ice ablation
                     if     ( TypeL .eq. 'none__' ) then
                       Melt(ii,jj) = - gT * meltfac * meltfac * (T0-Tf+DeltaT) * abs(TF_avg)
                     elseif ( TypeL .eq. 'locale' ) then
                       Melt(ii,jj) = - gT * meltfac * meltfac * (T0-Tf+DeltaT) * abs(TF_avg) * sin(alpha(ii,jj,3)) ! local slope
                     elseif ( TypeL .eq. 'cavity' ) then
                       Melt(ii,jj) = - gT * meltfac * meltfac * (T0-Tf+DeltaT) * abs(TF_avg) * sn  ! cavity slope between GL and front
                     endif
                     ! total melt (positive if melting) in Gt/yr
                     aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                     mmm_tot_p = mmm_tot_p - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                     mmm_avg_p = mmm_avg_p +               aaa * 1.d-9
                   endif
                enddo
                enddo
                !
                !write(*,939) kk_para, DeltaT, mmm_tot_p
                ! keep better DeltaT if we never find a satisfactory DeltaT value:
                reldif = ABS( ( mmm_tot_p - target_melt ) / target_melt )
                if ( reldif .lt. m_closest ) then
                  m_closest = reldif
                  DeltaT_closest = DeltaT
                  mmm_avg_p = mmm_tot_p / mmm_avg_p ! average melt rate (in m.w.e/yr)
                  mmm_min_p = MINVAL( -Melt )       ! minimum melt rate in m/yr (meters of ice per year) 
                  mmm_max_p = MAXVAL( -Melt )       ! maximum melt rate in m/yr (meters of ice per year) 
                endif
                !
                IF ( reldif .LT. 1.d-2 ) THEN ! This DeltaT value is satisfactory
                   EXIT
                ELSEIF ( ABS(DeltaT) .ge. maxDeltaT ) THEN  ! Bound DeltaT
                   DeltaT = DeltaT_closest
                   EXIT
                ELSEIF ( stun .eq. Ntun ) THEN ! no possibility of tunning DeltaT to get the correct melt rate
                   DeltaT = -99999.99
                ELSEIF ( mmm_tot_p .gt. target_melt ) THEN
                   ll_curr_over = .true.
                   if ( ll_curr_over .and. ll_prev_undr ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT - incr
                ELSEIF ( mmm_tot_p .lt. target_melt ) THEN
                   ll_curr_undr = .true.
                   if ( ll_curr_undr .and. ll_prev_over ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT + incr
                ENDIF
                ll_prev_over = ll_curr_over
                ll_prev_undr = ll_curr_undr
              ENDDO ! stun
              !
              if ( DeltaT .gt. -1.d2 ) then
                ! store mean melt pattern for each parametrization :
                !$OMP DO
                do ii=1,mlondim
                do jj=1,mlatdim
                  if( isfmask(ii,jj) .eq. kisf ) then
                    Pmelt(ii,jj,kk_para) = Pmelt(ii,jj,kk_para) - Melt(ii,jj)
                    Pnnnn(ii,jj,kk_para) = Pnnnn(ii,jj,kk_para) + 1  
                  endif
                enddo
                enddo
                !$OMP END DO
              endif

              !=== Future melt rates ===
              DO kk_CMIP5_anom = 1,mNmod  !== CMIP5 anomaly ==
      
                if     ( nn_TS .eq. 1 ) then  ! T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_1(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_1(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 2 ) then  ! T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_2(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_2(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 3 ) then  ! T,S taken offshore "nearby"
                   T_futu(:) = T_pres(:) + t_anom_3(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_3(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 4 ) then  ! sea-floor T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_4(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_4(kk_CMIP5_anom,kisf)
                elseif ( nn_TS .eq. 5 ) then  ! sea-floor T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_5(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_5(kk_CMIP5_anom,kisf)
                else
                   write(*,*) '~!@#$% error with T_futu and S_futu calculation >>>>>> STOP !!'
                   stop
                endif
 
                !-- Mean future Thermal Forcing
                TF_tot = 0.d0
                TF_avg = 0.d0
                do ii=1,mlondim
                do jj=1,mlatdim
                  zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                  zz2=-ice_base_topography(ii,jj) ! ice draft depth
                  if ( isfmask(ii,jj) .eq. kisf ) then
                    CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                    T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                    S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                    Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                    aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                    TF_tot = TF_tot + (T0-Tf+DeltaT) * aaa
                    TF_avg = TF_avg +                  aaa
                  endif
                enddo
                enddo
                TF_avg = TF_tot / TF_avg

                !-- Future melt rates
                if ( DeltaT .gt. -1.d2 ) then 
                  Melt(:,:) = 0.d0
                  mmm_tot_f = 0.d0
                  mmm_avg_f = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    if ( isfmask(ii,jj) .eq. kisf ) then
                      zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                      zz2=-ice_base_topography(ii,jj) ! ice draft depth
                      CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                      T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                      S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                      Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                      ! Melt in m/yr (meters of ice per year), positive if ice ablation
                      if     ( TypeL .eq. 'none__' ) then
                        Melt(ii,jj) = - gT * meltfac * meltfac * (T0-Tf+DeltaT) * abs(TF_avg)
                      elseif ( TypeL .eq. 'locale' ) then
                        Melt(ii,jj) = - gT * meltfac * meltfac * (T0-Tf+DeltaT) * abs(TF_avg) * sin(alpha(ii,jj,3)) ! local slope
                      elseif ( TypeL .eq. 'cavity' ) then
                        Melt(ii,jj) = - gT * meltfac * meltfac * (T0-Tf+DeltaT) * abs(TF_avg) * sn  ! cavity slope between GL and front
                      endif
                      ! total melt (positive if melting) in Gt/yr
                      aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                      mmm_tot_f = mmm_tot_f - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg_f = mmm_avg_f +               aaa * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg_f = mmm_tot_f / mmm_avg_f
                  ! saved quantities for offline statistical analysis:
                  kstat(kisf) = kstat(kisf) + 1
                  index_para     (kisf,kstat(kisf)) = kk_para
                  index_CMIP     (kisf,kstat(kisf)) = kk_CMIP5_anom
                  index_WOA      (kisf,kstat(kisf)) = kk_TS_pres
                  index_K        (kisf,kstat(kisf)) = kk_K 
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot_p
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot_f
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg_p
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg_f
                  min_melt_pres  (kisf,kstat(kisf)) = mmm_min_p
                  max_melt_pres  (kisf,kstat(kisf)) = mmm_max_p
                  DeltaT_out     (kisf,kstat(kisf)) = DeltaT
                  Kcoef          (kisf,kstat(kisf)) = gT
                endif

              ENDDO  !! DO kk_CMIP5_anom = 1,mNmod
 
            !*****************************************************************************************************************
            CASE (7:8) !== ORIGINAL BOX MODEL WITH UNIFORM MELT IN EACH BOX == 
                          ! kk_para= 7 -> BOX MODEL  2-box forced by T_bot (Olbers & Hellmer 2010; Reese et al. 2018)
                          ! kk_para= 8 -> BOX MODEL 10-box forced by T_bot
     
              gT  =  Kcoef_pct(kk_para,kk_K)   ! Effective Exchange Velocity PICO

              CC      = C_PICO        ! Circulation Parameter PICO
              alphap  = alpha_PICO    ! Thermal expansion coeff PICO
              beta    = beta_PICO     ! Salinity contraction coeff PICO
              rhostar = rhostar_PICO  ! EOS ref Density PICO
              nD      = 0             ! Number Boxes PICO
              if     ( kk_para .eq. 7 ) then; nD = nD1
              elseif ( kk_para .eq. 8 ) then; nD = nD3; endif

              ALLOCATE( Tbox(nD), Sbox(nD), Mbox(nD) )

              DeltaT = 0.d0
              incr = maxDeltaT / (Ntun-2)  ! initial increment for the dichotomy on DeltaT
              ll_curr_over = .false.
              ll_curr_undr = .false.
              ll_prev_over = .false.
              ll_prev_undr = .false.

              !=== Present-day Melt rates ===
              !- "Ambiant" present-day temperature and salinity
              zzz=front_bot_dep_max(kisf)   ! deepest entrence depth
              CALL find_z(mdepth,depth,zzz,kkinf,kksup,aainf,aasup)
              T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
              S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
              
              m_closest = 1.d20
              DeltaT_closest = -99999.99
              DO stun=1,Ntun  ! tuning iterations
                !- Temerature and salinity in Box #1 :
                T0c = T0 + DeltaT
                Tstar = lbd1*S0 + lbd2 + lbd3*Zbox(1,nD) - T0c  !NB: Tstar should be < 0
                g1 = Abox(1,nD) * gT
                tmp1 = g1 / (CC*rhostar*(beta*S0*meltfac-alphap))
                xbox = - 0.5*tmp1 + sqrt( max( (0.5*tmp1)**2 - tmp1*Tstar, 0.d0 ))
                Tbox(1) = T0c - xbox
                Sbox(1) = S0 - xbox*S0*meltfac
                qqq = CC*rhostar*(beta*(S0-Sbox(1))-alphap*(T0c-Tbox(1)))
                Mbox(1) = - gT * meltfac * ( lbd1*Sbox(1) + lbd2 + lbd3*Zbox(1,nD) - Tbox(1) )
                !
                !- Temperature and salinity in possible other boxes :
                DO kk=2,nD
                  Tstar = lbd1*Sbox(kk-1) + lbd2 + lbd3*Zbox(kk,nD) - Tbox(kk-1)
                  g1  = Abox(kk,nD) * gT
                  g2  = g1 * meltfac
                  xbox = - g1 * Tstar / ( qqq + g1 - g2*lbd1*Sbox(kk-1) )
                  Tbox(kk) = Tbox(kk-1) - xbox
                  Sbox(kk) = Sbox(kk-1) - xbox*Sbox(kk-1)*meltfac
                  Mbox(kk) = - gT * meltfac * ( lbd1*Sbox(kk) + lbd2 + lbd3*Zbox(kk,nD) - Tbox(kk) ) 
                ENDDO
                !  
                !- Attribute melt at each node to a box value :
                mmm_tot_p = 0.d0
                mmm_avg_p = 0.d0
                do ii=1,mlondim
                do jj=1,mlatdim
                  zzz=-ice_base_topography(ii,jj)
                  if ( isfmask(ii,jj) .eq. kisf ) then
                    rr = dGL(ii,jj) / ( dGL(ii,jj) + dIF(ii,jj) )
                    do kk=1,nD
                      if ( rr .ge. 1.0-sqrt(1.0*(nD-kk+1)/nD) .and. rr .le. 1.0-sqrt(1.0*(nD-kk)/nD) ) then
                        Melt(ii,jj) = - Mbox(kk)
                      endif
                    enddo
                    ! total melt (positive if melting) in Gt/yr
                    aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                    mmm_tot_p = mmm_tot_p - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                    mmm_avg_p = mmm_avg_p +               aaa * 1.d-9
                  endif
                enddo
                enddo
                !
                !write(*,939) kk_para, DeltaT, mmm_tot_p
                ! keep better DeltaT if we never find a satisfactory DeltaT value:
                reldif = ABS( ( mmm_tot_p - target_melt ) / target_melt )
                if ( reldif .lt. m_closest ) then
                  m_closest = reldif
                  DeltaT_closest = DeltaT
                  mmm_avg_p = mmm_tot_p / mmm_avg_p ! average melt rate (in m.w.e/yr)
                  mmm_min_p = MINVAL( -Melt )       ! minimum melt rate in m/yr (meters of ice per year) 
                  mmm_max_p = MAXVAL( -Melt )       ! maximum melt rate in m/yr (meters of ice per year) 
                endif
                !
                IF ( reldif .LT. 1.d-2 ) THEN ! This DeltaT value is satisfactory
                   EXIT
                ELSEIF ( ABS(DeltaT) .ge. maxDeltaT ) THEN  ! Bound DeltaT
                   DeltaT = DeltaT_closest
                   EXIT
                ELSEIF ( stun .eq. Ntun ) THEN ! no possibility of tunning DeltaT to get the correct melt rate
                   DeltaT = -99999.99
                ELSEIF ( mmm_tot_p .gt. target_melt ) THEN
                   ll_curr_over = .true.
                   if ( ll_curr_over .and. ll_prev_undr ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT - incr
                ELSEIF ( mmm_tot_p .lt. target_melt ) THEN
                   ll_curr_undr = .true.
                   if ( ll_curr_undr .and. ll_prev_over ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT + incr
                ENDIF
                ll_prev_over = ll_curr_over
                ll_prev_undr = ll_curr_undr
              ENDDO ! stun
              !
              if ( DeltaT .gt. -1.d2 .and. sum(Abox(:,nD)) .gt. 0.0 ) then
                ! store mean melt pattern for each parametrization :
                !$OMP DO
                do ii=1,mlondim
                do jj=1,mlatdim
                  if( isfmask(ii,jj) .eq. kisf ) then
                    Pmelt(ii,jj,kk_para) = Pmelt(ii,jj,kk_para) - Melt(ii,jj)
                    Pnnnn(ii,jj,kk_para) = Pnnnn(ii,jj,kk_para) + 1  
                  endif
                enddo
                enddo
                !$OMP END DO
              endif
             
              !=== Future melt rates ===
              DO kk_CMIP5_anom = 1,mNmod  !== CMIP5 anomaly ==
      
                if     ( nn_TS .eq. 1 ) then  ! T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_1(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_1(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 2 ) then  ! T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_2(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_2(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 3 ) then  ! T,S taken offshore "nearby"
                   T_futu(:) = T_pres(:) + t_anom_3(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_3(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 4 ) then  ! sea-floor T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_4(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_4(kk_CMIP5_anom,kisf)
                elseif ( nn_TS .eq. 5 ) then  ! sea-floor T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_5(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_5(kk_CMIP5_anom,kisf)
                else
                   write(*,*) '~!@#$% error with T_futu and S_futu calculation >>>>>> STOP !!'
                   stop
                endif
 
                !-- Future melt rates
                if ( DeltaT .gt. -1.d2 .and. sum(Abox(:,nD)) .gt. 0.0 ) then
                  !- "Ambiant" future temperature and salinity
                  zzz=front_bot_dep_max(kisf)   ! deepest entrence depth
                  CALL find_z(mdepth,depth,zzz,kkinf,kksup,aainf,aasup)
                  T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                  S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                  !
                  !- Temerature and salinity in Box #1 :
                  T0c = T0 + DeltaT
                  Tstar = lbd1*S0 + lbd2 + lbd3*Zbox(1,nD) - T0c  !NB: Tstar should be < 0
                  g1 = Abox(1,nD) * gT 
                  tmp1 = g1 / (CC*rhostar*(beta*S0*meltfac-alphap))
                  xbox = - 0.5*tmp1 + sqrt( max( (0.5*tmp1)**2 - tmp1*Tstar, 0.d0 ))
                  Tbox(1) = T0c - xbox
                  Sbox(1) = S0 - xbox*S0*meltfac
                  qqq = CC*rhostar*(beta*(S0-Sbox(1))-alphap*(T0c-Tbox(1)))
                  Mbox(1) = - gT * meltfac * ( lbd1*Sbox(1) + lbd2 + lbd3*Zbox(1,nD) - Tbox(1) )
                  !
                  !- Temperature and salinity in possible other boxes :
                  DO kk=2,nD
                    Tstar = lbd1*Sbox(kk-1) + lbd2 + lbd3*Zbox(kk,nD) - Tbox(kk-1)
                    g1  = Abox(kk,nD) * gT
                    g2  = g1 * meltfac
                    xbox = - g1 * Tstar / ( qqq + g1 - g2*lbd1*Sbox(kk-1) )
                    Tbox(kk) = Tbox(kk-1) - xbox
                    Sbox(kk) = Sbox(kk-1) - xbox*Sbox(kk-1)*meltfac
                    Mbox(kk) = - gT * meltfac * ( lbd1*Sbox(kk) + lbd2 + lbd3*Zbox(kk,nD) - Tbox(kk) ) 
                  ENDDO
                  !  
                  !- Attribute melt at each node to a box value :
                  mmm_tot_f = 0.d0
                  mmm_avg_f = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    if ( isfmask(ii,jj) .eq. kisf ) then
                      rr = dGL(ii,jj) / ( dGL(ii,jj) + dIF(ii,jj) )
                      do kk=1,nD
                        if ( rr .ge. 1.0-sqrt(1.0*(nD-kk+1)/nD) .and. rr .le. 1.0-sqrt(1.0*(nD-kk)/nD) ) then
                          Melt(ii,jj) = - Mbox(kk)
                        endif
                      enddo
                      ! total melt (positive if melting) in Gt/yr
                      aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                      mmm_tot_f = mmm_tot_f - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg_f = mmm_avg_f +               aaa * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg_f = mmm_tot_f / mmm_avg_f
                  ! saved quantities for offline statistical analysis:
                  kstat(kisf) = kstat(kisf) + 1
                  index_para     (kisf,kstat(kisf)) = kk_para
                  index_CMIP     (kisf,kstat(kisf)) = kk_CMIP5_anom
                  index_WOA      (kisf,kstat(kisf)) = kk_TS_pres
                  index_K        (kisf,kstat(kisf)) = kk_K 
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot_p
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot_f
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg_p
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg_f
                  min_melt_pres  (kisf,kstat(kisf)) = mmm_min_p
                  max_melt_pres  (kisf,kstat(kisf)) = mmm_max_p
                  DeltaT_out     (kisf,kstat(kisf)) = DeltaT
                  Kcoef          (kisf,kstat(kisf)) = gT
                  !Ccoef          (kisf,kstat(kisf)) = CC
                endif

              ENDDO  !! DO kk_CMIP5_anom = 1,mNmod

              DEALLOCATE( Tbox, Sbox, Mbox )

            !*****************************************************************************************************************
            CASE (9:12) !== BOX MODEL AS IMPLEMENTED IN PISM (section 2.4 of Reese et al. 2014) == 
                          ! kk_para= 9 -> BOX MODEL  5-box forced by T_bot, C = 1.e6
                          ! kk_para=10 -> BOX MODEL 10-box forced by T_bot, C = 1.e6
                          ! kk_para=11 -> BOX MODEL  5-box forced by T_bot, C = a . gammaT
                          ! kk_para=12 -> BOX MODEL 10-box forced by T_bot, C = a . gammaT

              gT  =  Kcoef_pct(kk_para,kk_K)   ! Effective Exchange Velocity PICO

              if ( kk_para .eq.  9 .or. kk_para .eq. 10 ) then
                CC = C_PICO   ! tuning with constant C value (best circulation parameter in Reese et al. 2018)
              else
                CC = ttt * gT ! tuning by varying C and gammaT together
              endif

              alphap  = alpha_PICO    ! Thermal expansion coeff PICO
              beta    = beta_PICO     ! Salinity contraction coeff PICO
              rhostar = rhostar_PICO  ! EOS ref Density PICO
              nD      = 0             ! Number Boxes PICO
              if     ( kk_para .eq.  9 .or. kk_para .eq. 11 ) then; nD = nD2
              elseif ( kk_para .eq. 10 .or. kk_para .eq. 12 ) then; nD = nD3; endif

              ALLOCATE( Tbox(nD), Sbox(nD) )

              DeltaT = 0.d0
              incr = maxDeltaT / (Ntun-2)  ! initial increment for the dichotomy on DeltaT
              ll_curr_over = .false.
              ll_curr_undr = .false.
              ll_prev_over = .false.
              ll_prev_undr = .false.

              !=== Present-day Melt rates ===
              !- "Ambiant" present-day temperature and salinity
              zzz=front_bot_dep_max(kisf)   ! deepest entrence depth
              CALL find_z(mdepth,depth,zzz,kkinf,kksup,aainf,aasup)
              T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
              S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
              m_closest = 1.d20
              DeltaT_closest = -99999.99
              DO stun=1,Ntun  ! tuning iterations
                mmm_tot_p = 0.d0
                mmm_avg_p = 0.d0
                Tbox(:)=0.d0 ; Sbox(:)=0.d0 ; qqq=0.d0
                !- Temerature and salinity in Box #1 :
                T0c = T0 + DeltaT
                do ii=1,mlondim
                do jj=1,mlatdim
                  rr = dGL(ii,jj) / ( dGL(ii,jj) + dIF(ii,jj) )
                  if ( isfmask(ii,jj) .eq. kisf .and. rr .le. 1.0-sqrt(1.0*(nD-1)/nD) ) then
                    zzz = -ice_base_topography(ii,jj)
                    aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                    Tstar = lbd1*S0 + lbd2 + lbd3*zzz - T0c  !NB: Tstar should be < 0
                    !g1 = gT * aaa
                    g1 = Abox(1,nD) * gT
                    tmp1 = g1 / (CC*rhostar*(beta*S0*meltfac-alphap))
                    xbox = - 0.5*tmp1 + sqrt( max( (0.5*tmp1)**2 - tmp1*Tstar, 0.d0 ))
                    TT = T0c - xbox
                    SS = S0 - xbox*S0*meltfac
                    Tbox(1) = Tbox(1) + TT * aaa
                    Sbox(1) = Sbox(1) + SS * aaa
                    qqq = qqq + CC*rhostar*(beta*(S0-SS)-alphap*(T0c-TT)) * aaa
                    Melt(ii,jj) = gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )
                    ! total melt (positive if melting) in Gt/yr
                    mmm_tot_p = mmm_tot_p - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                    mmm_avg_p = mmm_avg_p +               aaa * 1.d-9
                  endif
                enddo
                enddo
                Tbox(1) = Tbox(1) / Abox(1,nD)
                Sbox(1) = Sbox(1) / Abox(1,nD)
                qqq = qqq / Abox(1,nD)
                !
                !- Temperature and salinity in possible other boxes :
                DO kk=2,nD
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    rr = dGL(ii,jj) / ( dGL(ii,jj) + dIF(ii,jj) )
                    if ( isfmask(ii,jj) .eq. kisf .and. rr .ge. 1.0-sqrt(1.0*(nD-kk+1)/nD) .and. rr .le. 1.0-sqrt(1.0*(nD-kk)/nD) ) then
                      zzz = -ice_base_topography(ii,jj)
                      aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                      Tstar = lbd1*Sbox(kk-1) + lbd2 + lbd3*zzz - Tbox(kk-1)
                      !g1  = gT * aaa
                      g1  = Abox(kk,nD) * gT
                      g2  = g1 * meltfac
                      xbox = - g1 * Tstar / ( qqq + g1 - g2*lbd1*Sbox(kk-1) )
                      TT = Tbox(kk-1) - xbox
                      SS = Sbox(kk-1) - xbox*Sbox(kk-1)*meltfac
                      Tbox(kk) =  Tbox(kk) + TT * aaa
                      Sbox(kk) =  Sbox(kk) + SS * aaa
                      Melt(ii,jj) = gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )
                      ! total melt (positive if melting) in Gt/yr
                      mmm_tot_p = mmm_tot_p - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg_p = mmm_avg_p +               aaa * 1.d-9
                    endif
                  enddo
                  enddo
                  Tbox(kk) = Tbox(kk) / Abox(kk,nD)
                  Sbox(kk) = Sbox(kk) / Abox(kk,nD)
                ENDDO
                !  
                !write(*,939) kk_para, DeltaT, mmm_tot_p
                ! keep better DeltaT if we never find a satisfactory DeltaT value:
                reldif = ABS( ( mmm_tot_p - target_melt ) / target_melt )
                if ( reldif .lt. m_closest ) then
                  m_closest = reldif
                  DeltaT_closest = DeltaT
                  mmm_avg_p = mmm_tot_p / mmm_avg_p ! average melt rate (in m.w.e/yr)
                  mmm_min_p = MINVAL( -Melt )       ! minimum melt rate in m/yr (meters of ice per year) 
                  mmm_max_p = MAXVAL( -Melt )       ! maximum melt rate in m/yr (meters of ice per year) 
                endif
                !
                IF ( reldif .LT. 1.d-2 ) THEN ! This DeltaT value is satisfactory
                   EXIT
                ELSEIF ( ABS(DeltaT) .ge. maxDeltaT ) THEN  ! Bound DeltaT
                   DeltaT = DeltaT_closest
                   EXIT
                ELSEIF ( stun .eq. Ntun ) THEN ! no possibility of tunning DeltaT to get the correct melt rate
                   DeltaT = -99999.99
                ELSEIF ( mmm_tot_p .gt. target_melt ) THEN
                   ll_curr_over = .true.
                   if ( ll_curr_over .and. ll_prev_undr ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT - incr
                ELSEIF ( mmm_tot_p .lt. target_melt ) THEN
                   ll_curr_undr = .true.
                   if ( ll_curr_undr .and. ll_prev_over ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT + incr
                ENDIF
                ll_prev_over = ll_curr_over
                ll_prev_undr = ll_curr_undr
              ENDDO ! stun
              !
              if ( DeltaT .gt. -1.d2 .and. sum(Abox(:,nD)) .gt. 0.0 ) then
                ! store mean melt pattern for each parametrization :
                !$OMP DO
                do ii=1,mlondim
                do jj=1,mlatdim
                  if( isfmask(ii,jj) .eq. kisf ) then
                    Pmelt(ii,jj,kk_para) = Pmelt(ii,jj,kk_para) - Melt(ii,jj)
                    Pnnnn(ii,jj,kk_para) = Pnnnn(ii,jj,kk_para) + 1  
                  endif
                enddo
                enddo
                !$OMP END DO
              endif
             
              !=== Future melt rates ===
              DO kk_CMIP5_anom = 1,mNmod  !== CMIP5 anomaly ==
      
                if     ( nn_TS .eq. 1 ) then  ! T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_1(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_1(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 2 ) then  ! T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_2(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_2(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 3 ) then  ! T,S taken offshore "nearby"
                   T_futu(:) = T_pres(:) + t_anom_3(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_3(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 4 ) then  ! sea-floor T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_4(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_4(kk_CMIP5_anom,kisf)
                elseif ( nn_TS .eq. 5 ) then  ! sea-floor T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_5(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_5(kk_CMIP5_anom,kisf)
                else
                   write(*,*) '~!@#$% error with T_futu and S_futu calculation >>>>>> STOP !!'
                   stop
                endif
 
                !-- Future melt rates
                if ( DeltaT .gt. -1.d2 .and. sum(Abox(:,nD)) .gt. 0.0 ) then
                  !- "Ambiant" future temperature and salinity
                  zzz=front_bot_dep_max(kisf)   ! deepest entrence depth
                  CALL find_z(mdepth,depth,zzz,kkinf,kksup,aainf,aasup)
                  T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                  S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                  !
                  mmm_tot_f = 0.d0
                  mmm_avg_f = 0.d0
                  Tbox(:)=0.d0 ; Sbox(:)=0.d0 ; qqq=0.d0
                  !- Temerature and salinity in Box #1 :
                  T0c = T0 + DeltaT
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    rr = dGL(ii,jj) / ( dGL(ii,jj) + dIF(ii,jj) )
                    if ( isfmask(ii,jj) .eq. kisf .and. rr .le. 1.0-sqrt(1.0*(nD-1)/nD) ) then
                      zzz = -ice_base_topography(ii,jj)
                      aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                      Tstar = lbd1*S0 + lbd2 + lbd3*zzz - T0c  !NB: Tstar should be < 0
                      !g1 = gT * aaa
                      g1 = Abox(1,nD) * gT
                      tmp1 = g1 / (CC*rhostar*(beta*S0*meltfac-alphap))
                      xbox = - 0.5*tmp1 + sqrt( max( (0.5*tmp1)**2 - tmp1*Tstar, 0.d0 ))
                      TT = T0c - xbox
                      SS = S0 - xbox*S0*meltfac
                      Tbox(1) = Tbox(1) + TT * aaa
                      Sbox(1) = Sbox(1) + SS * aaa
                      qqq = qqq + CC*rhostar*(beta*(S0-SS)-alphap*(T0c-TT)) * aaa
                      Melt(ii,jj) = gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )
                      ! total melt (positive if melting) in Gt/yr
                      mmm_tot_f = mmm_tot_f - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg_f = mmm_avg_f +               aaa * 1.d-9
                    endif
                  enddo
                  enddo
                  Tbox(1) = Tbox(1) / Abox(1,nD)
                  Sbox(1) = Sbox(1) / Abox(1,nD)
                  qqq = qqq / Abox(1,nD)
                  !
                  !- Temperature and salinity in possible other boxes :
                  DO kk=2,nD
                    do ii=1,mlondim
                    do jj=1,mlatdim
                      rr = dGL(ii,jj) / ( dGL(ii,jj) + dIF(ii,jj) )
                      if ( isfmask(ii,jj) .eq. kisf .and. rr .ge. 1.0-sqrt(1.0*(nD-kk+1)/nD) .and. rr .le. 1.0-sqrt(1.0*(nD-kk)/nD) ) then
                        zzz = -ice_base_topography(ii,jj)
                        aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                        Tstar = lbd1*Sbox(kk-1) + lbd2 + lbd3*zzz - Tbox(kk-1)
                        !g1  = gT * aaa
                        g1 = Abox(kk,nD) * gT
                        g2  = g1 * meltfac
                        xbox = - g1 * Tstar / ( qqq + g1 - g2*lbd1*Sbox(kk-1) )
                        TT = Tbox(kk-1) - xbox
                        SS = Sbox(kk-1) - xbox*Sbox(kk-1)*meltfac
                        Tbox(kk) =  Tbox(kk) + TT * aaa
                        Sbox(kk) =  Sbox(kk) + SS * aaa
                        Melt(ii,jj) = gT * meltfac * ( lbd1*SS + lbd2 + lbd3*zzz - TT )
                        ! total melt (positive if melting) in Gt/yr
                        mmm_tot_f = mmm_tot_f - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                        mmm_avg_f = mmm_avg_f +               aaa * 1.d-9
                      endif
                    enddo
                    enddo
                    Tbox(kk) = Tbox(kk) / Abox(kk,nD)
                    Sbox(kk) = Sbox(kk) / Abox(kk,nD)
                  ENDDO ! DO kk=2,nD
                  !  
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg_f = mmm_tot_f / mmm_avg_f
                  !
                  if ( DeltaT .gt. -1.d2 .and. sum(Abox(:,nD)) .gt. 0.0 ) then
                    ! store mean melt pattern for each parametrization :
                    !$OMP DO
                    do ii=1,mlondim
                    do jj=1,mlatdim
                      if( isfmask(ii,jj) .eq. kisf ) then
                        Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                        Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1  
                      endif
                    enddo
                    enddo
                    !$OMP END DO
                  endif
                  ! saved quantities for offline statistical analysis:
                  kstat(kisf) = kstat(kisf) + 1
                  index_para     (kisf,kstat(kisf)) = kk_para
                  index_CMIP     (kisf,kstat(kisf)) = kk_CMIP5_anom
                  index_WOA      (kisf,kstat(kisf)) = kk_TS_pres
                  index_K        (kisf,kstat(kisf)) = kk_K 
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot_p
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot_f
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg_p
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg_f
                  min_melt_pres  (kisf,kstat(kisf)) = mmm_min_p
                  max_melt_pres  (kisf,kstat(kisf)) = mmm_max_p
                  DeltaT_out     (kisf,kstat(kisf)) = DeltaT
                  Kcoef          (kisf,kstat(kisf)) = gT
                  !Ccoef          (kisf,kstat(kisf)) = CC
                endif

              ENDDO  !! DO kk_CMIP5_anom = 1,mNmod

              DEALLOCATE( Tbox, Sbox )

            !*************************************************************************
            CASE(13:15) ! Plume model (Lazeroms et al. 2018)
                        ! Tuning both GamT and GefT

              E0      = E0_LAZER      ! Entrainment Coeff LAZER
              if     ( kk_para .eq. 13 ) then; TypeL='simple'; kkp=1        ! Zgl and Alpha are found between draft point and the central point
              elseif ( kk_para .eq. 14 ) then; TypeL='lazero'; kkp=2        ! original from Lazeroms et al.
              elseif ( kk_para .eq. 15 ) then; TypeL='appenB'; kkp=3; endif ! Appendix B in submitted TCD

              GamT = Kcoef_pct(kk_para,kk_K)       ! Thermal Stanton number LAZER
              GefT = GTS0_LAZER * GamT / GT_LAZER  ! Effective thermal Stanton number LAZER

              DeltaT = 0.d0
              incr = maxDeltaT / (Ntun-2)  ! initial increment for the dichotomy on DeltaT
              ll_curr_over = .false.
              ll_curr_undr = .false.
              ll_prev_over = .false.
              ll_prev_undr = .false.

              ! iterations to find present-day melting and tuning coeff. :
              K = 1.d0 ! old tuning factor (no more used for tuning)
              m_closest = 1.d20
              DeltaT_closest = -99999.99
              DO stun=1,Ntun  ! tuning iterations
                Melt(:,:) = 0.d0
                mmm_tot_p = 0.d0
                mmm_avg_p = 0.d0
                do ii=1,mlondim
                do jj=1,mlatdim
                  if ( isfmask(ii,jj) .eq. kisf ) then
                    zz1=-ice_base_topography(ii,jj)
                    zz2=MIN(zGL(ii,jj,kkp),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                    ! "Ambient" temperature and salinity at the effective grounding line depth :
                    CALL find_z(mdepth,depth,zz2,kkinf,kksup,aainf,aasup) 
                    T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                    S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                    ! Sea water freezing temperature at the effective grounding line depth :
                    Tf = lbd1*S0 + lbd2 + lbd3*zGL(ii,jj,kkp)  
                    ! Effective thermal Stanton number :
                    GTS = GamT * ( 0.545 + 3.5e-5 * (T0-Tf+DeltaT)/(-lbd3) * E0*sin(alpha(ii,jj,kkp))/(GefT+E0*sin(alpha(ii,jj,kkp))) )
                    ! Melt scale :
                    MM = 10. * (T0-Tf+DeltaT)**2 * sqrt( sin(alpha(ii,jj,kkp))/(Cd+E0*sin(alpha(ii,jj,kkp))) ) &
                    &    * sqrt( GTS/(GTS+E0*sin(alpha(ii,jj,kkp))) ) * E0*sin(alpha(ii,jj,kkp))/(GTS+E0*sin(alpha(ii,jj,kkp)))
                    ! Length scale :
                    ll = (T0-Tf+DeltaT)/lbd3 * (x0*GTS+E0*sin(alpha(ii,jj,kkp))) / (x0*(GTS+E0*sin(alpha(ii,jj,kkp))))
                    ! Dimensionless coordinate :
                    X_hat = MAX( 0.d0, MIN( 1.d0, ( zz1 - zGL(ii,jj,kkp) ) / ll ) )
                    ! Dimensionless melt curve :
                    M_hat = 0.d0
                    do kk=1,12
                      M_hat = M_hat + pp(kk) * X_hat**(kk-1)
                    enddo
                    ! Melt rate in m/yr:
                    Melt(ii,jj) = - K * MM * M_hat
                    ! total melt (positive if melting) in Gt/yr
                    aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                    mmm_tot_p = mmm_tot_p - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                    mmm_avg_p = mmm_avg_p +               aaa * 1.d-9
                  endif
                enddo
                enddo
                !
                !write(*,939) kk_para, DeltaT, mmm_tot_p
                ! keep better DeltaT if we never find a satisfactory DeltaT value:
                reldif = ABS( ( mmm_tot_p - target_melt ) / target_melt )
                if ( reldif .lt. m_closest ) then
                  m_closest = reldif
                  DeltaT_closest = DeltaT
                  mmm_avg_p = mmm_tot_p / mmm_avg_p ! average melt rate (in m.w.e/yr)
                  mmm_min_p = MINVAL( -Melt )       ! minimum melt rate in m/yr (meters of ice per year) 
                  mmm_max_p = MAXVAL( -Melt )       ! maximum melt rate in m/yr (meters of ice per year) 
                endif
                !
                IF ( reldif .LT. 1.d-2 ) THEN ! This DeltaT value is satisfactory
                   EXIT
                ELSEIF ( ABS(DeltaT) .ge. maxDeltaT ) THEN  ! Bound DeltaT
                   DeltaT = DeltaT_closest
                   EXIT
                ELSEIF ( stun .eq. Ntun ) THEN ! no possibility of tunning DeltaT to get the correct melt rate
                   DeltaT = -99999.99
                ELSEIF ( mmm_tot_p .gt. target_melt ) THEN
                   ll_curr_over = .true.
                   if ( ll_curr_over .and. ll_prev_undr ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT - incr
                ELSEIF ( mmm_tot_p .lt. target_melt ) THEN
                   ll_curr_undr = .true.
                   if ( ll_curr_undr .and. ll_prev_over ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT + incr
                ENDIF
                ll_prev_over = ll_curr_over
                ll_prev_undr = ll_curr_undr
              ENDDO ! stun
              !
              if ( DeltaT .gt. -1.d2 ) then
                ! store mean melt pattern for each parametrization :
                !$OMP DO
                do ii=1,mlondim
                do jj=1,mlatdim
                  if( isfmask(ii,jj) .eq. kisf ) then
                    Pmelt(ii,jj,kk_para) = Pmelt(ii,jj,kk_para) - Melt(ii,jj)
                    Pnnnn(ii,jj,kk_para) = Pnnnn(ii,jj,kk_para) + 1  
                  endif
                enddo
                enddo
                !$OMP END DO
              endif

              !=== Future melt rates ===
              DO kk_CMIP5_anom = 1,mNmod  !== CMIP5 anomaly ==
      
                if     ( nn_TS .eq. 1 ) then  ! T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_1(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_1(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 2 ) then  ! T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_2(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_2(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 3 ) then  ! T,S taken offshore "nearby"
                   T_futu(:) = T_pres(:) + t_anom_3(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_3(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 4 ) then  ! sea-floor T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_4(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_4(kk_CMIP5_anom,kisf)
                elseif ( nn_TS .eq. 5 ) then  ! sea-floor T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_5(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_5(kk_CMIP5_anom,kisf)
                else
                   write(*,*) '~!@#$% error with T_futu and S_futu calculation >>>>>> STOP !!'
                   stop
                endif

                !-- Future melt rates
                if ( DeltaT .gt. -1.d2 ) then
                  Melt(:,:) = 0.d0
                  mmm_tot_f = 0.d0
                  mmm_avg_f = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    if ( isfmask(ii,jj) .eq. kisf ) then
                      zz1=-ice_base_topography(ii,jj)
                      zz2=MIN(zGL(ii,jj,kkp),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                      ! "Ambient" future temperature and salinity at the effective grounding line depth :
                      CALL find_z(mdepth,depth,zz2,kkinf,kksup,aainf,aasup) 
                      T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                      S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                      ! Sea water freezing temperature at the effective grounding line depth :
                      Tf = lbd1*S0 + lbd2 + lbd3*zGL(ii,jj,kkp)  
                      ! Effective thermal Stanton number :
                      GTS = GamT * ( 0.545 + 3.5e-5 * (T0-Tf+DeltaT)/(-lbd3) * E0*sin(alpha(ii,jj,kkp))/(GefT+E0*sin(alpha(ii,jj,kkp))) )
                      ! Melt scale :
                      MM = 10. * (T0-Tf+DeltaT)**2 * sqrt( sin(alpha(ii,jj,kkp))/(Cd+E0*sin(alpha(ii,jj,kkp))) ) &
                      &    * sqrt( GTS/(GTS+E0*sin(alpha(ii,jj,kkp))) ) * E0*sin(alpha(ii,jj,kkp))/(GTS+E0*sin(alpha(ii,jj,kkp)))
                      ! Length scale :
                      ll = (T0-Tf+DeltaT)/lbd3 * (x0*GTS+E0*sin(alpha(ii,jj,kkp))) / (x0*(GTS+E0*sin(alpha(ii,jj,kkp))))
                      ! Dimensionless coordinate :
                      X_hat = MAX( 0.d0, MIN( 1.d0, ( zz1 - zGL(ii,jj,kkp) ) / ll ) )
                      ! Dimensionless melt curve :
                      M_hat = 0.d0
                      do kk=1,12
                        M_hat = M_hat + pp(kk) * X_hat**(kk-1)
                      enddo
                      ! Melt rate in m/yr:
                      Melt(ii,jj) = - K * MM * M_hat
                      ! total melt (positive if melting) in Gt/yr
                      aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                      mmm_tot_f = mmm_tot_f - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg_f = mmm_avg_f +               aaa * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg_f = mmm_tot_f / mmm_avg_f
                  ! saved quantities for offline statistical analysis:
                  kstat(kisf) = kstat(kisf) + 1
                  index_para     (kisf,kstat(kisf)) = kk_para
                  index_CMIP     (kisf,kstat(kisf)) = kk_CMIP5_anom
                  index_WOA      (kisf,kstat(kisf)) = kk_TS_pres
                  index_K        (kisf,kstat(kisf)) = kk_K 
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot_p
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot_f
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg_p
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg_f
                  min_melt_pres  (kisf,kstat(kisf)) = mmm_min_p
                  max_melt_pres  (kisf,kstat(kisf)) = mmm_max_p
                  DeltaT_out     (kisf,kstat(kisf)) = DeltaT
                  Kcoef          (kisf,kstat(kisf)) = GamT
                endif ! if ( DeltaT .gt. -1.d2 )

              ENDDO  !! DO kk_CMIP5_anom = 1,mNmod

            !*************************************************************************
            CASE(16:18) ! Plume model (Lazeroms et al. 2018)
                        ! Tuning both GamT and GefT

              E0      = E0_LAZER      ! Entrainment Coeff LAZER
              if     ( kk_para .eq. 16 ) then; TypeL='simple'; kkp=1        ! Zgl and Alpha are found between draft point and the central point
              elseif ( kk_para .eq. 17 ) then; TypeL='lazero'; kkp=2        ! original from Lazeroms et al.
              elseif ( kk_para .eq. 18 ) then; TypeL='appenB'; kkp=3; endif ! Appendix B in submitted TCD

              GefT = Kcoef_pct(kk_para,kk_K) 

              DeltaT = 0.d0
              incr = maxDeltaT / (Ntun-2)  ! initial increment for the dichotomy on DeltaT
              ll_curr_over = .false.
              ll_curr_undr = .false.
              ll_prev_over = .false.
              ll_prev_undr = .false.

              ! iterations to find present-day melting and tuning coeff. :
              K = 1.d0 ! old tuning factor (no more used for tuning)
              m_closest = 1.d20
              DeltaT_closest = -99999.99
              DO stun=1,Ntun  ! tuning iterations
                Melt(:,:) = 0.d0
                mmm_tot_p = 0.d0
                mmm_avg_p = 0.d0
                do ii=1,mlondim
                do jj=1,mlatdim
                  if ( isfmask(ii,jj) .eq. kisf ) then
                    zz1=-ice_base_topography(ii,jj)
                    zz2=MIN(zGL(ii,jj,kkp),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                    ! "Ambient" temperature and salinity at the effective grounding line depth :
                    CALL find_z(mdepth,depth,zz2,kkinf,kksup,aainf,aasup) 
                    T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                    S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                    ! Sea water freezing temperature at the effective grounding line depth :
                    Tf = lbd1*S0 + lbd2 + lbd3*zGL(ii,jj,kkp)  
                    ! here X_hat = x_tilda, equ. 28b :
                    ctau = (-lbd1*cpw*GefT*S0)/Lf 
                    X_hat = lbd3 * (zz1 - zGL(ii,jj,kkp)) / ( (T0-Tf+DeltaT) * ( 1.d0 + 6.d-1 * ((E0*sin(alpha(ii,jj,kkp)))/(GefT+ctau+E0*sin(alpha(ii,jj,kkp))))**0.75 ) )
                    X_hat = MAX( 0.d0, MIN( 1.d0, X_hat ) )
                    ! here M_hat = M0(X_tilde), equ. 26 :
                    M_hat = 1.d0/(2.d0*sqrt(2.d0)) * ( 3.d0*(1.d0-X_hat)**(4.d0/3.d0) -1.d0 ) * sqrt(1.d0 - (1.d0 - X_hat)**(4.d0/3.d0))
                    !
                    cr1 = (Lf*alpha_LAZER) / (cpw*GefT*beta_LAZER*S0)
                    MM = sqrt((beta_LAZER*S0*gg)/(-lbd3*(Lf/cpw)**3)) * sqrt((1-cr1*GefT)/(Cd+E0*sin(alpha(ii,jj,kkp)))) * ((GefT*E0*sin(alpha(ii,jj,kkp)))/(GefT+ctau+E0*sin(alpha(ii,jj,kkp))))**1.5 * (T0-Tf+DeltaT)**2
                    ! Melt rate in m/yr:
                    Melt(ii,jj) = - K * MM * M_hat
                    ! total melt (positive if melting) in Gt/yr
                    aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                    mmm_tot_p = mmm_tot_p - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                    mmm_avg_p = mmm_avg_p +               aaa * 1.d-9
                  endif
                enddo
                enddo
                !
                !write(*,939) kk_para, DeltaT, mmm_tot_p
                ! keep better DeltaT if we never find a satisfactory DeltaT value:
                reldif = ABS( ( mmm_tot_p - target_melt ) / target_melt )
                if ( reldif .lt. m_closest ) then
                  m_closest = reldif
                  DeltaT_closest = DeltaT
                  mmm_avg_p = mmm_tot_p / mmm_avg_p ! average melt rate (in m.w.e/yr)
                  mmm_min_p = MINVAL( -Melt )       ! minimum melt rate in m/yr (meters of ice per year) 
                  mmm_max_p = MAXVAL( -Melt )       ! maximum melt rate in m/yr (meters of ice per year) 
                endif
                !
                IF ( reldif .LT. 1.d-2 ) THEN ! This DeltaT value is satisfactory
                   EXIT
                ELSEIF ( ABS(DeltaT) .ge. maxDeltaT ) THEN  ! Bound DeltaT
                   DeltaT = DeltaT_closest
                   EXIT
                ELSEIF ( stun .eq. Ntun ) THEN ! no possibility of tunning DeltaT to get the correct melt rate
                   DeltaT = -99999.99
                ELSEIF ( mmm_tot_p .gt. target_melt ) THEN
                   ll_curr_over = .true.
                   if ( ll_curr_over .and. ll_prev_undr ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT - incr
                ELSEIF ( mmm_tot_p .lt. target_melt ) THEN
                   ll_curr_undr = .true.
                   if ( ll_curr_undr .and. ll_prev_over ) then
                     incr = incr * 1.d-1
                     ll_curr_over = .false.
                     ll_curr_undr = .false.
                     ll_prev_over = .false.
                     ll_prev_undr = .false.
                   endif
                   DeltaT = DeltaT + incr
                ENDIF
                ll_prev_over = ll_curr_over
                ll_prev_undr = ll_curr_undr
              ENDDO ! stun
              !
              if ( DeltaT .gt. -1.d2 ) then
                ! store mean melt pattern for each parametrization :
                !$OMP DO
                do ii=1,mlondim
                do jj=1,mlatdim
                  if( isfmask(ii,jj) .eq. kisf ) then
                    Pmelt(ii,jj,kk_para) = Pmelt(ii,jj,kk_para) - Melt(ii,jj)
                    Pnnnn(ii,jj,kk_para) = Pnnnn(ii,jj,kk_para) + 1  
                  endif
                enddo
                enddo
                !$OMP END DO
              endif

              !=== Future melt rates ===
              DO kk_CMIP5_anom = 1,mNmod  !== CMIP5 anomaly ==
      
                if     ( nn_TS .eq. 1 ) then  ! T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_1(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_1(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 2 ) then  ! T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_2(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_2(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 3 ) then  ! T,S taken offshore "nearby"
                   T_futu(:) = T_pres(:) + t_anom_3(kk_CMIP5_anom,kisf,:)
                   S_futu(:) = S_pres(:) + s_anom_3(kk_CMIP5_anom,kisf,:)
                elseif ( nn_TS .eq. 4 ) then  ! sea-floor T,S over the close continental shelf
                   T_futu(:) = T_pres(:) + t_anom_4(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_4(kk_CMIP5_anom,kisf)
                elseif ( nn_TS .eq. 5 ) then  ! sea-floor T,S over the "nearby" continental shelf
                   T_futu(:) = T_pres(:) + t_anom_5(kk_CMIP5_anom,kisf)
                   S_futu(:) = S_pres(:) + s_anom_5(kk_CMIP5_anom,kisf)
                else
                   write(*,*) '~!@#$% error with T_futu and S_futu calculation >>>>>> STOP !!'
                   stop
                endif

                !-- Future melt rates
                if ( DeltaT .gt. -1.d2 ) then
                  Melt(:,:) = 0.d0
                  mmm_tot_f = 0.d0
                  mmm_avg_f = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    if ( isfmask(ii,jj) .eq. kisf ) then
                      zz1=-ice_base_topography(ii,jj)
                      zz2=MIN(zGL(ii,jj,kkp),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                      ! "Ambient" future temperature and salinity at the effective grounding line depth :
                      CALL find_z(mdepth,depth,zz2,kkinf,kksup,aainf,aasup) 
                      T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                      S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                      ! Sea water freezing temperature at the effective grounding line depth :
                      Tf = lbd1*S0 + lbd2 + lbd3*zGL(ii,jj,kkp)  
                      ! here X_hat = x_tilda, equ. 28b :
                      ctau = (-lbd1*cpw*GefT*S0)/Lf 
                      X_hat = lbd3 * (zz1 - zGL(ii,jj,kkp)) / ( (T0-Tf+DeltaT) * ( 1.d0 + 6.d-1 * ((E0*sin(alpha(ii,jj,kkp)))/(GefT+ctau+E0*sin(alpha(ii,jj,kkp))))**0.75 ) )
                      X_hat = MAX( 0.d0, MIN( 1.d0, X_hat ) )
                      ! here M_hat = M0(X_tilde), equ. 26 :
                      M_hat = 1.d0/(2.d0*sqrt(2.d0)) * ( 3.d0*(1.d0-X_hat)**(4.d0/3.d0) -1.d0 ) * sqrt(1.d0 - (1.d0 - X_hat)**(4.d0/3.d0))
                      !
                      cr1 = (Lf*alpha_LAZER) / (cpw*GefT*beta_LAZER*S0)
                      MM = sqrt((beta_LAZER*S0*gg)/(-lbd3*(Lf/cpw)**3)) * sqrt((1-cr1*GefT)/(Cd+E0*sin(alpha(ii,jj,kkp)))) * ((GefT*E0*sin(alpha(ii,jj,kkp)))/(GefT+ctau+E0*sin(alpha(ii,jj,kkp))))**1.5 * (T0-Tf+DeltaT)**2
                      ! Melt rate in m/yr:
                      Melt(ii,jj) = - K * MM * M_hat
                      ! total melt (positive if melting) in Gt/yr
                      aaa = dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
                      mmm_tot_f = mmm_tot_f - Melt(ii,jj) * aaa * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg_f = mmm_avg_f +               aaa * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg_f = mmm_tot_f / mmm_avg_f
                  ! saved quantities for offline statistical analysis:
                  kstat(kisf) = kstat(kisf) + 1
                  index_para     (kisf,kstat(kisf)) = kk_para
                  index_CMIP     (kisf,kstat(kisf)) = kk_CMIP5_anom
                  index_WOA      (kisf,kstat(kisf)) = kk_TS_pres
                  index_K        (kisf,kstat(kisf)) = kk_K 
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot_p
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot_f
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg_p
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg_f
                  min_melt_pres  (kisf,kstat(kisf)) = mmm_min_p
                  max_melt_pres  (kisf,kstat(kisf)) = mmm_max_p
                  DeltaT_out     (kisf,kstat(kisf)) = DeltaT
                  Kcoef          (kisf,kstat(kisf)) = GefT
                endif ! if ( DeltaT .gt. -1.d2 )

              ENDDO  !! DO kk_CMIP5_anom = 1,mNmod


            !*********************************
            CASE DEFAULT

              write(*,*) 'WRONG PARAMETERIZATION CHOICE !!!!!!!!!'
              STOP

            !*********************************
            END SELECT

            !939 FORMAT('  para:',i2.2,' dT=',f6.3,' mm=',f6.1) 

          ENDDO  !! kk_para

        ENDDO  !! kk_K

      ENDDO !! kk_TS_pres

   ENDDO !! kk_tuning

  DEALLOCATE( bedrock_topography, ice_base_topography )
  DEALLOCATE( Melt, IF_mask, GL_mask )
  DEALLOCATE( dGL, dIF, zGL, alpha )

   !$OMP DO
   DO ii=1,mlondim
   DO jj=1,mlatdim
     if ( isfmask(ii,jj) .eq. kisf ) then
      do kk_para=1,nn_para
        if ( Pnnnn(ii,jj,kk_para) .gt. 0 ) then
          Pmelt(ii,jj,kk_para) = Pmelt(ii,jj,kk_para) / Pnnnn(ii,jj,kk_para)
          Pmelt_out(ii,jj,kk_para) = REAL(Pmelt(ii,jj,kk_para),4)
        else
          Pmelt_out(ii,jj,kk_para) = NF90_FILL_FLOAT
        endif
        if ( Fnnnn(ii,jj,kk_para) .gt. 0 ) then
          Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) / Fnnnn(ii,jj,kk_para)
          Fmelt_out(ii,jj,kk_para) = REAL(Fmelt(ii,jj,kk_para),4)
        else
          Fmelt_out(ii,jj,kk_para) = NF90_FILL_FLOAT
        endif
      enddo
     endif
   ENDDO
   ENDDO
   !$OMP END DO

   DEALLOCATE( isfmask, Pmelt, Fmelt )
   DEALLOCATE( Pnnnn, Fnnnn )
 
   !---------------------------------------------------------------------
   !---------------------------------------------------------------------
   ! Writing netcdf file with outputs :                                   
 
   write(file_out,839) kisf
   839 FORMAT('melt_ensemble_DeltaT_',i3.3,'.nc')   
   write(*,*) 'Creating ', TRIM(file_out)
   
   !status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM)
   status = NF90_CREATE(TRIM(file_out),or(NF90_CLOBBER,NF90_64BIT_OFFSET),fidM)
   call erreur(status,.TRUE.,'create output file')
  
   status = NF90_DEF_DIM(fidM,"lon",mlondim,dimID_lon); call erreur(status,.TRUE.,"def_dimID_lon")
   status = NF90_DEF_DIM(fidM,"lat",mlatdim,dimID_lat); call erreur(status,.TRUE.,"def_dimID_lat")
   status = NF90_DEF_DIM(fidM,"Nisf",mNisf,dimID_Nisf); call erreur(status,.TRUE.,"def_dimID_Nisf")
   status = NF90_DEF_DIM(fidM,"mstat",mtot,dimID_mstat); call erreur(status,.TRUE.,"def_dimID_mstat")
   status = NF90_DEF_DIM(fidM,"mpara",nn_para,dimID_mpara); call erreur(status,.TRUE.,"def_dimID_mpara")
   
   status = NF90_DEF_VAR(fidM,"lon",NF90_FLOAT,(/dimID_lon/),lon_ID)
   call erreur(status,.TRUE.,"def_var_lon_ID")
   status = NF90_DEF_VAR(fidM,"lat",NF90_FLOAT,(/dimID_lat/),lat_ID)
   call erreur(status,.TRUE.,"def_var_lat_ID")
   status = NF90_DEF_VAR(fidM,"Pres_Melt",NF90_FLOAT,(/dimID_lon,dimID_lat,dimID_mpara/),Pmelt_ID)
   call erreur(status,.TRUE.,"def_var_Pmelt_ID")
   status = NF90_DEF_VAR(fidM,"Futu_Melt",NF90_FLOAT,(/dimID_lon,dimID_lat,dimID_mpara/),Fmelt_ID)
   call erreur(status,.TRUE.,"def_var_Fmelt_ID")
   status = NF90_DEF_VAR(fidM,"total_melt_pres",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),total_melt_pres_ID)
   call erreur(status,.TRUE.,"def_var_total_melt_pres_ID")
   status = NF90_DEF_VAR(fidM,"mean_melt_pres",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),mean_melt_pres_ID)
   call erreur(status,.TRUE.,"def_var_mean_melt_pres_ID")
   status = NF90_DEF_VAR(fidM,"min_melt_pres",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),min_melt_pres_ID)
   call erreur(status,.TRUE.,"def_var_min_melt_pres_ID")
   status = NF90_DEF_VAR(fidM,"max_melt_pres",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),max_melt_pres_ID)
   call erreur(status,.TRUE.,"def_var_max_melt_pres_ID")
   status = NF90_DEF_VAR(fidM,"total_melt_futu",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),total_melt_futu_ID)
   call erreur(status,.TRUE.,"def_var_total_melt_futu_ID")
   status = NF90_DEF_VAR(fidM,"mean_melt_futu",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),mean_melt_futu_ID)
   call erreur(status,.TRUE.,"def_var_mean_melt_futu_ID")
   status = NF90_DEF_VAR(fidM,"DeltaT",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),DeltaT_ID)
   call erreur(status,.TRUE.,"def_var_DeltaT_ID")
   status = NF90_DEF_VAR(fidM,"Kcoef",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),Kcoef_ID)
   call erreur(status,.TRUE.,"def_var_Kcoef_ID")
   !status = NF90_DEF_VAR(fidM,"Ccoef",NF90_FLOAT,(/dimID_Nisf,dimID_mstat/),Ccoef_ID)
   !call erreur(status,.TRUE.,"def_var_Ccoef_ID")
   status = NF90_DEF_VAR(fidM,"index_para",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_para_ID)
   call erreur(status,.TRUE.,"def_var_index_para_ID")
   status = NF90_DEF_VAR(fidM,"index_CMIP",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_CMIP_ID)
   call erreur(status,.TRUE.,"def_var_index_CMIP_ID")
   status = NF90_DEF_VAR(fidM,"index_WOA",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_WOA_ID)
   call erreur(status,.TRUE.,"def_var_index_WOA_ID")
   status = NF90_DEF_VAR(fidM,"index_K",NF90_BYTE,(/dimID_Nisf,dimID_mstat/),index_K_ID)
   call erreur(status,.TRUE.,"def_var_index_K_ID")
   
   status = NF90_PUT_ATT(fidM,lon_ID,"standard_name","longitude")
   call erreur(status,.TRUE.,"put_att_lon_ID")
   status = NF90_PUT_ATT(fidM,lon_ID,"units","degrees_east")
   call erreur(status,.TRUE.,"put_att_lon_ID")
   status = NF90_PUT_ATT(fidM,lon_ID,"long_name","longitude coordinate")
   call erreur(status,.TRUE.,"put_att_lon_ID")
   status = NF90_PUT_ATT(fidM,lon_ID,"name","lon")
   call erreur(status,.TRUE.,"put_att_lon_ID")
   status = NF90_PUT_ATT(fidM,lat_ID,"standard_name","latitude")
   call erreur(status,.TRUE.,"put_att_lat_ID")
   status = NF90_PUT_ATT(fidM,lat_ID,"units","degrees_north")
   call erreur(status,.TRUE.,"put_att_lat_ID")
   status = NF90_PUT_ATT(fidM,lat_ID,"long_name","latitude coordinate")
   call erreur(status,.TRUE.,"put_att_lat_ID")
   status = NF90_PUT_ATT(fidM,lat_ID,"name","lat")
   call erreur(status,.TRUE.,"put_att_lat_ID")
   status = NF90_PUT_ATT(fidM,Pmelt_ID,"coordinates","lat lon Nisf")
   call erreur(status,.TRUE.,"put_att_Pmelt_ID")
   status = NF90_PUT_ATT(fidM,Pmelt_ID,"units","m/yr")
   call erreur(status,.TRUE.,"put_att_Pmelt_ID")
   status = NF90_PUT_ATT(fidM,Pmelt_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_Pmelt_ID")
   status = NF90_PUT_ATT(fidM,Pmelt_ID,"long_name","Present ice shelf melt rate")
   call erreur(status,.TRUE.,"put_att_Pmelt_ID")
   status = NF90_PUT_ATT(fidM,Pmelt_ID,"title","Present melt rate")
   call erreur(status,.TRUE.,"put_att_Pmelt_ID")
   status = NF90_PUT_ATT(fidM,Fmelt_ID,"coordinates","lat lon Nisf")
   call erreur(status,.TRUE.,"put_att_Fmelt_ID")
   status = NF90_PUT_ATT(fidM,Fmelt_ID,"units","m/yr")
   call erreur(status,.TRUE.,"put_att_Fmelt_ID")
   status = NF90_PUT_ATT(fidM,Fmelt_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_Fmelt_ID")
   status = NF90_PUT_ATT(fidM,Fmelt_ID,"long_name","Future ice shelf melt rate")
   call erreur(status,.TRUE.,"put_att_Fmelt_ID")
   status = NF90_PUT_ATT(fidM,Fmelt_ID,"title","Future melt rate")
   call erreur(status,.TRUE.,"put_att_Fmelt_ID")
   status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"units","Gt/yr")
   call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"long_name","Total present cavity melt flux")
   call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"title","Present Melt Flux")
   call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"units","Gt/yr")
   call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"long_name","Total future cavity melt flux")
   call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"title","Future Melt Flux")
   call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"units","m/yr")
   call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"long_name","Present mean melt rate over the cavity")
   call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"title","Present Melt Rate")
   call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,min_melt_pres_ID,"units","m/yr")
   call erreur(status,.TRUE.,"put_att_min_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,min_melt_pres_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_min_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,min_melt_pres_ID,"long_name","Present minimum melt rate over the cavity")
   call erreur(status,.TRUE.,"put_att_min_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,min_melt_pres_ID,"title","Minimum Melt Rate")
   call erreur(status,.TRUE.,"put_att_min_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,max_melt_pres_ID,"units","m/yr")
   call erreur(status,.TRUE.,"put_att_max_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,max_melt_pres_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_max_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,max_melt_pres_ID,"long_name","Present maximum melt rate over the cavity")
   call erreur(status,.TRUE.,"put_att_max_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,max_melt_pres_ID,"title","Maximum Melt Rate")
   call erreur(status,.TRUE.,"put_att_max_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"units","m/yr")
   call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"long_name","Future mean melt rate over the cavity")
   call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"title","Future Melt Rate")
   call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,DeltaT_ID,"units","K")
   call erreur(status,.TRUE.,"put_att_DeltaT_ID")
   status = NF90_PUT_ATT(fidM,DeltaT_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_DeltaT_ID")
   status = NF90_PUT_ATT(fidM,DeltaT_ID,"long_name","Thermal Forcing correction")
   call erreur(status,.TRUE.,"put_att_DeltaT_ID")
   status = NF90_PUT_ATT(fidM,DeltaT_ID,"title","DeltaT")
   call erreur(status,.TRUE.,"put_att_DeltaT_ID")
   status = NF90_PUT_ATT(fidM,Kcoef_ID,"units","m/s")
   call erreur(status,.TRUE.,"put_att_Kcoef_ID")
   status = NF90_PUT_ATT(fidM,Kcoef_ID,"_FillValue",NF90_FILL_FLOAT)
   call erreur(status,.TRUE.,"put_att_Kcoef_ID")
   status = NF90_PUT_ATT(fidM,Kcoef_ID,"long_name","K coefficient (tuning factor)")
   call erreur(status,.TRUE.,"put_att_Kcoef_ID")
   status = NF90_PUT_ATT(fidM,Kcoef_ID,"title","K coeff")
   call erreur(status,.TRUE.,"put_att_Kcoef_ID")
   !status = NF90_PUT_ATT(fidM,Ccoef_ID,"units","TBD")
   !call erreur(status,.TRUE.,"put_att_Ccoef_ID")
   !status = NF90_PUT_ATT(fidM,Ccoef_ID,"_FillValue",NF90_FILL_FLOAT)
   !call erreur(status,.TRUE.,"put_att_Ccoef_ID")
   !status = NF90_PUT_ATT(fidM,Ccoef_ID,"long_name","C coefficient (PICO's overturning)")
   !call erreur(status,.TRUE.,"put_att_Ccoef_ID")
   !status = NF90_PUT_ATT(fidM,Ccoef_ID,"title","C coeff")
   !call erreur(status,.TRUE.,"put_att_Ccoef_ID")
   status = NF90_PUT_ATT(fidM,index_para_ID,"units","-")
   call erreur(status,.TRUE.,"put_att_index_para_ID")
   status = NF90_PUT_ATT(fidM,index_para_ID,"long_name","index defining the melt parameterization")
   call erreur(status,.TRUE.,"put_att_index_para_ID")
   status = NF90_PUT_ATT(fidM,index_para_ID,"title","index_para")
   call erreur(status,.TRUE.,"put_att_index_para_ID")
   status = NF90_PUT_ATT(fidM,index_CMIP_ID,"units","-")
   call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
   status = NF90_PUT_ATT(fidM,index_CMIP_ID,"long_name","index defining the CMIP5 model used for the climate anomaly")
   call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
   status = NF90_PUT_ATT(fidM,index_CMIP_ID,"title","index_CMIP")
   call erreur(status,.TRUE.,"put_att_index_CMIP_ID")
   status = NF90_PUT_ATT(fidM,index_WOA_ID,"units","-")
   call erreur(status,.TRUE.,"put_att_index_WOA_ID")
   status = NF90_PUT_ATT(fidM,index_WOA_ID,"latg_name","index defining the sampling of WOA2013 range of uncertainty")
   call erreur(status,.TRUE.,"put_att_index_WOA_ID")
   status = NF90_PUT_ATT(fidM,index_WOA_ID,"title","index_WOA")
   call erreur(status,.TRUE.,"put_att_index_WOA_ID")
   status = NF90_PUT_ATT(fidM,index_K_ID,"units","-")
   call erreur(status,.TRUE.,"put_att_index_K_ID")
   status = NF90_PUT_ATT(fidM,index_K_ID,"long_name","index defining the sampling of Rignot (2013)'s range of uncertainty")
   call erreur(status,.TRUE.,"put_att_index_K_ID")
   status = NF90_PUT_ATT(fidM,index_K_ID,"title","index_K")
   call erreur(status,.TRUE.,"put_att_index_K_ID")
   
   status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using calculate_melt_rate_ensemble_DeltaT.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
   status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imin",imin(kisf));  call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
   status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"imax",imax(kisf));  call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
   status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmin",jmin(kisf));  call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
   status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"jmax",jmax(kisf));  call erreur(status,.TRUE.,"put_att_GLOBAL_ID")

   status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
   
   status = NF90_PUT_VAR(fidM,lon_ID,lon); call erreur(status,.TRUE.,"var_lon_ID")
   status = NF90_PUT_VAR(fidM,lat_ID,lat); call erreur(status,.TRUE.,"var_lat_ID")
   status = NF90_PUT_VAR(fidM,Pmelt_ID,Pmelt_out); call erreur(status,.TRUE.,"var_Pmelt_ID")
   status = NF90_PUT_VAR(fidM,Fmelt_ID,Fmelt_out); call erreur(status,.TRUE.,"var_Fmelt_ID")
   status = NF90_PUT_VAR(fidM,total_melt_pres_ID,total_melt_pres); call erreur(status,.TRUE.,"var_total_melt_pres_ID")
   status = NF90_PUT_VAR(fidM,mean_melt_pres_ID,mean_melt_pres); call erreur(status,.TRUE.,"var_mean_melt_pres_ID")
   status = NF90_PUT_VAR(fidM,max_melt_pres_ID,max_melt_pres); call erreur(status,.TRUE.,"var_max_melt_pres_ID")
   status = NF90_PUT_VAR(fidM,min_melt_pres_ID,min_melt_pres); call erreur(status,.TRUE.,"var_min_melt_pres_ID")
   status = NF90_PUT_VAR(fidM,total_melt_futu_ID,total_melt_futu); call erreur(status,.TRUE.,"var_total_melt_futu_ID")
   status = NF90_PUT_VAR(fidM,mean_melt_futu_ID,mean_melt_futu); call erreur(status,.TRUE.,"var_mean_melt_futu_ID")
   status = NF90_PUT_VAR(fidM,DeltaT_ID,DeltaT_out); call erreur(status,.TRUE.,"var_DeltaT_ID")
   status = NF90_PUT_VAR(fidM,Kcoef_ID,Kcoef); call erreur(status,.TRUE.,"var_Kcoef_ID")
   !status = NF90_PUT_VAR(fidM,Ccoef_ID,Ccoef); call erreur(status,.TRUE.,"var_Ccoef_ID")
   status = NF90_PUT_VAR(fidM,index_para_ID,index_para); call erreur(status,.TRUE.,"var_index_para_ID")
   status = NF90_PUT_VAR(fidM,index_CMIP_ID,index_CMIP); call erreur(status,.TRUE.,"var_index_CMIP_ID")
   status = NF90_PUT_VAR(fidM,index_WOA_ID,index_WOA); call erreur(status,.TRUE.,"var_index_WOA_ID")
   status = NF90_PUT_VAR(fidM,index_K_ID,index_K); call erreur(status,.TRUE.,"var_index_K_ID")
   
   status = NF90_CLOSE(fidM); call erreur(status,.TRUE.,"final")         
 
   !============================================
 
   DEALLOCATE( lon, lat ) 
   DEALLOCATE( Pmelt_out, Fmelt_out )

 ENDIF  !! ( melt_isf(kisf) .ne. 0.e0 )

ENDDO !! kisf

end program ensemble_melt

!==================================================================================================
!==================================================================================================

!--------------------------------------------------------------------
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

!--------------------------------------------------------------------
! calculate the determinant from a matrix
REAL(KIND=8) FUNCTION det(x1,y1,x2,y2)

REAL(KIND=8) :: x1,y1,x2,y2

det=x1*y2-x2*y1

return
END

!--------------------------------------------------------------------
!sort elements of a vector in ascending order
!inspired from "sort_shell" taken from "Numerical Recipes in Fortran, Vol 2"
SUBROUTINE sort_shell(arr,x,y,isfd,n)
!USE nrtype
IMPLICIT NONE
REAL(KIND=8), DIMENSION(n), INTENT(INOUT) :: arr,x,y,isfd
INTEGER, INTENT(IN) :: n
!Sorts an array arr into ascending numerical order by Shell’s method (diminishing increment
!sort). arr is replaced on output by its sorted rearrangement.
INTEGER :: i,j,inc
REAL(KIND=8) :: v,v2,v3,v4
inc=1
do !Determine the starting increment.
  inc=3*inc+1
  if (inc .gt. n) exit
end do
do !Loop over the partial sorts.
  inc=inc/3
  do i=inc+1,n !Outer loop of straight insertion.
    v=arr(i)
    v2=x(i)
    v3=y(i)
    v4=isfd(i)
    j=i
    do !Inner loop of straight insertion.
      if (arr(j-inc) .le. v) exit
      arr(j)=arr(j-inc)
      x(j)=x(j-inc)
      y(j)=y(j-inc)
      isfd(j)=isfd(j-inc)
      j=j-inc
      if (j .le. inc) exit
    end do
    arr(j)=v
    x(j)=v2
    y(j)=v3
    isfd(j)=v4
  end do
  if (inc .le. 1) exit
end do
END SUBROUTINE sort_shell


!--------------------------------------------------------------------
! To find the closest vertical levels and interpolation coefficients :
SUBROUTINE find_z(mdep,dep,z_out,kinf,ksup,ainf,asup)
IMPLICIT NONE
INTEGER, INTENT(IN) :: mdep
REAL(KIND=4), DIMENSION(mdep), INTENT(IN)  :: dep
REAL(KIND=8), INTENT(IN)  :: z_out ! z at which T_out and S_out are calculated
REAL(KIND=8), INTENT(OUT) :: ainf, asup
INTEGER, INTENT(OUT)      :: kinf, ksup
INTEGER :: k
REAL(KIND=8) :: ztmp
ztmp=1.e6
kinf=mdep-1
do k=1,mdep-1
  if ( z_out-dep(k) .lt. ztmp .and. z_out .ge. dep(k) ) then
    ztmp = z_out - dep(k)
    kinf = k
  endif
enddo
ksup = kinf + 1
ainf = ( dep(ksup) - z_out ) / ( dep(ksup) - dep(kinf) )
asup = ( z_out - dep(kinf) ) / ( dep(ksup) - dep(kinf) )
END SUBROUTINE find_z


!--------------------------------------------------------------------
! To sample normal disctributions of mean=0 and std=1 correctly :
SUBROUTINE norm_sample(is,ns,xout)
IMPLICIT NONE
INTEGER, INTENT(IN) ::       is, & ! sample id in [1,ns]
&                            ns    ! total number of samples
REAL(KIND=8), INTENT(OUT) :: xout  ! x value to sample
REAL(KIND=8) :: norm, x, deno, dx
INTEGER :: k, n

n=4001
deno=1.d0 / sqrt( 2 * dacos(-1.d0) )

IF ( is .LT. 1 .OR. is .GT. ns ) THEN
  write(*,*) '~!@#$%^* ERROR : WRONG INPUT VALUE GIVEN TO SUBROUTINE norm_sample  >>>> STOP'
  STOP
ENDIF

IF ( MOD(ns,2) .EQ. 0 ) THEN
  write(*,*) '~!@#$%^* ERROR in SUBROUTINE norm_sample : 2nd argt should be an odd number >>>> STOP'
  STOP
ENDIF

dx=20.d0 / (n-1)
IF ( ns .EQ. 1 ) THEN
  xout = 0.d0
ELSEIF ( is .eq. INT(ns/2)+1 ) THEN
  xout = 0.d0
ELSEIF ( is .lt. INT(ns/2)+1 ) THEN
  x=-10.d0
  norm=exp( - (x**2) / 2.d0 ) * deno * dx
  DO k=2,INT(n/2)+1
    x = x + dx
    norm = norm + exp( - (x**2) / 2.d0 ) * deno * dx
    IF ( norm .GE. 0.5d0*is/(INT(ns/2)+1) ) THEN
      xout = x
      EXIT
    ENDIF
  ENDDO
ELSE
  x=10.d0
  norm=1.d0 - exp( - (x**2) / 2.d0 ) * deno * dx
  DO k=2,INT(n/2)+1
    x = x - dx
    norm = norm - exp( - (x**2) / 2.d0 ) * deno * dx
    IF ( norm .LE. 1.d0-0.5d0*(ns-is+1)/(INT(ns/2)+1) ) THEN
      xout = x
      EXIT
    ENDIF
  ENDDO
ENDIF

END SUBROUTINE norm_sample

!-----------------------------
SUBROUTINE dist_sphere(lon1,lon2,lat1,lat2,distance)
IMPLICIT NONE
REAL*8, INTENT(IN) :: lon1,lon2,lat1,lat2
REAL*8, INTENT(OUT) :: distance
REAL*8 :: a, R, d2r

R = 6.371d6
d2r = dacos(-1.d0) / 180.d0

a =   dsin(5.d-1*(lat1-lat2)*d2r) * dsin(5.d-1*(lat1-lat2)*d2r)             &
         &       + dcos(lat2*d2r) * cos(lat1*d2r)                         &
         &       * dsin(5.d-1*(lon1-lon2)*d2r) * dsin(5.d-1*(lon1-lon2)*d2r)

distance = abs( 2.d0 * R * datan2( dsqrt(a), dsqrt(1.-a) ) )

END SUBROUTINE dist_sphere

!-----------------------------
subroutine hpsort_eps_epw (n, ra, ind, eps)
!========================================================================================
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
!                                                                            
! This file is distributed under the terms of the GNU General Public         
! License. See the file `LICENSE' in the root directory of the               
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
!                                                                            
! Adapted from flib/hpsort_eps
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! sort an array ra(1:n) into ascending order using heapsort algorithm,
! and considering two elements being equal if their values differ
! for less than "eps".
! n is input, ra is replaced on output by its sorted rearrangement.
! create an index table (ind) by making an exchange in the index array
! whenever an exchange is made on the sorted data array (ra).
! in case of equal values in the data array (ra) the values in the
! index array (ind) are used to order the entries.
! if on input ind(1)  = 0 then indices are initialized in the routine,
! if on input ind(1) != 0 then indices are assumed to have been
!                initialized before entering the routine and these
!                indices are carried around during the sorting process
!
! no work space needed !
! free us from machine-dependent sorting-routines !
!
! adapted from Numerical Recipes pg. 329 (new edition)
!
implicit none  
!-input/output variables
integer, intent(in)   :: n  
real(kind=4), intent(in)  :: eps
integer :: ind (n)  
real(kind=4) :: ra (n)
!-local variables
integer :: i, ir, j, l, iind  
real(kind=4) :: rra  
!
! initialize index array
IF (ind (1) .eq.0) then  
   DO i = 1, n  
      ind (i) = i  
   ENDDO
ENDIF
! nothing to order
IF (n.lt.2) return  
! initialize indices for hiring and retirement-promotion phase
l = n / 2 + 1  

ir = n  

sorting: do 

  ! still in hiring phase
  IF ( l .gt. 1 ) then  
     l    = l - 1  
     rra  = ra (l)  
     iind = ind (l)  
     ! in retirement-promotion phase.
  ELSE  
     ! clear a space at the end of the array
     rra  = ra (ir)  
     !
     iind = ind (ir)  
     ! retire the top of the heap into it
     ra (ir) = ra (1)  
     !
     ind (ir) = ind (1)  
     ! decrease the size of the corporation
     ir = ir - 1  
     ! done with the last promotion
     IF ( ir .eq. 1 ) then  
        ! the least competent worker at all !
        ra (1)  = rra  
        !
        ind (1) = iind  
        exit sorting  
     ENDIF
  ENDIF
  ! wheter in hiring or promotion phase, we
  i = l  
  ! set up to place rra in its proper level
  j = l + l  
  !
  DO while ( j .le. ir )  
     IF ( j .lt. ir ) then  
        ! compare to better underling
        IF ( hslt( ra (j),  ra (j + 1) ) ) then  
           j = j + 1  
        !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
           ! this means ra(j) == ra(j+1) within tolerance
         !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
        ENDIF
     ENDIF
     ! demote rra
     IF ( hslt( rra, ra (j) ) ) then  
        ra (i) = ra (j)  
        ind (i) = ind (j)  
        i = j  
        j = j + j  
     !else if ( .not. hslt ( ra(j) , rra ) ) then
        !this means rra == ra(j) within tolerance
        ! demote rra
       ! if (iind.lt.ind (j) ) then
       !    ra (i) = ra (j)
       !    ind (i) = ind (j)
       !    i = j
       !    j = j + j
       ! else
           ! set j to terminate do-while loop
       !    j = ir + 1
       ! endif
        ! this is the right place for rra
     ELSE
        ! set j to terminate do-while loop
        j = ir + 1  
     ENDIF
  ENDDO
  ra (i) = rra  
  ind (i) = iind  

END DO sorting    
contains 

!  internal function 
!  compare two real number and return the result

logical function hslt( a, b )
  REAL(kind=4) :: a, b
  IF( abs(a-b) <  eps ) then
    hslt = .false.
  ELSE
    hslt = ( a < b )
  end if
end function hslt

  !
end subroutine hpsort_eps_epw

