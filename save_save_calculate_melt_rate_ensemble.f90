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

INTEGER(KIND=4) :: nD, ii, jj, kk, kki, kkj, cpt, sizev, kkstart, max_nb_box, kkp, dimID_mstat, dimID_mpara, kndx,   &
&                  Pmelt_ID, Fmelt_ID, total_melt_pres_ID, mean_melt_pres_ID, total_melt_futu_ID, mean_melt_futu_ID, &
&                  index_para_ID, index_CMIP_ID, index_WOA_ID, index_tun_ID, nD1, nD2, nD3, im1, ip1, jm1, jp1

INTEGER(KIND=4) :: varXid, varYid, varid, dimid1, dimid2, dimid3, lenX, lenY, lenTime, status1, res, ierr, hwidth, &
&                cptglob, clL, imin_ID, jmin_ID, imax_ID, jmax_ID, cond12, cond23, cond13

REAL(KIND=8) :: x_NC_Init, x_NC_Fin, y_NC_Init, y_NC_Fin, x_NC_Res, y_NC_Res, localInteg, Integ, Integ_Reduced,    &
&                z1_0, z2_0, z1, z2, zzz, Tsrf0, Tbot0, Tsrf, Tbot, T0, Ssrf0, Sbot0, Ssrf, Sbot, S0, dz1, dz2,    &
&                dTsrf, dTbot, Tf, meltInt, time00, sealevel, lbd1, lbd2, lbd3, meltfac, K, gT, alphap,            &
&                dn, wn, div, x0, E0, Cd, GamT, GefT, GTS, MM, M_hat, X_hat, ll, AvgTmTf, localunity,              &
&                Area, Area_Reduced, Za, rhostar, CC, beta, g1, g2, rr, S_stratif, timsc, qqq, xerr,               &
&                mskcrit, time, Hcrit, Mcrit, isfslope_x, isfslope_y, tmp1, Tstar, xbox, lgl, closestgl,           &
&                locslope1, locslope2, TF_avg, TF_tot

REAL(KIND=8) :: yearinsec, rhoi_SI, rhosw_SI, rhofw_SI, Lf_SI, cpw_SI, grav_SI, lambda1, lambda2, lambda3,    &
&               rhoi, rhosw, rhofw, Lf, cpw, gravity, gammaT_SI, facPDC_SI, E0_LAZER, GT_LAZER, aaa, zz1, zz2,&
&               GTS0_LAZER, Zambiant_PICO, gammaT_PICO_SI, C_PICO_SI, alpha_PICO, beta_PICO, rhostar_SI_PICO, &
&               gammaT, facPDC, gammaT_PICO, C_PICO, rhostar_PICO, dist, mindist, lonGLmin, latGLmin, zGLmin, sn, snt

REAL(KIND=8) :: det, lonP2, latP2, lonP3, latP3, lonP12, latP12, lonP23, latP23, lonP31, latP31, radius, angle, angle1, angle2, &
&                zGLtmp, zGLloc, ltmp, zGLcp, zGLdp, xcp, ycp, xdp, ydp, x1, y1, x2, y2, zz, reflon, reflat, angle4,            &
&                locslope, mmm_tot, mmm_avg

REAL(KIND=8), DIMENSION(12) :: pp

REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::  Tbox, Sbox, Mbox, lonGLv, latGLv, isfdraftv, sinus, isfdraftvs

REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::  Zbox, Abox

REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: zGL, alpha

REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: Pmelt_out, Fmelt_out
REAL(KIND=16), DIMENSION(:,:,:), ALLOCATABLE :: Pmelt, Fmelt

INTEGER(KIND=4),ALLOCATABLE,DIMENSION(:,:,:) :: Pnnnn, Fnnnn

INTEGER(KIND=2),ALLOCATABLE,DIMENSION(:) :: indx, iin, jjn

INTEGER(KIND=4) :: fidA, status, dimID_y, dimID_x, my, mx, thickness_ID, bedrock_ID, yy_ID, xx_ID, Nbox_PICO, fidM, Melt_ID, &
&                  zGL_ID, XHat_ID, MHat_ID, MM_ID, alpha_ID, GM_ID, cpttmp, aa, nbDir, surface_ID

CHARACTER(LEN=150) :: file_msh, file_out                     

CHARACTER(LEN=10) ::  para_name, str_Hcrit, str_Integ, TYPE_LAZER, TypeL

!REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: xx, yy
!
!REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: thickness, bedrock, surface, lsrf, isfdraft
 
REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: dGL, dIF

REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: Melt

!INTEGER,ALLOCATABLE,DIMENSION(:,:) :: GM, IM, mask

LOGICAL :: llGL, cond1, cond2, cond3, llnD

!--------------------------------------------------

INTEGER(KIND=4) :: fidBED, dimID_londim, dimID_latdim, mlondim, mlatdim, lon_ID, lat_ID, bedrock_topography_ID,    &
&                  isfmask_ID, fidDFT, ice_base_topography_ID, surface_elevation_ID, dimID_Nisf, mlondim2,         &
&                  name_isf_ID, name_reg_ID, strlen, dimID_strlen, fidISF, err_melt_isf_ID, melt_isf_ID, mlatdim2, &
&                  IF_mask_ID, GL_mask_ID

CHARACTER(LEN=150) :: file_bed, file_isf, file_dft

INTEGER(KIND=2),ALLOCATABLE,DIMENSION(:,:) :: isfmask, &       ! 0 grounded; 1 ocean; >1 ice shelf number
&                                             IF_mask, GL_mask ! Ice shelf Front (IF) and Grounding Line (GL) mask

REAL(KIND=4),ALLOCATABLE,DIMENSION(:) :: lon, lat, err_melt_isf, melt_isf

REAL(KIND=4),ALLOCATABLE,DIMENSION(:,:) :: bedrock_topography, ice_base_topography, surface_elevation         

REAL(KIND=4) :: itmpmin, itmpmax, jtmpmin, jtmpmax, dlon, dlat

REAL(KIND=8) :: RT, deg2rad,  aainf, aasup, target_melt

!--------------------------------------------------
INTEGER(KIND=4) :: fidWOA, dimID_depth, mdepth, mdepth2, mNisf, mNisf2, s_er_3_ID, s_an_3_ID, t_er_3_ID, t_an_3_ID, s_er_2_ID, s_an_2_ID, t_er_2_ID, t_an_2_ID, s_er_1_ID, s_an_1_ID, t_er_1_ID, t_an_1_ID, fidCMIP, dimID_Nmod, mStrLen, mNmod, s_anom_3_ID, t_anom_3_ID, s_anom_2_ID, t_anom_2_ID, s_anom_1_ID, t_anom_1_ID, depth_ID, model_name_ID, front_max_lat_ID, front_min_lat_ID, front_max_lon_ID, front_min_lon_ID, front_ice_dep_avg_ID, front_ice_dep_min_ID, front_bot_dep_avg_ID, front_bot_dep_max_ID, kisf, kk_TS_pres, nn_TS_pres, kk_tuning, nn_tuning, kk_CMIP5_anom, kk_para, nn_para, mtot, kksup, kkinf, stun, nn_TS

CHARACTER(LEN=150) :: file_TS_WOA, file_CMIP5_anom

INTEGER(KIND=4),ALLOCATABLE,DIMENSION(:) :: kstat, imin, imax, jmin, jmax

INTEGER(KIND=2),ALLOCATABLE,DIMENSION(:,:) :: index_para, index_CMIP, index_WOA, index_tun
 
REAL*4,ALLOCATABLE,DIMENSION(:) :: front_max_lat, front_min_lat, front_max_lon, front_min_lon, front_ice_dep_avg,   &
&                                  front_ice_dep_min, front_bot_dep_avg, front_bot_dep_max, depth, T_pres, S_pres, T_futu, S_futu
 
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: s_er_3, s_an_3, t_er_3, t_an_3, s_er_2, s_an_2, t_er_2, t_an_2, s_er_1, s_an_1, t_er_1, t_an_1, &
&                                    total_melt_pres, mean_melt_pres, total_melt_futu, mean_melt_futu
 
REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: s_anom_3, t_anom_3, s_anom_2, t_anom_2, s_anom_1, t_anom_1

!--------------------------------------------------

file_dft         = '/store/njourd/DATA/DATA_BATHYMETRY/RTopo-2.0.1_30sec_ice_base_topography_ANT3.nc'
file_bed         = '/store/njourd/DATA/DATA_BATHYMETRY/RTopo-2.0.1_30sec_bedrock_topography_ANT3.nc'
file_isf         = 'RTopo-2.0.1_30sec_ice_shelf_mask_ANT3.nc'
file_TS_WOA      = 'WOA13_TS_per_ice_shelf.nc'
file_CMIP5_anom  = 'CMIP5_anom_TS_per_ice_shelf.nc'

!file_out = 'melt_ensemble_AMUNDSEN.nc'

!strlen=LEN(name_isf)

deg2rad = dacos(-1.d0) / 180.d0
RT = 6.371d6 ! Earth radius in m
write(*,*) 'pi = ', dacos(-1.d0)

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

!---------------------------------------
! Read TS in WOA 2013 :
 
write(*,*) 'Reading ', TRIM(file_TS_WOA)
 
status = NF90_OPEN(TRIM(file_TS_WOA),0,fidWOA); call erreur(status,.TRUE.,"read WOA T,S profiles")
 
! Read dimension IDs
status = NF90_INQ_DIMID(fidWOA,"depth",dimID_depth); call erreur(status,.TRUE.,"inq_dimID_depth")
status = NF90_INQ_DIMID(fidWOA,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
 
! Read values dimension values
status = NF90_INQUIRE_DIMENSION(fidWOA,dimID_depth,len=mdepth); call erreur(status,.TRUE.,"inq_dim_depth")
status = NF90_INQUIRE_DIMENSION(fidWOA,dimID_Nisf,len=mNisf); call erreur(status,.TRUE.,"inq_dim_Nisf")
  
! Allocation of arrays : 
ALLOCATE(  s_er_3(mNisf,mdepth)  ) 
ALLOCATE(  s_an_3(mNisf,mdepth)  ) 
ALLOCATE(  t_er_3(mNisf,mdepth)  ) 
ALLOCATE(  t_an_3(mNisf,mdepth)  ) 
ALLOCATE(  s_er_2(mNisf,mdepth)  ) 
ALLOCATE(  s_an_2(mNisf,mdepth)  ) 
ALLOCATE(  t_er_2(mNisf,mdepth)  ) 
ALLOCATE(  t_an_2(mNisf,mdepth)  ) 
!ALLOCATE(  s_er_1(mNisf,mdepth)  ) 
!ALLOCATE(  s_an_1(mNisf,mdepth)  ) 
!ALLOCATE(  t_er_1(mNisf,mdepth)  ) 
!ALLOCATE(  t_an_1(mNisf,mdepth)  ) 
ALLOCATE(  depth(mdepth)  ) 
 
! Read variable IDs
status = NF90_INQ_VARID(fidWOA,"s_er_3",s_er_3_ID); call erreur(status,.TRUE.,"inq_s_er_3_ID")
status = NF90_INQ_VARID(fidWOA,"s_an_3",s_an_3_ID); call erreur(status,.TRUE.,"inq_s_an_3_ID")
status = NF90_INQ_VARID(fidWOA,"t_er_3",t_er_3_ID); call erreur(status,.TRUE.,"inq_t_er_3_ID")
status = NF90_INQ_VARID(fidWOA,"t_an_3",t_an_3_ID); call erreur(status,.TRUE.,"inq_t_an_3_ID")
status = NF90_INQ_VARID(fidWOA,"s_er_2",s_er_2_ID); call erreur(status,.TRUE.,"inq_s_er_2_ID")
status = NF90_INQ_VARID(fidWOA,"s_an_2",s_an_2_ID); call erreur(status,.TRUE.,"inq_s_an_2_ID")
status = NF90_INQ_VARID(fidWOA,"t_er_2",t_er_2_ID); call erreur(status,.TRUE.,"inq_t_er_2_ID")
status = NF90_INQ_VARID(fidWOA,"t_an_2",t_an_2_ID); call erreur(status,.TRUE.,"inq_t_an_2_ID")
!status = NF90_INQ_VARID(fidWOA,"s_er_1",s_er_1_ID); call erreur(status,.TRUE.,"inq_s_er_1_ID")
!status = NF90_INQ_VARID(fidWOA,"s_an_1",s_an_1_ID); call erreur(status,.TRUE.,"inq_s_an_1_ID")
!status = NF90_INQ_VARID(fidWOA,"t_er_1",t_er_1_ID); call erreur(status,.TRUE.,"inq_t_er_1_ID")
!status = NF90_INQ_VARID(fidWOA,"t_an_1",t_an_1_ID); call erreur(status,.TRUE.,"inq_t_an_1_ID")
status = NF90_INQ_VARID(fidWOA,"depth",depth_ID); call erreur(status,.TRUE.,"inq_depth_ID")
 
! Read variable values
status = NF90_GET_VAR(fidWOA,s_er_3_ID,s_er_3); call erreur(status,.TRUE.,"getvar_s_er_3")
status = NF90_GET_VAR(fidWOA,s_an_3_ID,s_an_3); call erreur(status,.TRUE.,"getvar_s_an_3")
status = NF90_GET_VAR(fidWOA,t_er_3_ID,t_er_3); call erreur(status,.TRUE.,"getvar_t_er_3")
status = NF90_GET_VAR(fidWOA,t_an_3_ID,t_an_3); call erreur(status,.TRUE.,"getvar_t_an_3")
status = NF90_GET_VAR(fidWOA,s_er_2_ID,s_er_2); call erreur(status,.TRUE.,"getvar_s_er_2")
status = NF90_GET_VAR(fidWOA,s_an_2_ID,s_an_2); call erreur(status,.TRUE.,"getvar_s_an_2")
status = NF90_GET_VAR(fidWOA,t_er_2_ID,t_er_2); call erreur(status,.TRUE.,"getvar_t_er_2")
status = NF90_GET_VAR(fidWOA,t_an_2_ID,t_an_2); call erreur(status,.TRUE.,"getvar_t_an_2")
!status = NF90_GET_VAR(fidWOA,s_er_1_ID,s_er_1); call erreur(status,.TRUE.,"getvar_s_er_1")
!status = NF90_GET_VAR(fidWOA,s_an_1_ID,s_an_1); call erreur(status,.TRUE.,"getvar_s_an_1")
!status = NF90_GET_VAR(fidWOA,t_er_1_ID,t_er_1); call erreur(status,.TRUE.,"getvar_t_er_1")
!status = NF90_GET_VAR(fidWOA,t_an_1_ID,t_an_1); call erreur(status,.TRUE.,"getvar_t_an_1")
status = NF90_GET_VAR(fidWOA,depth_ID,depth); call erreur(status,.TRUE.,"getvar_depth")
 
! Close file
status = NF90_CLOSE(fidWOA); call erreur(status,.TRUE.,"close_file")

!---------------------------------------
! Read CMIP5 anomaly :
 
write(*,*) 'Reading ', TRIM(file_CMIP5_anom)
 
status = NF90_OPEN(TRIM(file_CMIP5_anom),0,fidCMIP); call erreur(status,.TRUE.,"read CMIP5 anomaly")
 
! Read dimension IDs
status = NF90_INQ_DIMID(fidCMIP,"StrLen",dimID_StrLen); call erreur(status,.TRUE.,"inq_dimID_StrLen")
status = NF90_INQ_DIMID(fidCMIP,"depth",dimID_depth); call erreur(status,.TRUE.,"inq_dimID_depth")
status = NF90_INQ_DIMID(fidCMIP,"Nmod",dimID_Nmod); call erreur(status,.TRUE.,"inq_dimID_Nmod")
status = NF90_INQ_DIMID(fidCMIP,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
 
! Read values dimension values
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

! Allocation of arrays : 
ALLOCATE(  s_anom_3(mNmod,mNisf,mdepth)  ) 
ALLOCATE(  t_anom_3(mNmod,mNisf,mdepth)  ) 
ALLOCATE(  s_anom_2(mNmod,mNisf,mdepth)  ) 
ALLOCATE(  t_anom_2(mNmod,mNisf,mdepth)  ) 
!ALLOCATE(  s_anom_1(mNmod,mNisf,mdepth)  ) 
!ALLOCATE(  t_anom_1(mNmod,mNisf,mdepth)  ) 
!ALLOCATE(  depth(mdepth)  ) 
 
! Read variable IDs
status = NF90_INQ_VARID(fidCMIP,"s_anom_3",s_anom_3_ID); call erreur(status,.TRUE.,"inq_s_anom_3_ID")
status = NF90_INQ_VARID(fidCMIP,"t_anom_3",t_anom_3_ID); call erreur(status,.TRUE.,"inq_t_anom_3_ID")
status = NF90_INQ_VARID(fidCMIP,"s_anom_2",s_anom_2_ID); call erreur(status,.TRUE.,"inq_s_anom_2_ID")
status = NF90_INQ_VARID(fidCMIP,"t_anom_2",t_anom_2_ID); call erreur(status,.TRUE.,"inq_t_anom_2_ID")
!status = NF90_INQ_VARID(fidCMIP,"s_anom_1",s_anom_1_ID); call erreur(status,.TRUE.,"inq_s_anom_1_ID")
!status = NF90_INQ_VARID(fidCMIP,"t_anom_1",t_anom_1_ID); call erreur(status,.TRUE.,"inq_t_anom_1_ID")
status = NF90_INQ_VARID(fidCMIP,"depth",depth_ID); call erreur(status,.TRUE.,"inq_depth_ID")
 
! Read variable values
status = NF90_GET_VAR(fidCMIP,s_anom_3_ID,s_anom_3); call erreur(status,.TRUE.,"getvar_s_anom_3")
status = NF90_GET_VAR(fidCMIP,t_anom_3_ID,t_anom_3); call erreur(status,.TRUE.,"getvar_t_anom_3")
status = NF90_GET_VAR(fidCMIP,s_anom_2_ID,s_anom_2); call erreur(status,.TRUE.,"getvar_s_anom_2")
status = NF90_GET_VAR(fidCMIP,t_anom_2_ID,t_anom_2); call erreur(status,.TRUE.,"getvar_t_anom_2")
!status = NF90_GET_VAR(fidCMIP,s_anom_1_ID,s_anom_1); call erreur(status,.TRUE.,"getvar_s_anom_1")
!status = NF90_GET_VAR(fidCMIP,t_anom_1_ID,t_anom_1); call erreur(status,.TRUE.,"getvar_t_anom_1")
status = NF90_GET_VAR(fidCMIP,depth_ID,depth); call erreur(status,.TRUE.,"getvar_depth")
 
! Close file
status = NF90_CLOSE(fidCMIP); call erreur(status,.TRUE.,"close_file")

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
E0_LAZER      = 3.6e-2                  ! Entrainment coeff. (no unit) in Lazeroms et al.
GT_LAZER      = 1.1e-5                  ! Turbulent heat exchange coeff. (no unit) in Lazeroms et al.
GTS0_LAZER    = 6.0e-4                  ! Effective heat exchange coeff. (no unit) in Lazeroms et al.
pp = (/  1.371330075095435e-1, &  ! p0
      &  5.527656234709359e1,  &  ! p1
      & -8.951812433987858e2,  &  ! p2
      &  8.927093637594877e3,  &  ! p3
      & -5.563863123811898e4,  &  ! p4 
      &  2.218596970948727e5,  &  ! p5
      & -5.820015295669482e5,  &  ! p6
      &  1.015475347943186e6,  &  ! p7
      & -1.166290429178556e6,  &  ! p8
      &  8.466870335320488e5,  &  ! p9
      & -3.520598035764990e5,  &  ! p10
      &  6.387953795485420e4   /) ! p11
x0   = 0.56    ! Dimensionless transition melting/freezing
Cd   = 2.5e-3  ! Drag coefficent
! 16 directions :
!nbDir = 16
nbDir = 32
ALLOCATE( iin(nbDir), jjn(nbDir) )
!iin=(/  0,  1,  1,  2,  1,  2,  1,  1,  0, -1, -1, -2, -1, -2, -1, -1 /)
!jjn=(/ -1, -2, -1, -1,  0,  1,  1,  2,  1,  2,  1,  1,  0, -1, -1, -2 /)
iin=(/  0,  1,  1,  2,  1,  2,  1,  1,  0, -1, -1, -2, -1, -2, -1, -1,  3,  3,  2,  1, -1, -2, -3, -3, -3, -3, -2, -1,  1,  2,  3,  3 /)
jjn=(/ -1, -2, -1, -1,  0,  1,  1,  2,  1,  2,  1,  1,  0, -1, -1, -2,  1,  2,  3,  3,  3,  3,  2,  1, -1, -2, -3, -3, -3, -3, -2, -1 /)

! Box models :
gammaT_PICO_SI  = 2.0e-5
C_PICO_SI       = 1.0e6                 ! Circulation parameter C (m^6/kg/s) in [0.1;9]*1.e6
alpha_PICO      = 7.5e-5                ! Thermal expansion coefficient for linear EOS (K^-1)
beta_PICO       = 7.7e-4                ! Salinity contraction coefficient for linear EOS (psu^-1)
rhostar_SI_PICO = 1027.51               ! In situ density for linear EOS (kg/m3)

gammaT = gammaT_SI * yearinsec
facPDC = facPDC_SI * yearinsec
gammaT_PICO = gammaT_PICO_SI * yearinsec
C_PICO = C_PICO_SI * (1.0e6*yearinsec**2) * yearinsec
rhostar_PICO = rhostar_SI_PICO / (1.0e6*yearinsec**2)

!--------------------
lbd1 = lambda1
lbd2 = lambda2
lbd3 = lambda3

!===================================================================================
!===================================================================================

nn_TS_pres = 2   ! in {2, 6, 10, 14, 18} [recommended=10]
nn_tuning  = 19  ! in {1, 3, 5, 9, ...}  [recommended>15]
nn_para    = 17

mtot=nn_TS_pres*mNmod*nn_tuning*nn_para

ALLOCATE( T_pres(mdepth), S_pres(mdepth), T_futu(mdepth), S_futu(mdepth) )
ALLOCATE( total_melt_pres(mNisf,mtot), mean_melt_pres(mNisf,mtot) )
ALLOCATE( total_melt_futu(mNisf,mtot), mean_melt_futu(mNisf,mtot) )
ALLOCATE( kstat(mNisf) )
ALLOCATE( index_para(mNisf,mtot), index_CMIP(mNisf,mtot), index_WOA(mNisf,mtot), index_tun(mNisf,mtot) )

max_nb_box = 10 ! maximum nb of boxes used in PICO
ALLOCATE( Zbox(max_nb_box,max_nb_box), Abox(max_nb_box,max_nb_box) )

total_melt_pres(:,:) = NF90_FILL_DOUBLE
total_melt_futu(:,:) = NF90_FILL_DOUBLE
mean_melt_pres(:,:)  = NF90_FILL_DOUBLE
mean_melt_futu(:,:)  = NF90_FILL_DOUBLE
index_para(:,:) = 0
index_CMIP(:,:) = 0
index_WOA (:,:) = 0
index_tun (:,:) = 0

kstat(:)=0

364 FORMAT(i2.2,', ',i8.8,', ',i2.2,', ',i2.2,', ',i2.2,', ',i2.2,', ',f9.2,', ',f9.3,', ',f9.2,', ',f9.3)

write(*,*) 'MAXVAL(isfmask) = ', MAXVAL(isfmask)

!DO kisf=2,mNisf
DO kisf=66,66

 IF ( melt_isf(kisf) .ne. 0.e0 ) THEN

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
   ALLOCATE( zGL(mlondim,mlatdim,3), alpha(mlondim,mlatdim,3) )

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

   !===============================================================================

   !=== Calculate BOXES parameters for PICO (once for all) ===
   write(*,*) '     > Identifying boxes characteristics:'
   ! Calculate local distance to grounding line and ice front :
   DO ii=1,mlondim
   DO jj=1,mlatdim
     if ( isfmask(ii,jj) .eq. kisf ) then
       do kki=1,mlondim
       do kkj=1,mlatdim
         if ( GL_mask(kki,kkj) .eq. kisf ) then
           CALL dist_sphere(DBLE(lon(ii)),DBLE(lon(kki)),DBLE(lat(jj)),DBLE(lat(kkj)),dist)
           if ( dist .le. dGL(ii,jj) )  dGL(ii,jj) = dist
         endif
         if ( IF_mask(kki,kkj) .eq. kisf ) then
           CALL dist_sphere(DBLE(lon(ii)),DBLE(lon(kki)),DBLE(lat(jj)),DBLE(lat(kkj)),dist)
           if ( dist .le. dIF(ii,jj) )  dIF(ii,jj) = dist
         endif
       enddo
       enddo
     endif
   ENDDO
   ENDDO
   !- Calculate total area and mean depth of each box for
   !  configurations of total nb of boxes from 1 to 10 :
   Zbox(:,:) = 0.e0
   Abox(:,:) = 0.e0
   DO nD=1,max_nb_box ! nb of boxes
      do ii=1,mlondim
      do jj=1,mlatdim
         !rr = dGL(ii,jj) / ( dGL(ii,jj) + dIF(ii,jj) )
         do kk=1,nD
           if ( isfmask(ii,jj) .eq. kisf .and. dGL(ii,jj)/(dGL(ii,jj)+dIF(ii,jj)) .ge. 1.0-sqrt(1.0*(nD-kk+1)/nD) .and. dGL(ii,jj)/(dGL(ii,jj)+dIF(ii,jj)) .le. 1.0-sqrt(1.0*(nD-kk)/nD) ) then
             !zzz=-ice_base_topography(ii,jj)
             Zbox(kk,nD) = Zbox(kk,nD) - ice_base_topography(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
             Abox(kk,nD) = Abox(kk,nD) +                              dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2
           endif
         enddo
      enddo
      enddo
      do kk=1,nD
        if ( Abox(kk,nD) .gt. 0.0 ) then
          Zbox(kk,nD) = Zbox(kk,nD)/Abox(kk,nD) ! mean depth of each box
        else
          Abox(:,nD) = 0.e0  ! 
          Zbox(:,nD) = 0.e0  ! If there are empty boxes (i.e. not enough points in the cavity),
          exit               ! we skip this configuration
        endif
      enddo
      write(*,*)   '       nD = ', nD
      write(*,563) '       Zbox(1:nD,nD) = ', Zbox(1:nD,nD)
   ENDDO ! nD
   !- Check boxes and choose 3 configurations that work (nD=nD1,nD2,nD3) :
   !  (ideally : nD3=max_nb_box, nD2=max_nb_box/2, nD1=2)
   do nD=max_nb_box,3,-1
     nD3=nD
     llnD = .true.
     do kk=1,nD3
       if ( Abox(kk,nD3) .le. 0.0 ) llnD = .false.  ! no empty box
     enddo
     do kk=2,nD3
       if ( Zbox(kk,nD3) .gt. Zbox(kk-1,nD3) ) llnD = .false.  ! cavity must be monotonic
     enddo
     if ( llnD ) then
       exit
     else
       Abox(:,nD3) = 0.d0
     endif
   enddo
   do nD=MIN(nD3-1,max_nb_box/2),2,-1
     nD2=nD
     llnD = .true.
     do kk=1,nD2
       if ( Abox(kk,nD2) .le. 0.0 ) llnD = .false.  ! no empty box
     enddo
     do kk=2,nD2
       if ( Zbox(kk,nD2) .gt. Zbox(kk-1,nD2) ) llnD = .false.  ! cavity must be monotonic
     enddo
     if ( llnD ) then
       exit
     else
       Abox(:,nD2) = 0.d0
     endif
   enddo
   do nD=MIN(nD2-1,2),1,-1
     nD1=nD
     llnD = .true.
     do kk=1,nD1
       if ( Abox(kk,nD1) .le. 0.0 ) llnD = .false.  ! no empty box
     enddo
     do kk=2,nD1
       if ( Zbox(kk,nD1) .gt. Zbox(kk-1,nD1) ) llnD = .false.  ! cavity must be monotonic
     enddo
     if ( llnD ) then
       exit
     else
       Abox(:,nD1) = 0.d0
     endif
   enddo
   write(*,*)   '       nD3 = ', nD3
   write(*,563) '       Zbox(1:nD3,nD3) = ', Zbox(1:nD3,nD3)
   write(*,*)   '       nD2 = ', nD2
   write(*,563) '       Zbox(1:nD2,nD2) = ', Zbox(1:nD2,nD2)
   write(*,*)   '       nD1 = ', nD1
   write(*,563) '       Zbox(1:nD1,nD1) = ', Zbox(1:nD1,nD1)
   563 FORMAT(a,10(f8.1))

   !=== Calculate slope and GL parameters for the plume model (once for all) ===
   TypeL   = TYPE_LAZER    ! type of parameterization

   DO kkp=1,3

     if     ( kkp .eq. 1 ) then; TypeL='simple'        ! Zgl and Alpha are found between draft point and the central point
     elseif ( kkp .eq. 2 ) then; TypeL='lazero2'        ! original from Lazeroms et al.
     elseif ( kkp .eq. 3 ) then; TypeL='appenB2'; endif ! Appendix B in submitted TCD

     !- Prelimin(kisf)ary calculations : 

     if (TypeL .eq. 'simple') then
       zGLmin = 0.d0
       DO ii=1,mlondim
       DO jj=1,mlatdim
         zzz=-ice_base_topography(ii,jj)
         if ( GL_mask(ii,jj) .eq. kisf .and. zzz .gt. zGLmin ) then
           zGLmin = zzz
           lonGLmin = DBLE(lon(ii))
           latGLmin = DBLE(lat(jj))
         endif
       ENDDO
       ENDDO
       write(*,*) 'Deepest GL (z,lon,lat) = ', zGLmin, lonGLmin, latGLmin
     endif 

     ! Main calculation
     DO ii=1,mlondim
     DO jj=1,mlatdim
  
       zzz=-ice_base_topography(ii,jj)
  
       if ( isfmask(ii,jj) .eq. kisf ) then
 
         if (TypeL .eq. 'simple') then
  
           ! for the 'simple' law
           zGL(ii,jj,kkp) = zGLmin
           CALL dist_sphere(lonGLmin,DBLE(lon(ii)),latGLmin,DBLE(lat(jj)),dist)
           !!-NJ dist = sqrt((xGLmin-xx(ii))**2+(yGLmin-yy(jj))**2)
           if (dist .eq. 0.d0) then
             alpha(ii,jj,kkp) = 0.d0
           else
             alpha(ii,jj,kkp) = atan(abs(zzz-zGL(ii,jj,kkp))/dist)
           endif         

         elseif (TypeL .eq. 'lazero') then

           cptglob = 0
           zGLtmp = 0.d0
           snt = 0.d0
           do aa=1,nbDir
             ip1=MAX(MIN(mlondim,ii+iin(aa)),1)
             jp1=MAX(MIN(mlatdim,jj+jjn(aa)),1)
             if ( isfmask(ip1,jp1) .eq. kisf ) then
               CALL dist_sphere(DBLE(lon(ip1)),DBLE(lon(ii)),DBLE(lat(jp1)),DBLE(lat(jj)),dist)
               sn = ( ice_base_topography(ii,jj) - ice_base_topography(ip1,jp1) ) / dist  ! local slope
             else
               sn = -99999.d0
             endif
             if ( sn .gt. 0.0 ) then ! if direction is valid (1st criterion)
               do kk=1,max(mlondim,mlatdim)
                 !im1=MAX(MIN(mlondim,ii+(kk-1)*iin(aa)),1)
                 !jm1=MAX(MIN(mlatdim,jj+(kk-1)*jjn(aa)),1)
                 ip1=MAX(MIN(mlondim,ii+ kk   *iin(aa)),1)
                 jp1=MAX(MIN(mlatdim,jj+ kk   *jjn(aa)),1)
                 if ( isfmask(ip1,jp1) .ne. kisf .or. GL_mask(ip1,jp1) .eq. kisf ) exit
               enddo 
               if ( GL_mask(ip1,jp1) .eq. kisf .and. ice_base_topography(ip1,jp1) .lt. ice_base_topography(ii,jj) ) then ! 2nd criterion
                 zGLtmp = zGLtmp - ice_base_topography(ip1,jp1)
                 snt = snt + sn
                 cptglob = cptglob + 1
               endif
             endif
           enddo  ! aa
           if ( cptglob .gt. 0 ) then
             zGL(ii,jj,kkp) = zGLtmp / cptglob
             alpha(ii,jj,kkp) = atan( snt / cptglob )
           else
             zGL(ii,jj,kkp) = - ice_base_topography(ii,jj)
             alpha(ii,jj,kkp) = 0.d0
           endif
           write(*,376) ii, jj, -ice_base_topography(ii,jj), zGL(ii,jj,kkp), alpha(ii,jj,kkp)
           376 FORMAT('lazero (',i4,',',i4,') z=',f6.1,' -> zGL=',f6.1,' alpha=',f9.6)

         elseif (TypeL .eq. 'appenB' ) then

           ! Local slope angle 
           im1=MAX(1,ii-1)  
           ip1=MIN(mlondim,ii+1)  
           jm1=MAX(1,jj-1)  
           jp1=MIN(mlatdim,jj+1)  
           if ( isfmask(im1,jj) .eq. kisf .and. isfmask(ip1,jj) .eq. kisf ) then
             CALL dist_sphere(DBLE(lon(ip1)),DBLE(lon(im1)),DBLE(lat(jj)),DBLE(lat(jj)),dist)
             isfslope_x = ( -ice_base_topography(ip1,jj) + ice_base_topography(im1,jj) ) / dist
           elseif ( isfmask(im1,jj) .eq. kisf .and. im1 .ne. ii ) then
             CALL dist_sphere(DBLE(lon(ii)),DBLE(lon(im1)),DBLE(lat(jj)),DBLE(lat(jj)),dist)
             isfslope_x = ( -ice_base_topography(ii  ,jj) + ice_base_topography(im1,jj) ) / dist
           elseif ( isfmask(ip1,jj) .eq. kisf .and. ip1 .ne. ii ) then
             CALL dist_sphere(DBLE(lon(ip1)),DBLE(lon(ii)),DBLE(lat(jj)),DBLE(lat(jj)),dist)
             isfslope_x = ( -ice_base_topography(ip1,jj) + ice_base_topography(ii  ,jj) ) / dist
           else
             isfslope_x = 0.d0
           endif
           if ( isfmask(ii,jm1) .eq. kisf .and. isfmask(ii,jp1) .eq. kisf ) then
             CALL dist_sphere(DBLE(lon(ii)),DBLE(lon(ii)),DBLE(lat(jp1)),DBLE(lat(jm1)),dist)
             isfslope_y = ( -ice_base_topography(ii,jp1) + ice_base_topography(ii,jm1) ) / dist
           elseif ( isfmask(ii,jm1) .eq. kisf .and. jm1 .ne. jj ) then
             CALL dist_sphere(DBLE(lon(ii)),DBLE(lon(ii)),DBLE(lat(jj)),DBLE(lat(jm1)),dist)
             isfslope_y = ( -ice_base_topography(ii  ,jj) + ice_base_topography(ii,jm1) ) / dist
           elseif ( isfmask(ii,jp1) .eq. kisf .and. jp1 .ne. jj ) then
             CALL dist_sphere(DBLE(lon(ii)),DBLE(lon(ii)),DBLE(lat(jp1)),DBLE(lat(jj)),dist)
             isfslope_y = ( -ice_base_topography(ii,jp1) + ice_base_topography(ii  ,jj) ) / dist
           else
             isfslope_y = 0.d0
           endif
           alpha(ii,jj,kkp) = atan( sqrt( isfslope_x**2 + isfslope_y**2 ) )

           ! Effective grounding line depth :
           zGLtmp = 0.d0
           snt = 0.d0
           do aa=1,nbDir
             ip1=MAX(MIN(mlondim,ii+iin(aa)),1)
             jp1=MAX(MIN(mlatdim,jj+jjn(aa)),1)
             if ( isfmask(ip1,jp1) .eq. kisf ) then
               CALL dist_sphere(DBLE(lon(ip1)),DBLE(lon(ii)),DBLE(lat(jp1)),DBLE(lat(jj)),dist)
               sn = ( ice_base_topography(ii,jj) - ice_base_topography(ip1,jp1) ) / dist
             else
               sn = -99999.d0
             endif
             if ( sn .gt. 0.0 ) then ! if direction is valid (1st criterion)
               do kk=2,max(mlondim,mlatdim)
                 !im1=MAX(MIN(mlondim,ii+(kk-1)*iin(aa)),1)
                 !jm1=MAX(MIN(mlatdim,jj+(kk-1)*jjn(aa)),1)
                 ip1=MAX(MIN(mlondim,ii+ kk   *iin(aa)),1)
                 jp1=MAX(MIN(mlatdim,jj+ kk   *jjn(aa)),1)
                 if ( isfmask(ip1,jp1) .ne. kisf .or. GL_mask(ip1,jp1) .eq. kisf ) exit
               enddo
               if ( GL_mask(ip1,jp1) .eq. kisf .and. ice_base_topography(ip1,jp1) .lt. ice_base_topography(ii,jj) ) then ! 2nd criterion
                 CALL dist_sphere(DBLE(lon(ii)),DBLE(lon(ip1)),DBLE(lat(jj)),DBLE(lat(jp1)),dist)
                 sn = MAX( 0.d0,  ( ice_base_topography(ii,jj)-ice_base_topography(ip1,jp1)) / dist ) ! weight factor (overall slope)
                 zGLtmp = zGLtmp - ice_base_topography(ip1,jp1) * sn 
                 snt = snt + sn
               endif
             endif
           enddo  ! aa
           if ( snt .gt. 0.d0 ) then
             zGL(ii,jj,kkp) = zGLtmp / snt
           else
             zGL(ii,jj,kkp) = - ice_base_topography(ii,jj)
           endif
           write(*,377) ii, jj, -ice_base_topography(ii,jj), zGL(ii,jj,kkp), alpha(ii,jj,kkp)
           377 FORMAT('appendixB (',i4,',',i4,') z=',f6.1,' -> zGL=',f6.1,' alpha=',f9.6)

         elseif (TypeL .eq. 'appenB2' .or. TypeL .eq. 'lazero2') then  ! Lionel's version with triangles in each direction
  
           ! Effective slope angle 
           ! tan(alpha) = sqrt( (dZb/dx)^2 + (dZb/dy)^2 ) (Lazeroms et al. 2017, Appendix B) :
           if ( isfmask(ii-1,jj) .eq. kisf .and. isfmask(ii+1,jj) .eq. kisf ) then
             CALL dist_sphere(lon(ii+1),lon(ii-1),lat(jj),lat(jj),dist)
             isfslope_x = ( -ice_base_topography(ii+1,jj) + ice_base_topography(ii-1,jj) ) / dist
           elseif ( isfmask(ii-1,jj) .eq. kisf ) then
             CALL dist_sphere(lon(ii),lon(ii-1),lat(jj),lat(jj),dist)
             isfslope_x = ( -ice_base_topography(ii  ,jj) + ice_base_topography(ii-1,jj) ) / dist
           elseif ( isfmask(ii+1,jj) .eq. kisf ) then
             CALL dist_sphere(lon(ii+1),lon(ii),lat(jj),lat(jj),dist)
             isfslope_x = ( -ice_base_topography(ii+1,jj) + ice_base_topography(ii  ,jj) ) / dist
           else
             isfslope_x = 0.d0
           endif
           if ( isfmask(ii,jj-1) .eq. kisf .and. isfmask(ii,jj+1) .eq. kisf ) then
             CALL dist_sphere(lon(ii),lon(ii),lat(jj+1),lat(jj-1),dist)
             isfslope_y = ( -ice_base_topography(ii,jj+1) + ice_base_topography(ii,jj-1) ) / dist
           elseif ( isfmask(ii,jj-1) .eq. kisf ) then
             CALL dist_sphere(lon(ii),lon(ii),lat(jj),lat(jj-1),dist)
             isfslope_y = ( -ice_base_topography(ii  ,jj) + ice_base_topography(ii,jj-1) ) / dist
           elseif ( isfmask(ii,jj+1) .eq. kisf ) then
             CALL dist_sphere(lon(ii),lon(ii),lat(jj+1),lat(jj),dist)
             isfslope_y = ( -ice_base_topography(ii,jj+1) + ice_base_topography(ii  ,jj) ) / dist
           else
             isfslope_y = 0.d0
           endif
 
           sn = 0.d0
           cptglob = 0
           ! define a radius (nj: here in the lon/lat space)
           !radius = 1.2 * sqrt( (lon(1)-lon(mlondim))**2 + (lat(1)-lat(mlatdim))**2 )
           radius = 1.1 * 360.d0
           nbDir = 16
           angle = (360.d0/nbdir)*4*atan(1.d0)/180.d0
           do aa=1,nbDir
             !
             ! calculate the two angles
             angle1 = (aa-1) * angle
             angle2 =  aa    * angle
             ! calculate 2 local slopes
             locslope1 = isfslope_x * COS(angle1) + isfslope_y * SIN(angle1)
             locslope2 = isfslope_x * COS(angle2) + isfslope_y * SIN(angle2)
             !locslope1 =   (   COS(ATAN(isfslope_x)) * SIN(ATAN(isfslope_y)) * SIN(angle1)   &
             !&               + COS(ATAN(isfslope_y)) * SIN(ATAN(isfslope_x)) * COS(angle1) ) &
             !&           / ( COS(ATAN(isfslope_x)) * COS(ATAN(isfslope_y)) )
             !locslope2 =   (   COS(ATAN(isfslope_x)) * SIN(ATAN(isfslope_y)) * SIN(angle2)   &
             !&               + COS(ATAN(isfslope_y)) * SIN(ATAN(isfslope_x)) * COS(angle2) ) &
             !&           / ( COS(ATAN(isfslope_x)) * COS(ATAN(isfslope_y)) )
             locslope = 0.5d0 * (locslope1+locslope2)
             !write(*,*) '%% ', aa, locslope, angle2
             ! if the average is negative, cycle
             IF (locslope .GE. 0.d0) CYCLE
             ! search the GL points inside a triangle (nj: here in the lon/lat space)
             ! define the two other points
             lonP2=lon(ii)+radius*cos(angle1)
             latP2=lat(jj)+radius*sin(angle1)
             lonP3=lon(ii)+radius*cos(angle2)
             latP3=lat(jj)+radius*sin(angle2)
             !
             !lonP12=lon(ii)-lonP2
             !latP12=lat(jj)-latP2
             !lonP23=lonP2-lonP3
             !latP23=latP2-latP3
             !lonP31=lonP3-lon(ii)
             !latP31=latP3-lat(jj)
             ! algorithm to find GL points inside the triangle could be two solutions
             ! - the 1st is to take all of them
             ! - the 2nd is to take the first GL point
             cpttmp=0
             zGLtmp=0.d0
             ltmp=0.d0
             div=0.d0
             !
             do kki=1,mlondim
             do kkj=1,mlatdim
               if (GL_mask(kki,kkj) .eq. kisf) then
                 ! The triangle is between points P2, P3 and (ii,jj).
                 ! Considering 1 segment, we look whether the opposite point of the triangle 
                 ! and the point under consideration (kki,kkj) are on the same side of the segment :
                 cond12 =   sign( 1.d0, det( lon(kki)-lonP2, lat(kkj)-latP2, lon(ii)-lonP2, lat(jj)-latP2 ) ) & ! segment [(ii,jj)-P2]
                 &        * sign( 1.d0, det( lonP3   -lonP2, latP3   -latP2, lon(ii)-lonP2, lat(jj)-latP2 ) )
                 cond13 =   sign( 1.d0, det( lon(kki)-lonP3, lat(kkj)-latP3, lon(ii)-lonP3, lat(jj)-latP3 ) ) & ! segment [(ii,jj)-P3]
                 &        * sign( 1.d0, det( lonP2   -lonP3, latP2   -latP3, lon(ii)-lonP3, lat(jj)-latP3 ) )
                 cond23 =   sign( 1.d0, det( lon(kki)-lonP3, lat(kkj)-latP3, lonP2  -lonP3, latP2  -latP3 ) ) & ! segment [P2-P3]
                 &        * sign( 1.d0, det( lon(ii) -lonP3, lat(jj) -latP3, lonP2  -lonP3, latP2  -latP3 ) )
                 !cond1=sign(1.d0,det(lonP31,latP31,lonP23,latP23))*sign(1.d0,det(lonP3-lon(kki),latP3-lat(kkj),lonP23,latP23)) .ge. 0.d0
                 !cond2=sign(1.d0,det(lonP12,latP12,lonP31,latP31))*sign(1.d0,det(lon(ii)-lon(kki),lat(jj)-lat(kkj),lonP31,latP31)) .ge. 0.d0
                 !cond3=sign(1.d0,det(lonP23,latP23,lonP12,latP12))*sign(1.d0,det(lonP2-lon(kki),latP2-lat(kkj),lonP12,latP12)) .ge. 0.d0
                 if ( cond12 .ge. 0.d0 .and. cond13 .ge. 0.d0 .and. cond23 .ge. 0.d0 ) then
                   zGLtmp=zGLtmp-ice_base_topography(kki,kkj)
                   cpttmp=cpttmp+1
                 endif
               endif
             enddo
             enddo
             !
             if (cpttmp .eq. 0) CYCLE
             zGLloc = zGLtmp / cpttmp
             !
             if (zGLloc .gt. zzz) then
               zGL(ii,jj,kkp) = zGL(ii,jj,kkp) + zGLloc
               sn = sn + abs(locslope)
               cptglob = cptglob + 1
             endif
             !
             !write(*,*) '@@', aa, zGL(ii,jj,kkp), cptglob, sn
             !write(*,*) '##', aa, zGLloc, ltmp, cpttmp, sqrt( isfslope_x**2 + isfslope_y**2 )
             IF ( TypeL .EQ. 'lazero2' ) THEN 
               IF (cptglob .GT. 0) THEN
                 zGL(ii,jj,kkp) = zGL(ii,jj,kkp) / cptglob
                 alpha(ii,jj,kkp) = ATAN(sn/cptglob)
               ELSE
                 zGL(ii,jj,kkp) = zzz
                 alpha(ii,jj,kkp) = 0.d0
               ENDIF
             ELSEIF (TypeL .EQ. 'appenB2') THEN
               alpha(ii,jj,kkp) = atan( sqrt( isfslope_x**2 + isfslope_y**2 ) )
               if (zGLloc .gt. zzz) then
                 ltmp=ltmp/cpttmp
                 wn=(zGLloc-zzz)/ltmp
                 zGL(ii,jj,kkp)=zGL(ii,jj,kkp)+wn*zGLloc
                 div=div+wn
               endif
             ENDIF
             !
           enddo  ! do aa=1,nbDir
           !
           IF (TypeL .EQ. 'appenB2') THEN
             if (div .gt. 0.d0) then
               zGL(ii,jj,kkp)=zGL(ii,jj,kkp)/div
             else
               zGL(ii,jj,kkp)=zzz
               alpha(ii,jj,kkp)=0.d0
             endif
           ENDIF

         endif  ! if ( TypeL .eq. 'lazero' .or. TypeL .eq. 'appenB' ) then
 
       endif ! if ( isfmask .eq. kisf )
 
     ENDDO
     ENDDO

   ENDDO !! kkp

   !==========

   DO kk_tuning=1,nn_tuning  !== Total cavity melt used for param calibration ==
   !DO kk_tuning=1,1

     ! Assuming that Rignot's uncertainties represent a normal distribution :
     ! Keeping a symmetric distribution to preserve the mean.
     CALL norm_sample(kk_tuning,nn_tuning,xerr)

     IF ( melt_isf(kisf) + xerr * err_melt_isf(kisf) .gt. 1.e-6 .AND. melt_isf(kisf) - xerr * err_melt_isf(kisf) .gt. 1.e-6 ) THEN

      target_melt = melt_isf(kisf) + xerr * err_melt_isf(kisf) ! Target cavity melt in Gt/yr
      write(*,*) '   target melt = ', target_melt

      DO kk_TS_pres=1,nn_TS_pres  !== Present-day T,S profiles ==
      !DO kk_TS_pres=1,1

        ! Assuming that WOA2013 uncertainties represent a normal distribution :
        CALL norm_sample(INT((nn_TS_pres-kk_TS_pres)/2)+1,INT(nn_TS_pres/2),xerr)
        if ( MOD(kk_TS_pres,2) .eq. 0 ) then
          ! WOA 2013 over the "nearby" continental shelf
          T_pres(:) = t_an_2(kisf,:) + xerr * t_er_2(kisf,:)
          S_pres(:) = s_an_2(kisf,:) + xerr * s_er_2(kisf,:)
          nn_TS = 2
        else
          ! WOA 2013 offshore "nearby"
          T_pres(:) = t_an_3(kisf,:) + xerr * t_er_3(kisf,:)
          S_pres(:) = s_an_3(kisf,:) + xerr * s_er_3(kisf,:)
          nn_TS = 3
        endif

        DO kk_CMIP5_anom = 1,mNmod  !== CMIP5 anomaly ==
        !DO kk_CMIP5_anom = 1,3
       
           if ( nn_TS .eq. 1 ) then
             T_futu(:) = T_pres(:) + t_anom_1(kk_CMIP5_anom,kisf,:)
             S_futu(:) = S_pres(:) + s_anom_1(kk_CMIP5_anom,kisf,:)
           elseif ( nn_TS .eq. 2 ) then
             ! T,S over the "nearby" continental shelf
             T_futu(:) = T_pres(:) + t_anom_2(kk_CMIP5_anom,kisf,:)
             S_futu(:) = S_pres(:) + s_anom_2(kk_CMIP5_anom,kisf,:)
           elseif ( nn_TS .eq. 3 ) then
             ! T,S taken offshore "nearby"
             T_futu(:) = T_pres(:) + t_anom_3(kk_CMIP5_anom,kisf,:)
             S_futu(:) = S_pres(:) + s_anom_3(kk_CMIP5_anom,kisf,:)
           else
             write(*,*) '~!@#$% error with T_futu and S_futu calculation >>>>>> STOP !!'
             stop
           endif
   
           !DO kk_para=1,nn_para  !== Melt parameterization ==
           DO kk_para=14,16  !== Melt parameterization ==
  
              kstat(kisf) = kstat(kisf)+1

              index_para(kisf,kstat(kisf)) = kk_para
              index_CMIP(kisf,kstat(kisf)) = kk_CMIP5_anom
              index_WOA (kisf,kstat(kisf)) = kk_TS_pres
              index_tun (kisf,kstat(kisf)) = kk_tuning 

              Melt(:,:) = 0.e0

              if ( MOD(kstat(kisf),1000) .eq. 0 ) write(*,*) '   ', kstat(kisf), kk_para, kk_CMIP5_anom, kk_TS_pres, kk_tuning

              SELECT CASE (kk_para)
   
              !*************************************************************************************
              CASE(1:4)     ! kk_para=1 -> Linear local               : m = k . TF(z)
                            ! kk_para=2 -> Quadratic local            : m = k . TF(z)  . |TF(z)| 
                            ! kk_para=3 -> Linear with bottom Temp    : m = k . TF_bot  
                            ! kk_para=4 -> Quadratic with bottom Temp : m = k . TF_bot . |TF_bot|
 
                gT  =  gammaT
  
                !=== Present-day melt rates === 
                DO stun=1,10 ! tuning iterations
                 mmm_tot = 0.d0
                 mmm_avg = 0.d0
                 do ii=1,mlondim
                 do jj=1,mlatdim
                   if     ( kk_para .eq. 1 .or. kk_para .eq. 2 ) then
                     zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                   elseif ( kk_para .eq. 3 .or. kk_para .eq. 4 ) then
                     zz1=front_bot_dep_max(kisf) ! deepest point of the entrance
                   endif
                   zz2=-ice_base_topography(ii,jj) ! ice draft depth
                   if ( isfmask(ii,jj) .eq. kisf ) then
                     CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                     T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                     S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                     Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                     ! Melt in m/yr (meters of ice per year), positive if ice ablation
                     if     ( kk_para .eq. 1 .or. kk_para .eq. 3 ) then
                       Melt(ii,jj) = - gT * meltfac * (T0-Tf)              ! Uniform exchange velocity
                     elseif ( kk_para .eq. 2 .or. kk_para .eq. 4 ) then
                       Melt(ii,jj) = - gT * meltfac * (T0-Tf) * abs(T0-Tf) ! Pollard & DeConto 2012
                     endif
                     ! total melt (positive if melting) in Gt/yr
                     mmm_tot = mmm_tot - Melt(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                     mmm_avg = mmm_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                   else
                     Melt(ii,jj) = 0.d0
                   endif
                 enddo
                 enddo
                 ! average melt rate (in m.w.e/yr) :
                 mmm_avg = mmm_tot / mmm_avg
                 !
                 IF ( ABS( ( mmm_tot - target_melt ) / target_melt ) .LT. 1.d-2 ) THEN ! This gT value is satisfactory
                   EXIT
                 ELSEIF ( NINT(SIGN(1.d0,target_melt)) .EQ. NINT(SIGN(1.d0,mmm_tot)) ) THEN ! tune gT and do another iteration
                   gT = gT * target_melt / mmm_tot 
                 ELSEIF ( NINT(SIGN(1.d0,target_melt)) .NE. NINT(SIGN(1.d0,mmm_tot)) .OR. stun .eq. 100 ) THEN ! no possibility of tunning gT to get the correct melt rate
                   gT = -99999.99
                   if ( kstat(kisf) .le. 10 ) write(*,*) 'WARNING1 : sign or non-convergence problem  >>> Fill_value'
                   EXIT
                 ENDIF
                ENDDO ! stun
                if ( gT .gt. 0.0 ) then
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg
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
                else
                  total_melt_pres(kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                  mean_melt_pres (kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                endif

                !=== Future melt rates ===
                if ( gT .gt. 0.0 ) then 
                  mmm_tot = 0.d0
                  mmm_avg = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                   if     ( kk_para .eq. 1 .or. kk_para .eq. 2 ) then
                     zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                   elseif ( kk_para .eq. 3 .or. kk_para .eq. 4 ) then
                     zz1=front_bot_dep_max(kisf) ! deepest point of the entrance
                   endif
                   zz2=-ice_base_topography(ii,jj) ! ice draft depth
                   if ( isfmask(ii,jj) .eq. kisf ) then
                      CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                      T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                      S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                      Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                      ! Melt in m/yr (meters of ice per year), positive if ice ablation
                      if     ( kk_para .eq. 1 .or. kk_para .eq. 3 ) then
                        Melt(ii,jj) = - gT * meltfac * (T0-Tf)              ! Uniform exchange velocity
                      elseif ( kk_para .eq. 2 .or. kk_para .eq. 4 ) then
                        Melt(ii,jj) = - gT * meltfac * (T0-Tf) * abs(T0-Tf) ! Pollard & DeConto 2012
                      endif
                      ! total melt (positive if melting) in Gt/yr
                      mmm_tot = mmm_tot - Melt(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg = mmm_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    else
                      Melt(ii,jj) = 0.d0
                    endif
                    !if ( Fnnnn(ii,jj,kk_para) .gt. 0 .and. ii .eq. 406 .and. jj .eq. 93 ) write(*,*) '@@ ', kk_para, Fmelt(ii,jj,kk_para), Fnnnn(ii,jj,kk_para)
                    !if ( Fnnnn(ii,jj,kk_para) .gt. 0 .and. ii .eq. 406 .and. jj .eq. 93 ) write(*,*) '^^ ', Melt(ii,jj), T0, S0
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg = mmm_tot / mmm_avg
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg
                else
                  total_melt_futu(kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                  mean_melt_futu (kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                endif

                !write(*,364) kisf, kstat(kisf), kk_para, kk_CMIP5_anom, kk_TS_pres, kk_tuning, &
                !& total_melt_pres(kisf,kstat(kisf)), mean_melt_pres(kisf,kstat(kisf)), total_melt_futu(kisf,kstat(kisf)), mean_melt_futu(kisf,kstat(kisf))

              !*****************************************************************************************************************
              CASE(5:7)   ! kk_para=5 -> Linear with mean TF                : m = k . <TF>           (avg over the ice draft)
                          ! kk_para=6 -> Quadratic with mean TF             : m = k . <TF> . |<TF>|  (avg over the ice draft) 
                          ! kk_para=7 -> Quadratic with mixed local/mean TF : m = k . TF(z) .|<TF>|  (avg over the ice draft)

                gT  =  gammaT
  
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
                    TF_tot = TF_tot + (T0-Tf) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                    TF_avg = TF_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                  endif
                enddo
                enddo
                TF_avg = TF_tot / TF_avg
         
                !=== Present-day melt rates === 
                DO stun=1,10 ! tuning iterations
                 mmm_tot = 0.d0
                 mmm_avg = 0.d0
                 do ii=1,mlondim
                 do jj=1,mlatdim
                   if ( isfmask(ii,jj) .eq. kisf ) then
                     ! Melt in m/yr (meters of ice per year), positive if ice ablation
                     if     ( kk_para .eq. 5 ) then
                       Melt(ii,jj) = - gT * meltfac * TF_avg
                     elseif ( kk_para .eq. 6 ) then
                       Melt(ii,jj) = - gT * meltfac * TF_avg * abs(TF_avg) 
                     elseif ( kk_para .eq. 7 ) then
                       zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                       zz2=-ice_base_topography(ii,jj) ! ice draft depth
                       CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                       T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                       S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                       Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                       Melt(ii,jj) = - gT * meltfac * (T0-Tf) * abs(TF_avg)
                     endif
                     ! total melt (positive if melting) in Gt/yr
                     mmm_tot = mmm_tot - Melt(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                     mmm_avg = mmm_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                   else
                     Melt(ii,jj) = 0.d0
                   endif
                 enddo
                 enddo
                 ! average melt rate (in m.w.e/yr) :
                 mmm_avg = mmm_tot / mmm_avg
                 !
                 IF ( ABS( ( mmm_tot - target_melt ) / target_melt ) .LT. 1.d-2 ) THEN ! This gT value is satisfactory
                   EXIT
                 ELSEIF ( NINT(SIGN(1.d0,target_melt)) .EQ. NINT(SIGN(1.d0,mmm_tot)) ) THEN ! tune gT and do another iteration
                   gT = gT * target_melt / mmm_tot 
                 ELSEIF ( NINT(SIGN(1.d0,target_melt)) .NE. NINT(SIGN(1.d0,mmm_tot)) .OR. stun .eq. 100 ) THEN ! no possibility of tunning gT to get the correct melt rate
                   gT = -99999.99
                   if ( kstat(kisf) .le. 10 ) write(*,*) 'WARNING2 : sign or non-convergence problem  >>> Fill_value'
                   EXIT
                 ENDIF
                ENDDO ! stun
                if ( gT .gt. 0.0 ) then
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg
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
                else
                  total_melt_pres(kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                  mean_melt_pres (kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                endif

                !== Mean future Thermal Forcing === 
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
                    TF_tot = TF_tot + (T0-Tf) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                    TF_avg = TF_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                  endif
                enddo
                enddo
                TF_avg = TF_tot / TF_avg

                !=== Future melt rates ===
                if ( gT .gt. 0.0 ) then 
                  mmm_tot = 0.d0
                  mmm_avg = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    if ( isfmask(ii,jj) .eq. kisf ) then
                      ! Melt in m/yr (meters of ice per year), positive if ice ablation
                      if     ( kk_para .eq. 5 ) then
                        Melt(ii,jj) = - gT * meltfac * TF_avg
                      elseif ( kk_para .eq. 6 ) then
                        Melt(ii,jj) = - gT * meltfac * TF_avg * abs(TF_avg) 
                      elseif ( kk_para .eq. 7 ) then
                        zz1=MIN(-ice_base_topography(ii,jj),front_bot_dep_max(kisf)) ! ice draft depth or deepest entrence depth
                        zz2=-ice_base_topography(ii,jj) ! ice draft depth
                        CALL find_z(mdepth,depth,zz1,kkinf,kksup,aainf,aasup)
                        T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                        S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                        Tf = lbd1*S0 + lbd2 + lbd3*zz2  ! Sea water freezing temperature at zz2 (ice draft depth)
                        Melt(ii,jj) = - gT * meltfac * (T0-Tf) * abs(TF_avg)
                      endif
                      ! total melt (positive if melting) in Gt/yr
                      mmm_tot = mmm_tot - Melt(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg = mmm_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    else
                      Melt(ii,jj) = 0.d0
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg = mmm_tot / mmm_avg
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg
                else
                  total_melt_futu(kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                  mean_melt_futu (kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                endif
 
                !write(*,364) kisf, kstat(kisf), kk_para, kk_CMIP5_anom, kk_TS_pres, kk_tuning, &
                !& total_melt_pres(kisf,kstat(kisf)), mean_melt_pres(kisf,kstat(kisf)), total_melt_futu(kisf,kstat(kisf)), mean_melt_futu(kisf,kstat(kisf))

              !*****************************************************************************************************************
              CASE (8:13)   ! kk_para= 8 -> BOX MODEL  2-box forced by T_bot (Olbers & Hellmer 2010; Reese et al. 2018)
                            ! kk_para= 9 -> BOX MODEL  5-box forced by T_bot
                            ! kk_para=10 -> BOX MODEL 10-box forced by T_bot
                            ! kk_para=11 -> BOX MODEL  2-box forced by mean temperature between max bot and avg bot depth at the entrance
                            ! kk_para=12 -> BOX MODEL  5-box forced by mean temperature between max bot and avg bot depth at the entrance
                            ! kk_para=13 -> BOX MODEL 10-box forced by mean temperature between max bot and avg bot depth at the entrance

                CC      = C_PICO        ! Circulation Parameter PICO
                gT      = gammaT_PICO   ! Effective Exchange Velocity PICO
                alphap  = alpha_PICO    ! Thermal expansion coeff PICO
                beta    = beta_PICO     ! Salinity contraction coeff PICO
                rhostar = rhostar_PICO  ! EOS ref Density PICO
                nD      = 0             ! Number Boxes PICO
                if     ( kk_para .eq.  8 .or. kk_para .eq. 11 ) then; nD = nD1
                elseif ( kk_para .eq.  9 .or. kk_para .eq. 12 ) then; nD = nD2
                elseif ( kk_para .eq. 10 .or. kk_para .eq. 13 ) then; nD = nD3; endif

                ALLOCATE( Tbox(nD), Sbox(nD), Mbox(nD) )

                !=== Present-day Melt rates ===
                !- "Ambiant" present-day temperature and salinity
                if ( kk_para .eq.  8 .or. kk_para .eq.  9 .or. kk_para .eq. 10 ) then
                  zzz=front_bot_dep_max(kisf)   ! deepest entrence depth
                  CALL find_z(mdepth,depth,zzz,kkinf,kksup,aainf,aasup)
                  T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                  S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                else
                  zz1=front_bot_dep_avg(kisf)            
                  zz2=front_bot_dep_max(kisf)
                  T0=0.d0; S0=0.d0; kkj=0 
                  do kki=INT(zz1),INT(zz2),5
                    CALL find_z(mdepth,depth,DBLE(kki),kkinf,kksup,aainf,aasup)
                    T0 = T0 + aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                    S0 = S0 + aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                    kkj = kkj + 1
                  enddo
                  T0 = T0 / kkj
                  S0 = S0 / kkj
                endif
                
                DO stun=1,10  ! tuning iterations
                  !- Temerature and salinity in Box #1 :
                  !write(*,*) '===== BOX 1 ====='
                  Tstar = lbd1*S0 + lbd2 + lbd3*Zbox(1,nD) - T0  !NB: Tstar should be < 0
                  g1 = Abox(1,nD) * gT
                  tmp1 = g1 / (CC*rhostar*(beta*S0*meltfac-alphap))
                  xbox = - 0.5*tmp1 + sqrt( (0.5*tmp1)**2 - tmp1*Tstar )
                  Tbox(1) = T0 - xbox
                  !write(*,*) ' Tbox(1) = ', Tbox(1)
                  Sbox(1) = S0 - xbox*S0*meltfac
                  !write(*,*) ' Sbox(1) = ', Sbox(1)
                  qqq = CC*rhostar*(beta*(S0-Sbox(1))-alphap*(T0-Tbox(1)))
                  !write(*,*) ' qqq (Sv) = ', qqq * 1.e-6 / yearinsec
                  Mbox(1) = - gT * meltfac * ( lbd1*Sbox(1) + lbd2 + lbd3*Zbox(1,nD) - Tbox(1) )
                  !write(*,*) ' Mbox(1) = ', Mbox(1), Zbox(1,nD)
                  !
                  !- Temperature and salinity in possible other boxes :
                  DO kk=2,nD
                    !write(*,*) '===== BOX ',kk,' ====='
                    Tstar = lbd1*Sbox(kk-1) + lbd2 + lbd3*Zbox(kk,nD) - Tbox(kk-1)
                    g1  = Abox(kk,nD) * gT
                    g2  = g1 * meltfac
                    xbox = - g1 * Tstar / ( qqq + g1 - g2*lbd1*Sbox(kk-1) )
                    Tbox(kk) = Tbox(kk-1) - xbox
                    !write(*,*) ' Tbox(kk) = ', Tbox(kk)
                    Sbox(kk) = Sbox(kk-1) - xbox*Sbox(kk-1)*meltfac
                    !write(*,*) ' Sbox(kk) = ', Sbox(kk)
                    Mbox(kk) = - gT * meltfac * ( lbd1*Sbox(kk) + lbd2 + lbd3*Zbox(kk,nD) - Tbox(kk) ) 
                    !write(*,*) ' Mbox(kk) = ', Mbox(kk), Zbox(kk,nD)
                  ENDDO
                  !  
                  !- Attribute melt at each node to a box value :
                  mmm_tot = 0.d0
                  mmm_avg = 0.d0
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
                      mmm_tot = mmm_tot - Melt(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg = mmm_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                    else
                      Melt(ii,jj) = 0.d0
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg = mmm_tot / mmm_avg
                  !
                  IF ( ABS( ( mmm_tot - target_melt ) / target_melt ) .LT. 1.d-2 ) THEN ! This gT value is satisfactory
                    EXIT
                  ELSEIF ( NINT(SIGN(1.d0,target_melt)) .EQ. NINT(SIGN(1.d0,mmm_tot)) ) THEN ! tune gT and do another iteration
                    gT = gT * target_melt / mmm_tot 
                  ELSEIF ( NINT(SIGN(1.d0,target_melt)) .NE. NINT(SIGN(1.d0,mmm_tot)) .OR. stun .eq. 100 ) THEN ! no possibility of tunning gT to get the correct melt rate
                    gT = -99999.99
                    if ( kstat(kisf) .le. 10 ) write(*,*) 'WARNING3 : sign or non-convergence problem  >>> Fill_value'
                    EXIT
                  ENDIF
                ENDDO ! stun
                !
                if ( gT .gt. 0.0 .and. sum(Abox(:,nD)) .gt. 0.0 ) then
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg
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
                else
                  total_melt_pres(kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                  mean_melt_pres (kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                endif
             
                !=== Future melt rates ===
                if ( gT .gt. 0.0 .and. sum(Abox(:,nD)) .gt. 0.0 ) then
                  !- "Ambiant" future temperature and salinity
                  if ( kk_para .eq.  8 .or. kk_para .eq.  9 .or. kk_para .eq. 10 ) then
                    zzz=front_bot_dep_max(kisf)   ! deepest entrence depth
                    CALL find_z(mdepth,depth,zzz,kkinf,kksup,aainf,aasup)
                    T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                    S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                  else
                    zz1=front_bot_dep_avg(kisf)            
                    zz2=front_bot_dep_max(kisf)
                    T0=0.d0; S0=0.d0; kkj=0 
                    do kki=INT(zz1),INT(zz2),5
                      CALL find_z(mdepth,depth,DBLE(kki),kkinf,kksup,aainf,aasup)
                      T0 = T0 + aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                      S0 = S0 + aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                      kkj = kkj + 1
                    enddo
                    T0 = T0 / kkj
                    S0 = S0 / kkj
                  endif
                  !
                  !- Temerature and salinity in Box #1 :
                  !write(*,*) '===== BOX 1 ====='
                  Tstar = lbd1*S0 + lbd2 + lbd3*Zbox(1,nD) - T0  !NB: Tstar should be < 0
                  g1 = Abox(1,nD) * gT 
                  tmp1 = g1 / (CC*rhostar*(beta*S0*meltfac-alphap))
                  xbox = - 0.5*tmp1 + sqrt( (0.5*tmp1)**2 - tmp1*Tstar )
                  Tbox(1) = T0 - xbox
                  !write(*,*) ' Tbox(1) = ', Tbox(1)
                  Sbox(1) = S0 - xbox*S0*meltfac
                  !write(*,*) ' Sbox(1) = ', Sbox(1)
                  qqq = CC*rhostar*(beta*(S0-Sbox(1))-alphap*(T0-Tbox(1)))
                  !write(*,*) ' qqq (Sv) = ', qqq * 1.e-6 / yearinsec
                  Mbox(1) = - gT * meltfac * ( lbd1*Sbox(1) + lbd2 + lbd3*Zbox(1,nD) - Tbox(1) )
                  !write(*,*) ' Mbox(1) = ', Mbox(1), Zbox(1,nD)
                  !
                  !- Temperature and salinity in possible other boxes :
                  DO kk=2,nD
                    !write(*,*) '===== BOX ',kk,' ====='
                    Tstar = lbd1*Sbox(kk-1) + lbd2 + lbd3*Zbox(kk,nD) - Tbox(kk-1)
                    g1  = Abox(kk,nD) * gT
                    g2  = g1 * meltfac
                    xbox = - g1 * Tstar / ( qqq + g1 - g2*lbd1*Sbox(kk-1) )
                    Tbox(kk) = Tbox(kk-1) - xbox
                    !write(*,*) ' Tbox(kk) = ', Tbox(kk)
                    Sbox(kk) = Sbox(kk-1) - xbox*Sbox(kk-1)*meltfac
                    !write(*,*) ' Sbox(kk) = ', Sbox(kk)
                    Mbox(kk) = - gT * meltfac * ( lbd1*Sbox(kk) + lbd2 + lbd3*Zbox(kk,nD) - Tbox(kk) ) 
                    !write(*,*) ' Mbox(kk) = ', Mbox(kk), Zbox(kk,nD)
                  ENDDO
                  !  
                  !- Attribute melt at each node to a box value :
                  mmm_tot = 0.d0
                  mmm_avg = 0.d0
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
                      mmm_tot = mmm_tot - Melt(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg = mmm_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    else
                      Melt(ii,jj) = 0.d0
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg = mmm_tot / mmm_avg
                  !
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg
                else
                  total_melt_futu(kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                  mean_melt_futu (kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                endif

                DEALLOCATE( Tbox, Sbox, Mbox )

                !write(*,364) kisf, kstat(kisf), kk_para, kk_CMIP5_anom, kk_TS_pres, kk_tuning, &
                !& total_melt_pres(kisf,kstat(kisf)), mean_melt_pres(kisf,kstat(kisf)), total_melt_futu(kisf,kstat(kisf)), mean_melt_futu(kisf,kstat(kisf))

              !*************************************************************************
              CASE(14:19) ! Plume model (Lazeroms et al. 2018)

                E0      = E0_LAZER      ! Entrainment Coeff LAZER
                GamT    = GT_LAZER      ! Heat Exchange Coeff LAZER
                GefT    = GTS0_LAZER    ! Effective Exchange Coeff LAZER
                if     ( kk_para .eq. 14 .or. kk_para .eq. 17 ) then; TypeL='simple'; kkp=1        ! Zgl and Alpha are found between draft point and the central point
                elseif ( kk_para .eq. 15 .or. kk_para .eq. 18 ) then; TypeL='lazero'; kkp=2        ! original from Lazeroms et al.
                elseif ( kk_para .eq. 16 .or. kk_para .eq. 19 ) then; TypeL='appenB'; kkp=3; endif ! Appendix B in submitted TCD

                ! iterations to find present-day melting and tuning coeff. :
                K = 1.0 ! tuning factor
                DO stun=1,10  ! tuning iterations
                  mmm_tot = 0.d0
                  mmm_avg = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    if ( isfmask(ii,jj) .eq. kisf ) then
                      zzz=-ice_base_topography(ii,jj)
                      ! "Ambient" temperature and salinity at the effective grounding line depth :
                      CALL find_z(mdepth,depth,zGL(ii,jj,kkp),kkinf,kksup,aainf,aasup) 
                      T0 = aainf*T_pres(kkinf)+aasup*T_pres(kksup)
                      S0 = aainf*S_pres(kkinf)+aasup*S_pres(kksup)
                      ! Sea water freezing temperature at the effective grounding line depth :
                      Tf = lbd1*S0 + lbd2 + lbd3*zGL(ii,jj,kkp)  
                      ! Effective heat exchange coefficient :
                      GTS = GamT * ( 0.545 + 3.5e-5 * (T0-Tf)/lbd3 * E0*sin(alpha(ii,jj,kkp))/(GefT+E0*sin(alpha(ii,jj,kkp))) )
                      ! Melt scale :
                      MM = 10. * (T0-Tf)**2 * sqrt( sin(alpha(ii,jj,kkp))/(Cd+E0*sin(alpha(ii,jj,kkp))) ) &
                      &    * sqrt( GTS/(GTS+E0*sin(alpha(ii,jj,kkp))) ) * E0*sin(alpha(ii,jj,kkp))/(GTS+E0*sin(alpha(ii,jj,kkp)))
                      ! Length scale :
                      ll = (T0-Tf)/lbd3 * (x0*GTS+E0*sin(alpha(ii,jj,kkp))) / (x0*(GTS+E0*sin(alpha(ii,jj,kkp))))
                      ! Dimensionless coordinate :
                      X_hat = MAX( 0.d0, MIN( 1.d0, ( zzz - zGL(ii,jj,kkp) ) / ll ) )
                      ! Dimensionless melt curve :
                      M_hat = 0.d0
                      do kk=1,12
                        M_hat = M_hat + pp(kk) * X_hat**(kk-1)
                      enddo
                      ! Melt rate in m/yr:
                      Melt(ii,jj) = - K * MM * M_hat
                      ! total melt (positive if melting) in Gt/yr
                      mmm_tot = mmm_tot - Melt(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg = mmm_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                    else
                      Melt(ii,jj) = 0.d0
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg = mmm_tot / mmm_avg
                  !
                  IF ( ABS( ( mmm_tot - target_melt ) / target_melt ) .LT. 1.d-2 ) THEN ! This gT value is satisfactory
                    EXIT
                  ELSEIF ( NINT(SIGN(1.d0,target_melt)) .EQ. NINT(SIGN(1.d0,mmm_tot)) ) THEN ! tune gT and do another iteration
                    K = K * target_melt / mmm_tot 
                  ELSEIF ( NINT(SIGN(1.d0,target_melt)) .NE. NINT(SIGN(1.d0,mmm_tot)) .OR. stun .eq. 100 ) THEN ! no possibility of tunning gT to get the correct melt rate
                    K = -99999.99
                    if ( kstat(kisf) .le. 10 ) write(*,*) 'WARNING4 : sign or non-convergence problem  >>> Fill_value'
                    EXIT
                  ENDIF
                ENDDO ! stun

                if ( K .gt. 0.0 ) then
                  total_melt_pres(kisf,kstat(kisf)) = mmm_tot
                  mean_melt_pres (kisf,kstat(kisf)) = mmm_avg
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
                else
                  total_melt_pres(kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                  mean_melt_pres (kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                endif

                !=== Future melt rates ===
                if ( K .gt. 0.0 ) then
                  mmm_tot = 0.d0
                  mmm_avg = 0.d0
                  do ii=1,mlondim
                  do jj=1,mlatdim
                    if ( isfmask(ii,jj) .eq. kisf ) then
                      zzz=-ice_base_topography(ii,jj)
                      ! "Ambient" future temperature and salinity at the effective grounding line depth :
                      CALL find_z(mdepth,depth,zGL(ii,jj,kkp),kkinf,kksup,aainf,aasup) 
                      T0 = aainf*T_futu(kkinf)+aasup*T_futu(kksup)
                      S0 = aainf*S_futu(kkinf)+aasup*S_futu(kksup)
                      ! Sea water freezing temperature at the effective grounding line depth :
                      Tf = lbd1*S0 + lbd2 + lbd3*zGL(ii,jj,kkp)  
                      ! Effective heat exchange coefficient :
                      GTS = GamT * ( 0.545 + 3.5e-5 * (T0-Tf)/lbd3 * E0*sin(alpha(ii,jj,kkp))/(GefT+E0*sin(alpha(ii,jj,kkp))) )
                      ! Melt scale :
                      MM = 10. * (T0-Tf)**2 * sqrt( sin(alpha(ii,jj,kkp))/(Cd+E0*sin(alpha(ii,jj,kkp))) ) &
                      &    * sqrt( GTS/(GTS+E0*sin(alpha(ii,jj,kkp))) ) * E0*sin(alpha(ii,jj,kkp))/(GTS+E0*sin(alpha(ii,jj,kkp)))
                      ! Length scale :
                      ll = (T0-Tf)/lbd3 * (x0*GTS+E0*sin(alpha(ii,jj,kkp))) / (x0*(GTS+E0*sin(alpha(ii,jj,kkp))))
                      ! Dimensionless coordinate :
                      X_hat = MAX( 0.d0, MIN( 1.d0, ( zzz - zGL(ii,jj,kkp) ) / ll ) )
                      ! Dimensionless melt curve :
                      M_hat = 0.d0
                      do kk=1,12
                        M_hat = M_hat + pp(kk) * X_hat**(kk-1)
                      enddo
                      ! Melt rate in m/yr:
                      Melt(ii,jj) = - K * MM * M_hat
                      ! total melt (positive if melting) in Gt/yr
                      mmm_tot = mmm_tot - Melt(ii,jj) * dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9  * rhoi_SI / rhofw_SI
                      mmm_avg = mmm_avg + dlon * dlat * dcos(lat(jj)*deg2rad) * RT**2 * deg2rad**2 * 1.d-9
                      ! store mean melt pattern for each parametrization :
                      Fmelt(ii,jj,kk_para) = Fmelt(ii,jj,kk_para) - Melt(ii,jj)
                      Fnnnn(ii,jj,kk_para) = Fnnnn(ii,jj,kk_para) + 1
                    else
                      Melt(ii,jj) = 0.d0
                    endif
                  enddo
                  enddo
                  ! average melt rate (in m.w.e/yr) :
                  mmm_avg = mmm_tot / mmm_avg
                  !
                  total_melt_futu(kisf,kstat(kisf)) = mmm_tot
                  mean_melt_futu (kisf,kstat(kisf)) = mmm_avg
                else
                  total_melt_futu(kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                  mean_melt_futu (kisf,kstat(kisf)) = NF90_FILL_DOUBLE
                endif ! if ( K .gt. 0.0 )

                !write(*,364) kisf, kstat(kisf), kk_para, kk_CMIP5_anom, kk_TS_pres, kk_tuning, &
                !& total_melt_pres(kisf,kstat(kisf)), mean_melt_pres(kisf,kstat(kisf)), total_melt_futu(kisf,kstat(kisf)), mean_melt_futu(kisf,kstat(kisf))

              !*********************************
              CASE DEFAULT

                write(*,*) 'WRONG PARAMETERIZATION CHOICE !!!!!!!!!'
                STOP

              !*********************************
              END SELECT

           ENDDO  !! kk_para

        ENDDO !! kk_CMIP5_anom

      ENDDO !! kk_TS_pres
     ENDIF !! IF ( target_melt .gt. 1.e-6 )

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
   839 FORMAT('melt_ensemble_',i3.3,'.nc')   
   write(*,*) 'Creating ', TRIM(file_out)
   
   !status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM)
   status = NF90_CREATE(TRIM(file_out),or(NF90_CLOBBER,NF90_64BIT_OFFSET),fidM)
   call erreur(status,.TRUE.,'create output file')
  
   status = NF90_DEF_DIM(fidM,"lon",mlondim,dimID_londim); call erreur(status,.TRUE.,"def_dimID_londim")
   status = NF90_DEF_DIM(fidM,"lat",mlatdim,dimID_latdim); call erreur(status,.TRUE.,"def_dimID_latdim")
   status = NF90_DEF_DIM(fidM,"Nisf",mNisf,dimID_Nisf); call erreur(status,.TRUE.,"def_dimID_Nisf")
   status = NF90_DEF_DIM(fidM,"mstat",mtot,dimID_mstat); call erreur(status,.TRUE.,"def_dimID_mstat")
   status = NF90_DEF_DIM(fidM,"mpara",nn_para,dimID_mpara); call erreur(status,.TRUE.,"def_dimID_mpara")
   
   status = NF90_DEF_VAR(fidM,"lon",NF90_FLOAT,(/dimID_londim/),lon_ID)
   call erreur(status,.TRUE.,"def_var_lon_ID")
   status = NF90_DEF_VAR(fidM,"lat",NF90_FLOAT,(/dimID_latdim/),lat_ID)
   call erreur(status,.TRUE.,"def_var_lat_ID")
   status = NF90_DEF_VAR(fidM,"Pres_Melt",NF90_FLOAT,(/dimID_londim,dimID_latdim,dimID_mpara/),Pmelt_ID)
   call erreur(status,.TRUE.,"def_var_Pmelt_ID")
   status = NF90_DEF_VAR(fidM,"Futu_Melt",NF90_FLOAT,(/dimID_londim,dimID_latdim,dimID_mpara/),Fmelt_ID)
   call erreur(status,.TRUE.,"def_var_Fmelt_ID")
   status = NF90_DEF_VAR(fidM,"total_melt_pres",NF90_DOUBLE,(/dimID_Nisf,dimID_mstat/),total_melt_pres_ID)
   call erreur(status,.TRUE.,"def_var_total_melt_pres_ID")
   status = NF90_DEF_VAR(fidM,"mean_melt_pres",NF90_DOUBLE,(/dimID_Nisf,dimID_mstat/),mean_melt_pres_ID)
   call erreur(status,.TRUE.,"def_var_mean_melt_pres_ID")
   status = NF90_DEF_VAR(fidM,"total_melt_futu",NF90_DOUBLE,(/dimID_Nisf,dimID_mstat/),total_melt_futu_ID)
   call erreur(status,.TRUE.,"def_var_total_melt_futu_ID")
   status = NF90_DEF_VAR(fidM,"mean_melt_futu",NF90_DOUBLE,(/dimID_Nisf,dimID_mstat/),mean_melt_futu_ID)
   call erreur(status,.TRUE.,"def_var_mean_melt_futu_ID")
   status = NF90_DEF_VAR(fidM,"index_para",NF90_SHORT,(/dimID_Nisf,dimID_mstat/),index_para_ID)
   call erreur(status,.TRUE.,"def_var_index_para_ID")
   status = NF90_DEF_VAR(fidM,"index_CMIP",NF90_SHORT,(/dimID_Nisf,dimID_mstat/),index_CMIP_ID)
   call erreur(status,.TRUE.,"def_var_index_CMIP_ID")
   status = NF90_DEF_VAR(fidM,"index_WOA",NF90_SHORT,(/dimID_Nisf,dimID_mstat/),index_WOA_ID)
   call erreur(status,.TRUE.,"def_var_index_WOA_ID")
   status = NF90_DEF_VAR(fidM,"index_tun",NF90_SHORT,(/dimID_Nisf,dimID_mstat/),index_tun_ID)
   call erreur(status,.TRUE.,"def_var_index_tun_ID")
   
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
   status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"_FillValue",NF90_FILL_DOUBLE)
   call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"long_name","Total present cavity melt flux")
   call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,total_melt_pres_ID,"title","Present Melt Flux")
   call erreur(status,.TRUE.,"put_att_total_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"units","Gt/yr")
   call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"_FillValue",NF90_FILL_DOUBLE)
   call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"long_name","Total future cavity melt flux")
   call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,total_melt_futu_ID,"title","Future Melt Flux")
   call erreur(status,.TRUE.,"put_att_total_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"units","m/yr")
   call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"_FillValue",NF90_FILL_DOUBLE)
   call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"long_name","Present mean melt rate over the cavity")
   call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_pres_ID,"title","Present Melt Rate")
   call erreur(status,.TRUE.,"put_att_mean_melt_pres_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"units","m/yr")
   call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"_FillValue",NF90_FILL_DOUBLE)
   call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"long_name","Future mean melt rate over the cavity")
   call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
   status = NF90_PUT_ATT(fidM,mean_melt_futu_ID,"title","Future Melt Rate")
   call erreur(status,.TRUE.,"put_att_mean_melt_futu_ID")
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
   status = NF90_PUT_ATT(fidM,index_tun_ID,"units","-")
   call erreur(status,.TRUE.,"put_att_index_tun_ID")
   status = NF90_PUT_ATT(fidM,index_tun_ID,"long_name","index defining the sampling of Rignot (2013)'s range of uncertainty")
   call erreur(status,.TRUE.,"put_att_index_tun_ID")
   status = NF90_PUT_ATT(fidM,index_tun_ID,"title","index_tun")
   call erreur(status,.TRUE.,"put_att_index_tun_ID")
   
   status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using calculate_melt_rate_ensemble.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
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
   status = NF90_PUT_VAR(fidM,total_melt_futu_ID,total_melt_futu); call erreur(status,.TRUE.,"var_total_melt_futu_ID")
   status = NF90_PUT_VAR(fidM,mean_melt_futu_ID,mean_melt_futu); call erreur(status,.TRUE.,"var_mean_melt_futu_ID")
   status = NF90_PUT_VAR(fidM,index_para_ID,index_para); call erreur(status,.TRUE.,"var_index_para_ID")
   status = NF90_PUT_VAR(fidM,index_CMIP_ID,index_CMIP); call erreur(status,.TRUE.,"var_index_CMIP_ID")
   status = NF90_PUT_VAR(fidM,index_WOA_ID,index_WOA); call erreur(status,.TRUE.,"var_index_WOA_ID")
   status = NF90_PUT_VAR(fidM,index_tun_ID,index_tun); call erreur(status,.TRUE.,"var_index_tun_ID")
   
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
