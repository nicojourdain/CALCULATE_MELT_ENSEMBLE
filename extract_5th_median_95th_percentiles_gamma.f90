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
&                  kk_K, nn_K, mstat, mpara2, DeltaT_ID

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
&                                          s_anom_5, t_anom_5, s_anom_4, t_anom_4, Kcoef_pct, DeltaT_out

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
&                maxDeltaT, ctau, cr1, alpha_LAZER, beta_LAZER, ttt, ttt_SI

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

file_in_stat        = 'melt_ensemble_statistics_ANT3.nc'

nn_K       = 3
nn_para    = 18

ALLOCATE( perct(nn_K), Kcoef_pct(nn_para,nn_K) )
perct = (/ 0.050, 0.500, 0.950 /)

!----------------------------------------------------------
! Read statistics from first calculation (with DeltaT = 0) :
 
write(*,*) 'Reading ', TRIM(file_in_stat)
 
status = NF90_OPEN(TRIM(file_in_stat),0,fidSTAT); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidSTAT,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
status = NF90_INQ_DIMID(fidSTAT,"mstat",dimID_mstat); call erreur(status,.TRUE.,"inq_dimID_mstat")
status = NF90_INQ_DIMID(fidSTAT,"mpara",dimID_mpara); call erreur(status,.TRUE.,"inq_dimID_mpara")
 
status = NF90_INQUIRE_DIMENSION(fidSTAT,dimID_Nisf,len=mNisf); call erreur(status,.TRUE.,"inq_dim_Nisf")
status = NF90_INQUIRE_DIMENSION(fidSTAT,dimID_mstat,len=mstat); call erreur(status,.TRUE.,"inq_dim_mstat")
  
ALLOCATE(  Kcoef(mNisf,mstat)  ) 
ALLOCATE(  index_para(mNisf,mstat)  ) 

status = NF90_INQ_VARID(fidSTAT,"Kcoef",Kcoef_ID); call erreur(status,.TRUE.,"inq_Kcoef_ID")
status = NF90_INQ_VARID(fidSTAT,"index_para",index_para_ID); call erreur(status,.TRUE.,"inq_index_para_ID")

status = NF90_GET_VAR(fidSTAT,Kcoef_ID,Kcoef); call erreur(status,.TRUE.,"getvar_Kcoef")
status = NF90_GET_VAR(fidSTAT,index_para_ID,index_para); call erreur(status,.TRUE.,"getvar_index_para")

status = NF90_CLOSE(fidSTAT); call erreur(status,.TRUE.,"close_file")

!----- calculate percentiles for each parameterization :

DO kk_para = 1, nn_para

  write(*,*) '    for param. ', kk_para  

  nnp=0
  do ii=10,75   !! NB: only take coeff of known ice shelves :
  do jj=1,mstat
    if ( index_para(ii,jj) .eq. kk_para .and. Kcoef(ii,jj) .ne. NF90_FILL_FLOAT ) nnp=nnp+1
  enddo
  enddo
  write(*,*) '          over ', nnp, ' values'
  
  ALLOCATE( Kcoef2(nnp), indx2(nnp) )
  nnp=0
  do ii=10,75   !! NB: only take coeff of known ice shelves :
  do jj=1,mstat
    if ( index_para(ii,jj) .eq. kk_para .and. Kcoef(ii,jj) .ne. NF90_FILL_FLOAT ) then
      nnp=nnp+1
      Kcoef2(nnp) = Kcoef(ii,jj)
      indx2(nnp) = nnp
    endif
  enddo
  enddo

  CALL hpsort_eps_epw(nnp,Kcoef2,indx2,1.d-20) ! this sorts Kcoef2
  do kk_K=1,nn_K
    Kcoef_pct(kk_para,kk_K) = Kcoef2(NINT(perct(kk_K)*nnp))
    write(*,*) 'kk_para, kk_K, NINT(perct(kk_K)*nnp), Kcoef_pct = ', kk_para, kk_K, NINT(perct(kk_K)*nnp), Kcoef_pct(kk_para,kk_K)
  enddo

  DEALLOCATE( Kcoef2, indx2 )

ENDDO
!--------------------------------------------------------

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
lambda3   =   -7.53e-8 * grav_SI * rhosw_SI ! Liquidus pressure coefficient (K/m)

!- conversions to Elmer units:
rhoi    = rhoi_SI  / (1.0e6*yearinsec**2)
rhosw   = rhosw_SI / (1.0e6*yearinsec**2)
rhofw   = rhofw_SI / (1.0e6*yearinsec**2)
Lf      = Lf_SI   * yearinsec**2
cpw     = cpw_SI  * yearinsec**2
gravity = grav_SI * yearinsec**2 * (-1.0)

meltfac    = rhosw_SI * cpw_SI / ( rhoi_SI * Lf_SI )

!--------------------------------------------------------
! Express the exchange velocities in m/yr (according to formulas in the paper):
!
!     - Linear already in m/yr
!     - Quadratic to correect because no meltfa^2 in my scripts...
!     - Box already in m/yr

Kcoef_pct(2,:) = Kcoef_pct(2,:) / meltfac
Kcoef_pct(4:6,:) = Kcoef_pct(4:6,:) / meltfac
  
101 FORMAT(a,' ',e10.3,' & ',e10.3,' & ',e10.3,' \\ ')

write(*,101) '                          & Linear local                               &', Kcoef_pct(1,1), Kcoef_pct(1,2), Kcoef_pct(1,3)
write(*,101) '                          & Quadratic local                            &', Kcoef_pct(2,1), Kcoef_pct(2,2), Kcoef_pct(2,3)
write(*,101) 'SIMPLE                    & Linear bottom                              &', Kcoef_pct(3,1), Kcoef_pct(3,2), Kcoef_pct(3,3)
write(*,101) '$\gamma_{\mbox{\tiny T}}$ & Quadratic bottom                           &', Kcoef_pct(4,1), Kcoef_pct(4,2), Kcoef_pct(4,3)
write(*,101) '                          & Quadratic mean                             &', Kcoef_pct(5,1), Kcoef_pct(5,2), Kcoef_pct(5,3)
write(*,101) '                          & Quadratic non-local                        &', Kcoef_pct(6,1), Kcoef_pct(6,2), Kcoef_pct(6,3)
write(*,*) '\hline'
write(*,101) '                          & 2-box discrete                             &', Kcoef_pct(7,1), Kcoef_pct(7,2), Kcoef_pct(7,3)
write(*,101) '                          & 10-box discrete                            &', Kcoef_pct(8,1), Kcoef_pct(8,2), Kcoef_pct(8,3)
write(*,101) 'BOX                       & 5-box PICO ~~$C=C_0$                       &', Kcoef_pct(9,1), Kcoef_pct(9,2), Kcoef_pct(9,3)
write(*,101) '$\gamma_{\mbox{\tiny T}}$ & 10-box PICO ~~$C=C_0$                      &', Kcoef_pct(10,1), Kcoef_pct(10,2), Kcoef_pct(10,3)
write(*,101) '                          & 5-box PICO ~~$C=a\gamma_{\mbox{\tiny T}}$  &', Kcoef_pct(11,1), Kcoef_pct(11,2), Kcoef_pct(11,3)
write(*,101) '                          & 10-box PICO ~~$C=a\gamma_{\mbox{\tiny T}}$ &', Kcoef_pct(12,1), Kcoef_pct(12,2), Kcoef_pct(12,3)
write(*,*) '\hline'
write(*,101) '                          & Polynomial, plume from deepest GL          &', Kcoef_pct(13,1), Kcoef_pct(13,2), Kcoef_pct(13,3)
write(*,101) '                          & Polynomial, multi-directional plumes       &', Kcoef_pct(14,1), Kcoef_pct(14,2), Kcoef_pct(14,3)
write(*,101) 'PLUME                     & Polynomial, multi-GL, local slope          &', Kcoef_pct(15,1), Kcoef_pct(15,2), Kcoef_pct(15,3)
write(*,101) '$\Gamma_{\mbox{\tiny T}}$ & Analytical, plume from deepest GL          &', Kcoef_pct(16,1), Kcoef_pct(16,2), Kcoef_pct(16,3)
write(*,101) '                          & Analytical, multi-directional plumes       &', Kcoef_pct(17,1), Kcoef_pct(17,2), Kcoef_pct(17,3)
write(*,101) '                          & Analytical, multi-GL, local slope          &', Kcoef_pct(18,1), Kcoef_pct(18,2), Kcoef_pct(18,3)

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

