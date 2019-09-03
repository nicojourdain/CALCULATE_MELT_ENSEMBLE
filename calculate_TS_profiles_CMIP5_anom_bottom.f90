program modif                                         

USE netcdf                                            

IMPLICIT NONE                                         

INTEGER :: fidT, dimID_depth, dimID_lon, dimID_lat, dimID_nbounds, mdepth, mnbounds, t_anom_ID, s_anom_ID, ks, kr, &
&          depth_ID, lon_ID, lat_ID, fidS, iCMIP, jCMIP, kCMIP, j50S, kisf, strlen, dimID_strlen, dimID_Nmod,&
&          mean_s_anom_4_ID, mean_t_anom_4_ID, mean_s_anom_5_ID, mean_t_anom_5_ID, kmod, Nmod, modnam_ID, dimID_x, kt, &
&          fidISF, status, dimID_Nisf, Nisf, front_max_lat_ID, fidM, dimID_y, fidB,                                    &
&          front_min_lat_ID, front_max_lon_ID, front_min_lon_ID, front_ice_dep_avg_ID, front_ice_dep_min_ID, fidA, k01,&
&          front_bot_dep_max_ID, front_bot_dep_avg_ID, mx, my, nav_lat_ID, nav_lon_ID, var_xxx_ID,  &
&          fidZ, dimID_z, melt_isf_ID, ii, jj, kk, NN

CHARACTER(LEN=150) :: file_isf, file_out, file_in_his, file_in_rcp, moddir

CHARACTER(LEN=4) :: varnam

CHARACTER(LEN=14),ALLOCATABLE,DIMENSION(:) :: modnam

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Ntot

REAL*4,ALLOCATABLE,DIMENSION(:) :: front_max_lat, front_min_lat, front_max_lon, front_min_lon, front_ice_dep_avg,    &
&                                  front_ice_dep_min, front_bot_dep_avg, front_bot_dep_max, melt_isf

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: nav_lat, nav_lon, tmp_lon

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: mean_t_anom_4, mean_s_anom_4, &
&                                    mean_t_anom_5, mean_s_anom_5, &
&                                    T_his, T_rcp, S_his, S_rcp
 
REAL*4,ALLOCATABLE,DIMENSION(:,:) :: var_xxx

REAL*4 :: Tmisval, Smisval, dif, mindif, front_inf, front_sup, &
&         lon_box_1, lon_box_2, lon_box_3, lat_box_1, lat_box_2, lat_box_3
 
file_isf = 'RTopo-2.0.1_30sec_ice_shelf_mask_ANT3.nc'
file_out = 'CMIP5_anom_TS_bottom_per_ice_shelf_ANT3.nc'

! Boxes :
! NB: 5.0 x 1.75 is the effective resolution at 70°S for a model of 1° resolution in longitude (assuming 5 delta X and a Mercator grid)
lon_box_1 =  5.0 ; lat_box_1 = 1.75 ! Continental shelf only
lon_box_2 = 10.0 ; lat_box_2 = 3.50 ! Continental shelf only

! CMIP5 models :
!Nmod=33
Nmod=31
ALLOCATE( modnam(Nmod) )

!moddir='/store/njourd/CMIP5_ANOM/REGRID_OCE'
moddir='/store/njourd/CMIP5_ANOM/SOUTH_EXTRACTION/LAST_LEVEL_REMOVED'

modnam=(/ 'ACCESS1-0     ', &
&         'ACCESS1-3     ', &
&         'BNU-ESM       ', &
&         'CCSM4         ', &
&         'CESM1-BGC     ', &
&         'CESM1-CAM5    ', &
&         'CESM1-WACCM   ', &
&         'CMCC-CESM     ', &
&         'CMCC-CM       ', &
&         'CMCC-CMS      ', &
&         'CNRM-CM5      ', &
&         'CSIRO-Mk3-6-0 ', &
&         'CanESM2       ', &
&         'FGOALS-g2     ', &
&         'FIO_ESM       ', &
&         'GFDL-CM3      ', &
&         'GFDL-ESM2G    ', &
&         'GFDL-ESM2M    ', &
&         'HadGEM2-CC    ', &
&         'HadGEM2-ES    ', &
&         'IPSL-CM5A-LR  ', &
&         'IPSL-CM5A-MR  ', &
&         'IPSL-CM5B-LR  ', &
&         'MIROC-ESM     ', &
&         'MIROC-ESM-CHEM', &
!&         'MIROC5        ', &
&         'MPI-ESM-LR    ', &
&         'MPI-ESM-MR    ', &
!&         'MRI-CGCM3     ', &
&         'NorESM1-M     ', &
&         'NorESM1-ME    ', &
&         'bcc-csm1-1    ', &
&         'inmcm4        '   /)

strlen=LEN(modnam)

!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
! Read model T,S 

ALLOCATE( mean_t_anom_4(Nmod,Nisf), mean_s_anom_4(Nmod,Nisf) )
ALLOCATE( mean_t_anom_5(Nmod,Nisf), mean_s_anom_5(Nmod,Nisf) )

mean_t_anom_4(:,:) = 0.e0
mean_s_anom_4(:,:) = 0.e0
mean_t_anom_5(:,:) = 0.e0
mean_s_anom_5(:,:) = 0.e0

DO kmod=1,Nmod
   
    write(*,*) '####################################################'
    write(*,*) '##           ', modnam(kmod), '           ##'
    write(*,*) '####################################################'
    
    write(file_in_his,201) TRIM(moddir), TRIM(modnam(kmod))
    201 FORMAT(a,'/SOUTH_',a,'_hist_1989-2009_bot.nc')

    write(file_in_rcp,202) TRIM(moddir), TRIM(modnam(kmod))
    202 FORMAT(a,'/SOUTH_',a,'_rcp85_2080-2100_bot.nc')
                                                           
    !-------------------------------------------                   
    ! Read historical temperature climatology : 
 
    varnam='Tbot'
 
    write(*,*) 'Reading ', TRIM(varnam), ' in ', TRIM(file_in_his)

    status = NF90_OPEN(TRIM(file_in_his),0,fidA)          
    call erreur(status,.TRUE.,"read historical climatology") 

    if ( kmod .eq. 1 ) then

         status = NF90_INQ_DIMID(fidA,"x",dimID_x)
         call erreur(status,.TRUE.,"inq_dimID_x")
         status = NF90_INQ_DIMID(fidA,"y",dimID_y)
         call erreur(status,.TRUE.,"inq_dimID_y")

         status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx)
         call erreur(status,.TRUE.,"inq_dim_x")
         status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my)
         call erreur(status,.TRUE.,"inq_dim_y")
 
         write(*,*) '(mx,my) = ', mx, my

         ALLOCATE(  nav_lat(mx,my)  ) 
         ALLOCATE(  nav_lon(mx,my)  ) 
         ALLOCATE(  tmp_lon(mx,my)  ) 

         status = NF90_INQ_VARID(fidA,"nav_lon",nav_lon_ID)
         call erreur(status,.TRUE.,"inq_nav_lon_ID")
         status = NF90_INQ_VARID(fidA,"nav_lat",nav_lat_ID)
         call erreur(status,.TRUE.,"inq_nav_lat_ID")

         status = NF90_GET_VAR(fidA,nav_lon_ID,nav_lon)
         call erreur(status,.TRUE.,"getvar_nav_lon")
         status = NF90_GET_VAR(fidA,nav_lat_ID,nav_lat)
         call erreur(status,.TRUE.,"getvar_nav_lat")

         ! find CMIP5's 50degS :
         mindif = 1.e4
         do jCMIP=1,my
           dif = abs( nav_lat(500,jCMIP) + 50.0 )
           if ( dif .lt. mindif ) then
             mindif = dif
             j50S = jCMIP
           endif
         enddo
         write(*,*) 'Index corresponding to 50degS : ', j50S, '      check: ', nav_lat(500,j50S), nav_lat(1,j50S)

         ALLOCATE(  T_his(mx,j50S)  )
         ALLOCATE(  T_rcp(mx,j50S)  )
         ALLOCATE(  S_his(mx,j50S)  )
         ALLOCATE(  S_rcp(mx,j50S)  )
         ALLOCATE(  Ntot(mx,j50S)  )

    endif

    ALLOCATE(  var_xxx (mx,j50S)  )

    status = NF90_INQ_VARID(fidA,TRIM(varnam),var_xxx_ID)
    call erreur(status,.TRUE.,"inq_var_xxx_ID")

    status = NF90_GET_VAR(fidA,var_xxx_ID,var_xxx,(/ 1, 1 /),(/mx, j50S /))
    call erreur(status,.TRUE.,"getvar_var_his")
                     
    status = NF90_GET_ATT(fidA,var_xxx_ID,"missing_value",Tmisval)
    call erreur(status,.TRUE.,"get_att_missval_ID")

    status = NF90_CLOSE(fidA)                      
    call erreur(status,.TRUE.,"fin_lecture_historical")     

    Ntot(:,:)  = 0
    T_his(:,:) = 0.e0
    do ii=1,mx
    do jj=1,j50S
        if ( var_xxx(ii,jj) .ge. 268.15 .and. var_xxx(ii,jj) .lt. 306.15 ) then
          T_his(ii,jj) = T_his(ii,jj) + var_xxx(ii,jj)
          Ntot(ii,jj) = Ntot(ii,jj) + 1
        endif
    enddo
    enddo
    do ii=1,mx
    do jj=1,j50S
      if ( Ntot(ii,jj) .gt. 0 ) then
        T_his(ii,jj) = T_his(ii,jj) / Ntot(ii,jj) - 273.15
      else
        T_his(ii,jj) = NF90_FILL_FLOAT
      endif
    enddo
    enddo
    write(*,*) '     T_his interval ', MINVAL(T_his), ' : ', MAXVAL(T_his) 
 
    DEALLOCATE( var_xxx )
                                  
    !---------------------------------------                   
    ! Read rcp85 temperature climatology : 

    write(*,*) 'Reading ', TRIM(varnam), ' in ', TRIM(file_in_rcp)

    status = NF90_OPEN(TRIM(file_in_rcp),0,fidB)
    call erreur(status,.TRUE.,"read rcp85 climatology")

    ALLOCATE(  var_xxx (mx,j50S)  )
                            
    status = NF90_INQ_VARID(fidB,TRIM(varnam),var_xxx_ID)
    call erreur(status,.TRUE.,"inq_var_xxx_ID")

    status = NF90_GET_VAR(fidB,var_xxx_ID,var_xxx,(/ 1, 1 /),(/mx, j50S /))
    call erreur(status,.TRUE.,"getvar_var_rcp")

    status = NF90_CLOSE(fidB)
    call erreur(status,.TRUE.,"fin_lecture_rcp85")

    Ntot(:,:)  = 0
    T_rcp(:,:) = 0.e0
    do ii=1,mx
    do jj=1,j50S
        if ( var_xxx(ii,jj) .ge. 268.15 .and. var_xxx(ii,jj) .lt. 306.15 ) then
          T_rcp(ii,jj) = T_rcp(ii,jj) + var_xxx(ii,jj)
          Ntot(ii,jj) = Ntot(ii,jj) + 1
        endif
    enddo
    enddo
    do ii=1,mx
    do jj=1,j50S
      if ( Ntot(ii,jj) .gt. 0 ) then
        T_rcp(ii,jj) = T_rcp(ii,jj) / Ntot(ii,jj) - 273.15
      else
        T_rcp(ii,jj) = NF90_FILL_FLOAT
      endif
    enddo
    enddo
    write(*,*) '     T_rcp interval ', MINVAL(T_rcp), ' : ', MAXVAL(T_rcp) 

    DEALLOCATE( var_xxx )

    !---------------------------------------                      
    ! Read historical salinity climatology : 
 
    varnam='Sbot'
 
    write(*,*) 'Reading ', TRIM(varnam), ' in ', TRIM(file_in_his)

    status = NF90_OPEN(TRIM(file_in_his),0,fidA)          
    call erreur(status,.TRUE.,"read historical climatology") 

    ALLOCATE(  var_xxx (mx,j50S)  )

    status = NF90_INQ_VARID(fidA,TRIM(varnam),var_xxx_ID)
    call erreur(status,.TRUE.,"inq_var_xxx_ID")

    status = NF90_GET_VAR(fidA,var_xxx_ID,var_xxx,(/ 1, 1 /),(/mx, j50S /))
    call erreur(status,.TRUE.,"getvar_var_his")
                                                  
    status = NF90_GET_ATT(fidA,var_xxx_ID,"missing_value",Smisval)
    call erreur(status,.TRUE.,"get_att_missval_ID")

    status = NF90_CLOSE(fidA)                      
    call erreur(status,.TRUE.,"fin_lecture_historical")     

    Ntot(:,:)  = 0
    S_his(:,:) = 0.e0
    do ii=1,mx
    do jj=1,j50S
        if ( var_xxx(ii,jj) .ge. 5.0 .and. var_xxx(ii,jj) .lt. 50.0 ) then
          S_his(ii,jj) = S_his(ii,jj) + var_xxx(ii,jj)
          Ntot(ii,jj) = Ntot(ii,jj) + 1
        endif
    enddo
    enddo
    do ii=1,mx
    do jj=1,j50S
      if ( Ntot(ii,jj) .gt. 0 ) then
        S_his(ii,jj) = S_his(ii,jj) / Ntot(ii,jj)
      else
        S_his(ii,jj) = NF90_FILL_FLOAT
      endif
    enddo
    enddo
    write(*,*) '     S_his interval ', MINVAL(S_his), ' : ', MAXVAL(S_his) 

    DEALLOCATE( var_xxx )
                                  
    !---------------------------------------                   
    ! Read rcp85 salinity climatology : 

    write(*,*) 'Reading ', TRIM(varnam), ' in ', TRIM(file_in_rcp)

    status = NF90_OPEN(TRIM(file_in_rcp),0,fidB)
    call erreur(status,.TRUE.,"read rcp85 climatology")

    ALLOCATE(  var_xxx (mx,j50S)  )
                            
    status = NF90_INQ_VARID(fidB,TRIM(varnam),var_xxx_ID)
    call erreur(status,.TRUE.,"inq_var_xxx_ID")

    status = NF90_GET_VAR(fidB,var_xxx_ID,var_xxx,(/ 1, 1 /),(/mx, j50S /))
    call erreur(status,.TRUE.,"getvar_var_rcp")

    status = NF90_CLOSE(fidB)
    call erreur(status,.TRUE.,"fin_lecture_rcp85")

    Ntot(:,:)  = 0
    S_rcp(:,:) = 0.e0
    do ii=1,mx
    do jj=1,j50S
        if ( var_xxx(ii,jj) .ge. 5.0 .and. var_xxx(ii,jj) .lt. 50.0 ) then
          S_rcp(ii,jj) = S_rcp(ii,jj) + var_xxx(ii,jj)
          Ntot(ii,jj) = Ntot(ii,jj) + 1
        endif
    enddo
    enddo
    do ii=1,mx
    do jj=1,j50S
      if ( Ntot(ii,jj) .gt. 0 ) then
        S_rcp(ii,jj) = S_rcp(ii,jj) / Ntot(ii,jj)
      else
        S_rcp(ii,jj) = NF90_FILL_FLOAT
      endif
    enddo
    enddo
    write(*,*) '     S_rcp interval ', MINVAL(S_rcp), ' : ', MAXVAL(S_rcp) 

    DEALLOCATE( var_xxx )

    !------------------------------------------------------
    ! Calculate T,S profiles for individual ice shelves :

    write(*,*) 'Calculating T,S profiles for individual ice shelves...'
   
    DO kisf=2,Nisf

     IF (       front_min_lon(kisf) .ne. 180.0 .and. front_max_lon(kisf) .ne. -180.0 &      ! to remove non-attributed ice shelves
     &    .and. front_min_lat(kisf) .ne. 90.00 .and. front_max_lat(kisf) .ne. -90.00 ) THEN

      if ( kisf .ne. 10 ) then ! not Ross IS
        tmp_lon(:,:) = nav_lon(:,:)
      else ! ROSS
        where ( nav_lon(:,:) .lt. 0.0 )
          tmp_lon(:,:) = nav_lon(:,:) + 360.0
        elsewhere
          tmp_lon(:,:) = nav_lon(:,:)
        endwhere
        front_max_lon(kisf) = front_max_lon(kisf) + 360.0
      endif
   
      write(*,*) '@@', kisf, front_max_lon(kisf)
 
      ! closest to the ice shelf front within 1 deg in lon and lat :
      NN = 0.e0
      do iCMIP=1,mx
      do jCMIP=1,j50S
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
        if (       tmp_lon(iCMIP,jCMIP) .ge. front_inf  &
        &    .and. tmp_lon(iCMIP,jCMIP) .le. front_sup  &
        &    .and. nav_lat(iCMIP,jCMIP) .ge. front_min_lat(kisf) - lat_box_1  &
        &    .and. nav_lat(iCMIP,jCMIP) .le. front_max_lat(kisf) + lat_box_1  ) then
           if (       T_his(iCMIP,jCMIP) .ge.  -5.0 .and. T_his(iCMIP,jCMIP) .lt.  33.0 &
           &    .and. S_his(iCMIP,jCMIP) .ge.   0.0 .and. S_his(iCMIP,jCMIP) .lt.  60.0 ) then
             mean_t_anom_4(kmod,kisf) = mean_t_anom_4(kmod,kisf) + T_rcp(iCMIP,jCMIP) - T_his(iCMIP,jCMIP)
             mean_s_anom_4(kmod,kisf) = mean_s_anom_4(kmod,kisf) + S_rcp(iCMIP,jCMIP) - S_his(iCMIP,jCMIP)
             NN = NN + 1
           endif
        endif
      enddo
      enddo
      !-
      if ( NN .ne. 0 ) then
        mean_t_anom_4(kmod,kisf) = mean_t_anom_4(kmod,kisf) / NN
        mean_s_anom_4(kmod,kisf) = mean_s_anom_4(kmod,kisf) / NN
      else
        mean_t_anom_4(kmod,kisf) = NF90_FILL_FLOAT
        mean_s_anom_4(kmod,kisf) = NF90_FILL_FLOAT
      endif
      write(*,*) '  4> ', mean_t_anom_4(kmod,kisf), NN

      ! within 5deg from the ice shelf front but only on the continental shelf (sea floor < 1500m)
      NN = 0.e0
      do iCMIP=1,mx
      do jCMIP=1,j50S
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
        if (       tmp_lon(iCMIP,jCMIP) .ge. front_inf  &
        &    .and. tmp_lon(iCMIP,jCMIP) .le. front_sup  &
        &    .and. nav_lat(iCMIP,jCMIP) .ge. front_min_lat(kisf) - lat_box_2  &
        &    .and. nav_lat(iCMIP,jCMIP) .le. front_max_lat(kisf) + lat_box_2  ) then
           if (         T_his(iCMIP,jCMIP) .ge.  -5.0 .and. T_his(iCMIP,jCMIP) .lt.  33.0   &
           &    .and.   S_his(iCMIP,jCMIP) .ge.   0.0 .and. S_his(iCMIP,jCMIP) .lt.  60.0   ) then
             mean_t_anom_5(kmod,kisf) = mean_t_anom_5(kmod,kisf) + T_rcp(iCMIP,jCMIP) - T_his(iCMIP,jCMIP)
             mean_s_anom_5(kmod,kisf) = mean_s_anom_5(kmod,kisf) + S_rcp(iCMIP,jCMIP) - S_his(iCMIP,jCMIP)
             NN = NN + 1
           endif
        endif
      enddo
      enddo
      !-
      if ( NN .ne. 0 ) then
        mean_t_anom_5(kmod,kisf) = mean_t_anom_5(kmod,kisf) / NN
        mean_s_anom_5(kmod,kisf) = mean_s_anom_5(kmod,kisf) / NN
      else
        mean_t_anom_5(kmod,kisf) = NF90_FILL_FLOAT
        mean_s_anom_5(kmod,kisf) = NF90_FILL_FLOAT
      endif
      write(*,*) '  5> ', mean_t_anom_5(kmod,kisf), NN
    
     ELSE
   
      mean_t_anom_4(kmod,kisf) = NF90_FILL_FLOAT
      mean_t_anom_5(kmod,kisf) = NF90_FILL_FLOAT
      mean_s_anom_4(kmod,kisf) = NF90_FILL_FLOAT
      mean_s_anom_5(kmod,kisf) = NF90_FILL_FLOAT

     ENDIF ! IF ( melt_isf(kisf) .NE. 0.0 )
 
    ENDDO ! kisf

ENDDO ! kmod

!---------------------------------------
! Writing new netcdf file :                                   

write(*,*) 'Creating ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create output file')

status = NF90_DEF_DIM(fidM,"Nisf",Nisf,dimID_Nisf); call erreur(status,.TRUE.,"def_dimID_Nisf")
status = NF90_DEF_DIM(fidM,"Nmod",Nmod,dimID_Nmod); call erreur(status,.TRUE.,"def_dimID_Nmod")
status = NF90_DEF_DIM(fidM,"StrLen",strlen,dimID_strlen); call erreur(status,.TRUE.,"def_dimID_strlen")

status = NF90_DEF_VAR(fidM,"model_name",NF90_CHAR,(/dimID_strlen,dimID_Nmod/),modnam_ID); call erreur(status,.TRUE.,"def_var_modnam_ID")
status = NF90_DEF_VAR(fidM,"t_anom_4",NF90_FLOAT,(/dimID_Nmod,dimID_Nisf/),mean_t_anom_4_ID); call erreur(status,.TRUE.,"def_var_mean_t_anom_4_ID")
status = NF90_DEF_VAR(fidM,"s_anom_4",NF90_FLOAT,(/dimID_Nmod,dimID_Nisf/),mean_s_anom_4_ID); call erreur(status,.TRUE.,"def_var_mean_s_anom_4_ID")
status = NF90_DEF_VAR(fidM,"t_anom_5",NF90_FLOAT,(/dimID_Nmod,dimID_Nisf/),mean_t_anom_5_ID); call erreur(status,.TRUE.,"def_var_mean_t_anom_5_ID")
status = NF90_DEF_VAR(fidM,"s_anom_5",NF90_FLOAT,(/dimID_Nmod,dimID_Nisf/),mean_s_anom_5_ID); call erreur(status,.TRUE.,"def_var_mean_s_anom_5_ID")

status = NF90_PUT_ATT(fidM,modnam_ID,"standard_name","model_name"); call erreur(status,.TRUE.,"put_att_modnam_ID")
status = NF90_PUT_ATT(fidM,modnam_ID,"long_name","CMIP5 model name"); call erreur(status,.TRUE.,"put_att_modnam_ID")

status = NF90_PUT_ATT(fidM,mean_t_anom_4_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_t_anom_4_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_t_anom_4_ID,"long_name","Mean seafloor temperature anomaly (futur-present) within box_1"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_t_anom_4_ID,"standard_name","sea_water_temperature"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_s_anom_4_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_s_anom_4_ID,"units","1"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_s_anom_4_ID,"long_name","Mean seafloor salinity anomaly (futur-present) within box_2"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_s_anom_4_ID,"standard_name","sea_water_salinity"); call erreur(status,.TRUE.,"put_att_t_anom_ID")

status = NF90_PUT_ATT(fidM,mean_t_anom_5_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_t_anom_5_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_t_anom_5_ID,"long_name","Mean sea water temperature anomaly (futur-present) withing 5deg over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_t_anom_5_ID,"standard_name","sea_water_temperature"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_s_anom_5_ID,"_FillValue",NF90_FILL_FLOAT); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_s_anom_5_ID,"units","1"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_s_anom_5_ID,"long_name","Mean sea water salinity anomaly (futur-present) withing 5deg over the cont. shelf"); call erreur(status,.TRUE.,"put_att_t_anom_ID")
status = NF90_PUT_ATT(fidM,mean_s_anom_5_ID,"standard_name","sea_water_salinity"); call erreur(status,.TRUE.,"put_att_t_anom_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","calculate_TS_profiles_CMIP5_anom_bottom.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"origin","CMIP5 data base (historical & rcp85)"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"end_definition") 

status = NF90_PUT_VAR(fidM,modnam_ID,modnam);                   call erreur(status,.TRUE.,"var_modnam_ID")
status = NF90_PUT_VAR(fidM,mean_t_anom_4_ID,mean_t_anom_4); call erreur(status,.TRUE.,"var_mean_t_anom_4_ID")
status = NF90_PUT_VAR(fidM,mean_s_anom_4_ID,mean_s_anom_4); call erreur(status,.TRUE.,"var_mean_s_anom_4_ID")
status = NF90_PUT_VAR(fidM,mean_t_anom_5_ID,mean_t_anom_5); call erreur(status,.TRUE.,"var_mean_t_anom_5_ID")
status = NF90_PUT_VAR(fidM,mean_s_anom_5_ID,mean_s_anom_5); call erreur(status,.TRUE.,"var_mean_s_anom_5_ID")

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
