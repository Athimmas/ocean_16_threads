!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module mix_submeso 

!BOP
! !MODULE: mix_submeso 

! !DESCRIPTION:
!  This module contains routines for computing submesoscale mixing
!  using the Fox-Kemper, Ferrari, and Hallberg parameterization
!  for restratification by mixed layer eddies.

! !REVISION HISTORY:
!  SVN:$Id: mix_submeso.F90

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks
   use domain
   use constants
   use broadcast
   use grid
   use io
   use vertical_mix
   use vmix_kpp
   use time_management
   use tavg
   use exit_mod
   use registry
   use communicate
   use hmix_gm_submeso_share
   use omp_lib
   
#ifdef CCSMCOUPLED
   use shr_sys_mod
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_submeso,   &
             submeso_sf,     &
	     submeso_flux 

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  variables to save from one call to next
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      ktp = 1, kbt = 2      ! refer to the top and bottom halves of a
                               !  grid cell, respectively

   !dir$ attributes offload:mic :: SF_SUBM_X
   !dir$ attributes offload:mic :: SF_SUBM_Y
   real (r8), dimension(:,:,:,:,:,:), allocatable, public :: &
      SF_SUBM_X,  &       ! components of the submesoscale 
      SF_SUBM_Y           !  streamfunction

   !dir$ attributes offload:mic :: FZTOP_SUBM
   real (r8), dimension(:,:,:,:), allocatable, public :: &
         FZTOP_SUBM 

!-----------------------------------------------------------------------
!
!  tavg ids for tavg diagnostics related to submesoscale mixing.
!  Zonal and meridional refer here to logical space only.
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_USUBM,        &   ! zonal      submeso velocity
      tavg_VSUBM,        &   ! meridional submeso velocity 
      tavg_WSUBM,        &   ! vertical   submeso velocity
      tavg_ADVT_SUBM,    &   ! vertically-integrated T submeso 
                             !  advection tendency
      tavg_ADVS_SUBM,    &   ! vertically-integrated S submeso 
                             !  advection tendency
      tavg_VNT_SUBM,     &   ! heat flux tendency in grid-y direction
                             !  due to submeso velocity
      tavg_VNS_SUBM,     &   ! salt flux tendency in grid-y direction
                             !  due to submeso velocity
      tavg_HLS_SUBM          ! horizontal length scale used in horizontal
                             !  buoyancy gradient scaling in submeso

   !dir$ attributes offload:mic :: TIME_SCALE
   real (r8), dimension(:,:,:), allocatable, public :: &
      TIME_SCALE             ! time scale used in horizontal length scale
                             !  calculation
   
   !dir$ attributes offload:mic :: max_hor_grid_scale
   real (r8), public :: &
      max_hor_grid_scale     ! maximum horizontal grid scale allowed

!-----------------------------------------------------------------------
!
!  namelist variables
!
!-----------------------------------------------------------------------

   !dir$ attributes offload:mic :: luse_const_horiz_len_scale
   logical (log_kind) ,public :: &
      luse_const_horiz_len_scale     ! if .true., then use a constant
                                     !  horizontal length scale given by
                                     !  hor_length_scale, otherwise the
                                     !  horizontal length scale varies both
                                     !  in space and time

   real (r8) :: &
      time_scale_constant            ! 1 day <= time scale constant <= 1 week

   !dir$ attributes offload:mic :: hor_length_scale
   !dir$ attributes offload:mic :: efficiency_factor
   real (r8),public :: &
      efficiency_factor,   &         ! 0.06 <= efficiency factor <= 0.08
      hor_length_scale               ! constant horizontal length scale used
                                     !  if luse_const_horiz_len_scale is true.
                                     !  if luse_const_horiz_len_scale is false,
                                     !  then hor_length_scale is used as the 
                                     !  lower limit.

!-----------------------------------------------------------------------
!
!  misc module variables
!
!-----------------------------------------------------------------------
   real (r8) :: &
      sqrt_grav              ! sqrt(grav)
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_submeso
! !INTERFACE:

   subroutine init_submeso

! !DESCRIPTION:
!  Initializes various submesoscale mixing options and allocates necessary
!  space.  Also computes some time-independent arrays.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables and namelist
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error,         &   ! error flag for namelist
      iblock                 ! block index

   namelist /mix_submeso_nml/ efficiency_factor,             &
                              time_scale_constant,           &
                              luse_const_horiz_len_scale,    &
                              hor_length_scale

!-----------------------------------------------------------------------
!
!  register init_submeso
!
!-----------------------------------------------------------------------

   call register_string ('init_submeso')

!-----------------------------------------------------------------------
!
!  default namelist values 
!
!-----------------------------------------------------------------------

   efficiency_factor             = 0.07_r8
   time_scale_constant           = 3.456e5_r8       ! 4 days, in seconds
   luse_const_horiz_len_scale    = .false.
   hor_length_scale              = 5.0e5_r8         ! 5 km

   max_hor_grid_scale            = 111.0e5_r8       ! about 1 degree

   if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old', iostat=nml_error)
     if (nml_error /= 0) then
       nml_error = -1
     else
       nml_error =  1
     endif
     do while (nml_error > 0)
       read(nml_in, nml=mix_submeso_nml, iostat=nml_error)
     enddo
     if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort, 'ERROR reading mix_submeso_nml')
   endif

   if (my_task == master_task) then

     write(stdout,*) ' '
     write(stdout,*) ' Document Namelist Parameters:'
     write(stdout,*) ' ============================ '
     write(stdout,*) ' '
     write(stdout,  mix_submeso_nml)
     write(stdout,*) ' '
     write(stdout,*) ' Submesoscale mixing options:'
     write(stdout,'(a21,1pe13.6)') ' efficiency factor = ', efficiency_factor
     write(stdout,'(a23,1pe13.6)') ' time scale constant = ', time_scale_constant
     if ( luse_const_horiz_len_scale ) then
       write(stdout,'(a45,1pe13.6)')  &
         ' using a constant horizontal length scale of ', hor_length_scale
     else
       write(stdout,'(a54)') ' horizontal length scale varies both in space and time'
     endif

     call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
   endif

   call broadcast_scalar (efficiency_factor,          master_task)
   call broadcast_scalar (time_scale_constant,        master_task)
   call broadcast_scalar (luse_const_horiz_len_scale, master_task)
   call broadcast_scalar (hor_length_scale,           master_task)

!-----------------------------------------------------------------------
!
!  define scalar constants
!
!-----------------------------------------------------------------------

   sqrt_grav  = sqrt(grav)

!-----------------------------------------------------------------------
!
!  allocate necessary arrays
!
!-----------------------------------------------------------------------

   allocate (SF_SUBM_X(nx_block,ny_block,2,2,km,nblocks_clinic),  &
             SF_SUBM_Y(nx_block,ny_block,2,2,km,nblocks_clinic))

   allocate (TIME_SCALE(nx_block,ny_block,nblocks_clinic))

   allocate (FZTOP_SUBM(nx_block,ny_block,nt,nblocks_clinic))

   SF_SUBM_X  = c0
   SF_SUBM_Y  = c0
   TIME_SCALE = c0
   FZTOP_SUBM = c0

!-----------------------------------------------------------------------
!
!  initialize various time-independent arrays
!
!-----------------------------------------------------------------------

   do iblock = 1,nblocks_clinic
     TIME_SCALE(:,:,iblock) = c1 / sqrt( FCORT(:,:,iblock)**2  &
                            + c1 / (time_scale_constant**2) )
   enddo

!-----------------------------------------------------------------------
!
!  define tavg fields related to the submesoscale parameterization 
!
!-----------------------------------------------------------------------

   call define_tavg_field (tavg_USUBM, 'USUBM', 3,                  &
    long_name='Submeso velocity in grid-x direction (diagnostic)',  &
                 units='cm/s', grid_loc='3211',                     &
                 coordinates='ULONG TLAT z_t time')

   call define_tavg_field (tavg_VSUBM, 'VSUBM', 3,                  &
    long_name='Submeso velocity in grid-y direction (diagnostic)',  &
                 units='cm/s', grid_loc='3121',                     &
                 coordinates='TLONG ULAT z_t time')

   call define_tavg_field (tavg_WSUBM, 'WSUBM', 3,                &
    long_name='Vertical submeso velocity (diagnostic)',           &
                 units='cm/s', grid_loc='3112',                   &
                 coordinates='TLONG TLAT z_w time')

   call define_tavg_field (tavg_HLS_SUBM, 'HLS_SUBM', 2,          &
    long_name='Horizontal length scale used in submeso',          &
                 units='cm', grid_loc='2110',                     &
                 coordinates='TLONG TLAT time')

   call define_tavg_field (tavg_ADVT_SUBM, 'ADVT_SUBM', 2,                       &
    long_name='Vertically-Integrated T submeso Advection Tendency (diagnostic)', &
                 units='cm degC/s', grid_loc='2110',                             &
                 coordinates='TLONG TLAT time')

   call define_tavg_field (tavg_ADVS_SUBM, 'ADVS_SUBM', 2,                       &
    long_name='Vertically-Integrated S submeso Advection Tendency (diagnostic)', &
                 scale_factor=1000.0_rtavg,                                      &
                 units='cm gram/kilogram/s', grid_loc='2110',                    &
                 coordinates='TLONG TLAT time')

   call define_tavg_field (tavg_VNT_SUBM, 'VNT_SUBM', 3,                          &
    long_name='Heat Flux Tendency in grid-y Dir due to submeso Vel (diagnostic)', &
                 units='degC/s', grid_loc='3121',                                 &
                 coordinates='TLONG ULAT z_t time')

   call define_tavg_field (tavg_VNS_SUBM, 'VNS_SUBM', 3,                          &
    long_name='Salt Flux Tendency in grid-y Dir due to submeso Vel (diagnostic)', &
                 scale_factor=1000.0_rtavg,                                       &
                 units='gram/kilogram/s', grid_loc='3121',                        &
                 coordinates='TLONG ULAT z_t time')

!-----------------------------------------------------------------------
!EOC

   end subroutine init_submeso

!***********************************************************************
!BOP
! !IROUTINE: submeso_sf 
! !INTERFACE:

   !dir$ attributes offload:mic :: submeso_sf
   subroutine submeso_sf ( TMIX, this_block )

! !DESCRIPTION:
!  The Fox-Kemper, Ferrari, and Hallberg [2008] submesoscale parameterization
!  for restratification by mixed layer eddies.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (block), intent(in) :: &
      this_block            ! block info for this sub block

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TMIX                  ! tracers at all vertical levels
                            !   at mixing time level
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j, k, kk,       &  ! dummy loop counters
      kp1,               &
      bid,temp              ! local block address for this sub block

   real (r8), dimension(nx_block,ny_block) :: &
      ML_DEPTH,          &  ! mixed layer depth
      HLS,               &  ! horizontal length scale
      WORK1, WORK2,      &  ! work arrays
      WORK3,             &
      USMT, VSMT,        &  ! arrays for submeso velocity computation
      USMB, VSMB,        &
      U_SUBM, V_SUBM,    &
      WBOT_SUBM,         &
      WTOP_SUBM

   real (r8), dimension(2) :: &
      reference_depth

   logical (log_kind), dimension(nx_block,ny_block) :: &
      CONTINUE_INTEGRAL     ! flag

   real (r8), dimension(nx_block,ny_block,2) :: &
      BX_VERT_AVG,       &  ! horizontal buoyancy differences vertically 
      BY_VERT_AVG           !  averaged within the mixed layer   

   real (r8) :: &
      zw_top, factor

   real (r8) :: &
      start_time, end_time

      
!-----------------------------------------------------------------------
!
!  initialize various quantities
!
!-----------------------------------------------------------------------
         
   bid = this_block%local_id
   

   WORK1     = c0
   WORK2     = c0
   WORK3     = c0
   HLS       = c0
   USMT      = c0
   USMB      = c0
   VSMT      = c0
   VSMB      = c0
   U_SUBM    = c0
   V_SUBM    = c0
   WTOP_SUBM = c0
   WBOT_SUBM = c0

   BX_VERT_AVG = c0
   BY_VERT_AVG = c0


     !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(k,kk,temp,j,i)num_threads(60)collapse(4)
     do k=1,km
        do kk=1,2
           do temp=1,2
              do j=1,ny_block
                 do i=1,nx_block

                    SF_SUBM_X(i,j,temp,kk,k,bid) = c0
                    SF_SUBM_Y(i,j,temp,kk,k,bid) = c0

                 enddo
              enddo
            enddo
          enddo
      enddo

   !print *,"base time is",end_time - start_time

   ML_DEPTH = zw(1)
   if ( vmix_itype == vmix_type_kpp )  &
     ML_DEPTH(:,:) = HMXL(:,:,bid) 
   

   do j=1,ny_block
      do i=1,nx_block
          CONTINUE_INTEGRAL(i,j) = .true.
          if( KMT(i,j,bid) == 0 ) then
           CONTINUE_INTEGRAL(i,j) = .false.
          endif
      enddo
   enddo
!-----------------------------------------------------------------------
!
!  compute vertical averages of horizontal buoyancy differences 
!  within the mixed layer 
!
!-----------------------------------------------------------------------

   !start_time = omp_get_wtime()
   do k=1,km
   
     zw_top = c0
     if ( k > 1 )  zw_top = zw(k-1)

    !!$OMP PARALLEL DO SHARED(CONTINUE_INTEGRAL,BX_VERT_AVG,RX,RY,ML_DEPTH)PRIVATE(i,WORK3)num_threads(60)SCHEDULE(dynamic,16)
    do j=1,ny_block
        do i=1,nx_block

            WORK3(i,j)=c0
            if( CONTINUE_INTEGRAL(i,j)  .and.  ML_DEPTH(i,j) > zw(k) )then
               WORK3(i,j) = dz(k)
           endif 

            if( CONTINUE_INTEGRAL(i,j)  .and.  ML_DEPTH(i,j) <= zw(k)  &
             .and.  ML_DEPTH(i,j) > zw_top ) then
                    WORK3(i,j) = ML_DEPTH(i,j) - zw_top
            endif

            if ( CONTINUE_INTEGRAL(i,j) ) then
                 BX_VERT_AVG(i,j,1) = BX_VERT_AVG(i,j,1)        &
                                      + RX(i,j,1,k,bid) * WORK3(i,j)
                 BX_VERT_AVG(i,j,2) = BX_VERT_AVG(i,j,2)        &
                                      + RX(i,j,2,k,bid) * WORK3(i,j)
                 BY_VERT_AVG(i,j,1) = BY_VERT_AVG(i,j,1)        &
                                      + RY(i,j,1,k,bid) * WORK3(i,j)
                 BY_VERT_AVG(i,j,2) = BY_VERT_AVG(i,j,2)        &
                           + RY(i,j,2,k,bid) * WORK3(i,j)
            endif

            if ( CONTINUE_INTEGRAL(i,j) .and.  ML_DEPTH(i,j) <= zw(k)  &
                  .and.  ML_DEPTH(i,j) > zw_top ) then
             CONTINUE_INTEGRAL(i,j) = .false.
            endif

      enddo
   enddo

  enddo

#ifdef CCSMCOUPLED
   if ( any(CONTINUE_INTEGRAL) ) then
     print *,'Incorrect mixed layer depth in submeso subroutine (I)'
   endif
#endif

    do j=1,ny_block
        do i=1,nx_block

           if ( KMT(i,j,bid) > 0 ) then
           BX_VERT_AVG(i,j,1) = - grav * BX_VERT_AVG(i,j,1) / ML_DEPTH(i,j)
           BX_VERT_AVG(i,j,2) = - grav * BX_VERT_AVG(i,j,2) / ML_DEPTH(i,j)
           BY_VERT_AVG(i,j,1) = - grav * BY_VERT_AVG(i,j,1) / ML_DEPTH(i,j)
           BY_VERT_AVG(i,j,2) = - grav * BY_VERT_AVG(i,j,2) / ML_DEPTH(i,j)
           endif
 
        enddo
     enddo
    
    !end_time = omp_get_wtime()
    ! print *,"Time at part1 is ",end_time - start_time

!-----------------------------------------------------------------------
!
!  compute horizontal length scale if necessary
!
!-----------------------------------------------------------------------
   !start_time = omp_get_wtime()


   if ( luse_const_horiz_len_scale ) then

    do j=1,ny_block
        do i=1,nx_block

           if ( KMT(i,j,bid) > 0 ) then
           HLS(i,j) = hor_length_scale
           endif

        enddo
     enddo
     

   else  

    do j=1,ny_block
        do i=1,nx_block
 
           WORK1(i,j)=c0
 
           if( KMT(i,j,bid) > 0 )then
                 WORK1(i,j) = sqrt( p5 * (                          &
                              ( BX_VERT_AVG(i,j,1)**2 + BX_VERT_AVG(i,j,2)**2 )  &
                                / DXT(i,j,bid)**2                                &
                              + ( BY_VERT_AVG(i,j,1)**2 + BY_VERT_AVG(i,j,2)**2 )  &
                                / DYT(i,j,bid)**2 ) )
                  WORK1(i,j) = WORK1(i,j) * ML_DEPTH(i,j) * (TIME_SCALE(i,j,bid)**2)
            endif


            CONTINUE_INTEGRAL(i,j) = .true.
            if( KMT(i,j,bid) == 0 ) then
                        CONTINUE_INTEGRAL(i,j) = .false.
            endif
 
            WORK2(i,j) = c0


        enddo
     enddo

          
     do k=2,km
        !!$OMP PARALLEL DO SHARED(WORK3,CONTINUE_INTEGRAL,WORK2,k,bid,dzw,zt,dzwr,RZ_SAVE,ML_DEPTH)PRIVATE(i,j)DEFAULT(NONE)num_threads(60)
        do j=1,ny_block
           do i=1,nx_block

              WORK3(i,j) = c0
              if ( CONTINUE_INTEGRAL(i,j)  .and.  ML_DEPTH(i,j) > zt(k) ) then
                   WORK3(i,j) = dzw(k-1) 
              endif

             if ( CONTINUE_INTEGRAL(i,j)  .and.  ML_DEPTH(i,j) <= zt(k)  &
                  .and.  ML_DEPTH(i,j) >= zt(k-1) ) then
                  WORK3(i,j) = ( (ML_DEPTH(i,j) - zt(k-1))**2 ) * dzwr(k-1)
             endif

             if ( CONTINUE_INTEGRAL(i,j) ) then
             WORK2(i,j) = WORK2(i,j) + sqrt(-RZ_SAVE(i,j,k,bid) * WORK3(i,j))
             endif

             if ( CONTINUE_INTEGRAL(i,j) .and.  ML_DEPTH(i,j) <= zt(k)  &
                     .and.  ML_DEPTH(i,j) >= zt(k-1) )then
                CONTINUE_INTEGRAL(i,j) = .false.
             endif

           enddo
        enddo

    enddo

#ifdef CCSMCOUPLED
     if ( any(CONTINUE_INTEGRAL) ) then
       print *,'Incorrect mixed layer depth in submeso subroutine (II)'
     endif
#endif

       do j=1,ny_block
           do i=1,nx_block

             if ( KMT(i,j,bid) > 0 ) then

                  WORK2(i,j) = sqrt_grav * WORK2(i,j) * TIME_SCALE(i,j,bid)

                  HLS(i,j) = max ( WORK1(i,j), WORK2(i,j), hor_length_scale )

             endif

           enddo
        enddo


   endif

   !end_time = omp_get_wtime()
   !print *,"Time at part2 is",end_time - start_time

!     if(master_task == my_task) then
!      open(unit=10,file="/home/aketh/ocn_correctness_data/changed.txt",status="unknown",position="append",action="write")
!       write(10,*),WORK1,WORK2,WORK3,HLS
!       close(10)
!   endif

  

!-----------------------------------------------------------------------
!
!  compute streamfunction due to submesoscale parameterization 
!
!-----------------------------------------------------------------------

   !start_time = omp_get_wtime() 

   do k=1,km

     reference_depth(ktp) = zt(k) - p25 * dz(k)
     reference_depth(kbt) = zt(k) + p25 * dz(k)

     do kk=ktp,kbt

       !!$OMP PARALLEL DO DEFAULT(NONE)PRIVATE(i,j)SHARED(reference_depth,ML_DEPTH,KMT,WORK3,WORK2,WORK1,TIME_SCALE,HLS,SF_SUBM_X,SF_SUBM_Y) &
       !!$OMP SHARED(BX_VERT_AVG,BY_VERT_AVG,DXT,DYT,max_hor_grid_scale,efficiency_factor,kk,k,bid)num_threads(60)  
       do j=1,ny_block
           do i=1,nx_block


               if ( reference_depth(kk) < ML_DEPTH(i,j)  .and.  &
                    KMT(i,j,bid) >= k ) then

                        WORK3(i,j) = ( c1 - ( c2 * reference_depth(kk) / ML_DEPTH(i,j) ) )**2
            
                        WORK2(i,j) = ( c1 - WORK3(i,j) )  &
                                    * ( c1 + ( 5.0_r8 / 21.0_r8 ) * WORK3(i,j) )

                        WORK1(i,j) = efficiency_factor * (ML_DEPTH(i,j)**2) * WORK2(i,j)  &
                                     * TIME_SCALE(i,j,bid) / HLS(i,j)

!     in the following negative sign is omitted to be consistent with
!     the GM implementation in hmix_gm subroutine. also, DXT and
!     DYT usage is approximate. 

                        SF_SUBM_X(i,j,1,kk,k,bid) = WORK1(i,j) * BX_VERT_AVG(i,j,1)  &
                                                          * min(DXT(i,j,bid),max_hor_grid_scale)
                        SF_SUBM_X(i,j,2,kk,k,bid) = WORK1(i,j) * BX_VERT_AVG(i,j,2)  &
                                                          * min(DXT(i,j,bid),max_hor_grid_scale)
                        SF_SUBM_Y(i,j,1,kk,k,bid) = WORK1(i,j) * BY_VERT_AVG(i,j,1)  &
                                                          * min(DYT(i,j,bid),max_hor_grid_scale)
                        SF_SUBM_Y(i,j,2,kk,k,bid) = WORK1(i,j) * BY_VERT_AVG(i,j,2)  &
                                                          * min(DYT(i,j,bid),max_hor_grid_scale)

               endif

           enddo
        enddo


     enddo

   enddo

   !end_time = omp_get_wtime()
   !print *,"Part 3 timings is",end_time - start_time

   USMT = c0
   VSMT = c0

   !start_time = omp_get_wtime()

   do k=1,km

!-----------------------------------------------------------------------
!
!  diagnostic computation of the submeso velocities
!
!-----------------------------------------------------------------------

     kp1 = k+1
     factor = c1
     if ( k == km ) then
       kp1 = k
       factor = c0
     endif

     !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60)
     do j=1,ny_block
       do i=1,nx_block

         if(j<=ny_block-1 .and. i<=nx_block-1) then

         WORK1(i,j) = (   SF_SUBM_X(i  ,j,1,kbt,k,  bid)    &
               + factor * SF_SUBM_X(i  ,j,1,ktp,kp1,bid)    &
                        + SF_SUBM_X(i+1,j,2,kbt,k,  bid)    &
               + factor * SF_SUBM_X(i+1,j,2,ktp,kp1,bid) )  &
                * p25 * HYX(i,j,bid)

         WORK2(i,j) = (   SF_SUBM_Y(i,j  ,1,kbt,k,  bid)    &
               + factor * SF_SUBM_Y(i,j  ,1,ktp,kp1,bid)    &
                        + SF_SUBM_Y(i,j+1,2,kbt,k,  bid)    &
               + factor * SF_SUBM_Y(i,j+1,2,ktp,kp1,bid) )  &
               * p25 * HXY(i,j,bid)

         endif


           USMB(i,j) = merge( WORK1(i,j), c0, k < KMT(i,j,bid) .and. k < KMTE(i,j,bid) )
           VSMB(i,j) = merge( WORK2(i,j), c0, k < KMT(i,j,bid) .and. k < KMTN(i,j,bid) )

            WORK1(i,j) = merge( USMT(i,j) - USMB(i,j), c0, k <= KMT(i,j,bid)  &
                                      .and. k <= KMTE(i,j,bid) )

            WORK2(i,j) = merge( VSMT(i,j) - VSMB(i,j), c0, k <= KMT(i,j,bid)  &
                                      .and. k <= KMTN(i,j,bid) )

            U_SUBM(i,j) = WORK1(i,j) * dzr(k) / HTE(i,j,bid)
            V_SUBM(i,j) = WORK2(i,j) * dzr(k) / HTN(i,j,bid)

       enddo
     enddo



     !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60)
     do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie

         if ( k < KMT(i,j,bid)  .and.  ( zw(k) < max( ML_DEPTH(i,j),  &
              ML_DEPTH(i+1,j), ML_DEPTH(i-1,j), ML_DEPTH(i,j+1),      &
              ML_DEPTH(i,j-1) ) ) ) then
           WBOT_SUBM(i,j) = WTOP_SUBM(i,j)             &
                           + TAREA_R(i,j,bid)          &
                       * ( WORK1(i,j) - WORK1(i-1,j)   &
                         + WORK2(i,j) - WORK2(i,j-1) )
         else
           WBOT_SUBM(i,j) = c0
         endif

       enddo
     enddo

!-----------------------------------------------------------------------
!
!  accumulate time average if necessary; checking is internal to routine
!
!-----------------------------------------------------------------------

     !if ( mix_pass /= 1 ) then

       !if ( k == 1 ) then
         !call accumulate_tavg_field (HLS, tavg_HLS_SUBM, bid, 1)  
       !endif

       !call accumulate_tavg_field (U_SUBM, tavg_USUBM, bid, k)
       !call accumulate_tavg_field (V_SUBM, tavg_VSUBM, bid, k)
       !call accumulate_tavg_field (WTOP_SUBM, tavg_WSUBM, bid, k) 

       !if (accumulate_tavg_now(tavg_ADVT_SUBM) ) then

         !WORK1 = p5 * HTE(:,:,bid) * U_SUBM * ( TMIX(:,:,k,1)  &
         !          + eoshift(TMIX(:,:,k,1), dim=1, shift=1) )
         !WORK2 = eoshift(WORK1, dim=1, shift=-1)
         !WORK3 = WORK1 - WORK2

         !WORK1 = p5 * HTN(:,:,bid) * V_SUBM * ( TMIX(:,:,k,1)  &
         !          + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )
         !WORK2 = eoshift(WORK1, dim=2, shift=-1)
         !WORK3 = WORK3 + WORK1 - WORK2

         !WORK1 = c0
         !do j=this_block%jb,this_block%je
           !do i=this_block%ib,this_block%ie
             !if ( k <= KMT(i,j,bid) ) then
               !WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
             !endif
           !enddo
         !enddo

         !call accumulate_tavg_field (WORK1, tavg_ADVT_SUBM, bid, k)

       !endif

       !if (accumulate_tavg_now(tavg_ADVS_SUBM) ) then

         !WORK1 = p5 * HTE(:,:,bid) * U_SUBM * ( TMIX(:,:,k,2)  &
         !          + eoshift(TMIX(:,:,k,2), dim=1, shift=1) )
         !WORK2 = eoshift(WORK1, dim=1, shift=-1)
         !WORK3 = WORK1 - WORK2

         !WORK1 = p5 * HTN(:,:,bid) * V_SUBM * ( TMIX(:,:,k,2)  &
         !          + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )
         !WORK2 = eoshift(WORK1, dim=2, shift=-1)
         !WORK3 = WORK3 + WORK1 - WORK2

         !WORK1 = c0
         !do j=this_block%jb,this_block%je
           !do i=this_block%ib,this_block%ie
             !if ( k <= KMT(i,j,bid) ) then
               !WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
             !endif
           !enddo
         !enddo

         !call accumulate_tavg_field (WORK1, tavg_ADVS_SUBM, bid, k)

       !endif

       !if ( accumulate_tavg_now(tavg_VNT_SUBM)  .or.  &
       !     accumulate_tavg_now(tavg_VNS_SUBM) ) then

         !WORK1 = p5 * V_SUBM * HTN(:,:,bid) * TAREA_R(:,:,bid)

         !if (accumulate_tavg_now(tavg_VNT_SUBM) ) then

          ! WORK2 = WORK1 * (    TMIX(:,:,k,1)  &
          !            + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )

       !   call accumulate_tavg_field (WORK2, tavg_VNT_SUBM, bid, k)

         !endif

         !if (accumulate_tavg_now(tavg_VNS_SUBM) ) then

           !WORK2 = WORK1 * (    TMIX(:,:,k,2)  &
           !           + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )

       !    call accumulate_tavg_field (WORK2, tavg_VNS_SUBM, bid, k)

         !endif

       !endif

     !endif ! mix_pass ne 1

!-----------------------------------------------------------------------
!
!  update remaining bottom-face fields to top-face fields for next pass
!
!-----------------------------------------------------------------------

     USMT = USMB
     VSMT = VSMB

     WTOP_SUBM = WBOT_SUBM

   enddo
  
   !end_time = omp_get_wtime()
   !print *,"Part 4 timings is",end_time - start_time


!-----------------------------------------------------------------------
!EOC

   end subroutine submeso_sf 

!***********************************************************************

! !IROUTINE: submeso_flux
! !INTERFACE:

   !dir$ attributes offload:mic :: submeso_flux
   subroutine submeso_flux (k, GTK, TMIX, tavg_HDIFE_TRACER, &
                     tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

! !DESCRIPTION:
!  Computes the fluxes due to submesoscale mixing of tracers
!
!INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! depth level index
   
   type (block), intent(in) :: &
      this_block            ! block info for this sub block
   
   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level   

   integer (int_kind), dimension(nt), intent(in) :: &
      tavg_HDIFE_TRACER, &! tavg id for east face diffusive flux of tracer
      tavg_HDIFN_TRACER, &! tavg id for north face diffusive flux of tracer
      tavg_HDIFB_TRACER   ! tavg id for bottom face diffusive flux of tracer
   
! !OUTPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,nt), intent(out) :: &
         GTK     ! submesoscale tracer flux at level k
      
!-----------------------------------------------------------------------
!
!      local variables
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ieast  = 1, iwest  = 2,       &
         jnorth = 1, jsouth = 2
      integer (int_kind)  :: &
         i,j,n,                &!dummy loop counters
         bid,                  &! local block address for this sub block
         kp1,kk
      
      real (r8) :: &
         fz, factor 
 
      real (r8), dimension(nx_block,ny_block) :: &
         CX, CY,                  &
         WORK1, WORK2,            &! local work space
         KMASK                     ! ocean mask
         
      real (r8), dimension(nx_block,ny_block,nt)  :: &
         FX, FY                    ! fluxes across east, north faces

      logical :: reg_match_init_gm

      real (r8) :: WORK1prev,WORK2prev,KMASKprev,fzprev,GTKmy

     bid = this_block%local_id

     !if ( k == 1) FZTOP_SUBM(:,:,:,bid) = c0        ! zero flux B.C. at the surface
     
      CX = merge(HYX(:,:,bid)*p25, c0, (k <= KMT (:,:,bid))   &
                                 .and. (k <= KMTE(:,:,bid)))
      CY = merge(HXY(:,:,bid)*p25, c0, (k <= KMT (:,:,bid))   &
                                 .and. (k <= KMTN(:,:,bid)))

      !if(my_task == master_task .and. nsteps_total == 3 .and. k == 45 .and. i == 45 .and. j == 45 )then

            !print *,"original"
            !print *,"HYX(45,45,bid)",HYX(i,j,bid)
            !print *,"KMT(45,45,bid)",KMT(i,j,bid) 
            !print *,"CX(45,45,bid)",CX(i,j)

      !endif
      
      KMASK = merge(c1, c0, k < KMT(:,:,bid))
            
      kp1 = k + 1
      if ( k == km )  kp1 = k

      if ( k < km ) then
        factor    = c1
      else
        factor    = c0
      endif
     
 
      do n = 1,nt
          do j=1,ny_block
            do i=1,nx_block

              if(i <= nx_block-1 ) then


            !if(my_task == master_task .and. nsteps_total == 3 .and. k == 45 .and. i == 45 - 1 .and. j == 45 .and. n == 1)then

                   !print *,"original in flux"
                   !print *,"SF_SUBM_X(i+1,j,iwest,kbt,k,bid) contribution is",SF_SUBM_X(i+1,j,iwest,kbt,k,bid)
                   !print *,"SF_SUBM_X(i+1,j,iwest,ktp,k,bid) contribution is",SF_SUBM_X(i+1,j,iwest,ktp,k,bid)
                   !print *,"SF_SUBM_X(i  ,j,ieast,kbt,k,bid) contribution is",SF_SUBM_X(i  ,j,ieast,kbt,k,bid)
                   !print *,"SF_SUBM_X(i  ,j,ieast,ktp,k,bid)  contribution is",SF_SUBM_X(i  ,j,ieast,ktp,k,bid)
                   !print *,"CX(i,j) contribution is",CX(i,j)
                   !print *,"TZ(i,j,k,n,bid) is",TZ(i,j,k,n,bid)
                   !print *,"TZ(i,j,kp1,n,bid) ",TZ(i,j,kp1,n,bid) 
                   !print *,"TZ(i+1,j,k,n,bid)",TZ(i+1,j,k,n,bid)  
                   !print *,"TZ(i+1,j,kp1,n,bid)",TZ(i,j,kp1,n,bid)

            !endif



              FX(i,j,n) = CX(i,j)                          &
               * ( SF_SUBM_X(i  ,j,ieast,ktp,k,bid) * TZ(i,j,k,n,bid)                        &
                 + SF_SUBM_X(i  ,j,ieast,kbt,k,bid) * TZ(i,j,kp1,n,bid)                    &
                 + SF_SUBM_X(i+1,j,iwest,ktp,k,bid) * TZ(i+1,j,k,n,bid)                    &
                 + SF_SUBM_X(i+1,j,iwest,kbt,k,bid) * TZ(i+1,j,kp1,n,bid) )

              endif  

              if(j <= ny_block -1 )then

              FY(i,j,n) =  CY(i,j)                          &
               * ( SF_SUBM_Y(i,j  ,jnorth,ktp,k,bid) * TZ(i,j,k,n,bid)                        &
                 + SF_SUBM_Y(i,j  ,jnorth,kbt,k,bid) * TZ(i,j,kp1,n,bid)                    &
                 + SF_SUBM_Y(i,j+1,jsouth,ktp,k,bid) * TZ(i,j+1,k,n,bid)                    &
                 + SF_SUBM_Y(i,j+1,jsouth,kbt,k,bid) * TZ(i,j+1,kp1,n,bid) )

              endif

              WORK1(i,j) = c0 ! zero halo regions so accumulate_tavg_field calls do not trap
              WORK2(i,j) = c0
              GTK(i,j,n) = c0
              

            enddo
          enddo

       end do

      do n = 1,nt

!-----------------------------------------------------------------------
!
!     calculate vertical submesoscale fluxes thru horizontal faces of T-cell
!
!-----------------------------------------------------------------------
 


        do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              
               if ( k < km ) then

                  WORK1(i,j) = SF_SUBM_X(i  ,j  ,ieast ,kbt,k  ,bid)     &
                             * HYX(i  ,j  ,bid) * TX(i  ,j  ,k  ,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jnorth,kbt,k  ,bid)     &
                             * HXY(i  ,j  ,bid) * TY(i  ,j  ,k  ,n,bid)  &
                             + SF_SUBM_X(i  ,j  ,iwest ,kbt,k  ,bid)     &
                             * HYX(i-1,j  ,bid) * TX(i-1,j  ,k  ,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jsouth,kbt,k  ,bid)     &
                             * HXY(i  ,j-1,bid) * TY(i  ,j-1,k  ,n,bid)


                  WORK2(i,j) = factor                                    &
                           * ( SF_SUBM_X(i  ,j  ,ieast ,ktp,kp1,bid)     &
                             * HYX(i  ,j  ,bid) * TX(i  ,j  ,kp1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jnorth,ktp,kp1,bid)     &
                             * HXY(i  ,j  ,bid) * TY(i  ,j  ,kp1,n,bid)  &
                             + SF_SUBM_X(i  ,j  ,iwest ,ktp,kp1,bid)     &
                             * HYX(i-1,j  ,bid) * TX(i-1,j  ,kp1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jsouth,ktp,kp1,bid)     &
                             * HXY(i  ,j-1,bid) * TY(i  ,j-1,kp1,n,bid) ) 
                   
                  if(k==1) then
 
                  fzprev = 0
    
                  else    

    
                  WORK1prev = SF_SUBM_X(i  ,j  ,ieast ,kbt,k-1 ,bid)     &
                             * HYX(i  ,j  ,bid) * TX(i  ,j  ,k-1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jnorth,kbt,k-1,bid)     &
                             * HXY(i  ,j  ,bid) * TY(i  ,j  ,k-1,n,bid)  &
                             + SF_SUBM_X(i  ,j  ,iwest ,kbt,k-1,bid)     &
                             * HYX(i-1,j  ,bid) * TX(i-1,j  ,k-1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jsouth,kbt,k-1,bid)     &
                             * HXY(i  ,j-1,bid) * TY(i  ,j-1,k-1,n,bid)

                  WORK2prev = factor &
                           * ( SF_SUBM_X(i  ,j  ,ieast ,ktp,kp1-1,bid)     &
                             * HYX(i  ,j  ,bid) * TX(i  ,j  ,kp1-1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jnorth,ktp,kp1-1,bid)     &
                             * HXY(i  ,j  ,bid) * TY(i  ,j  ,kp1-1,n,bid)  &
                             + SF_SUBM_X(i  ,j  ,iwest ,ktp,kp1-1,bid)     &
                             * HYX(i-1,j  ,bid) * TX(i-1,j  ,kp1-1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jsouth,ktp,kp1-1,bid)     &
                             * HXY(i  ,j-1,bid) * TY(i  ,j-1,kp1-1,n,bid) )

                  KMASKprev = merge(c1, c0, k-1 < KMT(i,j,bid))


                  fzprev = -KMASKprev * p25 &
                            * (WORK1prev + WORK2prev)

                  endif 

                  !if(fzprev /= FZTOP_SUBM(i,j,n,bid)) then 
                  !   print *,"wrong value OH NO",k
                  !else
                  !   print *,"its okay Yeah",k
                  !endif    

                  fz = -KMASK(i,j) * p25    &
                      * (WORK1(i,j) + WORK2(i,j))


            !if(my_task == master_task .and. nsteps_total == 3 .and. k == 45 .and. i == 45 .and. j == 45 .and. n == 1)then

                   !print *,"original in flux"
                   !print *,"FX(i,j,n) contribution is",FX(i,j,n)
                   !print *,"FX(i-1,j,n) contribution is",FX(i-1,j,n)
                   !print *,"FY(i,j,n) contribution is",FY(i,j,n)
                   !print *,"FY(i,j-1,n) contribution is",FY(i,j-1,n)

            !endif


                  GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                               + FY(i,j,n) - FY(i,j-1,n)  &
                        + fzprev - fz )*dzr(k)*TAREA_R(i,j,bid)

                  !FZTOP_SUBM(i,j,n,bid) = fz

               else  

                  WORK1prev = SF_SUBM_X(i  ,j  ,ieast ,kbt,k-1 ,bid)     &
                             * HYX(i  ,j  ,bid) * TX(i  ,j  ,k-1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jnorth,kbt,k-1,bid)     &
                             * HXY(i  ,j  ,bid) * TY(i  ,j  ,k-1,n,bid)  &
                             + SF_SUBM_X(i  ,j  ,iwest ,kbt,k-1,bid)     &
                             * HYX(i-1,j  ,bid) * TX(i-1,j  ,k-1,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jsouth,kbt,k-1,bid)     &
                             * HXY(i  ,j-1,bid) * TY(i  ,j-1,k-1,n,bid)

                  WORK2prev = factor &
                           * ( SF_SUBM_X(i  ,j  ,ieast ,ktp,km,bid)     &
                             * HYX(i  ,j  ,bid) * TX(i  ,j  ,km,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jnorth,ktp,km,bid)     &
                             * HXY(i  ,j  ,bid) * TY(i  ,j  ,km,n,bid)  &
                             + SF_SUBM_X(i  ,j  ,iwest ,ktp,km,bid)     &
                             * HYX(i-1,j  ,bid) * TX(i-1,j  ,km,n,bid)  &
                             + SF_SUBM_Y(i  ,j  ,jsouth,ktp,km,bid)     &
                             * HXY(i  ,j-1,bid) * TY(i  ,j-1,km,n,bid) )

                  KMASKprev = merge(c1, c0, k-1 < KMT(i,j,bid))


                  fzprev = -KMASKprev * p25 &
                            * (WORK1prev + WORK2prev)

                  GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                               + FY(i,j,n) - FY(i,j-1,n)  &
                     + fzprev )*dzr(k)*TAREA_R(i,j,bid)

                   !FZTOP_SUBM(i,j,n,bid) = c0
              
               endif   

            enddo
          enddo


!-----------------------------------------------------------------------
!
!     accumulate time average if necessary
!
!-----------------------------------------------------------------------

      if ( mix_pass /= 1 ) then
          if (accumulate_tavg_now(tavg_HDIFE_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FX(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            !call accumulate_tavg_field(WORK1,tavg_HDIFE_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFN_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FY(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            !call accumulate_tavg_field(WORK1,tavg_HDIFN_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFB_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FZTOP_SUBM(i,j,n,bid)*dzr(k)*TAREA_R(i,j,bid)
            enddo
            enddo
            !call accumulate_tavg_field(WORK1,tavg_HDIFB_TRACER(n),bid,k)
          endif
      endif   ! mix_pass ne 1  

!-----------------------------------------------------------------------
!
!     end of tracer loop
!
!-----------------------------------------------------------------------

      enddo !ends the do n loop 

!-----------------------------------------------------------------------

   end subroutine submeso_flux 

!***********************************************************************

 end module mix_submeso

