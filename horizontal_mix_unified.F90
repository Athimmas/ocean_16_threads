!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module horizontal_mix_unified

!BOP
! !MODULE: horizontal_mix_unified
!
! !DESCRIPTION:
!  This module contains driver routines for managing the individual
!  horizontal tracer and momentum mixing modules for unified grids.
!
! !REVISION HISTORY:
!  SVN:$Id: horizontal_mix.F90 47361 2013-05-21 20:54:30Z $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_ConstantsMod

   use kinds_mod
   use blocks, only: nx_block, ny_block, block , nx_block_unified , &
   nx_block_unified
   use distribution 
   use domain_size
   use domain, only: nblocks_clinic, distrb_clinic
   use constants
   use communicate, only: my_task, master_task
   use time_management, only: km, nt, mix_pass
   use broadcast, only: broadcast_scalar
   use io_types, only: nml_in, nml_filename, stdout
   use hmix_del2, only: init_del2u, init_del2t, hdiffu_del2, hdifft_del2
   use hmix_del4, only: init_del4u, init_del4t, hdiffu_del4, hdifft_del4
   use hmix_gm, only: init_gm, hdifft_gm,diag_gm_bolus
   use hmix_aniso, only: init_aniso, hdiffu_aniso
   use topostress, only: ltopostress
   use horizontal_mix, only:tavg_HDIFE_TRACER,tavg_HDIFN_TRACER,tavg_HDIFB_TRACER
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now, &
      tavg_in_which_stream, ltavg_on
   use timers, only: timer_start, timer_stop, get_timer
   use exit_mod, only: sigAbort, exit_pop, flushm
   use mix_submeso, only: init_submeso, submeso_flux, submeso_sf
   use hmix_gm_submeso_share, only: init_meso_mixing, tracer_diffs_and_isopyc_slopes
   use vertical_mix, only: implicit_vertical_mix
   use state_mod, only: pressz
   use omp_lib

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_horizontal_mix_unified, &
             hdifft_unified,tracer_diffs_and_isopyc_slopes_unified

!EOP
!BOC

! !PUBLIC VARIABLES

   !dir$ attributes offload:mic :: RX_UNIFIED
   !dir$ attributes offload:mic :: RY_UNIFIED
   !dir$ attributes offload:mic :: TX_UNIFIED
   !dir$ attributes offload:mic :: TY_UNIFIED
   !dir$ attributes offload:mic :: TZ_UNIFIED
   !dir$ attributes offload:mic :: SLX_UNIFIED
   !dir$ attributes offload:mic :: SLY_UNIFIED
   !dir$ attributes offload:mic :: RZ_SAVE_UNIFIED
   !dir$ attributes offload:mic :: HXY_UNIFIED
   !dir$ attributes offload:mic :: HYX_UNIFIED 

   real (r8), dimension(:,:,:,:,:), allocatable, public :: &
      RX_UNIFIED,RY_UNIFIED,            &     ! Dx(rho), Dy(rho)
      TX_UNIFIED,TY_UNIFIED,TZ_UNIFIED               ! tracer differences in each direction
   real (r8), dimension(:,:,:,:,:,:), allocatable, public :: &
      SLX_UNIFIED, SLY_UNIFIED                ! slope of isopycnal sfcs in x,y-direction
   real (r8), dimension(:,:,:,:), allocatable, public :: &
      RZ_SAVE_UNIFIED                 ! Dz(rho)
   real (r8), dimension(:,:,:), allocatable, public :: &
      HXY_UNIFIED,              &     ! dx/dy for y-z plane
      HYX_UNIFIED                     ! dy/dx for x-z plane

!  !State related global variables


   !dir$ attributes offload:mic :: tmin
   !dir$ attributes offload:mic :: tmax
   !dir$ attributes offload:mic :: smin
   !dir$ attributes offload:mic :: smax
   !dir$ attributes offload:mic :: pressz_unified
   real (r8), dimension(km) :: &
      tmin, tmax,        &! valid temperature range for level k
      smin, smax,        &! valid salinity    range for level k
      pressz_unified      ! ref pressure (bars) at each level

   !*** these constants will be used to construct the numerator
   !*** factor unit change (kg/m^3 -> g/cm^3) into numerator terms

   real (r8), parameter ::                     &
      mwjfnp0s0t0 =   9.99843699e+2_r8 * p001, &
      mwjfnp0s0t1 =   7.35212840e+0_r8 * p001, &
      mwjfnp0s0t2 =  -5.45928211e-2_r8 * p001, &
      mwjfnp0s0t3 =   3.98476704e-4_r8 * p001, &
      mwjfnp0s1t0 =   2.96938239e+0_r8 * p001, &
      mwjfnp0s1t1 =  -7.23268813e-3_r8 * p001, &
      mwjfnp0s2t0 =   2.12382341e-3_r8 * p001, &
      mwjfnp1s0t0 =   1.04004591e-2_r8 * p001, &
      mwjfnp1s0t2 =   1.03970529e-7_r8 * p001, &
      mwjfnp1s1t0 =   5.18761880e-6_r8 * p001, &
      mwjfnp2s0t0 =  -3.24041825e-8_r8 * p001, &
      mwjfnp2s0t2 =  -1.23869360e-11_r8* p001

   !*** these constants will be used to construct the denominator

   real (kind=r8), parameter ::          &
      mwjfdp0s0t0 =   1.0e+0_r8,         &
      mwjfdp0s0t1 =   7.28606739e-3_r8,  &
      mwjfdp0s0t2 =  -4.60835542e-5_r8,  &
      mwjfdp0s0t3 =   3.68390573e-7_r8,  &
      mwjfdp0s0t4 =   1.80809186e-10_r8, &
      mwjfdp0s1t0 =   2.14691708e-3_r8,  &
      mwjfdp0s1t1 =  -9.27062484e-6_r8,  &
      mwjfdp0s1t3 =  -1.78343643e-10_r8, &
      mwjfdp0sqt0 =   4.76534122e-6_r8,  &
      mwjfdp0sqt2 =   1.63410736e-9_r8,  &
      mwjfdp1s0t0 =   5.30848875e-6_r8,  &
      mwjfdp2s0t3 =  -3.03175128e-16_r8, &
      mwjfdp3s0t1 =  -1.27934137e-17_r8

!GRID related variables

   !dir$ attributes offload : mic :: KMT_unified
   integer (POP_i4), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      KMT_unified

   !dir$ attributes offload : mic :: KMTN_unified
   !dir$ attributes offload : mic :: KMTE_unified
   integer (POP_i4), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      KMTN_unified,KMTE_unified

   !dir$ attributes offload:mic :: DZT_unified
   real (POP_r8), dimension(:,:,:,:), allocatable, public :: &
      DZT_unified

   !dir$ attributes offload:mic :: DXT_UNIFIED
   !dir$ attributes offload:mic :: DYT_UNIFIED

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic), public :: &
      DYT_UNIFIED, DXT_UNIFIED

      !*** dimension(1:km)

   !dir$ attributes offload:mic :: dz_unified
   !dir$ attributes offload:mic :: zw_unified
   !dir$ attributes offload:mic :: dzr_unified
   !dir$ attributes offload:mic :: zt_unified
   real (POP_r8), dimension(km), public :: &
      dz_unified                ,&! thickness of layer k
      dzr_unified               ,&! reciprocals of dz, c2dz
      zt_unified                ,&! vert dist from sfc to midpoint of layer
      zw_unified                  ! vert dist from sfc to bottom of layer

   !*** dimension(0:km)

   !dir$ attributes offload:mic :: dzw_unified
   !dir$ attributes offload:mic :: dzwr_unified
   real (POP_r8), dimension(0:km), public :: &
      dzw_unified, dzwr_unified          ! midpoint of k to midpoint of k+1
                                      !   and its reciprocal

   !dir$ attributes offload:mic :: HTN_UNIFIED
   !dir$ attributes offload:mic :: HTE_UNIFIED
   !dir$ attributes offload:mic :: TAREA_R_UNIFIED
   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic), public :: &
   HTN_UNIFIED, HTE_UNIFIED,TAREA_R_UNIFIED


!mix_submeso related variables

   !dir$ attributes offload:mic :: SF_SUBM_X_UNIFIED
   !dir$ attributes offload:mic :: SF_SUBM_Y_UNIFIED
   real (r8), dimension(:,:,:,:,:,:), allocatable, public :: &
      SF_SUBM_X_UNIFIED,  &       ! components of the submesoscale 
      SF_SUBM_Y_UNIFIED           !  streamfunction

   !dir$ attributes offload:mic :: HMXL_UNIFIED
   !dir$ attributes offload:mic :: BOLUS_SP_UNIFIED
   !dir$ attributes offload:mic :: KPP_HBLT_UNIFIED
   real (r8), dimension(:,:,:), allocatable, public :: &
      HMXL_UNIFIED,KPP_HBLT_UNIFIED,BOLUS_SP_UNIFIED

   !dir$ attributes offload:mic :: TIME_SCALE_UNIFIED
   real (r8), dimension(:,:,:), allocatable, public :: &
      TIME_SCALE_UNIFIED     ! time scale used in horizontal length scale
                             !  calculation

   !dir$ attributes offload:mic :: hor_length_scale_unified
   !dir$ attributes offload:mic :: efficiency_factor_unified
   real (r8),public :: &
      efficiency_factor_unified,   &         ! 0.06 <= efficiency factor <= 0.08
      hor_length_scale_unified               ! constant horizontal length scale used
                                             !  if luse_const_horiz_len_scale is true.
                                             !  if luse_const_horiz_len_scale is false,
                                             !  then hor_length_scale is used as the 
                                             !  lower limit.

   real (r8) :: &
      sqrt_grav_unified              ! sqrt(grav)

   integer (int_kind), parameter :: &
         ktp = 1, kbt = 2     ! refer to the top and bottom halves of a 
                              ! grid cell, respectively

  !dir$ attributes offload:mic :: max_hor_grid_scale_unified
   real (r8), public :: &
      max_hor_grid_scale_unified     ! maximum horizontal grid scale allowed

   !dir$ attributes offload:mic :: FZTOP_SUBM_UNIFIED
   real (r8), dimension(:,:,:,:), allocatable, public :: &
         FZTOP_SUBM_UNIFIED



!vertical mix related variables


!variables related to hmix_gm

      !dir$ attributes offload:mic :: WTOP_ISOP_UNIFIED
      !dir$ attributes offload:mic :: WBOT_ISOP_UNIFIED
      !dir$ attributes offload:mic :: HXYS_UNIFIED
      !dir$ attributes offload:mic :: HYXW_UNIFIED
      !dir$ attributes offload:mic :: RB_UNIFIED
      !dir$ attributes offload:mic :: RBR_UNIFIED
      !dir$ attributes offload:mic :: BTP_UNIFIED
      !dir$ attributes offload:mic :: BL_DEPTH_UNIFIED
      !dir$ attributes offload:mic :: UIT_UNIFIED
      !dir$ attributes offload:mic :: VIT_UNIFIED   
      real (r8), dimension(:,:,:), allocatable ,public :: &
         HYXW_UNIFIED, HXYS_UNIFIED, &        ! west and south-shifted values of above
         RB_UNIFIED,                 &        ! Rossby radius
         RBR_UNIFIED,                &        ! inverse of Rossby radius
         BTP_UNIFIED,                &        ! beta plane approximation
         BL_DEPTH_UNIFIED,           &        ! boundary layer depth
         UIT_UNIFIED, VIT_UNIFIED             ! work arrays for isopycnal mixing velocities


      real (r8), dimension(:,:,:), allocatable, public :: WTOP_ISOP_UNIFIED,WBOT_ISOP_UNIFIED

      !dir$ attributes offload:mic :: HOR_DIFF_UNIFIED
      real (r8), dimension(:,:,:,:,:), public ,allocatable :: &
         HOR_DIFF_UNIFIED   ! 3D horizontal diffusion coefficient
                            !  for top and bottom half of a grid cell


      !dir$ attributes offload:mic :: ah_unified
      !dir$ attributes offload:mic :: ah_bolus_unified
      !dir$ attributes offload:mic :: ah_bkg_bottom_unified
      !dir$ attributes offload:mic :: ah_bkg_srfbl_unified                       
      !dir$ attributes offload:mic :: slm_r_unified
      !dir$ attributes offload:mic :: slm_b_unified 
      real (r8), public ::      &
         ah_unified,            &       ! isopycnal diffusivity
         ah_bolus_unified,      &       ! thickness (GM bolus) diffusivity
         ah_bkg_bottom_unified, &       ! backgroud horizontal diffusivity at k = KMT
         ah_bkg_srfbl_unified,  &       ! backgroud horizontal diffusivity within the surface boundary layer
         slm_r_unified,         &       ! max. slope allowed for redi diffusion
         slm_b_unified                  ! max. slope allowed for bolus transport

      type tlt_info_unified
        real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
           DIABATIC_DEPTH,  &   ! depth of the diabatic region at the
                                !  surface, i.e. mean mixed or boundary layer
                                !  depth
           THICKNESS,       &   ! transition layer thickness
           INTERIOR_DEPTH       ! depth at which the interior, adiabatic
                                !  region starts, i.e.
                                !   = TLT%DIABATIC_DEPTH + TLT%THICKNESS
        integer (int_kind), &
              dimension(nx_block,ny_block,max_blocks_clinic) :: &
           K_LEVEL,  &          ! k level at or below which the interior,
                                !  adiabatic region starts
           ZTW                  ! designates if the interior region
                                !  starts below depth zt or zw.
                                !  ( = 1 for zt, = 2 for zw )
      end type tlt_info_unified

      !dir$ attributes offload:mic :: TLT_UNIFIED
      type (tlt_info_unified),public :: &
         TLT_UNIFIED            ! transition layer thickness related fields 


!variables realted to vmix_kpp

   !dir$ attributes offload : mic :: zgrid_unified
   real (r8), dimension(:), allocatable, public :: &
      zgrid_unified           ! depth at cell interfaces


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_horizontal_mix
! !INTERFACE:

 subroutine init_horizontal_mix_unified()

 use mix_submeso, only: TIME_SCALE,efficiency_factor,hor_length_scale, &
                        max_hor_grid_scale
 use grid, only: KMT, dz, partial_bottom_cells, DZT, dzr, dzwr,KMTE, KMT,&
                   KMTN,zt,dzw,zw,DXT,DYT,HTN,HTE,TAREA_R
 use vmix_kpp, only:zgrid
 use prognostic
 use vertical_mix
 use omp_lib


   print *,"Intializing unified grids"

   allocate (HXY_UNIFIED (nx_block,ny_block,nblocks_clinic),    &
             HYX_UNIFIED (nx_block,ny_block,nblocks_clinic))
   allocate (SLX_UNIFIED   (nx_block,ny_block,2,2,km,nblocks_clinic),  &
             SLY_UNIFIED   (nx_block,ny_block,2,2,km,nblocks_clinic))

   allocate (TX_UNIFIED(nx_block,ny_block,km,nt,nblocks_clinic),  &
             TY_UNIFIED(nx_block,ny_block,km,nt,nblocks_clinic),  &
             TZ_UNIFIED(nx_block,ny_block,km,nt,nblocks_clinic))

   allocate (RX_UNIFIED(nx_block,ny_block,2,km,nblocks_clinic),  &
             RY_UNIFIED(nx_block,ny_block,2,km,nblocks_clinic))

   allocate (RZ_SAVE_UNIFIED(nx_block,ny_block,km,nblocks_clinic))

   allocate (SF_SUBM_X_UNIFIED(nx_block,ny_block,2,2,km,nblocks_clinic))
   allocate (SF_SUBM_Y_UNIFIED(nx_block,ny_block,2,2,km,nblocks_clinic))

   allocate (HMXL_UNIFIED(nx_block,ny_block,nblocks_clinic))
 
   allocate (HOR_DIFF_UNIFIED(nx_block,ny_block,2,km,nblocks_clinic))


   allocate(WTOP_ISOP_UNIFIED(nx_block,ny_block,nblocks_clinic), &
              WBOT_ISOP_UNIFIED(nx_block,ny_block,nblocks_clinic), &
                    UIT_UNIFIED(nx_block,ny_block,nblocks_clinic), &
                    VIT_UNIFIED(nx_block,ny_block,nblocks_clinic))

   allocate (KPP_HBLT_UNIFIED(nx_block,ny_block,nblocks_clinic), &
            BOLUS_SP_UNIFIED(nx_block,ny_block,nblocks_clinic))

   HXY_UNIFIED      = c0
   HYX_UNIFIED      = c0
   SLX_UNIFIED      = c0
   SLY_UNIFIED      = c0
   TX_UNIFIED       = c0
   TY_UNIFIED       = c0
   TZ_UNIFIED       = c0
   RX_UNIFIED       = c0
   RY_UNIFIED       = c0
   RZ_SAVE_UNIFIED  = c0

!! initialization for state 

   tmin =  -2.0_r8  ! limited   on the low  end
   tmax = 999.0_r8  ! unlimited on the high end
   smin =   0.0_r8  ! limited   on the low  end
   smax = 0.999_r8  ! unlimited on the high end

!! initialization for grid 

   allocate (DZT_unified(nx_block,ny_block,0:km+1,max_blocks_clinic))

   pressz_unified = pressz
   KMT_unified = KMT
   KMTE_unified = KMTE
   KMTN_unified = KMTN
   DZT_unified = DZT

   dz_unified = dz
   dzr_unified = dzr
   zt_unified = zt
   zw_unified = zw
   dzw_unified = dzw
   dzwr_unified = dzwr
   

   DXT_UNIFIED = DXT
   DYT_UNIFIED = DYT
   HTE_UNIFIED = HTE
   HTN_UNIFIED = HTN
   TAREA_R_UNIFIED = TAREA_R


!! initialization for mix_submeso

   allocate (TIME_SCALE_UNIFIED(nx_block,ny_block,nblocks_clinic))

   allocate (FZTOP_SUBM_UNIFIED(nx_block,ny_block,nt,nblocks_clinic))

   TIME_SCALE_UNIFIED = TIME_SCALE

      efficiency_factor_unified = efficiency_factor
      hor_length_scale_unified = hor_length_scale
      sqrt_grav_unified = sqrt(grav)
      max_hor_grid_scale_unified = max_hor_grid_scale

!! Variables related to vmix_kpp

  allocate  (zgrid_unified(0:km+1)) 

  zgrid_unified = zgrid

!-----------------------------------------------------------------------
!EOC

 end subroutine init_horizontal_mix_unified

!***********************************************************************


! !IROUTINE: hdifft_unified
! !INTERFACE:

 !dir$ attributes offload: mic :: hdifft_unified
 subroutine hdifft_unified(k, HDTK, TMIX, UMIX, VMIX, this_block)

! !DESCRIPTION:
!  This routine returns tendencies for horizontal diffusion of
!  tracers.  It is a driver routine which simply branches to the
!  proper horizontal mix routine based on the user choice of mixing
!  method.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: k   ! depth level index

   real (POP_r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TMIX     ! tracers at mix time level

   real (POP_r8), dimension(nx_block,ny_block,km), intent(in) :: &
      UMIX, VMIX   ! U,V velocities at mix time level

   type (block), intent(in) :: &
      this_block           ! block information for this subblock

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block,nt), intent(out) :: &
      HDTK                ! Hdiff(T) for nth tracer at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      bid,kk              ! local block id

   real (POP_r8), dimension(nx_block,ny_block) :: &
     WORK                 ! temporary to hold tavg field
   real (POP_r8), dimension(nx_block,ny_block,nt,km) :: &
      TDTK,HDTK_BUF       ! Hdiff(T) for nth tracer at level k from submeso_flux code

   real (POP_r8) :: &
      start_time,end_time ! Timers


!-----------------------------------------------------------------------
!
!  branch to the proper mix routine
!
!-----------------------------------------------------------------------

   !if(my_task == master_task)then
   !call flush(6)
   !print *,"hi"
   !call flush(6)
   !endif   

   bid = this_block%local_id


   HDTK = c0

      if (k == 1) then
        !start_time = omp_get_wtime() 
        call tracer_diffs_and_isopyc_slopes_unified(TMIX, this_block)
        !end_time = omp_get_wtime()
        !print *,"time at tracer_diffs 1 is ",end_time - start_time   
      endif

      !if (k == 1) then
         !start_time = omp_get_wtime()

         call hdifft_gm_unified(1, HDTK_BUF(:,:,:,1), TMIX, UMIX, VMIX, tavg_HDIFE_TRACER, &
                         tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

         !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)num_threads(59) 
         !do kk=2,km
         !call hdifft_gm(kk , HDTK_BUF(:,:,:,kk) , TMIX, UMIX,VMIX,tavg_HDIFE_TRACER, &
         !                        tavg_HDIFN_TRACER,tavg_HDIFB_TRACER,this_block)
         !enddo

         !end_time = omp_get_wtime()

         !print *,"time at hdifft_gm combined is ",end_time - start_time 
                    
      !endif

      !start_time = omp_get_wtime()  
      !HDTK = HDTK_BUF(:,:,:,k)
      !end_time = omp_get_wtime()

      !print *,"time at hdifft_gm is ",end_time - start_time
 
        if (k == 1) then
         !start_time = omp_get_wtime()
         call submeso_sf_unified(TMIX, this_block)
         !end_time = omp_get_wtime()
         !print *,"time at submeso_sf is ",end_time - start_time
        endif

        if(k==1) then
         !start_time = omp_get_wtime()
        !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)num_threads(60) 
        do kk=1,km
         call submeso_flux_unified(kk, TDTK(:,:,:,kk), TMIX, tavg_HDIFE_TRACER, &
                              tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)
        enddo
        !end_time = omp_get_wtime()
        !print *,"time at submeso_flux is ",end_time - start_time
        endif
        HDTK=HDTK+TDTK(:,:,:,k)
   

!-----------------------------------------------------------------------
!EOC

 end subroutine hdifft_unified

   !dir$ attributes offload : mic :: tracer_diffs_and_isopyc_slopes_unified
   subroutine tracer_diffs_and_isopyc_slopes_unified (TMIX, this_block)

! !DESCRIPTION:
!  Calculates common variables used in hmix_gm and mix_submeso.
!
! !INPUT PARAMETERS:


      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level
      type (block), intent(in) :: &
         this_block            ! block info for this sub block


!-----------------------------------------------------------------------
!
!      local variables
!/TXP,
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ieast  = 1, iwest  = 2,       &
         jnorth = 1, jsouth = 2
      integer (int_kind) :: &
         i,j,n,kk,k,        &! dummy loop counters
         ktmp,              &! array indices
         kn, ks,            &! cyclic pointers for 2-level local arrays
         bid                 ! local block address for this sub block
      !real (r8), dimension(nx_block,ny_block) :: &
      !   KMASKE, KMASKN    ! ocean mask
        !DRDT, DRDS              ! expansion coefficients d(rho)/dT,S
      real (r8), dimension(nx_block,ny_block,2) :: &
         TXP, TYP, TZP , TEMP

      real (r8), dimension(nx_block,ny_block) :: & 
         RZ                  ! Dz(rho)
      integer (int_kind), parameter :: &
         ktp = 1, kbt = 2     ! refer to the top and bottom halves of a 
                              ! grid cell, respectively

      real (r8), dimension(nx_block,ny_block,km) :: &
         DRDT, DRDS                ! expansion coefficients d(rho)/dT,S

      real (r8) :: tempi,tempip1,tempj,tempjp1
 
      real (r8) :: temp_ksi,temp_ksip1,temp_ksj,temp_ksjp1,kmask,kmaske,kmaskn

      real (r8) :: txpim1,kmaskeim1,temp_ksim1,txim1

      real (r8) :: typjm1,kmasknjm1,temp_ksjm1,tyjm1        

      real (r8) start_time,end_time

      logical(log_kind) :: match
      
!-----------------------------------------------------------------------
!
!  register tracer_diffs_and_isopyc_slopes
!
!-----------------------------------------------------------------------

      match = .true.

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      DRDT   = c0
      DRDS   = c0
      TXP    = c0
      TYP    = c0
      TZP    = c0
      TEMP   = c0


        kn = 1
        ks = 2

        !start_time = omp_get_wtime()  

        !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)NUM_THREADS(60) 
        do kk=1,km

        call state_unified (kk, kk, TMIX(:,:,kk,1), TMIX(:,:,kk,2),  &
                           this_block, DRHODT=DRDT(:,:,kk), DRHODS=DRDS(:,:,kk))

        enddo 

        !end_time = omp_get_wtime()

        !print *,"Time at first part is",end_time - start_time

        kk=1

            !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i,KMASKE,KMASKN,tempi,tempip1,tempj,tempjp1)num_threads(60)
            do j=1,ny_block
              do i=1,nx_block


!-----------------------------------------------------------------------
!
!     compute RX=Dx(rho) and RY=Dy(rho) for all vertical levels. 
!
!-----------------------------------------------------------------------

                if ( kk <= KMT_unified(i,j,bid) .and. kk <= KMTE_unified(i,j,bid) ) then
                  KMASKE = c1
                else
                  KMASKE = c0
                endif
                if ( kk <= KMT_unified(i,j,bid) .and. kk <= KMTN_unified(i,j,bid) ) then
                  KMASKN = c1
                else
                  KMASKN = c0
                endif
                TEMP(i,j,kn) = max(-c2, TMIX(i,j,kk,1))

                if(i <= nx_block-1) then 
                 tempi = max(-c2, TMIX(i,j,kk,1))
                 tempip1 = max(-c2, TMIX(i+1,j,kk,1))
                 TXP(i,j,kn) = KMASKE * (tempip1  &
                                            -tempi)
                 endif 

                if(j <= ny_block-1)then
                 tempjp1 = max(-c2, TMIX(i,j+1,kk,1))
                 tempj = max(-c2, TMIX(i,j,kk,1))                
                 TYP(i,j,kn) = KMASKN * (tempjp1  &
                                             -tempj )
                endif   

                do n=1,nt
                  if(i <= nx_block-1)&
                  TX_UNIFIED(i,j,kk,n,bid) = KMASKE  &
                              * (TMIX(i+1,j,kk,n) - TMIX(i,j,kk,n))


                   if(j <= ny_block-1)& 
                       TY_UNIFIED(i,j,kk,n,bid) = KMASKN  &
                              * (TMIX(i,j+1,kk,n) - TMIX(i,j,kk,n))
                enddo

                RX_UNIFIED(i,j,ieast ,kk,bid) = DRDT(i,j,kk) * TXP(i,j,kn)  &
                                         + DRDS(i,j,kk) * TX_UNIFIED(i,j,kk,2,bid) 

                RY_UNIFIED(i,j,jnorth,kk,bid) = DRDT(i,j,kk) * TYP(i,j,kn)  &
                                         + DRDS(i,j,kk) * TY_UNIFIED(i,j,kk,2,bid) 


                 if(i >= 2) then
                   RX_UNIFIED(i,j,iwest,kk,bid) = DRDT(i,j,kk) * TXP(i-1,j,kn)  &
                                     + DRDS(i,j,kk) * TX_UNIFIED (i-1,j,kk,2,bid)
                 endif 

                  if(j >= 2)then 
                     RY_UNIFIED(i,j,jsouth,kk,bid) = DRDT(i,j,kk) * TYP(i,j-1,kn)  &
                                      + DRDS(i,j,kk) * TY_UNIFIED (i,j-1,kk,2,bid)
                  endif 
              enddo
            enddo


!-------------------------------------------------------------------------
!
!
!      when kk is 1 ends
!
!
!-------------------------------------------------------------------------

        do kk=1,km

            if ( kk < km ) then

            !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i,temp_ksi,temp_ksip1,temp_ksj,temp_ksjp1,kmask,kmaske,kmaskn,temp_ksim1,kmaskeim1) &
            !!$OMP PRIVATE(txpim1,txim1,temp_ksjm1,kmasknjm1,typjm1,tyjm1)NUM_THREADS(60)
            do j=1,ny_block
              do i=1,nx_block
                 KMASK = merge(c1, c0, kk < KMT_unified(i,j,bid))


!-----------------------------------------------------------------------
!
!     compute RX=Dx(rho) and RY=Dy(rho) for all vertical levels. 
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     compute RZ=Dz(rho) and
!     SLX = RX / RZ = slope of isopycnal surfaces in x-direction
!     SLY = RY / RZ = slope of isopycnal surfaces in y-direction
!
!-----------------------------------------------------------------------


                 TEMP(i,j,ks) = max(-c2, TMIX(i,j,kk+1,1))
              
                 TZ_UNIFIED(i,j,kk+1,1,bid) = TMIX(i,j,kk  ,1) - TMIX(i,j,kk+1,1)
                 TZ_UNIFIED(i,j,kk+1,2,bid) = TMIX(i,j,kk  ,2) - TMIX(i,j,kk+1,2) 
                 TZP(i,j,ks) = TEMP(i,j,kn) - TEMP(i,j,ks)


!     RZ = Dz(rho) = DRDT*Dz(T) + DRDS*Dz(S)


                 RZ(i,j) = DRDT(i,j,kk) * TZP(i,j,ks) + DRDS(i,j,kk) * TZ_UNIFIED (i,j,kk+1,2,bid) 
                 RZ(i,j) = min(RZ(i,j),-eps2)

         
                 !if (match) then 


                    SLX_UNIFIED(i,j,ieast ,kbt,kk,bid) = KMASK * RX_UNIFIED(i,j,ieast ,kk,bid) / RZ(i,j)
                    SLX_UNIFIED(i,j,iwest ,kbt,kk,bid) = KMASK * RX_UNIFIED(i,j,iwest ,kk,bid) / RZ(i,j)
                    SLY_UNIFIED(i,j,jnorth,kbt,kk,bid) = KMASK * RY_UNIFIED(i,j,jnorth,kk,bid) / RZ(i,j)
                    SLY_UNIFIED(i,j,jsouth,kbt,kk,bid) = KMASK * RY_UNIFIED(i,j,jsouth,kk,bid) / RZ(i,j)


                 !endif

!-----------------------------------------------------------------------
!
!     compute Dx(rho), Dy(rho) at level kk+1
!
!-----------------------------------------------------------------------

                 KMASKE = merge(c1, c0, kk+1 <= KMT_unified(i,j,bid) .and.  &
                                kk+1 <= KMTE_unified(i,j,bid))

                 KMASKN = merge(c1, c0, kk+1 <= KMT_unified(i,j,bid) .and.  &
                                kk+1 <= KMTN_unified(i,j,bid))



                 temp_ksi = max(-c2, TMIX(i,j,kk+1,1))
                 temp_ksip1 = max(-c2, TMIX(i+1,j,kk+1,1))

                 temp_ksj = max(-c2, TMIX(i,j,kk+1,1))
                 temp_ksjp1 = max(-c2, TMIX(i,j+1,kk+1,1))
                 

                 if(i <= nx_block-1)then

                  TXP(i,j,ks) = KMASKE*(temp_ksip1  &
                                               - temp_ksi) 
                 endif 

                 if(j <= ny_block-1)then
                  TYP(i,j,ks) = KMASKN*(temp_ksjp1  &
                                            - temp_ksj)
                 
                 endif 

                 do n=1,nt
                  if(i <= nx_block-1)&
                  TX_UNIFIED(i,j,kk+1,n,bid) = KMASKE  &
                            * (TMIX(i+1,j,kk+1,n) - TMIX(i,j,kk+1,n))

                  if(j <= ny_block-1)&
                  TY_UNIFIED(i,j,kk+1,n,bid) = KMASKN  &
                            * (TMIX(i,j+1,kk+1,n) - TMIX(i,j,kk+1,n))
                 enddo

                 RX_UNIFIED(i,j,ieast ,kk+1,bid) = DRDT(i,j,kk+1) * TXP(i,j,ks)  &
                                         + DRDS(i,j,kk+1) * TX_UNIFIED(i,j,kk+1,2,bid) 
                 RY_UNIFIED(i,j,jnorth,kk+1,bid) = DRDT(i,j,kk+1) * TYP(i,j,ks)  &
                                         + DRDS(i,j,kk+1) * TY_UNIFIED(i,j,kk+1,2,bid) 

                 if(i >= 2)then

                 temp_ksim1 = max(-c2, TMIX(i-1,j,kk+1,1))

                 kmaskeim1 = merge(c1, c0, kk+1 <= KMT_unified(i-1,j,bid) .and.  &
                                kk+1 <= KMTE_unified(i-1,j,bid))

                 txpim1 = kmaskeim1 * (temp_ksi - temp_ksim1 ) 

                 txim1 = kmaskeim1 * (TMIX(i,j,kk+1,2) - TMIX(i-1,j,kk+1,2))                 
                   
                 RX_UNIFIED(i,j,iwest,kk+1,bid) = DRDT(i,j,kk+1) * txpim1  &
                                       + DRDS(i,j,kk+1) *  txim1

                 endif


                 if(j >= 2) then

                 temp_ksjm1 = max(-c2, TMIX(i,j-1,kk+1,1))

                 kmasknjm1 = merge(c1, c0, kk+1 <= KMT_unified(i,j-1,bid) .and.  &
                                kk+1 <= KMTN_unified(i,j-1,bid))

                 typjm1 = kmasknjm1 * (temp_ksj - temp_ksjm1 )

                 tyjm1 = kmasknjm1 * (TMIX(i,j,kk+1,2) - TMIX(i,j-1,kk+1,2))


                 RY_UNIFIED(i,j,jsouth,kk+1,bid) = DRDT(i,j,kk+1) * typjm1  &
                                        + DRDS(i,j,kk+1) * tyjm1

                 endif


                 RZ(i,j) = DRDT(i,j,kk+1) * TZP(i,j,ks) + DRDS(i,j,kk+1) * TZ_UNIFIED(i,j,kk+1,2,bid) 
                 RZ_SAVE_UNIFIED(i,j,kk+1,bid) = min(RZ(i,j),c0)
                 RZ(i,j) = min(RZ(i,j),-eps2)


            !if (match) then

!-----------------------------------------------------------------------
!
!     compute slope of isopycnal surfaces at level kk+1
!
!-----------------------------------------------------------------------

              if ( kk+1 <= KMT_unified(i,j,bid) ) then
                SLX_UNIFIED(i,j,ieast, ktp,kk+1,bid) = RX_UNIFIED(i,j,ieast ,kk+1,bid) / RZ(i,j)
                SLX_UNIFIED(i,j,iwest, ktp,kk+1,bid) = RX_UNIFIED(i,j,iwest ,kk+1,bid) / RZ(i,j)
                SLY_UNIFIED(i,j,jnorth,ktp,kk+1,bid) = RY_UNIFIED(i,j,jnorth,kk+1,bid) / RZ(i,j)
                SLY_UNIFIED(i,j,jsouth,ktp,kk+1,bid) = RY_UNIFIED(i,j,jsouth,kk+1,bid) / RZ(i,j)
              endif


           !endif 

              enddo
            enddo


!-----------------------------------------------------------------------
!
!     end of kk < km if block
!
!-----------------------------------------------------------------------

          endif

          ktmp = kn
          kn   = ks
          ks   = ktmp

        enddo   ! end of kk-loop

        if(k==1)then
          !if(.not. registry_match('init_gm')) then
          do n=3,nt 
           !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk,n,j,i)collapse(3)num_threads(60)
           do kk=1,km-1
             do j=1,ny_block
                do i=1,nx_block
                   TZ_UNIFIED(i,j,kk+1,n,bid) = TMIX(i,j,kk  ,n) - TMIX(i,j,kk+1,n)
                enddo
             enddo
            enddo
           enddo
          !endif
         endif



        !end_time = omp_get_wtime()

        !print *,"Time taken at second time",end_time - start_time

!-----------------------------------------------------------------------
!
     end subroutine tracer_diffs_and_isopyc_slopes_unified


 !dir$ attributes offload:mic :: state_unified
 subroutine state_unified(k, kk, TEMPK, SALTK, this_block, &
                         RHOOUT, RHOFULL, DRHODT, DRHODS)

! !DESCRIPTION:
!  Returns the density of water at level k from equation of state
!  $\rho = \rho(d,\theta,S)$ where $d$ is depth, $\theta$ is
!  potential temperature, and $S$ is salinity. the density can be
!  returned as a perturbation (RHOOUT) or as the full density
!  (RHOFULL). Note that only the polynomial EOS choice will return
!  a perturbation density; in other cases the full density is returned
!  regardless of which argument is requested.
!
!  This routine also computes derivatives of density with respect
!  to temperature and salinity at level k from equation of state
!  if requested (ie the optional arguments are present).
!
!  If $k = kk$ are equal the density for level k is returned.
!  If $k \neq kk$ the density returned is that for a parcel
!  adiabatically displaced from level k to level kk.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k,                    &! depth level index
      kk                     ! level to which water is adiabatically 
                            ! displaced

   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      TEMPK,             &! temperature at level k
      SALTK               ! salinity    at level k

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), optional, intent(out) :: & 
      RHOOUT,  &! perturbation density of water
      RHOFULL, &! full density of water
      DRHODT,  &! derivative of density with respect to temperature
      DRHODS    ! derivative of density with respect to salinity

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ib,ie,jb,je,       &! extent of physical domain
      bid,               &! local block index
      out_of_range        ! counter for out-of-range T,S values

   real (r8), dimension(nx_block,ny_block) :: &
      TQ,SQ,             &! adjusted T,S
      BULK_MOD,          &! Bulk modulus
      RHO_S,             &! density at the surface
      DRDT0,             &! d(density)/d(temperature), for surface
      DRDS0,             &! d(density)/d(salinity   ), for surface
      DKDT,              &! d(bulk modulus)/d(pot. temp.)
      DKDS,              &! d(bulk modulus)/d(salinity  )
      SQR,DENOMK,        &! work arrays
      WORK1, WORK2, WORK3, WORK4, T2

   real (r8) :: p, p2 ! temporary pressure scalars

   !*** MWJF numerator coefficients including pressure

   real (r8) ::                                                        &
      mwjfnums0t0, mwjfnums0t1, mwjfnums0t2, mwjfnums0t3,              &
      mwjfnums1t0, mwjfnums1t1, mwjfnums2t0,                           &
      mwjfdens0t0, mwjfdens0t1, mwjfdens0t2, mwjfdens0t3, mwjfdens0t4, &
      mwjfdens1t0, mwjfdens1t1, mwjfdens1t3,                           &
      mwjfdensqt0, mwjfdensqt2

!-----------------------------------------------------------------------
!
!  first check for valid range if requested
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

      TQ = TEMPK
      SQ = SALTK

      TQ = min(TEMPK,tmax(kk))
      TQ = max(TQ,tmin(kk))

      SQ = min(SALTK,smax(kk))
      SQ = max(SQ,smin(kk))

!-----------------------------------------------------------------------
!
!  now compute density or expansion coefficients
!
!-----------------------------------------------------------------------

      p   = c10*pressz_unified(kk)

      SQ  = c1000*SQ
#ifdef CCSMCOUPLED
      SQR = sqrt(SQ) 
#else
      SQR = sqrt(SQ)
#endif

      !***
      !*** first calculate numerator of MWJF density [P_1(S,T,p)]
      !***

      mwjfnums0t0 = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0)
      mwjfnums0t1 = mwjfnp0s0t1 
      mwjfnums0t2 = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2)
      mwjfnums0t3 = mwjfnp0s0t3
      mwjfnums1t0 = mwjfnp0s1t0 + p*mwjfnp1s1t0
      mwjfnums1t1 = mwjfnp0s1t1
      mwjfnums2t0 = mwjfnp0s2t0

      WORK1 = mwjfnums0t0 + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2 + &
              mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0 +              &
              mwjfnums1t1 * TQ + mwjfnums2t0 * SQ)

      !***
      !*** now calculate denominator of MWJF density [P_2(S,T,p)]
      !***

      mwjfdens0t0 = mwjfdp0s0t0 + p*mwjfdp1s0t0
      mwjfdens0t1 = mwjfdp0s0t1 + p**3 * mwjfdp3s0t1
      mwjfdens0t2 = mwjfdp0s0t2
      mwjfdens0t3 = mwjfdp0s0t3 + p**2 * mwjfdp2s0t3
      mwjfdens0t4 = mwjfdp0s0t4
      mwjfdens1t0 = mwjfdp0s1t0
      mwjfdens1t1 = mwjfdp0s1t1
      mwjfdens1t3 = mwjfdp0s1t3
      mwjfdensqt0 = mwjfdp0sqt0
      mwjfdensqt2 = mwjfdp0sqt2

      WORK2 = mwjfdens0t0 + TQ * (mwjfdens0t1 + TQ * (mwjfdens0t2 +    &
           TQ * (mwjfdens0t3 + mwjfdens0t4 * TQ))) +                   &
           SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+ &
           SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2))

      DENOMK = c1/WORK2

      if (present(RHOOUT)) then
         RHOOUT  = WORK1*DENOMK
      endif

      if (present(RHOFULL)) then
         RHOFULL = WORK1*DENOMK
      endif

      if (present(DRHODT)) then
         WORK3 = &! dP_1/dT
                 mwjfnums0t1 + TQ * (c2*mwjfnums0t2 +    &
                 c3*mwjfnums0t3 * TQ) + mwjfnums1t1 * SQ

         WORK4 = &! dP_2/dT
                 mwjfdens0t1 + SQ * mwjfdens1t1 +               &
                 TQ * (c2*(mwjfdens0t2 + SQ*SQR*mwjfdensqt2) +  &
                 TQ * (c3*(mwjfdens0t3 + SQ * mwjfdens1t3) +    &
                 TQ *  c4*mwjfdens0t4))

         DRHODT = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK
      endif

      if (present(DRHODS)) then
         WORK3 = &! dP_1/dS
                 mwjfnums1t0 + mwjfnums1t1 * TQ + c2*mwjfnums2t0 * SQ

         WORK4 = mwjfdens1t0 +   &! dP_2/dS
                 TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3) +   &
                 c1p5*SQR*(mwjfdensqt0 + TQ*TQ*mwjfdensqt2)

         DRHODS = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK * c1000
      endif

!-----------------------------------------------------------------------
!EOC

 end subroutine state_unified

   !dir$ attributes offload:mic :: submeso_sf_unified
   subroutine submeso_sf_unified ( TMIX, this_block )

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

                    SF_SUBM_X_unified(i,j,temp,kk,k,bid) = c0
                    SF_SUBM_Y_unified(i,j,temp,kk,k,bid) = c0

                 enddo
              enddo
            enddo
          enddo
      enddo

   !print *,"base time is",end_time - start_time

   ML_DEPTH(:,:) = HMXL_UNIFIED(:,:,bid) 
   

   do j=1,ny_block
      do i=1,nx_block
          CONTINUE_INTEGRAL(i,j) = .true.
          if( KMT_UNIFIED(i,j,bid) == 0 ) then
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
     if ( k > 1 )  zw_top = zw_unified(k-1)

    !!$OMP PARALLEL DO SHARED(CONTINUE_INTEGRAL,BX_VERT_AVG,RX_UNIFIED,RY_UNIFIED,ML_DEPTH)PRIVATE(i,WORK3)num_threads(60)SCHEDULE(dynamic,16)
    do j=1,ny_block
        do i=1,nx_block

            WORK3(i,j)=c0
            if( CONTINUE_INTEGRAL(i,j)  .and.  ML_DEPTH(i,j) > zw_unified(k) )then
               WORK3(i,j) = dz_unified(k)
           endif 

            if( CONTINUE_INTEGRAL(i,j)  .and.  ML_DEPTH(i,j) <= zw_unified(k)  &
             .and.  ML_DEPTH(i,j) > zw_top ) then
                    WORK3(i,j) = ML_DEPTH(i,j) - zw_top
            endif

            if ( CONTINUE_INTEGRAL(i,j) ) then
                 BX_VERT_AVG(i,j,1) = BX_VERT_AVG(i,j,1)        &
                                      + RX_UNIFIED(i,j,1,k,bid) * WORK3(i,j)
                 BX_VERT_AVG(i,j,2) = BX_VERT_AVG(i,j,2)        &
                                      + RX_UNIFIED(i,j,2,k,bid) * WORK3(i,j)
                 BY_VERT_AVG(i,j,1) = BY_VERT_AVG(i,j,1)        &
                                      + RY_UNIFIED(i,j,1,k,bid) * WORK3(i,j)
                 BY_VERT_AVG(i,j,2) = BY_VERT_AVG(i,j,2)        &
                                      + RY_UNIFIED(i,j,2,k,bid) * WORK3(i,j)
            endif

            if ( CONTINUE_INTEGRAL(i,j) .and.  ML_DEPTH(i,j) <= zw_unified(k)  &
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

           if ( KMT_UNIFIED(i,j,bid) > 0 ) then
           BX_VERT_AVG(i,j,1) = - grav * BX_VERT_AVG(i,j,1) / ML_DEPTH(i,j)
           BX_VERT_AVG(i,j,2) = - grav * BX_VERT_AVG(i,j,2) / ML_DEPTH(i,j)
           BY_VERT_AVG(i,j,1) = - grav * BY_VERT_AVG(i,j,1) / ML_DEPTH(i,j)
           BY_VERT_AVG(i,j,2) = - grav * BY_VERT_AVG(i,j,2) / ML_DEPTH(i,j)
           endif

        enddo
     enddo

!!!! if lsubmeso_scaling  false!!!!!!!

    do j=1,ny_block
        do i=1,nx_block
 
           WORK1(i,j)=c0
 
           if( KMT_UNIFIED(i,j,bid) > 0 )then
                 WORK1(i,j) = sqrt( p5 * (                          &
                              ( BX_VERT_AVG(i,j,1)**2 + BX_VERT_AVG(i,j,2)**2 )  &
                                / DXT_UNIFIED(i,j,bid)**2                                &
                              + ( BY_VERT_AVG(i,j,1)**2 + BY_VERT_AVG(i,j,2)**2 )  &
                                / DYT_UNIFIED(i,j,bid)**2 ) )
                  WORK1(i,j) = WORK1(i,j) * ML_DEPTH(i,j) * (TIME_SCALE_UNIFIED(i,j,bid)**2)
            endif


            CONTINUE_INTEGRAL(i,j) = .true.
            if( KMT_UNIFIED(i,j,bid) == 0 ) then
                        CONTINUE_INTEGRAL(i,j) = .false.
            endif
 
            WORK2(i,j) = c0


        enddo
     enddo

          
     do k=2,km
        !!$OMP PARALLEL DO SHARED(WORK3,CONTINUE_INTEGRAL,WORK2,k,bid,dzw_unified,zt_unified,dzwr_unified,RZ_SAVE_UNIFIED,ML_DEPTH)PRIVATE(i,j)DEFAULT(NONE)num_threads(60)
        do j=1,ny_block
           do i=1,nx_block

              WORK3(i,j) = c0
              if ( CONTINUE_INTEGRAL(i,j)  .and.  ML_DEPTH(i,j) > zt_unified(k) ) then
                   WORK3(i,j) = dzw_unified(k-1) 
              endif

             if ( CONTINUE_INTEGRAL(i,j)  .and.  ML_DEPTH(i,j) <= zt_unified(k)  &
                  .and.  ML_DEPTH(i,j) >= zt_unified(k-1) ) then
                  WORK3(i,j) = ( (ML_DEPTH(i,j) - zt_unified(k-1))**2 ) * dzwr_unified(k-1)
             endif

             if ( CONTINUE_INTEGRAL(i,j) ) then
             WORK2(i,j) = WORK2(i,j) + sqrt(-RZ_SAVE_UNIFIED(i,j,k,bid) * WORK3(i,j))
             endif

             if ( CONTINUE_INTEGRAL(i,j) .and.  ML_DEPTH(i,j) <= zt_unified(k)  &
                     .and.  ML_DEPTH(i,j) >= zt_unified(k-1) )then
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

             if ( KMT_unified(i,j,bid) > 0 ) then

                  WORK2(i,j) = sqrt_grav_unified * WORK2(i,j) * TIME_SCALE_unified(i,j,bid)

                  HLS(i,j) = max ( WORK1(i,j), WORK2(i,j), hor_length_scale_unified )

             endif

           enddo
        enddo

   do k=1,km

     reference_depth(ktp) = zt_unified(k) - p25 * dz_unified(k)
     reference_depth(kbt) = zt_unified(k) + p25 * dz_unified(k)

     do kk=ktp,kbt

       !!$OMP PARALLEL DO DEFAULT(NONE)PRIVATE(i,j)SHARED(reference_depth,ML_DEPTH,KMT,WORK3,WORK2,WORK1,TIME_SCALE_UNIFIED,HLS,SF_SUBM_X_UNIFIED,SF_SUBM_Y_UNIFIED) &
       !!$OMP SHARED(BX_VERT_AVG,BY_VERT_AVG,DXT_UNIFIED,DYT_UNIFIED,max_hor_grid_scale_unified,efficiency_factor_unified,kk,k,bid)num_threads(60)  
       do j=1,ny_block
           do i=1,nx_block


               if ( reference_depth(kk) < ML_DEPTH(i,j)  .and.  &
                    KMT_UNIFIED(i,j,bid) >= k ) then

                        WORK3(i,j) = ( c1 - ( c2 * reference_depth(kk) / ML_DEPTH(i,j) ) )**2
            
                        WORK2(i,j) = ( c1 - WORK3(i,j) )  &
                                    * ( c1 + ( 5.0_r8 / 21.0_r8 ) * WORK3(i,j) )

                        WORK1(i,j) = efficiency_factor_unified * (ML_DEPTH(i,j)**2) * WORK2(i,j)  &
                                     * TIME_SCALE_unified(i,j,bid) / HLS(i,j)

!     in the following negative sign is omitted to be consistent with
!     the GM implementation in hmix_gm subroutine. also, DXT and
!     DYT usage is approximate. 

                        SF_SUBM_X_UNIFIED(i,j,1,kk,k,bid) = WORK1(i,j) * BX_VERT_AVG(i,j,1)  &
                                                          * min(DXT_UNIFIED(i,j,bid),max_hor_grid_scale_unified)
                        SF_SUBM_X_UNIFIED(i,j,2,kk,k,bid) = WORK1(i,j) * BX_VERT_AVG(i,j,2)  &
                                                          * min(DXT_UNIFIED(i,j,bid),max_hor_grid_scale_unified)
                        SF_SUBM_Y_UNIFIED(i,j,1,kk,k,bid) = WORK1(i,j) * BY_VERT_AVG(i,j,1)  &
                                                          * min(DYT_UNIFIED(i,j,bid),max_hor_grid_scale_unified)
                        SF_SUBM_Y_UNIFIED(i,j,2,kk,k,bid) = WORK1(i,j) * BY_VERT_AVG(i,j,2)  &
                                                          * min(DYT_UNIFIED(i,j,bid),max_hor_grid_scale_unified)

               endif

           enddo
        enddo


     enddo

   enddo


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

         WORK1(i,j) = (   SF_SUBM_X_UNIFIED(i  ,j,1,kbt,k,  bid)    &
               + factor * SF_SUBM_X_UNIFIED(i  ,j,1,ktp,kp1,bid)    &
                        + SF_SUBM_X_UNIFIED(i+1,j,2,kbt,k,  bid)    &
               + factor * SF_SUBM_X_UNIFIED(i+1,j,2,ktp,kp1,bid) )  &
                * p25 * HYX_UNIFIED(i,j,bid)

         WORK2(i,j) = (   SF_SUBM_Y_UNIFIED(i,j  ,1,kbt,k,  bid)    &
               + factor * SF_SUBM_Y_UNIFIED(i,j  ,1,ktp,kp1,bid)    &
                        + SF_SUBM_Y_UNIFIED(i,j+1,2,kbt,k,  bid)    &
               + factor * SF_SUBM_Y_UNIFIED(i,j+1,2,ktp,kp1,bid) )  &
               * p25 * HXY_UNIFIED(i,j,bid)

         endif


           USMB(i,j) = merge( WORK1(i,j), c0, k < KMT_UNIFIED(i,j,bid) .and. k < KMTE_UNIFIED(i,j,bid) )
           VSMB(i,j) = merge( WORK2(i,j), c0, k < KMT_UNIFIED(i,j,bid) .and. k < KMTN_UNIFIED(i,j,bid) )

            WORK1(i,j) = merge( USMT(i,j) - USMB(i,j), c0, k <= KMT_UNIFIED(i,j,bid)  &
                                      .and. k <= KMTE_UNIFIED(i,j,bid) )

            WORK2(i,j) = merge( VSMT(i,j) - VSMB(i,j), c0, k <= KMT_UNIFIED(i,j,bid)  &
                                      .and. k <= KMTN_UNIFIED(i,j,bid) )

            U_SUBM(i,j) = WORK1(i,j) * dzr_unified(k) / HTE_UNIFIED(i,j,bid)
            V_SUBM(i,j) = WORK2(i,j) * dzr_unified(k) / HTN_UNIFIED(i,j,bid)

       enddo
     enddo



     !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60)
     do j=this_block%jb,this_block%je
       do i=this_block%ib,this_block%ie

         if ( k < KMT_UNIFIED(i,j,bid)  .and.  ( zw_unified(k) < max( ML_DEPTH(i,j),  &
              ML_DEPTH(i+1,j), ML_DEPTH(i-1,j), ML_DEPTH(i,j+1),      &
              ML_DEPTH(i,j-1) ) ) ) then
           WBOT_SUBM(i,j) = WTOP_SUBM(i,j)             &
                           + TAREA_R_UNIFIED(i,j,bid)          &
                       * ( WORK1(i,j) - WORK1(i-1,j)   &
                         + WORK2(i,j) - WORK2(i,j-1) )
         else
           WBOT_SUBM(i,j) = c0
         endif

       enddo
     enddo

     USMT = USMB
     VSMT = VSMB

     WTOP_SUBM = WBOT_SUBM

   enddo


 end subroutine submeso_sf_unified

 !dir$ attributes offload:mic :: submeso_flux_unified
 subroutine submeso_flux_unified (k, GTK, TMIX, tavg_HDIFE_TRACER, &
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

      real (r8) :: WORK1prev,WORK2prev,KMASKprev,fzprev

      CX = merge(HYX_UNIFIED(:,:,bid)*p25, c0, (k <= KMT_UNIFIED (:,:,bid))   &
                                 .and. (k <= KMTE_UNIFIED(:,:,bid)))
      CY = merge(HXY_UNIFIED(:,:,bid)*p25, c0, (k <= KMT_UNIFIED (:,:,bid))   &
                                 .and. (k <= KMTN_UNIFIED(:,:,bid)))
      
      KMASK = merge(c1, c0, k < KMT_UNIFIED(:,:,bid))
            
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

              FX(i,j,n) = CX(i,j)                          &
               * ( SF_SUBM_X_UNIFIED(i  ,j,ieast,ktp,k,bid) * TZ_UNIFIED(i,j,k,n,bid)                      &
                 + SF_SUBM_X_UNIFIED(i  ,j,ieast,kbt,k,bid) * TZ_UNIFIED(i,j,kp1,n,bid)                    &
                 + SF_SUBM_X_UNIFIED(i+1,j,iwest,ktp,k,bid) * TZ_UNIFIED(i+1,j,k,n,bid)                    &
                 + SF_SUBM_X_UNIFIED(i+1,j,iwest,kbt,k,bid) * TZ_UNIFIED(i+1,j,kp1,n,bid) )

              endif  

              if(j <= ny_block -1 )then

              FY(i,j,n) =  CY(i,j)                          &
               * ( SF_SUBM_Y_UNIFIED(i,j  ,jnorth,ktp,k,bid) * TZ_UNIFIED(i,j,k,n,bid)                      &
                 + SF_SUBM_Y_UNIFIED(i,j  ,jnorth,kbt,k,bid) * TZ_UNIFIED(i,j,kp1,n,bid)                    &
                 + SF_SUBM_Y_UNIFIED(i,j+1,jsouth,ktp,k,bid) * TZ_UNIFIED(i,j+1,k,n,bid)                    &
                 + SF_SUBM_Y_UNIFIED(i,j+1,jsouth,kbt,k,bid) * TZ_UNIFIED(i,j+1,kp1,n,bid) )

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

                  WORK1(i,j) = SF_SUBM_X_UNIFIED(i  ,j  ,ieast ,kbt,k  ,bid)     &
                             * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k  ,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jnorth,kbt,k  ,bid)     &
                             * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k  ,n,bid)  &
                             + SF_SUBM_X_UNIFIED(i  ,j  ,iwest ,kbt,k  ,bid)     &
                             * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k  ,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jsouth,kbt,k  ,bid)     &
                             * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k  ,n,bid)


                  WORK2(i,j) = factor                                    &
                           * ( SF_SUBM_X_UNIFIED(i  ,j  ,ieast ,ktp,kp1,bid)     &
                             * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,kp1,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jnorth,ktp,kp1,bid)     &
                             * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,kp1,n,bid)  &
                             + SF_SUBM_X_UNIFIED(i  ,j  ,iwest ,ktp,kp1,bid)     &
                             * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,kp1,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jsouth,ktp,kp1,bid)     &
                             * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,kp1,n,bid) ) 
                   
                  if(k==1) then
 
                  fzprev = 0
    
                  else    

    
                  WORK1prev = SF_SUBM_X_UNIFIED(i  ,j  ,ieast ,kbt,k-1 ,bid)     &
                             * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k-1,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jnorth,kbt,k-1,bid)     &
                             * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k-1,n,bid)  &
                             + SF_SUBM_X_UNIFIED(i  ,j  ,iwest ,kbt,k-1,bid)     &
                             * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k-1,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jsouth,kbt,k-1,bid)     &
                             * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k-1,n,bid)

                  WORK2prev = factor &
                           * ( SF_SUBM_X_UNIFIED(i  ,j  ,ieast ,ktp,kp1-1,bid)     &
                             * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,kp1-1,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jnorth,ktp,kp1-1,bid)     &
                             * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,kp1-1,n,bid)  &
                             + SF_SUBM_X_UNIFIED(i  ,j  ,iwest ,ktp,kp1-1,bid)     &
                             * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,kp1-1,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jsouth,ktp,kp1-1,bid)     &
                             * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,kp1-1,n,bid) )

                  KMASKprev = merge(c1, c0, k-1 < KMT_UNIFIED(i,j,bid))


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


                  GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                               + FY(i,j,n) - FY(i,j-1,n)  &
                        + fzprev - fz )*dzr_unified(k)*TAREA_R_UNIFIED(i,j,bid)

                  !FZTOP_SUBM(i,j,n,bid) = fz

               else  

                  WORK1prev = SF_SUBM_X_UNIFIED(i  ,j  ,ieast ,kbt,k-1 ,bid)     &
                             * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k-1,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jnorth,kbt,k-1,bid)     &
                             * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k-1,n,bid)  &
                             + SF_SUBM_X_UNIFIED(i  ,j  ,iwest ,kbt,k-1,bid)     &
                             * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k-1,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jsouth,kbt,k-1,bid)     &
                             * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k-1,n,bid)

                  WORK2prev = factor &
                           * ( SF_SUBM_X_UNIFIED(i  ,j  ,ieast ,ktp,km,bid)     &
                             * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,km,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jnorth,ktp,km,bid)     &
                             * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,km,n,bid)  &
                             + SF_SUBM_X_UNIFIED(i  ,j  ,iwest ,ktp,km,bid)     &
                             * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,km,n,bid)  &
                             + SF_SUBM_Y_UNIFIED(i  ,j  ,jsouth,ktp,km,bid)     &
                             * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,km,n,bid) )

                  KMASKprev = merge(c1, c0, k-1 < KMT_UNIFIED(i,j,bid))


                  fzprev = -KMASKprev * p25 &
                            * (WORK1prev + WORK2prev)

                  GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                               + FY(i,j,n) - FY(i,j-1,n)  &
                     + fzprev )*dzr_unified(k)*TAREA_R_UNIFIED(i,j,bid)

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
              WORK1(i,j) = FX(i,j,n)*dzr_unified(k)*TAREA_R_UNIFIED(i,j,bid)
            enddo
            enddo
            !call accumulate_tavg_field(WORK1,tavg_HDIFE_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFN_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FY(i,j,n)*dzr_unified(k)*TAREA_R_UNIFIED(i,j,bid)
            enddo
            enddo
            !call accumulate_tavg_field(WORK1,tavg_HDIFN_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFB_TRACER(n))) then
            do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie
              WORK1(i,j) = FZTOP_SUBM_UNIFIED(i,j,n,bid)*dzr_unified(k)*TAREA_R_UNIFIED(i,j,bid)
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
     

 end subroutine submeso_flux_unified

 !***********************************************************************
!BOP
! !IROUTINE: hdifft_gm_unified
! !INTERFACE:

      !dir$ attributes offload:mic :: hdifft_gm_unified
      subroutine hdifft_gm_unified (k, GTK, TMIX, UMIX, VMIX, tavg_HDIFE_TRACER, &
                                  tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

! !DESCRIPTION:
!  Gent-McWilliams eddy transport parameterization
!  and isopycnal diffusion.
!
!  This routine must be called successively with k = 1,2,3,...
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      integer (int_kind), intent(in) :: k  ! depth level index

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level

      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         UMIX, VMIX            ! U,V  at all vertical levels
                               !   at mixing time level

      integer (int_kind), dimension(nt), intent(in) :: &
         tavg_HDIFE_TRACER, &! tavg id for east face diffusive flux of tracer
         tavg_HDIFN_TRACER, &! tavg id for north face diffusive flux of tracer
         tavg_HDIFB_TRACER   ! tavg id for bottom face diffusive flux of tracer

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

! !OUTPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,nt), intent(out) :: &
         GTK     ! diffusion+bolus advection for nth tracer at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ieast  = 1, iwest  = 2,       &
         jnorth = 1, jsouth = 2

      integer (int_kind) :: &
         i,j,n,kk,          &! dummy loop counters
         kid,ktmp,          &! array indices
         kk_sub, kp1,       & 
         bid                 ! local block address for this sub block

      real (r8) :: &
         fz, dz_bottom, factor,fzprev, KMASKprev, WORK3prev, dzbottomprev

      real (r8), dimension(nx_block,ny_block) :: &
         CX, CY,                  &
         RZ,                      &! Dz(rho)
         SLA,                     &! absolute value of slope
         WORK1, WORK2,            &! local work space
         WORK3, WORK4,            &! local work space
         KMASK,                   &! ocean mask
         TAPER1, TAPER2, TAPER3,  &! tapering factors
         UIB, VIB,                &! work arrays for isopycnal mixing velocities
         U_ISOP, V_ISOP            ! horizontal components of isopycnal velocities

      real (r8), dimension(nx_block,ny_block,nt) :: &
         FX, FY                     ! fluxes across east, north faces

      real (r8), dimension(2) :: &
         reference_depth

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      U_ISOP = c0
      V_ISOP = c0
      WORK1  = c0
      WORK2  = c0
      WORK3  = c0
      WORK4  = c0

      if ( .not. implicit_vertical_mix )  print *, "Error in hmix_gm if ( .not. implicit_vertical_mix )"

      if ( k == 1 ) then

        if ( diag_gm_bolus ) then
          UIB = c0
          VIB = c0
          UIT_UNIFIED(:,:,bid) = c0
          VIT_UNIFIED(:,:,bid) = c0
          WBOT_ISOP_UNIFIED(:,:,bid) = c0
        endif

        HOR_DIFF_UNIFIED(:,:,ktp,k,bid) = ah_bkg_srfbl_unified

        BL_DEPTH_UNIFIED(:,:,bid) = KPP_HBLT_UNIFIED(:,:,bid)

        call smooth_hblt_unified ( .false., .true., bid,  &
                     SMOOTH_OUT=TLT_UNIFIED%DIABATIC_DEPTH(:,:,bid) )

        


     endif ! if k == 1 

 end subroutine hdifft_gm_unified 

!***********************************************************************
!BOP
! !IROUTINE: smooth_hblt
! !INTERFACE:

 !dir$ attributes offload:mic :: smooth_hblt_unified
 subroutine smooth_hblt_unified (overwrite_hblt, use_hmxl, &
                         bid, HBLT, KBL, SMOOTH_OUT)

! !DESCRIPTION:
!  This subroutine uses a 1-1-4-1-1 Laplacian filter one time
!  on HBLT or HMXL to reduce any horizontal two-grid-point noise.
!  If HBLT is overwritten, KBL is adjusted after smoothing.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in) :: &
      overwrite_hblt,   &    ! if .true.,  HBLT is overwritten
                             ! if .false., the result is returned in
                             !  a dummy array
      use_hmxl               ! if .true., smooth HMXL
                             ! if .false., smooth HBLT

   integer (int_kind), intent(in) :: &
      bid                    ! local block address

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), optional, intent(inout) :: &
      HBLT                   ! boundary layer depth

   integer (int_kind), dimension(nx_block,ny_block), optional, intent(inout) :: &
      KBL                    ! index of first lvl below hbl

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), optional, intent(out) ::  &
      SMOOTH_OUT              ! optional output array containing the
                              !  smoothened field if overwrite_hblt is false

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   character (char_len) ::  &
      message

   integer (int_kind) :: &
      i, j,              &  ! horizontal loop indices
      k                     ! vertical level index

   real (r8), dimension(nx_block,ny_block) ::  &
      WORK1, WORK2

   real (r8) ::  &
     cc, cw, ce, cn, cs, &  ! averaging weights
     ztmp                   ! temp for level depth

   real (r8) start_time,end_time

!-----------------------------------------------------------------------
!
!     consistency checks 
!
!-----------------------------------------------------------------------

   if ( overwrite_hblt  .and.  ( .not.present(KBL)  .or.        &
                                 .not.present(HBLT) ) ) then      
     message = 'incorrect subroutine arguments for smooth_hblt, error # 1'
     print *,message
     !call exit_POP (sigAbort, trim(message))
   endif

   if ( .not.overwrite_hblt  .and.  .not.present(SMOOTH_OUT) ) then 
     print *,message 
     message = 'incorrect subroutine arguments for smooth_hblt, error # 2'
     !call exit_POP (sigAbort, trim(message))
   endif

   if ( use_hmxl .and. .not.present(SMOOTH_OUT) ) then          
     message = 'incorrect subroutine arguments for smooth_hblt, error # 3'
     print *,message
     !call exit_POP (sigAbort, trim(message))
   endif

   if ( overwrite_hblt  .and.  use_hmxl ) then                  
     message = 'incorrect subroutine arguments for smooth_hblt, error # 4'
     print *,message 
     !call exit_POP (sigAbort, trim(message))
   endif

!-----------------------------------------------------------------------
!
!     perform one smoothing pass since we cannot do the necessary 
!     boundary updates for multiple passes.
!
!-----------------------------------------------------------------------


   if ( use_hmxl ) then
     WORK2 = HMXL_UNIFIED(:,:,bid)
   else
     WORK2 = HBLT
   endif

   WORK1 = WORK2

   do j=2,ny_block-1
     do i=2,nx_block-1
       if ( KMT_UNIFIED(i,j,bid) /= 0 ) then
         cw = p125
         ce = p125
         cn = p125
         cs = p125
         cc = p5
         if ( KMT_UNIFIED(i-1,j,bid) == 0 ) then
           cc = cc + cw
           cw = c0
         endif
         if ( KMT_UNIFIED(i+1,j,bid) == 0 ) then
           cc = cc + ce
           ce = c0
         endif
         if ( KMT_UNIFIED(i,j-1,bid) == 0 ) then
           cc = cc + cs
           cs = c0
         endif
         if ( KMT_UNIFIED(i,j+1,bid) == 0 ) then
           cc = cc + cn
           cn = c0
         endif
         WORK2(i,j) =  cw * WORK1(i-1,j)   &
                     + ce * WORK1(i+1,j)   &
                     + cs * WORK1(i,j-1)   &
                     + cn * WORK1(i,j+1)   &
                     + cc * WORK1(i,j)
       endif
     enddo
   enddo


   do k=1,km
     !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(ztmp,j,i)NUM_THREADS(60)
     do j=2,ny_block-1
       do i=2,nx_block-1

           ztmp = -zgrid_unified(k)

         if ( k == KMT_UNIFIED(i,j,bid)  .and.  WORK2(i,j) > ztmp ) then
           WORK2(i,j) = ztmp
         endif

       enddo
     enddo
     !!$OMP END PARALLEL DO
   enddo

   if ( overwrite_hblt  .and.  .not.use_hmxl ) then

     HBLT = WORK2

     do k=1,km
       do j=2,ny_block-1
         do i=2,nx_block-1

             ztmp = -zgrid_unified(k)

           if ( KMT_UNIFIED(i,j,bid) /= 0            .and.  &
                ( HBLT(i,j) >  -zgrid_unified(k-1) ) .and.  &
                ( HBLT(i,j) <= ztmp        ) ) KBL(i,j) = k
     
         enddo
       enddo
     enddo

   else

     SMOOTH_OUT = WORK2

   endif
 

!-----------------------------------------------------------------------

 end subroutine smooth_hblt_unified

 end module horizontal_mix_unified

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
