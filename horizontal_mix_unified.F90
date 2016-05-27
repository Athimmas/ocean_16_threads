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
   ny_block_unified
   use distribution 
   use domain_size
   use domain, only: nblocks_clinic, distrb_clinic
   use constants
   use communicate, only: my_task, master_task
   use time_management, only: km, nt, mix_pass,eod_last,nsteps_total
   use broadcast, only: broadcast_scalar
   use io_types, only: nml_in, nml_filename, stdout
   use hmix_del2, only: init_del2u, init_del2t, hdiffu_del2, hdifft_del2
   use hmix_del4, only: init_del4u, init_del4t, hdiffu_del4, hdifft_del4
   use hmix_gm, only: diag_gm_bolus,kappa_isop_type,kappa_thic_type, &
                      transition_layer_on,ah_bolus,slope_control,diff_tapering,&
                      slm_r,slm_b,use_const_ah_bkg_srfbl,ah_bkg_srfbl,cancellation_occurs, &
                      kappa_freq,ah_bkg_bottom
   use hmix_aniso, only: init_aniso, hdiffu_aniso
   use topostress, only: ltopostress
   use horizontal_mix, only:tavg_HDIFE_TRACER,tavg_HDIFN_TRACER,tavg_HDIFB_TRACER
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now, &
      tavg_in_which_stream, ltavg_on
   use timers, only: timer_start, timer_stop, get_timer
   use exit_mod, only: sigAbort, exit_pop, flushm
   use mix_submeso, only: init_submeso, submeso_flux, submeso_sf
   use hmix_gm_submeso_share, only: init_meso_mixing, tracer_diffs_and_isopyc_slopes
   use vertical_mix, only: implicit_vertical_mix,vmix_itype,vmix_type_kpp
   use state_mod, only: pressz
   use omp_lib

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_horizontal_mix_unified, &
             hdifft_unified,tracer_diffs_and_isopyc_slopes_unified,smooth_hblt_unified,state_unified, &
             submeso_flux_unified,transition_layer_unified,buoyancy_frequency_dependent_profile_unified,&
             merged_streamfunction_unified,apply_vertical_profile_to_isop_hor_diff_unified,submeso_sf_unified

!EOP
!BOC

! !PUBLIC VARIABLES

!variables for horizontal_mix

   real (POP_r8), dimension(nx_block,ny_block,nt,km,max_blocks_clinic) :: &
      TDTK,HDTK_BUF      ! Hdiff(T) for nth tracer at level k from submeso_flux code

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

   !dir$ attributes offload : mic :: KMT_UNIFIED
   integer (POP_i4), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      KMT_UNIFIED

   !dir$ attributes offload : mic :: KMTN_UNIFIED
   !dir$ attributes offload : mic :: KMTE_UNIFIED
   integer (POP_i4), dimension(nx_block,ny_block,max_blocks_clinic), &
      public :: &
      KMTN_UNIFIED,KMTE_UNIFIED

   !dir$ attributes offload:mic :: DZT_UNIFIED
   real (POP_r8), dimension(:,:,:,:), allocatable, public :: &
      DZT_unified

   !dir$ attributes offload:mic :: DXT_UNIFIED
   !dir$ attributes offload:mic :: DYT_UNIFIED
   !dir$ attributes offload:mic :: HUS_UNIFIED
   !dir$ attributes offload:mic :: HUW_UNIFIED
   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic), public :: &
      DYT_UNIFIED, DXT_UNIFIED,HUS_UNIFIED,HUW_UNIFIED

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

      !dir$ attributes offload:mic :: KAPPA_LATERAL_UNIFIED
      real (r8), dimension(:,:,:), allocatable ,public :: &
         KAPPA_LATERAL_UNIFIED      ! horizontal variation of KAPPA in cm^2/s

      !dir$ attributes offload:mic :: KAPPA_VERTICAL_UNIFIED 
      real (r8), dimension(:,:,:,:), allocatable ,public :: &
         KAPPA_VERTICAL_UNIFIED     ! vertical variation of KAPPA (unitless),
                            !  e.g. normalized buoyancy frequency dependent 
                            !  profiles at the tracer grid points
                            !  ( = N^2 / N_ref^2 ) OR a time-independent
                            !  user-specified function


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

      !dir$ attributes offload:mic :: SLA_SAVE_UNIFIED
      real (r8), dimension(:,:,:,:,:), allocatable , public :: &
         SLA_SAVE_UNIFIED             ! isopycnal slopes

      integer (int_kind), parameter ::   &
         kappa_type_const         = 1,   &
         kappa_type_depth         = 2,   &
         kappa_type_vmhs          = 3,   &
         kappa_type_hdgr          = 4,   &
         kappa_type_dradius       = 5,   &
         kappa_type_bfreq         = 6,   &
         kappa_type_bfreq_vmhs    = 7,   &
         kappa_type_bfreq_hdgr    = 8,   &
         kappa_type_bfreq_dradius = 9,   &
         kappa_type_eg            = 10,  &
         slope_control_tanh   = 1,       &
         slope_control_notanh = 2,       &
         slope_control_clip   = 3,       &
         slope_control_Gerd   = 4,       &
         kappa_freq_never           = 1, &
         kappa_freq_every_time_step = 2, &
         kappa_freq_once_a_day      = 3

      !dir$ attributes offload:mic :: compute_kappa_unified
      logical (log_kind), dimension(:), allocatable, public :: &
         compute_kappa_unified        ! compute spatially varying coefficients
                              !  this time step?

      !dir$ attributes offload:mic :: KAPPA_ISOP_UNIFIED
      !dir$ attributes offload:mic :: KAPPA_THIC_UNIFIED
      real (r8), dimension(:,:,:,:,:), allocatable, public :: &
       KAPPA_ISOP_UNIFIED, &      ! 3D isopycnal diffusion coefficient
                                  !  for top and bottom half of a grid cell
       KAPPA_THIC_UNIFIED         ! 3D thickness diffusion coefficient
                                  !  for top and bottom half of a grid cell

      !dir$ attributes offload:mic :: BUOY_FREQ_SQ_UNIFIED
      !dir$ attributes offload:mic :: SIGMA_TOPO_MASK_UNIFIED 
      real (r8), dimension(:,:,:,:), allocatable, public :: &
         BUOY_FREQ_SQ_UNIFIED,    & ! N^2 defined at level interfaces
         SIGMA_TOPO_MASK_UNIFIED    ! bottom topography mask used with kappa_type_eg


      !dir$ attributes offload:mic :: FZTOP_UNIFIED
      real (r8), dimension(:,:,:,:), allocatable,public :: &
         FZTOP_UNIFIED              ! vertical flux

      !dir$ attributes offload:mic :: SF_SLX_UNIFIED
      !dir$ attributes offload:mic :: SF_SLY_UNIFIED
      real (r8), dimension(:,:,:,:,:,:), allocatable ,public :: &
         SF_SLX_UNIFIED, SF_SLY_UNIFIED       ! components of the merged streamfunction

   !dir$ attributes offload : mic :: VDC_UNIFIED
   real (r8), dimension(:,:,:,:,:), allocatable, public, target :: &
      VDC_UNIFIED                 ! tracer diffusivity - public to allow
                          ! possible modification by Gent-McWilliams
                          ! horizontal mixing parameterization

   !dir$ attributes offload:mic :: VDC_GM_UNIFIED
   real (r8), dimension(:,:,:,:), allocatable, public, target :: &
      VDC_GM_UNIFIED              ! Gent-McWilliams contribution to VDC



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
 use grid, only: KMT, dz, partial_bottom_cells, DZT, dzr, dzwr,KMTE,&
                   KMTN,zt,dzw,zw,DXT,DYT,HTN,HTE,TAREA_R,HUS,HUW
 use vmix_kpp, only:zgrid,KPP_HBLT,HMXL
 use prognostic
 use vertical_mix
 use omp_lib
 use hmix_gm
 use hmix_gm_submeso_share, only: HXY,HYX
 use blocks
  
   integer :: iblock

  type (block) ::        &
      this_block           ! block information for current block


   if( .not. allocated (HXY_UNIFIED)) then

   print *,"Intializing unified grids"

  !------------------------done list--------------------------------!

   allocate (SLX_UNIFIED   (nx_block_unified,ny_block_unified,2,2,km,nblocks_clinic),  &
             SLY_UNIFIED   (nx_block_unified,ny_block_unified,2,2,km,nblocks_clinic))

   allocate (TX_UNIFIED(nx_block_unified,ny_block_unified,km,nt,nblocks_clinic),  &
             TY_UNIFIED(nx_block_unified,ny_block_unified,km,nt,nblocks_clinic),  &
             TZ_UNIFIED(nx_block_unified,ny_block_unified,km,nt,nblocks_clinic))

   allocate (RX_UNIFIED(nx_block_unified,ny_block_unified,2,km,nblocks_clinic),  &
             RY_UNIFIED(nx_block_unified,ny_block_unified,2,km,nblocks_clinic))

   allocate (RZ_SAVE_UNIFIED(nx_block_unified,ny_block_unified,km,nblocks_clinic))

   allocate (SF_SUBM_X_UNIFIED(nx_block_unified,ny_block_unified,2,2,km,nblocks_clinic))
   allocate (SF_SUBM_Y_UNIFIED(nx_block_unified,ny_block_unified,2,2,km,nblocks_clinic))

   allocate (HOR_DIFF_UNIFIED(nx_block_unified,ny_block_unified,2,km,nblocks_clinic))

   allocate(WTOP_ISOP_UNIFIED(nx_block_unified,ny_block_unified,nblocks_clinic), &
              WBOT_ISOP_UNIFIED(nx_block_unified,ny_block_unified,nblocks_clinic), &
                    UIT_UNIFIED(nx_block_unified,ny_block_unified,nblocks_clinic), &
                    VIT_UNIFIED(nx_block_unified,ny_block_unified,nblocks_clinic))

   allocate (HMXL_UNIFIED(nx_block_unified,ny_block_unified,1))
   allocate (KPP_HBLT_UNIFIED(nx_block_unified,ny_block_unified,1))

  allocate (HXY_UNIFIED (nx_block_unified,ny_block_unified,1),    &
            HYX_UNIFIED (nx_block_unified,ny_block_unified,1))


   !------------------------------------------------------------------!

   endif

   SLX_UNIFIED      = c0
   SLY_UNIFIED      = c0
   TX_UNIFIED       = c0
   TY_UNIFIED       = c0
   TZ_UNIFIED       = c0
   RX_UNIFIED       = c0
   RY_UNIFIED       = c0
   RZ_SAVE_UNIFIED  = c0

    do iblock = 1,nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       call merger( HMXL(:,:,iblock), HMXL_UNIFIED(:,:,1), iblock, this_block)
       call merger( KPP_HBLT(:,:,iblock), KPP_HBLT_UNIFIED(:,:,1), iblock, this_block)
       call merger( HXY(:,:,iblock), HXY_UNIFIED(:,:,1) , iblock, this_block)
       call merger( HYX(:,:,iblock), HYX_UNIFIED(:,:,1) , iblock, this_block)   
    enddo

   !if(my_task == master_task) print *,"HYX is",HYX_UNIFIED(45,45,1),HYX(45,45,1)

!! initialization for state 

   tmin =  -2.0_r8  ! limited   on the low  end
   tmax = 999.0_r8  ! unlimited on the high end
   smin =   0.0_r8  ! limited   on the low  end
   smax = 0.999_r8  ! unlimited on the high end

!! initialization for grid 

   if( .not. allocated(DZT_unified) ) then

   allocate (DZT_unified(nx_block,ny_block,0:km+1,max_blocks_clinic))

   endif

   pressz_unified = pressz
   KMT_UNIFIED = KMT
   KMTE_UNIFIED = KMTE
   KMTN_UNIFIED = KMTN
   DZT_UNIFIED = DZT

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
   HUS_UNIFIED = HUS
   HUW_UNIFIED = HUW

   !if(my_task == master_task) print *,"KMT is",KMT_UNIFIED(45,45,1),KMT(45,45,1)

!! initialization for mix_submeso

   if( .not. allocated(TIME_SCALE_UNIFIED)) then

   allocate (TIME_SCALE_UNIFIED(nx_block,ny_block,nblocks_clinic))

   allocate (FZTOP_SUBM_UNIFIED(nx_block,ny_block,nt,nblocks_clinic))

   endif

   TIME_SCALE_UNIFIED = TIME_SCALE

      efficiency_factor_unified = efficiency_factor
      hor_length_scale_unified = hor_length_scale
      sqrt_grav_unified = sqrt(grav)
      max_hor_grid_scale_unified = max_hor_grid_scale

!! Variables related to vmix_kpp

  if( .not. allocated (zgrid_unified )) then

  allocate  (zgrid_unified(0:km+1)) 

  endif

  zgrid_unified = zgrid

!! variables for allocatng hmix_gm variables

  if(.not. allocated(SLA_SAVE_UNIFIED )) then

  allocate (SLA_SAVE_UNIFIED(nx_block,ny_block,2,km,nblocks_clinic))
  allocate (RB_UNIFIED(nx_block,ny_block,nblocks_clinic))
  allocate (compute_kappa_unified(nblocks_clinic))
 
  allocate(KAPPA_ISOP_UNIFIED(nx_block,ny_block,2,km,nblocks_clinic),  &
           KAPPA_THIC_UNIFIED(nx_block,ny_block,2,km,nblocks_clinic))

  allocate (KAPPA_LATERAL_UNIFIED(nx_block,ny_block,nblocks_clinic),  &
            KAPPA_VERTICAL_UNIFIED(nx_block,ny_block,km,nblocks_clinic))

  allocate (BUOY_FREQ_SQ_UNIFIED(nx_block,ny_block,km,nblocks_clinic))

  allocate (FZTOP_UNIFIED(nx_block,ny_block,nt,nblocks_clinic))

  allocate (SF_SLX_UNIFIED(nx_block,ny_block,2,2,km,nblocks_clinic),  &
            SF_SLY_UNIFIED(nx_block,ny_block,2,2,km,nblocks_clinic))

  allocate (HYXW_unified(nx_block,ny_block,nblocks_clinic),    &
             HXYS_unified(nx_block,ny_block,nblocks_clinic))

  allocate (VDC_UNIFIED(nx_block,ny_block,0:km+1,2,nblocks_clinic))

  allocate (VDC_GM_UNIFIED(nx_block,ny_block,km,nblocks_clinic))
  
  allocate (RBR_UNIFIED (nx_block,ny_block,nblocks_clinic))

  endif

         KAPPA_ISOP_UNIFIED = KAPPA_ISOP
         KAPPA_THIC_UNIFIED = KAPPA_THIC
         KAPPA_LATERAL_UNIFIED = KAPPA_LATERAL
         KAPPA_VERTICAL_UNIFIED = KAPPA_VERTICAL
         HYXW_UNIFIED = HYXW
         HXYS_UNIFIED = HXYS 
         RB_UNIFIED = RB 
         RBR_UNIFIED = RBR

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

      if (k == 1) then
         !start_time = omp_get_wtime()

         call hdifft_gm_unified(1, HDTK_BUF(:,:,:,1,bid), TMIX, UMIX, VMIX, tavg_HDIFE_TRACER, &
                         tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

         !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)num_threads(59) 
         do kk=2,km
         call hdifft_gm_unified(kk , HDTK_BUF(:,:,:,kk,bid) , TMIX, UMIX,VMIX,tavg_HDIFE_TRACER, &
                                 tavg_HDIFN_TRACER,tavg_HDIFB_TRACER,this_block)
         enddo

         !end_time = omp_get_wtime()

         !print *,"time at hdifft_gm combined is ",end_time - start_time 
                    
      endif

     !if(nsteps_total == 6 .and. k == 1 .and. my_task == master_task) then

      !print *,"changed cont HDTK is ",HDTK_BUF(3,7,1,1)

      !endif

      !start_time = omp_get_wtime()  
      HDTK = HDTK_BUF(:,:,:,k,bid)
      !end_time = omp_get_wtime()

      !print *,"time at hdifft_gm is ",end_time - start_time
 
        if (k == 1) then
         !start_time = omp_get_wtime()
         call submeso_sf_unified(TMIX, this_block)
         !end_time = omp_get_wtime()
         !print *,"time at submeso_sf is ",end_time - start_time
        endif


        !if(my_task == master_task .and. nsteps_total == 3 .and. k == 45)then

           !print *,"KMT outside is",KMT_UNIFIED(45,45,bid)

           !print *,"HYX(45,45,bid)",HYX_UNIFIED(45,45,bid)

         !endif


        if(k==1) then
         !start_time = omp_get_wtime()
        !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)num_threads(60) 
        do kk=1,km
         call submeso_flux_unified(kk, TDTK(:,:,:,kk,bid), TMIX, tavg_HDIFE_TRACER, &
                              tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)
        enddo
        !end_time = omp_get_wtime()
        !print *,"time at submeso_flux is ",end_time - start_time
        endif
        HDTK=HDTK+TDTK(:,:,:,k,bid)
   
     !if( nsteps_total == 6 .and. k == 1 .and. my_task == master_task) then

      !print *,"changed cont TDTK is ",TDTK(3,7,1,1)

      !endif



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

      bid = this_block%local_id


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

            !if(my_task == master_task .and. nsteps_total == 3 .and. k == 45 .and. i == 45 - 1 .and. j == 45 .and. n == 1)then

                   !print *,"changed in flux"
                   !print *,"SF_SUBM_X(i+1,j,iwest,kbt,k,bid) contribution is",SF_SUBM_X_UNIFIED(i+1,j,iwest,kbt,k,bid)
                   !print *,"SF_SUBM_X(i+1,j,iwest,ktp,k,bid) contribution is",SF_SUBM_X_UNIFIED(i+1,j,iwest,ktp,k,bid)
                   !print *,"SF_SUBM_X(i  ,j,ieast,kbt,k,bid) contribution is",SF_SUBM_X_UNIFIED(i  ,j,ieast,kbt,k,bid)
                   !print *,"SF_SUBM_X(i  ,j,ieast,ktp,k,bid)  contributionis",SF_SUBM_X_UNIFIED(i  ,j,ieast,ktp,k,bid)
                   !print *,"CX(i,j) contribution is",CX(i,j)
                   !print *,"TZ(i,j,k,n,bid) is",TZ_UNIFIED(i,j,k,n,bid)
                   !print *,"TZ(i,j,kp1,n,bid) ",TZ_UNIFIED(i,j,kp1,n,bid)
                   !print *,"TZ(i+1,j,k,n,bid)",TZ_UNIFIED(i+1,j,k,n,bid)
                   !print *,"TZ(i+1,j,kp1,n,bid)",TZ_UNIFIED(i,j,kp1,n,bid)

            !endif


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

                  fz = -KMASK(i,j) * p25    &
                      * (WORK1(i,j) + WORK2(i,j))

             !if(my_task == master_task .and. nsteps_total == 3 .and. k == 45 .and. i == 45 .and. j == 45 .and. n == 1)then

                   !print *,"changed before in flux"
                   !print *,"FX(i,j,n) contribution is",FX(i,j,n)
                   !print *,"FX(i-1,j,n) contribution is",FX(i-1,j,n)
                   !print *,"FY(i,j,n) contribution is",FY(i,j,n)                   
                   !print *,"FY(i,j-1,n) contribution is",FY(i,j-1,n)

            !endif



                  GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                               + FY(i,j,n) - FY(i,j-1,n)  &
                        + fzprev - fz )*dzr_unified(k)*TAREA_R_UNIFIED(i,j,bid)

                  !FZTOP_SUBM(i,j,n,bid) = fz

               else  !k == km

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


      CX = merge(HYX_UNIFIED(:,:,bid)*p25, c0, (k <= KMT_UNIFIED (:,:,bid))   &
                                 .and. (k <= KMTE_UNIFIED(:,:,bid)))
      CY = merge(HXY_UNIFIED(:,:,bid)*p25, c0, (k <= KMT_UNIFIED (:,:,bid))   &
                                 .and. (k <= KMTN_UNIFIED(:,:,bid)))


     if ( k == 1 ) then

        if ( diag_gm_bolus ) then
          UIB = c0
          VIB = c0
          UIT_UNIFIED(:,:,bid) = c0
          VIT_UNIFIED(:,:,bid) = c0
          WBOT_ISOP_UNIFIED(:,:,bid) = c0
        endif

        HOR_DIFF_UNIFIED(:,:,ktp,k,bid) = ah_bkg_srfbl

        BL_DEPTH_UNIFIED(:,:,bid) = zw_unified(k)
        if ( vmix_itype == vmix_type_kpp )  &
                    BL_DEPTH_UNIFIED(:,:,bid) = KPP_HBLT_UNIFIED(:,:,bid)

        if ( transition_layer_on ) then

          if ( vmix_itype == vmix_type_kpp ) then

            !start_time = omp_get_wtime()

            call smooth_hblt_unified ( .false., .true., bid,  &
                               SMOOTH_OUT=TLT_UNIFIED%DIABATIC_DEPTH(:,:,bid) )

            !end_time = omp_get_wtime()

            !print *,"Smmoth takes ",end_time - start_time  

          else
            TLT_UNIFIED%DIABATIC_DEPTH(:,:,bid) = zw_unified(k)
          endif

          !!$OMP PARALLEL DO &
          !!$OMP DEFAULT(SHARED)PRIVATE(kid,i,j,kk_sub,kk)NUM_THREADS(60)COLLAPSE(3)
          do kk=1,km
            do kk_sub = ktp,kbt
                  do j=1,ny_block
                     do i=1,nx_block

                        kid = kk + kk_sub - 2

                        SLA_SAVE_UNIFIED(i,j,kk_sub,kk,bid) = dzw_unified(kid)*sqrt(p5*( &
                        (SLX_UNIFIED(i,j,1,kk_sub,kk,bid)**2                     &
                        + SLX_UNIFIED(i,j,2,kk_sub,kk,bid)**2)/DXT_UNIFIED(i,j,bid)**2   &
                        + (SLY_UNIFIED(i,j,1,kk_sub,kk,bid)**2                   &
                        + SLY_UNIFIED(i,j,2,kk_sub,kk,bid)**2)/DYT_UNIFIED(i,j,bid)**2)) &
                        + eps

                     enddo
                  enddo
            enddo
          enddo
          !!$OMP END PARALLEL DO
          call transition_layer_unified ( this_block ) 
 
        endif

      !if(my_task == master_task)then


         !print *,"Changed before"
         !print*,"KAPPA_THIC_UNIFIED(i,j,ktp,k,bid)",KAPPA_THIC_UNIFIED(3,7,ktp,4,bid),nsteps_total

      !endif
        

     if ( ( kappa_isop_type == kappa_type_vmhs           .or.    &
            kappa_thic_type == kappa_type_vmhs           .or.    &
            kappa_isop_type == kappa_type_hdgr           .or.    &
            kappa_thic_type == kappa_type_hdgr           .or.    &
            kappa_isop_type == kappa_type_dradius        .or.    &
            kappa_thic_type == kappa_type_dradius        .or.    &
            kappa_isop_type == kappa_type_bfreq          .or.    &
            kappa_thic_type == kappa_type_bfreq          .or.    &
            kappa_isop_type == kappa_type_bfreq_vmhs     .or.    &
            kappa_thic_type == kappa_type_bfreq_vmhs     .or.    &
            kappa_isop_type == kappa_type_bfreq_hdgr     .or.    &
            kappa_thic_type == kappa_type_bfreq_hdgr     .or.    &
            kappa_isop_type == kappa_type_bfreq_dradius  .or.    &
            kappa_thic_type == kappa_type_bfreq_dradius  .or.    &
            kappa_isop_type == kappa_type_eg             .or.    &
            kappa_thic_type == kappa_type_eg )            .and.  &
        ( ( kappa_freq == kappa_freq_every_time_step )           &
     .or. ( kappa_freq == kappa_freq_once_a_day .and. eod_last ) &
     .or. ( nsteps_total == 1 ) ) )  compute_kappa_unified(bid) = .true.

     if ( compute_kappa_unified(bid) ) then

               if ( kappa_isop_type == kappa_type_bfreq          .or.  &
               kappa_thic_type == kappa_type_bfreq          .or.  &
               kappa_isop_type == kappa_type_bfreq_vmhs     .or.  &
               kappa_thic_type == kappa_type_bfreq_vmhs     .or.  &
               kappa_isop_type == kappa_type_bfreq_hdgr     .or.  &
               kappa_thic_type == kappa_type_bfreq_hdgr     .or.  &
               kappa_isop_type == kappa_type_bfreq_dradius  .or.  &
               kappa_thic_type == kappa_type_bfreq_dradius )      &
            call buoyancy_frequency_dependent_profile_unified (TMIX, this_block)

            compute_kappa_unified(bid) = .false.

     endif

        !start_time = omp_get_wtime()

          !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk_sub,kk,j,i)NUM_THREADS(60)
          do kk_sub=ktp,kbt
            do kk=1,km
               do j=1,ny_block
                   do i=1,nx_block
                       KAPPA_ISOP_UNIFIED(i,j,kk_sub,kk,bid) =  KAPPA_LATERAL_UNIFIED(i,j,bid) &
                                                     *  KAPPA_VERTICAL_UNIFIED(i,j,kk,bid)
                   enddo
               enddo
            enddo
          enddo

          !start_time = omp_get_wtime()
          !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk_sub,kk,j,i)NUM_THREADS(60)collapse(3)schedule(dynamic,4)
           do kk_sub=ktp,kbt
            do kk=1,km
             do j=1,ny_block
              do i=1,nx_block
                 KAPPA_THIC_UNIFIED(i,j,kk_sub,kk,bid) =  ah_bolus  &
                                          * KAPPA_VERTICAL_UNIFIED(i,j,kk,bid)
              enddo
             enddo
            enddo
          enddo


      !if(my_task == master_task)then


         !print *,"Changed after"
         !print*,"KAPPA_THIC_UNIFIED(i,j,ktp,k,bid)",KAPPA_THIC_UNIFIED(3,7,ktp,4,bid),nsteps_total
         !print *,"-----------------------------------"

      !endif


    do kk=1,km

          kp1 = min(kk+1,km)
          reference_depth(ktp) = zt_unified(kp1)
          reference_depth(kbt) = zw_unified(kp1)
          if ( kk == km )  reference_depth(ktp) = zw_unified(kp1)

            do kk_sub = ktp,kbt 

             kid = kk + kk_sub - 2

             !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i,dzw,dz_bottom,zt)NUM_THREADS(60)
             do j=1,ny_block
                   do i=1,nx_block

                       if ( transition_layer_on ) then
                          SLA(i,j) = SLA_SAVE_UNIFIED(i,j,kk_sub,kk,bid)
                       else
                          SLA(i,j) = dzw_unified(kid)*sqrt(p5*( &
                                 (SLX_UNIFIED(i,j,1,kk_sub,kk,bid)**2 & 
                                + SLX_UNIFIED(i,j,2,kk_sub,kk,bid)**2)/DXT_UNIFIED(i,j,bid)**2 &
                                + (SLY_UNIFIED(i,j,1,kk_sub,kk,bid)**2 &
                                + SLY_UNIFIED(i,j,2,kk_sub,kk,bid)**2)/DYT_UNIFIED(i,j,bid)**2)) &
                                + eps
                        endif   

                        TAPER1(i,j) = c1 
                        if ( .not. transition_layer_on ) then

                        if ( kk == 1 ) then
                        dz_bottom = c0
                        else
                        dz_bottom = zt_unified(kk-1)
                        endif           


                        if (slope_control == slope_control_tanh) then

                        WORK1(i,j) = min(c1,zt_unified(kk)*RBR_UNIFIED(i,j,bid)/SLA(i,j))
                        TAPER1(i,j) = p5*(c1+sin(pi*(WORK1(i,j)-p5)))

!     use the Rossby deformation radius tapering
!     only within the boundary layer

                        TAPER1(i,j) = merge(TAPER1(i,j), c1,  &
                               dz_bottom <= BL_DEPTH_UNIFIED(i,j,bid))

                         else

!     sine function is replaced by
!     function = 4.*x*(1.-abs(x)) for |x|<0.5

                         WORK1(i,j) = min(c1,zt_unified(kk)*RBR_UNIFIED(i,j,bid)/SLA(i,j))
                         TAPER1(i,j) = (p5+c2*(WORK1(i,j)-p5)*(c1-abs(WORK1(i,j)-p5)))

                         TAPER1(i,j) = merge(TAPER1(i,j), c1,  &
                                  dz_bottom <= BL_DEPTH_UNIFIED(i,j,bid))

                         endif

                        endif ! not transition_layer
                        TAPER2(i,j) = c1
                        TAPER3(i,j) = c1


                         select case (slope_control)
                         case (slope_control_tanh)

!     method by Danabasoglu & Mcwilliams (1995)

                         TAPER2(i,j) = merge(p5*  &
                          (c1-tanh(c10*SLA(i,j)/slm_r-c4)), c0, SLA(i,j) < slm_r)

                         if ( diff_tapering ) then
                          TAPER3(i,j) = merge(p5*  &
                          (c1-tanh(c10*SLA(i,j)/slm_b-c4)), c0, SLA(i,j) < slm_b)
                         else
                          TAPER3(i,j) = TAPER2(i,j)
                         endif

                         case (slope_control_notanh)

!     similar to DM95 except replacing tanh by
!     function = x*(1.-0.25*abs(x)) for |x|<2
!              = sign(x)            for |x|>2
!     (faster than DM95)


                         if (SLA(i,j) > 0.2_r8*slm_r .and. &
                             SLA(i,j) < 0.6_r8*slm_r) then
                             TAPER2(i,j) = &
                             p5*(c1-(2.5_r8*SLA(i,j)/slm_r-c1)*  &
                             (c4-abs(c10*SLA(i,j)/slm_r-c4)))
                         else if (SLA(i,j) >= 0.6_r8*slm_r) then
                             TAPER2(i,j) = c0
                         endif

                         if ( diff_tapering ) then

                           if (SLA(i,j) > 0.2_r8*slm_b .and. &
                              SLA(i,j) < 0.6_r8*slm_b) then
                              TAPER3(i,j) = &
                              p5*(c1-(2.5_r8*SLA(i,j)/slm_b-c1)* &
                              (c4-abs(c10*SLA(i,j)/slm_b-c4)))
                           else if (SLA(i,j) >= 0.6_r8*slm_b) then
                              TAPER3(i,j) = c0
                         endif

                         else
                              TAPER3(i,j) = TAPER2(i,j)
                         endif

                         case (slope_control_clip)

!     slope clipping

                         do n=1,2

                         if (abs(SLX_UNIFIED(i,j,n,kk_sub,kk,bid)  &
                            * dzw_unified(kid) / HUS_UNIFIED(i,j,bid)) > slm_r) then
                            SLX_UNIFIED(i,j,n,kk_sub,kk,bid) =             &
                                        sign(slm_r * HUS_UNIFIED(i,j,bid)  &
                                           * dzwr_unified(kid),            &
                                        SLX_UNIFIED(i,j,n,kk_sub,kk,bid))
                         endif
                         enddo

                         do n=1,2

                         if (abs(SLY_UNIFIED(i,j,n,kk_sub,kk,bid)  &
                          * dzw_unified(kid) / HUW_UNIFIED(i,j,bid)) > slm_r) then
                          SLY_UNIFIED(i,j,n,kk_sub,kk,bid) =             &
                                  sign(slm_r * HUW_UNIFIED(i,j,bid)  &
                                     * dzwr_unified(kid),            &
                                  SLY_UNIFIED(i,j,n,kk_sub,kk,bid))
                          endif
                          enddo

                          case (slope_control_Gerd)

!     method by Gerdes et al (1991)


                          if (SLA(i,j) > slm_r)  &
                             TAPER2(i,j) = (slm_r/SLA(i,j))**2


                          if (diff_tapering) then

                             if (SLA(i,j) > slm_b)  &
                                TAPER3(i,j) = (slm_b/SLA(i,j))**2

                             else
                                TAPER3(i,j) = TAPER2(i,j)
                             endif

                         end select


                      if ( transition_layer_on ) then
                           TAPER2(i,j) = merge(c1, TAPER2(i,j), reference_depth(kk_sub) &
                                                         <= TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid))
                           TAPER3(i,j) = merge(c1, TAPER3(i,j), reference_depth(kk_sub) &
                                                         <= TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid))
                      endif
 
                       if ( transition_layer_on  .and.  use_const_ah_bkg_srfbl ) then

                           HOR_DIFF_UNIFIED(i,j,kk_sub,kk,bid) = ah_bkg_srfbl

                       else if ( transition_layer_on .and.               &
                               ( .not. use_const_ah_bkg_srfbl      .or.  &
                                 kappa_isop_type == kappa_type_eg  .or.  &
                                 kappa_thic_type == kappa_type_eg ) ) then

                           HOR_DIFF_UNIFIED(i,j,kk_sub,kk,bid) = KAPPA_ISOP_UNIFIED(i,j,kk_sub,kk,bid)

                       else

                       if ( .not. ( kk == 1 .and. kk_sub == ktp ) ) then

                          if ( use_const_ah_bkg_srfbl ) then
                              HOR_DIFF_UNIFIED(i,j,kk_sub,kk,bid) = &
                                   merge( ah_bkg_srfbl * (c1 - TAPER1(i,j) * TAPER2(i,j))   &
                                          * KAPPA_VERTICAL_UNIFIED(i,j,kk,bid), &
                                           c0, dz_bottom <= BL_DEPTH_UNIFIED(i,j,bid) )
                          else
                              HOR_DIFF_UNIFIED(i,j,kk_sub,kk,bid) =         &
                              merge( KAPPA_ISOP_UNIFIED(i,j,kk_sub,kk,bid)  &
                               * (c1 - TAPER1(i,j) * TAPER2(i,j)),  &
                              c0, dz_bottom <= BL_DEPTH_UNIFIED(i,j,bid) )
                          endif

                       endif

                       endif

      !if(my_task == master_task .and. nsteps_total == 6 .and. kk == 2+2 .and. i==3 .and. j==7 .and. kk_sub == ktp)then


         !print *,"Original"
         !print *,"KAPPA_THIC_UNIFIED(i,j,ktp,k,bid)",KAPPA_THIC_UNIFIED(i,j,kk_sub,kk,bid)
         !print *,"TAPER1(i,j)",TAPER1(i,j)
         !print *,"TAPER3(i,j)",TAPER3(i,j)

      !endif


                      KAPPA_ISOP_UNIFIED(i,j,kk_sub,kk,bid) =  &
                      TAPER1(i,j) * TAPER2(i,j) * KAPPA_ISOP_UNIFIED(i,j,kk_sub,kk,bid)

                      KAPPA_THIC_UNIFIED(i,j,kk_sub,kk,bid) =  &
                      TAPER1(i,j) * TAPER3(i,j) * KAPPA_THIC_UNIFIED(i,j,kk_sub,kk,bid)


                   enddo !j Loop
             enddo  !i loop

            enddo !kk_sub loop

           do j=1,ny_block
                   do i=1,nx_block

                       if (kk == KMT_UNIFIED(i,j,bid)) then
                          KAPPA_ISOP_UNIFIED(i,j,kbt,kk,bid) = c0
                          KAPPA_THIC_UNIFIED(i,j,kbt,kk,bid) = c0
                       endif

                   enddo
              enddo


    enddo !kk loop


    KAPPA_ISOP_UNIFIED(:,:,ktp,1,bid) = c0
    KAPPA_THIC_UNIFIED(:,:,ktp,1,bid) = c0

    FZTOP_UNIFIED(:,:,:,bid) = c0        ! zero flux B.C. at the surface

      if ( transition_layer_on ) then

          !start_time = omp_get_wtime()


          call merged_streamfunction_unified ( this_block )


          !end_time = omp_get_wtime()

          !print *,"Time taken at function1 is ",end_time - start_time

          !start_time = omp_get_wtime()


          call apply_vertical_profile_to_isop_hor_diff_unified ( this_block )

     endif 

  
    endif !k == 1 

    KMASK = merge(c1, c0, k < KMT_UNIFIED(:,:,bid))

     if ( k < km ) then

      do j=1,ny_block
         do i=1,nx_block

         WORK1(i,j) = dzw_unified(k)*KMASK(i,j)*TAREA_R_UNIFIED(i,j,bid)*      &
                 (dz_unified(k)*p25*KAPPA_ISOP_UNIFIED(i,j,kbt,k,  bid)*      &
               (HYX_UNIFIED (i,j,bid)*SLX_UNIFIED(i,j,ieast, kbt,k,  bid)**2   &
              + HYXW_UNIFIED(i,j,bid)*SLX_UNIFIED(i,j,iwest, kbt,k,  bid)**2   &
              + HXY_UNIFIED (i,j,bid)*SLY_UNIFIED(i,j,jnorth,kbt,k,  bid)**2   &
              + HXYS_UNIFIED(i,j,bid)*SLY_UNIFIED(i,j,jsouth,kbt,k,  bid)**2)  &
                 +dz_unified(k+1)*p25*KAPPA_ISOP_UNIFIED(i,j,ktp,k+1,bid)*     &
               (HYX_UNIFIED (i,j,bid)*SLX_UNIFIED(i,j,ieast, ktp,k+1,bid)**2   &
              + HYXW_UNIFIED(i,j,bid)*SLX_UNIFIED(i,j,iwest, ktp,k+1,bid)**2   &
              + HXY_UNIFIED (i,j,bid)*SLY_UNIFIED(i,j,jnorth,ktp,k+1,bid)**2   &
              + HXYS_UNIFIED(i,j,bid)*SLY_UNIFIED(i,j,jsouth,ktp,k+1,bid)**2))


              do n=1,size(VDC_UNIFIED,DIM=4)
                 VDC_GM_UNIFIED(i,j,k,bid) = WORK1(i,j)
                 VDC_UNIFIED(i,j,k,n,bid) = VDC_UNIFIED(i,j,k,n,bid) + WORK1(i,j)
              end do

          enddo
      enddo


      end if

     if ( ah_bkg_bottom /= c0 ) then

         do j=1,ny_block
            do i=1,nx_block

             if( k == KMT_UNIFIED(i,j,bid)) then
                  HOR_DIFF_UNIFIED(i,j,kbt,k,bid) = ah_bkg_bottom
             endif

            enddo
         enddo

      endif

     n = 1

     do j=1,ny_block
        do i=1,nx_block-1



         !if(my_task == master_task .and. nsteps_total == 1 .and. k == 1 .and. i == 3 .and. j == 5 .and. n == 1)then


             !print *,"Changed"
             !print *,"HOR_UNIFIED(i,  j,ktp,k,bid)",HOR_DIFF_UNIFIED(i,j,ktp,k,bid)
             !print *,"HOR_UNIFIED(i,j,kbt,k,bid)",HOR_DIFF_UNIFIED(i,j,kbt,k,bid)
             !print *,"HOR_UNIFIED(i+1,j,ktp,k,bid)",HOR_DIFF_UNIFIED(i+1,j,ktp,k,bid)
             !print *,"HOR_UNIFIED(i,j,ktp,k,bid)",HOR_DIFF_UNIFIED(i+1,j,kbt,k,bid)

         !endif

          WORK3(i,j) = KAPPA_ISOP_UNIFIED(i,  j,ktp,k,bid)  &
                     + HOR_DIFF_UNIFIED  (i,  j,ktp,k,bid)  &
                     + KAPPA_ISOP_UNIFIED(i,  j,kbt,k,bid)  &
                     + HOR_DIFF_UNIFIED  (i,  j,kbt,k,bid)  &
                     + KAPPA_ISOP_UNIFIED(i+1,j,ktp,k,bid)  &
                     + HOR_DIFF_UNIFIED  (i+1,j,ktp,k,bid)  &
                     + KAPPA_ISOP_UNIFIED(i+1,j,kbt,k,bid)  &
                     + HOR_DIFF_UNIFIED  (i+1,j,kbt,k,bid)
        enddo
      enddo

      do j=1,ny_block-1
        do i=1,nx_block
          WORK4(i,j) = KAPPA_ISOP_UNIFIED(i,j,  ktp,k,bid)  &
                     + HOR_DIFF_UNIFIED  (i,j,  ktp,k,bid)  &
                     + KAPPA_ISOP_UNIFIED(i,j,  kbt,k,bid)  &
                     + HOR_DIFF_UNIFIED  (i,j,  kbt,k,bid)  &
                     + KAPPA_ISOP_UNIFIED(i,j+1,ktp,k,bid)  &
                     + HOR_DIFF_UNIFIED  (i,j+1,ktp,k,bid)  &
                     + KAPPA_ISOP_UNIFIED(i,j+1,kbt,k,bid)  &
                     + HOR_DIFF_UNIFIED  (i,j+1,kbt,k,bid)
        enddo
      enddo


      kp1 = k + 1
      if ( k == km )  kp1 = k

      if ( k < km ) then
        dz_bottom = dz_unified(kp1)
        factor    = c1
      else
        dz_bottom = c0
        factor    = c0
      endif


     do n = 1,nt

!-----------------------------------------------------------------------
!
!     calculate horizontal fluxes thru vertical faces of T-cell
!     FX = dz*HYX*Ax(Az(KAPPA))*Dx(T) : flux in x-direction
!     FY = dz*HXY*Ay(Az(KAPPA))*Dy(T) : flux in y-direction
!
!-----------------------------------------------------------------------

        FX(:,:,n) = dz_unified(k) * CX * TX_UNIFIED(:,:,k,n,bid) * WORK3
        FY(:,:,n) = dz_unified(k) * CY * TY_UNIFIED(:,:,k,n,bid) * WORK4

           !if(my_task == master_task .and. nsteps_total == 1 .and. k == 1 .and. i == 3 .and. j == 5 .and. n == 1)then


             !print *,"Changed"
             !print *,"FX(i,j,n)",FX(i,j,n)
             !print *,"dz_unified",dz_unified(k)
             !print *,"CX",CX(i,j) 
             !print *,"TX_UNIFIED",TX_UNIFIED(i,j,k,n,bid)
             !print *,"WORK3",WORK3(i,j)
             !print *,"WORK4",WORK4(i,j)

             !endif


      end do


      if ( .not. cancellation_occurs ) then

        do j=1,ny_block
          do i=1,nx_block-1

           !if(my_task == master_task .and. nsteps_total == 6 .and. k == 1 .and. i == 3 .and. j == 7 )then


             !print *,"Changed"
             !print *,"dz_unified(k)",dz_unified(k)
             !print *,"KAPPA_ISOP_UNIFIED(i,j,ktp,k,bid)",KAPPA_ISOP_UNIFIED(i,j,ktp,k,bid)
             !print *,"SLX_UNIFIED(i,j,ieast,ktp,k,bid)",SLX_UNIFIED(i,j,ieast,ktp,k,bid)             
             !print *,"KAPPA_ISOP_UNIFIED(i,j,kbt,k,bid)",KAPPA_ISOP_UNIFIED(i,j,kbt,k,bid)
             !print *,"SLX_UNIFIED(i,j,ieast,kbt,k,bid)",SLX_UNIFIED(i,j,ieast,kbt,k,bid)          
             !print *,"SF_SLX_UNIFIED(i,j,ieast,ktp,k,bid)",SF_SLX_UNIFIED(i,j,ieast,ktp,k,bid)
             
             !endif


            WORK1(i,j) = KAPPA_ISOP_UNIFIED(i,j,ktp,k,bid)                     &
                         * SLX_UNIFIED(i,j,ieast,ktp,k,bid) * dz_unified(k)            &
                         - SF_SLX_UNIFIED(i,j,ieast,ktp,k,bid)
            WORK2(i,j) = KAPPA_ISOP_UNIFIED(i,j,kbt,k,bid)                     &
                         * SLX_UNIFIED(i,j,ieast,kbt,k,bid) * dz_unified(k)            &
                         - SF_SLX_UNIFIED(i,j,ieast,kbt,k,bid)
            WORK3(i,j) = KAPPA_ISOP_UNIFIED(i+1,j,ktp,k,bid)                   &
                         * SLX_UNIFIED(i+1,j,iwest,ktp,k,bid) * dz_unified(k)          &
                         - SF_SLX_UNIFIED(i+1,j,iwest,ktp,k,bid)
            WORK4(i,j) = KAPPA_ISOP_UNIFIED(i+1,j,kbt,k,bid)                   &
                         * SLX_UNIFIED(i+1,j,iwest,kbt,k,bid) * dz_unified(k)          &
                         - SF_SLX_UNIFIED(i+1,j,iwest,kbt,k,bid)
          enddo
        enddo

        do n = 1,nt
          do j=1,ny_block
            do i=1,nx_block-1


            !if(my_task == master_task .and. nsteps_total == 6 .and. k == 1 .and. i == 3 .and. j == 7 .and. n == 1)then


             !print *,"Changed"
             !print *,"FX(i,j,n)",FX(i,j,n)
             !print *,"CX(i,j)",CX(i,j)
             !print *,"WORK1(i,j)",WORK1(i,j)
             !print *,"WORK2(i,j)",WORK2(i,j)
             !print *,"WORK3(i,j)",WORK3(i,j)
             !print *,"WORK4(i,j)",WORK4(i,j)
             !print *,"TZ(i,j)",TZ_UNIFIED(i,j,k,n,bid)
             !print *,"TZ(i,j,kp1)",TZ_UNIFIED(i,j,kp1,n,bid)
             !print *,"TZ(i+1,j)",TZ_UNIFIED(i+1,j,k,n,bid)
             !print *,"TZ(i+1,j,kp1)",TZ_UNIFIED(i+1,j,kp1,n,bid)


             !endif

              FX(i,j,n) = FX(i,j,n) - CX(i,j)                          &
               * ( WORK1(i,j) * TZ_UNIFIED(i,j,k,n,bid)                        &
                   + WORK2(i,j) * TZ_UNIFIED(i,j,kp1,n,bid)                    &
                   + WORK3(i,j) * TZ_UNIFIED(i+1,j,k,n,bid)                    &
                   + WORK4(i,j) * TZ_UNIFIED(i+1,j,kp1,n,bid) )
            enddo
          enddo
        end do

        do j=1,ny_block-1
          do i=1,nx_block
            WORK1(i,j) = KAPPA_ISOP_UNIFIED(i,j,ktp,k,bid)                     &
                         * SLY_UNIFIED(i,j,jnorth,ktp,k,bid) * dz_unified(k)   &
                         - SF_SLY_UNIFIED(i,j,jnorth,ktp,k,bid)
            WORK2(i,j) = KAPPA_ISOP_UNIFIED(i,j,kbt,k,bid)                     &
                         * SLY_UNIFIED(i,j,jnorth,kbt,k,bid) * dz_unified(k)           &
                         - SF_SLY_UNIFIED(i,j,jnorth,kbt,k,bid)
            WORK3(i,j) = KAPPA_ISOP_UNIFIED(i,j+1,ktp,k,bid)                   &
                         * SLY_UNIFIED(i,j+1,jsouth,ktp,k,bid) * dz_unified(k)         &
                         - SF_SLY_UNIFIED(i,j+1,jsouth,ktp,k,bid)
            WORK4(i,j) = KAPPA_ISOP_UNIFIED(i,j+1,kbt,k,bid)                   &
                         * SLY_UNIFIED(i,j+1,jsouth,kbt,k,bid) * dz_unified(k)         &
                         - SF_SLY_UNIFIED(i,j+1,jsouth,kbt,k,bid)
          enddo
        enddo

       do n = 1,nt

          do j=1,ny_block-1
            do i=1,nx_block


             !if(my_task == master_task .and. nsteps_total == 6 .and. k == 1 .and. i == 3 .and. j == 7 .and. n == 1)then

                   !print *,"changed before"
                   !print *,"FY(i,j,n) contribution is",FY(i,j,n)
                   !print *,"CY(i,j,n) contribution is",CY(i,j)
                   !print *,"WORK1 contribution is",WORK1(i,j)
                   !print *,"WORK2 contribution is",WORK2(i,j)
                   !print *,"WORK3 contribution is",WORK3(i,j)
                   !print *,"WORK4 contribution is",WORK4(i,j)

                   !print *,"TZ contribution is",TZ_UNIFIED(i,j,k,n,bid)
                   !print *,"TZkp1 contribution is",TZ_UNIFIED(i,j,kp1,n,bid)
                   !print *,"TZjp1 contribution is",TZ_UNIFIED(i,j+1,k,n,bid)
                   !print *,"TZJp1kp1 contribution is",TZ_UNIFIED(i,j+1,kp1,n,bid)

           !endif


              FY(i,j,n) = FY(i,j,n) - CY(i,j)                          &
               * ( WORK1(i,j) * TZ_UNIFIED(i,j,k,n,bid)                        &
                   + WORK2(i,j) * TZ_UNIFIED(i,j,kp1,n,bid)                    &
                   + WORK3(i,j) * TZ_UNIFIED(i,j+1,k,n,bid)                    &
                   + WORK4(i,j) * TZ_UNIFIED(i,j+1,kp1,n,bid) )
            enddo
          enddo

        end do

      endif ! .not. cancellation_occurs

     do n = 1,nt

!-----------------------------------------------------------------------
!
!     calculate vertical fluxes thru horizontal faces of T-cell
!     - Az(dz*Ax(HYX*KAPPA*SLX*TX)) - Az(dz*Ay(HXY*KAPPA*SLY*TY))
!     calculate isopycnal diffusion from flux differences
!     DTK = (Dx(FX)+Dy(FY)+Dz(FZ)) / volume
!
!-----------------------------------------------------------------------

        GTK(:,:,n) = c0

        if ( k < km ) then

          WORK3 = c0

          if ( .not. cancellation_occurs ) then

!pw loop split to improve performance  -- 2

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                WORK3(i,j) = WORK3(i,j)                               &
                    + ( dz_unified(k) * KAPPA_ISOP_UNIFIED(i,j,kbt,k,bid)             &
                    * ( SLX_UNIFIED(i,j,ieast ,kbt,k  ,bid)                   &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k  ,n,bid)     &
                      + SLY_UNIFIED(i,j,jnorth,kbt,k  ,bid)                   &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k  ,n,bid)     &
                      + SLX_UNIFIED(i,j,iwest ,kbt,k  ,bid)                   &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k  ,n,bid)     &
                      + SLY_UNIFIED(i,j,jsouth,kbt,k  ,bid)                   &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k  ,n,bid) ) )

               enddo
             enddo


            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                WORK3(i,j) = WORK3(i,j)                               &
                    + ( SF_SLX_UNIFIED(i  ,j  ,ieast ,kbt,k  ,bid)            &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k  ,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jnorth,kbt,k  ,bid)            &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k  ,n,bid)     &
                      + SF_SLX_UNIFIED(i  ,j  ,iwest ,kbt,k  ,bid)            &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k  ,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jsouth,kbt,k  ,bid)            &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k  ,n,bid) )

               enddo
             enddo

            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                 WORK3(i,j) = WORK3(i,j)                              &
                    + ( factor                                        &
                    * ( SF_SLX_UNIFIED(i  ,j  ,ieast ,ktp,kp1,bid)            &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,kp1,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jnorth,ktp,kp1,bid)            &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,kp1,n,bid)     &
                      + SF_SLX_UNIFIED(i  ,j  ,iwest ,ktp,kp1,bid)            &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,kp1,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jsouth,ktp,kp1,bid)            &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,kp1,n,bid) ) )

               enddo
             enddo

           do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie

                WORK3(i,j) = WORK3(i,j)                               &
                    + ( dz_bottom * KAPPA_ISOP_UNIFIED(i,j,ktp,kp1,bid)       &
                    * ( SLX_UNIFIED(i  ,j  ,ieast ,ktp,kp1,bid)               &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,kp1,n,bid)     &
                      + SLY_UNIFIED(i  ,j  ,jnorth,ktp,kp1,bid)               &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,kp1,n,bid)     &
                      + SLX_UNIFIED(i  ,j  ,iwest ,ktp,kp1,bid)               &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,kp1,n,bid)     &
                      + SLY_UNIFIED(i  ,j  ,jsouth,ktp,kp1,bid)               &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,kp1,n,bid) ) )

               enddo
             enddo

           do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie


               if(k ==1) then

                 fzprev = c0
                 dzbottomprev = dz_unified(k)

               else ! k == 1

                WORK3prev = c0
                dzbottomprev = dz_unified(k)

                WORK3prev = WORK3prev                                 &  !done
                    + ( dz_unified(k-1) * KAPPA_ISOP_UNIFIED(i,j,kbt,k-1,bid)         &
                    * ( SLX_UNIFIED(i,j,ieast ,kbt,k-1 ,bid)                  &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k-1,n,bid)     &
                      + SLY_UNIFIED(i,j,jnorth,kbt,k-1,bid)                   &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k-1,n,bid)     &
                      + SLX_UNIFIED(i,j,iwest ,kbt,k-1,bid)                   &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k-1,n,bid)     &
                      + SLY_UNIFIED(i,j,jsouth,kbt,k-1,bid)                   &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k-1,n,bid) ) )


                WORK3prev = WORK3prev                                 &  !done
                    + ( SF_SLX_UNIFIED(i  ,j  ,ieast ,kbt,k-1,bid)            &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k-1,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jnorth,kbt,k-1,bid)            &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k-1,n,bid)     &
                      + SF_SLX_UNIFIED(i  ,j  ,iwest ,kbt,k-1,bid)            &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k-1,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jsouth,kbt,k-1,bid)            &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k-1,n,bid) )

                WORK3prev = WORK3prev                               &   !done
                    + ( dzbottomprev * KAPPA_ISOP_UNIFIED(i,j,ktp,k,bid)    &
                    * ( SLX_UNIFIED(i  ,j  ,ieast ,ktp,k,bid)               &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k,n,bid)     &
                      + SLY_UNIFIED(i  ,j  ,jnorth,ktp,k,bid)               &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k,n,bid)     &
                      + SLX_UNIFIED(i  ,j  ,iwest ,ktp,k,bid)               &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k,n,bid)     &
                      + SLY_UNIFIED(i  ,j  ,jsouth,ktp,k,bid)               &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k,n,bid) ) )

                 WORK3prev = WORK3prev                          &       !done
                    + ( c1                                            &
                    * ( SF_SLX_UNIFIED(i  ,j  ,ieast ,ktp,k  ,bid)            &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,k  ,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jnorth,ktp,k  ,bid)            &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,k  ,n,bid)     &
                      + SF_SLX_UNIFIED(i  ,j  ,iwest ,ktp,k  ,bid)            &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,k  ,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jsouth,ktp,k  ,bid)            &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,k,n,bid) ) )

                 KMASKprev = merge(c1, c0, k-1 < KMT_UNIFIED(i,j,bid))

                 fzprev = -KMASKprev * p25 * WORK3prev

             endif ! k == 1 

              fz = -KMASK(i,j) * p25 * WORK3(i,j)


            !if(my_task == master_task .and. nsteps_total == 6 .and. k == 1 .and. i == 3 .and. j == 7 .and. n == 1)then

                   !print *,"changed 1 at GTK write is"
                   !print *,"FX(i,j,n)",FX(i,j,n)
                   !print *,"FX(i-1,j,n)",FX(i-1,j,n)
                   !print *,"FY(i,j,n)", FY(i,j,n)
                   !print *,"FY(i,j-1,n)",FY(i,j-1,n)

              !endif


              GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &               !done
                           + FY(i,j,n) - FY(i,j-1,n)  &
                      + fzprev - fz )*dzr_unified(k)*TAREA_R_UNIFIED(i,j,bid)   

               enddo
            enddo 

          endif !cancellation_errors

       else ! k = km

           do j=this_block%jb,this_block%je
            do i=this_block%ib,this_block%ie

                WORK3prev = c0

                WORK3prev = WORK3prev                                 &   !done
                   + ( dz_unified(km-1) * KAPPA_ISOP_UNIFIED(i,j,kbt,km-1,bid)         &
                    * ( SLX_UNIFIED(i,j,ieast ,kbt,km-1 ,bid)                  &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,km-1,n,bid)     &
                      + SLY_UNIFIED(i,j,jnorth,kbt,km-1,bid)                   &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,km-1,n,bid)     &
                      + SLX_UNIFIED(i,j,iwest ,kbt,km-1,bid)                   &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,km-1,n,bid)     &
                      + SLY_UNIFIED(i,j,jsouth,kbt,km-1,bid)                   &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,km-1,n,bid) ) )

                WORK3prev = WORK3prev                                 &   !done
                    + ( SF_SLX_UNIFIED(i  ,j  ,ieast ,kbt,km-1,bid)            &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,km-1,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jnorth,kbt,km-1,bid)            &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,km-1,n,bid)     &
                      + SF_SLX_UNIFIED(i  ,j  ,iwest ,kbt,km-1,bid)            &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,km-1,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jsouth,kbt,km-1,bid)            &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,km-1,n,bid) )


                WORK3prev = WORK3prev                           &         !done
                    + ( dzbottomprev * KAPPA_ISOP_UNIFIED(i,j,ktp,km,bid)    &
                    * ( SLX_UNIFIED(i  ,j  ,ieast ,ktp,km,bid)               &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,km,n,bid)     &
                      + SLY_UNIFIED(i  ,j  ,jnorth,ktp,km,bid)               &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,km,n,bid)     &
                      + SLX_UNIFIED(i  ,j  ,iwest ,ktp,km,bid)               &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,km,n,bid)     &
                      + SLY_UNIFIED(i  ,j  ,jsouth,ktp,km,bid)               &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,km,n,bid) ) )


                 WORK3prev = WORK3prev                          &
                    + ( c1                                             &
                    * ( SF_SLX_UNIFIED(i  ,j  ,ieast ,ktp,km  ,bid)            &
                       * HYX_UNIFIED(i  ,j  ,bid) * TX_UNIFIED(i  ,j  ,km  ,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jnorth,ktp,km  ,bid)            &
                       * HXY_UNIFIED(i  ,j  ,bid) * TY_UNIFIED(i  ,j  ,km  ,n,bid)     &
                      + SF_SLX_UNIFIED(i  ,j  ,iwest ,ktp,km  ,bid)            &
                       * HYX_UNIFIED(i-1,j  ,bid) * TX_UNIFIED(i-1,j  ,km  ,n,bid)     &
                      + SF_SLY_UNIFIED(i  ,j  ,jsouth,ktp,km  ,bid)            &
                       * HXY_UNIFIED(i  ,j-1,bid) * TY_UNIFIED(i  ,j-1,km,n,bid) ) )


                 KMASKprev = merge(c1, c0, km-1 < KMT_UNIFIED(i,j,bid))

                 fzprev = -KMASKprev * p25 * WORK3prev


                 GTK(i,j,n) = ( FX(i,j,n) - FX(i-1,j,n)  &
                              + FY(i,j,n) - FY(i,j-1,n)  &
                        + fzprev )*dzr_unified(k)*TAREA_R_UNIFIED(i,j,bid)



           enddo
          enddo     


       endif    !k<km 


     enddo ! Tracer Loop


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


 !dir$ attributes offload:mic :: transition_layer_unified
 subroutine transition_layer_unified( this_block )

! !DESCRIPTION:
!  Compute transition layer related fields. the related algorithms
!  should work even with zero transition layer depth.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid,i,j             ! local block address for this sub block

      integer (int_kind), dimension(nx_block,ny_block) :: &
         K_START,   &        ! work arrays for TLT%K_LEVEL and 
         K_SUB               !  TLT%ZTW, respectively

      logical (log_kind), dimension(nx_block,ny_block) :: &
         COMPUTE_TLT         ! flag

      real (r8), dimension(nx_block,ny_block) :: &
         WORK                ! work space for TLT%THICKNESS

      real (r8), dimension(2) :: &
         reference_depth     ! zt or zw

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!----------------------------------------------------------------------- 

      bid = this_block%local_id

      K_START = 0
      K_SUB   = 0

      TLT_UNIFIED%THICKNESS(:,:,bid)      = c0
      TLT_UNIFIED%INTERIOR_DEPTH(:,:,bid) = c0
      TLT_UNIFIED%K_LEVEL(:,:,bid)        = 0
      TLT_UNIFIED%ZTW(:,:,bid)            = 0

      COMPUTE_TLT = merge(.true., .false., KMT_UNIFIED(:,:,bid) /= 0)

      !start_time = omp_get_wtime()
      
      do k=1,km

              !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60) 
              do j=1,ny_block
                   do i=1,nx_block

                      if ( COMPUTE_TLT(i,j)  .and.  &
                         TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid) < zw_unified(k) ) then

                         K_START(i,j) = k+1
                         K_SUB(i,j)   = ktp

                         TLT_UNIFIED%THICKNESS(i,j,bid) = zw_unified(k) - TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)
                         TLT_UNIFIED%K_LEVEL(i,j,bid)   = k
                         TLT_UNIFIED%ZTW(i,j,bid)       = 2

                         COMPUTE_TLT(i,j) = .false.

                      endif

          
                      if ( k /= 1  .and.  K_START(i,j) == k+1  .and. &
                         TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid) < zt_unified(k) ) then
                         K_START(i,j) = k
                         K_SUB(i,j)   = kbt

                         TLT_UNIFIED%THICKNESS(i,j,bid) = zt_unified(k) - TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)
                         TLT_UNIFIED%K_LEVEL(i,j,bid)   = k
                         TLT_UNIFIED%ZTW(i,j,bid)       = 1

                      endif

                   enddo
              enddo
      enddo
  
#ifdef CCSMCOUPLED
#ifndef _HIRES

      if ( any(COMPUTE_TLT) ) then
        print *,"Incorrect DIABATIC_DEPTH value in TLT computation"
      endif

#endif
#endif

             !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60)
             do j=1,ny_block
                   do i=1,nx_block

                      if ( KMT_UNIFIED(i,j,bid) == 0  .or.  K_START(i,j) > KMT_UNIFIED(i,j,bid)  .or.  &
                      ( K_START(i,j) == KMT_UNIFIED(i,j,bid)  .and.  K_SUB(i,j) == kbt ) ) then
                         COMPUTE_TLT(i,j) = .false.
                      else
                         COMPUTE_TLT(i,j) = .true.

                      endif

                   enddo
              enddo

             do k=1,km-1

             !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60)
             do j=1,ny_block
                   do i=1,nx_block

                     WORK(i,j) = c0

                      if ( COMPUTE_TLT(i,j)  .and.  K_SUB(i,j) == kbt  .and.  &
                         K_START(i,j) < KMT_UNIFIED(i,j,bid)  .and.  K_START(i,j) == k ) then

                         WORK(i,j) = max(SLA_SAVE_UNIFIED(i,j,kbt,k,bid), &
                         SLA_SAVE_UNIFIED(i,j,ktp,k+1,bid)) * RB_UNIFIED(i,j,bid)

                      endif

                      if ( WORK(i,j) /= c0  .and.  &
                         TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid) <  (zw_unified(k) - WORK(i,j)) ) then
                         COMPUTE_TLT(i,j) = .false.
                      endif


                      if ( WORK(i,j) /= c0  .and.  &
                         TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid) >= (zw_unified(k) - WORK(i,j)) ) then

                         K_START(i,j) = K_START(i,j) + 1
                         K_SUB(i,j)   = ktp

                         TLT_UNIFIED%THICKNESS(i,j,bid) = zw_unified(k) - TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)
                         TLT_UNIFIED%K_LEVEL(i,j,bid)   = k
                         TLT_UNIFIED%ZTW(i,j,bid)       = 2

                      endif
                   enddo
             enddo
             enddo

      do k=2,km

        reference_depth(ktp) = zt_unified(k)
        reference_depth(kbt) = zw_unified(k)

        do kk=ktp,kbt

              !!$OMP PARALLEL DO
              !DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60)SCHEDULE(DYNAMIC,16)
              do j=1,ny_block
                   do i=1,nx_block
                      if (kk == ktp) then
                   
                         WORK(i,j) = c0 
                         if ( COMPUTE_TLT(i,j)  .and.  K_START(i,j) <= KMT_UNIFIED(i,j,bid)  .and. &
                                                                   K_START(i,j) == k ) then
                         WORK(i,j) = max(SLA_SAVE_UNIFIED(i,j,ktp,k,bid), &
                         SLA_SAVE_UNIFIED(i,j,kbt,k,bid)) * RB_UNIFIED(i,j,bid)
                         endif


                      else
 
                         if ( COMPUTE_TLT(i,j)  .and.  K_START(i,j) < KMT_UNIFIED(i,j,bid)  .and. &
                            K_START(i,j) == k .and. k .lt. km) then
                            WORK(i,j) = max(SLA_SAVE_UNIFIED(i,j,kbt,k,bid), &
                            SLA_SAVE_UNIFIED(i,j,ktp,k+1,bid)) * RB_UNIFIED(i,j,bid)
                         endif


                         if ( COMPUTE_TLT(i,j)  .and.  K_START(i,j) == KMT_UNIFIED(i,j,bid)  .and. &
                               K_START(i,j) == k ) then
                               WORK(i,j) = SLA_SAVE_UNIFIED(i,j,kbt,k,bid) * RB_UNIFIED(i,j,bid)
                         endif

                      endif  

                      if ( WORK(i,j) /= c0  .and.  &
                               TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid) < (reference_depth(kk) - WORK(i,j)) )then 
                               COMPUTE_TLT(i,j) = .false.
                      endif
 
                      if ( WORK(i,j) /= c0  .and.  &
                        TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid) >= (reference_depth(kk) - WORK(i,j)) ) then

                                TLT_UNIFIED%THICKNESS(i,j,bid) = reference_depth(kk)  &
                                    - TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)
                                TLT_UNIFIED%K_LEVEL(i,j,bid)   = k
                                TLT_UNIFIED%ZTW(i,j,bid)       = kk

                       endif

                   enddo
              enddo
          enddo

             do j=1,ny_block
                   do i=1,nx_block

                       if ( COMPUTE_TLT(i,j)  .and.  K_START(i,j) == k ) then
                           K_START(i,j) = K_START(i,j) + 1
                       endif
                   enddo
              enddo  

      enddo

#ifdef CCSMCOUPLED
#ifndef _HIRES
      if ( any(COMPUTE_TLT) ) then
        print *,'Incorrect TLT computations'
      endif
#endif
#endif 

     do k=1,km

             !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60)
             do j=1,ny_block
                   do i=1,nx_block


                      if ( TLT_UNIFIED%K_LEVEL(i,j,bid) == k  .and.  &
                         TLT_UNIFIED%ZTW(i,j,bid) == 1 ) then

                            TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid) = zt_unified(k)
                      endif


                      if ( TLT_UNIFIED%K_LEVEL(i,j,bid) == k  .and.  &
                         TLT_UNIFIED%ZTW(i,j,bid) == 2 ) then

                            TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid) = zw_unified(k)
                      endif

                   enddo
              enddo
      enddo

             !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i)NUM_THREADS(60)
             do j=1,ny_block
                   do i=1,nx_block


                      COMPUTE_TLT(i,j) = .false.

                      if ( KMT_UNIFIED(i,j,bid) /= 0  .and.  &
                           TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid) == c0 ) then 
                           
                           COMPUTE_TLT(i,j) = .true.

                      endif    

                      if ( KMT_UNIFIED(i,j,bid) == 0  .and.  &
                            TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid) /= c0 )  then

                            COMPUTE_TLT(i,j) = .true.

                      endif  

                   enddo
              enddo


 end subroutine transition_layer_unified

 !dir$ attributes offload:mic :: buoyancy_frequency_dependent_profile_unified
 subroutine buoyancy_frequency_dependent_profile_unified (TMIX, this_block)

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k,        &          ! vertical loop index
         bid,i,j              ! local block address for this sub block

      integer (int_kind), dimension(nx_block,ny_block) :: &
         K_MIN                ! k index below SDL 

      real (r8), dimension(nx_block,ny_block) :: &
         TEMP_K,           &  ! temperature at level k
         TEMP_KP1,         &  ! temperature at level k+1
         RHOT,             &  ! dRHO/dT
         RHOS,             &  ! dRHO/dS
         BUOY_FREQ_SQ_REF, &  ! reference (normalization) value
                              !  of N^2
         SDL                  ! surface diabatic layer (see below)

      real (r8), dimension(nx_block,ny_block,km) :: &
         BUOY_FREQ_SQ_NORM    ! normalized N^2 defined at level interfaces


      bid = this_block%local_id

      BUOY_FREQ_SQ_NORM         = c0
      BUOY_FREQ_SQ_REF          = c0
      KAPPA_VERTICAL_UNIFIED(:,:,:,bid) = c1

      K_MIN = merge( km+1, 0, KMT_UNIFIED(:,:,bid) /= 0 )

      
      SDL = KPP_HBLT_UNIFIED(:,:,bid)
      if(transition_layer_on) SDL = TLT_UNIFIED%INTERIOR_DEPTH(:,:,bid)

        do k=1,km-1

          !if ( k == 1 ) TEMP_K = max( -c2, TMIX(:,:,k,1) )

          !TEMP_KP1 = max( -c2, TMIX(:,:,k+1,1) )

          call state_unified( k, k+1, TMIX(:,:,k,1), TMIX(:,:,k,2), &
                            this_block, DRHODT=RHOT, DRHODS=RHOS )

          !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(i,j)NUM_THREADS(60) 
          do j=1,ny_block
             do i=1,nx_block

                if ( k == 1 ) TEMP_K(i,j) = max( -c2, TMIX(i,j,k,1) )

                 TEMP_KP1(i,j) = max( -c2, TMIX(i,j,k+1,1) )



                 if ( k < KMT_UNIFIED(i,j,bid) ) then
                    BUOY_FREQ_SQ_UNIFIED(i,j,k,bid) = max( c0, - grav * dzwr_unified(k) &
                        * ( RHOT(i,j) * ( TEMP_K(i,j) - TEMP_KP1(i,j) ) &
                  + RHOS(i,j) * ( TMIX(i,j,k,  2) - TMIX(i,j,k+1,2) ) ) )
                 endif

                 TEMP_K(i,j) = TEMP_KP1(i,j)

             enddo
           enddo

        enddo

     do k=1,km-1
       !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(i,j)NUM_THREADS(60)
       do j=1,ny_block
        do i=1,nx_block

         if ( ( K_MIN(i,j) == km+1 ) .and. ( zw_unified(k) > SDL(i,j) ) .and.  &
                ( k <= KMT_UNIFIED(i,j,bid) )  .and.                 &
                ( BUOY_FREQ_SQ_UNIFIED(i,j,k,bid) > c0 ) ) then
          BUOY_FREQ_SQ_REF(i,j) = BUOY_FREQ_SQ_UNIFIED(i,j,k,bid)
          K_MIN(i,j) = k

         endif

         if ( ( k >= K_MIN(i,j) ) .and. ( k < KMT_UNIFIED(i,j,bid) ) .and. &
                  ( BUOY_FREQ_SQ_REF(i,j) /= c0 ) ) then
                      BUOY_FREQ_SQ_NORM(i,j,k) =  &
                      max( BUOY_FREQ_SQ_UNIFIED(i,j,k,bid) / BUOY_FREQ_SQ_REF(i,j),0.1_r8 )
                      BUOY_FREQ_SQ_NORM(i,j,k) =  &
                      min( BUOY_FREQ_SQ_NORM(i,j,k), c1 )
         else
             BUOY_FREQ_SQ_NORM(i,j,k) = c1
         endif

            enddo
         enddo
      enddo

      do k=1,km-1
       !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(i,j)NUM_THREADS(60)
       do j=1,ny_block
        do i=1,nx_block

        if ( k == KMT_UNIFIED(i,j,bid)-1 ) then
          BUOY_FREQ_SQ_NORM(i,j,k+1) = BUOY_FREQ_SQ_NORM(i,j,k)
        endif

        enddo
       enddo
      enddo

      do k=2,km
       !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(i,j)NUM_THREADS(60)
       do j=1,ny_block
        do i=1,nx_block

         if ( ( k > K_MIN(i,j) ) .and. ( k <= KMT_UNIFIED(i,j,bid) ) ) then
          KAPPA_VERTICAL_UNIFIED(i,j,k,bid) = BUOY_FREQ_SQ_NORM(i,j,k-1)
         endif

        enddo
       enddo
      enddo



 end subroutine buoyancy_frequency_dependent_profile_unified

      !dir$ attributes offload:mic :: merged_streamfunction_unified
 subroutine merged_streamfunction_unified ( this_block )

! !DESCRIPTION:
!  Construct a merged streamfunction that has the appropriate
!  behavior in the surface diabatic region, transition layer, and
!  adiabatic interior
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid,i,j,temp        ! local block address for this sub block
 
      real (r8), dimension(nx_block,ny_block,2) :: &
         WORK1, WORK2, WORK3, WORK4   ! work arrays

      real (r8), dimension(nx_block,ny_block) :: &
         WORK2_NEXT, WORK4_NEXT       ! WORK2 or WORK4 at next level

      real (r8), dimension(nx_block,ny_block) :: &
         WORK5, WORK6, WORK7          ! more work arrays

      logical (log_kind), dimension(nx_block,ny_block) :: &
         LMASK               ! flag

      real (r8), dimension(2) :: &
         reference_depth              ! zt or zw

      !real (r8) :: &
      !   start_time,end_time



!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------


      !start_time = omp_get_wtime() 
      bid = this_block%local_id

     !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(k,kk,temp,j,i)num_threads(60)collapse(4)
     do k=1,km
        do kk=1,2
           do temp=1,2 
              do j=1,ny_block
                 !!dir$ vector aligned
                 !!dir$ ivdep
                 do i=1,nx_block
 
                    SF_SLX_UNIFIED(i,j,temp,kk,k,bid) = c0
                    SF_SLY_UNIFIED(i,j,temp,kk,k,bid) = c0

                 enddo
              enddo
            enddo
          enddo
      enddo   

      !end_time = omp_get_wtime()
      !print *,end_time - start_time


      WORK1 = c0
      WORK2 = c0
      WORK3 = c0
      WORK4 = c0
      WORK5 = c0
      WORK6 = c0
      WORK7 = c0
      WORK2_NEXT = c0
      WORK4_NEXT = c0
 


!-----------------------------------------------------------------------
!
!     compute the interior streamfunction and its first derivative at the
!     INTERIOR_DEPTH level. WORK1 and WORK2 contain the streamfunction
!     and its first derivative, respectively, for the zonal component
!     for the east and west sides of a grid cell. WORK3 and WORK4 are
!     the corresponding fields for the meridional component for the
!     north and south sides of a grid cell. Note that these definitions
!     include a "dz". Also, the first derivative computations assume 
!     that the streamfunctions are located in the middle of the top or
!     bottom half of a grid cell, hence a factor of two in WORK2 and
!     WORK4 calculations. 
!
!-----------------------------------------------------------------------

      !start_time = omp_get_wtime()   
   
      do k=1,km-1

        do kk=1,2
 
          !!$OMP PARALLEL DO PRIVATE(I,J)DEFAULT(SHARED)NUM_THREADS(60) 
          do j=1,ny_block
           do i=1,nx_block

 
            LMASK(i,j) = TLT_UNIFIED%K_LEVEL(i,j,bid) == k  .and.            &
                       TLT_UNIFIED%K_LEVEL(i,j,bid) < KMT_UNIFIED(i,j,bid)  .and.  &
                       TLT_UNIFIED%ZTW(i,j,bid) == 1


            if ( LMASK(i,j) ) then 

             WORK1(i,j,kk) =  KAPPA_THIC_UNIFIED(i,j,kbt,k,bid)  &
                           * SLX_UNIFIED(i,j,kk,kbt,k,bid) * dz_unified(k)

             WORK2(i,j,kk) = c2 * dzwr_unified(k) * ( WORK1(i,j,kk)            &
              - KAPPA_THIC_UNIFIED(i,j,ktp,k+1,bid) * SLX_UNIFIED(i,j,kk,ktp,k+1,bid) &
                                            * dz_unified(k+1) )

      !if(my_task == master_task .and. nsteps_total == 6 .and. k == 2.and. i==3 .and. j==7 .and. kk == 1)then


         !print *,"CHanged1"
         !print *,"WORK2(i,j)",WORK2(3,7,1),k

      !endif


             WORK2_NEXT(i,j) = c2 * ( &
              KAPPA_THIC_UNIFIED(i,j,ktp,k+1,bid) * SLX_UNIFIED(i,j,kk,ktp,k+1,bid) - &
              KAPPA_THIC_UNIFIED(i,j,kbt,k+1,bid) * SLX_UNIFIED(i,j,kk,kbt,k+1,bid) )

             WORK3(i,j,kk) =  KAPPA_THIC_UNIFIED(i,j,kbt,k,bid)  &
                           * SLY_UNIFIED(i,j,kk,kbt,k,bid) * dz_unified(k)

             WORK4(i,j,kk) = c2 * dzwr_unified(k) * ( WORK3(i,j,kk)            &
              - KAPPA_THIC_UNIFIED(i,j,ktp,k+1,bid) * SLY_unified(i,j,kk,ktp,k+1,bid) &
                                            * dz_unified(k+1) )

             WORK4_NEXT(i,j) = c2 * ( &
              KAPPA_THIC_UNIFIED(i,j,ktp,k+1,bid) * SLY_UNIFIED(i,j,kk,ktp,k+1,bid) - &
              KAPPA_THIC_UNIFIED(i,j,kbt,k+1,bid) * SLY_UNIFIED(i,j,kk,kbt,k+1,bid) )

            endif

            if( LMASK(i,j) .and. abs( WORK2_NEXT(i,j) ) < abs( WORK2(i,j,kk) ) )then 

             WORK2(i,j,kk) = WORK2_NEXT(i,j)

      !if(my_task == master_task .and. nsteps_total == 6 .and. k == 2 .and. i==3 .and. j==7 .and. kk == 1)then


         !print *,"CHanged2"
         !print *,"WORK2(i,j)",WORK2(3,7,1),k

      !endif


            endif

           if ( LMASK(i,j) .and. abs( WORK4_NEXT(i,j) ) < abs( WORK4(i,j,kk ) )) then 
             WORK4(i,j,kk) = WORK4_NEXT(i,j)
           endif

          LMASK(i,j) = TLT_UNIFIED%K_LEVEL(i,j,bid) == k  .and.           &
                       TLT_UNIFIED%K_LEVEL(i,j,bid) < KMT_UNIFIED(i,j,bid)  .and. &
                       TLT_UNIFIED%ZTW(i,j,bid) == 2

          if ( LMASK(i,j) ) then

            WORK1(i,j,kk) =  KAPPA_THIC_UNIFIED(i,j,ktp,k+1,bid)     & 
                           * SLX_UNIFIED(i,j,kk,ktp,k+1,bid)

            WORK2(i,j,kk) =  c2 * ( WORK1(i,j,kk)                 &
                           - ( KAPPA_THIC_UNIFIED(i,j,kbt,k+1,bid)        &
                              * SLX_UNIFIED(i,j,kk,kbt,k+1,bid) ) )

      !if(my_task == master_task .and. nsteps_total == 6 .and. k == 2 .and. i==3 .and. j==7 .and. kk == 1)then


         !print *,"CHanged3"
         !print *,"WORK2(i,j)",WORK2(3,7,1),k

      !endif


            WORK1(i,j,kk) = WORK1(i,j,kk) * dz_unified(k+1)

            WORK3(i,j,kk) =  KAPPA_THIC_UNIFIED(i,j,ktp,k+1,bid)     &
                           * SLY_UNIFIED(i,j,kk,ktp,k+1,bid)

            WORK4(i,j,kk) =  c2 * ( WORK3(i,j,kk)                 &
                           - ( KAPPA_THIC_UNIFIED(i,j,kbt,k+1,bid)        &
                              * SLY_UNIFIED(i,j,kk,kbt,k+1,bid) ) )

            WORK3(i,j,kk) = WORK3(i,j,kk) * dz_unified(k+1)

            endif
 
          LMASK(i,j) = LMASK(i,j) .and. TLT_UNIFIED%K_LEVEL(i,j,bid) + 1 < KMT_UNIFIED(i,j,bid)

          if (k.lt.km-1) then ! added to avoid out of bounds access

            if( LMASK(i,j) ) then

         
      !if(my_task == master_task .and. nsteps_total == 6 .and. k == 2 .and. i==3 .and. j==7 .and. kk == 1)then


         !print *,"Changed"
         !print *,"WORK2(i,j)",WORK2(3,7,1),k
         !print *,"WORK2_NEXT(i,j)",WORK2_NEXT(i,j)
         !print *,"KAPPA_THIC_UNIFIED(i,j,kbt,k+1,bid)",KAPPA_THIC_UNIFIED(i,j,kbt,k+1,bid)
         !print *,"KAPPA_THIC_UNIFIED(i,j,ktp,k+2,bid)",KAPPA_THIC_UNIFIED(i,j,ktp,k+2,bid)

      !endif

               

              WORK2_NEXT(i,j) = c2 * dzwr_unified(k+1) * ( &
                KAPPA_THIC_UNIFIED(i,j,kbt,k+1,bid) * SLX_UNIFIED(i,j,kk,kbt,k+1,bid) * dz_unified(k+1)- &
                KAPPA_THIC_UNIFIED(i,j,ktp,k+2,bid) * SLX_UNIFIED(i,j,kk,ktp,k+2,bid) * dz_unified(k+2))

              WORK4_NEXT(i,j) = c2 * dzwr_unified(k+1) * ( &
                KAPPA_THIC_UNIFIED(i,j,kbt,k+1,bid) * SLY_UNIFIED(i,j,kk,kbt,k+1,bid) * dz_unified(k+1)- &
                KAPPA_THIC_UNIFIED(i,j,ktp,k+2,bid) * SLY_UNIFIED(i,j,kk,ktp,k+2,bid) * dz_unified(k+2))

              endif 

          end if
             
          if( LMASK(i,j) .and. abs( WORK2_NEXT(i,j) ) < abs( WORK2(i,j,kk) ) ) &
            WORK2(i,j,kk) = WORK2_NEXT(i,j)

      !if(my_task == master_task .and. nsteps_total == 6 .and. k == 2 .and. i==3 .and. j==7 .and. kk == 1)then


         !print *,"CHanged4"
         !print *,"WORK2(i,j)",WORK2(3,7,1),k

      !endif


          if( LMASK(i,j) .and. abs(WORK4_NEXT(i,j)) < abs(WORK4(i,j,kk)) ) &
            WORK4(i,j,kk) = WORK4_NEXT(i,j)

             enddo
          enddo
          !!$OMP END PARALLEL DO
        enddo
      enddo



      !if(my_task == master_task .and. nsteps_total == 6 )then


         !print *,"Changed1"
         !print *,"WORK2(i,j)",WORK2(3,7,1)

      !endif

!-----------------------------------------------------------------------
!
!     compute the depth independent interpolation factors used in the 
!     linear and quadratic interpolations within the diabatic and 
!     transition regions, respectively.
!
!-----------------------------------------------------------------------

          do j=1,ny_block
             do i=1,nx_block

                WORK5(i,j) = c0
                if (KMT_UNIFIED(i,j,bid) /= 0) then
                   WORK5(i,j) = c1 / ( c2 * TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid) &
                   + TLT_UNIFIED%THICKNESS(i,j,bid) )
 
                endif 

                WORK6(i,j) = c0
                if ((KMT_UNIFIED(i,j,bid) /= 0) .AND. (TLT_UNIFIED%THICKNESS(i,j,bid) > eps))then
                   WORK6(i,j) = WORK5(i,j) / TLT_UNIFIED%THICKNESS(i,j,bid)
      
                endif

             enddo
          enddo

      !start_time = omp_get_wtime()
      !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(i,j,k,kk,reference_depth)num_threads(60)SCHEDULE(DYNAMIC,6) 
      do k=1,km

        reference_depth(ktp) = zt_unified(k) - p25 * dz_unified(k)
        reference_depth(kbt) = zt_unified(k) + p25 * dz_unified(k)

        do kk=ktp,kbt

!-----------------------------------------------------------------------
!
!     diabatic region: use linear interpolation (in streamfunction) 
!
!-----------------------------------------------------------------------
     
          do j=1,ny_block
             do i=1,nx_block


                if ( reference_depth(kk) <= TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)  &
                       .and.  k <= KMT_UNIFIED(i,j,bid) ) then

                 !if(my_task == master_task .and. nsteps_total == 6 .and. k == 1 .and. i == 3 .and. j == 7 .and. kk == 1)then


                 !print *,"Changed1"
                 !print *,"reference_depth(kk)",reference_depth(kk)
                 !print *,"TLT%DIABATIC_DEPTH(i,j,bid)",TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)
                 !print *,"WORK1(i,j,1)",WORK1(i,j,1)
                 !print *,"WORK2(i,j,1)",WORK2(i,j,1)
                 !print *," KMT(i,j,bid)", KMT_UNIFIED(i,j,bid)
                 !print *,"WORK5(i,j)",WORK5(i,j)
                 !print *,"SF_SLX(i,j,1,kk,k,bid)",SF_SLX_UNIFIED(i,j,1,kk,k,bid)

                 !endif


                       SF_SLX_UNIFIED(i,j,1,kk,k,bid) = reference_depth(kk) * WORK5(i,j)  &
                             * ( c2 * WORK1(i,j,1) + TLT_UNIFIED%THICKNESS(i,j,bid)       &
                                * WORK2(i,j,1) )

                       SF_SLX_UNIFIED(i,j,2,kk,k,bid) = reference_depth(kk) * WORK5(i,j)  &
                              * ( c2 * WORK1(i,j,2) + TLT_UNIFIED%THICKNESS(i,j,bid)      &
                                * WORK2(i,j,2) )

                       SF_SLY_UNIFIED(i,j,1,kk,k,bid) = reference_depth(kk) * WORK5(i,j)  &
                             * ( c2 * WORK3(i,j,1) + TLT_UNIFIED%THICKNESS(i,j,bid)       &
                                * WORK4(i,j,1) )

                       SF_SLY_UNIFIED(i,j,2,kk,k,bid) = reference_depth(kk) * WORK5(i,j)  &
                             * ( c2 * WORK3(i,j,2) + TLT_UNIFIED%THICKNESS(i,j,bid)       &
                                * WORK4(i,j,2) )
     
                endif
      

!-----------------------------------------------------------------------
!
!     transition layer: use quadratic interpolation (in streamfunction) 
!
!-----------------------------------------------------------------------


                 if ( reference_depth(kk) > TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)   &
                .and.  reference_depth(kk) <= TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid) &
                .and.  k <= KMT_UNIFIED(i,j,bid) ) then


                 !if(my_task == master_task .and. nsteps_total == 6 .and. k == 1 .and. i == 3 .and. j == 7 .and. kk == 1 )then


                 !print *,"Changed2"
                 !print *,"reference_depth(kk)",reference_depth(kk)
                 !print *,"TLT%DIABATIC_DEPTH(i,j,bid)",TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)
                 !print *,"WORK1(i,j,1)",WORK1(i,j,1)
                 !print *,"WORK2(i,j,1)",WORK2(i,j,1)
                 !print *," KMT(i,j,bid)", KMT_UNIFIED(i,j,bid)
                 !print *,"WORK5(i,j)",WORK5(i,j)
                 !print *,"SF_SLX(i,j,1,kk,k,bid)",SF_SLX_UNIFIED(i,j,1,kk,k,bid)
                 !print *,"WORK7(i,j)",WORK7(i,j)

                 !endif


                       WORK7(i,j) = (TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)  &
                                       - reference_depth(kk))**2

                      SF_SLX_UNIFIED(i,j,1,kk,k,bid) = - WORK7(i,j) * WORK6(i,j)  &
                          * ( WORK1(i,j,1) + TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)  &
                             * WORK2(i,j,1) )                             &
                         + reference_depth(kk) * WORK5(i,j)               &
                          * ( c2 * WORK1(i,j,1) + TLT_UNIFIED%THICKNESS(i,j,bid)  &
                                     * WORK2(i,j,1) )

                      SF_SLX_UNIFIED(i,j,2,kk,k,bid) = - WORK7(i,j) * WORK6(i,j)  &
                          * ( WORK1(i,j,2) + TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)  &
                             * WORK2(i,j,2) )                             &
                         + reference_depth(kk) * WORK5(i,j)               &
                          * ( c2 * WORK1(i,j,2) + TLT_UNIFIED%THICKNESS(i,j,bid)  &
                                       * WORK2(i,j,2) )

                      SF_SLY_UNIFIED(i,j,1,kk,k,bid) = - WORK7(i,j) * WORK6(i,j)  &
                          * ( WORK3(i,j,1) + TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)  &
                             * WORK4(i,j,1) )                             &
                         + reference_depth(kk) * WORK5(i,j)               &
                          * ( c2 * WORK3(i,j,1) + TLT_UNIFIED%THICKNESS(i,j,bid)  &
                             * WORK4(i,j,1) )

                      SF_SLY_UNIFIED(i,j,2,kk,k,bid) = - WORK7(i,j) * WORK6(i,j)  &
                          * ( WORK3(i,j,2) + TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)  &
                             * WORK4(i,j,2) )                             &
                         + reference_depth(kk) * WORK5(i,j)               &
                          * ( c2 * WORK3(i,j,2) + TLT_UNIFIED%THICKNESS(i,j,bid)  &
                             * WORK4(i,j,2) )
 
                  endif


!-----------------------------------------------------------------------
!
!     interior, adiabatic region: no interpolation is needed. note that
!     "dzw" is introduced here, too, for consistency. 
!
!-----------------------------------------------------------------------

                 if ( reference_depth(kk) > TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)  & 
                       .and.  k <= KMT_UNIFIED(i,j,bid) ) then


               !if(my_task == master_task .and. nsteps_total == 6 .and. k == 1 .and. i == 3 .and. j == 7 .and. kk == 1)then


                 !print *,"changed3"
                 !print *,"reference_depth(kk)",reference_depth(kk)
                 !print *,"TLT%DIABATIC_DEPTH(i,j,bid)",TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)
                 !print *," KMT(i,j,bid)", KMT_UNIFIED(i,j,bid)
                 !print *,"KAPPA_THIC(i,j,kk,k,bid)",KAPPA_THIC_UNIFIED(i,j,kk,k,bid)
                 !print *,"dz(k)",dz_unified(k)
                 !print *,"SF_SLX(i,j,1,kk,k,bid)",SF_SLX_UNIFIED(i,j,1,kk,k,bid)

                !endif


                     SF_SLX_UNIFIED(i,j,1,kk,k,bid) =  KAPPA_THIC_UNIFIED(i,j,kk,k,bid)  &
                                       * SLX_UNIFIED(i,j,1,kk,k,bid) * dz_unified(k)

                     SF_SLX_UNIFIED(i,j,2,kk,k,bid) =  KAPPA_THIC_UNIFIED(i,j,kk,k,bid)  &
                                       * SLX_UNIFIED(i,j,2,kk,k,bid) * dz_unified(k)

                     SF_SLY_UNIFIED(i,j,1,kk,k,bid) =  KAPPA_THIC_UNIFIED(i,j,kk,k,bid)  &
                                       * SLY_UNIFIED(i,j,1,kk,k,bid) * dz_unified(k)

                     SF_SLY_UNIFIED(i,j,2,kk,k,bid) =  KAPPA_THIC_UNIFIED(i,j,kk,k,bid)  &
                                       * SLY_UNIFIED(i,j,2,kk,k,bid) * dz_unified(k)

                endif

             enddo
          enddo

        enddo  ! end of kk-loop

      enddo    ! end of k-loop



 end subroutine merged_streamfunction_unified

!***********************************************************************
!BOP
! !IROUTINE: apply_vertical_profile_to_isop_hor_diff 
! !INTERFACE:

      !dir$ attributes offload:mic :: apply_vertical_profile_to_isop_hor_diff_unified   
      subroutine apply_vertical_profile_to_isop_hor_diff_unified ( this_block ) 

! !DESCRIPTION:
!  Apply vertical tapers to KAPPA_ISOP and HOR_DIFF based on their
!  vertical location with respect to the diabatic, transition, and
!  adiabatic regions.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid,i,j             ! local block address for this sub block

      real (r8), dimension(2) :: &
         reference_depth

      bid = this_block%local_id

!-----------------------------------------------------------------------
!
!     start of tapering
!
!-----------------------------------------------------------------------
      !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(i,j,k,kk,reference_depth)num_threads(60)
      do k=1,km

        reference_depth(ktp) = zt_unified(k) - p25 * dz_unified(k)
        reference_depth(kbt) = zt_unified(k) + p25 * dz_unified(k)

        do kk=ktp,kbt

!-----------------------------------------------------------------------
!
!     diabatic region: no isopycnal diffusion 
!
!-----------------------------------------------------------------------
           do j=1,ny_block
              do i=1,nx_block

                 if ( reference_depth(kk) <= TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)  &
                 .and.  k <= KMT_UNIFIED(i,j,bid) ) then

                    KAPPA_ISOP_UNIFIED(i,j,kk,k,bid) = c0

                endif
  

!-----------------------------------------------------------------------
!
!      transition layer: a linear combination of isopcynal and horizontal
!      diffusion coefficients 
!
!-----------------------------------------------------------------------

                if( reference_depth(kk) > TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid)   &
              .and.  reference_depth(kk) <= TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid) &
              .and.  k <= KMT_UNIFIED(i,j,bid)  .and.                           &
                     TLT_UNIFIED%THICKNESS(i,j,bid) > eps ) then

                     HOR_DIFF_UNIFIED(i,j,kk,k,bid) = ( TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)  &
                           - reference_depth(kk) ) * HOR_DIFF_UNIFIED(i,j,kk,k,bid)  &
                                 / TLT_UNIFIED%THICKNESS(i,j,bid)

                     KAPPA_ISOP_UNIFIED(i,j,kk,k,bid) = ( reference_depth(kk)        &
                                            - TLT_UNIFIED%DIABATIC_DEPTH(i,j,bid) )  &
                          * KAPPA_ISOP_UNIFIED(i,j,kk,k,bid) / TLT_UNIFIED%THICKNESS(i,j,bid)

                 endif

!-----------------------------------------------------------------------
!
!     interior region: no horizontal diffusion
!
!-----------------------------------------------------------------------
                 !if(my_task == master_task .and. nsteps_total == 1 .and. k == 1 .and. i == 3 .and. j == 5 .and. kk == ktp)then

                 !print *,"Changed before",HOR_DIFF_UNIFIED(i,j,kk,k,bid) 

                 !endif


                 if ( reference_depth(kk) > TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)  &
                  .and.  k <= KMT_UNIFIED(i,j,bid) ) then

                          HOR_DIFF_UNIFIED(i,j,kk,k,bid) = c0
                 endif

                 !if(my_task == master_task .and. nsteps_total == 1 .and. k == 1 .and. i == 3 .and. j == 5 .and. kk == ktp)then


                 !print *,"Changed"
                 !print *,"HOR_UNIFIED(i,j,ktp,k,bid)",HOR_DIFF_UNIFIED(i,j,kk,k,bid)
                 !print *,"reference_depth(kk)",reference_depth(kk)
                 !print *,"TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)",TLT_UNIFIED%INTERIOR_DEPTH(i,j,bid)   
                 !print *,"KMT_UNIFIED(i,j,bid)",KMT_UNIFIED(i,j,bid)                 

                 !endif


              enddo
           enddo

        enddo  ! end of kk-loop

      enddo    ! end of k-loop


     !print *,KAPPA_ISOP(45,45,1,45,bid),HOR_DIFF(45,45,1,45,bid) 
     !if(my_task==master_task)then

       !open(unit=10,file="/home/aketh/ocn_correctness_data/changed.txt",status="unknown",position="append",action="write",form="unformatted")
       !write(10),HOR_DIFF,KAPPA_ISOP
       !close(10)

     !endif


!-----------------------------------------------------------------------
!EOC

      end subroutine apply_vertical_profile_to_isop_hor_diff_unified

 subroutine merger (TCUR , ARRAY , iblock , this_block )

 !-----------INPUT VARAIBLES-----------------------------------! 

 real (r8), dimension(nx_block,ny_block), intent(in) :: TCUR 

 integer (int_kind), intent(in) :: iblock

 type (block), intent(in) ::       &
      this_block           ! block information for current block

 !-----------OUTPUT VARIABLES----------------------------------!

 real (r8), dimension(nx_block_unified,ny_block_unified), intent(out) :: ARRAY

 !local variables

   integer (int_kind) :: k

   integer (int_kind) :: my_grid_blockno, block_row, &
   block_col,i_start,j_start,i_end,j_end,ib,ie,jb,je,i_index,j_index

   !logical (log_kind) :: written(164,196,60)

   !integer (int_kind) :: written_byi(164,196,60)

   !integer (int_kind) :: written_byj(164,196,60)

   !integer (int_kind) :: written_byk(164,196,60)

   !integer (int_kind) :: written_by_block(164,196,60)
  
   integer (int_kind) :: i,j  

 !-------------------------------------------------------------!


         my_grid_blockno = iblock - 1

         block_row = int( my_grid_blockno / 4  )

         block_col = mod(my_grid_blockno,4)

         i_start = block_col * (nx_block - 4) + 1 + 2

         j_start = block_row * (ny_block - 4) + 1 + 2

         i_end = i_start + (this_block%ie - this_block%ib) 

         j_end = j_start + (this_block%je - this_block%jb)

         ib = this_block%ib
 
         jb = this_block%jb

         ie = this_block%ie 

         je = this_block%je

         if(block_row == 0 ) then

         j_start = 1
         jb = 1 

         endif

         if(block_row == 3 ) then

         j_start = j_start
         je = this_block%je + 2

         endif
 
         if(block_col == 0 ) then

         i_start = 1
         ib = 1  

         endif 

         if(block_col == 3 ) then

         i_start = i_start
         ie = this_block%ie + 2

         endif
 

            j_index = j_start
             do j=jb,je
                   i_index = i_start
                     do i=ib,ie

                       ARRAY(i_index,j_index) = TCUR(i,j)

                       i_index = i_index + 1


                      end do
                j_index = j_index + 1
            end do
 
 end subroutine merger 



 subroutine splitter (SPLIT_ARRAY , MERGED_ARRAY , iblock , this_block )

 !-----------INPUT VARAIBLES-----------------------------------! 

 real (r8), dimension(nx_block_unified,ny_block_unified), intent(in) :: MERGED_ARRAY

 integer (int_kind), intent(in) :: iblock

 type (block), intent(in) ::       &
      this_block           ! block information for current block

 !-----------OUTPUT VARIABLES----------------------------------!

 real (r8), dimension(nx_block,ny_block), intent(out) :: SPLIT_ARRAY

 !local variables

   integer (int_kind) :: k

   integer (int_kind) :: my_grid_blockno, block_row, &
   block_col,i_start,j_start,i_end,j_end,ib,ie,jb,je,i_index,j_index

   !logical (log_kind) :: written(164,196,60)

   !integer (int_kind) :: written_byi(164,196,60)

   !integer (int_kind) :: written_byj(164,196,60)

   !integer (int_kind) :: written_byk(164,196,60)

   !integer (int_kind) :: written_by_block(164,196,60)
  
   integer (int_kind) :: i,j  

 !-------------------------------------------------------------!


         my_grid_blockno = iblock - 1

         block_row = int( my_grid_blockno / 4  )

         block_col = mod(my_grid_blockno,4)

         i_start = block_col * (nx_block - 4) + 1

         j_start = block_row * (ny_block - 4) + 1

         j_index = j_start
           do j=1,ny_block
               i_index = i_start
                 do i=1,nx_block

                     SPLIT_ARRAY(i,j) = MERGED_ARRAY(i_index,j_index)
                     i_index = i_index + 1


                  end do
               j_index = j_index + 1
           end do

 
 end subroutine splitter 

 end module horizontal_mix_unified

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
