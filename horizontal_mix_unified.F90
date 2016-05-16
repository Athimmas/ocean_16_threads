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
   use distribution, only: 
   use domain_size
   use domain, only: nblocks_clinic, distrb_clinic
   use constants, only: c0, blank_fmt, delim_fmt, ndelim_fmt
   use communicate, only: my_task, master_task
   use time_management, only: km, nt, mix_pass
   use broadcast, only: broadcast_scalar
   use grid, only: KMT, dz, partial_bottom_cells, DZT, dzr, dzwr,KMTE, KMT, KMTN
   use io_types, only: nml_in, nml_filename, stdout
   use hmix_del2, only: init_del2u, init_del2t, hdiffu_del2, hdifft_del2
   use hmix_del4, only: init_del4u, init_del4t, hdiffu_del4, hdifft_del4
   use hmix_gm, only: init_gm, hdifft_gm
   use hmix_aniso, only: init_aniso, hdiffu_aniso
   use topostress, only: ltopostress
   use horizontal_mix, only:tavg_HDIFE_TRACER,tavg_HDIFN_TRACER,tavg_HDIFB_TRACER
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now, &
      tavg_in_which_stream, ltavg_on
   use timers, only: timer_start, timer_stop, get_timer
   use exit_mod, only: sigAbort, exit_pop, flushm
   use mix_submeso, only: init_submeso, submeso_flux, submeso_sf
   use hmix_gm_submeso_share, only: init_meso_mixing, tracer_diffs_and_isopyc_slopes
   use prognostic
   use vertical_mix
   use omp_lib
   use state_mod, only: pressz

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

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_horizontal_mix
! !INTERFACE:

 subroutine init_horizontal_mix_unified()


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

   pressz_unified = pressz
   KMT_unified = KMT
   KMTE_unified = KMTE
   KMTN_unified = KMTN

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

         !call hdifft_gm(1, HDTK_BUF(:,:,:,1), TMIX, UMIX, VMIX, tavg_HDIFE_TRACER, &
         !                tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

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
 
        !if (k == 1) then
         !start_time = omp_get_wtime()
         !call submeso_sf(TMIX, this_block)
         !end_time = omp_get_wtime()
         !print *,"time at submeso_sf is ",end_time - start_time
        !endif

        !if(k==1) then
         !start_time = omp_get_wtime()
        !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)num_threads(60) 
        !do kk=1,km
         !call submeso_flux(kk, TDTK(:,:,:,kk), TMIX, tavg_HDIFE_TRACER, &
         !              tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)
        !enddo
        !end_time = omp_get_wtime()
        !print *,"time at submeso_flux is ",end_time - start_time
        !endif
        !HDTK=HDTK+TDTK(:,:,:,k)
   

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

 end module horizontal_mix_unified

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
