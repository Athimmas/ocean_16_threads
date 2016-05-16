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
   use grid, only: KMT, dz, partial_bottom_cells, DZT, dzr, dzwr
   use io_types, only: nml_in, nml_filename, stdout
   use hmix_del2, only: init_del2u, init_del2t, hdiffu_del2, hdifft_del2
   use hmix_del4, only: init_del4u, init_del4t, hdiffu_del4, hdifft_del4
   use hmix_gm, only: init_gm, hdifft_gm
   use hmix_aniso, only: init_aniso, hdiffu_aniso
   use topostress, only: ltopostress
   use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now, &
      tavg_in_which_stream, ltavg_on
   use timers, only: timer_start, timer_stop, get_timer
   use exit_mod, only: sigAbort, exit_pop, flushm
   use mix_submeso, only: init_submeso, submeso_flux, submeso_sf
   use hmix_gm_submeso_share, only: init_meso_mixing, tracer_diffs_and_isopyc_slopes
   use prognostic
   use vertical_mix
   use omp_lib

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_horizontal_mix_unified, &
             hdifft_unified

!EOP
!BOC

! !PUBLIC VARIABLES

!  !dir$ attributes offload:mic :: hmix_tracer_itype
!  integer (POP_i4) , public :: hmix_tracer_itype !users choice for type of mixing

!  !dir$ attributes offload:mic :: tavg_HDIFE_TRACER
!  !dir$ attributes offload:mic :: tavg_HDIFN_TRACER
!  !dir$ attributes offload:mic :: tavg_HDIFB_TRACER
!  integer (POP_i4) ,public ,dimension(nt) :: &
!      tavg_HDIFE_TRACER,            &! tavg id for east face diffusive flux of tracer
!      tavg_HDIFN_TRACER,            &! tavg id for north face diffusive flux of tracer
!      tavg_HDIFB_TRACER              ! tavg id for bottom face diffusive flux of tracer

!   !dir$ attributes offload:mic :: tavg_HDIFS
!   !dir$ attributes offload:mic :: tavg_HDIFT
!   integer (POP_i4),public ::            &
!      hmix_momentum_itype,          &! users choice for type of mixing
!      tavg_HDIFT,                   &! tavg id for horizontal diffusion
!      tavg_HDIFS                     ! tavg id for horizontal diffusion


!   !dir$ attributes offload:mic :: lsubmesoscale_mixing
!   logical (log_kind) ,public ::    &
!      lsubmesoscale_mixing           ! if true, submesoscale mixing is on
  
!-----------------------------------------------------------------------
!
!  horizontal mixing choices
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_horizontal_mix
! !INTERFACE:

 subroutine init_horizontal_mix_unified()


 print *,"Intializing unified grids" 

!-----------------------------------------------------------------------
!EOC

 end subroutine init_horizontal_mix_unified

!***********************************************************************
!BOP
! !IROUTINE: hdiffu
! !INTERFACE:

 subroutine hdiffu(k,HDUK,HDVK,UMIXK,VMIXK,this_block)

! !DESCRIPTION:
!  This routine returns tendencies for horizontal diffusion of
!  momentum.  It is a driver routine which simply branches to the
!  proper horizontal mix routine based on the user choice of mixing
!  method.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: k   ! depth level index

   real (POP_r8), dimension(nx_block,ny_block), intent(in) :: &
      UMIXK, VMIXK         ! U,V at level k and mix time level

   type (block), intent(in) :: &
      this_block           ! block information for this subblock

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block), intent(out) :: &
      HDUK,                   &! returned as Hdiff(U) at level k
      HDVK                     ! returned as Hdiff(V) at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  branch to the proper mix routine
!
!-----------------------------------------------------------------------

   call timer_start(timer_hdiffu, block_id=this_block%local_id)

   select case (hmix_momentum_itype)

   case (hmix_momentum_type_del2)
      call hdiffu_del2(k, HDUK, HDVK, UMIXK, VMIXK, this_block)
   case (hmix_momentum_type_del4)
      call hdiffu_del4(k, HDUK, HDVK, UMIXK, VMIXK, this_block)
   case (hmix_momentum_type_anis)
      call hdiffu_aniso(k, HDUK, HDVK, UMIXK, VMIXK, this_block)
   end select

   call timer_stop(timer_hdiffu, block_id=this_block%local_id)

!-----------------------------------------------------------------------
!EOC

 end subroutine hdiffu

!***********************************************************************
!BOP
! !IROUTINE: hdifft
! !INTERFACE:

 !dir$ attributes offload: mic :: hdifft 
 subroutine hdifft(k, HDTK, TMIX, UMIX, VMIX, this_block)

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

   select case (hmix_tracer_itype)
   !case (hmix_tracer_type_del2)
      !print *,"2 is called" 
      !call hdifft_del2(k, HDTK, TMIX, tavg_HDIFE_TRACER, tavg_HDIFN_TRACER, this_block)
   !case (hmix_tracer_type_del4)
      !print *,"4 is called"
      !call hdifft_del4(k, HDTK, TMIX, tavg_HDIFE_TRACER, tavg_HDIFN_TRACER, this_block)
   case (hmix_tracer_type_gm)
      if (k == 1) then
        !start_time = omp_get_wtime() 
        call tracer_diffs_and_isopyc_slopes(TMIX, this_block)
        !end_time = omp_get_wtime()
        !print *,"time at tracer_diffs 1 is ",end_time - start_time   
      endif

      if (k == 1) then
         !start_time = omp_get_wtime()

         call hdifft_gm(1, HDTK_BUF(:,:,:,1), TMIX, UMIX, VMIX, tavg_HDIFE_TRACER, &
                         tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

         !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)num_threads(59) 
         do kk=2,km
         call hdifft_gm(kk , HDTK_BUF(:,:,:,kk) , TMIX, UMIX,VMIX,tavg_HDIFE_TRACER, &
                                 tavg_HDIFN_TRACER,tavg_HDIFB_TRACER,this_block)
         enddo

         !end_time = omp_get_wtime()

         !print *,"time at hdifft_gm combined is ",end_time - start_time 
                    
      endif

      !start_time = omp_get_wtime()  
      HDTK = HDTK_BUF(:,:,:,k)
      !end_time = omp_get_wtime()

      !print *,"time at hdifft_gm is ",end_time - start_time
 
   end select
   
  

   if ( lsubmesoscale_mixing ) then
        if (.not. hmix_tracer_itype == hmix_tracer_type_gm) then
         if (k == 1) then
          !start_time = omp_get_wtime()
          call tracer_diffs_and_isopyc_slopes(TMIX, this_block)
          !end_time = omp_get_wtime()
          !print *,"time at tracer_diffs 2 is ",end_time - start_time
         endif
        endif
        if (k == 1) then
         !start_time = omp_get_wtime()
         call submeso_sf(TMIX, this_block)
         !end_time = omp_get_wtime()
         !print *,"time at submeso_sf is ",end_time - start_time
        endif
        if(k==1) then
         !start_time = omp_get_wtime()
        !!$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)num_threads(60) 
        do kk=1,km
         call submeso_flux(kk, TDTK(:,:,:,kk), TMIX, tavg_HDIFE_TRACER, &
                       tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)
        enddo
        end_time = omp_get_wtime()
        !print *,"time at submeso_flux is ",end_time - start_time
        endif
        HDTK=HDTK+TDTK(:,:,:,k)
   endif
   
  
   
   
!-----------------------------------------------------------------------
!
!  compute tavg diagnostic if requested
!
!-----------------------------------------------------------------------

   !if (accumulate_tavg_now(tavg_HDIFT)) then
     !WORK = c0
     !if (partial_bottom_cells) then
        !where (k <= KMT(:,:,bid)) WORK = DZT(:,:,k,bid)*HDTK(:,:,1)
     !else
        !where (k <= KMT(:,:,bid)) WORK = dz(k)*HDTK(:,:,1)
     !endif
     !call accumulate_tavg_field(WORK,tavg_HDIFT,bid,k)
   !endif

   !if (accumulate_tavg_now(tavg_HDIFS)) then
     !WORK = c0
     !if (partial_bottom_cells) then
        !where (k <= KMT(:,:,bid)) WORK = DZT(:,:,k,bid)*HDTK(:,:,2)
     !else
        !where (k <= KMT(:,:,bid)) WORK = dz(k)*HDTK(:,:,2)
     !endif
   !  call accumulate_tavg_field(WORK,tavg_HDIFS,bid,k)
   !endif

!-----------------------------------------------------------------------
!EOC

 end subroutine hdifft

!***********************************************************************

 subroutine iso_impvmixt_tavg(TNEW, bid)

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TNEW         ! on input, contains tracer to update from
                   ! on output, contains updated tracers at new time

   integer (int_kind), intent(in) :: &
      bid          ! local block address

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block) :: & 
      WORK1, WORK2

   integer (int_kind) ::  &
      k,n                  ! dummy loop indices

!-----------------------------------------------------------------------

   if (hmix_tracer_itype /= hmix_tracer_type_gm) return

   do n = 1,nt
      if (accumulate_tavg_now(tavg_HDIFB_TRACER(n))) then
         do k=1,km-1
            WORK1 = VDC_GM(:,:,k,bid)
            if (partial_bottom_cells) then
               WORK2 = merge(WORK1*(TNEW(:,:,k,n) - TNEW(:,:,k+1,n))/        &
                             (p5*(DZT(:,:,k,bid) + DZT(:,:,k+1,bid)))        &
                             ,c0, k < KMT(:,:,bid))*dzr(k)
            else
               WORK2 = merge(WORK1*(TNEW(:,:,k,n) - TNEW(:,:,k+1,n))*dzwr(k) &
                             ,c0, k < KMT(:,:,bid))*dzr(k)
            endif
            call accumulate_tavg_field(WORK2,tavg_HDIFB_TRACER(n),bid,k)
         end do
      endif
   end do

 end subroutine iso_impvmixt_tavg

!***********************************************************************

 end module horizontal_mix_unified

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
