!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module hmix_gm_submeso_share 

!BOP
! !MODULE: hmix_gm_submeso_share 

! !DESCRIPTION:
!  This module contains routines for computing tracer and density differences for
!  use in the hmix_gm and mix_submeso routines. In addition, isopycnal slopes are
!  computed if necessary.

! !REVISION HISTORY:
!  SVN:$Id: hmix_gm_submeso_share.F90

! !USES:

   use registry
   use blocks
   use kinds_mod
   use grid, only: KMTE, KMT, KMTN
   use constants
   use state_mod
   use time_management
   use domain_size, only: km, nt
   use domain, only: nblocks_clinic
   use omp_lib

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_meso_mixing,   &
             tracer_diffs_and_isopyc_slopes 

!-----------------------------------------------------------------------
!
!  variables to save from one call to next
!
!-----------------------------------------------------------------------
   !dir$ attributes offload:mic :: RX
   !dir$ attributes offload:mic :: RY
   !dir$ attributes offload:mic :: TX
   !dir$ attributes offload:mic :: TY
   !dir$ attributes offload:mic :: TZ
   !dir$ attributes offload:mic :: SLX
   !dir$ attributes offload:mic :: SLY
   !dir$ attributes offload:mic :: RZ_SAVE
   !dir$ attributes offload:mic :: HXY
   !dir$ attributes offload:mic :: HYX 
   
   real (r8), dimension(:,:,:,:,:), allocatable, public :: &
      RX,RY,            &     ! Dx(rho), Dy(rho)
      TX,TY,TZ               ! tracer differences in each direction
   real (r8), dimension(:,:,:,:,:,:), allocatable, public :: &
      SLX, SLY                ! slope of isopycnal sfcs in x,y-direction
   real (r8), dimension(:,:,:,:), allocatable, public :: &
      RZ_SAVE                 ! Dz(rho)
   real (r8), dimension(:,:,:), allocatable, public :: &
      HXY,              &     ! dx/dy for y-z plane
      HYX                     ! dy/dx for x-z plane



!***********************************************************************

 contains

!***********************************************************************

! !IROUTINE: init_meso_mixing
! !INTERFACE:

   subroutine init_meso_mixing(hmix_tracer_itype,hmix_tracer_type_gm)

! !DESCRIPTION:
!  Initializes various submesoscale and GM mixing options and allocates necessary
!  space. Also computes some time-independent arrays.
!

   integer (int_kind) :: &
      iblock,            &  ! block index
      hmix_tracer_itype,hmix_tracer_type_gm
 
!-----------------------------------------------------------------------
!
!  allocate GM and submeso_flux arrays
!
!-----------------------------------------------------------------------    
    
   allocate (HXY (nx_block,ny_block,nblocks_clinic),    &
             HYX (nx_block,ny_block,nblocks_clinic))
   if(hmix_tracer_itype == hmix_tracer_type_gm) then
    allocate (SLX   (nx_block,ny_block,2,2,km,nblocks_clinic),  &
              SLY   (nx_block,ny_block,2,2,km,nblocks_clinic))
   endif
   allocate (TX(nx_block,ny_block,km,nt,nblocks_clinic),  &
             TY(nx_block,ny_block,km,nt,nblocks_clinic),  &
             TZ(nx_block,ny_block,km,nt,nblocks_clinic))

   allocate (RX(nx_block,ny_block,2,km,nblocks_clinic),  &
             RY(nx_block,ny_block,2,km,nblocks_clinic))

   allocate (RZ_SAVE(nx_block,ny_block,km,nblocks_clinic))

 
   HXY      = c0
   HYX      = c0
   SLX      = c0
   SLY      = c0   
   TX       = c0
   TY       = c0
   TZ       = c0
   RX       = c0
   RY       = c0
   RZ_SAVE  = c0
   

!-----------------------------------------------------------------------
!
!  register init_meso_mixing
!
!-----------------------------------------------------------------------

   call register_string ('init_meso_mixing')
   

!-----------------------------------------------------------------------
!
!  initialize various time-independent arrays
!
!-----------------------------------------------------------------------

  do iblock = 1,nblocks_clinic
     
     !*** Hyx = dy/dx for x-z plane

     HYX(:,:,iblock) = HTE(:,:,iblock) / HUS(:,:,iblock)

     !*** Hxy = dx/dy for y-z plane

     HXY(:,:,iblock) = HTN(:,:,iblock) / HUW(:,:,iblock)

  enddo
  
   
!-----------------------------------------------------------------------

   end subroutine init_meso_mixing

!-----------------------------------------------------------------------


! !IROUTINE: tracer_diffs_and_isopyc_slopes
! !INTERFACE:

   !dir$ attributes offload : mic :: tracer_diffs_and_isopyc_slopes
   subroutine tracer_diffs_and_isopyc_slopes (TMIX, this_block)

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

   call register_string ('tracer_diffs_and_isopyc_slopes')

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

        !$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk)NUM_THREADS(60) 
        do kk=1,km

        call state (kk, kk, TMIX(:,:,kk,1), TMIX(:,:,kk,2),  &
                     this_block, DRHODT=DRDT(:,:,kk), DRHODS=DRDS(:,:,kk))

        enddo 

        !end_time = omp_get_wtime()

        !print *,"Time at first part is",end_time - start_time

        kk=1

            !$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i,KMASKE,KMASKN,tempi,tempip1,tempj,tempjp1)num_threads(60)
            do j=1,ny_block
              do i=1,nx_block


!-----------------------------------------------------------------------
!
!     compute RX=Dx(rho) and RY=Dy(rho) for all vertical levels. 
!
!-----------------------------------------------------------------------

                if ( kk <= KMT(i,j,bid) .and. kk <= KMTE(i,j,bid) ) then
                  KMASKE = c1
                else
                  KMASKE = c0
                endif
                if ( kk <= KMT(i,j,bid) .and. kk <= KMTN(i,j,bid) ) then
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
                  TX(i,j,kk,n,bid) = KMASKE  &
                              * (TMIX(i+1,j,kk,n) - TMIX(i,j,kk,n))


                   if(j <= ny_block-1)& 
                       TY(i,j,kk,n,bid) = KMASKN  &
                              * (TMIX(i,j+1,kk,n) - TMIX(i,j,kk,n))
                enddo

                RX(i,j,ieast ,kk,bid) = DRDT(i,j,kk) * TXP(i,j,kn)  &
                                         + DRDS(i,j,kk) * TX(i,j,kk,2,bid) 

                RY(i,j,jnorth,kk,bid) = DRDT(i,j,kk) * TYP(i,j,kn)  &
                                         + DRDS(i,j,kk) * TY(i,j,kk,2,bid) 


                 if(i >= 2) then
                   RX(i,j,iwest,kk,bid) = DRDT(i,j,kk) * TXP(i-1,j,kn)  &
                                     + DRDS(i,j,kk) * TX (i-1,j,kk,2,bid)
                 endif 

                  if(j >= 2)then 
                     RY(i,j,jsouth,kk,bid) = DRDT(i,j,kk) * TYP(i,j-1,kn)  &
                                      + DRDS(i,j,kk) * TY (i,j-1,kk,2,bid)
                  endif 
              enddo
            enddo



       !start_time = omp_get_wtime()
       match = registry_match('init_gm')
!-------------------------------------------------------------------------
!
!
!      when kk is 1 ends
!
!
!-------------------------------------------------------------------------

        do kk=1,km

            if ( kk < km ) then

            !$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(j,i,temp_ksi,temp_ksip1,temp_ksj,temp_ksjp1,kmask,kmaske,kmaskn,temp_ksim1,kmaskeim1) &
            !$OMP PRIVATE(txpim1,txim1,temp_ksjm1,kmasknjm1,typjm1,tyjm1)NUM_THREADS(60)
            do j=1,ny_block
              do i=1,nx_block
                 KMASK = merge(c1, c0, kk < KMT(i,j,bid))


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
              
                 TZ(i,j,kk+1,1,bid) = TMIX(i,j,kk  ,1) - TMIX(i,j,kk+1,1)
                 TZ(i,j,kk+1,2,bid) = TMIX(i,j,kk  ,2) - TMIX(i,j,kk+1,2) 
                 TZP(i,j,ks) = TEMP(i,j,kn) - TEMP(i,j,ks)


!     RZ = Dz(rho) = DRDT*Dz(T) + DRDS*Dz(S)


                 RZ(i,j) = DRDT(i,j,kk) * TZP(i,j,ks) + DRDS(i,j,kk) * TZ (i,j,kk+1,2,bid) 
                 RZ(i,j) = min(RZ(i,j),-eps2)

         
                 !if (match) then 


                    SLX(i,j,ieast ,kbt,kk,bid) = KMASK * RX(i,j,ieast ,kk,bid) / RZ(i,j)
                    SLX(i,j,iwest ,kbt,kk,bid) = KMASK * RX(i,j,iwest ,kk,bid) / RZ(i,j)
                    SLY(i,j,jnorth,kbt,kk,bid) = KMASK * RY(i,j,jnorth,kk,bid) / RZ(i,j)
                    SLY(i,j,jsouth,kbt,kk,bid) = KMASK * RY(i,j,jsouth,kk,bid) / RZ(i,j)


                 !endif

!-----------------------------------------------------------------------
!
!     compute Dx(rho), Dy(rho) at level kk+1
!
!-----------------------------------------------------------------------

                 KMASKE = merge(c1, c0, kk+1 <= KMT(i,j,bid) .and.  &
                                kk+1 <= KMTE(i,j,bid))

                 KMASKN = merge(c1, c0, kk+1 <= KMT(i,j,bid) .and.  &
                                kk+1 <= KMTN(i,j,bid))



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
                  TX(i,j,kk+1,n,bid) = KMASKE  &
                            * (TMIX(i+1,j,kk+1,n) - TMIX(i,j,kk+1,n))

                  if(j <= ny_block-1)&
                  TY(i,j,kk+1,n,bid) = KMASKN  &
                            * (TMIX(i,j+1,kk+1,n) - TMIX(i,j,kk+1,n))
                 enddo

                 RX(i,j,ieast ,kk+1,bid) = DRDT(i,j,kk+1) * TXP(i,j,ks)  &
                                         + DRDS(i,j,kk+1) * TX(i,j,kk+1,2,bid) 
                 RY(i,j,jnorth,kk+1,bid) = DRDT(i,j,kk+1) * TYP(i,j,ks)  &
                                         + DRDS(i,j,kk+1) * TY(i,j,kk+1,2,bid) 

                 if(i >= 2)then

                 temp_ksim1 = max(-c2, TMIX(i-1,j,kk+1,1))

                 kmaskeim1 = merge(c1, c0, kk+1 <= KMT(i-1,j,bid) .and.  &
                                kk+1 <= KMTE(i-1,j,bid))

                 txpim1 = kmaskeim1 * (temp_ksi - temp_ksim1 ) 

                 txim1 = kmaskeim1 * (TMIX(i,j,kk+1,2) - TMIX(i-1,j,kk+1,2))                 
                   
                 RX(i,j,iwest,kk+1,bid) = DRDT(i,j,kk+1) * txpim1  &
                                       + DRDS(i,j,kk+1) *  txim1

                 endif


                 if(j >= 2) then

                 temp_ksjm1 = max(-c2, TMIX(i,j-1,kk+1,1))

                 kmasknjm1 = merge(c1, c0, kk+1 <= KMT(i,j-1,bid) .and.  &
                                kk+1 <= KMTN(i,j-1,bid))

                 typjm1 = kmasknjm1 * (temp_ksj - temp_ksjm1 )

                 tyjm1 = kmasknjm1 * (TMIX(i,j,kk+1,2) - TMIX(i,j-1,kk+1,2))


                 RY(i,j,jsouth,kk+1,bid) = DRDT(i,j,kk+1) * typjm1  &
                                        + DRDS(i,j,kk+1) * tyjm1

                 endif


                 RZ(i,j) = DRDT(i,j,kk+1) * TZP(i,j,ks) + DRDS(i,j,kk+1) * TZ(i,j,kk+1,2,bid) 
                 RZ_SAVE(i,j,kk+1,bid) = min(RZ(i,j),c0)
                 RZ(i,j) = min(RZ(i,j),-eps2)


            !if (match) then

!-----------------------------------------------------------------------
!
!     compute slope of isopycnal surfaces at level kk+1
!
!-----------------------------------------------------------------------

              if ( kk+1 <= KMT(i,j,bid) ) then
                SLX(i,j,ieast, ktp,kk+1,bid) = RX(i,j,ieast ,kk+1,bid) / RZ(i,j)
                SLX(i,j,iwest, ktp,kk+1,bid) = RX(i,j,iwest ,kk+1,bid) / RZ(i,j)
                SLY(i,j,jnorth,ktp,kk+1,bid) = RY(i,j,jnorth,kk+1,bid) / RZ(i,j)
                SLY(i,j,jsouth,ktp,kk+1,bid) = RY(i,j,jsouth,kk+1,bid) / RZ(i,j)
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
           !$OMP PARALLEL DO DEFAULT(SHARED)PRIVATE(kk,n,j,i)collapse(3)num_threads(60)
           do kk=1,km-1
             do j=1,ny_block
                do i=1,nx_block
                   TZ(i,j,kk+1,n,bid) = TMIX(i,j,kk  ,n) - TMIX(i,j,kk+1,n)
                enddo
             enddo
            enddo
           enddo
          !endif
         endif



        !end_time = omp_get_wtime()

        !print *,"Time taken at second time",end_time - start_time

        !if(my_task == master_task)then
        !open(unit=10,file="/home/aketh/ocn_correctness_data/changed.txt",status="unknown",position="append",action="write",form="unformatted")
        !write(10),SLX,SLY,RX,RY,TX,TY
        !close(10)
        !endif

!-----------------------------------------------------------------------
!
     end subroutine tracer_diffs_and_isopyc_slopes
!
!***********************************************************************


 end module hmix_gm_submeso_share

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
