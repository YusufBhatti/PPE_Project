!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_setinv
!!
!! \brief
!!  This routine sets the chemically inert tracer concentrations of tracers (N2, O2, ...)
!!  This routine vertically redistributes condensed phase HNO3.
!!
!! \author Stacy Walters (NCAR)
!! \author Martin G. Schultz (FZ-Juelich)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - SW: original version (pre 2004)
!!  - MGS: merge of MOZ2 and MOZ3 code (2008-08-09)
!!
!! \limitations
!!  none
!!
!! \details
!!
!! \bibliographic_references
!!  - D.E. Kinnison et al., JGR 112, D20302, doi:10.1029/2006JD007879, 2007
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module mo_moz_setinv
!
! Note: invariants(:,:,mindex) is the air density (in molec. cm-3)
!
!  NOTE: should make use of inv_lst (1:nfs) defined in mo_moz_mods; species set inset_sim_dat

      save

      integer :: m_ndx, o2_ndx, n2_ndx, n2d_ndx
      integer :: stat_o2, stat_n2, stat_n2d                    ! 0 = not present
                                                               ! 1 = solution species
                                                               ! 2 = invariant
      integer, parameter :: max_invariants = 10


      private
      public :: invariants_inti
      public :: setinv


      contains

      subroutine invariants_inti

      use mo_moz_util,  ONLY : get_inv_ndx, get_spc_ndx
      use mo_exception, ONLY : finish, message, message_text, em_info

      implicit none

      stat_n2  = 0
      stat_n2d = 0
      stat_o2  = 0

      ! invariants
      m_ndx   = get_inv_ndx( 'M' )
      n2_ndx  = get_inv_ndx( 'N2' )
      n2d_ndx = get_inv_ndx( 'N2D' )
      o2_ndx  = get_inv_ndx( 'O2' )
      ! note: old MOZ2 code (troposphere only) also had H2O as invariant

      IF (m_ndx <= 0) THEN
          write(message_text,'(a)') 'Invalid chemistry code: Must have M as invariant!'
          CALL finish('init_invariants', message_text)
      END IF
      IF (n2_ndx > 0)  stat_n2  = 2
      IF (n2d_ndx > 0) stat_n2d = 2
      IF (o2_ndx > 0)  stat_o2  = 2

      ! other species needed for invariants calculation (WACCM/HAMMONIA)
      IF (stat_o2 == 0) THEN
          o2_ndx  = get_spc_ndx( 'O2' )
          IF (o2_ndx > 0)  stat_o2  = 1     ! O2 is a solution species
      END IF

      IF (stat_o2 == 2) THEN
        write(message_text,'(a,4i4)') 'inv_ndx(M, N2, N2D, O2) = ', m_ndx, n2_ndx, n2d_ndx, o2_ndx
        CALL message('init_invariants', message_text, level=em_info)
      ELSE
        write(message_text,'(a,3i4)') 'invariants: inv_ndx(M, N2, N2D) = ', m_ndx, n2_ndx, n2d_ndx
        CALL message('init_invariants', message_text, level=em_info)
        write(message_text,'(a,i4)') 'O2 is a solution species: spc_ndx(O2) = ', o2_ndx
        CALL message('init_invariants', message_text, level=em_info)
      END IF

      end subroutine invariants_inti


      subroutine setinv( invariants, tfld, h2ovmr, pmid, vmr )
!-----------------------------------------------------------------
!        ... set the invariant densities (molecules/cm**3)
!-----------------------------------------------------------------
      
      use mo_kind,       only : dp
      use mo_moz_mods,   only : plonl, plev, nfs, pcnstm1
      use mo_moz,        only : lhammoniagases
      use mo_exception,  only : message, message_text, em_debug
      implicit none

!-----------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------
      real(dp), intent(in)  ::      tfld(plonl,plev)           ! temperature
      real(dp), intent(in)  ::      h2ovmr(plonl,plev)         ! water vapor vmr
      real(dp), intent(in)  ::      pmid(plonl,plev)           ! pressure
      real(dp), intent(in)  ::      vmr(plonl,plev,pcnstm1)    ! volume mixing ratios
      real(dp), intent(out) ::      invariants(plonl,plev,max_invariants) ! invariant array
      
!-----------------------------------------------------------------
!        .. local variables
!-----------------------------------------------------------------
      real(dp), parameter ::  boltz = 1.38044e-16_dp   ! erg/k
      real(dp)            ::  n2m1(plonl,plev)         ! 1.-vmr(N2)
      integer :: k, ji, jt

!-----------------------------------------------------------------
!        note: invariants are in cgs density units.
!              the pmid array is in pascals and must be
!	       mutiplied by 10. to yield dynes/cm**2.
!-----------------------------------------------------------------

      invariants(:,:,:) = 0._dp

!-----------------------------------------------------------------
!       ... HAMMONIA case: Calculate mixing ratio of N2
!       Sum up volume mixing ratios of all other compounds and subtract 1.
!       This requires that O2, CO2, N2O and CH4 are contained in the chemical scheme
!-----------------------------------------------------------------
      if ( lhammoniagases ) then 
         n2m1(:,:) = 0.0_dp
         do k = 1,plev
!CDIR unroll=8
           do jt=1,pcnstm1
             n2m1(:,k)=n2m1(:,k) + vmr(:,k,jt)
           end do
         enddo
      else
         n2m1(:,:) = 0.205_dp         ! mole fraction of everything but N2
      end if
    
!-----------------------------------------------------------------
! ... set m, n2, and n2d densities
!-----------------------------------------------------------------
!        note: invariants are in cgs density units.
!              the pmid array is in pascals and must be
!        mutiplied by 10. to yield dynes/cm**2.
!-----------------------------------------------------------------
!+++sschr: commented out
!     WRITE(message_text,'(a,i0)')  'm_ndx = ',m_ndx
!     CALL message('setinv', message_text, level=em_debug)
!---sschr: commented out

      do k = 1,plev
         invariants(:,k,m_ndx) = 10._dp * pmid(:,k) / (boltz*tfld(:,k))
      end do

      if( stat_n2d == 2 ) then       !! N2D zero for now (## check for HAMMONIA ##)
!!       if( .not. is_ecmwf_sim ) then
!!          call invariants_interp( lat, ip, calday, invariants(1,1,n2d_ndx), pmid, plonl )
!!          do k = 1,plev
!!             invariants(:,k,n2d_ndx) = reduction_factor * invariants(:,k,n2d_ndx) * invariants(:,k,m_nd
!!          end do
      end if

      if( stat_n2 == 2 ) then
         do k = 1,plev
            invariants(:,k,n2_ndx) = (1._dp-n2m1(:,k)) * invariants(:,k,m_ndx)      ! N2
         end do
      end if

      if( stat_o2 == 2 ) then
         do k = 1,plev
            invariants(:,k,o2_ndx) = 0.205_dp * invariants(:,k,m_ndx)      ! O2
         end do
      end if

      end subroutine setinv

      end module mo_moz_setinv
