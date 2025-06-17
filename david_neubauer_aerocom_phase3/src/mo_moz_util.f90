!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_util
!!
!! \brief
!!  Provides a number of general purpose routines for atmospheric chemistry
!!  such as mmr to vmr conversion etc.
!!
!! \author Stacy Walters (NCAR)
!! \author Sebastian Rast (MPI-M)
!! \author Martin G. Schultz (FZ-Juelich)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - SW: original version (pre 2004)
!!  - SR and MGS: original source of chem_utils (2004)
!!  - MGS: cleaned merge with HAMMONIA code and chem_utils (2008-08-28)
!!  - MGS (2013-07-25) : general cleanup according to new MOZ mechanism
!!
!! \limitations
!!  none
!!
!! \details
!!
!! \bibliographic_references
!!  none
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

  MODULE mo_moz_util

  USE mo_kind,   ONLY: dp

  IMPLICIT NONE

  PRIVATE

  SAVE

!  subroutines from original chem_utls module:
  PUBLIC :: DOY,              &  ! day of year
            CALDAYR,          &  ! calendar day
            inti_mr_xform,    &  ! initialize mean air mass (dry mass) 
            mmr2vmr,          &  ! mass_mixing_ratio    -> volume mixing ratio
            vmr2mmr,          &  ! volume mixing ratio  -> mass_mixing_ratio 
            h2o_to_vmr,       &  ! H2O mass mixing ratio (q) to vmr
            negtrc,           &  ! avoid negative tracer concentrations
            get_spc_ndx,      &  ! obtain solution species index in MOZART chemistry
            get_inv_ndx,      &  ! obtain invariant species index in MOZART chemistry
            get_extfrc_ndx,   &  ! obtain index in "external forcing" (=aircraft emissions) list
            get_rxt_ndx          ! obtain index of named reaction 

!   physical constants
  REAL(dp), PARAMETER, PUBLIC      :: Rgas = 8.3144_dp    ! gas constant [ J mole-1 K-1 ]
  REAL(dp), PARAMETER, PUBLIC      :: mwair = 28.966_dp   ! molecular weight of air [ g/mole ]


  CONTAINS

  integer function DOY( mon, cday )
!-----------------------------------------------------------------------
!   ... Compute day of year ignoring leap years.
!           Returns values in the range [1,365].
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: mon, cday

      integer, save :: jdcon(12) = &
            (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)

      doy = jdcon(mon) + cday

  end function DOY


  real(dp) function CALDAYR( idate, isec )
!-----------------------------------------------------------------------
!   ... Calendar day with fractional part.  Returns values in the
!           range [1., 366.)
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
! ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: idate, isec

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------
      integer :: mon, day

      mon = MOD( idate,10000 ) / 100
      day = MOD( idate,100 )

      CALDAYR = DOY( mon, day ) + REAL( isec, dp )/86400._dp

  end function CALDAYR

  subroutine inti_mr_xform( amu, mbar )
!-----------------------------------------------------------------
!     Note: This routine was originally intended to compute the mass of
!           moist air, while we use the dry air mass for the MMR/VMR conversions.
!           In moz2 and moz3, this routine simply copies the dry air
!           mass field from ECHAM into the mbar argument.
!           In HAMMONIA the air mass is either taken from the msis-climatology or
!           computed from the sum of all tracers (the missing sum is assumed to be N2).
!-----------------------------------------------------------------

      use mo_moz_mods, only : plonl, plev

      implicit none

!-----------------------------------------------------------------
!	... dummy args
!-----------------------------------------------------------------
      real(dp), intent(in)  :: amu(plonl,plev)    ! mean atmos. dry mass (28.966_dp)
      real(dp), intent(out) :: mbar(plonl,plev)   ! mean atm mass output ( = amu )

!-----------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------
!     real(dp), parameter :: mfac = 1._dp / .622_dp

      integer :: k

      do k = 1,plev
! wet air mass calculation from original MOZART: no longer used
!         mbar(:,k) = amu(:,k) * (1._dp + sh(:,k)) / (1._dp + sh(:,k)*mfac)
         mbar(:,k) = amu(:,k)
      end do

  end subroutine inti_mr_xform

  subroutine mmr2vmr( vmr, mmr, mbar )
!-----------------------------------------------------------------
!	... xfrom from mass to volume mixing ratio
!-----------------------------------------------------------------

      use mo_moz_mods,   only : plonl, plev, pcnstm1, adv_mass

      implicit none

!-----------------------------------------------------------------
!	... dummy args
!-----------------------------------------------------------------
      real(dp), intent(in)  :: mbar(plonl,plev)
      real(dp), intent(in)  :: mmr(plonl,plev,pcnstm1)
      real(dp), intent(out) :: vmr(plonl,plev,pcnstm1)

!-----------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------
      integer :: k, m

      do m = 1,pcnstm1
         if( adv_mass(m) /= 0._dp ) then
            do k = 1,plev
               vmr(:,k,m) = mbar(:,k) * mmr(:,k,m) / adv_mass(m)
!++hs : baustelle    (limit vmr in mmr2vmr)
! Setting a minimum mixing ratio was necessary to avoid non converging
! chemistry. After some modifications to the model it should be tested if
! it is still necessary.
!!!sschr: ##### next line commented out
!              vmr(:,k,m) = MAX(vmr(:,k,m),1.E-100_dp)
!--hs

            end do
         end if
      end do

  end subroutine mmr2vmr

  subroutine vmr2mmr( vmr, mmr, mbar )
!-----------------------------------------------------------------
!	... xfrom from mass to volume mixing ratio
!-----------------------------------------------------------------

      use mo_moz_mods, only : plonl, plev, pcnstm1,  &
                              adv_mass

      implicit none

!-----------------------------------------------------------------
!	... dummy args
!-----------------------------------------------------------------
      real(dp), intent(in)  :: mbar(plonl,plev)
      real(dp), intent(in)  :: vmr(plonl,plev,pcnstm1)
      real(dp), intent(out) :: mmr(plonl,plev,pcnstm1)

!-----------------------------------------------------------------
!	... local variables
!-----------------------------------------------------------------
      integer :: k, m

!-----------------------------------------------------------------
!	... the non-group species
!-----------------------------------------------------------------
!++hs: some diagnostic print statements from HAMMONIA ... (removed here)
      do m = 1,pcnstm1
         if( adv_mass(m) /= 0._dp ) then
            do k = 1,plev
               mmr(:,k,m) = adv_mass(m) * vmr(:,k,m) / mbar(:,k)
            end do
         end if
      end do

  end subroutine vmr2mmr

!----------------------------------------------------------------------
! H2O conversion from mmr to vmr
! This routine also makes sure that H2O is positive and not unreasonably small
!----------------------------------------------------------------------

  subroutine h2o_to_vmr ( h2o_vmr, h2o_mmr, mbar, plonl )

      use mo_kind,     only : dp
      use mo_moz_mods, only : plev

      implicit none

!-----------------------------------------------------------------------
! ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  plonl                  ! long tile dim
      real(dp), dimension(plonl,plev), intent(out) :: h2o_vmr   ! water vapor vmr
      real(dp), dimension(plonl,plev), intent(in)  :: h2o_mmr   ! specific humidity ( mmr )
      real(dp), dimension(plonl,plev), intent(in)  :: mbar      ! atmos mean mass

!-----------------------------------------------------------------------
! ... local variables
!-----------------------------------------------------------------------
      real, parameter :: mh2oi  = 1._dp/18.01528_dp

      integer ::   k

      do k = 1,plev
         h2o_vmr(:,k) = mbar(:,k) * h2o_mmr(:,k) * mh2oi
         where(h2o_vmr(:,k) < 1.e-15_dp)
            h2o_vmr(:,k) = 0._dp
         end where
      end do

  end subroutine h2o_to_vmr

  subroutine negtrc( lat, header, fld )
!-----------------------------------------------------------------------
!  	... check for negative constituent values and
!	    replace with zero value
!-----------------------------------------------------------------------

      use mo_exception,   only : message, message_text, em_warn, em_info
      use mo_moz_mods,    only : plonl, plev, pcnstm1
      use mo_moz,         only : tracnam
      use mo_moz,         only : pdiags
      implicit none

!-----------------------------------------------------------------------
!  	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)          :: lat                      ! current latitude
      character(len=*), intent(in) :: header                   ! caller tag
      real(dp), intent(inout)      :: fld(plonl,plev,pcnstm1)  ! field to check

!-----------------------------------------------------------------------
!  	... local variables
!-----------------------------------------------------------------------
      integer :: m
      integer :: nneg                       ! flag counter
      integer :: iw, kw
      integer :: windex(2)
      real(dp)    :: worst

      do m  = 1,pcnstm1
         nneg = count( fld(:,:,m) < 0._dp )
         if( nneg > 0 ) then
            where( fld(:,:,m) < 0._dp )
               fld(:,:,m) = 0._dp
            endwhere
            if( pdiags%negtrc ) then
               worst     = minval( fld(:,:,m) )
               windex(:) = minloc( fld(:,:,m) )
               iw        = windex(1)
               kw        = windex(2)
               write(message_text,'(a,a,a,i0,a)') header(:len(header)), tracnam(m), ' has ',nneg,' neg values'
               CALL message('negtrc (MOZ)', message_text, level=em_warn)
               write(message_text,'(a,e25.15,a,i0,a,i0,a,i0)') ' worst =',worst,' @ long = ',iw,' lat = ',lat,' eta = ',kw
               CALL message('            ', message_text, level=em_info)
            end if
         end if
      end do

  end subroutine negtrc

  integer function get_spc_ndx( spc_name )
!-----------------------------------------------------------------------
!     ... return MOZART species index associated with spc_name
!     ATTENTION: this index is not necessarily equal to the ECHAM species
!     index (see get_tracer). Only a subset of the xt tracer field is 
!     passed to the chemistry driver CHEMDR and to the deposition routine
!     called from vdiff. 
!-----------------------------------------------------------------------

      use mo_moz_mods,  only : pcnstm1
      use mo_moz,       only : tracnam

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_spc_ndx = -1
      do m = 1,pcnstm1
         if( trim( spc_name ) == trim( tracnam(m) ) ) then
            get_spc_ndx = m
            exit
         end if
      end do

  end function get_spc_ndx


!++mgs : This routine is apparently obsolete. 
!     integer function get_het_ndx( het_name )
!-----------------------------------------------------------------------
!     ... return overall het process index associated with spc_name
!-----------------------------------------------------------------------
!
!     use mo_moz_mods,  only : hetcnt, het_lst
!
!     implicit none
!
!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
!     character(len=*), intent(in) :: het_name
!
!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
!     integer :: m
!
!     get_het_ndx = -1
!     do m = 1,max(1,hetcnt)
!        if( trim( het_name ) == trim( het_lst(m) ) ) then
!           get_het_ndx = m
!           exit
!        end if
!     end do
!
!     end function get_het_ndx
!--mgs

  integer function get_extfrc_ndx( frc_name )
!-----------------------------------------------------------------------
!     ... return overall external frcing index associated with spc_name
!-----------------------------------------------------------------------

      use mo_moz_mods,  only : extcnt, extfrc_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: frc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_extfrc_ndx = -1
      if( extcnt > 0 ) then
         do m = 1,max(1,extcnt)
            if( trim( frc_name ) == trim( extfrc_lst(m) ) ) then
               get_extfrc_ndx = m
               exit
            end if
         end do
      end if

  end function get_extfrc_ndx


  integer function get_rxt_ndx( rxt_tag )
!-----------------------------------------------------------------------
!     ... return overall external frcing index associated with spc_name
!-----------------------------------------------------------------------

      use mo_moz_mods,  only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: rxt_tag

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_rxt_ndx = -1
      do m = 1,rxt_tag_cnt
         if( trim( rxt_tag ) == trim( rxt_tag_lst(m) ) ) then
            get_rxt_ndx = rxt_tag_map(m)
            exit
         end if
      end do

  end function get_rxt_ndx


  integer function get_inv_ndx( invariant )
!-----------------------------------------------------------------------
!     ... return index of invariant species
!-----------------------------------------------------------------------

      use mo_moz_mods,  only : nfs, inv_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: invariant

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_inv_ndx = -1
      do m = 1,nfs
         if( trim( invariant ) == trim( inv_lst(m) ) ) then
            get_inv_ndx = m
            exit
         end if
      end do

  end function get_inv_ndx


  END MODULE mo_moz_util
