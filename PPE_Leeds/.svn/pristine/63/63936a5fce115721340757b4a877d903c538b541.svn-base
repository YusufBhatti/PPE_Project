!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE m_alloc_mods
  ! ------------------------------------------------------------------
  !
  ! module *mo_control* - control variables for model housekeeping.
  !
  ! u. schlese   dkrz-hamburg    dec-94
  !
  ! A. Rhodin    MPI-Hamburg     Jan-99:
  !   Subroutine m_control renamed to alloc_mods and moved from 
  !   module mo_control to module m_alloc_mods.
  !
  ! M.A. Giorgetta, MPI, May 2000
  !   alloc_mods allocates arrays for vertical columns in RRTM modules
  !   and diag in mo_tmp_buffer
  !
  ! L. Kornblueh, MPI, October 2003,
  !   removed diag from mo_tmp_buffer
  ! ------------------------------------------------------------------

  USE mo_control,     ONLY: nlev, nlevp1, nhgl
  USE mo_physc2,      ONLY: cevapcu
  USE mo_hdiff,       ONLY: diftcor, iq, ncdif
  USE mo_hyb,         ONLY: aktlrd, alpham, altrcp, ardprc, bb, ceta,      &
                            cetah, cpg, dela, delb, delpr, ralpha, ralphr, &
                            rddelb, rdelpr, rlnmar, rlnpr
  USE mo_semi_impl,   ONLY: vmax

  IMPLICIT NONE

  LOGICAL, SAVE :: lnot_used = .TRUE.

CONTAINS

  SUBROUTINE alloc_mods

    IMPLICIT NONE

    IF (lnot_used) THEN

       ! mo_physc2
       ALLOCATE (cevapcu(nlev))
       ! mo_hdiff
       ALLOCATE (diftcor(nlev))
       ! mo_hyp
       ALLOCATE (ralpha(nlev))
       ALLOCATE (rlnpr(nlev))
       ALLOCATE (dela(nlev))
       ALLOCATE (delb(nlev))
       ALLOCATE (rddelb(nlev))
       ALLOCATE (cpg(nlev))
       ALLOCATE (delpr(nlev))
       ALLOCATE (rdelpr(nlev))
       ALLOCATE (ralphr(nlev))
       ALLOCATE (alpham(nlev))
       ALLOCATE (ardprc(nlev))
       ALLOCATE (rlnmar(nlev))
       ALLOCATE (aktlrd(nlev))
       ALLOCATE (altrcp(nlev))
       ALLOCATE (ceta(nlev))
       ALLOCATE (cetah(nlevp1))
       ALLOCATE (bb(nlev,nlev))
       ! mo_semi_impl
       ALLOCATE (vmax(nlev))
       ! mo_diff
       ALLOCATE (ncdif(nlev))
       ALLOCATE (iq(nlev))

       lnot_used = .FALSE.

    ENDIF

  END SUBROUTINE alloc_mods
!------------------------------------------------------------------------------
  SUBROUTINE dealloc_mods
    !
    ! deallocate module variables
    !
    DEALLOCATE (cevapcu)
    DEALLOCATE (diftcor)
    DEALLOCATE (ralpha)
    DEALLOCATE (rlnpr)
    DEALLOCATE (dela)
    DEALLOCATE (delb)
    DEALLOCATE (rddelb)
    DEALLOCATE (cpg)
    DEALLOCATE (delpr)
    DEALLOCATE (rdelpr)
    DEALLOCATE (ralphr)
    DEALLOCATE (alpham)
    DEALLOCATE (ardprc)
    DEALLOCATE (rlnmar)
    DEALLOCATE (aktlrd)
    DEALLOCATE (altrcp)
    DEALLOCATE (ceta)
    DEALLOCATE (cetah)
    DEALLOCATE (bb)
    DEALLOCATE (vmax)
    DEALLOCATE (ncdif)
    DEALLOCATE (iq)

    lnot_used = .TRUE.

  END SUBROUTINE dealloc_mods
!------------------------------------------------------------------------------
END MODULE m_alloc_mods
