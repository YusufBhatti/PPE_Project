!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_scan_buffer
  !
  ! Authors:
  !  ?                               original source and changes
  !  A. Rhodin,      DWD, June 2002, new subroutine cleanup_scanbuf
  !

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc=>local_decomposition

  IMPLICIT NONE

  !                                        ! set in  > used in

  REAL(dp), ALLOCATABLE, TARGET :: dtm    (:,:,:)  ! ffti    ! dyn
  REAL(dp), ALLOCATABLE, TARGET :: dtl    (:,:,:)  ! ffti    ! dyn
  REAL(dp), ALLOCATABLE, TARGET :: dalpsl (:,:)    ! ffti    ! dyn
  REAL(dp), ALLOCATABLE, TARGET :: dalpsm (:,:)    ! ffti    ! dyn
  REAL(dp), ALLOCATABLE, TARGET :: dm     (:,:,:)  ! scan1sl ! fftd

  REAL(dp), ALLOCATABLE, TARGET :: vo     (:,:,:)  ! ffti    ! dyn,scan1sl,statd,tf2n
  REAL(dp), ALLOCATABLE, TARGET :: d      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE, TARGET :: t      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE, TARGET :: alps   (:,:)    ! ffti
  REAL(dp), ALLOCATABLE, TARGET :: u      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE, TARGET :: dudl   (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE, TARGET :: v      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE, TARGET :: dvdl   (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE, TARGET :: vol    (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE, TARGET :: vom    (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE, TARGET :: rh     (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE :: qte    (:,:,:)
  REAL(dp), ALLOCATABLE :: xlte   (:,:,:)
  REAL(dp), ALLOCATABLE :: xite   (:,:,:)
  REAL(dp), ALLOCATABLE :: xtte   (:,:,:,:)
  REAL(dp), ALLOCATABLE :: tte    (:,:,:)
  !--- Rebekka: included for adding only horizontal adv. ---
  REAL(dp), ALLOCATABLE :: vomh    (:,:,:)
  REAL(dp), ALLOCATABLE :: volh    (:,:,:)
  REAL(dp), ALLOCATABLE :: tteh    (:,:,:)
  REAL(dp), ALLOCATABLE :: qteh    (:,:,:)
  REAL(dp), ALLOCATABLE :: xlteh   (:,:,:)
  REAL(dp), ALLOCATABLE :: xiteh   (:,:,:)
  !--- end included for adding only horizontal adv. ---
  REAL(dp), ALLOCATABLE :: alpste (:,:)
  REAL(dp), ALLOCATABLE, TARGET :: u0     (:,:)              ! fftd
  REAL(dp), ALLOCATABLE, TARGET :: du0    (:,:)              ! fftd
  REAL(dp), ALLOCATABLE, TARGET :: ul     (:,:)              ! fftd
  REAL(dp), ALLOCATABLE, TARGET :: ulz    (:,:)              ! scan1
  REAL(dp), ALLOCATABLE :: alnpr  (:,:,:)
  REAL(dp), ALLOCATABLE :: alpha  (:,:,:)
  REAL(dp), ALLOCATABLE :: vervel (:,:,:)

  REAL(dp), ALLOCATABLE :: ugeo   (:,:,:)          !added geo. winds fxp 07-11
  REAL(dp), ALLOCATABLE :: vgeo   (:,:,:)

  LOGICAL, SAVE, PRIVATE :: lnot_used   = .TRUE.

CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE m_bufscan

    USE mo_tracdef,        ONLY: ntrac

    INTEGER :: ngpblks, nlev, nproma, ngl

    IF (lnot_used) THEN

       ngl     = ldc% nglat
       ngpblks = ldc% ngpblks
       nlev    = ldc% nlev
       nproma  = ldc% nproma
       ! zero for test_scan_buffer
       ALLOCATE (dtm    (nproma,nlev,ngpblks))       ;dtm    = 0.0_dp
       ALLOCATE (dtl    (nproma,nlev,ngpblks))       ;dtl    = 0.0_dp
       ALLOCATE (dalpsl (nproma,ngpblks))            ;dalpsl = 0.0_dp
       ALLOCATE (dalpsm (nproma,ngpblks))            ;dalpsm = 0.0_dp
       ALLOCATE (dm     (nproma,nlev,ngpblks))       ;dm     = 0.0_dp

       ALLOCATE (vo     (nproma,nlev,ngpblks))       ;vo     = 0.0_dp
       ALLOCATE (d      (nproma,nlev,ngpblks))       ;d      = 0.0_dp
       ALLOCATE (t      (nproma,nlev,ngpblks))       ;t      = 0.0_dp
       ALLOCATE (alps   (nproma,ngpblks))            ;alps   = 0.0_dp
       ALLOCATE (u      (nproma,nlev,ngpblks))       ;u      = 0.0_dp
       ALLOCATE (dudl   (nproma,nlev,ngpblks))       ;dudl   = 0.0_dp
       ALLOCATE (v      (nproma,nlev,ngpblks))       ;v      = 0.0_dp
       ALLOCATE (dvdl   (nproma,nlev,ngpblks))       ;dvdl   = 0.0_dp
       ALLOCATE (vol    (nproma,nlev,ngpblks))       ;vol    = 0.0_dp
       ALLOCATE (vom    (nproma,nlev,ngpblks))       ;vom    = 0.0_dp
       ALLOCATE (rh     (nproma,nlev,ngpblks))       ;rh     = 0.0_dp
       ALLOCATE (qte    (nproma,nlev,ngpblks))       ;qte    = 0.0_dp
       ALLOCATE (xlte   (nproma,nlev,ngpblks))       ;xlte   = 0.0_dp
       ALLOCATE (xite   (nproma,nlev,ngpblks))       ;xite   = 0.0_dp
       ALLOCATE (xtte   (nproma,nlev,ntrac,ngpblks)) ;xtte   = 0.0_dp
       ALLOCATE (tte    (nproma,nlev,ngpblks))       ;tte    = 0.0_dp
       !--- Rebekka: included for adding only horizontal adv. ---
       ALLOCATE (vomh   (nproma,nlev,ngpblks))       ;vomh   = 0.0_dp
       ALLOCATE (volh   (nproma,nlev,ngpblks))       ;volh   = 0.0_dp
       ALLOCATE (tteh   (nproma,nlev,ngpblks))       ;tteh   = 0.0_dp
       ALLOCATE (qteh   (nproma,nlev,ngpblks))       ;qteh   = 0.0_dp
       ALLOCATE (xlteh  (nproma,nlev,ngpblks))       ;xlteh  = 0.0_dp
       ALLOCATE (xiteh  (nproma,nlev,ngpblks))       ;xiteh  = 0.0_dp
       !--- end included for adding only horizontal adv. ---
       ALLOCATE (alpste (nproma,ngpblks))            ;alpste = 0.0_dp
       ALLOCATE (u0     (nlev,ngl))                  ;u0     = 0.0_dp
       ALLOCATE (du0    (nlev,ngl))                  ;du0    = 0.0_dp
       ALLOCATE (ul     (nlev,ngl))                  ;ul     = 0.0_dp
       ALLOCATE (ulz    (nlev,ngl))                  ;ulz    = 0.0_dp
       ALLOCATE (alnpr  (nproma,nlev,ngpblks))       ;alnpr  = 0.0_dp
       ALLOCATE (alpha  (nproma,nlev,ngpblks))       ;alpha  = 0.0_dp
       ALLOCATE (vervel (nproma,nlev,ngpblks))       ;vervel = 0.0_dp
 
       ALLOCATE (ugeo   (nproma,nlev,ngpblks))       ;ugeo   = 0.0_dp
       ALLOCATE (vgeo   (nproma,nlev,ngpblks))       ;vgeo   = 0.0_dp

       lnot_used = .FALSE.

    ENDIF

  END SUBROUTINE m_bufscan
  !------------------------------------------------------------------------------
  SUBROUTINE cleanup_scanbuffer
    !------------------------------------
    ! deallocate variables in this module
    !------------------------------------

    IF (.NOT. lnot_used) THEN

       DEALLOCATE (dtm    )
       DEALLOCATE (dtl    )
       DEALLOCATE (dalpsl )
       DEALLOCATE (dalpsm )
       DEALLOCATE (dm     )
       DEALLOCATE (vo     )
       DEALLOCATE (d      )
       DEALLOCATE (t      )
       DEALLOCATE (alps   )
       DEALLOCATE (u      )
       DEALLOCATE (dudl   )
       DEALLOCATE (v      )
       DEALLOCATE (dvdl   )
       DEALLOCATE (vol    )
       DEALLOCATE (vom    )
       DEALLOCATE (rh     )
       DEALLOCATE (qte    )
       DEALLOCATE (xlte   )
       DEALLOCATE (xite   )
       DEALLOCATE (xtte   )
       DEALLOCATE (tte    )
       !--- Rebekka: included for adding only horizontal adv. ---
       DEALLOCATE (vomh   )
       DEALLOCATE (volh   )
       DEALLOCATE (tteh   )
       DEALLOCATE (qteh   )
       DEALLOCATE (xlteh  )
       DEALLOCATE (xiteh  )
       !--- end included for adding only horizontal adv. ---
       DEALLOCATE (alpste )
       DEALLOCATE (u0     )
       DEALLOCATE (du0    )
       DEALLOCATE (ul     )
       DEALLOCATE (ulz    )
       DEALLOCATE (alnpr  )
       DEALLOCATE (alpha  )
       DEALLOCATE (vervel )
       
       DEALLOCATE (ugeo   )
       DEALLOCATE (vgeo   )

       lnot_used   = .TRUE.

    END IF

  END SUBROUTINE cleanup_scanbuffer
  !----------------------------------------------------------------------------

END MODULE mo_scan_buffer
