!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_input_arrays
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
USE mo_kind, ONLY: dp
IMPLICIT NONE
#else
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
#endif

CONTAINS

! Compresses packed data using integer mask. Destination (Dst) must be allocated and large enough to hold all data
! NOTICE: For performance reasons, mask must have even size - terminate with extra 0 if necessary
  SUBROUTINE MaskPack(Src,Dst,Mask,mul,rep)
  REAL(dp)          :: Src(*), Dst(*)
  INTEGER, POINTER :: Mask(:)
  INTEGER, OPTIONAL, INTENT(in) :: mul, rep

  INTEGER :: i, j, dsp, dep, ssp, sep, imul, irep, nx

    imul = 1
    irep = 1
    IF (PRESENT(mul)) imul = mul
    IF (PRESENT(rep)) irep = rep

    dsp = 1
    ssp = 1
    DO j=1,irep
      DO i=1,SIZE(Mask),2
        nx = Mask(i) * imul
        dep = dsp + nx - 1
        sep = ssp + nx - 1
        Dst(dsp:dep) = Src(ssp:sep)
        nx = Mask(i+1) * imul
        ssp = sep + nx + 1
        dsp = dep + 1
      ENDDO
    ENDDO

  END SUBROUTINE MaskPack

! Decompresses packed data using integer mask. Destination (Dst) must be allocated and large enough to hold all data
! NOTICE: For performance reasons, mask must have even size - terminate with extra 0 if necessary
  SUBROUTINE MaskUnpack(Src,Dst,Mask,mul,rep,Fill)
  REAL(dp)         :: Src(*), Dst(*)
  INTEGER, POINTER :: Mask(:)
  INTEGER, OPTIONAL, INTENT(in) :: mul, rep
  REAL(dp), OPTIONAL, INTENT(in) :: Fill

  INTEGER :: i, j, dsp, dep, ssp, sep, imul, irep, nx

    imul = 1
    irep = 1
    IF (PRESENT(mul)) imul = mul
    IF (PRESENT(rep)) irep = rep

    dsp = 1
    ssp = 1
    DO j=1,irep
      DO i=1,SIZE(Mask),2
        nx = Mask(i) * iMul
        dep = dsp + nx - 1
        sep = ssp + nx - 1
        Dst(dsp:dep) = Src(ssp:sep)
        nx = Mask(i+1) * iMul
        dsp = dep + nx + 1
        ssp = sep + 1
      ENDDO
    ENDDO

    IF (PRESENT(Fill)) THEN
      dsp = Mask(1) + 1
      ssp = dsp
      sep = SUM(Mask)
      DO j=1,irep
        dsp = ssp
        DO i=2,SIZE(Mask)-1,2
          nx = Mask(i) * iMul
          dep = dsp + nx - 1
          Dst(dsp:dep) = Fill
          nx = Mask(i+1) * iMul
          dsp = dep + nx + 1
        ENDDO
        IF (MOD(SIZE(Mask),2)==0) THEN
          nx = Mask(SIZE(Mask)) * iMul
          dep = dsp + nx - 1
          Dst(dsp:dep) = Fill
        ENDIF
        ssp = ssp + sep
      ENDDO
    ENDIF

  END SUBROUTINE MaskUnpack

! Creates a Fortran logical mask from an integer mask. Dst must have been allocated with sufficient size in advance
  SUBROUTINE MaskInt2Logical(Dst,Mask,lTrue,lFalse,lReverse,mul,rep)
  LOGICAL, INTENT(inout)        :: Dst(*)
  INTEGER, POINTER              :: Mask(:)
  LOGICAL, OPTIONAL, INTENT(in) :: lTrue, lFalse, lReverse
  INTEGER, OPTIONAL, INTENT(in) :: mul, rep

  INTEGER i, j, dsp, dep, r, m
  LOGICAL lCurr, True, False

    m = 1
    r = 1
    True  = .TRUE.
    False = .TRUE.
    IF (PRESENT(mul   )) m     = mul
    IF (PRESENT(rep   )) r     = rep
    IF (PRESENT(lTrue )) True  = lTrue
    IF (PRESENT(lFalse)) False = lFalse
    dsp = 1
    lCurr = .TRUE.
    IF (PRESENT(lReverse)) lCurr = .NOT. lReverse
    DO j=1,r
      DO i=1,SIZE(Mask)
        dep = dsp + Mask(i) * m
        IF ((True .AND. lCurr) .OR. (False .AND. .NOT. lCurr)) Dst(dsp:dep-1) = lCurr
        lCurr = .NOT. lCurr
        dsp = dep
      ENDDO
    ENDDO

  END SUBROUTINE MaskInt2Logical

! Creates an integer mask from a Fortran logical mask
  SUBROUTINE MaskLogical2Int(Dst,Mask,Sz)
  INTEGER,    POINTER :: Dst(:)
  LOGICAL, INTENT(in) :: Mask(*)
  INTEGER, INTENT(in) :: Sz

  INTEGER, POINTER :: Tmp(:)
  INTEGER i, n, dp
  LOGICAL lCurr

    ALLOCATE(Tmp(Sz+10)) ! 10 arbitrarily choosen - should ensure that the buffer is in any case large enough
    lCurr = .TRUE.
    n  = 0
    dp = 1
    DO i=1,Sz
      IF ((Mask(i) .AND. lCurr) .OR. .NOT. (Mask(i) .OR. lCurr)) THEN
        n = n + 1
      ELSE
        Tmp(dp) = n
        dp = dp + 1
        n  = 1
        lCurr = .NOT. lCurr
      ENDIF
    ENDDO
    IF (n>0) THEN
      Tmp(dp) = n
      dp = dp + 1
    ENDIF
    IF (lCurr) THEN
      Tmp(dp) = 0
      dp = dp + 1
    ENDIF
    ALLOCATE(Dst(dp-1))
    Dst(:) = Tmp(1:dp-1)
    DEALLOCATE(Tmp)

  END SUBROUTINE MaskLogical2Int

! Number of points for an unpacked field
  INTEGER FUNCTION MaskGetTotal(dta)
  INTEGER, INTENT(in) :: dta(:)

    MaskGetTotal = SUM(dta)

  END FUNCTION MaskGetTotal

! Number of points for a packed field
  INTEGER FUNCTION MaskGetValid(dta)
  INTEGER, INTENT(in) :: dta(:)

    MaskGetValid = SUM(dta(1:SIZE(dta):2))

  END FUNCTION MaskGetValid

! Tests if all valid data in mask are contiguous
  LOGICAL FUNCTION MaskIsCont(dta,FirstValid)
  INTEGER,           INTENT(in ) :: dta(:)
  INTEGER, OPTIONAL, INTENT(out) :: FirstValid

    MaskIsCont = .FALSE.
    IF ((SIZE(dta)==2) .OR. ((SIZE(dta)==4) .AND. (dta(1)==0))) THEN
      MaskIsCont = .TRUE.
      IF (PRESENT(FirstValid)) THEN
        IF (dta(1)==0) THEN
          FirstValid = dta(2)+1
        ELSE
          FirstValid = 1
        ENDIF
      ENDIF
    ENDIF

  END FUNCTION MaskIsCont

! Permutes the order of the dimensions of an array by swapping the data elements
  LOGICAL FUNCTION ArrayPermute(ArrSrc,ArrDst,Sz,Order)
  REAL(dp), DIMENSION(*), INTENT(in ) :: ArrSrc
  REAL(dp), DIMENSION(*), INTENT(out) :: ArrDst
  INTEGER,  DIMENSION(:), INTENT(in)  :: Sz,Order

  INTEGER :: SzSrc(7), SzDst(8), Ord(7), Curr(7), StSrc(8), StDst(8)
  INTEGER :: i, j, k, nDim, nDimLo, nDimHi, nBlk, SzBlk, SzRow, nCol, cpd, cps, ofs

    ! Check if order and size are internally consistent
    ArrayPermute = .FALSE.
    nDim = SIZE(Order)
    IF ((SIZE(Sz) /= nDim) .OR. ANY(Order>SIZE(Sz)) .OR. ANY(Order<1) .OR. (nDim>7)) RETURN
    j = 0
    DO i=1,nDim
      IF (IAND(j,ISHFT(1,Order(i))) /= 0) RETURN ! Each of the numbers [1,...,size(Sz)] must appear exactly once in Order
      j = IOR(j,ISHFT(1,Order(i)))
    ENDDO
    ArrayPermute = .TRUE.
    IF (ALL(Order(2:nDim)-Order(1:nDim-1)==1)) THEN ! Actually no reordering requested, just copy data
      SzRow = PRODUCT(Sz)
      ArrDst(1:SzRow) = ArrSrc(1:SzRow)
      RETURN
    ENDIF

    ! Reduce effective number of dimensions by joining adjecent dimensions to be moved together
    SzSrc(1:nDim) = Sz   (1:nDim)
    Ord  (1:nDim) = Order(1:nDim)
    DO WHILE (ANY(Ord(2:nDim)-Ord(1:nDim-1)==1)) ! As long as any neighbouring source dimensions are also adjecent in dst
      i = 1
      DO WHILE (Ord(i+1) /= Ord(i)+1) ! Find the first dimension where this is the case
        i = i + 1
      ENDDO
      ! Update size and order arrays
      SzSrc(Ord(i)) = SzSrc(Ord(i)) * SzSrc(Ord(i+1))
      SzSrc(Ord(i+1)) = 0
      SzSrc(1:nDim-1) = PACK(SzSrc(1:nDim),mask=SzSrc(1:nDim)>0)
      WHERE (Ord(1:nDim) > Ord(i)) Ord(1:nDim) = Ord(1:nDim) - 1
      nDim = nDim - 1
      Ord  (i+1:nDim) = Ord(i+2:nDim)
    ENDDO

    ! Determine number and size of independent blocks
    IF (Order(nDim) == nDim) THEN
      nDimHi = nDim - 1
      nBlk   = SzSrc(nDim)
    ELSE
      nDimHi = nDim
      nBlk   = 1
    ENDIF
    SzBlk = PRODUCT(SzSrc(1:nDimHi))

    ! Determine size of rows to move at once
    IF (Ord(1) == 1) THEN
      SzRow  = SzSrc(1)
      nDimLo = 2
    ELSE
      SzRow  = 1
      nDimLo = 1
    ENDIF
    nCol = PRODUCT(SzSrc(nDimLo:nDimHi))

    ! Calculate size of destination array and stride of data in source dimensions
    DO i=1,nDim
      SzDst(i) = SzSrc(Ord(i))
    ENDDO
    SzDst(nDim+1:8) = HUGE(i) 
    StDst(1) = 1
    DO i=1,nDim
      StDst(i+1) = StDst(i)*SzSrc(i)
    ENDDO
    StDst(nDim+2:8) = 0
    DO i=1,nDim
      StSrc(i) = StDst(Ord(i))
    ENDDO
    StSrc(nDim+1) = StDst(nDim+1)

    ! Real moving work starts here
    ofs = 0
    DO i=1,nBlk ! Loop over independent blocks
      curr(:) = 1
      cpd     = ofs
      cps     = ofs
      DO j=1,nCol ! Individual movements
        ArrDst(cpd+1:cpd+SzRow) = ArrSrc(cps+1:cps+SzRow) ! Move data
        cps = cps + StSrc(nDimLo)
        cpd = cpd + SzRow
        curr(nDimLo) = curr(nDimLo) + 1
        k = nDimLo
        DO WHILE (curr(k) > SzDst(k)) 
          cps     = cps - StSrc(k) * SzDst(k) + StSrc(k+1)
          curr(k) = 1
          k       = k + 1
          curr(k) = curr(k) + 1
        ENDDO
      ENDDO
      ofs = ofs + SzBlk
    ENDDO

  END FUNCTION ArrayPermute

! Copies a scalar or an array to fill an other array (multidimensional version of Fortran intrinsic "spread").
!  - Can operate along multiple dimensions at a time, but the order must be right
!  - The size of a destination dimension must be an integer multiplum of the corresponding source dimension
  LOGICAL FUNCTION ArrayCopy(ArrSrc,ArrDst,SzSrc,SzDst,CopyDims,nCopy)
  REAL(dp), DIMENSION(*), INTENT(in )           :: ArrSrc
  REAL(dp), DIMENSION(*), INTENT(out)           :: ArrDst
  INTEGER,  DIMENSION(:), INTENT(in)            :: SzSrc      
  INTEGER,  DIMENSION(:), INTENT(in), OPTIONAL  :: SzDst
  INTEGER,  DIMENSION(:), INTENT(in), OPTIONAL  :: CopyDims, nCopy

  INTEGER :: i, j, k, n, nBlk, SzBlk, SzS(7), SzD(7), nSDim, nDDim, ofS, ofD
  REAL(dp), POINTER :: Src(:), Dst(:), Swp(:)

    ! Either CopyDims and nCopy or SzDst must be present
    ArrayCopy = .FALSE.
    IF (.NOT. (PRESENT(SzDst) .OR. (PRESENT(CopyDims) .AND. PRESENT(nCopy)))) RETURN

    ! Make shape of 7D src array
    nSDim = SIZE(SzSrc)
    SzS(1:nSDim)   = SzSrc(:)
    SzS(nSDim+1:7) = 1
    ! Determine destination size
    ! CopyDims overrides SzDst if both are given
    IF (PRESENT(CopyDims)) THEN
      SzD(:) = SzS(:)
      DO i=1,size(CopyDims)
        DO j=7,CopyDims(i)+1,-1
          SzD(j) = SzD(j-1)
          SzS(j) = SzS(j-1)
        ENDDO
        SzD(CopyDims(i)) = nCopy(i)
        SzS(CopyDims(i)) = 1
      ENDDO
    ELSE
      nDDim = size(SzDst)
      SzD(1:nDDim)   = SzDst(:)
      SzD(nDDim+1:7) = 1
    ENDIF

    i = PRODUCT(SzD)

    ! Input is scalar, just copy easily
    IF (ALL(SzS==1)) THEN
      ArrayCopy = .TRUE.
      ArrDst(1:i) = ArrSrc(1)
      RETURN
    ENDIF

    k = COUNT(SzS-SzD/=0)
    IF (k==0) THEN ! Source and destination sizes are the same => simple copy
      ArrayCopy = .TRUE.
      ArrDst(1:i) = ArrSrc(1:i)
      RETURN
    ENDIF

    ! Do the necessary copying
    IF (k==1) THEN ! Only one dimension to copy => we can copy directly from ArrSrc to ArrDst which is faster
      DO WHILE (SzS(k) == SzD(k))
        k = k + 1
      ENDDO
      SzBlk = PRODUCT(SzS(1  :k))
      nBlk  = PRODUCT(SzS(k+1:7))
      n     = SzD(k) / SzS(k)
      IF (n*SzS(k) /= SzD(k)) RETURN ! Not an integer number of copies
      ofS = 0
      ofD = 0
      DO j=1,nBlk
        DO i=1,n
          ArrDst(ofD+1:ofD+SzBlk) = ArrSrc(ofS+1:ofS+SzBlk)
          ofD = ofD + SzBlk
        ENDDO
        ofS = ofS + SzBlk
      ENDDO
    ELSE ! Since Fortran don't allow pointers to assumed shaped arrays, we need to copy to tmp-buffers before actual work
      ALLOCATE(Src(i),Dst(i)) ! i still contains the total size of the final destination array, thus room enough
      j = PRODUCT(SzS)
      Dst(1:j) = ArrSrc(1:j) ! Copy to Dst rather than Src, since we swap buffers before doing the real work
      DO k=1,7
        IF (SzS(k) /= SzD(k)) THEN
          ! Swap working buffers
          Swp => Dst
          Dst => Src
          Src => Swp
          ! Do real data copies
          SzBlk = PRODUCT(SzS(1  :k))
          nBlk  = PRODUCT(SzS(k+1:7))
          n     = SzD(k) / SzS(k)
          IF (n*SzS(k) /= SzD(k)) RETURN ! Not an integer number of copies
          ofS = 0
          ofD = 0
          DO j=1,nBlk
            DO i=1,n
              Dst(ofD+1:ofD+SzBlk) = Src(ofS+1:ofS+SzBlk)
              ofD = ofD + SzBlk
            ENDDO
            ofS = ofS + SzBlk
          ENDDO
          SzS(k) = SzD(k) ! New source array now have the same size as the destination array along dimension k
        ENDIF
      ENDDO
      i = PRODUCT(SzD)
      ArrDst(1:i) = Dst(1:i)
      DEALLOCATE(Src,Dst)
    ENDIF

    ArrayCopy = .TRUE.

  END FUNCTION ArrayCopy

! Reverses the elements of an array along a given dimension
  SUBROUTINE ArrayReverse(Arr,Sz,RDim,Dst)
  REAL(dp), DIMENSION(*), INTENT(inout) :: Arr  ! Array data
  INTEGER,  DIMENSION(:), INTENT(in)    :: Sz   ! Shape of array (cannot be obtained from Arr due to 1D-access (*))
  INTEGER,                INTENT(in)    :: RDim ! Number of the dimension to reverse along
  REAL(dp), DIMENSION(*), OPTIONAL, INTENT(out) :: Dst

  INTEGER, ALLOCATABLE, DIMENSION(:) :: Tmp
  INTEGER :: i, j, s1, s2, s3, i_max, j_max, lo_lo, lo_hi, hi_lo, hi_hi

    s1 = PRODUCT(Sz(1:RDim-1)) ! Size of row
    s2 = s1*Sz(RDim)           ! Block size
    s3 = PRODUCT(Sz(RDim+1:))  ! #blocks
    IF (PRESENT(Dst)) THEN
      i_max = Sz(RDim)
    ELSE
      i_max = Sz(RDim)/2
      ALLOCATE(Tmp(s1))          ! Temporary space for one row
    ENDIF
    j_max = s3*s2-1            ! Very last  element in array - 1
    DO j=0,j_max,s2            ! Loop over blocks
      lo_lo = j+1              ! Very first element in this block (and in first row to move)
      hi_hi = j+s2             ! Very last  element in this block (and in first row to move)
      DO i=i_max,1,-1     ! Loop over rows (since we treat one from each end at a time, were done when the middle is reached)
        lo_hi = lo_lo + s1 - 1 ! Last  element of lower row to move
        hi_lo = hi_hi - s1 + 1 ! First element of upper row to move
        ! Swap lower and upper rows
        IF (PRESENT(Dst)) THEN
          Dst(lo_lo:lo_hi) = Arr(hi_lo:hi_hi)
        ELSE
          Tmp(     :     ) = Arr(lo_lo:lo_hi)
          Arr(lo_lo:lo_hi) = Arr(hi_lo:hi_hi)
          Arr(hi_lo:hi_hi) = Tmp(     :     )
        ENDIF
        lo_lo = lo_hi + 1      ! First element in next lower row is the last  in the current + 1 
        hi_hi = hi_lo - 1      ! Last  element in next upper row is the first in the current - 1
      ENDDO
    ENDDO
    IF (.NOT. PRESENT(Dst)) DEALLOCATE(Tmp)

  END SUBROUTINE ArrayReverse

! Parts the array in two groups and swaps them
  SUBROUTINE ArraySwap(Arr,Sz1,Sz2,Rep,Dst)
  REAL(dp), DIMENSION(*), INTENT(inout) :: Arr         ! Array data
  INTEGER,                INTENT(in)    :: Sz1, Sz2    ! Size of the two blocks to swap
  INTEGER,                INTENT(in)    :: Rep         ! Repeat for this many time two blocks (groups)
  REAL(dp), DIMENSION(*), OPTIONAL, INTENT(out) :: Dst

  INTEGER, ALLOCATABLE, DIMENSION(:) :: Tmp
  INTEGER :: i, ofs, Sz

    Sz = Sz1 + Sz2
    IF (.NOT. PRESENT(Dst)) THEN
      ALLOCATE(Tmp(Sz))        ! Temporary space for one group
    ENDIF
    ofs = 0
    DO i=1,Rep                 ! Loop over groups
      IF (PRESENT(Dst)) THEN
        Dst(ofs    +1:ofs+Sz2) = Arr(ofs+Sz1+1:ofs+Sz )
        Dst(ofs+Sz2+1:ofs+Sz ) = Arr(ofs    +1:ofs+Sz1)
      ELSE
        IF (Sz1 > Sz2) THEN ! Prevent unintended overwritings
          Tmp(        1:    Sz1) = Arr(ofs    +1:ofs+Sz1)
          Arr(ofs    +1:ofs+Sz2) = Arr(ofs+Sz1+1:ofs+Sz )
          Arr(ofs+Sz2+1:ofs+Sz ) = Tmp(        1:    Sz1)
        ELSE
          Tmp(        1:    Sz2) = Arr(ofs+Sz1+1:ofs+Sz )
          Arr(ofs+Sz2+1:ofs+Sz ) = Arr(ofs    +1:ofs+Sz1)
          Arr(ofs    +1:ofs+Sz2) = Tmp(        1:    Sz2)
        ENDIF
      ENDIF
      ofs = ofs + Sz
    ENDDO

    IF (.NOT. PRESENT(Dst)) DEALLOCATE(Tmp)

  END SUBROUTINE ArraySwap

END MODULE mo_input_arrays
