!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_nudging_sst
!BOP
  ! !MODULE: mo_nudging_sst (layer 2)

  ! !DESCRIPTION: 
  ! update SST field for nudging

  ! !REVISION HISTORY: 
  ! I. Kirchner, MPI Hamburg, April-2001
  ! R. Johanni, IPP Garching, May-2002, parallel version
  ! I. Kirchner, MPI Hamburg, Aug-2002, revision

  ! !USES:
  USE mo_kind,              ONLY: dp
  USE mo_decomposition,     ONLY: ldc => local_decomposition, &
                                  gdc => global_decomposition

!BOX
  IMPLICIT NONE

  PRIVATE
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: NudgingReadSST
  PUBLIC :: NudgingSSTnew
  PUBLIC :: NudgingSSTClose

  ! !PUBLIC DATA MEMBERS:

  !* namelist parameters
  CHARACTER(len=256), PUBLIC :: ndg_file_sst        ! template name of sst file

  REAL(kind=dp), PUBLIC      :: ndg_freez  = 271.65_dp ! correction of freezing point
                                                    ! with ERA15 (-1.5 C recommended)
  INTEGER, PUBLIC            :: nsstinc = 0         ! hours between new SST fields.
                                                    ! nsstinc = 0 means sea surface temperatures are not used
                                                    ! for nudging. This prevents too low sea ice coverage
                                                    ! in the Arctic when the temperatures used for nudging rise
                                                    ! above the freezing temperature of sea water, but remain
                                                    ! below the melting temperature of (salt-free) water ice.
  INTEGER, PUBLIC            :: nsstoff = 12        ! offset in hours to 00UTC
!EOP

  ! local variables

  REAL(kind=dp), ALLOCATABLE         :: zts(:,:)  ! local buffer
  REAL(kind=dp), ALLOCATABLE, TARGET :: sstn(:,:) ! global buffer

  INTEGER     :: sstblock = -1    ! number of sst file
  INTEGER     :: kfiles
  LOGICAL     :: lsstn            ! true if new sst needed
  INTEGER     :: iheads           ! YYYYMMDDHH sst time step
  INTEGER, SAVE   :: ipos_sst     ! position in sst file

  CHARACTER(len=256) :: mess

CONTAINS

  !======================================================================
!BOP
  !
  ! !IROUTINE:  NudgingReadSST
  ! !INTERFACE:

  SUBROUTINE NudgingReadSST

    ! !DESCRIPTION: 

    ! Update global SST field for nudging experiments. For first
    ! timestep, SST is taken from history files, whereas nudging
    ! fields will have to be read in immediately after model start.
    ! Update sst every time step in time step loop.

    ! !USES:
    !* general modules
    USE mo_mpi,               ONLY: p_parallel_io
    USE mo_time_control,      ONLY: &
         get_date_components, next_date, current_date, &
         out_convert_date, nudg_sst_time, add_date
    USE mo_exception,         ONLY: message, finish
    USE mo_transpose,         ONLY: scatter_gp
    USE mo_filename,          ONLY: str_filter

    !* nudging modules
    USE mo_nudging_constants, ONLY: lnudgcli, lnudg_run
#ifdef LITTLE_ENDIAN
    USE mo_nudging_utils,     ONLY: swap64, cpbread, WORD_LEN, HEAD_LEN
#else
    USE mo_nudging_utils,     ONLY: cpbread, WORD_LEN, HEAD_LEN
#endif
!EOP

    REAL(kind=dp), POINTER :: gl_sstn(:,:) !
    REAL(kind=dp), ALLOCATABLE :: zhbuf(:)
    INTEGER       :: &
         krtim, iday, isec, kret, iheads_need, iheads_new, &
         ihead(HEAD_LEN), iymd, ihms, yr, mo, dy, hr, mi, se
    LOGICAL            :: found
    CHARACTER (len=WORD_LEN) :: yhead(HEAD_LEN)
    CHARACTER(len=256) :: cfile

    INTRINSIC MOD

    IF (.NOT. lnudg_run .OR. nsstinc == 0) RETURN

    NULLIFY(gl_sstn)

    IF (.NOT. ALLOCATED(zts)) THEN
       ALLOCATE(zts(ldc%nglon,ldc%nglat))
       CALL message('',' Nudging SST local memory initialized')
    END IF

    IF (.NOT.ALLOCATED(sstn)) THEN
       ALLOCATE(sstn(ldc%nlon,ldc%nlat))
       CALL message('',' Nudging SST global memory initialized')
    END IF

    IF (p_parallel_io) THEN
       
       CALL get_date_components(next_date,yr,mo,dy,hr,mi,se)

       ! get actually date/time of NSTEP
       nudg_sst_time = current_date
       CALL out_convert_date(nudg_sst_time, iymd, ihms)

       ! find actually sst /date/time
       krtim = nsstoff + 24  ! in hours of the day
       DO
          krtim = krtim - nsstinc
          IF (krtim*10000 <= ihms) EXIT
       END DO
       IF (krtim < 0) THEN
          krtim = krtim + 24
          iday = -1
          isec = 0
          CALL add_date(iday, isec, nudg_sst_time)
          CALL out_convert_date(nudg_sst_time, iymd, ihms)
       END IF

       iheads_need = iymd*100 + krtim              ! YYYYMMDDHH
       IF (lnudgcli) iheads_need = MOD(iheads_need,1000000)    ! remove year

       lsstn = .FALSE.

       IF (sstblock < 0) THEN
          ! call at first time, no initialization before
          sstblock = 1

          cfile = str_filter(ndg_file_sst,yr,mo,dy,hr,mi,se,sstblock)
          WRITE(mess,*) 'use SST file : ',TRIM(cfile)
          CALL message('',mess)
          INQUIRE(file=cfile,exist=found)
          IF (.NOT.found) &
               CALL finish('NudgingReadSST','Nudging SST data file not found.')
          CALL pbopen(kfiles,cfile,'r',kret)
          IF (kret/=0) CALL finish('NudgingReadSST','nudging SST files available?')

          ipos_sst   = 0
          iheads     = -99
          lsstn      = .TRUE.
          
       ELSE IF (iheads < iheads_need) THEN
          ! the last record is older than needed, new sst will be read
          lsstn = .TRUE.

       ENDIF

       IF (lsstn) THEN

          DO
             ! search next sst data record

             CALL cpbread(kfiles,yhead,WORD_LEN*HEAD_LEN,kret)
             IF (kret/=WORD_LEN*HEAD_LEN) THEN
                ! open next SST data file

                CALL NudgingSSTClose
                sstblock = sstblock + 1
                IF ((sstblock > 3) .AND. lnudgcli) THEN
                   sstblock = 1
                ELSE IF (sstblock == 5) THEN
                   CALL finish('NudgingReadSST','Stop reading SST files.')
                END IF
                
                cfile = str_filter(ndg_file_sst,yr,mo,dy,hr,mi,se,sstblock)
                WRITE(mess,*) 'use SST file : ',TRIM(cfile)
                CALL message('',mess)
                INQUIRE(file=cfile,exist=found)
                IF (.NOT.found) &
                     CALL finish('NudgingReadSST','Nudging SST data file not found.')
                CALL pbopen(kfiles,cfile,'r',kret)

                ipos_sst = 0
                CYCLE
             ENDIF

             CALL util_i8toi4(yhead(1), ihead(1), HEAD_LEN)
             iheads_new = ihead(3)*100+ihead(4)
             IF (lnudgcli) iheads_new = MOD(iheads_new,1000000)  ! remove year

             IF (iheads_new == iheads_need) THEN                 ! next SST record found
                iheads = iheads_new
                EXIT

             ELSE IF ( (iheads_new > iheads_need) .AND. (iheads == -99)) THEN
                ! at initialization time the first SST may be not available
                CALL out_convert_date(current_date, iymd, ihms)
                iheads = iymd*100 + ihms/10000
                IF (lnudgcli) iheads = MOD(iheads,1000000)
                EXIT

             ELSE IF (iheads > iheads_need) THEN
                CALL finish('NudgingReadSST','sst date fault')

             END IF

             ! skip data set
             ipos_sst = ipos_sst + HEAD_LEN + ldc%nlon*ldc%nlat
             CALL pbseek(kfiles,ipos_sst*WORD_LEN,0,kret)
             WRITE (mess,*) 'skip record ',ihead
             CALL message('NudgingReadSST',mess)

          ENDDO

          ipos_sst = ipos_sst + HEAD_LEN + ldc%nlon*ldc%nlat
          WRITE (mess,*) 'USE SST record ',ihead
          CALL message('NudgingReadSST',mess)

          ! read new sst field
          ALLOCATE(zhbuf(ldc%nlat*ldc%nlon)); zhbuf(:) = 0.0_dp
          CALL message('NudgingReadSST',' Attention convert SST DATA to IEEE')

          CALL pbread(kfiles,zhbuf(1),ldc%nlon*ldc%nlat*WORD_LEN,kret)
#ifdef LITTLE_ENDIAN
          CALL swap64(zhbuf,ldc%nlon*ldc%nlat)
#endif
          CALL util_cray2ieee(zhbuf,sstn(1,1),ldc%nlon*ldc%nlat)

          DEALLOCATE (zhbuf)

       ENDIF

       gl_sstn => sstn

    ENDIF

    CALL scatter_gp(gl_sstn,zts,gdc)

  END SUBROUTINE NudgingReadSST

  !======================================================================
!BOP
  ! !IROUTINE:  NudgingSSTnew
  ! !INTERFACE:

  SUBROUTINE NudgingSSTnew(krow)

    ! !DESCRIPTION: 

    ! Replace the SST at every timestep by the field provided in the nudging 
    ! forcing files. Since no information about fractional seaice is provided 
    ! in this dataset, seaice=1 is set for SSTs lower or equal to the 
    ! freezing/melting temperature of seaice as it is done in earlier versions 
    ! of the model. Furthermore the new temperature replaces the ice temperature 
    ! tsi, while the water temperature is set to ctfreez (=271.38 K) (the latter 
    ! might not be necessary).
    ! 
    ! For SSTs larger than ctfreez, the water temperature tsw is replaced and 
    ! the seaice fraction is set to zero. The sea ice detection
    ! temperature can be modified using {\it NDG\_FREEZ}.

    ! !USES:
    USE mo_kind,              ONLY: dp
    USE mo_memory_g3b,        ONLY: tsi, tsw, slm, seaice
    USE mo_physc2,            ONLY: ctfreez
    USE mo_nudging_constants, ONLY: lnudg_run
    USE mo_exception,         ONLY: message
!EOP

    INTEGER                   :: krow
    INTEGER                   :: jrow, jn
    REAL(kind=dp), PARAMETER  :: zdt=0.01_dp  ! correct temperature near freezing point

    INTRINSIC MAX, MIN

    jrow = krow

    IF (nsstinc==0 .OR. .NOT. lnudg_run) THEN
      CALL clsst(jrow) ! use the standard sst for nudging
    ELSE               ! update sea surface temperatures at every time step

      IF (jrow == 1) CALL message('NudgingSSTnew','Update sst from nudging data set')
      DO jn = 1,ldc%nglon
        IF (slm(jn,jrow) < 1.0_dp) THEN             ! sea fraction present

          IF (zts(jn,jrow) <= ndg_freez) THEN
            ! below freezing level
            seaice(jn,jrow) = 1._dp
            tsi(jn,jrow)  = MIN(zts(jn,jrow),ctfreez-zdt)
            tsw(jn,jrow)  = ctfreez

          ELSE
            ! above freezing level
            seaice(jn,jrow) = 0._dp
            tsi(jn,jrow)  = ctfreez
            tsw(jn,jrow)  = MAX(zts(jn,jrow),ctfreez+zdt)

          ENDIF

        ENDIF
      ENDDO

    ENDIF

  END SUBROUTINE NudgingSSTnew

  !======================================================================
!BOP
  ! !IROUTINE:  NudgingSSTClose
  ! !INTERFACE:

  SUBROUTINE NudgingSSTClose(lmem)

    ! !DESCRIPTION: 
    ! close nudging sst file and deallocate memory (optional)

    ! !INPUT PARAMETERS: 
    USE mo_mpi,            ONLY: p_parallel_io
    LOGICAL, OPTIONAL, INTENT(in) :: lmem
!EOP

    INTEGER :: kret

    IF (p_parallel_io) CALL pbclose(kfiles,kret)
    IF (lmem) DEALLOCATE(zts,sstn)

  END SUBROUTINE NudgingSSTClose

  !======================================================================

END MODULE mo_nudging_sst
