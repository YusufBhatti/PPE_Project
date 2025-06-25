!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_nudging_init
!BOP
  ! !MODULE: mo_nudging_init (layer 1)

  ! !DESCRIPTION: 
  ! initialisation, function control and clean up of nudging procedure 

  ! !REVISION HISTORY: 
  ! I. Kirchner, MPI Hamburg, April-2001
  ! R. Johanni, IPP Garching, May-2002, parallel version
  ! I. Kirchner, MPI Hamburg, Aug-2002, revision
  ! J. S. Rast, MPI Hamburg, Aug-2008, bugfix: added nsstinc condition
  ! D. O'Donnell, MPI-M Hamburg, Feb 2009, bugfix. List-directed internal
  !               I/O to variable 'mess' caused attempt to write past end of 
  !               internal file. Replaced with formatted I/O.
  !               Corrections marked with '!>>dod' in the code.
  ! A. Voigt, MPI-M, Hamburg, March 2011, add ndg_file_nc
  !           ndg_file_nc: name of the netcdf-file containing the nudging data
  ! A. Voigt, MPI Hamburg, March 2011, add inudgformat to chose between CRAY 
  !           binary and netcdf format
  !                   0 --> CRAY binary format (default)
  !                   1 --> GRIB (not implemented)
  !                   2 --> netcdf format               
   
!BOX
  USE mo_nudging_constants, ONLY: inudgformat
         
  USE mo_nudging_io,   ONLY: ndg_file_nc, ndgfile

  IMPLICIT NONE

  PRIVATE
!EOX

  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC  :: NudgingInit

  ! !PUBLIC DATA MEMBERS:
  INTEGER, PARAMETER, PUBLIC ::   &!--- function switch for initialization
       NDG_INI_IO     = 0, &! initialize namelist parameters
       NDG_INI_STREAM = 1, &! initialize stream
       NDG_INI_MEM    = 2, &! allocate memory
       NDG_CLEAN_MEM  = 3, &! deallocate memory
       NDG_CLOSE      = 4   ! close nudging
!BOX
  CHARACTER(len=256) :: mess
!EOX
!EOP
CONTAINS

  !======================================================================
!BOP
  !
  ! !IROUTINE:  NudgingInit
  ! !INTERFACE:

  SUBROUTINE NudgingInit(ifunction)

    ! !DESCRIPTION: 
    ! setup nudging parameter and fields

    ! !USES:

    !* general modules
    USE mo_kind,          ONLY: dp
    USE mo_exception,     ONLY: message, finish
    USE mo_control,       ONLY: &
         lnudge, n2sp, nlev, nm, nlevp1, lnmi, nsp
    USE mo_filename,      ONLY: str_filter
    USE mo_physc2,        ONLY: ctfreez
    USE mo_namelist,      ONLY: open_nml, position_nml, POSITIONED
    USE mo_time_conversion, ONLY: OPERATOR(<)
    USE mo_time_control, ONLY: &
         dt_nudg_start, dt_nudg_stop, lstart, lresume, &
         init_nudgingtime_a, &
         init_nudgingtime_b, &
         init_nudgingtime_c, &
         lfirst_cycle, get_date_components, inp_convert_date, &
         nudg_start, next_date, nudg_stop, ndg_inp_date, &
         delta_time, ndg_date0, ndg_date1

    !* MPI functions
    USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_io, p_bcast
    USE mo_truncation,    ONLY: nmp, nnp
    USE mo_decomposition, ONLY: &
         ldc => local_decomposition, gdc => global_decomposition
    USE mo_transpose,     ONLY: scatter_sp

    !* nudging modules
    USE mo_nudging_constants, ONLY: &
         ndg_version, &
         lnudgdbx, lnudgini, lnudgimp, lnudgfrd, lnudgwobs, ldamplin, &
         ltintlin, nudgdamp, nudgdsize, lnudgpat, lnudgcli, &
         nudgd, nudgv, nudgt, nudgp, &
         nudgtrun, nudgsmin, nudgsmax, nudglmin, nudglmax, &
         lnudgstop, nmaxc, nfact, lsite, &
         linp_vor, linp_div, linp_tem, linp_lnp, &
         lnudg_run, nudgmin, nudgmax, &
         ilev_v_max, ilev_d_max, ilev_t_max, &
         ilev_v_min, ilev_d_min, ilev_t_min, &
         ino_v_lev,  ino_d_lev,  ino_t_lev, &
         NDG_WINDOW_ALL, NDG_WINDOW_CUT, NDG_WINDOW_CUT0 

    USE mo_nudging_sst, ONLY: &
         nsstinc, nsstoff, ndg_freez, ndg_file_sst, NudgingSSTClose
  
    USE mo_nudging_io, ONLY: &
         ndg_file_vor, ndg_file_div, ndg_file_stp, ndgblock, &
         lndgstep0, lndgstep1, lndgstep2, lndgstep3, CloseBlock

    USE mo_nudging_buffer, ONLY: &
         flagn, nudgdo, nudgvo, nudgto, nudgpo, &
         NudgingAllocMem, NudgingDeallocMem, &
         NudgingAllocRefMem, NudgingDeallocRefMem, &
         NudgingStreamInit, nudg, NdgInitCounter

    USE mo_nudging_utils, ONLY: cpbread, HEAD_LEN, WORD_LEN

    USE mo_memory_base,   ONLY: set_stream
    
    USE mo_io,            ONLY: io_read, io_open, IO_close                            !! AV
    USE mo_netCDF,        ONLY: IO_inq_dimid, IO_INQ_DIMLEN, &                        !! AV
                                IO_INQ_VARID,IO_GET_VAR_DOUBLE                        !! AV
    USE mo_filename,      ONLY: NETCDF                                                !! AV			

!EOP
!BOC
!BOX
    INTEGER, INTENT(in)  :: ifunction    ! function selector

    INTEGER  :: &
         i, kret, ibytes, snsp, &
         ihead(8), is1, is2, jm, ii, sum1, sum2, &
         yr, mo, dy, hr, mi, se, &
         ipos ! error return value from position_nml
    REAL(dp)      :: twodt
    LOGICAL            :: found
    CHARACTER(len=8)   :: yhead(HEAD_LEN)
    CHARACTER(len=256) :: cfile
    INTEGER            :: kfile
    INTEGER            :: inml, iunit   
    REAL(dp), POINTER :: flagn_global(:,:,:), flagn_local(:,:,:)
    
    INTEGER                :: ndimid, nts, nvarid               !! AV
    REAL(dp), ALLOCATABLE  :: timevals(:)                       !! AV
    INTEGER                :: ihead_nc(2)                       !! AV: header of the netcdf-file with nudging data
    

    INTRINSIC MAX, MIN

    INCLUDE 'ndgctl.inc'

    snsp = ldc%snsp

    IF (.NOT. lnudge) RETURN

    SELECT CASE(ifunction)

! ==============================================================================
!
!EOX
! ***** initialise constant parameters **** PART 1 nudging init

    CASE(NDG_INI_IO) !==========================================================
!BOX

      !*** read namelist *******************************************************

      IF (p_parallel_io) THEN
        inml = open_nml('namelist.echam')
        iunit = position_nml ('NDGCTL', inml, status=ipos)
        SELECT CASE (ipos)
        CASE (POSITIONED);  READ (iunit,ndgctl)
        END SELECT
      ENDIF

      IF (p_parallel) THEN
        CALL bc_nudging
        IF (lnudgpat) CALL finish('NudgingInit',&
             'Pattern nudging not allowed in parallel mode.')
      END IF

      CALL message('',ndg_version)



!EOX
      !*** define data channel names *******************************************
!BOX
      SELECT CASE(inudgformat)

        CASE(0)     ! nudging data in cray format
 
          CALL message(' File name template VOR     ',ndg_file_vor)
          CALL message(' File name template DIV     ',ndg_file_div)
          CALL message(' File name template STP     ',ndg_file_stp)

        CASE(2)     ! nudging data in netcdf-format
        
          CALL message(' File name template NetCDF  ',ndg_file_nc) 
           
        CASE DEFAULT

          CALL finish('NudgingInit','inudgformat set incorrectly')
    
      END SELECT

      CALL message('','')



!EOX
      !*** define SST dependent parameters *************************************

      !--- nsstinc equals 0 means no external sst is used
!BOX
      nsstinc = MAX(MIN(nsstinc,24),0)
      nsstoff = MAX(MIN(nsstoff,23),0)
      IF (ndg_freez < 0) THEN
         CALL message('Warning',&
              'Predefined NDG_FREEZ not correct ... reset to CTFREEZ')
         ndg_freez = ctfreez
      END IF

      !--- initialize sst memory

      IF (nsstinc==0) THEN
         CALL message('',' Use standard ECHAM SST')
      ELSE IF (lfirst_cycle) THEN
         !>>dod bugfix
         CALL message(' File name template SST : ',TRIM(ndg_file_sst))
         WRITE(mess,'(f7.2)') ndg_freez
         CALL message(' Use sea ice detection limit NDG_FREEZ = ',mess)
         !<<dod
      ENDIF



!EOX
      !*** general nudging options *********************************************
!BOX
      CALL message('','Nudging options:')

      IF (lnudgpat) THEN
         lnudgimp = .FALSE.
         CALL message('','Force explicit method.')
      END IF
      IF (lstart) THEN
         lnudgini = .TRUE.
         CALL message('','Force time synchronization with nudging data.')
      END IF
      !>>dod bugfix
      WRITE(mess,'(l1)') lnudgpat
      CALL message('  pattern assimilation            LNUDGPAT  = ',TRIM(mess))
      WRITE(mess,'(l1)') lnudgini
      CALL message('  Adjust date/time                LNUDGINI  = ',TRIM(mess))
      WRITE(mess,'(l1)') lnudgimp
      CALL message('  Implicit method                 LNUDGIMP  = ',TRIM(mess))
      WRITE(mess,'(l1)') lnudgdbx
      CALL message('  Debbugging                      LNUDGDBX  = ',TRIM(mess))
      WRITE(mess,'(l1)') lnudgwobs
      CALL message('  Store reference data in output  LNUDGWOBS = ',TRIM(mess))
      WRITE(mess,'(l1)') lsite
      CALL message('  Calculate SITE measure          LSITE     = ',TRIM(mess))
      WRITE(mess,'(l1)') lnudgcli
      CALL message('  Use nudging data cyclic         LNUDGCLI  = ',TRIM(mess))
      WRITE(mess,'(l1)') ltintlin
      CALL message('  Linear time interpolation       LTINTLIN  = ',TRIM(mess))

      CALL message('','   Time window options ....')
      WRITE(mess,'(l1)') ldamplin
      CALL message('    Linear damping in time        LDAMPLIN  = ',TRIM(mess))
      nudgdamp = MIN(MAX(nudgdamp,0.0_dp),1.0_dp)
      WRITE(mess,'(f8.6)') nudgdamp
      CALL message('    Reduce damping between reference times (fraction) ',TRIM(mess))
      nudgdsize = MIN(MAX(nudgdsize,0.0_dp),0.5_dp)
      WRITE(mess,'(f8.6)') nudgdsize
      CALL message('    Damping radius (fraction)             ',TRIM(mess))

      CALL message('','Nudging options for input data format:')
      WRITE(mess,'(i0)') inudgformat
      CALL message('                               INUDGFORMAT  = ',TRIM(mess))

!EOX
      !*** define nudging start and stop date **********************************
!BOX
      lnudgstop = init_nudgingtime_a()

!EOX
      !*** synchronize model date with nudging files ***************************
!BOX
      IF (lfirst_cycle) THEN

         DO
            !--- compose the file name used for time adjustment

            IF (lnudgini) THEN
               !--- use nudging start date for file name evaluation
               CALL message('','Nudging initial time adjustment ...')
               CALL get_date_components(nudg_start,yr,mo,dy,hr,mi,se)

            ELSE IF (lresume) THEN
               !--- use next model date for name evaluation
               CALL message('','Nudging rerun time adjustment ...')
               CALL get_date_components(next_date,yr,mo,dy,hr,mi,se)

            ELSE

               EXIT                                    ! no evaluation needed
            END IF

            !--- which data set must be available ?
            DO i = 1,nmaxc
               IF (nudgd(i) >= 0._dp) THEN; linp_div = .TRUE.; EXIT; END IF
               IF (nudgv(i) >= 0._dp) THEN; linp_vor = .TRUE.; EXIT; END IF
               IF (nudgt(i) >= 0._dp) THEN; linp_tem = .TRUE.; EXIT; END IF
            END DO
                   
            SELECT CASE (inudgformat)   

            CASE(0)       ! budging data is in cray-format
         
               IF (linp_div) THEN
                  cfile = str_filter(ndg_file_div,yr,mo,dy,hr,mi,se,ndgblock)
               ELSE IF (linp_vor) THEN
                  cfile = str_filter(ndg_file_vor,yr,mo,dy,hr,mi,se,ndgblock)
               ELSE IF (linp_tem) THEN
                  cfile = str_filter(ndg_file_stp,yr,mo,dy,hr,mi,se,ndgblock)
               ELSE
                  CALL finish('NudgingInit','No external data defined.')
               END IF

               !>>dod bugfix
               CALL message('   Adjust date using file: ',TRIM(cfile))
               !<<dod

               IF (p_parallel_io) THEN
                 INQUIRE(file=cfile,exist=found)
                 IF (.NOT.found) &
                    CALL finish('NudgingInit','Nudging data file not found.')

                 CALL pbopen  (kfile, cfile, 'r', kret)
                 CALL cpbread (kfile, yhead, HEAD_LEN*WORD_LEN, kret)
                 CALL util_i8toi4 (yhead(1), ihead(1),HEAD_LEN)
                 ! check dimension
                 IF (ihead(5)*ihead(6) /= n2sp*nlev) CALL finish('NudgingInit',&
                      'Dimension of nudging data mismatch with model.')

               ENDIF

               IF (p_parallel) CALL p_bcast(ihead, p_io)
               CALL inp_convert_date(ihead(3),ihead(4)*10000, ndg_date0)

               IF (p_parallel_io) THEN
                 ! skip first record and read second header
                 ibytes = (ihead(5)*ihead(6) + HEAD_LEN) * WORD_LEN
                 CALL pbseek  (kfile, ibytes, 0, kret)
                 CALL cpbread (kfile, yhead, HEAD_LEN*WORD_LEN, kret)
                 CALL util_i8toi4 (yhead(1), ihead(1), HEAD_LEN)
                 CALL pbclose(kfile,kret)
               ENDIF

               IF (p_parallel) CALL p_bcast(ihead, p_io)
               CALL inp_convert_date(ihead(3),ihead(4)*10000, ndg_date1)

            CASE(2)        ! nudging data is in netcdf-format
            
               cfile = str_filter(ndg_file_nc,yr,mo,dy,hr,mi,se,ndgblock)
               CALL message('   Adjust date using file: ',TRIM(cfile))

               IF (p_parallel_io) THEN
                 INQUIRE(file=cfile,exist=found)
                 IF (.NOT.found) &
                    CALL finish('NudgingInit','Nudging data file not found.')
                 ndgfile%format = NETCDF
                 CALL IO_open (cfile, ndgfile, IO_READ)
                 CALL IO_INQ_DIMID(ndgfile%file_id, 'time', ndimid)
                 CALL IO_INQ_DIMLEN(ndgfile%file_id, ndimid, nts)
                 CALL IO_INQ_VARID(ndgfile%file_id, 'time', nvarid)
                 ALLOCATE (timevals(nts))
                 CALL IO_GET_VAR_DOUBLE (ndgfile%file_id, nvarid, timevals)
                 ihead_nc(1) = FLOOR(timevals(1))                          ! ihead_nc(1) is YYYYMMDD
                 ihead_nc(2) = INT((timevals(1)-ihead_nc(1))*24._dp)       ! ihead_nc(2) is HH
                 DEALLOCATE (timevals)
               ENDIF

               IF (p_parallel) CALL p_bcast(ihead_nc, p_io)
               CALL inp_convert_date(ihead_nc(1),ihead_nc(2)*10000, ndg_date0)

               IF (p_parallel_io) THEN
                 ! skip first record and read second header
                 ALLOCATE (timevals(nts))
                 CALL IO_GET_VAR_DOUBLE (ndgfile%file_id, nvarid, timevals)
                 ihead_nc(1) = FLOOR(timevals(2))                          ! ihead_nc(1) is YYYYMMDD
                 ihead_nc(2) = INT((timevals(2)-ihead_nc(1))*24._dp)       ! ihead_nc(2) is HH
                 DEALLOCATE (timevals)
                 CALL IO_close(ndgfile) 
               ENDIF

               IF (p_parallel) CALL p_bcast(ihead_nc, p_io)
               CALL inp_convert_date(ihead_nc(1),ihead_nc(2)*10000, ndg_date1)

            CASE DEFAULT

               CALL finish('NudgingInit','inudgformat set incorrectly')
    
            END SELECT 


            IF (lnudgini) THEN
               ! ndg_inp_date defines the nudging data file name
               ! which will be used later for file name generation
               ndg_inp_date = ndg_date0

               IF (ltintlin) THEN
                  ! linear interpolation needs only two synoptic times
                  ! therefore nudging will start at the first time
                  IF (nudg_start < ndg_date0) nudg_start = ndg_date0
               ELSE
                  ! nonlinear time interpolation needs four time levels
                  ! therefore nudging will start at the second time
                  IF (nudg_start < ndg_date1) nudg_start = ndg_date1
               END IF

               ! calculate nudging start date and adjust echam time manager
               CALL init_nudgingtime_b

            ELSE
               ! for the rerun case the first nudging date needed must be
               ! recalculated otherwise the automatic evaluation of file 
               ! names will fail
               CALL init_nudgingtime_c(ltintlin)
 
            END IF

            EXIT
      
      
        END DO
      END IF


      CALL message('','Nudging basic initialization done.')
      CALL message('','')
!EOX
   CASE(NDG_INI_STREAM) ! =========================================================

      !*** define the output stream ********************************************
!BOX
      CALL NudgingStreamInit

      CALL message('','Nudging output stream initialization done.')
      CALL message('','')


! ==============================================================================
!
!EOX
! *****  Setup nudging coefficients and work space **** PART2 nudging init

   CASE(NDG_INI_MEM) ! =========================================================
!BOX
      IF (lnudg_run) THEN
         IF (lnudgstop) THEN
            IF (nudg_stop < next_date) CALL StopNudging
         END IF
         RETURN               ! initialisation was performed
      END IF

      IF (next_date < nudg_start) RETURN  ! no nudging for next time step

      CALL message('','***********************************************')
      CALL message('','')
      CALL message('','Start nudging mode.')
      CALL message('','')
      CALL message('','***********************************************')

      lnudg_run = .TRUE.

!EOX
      !*** setup spectral window for nudging ***********************************
!BOX
      !>>dod bugfix
      CALL message('','Nudging window options:')
      !<<dod

      nudgsmin = MIN(MAX(nudgsmin,-1),nm)
      nudgsmax = MIN(MAX(nudgsmax, 0),nm)
      nudgsmax = MAX(nudgsmax,nudgsmin)

      !>>dod bugfix
      WRITE(mess,'(a26,i0,a5,i0)') '  Modified wavenumber(s) ',nudgsmin,' ... ',nudgsmax
      CALL message('',mess)
      !<<dod

      IF (nudgsmin<0) THEN
         CALL message('','     global average included')
         nudgmin  = 1
         nudgsmin = 0

      ELSE IF (nudgsmin == 0) THEN
         CALL message('','     global average NOT included')
         nudgmin = 2

      ELSE
         nudgmin = nmp(nudgsmin+1)+1
      ENDIF

      nudgmax = nmp(nudgsmax+1)+nnp(nudgsmax+1)
      nudgmax = MAX(nudgmax,nudgmin)
      !>>ddo bugfix
      WRITE(mess,'(a36, i0, a5,i0)') '  Modified spectral index SPECC = ',nudgmin, ' ... ',nudgmax
      !<<dod
      CALL message('',mess)


!EOX
      !*** setup level window for nudging **************************************
!BOX
      nudglmin = MIN(MAX(nudglmin,1),nlev)
      nudglmax = MIN(MAX(nudglmax,1),nlev)
      nudglmax = MAX(nudglmax,nudglmin)


!EOX
      !*** compose nudging flag field for level-spectral domain ****************
!BOX      !     For parallel version: Allocate a global flag field, 
      !     set the global field and scatter it to the local field

      nudgtrun = MAX(MIN(nudgtrun,2),0)
      IF (.NOT.ALLOCATED(flagn)) ALLOCATE(flagn(nlevp1,2,snsp))
      flagn_local => flagn; flagn = 0._dp

      ALLOCATE(flagn_global(nlevp1,2,nsp)); flagn_global(:,:,:) = 0._dp

      SELECT CASE(nudgtrun)
      CASE (NDG_WINDOW_ALL)
         !>>dod bugfix
         WRITE(mess,'(a14,i0,a29)') '  Truncation (',nudgtrun,') ... USE all meridional parts'
         !<<dod
         flagn_global(nudglmin:nudglmax,:,nudgmin:nudgmax) = 1._dp
         flagn_global(nlevp1           ,:,nudgmin:nudgmax) = 1._dp

      CASE (NDG_WINDOW_CUT)
         !>>dod bugfix
         WRITE(mess,'(a14,i0,a29)') '  Truncation (',nudgtrun,')... cut off meridional parts'
         !<<dod
         DO jm=nudgsmin,nudgsmax
            is1 = nmp(jm+1) + 1
            is2 = nmp(jm+1) + 1 + nudgsmax - jm
            flagn_global(nudglmin:nudglmax,:,is1:is2) = 1._dp
            flagn_global(nlevp1           ,:,is1:is2) = 1._dp
         END DO

      CASE (NDG_WINDOW_CUT0)
         !>>dod bugfix
         WRITE(mess,'(a14,i0,a29)') '  Truncation (',nudgtrun,')... cut off except wave 0'
         !<<dod
         DO jm=nudgsmin,nudgsmax
            IF (jm==0) THEN
               is1 = nmp(1) + 1
               is2 = nmp(1) + nnp(1)
            ELSE
               is1 = nmp(jm+1) + 1
               is2 = nmp(jm+1) + 1 + nudgsmax - jm
            END IF
            flagn_global(nudglmin:nudglmax,:,is1:is2) = 1._dp
            flagn_global(nlevp1           ,:,is1:is2) = 1._dp
         END DO
      END SELECT

      CALL message('',mess)
      !>>dod bugfix
      WRITE(mess,'(a7,i0)') '      [',NDG_WINDOW_ALL
      CALL message(TRIM(mess),' ... all meridional parts used]')
      WRITE(mess,'(a7,i0)') '      [',NDG_WINDOW_CUT
      CALL message(mess,' ... triangular cut-off]')
      WRITE(mess,'(a7,i0)') '      [',NDG_WINDOW_CUT0
      CALL message(mess,' ... triangular cut-off except wave 0]')
      !<<dod

      !--- check flag field
      sum1 = 0; sum2 = 0
      DO i=1,nlevp1
         DO ii=1,nsp
            sum1 = INT(sum1 + flagn_global(i,1,ii)+flagn_global(i,2,ii))
            sum2 = INT(sum2 + (1-flagn_global(i,1,ii)) + (1-flagn_global(i,2,ii)))
         END DO
      END DO
      IF ((sum1+sum2)/=(n2sp*nlevp1)) THEN
         !>>dod bugfix
         WRITE(mess,'(a7,i0,a7,i0,a5,i0)') ' FLAGN= ',sum1,' FLAGO= ',sum2,' NO= ',n2sp*nlevp1
         !<<dod
         CALL message('',mess)
         CALL finish('NudgingInit','ERROR nudging flag field mismatch')
      END IF
      CALL message('','')

      CALL scatter_sp(flagn_global,flagn_local,gdc)

      DEALLOCATE(flagn_global) ! not needed any more


!EOX
      !*** which data channels are needed ? ************************************
!BOX  
      !--- find maximum level again, but using now correct number of levels
      linp_tem = .FALSE.
      linp_div = .FALSE.
      linp_vor = .FALSE.
      DO i = 1,nlev
         IF (nudgt(i) >= 0._dp) THEN; linp_tem = .TRUE.; ilev_t_max = i; END IF
         IF (nudgd(i) >= 0._dp) THEN; linp_div = .TRUE.; ilev_d_max = i; END IF
         IF (nudgv(i) >= 0._dp) THEN; linp_vor = .TRUE.; ilev_v_max = i; END IF
      END DO
           
      !--- clear upper part in coefficient arrays
      IF (linp_tem .AND. ilev_t_max < nlev) nudgt(ilev_t_max+1:nlev) = 0.0_dp
      IF (linp_div .AND. ilev_d_max < nlev) nudgd(ilev_d_max+1:nlev) = 0.0_dp
      IF (linp_vor .AND. ilev_v_max < nlev) nudgv(ilev_v_max+1:nlev) = 0.0_dp
 
      !--- initialize lower bound
      ilev_t_min = ilev_t_max
      ilev_d_min = ilev_d_max
      ilev_v_min = ilev_v_max
 
      !-- find minimum level
      DO i = ilev_t_max,1,-1; IF (nudgt(i) < 0._dp) EXIT; ilev_t_min = i
      END DO
      DO i = ilev_d_max,1,-1; IF (nudgd(i) < 0._dp) EXIT; ilev_d_min = i
      END DO
      DO i = ilev_v_max,1,-1; IF (nudgv(i) < 0._dp) EXIT; ilev_v_min = i
      END DO


!EOX
      !*** rescale nudging coefficients and merge with level window ************
!BOX
      !--- allocate weights of observations
      IF (nlev > nmaxc) CALL finish('NudgingInit','number of levels too large')
      IF (.NOT. ALLOCATED(nudgdo)) &
           ALLOCATE(nudgdo(nlev), nudgvo(nlev), nudgto(nlev))
      nudgto(:) = 0.0_dp ; nudgdo(:) = 0.0_dp ; nudgvo(:) = 0.0_dp

      IF (nudglmin > 1) THEN
         nudgt(1:nudglmin-1) = 0.0_dp
         nudgd(1:nudglmin-1) = 0.0_dp
         nudgv(1:nudglmin-1) = 0.0_dp
      ENDIF
      IF (nudglmax < nlev) THEN
         nudgt(nudglmax+1:nlev) = 0.0_dp
         nudgd(nudglmax+1:nlev) = 0.0_dp
         nudgv(nudglmax+1:nlev) = 0.0_dp
      ENDIF

      !--- clear lower part in coefficient arrays
      IF (ilev_t_min > 1) nudgt(1:ilev_t_min-1) = 0.0_dp
      IF (ilev_d_min > 1) nudgd(1:ilev_d_min-1) = 0.0_dp
      IF (ilev_v_min > 1) nudgv(1:ilev_v_min-1) = 0.0_dp
 
      IF (linp_tem) ino_t_lev = ilev_t_max - ilev_t_min + 1
      IF (linp_div) ino_d_lev = ilev_d_max - ilev_d_min + 1
      IF (linp_vor) ino_v_lev = ilev_v_max - ilev_v_min + 1
      IF (nudgp >= 0._dp) THEN
         linp_lnp = .TRUE.
         ino_t_lev = ino_t_lev + 1
         IF (.NOT. linp_tem) ilev_t_min = 1
      ELSE
         nudgp = 0.0_dp
      END IF

      CALL message('','Define external data amount:')
      !>>dod bugfix
      IF (linp_tem) THEN
         WRITE(mess,'(i0,a1,i0,a1)') ilev_t_min,':',ilev_t_max,')'
         CALL message('   Use temperature LEVELS = (',mess)
      END IF
      IF (linp_lnp) THEN
         CALL message('','   Use surface fields (log surface pressure)')
      END IF
      IF (linp_div) THEN
         WRITE(mess,'(i0,a1,i0,a1)') ilev_d_min,':',ilev_d_max,')'
         CALL message('   Use divergence  LEVELS = (',mess)
      END IF
      IF (linp_vor) THEN
         WRITE(mess,'(i0,a1,i0,a1)') ilev_v_min,':',ilev_v_max,')'
         CALL message('   Use vorticity   LEVELS = (',mess)
      END IF
      !<<dod

      !--- correct nudging coefficients
      IF (.NOT. linp_tem) nudgt(:) = 0.0_dp
      IF (.NOT. linp_lnp) nudgp    = 0.0_dp
      IF (.NOT. linp_div) nudgd(:) = 0.0_dp
      IF (.NOT. linp_vor) nudgv(:) = 0.0_dp

      !--- merge level intervalls in nudging coefficient field
      IF (linp_vor .AND. linp_div .AND. linp_tem) THEN
         nudglmin = MAX(nudglmin, MIN(ilev_v_min, ilev_d_min, ilev_t_min))
         nudglmax = MIN(nudglmax, MAX(ilev_v_max, ilev_d_max, ilev_t_max))

      ELSE IF (linp_vor .AND. linp_div ) THEN
         nudglmin = MAX(nudglmin, MIN(ilev_v_min, ilev_d_min))
         nudglmax = MIN(nudglmax, MAX(ilev_v_max, ilev_d_max))

      ELSE IF (linp_vor .AND. linp_tem) THEN
         nudglmin = MAX(nudglmin, MIN(ilev_v_min, ilev_t_min))
         nudglmax = MIN(nudglmax, MAX(ilev_v_max, ilev_t_max))

      ELSE IF (linp_div .AND. linp_tem) THEN
         nudglmin = MAX(nudglmin, MIN(ilev_d_min, ilev_t_min))
         nudglmax = MIN(nudglmax, MAX(ilev_d_max, ilev_t_max))

      ELSE IF (linp_vor) THEN
         nudglmin = MAX(nudglmin, ilev_v_min)
         nudglmax = MIN(nudglmax, ilev_v_max)

      ELSE IF (linp_div) THEN
         nudglmin = MAX(nudglmin, ilev_d_min)
         nudglmax = MIN(nudglmax, ilev_d_max)

      ELSE IF (linp_tem) THEN
         nudglmin = MAX(nudglmin, ilev_t_min)
         nudglmax = MIN(nudglmax, ilev_t_max)

      END IF

      !>>dod bugfix
      WRITE(mess,'(i0,a5,i0)') nudglmin,' ... ',nudglmax
      CALL message('  Assimilate external data into LEVEL(s) = ',mess)
      !<<dod


!EOX
      !*** prepare nudging weights level dependend *****************************

!BOX
      ! set weights for explicit nudging: new = N(read)*2*dt*nfact

      !--- rescale coefficients
      twodt = 2._dp*delta_time
      nudgpo    = nudgp         * nfact * twodt
      nudgto(:) = nudgt(1:nlev) * nfact * twodt
      nudgdo(:) = nudgd(1:nlev) * nfact * twodt
      nudgvo(:) = nudgv(1:nlev) * nfact * twodt

      !--- set weights for model data
      IF (lnudgpat) THEN
         nudgp          = 1._dp
         nudgt(1:nlev)  = 1._dp
         nudgd(1:nlev)  = 1._dp
         nudgv(1:nlev)  = 1._dp
      ELSE
         nudgp          = 1._dp-nudgpo
         nudgt(1:nlev)  = 1._dp-nudgto(:)
         nudgd(1:nlev)  = 1._dp-nudgdo(:)
         nudgv(1:nlev)  = 1._dp-nudgvo(:)
      END IF
             
      !--- convert coefficients for implicit nudging
      IF (lnudgimp) THEN
         !--- using implicit nudging reset weights
         nudgp         = 1._dp/(1._dp+nudgpo)
         nudgt(1:nlev) = 1._dp/(1._dp+nudgto(:))
         nudgd(1:nlev) = 1._dp/(1._dp+nudgdo(:))
         nudgv(1:nlev) = 1._dp/(1._dp+nudgvo(:))
             
         !--- reset weights for observations
         nudgpo    = nudgpo    * nudgp
         nudgto(:) = nudgto(:) * nudgt(1:nlev) 
         nudgdo(:) = nudgdo(:) * nudgd(1:nlev) 
         nudgvo(:) = nudgvo(:) * nudgv(1:nlev) 
      ENDIF

      !--- print out the time constant part of the nudging weights
      CALL message('',' ')
      !>>dod bugfix
      CALL message('','Nudging coefficients (1/sec) and weights for model and reference (%)')
      !<<dod
      WRITE(mess,'(a5,3a25,7x)') 'lev','TEM','DIV','VOR'
      CALL message('',mess)

      IF (lnudgimp) THEN
         CALL message('','    IMPLICIT nudging used')
         DO i=1,nlev
            WRITE(mess,'(i5,3(e12.3,2f10.3))') i,&
                 (1.0_dp-nudgt(i))/(nudgt(i)*twodt),nudgt(i)*100.0_dp,nudgto(i)*100.0_dp,&
                 (1.0_dp-nudgd(i))/(nudgd(i)*twodt),nudgd(i)*100.0_dp,nudgdo(i)*100.0_dp,&
                 (1.0_dp-nudgv(i))/(nudgv(i)*twodt),nudgv(i)*100.0_dp,nudgvo(i)*100.0_dp
            CALL message('',mess)
         ENDDO
         WRITE(mess,'(a,e12.3,2f10.3)') ' surface fields = ',&
              (1.0_dp-nudgp)/(nudgp*twodt),nudgp*100.0_dp,nudgpo*100._dp
         CALL message('',mess)

      ELSE
         CALL message('','    EXPLICIT nudging used')
         DO i=1,nlev
            WRITE(mess,'(i5,3(e12.3,2f10.3))') i,&
                 nudgto(i)/twodt,nudgt(i)*100.0_dp,nudgto(i)*100.0_dp,&
                 nudgdo(i)/twodt,nudgd(i)*100.0_dp,nudgdo(i)*100.0_dp,&
                 nudgvo(i)/twodt,nudgv(i)*100.0_dp,nudgvo(i)*100.0_dp
            CALL message('',mess)
         ENDDO
         WRITE(mess,'(a,e12.3,2f10.3)') ' surface fields = ',&
              nudgpo/twodt,nudgp*100.0_dp,nudgpo*100.0_dp
         CALL message('',mess)

      ENDIF

      CALL message('','')


!EOX
      !*** allocate memory *****************************************************
!BOX
      IF (.NOT. lndgstep3) CALL NudgingAllocRefMem

      CALL NudgingAllocMem

      CALL set_stream(nudg, lpout=.TRUE.)  ! switch output stream on

      CALL NdgInitCounter

      lnudgini = .FALSE.



! ==============================================================================
!
!EOX
!****** clean up work space memory

   CASE(NDG_CLEAN_MEM)
!BOX
      IF (lnudg_run) THEN
         CALL message('','Clean up nudging memory:')
         CALL NudgingDeallocMem
         lnudg_run = .FALSE.
      ELSE
         CALL message('','No nudging memory allocated.')
      END IF


! ==============================================================================
!
!EOX
!****** clean up i/o memory and close nudging

   CASE (NDG_CLOSE)
!BOX
      IF (lnmi) lnmi = .FALSE.

      CALL  NudgingDeallocRefMem

      lndgstep0 = .FALSE.
      lndgstep1 = .FALSE.
      lndgstep2 = .FALSE.
      lndgstep3 = .FALSE.

      CALL CloseBlock
!++jsr (bugfix: added nsstinc condition)
      IF (nsstinc > 0) THEN
         CALL NudgingSSTClose(.TRUE.)
      END IF

      CALL message('','Close nudging data files.')

      CALL set_stream(nudg, lpout=.FALSE., lrerun=.FALSE.)  ! switch output stream off

      lnudge = .FALSE.

      CALL message('','***********************************************')
      CALL message('','')
      CALL message('','Stop nudging mode.')
      CALL message('','')
      CALL message('','***********************************************')


! ==============================================================================

   CASE default

      CALL finish('NudgingInit','Wrong initialization mode.')
      
   END SELECT



 CONTAINS

   SUBROUTINE bc_nudging
     
     USE mo_mpi, ONLY: p_bcast, p_io

     CALL p_bcast (lnudgini,      p_io)
     CALL p_bcast (lnudgcli,      p_io)
     CALL p_bcast (dt_nudg_start, p_io)
     CALL p_bcast (dt_nudg_stop,  p_io)

     CALL p_bcast (ndg_file_vor,  p_io)
     CALL p_bcast (ndg_file_div,  p_io)
     CALL p_bcast (ndg_file_stp,  p_io)
     CALL p_bcast (ndg_file_nc,   p_io)

     CALL p_bcast (nsstinc,       p_io)
     CALL p_bcast (nsstoff,       p_io)
     CALL p_bcast (ndg_freez,     p_io)
     CALL p_bcast (ndg_file_sst,  p_io)

     CALL p_bcast (lnudgimp,      p_io)
     CALL p_bcast (lnudgpat,      p_io)
     CALL p_bcast (lnudgfrd,      p_io)

     CALL p_bcast (nudglmin,      p_io)
     CALL p_bcast (nudglmax,      p_io)

     CALL p_bcast (nudgtrun,      p_io)
     CALL p_bcast (nudgsmin,      p_io)
     CALL p_bcast (nudgsmax,      p_io)

     CALL p_bcast (nudgd,         p_io)
     CALL p_bcast (nudgv,         p_io)
     CALL p_bcast (nudgt,         p_io)
     CALL p_bcast (nudgp,         p_io)

     CALL p_bcast (ltintlin,      p_io)
     CALL p_bcast (ldamplin,      p_io)
     CALL p_bcast (nudgdamp,      p_io)
     CALL p_bcast (nudgdsize,     p_io)

     CALL p_bcast (lnudgdbx,      p_io)
     CALL p_bcast (lnudgwobs,     p_io)
     CALL p_bcast (lsite,         p_io)
     
     CALL p_bcast (inudgformat,   p_io)

   END SUBROUTINE bc_nudging

  END SUBROUTINE NudgingInit
!EOX
!EOC
! ==============================================================================
!BOP
  ! !IROUTINE:  StopNudging
  ! !INTERFACE:

  SUBROUTINE StopNudging

    ! !DESCRIPTION: 
    ! close nudging, clean up memory

!EOP
    CALL NudgingInit(NDG_CLEAN_MEM)
    CALL NudgingInit(NDG_CLOSE)

  END SUBROUTINE StopNudging
END MODULE mo_nudging_init
