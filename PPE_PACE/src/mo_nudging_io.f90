!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_nudging_io
!BOP
  ! !MODULE: mo_nudging_io (layer 2)

  ! !DESCRIPTION: 
  ! perform nudging I/O operations

  ! !REVISION HISTORY: 
  ! Ingo Kirchner, MPI Hamburg, Aug-2002, revision
  ! Aiko Voigt, MPI Hamburg, March 2011, implement netcdf-based reading of nudging data
  !   - the switch inudgformat (in ReadOneBlock) allows to specify whether the
  !     nudging data should be read in CRAY bindary format (inudgformat=0) or 
  !     netcdf-format (inudgformat=2) 

  ! !USES:
  USE mo_kind,          ONLY: dp
  USE mo_netCDF,        ONLY: FILE_INFO        !! AV: needed for variable ndgfile that
                                               !!     contains the netcdf-file with nudging data

!BOX
  IMPLICIT NONE
  PRIVATE
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: GetNudgeData
  PUBLIC :: CloseBlock

  ! !PUBLIC DATA MEMBERS:

  ! strings for file name generation
  CHARACTER(len=256), PUBLIC :: ndg_file_vor
  CHARACTER(len=256), PUBLIC :: ndg_file_div
  CHARACTER(len=256), PUBLIC :: ndg_file_stp
  CHARACTER(len=256), PUBLIC :: ndg_file_nc    !! AV: file name of the netcdf-file with nudging data
!EOP

  ! internal file ID of nudging data blocks and sst fields
  INTEGER, PUBLIC      :: kfiled, kfilev, kfilet
  INTEGER, PUBLIC      :: ndgblock =  0            ! block pointer nudging data
  
  CHARACTER(len=256) :: mess

  ! **** variables for I/O

  LOGICAL :: files_open = .FALSE.

  INTEGER, PUBLIC :: &
       ndgstep0  = 0,      ndgstep1  = 0,      ndgstep2  = 0,      ndgstep3  = 0
  LOGICAL, PUBLIC :: &
       lndgstep0 =.FALSE., lndgstep1 =.FALSE., lndgstep2 =.FALSE., lndgstep3 =.FALSE.
  
  TYPE (FILE_INFO), PUBLIC   :: ndgfile        !! AV: netcdf-file with nudging data 
       

CONTAINS

  !=============================================================================
!BOP
  ! !IROUTINE:  GetNudgeData
  ! !INTERFACE:

  SUBROUTINE GetNudgeData

    ! !DESCRIPTION: 
    ! read new data for nudging, dependend on the time step
    ! search next possible nudging data which cover the model time
    
    ! !USES:
    USE mo_time_control, ONLY: &
         ndg_date0, ndg_date1, ndg_date2, ndg_date3, &
         nudging_date_fit, write_date

    USE mo_nudging_buffer, ONLY: &
         buf_ref0, buf_ref1, buf_ref2, buf_ref3, &
         a_pat_div, a_pat_vor, a_pat_tem, a_pat_lnp, &
         a_nrm_div, a_nrm_vor, a_nrm_tem, a_nrm_lnp

    USE mo_nudging_constants, ONLY: &
         lnudgpat, ltintlin, lnudgdbx, &
         linp_vor, linp_div, linp_tem, linp_lnp
    USE mo_exception,        ONLY: message

!EOP
    read_data_loop : DO

      ! rotate buffer **********************************************************
      buf_ref0 = buf_ref1
      buf_ref1 = buf_ref2
      buf_ref2 = buf_ref3
      buf_ref3 = 0.0_dp

       IF (lndgstep1) THEN
          ! now step 1 available ***********************************************
          lndgstep0 = .TRUE.; ndg_date0 = ndg_date1; ndgstep0  = ndgstep1
       END IF

       IF (lndgstep2) THEN
          ! now step 2 available ***********************************************
          lndgstep1 = .TRUE.; ndg_date1 = ndg_date2; ndgstep1  = ndgstep2
          IF (lnudgpat) THEN
             IF(linp_vor) THEN
                a_pat_vor(:,2) = a_pat_vor(:,3); a_nrm_vor(:,2) = a_nrm_vor(:,3)
             END IF
             IF(linp_div) THEN
                a_pat_div(:,2) = a_pat_div(:,3); a_nrm_div(:,2) = a_nrm_div(:,3)
             END IF
             IF(linp_tem) THEN
                a_pat_tem(:,2) = a_pat_tem(:,3); a_nrm_tem(:,2) = a_nrm_tem(:,3)
             END IF
             IF(linp_lnp) THEN
                a_pat_lnp  (2) = a_pat_lnp  (3); a_nrm_lnp  (2) = a_nrm_lnp  (3)
             END IF
          END IF
       END IF

       IF (lndgstep3)  THEN
          ! now step 3 available ***********************************************
          lndgstep2 = .TRUE.; ndg_date2 = ndg_date3; ndgstep2  = ndgstep3
          IF (lnudgpat) THEN
             IF(linp_vor) THEN
                a_pat_vor(:,3) = a_pat_vor(:,4); a_nrm_vor(:,3) = a_nrm_vor(:,4)
             END IF
             IF(linp_div) THEN
                a_pat_div(:,3) = a_pat_div(:,4); a_nrm_div(:,3) = a_nrm_div(:,4)
             END IF
             IF(linp_tem) THEN
                a_pat_tem(:,3) = a_pat_tem(:,4); a_nrm_tem(:,3) = a_nrm_tem(:,4)
             END IF
             IF(linp_lnp) THEN
                a_pat_lnp  (3) = a_pat_lnp  (4); a_nrm_lnp  (3) = a_nrm_lnp  (4)
             END IF
          END IF
       END IF


       CALL ReadOneBlock ! read next data block ********************************

       IF ((lndgstep1 .AND. ltintlin) .OR. lndgstep0) THEN

          ! minimum of four steps 
          ! is necessary for non-linear time interpolation
          ! otherwise only two steps are needed

          IF (nudging_date_fit()) EXIT
       END IF

    END DO read_data_loop


    IF (lnudgdbx) THEN
       WRITE(mess,'(a32,i0,a5,i0)') 'Interpolation between ... steps ',ndgstep1,' ... ',ndgstep2
       CALL message('GetNudgeData',mess)
       CALL write_date(ndg_date1,' First point date : ')
       CALL write_date(ndg_date2,' Second point date: ')
    END IF

  END SUBROUTINE GetNudgeData

  !=============================================================================
!BOP
  ! !IROUTINE:  ReadOneBlock
  ! !INTERFACE:

  SUBROUTINE ReadOneBlock

    ! !DESCRIPTION: 
    ! read next nudging data time step

    ! !USES:
    USE mo_exception,     ONLY: finish, message
    USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_io, p_bcast
    USE mo_control,          ONLY: n2sp, lnmi, nlevp1
    USE mo_time_control,     ONLY: &
         get_step_from_header, get_date_components, &
         ndg_date3, next_date, write_date
    USE mo_nmi, ONLY: NMI_Make, NMI_MAKE_FOD

    USE mo_nudging_constants, ONLY: &
         linp_div, linp_vor, linp_tem, linp_lnp, &
         ino_d_lev,  ino_v_lev,  ino_t_lev, &
         ilev_d_min, ilev_v_min, ilev_t_min, &
         ilev_d_max, ilev_v_max,  &
         lnudgcli, lnudgpat, lnudgdbx, lnudgfrd, &
         inudgformat
    USE mo_nudging_buffer,   ONLY: &
         sdobs0, svobs0, stobs0, &
         sdobs3, svobs3, stobs3, &
         a_pat_div, a_pat_vor, a_pat_tem, a_pat_lnp
    USE mo_nudging_pattern,  ONLY: Nudg_Init_Alpha
#ifdef LITTLE_ENDIAN
    USE mo_nudging_utils,    ONLY: cpbread, WORD_LEN, HEAD_LEN, swap64
#else
    USE mo_nudging_utils,    ONLY: cpbread, WORD_LEN, HEAD_LEN
#endif
!EOP

    USE mo_netCDF,        ONLY: IO_inq_dimid, IO_INQ_DIMLEN, &                        !! AV
                                IO_INQ_VARID,IO_GET_VAR_DOUBLE,IO_GET_VARA_DOUBLE     !! AV
    USE mo_decomposition, ONLY: gdc => global_decomposition                           !! AV
    USE mo_control,       ONLY: nlev, nsp                                             !! AV
    USE mo_transpose,     ONLY: scatter_sp                                            !! AV

    INTEGER           :: kret, iilen
    INTEGER           :: ilen_d, ilen_v, ilen_t
    INTEGER           :: iheadd(HEAD_LEN), iheadv(HEAD_LEN), iheadt(HEAD_LEN), inph(HEAD_LEN)
    INTEGER           :: yr1,mo1,dy1,hh1,mm1,ss1
    INTEGER           :: yr2,mo2,dy2,hh2,mm2,ss2
    CHARACTER (len=WORD_LEN) :: yhead(HEAD_LEN)
    REAL(kind=dp), ALLOCATABLE :: phbuf(:), zhbuf(:)
    LOGICAL  :: lrd_check, endofblock
    
    INTEGER, SAVE          :: ndtsc                              !! AV: time step counter for netcdf-data  
    INTEGER                :: ndimid, nts, nvarid, kk, jj, kkmax !! AV
    REAL(dp), ALLOCATABLE  :: timevals(:)                        !! AV
    INTEGER                :: ihead_nc(2)                        !! AV
    INTEGER                :: start(4), COUNT(4)                 !! AV
    REAL(kind=dp), ALLOCATABLE         :: zin(:,:,:)             !! AV: buffer for reading netcdf-data
    REAL(kind=dp), POINTER             :: zin2(:,:,:)            !! AV: temporary buffer for netcdf-data
 
    INTRINSIC MAX, RESHAPE

    SELECT CASE (inudgformat)
    
    CASE(0)    ! read cray data

       ilen_d = n2sp*ino_d_lev
       ilen_v = n2sp*ino_v_lev
       ilen_t = n2sp*ino_t_lev
       iilen  = MAX(ilen_d,ilen_v,ilen_t)

       ALLOCATE(phbuf(iilen))
       ALLOCATE(zhbuf(iilen))
       zhbuf(:) = 0.0_dp

       IF (.NOT. files_open) CALL OpenOneBlock  ! open first data block ***********
       
       DO

          endofblock = .FALSE.
          IF (p_parallel_io) THEN
            ! read next nudging data set header
            IF (linp_div) THEN
               CALL cpbread(kfiled,yhead, HEAD_LEN*WORD_LEN,kret)
               IF (kret==HEAD_LEN*WORD_LEN) THEN
                  CALL util_i8toi4 (yhead(1), iheadd(1),HEAD_LEN)
               ELSE 
                  CALL message('ReadOneBlock',&
                       'nudging DATA inconsistent DIV or EOF'); endofblock = .TRUE.
               END IF
               inph = iheadd
            END IF

            IF (linp_vor .AND. .NOT. endofblock) THEN
               CALL cpbread(kfilev,yhead,HEAD_LEN*WORD_LEN,kret)
               IF (kret==HEAD_LEN*WORD_LEN) THEN
                  CALL util_i8toi4 (yhead(1), iheadv(1),HEAD_LEN)
               ELSE 
                  CALL message('ReadOneBlock',&
                       'nudging DATA inconsistent VOR or EOF'); endofblock = .TRUE.
               END IF
               inph = iheadv
            END IF

            IF ((linp_tem .OR. linp_lnp) .AND. .NOT. endofblock) THEN
               CALL cpbread(kfilet,yhead,HEAD_LEN*WORD_LEN,kret)
               IF (kret==HEAD_LEN*WORD_LEN) THEN
                  CALL util_i8toi4 (yhead(1), iheadt(1),HEAD_LEN)
               ELSE
                  CALL message('ReadOneBlock',&
                       'nudging DATA inconsistent TEM or EOF'); endofblock = .TRUE.
               END IF
               inph = iheadt
            END IF
          END IF

          IF (p_parallel) THEN
            CALL p_bcast(endofblock, p_io)
            CALL p_bcast(iheadd, p_io)
            CALL p_bcast(iheadv, p_io)
            CALL p_bcast(iheadt, p_io)
            CALL p_bcast(inph,   p_io)
          END IF

          IF (endofblock) THEN          ! end of nudging data block reached
             CALL OpenOneBlock
             CYCLE                      ! read all again from new files
          ENDIF


          ! check header consistency **********************************************
          lrd_check = .TRUE.
          IF ( (linp_div .AND. linp_vor) .AND. &
               (.NOT. (iheadd(3)==iheadv(3) .AND. iheadd(4)==iheadv(4))) )&
               lrd_check = .FALSE.

          IF ( (linp_div .AND. (linp_tem.OR.linp_lnp) ) .AND. &
               (.NOT. (iheadd(3)==iheadt(3) .AND. iheadd(4)==iheadt(4))) )&
               lrd_check = .FALSE.

          IF ( (linp_vor .AND. (linp_tem.OR.linp_lnp) ) .AND. &
               (.NOT. (iheadv(3)==iheadt(3) .AND. iheadv(4)==iheadt(4))) )&
               lrd_check = .FALSE.



          ! usage of data blocks **************************************************
          !
          ! lnudgpat = false --> block0,1,2,3 contains 4 time levels
          ! lnudgpat = true  --> block1,2 contains 2 time levels, 
          !                      block 0 used for climatology
           
          IF (lrd_check) THEN

             IF (lnudgcli) THEN
                ! set ndg_date3 preliminary and get components
                ndgstep3 = get_step_from_header(inph(3), inph(4)*10000)
                CALL get_date_components(ndg_date3,yr1,mo1,dy1,hh1,mm1,ss1)
                CALL get_date_components(next_date,yr2,mo2,dy2,hh2,mm2,ss2)
                IF (ndgblock == 1) THEN
                   yr1 = yr2 -1
                ELSE IF (ndgblock == 2) THEN
                   yr1 = yr2
                ELSE
                   yr1 = yr2 + 1
                END IF
                inph(3) = yr1*10000+mo1*100+dy1
             END IF
             ndgstep3 = get_step_from_header(inph(3), inph(4)*10000)

             IF (linp_div) THEN
                sdobs3(:,:,:) = 0.0_dp
             
                IF (iheadd(5)*iheadd(6) /= ilen_d) CALL finish('ReadOneBlock',&
                    'nudging data fault DIV, dimension mismatch')
                IF (p_parallel_io) THEN
                   CALL pbread(kfiled,phbuf,ilen_d*WORD_LEN,kret)
                   IF (kret/=ilen_d*WORD_LEN) &
                      CALL finish('ReadOneBlock','nudging data fault DIV')
                END IF
                CALL convert_input(sdobs3,ilev_d_min,ilev_d_max,ino_d_lev)

                IF (lnudgpat) THEN
                ! NB: If lnudgpat is set, the program runs NOT in parallel!
                   sdobs0(:,:,:) = 0.0_dp
                   ! read additional fields without header
                   CALL pbread(kfiled,phbuf,ilen_d*WORD_LEN,kret)
                   IF (kret/=ilen_d*WORD_LEN) &
                      CALL finish('ReadOneBlock','climate data fault DIV')

                   ! the climatology
                   CALL convert_input(sdobs0,ilev_d_min,ilev_d_max,ino_d_lev)
                
                   ! read constant factor
                   CALL pbread(kfiled,phbuf,ino_d_lev*WORD_LEN,kret)
                   IF (kret/=ino_d_lev*WORD_LEN) &
                       CALL finish('ReadOneBlock','alpha data fault DIV')
#ifdef LITTLE_ENDIAN
                   CALL swap64(phbuf,ino_d_lev)
#endif
                   CALL util_cray2ieee(phbuf,zhbuf,ino_d_lev)
                   phbuf(1:ino_d_lev) = zhbuf(1:ino_d_lev)
                   a_pat_div(ilev_d_min:ilev_d_max,4) = phbuf(1:ino_d_lev)
                   IF (lnudgdbx) THEN
                      WRITE(mess,'(a,40f12.5)') 'ALPHA(div)= ',a_pat_div(:,4)
                      CALL message('ReadOneBlock',mess)
                   END IF
                END IF
 
             END IF
 
             IF (linp_vor) THEN
                svobs3(:,:,:) = 0.0_dp

                IF (iheadv(5)*iheadv(6) /= ilen_v) CALL finish('ReadOneBlock',&
                   'nudging data fault VOR, dimension mismatch')
                IF (p_parallel_io) THEN
                   CALL pbread(kfilev,phbuf,ilen_v*WORD_LEN,kret)
                   IF (kret/=ilen_v*WORD_LEN) &
                      CALL finish('ReadOneBlock','nudging data fault VOR')
                END IF
                CALL convert_input(svobs3,ilev_v_min,ilev_v_max,ino_v_lev)

                IF (lnudgpat) THEN
                ! NB: If lnudgpat is set, the program runs NOT in parallel!
                   svobs0(:,:,:) = 0.0_dp
                   ! read additional fields without header
                   CALL pbread(kfilev,phbuf,ilen_v*WORD_LEN,kret)
                   IF (kret/=ilen_v*WORD_LEN) &
                      CALL finish('ReadOneBlock','climate data fault VOR')

                   ! the climatology
                   CALL convert_input(svobs0,ilev_v_min,ilev_v_max,ino_v_lev)

                   ! read constant factor
                   CALL pbread(kfilev,phbuf,ino_v_lev*WORD_LEN,kret)
                   IF (kret/=ino_v_lev*WORD_LEN) &
                      CALL finish('ReadOneBlock','alpha data fault VOR')
#ifdef LITTLE_ENDIAN
                   CALL swap64(phbuf,ino_v_lev)
#endif
                   CALL util_cray2ieee(phbuf,zhbuf,ino_v_lev)
                   phbuf(1:ino_v_lev) = zhbuf(1:ino_v_lev)
                   a_pat_vor(ilev_v_min:ilev_v_max,4) = phbuf(1:ino_v_lev)

                   IF (lnudgdbx) THEN
                      WRITE(mess,'(a,40f12.5)') 'ALPHA(vor)= ',a_pat_vor(:,4)
                      CALL message('ReadOneBlock',mess)
                   END IF
                END IF
             
             END IF
 
             IF (linp_tem .OR. linp_lnp) THEN
                stobs3(:,:,:) = 0.0_dp
 
                IF (iheadt(5)*iheadt(6) /= ilen_t) CALL finish('ReadOneBlock',&
                   'nudging data fault TEM/LNP, dimension mismatch')
                IF (p_parallel_io) THEN
                   CALL pbread(kfilet,phbuf,ilen_t*WORD_LEN,kret)
                   IF (kret/=ilen_t*WORD_LEN) &
                      CALL finish('ReadOneBlock','nudging data fault TEM/LNP')
                END IF
                CALL convert_input&
                     (stobs3,ilev_t_min,(ilev_t_min+ino_t_lev-1),ino_t_lev)
                IF (linp_lnp .AND. ino_t_lev < nlevp1) THEN
                   stobs3(nlevp1,:,:)                 = stobs3(ilev_t_min+ino_t_lev-1,:,:)
                   stobs3(ilev_t_min+ino_t_lev-1,:,:) = 0.0_dp
                END IF
 
                IF (lnudgpat) THEN
                ! NB: If lnudgpat is set, the program runs NOT in parallel!
                   stobs0(:,:,:) = 0.0_dp
                
                   ! read additional fields without header
                   CALL pbread(kfilet,phbuf,ilen_t*WORD_LEN,kret)
                   IF (kret/=ilen_t*WORD_LEN) &
                      CALL finish('ReadOneBlock','climate data fault TEM/LNP')
                
                   ! the climatology
                   CALL convert_input&
                        (stobs0,ilev_t_min,(ilev_t_min+ino_t_lev-1),ino_t_lev)
                
                   IF (linp_lnp .AND. ino_t_lev < nlevp1) THEN
                      stobs0(nlevp1,:,:)       = stobs0(ilev_t_min+ino_t_lev-1,:,:)
                      stobs0(ilev_t_min+ino_t_lev-1,:,:) = 0.0_dp
                   END IF
 
                   ! read constant factor
                   CALL pbread(kfilet,phbuf,ino_t_lev*WORD_LEN,kret)
                   IF (kret/=ino_t_lev*WORD_LEN) CALL finish('ReadOneBlock','alpha data fault TEM/LNP')
#ifdef LITTLE_ENDIAN
                   CALL swap64(phbuf,ino_t_lev)
#endif
                   CALL util_cray2ieee(phbuf,zhbuf,ino_t_lev)
                   phbuf(1:ino_t_lev) = zhbuf(1:ino_t_lev)
                   a_pat_tem(ilev_t_min:ilev_t_min+ino_t_lev-1,4) = phbuf(1:ino_t_lev)
                   IF (linp_lnp) THEN
                      a_pat_lnp(4) = a_pat_tem(ilev_t_min+ino_t_lev-1,4)
                      a_pat_tem(ilev_t_min+ino_t_lev-1,4) = 0.0_dp
                      IF (lnudgdbx) THEN
                         WRITE(mess,'(a,40f12.5)') 'ALPHA(lnp)= ',a_pat_lnp(4)
                         CALL message('ReadOneBlock',mess)
                      END IF
                   END IF
                   IF (linp_tem) THEN
                      IF (lnudgdbx) THEN
                         WRITE(mess,'(a,40f12.5)') 'ALPHA(tem)= ',a_pat_tem(:,4)
                         CALL message('ReadOneBlock',mess)
                      END IF
                   END IF
 
                END IF
 

                ! additional calulations for pattern nudging **********************
 
                IF (lnudgpat) THEN
                   CALL Nudg_Init_Alpha

                   IF (lnudgdbx) THEN
                      CALL message('ReadOneBlock','ALPHA corrected')
                      IF (linp_div) THEN
                         WRITE(mess,'(a,40f12.5)') 'ALPHA(div)= ',a_pat_div(:,4)
                         CALL message('ReadOneBlock',mess)
                      END IF
                      IF (linp_vor) THEN
                         WRITE(mess,'(a,40f12.5)') 'ALPHA(vor)= ',a_pat_vor(:,4)
                         CALL message('ReadOneBlock',mess)
                      END IF
                      IF (linp_tem) THEN
                         WRITE(mess,'(a,40f12.5)') 'ALPHA(tem)= ',a_pat_tem(:,4)
                         CALL message('ReadOneBlock',mess)
                      END IF
                      IF (linp_lnp) THEN
                         WRITE(mess,'(a,40f12.5)') 'ALPHA(lnp)= ',a_pat_lnp(4)
                         CALL message('ReadOneBlock',mess)
                      END IF

                   END IF
 
                END IF
 
             END IF
  
             lndgstep3 = .TRUE.

             WRITE(mess,'(a48,i0)') 'READ NEW nudging DATA, corresponding step is : ',ndgstep3
             CALL message('ReadOneBlock',mess)
             CALL write_date(ndg_date3,' Corresponding date: ')

             ! perform NMI filter now (optional) **********************************
             IF (lnmi .AND. lnudgfrd) CALL NMI_Make(NMI_MAKE_FOD)

             EXIT ! read block OBS3 HEAD3 success

          ELSE

             CALL finish('ReadOneBlock','Header synchronisation fault.')

          ENDIF

       ENDDO

       DEALLOCATE (zhbuf)
       DEALLOCATE (phbuf)

    CASE(2)                  !! read netcdf data

       IF (.NOT. files_open) THEN       
          CALL OpenOneBlock  !! open first data block ***********
          ndtsc=0            !! set time step counter ndtsc to zero when new nudging data file is opened
       END IF  
    
       DO
       
          endofblock = .FALSE.
          IF (p_parallel_io) THEN
             ! read next nudging data set header
             ndtsc=ndtsc + 1                                               ! increase counter by one
             CALL IO_INQ_DIMID(ndgfile%file_id, 'time', ndimid)
             CALL IO_INQ_DIMLEN(ndgfile%file_id, ndimid, nts)
             CALL IO_INQ_VARID(ndgfile%file_id, 'time', nvarid)
             ALLOCATE (timevals(nts))
             CALL IO_GET_VAR_DOUBLE (ndgfile%file_id, nvarid, timevals)
             ihead_nc(1) = FLOOR(timevals(ndtsc))                          ! ihead_nc(1) is YYYYMMDD
             ihead_nc(2) = INT((timevals(ndtsc)-ihead_nc(1))*24._dp)       ! ihead_nc(2) is HH 
             DEALLOCATE (timevals)
             IF (ndtsc.gt.nts) endofblock = .TRUE.                         ! end of netcdf-file is reached
          END IF
  
          IF (p_parallel) THEN  
             CALL p_bcast(endofblock, p_io)
             CALL p_bcast(ihead_nc, p_io)
          END IF

          IF (endofblock) THEN          ! end of nudging data block reached
             CALL OpenOneBlock
             ndtsc=0                    ! set time step counter ndtsc to zero when new nudging data file is opened
             CYCLE                      ! read all again from new files
          END IF
     
          IF (lnudgcli) THEN
             ! set ndg_date3 preliminary and get components
             ndgstep3 = get_step_from_header(ihead_nc(1), ihead_nc(2)*10000)
             CALL get_date_components(ndg_date3,yr1,mo1,dy1,hh1,mm1,ss1)
             CALL get_date_components(next_date,yr2,mo2,dy2,hh2,mm2,ss2)
             IF (ndgblock == 1) THEN
                yr1 = yr2 -1
             ELSE IF (ndgblock == 2) THEN
                yr1 = yr2
             ELSE
                yr1 = yr2 + 1
             END IF
             ihead_nc(1) = yr1*10000+mo1*100+dy1
          END IF                       ! end of lnudgcli
          ndgstep3 = get_step_from_header(ihead_nc(1), ihead_nc(2)*10000)

          sdobs3(:,:,:) = 0.0_dp
          svobs3(:,:,:) = 0.0_dp
          stobs3(:,:,:) = 0.0_dp
  
          !! read divergence       
          IF (linp_div) THEN
             IF (p_parallel_io) THEN
                CALL IO_INQ_VARID (ndgfile%file_id, 'sd', nvarid)
                COUNT(:) = (/ 2, nsp,nlev,1 /)
                start(:) = (/ 1, 1, 1,ndtsc /)
                ALLOCATE(zin(2,nsp,nlev))
                zin(:,:,:) = 0.0_dp
                CALL IO_GET_VARA_DOUBLE (ndgfile%file_id, nvarid, start, count, zin(:,:,:))
                ALLOCATE(zin2(nlev,2,nsp))
                zin2(:,:,:) = 0.0_dp
                DO kk=ilev_d_min,ilev_d_max           
                   DO jj=1,nsp
                      zin2(kk,1,jj) = zin(1,jj,kk)
                      zin2(kk,2,jj) = zin(2,jj,kk)
                   ENDDO
                ENDDO
                DEALLOCATE(zin)
             END IF
             CALL scatter_sp(zin2,sdobs3,gdc)
             IF (p_parallel_io) DEALLOCATE(zin2)
          END IF   ! end of linp_div  

          !! read vorticity
          IF (linp_vor) THEN
             IF (p_parallel_io) THEN
                CALL IO_INQ_VARID (ndgfile%file_id, 'svo', nvarid)
                COUNT(:) = (/ 2, nsp,nlev,1 /)
                start(:) = (/ 1, 1, 1,ndtsc /)
                ALLOCATE(zin(2,nsp,nlev))
                zin(:,:,:) = 0.0_dp
                CALL IO_GET_VARA_DOUBLE (ndgfile%file_id, nvarid, start, count, zin(:,:,:))
                ALLOCATE(zin2(nlev,2,nsp))
                zin2(:,:,:) = 0.0_dp
                DO kk=ilev_d_min,ilev_d_max           
                   DO jj=1,nsp
                      zin2(kk,1,jj) = zin(1,jj,kk)
                      zin2(kk,2,jj) = zin(2,jj,kk)
                   ENDDO
                ENDDO
                DEALLOCATE(zin)
             END IF
             CALL scatter_sp(zin2,svobs3,gdc)
             IF (p_parallel_io) DEALLOCATE(zin2)
          END IF   ! end of linp_vor    
         
          !! read temperature
          IF (linp_tem) THEN
             IF (p_parallel_io) THEN
                CALL IO_INQ_VARID (ndgfile%file_id, 't', nvarid)
                COUNT(:) = (/ 2, nsp,nlev,1 /)
                start(:) = (/ 1, 1, 1,ndtsc /)
                ALLOCATE(zin(2,nsp,nlev))
                zin(:,:,:) = 0.0_dp
                CALL IO_GET_VARA_DOUBLE (ndgfile%file_id, nvarid, start, count, zin(:,:,:))
                ALLOCATE(zin2(nlev,2,nsp))
                zin2(:,:,:) = 0.0_dp
                IF (linp_lnp)       kkmax=ilev_t_min+ino_t_lev-2 
                   !! if we nudge ln of surface pressure: ilev_t_max = ilev_t_min+ino_t_lev-2 
                IF (.NOT.linp_lnp)  kkmax=ilev_t_min+ino_t_lev-1 
                   !! if we don't nudge ln of surface pressure: ilev_t_max = ilev_t_min+ino_t_lev-1
                DO kk=ilev_t_min,kkmax
                   DO jj=1,nsp 
                      zin2(kk,1,jj) = zin(1,jj,kk)
                      zin2(kk,2,jj) = zin(2,jj,kk)
                   ENDDO
                ENDDO
                DEALLOCATE(zin)
             END IF
             CALL scatter_sp(zin2,stobs3(1:nlev,:,:),gdc)
             IF (p_parallel_io) DEALLOCATE(zin2)
          END IF   ! end of linp_tem 

          !! read logarithm of surface pressure
          IF (linp_lnp) THEN
             IF (p_parallel_io) THEN
                CALL IO_INQ_VARID (ndgfile%file_id, 'lsp', nvarid)
                COUNT(:) = (/ 2, nsp,1,1 /)
                start(:) = (/ 1, 1, 1,ndtsc /)
                ALLOCATE(zin(2,nsp,1))
                zin(:,:,:) = 0.0_dp
                CALL IO_GET_VARA_DOUBLE (ndgfile%file_id, nvarid, start, count, zin(:,:,:))
                ALLOCATE(zin2(1,2,nsp))
                zin2(:,:,:) = 0.0_dp
                DO jj=1,nsp 
                   zin2(1,1,jj) = zin(1,jj,1)
                   zin2(1,2,jj) = zin(2,jj,1)
                ENDDO
                DEALLOCATE(zin)
             END IF
             CALL scatter_sp(zin2,stobs3(nlevp1:nlevp1,:,:),gdc)
             IF (p_parallel_io) DEALLOCATE(zin2)
          END IF   ! end of linp_lnp 
  
          IF (lnudgpat) THEN
             CALL finish('ReadOneBlock','lnudgpat does not work yet for nudging data in netcdf-format')      
          END IF
  
          lndgstep3 = .TRUE.

          WRITE(mess,'(a48,i0)') 'READ NEW nudging DATA, corresponding step is : ',ndgstep3
          CALL message('ReadOneBlock',mess)
          CALL write_date(ndg_date3,' Corresponding date: ')

          ! perform NMI filter now (optional) **********************************
          IF (lnmi .AND. lnudgfrd) CALL NMI_Make(NMI_MAKE_FOD)

          EXIT ! read block OBS3 HEAD3 success
    
       ENDDO

    CASE DEFAULT

       CALL finish('ReadOneBlock','inudgformat set incorrectly')
    
    END SELECT 
    
  CONTAINS
 
    SUBROUTINE convert_input(feld,lmin,lmax,lno)

      USE mo_control,       ONLY: nsp
      USE mo_transpose,     ONLY: scatter_sp
      USE mo_decomposition, ONLY: gdc => global_decomposition

      INTEGER :: lmin, lmax, lno, i
      REAL(kind=dp)    :: feld(:,:,:)
      REAL(kind=dp), POINTER :: feld_global(:,:,:)
 
      IF (p_parallel_io) THEN

        ALLOCATE(feld_global(lmin:lmax,2,nsp))
        i = 2*nsp*lno
#ifdef LITTLE_ENDIAN
        CALL swap64(phbuf,i)
#endif
        CALL util_cray2ieee(phbuf,zhbuf,i)
        feld_global(lmin:lmax,:,:) = RESHAPE(zhbuf,(/lno,2,nsp/))
    ENDIF

    CALL scatter_sp(feld_global,feld(lmin:lmax,:,:),gdc)

    IF (p_parallel_io) DEALLOCATE(feld_global)
 
    END SUBROUTINE convert_input

  END SUBROUTINE ReadOneBlock

  !=============================================================================
!BOP
  ! !IROUTINE: OpenOneBlock  
  ! !INTERFACE:

  SUBROUTINE OpenOneBlock

    ! !DESCRIPTION: 
    ! open nudging data block

    ! !USES:
    USE mo_mpi,               ONLY: p_parallel_io
    USE mo_nudging_constants, ONLY: linp_div, linp_vor, linp_tem, &
                                    linp_lnp, lnudgcli, &
                                    inudgformat
    USE mo_filename,          ONLY: str_filter
    USE mo_exception,         ONLY: finish, message
    USE mo_time_control,      ONLY: get_date_components, ndg_inp_date, init_nudgingtime_d
!EOP
    USE mo_filename,         ONLY: NETCDF            !! AV
    USE mo_io,               ONLY: io_read, io_open  !! AV


    INTEGER :: kret, yr, mo, dy, hr, mi, se
    LOGICAL :: found
    CHARACTER(len=300) :: cfile
    
 
    IF (files_open) THEN
       CALL CloseBlock
       CALL message('','Nudging data blockes were closed in open.')
       ! set new input date for file evaluation
       IF (lndgstep2) THEN
          CALL init_nudgingtime_d
       ELSE
          CALL finish('OpenOneBlock',&
               'Minimum of two dates needed for nudging date evaluation.')
       END IF
    END IF

    ndgblock = ndgblock + 1    ! open next data block
    IF (lnudgcli .AND. ndgblock == 4) ndgblock = 1

    ! get time components for file name expansion
    CALL get_date_components(ndg_inp_date,yr,mo,dy,hr,mi,se)


    SELECT CASE (inudgformat)
    
      CASE(0)                  ! read cray-data
       
        IF (linp_div) THEN

           cfile = TRIM(str_filter(ndg_file_div,yr,mo,dy,hr,mi,se,ndgblock))
           WRITE(mess,*) 'use DIV file : ',TRIM(cfile)
           CALL message('',mess)

           IF (p_parallel_io) THEN
              INQUIRE(file=cfile,exist=found)
              IF (.NOT.found) THEN
                 CALL finish('OpenOneBlock','Nudging data file not found.')
              ELSE
                 CALL pbopen(kfiled,cfile,'r',kret)
                 IF(kret/=0) CALL finish('OpenOneBlock','DATA file DIV empty?')
              END IF
           END IF
           files_open = .TRUE.

        END IF
    
        IF (linp_vor) THEN

           cfile = TRIM(str_filter(ndg_file_vor,yr,mo,dy,hr,mi,se,ndgblock))
           WRITE(mess,*) 'use VOR file : ',TRIM(cfile)
           CALL message('',mess)

           IF (p_parallel_io) THEN
              INQUIRE(file=cfile,exist=found)
              IF (.NOT.found) THEN
                 CALL finish('OpenOneBlock','Nudging data file not found.')
              ELSE
                 CALL pbopen(kfilev,cfile,'r',kret)
                 IF(kret/=0) CALL finish('OpenOneBlock','DATA file VOR empty?')
              END IF
           END IF
           files_open = .TRUE.

        END IF

        IF (linp_tem .OR. linp_lnp) THEN

           cfile = TRIM(str_filter(ndg_file_stp,yr,mo,dy,hr,mi,se,ndgblock))
           WRITE(mess,*) 'use STP file : ',TRIM(cfile)
           CALL message('',mess)

           IF (p_parallel_io) THEN
              INQUIRE(file=cfile,exist=found)
              IF (.NOT.found) THEN
                 CALL finish('OpenOneBlock','Nudging data file not found.')
              ELSE
                 CALL pbopen(kfilet,cfile,'r',kret)
                 IF(kret/=0) CALL finish('OpenOneBlock','DATA file TEM empty?')
              END IF
           END IF
           files_open = .TRUE.

        END IF

        WRITE(mess,*) 'Attention : convert NUDGE DATA to IEEE'
        CALL message('OpenOneBlock',mess)
    
      CASE(2)                  ! read netcdf-data
   
       IF (p_parallel_io) THEN
           
           cfile = str_filter(ndg_file_nc,yr,mo,dy,hr,mi,se,ndgblock)
           WRITE(mess,*) 'use nudging data file : ',TRIM(cfile)
           CALL message('',mess)

           INQUIRE(file=cfile,exist=found)
           IF (.NOT.found) THEN
              CALL finish('OpenOneBlock','Nudging netcdf data file not found.')
           ELSE
              ndgfile%format = NETCDF
              CALL IO_open (cfile, ndgfile, IO_READ)
           END IF
        END IF
        files_open = .TRUE.

      CASE default
       
        IF (p_parallel_io) THEN
           WRITE (mess,*) 'inudgformat =', inudgformat, ' in ndgctl is not supported'
           CALL message('OpenOneBlock',mess)
           CALL finish('OpenOneBlock','Run terminated because of invalid value of inudgformat.')
        END IF  
     
    END SELECT
    
    
  END SUBROUTINE OpenOneBlock


  !=============================================================================
!BOP
  ! !IROUTINE:  CloseBlock
  ! !INTERFACE:

  SUBROUTINE CloseBlock

    ! !DESCRIPTION: 
    ! close nudging data block

    ! !USES:
    USE mo_nudging_constants, ONLY: linp_div, linp_vor, linp_tem, linp_lnp, &
                                    inudgformat                 !! AV
    USE mo_mpi,               ONLY: p_parallel_io
    USE mo_exception,         ONLY: message, &
                                    finish                      !! AV 
    
    USE mo_io,                ONLY: IO_close                    !! AV
!EOP

    INTEGER  :: kret

    IF (files_open) THEN
       IF (p_parallel_io) THEN
         SELECT CASE (inudgformat)
           CASE(0)                 ! cray-data
             IF (linp_div)               CALL pbclose(kfiled,kret)
             IF (linp_vor)               CALL pbclose(kfilev,kret)
             IF (linp_tem .OR. linp_lnp) CALL pbclose(kfilet,kret)
           CASE(2)                 ! netcdf-data
             CALL IO_close(ndgfile)  
           CASE DEFAULT
             CALL finish('CloseBlock','Run terminated because of invalid value of inudgformat')
         END SELECT
       END IF
       files_open = .FALSE.
    ELSE
       CALL message('','No Block open, close skipped.')
    END IF

  END SUBROUTINE CloseBlock

END MODULE mo_nudging_io
