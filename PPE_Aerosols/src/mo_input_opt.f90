#ifdef __xlC__
@PROCESS STRICT
@PROCESS OPT(2)
#endif 
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! This module contains procedures for mo_input dealing with namelist options 
! Never directly use this module from outside the mo_input framework - only via mo_input.f90
MODULE mo_input_opt
USE mo_input_strings
USE mo_input_types
IMPLICIT NONE

PUBLIC :: InputOptNew, InputOptDisp, InputOptImport

CONTAINS
! ----------------------------------------------------
!  Creating, adding, setting, modifying, displaying:
!     External settings for files and variables (opt)
! ----------------------------------------------------

! New Opt list element eventually added to global chain
  SUBROUTINE InputOptNew(Opt,lChain)
  TYPE (input_opt_list), POINTER :: Opt
  LOGICAL, OPTIONAL, INTENT(in)  :: lChain

    ALLOCATE(Opt)
    NULLIFY(Opt%next)
    IF (PRESENT(lChain)) THEN
      IF (lChain) THEN
        Opt%next     => AllInputOpts ! Reverses order, but order doesn't matter
        AllInputOpts => Opt
      ENDIF
    ELSE
      Opt%next     => AllInputOpts ! Reverses order, but order doesn't matter
      AllInputOpts => Opt
    ENDIF

    Opt%sub_model       = CurrSubModel
    Opt%flags           = 0
    Opt%file_name       = ''
    Opt%dim_rename      = ''
    Opt%weights         = ''
    Opt%def_vec         = ''
    Opt%var_name        = ''
    Opt%var_file_name   = ''
    Opt%depend          = ''
    Opt%dt_update_unit  = ''
    Opt%dt_file_unit    = ''
    Opt%action_var_miss = ''
    Opt%action_bad_val  = ''
    Opt%action_time_err = ''
    Opt%action_cycle    = ''
    Opt%action_dim_mis  = ''
    Opt%action_data     = ''
    Opt%data_start      = ''
    Opt%data_end        = ''
    Opt%valid_time      = ''
    Opt%labstime        = ''
    Opt%lnametime       = ''
    Opt%dt_update       = UNLIKELY_VAL
    Opt%dt_file         = UNLIKELY_VAL 
    Opt%offset_year     = UNLIKELY_VAL
    Opt%offset_day      = UNLIKELY_VAL
    Opt%offset_sec      = UNLIKELY_VAL
    Opt%offset_rec      = UNLIKELY_VAL
    Opt%init_rec        = UNLIKELY_VAL
    Opt%cycle_length    = UNLIKELY_VAL
    Opt%subset_index    = UNLIKELY_VAL
    Opt%file_weight     = UNLIKELY_RVAL
    Opt%mul             = UNLIKELY_RVAL
    Opt%add             = UNLIKELY_RVAL
    Opt%def_val         = UNLIKELY_RVAL
    Opt%fill_val        = UNLIKELY_RVAL
    Opt%missing_val     = UNLIKELY_RVAL
    Opt%valid_range(:)  = UNLIKELY_RVAL

  END SUBROUTINE InputOptNew

! Displays one or all options in list
  SUBROUTINE InputOptDisp(opt,cnt,unt,msg)
  TYPE (input_opt_list), OPTIONAL, POINTER :: opt
  INTEGER, OPTIONAL, INTENT(in) :: cnt
  INTEGER, OPTIONAL, INTENT(in) :: unt
  CHARACTER (len=*), INTENT(in), OPTIONAL :: msg

  TYPE (input_opt_list), POINTER :: curr
  INTEGER :: un, i, cr

    un = GetLogUnit(); IF (PRESENT(unt)) un = unt
    cr = -1;           IF (PRESENT(cnt)) cr = cnt

    IF (PRESENT(opt)) THEN
      curr => opt
    ELSE
      curr => AllInputOpts
    ENDIF
    i = 1
    IF (PRESENT(msg)) WRITE(un,'(a)') msg
    DO WHILE (ASSOCIATED(curr) .AND. (cr /= 0))
      WRITE(message_text,*) 'External_input_settings_(',i,')'
      CALL RmSpc(message_text)
      WRITE(un,'(a)') TRIM(message_text)
      WRITE(un,'(2a)') '  Setting type                     : ',TRIM(GetStringIndexed(flg_opt_type, &
                                                                    IAND(curr%flags,INPUT_OPT_NAMELIST)+1))
      IF (LEN_TRIM(curr%var_name) > 0) THEN
        WRITE(un,'(2a)') '  Variable/group name              : ',TRIM(curr%var_name)
        WRITE(un,'(2a)') '  Variable name(s) in file(s)      : ',TRIM(curr%var_file_name)
        WRITE(un,'(2a)') '  Variable dependent on variables  : ',TRIM(curr%depend)
        message_text = 'default'
        IF (curr%dt_update/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%dt_update
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Update time step                 : ',TRIM(message_text)
        WRITE(un,'(2a)') '  Update time step unit            : ',TRIM(curr%dt_update_unit)
        WRITE(un,'(2a)') '  Validity time                    : ',TRIM(curr%valid_time)
        WRITE(un,'(2a)') '  Match model and file data time   : ',TRIM(curr%labstime)
        WRITE(un,'(2a)') '  Extract time from file names     : ',TRIM(curr%lnametime)
        WRITE(un,'(2a)') '  Action for mismatching time      : ',TRIM(curr%action_time_err)
        WRITE(un,'(2a)') '  Action for data out of range     : ',TRIM(curr%action_bad_val)
        message_text = 'default'
        IF (curr%offset_day/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%offset_day
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Offset days                      : ',TRIM(message_text)
        message_text = 'default'
        IF (curr%offset_sec/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%offset_sec
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Offset seconds                   : ',TRIM(message_text)
        WRITE(un,'(2a)') '  Data cycling action              : ',TRIM(curr%action_cycle)
        message_text = 'default'
        IF (curr%cycle_length/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%cycle_length
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Cycle length                     : ',TRIM(message_text)
        WRITE(un,'(2a)') '  Action for missing variable      : ',TRIM(curr%action_var_miss)
        message_text = 'none'
        IF (curr%def_val/=UNLIKELY_RVAL) THEN
          WRITE(message_text,*) curr%def_val
          CALL RmSpc(message_text)
        ENDIF
        IF (LEN_TRIM(curr%def_vec) > 0) THEN
          message_text = TRIM(curr%def_vec)
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Default value for missing var    : ',TRIM(message_text)
        message_text = 'none'
        IF (curr%missing_val/=UNLIKELY_RVAL) THEN
          WRITE(message_text,*) curr%missing_val
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Value for missing data in file   : ',TRIM(message_text)
        message_text = 'none'
        IF (curr%valid_range(1)/=UNLIKELY_RVAL) THEN
          WRITE(message_text,*) curr%valid_range
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Valid range for file data (lo,hi): ',TRIM(message_text)
        message_text = 'none'
        IF (curr%fill_val/=UNLIKELY_RVAL) THEN
          WRITE(message_text,*) curr%fill_val
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Value to replace missing values  : ',TRIM(message_text)
        WRITE(un,'(2a)') '  Action for mismatching dimension : ',TRIM(curr%action_dim_mis)
        WRITE(un,'(2a)') '  Action for mismatching dimension : ',TRIM(curr%action_dim_mis)
        WRITE(un,'(2a)') '  Extra dimension weights          : ',TRIM(curr%weights)
        message_text = 'default'
        IF (curr%subset_index/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%subset_index
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Index for subsetting/-reading    : ',TRIM(message_text)
        message_text = 'default'
        IF (curr%mul/=UNLIKELY_RVAL) THEN
          WRITE(message_text,*) curr%mul
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Data multiplication scale        : ',TRIM(message_text)
        message_text = 'default'
        IF (curr%add/=UNLIKELY_RVAL) THEN
          WRITE(message_text,*) curr%add
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Data additional offset           : ',TRIM(message_text)
        WRITE(un,'(2a)') '  Input/model data combining action: ',TRIM(curr%action_data)
        message_text = 'default (1)'
        IF (curr%file_weight/=UNLIKELY_RVAL) THEN
          WRITE(message_text,*) curr%file_weight
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Weighting factor for file data   : ',TRIM(message_text)
      ENDIF
      IF (LEN_TRIM(curr%file_name) > 0) THEN
        WRITE(un,'(2a)') '  File name(s)                     : ',TRIM(curr%file_name)
        message_text = 'default'
        IF (curr%dt_file/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%dt_file
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  File time step                   : ',TRIM(message_text)
        WRITE(un,'(2a)') '  File time step unit              : ',TRIM(curr%dt_file_unit)
        WRITE(un,'(2a)') '  Data start time                  : ',TRIM(curr%data_start)
        WRITE(un,'(2a)') '  Data end   time                  : ',TRIM(curr%data_end)
        message_text = 'default'
        IF (curr%offset_year/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%offset_year
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Offset years                     : ',TRIM(message_text)
        message_text = 'default'
        IF (curr%offset_rec/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%offset_rec
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Offset records                   : ',TRIM(message_text)
        message_text = 'default'
        IF (curr%init_rec/=UNLIKELY_VAL) THEN
          WRITE(message_text,*) curr%init_rec
          CALL RmSpc(message_text)
        ENDIF
        WRITE(un,'(2a)') '  Record for initial variables     : ',TRIM(message_text)
        WRITE(un,'(2a)') '  Dimension rename                 : ',TRIM(curr%dim_rename)
      ENDIF

      i = i + 1
      curr => curr%next
      cr = cr - 1
    ENDDO

  END SUBROUTINE InputOptDisp

! Obtain settings for input variables/files from namelist file
  SUBROUTINE InputOptImport(Files,Unt)
  CHARACTER (len=*), INTENT(in) :: Files
  INTEGER, OPTIONAL, INTENT(in) :: Unt

  ! Namelist definition and variables
  CHARACTER (len=128) :: file_name, var_file_name, dim_rename, weights, def_vec, depend
  CHARACTER (len=64)  :: var_name
  CHARACTER (len=32)  :: action_dim_mis, action_data
  CHARACTER (len=16)  :: dt_update_unit, dt_file_unit, labstime, lnametime, valid_time, &
                         action_var_miss, action_time_err, action_cycle, action_bad_val
  CHARACTER (len=13)  :: data_start, data_end
  INTEGER             :: dt_update, dt_file, offset_year, offset_day, offset_sec, offset_rec, init_rec, cycle_length, subset_index
  REAL(dp)            :: mul, add, file_weight, def_val, fill_val, missing_val, valid_range(2)
  NAMELIST /input_ctl/                                                                                                         &
    file_name, dim_rename, var_name, var_file_name, dt_update_unit, dt_file_unit, action_var_miss, action_time_err,            &
    action_cycle,action_dim_mis,action_data,action_bad_val,data_start,data_end,valid_time,subset_index,dt_update,dt_file,      &
    offset_year, offset_day, offset_sec, offset_rec, init_rec, cycle_length, labstime, lnametime, mul, add, def_val, fill_val, &
    missing_val, valid_range, file_weight, weights, def_vec, depend

  ! Other local variables
  TYPE (input_opt_list), POINTER :: Opt, prev
  INTEGER :: ActStart, ActEnd, RUnit, iErr, un
  LOGICAL :: lCont, lCont2, lOpened

    RUnit = 10
    IF (PRESENT(Unt)) RUnit = Unt

    NULLIFY(prev, Opt)
    iErr = 0
    ActStart = 1
    lCont = .TRUE.
    DO WHILE (lCont)

      ! Determine file/unit from which to read namelist(s)
      IF (LEN_TRIM(Files) > 0) THEN
        ActEnd = INDEX(Files(ActStart:),',') + ActStart - 1
        IF (ActEnd < ActStart) THEN 
          ActEnd = LEN_TRIM(Files) + 1
          lCont = .FALSE.
        ENDIF
        iErr = 0
        lOpened = .FALSE.
        IF (lio) THEN
          INQUIRE(file=Files(ActStart:ActEnd-1),number=un,opened=lCont2) ! Test if file is already opened since reopening is an
          IF (lCont2) THEN                                               ! error on some systems
            RUnit = un
            REWIND(RUnit)
          ELSE
            OPEN(unit=RUnit,file=Files(ActStart:ActEnd-1),iostat=iErr,status='old',action='read',delim='apostrophe')
            lOpened = .TRUE.
          ENDIF
        ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
        CALL p_bcast(iErr,io_pe)
#endif
! ICON
        IF (iErr /= 0) CALL local_message('InputOptImport','Namelist file '//Files(ActStart:ActEnd-1)//' could not be opened')
      ELSE
        lCont = .FALSE.
      ENDIF

      IF (iErr == 0) THEN
        lCont2 = .TRUE.
        DO WHILE (lCont2) ! Loop over namelists

          CALL InputOptNew(Opt,.FALSE.)

          ! Set default values
          file_name       = TRIM(Opt%file_name)
          dim_rename      = TRIM(Opt%dim_rename)
          weights         = TRIM(Opt%weights)
          def_vec         = TRIM(Opt%def_vec)
          var_name        = TRIM(Opt%var_name)
          var_file_name   = TRIM(Opt%var_file_name)
          dt_update_unit  = TRIM(Opt%dt_update_unit)
          dt_file_unit    = TRIM(Opt%dt_file_unit)
          depend          = TRIM(Opt%depend)
          action_bad_val  = TRIM(Opt%action_bad_val)
          action_var_miss = TRIM(Opt%action_var_miss)
          action_dim_mis  = TRIM(Opt%action_dim_mis)
          action_time_err = TRIM(Opt%action_time_err)
          action_cycle    = TRIM(Opt%action_cycle)
          action_data     = TRIM(Opt%action_data)
          data_start      = TRIM(Opt%data_start)
          data_end        = TRIM(Opt%data_end)
          dt_update       = Opt%dt_update
          dt_file         = Opt%dt_file
          offset_year     = Opt%offset_year
          offset_day      = Opt%offset_day
          offset_sec      = Opt%offset_sec
          offset_rec      = Opt%offset_rec
          init_rec        = Opt%init_rec
          cycle_length    = Opt%cycle_length
          valid_time      = Opt%valid_time
          subset_index    = Opt%subset_index
          labstime        = Opt%labstime
          lnametime       = Opt%lnametime
          mul             = Opt%mul
          add             = Opt%add
          file_weight     = Opt%file_weight
          def_val         = Opt%def_val
          fill_val        = Opt%fill_val
          missing_val     = Opt%missing_val
          valid_range(:)  = Opt%valid_range(:) 

          ! Read namelist
          iErr = 1
          DO WHILE (iErr >= 0) ! Some compilers return iostat > 0 for "not a namelist/not requested namelist",
            IF (lio) &
              READ(Unit=RUnit,nml=input_ctl,iostat=iErr) ! all return -1 for EOF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
            CALL p_bcast(iErr,io_pe)
#endif            
            IF (iErr == 0) THEN ! Namelist found and read
              iErr = -100000 ! Get out of this loop to initialize new Opt

              IF ((IAND(mo_debug,INPUT_DBG_EXTERN) /= 0) .AND. lio) THEN
                IF (var_name /= '') THEN
                  CALL local_message('InputOptImport', &
                  'Imported input definitions from '//Files(ActStart:ActEnd-1)//' for variable(s) '//TRIM(var_name))
                ELSE
                  CALL local_message('InputOptImport', &
                  'Imported input definitions from '//Files(ActStart:ActEnd-1)//' for file(s) '//TRIM(file_name))
                ENDIF
              ENDIF
#if defined(INPUT_IN_ECHAM) || defined(USE_MPI)
              CALL p_bcast(file_name,io_pe)
              CALL p_bcast(dim_rename,io_pe)
              CALL p_bcast(weights,io_pe)
              CALL p_bcast(def_vec,io_pe)
              CALL p_bcast(depend,io_pe)
              CALL p_bcast(var_name,io_pe)
              CALL p_bcast(var_file_name,io_pe)
              CALL p_bcast(dt_update_unit ,io_pe)
              CALL p_bcast(dt_file_unit,io_pe)
              CALL p_bcast(action_bad_val,io_pe)
              CALL p_bcast(action_var_miss,io_pe)
              CALL p_bcast(action_dim_mis,io_pe)
              CALL p_bcast(action_time_err,io_pe)
              CALL p_bcast(action_cycle,io_pe)
              CALL p_bcast(action_data,io_pe)
              CALL p_bcast(data_start,io_pe)
              CALL p_bcast(data_end,io_pe)
              CALL p_bcast(dt_update,io_pe)
              CALL p_bcast(dt_file,io_pe)
              CALL p_bcast(offset_year,io_pe)
              CALL p_bcast(offset_day,io_pe)
              CALL p_bcast(offset_sec,io_pe)
              CALL p_bcast(offset_rec,io_pe)
              CALL p_bcast(init_rec,io_pe)
              CALL p_bcast(cycle_length,io_pe)
              CALL p_bcast(valid_time,io_pe)
              CALL p_bcast(subset_index,io_pe)
              CALL p_bcast(labstime,io_pe)
              CALL p_bcast(lnametime,io_pe)
              CALL p_bcast(mul,io_pe)
              CALL p_bcast(add,io_pe)
              CALL p_bcast(def_val,io_pe)
              CALL p_bcast(fill_val,io_pe)
              CALL p_bcast(missing_val,io_pe)
              CALL p_bcast(valid_range,io_pe)
              CALL p_bcast(file_weight,io_pe)
              CALL p_bcast(iErr,io_pe)
#endif
              Opt%file_name       = TRIM(file_name)
              Opt%dim_rename      = TRIM(dim_rename)
              Opt%weights         = TRIM(weights)
              Opt%def_vec         = TRIM(def_vec)
              Opt%var_name        = TRIM(var_name)
              Opt%var_file_name   = TRIM(var_file_name)
              Opt%depend          = TRIM(depend)
              Opt%dt_update_unit  = TRIM(dt_update_unit)
              Opt%dt_file_unit    = TRIM(dt_file_unit)
              Opt%action_bad_val  = TRIM(action_bad_val)
              Opt%action_var_miss = TRIM(action_var_miss)
              Opt%action_dim_mis  = TRIM(action_dim_mis)
              Opt%action_time_err = TRIM(action_time_err)
              Opt%action_cycle    = TRIM(action_cycle)
              Opt%action_data     = TRIM(action_data)
              Opt%data_start      = TRIM(data_start)
              Opt%data_end        = TRIM(data_end)
              Opt%dt_update       = dt_update
              Opt%dt_file         = dt_file
              Opt%offset_year     = offset_year
              Opt%offset_day      = offset_day
              Opt%offset_sec      = offset_sec
              Opt%offset_rec      = offset_rec
              Opt%init_rec        = init_rec
              Opt%cycle_length    = cycle_length
              Opt%valid_time      = valid_time
              Opt%subset_index    = subset_index
              Opt%labstime        = labstime
              Opt%lnametime       = lnametime
              Opt%mul             = mul
              Opt%add             = add
              Opt%file_weight     = file_weight
              Opt%def_val         = def_val
              Opt%fill_val        = fill_val
              Opt%missing_val     = missing_val
              Opt%valid_range(:)  = valid_range(:)
              Opt%flags           = INPUT_OPT_NAMELIST

              ! Save this VarOpt in global list with the read values
              IF (ASSOCIATED(prev)) THEN
                prev%next    => Opt
              ELSE
                AllInputOpts => Opt
              ENDIF
              prev => Opt
            ENDIF
          ENDDO
          !>>HK #532 better support for some platforms, to avoid unfinite loop
          !orig lCont2 = iErr /= -1
          lCont2 = iErr > -1
          !<<HK
        ENDDO ! Successful namelist read
        IF ((LEN_TRIM(Files) > 0) .AND. lio .AND. lOpened) CLOSE(Unit=RUnit)
      ENDIF
      ActStart = ActEnd + 1

    ENDDO ! File loop

    IF (ASSOCIATED(Opt)) DEALLOCATE(Opt) ! Get rid of the last one, which didn't get used

  END SUBROUTINE InputOptImport

END MODULE mo_input_opt
