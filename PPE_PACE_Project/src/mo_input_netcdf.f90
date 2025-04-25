!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! To be used from in mo_input.f90
! Do not include this file directly.
! Functions to specific to handle NetCDF files

! V.1.00, 2014-01-02, Stiig Wilkenskjeld, MPI-M, Extracted from mo_input itself
! V.2.00, 2014-08-29, Stiig Wilkenskjeld, MPI-M, Converted to a separate module

MODULE mo_input_netcdf
USE mo_input_calendar
USE mo_input_strings
USE mo_input_types
USE mo_input_base
USE mo_input_fileobj
IMPLICIT NONE

PUBLIC :: NCFileOpen, NCFileClose, NCFileGetNRec, NCFileRead
PUBLIC :: NCFileGetDimensions, NCFileGetVariables, NCFileGetRecTime, NCFileGetAttr

CONTAINS
! Open this file
  INTEGER FUNCTION NCFileOpen(file_name,message)
  INCLUDE 'netcdf.inc'
  CHARACTER (len=*), INTENT(in)            :: file_name
  CHARACTER (len=*), INTENT(out), OPTIONAL :: message

  INTEGER fid, st

    st = nf_open(file_name,nf_nowrite,fid)
    NCFileOpen = fid
    IF (st == NF_NOERR) RETURN
    NCFileOpen = -1
    IF (PRESENT(message)) THEN
      message = 'NetCDF error: '//TRIM(nf_strerror(st))
    ELSE
      message_text = 'NetCDF error: '//TRIM(nf_strerror(st))
    ENDIF

  END FUNCTION NCFileOpen

! Close a file
  LOGICAL FUNCTION NCFileClose(fid)
  INCLUDE 'netcdf.inc'
  INTEGER, INTENT(in) :: fid

  INTEGER :: st

    message_text = ''
    st = nf_close(fid)
    IF (st /= NF_NOERR) WRITE(message_text,*) '  Error closing file. Netcdf error: ',TRIM(nf_strerror(st))
    NCFileClose = st==NF_NOERR

  END FUNCTION NCFileClose

! Get the number of records in a file
  INTEGER FUNCTION NCFileGetNRec(IFile,udim)
  INCLUDE 'netcdf.inc'
  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_dim_list),  POINTER, OPTIONAL :: udim

  INTEGER st, res

    IF (PRESENT(udim)) THEN
      st = nf_inq_dimid(IFile%fid,TRIM(udim%dim_data%name_dim),IFile%time_fid(1))
      IF (st /= NF_NOERR) &
        st = nf_inq_dimid(IFile%fid,TRIM(udim%dim_data%name_alt),IFile%time_fid(1))
      IF (st /= NF_NOERR) THEN
        NCFileGetNRec = -1
        RETURN
      ENDIF
    ENDIF

    st = nf_inq_dimlen(IFile%fid,IFile%time_fid(1),res)
    IF (st /= NF_NOERR) THEN
      WRITE(message_text,*) '  Could not obtain dimension size for dimension: ',IFile%time_fid(1),'. Netcdf error: ', &
      TRIM(nf_strerror(st))
      CALL local_error('NCFileGetNRec',TRIM(IFile%name_file)//'. '//TRIM(message_text))
    ENDIF
    NCFileGetNRec = res

  END FUNCTION NCFileGetNRec

  ! Read variable data from a file
  LOGICAL FUNCTION NCFileRead(IFile,var,fid,ndim,Info,Ofs,RecNo,Buf)
  INCLUDE 'netcdf.inc'
  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_var_list),  POINTER :: var
  INTEGER,  INTENT(in)            :: fid 
  INTEGER,  INTENT(in)            :: ndim
  INTEGER,  INTENT(in)            :: Info(:)
  INTEGER,  INTENT(in)            :: Ofs
  INTEGER,  INTENT(in)            :: RecNo
  REAL(dp), INTENT(out)           :: Buf(*)

  INTEGER st, RSt(8), RSz(8), add

    NCFileRead = .TRUE.
    IF (.NOT. liodata) RETURN

    ! Find number of dimensions of variable in file (is this correct in all cases?)
    add = 0
    ! Starting reading point
    RSt(:) = 1
    RSt(1:ndim) = Info(1:ndim)
    IF (IAND(var%stat,INPUT_VAR_HAS_UNLIM) /= 0) THEN
      RSt(ndim+1) = RecNo
      add = 1
    ENDIF
    RSz(:) = 1
    RSz(1:ndim) = Info(ndim+1:2*ndim)
    ! Read!
    message_text = ''
    st = nf_get_vara_double(IFile%fid,fid,RSt(1:ndim+add),RSz(1:ndim+add),Buf(Ofs))
    IF (st /= NF_NOERR) THEN
      WRITE(message_text,*) '  Could not read data for variable: ',TRIM(var%name_var),', record: ',RecNo,'. Netcdf error: ', &
      TRIM(nf_strerror(st))
      NCFileRead = .FALSE.
    ENDIF

  END FUNCTION NCFileRead

! Check which dimensions are present in a file and set/check extends of these
  FUNCTION NCFileGetDimensions(IFile,Level)
  INCLUDE 'netcdf.inc'
  TYPE (input_fileobj_list), POINTER  :: NCFileGetDimensions
  TYPE (input_file_list),    POINTER  :: IFile
  INTEGER, INTENT(in),       OPTIONAL :: Level ! 1: names(def), 2: size(def), 4: direction, 8: data, 16: Id(def), 32: Unlim data

  TYPE (input_fileobj_list), POINTER  :: res, curr, new
  INTEGER :: flags, st, vid, n, i, unlim
  REAL(dp), DIMENSION(:), POINTER :: Dta

    NULLIFY(NCFileGetDimensions)
    message_text = ''

    flags = 19 ! 1+2+16
    IF (PRESENT(Level)) flags = Level

    st = nf_inq_ndims(IFile%fid,n)
    IF (st /= NF_NOERR) THEN
      WRITE(message_text,*) '  Could not obtain number of file dimensions. Netcdf error: '//TRIM(nf_strerror(st))
      CALL local_error('NCFileGetDimensions',TRIM(IFile%name_file)//'. '//TRIM(message_text))
    ENDIF
    st = nf_inq_unlimdim(IFile%fid,unlim)
    IF (st /= NF_NOERR) THEN
      WRITE(message_text,*) '  Could not obtain id of the unlimited dim. Netcdf error: '//TRIM(nf_strerror(st))
      CALL local_error('NCFileGetDimensions',TRIM(IFile%name_file)//'. '//TRIM(message_text))
    ENDIF

    NULLIFY(res,curr,new)    
    DO i=1,n

      ! Allocate, initialize and put in chain
      CALL InputFileObjNew(new)
      IF (ASSOCIATED(curr)) THEN
        curr%next => new
      ELSE
        res => new
      ENDIF
      curr => new
      new%flags = INPUT_FILEOBJ_DIM

      ! Get name/size if requested
      IF (IAND(flags,15) /= 0) THEN
        ! We use that the NetCDF_dim_ids are enumerated 1, 2 etc. Though ugly, this is a documented feature of the NetCDFs it
        ! therefore shouldn't be a big problem. Since NetCDF dimensions are principally unordered, we need either to scan the
        ! dimensions this way or search for specific dimension names, which we also don't want
        st = nf_inq_dim(IFile%fid,i,new%name_obj,new%size_global)
        IF (st /= NF_NOERR) THEN
          WRITE(message_text,*) '  Could not obtain dimension info for dimension: ',i,'. Netcdf error: '//TRIM(nf_strerror(st))
          CALL local_error('NCFileGetDimensions',TRIM(IFile%name_file)//'. '//TRIM(message_text))
        ENDIF
        new%Id = i
      ENDIF

      ! Get data/direction if required
      IF (IAND(flags,12) /= 0) THEN
        st = nf_inq_varid(IFile%fid,TRIM(new%name_obj),vid)
        IF (st /= NF_NOERR) THEN ! This may or may not be a fatal error thus ...
          WRITE(message_text,*) '  Dimension ',TRIM(new%name_obj),' has no associated data. Netcdf error: '//TRIM(nf_strerror(st))
          ! ... in this case we continue processing without data/direction, generating an error message, but no error
          new%flags = IOR(new%flags,INPUT_FILEDIM_DIR_ERROR)
        ELSE
          IF ((i == unlim) .AND. (IAND(flags,32) == 0)) THEN ! Don't get the strided data for unlimited dim unless expl. req.
            new%flags = IOR(new%flags,INPUT_FILEDIM_DIR_UNLIM)
          ELSE
            st = nf_inq_varndims(IFile%fid,vid,n)
            IF (st /= NF_NOERR) THEN
              WRITE(message_text,*) '  #Dimensions of dimension variable ',TRIM(new%name_obj),' could not be obtained.'
              new%flags = IOR(new%flags,INPUT_FILEDIM_DIR_ERROR)
            ELSE
              if (n > 1) THEN
                WRITE(message_text,*) '  Dimensions variable ',TRIM(new%name_obj),' has more than one dimension, reading skipped.'
                new%flags = IOR(new%flags,INPUT_FILEDIM_DIR_ERROR)
              ELSE
                ALLOCATE(Dta(new%size_global))
                st = nf_get_var_double(IFile%fid,vid,Dta)
                IF (st /= NF_NOERR) THEN ! This may or may not be a fatal error thus ...
                  WRITE(message_text,*) '  Dimension ',TRIM(new%name_obj),' data unreadable. Netcdf error: '//TRIM(nf_strerror(st))
                  ! ... in this case we continue processing without data/direction, generating an error message, but no error
                  new%flags = IOR(new%flags,INPUT_FILEDIM_DIR_ERROR)
                ELSE
                  ! Direction requested
                  IF (IAND(flags,4) /= 0) THEN ! We assume that data are monotonic and non-constant
                    IF (Dta(1) > Dta(2)) THEN
                      new%flags = IOR(new%flags,INPUT_FILEDIM_DIR_DEC)
                    ELSE
                      new%flags = IOR(new%flags,INPUT_FILEDIM_DIR_INC)
                    ENDIF
                  ENDIF
                  ! Values of the dimension variable requested
                  IF (IAND(flags,8) /= 0) THEN
                    new%points => Dta
                  ELSE
                    DEALLOCATE(Dta)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF

    ENDDO
    IFile%stat = IOR(IFile%stat,INPUT_FILE_DIMS_GOT)

    NCFileGetDimensions => Res

  END FUNCTION NCFileGetDimensions

! Check which variables are present in a file
  FUNCTION NCFileGetVariables(IFile)
  INCLUDE 'netcdf.inc'
  TYPE (input_fileobj_list), POINTER :: NCFileGetVariables
  TYPE (input_file_list),    POINTER :: IFile

  TYPE (input_fileobj_list), POINTER  :: res, curr, new
  INTEGER  :: st, typ, natt, n, i

    NULLIFY(NCFileGetVariables)
    message_text = ''

    IF (.NOT. liodata) RETURN

    st = nf_inq_nvars(IFile%fid,n)
    IF (st /= NF_NOERR) THEN
      WRITE(message_text,*) '  Could not obtain number of variables in file. Netcdf error: '//TRIM(nf_strerror(st))
      CALL local_error('NCFileGetVariables',TRIM(IFile%name_file)//'. '//TRIM(message_text))
    ENDIF

    NULLIFY(res,curr,new)
    DO i=1,n

      ! Allocate, initialize and put in chain
      CALL InputFileObjNew(new)
      IF (ASSOCIATED(curr)) THEN
        curr%next => new
      ELSE
        res => new
      ENDIF
      curr => new
      new%flags = INPUT_FILEOBJ_VAR

      ! Get name and dimension ids (also type and number of attributes is inquired to avoid multiple calls - not saved)
      !   We use that the NetCDF_var_ids are enumerated 1, 2 etc. Though ugly, this is a documented feature of the NetCDFs it
      !   therefore shouldn't be a big problem. Since NetCDF variables are principally unordered, we need either to scan the
      !   variables this way or search for specific variable names, which we also don't want
      st = nf_inq_var(IFile%fid,i,new%name_obj,typ,new%size_global,new%fids_dim,natt)
      IF (st /= NF_NOERR) THEN
        WRITE(message_text,*) '  Could not obtain info for variable nr.: ',i,'. Netcdf error: '//TRIM(nf_strerror(st))
        CALL local_error('NCFileGetVariables',TRIM(IFile%name_file)//'. '//TRIM(message_text))
      ENDIF
      new%Id = i

      ! Process relevant attributes of the variable
      ALLOCATE(new%points(5))
      new%points(:) = (/1._dp,0._dp,0._dp,0._dp,0._dp/)
      IF ((typ /= NF_FLOAT) .AND. (typ /= NF_DOUBLE)) THEN
        st = nf_get_att_double(IFile%fid,i,"scale_factor",new%points(1))
        IF (st /= NF_NOERR) THEN
          st = nf_get_att_double(IFile%fid,i,"mul_scale",new%points(1))
          IF (st /= NF_NOERR) new%points(1) = 1._dp
        ENDIF
        st = nf_get_att_double(IFile%fid,i,"add_offset",new%points(2))
        IF (st /= NF_NOERR) new%points(2) = 0._dp
      ENDIF
      st = nf_get_att_double(IFile%fid,i,"_FillValue",new%points(3))
      IF (st /= NF_NOERR) THEN
        st = nf_get_att_double(IFile%fid,i,"missing_value",new%points(3))
        IF (st /= NF_NOERR) new%points(3) = REAL(UNLIKELY_VAL,dp) / 1000._dp
      ENDIF
      st = nf_get_att_double(IFile%fid,i,"unpacked_valid_range",new%points(4))
      IF (st /= NF_NOERR) THEN
        st = nf_inq_atttype(IFile%fid,i,"valid_range",typ)
        IF ((st == NF_NOERR) .AND. ((typ /= NF_FLOAT) .AND. (typ /= NF_DOUBLE))) st = NF_NOERR - 1
        IF (st == NF_NOERR) st = nf_get_att_double(IFile%fid,i,"valid_range",new%points(4))
        IF (st /= NF_NOERR) new%points(4:5) = REAL(UNLIKELY_VAL,dp) / 1000._dp
      ENDIF

      IF ((new%points(1) /= 1._dp) .OR. (new%points(2) /= 0._dp)) new%flags = IOR(new%flags,INPUT_FILEVAR_HAS_SCALE)
      IF (new%points(3) /= REAL(UNLIKELY_VAL,dp) / 1000._dp) new%flags = IOR(new%flags,INPUT_FILEVAR_HAS_BAD)
      IF (new%points(4) /= REAL(UNLIKELY_VAL,dp) / 1000._dp) new%flags = IOR(new%flags,INPUT_FILEVAR_HAS_RANGE)

      IF (IAND(new%flags,INPUT_FILEVAR_MASK) == 0) DEALLOCATE(new%points)

    ENDDO
    IFile%stat = IOR(IFile%stat,INPUT_FILE_VARS_GOT)

    NCFileGetVariables => Res

  END FUNCTION NCFileGetVariables

! Get validity time of specified record in file
  FUNCTION NCFileGetRecTime(IFile,RecNo,lStep)
  INTEGER :: NCFileGetRecTime(4)
  TYPE (input_file_list), POINTER :: IFile
  INTEGER, INTENT(in)             :: RecNo
  LOGICAL, INTENT(in), OPTIONAL   :: lStep

  CHARACTER (len=64) :: Units
  INTEGER  :: nRec, st
  INTEGER  :: TmpTime(4,2)
  INTEGER(i8) :: df
  REAL(dp) :: Tmp(2)
  LOGICAL :: dStep

  INCLUDE 'netcdf.inc'

    ! Default error return
    NCFileGetRecTime = -1

    dStep = .FALSE.
    IF (PRESENT(lStep)) dStep = lStep

    IF (IAND(IFile%flags,INPUT_FILE_MULTI+INPUT_FILE_TIME_FROM_NAME) == INPUT_FILE_TIME_FROM_NAME) &
      CALL local_error('NCFileGetRecTime',TRIM(IFile%name_file)//'. Extracting time from file names only works for multi-files')

    IF ((IAND(IFile%flags,INPUT_FILE_TIME_FROM_NAME) /= 0) .AND. .NOT. dStep) THEN
      TmpTime(:,1) = InputFileTimeCurrGet(IFile,IFile%time_open)
      IF (RecNo > 1) THEN
        IF ((IFile%dt_file == 0) .OR. (IFile%dt_file == UNLIKELY_VAL)) &
          CALL local_error('NCFileGetRecTime',TRIM(IFile%name_file)//'. Determine record time step first')
        CALL CalendarTimeAdd(TmpTime(:,1),INT(IFile%dt_file,i8)*INT(RecNo-1,i8))
      ENDIF
      NCFileGetRecTime = TmpTime(:,1)
      RETURN
    ENDIF

    IF (IFile%time_fid(2) <= 0) CALL local_error('NCFileGetRecTime',TRIM(IFile%name_file)// &
      '. File does not contain the unlimited dimension. Was it intended as an initial file?')
 
    ! Process or reprocess unit string only if needed
    IF ((IAND(IFile%stat,INPUT_FILE_TIME_GOT)==0) .OR. (LEN_TRIM(IFile%time_unit) == 0)) THEN

      ! Get unit of time from file as string
      Units=''
      st = nf_get_att_text(IFile%fid,IFile%time_fid(2),'units',Units)
      IF (st /= NF_NOERR) THEN
        WRITE(message_text,*) '  Could not obtain unit for variable nr.: ',IFile%time_fid(2),'. Netcdf error: ', &
                              TRIM(nf_strerror(st))
        CALL local_error('NCFileGetRecTime',TRIM(IFile%name_file)//'. '//TRIM(message_text))
      ENDIF

      ! New string - reprocess!
      IF (TRIM(Units) /= TRIM(IFile%time_unit)) CALL GetTimeScaleFromStr(IFile,Units)
    ENDIF ! Get new unit string

    ! Get time data for requested record(s) from file
    nRec = 1
    IF (dStep) nRec = 2
    st = nf_get_vara_double(IFile%fid,IFile%time_fid(2),(/RecNo/),(/nRec/),Tmp)
    IF (st /= NF_NOERR) THEN
      WRITE(message_text,*) ' Could not read time data from record: ',RecNo,' of variable nr.: ',IFile%time_fid(2), &
                            '. Netcdf error: ',TRIM(nf_strerror(st))
      CALL local_error('NCFileGetRecTime',TRIM(IFile%name_file)//'. '//TRIM(message_text))
    ENDIF

    ! Process time(s)
    TmpTime(:,1) = GetRecTime(IFile,Tmp(1))

    IF (dStep) TmpTime(:,2) = GetRecTime(IFile,Tmp(2))

    ! Return value (record time or time between records)
    IF (nRec==2) THEN
      NCFileGetRecTime(1:3) = 0
      df = CalendarTimeDiff(TmpTime(:,1),TmpTime(:,2))
      IF (df > GREGORIAN_MIN_MONTH) df = -INT((REAL(df,dp) / REAL(GREGORIAN_MEAN_MONTH,dp)) + 0.5_dp,i8)
      NCFileGetRecTime(4) = INT(df)
    ELSE
      NCFileGetRecTime(:) = TmpTime(:,1)
    ENDIF

  END FUNCTION NCFileGetRecTime

  FUNCTION NCFileGetAttr(Attr,IFile,var,stat)
  CHARACTER(len=256)           :: NCFileGetAttr
  CHARACTER(len=*), INTENT(in) :: Attr
  TYPE (input_file_list), POINTER, OPTIONAL :: IFile
  TYPE (input_var_list),  POINTER, OPTIONAL :: var
  INTEGER,          INTENT(out),   OPTIONAL :: stat

  CHARACTER(len=256) :: tmp
  INTEGER fid, vid, st, tp, len

  INCLUDE 'netcdf.inc'

    NCFileGetAttr = ''
    IF (PRESENT(stat)) stat = 0
    vid = -1
    IF (PRESENT(var)) THEN
      IF (ASSOCIATED(var)) THEN
        vid = var%fid
        IF (.NOT. ASSOCIATED(var%parent_file)) RETURN
        fid = var%parent_file%fid
      ENDIF
    ENDIF
    IF (PRESENT(IFile) .AND. (vid < 0)) THEN
      fid = IFile%fid
      vid = NF_GLOBAL
    ENDIF
    st = nf_inq_att(fid,vid,Attr,tp,len)
    IF (st /= NF_NOERR) THEN
      IF (PRESENT(stat)) stat = st
      RETURN
    ENDIF
    IF (tp /= nf_char) THEN
      IF (PRESENT(stat)) stat = 1
      RETURN
    ENDIF
    IF (len > 256) THEN
      IF (PRESENT(stat)) stat = 2
      RETURN
    ENDIF
    tmp = ''
    st = nf_get_att_text(fid,vid,Attr,tmp)
    IF (st /= NF_NOERR) THEN
      IF (PRESENT(stat)) stat = st
      RETURN
    ENDIF
    st=256
    DO WHILE (ICHAR(tmp(st:st))==0)
      tmp(st:st) = ' '
      st = st - 1
    ENDDO
    NCFileGetAttr = TRIM(tmp)
    
  END FUNCTION NCFileGetAttr

#ifdef HAVE_F2003

! Register NetCDF files as a usable file type
  SUBROUTINE NCRegister()
  TYPE (input_file_type_fcns_list), POINTER :: new

    ALLOCATE(new)
    new%next                => AllInputFileTypes
    new%FileType            =  "NetCDF"
    new%Suffix              =  "nc"
    new%FcnOpen             => NCFileOpen
    new%FcnClose            => NCFileClose
    new%FcnGetNRec          => NCFileGetNRec
    new%FcnRead             => NCFileRead
    new%FcnGetDimensions    => NCFileGetDimensions
    new%FcnGetVariables     => NCFileGetVariables
    new%FcnGetRecTime       => NCFileGetRecTime
    new%FcnGetAttr          => NCFileGetAttr
    AllInputFileTypes       => new

  END SUBROUTINE NCRegister
#endif

END MODULE mo_input_netcdf
