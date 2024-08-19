!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!

! This file contains constants and type declaration for mo_input.f90
! This file should never be used directly from outside the mo_input framework - use via mo_input.f90

MODULE mo_input_types
#ifdef INPUT_IN_ECHAM
  USE mo_kind,          ONLY: dp, i8
  USE mo_exception,     ONLY: message_text
  USE mo_mpi
  USE mo_util_string,   ONLY: tolower
IMPLICIT NONE
#elif defined INPUT_IN_ICON
! To be filled in!
#elif defined USE_MPI
  USE mo_kind
  USE mo_mpi
  USE mo_input_strings, ONLY: message_text, tolower
IMPLICIT NONE
#else
  USE mo_input_strings, ONLY: message_text, tolower
IMPLICIT NONE
! From mo_kind.f90
  INTEGER, PARAMETER :: ps = 6
  INTEGER, PARAMETER :: rs = 37

  INTEGER, PARAMETER :: pd = 12
  INTEGER, PARAMETER :: rd = 307

  INTEGER, PARAMETER :: si =  9
  INTEGER, PARAMETER :: ri = 14

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(ps,rs)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd)
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(si)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(ri)
  PUBLIC dp, sp, i4, i8
#endif

PUBLIC

PUBLIC :: message_text, dp, i8, tolower

! General constants
INTEGER, PARAMETER :: INPUT_ALLFLAGS         =  2147483647      ! All possible flags: 0x7FFFFFFF = 2^31-1 
INTEGER, PARAMETER :: UNLIKELY_VAL           =  -452863529      ! Arbitrary value unlikely to be purposedly used
REAL(dp),PARAMETER :: UNLIKELY_RVAL          =  REAL(UNLIKELY_VAL,dp)/1000._dp  ! Floating point version of above
INTEGER, PARAMETER :: MAX_LINE_LEN           =         150      ! Write to new line if line is already as long as this
CHARACTER (len=*), PARAMETER :: opt_logical  = ".false. .true." ! Possible strings for logical variables

! Specialities which can only be applied if mo_input is compiled with a Fortran 2003 compiler
#ifdef HAVE_F2003
! List of procedure pointers for data processing
TYPE input_fcn_list
  TYPE (input_fcn_list),        POINTER :: next         ! Pointer to next procedure pointer
  PROCEDURE (DataProcTemplate), POINTER, NOPASS :: fcn
END TYPE input_fcn_list

! Type holding the functions to handle different file types generically
TYPE input_file_type_fcns_list
  TYPE (input_file_type_fcns_list), POINTER         :: Next                ! Next file type
  PROCEDURE (FileOpen),             POINTER, NOPASS :: FcnOpen             ! Function to open a file of the specified type
  PROCEDURE (FileClose),            POINTER, NOPASS :: FcnClose            ! Function to close a file of the specified type
  PROCEDURE (FileGetNRec),          POINTER, NOPASS :: FcnGetNRec          ! Function to get the number of records in a file
  PROCEDURE (FileRead),             POINTER, NOPASS :: FcnRead             ! Function to read a data chung from a file
  PROCEDURE (FileGetDimensions),    POINTER, NOPASS :: FcnGetDimensions    ! Fcn to get the list of dimensions present in file
  PROCEDURE (FileGetVariables),     POINTER, NOPASS :: FcnGetVariables     ! Fcn to get the list of variables present in file 
  PROCEDURE (FileGetRecTime),       POINTER, NOPASS :: FcnGetRecTime       ! Function to get the real time of a record in file
  PROCEDURE (FileGetAttr),          POINTER, NOPASS :: FcnGetAttr          ! Function to get a text attribute from var or file
  CHARACTER (len=10)                                :: FileType            ! Name of type of file handled by these functions
  CHARACTER (len=4)                                 :: Suffix              ! Files with this extention should be handled
END TYPE input_file_type_fcns_list
#endif

! Constants for identifying at which stage the user is asked to eventually insert additional actions
INTEGER, PARAMETER :: INPUT_ACTION_OBTAIN            =  1
INTEGER, PARAMETER :: INPUT_ACTION_AFTER_READ        =  2
INTEGER, PARAMETER :: INPUT_ACTION_AFTER_CONDENSE    =  3
INTEGER, PARAMETER :: INPUT_ACTION_AFTER_BROADCAST   =  4
INTEGER, PARAMETER :: INPUT_ACTION_AFTER_RESCALE     =  5
INTEGER, PARAMETER :: INPUT_ACTION_AFTER_REVERSE     =  6
INTEGER, PARAMETER :: INPUT_ACTION_AFTER_PACK        =  7
INTEGER, PARAMETER :: INPUT_ACTION_AFTER_DIM_CORRECT =  8
INTEGER, PARAMETER :: INPUT_ACTION_AFTER_INTERPOLATE =  9

! Constants for level parameter in the DataProcXXX subroutines
INTEGER, PARAMETER :: INPUT_PROC_DATA_ACTION       =  1 ! 0x00000001 : Perform the actual data processing
INTEGER, PARAMETER :: INPUT_PROC_PROGRESS          =  2 ! 0x00000002 : Progress acp
INTEGER, PARAMETER :: INPUT_PROC_DO_DATA_ACTION    =  3 ! 0x00000003 : Std. processing (ACTION+PROGRESS)
INTEGER, PARAMETER :: INPUT_PROC_NAME              =  4 ! 0x00000004 : Return name without processing data (msg must be present)
INTEGER, PARAMETER :: INPUT_PROC_MOD_DIM_FORWARD   =  8 ! 0x00000008 : Update dimensions according to the forward process(not yet)
INTEGER, PARAMETER :: INPUT_PROC_MOD_DIM_BACKWARD  = 16 ! 0x00000010 : Update dims according to the backward process (not yet)

! Used to hold a single data set
TYPE input_data_list
  TYPE (input_data_list), POINTER :: next               ! Pointer to next interpolation data
  TYPE (input_data_list), POINTER :: prev               ! Pointer to previous interpolation data
  REAL(dp), POINTER               :: dta(:,:,:,:,:,:,:) ! Pointer to the data themselves
  INTEGER                         :: flags              ! Flags for buffers (bitwise) INPUT_DATA_xxx constants below
  INTEGER                         :: time(4)            ! Time for which the data are valid
END TYPE input_data_list

! Constants for bitwise input_data_list%flags
INTEGER, PARAMETER :: INPUT_DATA_NREF_MASK   =      16383 ! 0x00003FFF : Number of different objects using this buffer
INTEGER, PARAMETER :: INPUT_DATA_NO_STEP     =      16384 ! 0x00004000 : Step for identification number
INTEGER, PARAMETER :: INPUT_DATA_NO_MASK     =     507904 ! 0x0007C000 : Identification number for this data type
INTEGER, PARAMETER :: INPUT_DATA_ALLOC       =     524288 ! 0x00080000 : Data of buffer has been allocated from the buffer system
INTEGER, PARAMETER :: INPUT_DATA_GLOBAL      =    1048576 ! 0x00100000 : This is a global buffer
INTEGER, PARAMETER :: INPUT_DATA_IN_USE      =    2097152 ! 0x00200000 : Data area currently taking part in an asyncroneous event
INTEGER, PARAMETER :: INPUT_DATA_DEFAULT     =    4194304 ! 0x00400000 : This datum is the default value to use
INTEGER, PARAMETER :: INPUT_DATA_FILL        =    8388608 ! 0x00800000 : Data (model,file) for _FillValue/missing_value
INTEGER, PARAMETER :: INPUT_DATA_RANGE       =   16777216 ! 0x01000000 : Data (model,file_low,file_high) for valid_range
INTEGER, PARAMETER :: INPUT_DATA_RESCALE     =   33554432 ! 0x02000000 : These data are used for scaling the data
INTEGER, PARAMETER :: INPUT_DATA_NUDGE_SCA   =   67108864 ! 0x04000000 : Datum is used for nudging with a scalar nudging factor
INTEGER, PARAMETER :: INPUT_DATA_NUDGE_FLD   =  134217728 ! 0x08000000 : Data is a field with nudging factors
INTEGER, PARAMETER :: INPUT_DATA_WEIGHTS     =  268435456 ! 0x10000000 : Data are used for weigthed entries along a dimension
INTEGER, PARAMETER :: INPUT_DATA_FUTURE_FRAC =  536870912 ! 0x20000000 : (Internal) future fraction for first interpolation
INTEGER, PARAMETER :: INPUT_DATA_USER        = 1073741824 ! 0x40000000 : This buffer holds some unspecified kind of user data
INTEGER, PARAMETER :: INPUT_DATA_TYPE_MASK   = 2143289344 ! 0x7FC00000 : Mask for type of data
INTEGER, PARAMETER :: flg_data_n             =         12 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_data     = &
  "alloc global in_use default fill_bad fill_range rescale nudge_scalar nudge_field weights future_frac user"

! Linked list to hold integer masks
TYPE input_mask_list
  TYPE (input_mask_list),    POINTER :: next            ! Next mask
  INTEGER,                   POINTER :: mask(:)         ! Mask data
  INTEGER                            :: npts            ! Number of points in total field for mask
  INTEGER                            :: nvalid          ! Number of "good" points in mask
  INTEGER                            :: flags           ! Bitwise flags for the mask INPUT_MASK_xxx
END TYPE input_mask_list

! Constants for bitwise flags for input_mask_list%flags
INTEGER, PARAMETER :: INPUT_MASK_FIRST_VALID =   16777215 ! 0x00FFFFFF : Start of valid data if continuous 
INTEGER, PARAMETER :: INPUT_MASK_TYPE_MODEL  =          0 ! 0x00000000 : Mask is a model mask describing points to use
INTEGER, PARAMETER :: INPUT_MASK_TYPE_DISTRIB=   16777216 ! 0x01000000 : Mask is a mask for PE-distribution of data
INTEGER, PARAMETER :: INPUT_MASK_TYPE_PE     =   33554432 ! 0x02000000 : Mask describes how dim(s) are distributed on PEs
INTEGER, PARAMETER :: INPUT_MASK_TYPE_EQUIV  =   50331648 ! 0x03000000 : Mask is used to map to model equivalent dimension(s)
INTEGER, PARAMETER :: INPUT_MASK_TYPE_MASK   =   50331648 ! 0x03000000 : Mask for mask types
INTEGER, PARAMETER :: INPUT_MASK_GLOBAL      =          0 ! 0x00000000 : Mask should be applied on global data
INTEGER, PARAMETER :: INPUT_MASK_CHUNK       =  134217728 ! 0x08000000 : Mask should be applied on chunk data
INTEGER, PARAMETER :: INPUT_MASK_LOCAL       =  268435456 ! 0x10000000 : Mask should be applied on PE-local data
INTEGER, PARAMETER :: INPUT_MASK_DOMAIN_MASK =  402653184 ! 0x18000000 : Mask for domain flags above
INTEGER, PARAMETER :: INPUT_MASK_CONT        =  536870912 ! 0x20000000 : All valid data are continuous
INTEGER, PARAMETER :: INPUT_MASK_PROCESSED   = 1073741824 ! 0x40000000 : Mask of PE has been processed (MaskCreateScatter internal)
INTEGER, PARAMETER :: flg_mask_n             =          2 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_mask  = 'continuous processed global chunk local'
CHARACTER (len=*), PARAMETER :: mask_type = 'model distributing PE-desc equivalence'

! Type to ease the construction and processing of the processing chain
TYPE input_curr
  INTEGER                          :: acp               ! "Actual Current pPinter"
  INTEGER                          :: acp_pre           ! Previous "Actual Current Pointer"
  TYPE (input_data_list),  POINTER :: buf               ! Actual buffer
  TYPE (input_data_list),  POINTER :: buf_pre           ! Previous buffer
  TYPE (input_mask_list),  POINTER :: mask              ! Actual mask
END TYPE input_curr

! Linked list to hold all "raw" dimensions/variables in a file
TYPE input_fileobj_list
  TYPE (input_fileobj_list), POINTER :: next            ! Next object
  CHARACTER (LEN=64)                 :: name_obj        ! Name of object
  INTEGER                            :: id              ! File Id of the object if appliable (eg. NetCDF files)
  INTEGER                            :: size_global     ! Size (#points) (dimensions, -ve=unlimited), ndims (vars)
  INTEGER                            :: flags           ! Bitwise flags for the object INPUT_FILEDIM_xxx/INPUT_FILEVAR_xxx
  REAL(dp), DIMENSION(:),    POINTER :: points          ! All data points of dimension/attribute values for variables
  INTEGER                            :: fids_dim(7)     ! Ids of the dimensions of a variable
END TYPE input_fileobj_list

! Constants for bitwise flags for input_fileobj_list%flags
INTEGER, PARAMETER :: INPUT_FILEOBJ_UNDET     =       0 ! 0x000000 : Object is of undefined type
INTEGER, PARAMETER :: INPUT_FILEOBJ_DIM       =       1 ! 0x000001 : Object is a file dimension description
INTEGER, PARAMETER :: INPUT_FILEOBJ_VAR       =       2 ! 0x000002 : Object is a file variable description
INTEGER, PARAMETER :: INPUT_FILEDIM_DIR_UNDET =       0 ! 0x000000 : No attempt to determine direction
INTEGER, PARAMETER :: INPUT_FILEDIM_DIR_ERROR =       4 ! 0x000004 : Direction indeterminable
INTEGER, PARAMETER :: INPUT_FILEDIM_DIR_UNLIM =       8 ! 0x000008 : Dimension is the unlimited one
INTEGER, PARAMETER :: INPUT_FILEDIM_DIR_DEC   =      16 ! 0x000010 : Values of the dimension are in decreasing order
INTEGER, PARAMETER :: INPUT_FILEDIM_DIR_INC   =      32 ! 0x000020 : Values of the dimension are in increasing order
INTEGER, PARAMETER :: INPUT_FILEDIM_DIR_MASK  =      60 ! 0x00003C : Mask for flags above
CHARACTER (len=*), PARAMETER :: flg_fileobj_dir = 'undetermined,undeterminable,unlimited,decreasing,increasing'
INTEGER, PARAMETER :: INPUT_FILEVAR_HAS_SCALE =      64 ! 0x000040 : Scaling attributes found
INTEGER, PARAMETER :: INPUT_FILEVAR_HAS_BAD   =     128 ! 0x000080 : Bad value specification found
INTEGER, PARAMETER :: INPUT_FILEVAR_HAS_RANGE =     256 ! 0x000100 : Specification of valid range found
INTEGER, PARAMETER :: INPUT_FILEVAR_MASK      =     448 ! 0x0001C0 : Mask for flags above
INTEGER, PARAMETER :: flg_filevar_n = 3                 ! #Flags above
CHARACTER (len=*), PARAMETER :: flg_fileobj_var = 'scale,missing_value,valid_range'

! Linked list to hold all dimensions needed by the model, including information only relevant to packed dimensions
! Splitted in two parts as workaround to Fortran not knowing about arrays of pointers to structures
TYPE input_dim_data
  TYPE (input_mask_list), POINTER :: mask               ! Mask to pack data for unpacked file and packed memory data
  TYPE (input_mask_list), POINTER :: local              ! Mask for points local to this PE (IODATA-PEs only)
  CHARACTER (LEN=64)              :: name_dim           ! Name of dimension
  CHARACTER (LEN=64)              :: name_alt           ! Alternative name of dimension
  REAL(dp)                        :: origin             ! Expected value of first entry of dimension values
  INTEGER                         :: sub_model          ! Sub-model to which this dimension belongs
  INTEGER                         :: size_global        ! Expected/actual global size of dimension
  INTEGER                         :: size_local         ! Local size of dimension
  INTEGER                         :: local_lo           ! Lower boundary of local data
  INTEGER                         :: local_hi           ! Upper boundary of local data
  INTEGER                         :: chunk_lo           ! Lower boundary of local read chunck (used on IODATA-PEs only)
  INTEGER                         :: chunk_hi           ! Upper boundary of local read chunck (used on IODATA-PEs only)
  INTEGER                         :: flags              ! Bitwise, see INPUT_DIM_xxx constants below
  INTEGER                         :: nDst               ! Number of destination ranks for this dimension
  INTEGER                         :: nsrc               ! Total number of points in unpacked dimensions of this packed dimension
END TYPE input_dim_data

TYPE input_dim_list
  TYPE (input_dim_list),  POINTER :: next               ! Next dimension
  TYPE (input_dim_list),  POINTER :: src_dims           ! List of dimensions which are packed into this dimension
  TYPE (input_dim_data),  POINTER :: dim_data           ! The information on the dimension itself
  CHARACTER (len=32)              :: file_dim_name      ! Renamed name of dimension, used within input_file_list only
  INTEGER                         :: lo                 ! Lower boundary of chunk/local when used as source for packed/equiv. dim.
  INTEGER                         :: hi                 ! Upper boundary of chunk/local when used as source for packed/equiv. dim.
  INTEGER                         :: flags              ! Bitwise. Used within input_file_list and input_var_list only
  INTEGER                         :: fid                ! Id in file
END TYPE input_dim_list

! Type for declaring one (list of) dimensions equivalent to another (list of) dimensions
TYPE input_eqdim_list
  TYPE (input_eqdim_list), POINTER :: next              ! Next set of dimensions which are equivalent
  TYPE (input_dim_list),   POINTER :: src_dims          ! Source      dimensions of equivalence
  TYPE (input_dim_list),   POINTER :: dst_dims          ! Destination dimensions of equivalence
  TYPE (input_mask_list),  POINTER :: mask              ! Mask to do the mapping
END TYPE input_eqdim_list

! Constants for bitwise input_dim_elm%flags
INTEGER, PARAMETER :: INPUT_DIM_UNLIM      =       1    ! 0x000001 : Split this dimension (user set)
INTEGER, PARAMETER :: INPUT_DIM_PACKED     =       2    ! 0x000002 : Packed dimension (indirect user)
INTEGER, PARAMETER :: INPUT_DIM_EXIST      =       4    ! 0x000004 : Dimension found in file (short term, system)
INTEGER, PARAMETER :: INPUT_DIM_REVERSE    =       8    ! 0x000008 : Dimension found, but in reverse order (system, var_list only)
INTEGER, PARAMETER :: INPUT_DIM_CYCLIC     =      16    ! 0x000010 : Dimension is cyclic and can be shifted (user)
INTEGER, PARAMETER :: INPUT_DIM_COLLAPSED  =      32    ! 0x000020 : Dimension not found - spread needed (system, var_list only)
INTEGER, PARAMETER :: INPUT_DIM_PARALLEL   =      64    ! 0x000040 : Dimension split parallel (system)
INTEGER, PARAMETER :: INPUT_DIM_DISTRIB    =     128    ! 0x000080 : Dimension parallelization requires scattering of data (sy)
INTEGER, PARAMETER :: INPUT_DIM_INC        =       0    ! 0x000000 : Dimension is expected to contain increasing values (user)
INTEGER, PARAMETER :: INPUT_DIM_DEC        =     INPUT_DIM_REVERSE !: Dimension is expected to contain decreasing values (user)
INTEGER, PARAMETER :: INPUT_DIM_NON_MODEL  =     256    ! 0x000100 : Dimension exist in file, but is unknown to the model
INTEGER, PARAMETER :: INPUT_DIM_READ       =     512    ! 0x000200 : Dimension were among those which were directly read (var only)
INTEGER, PARAMETER :: INPUT_DIM_LOCAL_DATA =    1024    ! 0x000400 : dim_data field is not associated with global dim list
INTEGER, PARAMETER :: INPUT_DIM_EXTRA      =    2048    ! 0x000800 : Extra dimension internally created for multi-variable
INTEGER, PARAMETER :: INPUT_DIM_FILE       =    4096    ! 0x001000 : Dimension header belongs to a file object
INTEGER, PARAMETER :: INPUT_DIM_VAR        =    8192    ! 0x002000 : Dimension header belongs to a variable object
INTEGER, PARAMETER :: INPUT_DIM_GLOBAL     =   16384    ! 0x004000 : Dimension header belongs to the global list
INTEGER, PARAMETER :: INPUT_DIM_MASK       =   32768    ! 0x008000 : Dimension header belongs to a mask object
INTEGER, PARAMETER :: INPUT_DIM_SOURCE     =   65536    ! 0x010000 : Dimension header is the source of a packed dimension
INTEGER, PARAMETER :: INPUT_DIM_CHUNKED    =  131072    ! 0x020000 : Dimension header is chunked and packed, try reduce src. chunks
INTEGER, PARAMETER :: INPUT_DIM_EQ_SRC     =  262144    ! 0x040000 : Dim. h. is part of the list of src dims of an equvalence decl.
INTEGER, PARAMETER :: INPUT_DIM_EQ_DST     =  524288    ! 0x080000 : Dim. h. is part of the list of dst dims of an equvalence decl.
INTEGER, PARAMETER :: INPUT_DIM_USEASGLOBAL= 1048576    ! 0x100000 : Use as global dimension despite of chunk definition
INTEGER, PARAMETER :: INPUT_DIM_TYPE       = &
  INPUT_DIM_VAR+INPUT_DIM_GLOBAL+INPUT_DIM_MASK+INPUT_DIM_SOURCE+INPUT_DIM_EXTRA+INPUT_DIM_EQ_SRC+INPUT_DIM_EQ_DST
INTEGER, PARAMETER :: flg_dim_n            =      21    ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_dim = "unlim packed exist reverse cyclic collapsed parallel distribute non_model " // &
"read local extra file variable global mask source chunked equiv_source equiv_dest use_as_global"

! Constants returned by InputDimGetRefEx
INTEGER, PARAMETER :: INPUT_DIM_FND_NOT      =     0    ! 0x00000 : Dimension not found
INTEGER, PARAMETER :: INPUT_DIM_FND_NORMAL   =     1    ! 0x00001 : Regular dimension found
INTEGER, PARAMETER :: INPUT_DIM_FND_SOURCE   =     2    ! 0x00002 : Found as a source dimension of a packed dimension
INTEGER, PARAMETER :: INPUT_DIM_FND_PACKED   =     3    ! 0x00003 : Found as a packed dimension to be unpacked
INTEGER, PARAMETER :: INPUT_DIM_FND_MASK     =    67    ! 0x00043 : Mask for dimension found
INTEGER, PARAMETER :: INPUT_DIM_FND_REVERSE  =     4    ! 0x00004 : Dimension is in opposite order than expected
INTEGER, PARAMETER :: INPUT_DIM_FND_SWAP     =     8    ! 0x00008 : Cyclic dimension origin does not match expectations
INTEGER, PARAMETER :: INPUT_DIM_FND_SCATTER  =    16    ! 0x00010 : Dimension is to be scattered
INTEGER, PARAMETER :: INPUT_DIM_FND_DISTRIB  =    32    ! 0x00020 : Dimension "participates" in a dimension to be scattered
INTEGER, PARAMETER :: INPUT_DIM_FND_EQ_SRC   =    64    ! 0x00040 : Dimension is a source dimension of a dimension equivalence
INTEGER, PARAMETER :: INPUT_DIM_FND_EQ_DST   =   128    ! 0x00080 : Dimension is a destination dimension of a dimension equivalence

! Linked list of files to be read from mo_input
TYPE input_file_list
  TYPE (input_file_list),    POINTER :: next           ! Next input file
  TYPE (input_fileobj_list), POINTER :: fdims          ! Pointer to all dimensions found in current file
  TYPE (input_dim_list),     POINTER :: dims           ! Pointer to model dimensions found in current file
  TYPE (input_fileobj_list), POINTER :: fvars          ! Pointer to all variables found in current file
#ifdef HAVE_F2003
  TYPE (input_file_type_fcns_list), POINTER :: ft      ! Functions to operate on the specific file type
#endif
  CHARACTER (len=64)              :: name_pre          ! Prefix of filename (before time specification)
  CHARACTER (len=64)              :: name_suf          ! Suffix of filename (after time specification)
  CHARACTER (len=128)             :: name_file         ! Exact complete name of currently opened file
  CHARACTER (len=64)              :: time_unit         ! String with the unit of the time attribute (NetCDF specific)
  INTEGER                         :: time_fid(2)       ! File obj identifieres for the splitted dimension as dim and var resp.
  INTEGER                         :: time_scale(5)     ! Scaling of file to model time for linear file times
                                                       ! mod_t = CalendarTimeAdd(time_scale(1:4),time_scale(5)*read_val)
  INTEGER                         :: sub_model         ! Sub-model to which this file belongs
  INTEGER                         :: fid               ! NetCDF-Id of the file (could also be used as unit number)
  INTEGER                         :: nUsed             ! #Times the file has been opened
  INTEGER                         :: flags             ! Bitwise, see INPUT_FILE_xxx constants below. File properties
  INTEGER                         :: stat              ! Bitwise, see INPUT_FILE_xxx constants below. File status
  INTEGER                         :: nVar              ! Number of variables associated with file
  INTEGER                         :: nRead             ! Number of variables read from current record
  INTEGER                         :: nRead_dbg         ! Number of variables read in current model time step
  INTEGER                         :: nRec              ! Number of records (time steps) in file not yet read
  INTEGER                         :: curr_rec          ! Number of next record to read
  INTEGER                         :: init_rec          ! Number of record to read for initial variables with time dimension
  INTEGER                         :: rec_step          ! Number of file records between reads
  INTEGER                         :: time_offset(4)    ! Time offset (years,days,seconds,records) in file for var. validity time
  INTEGER                         :: dt_file           ! Time step in file (+ve: seconds) or (-ve: months)
  INTEGER                         :: time_open(4)      ! Time of currently opened file
  INTEGER                         :: time_file(4)      ! Time of next file to open file
  INTEGER                         :: time_start(5)     ! Start file time + start record
  INTEGER                         :: time_cycle(7)     ! File time to cycle to @time_enddata  +rec number+#rec in cycle offset
  INTEGER                         :: time_enddata(6)   ! First time + rec of unavailable data + #rec available
END TYPE input_file_list

! Flags for input_file_list%flags
!   Time validity
INTEGER, PARAMETER :: INPUT_FILE_INITIAL         =          1 ! 0x00000001 : This file is an initial file (user)
INTEGER, PARAMETER :: INPUT_FILE_TIME_ABSOLUTE   =          2 ! 0x00000002 : Use absolute time if possible (user set)
INTEGER, PARAMETER :: INPUT_FILE_TIME_FROM_NAME  =          4 ! 0x00000004 : Assume that the time file name is that of the 1st rec
INTEGER, PARAMETER :: INPUT_FILE_NOCYCLE         =          8 ! 0x00000008 : Keep end (start) data instead of cycling file (user)
INTEGER, PARAMETER :: flg_file_time_n            =          4 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_file_time = "initial absolute name_time no_cycle"
!   Name construction
INTEGER, PARAMETER :: INPUT_FILE_NAME_YEAR_SIGN  =         16 ! 0x00000010 : Construct name with signed + 6 digit year (ind. user)
INTEGER, PARAMETER :: INPUT_FILE_NAME_YEAR_LONG  =         32 ! 0x00000020 : Construct name with 6 digit year (indirectly user)
INTEGER, PARAMETER :: INPUT_FILE_NAME_YEAR       =         48 ! 0x00000030 : Construct name with 4 digit year (indirectly user)
INTEGER, PARAMETER :: INPUT_FILE_NAME_MONTH      =         64 ! 0x00000040 : Construct name with 2 digit month of year (ind. user)
INTEGER, PARAMETER :: INPUT_FILE_NAME_DAY        =        128 ! 0x00000080 : Construct name with 2 digit day of month (ind. user)
INTEGER, PARAMETER :: INPUT_FILE_NAME_HOUR       =        256 ! 0x00000100 : Construct name with 2 digit hour of day (ind. user)
INTEGER, PARAMETER :: INPUT_FILE_NAME_MINUTE     =        512 ! 0x00000200 : Construct name with 2 digit min. of hour (ind. user)
INTEGER, PARAMETER :: INPUT_FILE_NAME_SECOND     =       1024 ! 0x00000400 : Construct name with 2 digit sec. of min.(ind. user)
INTEGER, PARAMETER :: INPUT_FILE_MULTI           =       2032 ! 0x000007F0 : File cont.time dependent parts (sum of 6 flags above)
!   Time error handling
INTEGER, PARAMETER :: INPUT_FILE_TIME_ERR        =          0 ! 0x00000000 : Time check: Mismatch is an error (user set, default).
INTEGER, PARAMETER :: INPUT_FILE_TIME_WARN       =       4096 ! 0x00001000 : Time check: Warn if mismatch (user set).
INTEGER, PARAMETER :: INPUT_FILE_TIME_IGNORE     =       8096 ! 0x00002000 : Ignore if absolute times mismatch (user set).
INTEGER, PARAMETER :: INPUT_FILE_TIME_MASK       =      12288 ! 0x00003000 : Time check mask (combines 3 flags above)
CHARACTER (len=*), PARAMETER :: opt_time_err = "err warn ignore"     ! Possible strings for action_time_err
! Free bits 0x7FFFC800

! Status flags stored input_file_list%stat
!   Actual open/initialization status
INTEGER, PARAMETER :: INPUT_FILE_OPEN            =          1 ! 0x00000001 : File currently open (for reading) (system)
INTEGER, PARAMETER :: INPUT_FILE_MULTI_OPENED    =          2 ! 0x00000002 : At least one file of a multi-file has been found(sys)
INTEGER, PARAMETER :: INPUT_FILE_DIMS_GOT        =          4 ! 0x00000004 : File opened and file dimensions obtained(sys)
INTEGER, PARAMETER :: INPUT_FILE_VARS_GOT        =          8 ! 0x00000008 : File opened and file variables  obtained(sys)
INTEGER, PARAMETER :: INPUT_FILE_TIME_GOT        =         16 ! 0x00000010 : File opened and file time info  obtained(sys)
INTEGER, PARAMETER :: INPUT_FILE_FIRST_CYCLE     =         32 ! 0x00000020 : Data before start time - cycle to start (system)
INTEGER, PARAMETER :: INPUT_FILE_END_REACHED     =         64 ! 0x00000040 : No more data should be read (system)
INTEGER, PARAMETER :: INPUT_FILE_CYCLED          =        128 ! 0x00000080 : Has Returned to cycle time (system)
INTEGER, PARAMETER :: INPUT_FILE_READ            =        256 ! 0x00000100 : Read by last update call (system)
INTEGER, PARAMETER :: INPUT_FILE_ENDDATE_GOT     =        512 ! 0x00000200 : First date without data obtained (system)
INTEGER, PARAMETER :: INPUT_FILE_EXTERN_APPLIED  =       1024 ! 0x00000400 : Eventual external parameters has been applied (sys)
INTEGER, PARAMETER :: flg_file_stat_n            =         11 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_file_stat    = &
  "open multi_present has_dims has_vars has_time first_cycle end_reached cycled read_now enddate_correct settings_applied"
!   Flags for (to be) performed checks
INTEGER, PARAMETER :: INPUT_FILE_CHECK_DIM_SIZE  =       2048 ! 0x00000800 : File dimensions
INTEGER, PARAMETER :: INPUT_FILE_CHECK_NRECS     =       4096 ! 0x00001000 : Number of records in file
INTEGER, PARAMETER :: INPUT_FILE_CHECK_TIME_STEP =       8192 ! 0x00002000 : Time step in file
INTEGER, PARAMETER :: INPUT_FILE_CHECK_TIME_CYCLE=      16384 ! 0x00004000 : Cyclicity time
INTEGER, PARAMETER :: INPUT_FILE_CHECK_VAR       =      32768 ! 0x00008000 : Variables
INTEGER, PARAMETER :: INPUT_FILE_CHECK_VAR_DIMS  =      65536 ! 0x00010000 : Dimensionality of file variables
INTEGER, PARAMETER :: INPUT_FILE_CHECK_INITIAL   = INPUT_FILE_CHECK_DIM_SIZE + INPUT_FILE_CHECK_VAR + INPUT_FILE_CHECK_VAR_DIMS
INTEGER, PARAMETER :: INPUT_FILE_CHECK_MASK      =     129024 ! 0x0001F800 : Mask for check flags
INTEGER, PARAMETER :: flg_file_check_n           =          6 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_file_check = "dim_size nrec time_step time_cycle var var_dims"
! Free bits 0x7FFFE0000

! For passing to InputFileAdd%cycle_act
INTEGER, PARAMETER :: INPUT_FILE_REPEAT_NONE     = HUGE(i8)-1 ! Keep at end value
INTEGER, PARAMETER :: INPUT_FILE_REPEAT_ALL      = HUGE(i8)-2 ! Repeat entire forcing period
INTEGER, PARAMETER :: INPUT_FILE_REPEAT_YEAR     = HUGE(i8)-3 ! Repeat last year of forcing period
INTEGER, PARAMETER :: INPUT_FILE_REPEAT_DAY      = HUGE(i8)-4 ! Repeat last day of forcing period
INTEGER, PARAMETER :: INPUT_FILE_REPEAT_RECS     = HUGE(i8)-5 ! Repeat last n records of forcing period
INTEGER, PARAMETER :: INPUT_FILE_REPEAT_SECS     = HUGE(i8)-6 ! Repeat last n seconds of forcing period
INTEGER, PARAMETER :: INPUT_FILE_REPEAT_MONTHS   = HUGE(i8)-7 ! Repeat last n months of forcing period
CHARACTER (len=*), PARAMETER :: opt_cycle    = "none all years days records seconds months" ! Possible strings for action_cycle

! Linked list of variable groups
TYPE input_group_list
  TYPE (input_group_list), POINTER :: next
  TYPE (input_file_list),  POINTER :: group_file        ! The group file (or file of first variable if multiple)
  CHARACTER (len=64)               :: name_group        ! Name of the group
  INTEGER                          :: flags             ! Bitwise, INPUT_GROUP_xxx below and INPUT_VAR_VALID_xxx from ..._var_list
  INTEGER                          :: dt_group          ! Update time step of the group
  INTEGER                          :: nIntrp, nAtOnce, nFuture ! Interpolation settings of group
  INTEGER                          :: time(4)
  INTEGER                          :: curr_rec
END TYPE input_group_list

! Flags for input_group_list%flags
INTEGER, PARAMETER :: INPUT_GROUP_MULTI_FILE     =          1 ! 0x00000001 : Group variables are in multiple files
INTEGER, PARAMETER :: INPUT_GROUP_MSG_DONE       =          2 ! 0x00000002 : Message from group has been displayed
INTEGER, PARAMETER :: flg_group_n                =          2 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: opt_group        = "multi_file msg_done" ! Possible strings for action_cycle

! Linked list of variables managed by mo_input
TYPE input_var_list
  TYPE (input_var_list),  POINTER :: next               ! Next input variable
  TYPE (input_var_list),  POINTER :: depend             ! Variable is dependent on this variable
  TYPE (input_file_list), POINTER :: parent_file        ! (Multi-)file to read variable from
  TYPE (input_dim_list),  POINTER :: dims               ! Dimensions of variable in model
  TYPE (input_dim_list),  POINTER :: dims_at_read       ! Dimensions of variable in file 
  TYPE (input_data_list), POINTER :: intrp              ! Fields from which to interpolate (newest first)
  TYPE (input_data_list), POINTER :: buffers            ! Temporary buffers for reshape/permute/distribute etc.
  TYPE (input_data_list), POINTER :: buffers_intrp      ! Shortcut to the buffers to use from time interpolation onward
  TYPE (input_mask_list), POINTER :: masks              ! Masks to (un)pack and distribute data 
  TYPE (input_mask_list), POINTER :: masks_intrp        ! Shortcut to masks to use from time interpolation onward
  TYPE (input_data_list), POINTER :: aux_data           ! Any extra data needed for the processing
  TYPE (input_group_list),POINTER :: group              ! Group to which this variable belongs
#ifdef HAVE_F2003
  TYPE (input_fcn_list),  POINTER :: DataFcns           ! List of procedures to be executed when reading the data
  TYPE (input_fcn_list),  POINTER :: ProcInterpol       ! List of procedures to be executed when interpolating the data
  PROCEDURE (DataProcTemplate), POINTER, NOPASS :: TimeIntrpFcn ! User function to interpolate in time
  PROCEDURE (DataActionAddNone),POINTER, NOPASS :: AddActionFcn ! User function to add additional processing actions
#endif
  INTEGER,                POINTER :: saction(:)         ! List of actions to be performed to obtain data 
                                                        ! (series of INPUT_VAR_ACTION_xxx below + their data)
  REAL (dp),              POINTER :: dta(:,:,:,:,:,:,:) ! Pointer to actual data
  CHARACTER (LEN=64)              :: name_var           ! Name of variable
  CHARACTER (LEN=128)             :: name_alt           ! Alternative name(s) of variable
  INTEGER                         :: acp_intrp          ! Pointer into saction to use for interpolation
  INTEGER                         :: nIntrp             ! #Fields from which to interpolate
  INTEGER                         :: nAtOnce            ! #Time steps to read at once from file
  INTEGER                         :: nFuture            ! Minimum #time steps for interpolation in the future
  INTEGER                         :: sub_model          ! Sub-model to which this variable belongs
  INTEGER                         :: fid                ! File identification number of variable
  INTEGER                         :: flags              ! Bitwise, see INPUT_VAR_xxx constants below. Properties of variable
  INTEGER                         :: stat               ! Bitwise, see INPUT_VAR_xxx constants below. Status of variable
  INTEGER                         :: time(4)            ! Time of most recently read field
  INTEGER                         :: nToRead            ! #Records to be read before exiting "InputUpdate"
  INTEGER                         :: subset_index       ! Index into missing, non-expanding dimension
  INTEGER                         :: ndim               ! Number of dimensions of variable as seen in model
  INTEGER                         :: edims(7)           ! Expected dimensions sizes for local memory variable
  INTEGER                         :: nLastUpdate        ! Number of time-steps since last update
  INTEGER                         :: nLastRead          ! Number of time-steps since last read
  INTEGER                         :: nupdate(3)         ! Number of time steps between updates of the linear interpolation
  INTEGER                         :: dt_update          ! Time between updates (-ve: months, +ve: seconds)
  INTEGER                         :: next_update        ! Number of time steps until next update
  INTEGER                         :: nRead              ! Number of times new data were filled into the variable
END TYPE input_var_list

! Constants for bitwise input_var_list%flags. Bitwise relations between all INPUT_VAR_MISS_xx is important, don't change
!   Error handling if missing
INTEGER, PARAMETER :: INPUT_VAR_MISS_DEFAULT     =         0 ! 0x00000000 : Missing: To be specified dependent on def_val avail
INTEGER, PARAMETER :: INPUT_VAR_MISS_IGNORE      =         1 ! 0x00000001 : Missing: Ignore (neither warn nor error, user set)
INTEGER, PARAMETER :: INPUT_VAR_MISS_WARN        =         2 ! 0x00000002 : Missing: Warn  (user set)
INTEGER, PARAMETER :: INPUT_VAR_MISS_WARN_IGNORE =         3 ! 0x00000003 : Missing: Warn only when first = warn+ignore (user)
INTEGER, PARAMETER :: INPUT_VAR_MISS_ERR         =         4 ! 0x00000004 : Missing: Error (user set, default)
INTEGER, PARAMETER :: INPUT_VAR_MISS_ERR_IGNORE  =         5 ! 0x00000005 : Missing: Error only when first = err + ignore (user)
INTEGER, PARAMETER :: INPUT_VAR_MISS_ERR_WARN    =         6 ! 0x00000006 : Missing: Err when first, else warn =err+warn(user)
INTEGER, PARAMETER :: INPUT_VAR_MISS_MASK        =         7 ! 0x00000007 : Mask for error handling types above
CHARACTER (len=*), PARAMETER :: opt_var_miss = "ignore warn warn_ignore err err_ignore err_warn" ! Pos. strs for action_var_miss
!   Error handling if bad/out of range values are about to enter the model
INTEGER, PARAMETER :: INPUT_VAR_BAD_VAL_DEFAULT  =         0 ! 0x00000000 : Bad val: Ignore (neither warn nor error, user set)
INTEGER, PARAMETER :: INPUT_VAR_BAD_VAL_REPLACE  =  33554432 ! 0x02000000 : Bad val: Replace with fill-/def-value (def if given)
INTEGER, PARAMETER :: INPUT_VAR_BAD_VAL_WARN     =  67108864 ! 0x04000000 : Bad val: Warn about bad values (def if no fill/def)
INTEGER, PARAMETER :: INPUT_VAR_BAD_VAL_WARN_REPLACE=100663296!0x06000000 : Bad val: Warn about bad values and replace them
INTEGER, PARAMETER :: INPUT_VAR_BAD_VAL_ERR      = 134217728 ! 0x08000000 : Bad val: Bad values causes error termination
INTEGER, PARAMETER :: INPUT_VAR_BAD_VAL_IGNORE   = 167772160 ! 0x0A000000 : Bad val: Ignore (neither warn nor error, user set)
INTEGER, PARAMETER :: INPUT_VAR_BAD_VAL_MASK     = 234881024 ! 0x0E000000 : Mask for error handling types above
CHARACTER (len=*), PARAMETER :: opt_bad_val = "default replace warn warn_replace err ignore" ! Pos. strs for action_var_bad_val
!   Variable type and treatement                   
INTEGER, PARAMETER :: INPUT_VAR_INTERPOL_USER    =         8 ! 0x00000008 : Variable interp is to be done by the user(ind.user)
INTEGER, PARAMETER :: INPUT_VAR_GLOBAL           =        16 ! 0x00000010 : Local and global size are the same (user set)
INTEGER, PARAMETER :: INPUT_VAR_AUTORESET_OFF    =        32 ! 0x00000020 : Don't autom. reset UPDATED and READ flags (user)
INTEGER, PARAMETER :: INPUT_VAR_GROUP            =        64 ! 0x00000040 : Pseudo-var. - placeholder for group parameters
INTEGER, PARAMETER :: INPUT_VAR_MODEL            =       128 ! 0x00000080 : Variable normally updated by the model
INTEGER, PARAMETER :: flg_var_type_n             =         5 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_var_type = "user_interpol global autoreset_off group model"
!   Temporal validity constants (all user)
INTEGER, PARAMETER :: INPUT_VAR_VALID_ACTUAL     =         0 ! 0x00000000 : Variable is valid at indicated time
INTEGER, PARAMETER :: INPUT_VAR_VALID_BEFORE     =       256 ! 0x00000100 : Variable is valid in the middle of interval (before)
INTEGER, PARAMETER :: INPUT_VAR_VALID_START      =       512 ! 0x00000200 : Variable is valid at interval start
INTEGER, PARAMETER :: INPUT_VAR_VALID_MIDPOINT   =       768 ! 0x00000300 : Variable is valid in the middle of interval (after)
INTEGER, PARAMETER :: INPUT_VAR_VALID_END        =      1024 ! 0x00000400 : Variable is valid at interval end
INTEGER, PARAMETER :: INPUT_VAR_VALID_MIDDAY     =      1280 ! 0x00000500 : Variable is valid at the middel of the day
INTEGER, PARAMETER :: INPUT_VAR_VALID_MIDMONTH   =      1536 ! 0x00000600 : Variable is valid at the middel of the month
INTEGER, PARAMETER :: INPUT_VAR_VALID_MIDYEAR    =      1792 ! 0x00000700 : Variable is valid at the middel of the year
INTEGER, PARAMETER :: INPUT_VAR_VALID_MASK       =      1792 ! 0x00000700 : Mask for flags above
CHARACTER (len=*), PARAMETER :: opt_time_val = "actual before start midpoint end midday midmonth midyear" ! Record val. strings
!   Scale data with time options
INTEGER, PARAMETER :: INPUT_VAR_DIV_MONTH        =      2048 ! 0x00000800 : Divide data with length of month
INTEGER, PARAMETER :: INPUT_VAR_DIV_YEAR         =      4096 ! 0x00001000 : Divide data with length of year
INTEGER, PARAMETER :: INPUT_VAR_DIV_MASK         =      6144 ! 0x00001800 : Mask for above 
CHARACTER (len=*), PARAMETER :: opt_timediv = "div_month div_year"
!   Flags specifying how to combine file and model data
INTEGER, PARAMETER :: INPUT_VAR_ACTION_REPLACE   =         0 ! 0x00000000 : Replace model data with file data
INTEGER, PARAMETER :: INPUT_VAR_ACTION_ADD       =      8192 ! 0x00002000 : Add file to model data. Must be first non-z in flg grp
INTEGER, PARAMETER :: INPUT_VAR_ACTION_SUB       =     16384 ! 0x00004000 : Subtract file data from model data
INTEGER, PARAMETER :: INPUT_VAR_ACTION_SUBR      =     24576 ! 0x00006000 : Subtract model data from file data
INTEGER, PARAMETER :: INPUT_VAR_ACTION_MUL       =     32768 ! 0x00008000 : Multiply model and file data
INTEGER, PARAMETER :: INPUT_VAR_ACTION_DIV       =     40960 ! 0x0000A000 : Divide file data with model data
INTEGER, PARAMETER :: INPUT_VAR_ACTION_DIVR      =     49152 ! 0x0000C000 : Divide model data with file data
INTEGER, PARAMETER :: INPUT_VAR_ACTION_ADD_FLUX  =     57344 ! 0x0000E000 : Multiply file data with dt_model and add to model data
INTEGER, PARAMETER :: INPUT_VAR_ACTION_NUDGE     =     65536 ! 0x00010000 : File data * weight + model data * (1 - weight)
INTEGER, PARAMETER :: INPUT_VAR_ACTION_NUDGE_FLD =     73728 ! 0x00012000 : As nudge, but indv. weight for each pt. (var%userdata)
INTEGER, PARAMETER :: INPUT_VAR_ACTION_MASK      =    122880 ! 0x0001E000 : Mask for all INPUT_VAR_ACTION flags
CHARACTER (len=*), PARAMETER :: opt_data = "replace add sub subr mul div divr flux_add nudge nudge_fld"
!   Flags specifying actions on mismatching dimensions 
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_ERR      =         0 ! 0x00000000 : It is an error if not all dims are present (user)
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_WARN     =    262144 ! 0x00040000 : Warning if not all dim. are pres., combinable (user)
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_VAR_MISS =    524288 ! 0x00080000 : Treat as if variable was missing (user)
! 0x000C0000: Warn_var_miss
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_RESOLVE  =   1048576 ! 0x00100000 : User defined action will resolve dim mismatches (F2003)
! 0x00140000: Warn_resolve, 0x00180000, 0x001C0000: Crap
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_SPREAD   =   2097152 ! 0x00200000 : Spread data over missing dimensions (user)
! 0x00240000: Warn_spread, 0x00280000, 0x002C0000, 0x00300000, 0x00340000, 0x00380000, 0x003C0000: Crap
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_SUBSET   =   4194304 ! 0x00400000 : Set only a subset of the var. (max 1 miss dim.allowed)
! 0x00440000: Warn_subset,        0x00480000, 0x004C0000, 0x00500000, 0x00540000, 0x00580000 0x005C0000: Crap, 
! 0x00600000: Spread_subset,      0x00640000, Warn_spread_subset, 
! 0x00680000, 0x006C0000, 0x00700000, 0x00740000, 0x00780000, 0x007C0000: Crap
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_SUM      =   8388608 ! 0x00800000 : Sum over additional dimension(max 1 extra dim.allowed)
! 0x00840000: Warn_sum,           0x00880000, 0x008C0000, 0x00900000, 0x00940000, 0x00980000, 0x009C0000: Crap, 
! 0x00A00000: spread_sum          0x00A40000: Warn_spread_sum,    
! 0x00A80000, 0x00AC0000, 0x00B00000, 0x00B40000, 0x00B80000, 0x00BC0000: Crap, 
! 0x00C00000: subset_sum          0x00C40000: Warn_subset_sum,    
! 0x00C80000, 0x00CC0000, 0x00D00000, 0x00D40000, 0x00D80000, 0x00DC0000, 0x00E00000, 0x00E40000: Crap
! 0x00E80000, 0x00EC0000, 0x00F00000, 0x00F40000, 0x00F80000, 0x00FC0000: Crap
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_AVG      =  16777216 ! 0x01000000 : Average additional dimension (max 1 extra dim.allowed)
! 0x01040000: Warn_avg,           0x01080000, 0x010C0000, 0x01100000, 0x01140000 0x01180000, 0x011C0000, Crap, 
! 0x01200000: spread_avg          0x01240000: Warn_spread_avg, 
! 0x01280000, 0x012C0000,0x01300000, 0x01340000, 0x01380000, 0x013C0000: Crap, 
! 0x01400000: subset_avg          0x01440000: Warn_subset_avg,    
! 0x01480000, 0x014C0000, 0x01500000, 0x01540000, 0x01580000, 0x015C0000, 0x01600000, 0x01640000: Crap
! 0x01680000, 0x016C0000, 0x01700000, 0x01740000, 0x01780000, 0x017C0000: Crap
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_NORM     =  25165824 ! 0x01800000 : Normalize additional dimension (max 1 extra dim.)
! 0x01840000: Warn_norm,          0x01880000, 0x018C0000, 0x01900000, 0x01940000, 0x01980000, 0x019C0000: Crap, 
! 0x01A00000: spread_norm         0x01A40000: Warn_spread_norm,
! 0x01A80000, 0x01AC0000, 0x01B00000, 0x01B40000, 0x01B80000 0x01BC0000: Crap, 
! 0x01C00000: subset_norm         0x01C40000: Warn_subset_norm,
! 0x01D00000, 0x01D40000, 0x01D80000, 0x01DC0000, 0x01E00000, 0x01E40000, 0x01E80000, 0x01EC000: Crap
! 0x01F00000, 0x01F40000, 0x01F80000, 0x01FC0000: Crap
INTEGER, PARAMETER :: INPUT_VAR_DIM_MIS_MASK     =  33292288 ! 0x01FC0000 : Mask for all INPUT_VAR_DIM_MIS_MASK flags
CHARACTER (len=*), PARAMETER :: opt_dim_mis = &
"err,warn,var_miss,warn_var_miss,resolve,warn_resolve,,,spread,warn_spread,,,,,,,subset,warn_subset,,,,,,,"                     //&
"spread subset,w(spread subset),,,,,,,sum,warn_sum,,,,,,,spread sum,w(spread sum),,,,,,,subset sum,w(subset sum),,,,,,,,,,,,,,,"//&
"avg,warn_avg,,,,,,,spread avg,w(spread avg),,,,,,,subset avg,w(subset avg),,,,,,,,,,,,,,norm,warn_norm,,,,,,,"                 //&
"spread norm,w(spread norm),,,,,,,subset norm,w(subset norm)"
! Free bits 0x70020000

! Status flags (all system) stored in input_var_list%stat for normal variables
INTEGER, PARAMETER :: INPUT_VAR_INITIALIZED      =         1 ! 0x00000001 : Variable has been initialized
INTEGER, PARAMETER :: INPUT_VAR_EXIST            =         2 ! 0x00000002 : Found in current file
INTEGER, PARAMETER :: INPUT_VAR_PROC_IN_FILE     =         4 ! 0x00000004 : Process chain has been defined for this variable
INTEGER, PARAMETER :: INPUT_VAR_UPDATED          =         8 ! 0x00000008 : Updated by last update call
INTEGER, PARAMETER :: INPUT_VAR_READ             =        16 ! 0x00000010 : Read by last update or read var from initial file
INTEGER, PARAMETER :: INPUT_VAR_HAS_UNLIM        =        32 ! 0x00000020 : Variable contains the unlimited dimension in file
INTEGER, PARAMETER :: INPUT_VAR_MULTI_SPREAD     =        64 ! 0x00000040 : INPUT_VAR_MULTI variable is spread in different files
INTEGER, PARAMETER :: INPUT_VAR_MULTI            =       128 ! 0x00000080 : Variable is composed of more sub-vars
INTEGER, PARAMETER :: INPUT_VAR_CONT_SEND        =       256 ! 0x00000100 : Scatter sends contains only continuous data
INTEGER, PARAMETER :: INPUT_VAR_BEFORE_FIRST     =       512 ! 0x00000200 : Postpone extra reads until data available
INTEGER, PARAMETER :: INPUT_VAR_DT_SET           =      1024 ! 0x00000400 : Time step of variable has been set
INTEGER, PARAMETER :: INPUT_VAR_DEPENDENT        =      2048 ! 0x00000800 : Variable is part of a dependency chain
INTEGER, PARAMETER :: flg_var_stat_n             =        12 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_var_stat = &
  "init exist proc_chain updated read unlim multi-file multi-var continuous before dt_set depend"
! Free bits 0x7FFFF000

! Status flags (all system) stored in input_var_list%stat for group description variables
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_IFILE     =         1 ! 0x00000001 : Group has IFile specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_DT        =         2 ! 0x00000002 : Group has dt_update and/or dt_unit specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_MISS_ACT  =         4 ! 0x00000004 : Group has miss parameter specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_MIS_DIM   =         8 ! 0x00000008 : Group has mis_dim parameter specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_DATA_ACT  =        16 ! 0x00000010 : Group has data_action parameter specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_VALID_TIME=        32 ! 0x00000020 : Group has valid_time parameter specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_DEPEND    =        64 ! 0x00000040 : Group has a dependency (depend/depend_ref) specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_LAUTO     =       128 ! 0x00000080 : Group has lauto specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_NINTRP    =       256 ! 0x00000100 : Group has nIntrp specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_NATONCE   =       512 ! 0x00000200 : Group has nAtOnce specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_NFUTURE   =      1024 ! 0x00000400 : Group has nFuture specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_INDEX     =      2048 ! 0x00000800 : Group has subset_index specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_MISS_VAL  =      4096 ! 0x00001000 : Group has miss_val specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_WEIGHT    =      8192 ! 0x00002000 : Group has weight specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_WEIGHTS   =     16384 ! 0x00004000 : Group has weights specified
INTEGER, PARAMETER :: INPUT_GRPVAR_SET_SUB_MODEL =     43768 ! 0x00008000 : Group has sub_model specified
INTEGER, PARAMETER :: flg_grpvar_stat_n          =        16 ! Number of flags to display
CHARACTER (len=*), PARAMETER :: flg_grpvar_stat = "ifile dt miss_act mis_dim_act data_act valid_time depend lauto "// &
                                                  "nintrp natonce nfuture subset_index miss_val weight weights sub_model"

INTEGER, PARAMETER :: INPUT_VAR_FILE_TIME_STEP   =  UNLIKELY_VAL-1  ! Unlikely value to use for special var dt_update

! Constants for action specification in input_var_list%saction (all system set). In [] are the data to follow.
INTEGER, PARAMETER :: INPUT_VAR_DO_NOACTION         =      0 ! Dummy value indicating to terminate processing
INTEGER, PARAMETER :: INPUT_VAR_DO_MESSAGE          =      1 ! Create debug update message with time and variable name []
INTEGER, PARAMETER :: INPUT_VAR_DO_TOREAD           =      2 ! Read variable [fid,ndim in file,read_lo(1:ndim),read_sz(1:ndim)]
INTEGER, PARAMETER :: INPUT_VAR_DO_MULTI_READ       =      3 ! Hold var-fids for multi-var reads [nvars,fid(1:nvar)]
INTEGER, PARAMETER :: INPUT_VAR_DO_SUM              =      4 ! Sum over additional dim evt using weights [sum_dim,sz(sum_dim)]
INTEGER, PARAMETER :: INPUT_VAR_DO_AVG              =      5 ! Average over add dim using weights if given [sum_dim,sz(sum_dim)]
INTEGER, PARAMETER :: INPUT_VAR_DO_NORM             =      6 ! Normalize over add dim evt. using weights [sum_dim,sz(sum_dim)]
INTEGER, PARAMETER :: INPUT_VAR_DO_SPREAD           =      7 ! Spread along dims[ndim,sz(1:ndim),#sp_dim,dims(1:nsp),ncp(1:nsp)]
INTEGER, PARAMETER :: INPUT_VAR_DO_SUBSET           =      8 ! Sets a single index along mismatching dim [set/read dim,index]
INTEGER, PARAMETER :: INPUT_VAR_DO_REVERSE          =      9 ! Reverse a dimension [ndim,size(1:ndim),rev_dim]
INTEGER, PARAMETER :: INPUT_VAR_DO_SWAP             =     10 ! Swap two data blocks [blk1_sz,blk2_sz,rep]
INTEGER, PARAMETER :: INPUT_VAR_DO_PERMUTE          =     11 ! Permute variable dimensions [ndim,size(1:ndim),order(1:ndim)]
INTEGER, PARAMETER :: INPUT_VAR_DO_PACK             =     12 ! Pack variable [mul,rep]
INTEGER, PARAMETER :: INPUT_VAR_DO_UNPACK           =     13 ! Unpack variable [mul,rep]
INTEGER, PARAMETER :: INPUT_VAR_DO_EQUIV            =     14 ! Replace equivalent dimensions [mul,rep]
INTEGER, PARAMETER :: INPUT_VAR_DO_FILL_BAD         =     15 ! Replace bad values with user selected data []
INTEGER, PARAMETER :: INPUT_VAR_DO_FILL_RANGE       =     16 ! Replaces values outside valid range with selected data []
INTEGER, PARAMETER :: INPUT_VAR_DO_RESCALE          =     17 ! Rescale variable data using multiplicative scale and add offset [n]
INTEGER, PARAMETER :: INPUT_VAR_DO_DIVTIME          =     18 ! Divide cumulated data to be "flux style data" []
INTEGER, PARAMETER :: INPUT_VAR_DO_INTERPOLATE      =     19 ! Temporal interpolation []
INTEGER, PARAMETER :: INPUT_VAR_DO_BROADCAST        =     20 ! Broadcast/recieve variable []
INTEGER, PARAMETER :: INPUT_VAR_DO_SCATTER          =     21 ! Scatter/recieve variable []
INTEGER, PARAMETER :: INPUT_VAR_DO_COPY             =     22 ! Copy data between two identical buffers []
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_ADD      =     23 ! model_data = model_data + read_data
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_SUB      =     24 ! model_data = model_data - read_data
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_SUBR     =     25 ! model_data = read_data  - model_data
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_MUL      =     26 ! model_data = model_data * read_data
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_DIV      =     27 ! model_data = model_data / read_data
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_DIVR     =     28 ! model_data = read_data  / model_data
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_ADD_FLUX =     29 ! model_data = model_data + read_data * dt_model
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_NUDGE    =     30 ! model_data(i) = model_data(i)*(1 - weight) + read_data(i)*weight
INTEGER, PARAMETER :: INPUT_VAR_DO_COMBINE_NUDGE_FLD=     31 ! model_data(i) = model_data(i)*(1-weight(i))+read_data(i)*weight(i)
INTEGER, PARAMETER :: INPUT_VAR_DO_USER_ACTION      =     32 ! User action invoked [n] n=#total acp increase
CHARACTER (len=*), PARAMETER :: flg_var_spatial  = & 
     "term/skip message read multi sum avg norm spread subset reverse swap permute pack unpack equiv fill_bad fill_range " &
  // "rescale divtime interpolate broadcast scatter copy add sub subr mul div divr add_flux nudge nudge_fld user_action"
! These actions adds a buffer but do not care about the shape of this buffer. Some actions may save a buffer+data copy if the
! previous action is one of these. Care should be taken about INPUT_VAR_DO_BROADCAST which only allocates a buffer on some PEs.
! TODO: Do not also the rescaling operators (RESCALE and TIMEDIV) fit in here?
INTEGER, PARAMETER :: INPUT_VAR_1D_SET(9) = (/INPUT_VAR_DO_BROADCAST, INPUT_VAR_DO_SCATTER, INPUT_VAR_DO_TOREAD,  & !Parallel/rd
                                              INPUT_VAR_DO_REVERSE,   INPUT_VAR_DO_SWAP,    INPUT_VAR_DO_PERMUTE, & !Array manipul
                                              INPUT_VAR_DO_PACK,      INPUT_VAR_DO_UNPACK,  INPUT_VAR_DO_EQUIV/)    !Data copy

! Constants for debug output
INTEGER, PARAMETER :: INPUT_DBG_TIME_YMDHMS      =         0 ! 0x00000000 : [[+-]yy]yyyy:mm:dd hh:mm:ss
INTEGER, PARAMETER :: INPUT_DBG_TIME_YMDHM       =         1 ! 0x00000001 : [[+-]yy]yyyy:mm:dd hh:mm
INTEGER, PARAMETER :: INPUT_DBG_TIME_YMDH        =         2 ! 0x00000002 : [[+-]yy]yyyy:mm:dd hh
INTEGER, PARAMETER :: INPUT_DBG_TIME_YMD         =         3 ! 0x00000003 : [[+-]yy]yyyy:mm:dd
INTEGER, PARAMETER :: INPUT_DBG_TIME_YM          =         4 ! 0x00000004 : [[+-]yy]yyyy:mm
INTEGER, PARAMETER :: INPUT_DBG_TIME_Y           =         5 ! 0x00000005 : [[+-]yy]yyyy
INTEGER, PARAMETER :: INPUT_DBG_TIME_MASK        =         7 ! 0x00000007 : Accuracy of printed date/times
INTEGER, PARAMETER :: INPUT_DBG_AVG              =         8 ! 0x00000008 : Calculate and print averages of data passed to the model
INTEGER, PARAMETER :: INPUT_DBG_EXTERN           =        16 ! 0x00000010 : Disp files/vars for which external settings are applied
INTEGER, PARAMETER :: INPUT_DBG_DBG_FLUSH        =        32 ! 0x00000020 : Flush debug output unit after each message
INTEGER, PARAMETER :: INPUT_DBG_DBG_IO           =         0 ! 0x00000000 : Debug output is done only from the IO-PE (default)
INTEGER, PARAMETER :: INPUT_DBG_DBG_IODATA       =        64 ! 0x00000040 : Debug output is done from all IOData-PEs
INTEGER, PARAMETER :: INPUT_DBG_DBG_ALL          =       128 ! 0x00000080 : Debug output is done from all PEs
INTEGER, PARAMETER :: INPUT_DBG_QUIET            =         0 ! 0x00000000 : Don't display anything - not even warnings
INTEGER, PARAMETER :: INPUT_DBG_WARN             =       256 ! 0x00000100 : Display warnings only
INTEGER, PARAMETER :: INPUT_DBG_OPENCLOSE        =       512 ! 0x00000200 : + Opening and closing of files
INTEGER, PARAMETER :: INPUT_DBG_READ_SUMMARY     =       768 ! 0x00000300 : + Time steps where read (one line per time step)
INTEGER, PARAMETER :: INPUT_DBG_READ_GROUP       =      1024 ! 0x00000400 : + Time steps where read (one line per var group)
INTEGER, PARAMETER :: INPUT_DBG_READ_FILE        =      1280 ! 0x00000500 : + Time steps where read (one line per file)
INTEGER, PARAMETER :: INPUT_DBG_READ_VAR         =      1536 ! 0x00000600 : + Variable read
INTEGER, PARAMETER :: INPUT_DBG_UPDATE_SUMMARY   =      1792 ! 0x00000700 : + Time steps where updated (one line per time step)
INTEGER, PARAMETER :: INPUT_DBG_UPDATE_GROUP     =      2048 ! 0x00000800 : + Time steps where updated (one line per var group)
INTEGER, PARAMETER :: INPUT_DBG_UPDATE_VAR       =      2304 ! 0x00000900 : + Variable updated (interpolated)
INTEGER, PARAMETER :: INPUT_DBG_MASK             =      3840 ! 0x00000F00 : Mask for dependent flags above
INTEGER, PARAMETER :: INPUT_DBG_DEFAULT          = INPUT_DBG_READ_GROUP + INPUT_DBG_EXTERN + INPUT_DBG_AVG + INPUT_DBG_TIME_YMDH
CHARACTER (len=*), PARAMETER :: flg_dbg = &
  "time_ymdhms time_ymdhm time_ymdh time_ymd time_ym time_y avg extern flush dbg_iodata dbg_all "// &
  "quiet warn openclose read_summary read_group read_file read_var update_summary update_group update_var"

! Type for namelist specification of input variables/files
TYPE input_opt_list
  TYPE (input_opt_list), POINTER :: next                ! Next variable input entry (system)
  CHARACTER (len=128)            :: file_name           ! Name of file evt. including time macros (file)
  CHARACTER (len=128)            :: var_file_name       ! (Second) name(s) of variable(s) in file (var)
  CHARACTER (len=128)            :: dim_rename          ! Comma sep. list of dimensions to rename (file)
  CHARACTER (len=128)            :: weights             ! String with weights for sub-dimension scaling (var)
  CHARACTER (len=128)            :: def_vec             ! String with values for vector defaults (var)
  CHARACTER (len=128)            :: depend              ! Name of variables on which this variable depends (var)
  CHARACTER (len=64)             :: var_name            ! Name of variable (var)
  CHARACTER (len=32)             :: action_dim_mis      ! How to treat data with mismatching dim. (INPUT_VAR_DIM_MIS_???) (var)
  CHARACTER (len=32)             :: action_data         ! How to combine model and file data (INPUT_VAR_ACTION_???) (var)
  CHARACTER (len=16)             :: dt_update_unit      ! Unit of dt_update (var)
  CHARACTER (len=16)             :: dt_file_unit        ! Unit of dt_file (file)
  CHARACTER (len=16)             :: action_var_miss     ! Action when var missing (INPUT_VAR_MISS_???) (var)
  CHARACTER (len=16)             :: action_bad_val      ! Action when bad values are about to enter the model (var)
  CHARACTER (len=16)             :: action_time_err     ! Action when record time mismatches (INPUT_FILE_TIME_???) (file)
  CHARACTER (len=16)             :: action_cycle        ! How to cycle data (INPUT_FILE_REPEAT_???) (file)
  CHARACTER (len=16)             :: valid_time          ! At which time of month is data valid (var) 
  CHARACTER (len=16)             :: labstime            ! (Logical) Check time of record and calculate start record (file)
  CHARACTER (len=16)             :: lnametime           ! (Logical) Extract time of data from file name (file)
  CHARACTER (len=13)             :: data_start          ! Date of first existing data ([[+-]yy]yyyymmdd[hh]) (file)
  CHARACTER (len=13)             :: data_end            ! Date of last existing data ([[+-]yy]yyyymmdd[hh]) (file)
  INTEGER                        :: sub_model           ! Submodel to which this belongs (system)
  INTEGER                        :: flags               ! Bitwise flags for option object
  INTEGER                        :: dt_update           ! Time between updates of model var (var)
  INTEGER                        :: dt_file             ! Time between records in file (file)
  INTEGER                        :: offset_year         ! Year difference between model and file time (file)
  INTEGER                        :: offset_day          ! Day difference between model and file time (file)
  INTEGER                        :: offset_sec          ! Second difference between model and file time (file)
  INTEGER                        :: offset_rec          ! Record offset when calculating first record to read (file)
  INTEGER                        :: init_rec            ! Record to read for initial variables with time dimension (file)
  INTEGER                        :: cycle_length        ! Time for cycling, unit given by "action_cycle" (file)
  INTEGER                        :: subset_index        ! Index into one dimension which may be missing (var)
  REAL(dp)                       :: mul                 ! Multiplicative scale for data unit conversion (var)
  REAL(dp)                       :: add                 ! Additional offset for data unit conversion (var)
  REAL(dp)                       :: def_val             ! Default value to use if var is missing (var)
  REAL(dp)                       :: fill_val            ! Value to insert for _FillValue/missing_value (var)
  REAL(dp)                       :: file_weight         ! Weight of file data for nudging (0-1, var)
  REAL(dp)                       :: missing_val         ! Missing value to replace with fill_val/cause warning/error
  REAL(dp)                       :: valid_range(2)      ! Range of values to accept without replacing/warning/error
END TYPE input_opt_list

! Flags (all system) stored in input_opt_list%flags
INTEGER, PARAMETER :: INPUT_OPT_SPECIFICATION    =      0    ! 0x00000000 : Option hardcoded
INTEGER, PARAMETER :: INPUT_OPT_NAMELIST         =      1    ! 0x00000001 : Option specified via namelist
CHARACTER (len=*), PARAMETER :: flg_opt_type = "spec namelist"

! --- Private variables ---
! Global list type variables
TYPE (input_file_list),  POINTER :: AllInputFiles, LastFile, DummyFile      ! Global list of files, end of g.list and "dummy file"
TYPE (input_dim_list),   POINTER :: AllInputDims                            ! Global list of dimensions
TYPE (input_eqdim_list), POINTER :: AllInputEqDims                          ! Global list of equivalent dimensions
TYPE (input_eqdim_list), POINTER :: AllInputColDimDistrib                   ! Global list of dim remaps for collective dim distrb
TYPE (input_var_list),   POINTER :: AllInputVars, AllInputGrpVars, DummyVar ! Global list of vars, grp desc. vars, inactive var 
TYPE (input_data_list),  POINTER :: AllInputData                            ! Global list of shared data buffers
TYPE (input_opt_list),   POINTER :: AllInputOpts                            ! Global list of external (namelist) specifications
TYPE (input_mask_list),  POINTER :: AllInputMasks                           ! Global list of shared masks
TYPE (input_group_list), POINTER :: AllInputGroups, CurrGroup               ! Global list of variable groups
#ifdef HAVE_F2003
TYPE (input_file_type_fcns_list), POINTER :: AllInputFileTypes => NULL() ! Global list to all file types
#endif
LOGICAL                          :: ColDimAssociated = .FALSE.         ! AllInputColDimDistrib linked to end of AllInputEqDims
! Parallelisation variables
INTEGER :: n_pe                 ! Total number of PEs running the program
INTEGER :: n_iope               ! Total number of PEs reading data
INTEGER, PUBLIC :: my_pe        ! Number of this PE (0 indexed, -1 = non-parallel run)
INTEGER :: io_pe                ! Number of the master (namelists, file context checks) IO PE
LOGICAL :: lio                  ! Is this the PE which does the master IO?
LOGICAL :: liodata              ! Should this PE read data?
LOGICAL :: lioserver            ! Should this PE distribute read data?
LOGICAL :: ldbg                 ! Should this PE produce debug output
INTEGER, POINTER :: src_pe(:)   ! Which PEs delivers data to the PEs (same on all PEs, length = n_pe)
! Time variables
INTEGER :: exp_start(4)  = UNLIKELY_VAL ! Start time of experiment (i.e. start of first run in case of multiple runs) 
INTEGER :: model_time(4) = UNLIKELY_VAL ! Current time of current sub-model (after InputUpdate call, one dt ahead of main model)
INTEGER :: model_dt      = UNLIKELY_VAL ! Time step of current sub-model (in seconds)
INTEGER, POINTER :: submodel_settings(:,:) => NULL() ! Current time, time step, debug setting and log unit for each sub-model
LOGICAL :: lInLoop = .FALSE.    ! Has InputUpdate been called yet?
! Submodel variables
INTEGER :: CurrSubModel = 1     ! Number of the current sub-model
INTEGER :: NextSubModel = 0     ! Number of next new sub-model minus one
! Various
INTEGER :: mo_debug             ! Flags for debug output (INPUT_DBG_xxx above)
LOGICAL :: LogOpened = .FALSE.  ! Has mo_input opened a special file for its log output?
LOGICAL :: moInputInitialized = .FALSE. ! Has mo_input been initialized?

CONTAINS

! Dummy procedures for handling file types (parameters are attempted used to avoid compiler errors)
  INTEGER FUNCTION FileOpen(file_name,message)
  CHARACTER (len=*), INTENT(in)            :: file_name
  CHARACTER (len=*), INTENT(out), OPTIONAL :: message

    IF (PRESENT(message)) message=TRIM(file_name)
    FileOpen = 0

  END FUNCTION FileOpen

  LOGICAL FUNCTION FileClose(fid)
  INTEGER, INTENT(in) :: fid

    FileClose = fid==1

  END FUNCTION FileClose

  INTEGER FUNCTION FileGetNRec(IFile,udim)
  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_dim_list),  POINTER, OPTIONAL :: udim

    IF (PRESENT(udim)) udim => IFile%dims
    FileGetNRec = -1

  END FUNCTION FileGetNRec

  LOGICAL FUNCTION FileRead(IFile,var,fid,ndim,Info,Ofs,RecNo,Buf)
  TYPE (input_file_list), POINTER :: IFile
  TYPE (input_var_list),  POINTER :: var
  INTEGER,  INTENT(in)            :: fid
  INTEGER,  INTENT(in)            :: ndim
  INTEGER,  INTENT(in)            :: Info(:)
  INTEGER,  INTENT(in)            :: Ofs
  INTEGER,  INTENT(in)            :: RecNo
  REAL(dp), INTENT(out)           :: Buf(*)

    var%parent_file => IFile
    Buf(1:1) = fid
    FileRead = RecNo+Ofs+SUM(Info)-ndim==0

  END FUNCTION FileRead

  FUNCTION FileGetDimensions(IFile,Level)
  TYPE (input_fileobj_list), POINTER  :: FileGetDimensions
  TYPE (input_file_list),    POINTER  :: IFile
  INTEGER, INTENT(in),       OPTIONAL :: Level

    FileGetDimensions => IFile%fdims
    IF (PRESENT(Level)) IFile%fid = Level

  END FUNCTION FileGetDimensions

  FUNCTION FileGetVariables(IFile)
  TYPE (input_fileobj_list), POINTER :: FileGetVariables
  TYPE (input_file_list),    POINTER :: IFile

    FileGetVariables => IFile%fvars

  END FUNCTION FileGetVariables

  FUNCTION FileGetRecTime(IFile,RecNo,lStep)
  INTEGER :: FileGetRecTime(4)
  TYPE (input_file_list), POINTER :: IFILE
  INTEGER, INTENT(in)             :: RecNo
  LOGICAL, INTENT(in), OPTIONAL   :: lStep

    IF (PRESENT(lStep)) IFile%fid = 1 
    FileGetRecTime(:) = RecNo

  END FUNCTION FileGetRecTime

  FUNCTION FileGetAttr(Attr,IFile,var,stat)
  CHARACTER(len=256)           :: FileGetAttr
  CHARACTER(len=*), INTENT(in) :: Attr
  TYPE (input_file_list), POINTER, OPTIONAL :: IFile
  TYPE (input_var_list),  POINTER, OPTIONAL :: var
  INTEGER,          INTENT(out),   OPTIONAL :: stat

    FileGetAttr = TRIM(Attr)
    IF (PRESENT(IFile)) NULLIFY(IFile)
    IF (PRESENT(var  )) NULLIFY(var)
    IF (PRESENT(stat )) stat = 0

  END FUNCTION FileGetAttr

  ! Other template procedures

  ! Dummy function doing nothing to the data
  SUBROUTINE DataProcTemplate(level,var,curr,msg)
  INTEGER,                    INTENT(in   ) :: level
  TYPE (input_var_list),            POINTER :: var
  TYPE (input_curr),          INTENT(inout) :: curr
  CHARACTER (len=32), OPTIONAL, INTENT(out) :: msg

    IF ( IAND(level,INPUT_PROC_PROGRESS) /= 0) curr%acp = SIZE(var%saction) + 1 ! Skip over all actions
    IF ((IAND(level,INPUT_PROC_NAME) /= 0) .AND. PRESENT(msg)) msg = 'NoAction'

  END SUBROUTINE DataProcTemplate

#ifdef HAVE_F2003
  ! Although never used, this subroutine cannot be declared as an abstract interface, 
  ! since for some reason these do not accept derived types
  SUBROUTINE DataActionAddNone(Level,var,Sz,dims,curr,DimIds)
  INTEGER, INTENT(in)              :: Level ! Reason why the procedure is called
  TYPE (input_var_list),   POINTER :: var
  INTEGER,           INTENT(inout) :: Sz(:)
  TYPE (input_dim_list),   POINTER :: dims
  TYPE (input_curr), INTENT(inout) :: curr
  INTEGER,           INTENT(inout) :: DimIds(:)

    ! Reduce number of compiler warnings
    curr%acp=curr%acp
    Sz(1)=Sz(1)
    dims=>dims
    var=>var
    DimIds(1)=Level

  END SUBROUTINE DataActionAddNone
#endif

END MODULE mo_input_types
