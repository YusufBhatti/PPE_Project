!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! handling a matrix of emissioning conditions
!!
!! concepts of the routines:
!! see document: "New ECHAM5 Boundary condition scheme" by M. G. Schultz, H. Schmidt, O. Stein,
!!                                                         S. Schroeder - June 2008
!! see also: http://hammoz.icg.fz-juelich.de/data/BoundaryConditions
!!
!! @author S. Schroeder, FZ-Juelich
!!
!! $Id: 1423$
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-11)
!!
!! @par Copyright
!! 2009 by MPI-M and FZJ
!! This software is provided for non-commercial use only.
!!
MODULE mo_emi_matrix
  USE mo_kind,                     ONLY: dp
  USE mo_boundary_condition,       ONLY: BC_EVERYWHERE, BC_BOTTOM, BC_LEVEL
  USE mo_external_field_processor, ONLY: EF_INACTIVE, EF_VALUE, EF_FILE, EF_MODULE
  USE mo_external_field_processor, ONLY: EF_3D, EF_LONLAT, EF_LATLEV, EF_LEV, EF_LAT, EF_SINGLE
  USE mo_external_field_processor, ONLY: EF_TIMERESOLVED, EF_IGNOREYEAR, EF_CONSTANT, EF_NOINTER, EF_LINEAR, EF_CUBIC

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: em_read
  PUBLIC :: em_get_sector_info
  PUBLIC :: em_get_bc_from_matrix
  PUBLIC :: em_add_spec_to_sector
  PUBLIC :: em_get_SectorIndex
  PUBLIC :: maxsectors
  PUBLIC :: maxvars
!!   PUBLIC :: ematrix

  PUBLIC :: EM_NONE, EM_SURFACE, EM_VOLUME, EM_LEVEL50M, EM_FIRE

  INTEGER, PARAMETER :: maxsectors   = 30    ! max number of sectors
  INTEGER, PARAMETER :: maxvars      = 50    ! max number of variables in one sector

  INTEGER, PARAMETER :: EM_NONE     = 0
  INTEGER, PARAMETER :: EM_SURFACE  = 1
  INTEGER, PARAMETER :: EM_VOLUME   = 2
  INTEGER, PARAMETER :: EM_LEVEL50M = 3
  INTEGER, PARAMETER :: EM_FIRE     = 4

  TYPE, PUBLIC :: emi_var
  ! type of emission variable
    CHARACTER(LEN=512) :: ev_varname
    CHARACTER(LEN=512) :: ev_translation
    REAL(KIND=dp)      :: ev_factor
  END TYPE emi_var

  TYPE, PUBLIC :: emi_sector
  ! type of emission sector
    CHARACTER(LEN=64)  :: es_sectorname
    CHARACTER(LEN=512) :: es_file
    CHARACTER(LEN=512) :: es_varname
    INTEGER            :: es_nvars
    TYPE(emi_var)      :: es_variables(maxvars)
    INTEGER            :: es_eftype       ! external field type
    INTEGER            :: es_efgeometry   ! external field geometry
    INTEGER            :: es_eftimedef    ! external field time axis
    INTEGER            :: es_efinterpolate! external field time interpolation
    REAL(dp)           :: es_eftimeoffset ! external field time offset
    REAL(dp)           :: es_efvalue  
    INTEGER            :: es_emtype       ! emission (application) type
  END TYPE emi_sector

  TYPE, PUBLIC :: emi_matrix
  ! type of emission matrix
    INTEGER            :: em_nsectors
    TYPE(emi_sector)   :: em_sectors(maxsectors)
  END TYPE emi_matrix

  TYPE(emi_matrix), SAVE :: ematrix

  ! subprograms

  CONTAINS

!-----------------------------------------------------------------------
!>
!! get next data line of input file
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-11)
!!
  SUBROUTINE get_next_dataline (iunit, inline, eof)

  INTEGER, INTENT(IN)             :: iunit
  CHARACTER(LEN=512), INTENT(OUT) :: inline
  LOGICAL, INTENT(OUT)            :: eof

  INTEGER                         :: status

    status = 0
    inline ='#'
    DO WHILE ((inline(1:1) == '#') .AND. (status .GE. 0))
      read(iunit,'(a)',iostat=status) inline
      inline = ADJUSTL(inline)
      IF (LEN_TRIM(inline) == 0) inline(1:1) = '#'     ! skip empty lines
    ENDDO
    eof = (status < 0)
  END SUBROUTINE get_next_dataline

!>
!! read file with description and emission matrix
!!
!! Description?
!!
!! @par Revision History
!! code implementation by S. Schroeder (2008-11)
!!
  SUBROUTINE em_read (nsectors)
  ! should be called by p_io
  ! values should be broadcasted to other PEs after the call

  USE mo_filename,      ONLY: find_next_free_unit
  USE mo_util_string,   ONLY: toupper
  USE mo_exception,     ONLY: finish, message, em_error, em_warn, em_info, em_param, message_text
  USE mo_mpi,           ONLY: p_parallel, p_parallel_io, p_io, p_bcast
  USE mo_submodel,      ONLY: emi_basepath
  USE mo_util_string,   ONLY: separator

  INTEGER, INTENT(out)          :: nsectors

  INTEGER              :: iunit, start_ndx, tokenlen, stop_ndx, i, j, ierr
  INTEGER              :: eftype, efgeometry, eftimedef
  INTEGER              :: emtype, nvars, matvec(maxsectors), nheadersec, efinterpolate
  REAL(KIND=dp)        :: factor, efvalue, eftimeoffset
  CHARACTER(LEN=512)   :: inline, pathname, filename, cfile, ctoken
  CHARACTER(LEN=256)   :: eftypestring, domainstring, errmsg
  CHARACTER(LEN=64)    :: secname, varname, help(maxsectors+1), translation
  LOGICAL              :: eof, lfound, llast
  LOGICAL              :: lcorrect  ! try to avoid calling finish only for p_parallel_io
                                    ! (because of bunch of error messages from all other processors when being killed)

  lcorrect = .TRUE.
  CALL message('',separator)
  CALL message('','',level=em_param)
  !>>SF: add security to avoid empty basepath:
  IF (LEN(TRIM(emi_basepath)) /= 0 ) THEN
     CALL message('em_read', 'Basepath for emission data files: '//TRIM(emi_basepath), level=em_info)
  ELSE
     WRITE(message_text,'(a)') &
          'Basepath for emission data files is empty! Please set it in your submodelctl namelist!'
     CALL message('em_read', message_text, level=em_error)
     lcorrect = .FALSE.
  ENDIF
  !<<SF
  CALL message('','',level=em_param)

  ! Initialisation
  nsectors = 0
  ematrix%em_nsectors = 0
  DO i=1, maxsectors
    ematrix%em_sectors(i)%es_nvars  = 0
    ematrix%em_sectors(i)%es_emtype = -1
    ematrix%em_sectors(i)%es_eftype = -1
    ematrix%em_sectors(i)%es_efgeometry = -1
    ematrix%em_sectors(i)%es_eftimedef = -1
    ematrix%em_sectors(i)%es_eftimeoffset = 0.0_dp
    ematrix%em_sectors(i)%es_efinterpolate= -1
    ematrix%em_sectors(i)%es_efvalue = 0._dp
  END DO
  matvec(:) = -1
  eftimeoffset = 0.0_dp
  efinterpolate= EF_NOINTER

  ! read the emissions matrix on the I/O processor
  IF (p_parallel_io) THEN
    iunit = find_next_free_unit (30, 100)
    OPEN(iunit,file="emi_spec.txt",status='OLD',iostat=ierr)
    IF (ierr /= 0) THEN
      CALL message('em_read', 'emi_spec not found!',level=em_error)
      lcorrect = .FALSE.
    ELSE

    ! first section: SECTORS
      CALL get_next_dataline (iunit, inline, eof)
      IF (eof) THEN
        CALL message('em_read', 'no information found in emi_spec!',level=em_error)
        lcorrect = .FALSE.
      ELSE
        DO WHILE ((toupper(inline(1:6)) /= 'MATRIX') .AND. (.NOT. eof))
          nsectors = nsectors + 1

          ! read emission sector (short) name
          start_ndx = 1
          CALL extract_token(inline, '=', 64, start_ndx, tokenlen, ctoken)
          secname = ctoken(1:tokenlen)

          ! read boundary condition type
          CALL extract_token(inline, ',', 64, start_ndx, tokenlen, ctoken)
          eftypestring = toupper(ctoken(1:tokenlen))

          !! test validity of eftypestring
          eftype = -1
          IF (TRIM(eftypestring) == 'EF_INACTIVE') THEN
            eftype = EF_INACTIVE
            CALL message('em_read', 'Detected inactive boundary condition for sector '  &
                         //TRIM(secname)//' in emi_spec.txt!', level=em_warn)
          END IF
          IF (TRIM(eftypestring) == 'EF_VALUE' ) eftype = EF_VALUE
          IF (TRIM(eftypestring) == 'EF_MODULE') eftype = EF_MODULE
          IF (TRIM(eftypestring) == 'EF_FILE'  ) eftype = EF_FILE
          IF (TRIM(eftypestring) == '') THEN
            CALL message('em_read', 'No processing information for sector '  &
                                 //TRIM(secname)//' in emi_spec.txt! Will use EF_INACTIVE.', level=em_warn)
            eftype = EF_INACTIVE
          END IF
          IF (eftype < 0) THEN
            CALL message ('em_read', 'Invalid EF type for sector '//TRIM(secname)   &
                          //' in emi_spec.txt: "'//TRIM(eftypestring)//'"!', &
                           level=em_error)
            lcorrect = .FALSE.
            CALL get_next_dataline (iunit, inline,eof)
            CYCLE     ! ignore rest of input. Program will abort anyhow
          END IF

          !! for eftype == EF_FILE must read filename, varname and domainstring
          !! all others except EF_INACTIVE need only domainstring 
          !! The domainstring contains information related to the definition of a 
          !! boundary condition plus a keyword to identify the way how emissions
          !! from this sector shall be applied (surface, volume, level50m, etc.)
          filename=''
          varname=''
          efgeometry = -1
          eftimedef = -1
          efvalue = 0._dp
          emtype = -1

          IF (eftype == EF_FILE) THEN
            CALL extract_token(inline, ',', 512, start_ndx, tokenlen, ctoken)
            filename = ctoken(1:tokenlen)
            IF (TRIM(filename) == '') THEN
              CALL message('em_read', 'Empty filename for sector '//TRIM(secname)  &
                           //' in spite of EF_FILE!', level=em_error)
              lcorrect = .FALSE.
            END IF
            CALL extract_token(inline, ',', 64, start_ndx, tokenlen, ctoken)
            varname = ctoken(1:tokenlen)
            IF (TRIM(varname) == '') THEN
              CALL message('em_read', 'Empty varname for sector '//TRIM(secname)  &
                         //' in spite of EF_FILE!', level=em_error)
              lcorrect = .FALSE.
            END IF
          END IF

          CALL extract_token(inline, '#', 256, start_ndx, tokenlen, ctoken)
          domainstring = toupper(ctoken(1:tokenlen))
          ierr = 0
          SELECT CASE (eftype)
          CASE (EF_INACTIVE)     ! nothing to be done
          CASE (EF_FILE)
             CALL parse_ds(domainstring, ierr, errmsg, emtype, &
                           efgeometry=efgeometry, eftimedef=eftimedef, eftimeoffset=eftimeoffset, &
                           efinterpolate=efinterpolate)
          CASE (EF_VALUE)
             CALL parse_ds(domainstring, ierr, errmsg, emtype, efvalue=efvalue)
          CASE (EF_MODULE)
             CALL parse_ds(domainstring, ierr, errmsg, emtype)
          CASE DEFAULT                ! error branch: message already printed above
          END SELECT

          IF (ierr /= 0) THEN
            CALL message('em_read', errmsg, level=em_error)
            CALL message('', 'sector = '//TRIM(secname)//', domainstring = '//TRIM(domainstring))
            lcorrect = .FALSE.
          END IF

          IF (eftype >= 0) THEN
            IF (nsectors <= maxsectors) THEN
              ematrix%em_sectors(nsectors)%es_sectorname   = secname
              ematrix%em_sectors(nsectors)%es_file         = filename
              ematrix%em_sectors(nsectors)%es_varname      = varname
              ematrix%em_sectors(nsectors)%es_eftype       = eftype   
              ematrix%em_sectors(nsectors)%es_efgeometry   = efgeometry
              ematrix%em_sectors(nsectors)%es_eftimedef    = eftimedef
              ematrix%em_sectors(nsectors)%es_eftimeoffset = eftimeoffset
              ematrix%em_sectors(nsectors)%es_efinterpolate= efinterpolate
              ematrix%em_sectors(nsectors)%es_efvalue      = efvalue  
              ematrix%em_sectors(nsectors)%es_emtype       = emtype   
            ELSE
              CALL message('em_read', 'Too many sectors defined in emi_spec.txt!',level=em_error)
              lcorrect = .FALSE.
              ! nevertheless overread all following sector definitions to check the rest of emi_spec
              ! ==> don't leave loop but read until "MATRIX" or eof is read (don't end loop if (nsectors < maxsectors))
              !     (this might lead to more than one message of the above type)
            END IF
          END IF

          CALL get_next_dataline (iunit, inline,eof)
        ENDDO
        nsectors = MIN(nsectors,maxsectors)
        ematrix%em_nsectors = nsectors
        IF (ematrix%em_nsectors == 0) THEN
          CALL message('em_read', 'No sectors defined in emi_spec.txt!',level=em_error)
          lcorrect = .FALSE.
        END IF

! second section: Species-sector-matrix

        IF (eof) THEN
          CALL message('em_read', 'MATRIX keyword not found in emi_spec!',level=em_error)
          lcorrect = .FALSE.
        END IF

      ! first get the matrix header line

        CALL get_next_dataline (iunit, inline, eof)
        IF (.NOT. eof) THEN
          start_ndx = 1
          CALL condense_string(inline, tokenlen)
          i = 1
          ! get dummy string (e.g. 'SPEC' heading)
          CALL extract_token(inline, ' ', 32, start_ndx, tokenlen, ctoken, llast=llast)
          ! extract sector names and associate them with sector definition above
          DO WHILE (.NOT. llast)
            CALL extract_token(inline, ' ', 32, start_ndx, tokenlen, ctoken, llast=llast)
            DO j=1, ematrix%em_nsectors
              IF (ctoken(1:tokenlen) == TRIM(ematrix%em_sectors(j)%es_sectorname)) matvec(i)=j
            END DO
            IF (matvec(i) == -1 .AND. .NOT. llast) THEN
              CALL message('em_read', 'No sector definition for sector '//TRIM(ctoken)//'!', &
                           level=em_error)
            END IF
            i = i + 1
          END DO
          ! nheadersec = number of sectors defined in matrix header line
          ! i-2: last token with tokenlen=0 is read and additionally i is incremented
          nheadersec = i - 2
        ELSE
          CALL message('em_read', 'matrix header line (and matrix) missing in emi_spec!',level=em_error)
          lcorrect = .FALSE.
          nheadersec = 0
        END IF

        IF(ANY(matvec(1:nheadersec)==-1)) THEN
          IF (nheadersec > 0) CALL message('em_read', 'undefined sector(s) found in matrix header line',level=em_error)
          lcorrect = .FALSE.
        END IF
  
      ! now read matrix
        DO WHILE ((toupper(inline(1:5)) /= 'ALIAS') .AND. (.NOT. eof))
          CALL get_next_dataline (iunit, inline,eof)
          IF (.NOT. eof) THEN
            start_ndx = 1
            CALL condense_string(inline, tokenlen)
            i = 1
            ! get variable (i.e. species) name
            CALL extract_token(inline, ' ', 32, start_ndx, tokenlen, ctoken, llast=llast)
            varname = toupper(ctoken(1:tokenlen))
            DO WHILE (.NOT. llast)
              CALL extract_token(inline, ' ', 32, start_ndx, tokenlen, ctoken, llast=llast)
              IF (matvec(i) > 0 .AND. .NOT. llast) THEN
                IF (ematrix%em_sectors(matvec(i))%es_eftype /= EF_INACTIVE) THEN
                ! check for "not available" character(s)
                  IF (ctoken(1:1) /= '-' .AND. ctoken(1:1) /= 'X' .AND. ctoken(1:1) /= 'x') THEN
                    factor = -1.0_dp
                    READ(ctoken, *,iostat=ierr) factor
                    IF (ierr /= 0) THEN
                      CALL message('em_read', 'factor for sector '//TRIM(ematrix%em_sectors(matvec(i))%es_sectorname)//&
                                              ' of species '//TRIM(varname)//' is not a number!', level=em_error)
                      lcorrect = .FALSE.
                    END IF
                    nvars = ematrix%em_sectors(matvec(i))%es_nvars + 1
                    ematrix%em_sectors(matvec(i))%es_nvars = nvars
                    ematrix%em_sectors(matvec(i))%es_variables(nvars)%ev_varname = varname
                    ematrix%em_sectors(matvec(i))%es_variables(nvars)%ev_factor = factor
                    ematrix%em_sectors(matvec(i))%es_variables(nvars)%ev_translation = "none"
                  END IF
                END IF
              END IF
              i = i + 1
            END DO
            IF (((i-2) .lt. nheadersec) .AND. (TRIM(varname) /= 'ALIAS')) THEN
              CALL message('em_read', 'incomplete matrix line for species '//TRIM(varname), level=em_error)
              lcorrect = .FALSE.
            ENDIF
          END IF
        END DO

! third section translation table (optional!)

        DO WHILE (.NOT. eof)
          CALL get_next_dataline (iunit, inline,eof)
          IF (.NOT. eof) THEN
            start_ndx = 1
            CALL condense_string(inline, tokenlen)
            CALL extract_token(inline, ' ', 32, start_ndx, tokenlen, ctoken)
            varname = toupper(ctoken(1:tokenlen))
            CALL extract_token(inline, ' ', 32, start_ndx, tokenlen, ctoken)
            translation = ctoken(1:tokenlen)
            ! read aliassectors
            CALL extract_token(inline, '#', 256, start_ndx, tokenlen, ctoken)
            IF (tokenlen==0) THEN
              ! apply aliasname to ALL sectors of that variable
              DO i=1,nheadersec
                IF (matvec(i) > 0) THEN
                  DO j= 1,ematrix%em_sectors(matvec(i))%es_nvars
                    IF (TRIM(ematrix%em_sectors(matvec(i))%es_variables(j)%ev_varname) == TRIM(varname)) &
                      ematrix%em_sectors(matvec(i))%es_variables(j)%ev_translation=translation
                  END DO
                END IF
              END DO
            ELSE
              ! apply aliasname of that variable only to given sectors
              llast = .false.
              inline=ctoken
              start_ndx=1
              CALL condense_string(inline, tokenlen)
              DO WHILE (.NOT. llast)
                CALL extract_token(inline, ' ', 32, start_ndx, tokenlen, ctoken, llast=llast)
                DO i=1,nheadersec
                  IF (matvec(i) > 0) THEN
                    IF (TRIM(ematrix%em_sectors(matvec(i))%es_sectorname) == ctoken(1:tokenlen)) THEN
                      DO j= 1,ematrix%em_sectors(matvec(i))%es_nvars
                        IF (TRIM(ematrix%em_sectors(matvec(i))%es_variables(j)%ev_varname) == TRIM(varname)) &
                          ematrix%em_sectors(matvec(i))%es_variables(j)%ev_translation=translation
                      END DO
                    END IF
                  END IF
                END DO
              END DO
            END IF
          END IF
        END DO
      END IF  ! eof (no information)
    END IF  ! iostat (emi_spec) /= 0
  END IF  ! p_parallel_io

  IF (p_parallel) CALL p_bcast(lcorrect, p_io)
  IF (.NOT. lcorrect) CALL finish('em_read','program stopped due to errors!')
  
  !! distribute across all processors
  IF (p_parallel) THEN
    CALL p_bcast_em(ematrix, p_io)
    CALL p_bcast(nsectors,p_io)
  END IF

  END SUBROUTINE em_read

  !-----------------------------------------------------------------------------------
  ! extract a substring from an input line
  ! this routine handles white spaces
  ! start_ndx is modified to point to the next character after the delimiter or to the string end
  ! if delimiter is not found, return the remainder of the string
  SUBROUTINE extract_token(cstring, cdelim, maxlen, start_ndx, tokenlen, ctoken, lfound, llast)
  
  CHARACTER(LEN=*), INTENT(in)     :: cstring
  CHARACTER(LEN=1), INTENT(in)     :: cdelim
  INTEGER, INTENT(in)              :: maxlen
  INTEGER, INTENT(inout)           :: start_ndx
  INTEGER, INTENT(out)             :: tokenlen
  CHARACTER(LEN=maxlen), INTENT(out) :: ctoken
  LOGICAL, INTENT(out), OPTIONAL   :: lfound, llast

  CHARACTER(LEN=512)         :: ctmp
  INTEGER                    :: istop
  LOGICAL                    :: lofound, lolast

  lofound = .FALSE.
  lolast = .FALSE.
  IF (start_ndx <= 0) start_ndx = 1     ! make sure it works
  ! check for end of input string
  IF (start_ndx > LEN_TRIM(cstring)) THEN
    start_ndx = LEN_TRIM(cstring)+1
    lolast = .TRUE.
    ctoken = ''
  ELSE
    ! check for delimiter
    istop = INDEX(cstring(start_ndx:), cdelim)
    IF (istop == 0) THEN        ! delimiter not found - assume here comes the last token
      lolast = .TRUE.
      ctmp = TRIM(ADJUSTL(cstring(start_ndx:)))   
      start_ndx = LEN_TRIM(cstring)+1
    ELSE                        ! delimiter found
      lofound = .TRUE. 
      ctmp = TRIM(ADJUSTL(cstring(start_ndx:start_ndx+istop-2)))
      start_ndx = start_ndx + istop  ! points to one character behind delimiter
    END IF
    ! return token
    IF (LEN_TRIM(ctmp) > maxlen) THEN
      ctoken = ctmp(1:maxlen)
    ELSE
      ctoken = TRIM(ctmp)
    END IF
  END IF
  tokenlen = LEN_TRIM(ctoken)
  IF (PRESENT(lfound)) THEN
    lfound = lofound
  END IF
  IF (PRESENT(llast)) THEN
    llast = lolast
  END IF

  END SUBROUTINE extract_token
   
  !-----------------------------------------------------------------------------------
  ! remove duplicate blanks from a string
  ! This is needed for parsing the MATRIX section
  SUBROUTINE condense_string(cstring, len)
  
  CHARACTER(LEN=*), INTENT(inout)  :: cstring
  INTEGER, INTENT(out)             :: len

  INTEGER                    :: i, j
  LOGICAL                    :: lblank

  lblank = .FALSE.
  len = LEN_TRIM(cstring)
  IF (len == 0) RETURN

  j = 0
  ! ignore leading blanks
  i = 1
  DO WHILE (cstring(i:i) == ' ')
    i = i+1
  END DO
  DO i=1,len
    IF (lblank) THEN                 ! a blank has been found at the last step
      IF (cstring(i:i) /= ' ') THEN  ! ignore another blank
        lblank = .FALSE.
        j = j+1
        cstring(j:j) = cstring(i:i)
      END IF
    ELSE                             ! last character was not a blank -> copy
      j = j+1
      cstring(j:j) = cstring(i:i)
      IF (cstring(i:i) == ' ') lblank = .TRUE.
    END IF
  END DO
  cstring(j+1:) = ' '

  len = j
  END SUBROUTINE condense_string
  
  !-----------------------------------------------------------------------------------
  ! subroutines for string parsing (see em_read)
  SUBROUTINE parse_ds(cstring, ierr, errmsg, emtype, efgeometry, eftimedef, eftimeoffset, efinterpolate, efvalue)

  CHARACTER(LEN=*),   INTENT(in)  :: cstring
  INTEGER,            INTENT(out) :: ierr
  CHARACTER(LEN=256), INTENT(out) :: errmsg
  INTEGER,            INTENT(out) :: emtype
  INTEGER, OPTIONAL,  INTENT(out) :: efgeometry, eftimedef, efinterpolate
  REAL(dp), OPTIONAL, INTENT(out) :: efvalue, eftimeoffset

  INTEGER           :: start_ndx, hstart_ndx, tokenlen, htokenlen, i
  LOGICAL           :: llast, lfound, lnum
  CHARACTER(LEN=32) :: ctoken, chtoken

  !-- set default values
  ierr = 0
  errmsg = ''
  emtype     = EM_SURFACE
  IF (PRESENT(efgeometry)) THEN
    efgeometry = EF_LONLAT    ! this is typical for emission files
  END IF
  IF (PRESENT(eftimedef)) THEN
    eftimedef  = EF_TIMERESOLVED
  END IF
  IF (PRESENT(eftimeoffset)) THEN
    eftimeoffset  = 0._dp
  END IF
  IF (PRESENT(efvalue)) THEN
    efvalue = 0._dp
  END IF
  IF (PRESENT(efinterpolate)) THEN
    efinterpolate = EF_NOINTER
  END IF

  !-- split domainstring at commas and analyze each substring
  ! note: at least one sub string must be given
  IF (LEN_TRIM(cstring) == 0) THEN
    ierr = 1
    errmsg = 'Empty attribute string! You must at least provide a value for '  &
             //'the application mode (surface, volume, etc.)'
    RETURN
  END IF

  start_ndx = 1    ! absolute position in string
  llast = .FALSE.
  DO WHILE(.NOT. llast)
    CALL extract_token(cstring, ',', 32, start_ndx, tokenlen, ctoken, llast=llast)
    lfound = .FALSE.
    ! look for known tokens (token is uppercase because domainstring is uppercase)
    IF (PRESENT(efgeometry)) THEN
      SELECT CASE (TRIM(ctoken))
      CASE ('EF_LONLAT')
        efgeometry = EF_LONLAT
        lfound = .TRUE.
      CASE ('EF_3D')
        efgeometry = EF_3D
        lfound = .TRUE.
      CASE ('EF_LATLEV')
        efgeometry = EF_LATLEV
        lfound = .TRUE.
      CASE ('EF_LAT')
        efgeometry = EF_LAT
        lfound = .TRUE.
      CASE ('EF_LEV')
        efgeometry = EF_LEV
        lfound = .TRUE.
      CASE ('EF_SINGLE')
        efgeometry = EF_SINGLE
        lfound = .TRUE.
      CASE DEFAULT       !! nothing, but valid
      END SELECT
    END IF
    IF (PRESENT(eftimedef)) THEN
      SELECT CASE (TRIM(ctoken))
      CASE ('EF_TIMERESOLVED')
        eftimedef = EF_TIMERESOLVED
        lfound = .TRUE.
      CASE ('EF_IGNOREYEAR')
        eftimedef = EF_IGNOREYEAR
        lfound = .TRUE.
      CASE ('EF_CONSTANT')
        eftimedef = EF_CONSTANT
        lfound = .TRUE.
      CASE DEFAULT
      END SELECT
    END IF
!! skip interpolation flags for now...
!!      CASE ('EF_NOINTER')
!!      CASE ('EF_LINEAR')
!!      CASE ('EF_CUBIC')
!! we could also parse EF_TIMEINDEX (split at '=' and extract value)
!! do nothing for now
    IF (.NOT. lfound) THEN
      SELECT CASE (TRIM(ctoken))
      CASE ('SURFACE')
        emtype = EM_SURFACE
        lfound = .TRUE.
      CASE ('VOLUME')
        emtype = EM_VOLUME
        lfound = .TRUE.
      CASE ('LEVEL50M')
        emtype = EM_LEVEL50M
        lfound = .TRUE.
      CASE ('FIRE')
        emtype = EM_FIRE
        lfound = .TRUE.
      CASE DEFAULT
      END SELECT
    END IF
    IF (PRESENT(eftimeoffset)) THEN
      hstart_ndx = 1
      CALL extract_token(ctoken, '=', 32, hstart_ndx, htokenlen, chtoken)
      IF (TRIM(chtoken) == 'EF_TIMEOFFSET') THEN
        CALL extract_token(ctoken, '=', 32, hstart_ndx, htokenlen,chtoken)
        lnum = testnum(chtoken)
        IF (lnum) THEN
          READ(chtoken, *) eftimeoffset
          lfound = .TRUE.
        ELSE
          ierr = 2
          errmsg = 'Invalid floating point number! '//TRIM(chtoken)
        END IF
      END IF
    END IF
    IF (PRESENT(efinterpolate)) THEN
      hstart_ndx = 1
      CALL extract_token(ctoken, '=', 32, hstart_ndx, htokenlen, chtoken)
      IF (TRIM(chtoken) == 'EF_INTERPOLATE') THEN
        CALL extract_token(ctoken, '=', 32, hstart_ndx, htokenlen,chtoken)
        SELECT CASE (TRIM(chtoken))
        CASE ('EF_NOINTER')
          lfound = .TRUE.
        CASE ('EF_LINEAR')
          efinterpolate = EF_LINEAR
          lfound = .TRUE.
        CASE ('EF_CUBIC')
          ! still to be done ==> for the moment set to EF_LINEAR
          efinterpolate = EF_LINEAR
          lfound = .TRUE.
        CASE DEFAULT
        END SELECT
      END IF
    END IF
    IF (PRESENT(efvalue) .AND. .NOT. lfound) THEN
    ! test if token is numerical value (use better test from library??)
      lnum = testnum(ctoken)
      IF (lnum) THEN
        READ(ctoken, *) efvalue
        lfound = .TRUE.
      ELSE
        ierr = 2
        errmsg = 'Invalid floating point number! '//TRIM(ctoken)
      END IF
    END IF
    ! if nothing has been found until now, we have an invalid token string
    IF (.NOT. lfound) THEN
      ierr = 3
      errmsg = 'Undefined token in domainstring: '//TRIM(ctoken)
      RETURN
    END IF
  END DO

  !! consistency checks
  IF (PRESENT(efgeometry)) THEN
    IF (emtype == EM_SURFACE .AND. efgeometry == EF_3D) THEN
      ierr = 4
      errmsg = 'Inconsistent setting of EM_TYPE (surface) and EF_GEOMETRY (3D)!'
    END IF
  END IF
  !! more of these ..???
  END SUBROUTINE parse_ds

  FUNCTION testnum(string)    RESULT(lnum)

  CHARACTER(len=*), INTENT(in)    :: string
  LOGICAL                         :: lnum, lsign
  CHARACTER(len=16)               :: ctest
  INTEGER                         :: i, j, jmax

  lnum = .TRUE.
  ctest = '0123456789+-.'     ! must not contain exponential sign in position 1
  j = -1
  jmax = 3
  DO i=1,LEN_TRIM(string)
    lsign = (string(1:1) == '+' .OR. string(1:1) == '-')
    IF (INDEX(TRIM(ctest), string(i:i)) == 0) lnum = .FALSE.
    IF (i == 1) ctest = '0123456789.eE'      ! sign not allowed again until after 'e', 'E'
    IF (string(i:i) == '.') THEN
      ctest = '0123456789eE'    ! only one decimal point allowed
      IF (i == 1 .AND. LEN_TRIM(string) == 1) lnum = .FALSE.   ! avoid dot-only
      IF (lsign .AND. i == 2 .AND. LEN_TRIM(string) == 2) lnum = .FALSE.
    END IF
    IF (string(i:i) == 'e' .OR. string(i:i) == 'E') THEN
      ctest = '0123456789+-'
      j = 0
    END IF
    IF (j == 1 .AND. (string(i:i) == '+' .OR. string(i:i) == '-')) jmax = 4
    IF (j > 0) ctest = '0123456789'
    IF (j >= 0) j = j+1              ! count digits of exponent
    IF (j > jmax) lnum = .FALSE.        ! overflow
  END DO

  END FUNCTION testnum

  !! broadcast emission matrix structure
  SUBROUTINE p_bcast_em (em_struc, p_source, comm)
  USE mo_mpi, ONLY: p_bcast, p_all_comm

  TYPE(emi_matrix),  INTENT(INOUT) :: em_struc
  INTEGER,           INTENT(in)    :: p_source
  INTEGER, OPTIONAL, INTENT(in)    :: comm

  INTEGER :: p_comm, jsec, jvar

  IF (PRESENT(comm)) THEN
     p_comm = comm
  ELSE
     p_comm = p_all_comm
  ENDIF

  CALL p_bcast(em_struc%em_nsectors,    p_source, p_comm)
  DO jsec=1, em_struc%em_nsectors
    CALL p_bcast(em_struc%em_sectors(jsec)%es_sectorname,   p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_file,         p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_varname,      p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_nvars,        p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_eftype,       p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_efgeometry,   p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_eftimedef,    p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_eftimeoffset, p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_efinterpolate,p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_efvalue,      p_source, p_comm)
    CALL p_bcast(em_struc%em_sectors(jsec)%es_emtype,       p_source, p_comm)
    DO jvar=1, em_struc%em_sectors(jsec)%es_nvars
      CALL p_bcast(em_struc%em_sectors(jsec)%es_variables(jvar)%ev_varname,     p_source, p_comm)
      CALL p_bcast(em_struc%em_sectors(jsec)%es_variables(jvar)%ev_translation, p_source, p_comm)
      CALL p_bcast(em_struc%em_sectors(jsec)%es_variables(jvar)%ev_factor,      p_source, p_comm)
    END DO
  END DO

  END SUBROUTINE p_bcast_em


  SUBROUTINE em_get_sector_info(i, sectorname, nvars)

  USE mo_exception,   ONLY: message, em_error, message_text

  INTEGER, INTENT(in)             :: i
  CHARACTER(LEN=64), INTENT(out)  :: sectorname
  INTEGER, INTENT(out)            :: nvars

  IF (i > ematrix%em_nsectors) THEN
     WRITE(message_text,'(2(a,i0),a)') 'Sector index',i,' larger than',ematrix%em_nsectors,'! Should not happen!!'
     CALL message('em_get_sector_info', message_text, level=em_error)
  END IF

  sectorname = ematrix%em_sectors(i)%es_sectorname
  nvars      = ematrix%em_sectors(i)%es_nvars

  END SUBROUTINE

  SUBROUTINE em_get_bc_from_matrix( i, j, sectorname, bc_struc, emtype, emfactor, varname)

  USE mo_control,                  ONLY: nlev
  USE mo_exception,                ONLY: message, em_error
  USE mo_boundary_condition,       ONLY: bc_nml, BC_EVERYWHERE, BC_BOTTOM, BC_LEVEL
  USE mo_boundary_condition,       ONLY: BC_REPLACE
  USE mo_filename,                 ONLY: str_filter
  USE mo_external_field_processor, ONLY: EF_NOINTER
  USE mo_submodel,                 ONLY: emi_basepath

  INTEGER,          INTENT(IN)   :: i, j
  CHARACTER(LEN=*), INTENT(IN)   :: sectorname
  TYPE(bc_nml),     INTENT(OUT)  :: bc_struc
  INTEGER,          INTENT(out)  :: emtype
  REAL(dp),         INTENT(out)  :: emfactor  ! emission scaling factor
  CHARACTER(LEN=64), INTENT(out) :: varname   ! actual name of a species

  CHARACTER(LEN=512) :: cfile, help
  INTEGER            :: ndx, ndx_c, start_ndx, stop_ndx
  LOGICAL            :: found

  !-- set default values
  emfactor = 1._dp
  bc_struc%ef_interpolate = EF_NOINTER      !! add to parse_ds later (?)
  bc_struc%ef_actual_unit= ''               !! add to parse_ds later (?)
  bc_struc%bc_mode = BC_REPLACE             !! needed for extended emi diagnostics

  !-- fill easy stuff
  emtype = ematrix%em_sectors(i)%es_emtype
  varname = ematrix%em_sectors(i)%es_variables(j)%ev_varname
  bc_struc%ef_type       = ematrix%em_sectors(i)%es_eftype
  bc_struc%ef_geometry   = ematrix%em_sectors(i)%es_efgeometry
  bc_struc%ef_timedef    = ematrix%em_sectors(i)%es_eftimedef
  bc_struc%ef_timeoffset = ematrix%em_sectors(i)%es_eftimeoffset
  bc_struc%ef_interpolate= ematrix%em_sectors(i)%es_efinterpolate
  bc_struc%ef_value      = ematrix%em_sectors(i)%es_efvalue

  !-- construct correct file name template
  IF (bc_struc%ef_type == EF_FILE) THEN
    IF (LEN_TRIM(ematrix%em_sectors(i)%es_file) == 0) THEN
      CALL message('em_get_bc_from_matrix', 'Empty file name for EF_FILE in ' &
                   //TRIM(ematrix%em_sectors(i)%es_sectorname), level=em_error)
    END IF
    ndx_c = INDEX(ematrix%em_sectors(i)%es_file,'%C')
    IF (ndx_c > 0) THEN
      help = ematrix%em_sectors(i)%es_file
      help(ndx_c:ndx_c)='$'
      ndx=INDEX(help,'%')
      DO WHILE (ndx > 0)
        help(ndx:ndx)='$'
        ndx=INDEX(help,'%')
      ENDDO
      help(ndx_c:ndx_c)='%'
      IF (TRIM(ematrix%em_sectors(i)%es_variables(j)%ev_translation) /= "none") THEN
        cfile = str_filter(help,0, 0, 0, 0, 0, 0,0,TRIM(ematrix%em_sectors(i)%es_variables(j)%ev_translation),"T00","L00")
      ELSE
        cfile = str_filter(help,0, 0, 0, 0, 0, 0,0,TRIM(ematrix%em_sectors(i)%es_variables(j)%ev_varname),"T00","L00")
      ENDIF
      ndx=INDEX(cfile,'$')
      DO WHILE (ndx > 0)
        cfile(ndx:ndx)='%'
        ndx=INDEX(cfile,'$')
      ENDDO

      IF (cfile(1:1) == '/') THEN
        bc_struc%ef_template = TRIM(cfile)
      ELSE
        bc_struc%ef_template = TRIM(emi_basepath)//'/'//TRIM(cfile)
      ENDIF
    ELSE
      IF (ematrix%em_sectors(i)%es_file(1:1) == '/') THEN
        bc_struc%ef_template = TRIM(ematrix%em_sectors(i)%es_file)
      ELSE
        bc_struc%ef_template = TRIM(emi_basepath)//'/'//TRIM(ematrix%em_sectors(i)%es_file)
      ENDIF
    ENDIF
  ENDIF

  !-- construct bc_domain value
  ! note: except for "volume" emissions (i.e. aircraft, volcanoes) all boundary conditions
  ! are defined as BC_BOTTOM, i.e. as 2D fields. Vertical distribution takes place where
  ! necessary in emi_interface in mo_emi_interface.
  SELECT CASE (ematrix%em_sectors(i)%es_emtype)
  CASE (EM_SURFACE)
    bc_struc%bc_domain = BC_BOTTOM
  CASE (EM_VOLUME)
    bc_struc%bc_domain = BC_EVERYWHERE
    !! WARNING: this might be wrong (e.g. aircraft)
    !! ==> if read from file: BC_ALTITUDE || BC_PRESSURE is set accordingly after reading
  CASE (EM_LEVEL50M)
    bc_struc%bc_domain = BC_BOTTOM
!   bc_struc%bc_domain = BC_LEVEL
!   bc_struc%bc_minlev = nlev-1
!   bc_struc%bc_maxlev = nlev-1
  CASE (EM_FIRE)
     !>>dod allow 3-D fire emissions
     IF (bc_struc%ef_geometry == EF_3D) THEN
        bc_struc%bc_domain = BC_EVERYWHERE
     ELSE
        bc_struc%bc_domain = BC_BOTTOM
     END IF
  END SELECT

  !!--- replace %C0 in varname field by species name
  !! must be made more flexible to allow for "emis_%C0" for example  ### mgs
  IF (index(ematrix%em_sectors(i)%es_varname,'%') > 0) THEN
    IF (trim(ematrix%em_sectors(i)%es_variables(j)%ev_translation) == "none") THEN
      bc_struc%ef_varname = trim(ematrix%em_sectors(i)%es_variables(j)%ev_varname)
    ELSE
      bc_struc%ef_varname = trim(ematrix%em_sectors(i)%es_variables(j)%ev_translation)
    ENDIF
  ELSE
    bc_struc%ef_varname = trim(ematrix%em_sectors(i)%es_varname)
  ENDIF
  emfactor = ematrix%em_sectors(i)%es_variables(j)%ev_factor
  END SUBROUTINE em_get_bc_from_matrix

  SUBROUTINE em_add_spec_to_sector(i, factor, varname)
  INTEGER,           INTENT(IN) :: i
  REAL(dp),          INTENT(IN) :: factor  ! emission scaling factor
  CHARACTER(LEN=64), INTENT(IN) :: varname   ! actually name of a species

  INTEGER :: nvars

  nvars = ematrix%em_sectors(i)%es_nvars + 1
  ematrix%em_sectors(i)%es_nvars  = nvars

  ematrix%em_sectors(i)%es_variables(nvars)%ev_varname = varname
  ematrix%em_sectors(i)%es_variables(nvars)%ev_factor = factor

  END SUBROUTINE em_add_spec_to_sector

  FUNCTION em_get_SectorIndex (name) RESULT(index)
  USE mo_exception,                ONLY: message, em_error
  CHARACTER(LEN=*), INTENT(IN) :: name
  INTEGER :: index

  LOGICAL :: found
  INTEGER :: i

! Determine index of sector in matrix

  i = 0
  found = .false.
  DO WHILE ((.NOT. found) .AND. (i <= ematrix%em_nsectors))
    i = i + 1
    IF (ematrix%em_sectors(i)%es_sectorname == name) found = .true.
  ENDDO
  IF (i .gt. ematrix%em_nsectors) CALL message('em_get_bc_from_matrix', 'no sector '//trim(name)//' defined!',level=em_error)
  index = i
  END FUNCTION em_get_SectorIndex

END MODULE mo_emi_matrix
