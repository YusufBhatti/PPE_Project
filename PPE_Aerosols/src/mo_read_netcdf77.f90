!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_read_netcdf77
  !
  ! Description: 
  !
  ! subroutines for reading (hyperslabs) of variables from netcdf files
  ! order of dimensions in netcdf files and fortran variable can be different
  !
  ! Authors:
  !
  ! J.S. Rast, MPI, May 2005, change all netcdf routines to fortran 77 calls
  ! D.O'Donnell, MPI-M, Feb 2008, added read 4-D variable

  USE mo_kind,      ONLY: dp
  USE mo_exception, ONLY: finish  
  USE mo_netcdf,    ONLY: nf_check, nf_max_var_dims, nf_open, nf__open,  &
                          nf_close, chunksize, nf_nowrite,               &
                          nf_inq_dimlen, nf_inq_dimid, nf_inq_varid,     &
                          nf_get_vara_double, nf_inq_varndims,           &
                          nf_inq_nvars, nf_inq_vardimid, nf_inq_varname, &
                          nf_get_var_double

  IMPLICIT NONE

  PRIVATE

  PUBLIC    :: read_var_nf77_0d          ! read 0D-variable
  PUBLIC    :: read_var_nf77_1d          ! read 1D-variable 
  PUBLIC    :: read_var_nf77_2d          ! read 2D-variable 
  PUBLIC    :: read_var_nf77_3d          ! read 3D-variable 
  PUBLIC    :: read_var_nf77_4d          ! read 4D-variable 
  PUBLIC    :: read_sumvar_nf77_2d       ! read sum of hyperslab of 3D-variables (third dim. has to 
  ! have length 12, sum of all such 3D-variables is calc.
  PUBLIC    :: read_var_hs_nf77_0d
  PUBLIC    :: read_var_hs_nf77_1d       ! read hyperslab of 2D-variable (second dim. is time)
  PUBLIC    :: read_var_hs_nf77_2d       ! read hyperslab of 3D-variable (third dim. is time)
  PUBLIC    :: read_var_hs_nf77_3d       ! read hyperslab of 4D-variable (fourth dim. is time)
  PUBLIC    :: search_var_nf77           ! searches for variable in file
  PUBLIC    :: read_diml_nf77            ! gives length of a dimension   

CONTAINS
!===============================================================================
  SUBROUTINE read_var_nf77_0d (file_name, varname, var, ierr)

    ! Description:
    ! 
    ! read variable varname from file file_name
    !
    ! Method:
    ! 
    ! Open the netCDF file named file_name and read varname into var.
    ! The routine is to be intended to run on the I/O-processor
    ! 
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: varname   ! name of variable to be read
    ! REAL(dp), INTENT(out)            :: var       ! variable value on exit
    ! INTEGER, INTENT(out)             :: ierr      ! ierr is 0 if variable varname 
    !    was sucessfully read, 1 otherwise
    !
    ! Authors:
    ! 
    ! J.S. Rast, MPI, January 2010, original source
    !
    ! Arguments:
    !
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: varname 
    REAL(dp), INTENT(out)                 :: var
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, znvar
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid), fname=TRIM(file_name))
! inquire variable and number of dimensions
    ierr=nf_inq_varid(zncid, TRIM(varname), znvar)
    IF (ierr == 0) THEN
       CALL nf_check(nf_get_var_double(zncid, znvar, var), fname=TRIM(file_name))
    ELSE
       CALL finish('read_var_nf77_0d:', &
               'variable '//TRIM(varname)//' not found in '//TRIM(file_name))
    END IF
    CALL nf_check(nf_close(zncid))
  END SUBROUTINE read_var_nf77_0d
!===============================================================================
!===============================================================================
  SUBROUTINE read_var_nf77_1d (file_name, dimname, varname, varptr, ierr)

    ! Description:
    !
    ! read in from a netCDF file variable varname(dimname).
    !
    ! Method:
    ! 
    ! Open the netCDF file named file_name and read in varptr varname(dimname).
    ! The routine is to be intended to run on the I/O-processor
    ! 
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname ! name of dimension in the netCDF
    !    file which shall appear as the dimension of varptr.
    ! CHARACTER (LEN = *), INTENT (in) :: varname   ! name of variable to be read
    ! REAL(dp), INTENT(out)            :: varptr(:) ! variable values on exit, has to have
    !    dimension exactly as the lengths of dimname in netCDF file
    ! INTEGER, INTENT(out)             :: ierr ! ierr is 0 if variable varname is successfully 
    !    read, 1 otherwise
    !
    ! Authors:
    ! 
    ! J.S. Rast, MPI, August 2004, original source
    !
    ! Arguments:
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname
    CHARACTER (LEN = *), INTENT (in)      :: varname 
    REAL(dp), INTENT(out)                 :: varptr(:)
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, zdimlength, znvar, zvardims
    INTEGER                               :: zdims, zdimid
    zdims = SIZE (varptr)
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid), fname=TRIM(file_name))
! verify lenght of dimension
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname), zdimid), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid, zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims) THEN
       CALL finish('read_var_nf77_1d:', 'wrong length of dim in file '//TRIM(file_name))
    END IF
! inquire variable and number of dimensions
    ierr=nf_inq_varid(zncid, TRIM(varname), znvar)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_varndims(zncid, znvar, zvardims), fname=TRIM(file_name))
       IF (zvardims /= 1) THEN
          CALL finish('read_var_nf77_1d:', &
               'wrong number of dimension of variable '//TRIM(varname))
       END IF
       CALL nf_check(nf_get_vara_double(zncid, znvar, (/1/), (/zdims/), varptr), fname=TRIM(file_name))
    ELSE
       CALL finish('read_var_nf77_1d:', &
               'variable '//TRIM(varname)//' not found in '//TRIM(file_name))
    END IF
  END SUBROUTINE read_var_nf77_1d
!===============================================================================
!===============================================================================
  SUBROUTINE read_var_nf77_2d (file_name, dimname1, dimname2, var_name, varptr, ierr)

    ! Description:
    !
    ! read in a variable from a netCDF file ******************(not yet tested)*******************
    !
    ! Method:
    ! 
    ! Open the netCDF file named file_name and read in varptr the values of var_name
    ! dimname1,dimname2 are the names of the dimension in the order as they will
    ! be in varptr. The order of the dimensions in the netCDF file is irrelevant.
    ! The routine is to be intended to run on the I/O-processor
    ! 
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname1 ! name of dimension in the netCDF
    !    file which shall appear as the first dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname2 ! name of dimension in the netCDF
    !    file which shall appear as the second dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: var_name ! name of the variable to be read
    ! REAL(dp), INTENT(out)            :: varptr(:,:,:) ! variable values on exit, has to
    !    dimension exactly as the lengths of dimname1,...,dimname3 in netCDF file
    ! INTEGER, INTENT(out) :: ierr ! ierr is 0 if variable var_name is successfully 
    !    read, 1 if it was not found in file file_name
    !
    ! Authors:
    ! 
    ! J.S. Rast, MPI, December 2002, original source
    !
    ! Arguments:
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname1
    CHARACTER (LEN = *), INTENT (in)      :: dimname2
    CHARACTER (LEN = *), INTENT (in)      :: var_name
    REAL(dp), INTENT(out)                 :: varptr(:,:)
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, zdimlength, zvarid, zvardims, i, j
    INTEGER, DIMENSION(2)                 :: zdims, zdimid, zorder, zcountvar
    INTEGER, DIMENSION(NF_MAX_VAR_DIMS)   :: zvardimids
    REAL(dp), POINTER                     :: zin(:,:)

! define dimension length vector
    zdims(1) = SIZE (varptr, DIM = 1)
    zdims(2) = SIZE (varptr, DIM = 2)
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid), fname=TRIM(file_name))
! verify lenghts of dimensions
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname1), zdimid(1)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(1), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(1)) THEN
       CALL finish('read_var_nf77_2d:', 'wrong length of dim1 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname2), zdimid(2)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(2), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(2)) THEN
       CALL finish('read_var_nf77_2d:', 'wrong length of dim2 in initial file '//TRIM(file_name))
    END IF
    ierr=nf_inq_varid(zncid, TRIM(var_name), zvarid)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_varndims(zncid, zvarid, zvardims), fname=TRIM(file_name))
       CALL nf_check(nf_inq_vardimid(zncid, zvarid, zvardimids), fname=TRIM(file_name))
       IF (zvardims /= 2) THEN
          CALL finish('read_var_nf77_2d:', &
               'wrong number of dimension of variable '//TRIM(var_name))
       END IF
       DO i = 1, zvardims
          DO j = 1, zvardims
             IF (zvardimids(j) == zdimid(i)) THEN
                zorder(i)=j
                EXIT
             END IF
          END DO
       END DO
       ALLOCATE(zin(zdims(zorder(1)), zdims(zorder(2))))
       zcountvar(1)=zdims(zorder(1))
       zcountvar(2)=zdims(zorder(2))
       CALL nf_check(nf_get_vara_double(zncid, zvarid, (/1,1/), zcountvar, zin), fname=TRIM(file_name))
       CALL nf_check(nf_close(zncid), fname=TRIM(file_name))
       IF (ALL(zorder == (/ 1, 2 /))) THEN
         varptr = zin         
       ELSE
         varptr = RESHAPE(zin, zdims, order=zorder)
       ENDIF
       IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    ELSE
       CALL finish('read_var_nf77_2d:', &
               'variable '//TRIM(var_name)//' not found in '//TRIM(file_name))
    END IF
    END SUBROUTINE read_var_nf77_2d
!===============================================================================
!===============================================================================
  SUBROUTINE read_var_nf77_3d (file_name, dimname1, dimname2, dimname3, var_name, varptr, ierr)

    ! Description:
    !
    ! read in a variable from a netCDF file
    !
    ! Method:
    ! 
    ! Open the netCDF file named file_name and read in varptr the values of var_name
    ! dimname1,...,dimname3 are the names of the dimension in the order as they will
    ! be in varptr. The order of the dimensions in the netCDF file is irrelevant.
    ! The routine is to be intended to run on the I/O-processor
    ! 
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname1 ! name of dimension in the netCDF
    !    file which shall appear as the first dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname2 ! name of dimension in the netCDF
    !    file which shall appear as the second dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname3 ! name of dimension in the netCDF
    !    file which shall appear as the third dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: var_name ! name of the variable to be read
    ! REAL(dp), INTENT(out) :: varptr(:,:,:) ! variable values on exit, has to
    !    dimension exactly as the lengths of dimname1,...,dimname3 in netCDF file
    ! INTEGER, INTENT(out) :: ierr ! ierr is 0 if variable var_name is successfully 
    !    read, 1 if it was not found in file file_name
    !
    ! Authors:
    ! 
    ! J.S. Rast, MPI, December 2002, original source
    !
    ! Arguments:
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname1
    CHARACTER (LEN = *), INTENT (in)      :: dimname2
    CHARACTER (LEN = *), INTENT (in)      :: dimname3
    CHARACTER (LEN = *), INTENT (in)      :: var_name
    REAL(dp), INTENT(out)                 :: varptr(:,:,:)
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, zdimlength, zvarid, zvardims, i, j
    INTEGER, DIMENSION(3)                 :: zdims, zdimid, zorder, zcountvar
    INTEGER, DIMENSION(NF_MAX_VAR_DIMS)   :: zvardimids
    REAL(dp), POINTER                     :: zin(:,:,:)

! define dimension length vector
    zdims(1) = SIZE (varptr, DIM = 1)
    zdims(2) = SIZE (varptr, DIM = 2)
    zdims(3) = SIZE (varptr, DIM = 3)
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid), fname=TRIM(file_name))
! verify lenghts of dimensions
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname1), zdimid(1)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(1), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(1)) THEN
       CALL finish('rad_var_nf77_3d:', 'wrong length of dim1 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname2), zdimid(2)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(2), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(2)) THEN
       CALL finish('rad_var_nf77_3d:', 'wrong length of dim2 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname3), zdimid(3)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(3), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(3)) THEN
       CALL finish('rad_var_nf77_3d:', 'wrong length of dim3 in initial file '//TRIM(file_name))
    END IF
    ierr=nf_inq_varid(zncid, TRIM(var_name), zvarid)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_varndims(zncid, zvarid, zvardims), fname=TRIM(file_name))
       CALL nf_check(nf_inq_vardimid(zncid, zvarid, zvardimids), fname=TRIM(file_name))
       IF (zvardims /= 3) THEN
          CALL finish('rad_var_nf77_3d:', &
               'wrong number of dimension of variable '//TRIM(var_name))
       END IF
       DO i = 1, zvardims
          DO j = 1, zvardims
             IF (zvardimids(j) == zdimid(i)) THEN
                zorder(i)=j
                EXIT
             END IF
          END DO
       END DO
       ALLOCATE(zin(zdims(zorder(1)), zdims(zorder(2)), zdims(zorder(3))))
       zcountvar(1)=zdims(zorder(1))
       zcountvar(2)=zdims(zorder(2))
       zcountvar(3)=zdims(zorder(3))
       CALL nf_check(nf_get_vara_double(zncid, zvarid, (/1,1,1/), zcountvar, zin), fname=TRIM(file_name))
       CALL nf_check(nf_close(zncid), fname=TRIM(file_name))
       IF (ALL(zorder == (/ 1, 2, 3 /))) THEN
         varptr = zin         
       ELSE
         varptr = RESHAPE(zin, zdims, order=zorder)
       ENDIF
       IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    ELSE
       CALL finish('read_var_nf77_3d:', &
               'variable '//TRIM(var_name)//' not found in '//TRIM(file_name))
    END IF
    END SUBROUTINE read_var_nf77_3d
!===============================================================================
!===============================================================================
    SUBROUTINE read_var_nf77_4d (file_name, dimname1, dimname2, dimname3, dimname4, &
                               var_name, varptr, ierr)

    ! Description:
    !
    ! read in a variable from a netCDF file
    !
    ! Method:
    ! 
    ! Open the netCDF file named file_name and read in varptr the values of var_name
    ! dimname1,...,dimname4 are the names of the dimension in the order as they will
    ! be in varptr. The order of the dimensions in the netCDF file is irrelevant.
    ! The routine is to be intended to run on the I/O-processor
    ! 
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname1 ! name of dimension in the netCDF
    !    file which shall appear as the first dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname2 ! name of dimension in the netCDF
    !    file which shall appear as the second dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname3 ! name of dimension in the netCDF
    !    file which shall appear as the third dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname4 ! name of dimension in the netCDF
    !    file which shall appear as the fourth dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: var_name ! name of the variable to be read
    ! REAL(dp), INTENT(out) :: varptr(:,:,:,:) ! variable values on exit, has to
    !    dimension exactly as the lengths of dimname1,...,dimname4 in netCDF file
    ! INTEGER, INTENT(out) :: ierr ! ierr is 0 if variable var_name is successfully 
    !    read, 1 if it was not found in file file_name
    !
    ! Authors:
    ! 
    ! D. O'Donnell, MPI-M, Feb 2008, based on other read_var_... subroutines by J.S.Rast
    !
    !
    ! Arguments:
    !
    
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname1
    CHARACTER (LEN = *), INTENT (in)      :: dimname2
    CHARACTER (LEN = *), INTENT (in)      :: dimname3
    CHARACTER (LEN = *), INTENT (in)      :: dimname4
    CHARACTER (LEN = *), INTENT (in)      :: var_name
    REAL(dp), INTENT(out)                 :: varptr(:,:,:,:)
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, zdimlength, zvarid, zvardims, i, j
    INTEGER, DIMENSION(4)                 :: zdims, zdimid, zorder, zcountvar
    INTEGER, DIMENSION(NF_MAX_VAR_DIMS)   :: zvardimids
    REAL(dp), POINTER                     :: zin(:,:,:,:)

! define dimension length vector
    zdims(1) = SIZE (varptr, DIM = 1)
    zdims(2) = SIZE (varptr, DIM = 2)
    zdims(3) = SIZE (varptr, DIM = 3)
    zdims(4) = SIZE (varptr, DIM = 4)
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid), fname=TRIM(file_name))
! verify lenghts of dimensions
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname1), zdimid(1)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(1), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(1)) THEN
       CALL finish('read_var_nf77_4d:', 'wrong length of dim1 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname2), zdimid(2)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(2), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(2)) THEN
       CALL finish('read_var_nf77_4d:', 'wrong length of dim2 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname3), zdimid(3)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(3), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(3)) THEN
       CALL finish('read_var_nf77_4d:', 'wrong length of dim3 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname4), zdimid(4)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(4), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(4)) THEN
       CALL finish('read_var_nf77_4d:', 'wrong length of dim4 in initial file '//TRIM(file_name))
    END IF
    ierr=nf_inq_varid(zncid, TRIM(var_name), zvarid)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_varndims(zncid, zvarid, zvardims), fname=TRIM(file_name))
       CALL nf_check(nf_inq_vardimid(zncid, zvarid, zvardimids), fname=TRIM(file_name))
       IF (zvardims /= 4) THEN
          CALL finish('read_var_nf77_4d:', &
               'wrong number of dimension of variable '//TRIM(var_name))
       END IF
       DO i = 1, zvardims
          DO j = 1, zvardims
             IF (zvardimids(j) == zdimid(i)) THEN
                zorder(i)=j
                EXIT
             END IF
          END DO
       END DO
       ALLOCATE(zin(zdims(zorder(1)), zdims(zorder(2)), zdims(zorder(3)), zdims(zorder(4)) ))
       zcountvar(1)=zdims(zorder(1))
       zcountvar(2)=zdims(zorder(2))
       zcountvar(3)=zdims(zorder(3))
       zcountvar(4)=zdims(zorder(4))
       CALL nf_check(nf_get_vara_double(zncid, zvarid, (/1,1,1,1/), zcountvar, zin), fname=TRIM(file_name))
       CALL nf_check(nf_close(zncid), fname=TRIM(file_name))
       IF (ALL(zorder == (/ 1, 2, 3, 4 /))) THEN
         varptr = zin         
       ELSE
         varptr = RESHAPE(zin, zdims, order=zorder)
       ENDIF
       IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    ELSE
       CALL finish('read_var_nf77_4d:', &
               'variable '//TRIM(var_name)//' not found in '//TRIM(file_name))
    END IF
  END SUBROUTINE read_var_nf77_4d
!===============================================================================
!===============================================================================
  SUBROUTINE read_sumvar_nf77_2d (file_name, dimname1, dimname2, dimname3, imonth, varptr, ierr)

    ! Description:
    !
    ! read in from a netCDF file the sum of all variables containing dimname1 
    ! and dimname2 and a time dimension of size 12
    !
    ! Method:
    ! 
    ! Open the netCDF file named file_name and read in varptr the sum of the
    ! values of all variables containing the dimensions dimname1, dimname2 
    ! which are the names of the dimensions in the order 
    ! as they will be in varptr. The order of the dimensions in the netCDF file is 
    ! irrelevant.
    ! The routine is to be intended to run on the I/O-processor
    ! 
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname1 ! name of dimension in the netCDF
    !    file which shall appear as the first dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname2 ! name of dimension in the netCDF
    !    file which shall appear as the second dimension of varptr. var_name can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname3 ! name of the time dimension
    !    currently, there must be 12 values for the 12 months
    ! INTEGER, INTENT (in)             :: imonth   ! number of month in the
    !    year (range 1..12)
    ! REAL(dp), INTENT(out) :: varptr(:,:) ! variable values on exit, has to
    !    dimension exactly as the lengths of dimname1,dimname2 in netCDF file
    ! INTEGER, INTENT(out) :: ierr ! ierr is 0 if variable var_name is successfully 
    !    read, 1 otherwise
    !
    ! Authors:
    ! 
    ! J.S. Rast, MPI, January 2003, original source
    !
    ! Arguments:
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname1
    CHARACTER (LEN = *), INTENT (in)      :: dimname2
    CHARACTER (LEN = *), INTENT (in)      :: dimname3
    INTEGER, INTENT(in)                   :: imonth
    REAL(dp), INTENT(out)                 :: varptr(:,:)
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, znvar, zdimlength, zvardims, i, j, ivar
    INTEGER, DIMENSION(3)                 :: zdims, zdimid, zvardimids, zorder, zstart, zcount, &
                                             zstartvar, zcountvar
    LOGICAL, DIMENSION(3)                 :: zorderlog 
    REAL(dp), POINTER                     :: zin(:,:,:)
!,zinre(:,:,:)
    character (len = 256)                 :: zvarname

! define dimension length vector
    zdims(1) = SIZE (varptr, DIM = 1)
    zdims(2) = SIZE (varptr, DIM = 2)
    zdims(3) = 1 ! the dimension length is 12, but we only read one month
! define start vector of dimensions and number of indices in each dimension to be read
    zstart(1) = 1
    zstart(2) = 1
    zstart(3) = imonth
    zcount(1) = zdims(1)
    zcount(2) = zdims(2)
    zcount(3) = 1
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid),fname=TRIM(file_name))
! verify lenghts of dimensions
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname1), zdimid(1)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(1), zdimlength),fname=TRIM(file_name))
    IF (zdimlength /= zdims(1)) THEN
       CALL finish('read_sumvar_nf77_2d:', 'wrong length of dim1 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname2), zdimid(2)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(2), zdimlength),fname=TRIM(file_name))
    IF (zdimlength /= zdims(2)) THEN
       CALL finish('read_sumvar_nf77_2d:', 'wrong length of dim2 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname3), zdimid(3)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(3), zdimlength),fname=TRIM(file_name))
    IF (zdimlength /= 12) THEN
       CALL finish('read_sumvar_nf77_2d:', 'wrong length of dim3 in initial file '//TRIM(file_name))
    END IF
! determine number of variables in netCDF file
    CALL nf_check(nf_inq_nvars(zncid, znvar),fname=TRIM(file_name))
    IF (znvar == 0) THEN
       ierr=1
       RETURN
    END IF
    varptr = 0._dp
    DO ivar=1,znvar
       CALL nf_check(nf_inq_varndims(zncid,ivar,zvardims),fname=TRIM(file_name))
       CALL nf_check(nf_inq_varname(zncid,ivar,zvarname),fname=TRIM(file_name))
       IF (zvardims == 3) THEN
          CALL nf_check(nf_inq_vardimid(zncid,ivar,zvardimids),fname=TRIM(file_name))
          zorderlog(:)=.false.
          DO i = 1, zvardims
             DO j = 1, zvardims
                IF (zvardimids(j) == zdimid(i)) THEN
                   zorder(i)=j
                   zorderlog(i)=.true.
                   EXIT
                END IF
             END DO
          END DO
          IF ( ANY (zorderlog .eqv. .false. )) THEN
             CALL finish ('read_sumvar_nf77_2d:','dimensions of variable '//TRIM(zvarname)//' of file '// &
                  TRIM(file_name)//' unknown')
          END IF
          ALLOCATE(zin(zdims(zorder(1)), zdims(zorder(2)), zdims(zorder(3))))
!          ALLOCATE(zinre(zdims(1),zdims(2),zdims(3)))
          zstartvar(1)=zstart(zorder(1))
          zstartvar(2)=zstart(zorder(2))
          zstartvar(3)=zstart(zorder(3))
          zcountvar(1)=zcount(zorder(1))
          zcountvar(2)=zcount(zorder(2))
          zcountvar(3)=zcount(zorder(3))
          CALL nf_check(nf_get_vara_double(zncid, ivar, zstartvar, zcountvar, zin),fname=TRIM(file_name))
!          zinre=RESHAPE(zin, zdims, order=zorder)
          varptr = varptr + RESHAPE(RESHAPE(zin, zdims, order=zorder),zdims(1:2))
          IF (ASSOCIATED(zin)) DEALLOCATE(zin)
       END IF
    END DO
    CALL nf_check(nf_close(zncid),fname=TRIM(file_name))
    END SUBROUTINE read_sumvar_nf77_2d
!===============================================================================
!===============================================================================
  SUBROUTINE read_var_hs_nf77_0d (file_name, dimname, &
                                 idx, varname, varptr, ierr)

    ! Description:
    !
    ! read in from a netCDF file variable varname.
    ! dimname2 has to be the unlimited time dimension of arbitrary length.
    ! The hyperslab dimname2(idx) will be read.
    !
    ! Method:
    !
    ! Open the netCDF file named file_name and read in varptr varname.
    ! The routine is intended to be run on the I/O-processor
    !
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname ! name of the unlimited (time) dimension
    !    can have an arbitrary length, but idx <= length of this dimension constrained
    ! INTEGER, INTENT (in)             :: idx   ! index number of unlimited (time)
    !    dimension2 (<= length of unlimited (time) dimension)
    ! CHARACTER (LEN = *), INTENT (in) :: varname   ! name of variable to be read
    ! REAL(dp), INTENT(out)            :: varptr ! variable value on exit
    ! INTEGER, INTENT(out)             :: ierr ! ierr is 0 if variable varname is successfully
    !    read, 1 otherwise
    !
    ! Authors:
    !
    ! J.S. Rast, MPI, March 2010, original source
    !
    ! Arguments:
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname
    INTEGER, INTENT(in)                   :: idx
    CHARACTER (LEN = *), INTENT (in)      :: varname
    REAL(dp), INTENT(out)                 :: varptr
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, znvar, zdimlength, zvardims
    INTEGER                               :: zdimid
    INTEGER, DIMENSION(1)                 :: zstart, zcount, zvardimid
    REAL(dp), POINTER                     :: zin(:)

    character (len = 256)                 :: zerrstrg

! define dimension length vector
!    zdims(1) = SIZE (varptr, DIM = 1)
!    zdims(2) = 1 ! the dimension length is arbitrary, but we only read one
! month define start vector of dimensions and number of indices in each
! dimension to be read
    zstart(1) = idx
    zcount(1) = 1
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid),fname=TRIM(file_name))
! verify lenght of time dimension
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname), zdimid),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid, zdimlength),fname=TRIM(file_name))
    IF (zdimlength < idx) THEN
       WRITE (zerrstrg,*) idx
       CALL finish('read_var_hs_nf77_0d:', 'actual length of unlimited dimension is less than index ' &
                   //TRIM(zerrstrg))
    END IF
! inquire variable and number of dimension including the unlimited (time) dimension
    ierr=nf_inq_varid(zncid, TRIM(varname), znvar)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_varndims(zncid, znvar, zvardims),fname=TRIM(file_name))
       CALL nf_check(nf_inq_vardimid(zncid, znvar, zvardimid),fname=TRIM(file_name))
       IF (zvardims /= 1) THEN
          CALL finish('read_var_hs_nf77_0d:', &
               'wrong number of dimension of variable '//TRIM(varname))
       END IF
       ALLOCATE(zin(1))
       CALL nf_check(nf_get_vara_double(zncid, znvar, zstart, zcount, zin),fname=TRIM(file_name))
       varptr = zin(1)
       IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    ELSE
       CALL finish('read_var_hs_nf77_1d:', &
               'variable '//TRIM(varname)//' not found in '//TRIM(file_name))
    END IF
    CALL nf_check(nf_close(zncid),fname=TRIM(file_name))
  END SUBROUTINE read_var_hs_nf77_0d
!===============================================================================
!===============================================================================
  SUBROUTINE read_var_hs_nf77_1d (file_name, dimname1, dimname2, &
                                 idx, varname, varptr, ierr)

    ! Description:
    !
    ! read in from a netCDF file variable varname(dimname1).
    ! dimname2 has to be the unlimited time dimension of arbitrary length.
    ! The hyperslab dimname2(idx) will be read.
    !
    ! Method:
    !
    ! Open the netCDF file named file_name and read in varptr varname(dimname1).
    ! The routine is intended to be run on the I/O-processor
    !
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname1 ! name of dimension in the netCDF
    !    file which shall appear as the first dimension of varptr. varname can
    !    have this dimension at any of the four possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname2 ! name of the unlimited (time) dimension
    !    can have an arbitrary length, but idx <= length of this dimension constrained
    ! INTEGER, INTENT (in)             :: idx   ! index number of unlimited (time)
    !    dimension2 (<= length of unlimited (time) dimension)
    ! CHARACTER (LEN = *), INTENT (in) :: varname   ! name of variable to be read
    ! REAL(dp), INTENT(out)                :: varptr(:) ! variable values on exit, has to have
    !    dimension exactly as the lengths of dimname1 in netCDF file
    ! INTEGER, INTENT(out)             :: ierr ! ierr is 0 if variable varname is successfully
    !    read, 1 otherwise
    !
    ! Authors:
    !
    ! J.S. Rast, MPI, August 2004, original source
    !
    ! Arguments:
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname1
    CHARACTER (LEN = *), INTENT (in)      :: dimname2
    INTEGER, INTENT(in)                   :: idx
    CHARACTER (LEN = *), INTENT (in)      :: varname
    REAL(dp), INTENT(out)                 :: varptr(:)
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, znvar, zdimlength, zvardims
    INTEGER, DIMENSION(2)                 :: zdims, zdimid, zvardimids, &
                                             zstart, zcount
    REAL(dp), POINTER                     :: zin(:,:)

    character (len = 256)                 :: zerrstrg

! define dimension length vector
    zdims(1) = SIZE (varptr, DIM = 1)
    zdims(2) = 1 ! the dimension length is arbitrary, but we only read one
! month define start vector of dimensions and number of indices in each
! dimension to be read
    zstart(1) = 1
    zstart(2) = idx
    zcount(1) = zdims(1)
    zcount(2) = 1
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid),fname=TRIM(file_name))
! verify lenghts of dimensions
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname1), zdimid(1)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(1), zdimlength),fname=TRIM(file_name))
    IF (zdimlength /= zdims(1)) THEN
       CALL finish('read_var_hs_nf77_1d:', 'wrong length of dim1 in file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname2), zdimid(2)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(2), zdimlength),fname=TRIM(file_name))
    IF (zdimlength < idx) THEN
       WRITE (zerrstrg,*) idx
       CALL finish('read_var_hs_nf77_1d:', 'actual length of unlimited dimension is less than index ' &
                   //TRIM(zerrstrg))
    END IF
! inquire variable and number of dimension including the unlimited (time) dimension
    ierr=nf_inq_varid(zncid, TRIM(varname), znvar)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_varndims(zncid, znvar, zvardims),fname=TRIM(file_name))
       CALL nf_check(nf_inq_vardimid(zncid, znvar, zvardimids),fname=TRIM(file_name))
       IF (zvardims /= 2) THEN
          CALL finish('read_var_hs_nf77_1d:', &
               'wrong number of dimension of variable '//TRIM(varname))
       END IF
       ALLOCATE(zin(zdims(1), zdims(2)))
       CALL nf_check(nf_get_vara_double(zncid, znvar, zstart, zcount, zin),fname=TRIM(file_name))
       varptr = zin(:,1)
       IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    ELSE
       CALL finish('read_var_hs_nf77_1d:', &
               'variable '//TRIM(varname)//' not found in '//TRIM(file_name))
    END IF
    CALL nf_check(nf_close(zncid),fname=TRIM(file_name))
    END SUBROUTINE read_var_hs_nf77_1d
!===============================================================================
!===============================================================================
  SUBROUTINE read_var_hs_nf77_2d (file_name, dimname1, dimname2, dimname3, idx, varname, varptr, ierr)

    ! Description:
    !
    ! read in from a netCDF file variable varname(dimname1, dimname2).
    ! dimname3 has to be the unlimited time dimension of arbitrary length.
    ! The hyperslab dimname3(idx) will be read.
    !
    ! Method:
    ! 
    ! Open the netCDF file named file_name and read in varptr varname(dimname1, dimname2).
    ! The order of the dimensions in the netCDF file is irrelevant.
    ! The routine is to be intended to run on the I/O-processor
    ! 
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname1 ! name of dimension in the netCDF
    !    file which shall appear as the first dimension of varptr. varname can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname2 ! name of dimension in the netCDF
    !    file which shall appear as the second dimension of varptr. varname can 
    !    have this dimension at any of the three possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname3 ! name of the unlimited (time) dimension
    !    can have an arbitrary length, but idx <= length of this dimension constrained
    ! INTEGER, INTENT (in)             :: idx   ! index number of unlimited (time)
    !    dimension3
    ! CHARACTER (LEN = *), INTENT (in) :: varname   ! name of variable to be read
    ! REAL(dp), INTENT(out)            :: varptr(:,:) ! variable values on exit, has to
    !    dimension exactly as the lengths of dimname1,dimname2 in netCDF file
    ! INTEGER, INTENT(out)             :: ierr ! ierr is 0 if variable varname is successfully 
    !    read, 1 otherwise
    !
    ! Authors:
    ! 
    ! J.S. Rast, MPI, January 2003, original source
    !
    ! Arguments:
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname1
    CHARACTER (LEN = *), INTENT (in)      :: dimname2
    CHARACTER (LEN = *), INTENT (in)      :: dimname3
    INTEGER, INTENT(in)                   :: idx
    CHARACTER (LEN = *), INTENT (in)      :: varname 
    REAL(dp), INTENT(out)                 :: varptr(:,:)
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, znvar, zdimlength, zvardims, i, j
    INTEGER, DIMENSION(3)                 :: zdims, zdimid, zvardimids, zorder, zstart, zcount, &
                                             zstartvar, zcountvar
    REAL(dp), POINTER                     :: zin(:,:,:)
!,zinre(:,:,:)
    character (len = 256)                 :: zerrstrg

! define dimension length vector
    zdims(1) = SIZE (varptr, DIM = 1)
    zdims(2) = SIZE (varptr, DIM = 2)
    zdims(3) = 1 ! the dimension length normally is arbitrary, but we only read one month
! define start vector of dimensions and number of indices in each dimension to be read
    zstart(1) = 1
    zstart(2) = 1
    zstart(3) = idx
    zcount(1) = zdims(1)
    zcount(2) = zdims(2)
    zcount(3) = 1
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid), fname=TRIM(file_name))
! verify lenghts of dimensions
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname1), zdimid(1)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(1), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(1)) THEN
       CALL finish('read_var_hs_nf77_2d:', 'wrong length of dim1 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname2), zdimid(2)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(2), zdimlength), fname=TRIM(file_name))
    IF (zdimlength /= zdims(2)) THEN
       CALL finish('read_var_hs_nf77_2d:', 'wrong length of dim2 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname3), zdimid(3)), fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(3), zdimlength), fname=TRIM(file_name))
    IF (zdimlength < idx) THEN
       WRITE (zerrstrg,*) idx
       CALL finish('read_var_hs_nf77_2d:', 'actual length of unlimited dimension is less than index ' &
                   //TRIM(zerrstrg))
    END IF
! inquire variable and number of dimension including the unlimited (time) dimension
    ierr=nf_inq_varid(zncid, TRIM(varname), znvar)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_varndims(zncid, znvar, zvardims), fname=TRIM(file_name))
       CALL nf_check(nf_inq_vardimid(zncid, znvar, zvardimids),fname=TRIM(file_name))
       IF (zvardims /= 3) THEN
          CALL finish('read_var_hs_nf77_2d:', &
               'wrong number of dimension of variable '//TRIM(varname))
       END IF
       DO i = 1, zvardims
          DO j = 1, zvardims
             IF (zvardimids(j) == zdimid(i)) THEN
                zorder(i)=j
                EXIT
             END IF
          END DO
       END DO
       ALLOCATE(zin(zdims(zorder(1)), zdims(zorder(2)), zdims(zorder(3))))
       zstartvar(1)=zstart(zorder(1))
       zstartvar(2)=zstart(zorder(2))
       zstartvar(3)=zstart(zorder(3))
       zcountvar(1)=zcount(zorder(1))
       zcountvar(2)=zcount(zorder(2))
       zcountvar(3)=zcount(zorder(3))
       CALL nf_check(nf_get_vara_double(zncid, znvar, zstartvar, zcountvar, zin), fname=TRIM(file_name))
       varptr = RESHAPE(RESHAPE(zin, zdims, order=zorder),zdims(1:2))
       IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    ELSE
       CALL finish('read_var_hs_nf77_2d:', &
               'variable '//TRIM(varname)//' not found in '//TRIM(file_name))
    END IF
    CALL nf_check(nf_close(zncid),fname=TRIM(file_name))
    END SUBROUTINE read_var_hs_nf77_2d
!===============================================================================
!===============================================================================
  SUBROUTINE read_var_hs_nf77_3d (file_name, dimname1, dimname2, dimname3, dimname4, &
                                 idx, varname, varptr, ierr)

    ! Description:
    !
    ! read in from a netCDF file variable varname(dimname1, dimname2, dimname3).
    ! dimname4 has to be the unlimited time dimension of arbitrary length.
    ! The hyperslab dimname4(idx) will be read.
    !
    ! Method:
    ! 
    ! Open the netCDF file named file_name and read in varptr varname(dimname1, dimname2, dimname3).
    ! The order of the dimensions in the netCDF file is irrelevant.
    ! The routine is to be intended to run on the I/O-processor
    ! 
    ! Parameter list:
    !
    ! CHARACTER (LEN = *), INTENT (in) :: file_name ! file name of the netCDF file
    !    from which variable should be read
    ! CHARACTER (LEN = *), INTENT (in) :: dimname1 ! name of dimension in the netCDF
    !    file which shall appear as the first dimension of varptr. varname can 
    !    have this dimension at any of the four possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname2 ! name of dimension in the netCDF
    !    file which shall appear as the second dimension of varptr. varname can 
    !    have this dimension at any of the four possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname3 ! name of dimension in the netCDF
    !    file which shall appear as the third dimension of varptr. varname can
    !    have this dimension at any of the four possible positions
    ! CHARACTER (LEN = *), INTENT (in) :: dimname4 ! name of the unlimited (time) dimension
    !    can have an arbitrary length, but idx <= length of this dimension constrained
    ! INTEGER, INTENT (in)             :: idx   ! index number of unlimited (time)
    !    dimension4 (<= length of unlimited (time) dimension)
    ! CHARACTER (LEN = *), INTENT (in) :: varname   ! name of variable to be read
    ! REAL(dp), INTENT(out)            :: varptr(:,:,:) ! variable values on exit, has to
    !    dimension exactly as the lengths of dimname1,dimname2,dimname3 in netCDF file
    ! INTEGER, INTENT(out)             :: ierr ! ierr is 0 if variable varname is successfully 
    !    read, 1 otherwise
    !
    ! Authors:
    ! 
    ! J.S. Rast, MPI, January 2003, original source
    !
    ! Arguments:
    !
    CHARACTER (LEN = *), INTENT (in)      :: file_name
    CHARACTER (LEN = *), INTENT (in)      :: dimname1
    CHARACTER (LEN = *), INTENT (in)      :: dimname2
    CHARACTER (LEN = *), INTENT (in)      :: dimname3
    CHARACTER (LEN = *), INTENT (in)      :: dimname4
    INTEGER, INTENT(in)                   :: idx
    CHARACTER (LEN = *), INTENT (in)      :: varname 
    REAL(dp), INTENT(out)                 :: varptr(:,:,:)
    INTEGER, INTENT (out)                 :: ierr
    !
    ! Local variables:
    !
    INTEGER                               :: zncid, znvar, zdimlength, zvardims, i, j
    INTEGER, DIMENSION(4)                 :: zdims, zdimid, zvardimids, zorder, zstart, zcount, &
                                             zstartvar, zcountvar
    REAL(dp), POINTER                     :: zin(:,:,:,:)


    character (len = 256)                 :: zerrstrg

! define dimension length vector
    zdims(1) = SIZE (varptr, DIM = 1)
    zdims(2) = SIZE (varptr, DIM = 2)
    zdims(3) = SIZE (varptr, DIM = 3)
    zdims(4) = 1 ! the dimension length is arbitrary, but we only read one 
! month define start vector of dimensions and number of indices in each 
! dimension to be read
    zstart(1) = 1
    zstart(2) = 1
    zstart(3) = 1

    zstart(4) = idx
    zcount(1) = zdims(1)
    zcount(2) = zdims(2)
    zcount(3) = zdims(3)
    zcount(4) = 1
! open netCDF file
    CALL nf_check(nf__open(TRIM(file_name), nf_nowrite, chunksize, zncid),fname=TRIM(file_name))
! verify lenghts of dimensions
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname1), zdimid(1)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(1), zdimlength),fname=TRIM(file_name))
    IF (zdimlength /= zdims(1)) THEN
       CALL finish('read_var_hs_nf77_3d:', 'wrong length of dim1 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname2), zdimid(2)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(2), zdimlength),fname=TRIM(file_name))
    IF (zdimlength /= zdims(2)) THEN
       CALL finish('read_var_hs_nf77_3d:', 'wrong length of dim2 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname3), zdimid(3)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(3), zdimlength),fname=TRIM(file_name))
    IF (zdimlength /= zdims(3)) THEN
       CALL finish('read_var_hs_nf77_3d:', 'wrong length of dim3 in initial file '//TRIM(file_name))
    END IF
    CALL nf_check(nf_inq_dimid(zncid, TRIM(dimname4), zdimid(4)),fname=TRIM(file_name))
    CALL nf_check(nf_inq_dimlen(zncid, zdimid(4), zdimlength),fname=TRIM(file_name))
    IF (zdimlength < idx) THEN
       WRITE (zerrstrg,*) idx
       CALL finish('read_var_hs_nf77_3d:', 'actual length of unlimited dimension is less than index ' &
                   //TRIM(zerrstrg))
    END IF
! inquire variable and number of dimension including the unlimited (time) dimension
    ierr=nf_inq_varid(zncid, TRIM(varname), znvar)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_varndims(zncid, znvar, zvardims),fname=TRIM(file_name))
       CALL nf_check(nf_inq_vardimid(zncid, znvar, zvardimids),fname=TRIM(file_name))
       IF (zvardims /= 4) THEN
          CALL finish('read_var_hs_nf77_3d:', &
               'wrong number of dimension of variable '//TRIM(varname))
       END IF
       DO i = 1, zvardims
          DO j = 1, zvardims
             IF (zvardimids(j) == zdimid(i)) THEN
                zorder(i)=j
                EXIT
             END IF
          END DO
       END DO
       ALLOCATE(zin(zdims(zorder(1)), zdims(zorder(2)), zdims(zorder(3)), zdims(zorder(4))))
       zstartvar(1)=zstart(zorder(1))
       zstartvar(2)=zstart(zorder(2))
       zstartvar(3)=zstart(zorder(3))
       zstartvar(4)=zstart(zorder(4))
       zcountvar(1)=zcount(zorder(1))
       zcountvar(2)=zcount(zorder(2))
       zcountvar(3)=zcount(zorder(3))
       zcountvar(4)=zcount(zorder(4))
       CALL nf_check(nf_get_vara_double(zncid, znvar, zstartvar, zcountvar, zin),fname=TRIM(file_name))
       varptr = RESHAPE(RESHAPE(zin, zdims, order=zorder),zdims(1:3))
       IF (ASSOCIATED(zin)) DEALLOCATE(zin)
    ELSE
       CALL finish('read_var_hs_nf77_3d:', &
               'variable '//TRIM(varname)//' not found in '//TRIM(file_name))
    END IF
    CALL nf_check(nf_close(zncid),fname=TRIM(file_name))
  END SUBROUTINE read_var_hs_nf77_3d
!===============================================================================
!=============================================================================================
  LOGICAL FUNCTION search_var_nf77(file_name, varname)

    ! Description:
    !
    ! search for the variable varname in file file_name
    !
    ! LOGICAL             :: search_var_nf77 ! if varname is in file_name, then
    !    search_var_nf77=.true., else it is .false.
    ! CHARACTER (LEN = *) :: file_name ! file name of the netCDF file
    !    in which the variable varname should be searched
    ! CHARACTER (LEN = *) :: varname ! name of variable you are searching
    !    for
    !
    ! Authors:
    ! 
    ! J.S. Rast, MPI, July 2006, original source
    !

    CHARACTER (LEN = *) :: file_name
    CHARACTER (LEN = *) :: varname

    ! local variables
    INTEGER             :: zncid, znvar, ierr

    search_var_nf77=.false.
! open netCDF file
    CALL nf_check(nf_open(TRIM(file_name), nf_nowrite, zncid),fname=TRIM(file_name))
    ierr=nf_inq_varid(zncid, TRIM(varname), znvar)
    IF (ierr == 0) THEN
       search_var_nf77=.true.
    END IF
  END FUNCTION search_var_nf77
!===============================================================================
!===============================================================================
  INTEGER FUNCTION read_diml_nf77(file_name, dimname)

    ! Description:
    !
    ! gives length of dimension dimname in file file_name
    ! 
    ! INTEGER              :: read_diml_nf77  ! length of dimension dimname,
    !    if dimname is not present in file_name, read_diml_nf77=-1
    ! CHARACTER (LEN = *)  :: file_name       ! file name of the netCDF file 
    !    containing dimension dimname
    ! CHARACTER (LEN = *)  :: dimname         ! name of dimension the length
    !    of which is searched
    !
    ! Authors:
    !
    ! J.S. Rast, MPI, January 2010, original source

    CHARACTER (LEN = *)  :: file_name
    CHARACTER (LEN = *)  :: dimname

    ! local variables
    INTEGER              :: zncid, zdimid, ierr
    read_diml_nf77=-1
    CALL nf_check(nf_open(TRIM(file_name), nf_nowrite, zncid),fname=TRIM(file_name))
    ierr=nf_inq_dimid(zncid, TRIM(dimname), zdimid)
    IF (ierr == 0) THEN
       CALL nf_check(nf_inq_dimlen(zncid, zdimid, read_diml_nf77),fname=TRIM(file_name))
    END IF
    CALL nf_check(nf_close(zncid))
    END FUNCTION read_diml_nf77
!===============================================================================
  END MODULE mo_read_netcdf77
