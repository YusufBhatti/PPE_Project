!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_rad_forcing_diag
  !
  ! Authors:
  ! -------
  ! C. Timmreck, M.A. Thomas, S. Lorenz, S. Rast, M. Schultz, MPI-MET 2001-2009

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream, GRIB, HYBRID_H
  USE mo_memory_base,   ONLY: new_stream, add_stream_element,       &
                              default_stream_setting, add_stream_reference, AUTO
  USE mo_netCDF,        ONLY: max_dim_name
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_control,       ONLY: nlon, nlev, ngl  ! size of global arrays

  IMPLICIT NONE

  PRIVATE

  !--- Service routines: ----------------------------------------------------------------

  PUBLIC :: init_stream_rad_forcing            ! construct the radf_stream stream

  !--- Declarations for stream radf_stream: ----------------------------------------------------

  TYPE (t_stream), PUBLIC, POINTER :: radf_stream
  
  !
  !  variables for double radiation
  !
  REAL(dp), POINTER, PUBLIC :: emtef01(:,:,:)
  REAL(dp), POINTER, PUBLIC :: emtef1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: emter1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: netht_lw(:,:,:)
  REAL(dp), POINTER, PUBLIC :: netht_sw(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof01(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsol1(:,:,:)
  !
  REAL(dp), POINTER, PUBLIC :: emtef01m(:,:,:)
  REAL(dp), POINTER, PUBLIC :: emtef1m(:,:,:)
  REAL(dp), POINTER, PUBLIC :: emter1m(:,:,:)
  REAL(dp), POINTER, PUBLIC :: netht_lwm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: netht_swm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof01m(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof1m(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsol1m(:,:,:)
  !
  REAL(dp), PUBLIC, POINTER :: d_aflx_sw(:,:,:)
  REAL(dp), PUBLIC, POINTER :: d_aflx_lw(:,:,:)
  REAL(dp), PUBLIC, POINTER :: d_aflx_swc(:,:,:)
  REAL(dp), PUBLIC, POINTER :: d_aflx_lwc(:,:,:)
  
  !--- Declare dimensions for the streams: ----------------------------------------------------------------------

  INTEGER :: lnlon,   lnlev,  lngl ! size of local arrays on PE
  INTEGER :: nlevp1
  INTEGER :: lnlevp1
  INTEGER :: dim3(3), dim3p(3)
  CHARACTER (max_dim_name) :: dim3n(3)
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE init_stream_rad_forcing

    !--- Set dimensions to be used in the declarations:

    lnlon=lc%nproma
    lnlev=lc%nlev
    lngl =lc%ngpblks

    nlevp1  = nlev  + 1
    lnlevp1 = lnlev + 1

    dim3p = (/ lnlon,  lnlevp1, lngl  /)
    dim3  = (/  nlon,   nlevp1,  ngl  /)
    dim3n = (/  "lon ","ilev","lat "/)
    
    
    !--- 1) Construct the radf stream: ------------------------------------------------------------------------------

    CALL new_stream (radf_stream,'radf',filetype=GRIB,lrerun=.true.)

    !--- Add standard fields for post-processing:

    CALL add_stream_reference (radf_stream, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (radf_stream, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (radf_stream, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (radf_stream, 'gboxarea','geoloc',lpost=.TRUE.)

    !--- 2) Add stream elements: ------------------------------------------------------------------------------------

    CALL default_stream_setting (radf_stream, units = 'Wm**-2',    &
                                      lrerun   = .TRUE. ,     &
                                      contnorest=.TRUE. ,    &
                                      laccu    = .TRUE. ,     &
                                      lpost    = .TRUE. ,     &
                                      gdims    = dim3 ,       &
                                      dimnames = dim3n ,      &
                                      leveltype = HYBRID_H,   &
                                      table    = 199,         &
                                      code     = AUTO         )

    !--- 2.1) flux anomalies, perturbed minus unperturbed
    
    CALL add_stream_element (radf_stream, 'd_aflx_sw', d_aflx_sw, dim3p, code=73,  &
                             longname='Accumulated SW flux anomalies - all sky ')
    CALL add_stream_element (radf_stream, 'd_aflx_lw', d_aflx_lw, dim3p, code=74,  &
                             longname='Accumulated LW flux anomalies - all sky')
    CALL add_stream_element (radf_stream, 'd_aflx_swc', d_aflx_swc, dim3p, code=75,  &
                             longname='Accumulated SW flux anomalies - clear')
    CALL add_stream_element (radf_stream, 'd_aflx_lwc', d_aflx_lwc, dim3p, code=76,  &
                             longname='Accumulated LW flux anomalies - clear')
    !
    !--- 2.2) heating rates
    !
    CALL add_stream_element (radf_stream, 'netht_sw', netht_sw, code=77,  &
                             longname='Net SW heating rates', units='K/s')
    CALL add_stream_element (radf_stream, 'netht_lw', netht_lw, code=78,  &
                             longname='Net LW heating rates', units='K/s')

    !
    !--- 2.3) flux variables of second call of radiation not written per default
    !
    CALL add_stream_element (radf_stream,'emtef01', emtef01, lpost=.FALSE.)
    CALL add_stream_element (radf_stream,'emtef1',  emtef1,  lpost=.FALSE.)
    CALL add_stream_element (radf_stream,'emter1',  emter1,  lpost=.FALSE.)
    CALL add_stream_element (radf_stream,'trsof01', trsof01, lpost=.FALSE.)
    CALL add_stream_element (radf_stream,'trsof1',  trsof1,  lpost=.FALSE.)
    CALL add_stream_element (radf_stream,'trsol1',  trsol1,  lpost=.FALSE.)
    

    netht_lwm  => netht_lw
    netht_swm  => netht_sw

    emtef01m   => emtef01
    emtef1m    => emtef1
    emter1m    => emter1
    trsof01m   => trsof01
    trsof1m    => trsof1
    trsol1m    => trsol1

  END SUBROUTINE init_stream_rad_forcing
  

!---------------------------------------------------------------------------------------
END MODULE mo_rad_forcing_diag
