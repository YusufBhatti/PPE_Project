!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
! (c) COPYRIGHT British Crown / Met Office 2008  
! Please refer to Met_Office_license.txt for details.
!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Jul 2008 - A. Bodas-Salcedo - Added definitions of ISCCP axes
! Oct 2008 - H. Chepfer       - Added PARASOL_NREFL

MODULE mo_cosp_constants

    USE mo_kind,  ONLY: dp

    IMPLICIT NONE
    
    ! Indices to address arrays of LS and CONV hydrometeors
    INTEGER, PARAMETER :: I_LSCLIQ = 1
    INTEGER, PARAMETER :: I_LSCICE = 2
    
    ! Missing value
    REAL(dp), PARAMETER :: R_UNDEF =   -1.0E30_dp
    ! Number of possible output variables
    INTEGER, PARAMETER :: N_OUT_LIST = 27
    ! Value for forward model result from a level that is under the ground
    REAL(dp), PARAMETER :: R_GROUND = -1.0E20_dp

!--- Lidar constants
    ! CFAD constants
    INTEGER, PARAMETER :: SR_BINS       =   15
    INTEGER, PARAMETER :: DPOL_BINS     =   6
    REAL(dp), PARAMETER    :: LIDAR_UNDEF   =   999.999_dp
    ! Other constants
    INTEGER, PARAMETER :: LIDAR_NCAT    =   4
    INTEGER, PARAMETER :: PARASOL_NREFL =   5 ! parasol
    REAL(dp), PARAMETER, DIMENSION(PARASOL_NREFL) :: PARASOL_SZA = (/0.0_dp, 20.0_dp, 40.0_dp, 60.0_dp, 80.0_dp/)
    REAL(dp), PARAMETER :: DEFAULT_LIDAR_REFF = 30.0e-6_dp ! Default lidar effective radius
    
    !  number of hydrometeor classes
    INTEGER, PARAMETER :: N_HYDRO = 2

    !levels for vgrid
    INTEGER :: Nlr = 40

    LOGICAL :: use_vgrid = .TRUE. ! DON'T CHANGE output on 40 level satellite grid?

END MODULE mo_cosp_constants
