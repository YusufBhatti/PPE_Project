!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_jsbach_comm_to_echam5mods
  ! 
  ! Description: 
  !   <Say what this module is for> 
  ! 
  ! Current Code Owner: <Name of person responsible for this code> 
  ! 
  ! History: 
  !  
  ! Version   Date     Comment 
  ! -------   ----     ------- 
  ! <version> <date>   Original code. <Your name> 
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 

  USE mo_kind, ONLY: dp

  IMPLICIT NONE 

  INTEGER           :: nlon, nlat, nland
  INTEGER,  POINTER :: kpoints(:)
  REAL(dp), POINTER :: lon(:), lat(:)
  LOGICAL,  POINTER :: mask(:,:)

  INTEGER           :: domain_nlon, domain_nlat, domain_nland
  LOGICAL,  POINTER :: domain_mask(:,:)

  INTEGER           :: jsb_ntiles, jsb_nsoil, jsb_ntsoil

END MODULE mo_jsbach_comm_to_echam5mods

!- End of module header
