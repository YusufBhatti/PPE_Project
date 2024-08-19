!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_jsbach_time
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

  IMPLICIT NONE 


CONTAINS 
  ! Define procedures contained in this module. 

  SUBROUTINE jsbach_start_timestep

    ! Do initialization for time step, i.e. this is the place to put computations that
    ! have to be done once at the start of a new time step

    ! Called either from inside the time loop of the jsbach driver or, in echam coupled
    ! mode, either from *stepon* before the call to *scan1* or from the beginning of *scan1*

  END SUBROUTINE jsbach_start_timestep


END MODULE mo_jsbach_time

!- End of module header
