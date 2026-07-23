!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_utils

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


CONTAINS 

  FUNCTION get_month(date) RESULT(imon)

    USE mo_time_conversion, ONLY : time_days, time_native, TC_get, TC_convert
  
    TYPE(time_days)   :: date
    INTEGER           :: imon

    TYPE(time_native) :: date_nat
    INTEGER           :: yr, mo, dy, hr, mn, se

    CALL TC_convert(date,date_nat)
    CALL TC_get(date_nat,yr,mo,dy,hr,mn,se)

    ! Time 00:00 belongs to the previous day, therefore 00Z of the first of
    ! each month belongs to the previous month
    ! CHECK THIS
    IF (dy == 1 .AND. hr == 0) THEN
       mo = mo-1
       IF (mo < 1) mo = 12
    ENDIF
    imon = mo

  END FUNCTION get_month

  FUNCTION get_current_month() RESULT(imon)

    USE mo_time_control, Only: current_date

    INTEGER :: imon

    imon = get_month(current_date)

  END FUNCTION get_current_month

  SUBROUTINE average_tiles(in, mask, fract, out)

    REAL(dp),    INTENT(in)  :: in(:,:)
    LOGICAL, INTENT(in)  :: mask(:,:)
    REAL(dp),    INTENT(in)  :: fract(:,:)
    REAL(dp),    INTENT(out) :: out(:)

    REAL(dp) :: f(SIZE(fract,1),SIZE(fract,2))
    REAL(dp) :: s(SIZE(in,1))

    IF (SIZE(in,2) == 1) THEN
       out(:) = in(:,1)
       RETURN
    END IF

    f = MERGE(fract,0.0_dp,mask)
    s = SUM(f, DIM=2)
    out = SUM(in * f, DIM=2)

    WHERE (s(:) > 0._dp)
       out(:) = out(:) / s(:)
    END WHERE

  END SUBROUTINE average_tiles

END MODULE mo_utils

!- End of module header
