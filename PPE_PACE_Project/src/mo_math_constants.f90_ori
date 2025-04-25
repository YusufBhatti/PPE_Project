!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!!  Module defines main mathematical constants
!!
!!  These constants are provided to more significant digits than is necessary for a 64-bit
!!  double precision number; they may be used for other purposes where the extra precision 
!!  is necessary or useful.
!! 
!!  The variables are calculated by a multi-precision floating point package 
!!  (mpmath in python). 
!!
!! @par Revision History
!!  Developed  by Luis Kornblueh (2004)
!!  Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!!  Modified according to style guide by Thomas Heinze (2005-06-24):
!!   - module renamed from mo_math to mo_math_constants
!!   - eps moved from mo_physical_constants
!!   - pid180i renamed to rad2deg
!!  Including some more constants from math.h by Thomas Heinze (2005-07-18)
!!  pid5 renamed to pi_5 by Thomas Heinze (2005-07-26)
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Thomas Heinze (2006-05-18):
!!  - introduced dbl_eps
!!  Modified by Luis Kornblueh (2012)
!!  - replaced arbitrary number of digits by 42 significant decimal digits calculated by 
!!    multi-precision algorithms
!!
MODULE mo_math_constants
  
  USE mo_kind,      ONLY:  wp

  IMPLICIT NONE
  
  PUBLIC

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  ! Mathematical constants defined:
  !
  !--------------------------------------------------------------
  ! Fortran name | C name       | meaning                       |
  !--------------------------------------------------------------
  ! euler        | M_E          |  e                            | 
  ! log2e        | M_LOG2E      |  log2(e)                      | 
  ! log10e       | M_LOG10E     |  log10(e)                     | 
  ! ln2          | M_LN2        |  ln(2)                        | 
  ! ln10         | M_LN10       |  ln(10)                       | 
  ! pi           | M_PI         |  pi                           | 
  ! pi_2         | M_PI_2       |  pi/2                         | 
  ! pi_4         | M_PI_4       |  pi/4                         | 
  ! rpi          | M_1_PI       |  1/pi                         | 
  ! rpi_2        | M_2_PI       |  2/pi                         | 
  ! rsqrtpi_2    | M_2_SQRTPI   |  2/(sqrt(pi))                 | 
  ! sqrt2        | M_SQRT2      |  sqrt(2)                      | 
  ! sqrt1_2      | M_SQRT1_2    |  1/sqrt(2)                    | 
  ! sqrt3        |              |  sqrt(3)                      | 
  ! sqrt1_3      |              |  1/sqrt(3)                    |
  ! half angle of pentagon                                      | 
  ! pi_5         |              |  pi/5                         |
  ! latitude of the lowest major  triangle corner               |
  ! and latitude of the major hexagonal faces centers           |
  ! phi0         |              |  pi/2 -2acos(1/(2*sin(pi/5))) |
  ! conversion factor from radians to degree                    |
  ! rad2deg      |              |  180/pi                       |
  ! conversion factor from degree to radians                    | 
  ! deg2rad      |              |  pi/180                       |
  ! one_third    |              |  1/3                          |
  !-------------------------------------------------------------|
  
  REAL(wp), PARAMETER :: euler     =  2.71828182845904523536028747135266249775725_wp
  REAL(wp), PARAMETER :: log2e     =  1.44269504088896340735992468100189213742665_wp
  REAL(wp), PARAMETER :: log10e    =  0.434294481903251827651128918916605082294397_wp
  REAL(wp), PARAMETER :: ln2       =  0.6931471805599453094172321214581765680755_wp
  REAL(wp), PARAMETER :: ln10      =  2.3025850929940456840179914546843642076011_wp
  REAL(wp), PARAMETER :: pi        =  3.14159265358979323846264338327950288419717_wp
  REAL(wp), PARAMETER :: pi_2      =  1.57079632679489661923132169163975144209858_wp
  REAL(wp), PARAMETER :: pi_4      =  0.785398163397448309615660845819875721049292_wp
  REAL(wp), PARAMETER :: rpi       =  0.318309886183790671537767526745028724068919_wp
  REAL(wp), PARAMETER :: rpi_2     =  0.636619772367581343075535053490057448137839_wp
  REAL(wp), PARAMETER :: rsqrtpi_2 =  1.1283791670955125738961589031215451716881_wp
  REAL(wp), PARAMETER :: sqrt2     =  1.41421356237309504880168872420969807856967_wp
  REAL(wp), PARAMETER :: sqrt1_2   =  0.707106781186547524400844362104849039284836_wp
  REAL(wp), PARAMETER :: sqrt3     =  1.73205080756887729352744634150587236694281_wp
  REAL(wp), PARAMETER :: sqrt1_3   =  0.577350269189625764509148780501957455647602_wp
  REAL(wp), PARAMETER :: cos45     =  0.707106781186547524400844362104849039284836_wp
  REAL(wp), PARAMETER :: pi_5      =  0.628318530717958647692528676655900576839434_wp
  REAL(wp), PARAMETER :: pi2       =  6.28318530717958647692528676655900576839434_wp
  REAL(wp), PARAMETER :: phi0      =  0.463647609000806116214256231461214402028537_wp
  REAL(wp), PARAMETER :: rad2deg   =  57.2957795130823208767981548141051703324055_wp
  REAL(wp), PARAMETER :: deg2rad   =  0.0174532925199432957692369076848861271344287_wp

  REAL(wp), PARAMETER :: one_third =  1.0_wp/3.0_wp

END MODULE mo_math_constants
