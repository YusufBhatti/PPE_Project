!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE baresoildiff(kidia, kfdia, klon, zdiagt,  &
    ws1, wsbs, cveg, cdel1, fksat, vpor, bclapp, fmpot)
!
!   *** Diffusion of water between vegetated and bare soil part of gridbox
!   *** Numerical recipes/Richtmyer & -Morton diffusion  *************
!
!   Version 1.0 - July 2012, Stefan Hagemann 
!
!   IAVG: Averaging of diffusivities between bare soil and vegetated part.
!     1  zdiff(jl) = zdum2 * (1._dp-cveg(jl)) + zdum3 * cveg(jl) 
!     2  zdiff(jl) = (zdum2 + zdum3) * 0.5 
!     3  zdiff(jl) = (zdum2 + zdum3) * 0.5 * MIN(1._dp-cveg(jl), cveg(jl))
!        The smaller part determines the diffusive exchange 
!        --> reduce diffusivity 
!     4  zdiff(jl) = (zdum2 + zdum3) * 0.5 * 0.05
!        The diffusive exchange takes only part in a small area of the grid box,
!        e.g. 5%. In principle this should depend on the mixing ratio of vegetation
!        and bare soil within a gridbox.
!     5  zdiff(jl) = (zdum2 + zdum3) * 0.5 * 0.01 as 4, just with 1%

USE mo_exception, ONLY: finish
USE mo_kind,      ONLY: dp
IMPLICIT NONE

INTEGER,  INTENT(IN)    :: kidia
INTEGER,  INTENT(IN)    :: kfdia
INTEGER,  INTENT(IN)    :: klon
REAL(dp), INTENT(IN)    :: zdiagt           ! [ZDIAGT] = s
REAL(dp), INTENT(IN)    :: ws1(klon)
REAL(dp), INTENT(INOUT) :: wsbs(klon)
REAL(dp), INTENT(IN)    :: cveg(klon)
REAL(dp), INTENT(IN)    :: cdel1
REAL(dp), INTENT(IN)    :: fksat(klon)
REAL(dp), INTENT(IN)    :: vpor(klon)
REAL(dp), INTENT(IN)    :: bclapp(klon)
REAL(dp), INTENT(IN)    :: fmpot(klon)

! Parameters
INTEGER,  PARAMETER :: iavg = 5
REAL(dp), PARAMETER :: zeps = 1.e-15_dp

! Diffusion
REAL(dp) :: zda (klon), zdb(klon), zdiff(klon)
REAL(dp) :: ztri(klon), zdc(klon), wsveg(klon)

! Other locals 
INTEGER  :: jl
REAL(dp) :: zdum, zdum2, zdum3, zdum4

  DO jl=kidia, kfdia
    IF (cdel1 > 0 .AND. cveg(jl) > 1.e-10_dp) THEN
        
!     Soil moisture diffusivity [m^2/day]
!     Diffusivity of bare soil part and vegetated part is averaged to calculate diffusivity 
      IF (wsbs(jl) > 0) THEN
        zdum = wsbs(jl)/cdel1
          
!       Calculating the diffusivity of bare soil part
        zdum4 = zdum / vpor(jl)
        IF (zdum4-1.e-10_dp > zeps) THEN
          zdum2 = bclapp(jl) * fksat(jl) * fmpot(jl) /  &
                  zdum * zdum4 ** (bclapp(jl)+3.)
        ELSE
          zdum2 = 0._dp
        END IF
      ELSE
        zdum2=0._dp
      END IF

      IF (ws1(jl) > 0._dp) THEN

        wsveg(jl) = MAX( (ws1(jl) - wsbs(jl) *(1._dp-cveg(jl))) / cveg(jl), 0._dp)
        zdum = wsveg(jl) / cdel1
          
!       Calculate the diffusivity of vegetated part
        zdum4 = zdum / vpor(jl)
        IF (zdum4-1.e-10_dp > zeps) THEN
          zdum3 = bclapp(jl) * fksat(jl) * fmpot(jl) /  &
                  zdum *  zdum4 ** (bclapp(jl)+3._dp)
        ELSE
          zdum3 = 0._dp
        END IF
      ELSE
        wsveg(jl)  =  0._dp
        zdum3=0._dp
      END IF
        
!     Calculate the diffusivity weighted by the fractions
      SELECT CASE (iavg)
        CASE (1)
          zdiff(jl) = zdum2 * (1._dp-cveg(jl)) + zdum3 * cveg(jl) 
        CASE (2)
          zdiff(jl) = (zdum2 + zdum3) * 0.5_dp 
        CASE (3)
          zdiff(jl) = (zdum2 + zdum3) * 0.5_dp * MIN(1._dp-cveg(jl), cveg(jl))
        CASE (4)
          zdiff(jl) = (zdum2 + zdum3) * 0.5_dp * 0.05_dp
        CASE (5)
          zdiff(jl) = (zdum2 + zdum3) * 0.5_dp * 0.01_dp
        CASE DEFAULT 
          CALL finish('baresoildiff','iavg < 1 or > 5 ==> STOP')
      END SELECT
    ELSE IF (cveg(jl) <= 1.e-10_dp) THEN
      zdiff(jl) = 0._dp
      wsveg(jl) = 0._dp
    END IF
  END DO

! Calculation of diffusion coefficients
  DO jl=kidia, kfdia
    IF (cdel1 > 0) THEN
      zda  (jl) = zdiff(jl) * zdiagt / cdel1 / cdel1 
      zdc  (jl) = zda(jl) + 1._dp
      wsbs (jl) = wsbs(jl) / zdc(jl)

      ztri (jl) = -zda(jl) / zdc(jl)
      zdb  (jl) = zda(jl)+1._dp + zda(jl)*ztri(jl)
      wsveg(jl) = (wsveg(jl)/cdel1 + zda(jl)  &
            * wsbs(jl)/cdel1 ) / zdb(jl) * cdel1
      wsbs (jl) =  MAX(0.0_dp, wsbs(jl) - ztri(jl)*wsveg(jl))
      wsbs (jl) =  MIN(wsbs(jl), ws1(jl))
    END IF
  END DO

END SUBROUTINE baresoildiff 

