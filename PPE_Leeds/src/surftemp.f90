#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE surftemp(klon, pdt,                                         &
!*  CONSTANTS
        pemi, pboltz, pcp, pc16, platev, platsu,                       &
!*  COEFFICIENTS FROM THE ELIMINATION
               pfscoe, pescoe, pfqcoe, peqcoe,                         &
!*  OLD VALUES AT THE SURFACE
               psold, pqsold, pdqsold,                                 &
!*  OTHER FLUXES
               pnetrad, pgrdfl,                                        &
!*  DIFFUSION COEFFICIENTS, CAIR AND CSAT
!*                    FOR EVAP AND SOIL HEAT CAPACITY
               pcfh, pcair, pcsat, pfracsu, pgrdcap,                   &
!   Logical land mask
               lpland,                                                 &
!*  OUTPUT
               psnew, pqsnew)
!
!
!    COMPUTES THE ENERGY BALANCE AT THE SURFACE WITH AN IMPLICIT SCHEME
!    THAT IS CONNECTED TO THE RICHTMYER AND MORTON ALGORITHM OF THE PBL.
!
!     INPUT
!     -----
!     KLON     : LENGTH OF ARRAYS TO BE USED
!     PDT      : TIMESTEP IN SECONDS
!
!     PEMI     : SURFACE EMISSIVITY
!     PBOLTZ   : STEFAN-BOLTZMANN CONSTANT
!     PCP      : SPECIFIC HEAT OF AIR
!     PC16     : CPD*VTMPC2=CPD*(DELTA-1) (CF. VDIFF),FOR SENS.HEAT FL.
!     PLATEV   : LATENT HEAT OF EVAPORATION
!     PLATSU   : LATENT HEAT OF SUBLIMATION
!
!     PFSCOE, PESCOE : COEFFICIENTS OF THE RICHTMYER AND MORTON SCHEME
!                      FOR DRY STATIC ENERGY
!     PFQCOE, PEQCOE : AS ABOVE BUT FOR SPECIFIC HUMIDITY
!
!     PSOLD   : OLD SURFACE DRY STATIC ENERGY (TS * CP)
!     PQSOLD  : SATURATED  SPECIFIC HUMIDITY FOR OLD TEMPERATURE
!     PDQSOLD : DERIVATIVE OF SATURATED  SPECIFIC HUMIDITY AT THE
!               OLD TEMPERATURE
!
!     PNETRAD : NET RADIATION AT THE SURFACE (UPWARD LONGWAVE IS
!               INCLUDED BUT FOR THE OLD SURFACE TEMPERATURE)
!     PGRDFL  : GROUND HEAT FLUX
!
!     PCFH    : DIFFUSION COEFFICIENT FOR STATIC ENERGY AND MOISTURE
!     PCAIR   : COEFFICIENT IN LATENT HEAT FLUX FORMULA (SEE VDIFF)
!     PCSAT   : COEFFICIENT IN LATENT HEAT FLUX FORMULA (SEE VDIFF)
!     PFRACSU : FRACTION OF SURFACE FOR SUBLIMATION
!     PGRDCAP : SURFACE HEAT CAPACITY
!
!     OUTPUT
!     ------
!     PSNEW   : NEW SURFACE STATIC ENERGY
!     PQSNEW  : NEW SATURATED SURFACE AIR MOISTURE
!
!
!     AUTHOR.
!     -------
!
!     J. POLCHER  *LMD*  AND  J.-P. SCHULZ  *MPI*,  MAY 1995
!
!
!     MODIFICATIONS.
!     --------------
!
!     J.-P. SCHULZ  *MPI*,  OCTOBER 1997:
!        MODIFY ACCORDING TO LATENT HEAT FLUX FORMULATION IN VDIFF
!        USING ZCAIR AND ZCSAT COEFFICIENTS.
!
!     J.-P. SCHULZ  *MPI*,  AUGUST 1998:
!        MODIFY ACCORDING TO SENSIBLE HEAT FLUX FORMULATION IN VDIFF.
!
  USE mo_kind,  ONLY: dp
!
  IMPLICIT NONE
!
!* ARGUMENTS
!
  INTEGER :: klon
  REAL(dp):: pdt, pemi, pboltz, pcp(klon), pc16, platev, platsu
  REAL(dp):: pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
  REAL(dp):: psold(klon), pqsold(klon), pdqsold(klon)
  REAL(dp):: pnetrad(klon), pgrdfl(klon)
  REAL(dp):: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
  REAL(dp):: pgrdcap(klon)
  REAL(dp):: psnew(klon), pqsnew(klon)
  LOGICAL :: lpland(klon)
!
!
  INTEGER :: jl
  REAL(dp):: zcolin, zcohfl, zcoind, zicp, zca, zcs
!
!****************************************************************
!     MAIN PROGRAM
!****************************************************************
!
 DO jl = 1,klon
!
   IF (lpland(jl)) THEN
     zicp = 1._dp/pcp(jl)
!
     zca    = platsu*pfracsu(jl) +  platev*(pcair(jl) - pfracsu(jl))
     zcs    = platsu*pfracsu(jl) +  platev*(pcsat(jl) - pfracsu(jl))
!
     zcolin = pgrdcap(jl)*zicp +                                       &
                       pdt*(zicp*4._dp*pemi*pboltz*                    &
                       ((zicp*psold(jl))**3._dp) -                     &
                       pcfh(jl)*(zca*peqcoe(jl) - zcs -                &
                                 zicp*pc16*psold(jl)*                  &
                                 (pcair(jl)*peqcoe(jl) - pcsat(jl)))*  &
                            zicp*pdqsold(jl))
 
     zcohfl = -pdt*pcfh(jl)*(pescoe(jl)-1._dp)
!
     zcoind = pdt * (pnetrad(jl) + pcfh(jl)*pfscoe(jl) +  pcfh(jl)*    &
                       ((zca*peqcoe(jl)-zcs)*pqsold(jl) +              &
                       zca*pfqcoe(jl) - zicp*pc16*psold(jl)*           &
                      ((pcair(jl)*peqcoe(jl) - pcsat(jl))*pqsold(jl) + &
                          pcair(jl)*pfqcoe(jl))) + pgrdfl(jl))
!
!
     psnew(jl) = (zcolin * psold(jl) + zcoind) / (zcolin + zcohfl)
     pqsnew(jl) = pqsold(jl) + zicp * pdqsold(jl)*(psnew(jl) -         &
                                                           psold(jl))
   ELSE
     psnew(jl) = psold(jl)
     pqsnew(jl) = pqsold(jl)
   END IF
!
 END DO
!
  RETURN
END SUBROUTINE surftemp
