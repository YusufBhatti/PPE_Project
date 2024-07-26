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
SUBROUTINE collect (kproma                                             &
            , pahflw,  pahfsw,  pahfice                                &
            , ptrflw,  psoflw                                          &
            , pqres,   pevapw,  pevapi                                 &
            , pustrw,  pvstrw,  pustri,  pvstri                        &
            , palake,  pslf,    pseaice                                &
            , pwind10w,  pco2,   pco2_flux_ocean                       &
            , pawhea,  pawsol,  pawust                                 &
            , pawvst,  paicon,  paiqre,  paifre                        &
            , paiust,  paivst,  pawsta,  pco2atmos, pco2flux           &
            , prsfc,   pssfc,   prsfl,   pssfl          )
!
!  ---------------------------------------------------------------------
!
!  Collects surface values for input to the ocean model
!
!  *collect* is called from *physc*
!
!     Authors:
!
!     R. Voss, DKRZ,  August 1997, origianl source
!     R. Voss, U. Schlese, DKRZ, December 1999, modified for echam5
!     S. Legutke, DKRZ, July 2000, modifications for coupling with C-HOPE
!     I. Kirchner, MPI Hamburg, December 2000, time control
!     S. Legutke,  MPI,M&D, Jan  2002, time control for coupling revisited;
!                  exchange-field accumulation with timestep length weights
!     U. Schlese, M. Esch MPI, Jan 2003, 10m wind instead of ustar3
!     S. Legutke,  MPI,M&D, Apr  2003, coupling revisited;
!        - #ifdef cpl_hope   for coupling with HOPE
!        - #ifdef cpl_mpiom for coupling with MPIOM
!     M. Esch, MPI, Jan 2010, substitute cpl_co2 flag with namelist switch
!     V. Gayler, MPI, Mar 2014, accumulation on pawfre now in
!                               collect_hydrology_lake (mo_hydrology)
!                           
!     Method:
!
!     Fluxes calculated by ECHAM are accumulated and stored in arrays
!     which serve to transfer these data to *OASIS*.
!
USE mo_kind,         ONLY:  wp
USE mo_physical_constants,    ONLY:  alf, rhoh2o
USE mo_time_control, ONLY:  l_putocean, delta_time
USE mo_couple,       ONLY:  couple_a2o_time
USE mo_control,      ONLY:  lcouple_co2
!
IMPLICIT NONE
!
!  scalar arguments
!
  INTEGER, INTENT (IN) :: kproma
!
! array arguments
!
  REAL(wp) ::  pahflw(kproma),  pahfsw(kproma)                         &
            , ptrflw(kproma),  psoflw(kproma)                          &
            , pahfice(kproma), pqres(kproma)                           &
            , pevapw(kproma),  pevapi(kproma)                          &
            , pustrw(kproma),  pvstrw(kproma)                          &
            , pustri(kproma),  pvstri(kproma)                          &
            , palake(kproma),  pslf(kproma),    pseaice(kproma)        &
            , pwind10w(kproma), pco2(kproma), pco2_flux_ocean(kproma)  &
            , pawhea(kproma),  pawsol(kproma)                          &
            , pawust(kproma),  pawvst(kproma),  pawsta(kproma)         &
            , paicon(kproma),  paiqre(kproma),  paifre(kproma)         &
            , paiust(kproma),  paivst(kproma),  pco2atmos(kproma)      &
            , pco2flux(kproma)                                         &
            , prsfc(kproma),   pssfc(kproma),   prsfl(kproma)          &
            , pssfl(kproma)
!
! local scalars
!
  INTEGER :: jl

  REAL(wp) ::  zzf1, zzf2, zrcouple

!     accumulate variables for coupling
!
   DO jl=1,kproma
      IF (pslf(jl).LT.1._wp) THEN
        zzf1=1._wp-pseaice(jl)
        zzf2=      pseaice(jl)

#if defined  __cpl_mpiom
          pawhea(jl) = pawhea(jl)+( pahflw(jl)+pahfsw(jl)              &
                                   +ptrflw(jl)+psoflw(jl)              &
                                  -(pssfl(jl)+pssfc(jl))*alf*zzf1)     &
                                  *delta_time
          pawust(jl) = pawust(jl)+pustrw(jl)*delta_time
          pawvst(jl) = pawvst(jl)+pvstrw(jl)*delta_time
#endif
          pawsol(jl) = pawsol(jl)+psoflw(jl)*delta_time
          pawsta(jl) = pawsta(jl)+pwind10w(jl)*delta_time

#if defined __cpl_mpiom
          paicon(jl) = paicon(jl)+ pahfice(jl)*delta_time
          paiqre(jl) = paiqre(jl)+ pqres(jl)*delta_time
          paiust(jl) = paiust(jl)+ pustri(jl)*delta_time
          paivst(jl) = paivst(jl)+ pvstri(jl)*delta_time
#endif
        paifre(jl) = paifre(jl)+(pssfc(jl)+pssfl(jl)+pevapi(jl))       &
                               *zzf2*delta_time
        IF(lcouple_co2) THEN
          pco2atmos(jl) = pco2atmos(jl)+pco2(jl)*delta_time
          pco2flux(jl)  = pco2flux(jl)+pco2_flux_ocean(jl)*delta_time
        END IF
      END IF
    END DO
!
!    prepare coupling fields before transfer:
!            average and convert freshwater fluxes to [m/s]
!
   IF (l_putocean) THEN

     zrcouple = 1.0_wp/couple_a2o_time

     DO jl = 1,kproma
       pawhea(jl) = pawhea(jl)*zrcouple
       pawust(jl) = pawust(jl)*zrcouple
       pawvst(jl) = pawvst(jl)*zrcouple
       paicon(jl) = paicon(jl)*zrcouple
       paiqre(jl) = paiqre(jl)*zrcouple
       paifre(jl) = paifre(jl)*zrcouple/rhoh2o
       paiust(jl) = paiust(jl)*zrcouple
       paivst(jl) = paivst(jl)*zrcouple
       pawsol(jl) = pawsol(jl)*zrcouple
       pawsta(jl) = pawsta(jl)*zrcouple

       IF(lcouple_co2) THEN
         pco2atmos(jl) = pco2atmos(jl)*zrcouple
         pco2flux(jl)  = pco2flux(jl)*zrcouple
       END IF

     END DO

   END IF

   RETURN
   END
