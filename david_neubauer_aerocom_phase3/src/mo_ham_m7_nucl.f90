!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7_nucl.f90
!!
!! \brief
!! mo_ham_m7_nucl provides routines for calculating aerosol nucleation rates.
!!
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Philip Stier, philip.stier@physics.ox.ac.uk
!!
!! \revision_history
!!   -# Philip Stier (MPI-Met) - original code (2003-01)
!!   -# Jan Kazil (MPI-Met) (2007-10-07)
!!
!! \limitations
!! None
!!
!! \details
!! This module provides routines for calculating nucleation rates for
!! neutral and charged H2SO4/H2O and neutral H2SO4/organic clusters.
!!
!! \bibliographic_references
!!   - Vignati, E. (1999), Modelling Interactions between Aerosols and 
!!     Gaseous Compounds in the Polluted Marine Atmosphere. PhD-Thesis,
!!     RISO National Laborartory Copenhagen, Riso-R-1163(EN)
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_m7_nucl
  
  USE mo_kind, ONLY: dp
  USE mo_netcdf
  
  IMPLICIT NONE
  
  ! Sulfate aerosol formation rates (cm-3 s-1) lookup table for the Kazil and
  ! Lovejoy nucleation scheme:
  
  INTEGER, PRIVATE :: t_n,rh_n,h2so4_n,ipr_n,cs_n
  
  REAL(dp), PRIVATE, ALLOCATABLE :: t_table(:) ! Temperature
  REAL(dp), PRIVATE, ALLOCATABLE :: rh_table(:) ! Relative humidity
  REAL(dp), PRIVATE, ALLOCATABLE :: h2so4_table(:) ! [H2SO4(g)]
  REAL(dp), PRIVATE, ALLOCATABLE :: ipr_table(:) ! Ionization rate
  REAL(dp), PRIVATE, ALLOCATABLE :: cs_table(:) ! H2SO4 condensation sink
      
  REAL(dp), PRIVATE, ALLOCATABLE :: pfr_table(:,:,:,:,:)
  
!!  INCLUDE 'netcdf.inc'
  
CONTAINS

SUBROUTINE nucl_vehkamaeki(kproma, kbdim,   klev,        & ! ECHAM5 dimensions
                         ptp1,     psatrat, pmolecH2SO4, & ! ECHAM5 temperature, water vapor saturaion ratio, [H2SO4]g
                         pxtrnucr, pntot               )   ! nucleation rate, number of molecules in the critical cluster
  !
  !   Authors:
  !   ---------
  !   C. TIMMRECK, MPI HAMBURG                                                    2002
  !
  !   Purpose
  !   ---------
  !   Calculation of classical nucleation rate
  !               
  !   calculation of the nucleation rate after Vehkamaeki et al. (2002)
  !   The calculation of the nucrate ZKNH2SO4 is in cm^-3 s^-1
  !   and a coarse approxmation for the first class
  !
  !   Modifications:
  !   --------------
  !   R. Hommel; rewrite in f90, adopted to ECHAM5;        MPI HAMBURG;      Dec. 2002
  !   P. Stier; bugfixes, modularisation and optimization; MPI HAMBURG;      2003-2004
  !
  !   H2SO4 still fixed to xxx molc/cm3, no sulfur cycle coupling yet
  !
  !   References:
  !   -----------
  !   Vehkamaeki et al. (2002), An improved parameterization for sulfuric
  !      acid/water nucleation rates for tropospheric and stratospheric
  !      conditions, J. Geophys. Res, 107, D22, 4622
  !
  !   Parameters
  !   ----------
  !   prho = prhop_neu in *sam*
  !
  !   psatrat = water vapor saturation ratio
  !
  !   pxtrnucr = nucleation rate in [1/m3s]
  !   xrhoc    = density of the critical nucleus in kg/m^3
  !   zrxc = ?

  !----------------------------------------------------

  IMPLICIT NONE

  !----------------------------------------------------
  
  INTEGER :: kproma, kbdim, klev
  
  INTEGER :: jk, jl
  
  !----------------------------------------------------
  !
  
  REAL(dp)::   ptp1(kbdim,klev), psatrat(kbdim,klev), &
               pxtrnucr(kbdim,klev),  &
               pmolecH2SO4(kbdim,klev), &            ! revisited, ok
               pntot(kbdim,klev)

  !----------------------------------------------------  
  ! Local Arrays
  
  REAL(dp)::   zrxc(kbdim) 

  REAL(dp)::   zrhoa, zrh, zt, x, zjnuc, zrc, zxmole, zntot

  REAL(dp)::   zlogrh, zlogrh2, zlogrh3, zlogrhoa, zlogrhoa2, zlogrhoa3, &
               zix, zt2, zt3

  !--- 0) Initializations:

  DO jk=1, klev
     DO jl=1,kproma

  !----1.) Parameterization of  nucleation rate after Vehkamaeki et al. (2002)

        ! t: temperature in K (190.15-300.15K)                                  
        ! zrh: saturatio ratio of water (0.0001-1)                               
        ! zrhoa: sulfuric acid concentration in 1/cm3 (10^4-10^11 1/cm3)         
        ! jnuc: nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)                  
        ! ntot: total number of molecules in the critical cluster (ntot>4)      
        ! x: molefraction of H2SO4 in the critical cluster                      
        ! rc: radius of the critical cluster in nm                              

        ! Calculate nucleation only for valid thermodynamic conditions:

        zrhoa = max(pmolecH2SO4(jl,jk),1.E+4_dp)
        zrhoa = min(zrhoa,1.E11_dp)
        
        zrh   = max(psatrat(jl,jk),1.E-4_dp)
        zrh   = min(zrh,1.0_dp)
        
        zt    = max(ptp1(jl,jk),190.15_dp)
        zt    = min(zt,300.15_dp)

        zt2 = zt*zt
        zt3 = zt2*zt

        ! Equation (11) - molefraction of H2SO4 in the critical cluster

        zlogrh  = LOG(zrh)
        zlogrh2 = zlogrh*zlogrh
        zlogrh3 = zlogrh2*zlogrh

        zlogrhoa  = LOG(zrhoa)
        zlogrhoa2 = zlogrhoa*zlogrhoa
        zlogrhoa3 = zlogrhoa2*zlogrhoa

        x=0.7409967177282139_dp - 0.002663785665140117_dp*zt   &
          + 0.002010478847383187_dp*zlogrh    &
          - 0.0001832894131464668_dp*zt*zlogrh    &
          + 0.001574072538464286_dp*zlogrh2        &
          - 0.00001790589121766952_dp*zt*zlogrh2    &
          + 0.0001844027436573778_dp*zlogrh3     &
          -  1.503452308794887e-6_dp*zt*zlogrh3    &
          - 0.003499978417957668_dp*zlogrhoa   &
          + 0.0000504021689382576_dp*zt*zlogrhoa

        zxmole=x

        zix = 1.0_dp/x

        ! Equation (12) - nucleation rate in 1/cm3s

        zjnuc=0.1430901615568665_dp + 2.219563673425199_dp*zt -   &
              0.02739106114964264_dp*zt2 +     &
              0.00007228107239317088_dp*zt3 + 5.91822263375044_dp*zix +     &
              0.1174886643003278_dp*zlogrh + 0.4625315047693772_dp*zt*zlogrh -     &
              0.01180591129059253_dp*zt2*zlogrh +     &
              0.0000404196487152575_dp*zt3*zlogrh +    &
              (15.79628615047088_dp*zlogrh)*zix -     &
              0.215553951893509_dp*zlogrh2 -    &
              0.0810269192332194_dp*zt*zlogrh2 +     &
              0.001435808434184642_dp*zt2*zlogrh2 -    &
              4.775796947178588e-6_dp*zt3*zlogrh2 -     &
              (2.912974063702185_dp*zlogrh2)*zix -   &
              3.588557942822751_dp*zlogrh3 +     &
              0.04950795302831703_dp*zt*zlogrh3 -     &
              0.0002138195118737068_dp*zt2*zlogrh3 +    &
              3.108005107949533e-7_dp*zt3*zlogrh3 -     &
              (0.02933332747098296_dp*zlogrh3)*zix +     &
              1.145983818561277_dp*zlogrhoa -    &
              0.6007956227856778_dp*zt*zlogrhoa +    &
              0.00864244733283759_dp*zt2*zlogrhoa -    &
              0.00002289467254710888_dp*zt3*zlogrhoa -    &
              (8.44984513869014_dp*zlogrhoa)*zix +    &
              2.158548369286559_dp*zlogrh*zlogrhoa +   &
              0.0808121412840917_dp*zt*zlogrh*zlogrhoa -    &
              0.0004073815255395214_dp*zt2*zlogrh*zlogrhoa -   &
              4.019572560156515e-7_dp*zt3*zlogrh*zlogrhoa +    &
              (0.7213255852557236_dp*zlogrh*zlogrhoa)*zix +    &
              1.62409850488771_dp*zlogrh2*zlogrhoa -    &
              0.01601062035325362_dp*zt*zlogrh2*zlogrhoa +   &
              0.00003771238979714162_dp*zt2*zlogrh2*zlogrhoa +    &
              3.217942606371182e-8_dp*zt3*zlogrh2*zlogrhoa -    &
              (0.01132550810022116_dp*zlogrh2*zlogrhoa)*zix +    &
              9.71681713056504_dp*zlogrhoa2 -    &
              0.1150478558347306_dp*zt*zlogrhoa2 +    &
              0.0001570982486038294_dp*zt2*zlogrhoa2 +    &
              4.009144680125015e-7_dp*zt3*zlogrhoa2 +    &
              (0.7118597859976135_dp*zlogrhoa2)*zix -    &
              1.056105824379897_dp*zlogrh*zlogrhoa2 +    &
              0.00903377584628419_dp*zt*zlogrh*zlogrhoa2 -    &
              0.00001984167387090606_dp*zt2*zlogrh*zlogrhoa2 +    &
              2.460478196482179e-8_dp*zt3*zlogrh*zlogrhoa2 -    &
              (0.05790872906645181_dp*zlogrh*zlogrhoa2)*zix -    &
              0.1487119673397459_dp*zlogrhoa3 +    &
              0.002835082097822667_dp*zt*zlogrhoa3 -    &
              9.24618825471694e-6_dp*zt2*zlogrhoa3 +    &
              5.004267665960894e-9_dp*zt3*zlogrhoa3 -    &
              (0.01270805101481648_dp*zlogrhoa3)*zix
        
        zjnuc=EXP(zjnuc)      !   add. Eq. (12) [1/(cm^3s)]      


        ! Equation (13) - total number of molecules in the critical cluster

        zntot=-0.002954125078716302_dp - 0.0976834264241286_dp*zt +   &
               0.001024847927067835_dp*zt2 - 2.186459697726116e-6_dp*zt3 -    &
               0.1017165718716887_dp*zix - 0.002050640345231486_dp*zlogrh -   &
               0.007585041382707174_dp*zt*zlogrh +    &
               0.0001926539658089536_dp*zt2*zlogrh -   &
               6.70429719683894e-7_dp*zt3*zlogrh -    &
               (0.2557744774673163_dp*zlogrh)*zix +   &
               0.003223076552477191_dp*zlogrh2 +   &
               0.000852636632240633_dp*zt*zlogrh2 -    &
               0.00001547571354871789_dp*zt2*zlogrh2 +   &
               5.666608424980593e-8_dp*zt3*zlogrh2 +    &
               (0.03384437400744206_dp*zlogrh2)*zix +   &
               0.04743226764572505_dp*zlogrh3 -    &
               0.0006251042204583412_dp*zt*zlogrh3 +   &
               2.650663328519478e-6_dp*zt2*zlogrh3 -    &
               3.674710848763778e-9_dp*zt3*zlogrh3 -   &
               (0.0002672510825259393_dp*zlogrh3)*zix -    &
               0.01252108546759328_dp*zlogrhoa +   &
               0.005806550506277202_dp*zt*zlogrhoa -    &
               0.0001016735312443444_dp*zt2*zlogrhoa +   &
               2.881946187214505e-7_dp*zt3*zlogrhoa +    &
               (0.0942243379396279_dp*zlogrhoa)*zix -   &
               0.0385459592773097_dp*zlogrh*zlogrhoa -   &
               0.0006723156277391984_dp*zt*zlogrh*zlogrhoa +   &
               2.602884877659698e-6_dp*zt2*zlogrh*zlogrhoa +    &
               1.194163699688297e-8_dp*zt3*zlogrh*zlogrhoa -   &
               (0.00851515345806281_dp*zlogrh*zlogrhoa)*zix -    &
               0.01837488495738111_dp*zlogrh2*zlogrhoa +   &
               0.0001720723574407498_dp*zt*zlogrh2*zlogrhoa -   &
               3.717657974086814e-7_dp*zt2*zlogrh2*zlogrhoa -    &
               5.148746022615196e-10_dp*zt3*zlogrh2*zlogrhoa +    &
               (0.0002686602132926594_dp*zlogrh2*zlogrhoa)*zix -   &
               0.06199739728812199_dp*zlogrhoa2 +    &
               0.000906958053583576_dp*zt*zlogrhoa2 -   &
               9.11727926129757e-7_dp*zt2*zlogrhoa2 -    &
               5.367963396508457e-9_dp*zt3*zlogrhoa2 -   &
               (0.007742343393937707_dp*zlogrhoa2)*zix +    &
               0.0121827103101659_dp*zlogrh*zlogrhoa2 -   &
               0.0001066499571188091_dp*zt*zlogrh*zlogrhoa2 +    &
               2.534598655067518e-7_dp*zt2*zlogrh*zlogrhoa2 -    &
               3.635186504599571e-10_dp*zt3*zlogrh*zlogrhoa2 +    &
               (0.0006100650851863252_dp*zlogrh*zlogrhoa2)*zix +   &
               0.0003201836700403512_dp*zlogrhoa3 -    &
               0.0000174761713262546_dp*zt*zlogrhoa3 +   &
               6.065037668052182e-8_dp*zt2*zlogrhoa3 -    &
               1.421771723004557e-11_dp*zt3*zlogrhoa3 +   &
               (0.0001357509859501723_dp*zlogrhoa3)*zix

        zntot=EXP(zntot)  !  add. Eq. (13)

          
        ! Equation (14) - radius of the critical cluster in nm

        zrc=EXP(-1.6524245_dp+0.42316402_dp*x+0.33466487_dp*LOG(zntot))    ! [nm]

        ! Conversion [nm -> m]

        zrxc(jl)=zrc*1e-9_dp

        !----1.2) Limiter

        IF(zjnuc<1.e-7_dp .OR. zntot<4.0_dp) zjnuc=0.0_dp

        ! limitation to 1E+10 [1/cm3s]
      
        zjnuc=MIN(zjnuc,1.e10_dp)

        pxtrnucr(jl,jk) = zjnuc

        ! convert total number of molecules in the critical cluster
        ! to number of sulfate molecules:

        pntot(jl,jk)=zntot*zxmole

      ENDDO ! kproma
      
    ENDDO ! klev
  
  END SUBROUTINE nucl_vehkamaeki
  
  !=============================================================================
  
  SUBROUTINE nucl_kazil_lovejoy(kproma,kbdim,klev, &
                                ptemp,prh,ph2so4,pcs,pipr, &
                                ppfr,pns)
    
    ! *nucl_kazil_lovejoy* returns the formation rate of sulfate aerosol
    ! particles containing a fixed number of H2SO4 molecules. The particle
    ! formation rate is determined by interpolating a lookup table that was
    ! generated with the code PARNUC, which implements the method described by
    ! Kazil and Lovejoy (2007). Currently. lookup tables for different
    ! nucleation processes (neutral and/or negative nucleation of H2SO4 nad H2O)
    ! and with different resolutions are available.
    !
    ! References:
    !
    ! Kazil J., and Lovejoy, E. R., A semi-analytical method for
    ! calculating rates of new sulfate aerosol formation from the gas phase,
    ! Atmos. Chem. Phys., 7, 3447-3459, 2007
  
    USE mo_kind, ONLY: dp
    
    IMPLICIT NONE
    
    !
    ! Input variables:
    !
    
    INTEGER :: kproma,             & !   "    geographic block number of locations
               kbdim,              & !   "    geographic block maximum number of locations
               klev                  !   "    numer of levels
    
    REAL(dp):: ptemp(kbdim,klev),  & ! Temperature (K)
               prh(kbdim,klev),    & ! Relative humidity (%)
               ph2so4(kbdim,klev), & ! H2SO4 gas phase concentration (cm-3)
               pcs(kbdim,klev),    & ! H2SO4 condensation sink (s-1)
               pipr(kbdim,klev)      ! Ion pair production rate (cm-3 s-1)
    
    !
    ! Output variables:
    !
    
    REAL(dp) :: ppfr(kbdim,klev)     ! Sulfate aerosol formation rate (cm-3 s-1)
    REAL(dp) :: pns(kbdim,klev)      ! Number of H2SO4 molecules in the newly formed aerosol particles
    
    !
    ! Local variables:
    !
    
    REAL(dp) :: ztemp,  & ! Temperature (K)
                zrh,    & ! Relative humidity (%)
                zh2so4, & ! H2SO4 gas phase concentration (cm-3)
                zipr,   & ! Ionization rate (cm-3 s-1)
                zcs       ! H2SO4 condensation sink (s-1)
    
    ! Interpolation coefficients:
    
    REAL(dp) :: zc(32)
    
    ! Temporary values:
    
    REAL(dp) :: zv,zw,zx,zy,zz, &
                zv_,zw_,zx_,zy_,zz_, &
                zvw,zxy,zv_w,zx_y,zvw_,zxy_,zv_w_,zx_y_, &
                zxyz,zxyz_,zx_yz,zx_yz_,zxy_z,zxy_z_,zx_y_z,zx_y_z_
    
    ! Array indices:
    
    INTEGER :: ik
    
    INTEGER :: i_temp_0,i_rh_0,i_h2so4_0,i_ipr_0,i_cs_0
    INTEGER :: i_temp_1,i_rh_1,i_h2so4_1,i_ipr_1,i_cs_1
    
    ! Loop indices:
    
    INTEGER :: ji,jj
    
    ! Logical variables:
    
    LOGICAL :: lerror,lset_zero
    
    ! Number of H2SO4 molecules in the newly formed aerosol:
    pns(1:kproma,:) = 15.0_dp
    
    !
    ! Calculate the nucleation rate by interpolating the lookup table:
    !
    
    DO ji=1,klev ! Levels loop
      DO jj=1,kproma ! Geographic locations loop
        
        ztemp  = ptemp(jj,ji)
        zrh    = prh(jj,ji)
        zh2so4 = ph2so4(jj,ji)
        zipr   = pipr(jj,ji)
        zcs    = pcs(jj,ji)
        
        ! Check if the current variables which particle formation depends upon
        ! are within the limits of the particle formation lookup table, and if
        ! not, take appropriate action:
        
        lerror    = .false.
        lset_zero = .false.
        
        IF (ztemp < t_table(1)) THEN
          lerror    = .true.
          ztemp    = t_table(1)
        ENDIF
        
        IF (ztemp > t_table(t_n)) THEN
          lerror    = .true.
          ztemp     = t_table(t_n)
          lset_zero = .true.
        ENDIF
        
        IF (zrh < rh_table(1)) THEN
          lerror    = .true.
          zrh       = rh_table(1)
          lset_zero = .true.
        ENDIF
        
        IF (zrh > rh_table(rh_n)) THEN
          lerror    = .true.
          zrh       = rh_table(rh_n)
        ENDIF
        
        IF (zh2so4 < h2so4_table(1)) THEN
          lerror    = .true.
          zh2so4    = h2so4_table(1)
          lset_zero = .true.
        ENDIF
        
        IF (zh2so4 > h2so4_table(h2so4_n)) THEN
          lerror    = .true.
          zh2so4    = h2so4_table(h2so4_n)
        ENDIF
        
        IF (zipr < ipr_table(1)) THEN
          lerror    = .true.
          zipr      = ipr_table(1)
        ENDIF
        
        IF (zipr > ipr_table(ipr_n)) THEN
          lerror    = .true.
          zipr      = ipr_table(ipr_n)
        ENDIF
        
        IF (zcs < cs_table(1)) THEN
          lerror    = .true.
          zcs       = cs_table(1)
        ENDIF
        
        IF (zcs > cs_table(cs_n)) THEN
          lerror    = .true.
          zcs       = cs_table(cs_n)
        ENDIF
        
        ! If needed, set the particle formation rate to zero:
        
        IF (lset_zero) THEN
          ppfr(jj,ji) = 0.0_dp
          CYCLE
        ENDIF
        
        ! Identify the intervals in the ambient conditions arrays where the
        ! current ambient conditions are located:
        
        i_temp_0 = 1
        i_temp_1 = t_n
        
        DO WHILE(i_temp_1-i_temp_0.gt.1)
          
          ik = (i_temp_1+i_temp_0)/2
          
          IF (t_table(ik).gt.ztemp) THEN
            i_temp_1 = ik
          ELSE
            i_temp_0 = ik
          ENDIF
          
        ENDDO
        
        i_rh_0 = 1
        i_rh_1 = rh_n
        
        DO WHILE(i_rh_1-i_rh_0.gt.1)
          
          ik = (i_rh_1+i_rh_0)/2
          
          IF (rh_table(ik).gt.zrh) THEN
            i_rh_1 = ik
          ELSE
            i_rh_0 = ik
          ENDIF
          
        ENDDO
        
        i_h2so4_0 = 1
        i_h2so4_1 = h2so4_n
        
        DO WHILE(i_h2so4_1-i_h2so4_0.gt.1)
          
          ik = (i_h2so4_1+i_h2so4_0)/2
          
          IF (h2so4_table(ik).gt.zh2so4) THEN
            i_h2so4_1 = ik
          ELSE
            i_h2so4_0 = ik
          ENDIF
          
        ENDDO
        
        i_ipr_0 = 1
        i_ipr_1 = ipr_n
        
        DO WHILE(i_ipr_1-i_ipr_0.gt.1)
          
          ik = (i_ipr_1+i_ipr_0)/2
          
          IF (ipr_table(ik).gt.zipr) THEN
            i_ipr_1 = ik
          ELSE
            i_ipr_0 = ik
          ENDIF
          
        ENDDO
        
        i_cs_0 = 1
        i_cs_1 = cs_n
        
        DO WHILE(i_cs_1-i_cs_0.gt.1)
          
          ik = (i_cs_1+i_cs_0)/2
          
          IF (cs_table(ik).gt.zcs) THEN
            i_cs_1 = ik
          ELSE
            i_cs_0 = ik
          ENDIF
          
        ENDDO
        
        ! Calculate coefficients for the interpolation of the lookup tables:
        
        zv = (ztemp-t_table(i_temp_0))/(t_table(i_temp_1)-t_table(i_temp_0))
        zw = (zrh-rh_table(i_rh_0))/(rh_table(i_rh_1)-rh_table(i_rh_0))
        zx = (zh2so4-h2so4_table(i_h2so4_0))/(h2so4_table(i_h2so4_1)-h2so4_table(i_h2so4_0))
        zy = (zipr-ipr_table(i_ipr_0))/(ipr_table(i_ipr_1)-ipr_table(i_ipr_0))
        zz = (zcs-cs_table(i_cs_0))/(cs_table(i_cs_1)-cs_table(i_cs_0))
        
        zv_ = zv - 1.0_dp
        zw_ = zw - 1.0_dp
        zx_ = zx - 1.0_dp
        zy_ = zy - 1.0_dp
        zz_ = zz - 1.0_dp
        
        zvw   = zv*zw
        zvw_  = zv*zw_
        zv_w  = zv_*zw
        zv_w_ = zv_*zw_
        
        zxy   = zx*zy
        zxy_  = zx*zy_
        zx_y  = zx_*zy
        zx_y_ = zx_*zy_
        
        zxyz    = zxy*zz
        zxyz_   = zxy*zz_
        zxy_z   = zxy_*zz
        zx_yz   = zx_y*zz
        zxy_z_  = zxy_*zz_
        zx_yz_  = zx_y*zz_
        zx_y_z  = zx_y_*zz
        zx_y_z_ = zx_y_*zz_
        
        zc(1)  = zv_w_*zx_y_z_
        zc(2)  = zv_w_*zx_y_z
        zc(3)  = zv_w_*zx_yz_
        zc(4)  = zv_w_*zx_yz
        zc(5)  = zv_w_*zxy_z_
        zc(6)  = zv_w_*zxy_z
        zc(7)  = zv_w_*zxyz_
        zc(8)  = zv_w_*zxyz
        zc(9)  = zv_w*zx_y_z_
        zc(10) = zv_w*zx_y_z
        zc(11) = zv_w*zx_yz_
        zc(12) = zv_w*zx_yz
        zc(13) = zv_w*zxy_z_
        zc(14) = zv_w*zxy_z
        zc(15) = zv_w*zxyz_
        zc(16) = zv_w*zxyz
        zc(17) = zvw_*zx_y_z_
        zc(18) = zvw_*zx_y_z
        zc(19) = zvw_*zx_yz_
        zc(20) = zvw_*zx_yz
        zc(21) = zvw_*zxy_z_
        zc(22) = zvw_*zxy_z
        zc(23) = zvw_*zxyz_
        zc(24) = zvw_*zxyz
        zc(25) = zvw*zx_y_z_
        zc(26) = zvw*zx_y_z
        zc(27) = zvw*zx_yz_
        zc(28) = zvw*zx_yz
        zc(29) = zvw*zxy_z_
        zc(30) = zvw*zxy_z
        zc(31) = zvw*zxyz_
        zc(32) = zvw*zxyz
         
        ! Interpolation of the particle formation rates table
        ! within the current 5D hypercuboid:
        
        ppfr(jj,ji) = &
         -  zc(1)*pfr_table(i_cs_0,i_ipr_0,i_h2so4_0,i_rh_0,i_temp_0) &
         +  zc(2)*pfr_table(i_cs_1,i_ipr_0,i_h2so4_0,i_rh_0,i_temp_0) &
         +  zc(3)*pfr_table(i_cs_0,i_ipr_1,i_h2so4_0,i_rh_0,i_temp_0) &
         -  zc(4)*pfr_table(i_cs_1,i_ipr_1,i_h2so4_0,i_rh_0,i_temp_0) &
         +  zc(5)*pfr_table(i_cs_0,i_ipr_0,i_h2so4_1,i_rh_0,i_temp_0) &
         -  zc(6)*pfr_table(i_cs_1,i_ipr_0,i_h2so4_1,i_rh_0,i_temp_0) &
         -  zc(7)*pfr_table(i_cs_0,i_ipr_1,i_h2so4_1,i_rh_0,i_temp_0) &
         +  zc(8)*pfr_table(i_cs_1,i_ipr_1,i_h2so4_1,i_rh_0,i_temp_0) &
         +  zc(9)*pfr_table(i_cs_0,i_ipr_0,i_h2so4_0,i_rh_1,i_temp_0) &
         - zc(10)*pfr_table(i_cs_1,i_ipr_0,i_h2so4_0,i_rh_1,i_temp_0) &
         - zc(11)*pfr_table(i_cs_0,i_ipr_1,i_h2so4_0,i_rh_1,i_temp_0) &
         + zc(12)*pfr_table(i_cs_1,i_ipr_1,i_h2so4_0,i_rh_1,i_temp_0) &
         - zc(13)*pfr_table(i_cs_0,i_ipr_0,i_h2so4_1,i_rh_1,i_temp_0) &
         + zc(14)*pfr_table(i_cs_1,i_ipr_0,i_h2so4_1,i_rh_1,i_temp_0) &
         + zc(15)*pfr_table(i_cs_0,i_ipr_1,i_h2so4_1,i_rh_1,i_temp_0) &
         - zc(16)*pfr_table(i_cs_1,i_ipr_1,i_h2so4_1,i_rh_1,i_temp_0) &
         + zc(17)*pfr_table(i_cs_0,i_ipr_0,i_h2so4_0,i_rh_0,i_temp_1) &
         - zc(18)*pfr_table(i_cs_1,i_ipr_0,i_h2so4_0,i_rh_0,i_temp_1) &
         - zc(19)*pfr_table(i_cs_0,i_ipr_1,i_h2so4_0,i_rh_0,i_temp_1) &
         + zc(20)*pfr_table(i_cs_1,i_ipr_1,i_h2so4_0,i_rh_0,i_temp_1) &
         - zc(21)*pfr_table(i_cs_0,i_ipr_0,i_h2so4_1,i_rh_0,i_temp_1) &
         + zc(22)*pfr_table(i_cs_1,i_ipr_0,i_h2so4_1,i_rh_0,i_temp_1) &
         + zc(23)*pfr_table(i_cs_0,i_ipr_1,i_h2so4_1,i_rh_0,i_temp_1) &
         - zc(24)*pfr_table(i_cs_1,i_ipr_1,i_h2so4_1,i_rh_0,i_temp_1) &
         - zc(25)*pfr_table(i_cs_0,i_ipr_0,i_h2so4_0,i_rh_1,i_temp_1) &
         + zc(26)*pfr_table(i_cs_1,i_ipr_0,i_h2so4_0,i_rh_1,i_temp_1) &
         + zc(27)*pfr_table(i_cs_0,i_ipr_1,i_h2so4_0,i_rh_1,i_temp_1) &
         - zc(28)*pfr_table(i_cs_1,i_ipr_1,i_h2so4_0,i_rh_1,i_temp_1) &
         + zc(29)*pfr_table(i_cs_0,i_ipr_0,i_h2so4_1,i_rh_1,i_temp_1) &
         - zc(30)*pfr_table(i_cs_1,i_ipr_0,i_h2so4_1,i_rh_1,i_temp_1) &
         - zc(31)*pfr_table(i_cs_0,i_ipr_1,i_h2so4_1,i_rh_1,i_temp_1) &
         + zc(32)*pfr_table(i_cs_1,i_ipr_1,i_h2so4_1,i_rh_1,i_temp_1)
        
        ppfr(jj,ji) = exp(ppfr(jj,ji))
        
      ENDDO
    ENDDO
    
  END SUBROUTINE nucl_kazil_lovejoy
  
  !=============================================================================
  
  SUBROUTINE nucl_activation(kproma,kbdim,klev,ph2so4,pforest,ppbl,ppfr,pns)
    
    ! *nucl_activation* returns the formation rate of aerosol particles from
    ! organic nucleation (activation type).
    !
    ! References:
    !
    ! Kulmala, M., Lehtinen, K. E. J., and Laaksonen, A.: Cluster activation
    ! theory as an explanation of the linear dependence between formation rate of
    ! 3nm particles and sulphuric acid concentration, Atmos. Chem. Phys., 6,
    ! 787-793, 2006
    ! 
    ! Riipinen, I., Sihto, S.-L., Kulmala, M., Arnold, F., Dal Maso, M., Birmili,
    ! W., Saarnio, K., Teinilä, K., Kerminen, V.-M., Laaksonen, A., and Lehtinen,
    ! K. E. J.: Connections between atmospheric sulphuric acid and new particle
    ! formation during QUEST III-IV campaigns in Heidelberg and Hyytiälä, 
    ! Atmos. Chem. Phys., 7, 1899-1914, 2007. 
  
    IMPLICIT NONE
    
    !
    ! Input variables:
    !
    
    INTEGER :: kproma,              & !   "    geographic block number of locations
               kbdim,               & !   "    geographic block maximum number of locations
               klev                   !   "    numer of levels
    
    REAL(dp):: ph2so4(kbdim,klev),  & ! Gas phase H2SO4 concentration [molec. cm-3]
               ppbl(kbdim),         & ! Planetary boundary layer top level
               pforest(kbdim)         ! Forest fraction
    
    !
    ! Output variables:
    !
    
    REAL(dp) :: ppfr(kbdim,klev)     ! Aerosol formation rate (cm-3 s-1)
    REAL(dp) :: pns(kbdim,klev)      ! Number of H2SO4 molecules in the newly formed aerosol particles
    
    !
    ! Local variables:
    !
    
    INTEGER :: jk,jl
    
    ! Activation nucleation (Kulmala et al., ACP 2006), parameterized by
    ! Riipinen et al., ACP 2007, using Hyytiälä QUEST II campaign median.
    ! Proceeds only over forests, in proportion to the forest fraction,
    ! and only in the boundary layer.
    
    DO jl = 1, kproma
      ! Above the boundary layer:
      DO jk = 1, int(ppbl(jl)) - 1
        pns(jl,jk)  = 0.0_dp
        ppfr(jl,jk) = 0.0_dp
      ENDDO
      ! In the boundary layer:
      DO jk = int(ppbl(jl)), klev
        pns(jl,jk)  = 1.0_dp
        ppfr(jl,jk) = 1.0E-6_dp*ph2so4(jl,jk)*pforest(jl)
      ENDDO
    ENDDO
    
  END SUBROUTINE nucl_activation
  
  !=============================================================================
  
  SUBROUTINE nucl_kinetic(kproma,kbdim,klev,ph2so4,pforest,ppbl,ppfr,pns)
    
    ! *nucl_kinetic* returns the formation rate of aerosol particles from
    ! organic nucleation (kinetic type).
    !
    ! References:
    !
    ! Kuang, C., P. H. McMurry, A. V. McCormick, and F. L. Eisele,
    ! Dependence of nucleation rates on sulfuric acid vapor concentration in
    ! diverse atmospheric locations, J. Geophys. Res., 113, D10209,
    ! doi:10.1029/2007JD009253, 2008
    ! 
    ! Laakso, L., Anttila, T., Lehtinen, K. E. J., Aalto, P. P., Kulmala,
    ! M., Hõrrak, U., Paatero, J., Hanke, M., and Arnold, F.: Kinetic nucleation
    ! and ions in boreal forest particle formation events, Atmos. Chem. Phys., 4,
    ! 2353-2366, 2004
    !
    ! Sihto, S.-L., Kulmala, M., Kerminen, V.-M., Dal Maso, M., Petäjä, T.,
    ! Riipinen, I., Korhonen, H., Arnold, F., Janson, R., Boy, M., Laaksonen,
    ! A., and Lehtinen, K. E. J.: Atmospheric sulphuric acid and aerosol
    ! formation: implications from atmospheric measurements for nucleation and
    ! early growth mechanisms, Atmos. Chem. Phys., 6, 4079-4091, 2006
    
    IMPLICIT NONE
    
    !
    ! Input variables:
    !
    
    INTEGER :: kproma,              & !   "    geographic block number of locations
               kbdim,               & !   "    geographic block maximum number of locations
               klev                   !   "    numer of levels
    
    REAL(dp):: ph2so4(kbdim,klev),  & ! Gas phase H2SO4 concentration [molec. cm-3]
               ppbl(kbdim),         & ! Planetary boundary layer top level
               pforest(kbdim)         ! Forest fraction
    
    !
    ! Output variables:
    !
    
    REAL(dp) :: ppfr(kbdim,klev)      ! Aerosol formation rate (cm-3 s-1)
    REAL(dp) :: pns(kbdim,klev)       ! Number of H2SO4 molecules in the newly formed aerosol particles
    
    !
    ! Local variables:
    !
    
    INTEGER :: jk,jl
    
    ! Kinetic nucleation (Laakso et al., ACP 2004), parameterized by Kuang et
    ! al., JGR 2008, using Hyytiälä QUEST II campaign data (Sihto et al., ACP
    ! 2006). Proceeds only over forests, in proportion to the forest fraction,
    ! and only in the boundary layer.
    
    DO jl = 1, kproma
      ! Above the boundary layer:
      DO jk = 1, int(ppbl(jl)) - 1
        pns(jl,jk)  = 0.0_dp
        ppfr(jl,jk) = 0.0_dp
      ENDDO
      ! In the boundary layer:
      DO jk = int(ppbl(jl)), klev
        pns(jl,jk)  = 2.0_dp
        ppfr(jl,jk) = 3.86E-13_dp*ph2so4(jl,jk)*ph2so4(jl,jk)*pforest(jl)
      ENDDO
    ENDDO
    
  END SUBROUTINE nucl_kinetic
  
!=============================================================================
  
  SUBROUTINE ham_nucl_initialize()
    
    ! *ham_nucl_initialize* reads data needed to calculate formation rates of
    ! aerosol particles from the gas phase.
    
    USE mo_ham_m7ctl, ONLY: nsnucl
    USE mo_netcdf
    USE mo_mpi, ONLY: p_parallel,p_parallel_io,p_io,p_bcast
    
    IMPLICIT NONE
    
    !
    ! Local variables:
    !
    
    CHARACTER(LEN=32) :: file
    
    INTEGER :: i_status,i_nc_id,i_dim_id,i_var_id
    
    INTEGER :: ji
               
    REAL(dp) :: zvalue
    
    IF (nsnucl==2) THEN ! Kazil and Lovejoy nucleation scheme
      
      ! Read tabulated formation rates of sulfate aerosol particles containing
      ! at least 15 H2SO4 molecules from the netCDF archive:
      
      file = 'parnuc.15H2SO4.nc'
      
      IF (p_parallel) THEN ! In parallel mode
        
        IF (p_parallel_io) THEN ! We are on the I/O processor
          
          i_status = nf_open(file,NF_NOWRITE,i_nc_id)
          CALL nf_check(i_status,file)
          
          !
          ! Read the lenghts of the dimensions of the particle formation rate tables:
          !
          
          CALL IO_INQ_DIMID(i_nc_id,'temperature',i_dim_id)
        	CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,t_n)
        	
          CALL IO_INQ_DIMID(i_nc_id,'RH',i_dim_id)
        	CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,rh_n)
        	
          ! As the dimension/variable name used for the H2SO4 gas phase
          ! concentration we try both '[H2SO4]' and 'H2SO4'. '[H2SO4]' is
          ! obsolete, as newer versions of the netCDF library do not allow
          ! special characters such as '[' as first character of a dimension/
          ! variable name, but lookup tables using '[H2SO4]' may still be in
          ! use.
          
          i_status = NF_INQ_DIMID(i_nc_id,'H2SO4',i_dim_id)
          IF (i_status /= nf_noerr) i_status = NF_INQ_DIMID(i_nc_id,'[H2SO4]',i_dim_id)
          CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,h2so4_n)
        	
          CALL IO_INQ_DIMID(i_nc_id,'ionization',i_dim_id)
        	CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,ipr_n)
        	
          CALL IO_INQ_DIMID(i_nc_id,'condensation_sink',i_dim_id)
          CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,cs_n)
          
          i_status = nf_close(i_nc_id)
          
        ENDIF
        
        ! Broadcast the dimensions:
        
        CALL p_bcast(t_n,p_io)
        CALL p_bcast(rh_n,p_io)
        CALL p_bcast(h2so4_n,p_io)
        CALL p_bcast(ipr_n,p_io)
        CALL p_bcast(cs_n,p_io)
        
      ELSE ! In serial mode
        
        i_status = nf_open(file,NF_NOWRITE,i_nc_id)
        CALL nf_check(i_status,file)
        
        !
        ! Read the lenghts of the dimensions of the particle formation rate tables:
        !
        
        CALL IO_INQ_DIMID(i_nc_id,'temperature',i_dim_id)
        CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,t_n)
        
        CALL IO_INQ_DIMID(i_nc_id,'RH',i_dim_id)
        CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,rh_n)
        
        ! As the dimension/variable name used for the H2SO4 gas phase
        ! concentration we try both '[H2SO4]' and 'H2SO4'. '[H2SO4]' is
        ! obsolete, as newer versions of the netCDF library do not allow
        ! special characters such as '[' as first character of a dimension/
        ! variable name, but lookup tables using '[H2SO4]' may still be in
        ! use.
        
        i_status = NF_INQ_DIMID(i_nc_id,'H2SO4',i_dim_id)
        IF (i_status /= nf_noerr) i_status = NF_INQ_DIMID(i_nc_id,'[H2SO4]',i_dim_id)
        CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,h2so4_n)
        
        CALL IO_INQ_DIMID(i_nc_id,'ionization',i_dim_id)
        CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,ipr_n)
        
        CALL IO_INQ_DIMID(i_nc_id,'condensation_sink',i_dim_id)
        CALL IO_INQ_DIMLEN(i_nc_id,i_dim_id,cs_n)
        
        i_status = nf_close(i_nc_id)
        
      ENDIF
      
      !
      ! Allocate memory for the ambient conditions and the particle
      ! formation rate arrays:
      !
      
      ALLOCATE(t_table(t_n))
      ALLOCATE(rh_table(rh_n))
      ALLOCATE(h2so4_table(h2so4_n))
      ALLOCATE(ipr_table(ipr_n))
      ALLOCATE(cs_table(cs_n))
      
      ALLOCATE(pfr_table(cs_n,ipr_n,h2so4_n,rh_n,t_n))
      
      IF (p_parallel) THEN ! In parallel mode
        
        IF (p_parallel_io) THEN ! We are on the I/O processor
          
          i_status = nf_open(file,NF_NOWRITE,i_nc_id)
          CALL nf_check(i_status,file)
          
          !
          ! Read the ambient conditions:
          !
          
          ! Temperature:
          
          CALL IO_INQ_VARID(i_nc_id,'temperature',i_var_id)
          
          DO ji = 1, t_n
            i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
            CALL nf_check(i_status,file)
            t_table(ji) = zvalue
          ENDDO
          
          ! RH:
          
          CALL IO_INQ_VARID(i_nc_id,'RH',i_var_id)
          
          DO ji = 1, rh_n
            i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
            CALL nf_check(i_status,file)
            rh_table(ji) = zvalue
          ENDDO
          
          ! [H2SO4]g:
          
          ! As the dimension/variable name used for the H2SO4 gas phase
          ! concentration we try both '[H2SO4]' and 'H2SO4'. '[H2SO4]' is
          ! obsolete, as newer versions of the netCDF library do not allow
          ! special characters such as '[' as first character of a dimension/
          ! variable name, but lookup tables using '[H2SO4]' may still be in
          ! use.
          
          i_status = NF_INQ_VARID(i_nc_id,'H2SO4',i_var_id)
          IF (i_status /= nf_noerr) i_status = NF_INQ_DIMID(i_nc_id,'[H2SO4]',i_var_id)
          
          DO ji = 1, h2so4_n
            i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
            CALL nf_check(i_status,file)
            h2so4_table(ji) = zvalue
          ENDDO
          
          ! Ionization rate:
          
          CALL IO_INQ_VARID(i_nc_id,'ionization',i_var_id)
          
          DO ji = 1, ipr_n
            i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
            CALL nf_check(i_status,file)
            ipr_table(ji) = zvalue
          ENDDO
          
          ! H2SO4 condensation sink:
          
          CALL IO_INQ_VARID(i_nc_id,'condensation_sink',i_var_id)
          
          DO ji = 1, cs_n
            i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
            CALL nf_check(i_status,file)
            cs_table(ji) = zvalue
          ENDDO
          
          !
          ! Read the particle formation rate:
          !
          
          CALL IO_INQ_VARID(i_nc_id,'pfr',i_var_id)
        	
          i_status = nf_get_var_double(i_nc_id,i_var_id,pfr_table)
          CALL nf_check(i_status,file)
          
          i_status = nf_close(i_nc_id)
          
        ENDIF
        
        ! Broadcast the ambient conditions:
        
        CALL p_bcast(t_table,p_io)
        CALL p_bcast(rh_table,p_io)
        CALL p_bcast(h2so4_table,p_io)
        CALL p_bcast(ipr_table,p_io)
        CALL p_bcast(cs_table,p_io)
        
        ! Broadcast the particle formation rates:
        
        CALL p_bcast(pfr_table,p_io)
        
      ELSE ! In serial mode
        
        i_status = nf_open(file,NF_NOWRITE,i_nc_id)
        CALL nf_check(i_status,file)
        
        !
        ! Read the ambient conditions:
        !
        
        ! Temperature:
        
        CALL IO_INQ_VARID(i_nc_id,'temperature',i_var_id)
        
        DO ji = 1, t_n
          i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
          CALL nf_check(i_status,file)
          t_table(ji) = zvalue
        ENDDO
        
        ! RH:
        
        CALL IO_INQ_VARID(i_nc_id,'RH',i_var_id)
        
        DO ji = 1, rh_n
          i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
          CALL nf_check(i_status,file)
          rh_table(ji) = zvalue
        ENDDO
        
        ! [H2SO4]g:
        
        ! As the dimension/variable name used for the H2SO4 gas phase
        ! concentration we try both '[H2SO4]' and 'H2SO4'. '[H2SO4]' is
        ! obsolete, as newer versions of the netCDF library do not allow
        ! special characters such as '[' as first character of a dimension/
        ! variable name, but lookup tables using '[H2SO4]' may still be in
        ! use.
        
        i_status = NF_INQ_VARID(i_nc_id,'H2SO4',i_var_id)
        IF (i_status /= nf_noerr) i_status = NF_INQ_DIMID(i_nc_id,'[H2SO4]',i_var_id)
        
        DO ji = 1, h2so4_n
          i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
          CALL nf_check(i_status,file)
          h2so4_table(ji) = zvalue
        ENDDO
        
        ! Ionization rate:
        
        CALL IO_INQ_VARID(i_nc_id,'ionization',i_var_id)
        
        DO ji = 1, ipr_n
          i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
          CALL nf_check(i_status,file)
          ipr_table(ji) = zvalue
        ENDDO
        
        ! H2SO4 condensation sink:
        
        CALL IO_INQ_VARID(i_nc_id,'condensation_sink',i_var_id)
        
        DO ji = 1, cs_n
          i_status = nf_get_var1_double(i_nc_id,i_var_id,ji,zvalue)
          CALL nf_check(i_status,file)
          cs_table(ji) = zvalue
        ENDDO
        
        !
        ! Read the particle formation rate:
        !
        
        CALL IO_INQ_VARID(i_nc_id,'pfr',i_var_id)
        
        i_status = nf_get_var_double(i_nc_id,i_var_id,pfr_table)
        CALL nf_check(i_status,file)
        
        i_status = nf_close(i_nc_id)
        
      ENDIF
      
    ENDIF
    
  END SUBROUTINE ham_nucl_initialize
  
  !=============================================================================
  
  SUBROUTINE ham_nucl_cleanup()
    
    ! *ham_nucl_cleanup* deallocates the memory of the module mo_ham_m7_nucl.
    
    USE mo_ham_m7ctl,   ONLY: nsnucl
    
    IMPLICIT NONE
    
    IF (nsnucl==2) THEN ! Kazil and Lovejoy nucleation scheme
      
      DEALLOCATE(t_table)
      DEALLOCATE(rh_table)
      DEALLOCATE(h2so4_table)
      DEALLOCATE(ipr_table)
      DEALLOCATE(cs_table)
      
      DEALLOCATE(pfr_table)
      
    ENDIF
    
  END SUBROUTINE ham_nucl_cleanup
  
END MODULE mo_ham_m7_nucl
