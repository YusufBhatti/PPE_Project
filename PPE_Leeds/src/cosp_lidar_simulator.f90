!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
subroutine cosp_lidar_simulator(npoints,nlev,npart,nrefl & ! parasol
                , pres, presf, temp &
                , q_lsliq, q_lsice  &
                , ls_radliq, ls_radice  &
!!                , frac_out, ice_type &
                , ice_type &
                , pmol, pnorm, tautot, refl ) ! parasol

!---------------------------------------------------------------------------------
! Purpose: To compute lidar signal from model-simulated profiles of cloud water
!          and cloud fraction in each sub-column of each model gridbox.
!
! References: 
! Chepfer H., S. Bony, D. Winker, M. Chiriaco, J.-L. Dufresne, G. Seze (2008),
! Use of CALIPSO lidar observations to evaluate the cloudiness simulated by a 
! climate model, Geophys. Res. Lett., 35, L15704, doi:10.1029/2008GL034207.
!
! Previous references:
! Chiriaco et al, MWR, 2006; Chepfer et al., MWR, 2007
!
! Contacts: Helene Chepfer (chepfer@lmd.polytechnique.fr), Sandrine Bony (bony@lmd.jussieu.fr)
!
! May 2007: ActSim code of M. Chiriaco and H. Chepfer rewritten by S. Bony
!
! May 2008, H. Chepfer:
! - Units of pressure inputs: Pa 
! - Non Spherical particles : LS Ice NS coefficients, CONV Ice NS coefficients
! - New input: ice_type (0=ice-spheres ; 1=ice-non-spherical)
!
! June 2008, A. Bodas-Salcedo:
! - Ported to Fortran 90 and optimisation changes
!
! August 2008, J-L Dufresne:
! - Optimisation changes (sum instructions suppressed)
!
! October 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - Interface with COSP v2.0:
!      cloud fraction removed from inputs
!      in-cloud condensed water now in input (instead of grid-averaged value)
!      depolarisation diagnostic removed
!      parasol (polder) reflectances (for 5 different solar zenith angles) added
!
! Inputs:
!  npoints  : number of horizontal points
!  nlev : number of vertical levels
!  npart: numberb of cloud meteors (stratiform_liq, stratiform_ice, conv_liq, conv_ice). 
!        Currently npart must be 4
!  nrefl: number of solar zenith angles for parasol reflectances ! parasol
!  pres : pressure in the middle of atmospheric layers (full levels): Pa
!  presf: pressure in the interface of atmospheric layers (half levels): Pa
!     presf(..,1) : surface pressure ; presf(..,nlev+1)= TOA pressure
!  temp : temperature of atmospheric layers: K
!  q_lsliq: LS sub-column liquid water mixing ratio (kg/kg)
!  q_lsice: LS sub-column ice water mixing ratio (kg/kg)
!  ls_radliq: effective radius of LS liquid particles (meters)
!  ls_radice: effective radius of LS ice particles (meters)
!  frac_out: cloud cover in each sub-column of the gridbox (output from scops)
!  ice_type : ice particle shape hypothesis (ice_type=0 for spheres, ice_type=1 for non spherical particles)
!
! Outputs:
!  pmol : molecular attenuated backscatter lidar signal power (m^-1.sr^-1)
!  pnorm: total attenuated backscatter lidar signal power (m^-1.sr^-1)
!  tautot: optical thickess integrated from top to level z
!  refl : parasol(polder) reflectance ! parasol
!
! Version 1.0 (June 2007)
! Version 1.1 (May 2008)
! Version 1.2 (June 2008)
! Version 2.0 (October 2008)
! Version 2.1 (December 2008)
!---------------------------------------------------------------------------------

   USE mo_math_constants, ONLY: pi       ! pi
   USE mo_kind,           ONLY: wp

      IMPLICIT NONE

      LOGICAL ok_parasol
      PARAMETER (ok_parasol=.true.)  ! set to .true. if you want to activate parasol reflectances
                                     ! (caution: routine "parasol" not vectorized for the moment)

      INTEGER i, k
      
      INTEGER INDX_LSLIQ,INDX_LSICE
      PARAMETER (INDX_LSLIQ=1,INDX_LSICE=2)
! inputs:
      INTEGER npoints,nlev,npart,ice_type
      INTEGER nrefl ! parasol
      REAL(wp) ::  pres(npoints,nlev)    ! pressure full levels
      REAL(wp) ::  presf(npoints,nlev+1) ! pressure half levels
      REAL(wp) ::  temp(npoints,nlev)
      REAL(wp) ::  q_lsliq(npoints,nlev), q_lsice(npoints,nlev)
      REAL(wp) ::  ls_radliq(npoints,nlev), ls_radice(npoints,nlev)
!!      REAL(wp) ::  frac_out(npoints,nlev)

! outputs (for each subcolumn):

      REAL(wp) ::  pmol(npoints,nlev)  ! molecular backscatter signal power (m^-1.sr^-1)
      REAL(wp) ::  pnorm(npoints,nlev) ! total lidar backscatter signal power (m^-1.sr^-1)
      REAL(wp) ::  tautot(npoints,nlev)! optical thickess integrated from top
      REAL(wp) ::  refl(npoints,nrefl)! parasol reflectance ! parasol

! actsim variables:

      REAL(wp) ::  km, rdiffm, Qscat, Cmol
      PARAMETER (Cmol = 6.2446e-32_wp) ! depends on wavelength
      PARAMETER (km = 1.38e-23_wp)     ! Boltzmann constant (J/K)

      PARAMETER (rdiffm = 0.7_wp)      ! multiple scattering correction parameter
      PARAMETER (Qscat = 2.0_wp)       ! particle scattering efficiency at 532 nm

      REAL(wp) ::  rholiq, rhoice
      PARAMETER (rholiq=1.0e+03_wp)     ! liquid water (kg/m3)
      PARAMETER (rhoice=0.5e+03_wp)     ! ice (kg/m3)

      REAL(wp) ::  rhopart(npart)
      REAL(wp) ::  polpart(npart,5)  ! polynomial coefficients derived from Mie theory

!   grid-box variables:
      REAL(wp) ::  rad_part(npoints,nlev,npart)
      REAL(wp) ::  rhoair(npoints,nlev), zheight(npoints,nlev+1)
      REAL(wp) ::  beta_mol(npoints,nlev), alpha_mol(npoints,nlev)
      REAL(wp) ::  kp_part(npoints,nlev,npart)

!   sub-column variables:
!!      REAL(wp) ::  frac_sub(npoints,nlev)
      REAL(wp) ::  qpart(npoints,nlev,npart) ! mixing ratio particles in each subcolumn
      REAL(wp) ::  alpha_part(npoints,nlev,npart)
      REAL(wp) ::  tau_mol_lay(npoints)      ! temporary variable, moL. opt. thickness of layer k
      REAL(wp) ::  tau_mol(npoints,nlev)     ! optical thickness between TOA and bottom of layer k
      REAL(wp) ::  tau_part(npoints,nlev,npart)
      REAL(wp) ::  betatot(npoints,nlev)
      REAL(wp) ::  tautot_lay(npoints)   ! temporary variable, total opt. thickness of layer k
      REAL(wp) ::  tautot_S_liq(npoints),tautot_S_ice(npoints)     ! parasol


!------------------------------------------------------------
!---- 1. Preliminary definitions and calculations :
!------------------------------------------------------------

      if ( npart .ne. 2 ) then
        print *,'Error in lidar_simulator, npart should be 2, not',npart
        stop
      endif

! Polynomial coefficients for spherical liq/ice particles derived from Mie theory.
! Polynomial coefficients for non spherical particles derived from a composite of
! Ray-tracing theory for large particles (e.g. Noel et al., Appl. Opt., 2001)
! and FDTD theory for very small particles (Yang et al., JQSRT, 2003).

! We repeat the same coefficients for LS and CONV cloud to make code more readable
!*     LS Liquid water coefficients:
         polpart(INDX_LSLIQ,1) =  2.6980e-8_wp     
         polpart(INDX_LSLIQ,2) = -3.7701e-6_wp
         polpart(INDX_LSLIQ,3) =  1.6594e-4_wp
         polpart(INDX_LSLIQ,4) = -0.0024_wp
         polpart(INDX_LSLIQ,5) =  0.0626_wp
!*     LS Ice coefficients: 
      if (ice_type.eq.0) then     
         polpart(INDX_LSICE,1) = -1.0176e-8_wp   
         polpart(INDX_LSICE,2) =  1.7615e-6_wp
         polpart(INDX_LSICE,3) = -1.0480e-4_wp
         polpart(INDX_LSICE,4) =  0.0019_wp
         polpart(INDX_LSICE,5) =  0.0460_wp
      endif
!*     LS Ice NS coefficients: 
      if (ice_type.eq.1) then 
         polpart(INDX_LSICE,1) = 1.3615e-8_wp  
         polpart(INDX_LSICE,2) = -2.04206e-6_wp
         polpart(INDX_LSICE,3) = 7.51799e-5_wp
         polpart(INDX_LSICE,4) = 0.00078213_wp
         polpart(INDX_LSICE,5) = 0.0182131_wp
      endif


! density:
!*    clear-sky air:
      rhoair = pres/(287.04_wp*temp)

!*    liquid/ice particules:
      rhopart(INDX_LSLIQ) = rholiq
      rhopart(INDX_LSICE) = rhoice

! effective radius particles:
      rad_part(:,:,INDX_LSLIQ) = ls_radliq(:,:)
      rad_part(:,:,INDX_LSICE) = ls_radice(:,:)
      rad_part(:,:,:)=MAX(rad_part(:,:,:),0._wp)
      
! altitude at half pressure levels:
      zheight(:,1) = 0.0_wp
      do k = 2, nlev+1 ! Unvectorized loop : unimportant
        zheight(:,k) = zheight(:,k-1) &
                  -(presf(:,k)-presf(:,k-1))/(rhoair(:,k-1)*9.81_wp)
      enddo

! cloud fraction (0 or 1) in each sub-column:
! (if frac_out=1or2 -> frac_sub=1; if frac_out=0 -> frac_sub=0)
!!      frac_sub = MIN( frac_out, 1.0_wp )

!------------------------------------------------------------
!---- 2. Molecular alpha and beta:
!------------------------------------------------------------


      beta_mol = pres/km/temp * Cmol
      alpha_mol = 8.0_wp*pi/3.0_wp * beta_mol

!------------------------------------------------------------
!---- 3. Particles alpha and beta:
!------------------------------------------------------------

! polynomes kp_lidar derived from Mie theory:
      do i = 1, npart ! Unvectorized loop : unimportant
       where ( rad_part(:,:,i).gt.0.0_wp)
         kp_part(:,:,i) = &
            polpart(i,1)*(rad_part(:,:,i)*1e6_wp)**4 &
          + polpart(i,2)*(rad_part(:,:,i)*1e6_wp)**3 &
          + polpart(i,3)*(rad_part(:,:,i)*1e6_wp)**2 &
          + polpart(i,4)*(rad_part(:,:,i)*1e6_wp) &
          + polpart(i,5)
        elsewhere
         kp_part(:,:,i) = 0._wp
        endwhere
      enddo
      
! mixing ratio particules in each subcolumn:
          qpart(:,:,INDX_LSLIQ) = q_lsliq(:,:) ! oct08
          qpart(:,:,INDX_LSICE) = q_lsice(:,:) ! oct08


! alpha of particles in each subcolumn:
      do i = 1, npart ! Unvectorized loop : unimportant
        where ( rad_part(:,:,i).gt.0.0_wp)
          alpha_part(:,:,i) = 3.0_wp/4.0_wp * Qscat &
                 * rhoair(:,:) * qpart(:,:,i) &
                 / (rhopart(i) * rad_part(:,:,i) )
        elsewhere
          alpha_part(:,:,i) = 0._wp
        endwhere
      enddo

!------------------------------------------------------------
!---- 4. Backscatter signal:
!------------------------------------------------------------
! optical thickness (molecular):
!     opt. thick of each layer
      tau_mol(:,1:nlev) = alpha_mol(:,1:nlev) &
         & *(zheight(:,2:nlev+1)-zheight(:,1:nlev))
!     opt. thick from TOA
      DO k = nlev-1, 1, -1
        tau_mol(:,k) = tau_mol(:,k) + tau_mol(:,k+1)
      ENDDO

! optical thickness (particles):

      tau_part = rdiffm * alpha_part
      DO i = 1, npart
!       opt. thick of each layer
        tau_part(:,:,i) = tau_part(:,:,i) &
           & * (zheight(:,2:nlev+1)-zheight(:,1:nlev) )
!       opt. thick from TOA
        DO k = nlev-1, 1, -1 
          tau_part(:,k,i) = tau_part(:,k,i) + tau_part(:,k+1,i)
        ENDDO
      ENDDO

! molecular signal:
!      Upper layer 
       pmol(:,nlev) = beta_mol(:,nlev) / (2._wp*tau_mol(:,nlev)) &
            & * (1._wp-exp(-2.0_wp*tau_mol(:,nlev)))
!      Other layers
       DO k= nlev-1, 1, -1
        tau_mol_lay(:) = tau_mol(:,k)-tau_mol(:,k+1) ! opt. thick. of layer k
        WHERE (tau_mol_lay(:).GT.0._wp)
          pmol(:,k) = beta_mol(:,k) * EXP(-2.0_wp*tau_mol(:,k+1)) / (2._wp*tau_mol_lay(:)) &
            & * (1._wp-exp(-2.0_wp*tau_mol_lay(:)))
        ELSEWHERE
!         This must never happend, but just in case, to avoid div. by 0
          pmol(:,k) = beta_mol(:,k) * EXP(-2.0_wp*tau_mol(:,k+1))
        END WHERE
      END DO
!
! Total signal (molecular + particules):
!
! For performance reason on vector computers, the 2 following lines should not be used
! and should be replace by the later one.
!      betatot(:,:) = beta_mol(:,:) + sum(kp_part*alpha_part,dim=3)
!      tautot(:,:)  = tau_mol(:,:)  + sum(tau_part,dim=3)
      betatot(:,:) = beta_mol(:,:)
      tautot(:,:)  = tau_mol(:,:)
      DO i = 1, npart
           betatot(:,:) = betatot(:,:) + kp_part(:,:,i)*alpha_part(:,:,i)
           tautot(:,:) = tautot(:,:)  + tau_part(:,:,i)
      ENDDO ! i
!
!     Upper layer 
      pnorm(:,nlev) = betatot(:,nlev) / (2._wp*tautot(:,nlev)) &
            & * (1._wp-exp(-2.0_wp*tautot(:,nlev)))
!     Other layers
      DO k= nlev-1, 1, -1
        tautot_lay(:) = tautot(:,k)-tautot(:,k+1) ! optical thickness of layer k
        WHERE (tautot_lay(:).GT.0._wp)
          pnorm(:,k) = betatot(:,k) * EXP(-2.0_wp*tautot(:,k+1)) / (2._wp*tautot_lay(:)) &
               & * (1._wp-EXP(-2.0_wp*tautot_lay(:)))
        ELSEWHERE
!         This must never happend, but just in case, to avoid div. by 0
          pnorm(:,k) = betatot(:,k) * EXP(-2.0_wp*tautot(:,k+1))
        END WHERE
      END DO

!-------- End computation Lidar --------------------------
!---------------------------------------------------------
!  Parasol/Polder module
!
!  Purpose : Compute reflectance for one particular viewing direction
!  and 5 solar zenith angles (calculation valid only over ocean)
! ---------------------------------------------------------

! initialization:
    refl(:,:) = 0.0_wp

! activate parasol calculations:
    if (ok_parasol) then

!     Optical thickness from TOA to surface
      tautot_S_liq = 0._wp
      tautot_S_ice = 0._wp
      tautot_S_liq(:) = tautot_S_liq(:) &
         + tau_part(:,1,1)
      tautot_S_ice(:) = tautot_S_ice(:) &
         + tau_part(:,1,2)

      call parasol(npoints,nrefl        &
                 ,tautot_S_liq,tautot_S_ice &
                 ,refl)

    endif ! ok_parasol

END SUBROUTINE cosp_lidar_simulator
!
!---------------------------------------------------------------------------------
!
SUBROUTINE parasol(npoints,nrefl      &
                       ,tautot_S_liq,tautot_S_ice  &
                       ,refl)
!---------------------------------------------------------------------------------
! Purpose: To compute Parasol reflectance signal from model-simulated profiles 
!          of cloud water and cloud fraction in each sub-column of each model 
!          gridbox.
!
!
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - optimization for vectorization
!
! Version 2.0 (October 2008)
! Version 2.1 (December 2008)
!---------------------------------------------------------------------------------
    USE mo_kind,                ONLY: wp
    USE mo_math_constants,      ONLY: pi !pi 

    IMPLICIT NONE

! inputs
    INTEGER npoints              ! Number of horizontal gridpoints
    INTEGER nrefl                ! Number of angles for which the reflectance 
                                 ! is computed. Can not be greater then ntetas
    REAL(wp)  tautot_S_liq(npoints)   ! liquid water cloud optical thickness, 
                                   ! integrated from TOA to surface
    REAL(wp)  tautot_S_ice(npoints)   ! same for ice water clouds only
! outputs
    REAL(wp)  refl(npoints,nrefl)     ! Parasol reflectances
!
! Local variables
    REAL(wp)  tautot_S(npoints)       ! cloud optical thickness, from TOA to surface
    REAL(wp)  frac_taucol_liq(npoints), frac_taucol_ice(npoints)


!   look up table variables:
    INTEGER ny, it
    INTEGER ntetas, nbtau        ! number of angle and of optical thickness
                                   ! of the look-up table
    PARAMETER (ntetas=5, nbtau=7)
    REAL(wp)  aa(ntetas,nbtau-1), ab(ntetas,nbtau-1)
    REAL(wp)  ba(ntetas,nbtau-1), bb(ntetas,nbtau-1)  
    REAL(wp)  tetas(ntetas),optau(nbtau)                        
    REAL(wp)  r_norm(ntetas)
    REAL(wp)  rlumA(ntetas,nbtau), rlumB(ntetas,nbtau)       
    REAL(wp)  rlumA_mod(npoints,5), rlumB_mod(npoints,5) 

    DATA optau   /0._wp, 1._wp, 5._wp, 10._wp, 20._wp, 50._wp, 100._wp/
    DATA tetas /0._wp, 20._wp, 40._wp, 60._wp, 80._wp/
    
! Look-up table for spherical liquid particles:
    DATA (rlumA(1,ny),ny=1,nbtau) /0.03_wp, 0.090886_wp, 0.283965_wp, &
     0.480587_wp, 0.695235_wp, 0.908229_wp, 1.0_wp /
    DATA (rlumA(2,ny),ny=1,nbtau) /0.03_wp, 0.072185_wp, 0.252596_wp, &
      0.436401_wp,  0.631352_wp, 0.823924_wp, 0.909013_wp /
    DATA (rlumA(3,ny),ny=1,nbtau) /0.03_wp, 0.058410_wp, 0.224707_wp, &
      0.367451_wp,  0.509180_wp, 0.648152_wp, 0.709554_wp /
    DATA (rlumA(4,ny),ny=1,nbtau) /0.03_wp, 0.052498_wp, 0.175844_wp, &
      0.252916_wp,  0.326551_wp, 0.398581_wp, 0.430405_wp /
    DATA (rlumA(5,ny),ny=1,nbtau) /0.03_wp, 0.034730_wp, 0.064488_wp, &
      0.081667_wp,  0.098215_wp, 0.114411_wp, 0.121567_wp /

! Look-up table for ice particles:
    DATA (rlumB(1,ny),ny=1,nbtau) /0.03_wp, 0.092170_wp, 0.311941_wp, &
       0.511298_wp, 0.712079_wp , 0.898243_wp , 0.976646_wp  /
    DATA (rlumB(2,ny),ny=1,nbtau) /0.03_wp, 0.087082_wp, 0.304293_wp, &
       0.490879_wp,  0.673565_wp, 0.842026_wp, 0.912966_wp /
    DATA (rlumB(3,ny),ny=1,nbtau) /0.03_wp, 0.083325_wp, 0.285193_wp, &
      0.430266_wp,  0.563747_wp, 0.685773_wp,  0.737154_wp /
    DATA (rlumB(4,ny),ny=1,nbtau) /0.03_wp, 0.084935_wp, 0.233450_wp, &
      0.312280_wp, 0.382376_wp, 0.446371_wp, 0.473317_wp /
    DATA (rlumB(5,ny),ny=1,nbtau) /0.03_wp, 0.054157_wp, 0.089911_wp, &
      0.107854_wp, 0.124127_wp, 0.139004_wp, 0.145269_wp /

!--------------------------------------------------------------------------------
! Lum_norm=f(tetaS,tau_cloud) derived from adding-doubling calculations
!        valid ONLY ABOVE OCEAN (albedo_sfce=5%)
!        valid only in one viewing direction (theta_v=30 G� @, phi_s-phi_v=320 G� @)
!        based on adding-doubling radiative transfer computation
!        for tau values (0 to 100) and for tetas values (0 to 80)
!        for 2 scattering phase functions: liquid spherical, ice non spherical

    IF ( nrefl.GT. ntetas ) THEN
        PRINT *,'Error in lidar_simulator, nrefl should be less then ',ntetas,' not',nrefl
        STOP
    ENDIF

    rlumA_mod=0._wp
    rlumB_mod=0._wp
!
   
    r_norm(:)=1._wp/ COS(pi/180._wp*tetas(:))
!
    tautot_S_liq(:)=MAX(tautot_S_liq(:),optau(1))
    tautot_S_ice(:)=MAX(tautot_S_ice(:),optau(1))
    tautot_S(:) = tautot_S_ice(:) + tautot_S_liq(:)
!
! relative fraction of the opt. thick due to liquid or ice clouds
    WHERE (tautot_S(:) .GT. 0._wp)
        frac_taucol_liq(:) = tautot_S_liq(:) / tautot_S(:)
        frac_taucol_ice(:) = tautot_S_ice(:) / tautot_S(:)
    ELSEWHERE
        frac_taucol_liq(:) = 1._wp
        frac_taucol_ice(:) = 0._wp
    END WHERE
    tautot_S(:)=MIN(tautot_S(:),optau(nbtau))
!
! Linear interpolation :

    DO ny=1,nbtau-1
! microphysics A (liquid clouds) 
      aA(:,ny) = (rlumA(:,ny+1)-rlumA(:,ny))/(optau(ny+1)-optau(ny))
      bA(:,ny) = rlumA(:,ny) - aA(:,ny)*optau(ny)
! microphysics B (ice clouds)
      aB(:,ny) = (rlumB(:,ny+1)-rlumB(:,ny))/(optau(ny+1)-optau(ny))
      bB(:,ny) = rlumB(:,ny) - aB(:,ny)*optau(ny)
    ENDDO
!
    DO it=1,ntetas
      DO ny=1,nbtau-1
        WHERE (tautot_S(:).GE.optau(ny).AND.tautot_S(:).LE.optau(ny+1))
            rlumA_mod(:,it) = aA(it,ny)*tautot_S(:) + bA(it,ny)
            rlumB_mod(:,it) = aB(it,ny)*tautot_S(:) + bB(it,ny)
        END WHERE
      END DO
    END DO
!
    DO it=1,ntetas
      refl(:,it) = frac_taucol_liq(:) * rlumA_mod(:,it) &
         + frac_taucol_ice(:) * rlumB_mod(:,it)
! normalized radiance -> reflectance: 
      refl(:,it) = refl(:,it) * r_norm(it)
    ENDDO

    RETURN
  END SUBROUTINE parasol

