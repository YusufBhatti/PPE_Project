!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_perturbations.f90
!!
!! \brief
!! mo_hammoz_perturbations provides a namelist interfact to parameters used in
!! parametric perturbation experiments.
!!
!! \author Philip Stier (Oxford)
!!
!! \responsible_coder
!! Philip Stier (Oxford), philip.stier@physics.ox.ac.uk
!!
!! \revision_history
!!   -# P. Stier (Oxford) - original version - (2013-09)
!!
!! \limitations
!! None
!!
!! \details
!! Perturbation design is largely based on Lee et al. (2013)
!!
!! \bibliographic_references
!! Lee, L. A., K. J. Pringle, C. L. Reddington, G. W. Mann, P. Stier, D. V. Spracklen,
!!   J. R. Pierce and K. S. Carslaw, The magnitude and causes of uncertainty in global
!!   model simulations of cloud condensation nuclei, Atmos. Chem. Phys., 13, 8879–8914,
!!   2013, www.atmos-chem-phys.net/13/8879/2013/ doi:10.5194/acp-13-8879-2013.
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

MODULE mo_hammoz_perturbations

  USE mo_kind,                ONLY: dp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init_hammoz_perturbations, lo_hammoz_perturbations, &
            init_hammoz_emi_perturbations, &
            scale_nuc_ft, &
            scale_emi_cmr_ff, scale_emi_cmr_bb, &
	    scale_drydep_acc,  &
            scale_wetdep_ic, scale_wetdep_bc, & 
            bc_rad_ni, du_rad_ni, oc_rad_ni, &
            scale_so4_coating, &
	    scale_so2_reactions, &
            scale_kappa_so4, scale_kappa_oc, &
	    scale_vertical_velocity, scale_dms_sc, &
	    scale_seasalt_expo


  LOGICAL :: lo_hammoz_perturbations=.TRUE.

  !--- Parameter list to perturb with short description:

  REAL(dp)    :: ft_nuc, & ! scalable - done
                 ageing, & ! scalable - done
                 act_diam, & ! tricky as not in HAM
                 drydep_aer_acc, & ! scalable - done
                 acc_width, & !tricky with radiation, fine otherwise? Check what they have done!
                 ait_width, & !tricky with radiation, fine otherwise? Check what they have done!
                 nucait_width, & !tricky with radiation, fine otherwise? Check what they have done!
                 aitacc_width, & !tricky with radiation, fine otherwise? Check what they have done!
                 ff_ems, & ! scalable - done through scaling of emission scale factors in init_hammoz_emi_perturbations
                 bb_ems, & ! scalable - done through scaling of emission scale factors in init_hammoz_emi_perturbations
                 bf_ems, & ! scalable - done through scaling of emission scale factors in init_hammoz_emi_perturbations
                 ff_diam, & ! scalable - done in mo_ham_m7_emissions
                 bb_diam, & ! scalable - done in mo_ham_m7_emissions
                 bf_diam, & ! scalable - done in mo_ham_m7_emissions
                 ss_acc, & ! scalable
                 anth_s02, & ! scalable
                 volc_so2, & ! scalable
                 dms_flux, & ! scalable
                 bio_soa, & ! scalable
                 anth_soa ! ?

  REAL(dp)    :: scale_nuc_ft = 1.0_dp, &     ! Scale factor for free tropospheric nucleation scheme (Vehkamaeki or Kazil & Lovejoy)
                 scale_emi_cmr_ff = 1.0_dp, & ! Scale factor for emission count median raduius for fossil fules
                 scale_emi_cmr_bb = 1.0_dp, & ! Scale factor for emission count median raduius for wildfires
                 scale_emi_ff = 1.0_dp,     & ! Scale factor for fossil fule emissions
                 scale_emi_bb = 1.0_dp,     & ! Scale factor for wildfire emissions
                 scale_emi_bf = 1.0_dp,     & ! Scale factor for biofuel emissions
                 scale_emi_dms = 1.0_dp,    &   ! Scale factor for dms emissions (terrestrial + oceanic)
                 scale_emi_ssa = 1.0_dp,     &! Scale factor for SSA emissions
                 scale_emi_du = 1.0_dp,    &! Scale factor for dust emissions
                 scale_emi_so2 = 1.0_dp,   & ! Scale factor for ANTH SO2 emissions
                 scale_so2_reactions = 1.0_dp, &! Scale factor for all SO2 reactions
                 scale_dms_sc = 1.0_dp,  &    ! Scale factor schmidt number ratio of DMS
                 scale_seasalt_expo = 1.0_dp    ! Scale factor schmidt number ratio of DMS

  REAL(dp)    :: scale_drydep_acc = 1.0_dp ! Scale factor for dry deposition of accumulation modes

  REAL(dp)    :: scale_wetdep_ic = 1.0_dp, & ! Scale factor for in-cloud wet deposition
                 scale_wetdep_bc = 1.0_dp  ! Scale factor for below-cloud wet deposition

  REAL(dp)    :: scale_vertical_velocity = 1.0_dp ! Scale (total) vertical velocity - ONLY for nactivpdf == 0

  REAL(dp)    :: kappa_ss  = 1.0_dp, & ! Sea salt kappa
                 kappa_so4 = 0.6_dp, & ! Sulfate kappa
                 kappa_oc  = 0.06_dp   ! OC kappa

  REAL(dp)    :: bc_rad_ni = 0.71_dp, & ! Absolute value of imaginary part of the refractive index (for BC) at 550 nm (default)
                 oc_rad_ni = 0.0055_dp, & ! Absolute value of imaginary part of the refractive index (for OC) at 550 nm (default)
                 du_rad_ni = 0.001_dp   ! Absolute value of imaginary part of the refractive index (for DUST) at 550 nm (default)


  REAL(dp)    :: scale_so4_coating = 1.0_dp  ! Scale the coating thickness of SO4 required to 'age' particles (move them
                                             !  from insoluble to soluble modes)
  REAL(dp)    :: scale_kappa_so4 = 0.6_dp


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! setham modifies pre-set switches of the hamctl namelist for the
!! configuration of the ECHAM/HAM aerosol model
!!
!! @author see above
!!
!! $Id: ????$
!!
!! @par Revision History
!! see above
!!
!! @par This subroutine is called by
!! init_ham
!!
!! @par Externals:
!! <ol>
!! <li>None
!! </ol>
!!
!! @par Notes
!!
!!
!! @par Responsible coder
!! philip.stier@physics.ox.ac.uk
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE init_hammoz_perturbations

    ! *init_ham_perturbations* reads the perturbation namelist and populates the
    !           perturbation scaling factors
    !
    ! Authors:
    ! --------
    ! Philip Stier, Oxford                        09/2013

    USE mo_mpi,                 ONLY: p_parallel, p_parallel_io, p_bcast, p_io
    USE mo_namelist,            ONLY: open_nml, position_nml, POSITIONED
    USE mo_exception,           ONLY: finish, message, message_text, em_info, em_error, em_param,  &
                                      em_info, em_warn
    USE mo_util_string,         ONLY: separator
    USE mo_submodel,            ONLY: print_value
    USE mo_submodel,            ONLY: lham, lhammoz
    USE mo_util_string,         ONLY: separator
    USE mo_srtm_setup,          ONLY: ssi_amip

    IMPLICIT NONE


    INCLUDE 'hammoz_perturbations.inc'

!!$    INTEGER, INTENT(in)       :: nmod         ! number of modes/bins

    !--- Local variables

    CHARACTER(len=24)         :: csubmname    ! name of aerosol sub model
    INTEGER                   :: jj, ierr, inml, iunit
    REAL(dp)                  :: fire_prop_tot


    !--- 1) Read namelist:

    CALL message('',separator)
    CALL message('init_hammoz_perturbations', 'Reading namelist hammoz_perturbation.inc...', level=em_info)

    IF (p_parallel_io) THEN
       inml = open_nml('namelist.echam')
       iunit = position_nml ('HAMMOZ_PERTURBATIONS', inml, status=ierr)
       SELECT CASE (ierr)
       CASE (POSITIONED)
          READ (iunit, hammoz_perturbations)
       CASE DEFAULT
          WRITE(message_text,'(a,i0)') 'Namelist hammoz_perturbations not correctly read! ierr = ', ierr
          CALL finish('setham', message_text)
       END SELECT
    ENDIF

    !--- 2) Broadcast over processors:
    IF (p_parallel) THEN
!       CALL p_bcast (scale_nuc_bl,         p_io)
!       CALL p_bcast (scale_nuc_ft,         p_io)
!       CALL p_bcast (scale_emi_cmr_ff,     p_io)
!       CALL p_bcast (scale_emi_cmr_bb,     p_io)
!       CALL p_bcast (scale_emi_cmr_bf,     p_io)
!       CALL p_bcast (scale_emi_ff,         p_io)
!       CALL p_bcast (scale_emi_bb,         p_io)
!       CALL p_bcast (scale_emi_bf,         p_io)
!       CALL p_bcast (scale_emi_bc,        p_io)
!       CALL p_bcast (scale_emi_dms,        p_io)
!       CALL p_bcast (scale_emi_ssa,        p_io)
!       CALL p_bcast (scale_emi_du,        p_io)
!       CALL p_bcast (scale_emi_so2,        p_io)
       CALL p_bcast (scale_drydep_acc,     p_io)
!       CALL p_bcast (scale_wetdep_ic,      p_io)
     !  CALL p_bcast (scale_wetdep_ic_bc_only,    p_io)
     !  CALL p_bcast (scale_wetdep_bc_bc_only,     p_io)
!       CALL p_bcast (bc_rad_ni,      p_io)
!       CALL p_bcast (oc_rad_ni,      p_io) 
!       CALL p_bcast (du_rad_ni,      p_io) 
     !  CALL p_bcast (prop_fire_in_pbl_p1,  p_io)
     !  CALL p_bcast (prop_fire_in_pbl_p2,  p_io)
     !  CALL p_bcast (prop_fire_in_pbl,     p_io)
!       CALL p_bcast (scale_so4_coating,    p_io)
!       CALL p_bcast (pH_pert,      p_io)       
!    !   CALL p_bcast (scale_tr_entrainment, p_io)
!       CALL p_bcast (scale_vertical_velocity, p_io)
!    !   CALL p_bcast (scale_intra_mode_coagulation, p_io)
!    !   CALL p_bcast (scale_inter_mode_coagulation, p_io)
!       ! Solar constant
!       CALL p_bcast (scale_solar_const, p_io)
       ! Cloud params
   !    CALL p_bcast (KK_exponent, p_io)
   !    CALL p_bcast (KK_LWP_exponent, p_io)
   !    CALL p_bcast (scale_activation, p_io)
   !    CALL p_bcast (scale_accretion, p_io)
   !    CALL p_bcast (scale_water, p_io)
 !      CALL p_bcast (scale_RH, p_io)
!       CALL p_bcast (scale_kappa_ss, p_io)
!       CALL p_bcast (scale_kappa_so4, p_io)
!       ! SO2 chemistry
!       CALL p_bcast (scale_so2_reactions,        p_io)
!       CALL p_bcast (scale_dms_sc,        p_io)
!       CALL p_bcast (scale_seasalt_expo,        p_io)


    END IF

    !---------------------------------------------------------------------------------------------------
    !--- 3) Consistency and dependency checks:

!    fire_prop_tot = prop_fire_in_pbl_p1 + prop_fire_in_pbl_p2 + prop_fire_in_pbl
!    IF (ABS(fire_prop_tot - 1.0_dp) > 0.01) THEN
!        WRITE(message_text,'(a,f6.3)') 'Total proportion of fire emissions emitted /= 1.0, total = ', fire_prop_tot
!        CALL finish('setham', message_text)
!    ENDIF
    !--- 4) Write out parameter settings

    csubmname = 'UNKNOWN'
    IF (lham) csubmname = 'HAM'
    IF (lhammoz) csubmname = 'HAMMOZ'

    CALL message('','')
    CALL message('', separator)
    CALL message('hammoz_perturbations','Parameter settings for the ECHAM-'//TRIM(csubmname)//' model', &
                 level=em_info)

    CALL message('','---')
    CALL print_value('Nucleation scaling factor (scale_nuc_bl)', scale_nuc_bl)
    CALL print_value('Nucleation scaling factor (scale_nuc_ft)', scale_nuc_ft)

    CALL message('','---')
    CALL print_value('Emission radius scaling factor (scale_emi_cmr_ff)', scale_emi_cmr_ff)
    CALL print_value('Emission radius scaling factor (scale_emi_cmr_bb)', scale_emi_cmr_bb)
    CALL print_value('Emission radius scaling factor (scale_emi_cmr_bf)', scale_emi_cmr_bf)

    CALL message('','---')
    CALL print_value('Emission of primary SO4 - fraction of SO2 emitted as SO4', emi_prim_so4_frac)
    CALL print_value('Emission of primary SO4 - radius of emitted SO4         ', emi_prim_so4_cmr)

    CALL message('','---')
    CALL print_value('Emission scaling factor (scale_emi_ff)', scale_emi_ff)
    CALL print_value('Emission scaling factor (scale_emi_bb)', scale_emi_bb)
    CALL print_value('Emission scaling factor (scale_emi_bf)', scale_emi_bf)
    CALL print_value('Emission scaling factor (scale_emi_bc)', scale_emi_bc)
    CALL print_value('Emission scaling factor (scale_emi_dms)', scale_emi_dms)
    CALL print_value('Emission scaling factor (scale_emi_ssa)', scale_emi_ssa)
    CALL print_value('Emission scaling factor (scale_emi_du)', scale_emi_du)
    CALL print_value('Emission scaling factor (scale_emi_so2)', scale_emi_so2)

    !dwp NOTE that the scale_emi_bc happens *before* the other scalings to avoid circular dependencies. 
    !         I strongly suggest to set emi_matrix values to 1. and only use either the sector 
    !         scalings OR the species scalings above.

    CALL message('','---')
    CALL print_value('Dry deposition velocity scaling factor (scale_drydep_acc)', scale_drydep_acc)

    CALL message('','---')
    CALL print_value('Wet deposition scaling factor, in-cloud (scale_wetdet_ic)', scale_wetdep_ic)
    CALL print_value('Wet deposition scaling factor, below-cloud (scale_wetdep_bc)', scale_wetdep_bc)
   ! CALL print_value('Wet deposition scaling factor, in-cloud BC (scale_wetdep_ic_bc_only)', scale_wetdep_ic_bc_only)
   ! CALL print_value('Wet deposition scaling factor, below-cloud BC (scale_wetdep_bc_bc_only)', scale_wetdep_bc_bc_only)

  !  CALL message('','---')
  !  CALL print_value('Entrainment of tracer mass flux scaling factor (scale_tr_entrainment)', scale_tr_entrainment)
    CALL print_value('Vertical velocity scaling factor (scale_vertical_velocity)', scale_vertical_velocity)

    CALL message('', '---')
    CALL print_value('BC imaginary refractive index (bc_rad_ni)', bc_rad_ni)
    CALL print_value('BC imaginary refractive index (oc_rad_ni)', oc_rad_ni)
    CALL print_value('BC imaginary refractive index (du_rad_ni)', du_rad_ni)

!    CALL message('', '---')
!    CALL print_value('The proportion of (G/F)FIRE tracer emitted into the layer &
 !   &imediately above the PBL (PBL+1) (prop_fire_in_pbl_p1)', prop_fire_in_pbl_p1)
 !   CALL print_value('The proportion of (G/F)FIRE tracer emitted into PBL+2 (prop_fire_in_pbl_p2)', prop_fire_in_pbl_p2)
  !  CALL print_value('The proportion of (G/F)FIRE tracer emitted into the PBL (prop_fire_in_pbl)', prop_fire_in_pbl)

    CALL message('', '---')
    CALL print_value('The scaling of SO4 layer thickness cutoff (scale_so4_coating)', scale_so4_coating)

!    CALL print_value('The scaling of intra-mode coagulation (scale_intra_mode_coagulation)', scale_intra_mode_coagulation)
 
!CALL print_value('The scaling of inter-mode coagulation (scale_inter_mode_coagulation)', scale_inter_mode_coagulation)

    CALL message('', '---')
 !   CALL print_value('The scaling of ARG activation', scale_activation)
 !   CALL message('hammoz_perturbations', 'Be careful - the below perturbations are only applied &
 !   &for the KK scheme (nauto=2)')
 !   CALL print_value('The KK exponent', KK_exponent)
 !   CALL print_value('The LWP KK exponent', KK_LWP_exponent)
 !   CALL print_value('The scaling of accretion', scale_accretion)
 !   CALL print_value('The aerosol water uptake', scale_RH)

    CALL message('', '---')
    CALL print_value('The scaling of the (AMIP) solar constant', scale_solar_const)
    ! Actually do the scaling here
    ssi_amip = ssi_amip * scale_solar_const

  END SUBROUTINE init_hammoz_perturbations


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! init_hammoz_emi_perturbations scales the default scaling factors for emissions
!! in the emission matrix ematrix
!!
!! @author Philip Stier (University of Oxford)
!!
!! @par Revision History
!! - Philip Stier - original code (04/2015)
!!                - revised (10/2016)
!!
!! @par This subroutine is called by
!! init_subm in mo_submodel_interface
!!
!! @par Externals:
!! <ol>
!! <li>em_get_SectorIndex
!! <li>em_add_spec_to_sector
!! </ol>
!!
!! @par Notes
!!
!!
!! @par Responsible coder
!! philip.stier@physics.ox.ac.uk
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE init_hammoz_emi_perturbations

    USE mo_emi_matrix,  ONLY: ematrix, em_get_SectorIndex, maxsectors, EM_FIRE
    USE mo_submodel,    ONLY: print_value
    USE mo_exception,   ONLY: message, message_text
    USE mo_util_string, ONLY: separator


    INTEGER :: i, ind, nvars, icat, isec, nsec


    !CHARACTER(LEN=64), DIMENSION(1) :: sectors_bf=(/ 'AWB' /)
   ! CHARACTER(LEN=64), DIMENSION(6) :: sectors_ff=(/ 'AIRC', 'DOM', 'ENE', 'IND', 'SHIPS', 'TRA' /)

!    CHARACTER(LEN=64), DIMENSION(6) :: sectors_ff ! YAB 29/04/24 Needed to make new sectors_ff string for compilier
!
!    sectors_ff(1) = 'AIRC'
!    sectors_ff(2) = 'DOM'
!    sectors_ff(3) = 'ENE'
!    sectors_ff(4) = 'IND'
!    sectors_ff(5) = 'SHIPS'
!    sectors_ff(6) = 'TRA'
!        
!    CHARACTER(LEN=64), DIMENSION(8) :: sectors_sulf!  YAB 25/06/24 New Anth SO2 emissions
!
!    sectors_sulf(1) = 'AIRC'
!    sectors_sulf(2) = 'DOM'
!    sectors_sulf(3) = 'ENE'
!    sectors_sulf(4) = 'IND'
!    sectors_sulf(5) = 'SHIPS'
!    sectors_sulf(6) = 'TRA'
!    sectors_sulf(7) = 'AWB'
!    sectors_sulf(8) = 'WST'
    CHARACTER(LEN=64), DIMENSION(1) :: sectors_bf
    CHARACTER(LEN=64), DIMENSION(6) :: sectors_ff
    CHARACTER(LEN=64), DIMENSION(8) :: sectors_so2

    ! Initialize arrays
    sectors_bf(1) = 'AWB'
    sectors_ff(1) = 'AIRC'
    sectors_ff(2) = 'DOM'
    sectors_ff(3) = 'ENE'
    sectors_ff(4) = 'IND'
    sectors_ff(5) = 'SHIPS'
    sectors_ff(6) = 'TRA'
    sectors_so2(1) = 'AIRC'
    sectors_so2(2) = 'DOM'
    sectors_so2(3) = 'ENE'
    sectors_so2(4) = 'IND'
    sectors_so2(5) = 'SHIPS'
    sectors_so2(6) = 'TRA'
    sectors_so2(7) = 'AWB'
    sectors_so2(8) = 'WST'


    CALL message('',separator)

    !--- Left emission sectors: AGR, BIOGENIC SLV TERR DUST OCEANI SEASALT VOLCE VOLCC
    !--- Left emission sectors: AGR, BIOGENIC SLV TERR DUST OCEANI SEASALT VOLCE VOLCC

    !--- 1) Map broad sectors to HAMMOZ sectors:

    !--- Biofuel:

    CALL message('',separator)

    nsec=SIZE(sectors_bf)
    DO isec=1, nsec
      ind=em_get_sectorindex(TRIM(sectors_bf(isec)))
      nvars=ematrix%em_sectors(ind)%es_nvars
      ematrix%em_sectors(ind)%es_variables(1:nvars)%ev_factor=ematrix%em_sectors(ind)%es_variables(1:nvars)%ev_factor*scale_emi_bf
      CALL message('hammoz_perturbations','BF Emissions scale factor for sector:')
      CALL message('', ematrix%em_sectors(ind)%es_sectorname )
      DO i=1, nvars
        CALL print_value(ematrix%em_sectors(ind)%es_variables(i)%ev_varname, ematrix%em_sectors(ind)%es_variables(i)%ev_factor)
      ENDDO
    ENDDO

    !    --- Fossil fuel:

    nsec=SIZE(sectors_ff)
    DO isec=1, nsec
      ind=em_get_sectorindex(TRIM(sectors_ff(isec)))
      nvars=ematrix%em_sectors(ind)%es_nvars
      ematrix%em_sectors(ind)%es_variables(1:nvars)%ev_factor=ematrix%em_sectors(ind)%es_variables(1:nvars)%ev_factor*scale_emi_ff
      CALL message('hammoz_perturbations','FF Emissions scale factor for sector:')
      CALL message('', ematrix%em_sectors(ind)%es_sectorname )
      DO i=1, nvars
        CALL print_value(ematrix%em_sectors(ind)%es_variables(i)%ev_varname, ematrix%em_sectors(ind)%es_variables(i)%ev_factor)
      ENDDO
    ENDDO

!    --- Biomass burning:

    DO ind=1, maxsectors
    !Find any sectors that are type FIRE since it can change for different emissions datasets
      IF (ematrix%em_sectors(ind)%es_emtype == EM_FIRE) THEN
        nvars=ematrix%em_sectors(ind)%es_nvars
        ematrix%em_sectors(ind)%es_variables(1:nvars)%ev_factor=ematrix%em_sectors(ind)%es_variables(1:nvars)%ev_factor*scale_emi_bb
        CALL message('hammoz_perturbations','Bio Emissions scale factor for sector:')
        CALL message('', ematrix%em_sectors(ind)%es_sectorname )
        DO i=1, nvars
          CALL print_value(ematrix%em_sectors(ind)%es_variables(i)%ev_varname, ematrix%em_sectors(ind)%es_variables(i)%ev_factor)
        ENDDO
      END IF
    ENDDO

!    --- DMS:

!    DO ind=1, maxsectors
!     !  Find any sectors that are type DMS since it can change for different emissions datasets
!      nvars=ematrix%em_sectors(ind)%es_nvars
!      DO i=1, nvars
!        IF (TRIM(ematrix%em_sectors(ind)%es_variables(i)%ev_varname) == "DMS") THEN
!          ematrix%em_sectors(ind)%es_variables(i)%ev_factor=ematrix%em_sectors(ind)%es_variables(i)%ev_factor*scale_emi_dms
!          CALL message('hammoz_perturbations','DMS Emissions scale factor for sector:')
!          CALL message('', ematrix%em_sectors(ind)%es_sectorname )
!          CALL print_value(ematrix%em_sectors(ind)%es_variables(i)%ev_varname, ematrix%em_sectors(ind)%es_variables(i)%ev_factor)
!        END IF
!      ENDDO
!    ENDDO

!    --- DUST:

    DO ind=1, maxsectors
!       Find any sectors that are type DUST since it can change for different emissions datasets
      nvars=ematrix%em_sectors(ind)%es_nvars
      DO i=1, nvars
        IF (TRIM(ematrix%em_sectors(ind)%es_variables(i)%ev_varname) == "DU") THEN
          ematrix%em_sectors(ind)%es_variables(i)%ev_factor=ematrix%em_sectors(ind)%es_variables(i)%ev_factor*scale_emi_du
          CALL message('hammoz_perturbations','DUST Emissions scale factor for sector:')
          CALL message('', ematrix%em_sectors(ind)%es_sectorname )
          CALL print_value(ematrix%em_sectors(ind)%es_variables(i)%ev_varname, ematrix%em_sectors(ind)%es_variables(i)%ev_factor)
        END IF
      ENDDO
    ENDDO

!    -- Anthropogenic SO2: YAB

    nsec=SIZE(sectors_so2)
    DO isec=1, nsec
      ind=em_get_sectorindex(TRIM(sectors_so2(isec)))
      nvars=ematrix%em_sectors(ind)%es_nvars
      ematrix%em_sectors(ind)%es_variables(1:nvars)%ev_factor=ematrix%em_sectors(ind)%es_variables(1:nvars)%ev_factor*scale_emi_so2
      CALL message('hammoz_perturbations','FF Emissions scale factor for sector:')
      CALL message('', ematrix%em_sectors(ind)%es_sectorname )
      DO i=1, nvars
        CALL print_value(ematrix%em_sectors(ind)%es_variables(i)%ev_varname, ematrix%em_sectors(ind)%es_variables(i)%ev_factor)
      ENDDO
    ENDDO

   !    --- Sea salt Aerosol:

    DO ind=1, maxsectors
!       Find any sectors that are type SSA since it can change for different emissions datasets
      nvars=ematrix%em_sectors(ind)%es_nvars
      DO i=1, nvars
        IF (TRIM(ematrix%em_sectors(ind)%es_variables(i)%ev_varname) == "SS") THEN
          ematrix%em_sectors(ind)%es_variables(i)%ev_factor=ematrix%em_sectors(ind)%es_variables(i)%ev_factor*scale_emi_ssa
          CALL message('hammoz_perturbations','SEASALT Emissions scale factor for sector:')
          CALL message('', ematrix%em_sectors(ind)%es_sectorname )
          CALL print_value(ematrix%em_sectors(ind)%es_variables(i)%ev_varname, ematrix%em_sectors(ind)%es_variables(i)%ev_factor)
        END IF
      ENDDO
    ENDDO


    CALL message('',separator)

  END SUBROUTINE init_hammoz_emi_perturbations

END MODULE mo_hammoz_perturbations
