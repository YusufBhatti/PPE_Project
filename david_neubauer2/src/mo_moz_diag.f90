!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_diag
!!
!! \brief
!!  This module defines the diagnostic stream and pointers for MOZ diagnostics
!!
!! \author Martin G. Schultz (FZ-Juelich)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - draft ECHAM6-HAMMOZ version (m.schultz, 2010-05-03)
!!  - cleanup and addition of diagnostics (m.schultz, 2014-03-19)
!!
!! \limitations
!!  This module contains a lot of information that is specific for the chemical mechanism
!!  of the model. If you simplify the mechanism, it should still work (no guarantee!), but
!!  if you extend the mechanism or change speciesd or reaction names, you will have to
!!  edit the code in this module.
!!  Be aware of:
!!  - species names (names_... variables)
!!  - reaction tag names (used to build the ndx... fields)
!!  - number of reactions for a given process (only partly flexible - don't forget to adapt 
!!    this also in subroutine moz_rates_diag!)
!!
!! \details
!!  Chemistry specific diagnostics include:
!!  - OH concentration (in molec cm-3), also as "true mean value"
!!  - chemical lifetimes of methane (CH4) and methyl chloroform (CH3CCl3)
!!  - surface area densities of stratospheric aerosol (sulfate, NAT, ice)
!!  - condensed HNO3 (stratosphere)
!!  - UV albedo
!!  - species families (NOy, Cly, Bry)
!!  - various reaction rates (turnovers)
!!
!! \bibliographic_references
!!  none
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_moz_diag

!! special diagnostics for MOZART chemistry

!! @par Current responsible coder
!! m.schultz@fz-juelich.de
!!
!## TODO:  - add local time diagnostics, aircraft flight track diag -- fix diags.

  USE mo_kind,               ONLY: dp
  USE mo_submodel_diag,      ONLY: vmem2d, vmem3d
  USE mo_moz,                ONLY: tracnam
  USE mo_moz_mods,           ONLY: pcnstm1, rxt_tag_lst     ! n) umber of MOZART tracers and reaction tags
  USE mo_time_event,         ONLY: io_time_event, TRIG_FIRST, TIME_INC_HOURS !s.stadtler

  IMPLICIT NONE

  PRIVATE

  ! -- public module routines
  PUBLIC  :: init_moz_streams, moz_family_diag, moz_rates_diag, mozh_diag

  ! -- public module variables: diagnostic pointers
  REAL(dp), PUBLIC, POINTER :: dpalb(:,:) => NULL()           ! UV albedo                                          (mo_moz_uvalbedo)
  REAL(dp), PUBLIC, POINTER :: dpohconc(:,:,:) => NULL()      ! OH in concentration units                          (mo_moz_chemdr)
  REAL(dp), PUBLIC, POINTER :: dpratech4(:,:,:) => NULL()     ! Methane rate (for calculating lifetime )           (mo_moz_chemdr)
  REAL(dp), PUBLIC, POINTER :: dpratemcl(:,:,:) => NULL()     ! Methyl chloroform rate (for calculating lifetime)  (mo_moz_chemdr)
  REAL(dp), PUBLIC, POINTER :: dpsad1(:,:,:) => NULL()        ! Stratospheric aerosol surface area density/sulfate (mo_moz_chemdr)
  REAL(dp), PUBLIC, POINTER :: dpsad2(:,:,:) => NULL()        ! Stratospheric aerosol surface area density/NAT     (mo_moz_chemdr)
  REAL(dp), PUBLIC, POINTER :: dpsad3(:,:,:) => NULL()        ! Stratospheric aerosol surface area density/ice     (mo_moz_chemdr)
  REAL(dp), PUBLIC, POINTER :: dphno3_cond(:,:,:) => NULL()   ! Stratospheric condensed nitric acid                (mo_moz_chemdr)

  REAL(dp), PUBLIC, POINTER :: dpnoy(:,:,:) => NULL()         ! (Inorganic) NOy family
                                                              ! NOy = N + NO + NO2 + NO3 + 2*N2O5 + HNO3 + HO2NO2 
                                                              !       + BrONO2 + ClONO2
  REAL(dp), PUBLIC, POINTER :: dpcly(:,:,:) => NULL()         ! (Inorganic) Cly family
                                                              ! Cly = HCl + ClONO2 + HOCl + ClO + Cl + 2*Cl2O2 + 2*Cl2 
                                                              !       + OClO + BrCl
  REAL(dp), PUBLIC, POINTER :: dpbry(:,:,:) => NULL()         ! (Inorganic) Bry family
                                                              ! Bry = Br + BrO + HOBr + HBr + BrONO2 + BrCl
  REAL(dp), PUBLIC, POINTER :: dpro2(:,:,:) => NULL()         ! Sum of organic peroxy radicals
                                                              ! RO2 = CH3O2 + ...
  REAL(dp), PUBLIC, POINTER :: dprooh(:,:,:) => NULL()        ! Sum of organic hydrogen peroxides
                                                              ! ROOH = CH3OOH + ...
  REAL(dp), PUBLIC, POINTER :: dprcoo2(:,:,:) => NULL()       ! RC(O)O2 = CH3CO3 + HOCH2CO3 + HCOCO3 + MCO3 + CO2H3CO3 + LIECO3 + LHC4ACCO3

  REAL(dp), PUBLIC, POINTER :: dppo3_ho2(:,:,:) => NULL()     ! Rate of ozone production due to NO+HO2
  REAL(dp), PUBLIC, POINTER :: dppo3_ch3o2(:,:,:) => NULL()   ! Rate of ozone production due to NO+CH3O2
  REAL(dp), PUBLIC, POINTER :: dppo3_ro2(:,:,:) => NULL()     ! Rate of ozone production due to NO+RO2 (incl. CH3O2!)
  REAL(dp), PUBLIC, POINTER :: dppo3_hno3(:,:,:) => NULL()    ! Rate of Ox production due to HNO3+OH
  REAL(dp), PUBLIC, POINTER :: dppo3_hono(:,:,:) => NULL()    ! Rate of Ox production due to HONO+OH
  REAL(dp), PUBLIC, POINTER :: dppo3_pan(:,:,:) => NULL()     ! Rate of Ox production due to PAN+OH
  REAL(dp), PUBLIC, POINTER :: dppo3_jpan(:,:,:) => NULL()    ! Rate of Ox production due to PAN+hv
  REAL(dp), PUBLIC, POINTER :: dppo3_janit(:,:,:) => NULL()   ! Rate of Ox production due to AN+hv
  REAL(dp), PUBLIC, POINTER :: dppo3_jho2no2(:,:,:) => NULL() ! Rate of Ox production due to HO2NO2+hv
  REAL(dp), PUBLIC, POINTER :: dplo3_ho2(:,:,:) => NULL()     ! Rate of ozone loss due to O3+HO2
  REAL(dp), PUBLIC, POINTER :: dplo3_oh(:,:,:) => NULL()      ! Rate of ozone loss due to O3+OH
  REAL(dp), PUBLIC, POINTER :: dplo3_o1d(:,:,:) => NULL()     ! Rate of ozone loss due to O1D+H2O (after O3+hv)
  REAL(dp), PUBLIC, POINTER :: dplo3_ene(:,:,:) => NULL()     ! Rate of ozone loss due to alkenes
  REAL(dp), PUBLIC, POINTER :: dplo3_c6h5o(:,:,:) => NULL()   ! Rate of ozone loss due to reaction of O3 with phenolate
  REAL(dp), PUBLIC, POINTER :: dplo3_c6h5ono2(:,:,:) => NULL()! Rate of ozone loss due to reaction of NO2 with phenolate
  REAL(dp), PUBLIC, POINTER :: dplo3_no(:,:,:) => NULL()   ! Rate of ozone loss due to reaction of O3 with NO
  REAL(dp), PUBLIC, POINTER :: dplo3_no2(:,:,:) => NULL()   ! Rate of ozone loss due to reaction of O3 with NO2
  REAL(dp), PUBLIC, POINTER :: dplox_het(:,:,:) => NULL()     ! Rate of Ox loss due to heterogeneous losses
  REAL(dp), PUBLIC, POINTER :: dplno3_jno3a(:,:,:) => NULL()   ! Rate of Ox production due to NO3 photolysis yielding NO2 + O
  REAL(dp), PUBLIC, POINTER :: dplno3_jno3b(:,:,:) => NULL()    ! Rate of Ox loss due to NO3 photolysis
  REAL(dp), PUBLIC, POINTER :: dplno3_ho2(:,:,:) => NULL()     ! Rate of Ox loss due to NO3+HO2
  REAL(dp), PUBLIC, POINTER :: dplno3_rcho(:,:,:) => NULL()    ! Rate of Ox loss due to NO3+RCHO and phenols (PHENOL, CATECHOL, XYLOL, CRESOL)
  REAL(dp), PUBLIC, POINTER :: dplno3_rene(:,:,:) => NULL()    ! Rate of Ox loss due to NO3+RENE (alkenes)
  REAL(dp), PUBLIC, POINTER :: dplno3_ro2no3(:,:,:) => NULL()  ! Rate of Ox loss due to NO3+RO2
  REAL(dp), PUBLIC, POINTER :: dplo3_xo(:,:,:) => NULL()      ! Rate of Ox loss due to XO+HO2
  REAL(dp), PUBLIC, POINTER :: dpo3netchem(:,:,:) => NULL()   ! Rate of ozone netto change due to chemical processes
  REAL(dp), PUBLIC, POINTER :: dprco_oh(:,:,:) => NULL()      ! Rate of CO loss due to CO+OH
  REAL(dp), PUBLIC, POINTER :: dprch4_oh(:,:,:) => NULL()     ! Rate of CH4 loss due to CH4+OH
  REAL(dp), PUBLIC, POINTER :: dplch4(:,:,:) => NULL()        ! Rate of CH4 loss due to all reactions (OH, O1D, photolysis)
  REAL(dp), PUBLIC, POINTER :: dprno2_hv(:,:,:) => NULL()     ! Rate of photolysis reaction of NO2
  REAL(dp), PUBLIC, POINTER :: dprno2_oh(:,:,:) => NULL()     ! Rate of NO2+OH
  REAL(dp), PUBLIC, POINTER :: dprno2_ohreact(:,:,:) => NULL()! OH-reactivity due to NO2+OH
  REAL(dp), PUBLIC, POINTER :: dprho2_ho2(:,:,:) => NULL()    ! Rate of reaction HO2+HO2 -> H2O2
  REAL(dp), PUBLIC, POINTER :: dplpan(:,:,:) => NULL()        ! Rate of PAN loss due to thermal decomposition, hv, and OH
  REAL(dp), PUBLIC, POINTER :: dploh(:,:,:) => NULL()         ! Rate of OH loss due to all reactions
  REAL(dp), PUBLIC, POINTER :: dpohreact(:,:,:) => NULL()     ! OH reactivity
  REAL(dp), PUBLIC, POINTER :: dplrcoo2_no2(:,:,:) => NULL()  ! Rate of RC(O)O2 loss due to NO2
!!++mgs20140320: temporary addition - possibly irrelevant later
  REAL(dp), PUBLIC, POINTER :: dpln2o_j(:,:,:) => NULL()      ! Rate of N2O loss due to photolysis
  REAL(dp), PUBLIC, POINTER :: dpln2o_o1d(:,:,:) => NULL()    ! Rate of N2O loss due to reaction with O1D
!!--mgs
!! s.stadtler
  REAL(dp), PUBLIC, POINTER :: dpk_n2o5(:,:,:) => NULL()      ! Ratecoefficient of N2O5 loss due to heterogeneous reaction on aerosol
  REAL(dp), PUBLIC, POINTER :: dpk_no3(:,:,:) => NULL()       ! Ratecoefficient of NO3  loss due to heterogeneous reaction on aerosol
  REAL(dp), PUBLIC, POINTER :: dpk_no2(:,:,:) => NULL()       ! Ratecoefficient of NO2  loss due to heterogeneous reaction on aerosol
  REAL(dp), PUBLIC, POINTER :: dpk_ho2(:,:,:) => NULL()       ! Ratecoefficient of HO2  loss due to heterogeneous reaction on aerosol
  REAL(dp), PUBLIC, POINTER :: dpk_hno3(:,:,:) => NULL()      ! Ratecoefficient of HNO3  loss due to heterogeneous reaction on aerosol
  REAL(dp), PUBLIC, POINTER :: dpk_o3(:,:,:) => NULL()        ! Ratecoefficient of O3  loss due to heterogeneous reaction on aerosol
  REAL(dp), PUBLIC, POINTER :: dpgamma_n2o5(:,:,:) => NULL()  ! Reactionprobability of N2O5 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpgamma_n2o5_a(:,:,:) => NULL()! Reactionprobability of N2O5 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpgamma_ho2(:,:,:) => NULL()   ! Reactionprobability of HO2 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpgamma_ho2_a(:,:,:) => NULL() ! Reactionprobability of HO2 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpgamma_ho2_cloud(:,:,:) => NULL()   ! Reactionprobability of HO2 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpgamma_ho2_tot(:,:,:) => NULL()   ! Overall reactionprobability of HO2 on aerosols and clouds
  REAL(dp), PUBLIC, POINTER :: dpgamma_hno3(:,:,:) => NULL()   ! Reactionprobability of HNO3 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpgamma_hno3_a(:,:,:) => NULL() ! Reactionprobability of HNO3 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpgamma_o3(:,:,:) => NULL()   ! Reactionprobability of O3 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpgamma_o3_a(:,:,:) => NULL() ! Reactionprobability of O3 on aerosols
  REAL(dp), PUBLIC, POINTER :: dpsa_aerosol(:,:,:) => NULL()  ! surfaceareadensity of all aerosols
  REAL(dp), PUBLIC, POINTER :: dpsa_aerosol_wet(:,:,:) => NULL()  ! surface area density of wet aerosols
  REAL(dp), PUBLIC, POINTER :: dpsa_cloud(:,:,:) => NULL()    ! surface area density of cloud droplets
  REAL(dp), PUBLIC, POINTER :: dpppan(:,:,:) => NULL()         ! Rate of PAN production
  REAL(dp), PUBLIC, POINTER :: dplghno3(:,:,:) => NULL()      ! Rate of HNO3 loss via gasphase
  REAL(dp), PUBLIC, POINTER :: dpn2o5(:,:,:) => NULL()        ! Rate of N2O5 Photolysis
  REAL(dp), PUBLIC, POINTER :: dplgn2o5(:,:,:) => NULL()      ! Rate of N2O5 thermal loss
  REAL(dp), PUBLIC, POINTER :: dpno2_no3(:,:,:) => NULL()     ! Rate of N2O5 production
  REAL(dp), PUBLIC, POINTER :: dppno3(:,:,:) => NULL()        ! Rate of NO3 production
  REAL(dp), PUBLIC, POINTER :: dplgno3(:,:,:) => NULL()       ! Rate of NO3 loss via gasphase
  REAL(dp), PUBLIC, POINTER :: dplnox_het(:,:,:) => NULL()       ! Rate of NOx heterogeneous loss, rate of HNO3 heterogeneous production
  REAL(dp), PUBLIC, POINTER :: dppnox(:,:,:) => NULL()        ! Rate of NOx production
  REAL(dp), PUBLIC, POINTER :: dplnox(:,:,:) => NULL()        ! Rate of NOx loss
  REAL(dp), PUBLIC, POINTER :: dpvmro3(:,:) => NULL()         ! Hourly Ozone mixing ratio
  REAL(dp), PUBLIC, POINTER :: dpvmrno(:,:) => NULL()         ! Hourly NO mixing ratio
  REAL(dp), PUBLIC, POINTER :: dpvmrno2(:,:) => NULL()        ! Hourly NO2 mixing ratio
  REAL(dp), PUBLIC, POINTER :: dpvmrco(:,:) => NULL()         ! Hourly CO mixing ratio
!! s.stadtler

!! templates for pointer definitions...
!!  REAL(dp), PUBLIC, POINTER :: diagfield(:,:,:)       
!!  TYPE (vmem3d), PUBLIC, ALLOCATABLE :: diaglist(:)  
!!

  ! -- private module variables
!## TODO: let user control detail of MOZ diagnostics via MOZCTL namelist
  LOGICAL                   :: ldetail                     ! turn on detailed diagnostics
  LOGICAL                   :: ldetail2                    ! turn on detailed diagnostics for Scarlet's work
  TYPE(io_time_event), SAVE :: mozh_putdata                ! Output interval for hourly moz stream

  INTEGER, PARAMETER        :: nmaxtpf = 200               ! maximum number of tracers per family
  INTEGER, PARAMETER        :: ntracene = 11               ! number of alkenes reacting with O3
  INTEGER, PARAMETER        :: ntracrene = 8               ! number of alkenes reacting with O3
  INTEGER, PARAMETER        :: ntracrcho = 15              ! number of aldehydes and phenols reacting with NO3
  INTEGER, PARAMETER        :: ntracphenoxy = 2            ! number of phenoxy radicals reacting with O3
  INTEGER, PARAMETER        :: ntracjanit = 14             ! number of alkyl nitrates yielding NO2 by photolysis
  INTEGER                   :: ntracnoy, ntraccly, ntracbry, ntracro2, ntracrooh, ntracrcoo2, ntracoh, ntracohreact  ! number of tracers in families
  INTEGER                   :: nreacch4                    ! number of reactions with CH4
  INTEGER                   :: nreacpan                    ! number of reactions with PAN
  INTEGER                   :: nreacno3                    ! number of reactions with NO3 
  INTEGER                   :: nreachno3                   ! number of reactions with HNO3 
  INTEGER                   :: nreacn2o5                   ! number of reactions with N2O5
  INTEGER                   :: nreacpno3                   ! number of reactions producing NO3 
  INTEGER                   :: nreacpnox                   ! number of reactions producing NOx 
  INTEGER                   :: nreaclnox                   ! number of reactions leading to loss of NOx 
  INTEGER                   :: nreacxo                     ! number of halogen reactions (Ox) 
  INTEGER                   :: ndxnoy(nmaxtpf), ndxcly(nmaxtpf), ndxbry(nmaxtpf), &
                               ndxro2(nmaxtpf), ndxrooh(nmaxtpf), ndxrcoo2(nmaxtpf), &
                               ndxmozh(nmaxtpf)
  INTEGER                   :: ntracmozh                   !number of tracers in moz hourly output
  ! indices for reaction rates (turnovers) 
  ! These are 3-tuples, where the first element is the index for the rate coefficient (rxt_tag_lst)
  ! and the second and third contain the species indices. If there is only one species, the third
  ! element should be -1.
  ! The first set of variables is for "single reaction" diagnostics, the second set for cases where several 
  ! reactions must be summed.
  INTEGER, DIMENSION(3)     :: ndxpo3_ho2, ndxpo3_ch3o2, ndxpo3_hno3, ndxpo3_pan,    &  ! Ox production
                               ndxpo3_jpan, ndxpo3_jho2no2, ndxlhhno3, ndxlno3_jno3a,&  ! Ox production
                               ndxpo3_hono,                                          &  ! Ox production
                               ndxlo3_ho2, ndxlo3_oh, ndxlno3_jno3b,                 &  ! Ox loss
                               ndxlo3_o1d, ndxlno3_ho2,                              &  ! Ox loss 
                               ndxlo3_no, ndxlo3_no2,                                &  ! O3 loss 
                               ndxch4_oh,                                            &  ! CH4+OH
                               ndxno2_oh,                                            &  ! HNO3 formation from NO2+OH
                               ndxno2_ohreact,                                       &  ! OH reactivity due NO2+OH
                               ndxno2_hv,                                            &  ! NO2 Photolysis
                               ndxho2_ho2,                                           &  ! HO2+HO2
                               ndxppan,                                              &  ! PAN production
                               ndxlgn2o5,                                            &  ! N2O5 thermal loss
                               ndxno2_no3,                                           &  ! N2O5 production due to NO2 and NO3
                               ndxph2so4,                                            &  ! H2SO4 production
                               ndxpso2dmsno3                                            ! DMS reaction with NO3


!!++mgs20140320: temporary addition - possibly irrelevant later
  INTEGER, DIMENSION(3)     :: ndxn2o_j                                                 ! N2O loss reactions
  INTEGER, DIMENSION(2,3)   :: ndxn2o_o1d                                               ! N2O loss reactions
!!--mgs

!!###PRELIMINARY - NEED TO FIGURE OUT DIMENSIONS...
  INTEGER, DIMENSION(nmaxtpf,3) :: ndxpo3_ro2,                                       &  ! ozone production due to RO2
                               ndxpo3_janit,                                         &  ! ozone production due to AN + hv
                               ndxlo3_ene,                                           &  ! ozone loss due to alkenes
                               ndxlno3_rcho,                                         &  ! ozone loss due to aldehydes and phenols 
                               ndxlno3_ro2no3,                                       &  ! ozone loss due to RO2+NO3 
                               ndxlno3_rene,                                         &  ! ozone loss due to NO3+alkene
                               ndxlo3_c6h5ono2,                                      &  ! Ox loss due to phenoxy + NO2  
                               ndxlo3_c6h5o,                                         &  ! Ox loss due to phenoxy + O3  
                               ndxlo3_xo,                                            &  ! Ox loss due to ClO and BrO 
                               ndxlox_het,                                           &  ! Ox loss due to heterogeneous reactions 
                               ndxco_oh,                                             &  ! CO+OH
                               ndxlch4,                                              &  ! CH4 loss
                               ndxlpan,                                              &  ! PAN loss
                               ndxloh,                                               &  ! OH loss
                               ndxohreact,                                           &  ! OH reactivity
                               ndxlrcoo2_no2,                                        &  ! loss RCOO2 due to NO2
                               ndxn2o5,                                              &  ! N2O5 Photolysis
                               ndxpno3,                                              &  ! NO3 production
                               ndxlgno3,                                             &  ! NO3 loss via gasphase
                               ndxlnox_het,                                          &  ! NOx heterogeneous loss via HNO3
                               ndxpnox,                                              &  ! NOx production
                               ndxlnox,                                              &  ! NOx loss
                               ndxlghno3                                          ! HNO3 loss via gasphase
  CONTAINS


  SUBROUTINE init_moz_streams

    USE mo_exception,           ONLY: finish
    USE mo_linked_list,         ONLY: t_stream, SURFACE, HYBRID
    USE mo_filename,            ONLY: trac_filetype
    USE mo_memory_base,         ONLY: AUTO, new_stream, get_stream,                  &
                                      default_stream_setting,                        &
                                      add_stream_reference, add_stream_element
    USE mo_time_control,        ONLY: putdata
    USE mo_moz,                 ONLY: lphotolysis
    USE mo_moz_photo,           ONLY: init_photo_stream
    USE mo_moz_util,            ONLY: get_spc_ndx, get_rxt_ndx
    USE mo_moz_mods,            ONLY: pcnstm1
    USE mo_util_string,         ONLY: tolower

    TYPE (t_stream), POINTER   :: stream_moz, stream_mozh, stream_tracer

    INTEGER                    :: i, jt

    ! family member names: these will be compared to the tracnam definitions in mo_moz_subs.
    ! the number of names must be adapted if new species are added to the mechanism. 
    ! too many names are not a problem, they will be ignored.
    ! if a species is counted twice, its names must appear twice.
    INTEGER, PARAMETER         :: ndefnoy = 10, ndefcly = 11, ndefbry = 8, ndefro2 = 44, ndefrooh = 28,  &
                                  ndefrcoo2 = 9, ndefmozh = 4
    CHARACTER(len=24)          :: names_noy(ndefnoy), names_cly(ndefcly), names_bry(ndefbry),    &
                                  names_ro2(ndefro2), names_rooh(ndefrooh), names_rcoo2(ndefrcoo2), &
                                  names_mozh(ndefmozh) ! s.stadtler
    CHARACTER(len=12), PARAMETER  :: names_ene(11)     = (/'C2H4        ', 'C3H6        ', 'MACR        ', 'MVK         ',      &
                                                           'C5H8        ', 'MBO         ', 'APIN        ', 'BPIN        ',      &
                                                           'LIMON       ', 'MYRC        ', 'BCARY       '   /)                 ,&
                                    names_rene(8)     = (/ 'C3H6        ', 'C5H8        ', 'MBO         ', 'APIN        ',      &
                                                           'BPIN        ', 'LIMON       ', 'MYRC        ', 'BCARY       ' /)   ,&
                                    names_ene_no3(13) = (/ 'CH2O        ', 'CH3CHO      ', 'C3H6        ', 'MACR        ',      &
                                                           'GLYALD      ', 'CH3COCHO    ',                                      &
                                                           'C5H8        ', 'MBO         ', 'APIN        ', 'BPIN        ',      &
                                                           'LIMON       ', 'MYRC        ', 'BCARY       ' /)                   ,&
                                    names_reac_no3(3) = (/ 'OH          ', 'HO2         ', 'NO          ' /)                   ,&
                                    names_rcho(15)    = (/ 'HCHO        ', 'CH3CHO      ', 'GLYOXAL     ', 'GLYALD      ',      &
                                                           'CH3COCHO    ', 'MACR        ', 'CO2H3CHO    ', 'LIECHO      ',      &
                                                           'LHC4ACCHO   ', 'TERPROD1    ', 'NC4CHO      ', 'PHENOL      ',      &
                                                           'CRESOL      ', 'CATECHOL    ', 'XYLOL       '/)                    ,&
                                    names_phenoxy(2)  = (/ 'C6H5O       ', 'CATEC1O     ' /)                                   ,&
                                    names_janit(14)    = (/ 'NOA         ', 'MEKNO3      ', 'ALKNO3      ', 'MVKN        ',     &
                                                            'MACRN       ', 'NC4CHO      ', 'ISOPDNO3    ', 'ISOPBNO3    ',     &
                                                            'LISOPACNO3  ', 'C59OOH      ', 'LISOPNO3OOH ', 'LISOPNO3NO3 ',     &
                                                            'TERPNO3     ', 'NTERPNO3    ' /)

    CHARACTER(len=12) ::  namestemp
 
    ldetail  = .TRUE.                 ! detailed diagnostics
    ldetail2 = .FALSE.               ! detailed diagnostics for Scarlet's work: #temporary

    ! Initialize independent streams
    IF (lphotolysis) CALL init_photo_stream
    ! Note: lightning can run independently of lmoz and therefore manages its own stream

    ! Define moz diagnostic stream
    CALL new_stream(stream_moz,'moz',filetype=trac_filetype,interval=putdata,lrerun=.FALSE.)
    CALL default_stream_setting (stream_moz,               &
                                 lrerun    = .FALSE.,      &
                                 laccu     = .FALSE.,      &
                                 lpost     =  .TRUE.,      &
                                 leveltype =  SURFACE,     &
                                 table     =  199,         &
                                 code      =  AUTO         )

    CALL add_stream_reference (stream_moz, 'geosp'   ,'g3b'   )
    CALL add_stream_reference (stream_moz, 'aps'     ,'g3b'   )
    CALL add_stream_reference (stream_moz, 'gboxarea','geoloc')

    ! additional tracers (species families)
    ! Add instantaneous diagnostics (output depends on ldetail)
!+++sschr: 20140321: add tracers to moz stream (and not to tracer stream)
!+++                 (for diagnostic reasons theses tracers should go to output as "vmr"
!+++                 -- that can't be recalculated if once written to the tracer stream
!+++                 as "mmr")

    CALL default_stream_setting (stream_moz, laccu = .FALSE., leveltype =  HYBRID)

    CALL add_stream_element (stream_moz,'NOy', dpnoy,                            &
                             longname='NOy volume mixing ratio',                    &
                             units='mole mole-1', lpost=ldetail)
    CALL add_stream_element (stream_moz,'Cly', dpcly,                            &
                             longname='Cly volume mixing ratio',                    &
                             units='mole mole-1', lpost=ldetail)
    CALL add_stream_element (stream_moz,'Bry', dpbry,                            &
                             longname='Bry volume mixing ratio',                    &
                             units='mole mole-1', lpost=ldetail)
    CALL add_stream_element (stream_moz,'RO2', dpro2,                            &
                             longname='Organic peroxy radicals volume mixing ratio',   &
                             units='mole mole-1', lpost=ldetail)
    CALL add_stream_element (stream_moz,'ROOH', dprooh,                          &
                             longname='Organic hydrogen peroxides volume mixing ratio',&
                             units='mole mole-1', lpost=ldetail)
    CALL add_stream_element (stream_moz,'RCOO2', dprcoo2,                          &
                             longname='Organic carbonate radicals volume mixing ratio',&
                             units='mole mole-1', lpost=ldetail)
!---sschr: 20140321: add tracers to moz stream

    ! OH related
    CALL add_stream_element (stream_moz,'OH_conc', dpohconc,                        &
                             longname='OH concentration',                           &
                             units='molecules cm-3', lpost=.true.)
    CALL add_stream_element (stream_moz, 'rate_CH4', dpratech4,                     &
                             units = 'years',                                       &
                             longname = 'CH4 rate by OH', lpost=.true.) !s.stadtler
    CALL add_stream_element (stream_moz, 'rate_CH3CCl3', dpratemcl,                 &
                             units = 'years',                                       &
                             longname = 'CH3CCl3 rate by OH', lpost=.true.) !s.stadtler
    ! Stratospheric aerosol
    CALL add_stream_element (stream_moz, 'sad_strat_sulfate', dpsad1,               &
                             units = 'cm2 cm-3',                                    &
                             longname = 'surface area density of stratospheric sulfate aerosol', lpost=ldetail)
    CALL add_stream_element (stream_moz, 'sad_strat_NAT', dpsad2,                   &
                             units = 'cm2 cm-3',                                    &
                             longname = 'surface area density of stratospheric nitric acid ternary aerosol', lpost=ldetail)
    CALL add_stream_element (stream_moz, 'sad_strat_ice', dpsad3,                   &
                             units = 'cm2 cm-3',                                    &
                             longname = 'surface area density of stratospheric ice aerosol', lpost=ldetail)
    CALL add_stream_element (stream_moz, 'hno3_cond', dphno3_cond,                  &
                             units = 'mole mole-1',                                 &
                             longname = 'volume mixing ratio of condensed nitric acid', lpost=ldetail)

     ! Chemical production and loss rates
     ! s.stadtler ozone budget production terms
     CALL add_stream_element (stream_moz, 'prodo3viaho2', dppo3_ho2,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox production due to NO+HO2', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viach3o2', dppo3_ch3o2,            &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox production due to NO+CH3O2', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viaro2', dppo3_ro2,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox production due to NO+RO2 (incl. CH3O2!)', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viahno3', dppo3_hno3,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox production due to HNO3+OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viapan', dppo3_pan,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox production due to PAN+OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viajpan', dppo3_jpan,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox production due to PAN+hv', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viajho2no2', dppo3_jho2no2,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox production due to HO2NO2+hv', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viajno3', dplno3_jno3a,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox loss due to NO3+hv', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viahono', dppo3_hono,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox production due to HONO+OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodo3viajanit', dppo3_janit,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of NO2 production due to AN+hv', lpost=ldetail)

     ! s.stadtler ozone budget loss terms
     CALL add_stream_element (stream_moz, 'losso3viaho2', dplo3_ho2,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of ozone loss due to O3+HO2', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viaoh', dplo3_oh,                  &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of ozone loss due to O3+OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viaalkenes', dplo3_ene,            &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of ozone loss due to alkenes', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viaphenoxy', dplo3_c6h5o,            &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of ozone loss due to C6H5O+O3', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'o3netchem', dpo3netchem,                 &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of ozone netto change due to chemical processes', lpost=.False.)
     CALL add_stream_element (stream_moz, 'losso1dviah2o', dplo3_o1d,               &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of ozone loss due to O1D+H2O (after O3+hv)', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossoxhetchem', dplox_het,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of Ox loss due to heterogeneous reactions', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viajno3', dplno3_jno3b,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox loss due to NO3+hv', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viano3', dplno3_ho2,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of ozone loss due to NO3+HO2', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viarcho', dplno3_rcho,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox loss due to NO3+RCHO and phenols', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossno3viarene', dplno3_rene,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox loss due to RENE+NO3', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossno3viaro2', dplno3_ro2no3,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox loss due to RO2+NO3 (incl. CH3O2)', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossno2viaphenoxy', dplo3_c6h5ono2,                  &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox loss due to C6H5O+NO2', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viaxo', dplo3_xo,                  &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of Ox loss due to XO+HO2 (X=Cl,Br)', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viano', dplo3_no,                  &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of O3 loss due to O3+NO', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losso3viano2', dplo3_no2,                  &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of O3 loss due to O3+NO2', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossco', dprco_oh,                       &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of CO loss due to CO+OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossch4', dprch4_oh,                     &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of CH4 loss due to CH4+OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'prodh2o2viaho2', dprho2_ho2,             &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of reaction HO2+HO2 -> H2O2', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossno2viaoh', dprno2_oh,                &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of NO2+OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'ohreactviano2', dprno2_ohreact,   &
                              units = 's-1',                            &
                              longname = 'OH-reactivity due to NO2+OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'losspan', dplpan,                        &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of PAN loss due to thermal decomposition, hv, and OH', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossoh', dploh,                          &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of OH loss due to all reactions', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'oh_reactivity', dpohreact,               &
                              units = 's-1',                                        &
                              longname = 'OH reactivity', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossrcoo2viano2', dplrcoo2_no2,          &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of RCOO2 loss due to NO2', lpost=ldetail)
!!++mgs20140320: temporary addition - possibly irrelevant later
     CALL add_stream_element (stream_moz, 'lossno2viaphoto', dprno2_hv,            &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of photolysis reaction of NO2', lpost=ldetail)
     CALL add_stream_element (stream_moz, 'lossn2oviaphoto', dpln2o_j,              &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of N2O loss due to photolysis', lpost=.False.)
     CALL add_stream_element (stream_moz, 'lossn2oviao1d', dpln2o_o1d,              &
                              units = 'mole mole-1 s-1',                            &
                              longname = 'Rate of N2O loss due to reaction with O1D', lpost=.False.)
!!--mgs
!! s.stadtler
     CALL add_stream_element (stream_moz, 'k_het_no3', dpk_no3,                     &
                              units = 's-1',                                        &
                              longname = 'Ratecoefficient of NO3 loss due to heterogeneous reaction on aerosols', lpost=.False.)
     CALL add_stream_element (stream_moz, 'k_het_no2', dpk_no2,                     &
                              units = 's-1',                                        &
                              longname = 'Ratecoefficient of NO2 loss due to heterogeneous reaction on aerosols', lpost=.False.)
     CALL add_stream_element (stream_moz, 'k_het_n2o5', dpk_n2o5,                   &
                              units = 's-1',                                        &
                              longname = 'Ratecoefficient of N2O5 loss due to heterogeneous reaction on aerosols', lpost=.False.)
     CALL add_stream_element (stream_moz, 'k_het_ho2', dpk_ho2,                   &
                              units = 's-1',                                        &
                              longname = 'Ratecoefficient of HO2 loss due to heterogeneous reaction on aerosols', lpost=.False.)
     CALL add_stream_element (stream_moz, 'k_het_hno3', dpk_hno3,                   &
                              units = 's-1',                                        &
                              longname = 'Ratecoefficient of HNO3 loss due to heterogeneous reaction on aerosols', lpost=.False.)
     CALL add_stream_element (stream_moz, 'k_het_o3', dpk_o3,                   &
                              units = 's-1',                                        &
                              longname = 'Ratecoefficient of O3 loss due to heterogeneous reaction on aerosols', lpost=.False.)
     CALL add_stream_element (stream_moz, 'gamma_n2o5', dpgamma_n2o5,               &
                              units = '1',                                          &
                              longname = 'Reactionprobability of N2O5 on aerosols weighted by S_a', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'gamma_n2o5_a', dpgamma_n2o5_a,            &
                              units = '1',                                          &
                              longname = 'Reactionprobability of N2O5 on aerosols calculated by k', lpost=.False.)
     CALL add_stream_element (stream_moz, 'gamma_ho2', dpgamma_ho2,               &
                              units = '1',                                          &
                              longname = 'Reactionprobability of HO2 on aerosols weighted by S_a_wet', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'gamma_ho2_a', dpgamma_ho2_a,            &
                              units = '1',                                          &
                              longname = 'Reactionprobability of HO2 on aerosols calculated by k', lpost=.False.)
     CALL add_stream_element (stream_moz, 'gamma_ho2_cloud', dpgamma_ho2_cloud,            &
                              units = '1',                                          &
                              longname = 'Reactionprobability of HO2 on cloud droplets', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'gamma_ho2_tot', dpgamma_ho2_tot,            &
                              units = '1',                                          &
                              longname = 'Overall reactionprobability of HO2 on cloud droplets and aerosols', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'gamma_hno3', dpgamma_hno3,               &
                              units = '1',                                          &
                              longname = 'Reactionprobability of HNO3 on aerosols weighted by S_a', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'gamma_hno3_a', dpgamma_hno3_a,            &
                              units = '1',                                          &
                              longname = 'Reactionprobability of HNO3 on aerosols calculated by k', lpost=.False.)
     CALL add_stream_element (stream_moz, 'gamma_o3', dpgamma_o3,               &
                              units = '1',                                          &
                              longname = 'Reactionprobability of O3 on aerosols weighted by S_a', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'gamma_o3_a', dpgamma_o3_a,            &
                              units = '1',                                          &
                              longname = 'Reactionprobability of O3 on aerosols calculated by k', lpost=.False.)
     CALL add_stream_element (stream_moz, 'S_a', dpsa_aerosol,                      &
                              units = 'cm2 cm-3',                                   &
                              longname = 'Total aerosol surface area density', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'S_a_wet', dpsa_aerosol_wet,                      &
                              units = 'cm2 cm-3',                                   &
                              longname = 'Total wet aerosol surface area density', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'S_a_cloud', dpsa_cloud,                      &
                              units = 'cm2 cm-3',                                   &
                              longname = 'Total cloud droplet surface area density', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'prodpan', dpppan,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of PAN production', lpost=.False.)
     CALL add_stream_element (stream_moz, 'losshno3', dplghno3,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of HNO3 loss via gasphase', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'photn2o5', dpn2o5,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of N2O5 photolysis', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'thlossn2o5', dplgn2o5,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of N2O5 thermal loss', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'prodn2o5', dpno2_no3,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of N2O5 production NO2 + NO3', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'prodno3', dppno3,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of NO3 production', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'lossno3', dplgno3,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of NO3 loss via gasphase', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'hetlossnox', dplnox_het,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of NOx heterogeneous loss', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'prodnox', dppnox,                      &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of NOx production', lpost=ldetail2)
     CALL add_stream_element (stream_moz, 'lossnox', dplnox,                               &
                              units = 'mole mole-1 s-1',                                   &
                              longname = 'Rate of NOx loss', lpost=ldetail2)

    CALL add_stream_element (stream_moz, 'uvalbedo',  dpalb,                               &
                             units='1',                                                    &
                             longname='mean UV albedo', lpost=ldetail)


    ! Define hourly moz diagnostic stream

    mozh_putdata = io_time_event(1, TIME_INC_HOURS, TRIG_FIRST, 0)

    !CALL new_stream(stream_mozh,'mozh',filetype=trac_filetype, &
    !                interval=mozh_putdata,lrerun=.FALSE.)
    CALL new_stream(stream_mozh,'mozh',filetype=trac_filetype, &
                    interval=mozh_putdata,lrerun=.FALSE.)
    CALL default_stream_setting (stream_mozh,              &
                                 lrerun    = .FALSE.,      &
                                 laccu     = .FALSE.,      &
                                 lpost     =  .TRUE.,      &
                                 leveltype =  SURFACE,     &
                                 table     =  199,         &
                                 code      =  AUTO         )
    CALL add_stream_reference (stream_mozh, 'geosp'   ,'g3b'   )
    CALL add_stream_reference (stream_mozh, 'aps'     ,'g3b'   )
    CALL add_stream_reference (stream_mozh, 'gboxarea','geoloc')

    CALL default_stream_setting (stream_mozh, laccu = .FALSE., leveltype=SURFACE)

    CALL add_stream_element (stream_mozh,'O3', dpvmro3,                         &
                             longname='O3 volume mixing ratio at the surface',  &
                             units='mole mole-1', lpost=ldetail2)
    CALL add_stream_element (stream_mozh,'NO', dpvmrno,                         &
                             longname='NO volume mixing ratio at the surface',  &
                             units='mole mole-1', lpost=ldetail2)
    CALL add_stream_element (stream_mozh,'NO2', dpvmrno2,                       &
                             longname='NO2 volume mixing ratio at the surface', &
                             units='mole mole-1', lpost=ldetail2)
    CALL add_stream_element (stream_mozh,'CO', dpvmrco,                       &
                             longname='CO volume mixing ratio at the surface', &
                             units='mole mole-1', lpost=ldetail2)
    ! Initialize MOZ species indices for hourly output at surface

    names_mozh(:) = (/ 'O3                      ', &
                       'NO                      ', &
                       'NO2                     ', &
                       'CO                      '/)

    ntracmozh = 0
    DO i = 1, ndefmozh
       jt = get_spc_ndx(names_mozh(i))
       IF (jt > 0) THEN
          ntracmozh = ntracmozh +1
          ndxmozh(ntracmozh) = jt
       END IF
    END DO

!! s.stadtler

    ! ---- Initialize MOZ species indices for families ----
    names_noy(:) = (/ 'N                       ', 'NO                      ', 'NO2                     ',  &
                      'NO3                     ', 'N2O5                    ', 'N2O5                    ',  &
                      'HNO3                    ', 'HO2NO2                  ', 'BRONO2                  ',  &
                      'CLONO2                  '  /)
  
    names_cly(:) = (/ 'HCL                     ', 'CLONO2                  ', 'HOCL                    ',  &
                      'CLO                     ', 'CL                      ', 'CL2O2                   ',  &
                      'CL2O2                   ', 'CL2                     ', 'CL2                     ',  &
                      'OCLO                    ', 'BRCL                    '  /)
  
    names_bry(:) = (/ 'BR                      ', 'BRO                     ', 'HOBR                    ',  &
                      'HBR                     ', 'BRONO2                  ', 'BRCL                    ',  &
                      'OCLO                    ', 'BRCL                    '  /)

    names_ro2(:) = (/ 'CH3O2                   ', 'HOCH2OO                 ', 'C2H5O2                  ',  &
                      'EO2                     ', 'CH3CO3                  ', 'HOCH2CO3                ',  &
                      'HCOCO3                  ', 'C3H7O2                  ', 'PO2                     ',  &
                      'PRONO3BO2               ', 'CH3COCH2O2              ', 'MEKO2                   ',  &
                      'ENEO2                   ', 'MCO3                    ', 'MACRO2                  ',  &
                      'LHMVKABO2               ', 'IBUTALOHO2              ', 'CO2H3CO3                ',  &
                      'MALO2                   ', 'ALKO2                   ', 'MBOO2                   ',  &
                      'ISOPBO2                 ', 'ISOPDO2                 ', 'LISOPACO2               ',  &
                      'NISOPO2                 ', 'LIECO3                  ', 'IEC1O2                  ',  &
                      'LNISO3                  ', 'LHC4ACCO3               ', 'LC578O2                 ',  &
                      'C59O2                   ', 'MDIALO2                 ', 'BENZO2                  ',  &
                      'PHENO2                  ', 'C6H5O2                  ', 'TOLO2                   ',  &
                      'BZOO                    ', 'ACBZO2                  ', 'XYLOLO2                 ',  &
                      'XYLENO2                 ', 'TERPO2                  ', 'NTERPO2                 ',  &
                      'TERP2O2                 ', 'CATEC1O2                '   /)

    names_rooh(:) = (/'CH3OOH                  ', 'C2H5OOH                 ', 'EOOH                    ',  &
                      'C3H7OOH                 ', 'MEKOOH                  ', 'POOH                    ',  &
                      'ROOH                    ', 'PR2O2HNO3               ', 'MACROOH                 ',  &
                      'LHMVKABOOH              ', 'ALKOOH                  ', 'MBOOOH                  ',  &
                      'LISOPACOOH              ', 'ISOPBOOH                ', 'ISOPDOOH                ',  &
                      'NISOPOOH                ', 'LNISOOH                 ', 'LC578OOH                ',  &
                      'C59OOH                  ', 'PHENOOH                 ', 'BENZOOH                 ',  &
                      'TOLOOH                  ', 'C6H5OOH                 ', 'BZOOH                   ',  &
                      'XYLOLOOH                ', 'XYLENOOH                ', 'TERPOOH                 ',  &
                      'TERP2OOH                '   /)

    names_rcoo2(:)= (/'CH3CO3                  ', 'HOCH2CO3                ', 'HCOCO3                  ',  &
                      'MALO2                   ', 'LIECO3                  ', 'LHC4ACCO3               ',  &
                      'MDIALO2                 ', 'MCO3                    ', 'CO2H3CO3                ' /)

    ntracnoy = 0
    DO i = 1, ndefnoy
       jt = get_spc_ndx(names_noy(i))
       IF (jt > 0) THEN
          ntracnoy = ntracnoy + 1
          ndxnoy(ntracnoy) = jt
       END IF
    END DO
    !##TODO: report on species found for NOy definition (and others, of course)

    ntraccly = 0
    DO i = 1, ndefcly
       jt = get_spc_ndx(names_cly(i))
       IF (jt > 0) THEN
          ntraccly = ntraccly + 1
          ndxcly(ntraccly) = jt
       END IF
    END DO
 
    ntracbry = 0
    DO i = 1, ndefbry
       jt = get_spc_ndx(names_bry(i))
       IF (jt > 0) THEN
          ntracbry = ntracbry + 1
          ndxbry(ntracbry) = jt
       END IF
    END DO
 
    ntracro2 = 0
    DO i = 1, ndefro2
       jt = get_spc_ndx(names_ro2(i))
       IF (jt > 0) THEN
          ntracro2 = ntracro2 + 1
          ndxro2(ntracro2) = jt
       END IF
    END DO
 
    ntracrooh = 0
    DO i = 1, ndefrooh
       jt = get_spc_ndx(names_rooh(i))
       IF (jt > 0) THEN
          ntracrooh = ntracrooh + 1
          ndxrooh(ntracrooh) = jt
       END IF
    END DO

    ntracrcoo2 = 0
    DO i = 1, ndefrcoo2
       jt = get_spc_ndx(names_rcoo2(i))
       IF (jt > 0) THEN
          ntracrcoo2 = ntracrcoo2 + 1
          ndxrcoo2(ntracrcoo2) = jt
       END IF
    END DO

    ! --- Set indices for reaction rate diagnostics
    ! Ox production
    ndxpo3_ho2(1) = get_rxt_ndx('NO_HO2')
    ndxpo3_ho2(2) = get_spc_ndx('NO')
    ndxpo3_ho2(3) = get_spc_ndx('HO2')

    ndxpo3_ch3o2(1) = get_rxt_ndx('CH3O2_NO')
    ndxpo3_ch3o2(2) = get_spc_ndx('CH3O2')
    ndxpo3_ch3o2(3) = get_spc_ndx('NO')

    ndxpo3_hno3(1) = get_rxt_ndx('HNO3_OH')
    ndxpo3_hno3(2) = get_spc_ndx('HNO3')
    ndxpo3_hno3(3) = get_spc_ndx('OH')

    ndxpo3_hono(1) = get_rxt_ndx('HONO_OH')
    ndxpo3_hono(2) = get_spc_ndx('HONO')
    ndxpo3_hono(3) = get_spc_ndx('OH')
    
    ndxpo3_pan(1) = get_rxt_ndx('PAN_OH')
    ndxpo3_pan(2) = get_spc_ndx('PAN')
    ndxpo3_pan(3) = get_spc_ndx('OH')

    ndxpo3_jpan(1) = get_rxt_ndx('jpan')
    ndxpo3_jpan(2) = get_spc_ndx('PAN')
    ndxpo3_jpan(3) = -1

    ndxpo3_jho2no2(1) = get_rxt_ndx('jho2no2_a')
    ndxpo3_jho2no2(2) = get_spc_ndx('HO2NO2')
    ndxpo3_jho2no2(3) = -1

    DO i = 1, ntracjanit
       namestemp=tolower(names_janit(i)) !HK #531 circumvent Cray bug
       ndxpo3_janit(i,1) = get_rxt_ndx('j'//TRIM(namestemp)) ! lowercase 
       ndxpo3_janit(i,2) = get_spc_ndx(names_janit(i))
       ndxpo3_janit(i,3) = -1
    END DO

    ndxlhhno3(1) = get_rxt_ndx('het_HNO3_tropo')
    ndxlhhno3(2) = get_spc_ndx('HNO3')
    ndxlhhno3(3) = -1
    
    DO i = 1, ntracro2
       ndxpo3_ro2(i,1) = get_rxt_ndx(TRIM(tracnam(ndxro2(i)))//'_NO')
       ndxpo3_ro2(i,2) = get_spc_ndx(tracnam(ndxro2(i)))
       ndxpo3_ro2(i,3) = get_spc_ndx('NO')
    END DO

    ! O3 loss
    ndxlo3_ho2(1) = get_rxt_ndx('HO2_O3')
    ndxlo3_ho2(2) = get_spc_ndx('O3')
    ndxlo3_ho2(3) = get_spc_ndx('HO2')

    ndxlo3_oh(1) = get_rxt_ndx('OH_O3')
    ndxlo3_oh(2) = get_spc_ndx('O3')
    ndxlo3_oh(3) = get_spc_ndx('OH')

    ndxlo3_o1d(1) = get_rxt_ndx('O1D_H2O')
    ndxlo3_o1d(2) = get_spc_ndx('O1D')
    ndxlo3_o1d(3) = get_spc_ndx('H2O')

    ndxlo3_no(1) = get_rxt_ndx('NO_O3')
    ndxlo3_no(2) = get_spc_ndx('O3')
    ndxlo3_no(3) = get_spc_ndx('NO')

    ndxlo3_no2(1) = get_rxt_ndx('NO2_O3')
    ndxlo3_no2(2) = get_spc_ndx('O3')
    ndxlo3_no2(3) = get_spc_ndx('NO2')

    DO i = 1, ntracene
       ndxlo3_ene(i,1) = get_rxt_ndx(TRIM(names_ene(i))//'_O3')
       ndxlo3_ene(i,2) = get_spc_ndx(names_ene(i))
       ndxlo3_ene(i,3) = get_spc_ndx('O3')
    END DO
 
    DO i = 1, ntracphenoxy
       ndxlo3_c6h5o(i,1) = get_rxt_ndx(TRIM(names_phenoxy(i))//'_O3')
       ndxlo3_c6h5o(i,2) = get_spc_ndx(names_phenoxy(i))
       ndxlo3_c6h5o(i,3) = get_spc_ndx('O3')
    END DO
    
    DO i = 1, ntracphenoxy
       ndxlo3_c6h5ono2(i,1) = get_rxt_ndx(TRIM(names_phenoxy(i))//'_NO2')
       ndxlo3_c6h5ono2(i,2) = get_spc_ndx(names_phenoxy(i))
       ndxlo3_c6h5ono2(i,3) = get_spc_ndx('NO2')
    END DO
    
    ! Heterogeneous Ox loss
    ndxlox_het(1,1) = get_rxt_ndx('het_N2O5_tropo')
    ndxlox_het(1,2) = get_spc_ndx('N2O5')
    ndxlox_het(1,3) = -1
    ndxlox_het(2,1) = get_rxt_ndx('het_O3_tropo') 
    ndxlox_het(2,2) = get_spc_ndx('O3')
    ndxlox_het(2,3) = -1
    ndxlox_het(3,1) = get_rxt_ndx('het_NO2_tropo') 
    ndxlox_het(3,2) = get_spc_ndx('NO2')
    ndxlox_het(3,3) = -1
    ndxlox_het(4,1) = get_rxt_ndx('het_HNO3_tropo') 
    ndxlox_het(4,2) = get_spc_ndx('HNO3')
    ndxlox_het(4,3) = -1
    ndxlox_het(5,1) = get_rxt_ndx('het_NO3_tropo') 
    ndxlox_het(5,2) = get_spc_ndx('NO3')
    ndxlox_het(5,3) = -1

    ! Loss due to NO3
    ndxlno3_jno3a(1) = get_rxt_ndx('jno3_a') ! NO3 + hv -> NO2 + O 
    ndxlno3_jno3a(2) = get_spc_ndx('NO3')
    ndxlno3_jno3a(3) = -1

    ndxlno3_jno3b(1) = get_rxt_ndx('jno3_b') ! NO3 + hv -> NO + O2
    ndxlno3_jno3b(2) = get_spc_ndx('NO3')
    ndxlno3_jno3b(3) = -1

    ndxlno3_ho2(1) = get_rxt_ndx('NO3_HO2') ! heterogeneous HNO3 production
    ndxlno3_ho2(2) = get_spc_ndx('NO3')
    ndxlno3_ho2(3) = get_spc_ndx('HO2')

    DO i = 1, ntracrcho
       ndxlno3_rcho(i,1) = get_rxt_ndx(TRIM(names_rcho(i))//'_NO3')
       ndxlno3_rcho(i,2) = get_spc_ndx(names_rcho(i))
       ndxlno3_rcho(i,3) = get_spc_ndx('NO3')
    END DO

    DO i = 1, ntracro2
       ndxlno3_ro2no3(i,1) = get_rxt_ndx(TRIM(tracnam(ndxro2(i)))//'_NO3')
       ndxlno3_ro2no3(i,2) = get_spc_ndx(tracnam(ndxro2(i)))
       ndxlno3_ro2no3(i,3) = get_spc_ndx('NO3')
    END DO

    DO i = 1, ntracrene
       ndxlno3_rene(i,1) = get_rxt_ndx(TRIM(names_rene(i))//'_NO3')
       ndxlno3_rene(i,2) = get_spc_ndx(names_rene(i))
       ndxlno3_rene(i,3) = get_spc_ndx('NO3')
    END DO

    ndxlo3_xo(1,1) = get_rxt_ndx('BRO_HO2')   ! *1
    ndxlo3_xo(1,2) = get_spc_ndx('BRO')
    ndxlo3_xo(1,3) = get_spc_ndx('HO2')
    ndxlo3_xo(2,1) = get_rxt_ndx('CLO_HO2')   ! *1
    ndxlo3_xo(2,2) = get_spc_ndx('CLO') 
    ndxlo3_xo(2,3) = get_spc_ndx('HO2') 
    ndxlo3_xo(3,1) = get_rxt_ndx('CLO_OH_a')  ! *1
    ndxlo3_xo(3,2) = get_spc_ndx('CLO') 
    ndxlo3_xo(3,3) = get_spc_ndx('OH')
    ndxlo3_xo(4,1) = get_rxt_ndx('CLO_OH_b')  ! *1
    ndxlo3_xo(4,2) = get_spc_ndx('CLO') 
    ndxlo3_xo(4,3) = get_spc_ndx('OH')
    ndxlo3_xo(5,1) = get_rxt_ndx('CLO_CH3O2') ! *1
    ndxlo3_xo(5,2) = get_spc_ndx('CLO') 
    ndxlo3_xo(5,3) = get_spc_ndx('CH3O2')
    ndxlo3_xo(6,1) = get_rxt_ndx('BRO_CLO_b') ! *2
    ndxlo3_xo(6,2) = get_spc_ndx('BRO') 
    ndxlo3_xo(6,3) = get_spc_ndx('CLO')
    ndxlo3_xo(7,1) = get_rxt_ndx('BRO_CLO_c') ! *2
    ndxlo3_xo(7,2) = get_spc_ndx('BRO') 
    ndxlo3_xo(7,3) = get_spc_ndx('CLO')
    ndxlo3_xo(8,1) = get_rxt_ndx('BRO_BRO')   ! *2
    ndxlo3_xo(8,2) = get_spc_ndx('BRO') 
    ndxlo3_xo(8,3) = get_spc_ndx('BRO')
    ndxlo3_xo(9,1) = get_rxt_ndx('CLO_CLO_a') ! *2
    ndxlo3_xo(9,2) = get_spc_ndx('CLO') 
    ndxlo3_xo(9,3) = get_spc_ndx('CLO')
    ndxlo3_xo(10,1) = get_rxt_ndx('CLO_CLO_b')! *2
    ndxlo3_xo(10,2) = get_spc_ndx('CLO') 
    ndxlo3_xo(10,3) = get_spc_ndx('CLO')
    nreacxo = 5 ! *1
    ! + 5 count double

 
    ! CO
    ndxco_oh(1,1) = get_rxt_ndx('CO_OH_a')
    ndxco_oh(1,2) = get_spc_ndx('CO')
    ndxco_oh(1,3) = get_spc_ndx('OH')
    ndxco_oh(2,1) = get_rxt_ndx('CO_OH_b')
    ndxco_oh(2,2) = get_spc_ndx('CO')
    ndxco_oh(2,3) = get_spc_ndx('OH')

    ! CH4
    ndxch4_oh(1) = get_rxt_ndx('CH4_OH')
    ndxch4_oh(2) = get_spc_ndx('CH4')
    ndxch4_oh(3) = get_spc_ndx('OH')

    ! CH4 loss ## careful about mechanism changes!
    ndxlch4(1,:) = ndxch4_oh(:)
    ndxlch4(2,1) = get_rxt_ndx('jch4_a')
    ndxlch4(2,2) = get_spc_ndx('CH4')
    ndxlch4(2,3) = -1
    ndxlch4(3,1) = get_rxt_ndx('jch4_b')
    ndxlch4(3,2) = get_spc_ndx('CH4')
    ndxlch4(3,3) = -1
    ndxlch4(4,1) = get_rxt_ndx('O1D_CH4_a')
    ndxlch4(4,2) = get_spc_ndx('O1D')
    ndxlch4(4,3) = get_spc_ndx('CH4')
    ndxlch4(5,1) = get_rxt_ndx('O1D_CH4_b')
    ndxlch4(5,2) = get_spc_ndx('O1D')
    ndxlch4(5,3) = get_spc_ndx('CH4')
    ndxlch4(6,1) = get_rxt_ndx('O1D_CH4_c')
    ndxlch4(6,2) = get_spc_ndx('O1D')
    ndxlch4(6,3) = get_spc_ndx('CH4')
    ndxlch4(7,1) = get_rxt_ndx('CL_CH4')
    ndxlch4(7,2) = get_spc_ndx('CL')
    ndxlch4(7,3) = get_spc_ndx('CH4')
    ndxlch4(8,1) = get_rxt_ndx('F_CH4')
    ndxlch4(8,2) = get_spc_ndx('F')
    ndxlch4(8,3) = get_spc_ndx('CH4')
    nreacch4 = 8

    ! PAN
    ndxlpan(1,1) = get_rxt_ndx('jpan')
    ndxlpan(1,2) = get_spc_ndx('PAN')
    ndxlpan(1,3) = -1
    ndxlpan(2,1) = get_rxt_ndx('PAN')
    ndxlpan(2,2) = get_spc_ndx('PAN')
    ndxlpan(2,3) = -1
    ndxlpan(3,1) = get_rxt_ndx('PAN_OH')
    ndxlpan(3,2) = get_spc_ndx('PAN')
    ndxlpan(3,3) = get_spc_ndx('OH')
    nreacpan = 3
    ! s.stadtler Production PAN
    ndxppan(1) = get_rxt_ndx('CH3CO3_NO2')
    ndxppan(2) = get_spc_ndx('CH3CO3')
    ndxppan(3) = get_spc_ndx('NO2')
    ! s.stadtler
 
    ! NO2+OH
    ndxno2_oh(1) = get_rxt_ndx('NO2_OH') !HNO3 prodcution (for heterogeneous production see n2o5 and no3)
    ndxno2_oh(2) = get_spc_ndx('NO2')
    ndxno2_oh(3) = get_spc_ndx('OH')

    ndxno2_ohreact(1) = get_rxt_ndx('NO2_OH') 
    ndxno2_ohreact(2) = get_spc_ndx('NO2')
    ndxno2_ohreact(3) = -1

    ! NO2 Photolysis
    ndxno2_hv(1) = get_rxt_ndx('jno2') 
    ndxno2_hv(2) = get_spc_ndx('NO2')
    ndxno2_hv(3) = -1

    !HNO3 loss via Photolysis and OH s.stadtler
    ndxlghno3(1,1) = get_rxt_ndx('jhno3')
    ndxlghno3(1,2) = get_spc_ndx('HNO3')
    ndxlghno3(1,3) = -1
    ndxlghno3(2,1) = get_rxt_ndx('HNO3_OH')
    ndxlghno3(2,2) = get_spc_ndx('HNO3')
    ndxlghno3(2,3) = get_spc_ndx('OH')
    nreachno3 = 2

    ! N2O5 s.stadtler
    ! N2O5 Photolysis Muss noch verändert werden!
    ndxn2o5(1,1) = get_rxt_ndx('jn2o5_a') ! Das haette ich doch in eine Schleife schreiben können?
    ndxn2o5(1,2) = get_spc_ndx('N2O5')    ! Ist nur eh bald anders...
    ndxn2o5(1,3) = -1
    ndxn2o5(2,1) = get_rxt_ndx('jn2o5_b')
    ndxn2o5(2,2) = get_spc_ndx('N2O5')
    ndxn2o5(2,3) = -1
    nreacn2o5 = 2
    ! Thermal Loss
    ndxlgn2o5(1) = get_rxt_ndx('N2O5')
    ndxlgn2o5(2) = get_spc_ndx('N2O5')
    ndxlgn2o5(3) = -1
    ! NO2+NO3
    ndxno2_no3(1) = get_rxt_ndx('NO2_NO3')
    ndxno2_no3(2) = get_spc_ndx('NO2')
    ndxno2_no3(3) = get_spc_ndx('NO3')

    ! NO3
    ! Production Schleife nicht moeglich 
    ndxpno3(1,1)= get_rxt_ndx('jho2no2_a')
    ndxpno3(1,2)= get_spc_ndx('HO2NO2')
    ndxpno3(1,3)= -1
    ndxpno3(2,1)= get_rxt_ndx('jn2o5_a')
    ndxpno3(2,2)= get_spc_ndx('N2O5')
    ndxpno3(2,3)= -1
    ndxpno3(3,1)= get_rxt_ndx('jn2o5_b')
    ndxpno3(3,2)= get_spc_ndx('N2O5')
    ndxpno3(3,3)= -1
    ndxpno3(4,1)= get_rxt_ndx('jpan') !Nur 0.4 NO3 werden produziert
    ndxpno3(4,2)= get_spc_ndx('PAN')
    ndxpno3(4,3)= -1
    ndxpno3(5,1)= get_rxt_ndx('HNO3_OH')
    ndxpno3(5,2)= get_spc_ndx('HNO3')
    ndxpno3(5,3)= get_spc_ndx('OH')
    ndxpno3(6,1)= get_rxt_ndx('PAN_OH')
    ndxpno3(6,2)= get_spc_ndx('PAN')
    ndxpno3(6,3)= get_spc_ndx('OH')
    ndxpno3(7,1)= get_rxt_ndx('N2O5') 
    ndxpno3(7,2)= get_spc_ndx('N2O5')
    ndxpno3(7,3)= -1
    ndxpno3(8,1)= get_rxt_ndx('NO2_O3') 
    ndxpno3(8,2)= get_spc_ndx('NO2')
    ndxpno3(8,3)= get_spc_ndx('O3')
    nreacpno3 = 8
    ! Gasphase Loss
    nreacno3 = 3
    DO i = 1, nreacno3
       ndxlgno3(i,1) = get_rxt_ndx('NO3_'//TRIM(names_reac_no3(i)))
       ndxlgno3(i,2) = get_spc_ndx(names_reac_no3(i))
       ndxlgno3(i,3) = get_spc_ndx('NO3')
    END DO
    nreacno3 = 13 
    DO i = 1, nreacno3
       ndxlgno3(i+3,1) = get_rxt_ndx(TRIM(names_ene_no3(i))//'_NO3')
       ndxlgno3(i+3,2) = get_spc_ndx(names_ene_no3(i))
       ndxlgno3(i+3,3) = get_spc_ndx('NO3')
    END DO
    ! Heterogeneous NOx loss via HNO3
    ! Actually heterogeneous HNO3 production, which is lost instantaneously (see jam003a)
    ndxlnox_het(1,1) = get_rxt_ndx('het_NO3_tropo') 
    ndxlnox_het(1,2) = get_spc_ndx('NO3')
    ndxlnox_het(1,3) = -1
    ndxlnox_het(2,1) = get_rxt_ndx('het_N2O5_tropo')
    ndxlnox_het(2,2) = get_spc_ndx('N2O5')
    ndxlnox_het(2,3) = -1
    ndxlnox_het(3,1) = get_rxt_ndx('het_HNO3_tropo') 
    ndxlnox_het(3,2) = get_spc_ndx('HNO3')
    ndxlnox_het(3,3) = -1
    ndxlnox_het(4,1) = get_rxt_ndx('het_NO2_tropo') ! scale by 0.5 
    ndxlnox_het(4,2) = get_spc_ndx('NO2')
    ndxlnox_het(4,3) = -1

    ! NOx
    ! Production
    ndxpnox(1,1)= get_rxt_ndx('jno3_a')
    ndxpnox(1,2)= get_spc_ndx('NO3')
    ndxpnox(1,3)= -1
    ndxpnox(2,1)= get_rxt_ndx('jno3_b')
    ndxpnox(2,2)= get_spc_ndx('NO3')
    ndxpnox(2,3)= -1
    ndxpnox(3,1)= get_rxt_ndx('NO3_OH')
    ndxpnox(3,2)= get_spc_ndx('NO3')
    ndxpnox(3,3)= get_spc_ndx('OH')
    ndxpnox(4,1)= get_rxt_ndx('NO3_HO2') 
    ndxpnox(4,2)= get_spc_ndx('NO3')
    ndxpnox(4,3)= get_spc_ndx('HO2')
    ndxpnox(5,1)= get_rxt_ndx('NO3_NO')
    ndxpnox(5,2)= get_spc_ndx('NO3')
    ndxpnox(5,3)= get_spc_ndx('NO')
    ndxpnox(6,1)= get_rxt_ndx('jhno3')
    ndxpnox(6,2)= get_spc_ndx('HNO3')
    ndxpnox(6,3)= -1
    ndxpnox(7,1)= get_rxt_ndx('jn2o5_a')! Weil das zwei mal vorkommt
    ndxpnox(7,2)= get_spc_ndx('N2O5')
    ndxpnox(7,3)= -1
    ndxpnox(8,1)= get_rxt_ndx('jn2o5_a')
    ndxpnox(8,2)= get_spc_ndx('N2O5')
    ndxpnox(8,3)= -1
    ndxpnox(9,1)= get_rxt_ndx('jn2o5_b')
    ndxpnox(9,2)= get_spc_ndx('N2O5')
    ndxpnox(9,3)= -1
    ndxpnox(10,1)= get_rxt_ndx('jn2o5_b')
    ndxpnox(10,2)= get_spc_ndx('N2O5')
    ndxpnox(10,3)= -1
    ndxpnox(11,1)= get_rxt_ndx('N2O5')
    ndxpnox(11,2)= get_spc_ndx('N2O5')
    ndxpnox(11,3)= -1
    ndxpnox(12,1)= get_rxt_ndx('jhono')
    ndxpnox(12,2)= get_spc_ndx('HONO')
    ndxpnox(12,3)= -1
    ndxpnox(13,1)= get_rxt_ndx('HONO_OH')
    ndxpnox(13,2)= get_spc_ndx('HONO')
    ndxpnox(13,3)= get_spc_ndx('OH')
    ndxpnox(14,1)= get_rxt_ndx('jpan')! nur 0.6 NOx werden produziert
    ndxpnox(14,2)= get_spc_ndx('PAN')
    ndxpnox(14,3)= -1
    ndxpnox(15,1)= get_rxt_ndx('PAN')
    ndxpnox(15,2)= get_spc_ndx('PAN')
    ndxpnox(15,3)= -1
    ndxpnox(16,1)= get_rxt_ndx('jho2no2_b')
    ndxpnox(16,2)= get_spc_ndx('HO2NO2')
    ndxpnox(16,3)= -1
    ndxpnox(17,1)= get_rxt_ndx('HO2NO2')
    ndxpnox(17,2)= get_spc_ndx('HO2NO2')
    ndxpnox(17,3)= -1
    ndxpnox(18,1)= get_rxt_ndx('HO2NO2_OH')
    ndxpnox(18,2)= get_spc_ndx('HO2NO2')
    ndxpnox(18,3)= get_spc_ndx('OH')
    nreacpnox = 18
    ! NOx loss
    ndxlnox(1,1)= get_rxt_ndx('NO2_OH')
    ndxlnox(1,2)= get_spc_ndx('NO2')
    ndxlnox(1,3)= get_spc_ndx('OH')
    ndxlnox(2,1)= get_rxt_ndx('NO2_HO2')
    ndxlnox(2,2)= get_spc_ndx('NO2')
    ndxlnox(2,3)= get_spc_ndx('HO2')
    ndxlnox(3,1)= get_rxt_ndx('NO2_O3')
    ndxlnox(3,2)= get_spc_ndx('NO2')
    ndxlnox(3,3)= get_spc_ndx('O3')
    ndxlnox(4,1)= get_rxt_ndx('NO2_NO3')
    ndxlnox(4,2)= get_spc_ndx('NO2')
    ndxlnox(4,3)= get_spc_ndx('NO3')
    ndxlnox(5,1)= get_rxt_ndx('CH3CO3_NO2')
    ndxlnox(5,2)= get_spc_ndx('CH3CO3')
    ndxlnox(5,3)= get_spc_ndx('NO2')
    ndxlnox(6,1)= get_rxt_ndx('NO_OH')
    ndxlnox(6,2)= get_spc_ndx('NO')
    ndxlnox(6,3)= get_spc_ndx('OH')
    nreaclnox = 6

    ! HO2+HO2
    ndxho2_ho2(1) = get_rxt_ndx('HO2_HO2')
    ndxho2_ho2(2) = get_spc_ndx('HO2')
    ndxho2_ho2(3) = get_spc_ndx('HO2')

    ! OH reactivity
    ntracohreact = 0
    DO jt=1,pcnstm1
      IF ((get_rxt_ndx(TRIM(tracnam(jt))//'_OH') > 0) .OR. &
          (get_rxt_ndx('OH_'//TRIM(tracnam(jt))) > 0)) THEN
        ntracohreact = ntracohreact + 1
        i = get_rxt_ndx(TRIM(tracnam(jt))//'_OH')
        IF (i > 0) THEN
          ndxohreact(ntracohreact,1) = i
        ELSE
          ndxohreact(ntracohreact,1) = get_rxt_ndx('OH_'//TRIM(tracnam(jt)))
        END IF
        ndxohreact(ntracohreact,2)   = get_spc_ndx(tracnam(jt))
        ndxohreact(ntracohreact,3)   = -1
      END IF
    END DO
    ! with the above loop reactions labeled _a/_b/_c/_d (OH_OH, CO_OH, CLO_OH, DMS_OH)
    ! are not found ==> do this explicitly
    ntracohreact = ntracohreact + 1
    ndxohreact(ntracohreact,1) = get_rxt_ndx('OH_OH_a')
    ndxohreact(ntracohreact,2) = get_spc_ndx('OH')
    ndxohreact(ntracohreact,3) = -1
    ntracohreact = ntracohreact + 1
    ndxohreact(ntracohreact,1) = get_rxt_ndx('OH_OH_b')
    ndxohreact(ntracohreact,2) = get_spc_ndx('OH')
    ndxohreact(ntracohreact,3) = -1
    ntracohreact = ntracohreact + 1
    ndxohreact(ntracohreact,1) = get_rxt_ndx('CO_OH_a')
    ndxohreact(ntracohreact,2) = get_spc_ndx('CO')
    ndxohreact(ntracohreact,3) = -1
    ntracohreact = ntracohreact + 1
    ndxohreact(ntracohreact,1) = get_rxt_ndx('CO_OH_b')
    ndxohreact(ntracohreact,2) = get_spc_ndx('CO')
    ndxohreact(ntracohreact,3) = -1
    ntracohreact = ntracohreact + 1
    ndxohreact(ntracohreact,1) = get_rxt_ndx('CLO_OH_a')
    ndxohreact(ntracohreact,2) = get_spc_ndx('CLO')
    ndxohreact(ntracohreact,3) = -1
    ntracohreact = ntracohreact + 1
    ndxohreact(ntracohreact,1) = get_rxt_ndx('CLO_OH_b')
    ndxohreact(ntracohreact,2) = get_spc_ndx('CLO')
    ndxohreact(ntracohreact,3) = -1 
    ntracohreact = ntracohreact + 1
    ndxohreact(ntracohreact,1) = get_rxt_ndx('DMS_OH_a')
    ndxohreact(ntracohreact,2) = get_spc_ndx('DMS')
    ndxohreact(ntracohreact,3) = -1
    ntracohreact = ntracohreact + 1
    ndxohreact(ntracohreact,1) = get_rxt_ndx('DMS_OH_b')
    ndxohreact(ntracohreact,2) = get_spc_ndx('DMS')
    ndxohreact(ntracohreact,3) = -1
    ! Total OH loss
    ! count twice for OH+OH
    ! not right when ndxohreact gets used !!!!
    ntracoh = ntracohreact
    DO i = 1, ntracoh
      DO jt = 1, 2
        ndxloh(i,jt) = ndxohreact(i,jt)
      END DO
        ndxloh(i,3)  = get_spc_ndx('OH')
    END DO
    ntracoh = ntracoh + 1
    ndxloh(ntracoh,1) = get_rxt_ndx('OH_OH_a')
    ndxloh(ntracoh,2) = get_spc_ndx('OH')
    ndxloh(ntracoh,3) = get_spc_ndx('OH')
    ntracoh = ntracoh + 1
    ndxloh(ntracoh,1) = get_rxt_ndx('OH_OH_b')
    ndxloh(ntracoh,2) = get_spc_ndx('OH')
    ndxloh(ntracoh,3) = get_spc_ndx('OH')

    DO i = 1, ntracrcoo2
       ndxlrcoo2_no2(i,1) = get_rxt_ndx(TRIM(tracnam(ndxrcoo2(i)))//'_NO2')
       ndxlrcoo2_no2(i,2) = get_spc_ndx(tracnam(ndxrcoo2(i)))
       ndxlrcoo2_no2(i,3) = get_spc_ndx('NO2')
    END DO

!!++mgs20140320: temporary addition - possibly irrelevant later
    ndxn2o_j(1) = get_rxt_ndx('jn2o')
    ndxn2o_j(2) = get_spc_ndx('N2O')
    ndxn2o_j(3) = -1

    ndxn2o_o1d(1,1) = get_rxt_ndx('O1D_N2O_a')
    ndxn2o_o1d(1,2) = get_spc_ndx('N2O')
    ndxn2o_o1d(1,3) = get_spc_ndx('O1D')
    ndxn2o_o1d(2,1) = get_rxt_ndx('O1D_N2O_b')
    ndxn2o_o1d(2,2) = get_spc_ndx('N2O')
    ndxn2o_o1d(2,3) = get_spc_ndx('O1D')
!!--mgs


  END SUBROUTINE init_moz_streams


!!mgs!!  ! stream _tmdiag:
!!mgs!!  ! tracers for which the nbudget flag is set
!!mgs!!  !--------------------------------------------------------------------------
!!mgs!!  USE mo_tracdef,        ONLY: trlist
!!mgs!!  USE mo_tracer,         ONLY: ntrac
!!mgs!!  USE mo_decomposition,  ONLY: lc => local_decomposition
!!mgs!!  !create tmdiag stream and add variables
!!mgs!!  IF (ANY(trlist%ti(1:ntrac)%nbudget == 1)) THEN
!!mgs!!     CALL new_stream (tmdiag,'tmdiag')
!!mgs!!     CALL default_stream_setting (tmdiag,lrerun = .FALSE., &
!!mgs!!                                  contnorest = .TRUE., &
!!mgs!!                                  laccu = .FALSE., &
!!mgs!!                                  lpost = .TRUE., &
!!mgs!!                                  table = 199, &
!!mgs!!                                  code = AUTO )
!!mgs!!     DO itrac = 1,ntrac
!!mgs!!        IF (trlist%ti(itrac)%nbudget == 1) THEN
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_emission', &
!!mgs!!                tmdiag_p, &
!!mgs!!                longname='mass budget of '//TRIM(trlist%ti(itrac)%basename)// &
!!mgs!!                ' emissions', units='1/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_drydep', &
!!mgs!!                tmdiag_sur, &
!!mgs!!                longname='mass budget of '//TRIM(trlist%ti(itrac)%basename)// &
!!mgs!!                ' dry deposition', units='1/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_ddepvel', &
!!mgs!!                tmdiag_sur, &
!!mgs!!                longname='dry deposition velocity of ' &
!!mgs!!                //TRIM(trlist%ti(itrac)%basename), &
!!mgs!!                units='m/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_chem', &
!!mgs!!                tmdiag_p, &
!!mgs!!                longname='mass budget of '//TRIM(trlist%ti(itrac)%basename)// &
!!mgs!!                ' chemical reaction', units='1/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_diff', &
!!mgs!!                tmdiag_p, &
!!mgs!!                longname='mass budget of '//TRIM(trlist%ti(itrac)%basename)// &
!!mgs!!                ' diffusion', units='1/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_conv', &
!!mgs!!                tmdiag_p, &
!!mgs!!                longname='mass budget of '//TRIM(trlist%ti(itrac)%basename)// &
!!mgs!!                ' convection', units='1/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_adv', &
!!mgs!!                tmdiag_p, &
!!mgs!!                longname='mass budget of '//TRIM(trlist%ti(itrac)%basename)// &
!!mgs!!                ' advection', units='1/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_wetdepc', &
!!mgs!!                tmdiag_p, longname='mass budget of ' &
!!mgs!!                //TRIM(trlist%ti(itrac)%basename)//' conv wet deposition', &
!!mgs!!                units='1/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_wetdeps', &
!!mgs!!                tmdiag_p, longname='mass budget of ' &
!!mgs!!                //TRIM(trlist%ti(itrac)%basename)//' strat wet deposition', &
!!mgs!!                units='1/s')
!!mgs!!           CALL add_stream_element (tmdiag, &
!!mgs!!                TRIM(trlist%ti(itrac)%basename)//'_total', &
!!mgs!!                tmdiag_p, longname='total mass budget of ' &
!!mgs!!                //TRIM(trlist%ti(itrac)%basename), units='1/s')
!!mgs!!        END IF
!!mgs!!     END DO
!!mgs!!  END IF


!! extra drydep diagnostics (integrated mass flux ...)
!! make use of ddep diag_list in mo_hammoz_drydep !!!

!!baustelle!!   moz diagnostics copied from now obsolete moz_drydep routine - should be re-activated
!!   !++jsr diagnostics
!! !!baustelle!!  Needs to be adapted for MOZ3 chemistry without OX definition (plain O3)
!!   ! the OX group has been renamed in order to have the dry deposition inherent
!!   ! names. The tmdiag stream still contains the name OX. So, connect OX with
!!   ! O3 here.
!!       IF (do_ox_pl) THEN
!!          CALL get_tracer ('O3', idx=idx_ox)
!!          IF (trlist%ti(idx_ox)%ndrydep == 2) THEN
!!             CALL set_spec_ddep_flux ( zdrydepflux(1:kproma,idx_ox), kproma, krow, 'OX')
!!             CALL set_spec_ddep_flux ( pxtm1(1:kproma,klev,idx_ox)*pdensair(1:kproma)* &
!!                  ox_veg_vddep(1:kproma,krow), kproma, krow, 'OX_VEG')
!!          END IF
!!       END IF
!!       IF (do_co_pl) THEN
!!          CALL get_tracer ('CO', idx=idx_co)
!!          IF (trlist%ti(idx_co)%ndrydep == 2) THEN
!!             CALL set_spec_ddep_flux ( zdrydepflux(1:kproma,idx_co), kproma, krow, 'CO')
!!          END IF
!!       END IF
!!       IF (do_diag_dep) THEN
!!          CALL set_ddepstream (zdrydepflux, vdrydep(:,:,krow), kbdim, kproma, krow, ntrac)
!!       END IF
!!       IF (ANY(trlist%ti(1:ntrac)%nbudget == 1)) THEN
!!        CALL get_stream(tmdiag,'tmdiag')
!!        DO jt=1,ntrac
!!           IF ( trlist%ti(jt)%nbudget == 1 ) THEN
!!              specname = trlist%ti(jt)%basename
!!              IF (specname == 'O3') THEN
!!                 specname = 'OX'
!!              END IF
!!              CALL get_stream_element(tmdiag,TRIM(specname)//'_drydep',tmdiag_sur)
!!              IF ( trlist%ti(jt)%ndrydep == 2 ) THEN
!! !++mgs: replaced call to mflux_to_mmrt by direct computation
!! !!!             CALL mflux_to_mmrt (zdrydepflux(1:kproma,jt), zbuf(1:kproma), klev)
!!                 ! Compute air mass [kg]: dp/g * surface_area
!!                 zma(1:kproma) = ( paphp1(1:kproma,klev+1) - paphp1(1:kproma,klev) ) / g
!!                 ! Get the result
!!                 zbuf(1:kproma) = (zdrydepflux(1:kproma,jt)/zma(1:kproma))
!! !--mgs
!!                 tmdiag_sur(1:kproma,krow)=zbuf(1:kproma)
!!              ELSE
!!                 tmdiag_sur(1:kproma,krow)=0._dp
!!              END IF
!!           END IF
!!        ENDDO
!!     END IF
  ! s.stadtler
  SUBROUTINE mozh_diag(kproma, plonl, plev, lat, pcnstm1, vmr)

    INTEGER, INTENT(in)     :: kproma, plonl, plev, lat, pcnstm1
    REAL(dp), INTENT(in)    :: vmr(plonl, plev, pcnstm1)      ! species mixing ratios (from chemdr)

    INTEGER                 :: jt

    !
    dpvmro3(1:kproma, lat)  = 0._dp
    dpvmrno(1:kproma, lat)  = 0._dp
    dpvmrno2(1:kproma, lat) = 0._dp
    dpvmrco(1:kproma, lat)  = 0._dp
    DO jt = 1, ntracnoy
        IF (jt == 1) dpvmro3(1:kproma, lat)  = vmr(1:kproma, plev, ndxmozh(jt))
        IF (jt == 2) dpvmrno(1:kproma, lat)  = vmr(1:kproma, plev, ndxmozh(jt))
        IF (jt == 3) dpvmrno2(1:kproma, lat) = vmr(1:kproma, plev, ndxmozh(jt))
        IF (jt == 4) dpvmrco(1:kproma, lat)  = vmr(1:kproma, plev, ndxmozh(jt))
    END DO

  END SUBROUTINE mozh_diag
  ! s.stadtler


  ! ---- calculate tracer family mixing ratios ----

  SUBROUTINE moz_family_diag(kproma, plonl, plev, lat, pcnstm1, vmr)

    INTEGER, INTENT(in)     :: kproma, plonl, plev, lat, pcnstm1
    REAL(dp), INTENT(in)    :: vmr(plonl, plev, pcnstm1)      ! species mixing ratios (from chemdr)

    INTEGER                 :: jt

    ! NOy
    dpnoy(1:kproma, 1:plev, lat) = 0._dp
    DO jt = 1, ntracnoy
        dpnoy(1:kproma, 1:plev, lat) =  dpnoy(1:kproma, 1:plev, lat) + vmr(1:kproma, 1:plev, ndxnoy(jt))
    END DO
    ! Cly
    dpcly(1:kproma, 1:plev, lat) = 0._dp
    DO jt = 1, ntraccly
        dpcly(1:kproma, 1:plev, lat) =  dpcly(1:kproma, 1:plev, lat) + vmr(1:kproma, 1:plev, ndxcly(jt))
    END DO
    ! Bry
    dpbry(1:kproma, 1:plev, lat) = 0._dp
    DO jt = 1, ntracbry
        dpbry(1:kproma, 1:plev, lat) =  dpbry(1:kproma, 1:plev, lat) + vmr(1:kproma, 1:plev, ndxbry(jt))
    END DO
    ! RO2
    dpro2(1:kproma, 1:plev, lat) = 0._dp
    DO jt = 1, ntracro2
        dpro2(1:kproma, 1:plev, lat) =  dpro2(1:kproma, 1:plev, lat) + vmr(1:kproma, 1:plev, ndxro2(jt))
    END DO
    ! ROOH
    dprooh(1:kproma, 1:plev, lat) = 0._dp
    DO jt = 1, ntracrooh
        dprooh(1:kproma, 1:plev, lat) =  dprooh(1:kproma, 1:plev, lat) + vmr(1:kproma, 1:plev, ndxrooh(jt))
    END DO
    ! RCOO2
    dprcoo2(1:kproma, 1:plev, lat) = 0._dp
    DO jt = 1, ntracrcoo2
        dprcoo2(1:kproma, 1:plev, lat) =  dprcoo2(1:kproma, 1:plev, lat) + vmr(1:kproma, 1:plev, ndxrcoo2(jt))
    END DO


  END SUBROUTINE moz_family_diag

  ! ---- calculate reaction rate (turnover) of one specific reaction
  ! Note: this is the working routine; the interface to chemdr is moz_rates_diag below
  SUBROUTINE rate_diag(ndxtuple, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, turnover)

    USE mo_exception,       ONLY: finish, message_text

    INTEGER, INTENT(in)     :: ndxtuple(3), plonl, plev, pcnstm1, rxntot
    REAL(dp), INTENT(in)    :: vmr(plonl, plev, pcnstm1)               ! species mixing ratios (from chemdr)
    REAL(dp), INTENT(in)    :: reaction_rates(plonl, plev, rxntot)     ! reaction rates (from chemdr)
    REAL(dp), INTENT(out)   :: turnover(plonl, plev)                   ! k*A*B
    
    !##TODO: unit conversion (?)
    turnover(:,:) = 0._dp
    IF (ndxtuple(1) > 0) THEN
       turnover(:,:) = reaction_rates(:, :, ndxtuple(1))
       IF (ndxtuple(2) > 0) THEN
          turnover(:,:) = turnover(:,:) * vmr(:, :, ndxtuple(2))
       ELSE
          WRITE(message_text, '(a,3i4)') 'Invalid indices for reaction rate diagnostics: ', ndxtuple(:)
          CALL finish('moz_rates_diag', message_text)
       END IF
       IF (ndxtuple(3) > 0) THEN
          turnover(:,:) = turnover(:,:) * vmr(:, :, ndxtuple(3))
       END IF
    END IF

  END SUBROUTINE rate_diag
 
  SUBROUTINE moz_rates_diag(kproma, plonl, plev, lat, pcnstm1, rxntot, vmr, reaction_rates)

    INTEGER, INTENT(in)     :: kproma, plonl, plev, lat, pcnstm1, rxntot
    REAL(dp), INTENT(in)    :: vmr(plonl, plev, pcnstm1)               ! species mixing ratios (from chemdr)
    REAL(dp), INTENT(in)    :: reaction_rates(plonl, plev, rxntot)     ! reaction rates (from chemdr)

    INTEGER                 :: i
    REAL(dp)                :: zturnover(plonl, plev)

    ! s.stadtler Ox production
    IF (associated(dppo3_ho2)) THEN
       CALL rate_diag(ndxpo3_ho2, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dppo3_ho2(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dppo3_ch3o2)) THEN
       CALL rate_diag(ndxpo3_ch3o2, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dppo3_ch3o2(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dppo3_ro2)) THEN
       dppo3_ro2(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracro2
          CALL rate_diag(ndxpo3_ro2(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dppo3_ro2(1:kproma, 1:plev, lat) = dppo3_ro2(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dppo3_hno3)) THEN
       CALL rate_diag(ndxpo3_hno3, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dppo3_hno3(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dppo3_hono)) THEN
       CALL rate_diag(ndxpo3_hono, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dppo3_hono(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dppo3_pan)) THEN
       CALL rate_diag(ndxpo3_pan, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dppo3_pan(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dppo3_jpan)) THEN
       CALL rate_diag(ndxpo3_jpan, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dppo3_jpan(1:kproma, 1:plev, lat) = 0.4_dp*zturnover(1:kproma, 1:plev) ! s.stadtler I am not sure if this works!!
    END IF

    IF (associated(dppo3_jho2no2)) THEN
       CALL rate_diag(ndxpo3_jho2no2, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dppo3_jho2no2(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev) 
    END IF   
 
    IF (associated(dppo3_janit)) THEN
       dppo3_janit(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracjanit
          CALL rate_diag(ndxpo3_janit(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dppo3_janit(1:kproma, 1:plev, lat) = dppo3_janit(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO  
    END IF 
 
    ! s.stadtler Ox loss
    IF (associated(dplo3_ho2)) THEN
       CALL rate_diag(ndxlo3_ho2, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplo3_ho2(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dplo3_oh)) THEN
       CALL rate_diag(ndxlo3_oh, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplo3_oh(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dplo3_o1d)) THEN
       CALL rate_diag(ndxlo3_o1d, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplo3_o1d(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dplo3_ene)) THEN
       dplo3_ene(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracene
          CALL rate_diag(ndxlo3_ene(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplo3_ene(1:kproma, 1:plev, lat) = dplo3_ene(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO  
    END IF 

    IF (associated(dplo3_c6h5o)) THEN
       dplo3_c6h5o(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracphenoxy
         CALL rate_diag(ndxlo3_c6h5o(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dplo3_c6h5o(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dplo3_c6h5ono2)) THEN
       dplo3_c6h5ono2(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracphenoxy
         CALL rate_diag(ndxlo3_c6h5ono2(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dplo3_c6h5ono2(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dplox_het)) THEN
       dplox_het(1:kproma, 1:plev, lat) = 0._dp
       CALL rate_diag(ndxlox_het(1,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover) ! N2O5
       dplox_het(1:kproma, 1:plev, lat) = dplox_het(1:kproma, 1:plev, lat) + 3.0_dp*zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxlox_het(2,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover) !O3 
       dplox_het(1:kproma, 1:plev, lat) = dplox_het(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxlox_het(3,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover) ! NO2
       dplox_het(1:kproma, 1:plev, lat) = dplox_het(1:kproma, 1:plev, lat) + 0.5_dp*zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxlox_het(4,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover) ! HNO3
       dplox_het(1:kproma, 1:plev, lat) = dplox_het(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxlox_het(5,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover) ! NO3
       dplox_het(1:kproma, 1:plev, lat) = dplox_het(1:kproma, 1:plev, lat) + 2.0_dp*zturnover(1:kproma, 1:plev)
    END IF 

    IF (associated(dplno3_jno3a)) THEN
       CALL rate_diag(ndxlno3_jno3a, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplno3_jno3a(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev) 
    END IF

    IF (associated(dplno3_jno3b)) THEN
       CALL rate_diag(ndxlno3_jno3b, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplno3_jno3b(1:kproma, 1:plev, lat) = 2.0_dp*zturnover(1:kproma, 1:plev) ! s.stadtler I am not sure if this works!!
    END IF

    IF (associated(dplno3_ho2)) THEN
       CALL rate_diag(ndxlno3_ho2, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplno3_ho2(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev) 
    END IF

    IF (associated(dplno3_rcho)) THEN
       dplno3_rcho(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracrcho
          CALL rate_diag(ndxlno3_rcho(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplno3_rcho(1:kproma, 1:plev, lat) = dplno3_rcho(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO  
    END IF 

    IF (associated(dplno3_ro2no3)) THEN
       dplno3_ro2no3(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracro2
          CALL rate_diag(ndxlno3_ro2no3(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplno3_ro2no3(1:kproma, 1:plev, lat) = dplno3_ro2no3(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF 

    IF (associated(dplno3_rene)) THEN
       dplno3_rene(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracrene
          CALL rate_diag(ndxlno3_rene(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplno3_rene(1:kproma, 1:plev, lat) = dplno3_rene(1:kproma, 1:plev, lat) + 2.0_dp*zturnover(1:kproma, 1:plev)!Do not know if it is correct
       END DO  
    END IF 

    IF (associated(dplo3_xo)) THEN
       dplo3_xo(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreacxo ! *1 counting reactions
          CALL rate_diag(ndxlo3_xo(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplo3_xo(1:kproma, 1:plev, lat) = dplo3_xo(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO  
       DO i = nreacxo, nreacxo+5 ! *2 counting reactions
          CALL rate_diag(ndxlo3_xo(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplo3_xo(1:kproma, 1:plev, lat) = dplo3_xo(1:kproma, 1:plev, lat) + 2.0_dp*zturnover(1:kproma, 1:plev)
       END DO  
    END IF

    IF (associated(dplo3_no)) THEN
       CALL rate_diag(ndxlo3_no, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplo3_no(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dplo3_no2)) THEN
       CALL rate_diag(ndxlo3_no2, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplo3_no2(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF
 
    !! --- s.stadtler
    IF (associated(dprco_oh)) THEN
       CALL rate_diag(ndxco_oh(1,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dprco_oh(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxco_oh(2,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dprco_oh(1:kproma, 1:plev, lat) = dprco_oh(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dprch4_oh)) THEN
       CALL rate_diag(ndxch4_oh, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dprch4_oh(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dplch4)) THEN
       dplch4(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreacch4
          CALL rate_diag(ndxlch4(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplch4(1:kproma, 1:plev, lat) = dplch4(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dprno2_hv)) THEN
       CALL rate_diag(ndxno2_hv, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dprno2_hv(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dprho2_ho2)) THEN
       CALL rate_diag(ndxho2_ho2, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dprho2_ho2(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dplpan)) THEN
       dplpan(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreacpan
          CALL rate_diag(ndxlpan(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplpan(1:kproma, 1:plev, lat) = dplpan(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dplrcoo2_no2)) THEN
       dplrcoo2_no2(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracrcoo2
          CALL rate_diag(ndxlrcoo2_no2(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
          dplrcoo2_no2(1:kproma, 1:plev, lat) = dplrcoo2_no2(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO  
    END IF

    IF (associated(dploh)) THEN
       dploh(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracoh
         CALL rate_diag(ndxloh(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dploh(1:kproma, 1:plev, lat) = dploh(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dpohreact)) THEN
       dpohreact(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, ntracohreact
         CALL rate_diag(ndxohreact(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dpohreact(1:kproma, 1:plev, lat) = dpohreact(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

!!++mgs20140320: temporary addition - possibly irrelevant later
    IF (associated(dpln2o_j)) THEN
       CALL rate_diag(ndxn2o_j, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dpln2o_j(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dpln2o_o1d)) THEN
       CALL rate_diag(ndxn2o_o1d(1,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dpln2o_o1d(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxn2o_o1d(2,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dpln2o_o1d(1:kproma, 1:plev, lat) = dpln2o_o1d(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
    END IF
! s.stadtler
    IF (associated(dpppan)) THEN
       CALL rate_diag(ndxppan, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dpppan(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dprno2_oh)) THEN
       CALL rate_diag(ndxno2_oh, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dprno2_oh(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dprno2_ohreact)) THEN
       CALL rate_diag(ndxno2_ohreact, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dprno2_ohreact(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dplghno3)) THEN
       dplghno3(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreachno3
         CALL rate_diag(ndxlghno3(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dplghno3(1:kproma, 1:plev, lat) = dplghno3(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dpn2o5)) THEN
       dpn2o5(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreacn2o5
         CALL rate_diag(ndxn2o5(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dpn2o5(1:kproma, 1:plev, lat) = dpn2o5(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dplgn2o5)) THEN
       CALL rate_diag(ndxlgn2o5, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplgn2o5(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    
    IF (associated(dpno2_no3)) THEN
       CALL rate_diag(ndxno2_no3, plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dpno2_no3(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
    END IF

    IF (associated(dppno3)) THEN
       dppno3(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreacpno3 
         CALL rate_diag(ndxpno3(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dppno3(1:kproma, 1:plev, lat) = dppno3(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dplgno3)) THEN
       dplgno3(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreacno3 
         CALL rate_diag(ndxlgno3(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dplgno3(1:kproma, 1:plev, lat) = dplgno3(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dplnox_het)) THEN
       CALL rate_diag(ndxlnox_het(1,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplnox_het(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxlnox_het(2,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplnox_het(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxlnox_het(3,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplnox_het(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)
       CALL rate_diag(ndxlnox_het(4,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
       dplnox_het(1:kproma, 1:plev, lat) = zturnover(1:kproma, 1:plev)*0.5_dp
    END IF

    IF (associated(dppnox)) THEN
       dppnox(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreacpnox 
         CALL rate_diag(ndxpnox(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dppnox(1:kproma, 1:plev, lat) =  dppnox(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF

    IF (associated(dplnox)) THEN
       dplnox(1:kproma, 1:plev, lat) = 0._dp
       DO i = 1, nreaclnox 
         CALL rate_diag(ndxlnox(i,:), plonl, plev, pcnstm1, rxntot, vmr, reaction_rates, zturnover)
         dplnox(1:kproma, 1:plev, lat) = dplnox(1:kproma, 1:plev, lat) + zturnover(1:kproma, 1:plev)
       END DO
    END IF
! s.stadtler
  END SUBROUTINE moz_rates_diag

END MODULE mo_moz_diag
