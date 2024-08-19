!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! ham_salsactl
!!
!! \brief
!! namelist for the SALSA aerosol model
!!
!! \author Anton Laakso (FMI)
!!
!! \responsible_coder
!! Anton Laakso, anton.laakso@fmi.fi
!!
!! \revision_history
!!   -# Anton Laakso (FMI) +SALSA (2013-6)
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
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

NAMELIST /ham_salsactl/ nsnucl,     & !< Choice of the sulfate aerosol nucleation scheme:
                                      !!  nsnucl = 1  Binary
                                      !!         = 2  activation type nucleation
                                      !!         = 3  Kinetic
                                      !!         = 4  Ternary
                                      !!         = 5  nucleation with ORGANICs
                                      !!         = 6  activation type of nucleation with H2SO4+ORG
                                      !!         = 7  heteromolecular nucleation with H2SO4*ORG
                                      !!         = 8  homomolecular nucleation of  H2SO4 + 
                                      !!              heteromolecular nucleation with H2SO4*ORG
                                      !!         = 9  homomolecular nucleation of  H2SO4 and ORG + 
                                      !!              heteromolecular nucleation with H2SO4*ORG
                        locgas,     & !< Organic carbon emisison in gas phase 
                        lsol2b,     & !< Repartition soluble material from b-regions to a-regions
                        nj3,        & !< Choice of the particle formation scheme
                                      !!
                                      !!    j3_scheme = 1 Kerminen and Kulmala
                                      !!              = 2 Lehtinen et al. (2007)
                                      !!              = 3 Anttila et al. (2010)
                        act_coeff     !< Activation coefficient [unit?] 
