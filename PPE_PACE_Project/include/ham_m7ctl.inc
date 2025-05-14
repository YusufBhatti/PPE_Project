!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! ham_m7ctl
!!
!! \brief
!! namelist for the M7 aerosol model
!!
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Philip Stier, philip.stier@physics.ox.ac.uk
!!
!! \revision_history
!!   -# Philip Stier (MPI-Met) - original code (2003-01)
!!   -# Jan Kazil (MPI-M) (2008-05)
!!
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

NAMELIST /ham_m7ctl/  nwater,     & !< Aerosol water uptake scheme:
                                    !!
                                    !! nwater = 0 Jacobson et al., JGR 1996
                                    !!        = 1 Kappa-Koehler theory based approach (Petters and Kreidenweis, ACP 2007)
                      nsnucl,     & !< Choice of the sulfate aerosol nucleation scheme:
                                    !! 
                                    !!  nsnucl = 0 off
                                    !!         = 1 Vehkamaeki et al., JGR 2002
                                    !!         = 2 Kazil and Lovejoy, ACP 2007
                                    !!
                      nonucl,     & !< Choice of the organic aerosol nucleation scheme:
                                    !! 
                                    !!  nonucl = 0 off
                                    !!         = 1 Activation nucleation, Kulmala et al., ACP 2006
                                    !!         = 2 Kinetic nucleation, Laakso et al., ACP 2004
                                    !!
                      lnucl_stat    !< Sample the cloud-free volume as function of T, RH, [H2SO4(g)],
                                    !! H2SO4 condensation sink, and ionization rate (memory intensive)
