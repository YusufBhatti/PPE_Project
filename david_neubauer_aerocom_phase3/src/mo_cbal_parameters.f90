!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cbal_parameters
  !
  !! Description:
  !!   Declaration of parameters used in cbalance
  ! 
  USE mo_kind, ONLY : dp

  IMPLICIT NONE 

  !! Public Parameters: 

  !! Parameters for landuse change

  REAL(dp), PARAMETER, private :: def_frac_wood_2_atmos    = 0.8_dp !! Fraction of carbon and nitrogen from the wood pool to be 
                                                                    !!   released into the atmosphere on human landcover change
  REAL(dp), PARAMETER, private :: def_frac_green_2_atmos   = 0.8_dp !! Fraction of carbon and nitrogen from the green pool to be 
                                                                    !!   released into the atmosphere on human landcover change 
  REAL(dp), PARAMETER, private :: def_frac_reserve_2_atmos = 0.8_dp !! Fraction of carbon from the reserve pool to be released 
                                                                    !!   into the atmosphere of human landcover change
  REAL(dp), PARAMETER, private :: def_frac_mobile_2_atmos  = 0.8_dp !! Fraction of nitrogen from the plant mobile N poolto be 
                                                                    !!   released into the atmosphere on human landcover change 
  REAL(dp), PARAMETER, private :: def_frac_harvest_2_atmos = 0.2_dp !! Fraction of harvested carbon released immediately to atm.


  !! Carbon-to-Nitrogen ratios for all carbon pools
 
  REAL(dp), PARAMETER, private :: def_cn_green = 35._dp !! For green pool (Ref: CLASS(50); CLM(av.CN=35); DGVM-LSM,Bonan,2003 (29)
  REAL(dp), PARAMETER, private :: def_cn_woods = 150._dp !! For wood pool (Ref: NCAR, CLM(CN = 50); 
                                                         !!                set to 150 by Daniel Goll and Daniela Kracher)
  REAL(dp), PARAMETER, private :: def_cn_litter_green = 55.0_dp  !! For green-litter (Ref: CASA; Biome-BGC)CN=150; CLM(av.CN=55)  
                                                        !!    Note that this is NOT the C/N-ratio of the green litter POOLS,
                                                        !!    but the C/N-ratio of FALLING LEAVES and FINE ROOTS. It determines
                                                        !!    only partly the C/N-ratio of the green litter pools, because they
                                                        !!    gain inputs from 2 pools (green,reserve) with different C/N-ratios
  REAL(dp), PARAMETER, private :: def_cn_litter_wood = 330.0_dp  !! For wood-litter pool (Ref: CASA) CN=338; CLM(CN =500); 
                                                                 !!   DGVM-LSM,Bonan,2003(330)
  REAL(dp), PARAMETER, private :: def_cn_slow = 10.0_dp          !! For slow pool  (Ref: Century=7-10);  CLM CN = 10

  !! 10% decay times [days] for anthropogenic c-pools (standard mode)
  REAL(dp), PARAMETER, private :: def_tau_onSite          =   1.0_dp*365._dp  ! anthro annual pool time scale
  REAL(dp), PARAMETER, private :: def_tau_paper           =  10.0_dp*365._dp  ! anthro deca.  pool time scale
  REAL(dp), PARAMETER, private :: def_tau_construction    = 100.0_dp*365._dp  ! anthro cent.  pool time scale
  
  !! Fractions for carbon to atmosphere

  REAL(dp) :: frac_wood_2_atmos
  REAL(dp) :: frac_green_2_atmos
  REAL(dp) :: frac_reserve_2_atmos
  REAL(dp) :: frac_mobile_2_atmos
  REAL(dp) :: frac_harvest_2_atmos

  !! Carbon-to-Nitrogen ratios for all carbon pools
 
  REAL(dp) :: cn_green
  REAL(dp) :: cn_woods
  REAL(dp) :: cn_litter_green
  REAL(dp) :: cn_litter_wood
  REAL(dp) :: cn_slow

  !! Constants for anthropogenic pools (tau_XXX is rescaled to division constants with unit [days]
  REAL(dp) :: frac_wood_2_onSite
  REAL(dp) :: frac_wood_2_paper
  REAL(dp) :: frac_wood_2_construction
  REAL(dp) :: tau_onSite
  REAL(dp) :: tau_paper
  REAL(dp) :: tau_construction
  REAL(dp) :: tau_max_green
  REAL(dp) :: tau_max_wood
  REAL(dp) :: tau_min_green
  REAL(dp) :: tau_min_wood

CONTAINS

  SUBROUTINE config_cbal_parameters(nml_unit)
  USE mo_mpi,       ONLY: p_parallel_io, p_parallel, p_io, p_bcast
  USE mo_namelist,  ONLY: position_nml, POSITIONED, LENGTH_ERROR, READ_ERROR
  USE mo_io_units,  ONLY: nout
  USE mo_exception, ONLY: finish, message
  INTEGER nml_unit

  INTEGER f_unit, read_status
  REAL (dp) :: rBuf(13)

  INCLUDE 'cbal_parameters_ctl.inc'

! Set default values

    frac_wood_2_atmos             = def_frac_wood_2_atmos    
    frac_green_2_atmos            = def_frac_green_2_atmos
    frac_reserve_2_atmos          = def_frac_reserve_2_atmos
    frac_mobile_2_atmos           = def_frac_mobile_2_atmos
    frac_harvest_2_atmos          = def_frac_harvest_2_atmos
    cn_green                      = def_cn_green
    cn_woods                      = def_cn_woods
    cn_litter_green               = def_cn_litter_green
    cn_litter_wood                = def_cn_litter_wood
    cn_slow                       = def_cn_slow             
    tau_onSite                    = def_tau_onSite
    tau_paper                     = def_tau_paper
    tau_construction              = def_tau_construction

! Read namelist if available
    IF (p_parallel_io) THEN 
      f_unit = position_nml ('CBAL_PARAMETERS_CTL', nml_unit,status=read_status)
      SELECT CASE (read_status)
        CASE (POSITIONED)
          READ (f_unit, cbal_parameters_ctl)
          CALL message('config_cbal_parameters','Namelist cbal_parameters_ctl: ')
          WRITE(nout,cbal_parameters_ctl)
        CASE (LENGTH_ERROR)
          CALL finish('config_cbal_parameters','Length error in namelist cbal_parameters_ctl')
        CASE (READ_ERROR)
          CALL finish('config_cbal_parameters','Error reading namelist cbal_parameters_ctl')
      END SELECT

    ENDIF

! Distribute parameters
   IF (p_parallel) THEN
     rbuf = (/frac_wood_2_atmos,frac_green_2_atmos,frac_reserve_2_atmos,frac_mobile_2_atmos,frac_harvest_2_atmos, &
              cn_green,cn_woods,cn_litter_green,cn_litter_wood,cn_slow,    &
              tau_onSite,tau_paper,tau_construction/)
     CALL p_bcast(rbuf,p_io)
     frac_wood_2_atmos             = rbuf( 1) 
     frac_green_2_atmos            = rbuf( 2)   
     frac_reserve_2_atmos          = rbuf( 3) 
     frac_mobile_2_atmos           = rbuf( 4)
     frac_harvest_2_atmos          = rbuf( 5) 
     cn_green                      = rbuf( 6) 
     cn_woods                      = rbuf( 7) 
     cn_litter_green               = rbuf( 8) 
     cn_litter_wood                = rbuf( 9) 
     cn_slow                       = rbuf(10) 
     tau_onSite                    = rbuf(11)
     tau_paper                     = rbuf(12)
     tau_construction              = rbuf(13)
   ENDIF

! Rescale anthropogenic tau-constants to 10% decay times
   tau_onSite       = 1._dp/(1._dp - 0.1_dp**(1._dp/tau_onSite      ))
   tau_paper        = 1._dp/(1._dp - 0.1_dp**(1._dp/tau_paper       ))
   tau_construction = 1._dp/(1._dp - 0.1_dp**(1._dp/tau_construction))

  END SUBROUTINE config_cbal_parameters

END MODULE mo_cbal_parameters
