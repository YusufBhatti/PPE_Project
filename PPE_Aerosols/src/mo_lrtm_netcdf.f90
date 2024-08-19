!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_lrtm_netcdf

  USE mo_mpi,       ONLY: p_parallel_io, p_io, p_bcast
  USE mo_netcdf,    ONLY: io_inq_varid, io_get_vara_double
  USE mo_io,        ONLY: io_open, io_close, io_read, file_info  
  USE mo_exception, ONLY: finish

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: lrtm_read

  INTEGER, PARAMETER :: maxAbsorberNameLength =  5, &
       Absorber              = 12

  CHARACTER(len = maxAbsorberNameLength), PARAMETER :: &
       AbsorberNames(Absorber) = (/        &
       'N2   ',  &
       'CCL4 ',  &
       'CFC11',  &
       'CFC12',  &
       'CFC22',  &
       'H2O  ',  &
       'CO2  ',  &
       'O3   ',  &
       'N2O  ',  & 
       'CO   ',  &
       'CH4  ',  &
       'O2   '  /)

  INTEGER, PARAMETER :: &
       keylower  = 9,   &
       keyupper  = 5,   &
       Tdiff     = 5,   &
       ps        = 59,  &
       plower    = 13,  &
       pupper    = 47,  &
       Tself     = 10,  &
       Tforeign  = 4,   &
       pforeign  = 4,   &
       T         = 19,  &
       Tplanck   = 181, &
       band      = 16,  &
       GPoint    = 16,  &
       GPointSet = 2

  INTEGER, PARAMETER :: gPointSetNumber = 1
  INTEGER :: varid

  TYPE(file_info) :: rrtmg_lw

CONTAINS 

  !=============================================================================

  INTEGER FUNCTION AbsorberIndex(AbsorberName)
    CHARACTER(len=*), INTENT(in)  :: AbsorberName

    INTEGER :: m

    AbsorberIndex = -1
    DO m = 1, Absorber
      IF (TRIM(AbsorberNames(m)) == TRIM(AbsorberName)) THEN
        AbsorberIndex = m
      END IF
    END DO

    IF (AbsorberIndex == -1) THEN
      CALL finish('Absorber name index lookup failed.')
    END IF

  END FUNCTION AbsorberIndex

  !=============================================================================

  SUBROUTINE lrtm_read

    USE rrlw_planck, ONLY: chi_mls, totplanck, totplanck16

  
    IF (p_parallel_io) THEN
      CALL io_open ('rrtmg_lw.nc', rrtmg_lw, io_read)
    ENDIF

	IF (p_parallel_io) THEN
 
	  CALL io_inq_varid(rrtmg_lw%file_id, 'AbsorberAmountMLS', varid)
	  CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
		   (/6, 1/), &
		   (/SIZE(chi_mls, 1), SIZE(chi_mls, 2) /), &
		   chi_mls)
	  CALL io_inq_varid(rrtmg_lw%file_id, 'IntegratedPlanckFunction', varid)
	  CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
		   (/1, 1/), &
		   (/SIZE(totplanck,1), SIZE(totplanck,2)/), &
		   totplanck)
	  CALL io_inq_varid(rrtmg_lw%file_id, 'IntegratedPlanckFunctionBand16', varid)
	  CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
		   (/ 1 /), &
		   (/SIZE(totplanck16)/), &
		   totplanck16)
	END IF 
    CALL p_bcast(chi_mls,     p_io) 
    CALL p_bcast(totplanck,   p_io) 
    CALL p_bcast(totplanck16, p_io) 
	
    CALL lw_kgb01  ! molecular absorption coefficients
    CALL lw_kgb02
    CALL lw_kgb03
    CALL lw_kgb04
    CALL lw_kgb05
    CALL lw_kgb06
    CALL lw_kgb07
    CALL lw_kgb08
    CALL lw_kgb09
    CALL lw_kgb10
    CALL lw_kgb11
    CALL lw_kgb12
    CALL lw_kgb13
    CALL lw_kgb14
    CALL lw_kgb15
    CALL lw_kgb16

    IF (p_parallel_io) THEN
      CALL io_close(rrtmg_lw)
    ENDIF

  END SUBROUTINE lrtm_read 

    SUBROUTINE lw_kgb01

      USE rrlw_kg01, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, selfrefo, forrefo, no1

      INTEGER, PARAMETER :: bandNumber = 1
      INTEGER, PARAMETER :: numGPoints = no1

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id, 'PlanckFractionLowerAtmos', varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kao_mn2)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kbo_mn2)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(fracrefbo, p_io)
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kbo,       p_io)
      CALL p_bcast(kao_mn2,   p_io)
      CALL p_bcast(kbo_mn2,   p_io) 
      CALL p_bcast(selfrefo,  p_io)
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb01

    SUBROUTINE lw_kgb02

      USE rrlw_kg02, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no2

      INTEGER, PARAMETER :: bandNumber = 2
      INTEGER, PARAMETER :: numGPoints = no2

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

      ENDIF

      CALL p_bcast(fracrefao, p_io)
      CALL p_bcast(fracrefbo, p_io)
      CALL p_bcast(kao,       p_io)
      CALL p_bcast(kbo,       p_io)
      CALL p_bcast(selfrefo,  p_io)
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb02

    SUBROUTINE lw_kgb03

      USE rrlw_kg03, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, selfrefo, forrefo, no3

      INTEGER, PARAMETER :: bandNumber = 3
      INTEGER, PARAMETER :: numGPoints = no3

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keyupper,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keyupper,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
             (/keylower,T,numGPoints,1,1,1/), &
             kao_mn2o)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
             (/keyupper,T,numGPoints,1,1,1/), &
             kbo_mn2o)

      ENDIF

      CALL p_bcast(fracrefao, p_io)
      CALL p_bcast(fracrefbo, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kbo,       p_io) 
      CALL p_bcast(kao_mn2o,  p_io)
      CALL p_bcast(kbo_mn2o,  p_io) 
      CALL p_bcast(selfrefo,  p_io)
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb03

    SUBROUTINE lw_kgb04

      USE rrlw_kg04, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no4

      INTEGER, PARAMETER :: bandNumber = 4
      INTEGER, PARAMETER :: numGPoints = no4

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keyupper,1,1/), &
             fracrefbo(:,1:5))

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keyupper,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

      ENDIF

      CALL p_bcast(fracrefao, p_io)
      CALL p_bcast(fracrefbo, p_io)
      CALL p_bcast(kao,       p_io)
      CALL p_bcast(kbo,       p_io)
      CALL p_bcast(selfrefo,  p_io)
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb04

    SUBROUTINE lw_kgb05

      USE rrlw_kg05, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mo3, selfrefo, forrefo, ccl4o, no5

      INTEGER, PARAMETER :: bandNumber = 5
      INTEGER, PARAMETER :: numGPoints = no5

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keyupper,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, & 
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keyupper,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('O3'),bandNumber,gPointSetNumber/), &
             (/keylower,T,numGPoints,1,1,1/), &
             kao_mo3)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CCL4'),bandNumber,gPointSetNumber/), &
             (/1,1,numGPoints,1,1,1/), &
             ccl4o)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(fracrefbo, p_io)
      CALL p_bcast(kao,       p_io)
      CALL p_bcast(kbo,       p_io)
      CALL p_bcast(kao_mo3,   p_io)
      CALL p_bcast(selfrefo,  p_io)
      CALL p_bcast(forrefo,   p_io)
      CALL p_bcast(ccl4o,     p_io)

    END SUBROUTINE lw_kgb05

    SUBROUTINE lw_kgb06

      USE rrlw_kg06, ONLY : fracrefao, kao, kao_mco2, selfrefo, forrefo, cfc11adjo, cfc12o, no6

      INTEGER, PARAMETER :: bandNumber = 6
      INTEGER, PARAMETER :: numGPoints = no6

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kao_mco2)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CFC11'),bandNumber,gPointSetNumber/), &
             (/1,1,numGPoints,1,1,1/), &
             cfc11adjo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CFC12'),bandNumber,gPointSetNumber/), &
             (/1,1,numGPoints,1,1,1/), &
             cfc12o)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kao_mco2,  p_io) 
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io) 
      CALL p_bcast(cfc11adjo, p_io) 
      CALL p_bcast(cfc12o,    p_io)

    END SUBROUTINE lw_kgb06

    SUBROUTINE lw_kgb07 

      USE rrlw_kg07, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, selfrefo, forrefo, no7

      INTEGER, PARAMETER :: bandNumber = 7
      INTEGER, PARAMETER :: numGPoints = no7

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, & 
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
             (/keylower,T,numGPoints,1,1,1/), &
             kao_mco2)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kbo_mco2)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(fracrefbo, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kbo,       p_io) 
      CALL p_bcast(kao_mco2,  p_io) 
      CALL p_bcast(kbo_mco2,  p_io)
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb07

    SUBROUTINE lw_kgb08         

      USE rrlw_kg08, ONLY : fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, kao_mo3, kbo, kbo_mco2, kbo_mn2o, &
           selfrefo, forrefo, cfc12o, cfc22adjo, no8

      INTEGER, PARAMETER :: bandNumber = 8
      INTEGER, PARAMETER :: numGPoints = no8

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('O3'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kao_mo3)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kao_mco2)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kbo_mco2)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kao_mn2o)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kbo_mn2o)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CFC12'),bandNumber,gPointSetNumber/), &
             (/1,1,numGPoints,1,1,1/), &
             cfc12o)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CFC22'),bandNumber,gPointSetNumber/), &
             (/1,1,numGPoints,1,1,1/), &
             cfc22adjo)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(fracrefbo, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kao_mco2,  p_io) 
      CALL p_bcast(kao_mn2o,  p_io) 
      CALL p_bcast(kao_mo3,   p_io)
      CALL p_bcast(kbo,       p_io) 
      CALL p_bcast(kbo_mco2,  p_io)
      CALL p_bcast(kbo_mn2o,  p_io) 
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io) 
      CALL p_bcast(cfc12o,    p_io) 
      CALL p_bcast(cfc22adjo, p_io)

    END SUBROUTINE lw_kgb08

    SUBROUTINE lw_kgb09 

      USE rrlw_kg09, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, selfrefo, forrefo, no9

      INTEGER, PARAMETER :: bandNumber = 9
      INTEGER, PARAMETER :: numGPoints = no9

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
             (/keylower,T,numGPoints,1,1,1/), &
             kao_mn2o)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2O'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kbo_mn2o)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(fracrefbo, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kbo,       p_io) 
      CALL p_bcast(kao_mn2o,  p_io) 
      CALL p_bcast(kbo_mn2o,  p_io) 
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb09

    SUBROUTINE lw_kgb10         

      USE rrlw_kg10, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no10

      INTEGER, PARAMETER :: bandNumber = 10
      INTEGER, PARAMETER :: numGPoints = no10

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(fracrefbo, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kbo,       p_io) 
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb10

    SUBROUTINE lw_kgb11
      USE rrlw_kg11, ONLY : fracrefao, fracrefbo, kao, kbo, kao_mo2, kbo_mo2, selfrefo, forrefo, no11

      INTEGER, PARAMETER :: bandNumber = 11
      INTEGER, PARAMETER :: numGPoints = no11

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)


        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('O2'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kao_mo2)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('O2'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kbo_mo2)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(fracrefbo, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kbo,       p_io) 
      CALL p_bcast(kao_mo2,   p_io) 
      CALL p_bcast(kbo_mo2,   p_io) 
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb11

    SUBROUTINE lw_kgb12

      USE rrlw_kg12, ONLY : fracrefao, kao, selfrefo, forrefo, no12

      INTEGER, PARAMETER :: bandNumber = 12
      INTEGER, PARAMETER :: numGPoints = no12

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb12

    SUBROUTINE lw_kgb13         

      USE rrlw_kg13, ONLY : fracrefao, fracrefbo, kao, kao_mco2, kao_mco, kbo_mo3, selfrefo, forrefo, no13

      INTEGER, PARAMETER :: bandNumber = 13
      INTEGER, PARAMETER :: numGPoints = no13  

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/),  &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('O3'),bandNumber,gPointSetNumber/), &
             (/1,T,numGPoints,1,1,1/), &
             kbo_mo3)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CO2'),bandNumber,gPointSetNumber/), &
             (/keylower,T,numGPoints,1,1,1/), &
             kao_mco2)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('CO'),bandNumber,gPointSetNumber/), &
             (/keylower,T,numGPoints,1,1,1/), &
             kao_mco)

      ENDIF

      CALL p_bcast(fracrefao, p_io)
      CALL p_bcast(fracrefbo, p_io)
      CALL p_bcast(kao,       p_io)
      CALL p_bcast(kao_mco2,  p_io)
      CALL p_bcast(kao_mco,   p_io)
      CALL p_bcast(kbo_mo3,   p_io) 
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb13

    SUBROUTINE lw_kgb14 

      USE rrlw_kg14, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no14

      INTEGER, PARAMETER :: bandNumber = 14
      INTEGER, PARAMETER :: numGPoints = no14

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

      ENDIF

      CALL p_bcast(fracrefao, p_io)
      CALL p_bcast(fracrefbo, p_io)
      CALL p_bcast(kao,       p_io)
      CALL p_bcast(kbo,       p_io)
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb14

    SUBROUTINE lw_kgb15 

      USE rrlw_kg15, ONLY : fracrefao, kao, kao_mn2, selfrefo, forrefo, no15

      INTEGER, PARAMETER :: bandNumber = 15
      INTEGER, PARAMETER :: numGPoints = no15

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'AbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,AbsorberIndex('N2'),bandNumber,gPointSetNumber/), &
             (/keylower,T,numGPoints,1,1,1/), &
             kao_mn2)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kao_mn2,   p_io)
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb15

    SUBROUTINE lw_kgb16         

      USE rrlw_kg16, ONLY : fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, no16

      INTEGER, PARAMETER :: bandNumber = 16
      INTEGER, PARAMETER :: numGPoints = no16

      IF (p_parallel_io) THEN

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,keylower,1,1/), &
             fracrefao)

        CALL io_inq_varid(rrtmg_lw%file_id,'PlanckFractionUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/numGPoints,1,1,1/), &
             fracrefbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsLowerAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/keylower,Tdiff,plower,numGPoints,1,1/), &
             kao)

        CALL io_inq_varid(rrtmg_lw%file_id,'KeySpeciesAbsorptionCoefficientsUpperAtmos',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,1,1,bandNumber,gPointSetNumber/), &
             (/1,Tdiff,pupper,numGPoints,1,1/), &
             kbo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20SelfAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tself,numGPoints,1,1/), &
             selfrefo)

        CALL io_inq_varid(rrtmg_lw%file_id,'H20ForeignAbsorptionCoefficients',varid)
        CALL io_get_vara_double(rrtmg_lw%file_id, varid, &
             (/1,1,bandNumber,gPointSetNumber/), &
             (/Tforeign,numGPoints,1,1/), &
             forrefo)

      ENDIF

      CALL p_bcast(fracrefao, p_io) 
      CALL p_bcast(fracrefbo, p_io) 
      CALL p_bcast(kao,       p_io) 
      CALL p_bcast(kbo,       p_io) 
      CALL p_bcast(selfrefo,  p_io) 
      CALL p_bcast(forrefo,   p_io)

    END SUBROUTINE lw_kgb16

END MODULE mo_lrtm_netcdf
