!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module containing the submodel for interactive CO2 calculation
!!
!! @remarks This module contains all subroutines to define and transport
!!   interactive CO2 as a tracer in echam
!! @author Reiner Schnur, MPI-M, Hamburg (2004-12)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   MPI, November-December 2004: original implementation (L. Kornblueh,
!!   R. Schnur)
!!   November 2004, added some diagnostics and new stream variables for 
!!   land/ocean fluxes (R. Schnur)
!!   December 2004, added routines to use anthropogenic CO2 
!!   emissions (R. Schnur)
!!   December 2009, adapted to echam6 (S. Rast)
!!   January  2010, substitute cpl_co2-flags by namelist switch 
!!                  lcouple_co2 (M. Esch)
!!
!

MODULE mo_co2


  USE mo_kind,          ONLY: wp
  USE mo_tracer,        ONLY: new_tracer, get_tracer
  USE mo_tracdef,       ONLY: OFF, ON, RESTART, CONSTANT
  USE mo_memory_base,   ONLY: new_stream, default_stream_setting,   &
                              add_stream_element, delete_stream,    &
                              GAUSSIAN
  USE mo_linked_list,   ONLY: t_stream
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_mpi,           ONLY: p_parallel_io, p_io, p_bcast
  USE mo_exception,     ONLY: message_text, message, finish
  USE mo_control,       ONLY: lcouple, lcouple_co2
  USE mo_submodel,      ONLY: lco2, print_value
  USE mo_physical_constants, ONLY: amco2, amd      ! molecular weight of CO2 [g mol-1],
                                                   ! molecular weight of dry air [g mol-1]

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_co2,             init_co2_memory,      cleanup_co2_memory,   &
            reference_co2,        co2_mbalance,         co2_te_check,         &
            diag_co2,             read_co2_emission,    co2_exchange,         &
            co2_flux_correction
  PUBLIC :: co2_flux_atmosphere_ocean

  LOGICAL, PUBLIC           :: lco2_flxcor=.TRUE.    !flux correction
  LOGICAL, PUBLIC           :: lco2_mixpbl=.TRUE.    !ideal mixing of CO2 in PBL
  LOGICAL, PUBLIC           :: lco2_2perc=.FALSE.    !2% limitation of tendency
  LOGICAL, PUBLIC           :: lco2_emis=.FALSE.     !emissions from file
  LOGICAL, PUBLIC           :: lco2_clim=.FALSE.     !use transported CO2
  LOGICAL, PUBLIC           :: lco2_scenario=.FALSE. !use CO2 concentrations from external GHG 
                                                     ! scenario (if lco2=false and ighg=1)
  REAL(wp), POINTER, PUBLIC :: co2(:,:,:)
  REAL(wp), POINTER, PUBLIC :: co2m1(:,:,:)
  REAL(wp), POINTER, PUBLIC :: co2atmos(:,:)         ! CO2 concentration in lowest level, accu-
                                                     ! mulated for coupling (collect.f90)
  REAL(wp), POINTER, PUBLIC :: co2flux_cpl(:,:)      ! upward ocean CO2 flux, accumulated for 
                                                     ! coupling (collect.f90)
  REAL(wp), POINTER, PUBLIC :: co2trans(:,:)         ! CO2 ocean-atmosphere transfer velocity (not
                                                     ! including wind speed)
  REAL(wp), POINTER, PUBLIC :: co2ocean(:,:)         ! CO2 concentration of the ocean
  REAL(wp), POINTER, PUBLIC :: co2_flux_ocean(:,:)   ! CO2 flux from ocean
  REAL(wp), POINTER, PUBLIC :: co2_flux(:,:)         ! Total natural CO2 flux from ocean and land 
                                                     !  (all except emissions and flux correction)
  REAL(wp), POINTER, PUBLIC :: co2_flux_land(:,:)    ! Total natural CO2 flux from land (all except
                                                     !   emissions and flux corr.)
  REAL(wp), POINTER, PUBLIC :: co2_flux_npp(:,:)     ! CO2 flux from NPP
  REAL(wp), POINTER, PUBLIC :: co2_flux_soilresp(:,:)  ! CO2 flux from soil respiration
  REAL(wp), POINTER, PUBLIC :: co2_flux_herbivory(:,:) ! CO2 flux from herbivory
  REAL(wp), POINTER, PUBLIC :: co2_flux_dynveg(:,:)    ! CO2 flux from fire (dynveg)
  REAL(wp), POINTER, PUBLIC :: co2_flux_acc(:,:)       ! Same as above but accumulated
  REAL(wp), POINTER, PUBLIC :: co2_flux_ocean_acc(:,:) ! ...
  REAL(wp), POINTER, PUBLIC :: co2_flux_land_acc(:,:)  ! ...
  REAL(wp), POINTER, PUBLIC :: co2_flux_npp_acc(:,:)   ! ...
  REAL(wp), POINTER, PUBLIC :: co2_flux_soilresp_acc(:,:)  ! ...
  REAL(wp), POINTER, PUBLIC :: co2_flux_herbivory_acc(:,:) ! ...
  REAL(wp), POINTER, PUBLIC :: co2_flux_dynveg_acc(:,:)    ! ...
  REAL(wp), POINTER, PUBLIC :: co2_burden(:,:)
  REAL(wp), POINTER, PUBLIC :: co2_burden_old(:,:)
  REAL(wp), POINTER, PUBLIC :: co2_burden_acc(:,:)
  REAL(wp), POINTER, PUBLIC :: co2_flux_corr_acc(:,:)      ! CO2 flux correction, accumulated for 
                                                           !  output time step
  REAL(wp), POINTER, PUBLIC :: co2_burden_corr_acc1(:,:)   ! CO2 burden correction, accumulated for
                                                           !   daily correction time step
  ! Note that co2_flux_corr_acc1 and co2_flux_corr_acc1_old are globally constant, 
  !  but are kept in spatial fields because they need to be written/read to/from the restart files.
  REAL(wp), POINTER, PUBLIC :: co2_flux_corr_acc1(:,:)     ! CO2 flux correction, accumulated for 
                                                           !  daily correction time step
  REAL(wp), POINTER, PUBLIC :: co2_flux_corr_acc1_old(:,:) ! CO2 flux correction, accumulated for 
                                                           !  daily correction time step
  REAL(wp), POINTER, PUBLIC :: co2_burden_corr_acc2(:,:)   ! CO2 burden correction in tendencies, 
                                                           !  accumulated for output time step

  ! co2_emission is read from input file for scenario simulations with prescribed emissions
  REAL(wp), POINTER, PUBLIC :: co2_emission(:,:)           ! CO2 emissions from fossil fuel 
                                                           !  combustion (including sources like 
                                                           !  cement
  REAL(wp), POINTER, PUBLIC :: co2_emission_acc(:,:)       ! production, flaring, as specified by 
                                                           !  scenario) .. as read from input file

  ! co2_emission_lcc and co2_emission_harvest are computed in JSBACH from scenarios 
  !  of land cover change
  REAL(wp), POINTER, PUBLIC :: co2_emission_lcc(:,:)       ! CO2 emission from land use change
  REAL(wp), POINTER, PUBLIC :: co2_emission_harvest(:,:)   ! CO2 emission from harvest
  REAL(wp), POINTER, PUBLIC :: co2_emission_lcc_acc(:,:)   ! Same as above but accumulated
  REAL(wp), POINTER, PUBLIC :: co2_emission_harvest_acc(:,:) ! ...

  ! The following emissions are diagnosed (also for concentration-driven scenario simulations)
  REAL(wp), POINTER, PUBLIC :: co2_flux_anthro_acc(:,:)    ! Anthropogenic CO2 emissions (fossil 
                                                           !  and land use change/harvest)

  REAL(wp), POINTER, PUBLIC :: co2_flux_total(:,:)         ! Total CO2 Flux (incl. emissions, but 
                                                           !  not flux correction)

  INTEGER :: emis_no_years = -1
  REAL(wp) :: emis_base_year
  REAL(wp), ALLOCATABLE :: emis_years(:)
  INTEGER :: emis_id = -1  ! NetCDF file id

  INTEGER, PUBLIC :: ico2idx
  TYPE(t_stream), POINTER :: xtco2

  LOGICAL, SAVE, PUBLIC :: l_co2flxcorr = .FALSE.

CONTAINS

  !---------------------------------------------------------------------------
  !>
  !! @brief Initializes the interactive CO2 submodel
  !! 
  !! @remarks This routine reads the CO2 name list and 
  !! initializes the CO2 tracer
  !

  SUBROUTINE init_co2(co2mmr)

    USE mo_mpi,                  ONLY: p_parallel_io, p_parallel, p_bcast
    USE mo_namelist,             ONLY: open_nml, position_nml,  &
                                       POSITIONED,              &
                                       LENGTH_ERROR, READ_ERROR
    USE mo_time_control,         ONLY: lstart
    USE mo_radiation_parameters, ONLY: ighg

    REAL(wp), INTENT(in) :: co2mmr

    INTEGER :: ierr
    CHARACTER(LEN=32)    :: trsubname      ! subname ('_clim') for climatological CO2 tracer
    CHARACTER(LEN=256)   :: trlongname     ! longname for CO2 tracer
    CHARACTER(LEN=256)   :: cmodule        !< module name

    INTEGER :: inml, iunit

    INCLUDE 'co2ctl.inc'

    IF (.NOT. lco2 .AND. ighg == 1) lco2_scenario = .TRUE.

    IF (p_parallel_io) THEN

      inml = open_nml('namelist.echam')
      iunit =  position_nml ('CO2CTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
        READ (iunit, co2ctl)
      CASE (LENGTH_ERROR)
        CALL finish ('init_submodel_co2', &
             'length error in namelist co2ctl')
      CASE (READ_ERROR)
        CALL finish ('init_submodel_co2', &   
             'read error in namelist.echam')
      END SELECT
    END IF

    IF (.NOT. lco2) THEN
       IF (p_parallel) CALL p_bcast (lco2_scenario, p_io)
       IF (lco2_scenario) THEN
          IF (ighg == 1) THEN
             CALL message('init_co2', 'using CO2 from external GHG scenario.')
          ELSE
             CALL finish('init_co2', 'Need ighg=1 to use CO2 '// &
                         'concentrations from external scenario')
          END IF
       ELSE
          CALL print_value('init_co2, using fixed global CO2 mass mixing ratio'// &
                           ' co2mmr=',co2mmr)
       END IF
       RETURN
    END IF
    CALL message('init_co2', 'Using interactive CO2 transport')
    IF (lco2_scenario) THEN
       CALL message('init_co2', 'Setting lco2_scenario to .false. since interactive '// &
                                'CO2 is used (lco2=.true.)')
       lco2_scenario = .FALSE.
    END IF

    IF (p_parallel) THEN
       CALL p_bcast (lco2_flxcor, p_io)
       CALL p_bcast (lco2_mixpbl, p_io)
       CALL p_bcast (lco2_2perc, p_io)
       CALL p_bcast (lco2_emis, p_io)
       CALL p_bcast (lco2_clim, p_io)
    END IF
    IF (p_parallel_io) THEN
       CALL print_value('init_submodel_co2, flux correction, lco2_flxcor', &
                        lco2_flxcor)
       CALL print_value('init_submodel_co2, mix CO2 in PBL, lco2_mixpbl', &
                        lco2_mixpbl)
       CALL print_value('init_submodel_co2, limit tendencies, lco2_2perc', &
                        lco2_2perc)
       CALL print_value('init_submodel_co2, use external co2 emissions,'// &
                        ' lco2_emis',lco2_emis)
       CALL print_value('init_submodel_co2, use climatological CO2,'// &
                        ' lco2_clim',lco2_clim)
    END IF

    IF (lco2_clim) THEN
      trsubname = 'clim'
      trlongname = 'mass fraction of CO2 in air (climatology)'
    ELSE
      trsubname = ''
      trlongname = 'mass fraction of CO2 in air'
    END IF
    CALL get_tracer ('CO2',subname=TRIM(trsubname), &
                     modulename=cmodule, idx=ico2idx, ierr=ierr)
    IF (ierr==1) THEN
      CALL new_tracer ('CO2', 'CO2',                 &
           subname = TRIM(trsubname),                &
           longname = TRIM(trlongname),              &
           idx = ico2idx,                            &
           units = 'kg kg-1',                        &
           code = 1, table = 191, bits = 16,         &
           ninit = RESTART+CONSTANT, vini = co2mmr,  &
           nwrite = ON, nrerun = ON, nemis = OFF,    &
           ierr = ierr)
    ELSE
      CALL finish ('mo_co2 ( init_co2 )','Tracer ''CO2'' already exists '//&
                   'defined by module '//TRIM(cmodule))
    END IF

    IF (p_parallel_io .AND. lstart) THEN
       CALL print_value('init_co2, global CO2 mass mixing ratio'// &
                        ' initialized to co2mmr=',co2mmr)
    END IF

  END SUBROUTINE init_co2

  !---------------------------------------------------------------------------
  !>
  !! @brief Initializes the memory for the CO2 submodel
  !! 
  !! @remarks The CO2 stream is defined and written for interactive CO2 as well
  !! as for the case that CO2 is not transported. In the latter case, one still
  !! wants to have the diagnostics to look at the CO2 cycle.
  !

  SUBROUTINE init_co2_memory(co2mmr)

    USE mo_filename,      ONLY: out_filetype                                            

    REAL(wp), INTENT(in) :: co2mmr

    INTEGER :: nlon, nglat, lnlon, lnglat 

    LOGICAL :: lpost = .TRUE.
    LOGICAL :: lrerun = .TRUE.

    nlon   = lc%nlon
    nglat  = lc%nlat
    lnlon  = lc%nproma
    lnglat = lc%ngpblks

    CALL new_stream (xtco2, 'co2', filetype=out_filetype,  &
         post_suf='_co2', rest_suf='_co2') 

    CALL default_stream_setting (xtco2,     &
         table=191, bits=16, repr=GAUSSIAN, &
         lrerun=lrerun, lpost=lpost)

    IF (.NOT. lco2) THEN
       ALLOCATE(co2m1(lnlon, lc%nlev, lnglat))
       co2m1 = co2mmr ! If lco2_scenario=.true. this will be overwritten each
                      ! time step by values from greenhouse gas scenario
    END IF

    ! Instantaneous CO2 flux from the surface (ocean+land) into the atmosphere
    ! This will be updated every time step for the land surface, but only
    ! at each coupling time step for the ocean
    IF (lcouple) THEN
       CALL add_stream_element (xtco2, 'co2atmos', co2atmos,                 &
            (/ lnlon, lnglat /), (/ nlon, nglat /),                          &
            code=16,                                                         &
            longname='CO2 concentration at surface for ocean coupling',      &
            units='kg kg-1', laccu=.FALSE., lpost=.FALSE.,                   &
            ! is accumulated in collect at coupling interval
            lrerun=lrerun, contnorest=.TRUE.) 
    END IF
    CALL add_stream_element (xtco2, 'co2_flx_inst', co2_flux,                &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=2, bits=24,                                                    &
         longname='upward surface CO2 flux',                                 &
         units='kg m-2 s-1', laccu=.FALSE., lrerun=.FALSE.)
    CALL add_stream_element (xtco2, 'co2_flx_l_inst', co2_flux_land,         &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=3, bits=24,                                                    &
         longname='upward land CO2 flux',                                    &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
    CALL add_stream_element (xtco2, 'co2_flx_npp_inst', co2_flux_npp,        &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=27, bits=24,                                                   &
         longname='upward CO2 flux - NPP',                                   &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
    CALL add_stream_element (xtco2, 'co2_flx_resp_inst', co2_flux_soilresp,  &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=28, bits=24,                                                   &
         longname='upward CO2 flux - soil respiration',                      &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
    CALL add_stream_element (xtco2, 'co2_flx_herb_inst', co2_flux_herbivory, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=29, bits=24,                                                   &
         longname='upward CO2 flux - herbivory',                             &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
    CALL add_stream_element (xtco2, 'co2_flx_lcc_inst', co2_emission_lcc,    &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=30, bits=24,                                                   &
         longname='upward CO2 flux - land use change',                       &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
    CALL add_stream_element (xtco2, 'co2_flx_harv_inst', co2_emission_harvest,&
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=31, bits=24,                                                   &
         longname='upward CO2 flux - harvest',                               &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
    CALL add_stream_element (xtco2, 'co2_flx_fire_inst', co2_flux_dynveg,    &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=32, bits=24,                                                   &
         longname='upward CO2 flux - fire',                                  &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
    CALL add_stream_element (xtco2, 'co2_flx_o_inst', co2_flux_ocean,        &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=4, bits=24,                                                    &
         longname='upward ocean CO2 flux',                                   &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE., lrerun=.FALSE.)
   CALL add_stream_element (xtco2, 'co2flx_cpl', co2flux_cpl,                &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=19,                                                            &
         longname='upward ocean CO2 flux for ocean coupling',                &
         units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE.,                   &
         ! is accumulated in collect at coupling interval
         lrerun=lrerun, contnorest=.TRUE.)
   IF (lcouple_co2) THEN 
      CALL add_stream_element (xtco2, 'co2trans', co2trans,                  &
           (/ lnlon, lnglat /), (/ nlon, nglat /),                           &
           code=18, bits=24,                                                 &
           longname='transfer velocity of ocean/atmosphere CO2 flux',        &
           units='10-9 mol s m-4', laccu=.FALSE., lpost=.TRUE.)
      CALL add_stream_element (xtco2, 'co2ocean', co2ocean,                  &
           (/ lnlon, lnglat /), (/ nlon, nglat /),                           &
           code=17,bits=24,                                                  &
           longname='pco2 of surface ocean water',                           &
           units='ppm CO2', laccu=.FALSE., lpost=.TRUE.)
   ENDIF
   CALL add_stream_element (xtco2, 'co2_flx_total_inst', co2_flux_total,     &
        (/ lnlon, lnglat /), (/ nlon, nglat /),                              &
        code=15, bits=24,                                                    &
        longname='total upward surface CO2 flux',                            &
        units='kg m-2 s-1', laccu=.FALSE., lpost=lpost .AND. lco2,           &
        contnorest=.TRUE.) 

   CALL add_stream_element (xtco2, 'co2_burden_inst', co2_burden,            &
        (/ lnlon, lnglat /), (/ nlon, nglat /),                              &
        code=14,                                                             &
        longname='CO2 content',                                              &
        units='kg m-2', laccu=.FALSE., lrerun=lrerun, contnorest=.TRUE.)
   CALL add_stream_element (xtco2, 'co2_burden_old', co2_burden_old,         &
        (/ lnlon, lnglat /), (/ nlon, nglat /),                              &
        longname='old CO2 content',                                          &
        units='kg m-2', laccu=.FALSE., lpost=.FALSE., lrerun=lrerun,         &
        contnorest=.TRUE.)

   IF (lco2_emis) THEN
      CALL add_stream_element (xtco2, 'co2_emis_inst', co2_emission,         &
           (/ lnlon, lnglat /), (/ nlon, nglat /),                           &
           code=12, bits=24,                                                 &
           longname='upward CO2 emission',                                   &
           units='kg m-2 s-1', laccu=.FALSE., lpost=.FALSE.,                 &
           lrerun=lrerun, contnorest=.TRUE.)
      ! lrerun must be true because emissions are read only at beginning of new year
   END IF

    ! Accumulated CO2 flux
    CALL add_stream_element (xtco2, 'co2_flux', co2_flux_acc,                &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=5, bits=24,                                                    &
         longname='total upward surface CO2 flux (acc.)',                    &
         units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flx_land', co2_flux_land_acc,       &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=6, bits=24,                                                    &
         longname='total upward land CO2 flux (acc.)',                       &
         units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flx_npp', co2_flux_npp_acc,         &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=21, bits=24,                                                   &
         longname='upward CO2 flux - NPP (acc.)',                            &
         units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flx_resp', co2_flux_soilresp_acc,   &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=22, bits=24,                                                   &
         longname='upward CO2 flux - soil respiration (acc.)',               &
         units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flx_herb', co2_flux_herbivory_acc,  &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=23, bits=24,                                                   &
         longname='upward CO2 flux - herbivory (acc.)',                      &
         units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flx_lcc', co2_emission_lcc_acc,     &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=24, bits=24,                                                   &
         longname='upward CO2 flux - land use change (acc.)',                &
         units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flx_harvest', co2_emission_harvest_acc, &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=25, bits=24,                                                   &
         longname='upward CO2 flux - harvest (acc.)',                        &
         units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flx_fire', co2_flux_dynveg_acc,     &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=26, bits=24,                                                   &
         longname='upward CO2 flux - fire (acc.)',                           &
         units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flx_ocean', co2_flux_ocean_acc,     &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=7, bits=24,                                                    &
         longname='upward ocean CO2 flux (acc.)',                            &
         units='kg m-2 s-1', laccu=.TRUE., lpost=lpost.AND.lcouple,          &
         lrerun=lrerun, contnorest=.TRUE.) 
    IF (lco2_emis) THEN
       CALL add_stream_element (xtco2, 'co2_emis', co2_emission_acc,         &
            (/ lnlon, lnglat /), (/ nlon, nglat /),                          &
            code=11, bits=24,                                                &
            longname='upward CO2 emission (acc.)',                           &
            units='kg m-2 s-1', laccu=.TRUE., lrerun=lrerun,                 &
            contnorest=.TRUE.) 
    END IF

    CALL add_stream_element (xtco2, 'co2_flx_anthro', co2_flux_anthro_acc,   &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=20, bits=24,                                                   &
         longname='upward anthropogenic CO2 flux (acc.)',                    &
         units='kg m-2 s-1', laccu=.TRUE., lpost=lpost.AND.lcouple,          &
         lrerun=lrerun, contnorest=.TRUE.) 

    ! Accumulated vertical integral of CO2 concentration
    CALL add_stream_element (xtco2, 'co2_burden', co2_burden_acc,            &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=8,                                                             &
         longname='CO2 content (acc.)',                                      &
         units='kg m-2', laccu=.TRUE., lpost=.TRUE., lrerun=lrerun,          &
         contnorest=.TRUE.) 
    CALL add_stream_element (xtco2, 'co2_flux_corr', co2_flux_corr_acc,      &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=10,                                                            &
         longname='CO2 flux correction (acc.)',                              &
         units='kg m-2 s-1', laccu=.TRUE., lpost=lpost .AND. lco2,           &
         lrerun=lrerun .AND. lco2, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_burden_corr_acc1',                  &
         co2_burden_corr_acc1,                                               &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         units='kg m-2', laccu=.FALSE., lpost=.FALSE.,                       &
         ! Is accumulated manually for daily correction time step
         lrerun=lrerun .AND. lco2, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flux_corr_acc1', co2_flux_corr_acc1,&
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         units='kg m-2', laccu=.FALSE., lpost=.FALSE.,                       &
         lrerun=lrerun .AND. lco2, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_flux_corr_acc1_old',                &
         co2_flux_corr_acc1_old,                                             &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         units='kg m-2', laccu=.FALSE., lpost=.FALSE.,                       &
         lrerun=lrerun .AND. lco2, contnorest=.TRUE.)
    CALL add_stream_element (xtco2, 'co2_burden_corr_acc2',                  &
         co2_burden_corr_acc2,                                               &
         (/ lnlon, lnglat /), (/ nlon, nglat /),                             &
         code=9,                                                             &
         longname='CO2 correction in tendencies (acc.)',                     &
         units='kg m-2', laccu=.TRUE., lpost=lpost .AND. lco2 .AND. lco2_2perc, &
         lrerun=.FALSE., contnorest=.TRUE.)

    IF (lcouple) co2atmos = 0._wp
    IF (lcouple_co2) THEN
       co2trans = 0._wp
       co2ocean = 0._wp
    END IF
    co2flux_cpl = 0._wp
    co2_flux = 0._wp
    co2_flux_land = 0._wp
    co2_flux_npp = 0._wp
    co2_flux_soilresp = 0._wp
    co2_flux_herbivory = 0._wp
    co2_emission_lcc = 0._wp
    co2_emission_harvest = 0._wp
    co2_flux_dynveg = 0._wp
    co2_flux_ocean = 0._wp
    co2_flux_acc = 0._wp
    co2_flux_land_acc = 0._wp
    co2_flux_npp_acc = 0._wp
    co2_flux_soilresp_acc = 0._wp
    co2_flux_herbivory_acc = 0._wp
    co2_emission_lcc_acc = 0._wp
    co2_emission_harvest_acc = 0._wp
    co2_flux_dynveg_acc = 0._wp
    co2_flux_ocean_acc = 0._wp
    co2_flux_anthro_acc = 0._wp
    co2_flux_total = 0._wp
    co2_burden = -1._wp
    co2_burden_old = 0._wp
    co2_burden_acc = 0._wp
    co2_flux_corr_acc = 0._wp
    co2_burden_corr_acc1 = 0._wp
    co2_flux_corr_acc1 = 0._wp
    co2_flux_corr_acc1_old = 0._wp
    co2_burden_corr_acc2 = 0._wp
    IF (lco2_emis) THEN
       co2_emission = 0._wp
       co2_emission_acc = 0._wp
    END IF

  END SUBROUTINE init_co2_memory

  SUBROUTINE cleanup_co2_memory

    CALL delete_stream (xtco2)

  END SUBROUTINE cleanup_co2_memory

  !---------------------------------------------------------------------------
  !>
  !! @brief Pointer co2 and co2m1
  !! 
  !! @remarks Get pointers to the 3d-field of CO2 at time t and t-dt.
  !! Using the pointers instead of "get_tracer" calls speeds up the program
  !
  SUBROUTINE reference_co2
    
    INTEGER :: ierr
    CHARACTER (LEN=32)      :: trname

    IF (.NOT. lco2) RETURN

    IF (lco2_clim) THEN
      trname = 'CO2_clim'
    ELSE
      trname = 'CO2'
    END IF

    CALL get_tracer (TRIM(trname), idx = ico2idx, pxt = co2, pxtm1 = co2m1, &
         ierr = ierr) 

  END SUBROUTINE reference_co2

  !---------------------------------------------------------------------------
  !>
  !! @brief Read CO2 emissions from files
  !! 
  !! @remarks Get pointers to the 3d-field of CO2 at time t and t-dt.
  !! Using the pointers instead of "get_tracer" calls speeds up the program
  !
  SUBROUTINE read_co2_emission

    ! Read annual carbon emissions from netCDF file and convert to CO2 

    USE mo_time_conversion, ONLY: time_native, &
                                  TC_get, TC_convert
    USE mo_time_control, ONLY: previous_date, current_date, lstart, lresume
    USE mo_transpose, ONLY: scatter_gp
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    USE mo_control, ONLY: nlon, ngl
    USE mo_netcdf

    REAL(wp), PARAMETER :: amc = 12.0107_wp   ! Molecular weight of C [g mol-1]

    INTEGER :: nvarid, ndimid
    INTEGER :: zcount(3), zstart(3)
    REAL(wp), POINTER :: ztmp1(:,:)

    TYPE(time_native) :: emis_date
    INTEGER :: iyear
    INTEGER :: yr, mo, dy, hr, mn, se, yrm

    IF (.NOT.lco2_emis) RETURN
    IF (emis_no_years < 0) THEN

       IF (p_parallel_io) THEN
          CALL nf_check(nf_open('carbon_emissions.nc', nf_nowrite, emis_id))
          CALL nf_check(nf_inq_dimid(emis_id, 'time', ndimid))
          CALL nf_check(nf_inq_dimlen(emis_id, ndimid, emis_no_years))
       END IF
       CALL p_bcast(emis_no_years, p_io)
       ALLOCATE(emis_years(emis_no_years))
       IF (p_parallel_io) THEN
          CALL nf_check(nf_inq_varid(emis_id, 'time', nvarid))
          CALL nf_check(nf_get_var_double(emis_id, nvarid, emis_years))
       ENDIF

       CALL p_bcast(emis_years, p_io)
       emis_base_year = emis_years(1)
       
    END IF

    CALL TC_convert(previous_date, emis_date) 
    CALL TC_get (emis_date, yrm, mo, dy, hr, mn, se)
    CALL TC_convert(current_date, emis_date) 
    CALL TC_get (emis_date, yr, mo, dy, hr, mn, se)

!!$    IF (yrm == yr) RETURN      ! same year, nothing to do
    IF (yrm == yr .AND. .NOT. (lstart .OR. lresume)) RETURN      ! same year, nothing to do
    ! We're either at the beginning of a new year or this is an initial or resumed run.
    ! The emission fields for yr-1,yr,yr+1 have therefore to be read from the netCDF file.

    iyear = yr-INT(emis_base_year)+1   ! set right index to access in emission fields

    IF (yr < emis_base_year .OR. iyear > emis_no_years) THEN
       CALL finish('read_co2_emission', 'Year not present in emission file')
    END IF

    co2_emission = 0._wp

    IF (p_parallel_io) THEN
       ALLOCATE(ztmp1(nlon,ngl))
       zstart(:) = (/ 1,1,iyear /)
       zcount(:) = (/ nlon,ngl,1 /)
       CALL nf_check(nf_inq_varid(emis_id, 'carbon_emission', nvarid))
       CALL nf_check(nf_get_vara_double(emis_id, nvarid, zstart, zcount, ztmp1))
    END IF
    CALL scatter_gp(ztmp1, co2_emission, gl_dc)

    IF (p_parallel_io) THEN
       DEALLOCATE(ztmp1)
       CALL print_value (' read_co2_emission for year: ', yr)
    END IF
    
    ! Convert carbon emission in g m2-1 s-1 to CO2 emission in kg m2-1 s-1
    co2_emission   = co2_emission   * amco2 / (amc * 1000._wp)

  END SUBROUTINE read_co2_emission

  !---------------------------------------------------------------------------
  !>
  !! @brief Calculate burden and burden correction if necessary
  !! 
  !! @remarks This routine calculates the actual CO2 burden and flux 
  !!   correction in order to close the CO2 budget in the simulation
  !
  SUBROUTINE co2_mbalance(kproma,          kbdim,             klev,           &
                          klevp1,          krow,            paphm1)

    USE mo_physical_constants,      ONLY: rgrav
    USE mo_time_control,   ONLY: lstart, delta_time
    
    INTEGER, INTENT(in)       :: kproma, kbdim, klev, klevp1, krow
    REAL(wp), INTENT(in)      :: paphm1(kbdim,klevp1)

    ! Local variables

    INTEGER                   :: jk
    REAL(wp)                  :: zflux_corr_acc(kproma)

    co2_burden_old(1:kproma,krow) = co2_burden(1:kproma,krow)
    co2_burden(1:kproma,krow) = 0._wp
    DO jk=1,klev
       co2_burden(1:kproma,krow) = co2_burden(1:kproma,krow) + &
           co2m1(1:kproma,jk,krow)*(paphm1(1:kproma,jk+1)-paphm1(1:kproma,jk))*rgrav
    END DO

    IF (.NOT. lco2) RETURN

    IF (ico2idx > 0 .AND. lco2_flxcor) THEN
       IF (ALL(co2_burden_old(1:kproma,krow) < 0._wp)) THEN
          IF (krow == 1) THEN
             CALL message('co2_mbalance','co2_burden < 0')
             CALL message('co2_mbalance', &
                          '   assuming that it was not present in restart file,')
             CALL message('co2_mbalance', &
                          '   setting correction in CO2 burden to zero')
             CALL message('co2_mbalance', &
                          '   (this should only happen on the first time step')
             CALL message('co2_mbalance', &
                          '    of the first restart after introduction of')
             CALL message('co2_mbalance', &
                          '    flux correction, or for lstart=true)')
          END IF
          IF (ANY(co2_flux_total(1:kproma,krow) /= 0._wp)) THEN
             CALL finish('co2_mbalance', &
                         'co2_flux_total should also be zero in this case!')
          END IF
          IF (ANY(co2_flux_corr_acc1_old /= 0._wp)) THEN
             CALL finish('co2_mbalance', &
                         'co2_flux_corr_acc1_old should also be zero in this case!')
          END IF
          co2_burden_old(1:kproma,krow) = co2_burden(1:kproma,krow)
       END IF

       IF (l_co2flxcorr .OR. lstart) THEN
          co2_burden_corr_acc1(1:kproma,krow) = 0._wp! set to zero after coupling
          IF (krow == 1) CALL message('co2_mbalance', &
                         'Resetting acc. CO2 burden correction')
       ELSE
          co2_flux_corr_acc1_old(1:kproma,krow)=co2_flux_corr_acc1(1:kproma,krow)
       END IF

       ! Compute flux correction field from 
       !    (old burden + fluxes + global mean flux corr.) - (new burden).
       ! This is the instant. flux correction multiplied by zdtime, i.e. it is an 
       ! accumulated value
       ! Uses old co2_flux_total which is only updated in vdiff
       zflux_corr_acc = 0.0_wp
       zflux_corr_acc(1:kproma) = &
            (co2_burden_old(1:kproma,krow) + (co2_flux_total(1:kproma,krow)+    &
             co2_flux_corr_acc1_old(1:kproma,krow))*delta_time)                 &
           - co2_burden(1:kproma,krow)
       ! Accumulate flux correction field for output interval
       ! (co2_flux_corr_acc has laccu=.TRUE., i.e. is divided by interval 
       !  length before output)
       co2_flux_corr_acc(1:kproma,krow) = co2_flux_corr_acc(1:kproma,krow)        &
                                           + zflux_corr_acc(1:kproma)

       ! Accumulate CO2 burden correction for daily correction update
       ! (co2_burden_corr_acc1 has laccu=.FALSE. and is reset above to zero at 
       !  correction time step, but not divided by interval length)

       co2_burden_corr_acc1(1:kproma,krow) = co2_burden_corr_acc1(1:kproma,krow)&
                                              + zflux_corr_acc(1:kproma)
    END IF
  END SUBROUTINE co2_mbalance


  !---------------------------------------------------------------------------
  !>
  !! @brief Calculate net CO2 flux at surface
  !! 
  !! @remarks This routine adds all fluxes and emissions at surface to get net CO2 flux as
  !! needed as lower boundary condition in vdiff
  !
  SUBROUTINE co2_exchange(kproma, kbdim, ktrac, pfrl, pxtems, krow)

    INTEGER,  INTENT(in)    :: kproma, kbdim, ktrac, krow
    REAL(wp), INTENT(in)    :: pfrl(kbdim)                   ! Land fraction
    REAL(wp), INTENT(inout) :: pxtems(kbdim,ktrac)

    ! Copy natural CO2 fluxes (from ocean and from land) into emissions field
    pxtems(1:kproma,ico2idx) = co2_flux(1:kproma,krow)

    ! Copy CO2 emissions from landcover change and harvest into emissions field
    pxtems(1:kproma,ico2idx) = pxtems(1:kproma,ico2idx) + &
         pfrl(1:kproma) * (co2_emission_lcc(1:kproma,krow) + co2_emission_harvest(1:kproma,krow))

    ! Copy CO2 emissions from fossil fuel combustion into emissions field
    IF (lco2_emis) THEN
       pxtems(1:kproma, ico2idx) = pxtems(1:kproma,ico2idx) + &
                                   co2_emission(1:kproma,krow)
    END IF

    co2_flux_total(1:kproma,krow) = pxtems(1:kproma,ico2idx)

    ! Add CO2 flux correction
    pxtems(1:kproma,ico2idx) =                    &
         pxtems(1:kproma,ico2idx) + co2_flux_corr_acc1(1:kproma,krow)
    
  END SUBROUTINE co2_exchange

  !---------------------------------------------------------------------------
  !>
  !! @brief ! Calculation of the CO2 flux between atmosphere and 
  !!          ocean [kg m-2 s-1]
  !! 
  !! @remarks This routine is active only if ocean coupling is enabled
  !
  SUBROUTINE co2_flux_atmosphere_ocean(jrow, kproma, klev)

    USE mo_memory_g3b,   ONLY: wind10w

    INTEGER, INTENT(in) :: jrow, kproma, klev

    ! The range (Max-Min) of co2trans is smaller than 1.e-14. Array with a range smaller 1.e-12 
    !  are set to zero in grib output. To avoid this problem the co2trans was multiplied by 1.e6 
    !  (mo_couple).
    
    co2_flux_ocean(1:kproma,jrow) = amco2 * co2trans(1:kproma,jrow) * wind10w(1:kproma,jrow)**2 * &
                                    (co2ocean(1:kproma,jrow) - (co2m1(1:kproma,klev,jrow) * &
                                    1.0E+06_wp * amd / amco2)) * 1.e-6

  END SUBROUTINE co2_flux_atmosphere_ocean

  !---------------------------------------------------------------------------
  !>
  !! @brief ! Check CO2 tendencies.
  !! 
  !! @remarks Check CO2 mixing ratios for non-negativity and limit the 
  !!   tendencies if lco2_2perc is set
  !
  SUBROUTINE co2_te_check(kproma,           kbdim,           klev,           &
                          klevp1,           ktrac,           krow,           &
                          paphm1,           loland,          pxtm1,          &
                          pxtte)

    USE mo_time_control,         ONLY: time_step_len, delta_time
    USE mo_geoloc,               ONLY: philat_2d, philon_2d
    USE mo_physical_constants,            ONLY: rgrav

    INTEGER, INTENT(in)     :: kproma, kbdim, klev, klevp1, ktrac, krow
    REAL(wp), INTENT(in)    :: paphm1(kbdim,klevp1), &  ! p at half levels t-dt
                               pxtm1(kbdim,klev,ktrac)   ! tracer mixing ratios at t-dt
    REAL(wp), INTENT(inout) :: pxtte(kbdim,klev,ktrac)   ! tracer tendencies
    LOGICAL, INTENT(in)     :: loland(kproma)         ! logical land mask

    INTEGER                 :: jk,jl    
    REAL(wp)                :: zco2eps, zco2tend, dt

    dt=time_step_len
    ! Check for negative CO2 concentration
    IF (MINVAL(pxtm1(1:kproma,1:klev,ico2idx)) <= 0._wp) THEN
       CALL finish('co2_te_check','Negative CO2 concentration. Stop.')
    END IF
    ! Limit change in CO2 concentration to +/- 2% of absolute 
    !  concentration over ocean
    IF (lco2_2perc) THEN
      DO jk=1,klev
        DO jl=1,kproma
          IF (.NOT. loland(jl)) THEN
            zco2eps = 0.02_wp * pxtm1(jl,jk,ico2idx) / dt
            zco2tend = pxtte(jl,jk,ico2idx)
            pxtte(jl,jk,ico2idx) =                                            &
                     MIN(MAX(pxtte(jl,jk,ico2idx), -zco2eps), zco2eps)
            IF (pxtte(jl,jk,ico2idx) /= zco2tend) THEN
              CALL print_value('co2_te_check: Limit change in CO2 '// &
                               'mixing ratio, lat',philat_2d(jl,krow))
              CALL print_value('co2_te_check: Limit change in CO2 '// &
                               'mixing ratio, lon',philon_2d(jl,krow))
              CALL print_value('co2_te_check: Limit change in CO2 '// &
                               'mixing ratio, lev',klev)
            END IF
            co2_burden_corr_acc2(jl,krow) = co2_burden_corr_acc2(jl,krow)     &
                                 + (pxtte(jl,jk,ico2idx)-zco2tend)*dt &
                                 * delta_time*(paphm1(jl,jk+1)-paphm1(jl,jk))*rgrav
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE co2_te_check


  !---------------------------------------------------------------------------
  !>
  !! @brief ! Diagnose CO2 fluxes.
  !! 
  !! @remarks Called in physc_subm_4 
  !
  SUBROUTINE diag_co2(kproma,                                krow,                &
                      pfrl,             pfrw,                pfri)

    USE mo_time_control, ONLY: delta_time, lstart

    ! Scalar arguments
    INTEGER, INTENT(in) :: kproma, krow

    ! Array arguments
    REAL(wp), INTENT(in) :: pfrl(kproma)          ! Land fraction
    REAL(wp), INTENT(in) :: pfrw(kproma)          ! Ocean fraction
    REAL(wp), INTENT(in) :: pfri(kproma)          ! Sea ice fraction

    co2_flux_acc(1:kproma,krow)             = co2_flux_acc(1:kproma,krow)           &
                + co2_flux(1:kproma,krow) * delta_time
    co2_flux_land_acc(1:kproma,krow)        = co2_flux_land_acc(1:kproma,krow)      &
                + pfrl(1:kproma) * co2_flux_land(1:kproma,krow) * delta_time
    co2_flux_npp_acc(1:kproma,krow)         = co2_flux_npp_acc(1:kproma,krow)       &
                + pfrl(1:kproma) * co2_flux_npp(1:kproma,krow) * delta_time
    co2_flux_soilresp_acc(1:kproma,krow)    = co2_flux_soilresp_acc(1:kproma,krow)  &
                + pfrl(1:kproma) * co2_flux_soilresp(1:kproma,krow) * delta_time
    co2_flux_herbivory_acc(1:kproma,krow)   = co2_flux_herbivory_acc(1:kproma,krow) &
                + pfrl(1:kproma) * co2_flux_herbivory(1:kproma,krow) * delta_time
    co2_flux_dynveg_acc(1:kproma,krow)      = co2_flux_dynveg_acc(1:kproma,krow)    &
                + pfrl(1:kproma) * co2_flux_dynveg(1:kproma,krow) * delta_time
    co2_emission_lcc_acc(1:kproma,krow)     = co2_emission_lcc_acc(1:kproma,krow)   &
                + pfrl(1:kproma) * co2_emission_lcc(1:kproma,krow) * delta_time
    co2_emission_harvest_acc(1:kproma,krow) = co2_emission_harvest_acc(1:kproma,krow) &
                + pfrl(1:kproma) * co2_emission_harvest(1:kproma,krow) * delta_time
    co2_flux_ocean_acc(1:kproma,krow)       = co2_flux_ocean_acc(1:kproma,krow)     &
                + (pfrw(1:kproma)+pfri(1:kproma))                                   &
                * co2_flux_ocean(1:kproma,krow) * delta_time

    ! dCO2/dt = E_fossil + E_lcc + F_ocean + F_land + CO2_flux_correction
    ! where E are anthropogenic emissions and F are natural fluxes
    ! E_lcc includes the emissions due to land cover change and, for LU transitions, harvest 
    ! (excluding natural fire)
    IF (.NOT. lstart) THEN
       co2_flux_anthro_acc(1:kproma,krow) = co2_flux_anthro_acc(1:kproma,krow)  &
                   + (co2_burden(1:kproma,krow)-co2_burden_old(1:kproma,krow))  &
                   - co2_flux_corr_acc1(1:kproma,krow) * delta_time             & 
                                                                              ! CO2 flux correction
                   - co2_flux(1:kproma,krow) * delta_time                     ! F_ocean + F_land
    END IF
    IF (lco2_emis) THEN
       co2_emission_acc(1:kproma,krow) = co2_emission_acc(1:kproma,krow)   &
                                    + co2_emission(1:kproma,krow) * delta_time
    END IF
    
    co2_burden_acc(1:kproma,krow) = co2_burden_acc(1:kproma,krow)          &
                                    + co2_burden(1:kproma,krow) * delta_time
  END SUBROUTINE diag_co2

  SUBROUTINE co2_flux_correction

    USE mo_control,       ONLY: ngl, nlon
    USE mo_decomposition, ONLY: gl_dc => global_decomposition
    USE mo_transpose,     ONLY: gather_gp
    USE mo_gaussgrid,     ONLY: gl_budw
    USE mo_mpi,           ONLY: p_pe, p_io, p_bcast
    USE mo_time_control,  ONLY: get_date_components, current_date, next_date

    IMPLICIT NONE

    ! local arrays
    REAL(wp), POINTER :: zco2_burden_corr_glob(:,:)
    REAL(wp)          :: ztime = 0.0_wp
    INTEGER           :: day_of_month, day_of_month_at_prev_ts
    REAL(wp)          :: ztemp(ngl)

    REAL(wp)     :: zco2_burden_corr_mean, zco2_flux_corr_mean !, zco2_flux_corr

    INTEGER :: jlat

    IF (.NOT. lco2) RETURN

    ! Set values to zero if CO2 flux correction not requested (e.g. in case of CO2 nudging)
    IF (.NOT. lco2_flxcor) THEN
       co2_flux_corr_acc1(:,:) = 0._wp
       co2_flux_corr_acc1_old(:,:) = 0._wp
       RETURN
    END IF

    CALL get_date_components(next_date,DAY=day_of_month)
    CALL get_date_components(current_date, DAY=day_of_month_at_prev_ts)
    l_co2flxcorr = day_of_month /= day_of_month_at_prev_ts

    ztime = 86400._wp

    IF (.NOT. l_co2flxcorr) THEN
       RETURN
    END IF

    CALL message('co2_flux_correction','Computing global mean CO2 flux correction')

    IF (p_pe == p_io) THEN
       ALLOCATE(zco2_burden_corr_glob(nlon,ngl))
    END IF

!!$    co2_flux_corr_old = co2_flux_corr
    co2_flux_corr_acc1_old(:,:) = co2_flux_corr_acc1(:,:)

    CALL gather_gp (zco2_burden_corr_glob, co2_burden_corr_acc1, gl_dc)

    ! Global mean
    IF (p_pe == p_io) THEN
       DO jlat=1,ngl
          ztemp(jlat) = SUM(zco2_burden_corr_glob(1:nlon,jlat)) * gl_budw(jlat)
       END DO
       zco2_burden_corr_mean = SUM(ztemp)

       ! Apply global CO2 burden correction to CO2 flux
       zco2_flux_corr_mean = zco2_burden_corr_mean / ztime
       WRITE(message_text,*) 'CO2 flux correction (global): ',zco2_flux_corr_mean, ' kg m-2 s-1'
       CALL message('', message_text)

    END IF

    CALL p_bcast(zco2_flux_corr_mean, p_io)
    co2_flux_corr_acc1(:,:) = zco2_flux_corr_mean

    IF (p_pe == p_io) THEN
       DEALLOCATE(zco2_burden_corr_glob)
    END IF

    RETURN

  END SUBROUTINE co2_flux_correction


END MODULE mo_co2
