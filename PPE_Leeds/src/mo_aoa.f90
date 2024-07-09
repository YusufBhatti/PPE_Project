!>
!! @brief Module to provide tracer concentrations for age-of-air calculation
!!
!! @remarks
!!   This module was implemented for use as a submodel for echam6. It provides
!!   concentrations of passive tracers, which can be used to calculate mean
!!   age-of-air (mean_age) as well as age spectra (spec_winter, spec_summer).
!!   
!!
!! @author F.Bunzel, MPI-M, Hamburg (2010-12-23)
!!
!! $ID: n/a$
!!
!! @par Origin Module was derived from a transport module compiled by J.S.Rast
!!
!! @par Copyright
!!   2002-2010 by the Deutsche Wetterdienst (DWD) and the Max-Planck-Institut
!!   for Meteorology (MPI-M).  This software is provided for non-commerical 
!!   use only.  See the LICENSE and the WARRANTY conditions
!! 
!! @par License
!!   The use of ICON is hereby granted free of charge for an unlimited time,
!!   provided:
!!   <ol>
!!    <li> Its use is limited to own non-commercial and non-violent purposes;
!!    <li> The code is not re-distributed without the consent of DWD and MPI-M;
!!    <li> This header appears in all copies of the code;
!!    <li> You accept the warranty conditions (see WARRANTY).
!!   </ol>
!!   Commericial use of the code is allowed subject to a separate licensing 
!!   agreement with the DWD and MPI-M
!!
!! @par Warranty
!!   This code is distributed in the home that it will be useful, but WITHOUT
!!   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
!!   FITNESS FOR A PARTICULAR PURPOSE.  
!
MODULE mo_aoa

  USE mo_kind,                 ONLY: dp
  USE mo_time_conversion,      ONLY: time_days

  IMPLICIT NONE

  PRIVATE
  PUBLIC                           :: init_aoa, bcond_aoa, get_pointer2trac, &
                                      tracer_reset, tf_reset

  CHARACTER(LEN=3), PARAMETER      :: c_mod_name='aoa'
  INTEGER                          :: dt_start_emi_winter_1(6)=(/1978,1,2,0,0,0/)
  INTEGER                          :: dt_start_emi_summer_1(6)=(/1978,7,2,0,0,0/)
  INTEGER                          :: dt_start_emi_winter_2(6)=(/1990,1,2,0,0,0/)
  INTEGER                          :: dt_start_emi_summer_2(6)=(/1990,7,2,0,0,0/)
  INTEGER                          :: dt_start_emi_winter_3(6)=(/2002,1,2,0,0,0/)
  INTEGER                          :: dt_start_emi_summer_3(6)=(/2002,7,2,0,0,0/)
  INTEGER                          :: dt_start_emi_winter_4(6)=(/2014,1,2,0,0,0/)
  INTEGER                          :: dt_start_emi_summer_4(6)=(/2014,7,2,0,0,0/)
  INTEGER                          :: idx_lev_emission, idx_xtrans, idx_xtrans_winter, idx_xtrans_summer
  REAL(dp)                         :: conc_increase=1.0e-6_dp ! emission mass flux in kg/s
  REAL(dp)                         :: emission_plev=1013.25_dp  ! pressure level in hPa
  REAL(dp)                         :: emission_lat_s=-5._dp, emission_lat_n=5._dp
  REAL(dp), POINTER                :: emission_mask(:,:)
  REAL(dp), POINTER                :: pxt_winter(:,:,:), pxtm1_winter(:,:,:)
  REAL(dp), POINTER                :: pxt_summer(:,:,:), pxtm1_summer(:,:,:)
  TYPE(time_days), SAVE            :: start_emi_winter_1,start_emi_summer_1
  TYPE(time_days), SAVE            :: start_emi_winter_2,start_emi_summer_2
  TYPE(time_days), SAVE            :: start_emi_winter_3,start_emi_summer_3
  TYPE(time_days), SAVE            :: start_emi_winter_4,start_emi_summer_4

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief Initialization of age-of-air submodel
  !
  SUBROUTINE init_aoa
    USE mo_exception,             ONLY: finish, message
    USE mo_mpi,                   ONLY: p_parallel_io, p_parallel, p_bcast, p_io
    USE mo_namelist,              ONLY: open_nml, position_nml, close_nml, POSITIONED, &
                                        MISSING, LENGTH_ERROR, READ_ERROR
!    USE mo_io_units,              ONLY: nnml
    USE mo_tracer,                ONLY: new_tracer, get_tracer
    USE mo_tracdef,               ONLY: ON, OFF, GAS, RESTART, CONSTANT
    USE mo_time_conversion,       ONLY: TC_set, TC_convert, time_native
    USE mo_memory_base,           ONLY: new_stream, add_stream_element, t_stream
    USE mo_gaussgrid,             ONLY: philat
    USE mo_transpose,             ONLY: scatter_gp
    USE mo_decomposition,         ONLY: dc => local_decomposition, global_decomposition
    USE mo_control,               ONLY: nlev, nvclev, vct ! for hyai, hybi coefficients
    USE mo_util_string,           ONLY: int2string

    INTEGER                          :: jlev, jlon, jlat
    INTEGER                          :: inml,iunit,ierr
    REAL(dp), POINTER                :: zmask(:,:)                                        
    REAL(dp), DIMENSION(nvclev)      :: hyai, hybi
    REAL(dp)                         :: ztp
    TYPE(t_stream), POINTER          :: trans
    TYPE(time_native)                :: date_nat

  include 'aoactl.inc'

  IF (p_parallel_io) THEN
     inml = open_nml ('namelist.echam')
     iunit = position_nml ('AOACTL', inml, status=ierr)
!     CALL position_nml ('AOACTL', status=ierr)
     SELECT CASE (ierr)
     CASE (POSITIONED)
        READ(iunit, aoactl)
     CASE (MISSING)
        CALL finish ('init_aoa','namelist aoactl not found in namelist.echam')
     CASE (LENGTH_ERROR)
        CALL finish ('init_aoa','namelist aoactl has wrong length')
     CASE (READ_ERROR)
        CALL finish ('init_aoa','cannot read namelist aoactl')
     END SELECT
     CALL close_nml(inml)
  END IF
  IF (p_parallel) THEN
     CALL p_bcast (dt_start_emi_winter_1, p_io)
     CALL p_bcast (dt_start_emi_summer_1, p_io)
     CALL p_bcast (dt_start_emi_winter_2, p_io)
     CALL p_bcast (dt_start_emi_summer_2, p_io)
     CALL p_bcast (dt_start_emi_winter_3, p_io)
     CALL p_bcast (dt_start_emi_summer_3, p_io)
     CALL p_bcast (dt_start_emi_winter_4, p_io)
     CALL p_bcast (dt_start_emi_summer_4, p_io)
     CALL p_bcast (conc_increase, p_io)
     CALL p_bcast (emission_plev, p_io)
     CALL p_bcast (emission_lat_s, p_io)
     CALL p_bcast (emission_lat_n, p_io)
  END IF
  CALL get_tracer('mean_age',idx=idx_xtrans,ierr=ierr)
  IF (ierr == 1) THEN
     CALL new_tracer('mean_age',              c_mod_name,                &
                     idx=idx_xtrans,          units='kg/kg',             &
                     ndrydep=OFF,             nwetdep=OFF,               &
                     nphase=GAS,              moleweight=1._dp,          &
                     longname='tracer mass mixing ratio in kg/kg',       &
                     table=131,               code=182,                  &
                     nwrite=ON,               ninit=RESTART+CONSTANT,    &
                     vini=0.0_dp,             nrerun=ON)
  END IF

  CALL get_tracer('spec_winter',idx=idx_xtrans_winter,ierr=ierr)
  IF (ierr == 1) THEN
     CALL new_tracer('spec_winter',           c_mod_name,                &
                     idx=idx_xtrans_winter,    units='kg/kg',            &
                     ndrydep=OFF,             nwetdep=OFF,               &
                     nphase=GAS,              moleweight=1._dp,          &
                     longname='tracer mass mixing ratio in kg/kg',       &
                     table=131,               code=183,                  &
                     nwrite=ON,               ninit=RESTART+CONSTANT,    &
                     vini=0.0_dp,             nrerun=ON)
  END IF

  CALL get_tracer('spec_summer',idx=idx_xtrans_summer,ierr=ierr)
  IF (ierr == 1) THEN
     CALL new_tracer('spec_summer',           c_mod_name,                &
                     idx=idx_xtrans_summer,    units='kg/kg',            &
                     ndrydep=OFF,             nwetdep=OFF,               &
                     nphase=GAS,              moleweight=1._dp,          &
                     longname='tracer mass mixing ratio in kg/kg',       &
                     table=131,               code=184,                  &
                     nwrite=ON,               ninit=RESTART+CONSTANT,    &
                     vini=0.0_dp,             nrerun=ON)
  END IF

  ! Convert dt_start_1 in internal time format
  CALL TC_set(&
       dt_start_emi_winter_1(1), dt_start_emi_winter_1(2), dt_start_emi_winter_1(3), &
       dt_start_emi_winter_1(4), dt_start_emi_winter_1(5), dt_start_emi_winter_1(6), date_nat)
  CALL TC_convert(date_nat, start_emi_winter_1)

  ! Convert dt_start_ref_1 in internal time format
  CALL TC_set(&
       dt_start_emi_summer_1(1), dt_start_emi_summer_1(2), dt_start_emi_summer_1(3), &
       dt_start_emi_summer_1(4), dt_start_emi_summer_1(5), dt_start_emi_summer_1(6), date_nat)
  CALL TC_convert(date_nat, start_emi_summer_1)

  ! Convert dt_start_2 in internal time format
  CALL TC_set(&
       dt_start_emi_winter_2(1), dt_start_emi_winter_2(2), dt_start_emi_winter_2(3), &
       dt_start_emi_winter_2(4), dt_start_emi_winter_2(5), dt_start_emi_winter_2(6), date_nat)
  CALL TC_convert(date_nat, start_emi_winter_2)

  ! Convert dt_start_ref_2 in internal time format
  CALL TC_set(&
       dt_start_emi_summer_2(1), dt_start_emi_summer_2(2), dt_start_emi_summer_2(3), &
       dt_start_emi_summer_2(4), dt_start_emi_summer_2(5), dt_start_emi_summer_2(6), date_nat)
  CALL TC_convert(date_nat, start_emi_summer_2)

  ! Convert dt_start_3 in internal time format
  CALL TC_set(&
       dt_start_emi_winter_3(1), dt_start_emi_winter_3(2), dt_start_emi_winter_3(3), &
       dt_start_emi_winter_3(4), dt_start_emi_winter_3(5), dt_start_emi_winter_3(6), date_nat)
  CALL TC_convert(date_nat, start_emi_winter_3)

  ! Convert dt_start_ref_3 in internal time format
  CALL TC_set(&
       dt_start_emi_summer_3(1), dt_start_emi_summer_3(2), dt_start_emi_summer_3(3), &
       dt_start_emi_summer_3(4), dt_start_emi_summer_3(5), dt_start_emi_summer_3(6), date_nat)
  CALL TC_convert(date_nat, start_emi_summer_3)

  ! Convert dt_start_4 in internal time format
  CALL TC_set(&
       dt_start_emi_winter_4(1), dt_start_emi_winter_4(2), dt_start_emi_winter_4(3), &
       dt_start_emi_winter_4(4), dt_start_emi_winter_4(5), dt_start_emi_winter_4(6), date_nat)
  CALL TC_convert(date_nat, start_emi_winter_4)

  ! Convert dt_start_ref_4 in internal time format
  CALL TC_set(&
       dt_start_emi_summer_4(1), dt_start_emi_summer_4(2), dt_start_emi_summer_4(3), &
       dt_start_emi_summer_4(4), dt_start_emi_summer_4(5), dt_start_emi_summer_4(6), date_nat)
  CALL TC_convert(date_nat, start_emi_summer_4)

  ! Calculate indices of location of emissions
  ! Create a 2d mask field for the location, levels are treated separately
  CALL new_stream(trans,'trans',lpost=.true.,lrerun=.false.)
  CALL add_stream_element(trans,'emission_mask',emission_mask,lpost=.true.,lrerun=.false.)

  ALLOCATE (zmask(dc%nlon,dc%nlat))
  zmask=0._dp
  DO jlat=1,dc%nlat
     IF (philat(jlat) > emission_lat_s .AND. philat(jlat) < emission_lat_n) THEN
        DO jlon = 1,dc%nlon
           zmask(jlon,jlat)=1.0_dp
        END DO
     END IF
  END DO
  CALL scatter_gp(zmask,emission_mask,global_decomposition)
  IF (p_parallel_io) THEN
     DEALLOCATE(zmask)
  END IF

  ! hyam,hybm
   hyai = vct(1:nvclev)
   hybi = vct(nvclev+1:2*nvclev)
   idx_lev_emission=0
   DO jlev = 1, nlev
      ztp = hyai(jlev+1) + hybi(jlev+1) * 101325.
      IF ( ztp >= emission_plev*100._dp ) THEN
         idx_lev_emission = jlev
         EXIT
      END IF
   END DO

   IF (p_parallel_io) THEN
      CALL message ('init_aoa','index of level for emissions: '&
                     //int2string(idx_lev_emission))
   END IF

END SUBROUTINE init_aoa

  !--------------------------------------------------------------------------
  !>
  !! @brief Apply emissions
  !
SUBROUTINE bcond_aoa(kproma, kbdim, klev, krow, pxtte, pxtm1)
  USE mo_time_control,         ONLY: current_date, next_date, time_step_len
  USE mo_time_conversion,      ONLY: OPERATOR(>), OPERATOR(<), OPERATOR(==), TC_get
  USE mo_tracdef,              ONLY: ntrac

! !INPUT PARAMETERS:
  INTEGER, INTENT(IN)             :: kproma, kbdim, klev, krow
  REAL(dp), INTENT(IN)            :: pxtm1(kbdim,klev,ntrac)

! !INPUT/OUTPUT PARAMETERS:
  REAL(dp), INTENT(INOUT)         :: pxtte(kbdim,klev,ntrac)
  REAL(dp)                        :: conc

! !LOCAL VARIABLES
  LOGICAL                         :: linject,linject_winter,linject_summer
  INTEGER                         :: jk, current_day, current_second, start_day, start_second

  jk=idx_lev_emission
  linject=.false.
  linject_winter=.false.
  linject_summer=.false.
  IF (next_date > start_emi_winter_1) THEN
     linject=.true.
  END IF
  IF (current_date == start_emi_winter_1   .OR.                                         &
      ((current_date < start_emi_winter_1) .AND. (next_date > start_emi_winter_1)) .OR. &
      current_date == start_emi_winter_2   .OR.                                         &
      ((current_date < start_emi_winter_2) .AND. (next_date > start_emi_winter_2)) .OR. &
      current_date == start_emi_winter_3   .OR.                                         &
      ((current_date < start_emi_winter_3) .AND. (next_date > start_emi_winter_3)) .OR. &
      current_date == start_emi_winter_4   .OR.                                         &
      ((current_date < start_emi_winter_4) .AND. (next_date > start_emi_winter_4))      ) THEN
     linject_winter=.true.
  END IF
  IF (current_date == start_emi_summer_1   .OR.                                         &
      ((current_date < start_emi_summer_1) .AND. (next_date > start_emi_summer_1)) .OR. & 
      current_date == start_emi_summer_2   .OR.                                         &
      ((current_date < start_emi_summer_2) .AND. (next_date > start_emi_summer_2)) .OR. & 
      current_date == start_emi_summer_3   .OR.                                         &
      ((current_date < start_emi_summer_3) .AND. (next_date > start_emi_summer_3)) .OR. & 
      current_date == start_emi_summer_4   .OR.                                         &
      ((current_date < start_emi_summer_4) .AND. (next_date > start_emi_summer_4))      ) THEN
     linject_summer=.true.
  END IF

  IF (linject) THEN
     call TC_get(start_emi_winter_1, start_day, start_second)
     call TC_get(current_date, current_day, current_second)
     conc=conc_increase*((time_step_len/86400._dp)+(current_day-start_day) &
          +((current_second-start_second)/86400._dp))
     WHERE (emission_mask(1:kproma,krow)>0._dp)
       pxtte(1:kproma,jk,idx_xtrans)=(conc-pxtm1(1:kproma,jk,idx_xtrans))/time_step_len
     END WHERE
  END IF

  IF (linject_winter) THEN
     WHERE (emission_mask(1:kproma,krow)>0._dp)
       pxtte(1:kproma,jk,idx_xtrans_winter)=(1._dp-pxtm1(1:kproma,jk,idx_xtrans_winter))/time_step_len
     END WHERE
  ELSE
     WHERE (emission_mask(1:kproma,krow)>0._dp)
       pxtte(1:kproma,jk,idx_xtrans_winter)=(0._dp-pxtm1(1:kproma,jk,idx_xtrans_winter))/time_step_len
     END WHERE
  END IF

  IF (linject_summer) THEN
     WHERE (emission_mask(1:kproma,krow)>0._dp)
       pxtte(1:kproma,jk,idx_xtrans_summer)=(1._dp-pxtm1(1:kproma,jk,idx_xtrans_summer))/time_step_len
     END WHERE
  ELSE
     WHERE (emission_mask(1:kproma,krow)>0._dp)
       pxtte(1:kproma,jk,idx_xtrans_summer)=(0._dp-pxtm1(1:kproma,jk,idx_xtrans_summer))/time_step_len
     END WHERE
  END IF
END SUBROUTINE bcond_aoa

  !-----------------------------------------------------------------------------
  !>
  !! @brief get pointer to tracer concentrations
  !
SUBROUTINE get_pointer2trac
!
! !DESCRIPTION:
! get pointer to tracer concentrations
!
! !REVISION HISTORY:
! original source by F.Bunzel (2010-10-04)
!
! !USES:
  USE mo_exception,            ONLY: message
  USE mo_tracer,               ONLY: get_tracer

! !LOCAL VARIABLES:
  INTEGER                      :: ierr

  CALL get_tracer('spec_winter',idx=idx_xtrans_winter,pxt=pxt_winter,pxtm1=pxtm1_winter,ierr=ierr)
  IF (ierr == 0) CALL message ('get_pointer2trac','winter tracer concentrations submitted!')

  CALL get_tracer('spec_summer',idx=idx_xtrans_summer,pxt=pxt_summer,pxtm1=pxtm1_summer,ierr=ierr)
  IF (ierr == 0) CALL message ('get_pointer2trac','summer tracer concentrations submitted!')

END SUBROUTINE get_pointer2trac

  !-----------------------------------------------------------------------------
  !>
  !! @brief reset tracer concentrations, if specified (via namelist) date is reached
  !
SUBROUTINE tracer_reset(kproma,kbdim,klev,krow,pxtte)
!
! !DESCRIPTION:
! reset tracer concentrations, if specified (via namelist) date is reached
!
! !REVISION HISTORY:
! original source by F. Bunzel (2010-10-04)
!
! !USES:
  USE mo_exception,            ONLY: message
  USE mo_time_control,         ONLY: next_date, echam_time, get_time_step
  USE mo_time_manager,         ONLY: manager_state
  USE mo_time_conversion,      ONLY: OPERATOR(==), OPERATOR(<), OPERATOR(>)
  USE mo_tracdef,              ONLY: ntrac

! !INPUT PARAMETERS:
  INTEGER, INTENT(IN)             :: kproma, kbdim, klev, krow

! !INPUT/OUTPUT PARAMETERS:
  REAL(dp), INTENT(INOUT)         :: pxtte(kbdim,klev,ntrac)

  ! LOCAL VARIABLES:
  TYPE(TIME_DAYS)              :: nnext_date
  INTEGER                      :: istep
  
  istep = get_time_step()
  CALL manager_state(echam_time,nnext_date,istep+2)
  IF (next_date == start_emi_winter_1   .OR.                                          &
      ((next_date < start_emi_winter_1) .AND. (nnext_date > start_emi_winter_1)) .OR. &
      next_date == start_emi_winter_2   .OR.                                          &
      ((next_date < start_emi_winter_2) .AND. (nnext_date > start_emi_winter_2)) .OR. &
      next_date == start_emi_winter_3   .OR.                                          &
      ((next_date < start_emi_winter_3) .AND. (nnext_date > start_emi_winter_3)) .OR. &
      next_date == start_emi_winter_4   .OR.                                          &
      ((next_date < start_emi_winter_4) .AND. (nnext_date > start_emi_winter_4))      )THEN
     pxt_winter(1:kproma,:,krow)=0._dp
     pxtm1_winter(1:kproma,:,krow)=0._dp
     pxtte(1:kproma,:,idx_xtrans_winter)=0._dp
     CALL message ('tracer_reset','winter tracer reset!')
  END IF

  IF (next_date == start_emi_summer_1 .OR.                                            &
      ((next_date < start_emi_summer_1) .AND. (nnext_date > start_emi_winter_1)) .OR. &
      next_date == start_emi_summer_2 .OR.                                            &
      ((next_date < start_emi_summer_2) .AND. (nnext_date > start_emi_winter_2)) .OR. &
      next_date == start_emi_summer_3 .OR.                                            &
      ((next_date < start_emi_summer_3) .AND. (nnext_date > start_emi_winter_3)) .OR. &
      next_date == start_emi_summer_4 .OR.                                            &
      ((next_date < start_emi_summer_4) .AND. (nnext_date > start_emi_winter_4))      ) THEN
     pxt_summer(1:kproma,:,krow)=0._dp
     pxtm1_summer(1:kproma,:,krow)=0._dp
     pxtte(1:kproma,:,idx_xtrans_summer)=0._dp
     CALL message ('tracer_reset','summer tracer reset!')
  END IF

END SUBROUTINE tracer_reset

  !-----------------------------------------------------------------------------
  !>
  !! @brief reset time filtered tracer concentrations, if specified (via namelist) date is reached
  !
SUBROUTINE tf_reset
!
! !DESCRIPTION:
! reset time filtered tracer concentrations, if specified (via namelist) date is reached
!
! !REVISION HISTORY:
! original source by F. Bunzel (2010-10-04)
!
! !USES:
  USE mo_memory_g1b,           ONLY: xtf
  USE mo_time_control,         ONLY: next_date, echam_time, get_time_step
  USE mo_time_manager,         ONLY: manager_state
  USE mo_time_conversion,      ONLY: OPERATOR(==), OPERATOR(<), OPERATOR(>)

  ! LOCAL VARIABLES:
  TYPE(TIME_DAYS)              :: nnext_date
  INTEGER                      :: istep
  
  istep = get_time_step()
  CALL manager_state(echam_time,nnext_date,istep+2)
  IF (next_date == start_emi_winter_1   .OR.                                          &
      ((next_date < start_emi_winter_1) .AND. (nnext_date > start_emi_winter_1)) .OR. &
      next_date == start_emi_winter_2   .OR.                                          &
      ((next_date < start_emi_winter_2) .AND. (nnext_date > start_emi_winter_2)) .OR. &
      next_date == start_emi_winter_3   .OR.                                          &
      ((next_date < start_emi_winter_3) .AND. (nnext_date > start_emi_winter_3)) .OR. &
      next_date == start_emi_winter_4   .OR.                                          &
      ((next_date < start_emi_winter_4) .AND. (nnext_date > start_emi_winter_4))      )THEN
     xtf(:,:,idx_xtrans_winter,:) = 0._dp
  END IF

  IF (next_date == start_emi_summer_1 .OR.                                            &
      ((next_date < start_emi_summer_1) .AND. (nnext_date > start_emi_winter_1)) .OR. &
      next_date == start_emi_summer_2 .OR.                                            &
      ((next_date < start_emi_summer_2) .AND. (nnext_date > start_emi_winter_2)) .OR. &
      next_date == start_emi_summer_3 .OR.                                            &
      ((next_date < start_emi_summer_3) .AND. (nnext_date > start_emi_winter_3)) .OR. &
      next_date == start_emi_summer_4 .OR.                                            &
      ((next_date < start_emi_summer_4) .AND. (nnext_date > start_emi_winter_4))      ) THEN
     xtf(:,:,idx_xtrans_summer,:) = 0._dp
  END IF

END SUBROUTINE tf_reset

END MODULE mo_aoa
