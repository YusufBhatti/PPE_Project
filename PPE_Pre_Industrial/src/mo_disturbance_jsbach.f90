!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_disturbance_jsbach
!
! Calculating disturbances (fire and windbreak) using the jsbach algorithms
!
! History:
! 2011-04-29 - StW - Original code obtained by cut'n'paste from mo_dynveg and mo_cbal_cpools

! Control parameters for jsbach disturbance algorithms (to avoid circular
! dependencies between mo_disturbance and mo_disturbance_jsbach).

  USE mo_land_surface,  ONLY: fract_small
  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: finish, message
  USE mo_jsbach_lctlib, ONLY: lctlib_type
IMPLICIT NONE

! Fire namelist parameters
  REAL(dp), PARAMETER, PRIVATE :: def_fire_litter_threshold  = 16.67_dp   !! minimal amount of litter [mol(C)/m^2(grid box)] for 
                                                                          !! fire
  REAL(dp), PARAMETER, PRIVATE :: def_fire_rel_hum_threshold = 70._dp     !! maximal relative humidity for fire
  REAL(dp), PARAMETER, PRIVATE :: def_fire_minimum_woody     =  0.002_dp  !! minimal fraction of act_fpc of woody PFT to be burned
                                                                          !! each year
  REAL(dp), PARAMETER, PRIVATE :: def_fire_minimum_grass     =  0.006_dp  !! minimal fraction of act_fpc of grass PFT to be burned
                                                                          !! each year
  REAL(dp), PARAMETER, PRIVATE :: def_fire_tau_woody         =  6.0_dp    !! return period of fire for woody PFT [year] at 0% 
                                                                          !! relative humidity
  REAL(dp), PARAMETER, PRIVATE :: def_fire_tau_grass         =  2.0_dp    !! return period of fire for grass PFT [year] at 0% 
                                                                          !! relative humidity
! Wind break namelist parameters
  REAL(dp), PARAMETER, PRIVATE :: def_wind_threshold         =  2.25_dp   !! factor by which the previous day maximum wind speed
                                                                          !! must be larger than the climatological daily maximum
                                                                          !! wind speed to allow any wind damage
  REAL(dp), PARAMETER, PRIVATE :: def_wind_damage_scale      =  5.e-03_dp !! scaling factor for wind damage

! Actual fire parameters
  REAL(dp), PRIVATE :: litter_threshold
  REAL(dp), PRIVATE :: rel_hum_threshold
  REAL(dp), PRIVATE, DIMENSION(:), ALLOCATABLE :: fire_minimum
  REAL(dp), PRIVATE, DIMENSION(:), ALLOCATABLE :: fire_tau

! Actual wind parameters
  REAL(dp), PRIVATE :: wnd_threshold
  REAL(dp), PRIVATE :: wnd_damage_scale

! Other globals
  INTEGER , PRIVATE :: module_status = 0

! Function declarations

PUBLIC config_fire_jsbach
PUBLIC config_windbreak_jsbach
PUBLIC init_fire_jsbach
PUBLIC init_windbreak_jsbach
PUBLIC burned_frac_jsbach
PUBLIC broken_woody_frac_jsbach

CONTAINS

  SUBROUTINE config_fire_jsbach(lctlib,nml_unit)
    USE mo_namelist,         ONLY: position_nml, POSITIONED, LENGTH_ERROR, READ_ERROR
    USE mo_io_units,         ONLY: nout
    USE mo_mpi,              ONLY: p_parallel_io, p_io, p_parallel, p_bcast

! Input parameters
    TYPE(lctlib_type), INTENT(IN) :: lctlib
    INTEGER,           INTENT(IN) :: nml_unit

! Local namelist parameters
    REAL(dp)      :: fire_litter_threshold 
    REAL(dp)      :: fire_rel_hum_threshold
    REAL(dp)      :: fire_minimum_woody
    REAL(dp)      :: fire_minimum_grass
    REAL(dp)      :: fire_tau_woody
    REAL(dp)      :: fire_tau_grass

! Other locals
    INTEGER  :: read_status
    REAL(dp) :: rbuf(6)
    INTEGER  :: f_unit

    INCLUDE 'fire_jsbach_ctl.inc'

    IF (IAND(module_status,1) /= 0) THEN
      CALL message('config_fire_jsbach','Warning: Re-initializing module')
      DEALLOCATE(fire_minimum,fire_tau)
    ENDIF

    IF (p_parallel_io) THEN
! Map default values
      fire_litter_threshold  = def_fire_litter_threshold          
      fire_rel_hum_threshold = def_fire_rel_hum_threshold
      fire_minimum_woody     = def_fire_minimum_woody
      fire_minimum_grass     = def_fire_minimum_grass
      fire_tau_woody         = def_fire_tau_woody
      fire_tau_grass         = def_fire_tau_grass

! Read namelist if available
      f_unit = position_nml ('FIRE_JSBACH_CTL', nml_unit, status=read_status)
      SELECT CASE (read_status)
        CASE (POSITIONED)
          READ (f_unit, fire_jsbach_ctl)
          CALL message('config_fire_jsbach','Namelist FIRE_JSBACH_CTL: ')
          WRITE(nout, fire_jsbach_ctl)
        CASE (LENGTH_ERROR)
          CALL finish ('config_fire_jsbach','Length error in namelist fire_jsbach_ctl')
        CASE (READ_ERROR)
          CALL finish ('config_fire_jsbach','Error reading namelist fire_jsbach_ctl ')
        END SELECT

! Distribute namelist parameters
      rbuf = (/fire_litter_threshold,                                    &
               fire_rel_hum_threshold,fire_minimum_woody,                &
               fire_minimum_grass,fire_tau_woody,fire_tau_grass/)
    ENDIF

    IF (p_parallel) THEN
      CALL p_bcast(rbuf,p_io)
    ENDIF
     
    litter_threshold       = rbuf( 1)
    rel_hum_threshold      = rbuf( 2)

    ALLOCATE(fire_minimum(lctlib%nlct), &
             fire_tau    (lctlib%nlct))

    WHERE(lctlib%woody_pft(:))
      fire_minimum(:) = rbuf(3)
      fire_tau(:)     = rbuf(5)
    ELSEWHERE
      fire_minimum(:) = rbuf(4)
      fire_tau(:)     = rbuf(6)
    ENDWHERE

  END SUBROUTINE config_fire_jsbach

  SUBROUTINE config_windbreak_jsbach(nml_unit)
    USE mo_namelist,         ONLY: position_nml, POSITIONED, LENGTH_ERROR, READ_ERROR
    USE mo_io_units,         ONLY: nout
    USE mo_mpi,              ONLY: p_parallel_io, p_io, p_parallel, p_bcast

! Input parameters
    INTEGER, INTENT(IN) :: nml_unit

! Local namelist parameters
    REAL(dp) :: wind_threshold
    REAL(dp) :: wind_damage_scale

! Other locals
    INTEGER  :: read_status
    INTEGER  :: f_unit

    INCLUDE 'windbreak_jsbach_ctl.inc'

    IF (IAND(module_status,2) /= 0) CALL message('config_windbreak_jsbach','Re-initializing module')

    IF (p_parallel_io) THEN
! Map default values
      wind_threshold    = def_wind_threshold
      wind_damage_scale = def_wind_damage_scale

! Read namelist if available
      f_unit = position_nml ('WINDBREAK_JSBACH_CTL', nml_unit, status=read_status)
      SELECT CASE (read_status)
        CASE (POSITIONED)
          READ (f_unit, windbreak_jsbach_ctl)
          CALL message('config_windbreak_jsbach','Namelist WINDBREAK_JSBACH_CTL: ')
          WRITE(nout, windbreak_jsbach_ctl)
       CASE (LENGTH_ERROR)
          CALL finish ('config_windbreak_jsbach','Length error in namelist windbreak_jsbach_ctl')
       CASE (READ_ERROR)
          CALL finish ('config_windbreak_jsbach','Error reading namelist windbreak_jsbach_ctl ')
       END SELECT
    ENDIF

! Distribute namelist parameters
    IF (p_parallel) THEN
      CALL p_bcast(wind_threshold,p_io)
      CALL p_bcast(wind_damage_scale,p_io)
    ENDIF
    
    wnd_threshold    = wind_threshold
    wnd_damage_scale = wind_damage_scale

  END SUBROUTINE config_windbreak_jsbach

  SUBROUTINE init_fire_jsbach(lctlib,nml_unit)
    TYPE(lctlib_type), INTENT(IN) :: lctlib
    INTEGER,           INTENT(IN) :: nml_unit

    CALL config_fire_jsbach(lctlib,nml_unit)

    module_status = IOR(module_status,1)

  END SUBROUTINE init_fire_jsbach

  SUBROUTINE init_windbreak_jsbach(nml_unit)
    INTEGER, INTENT(IN) :: nml_unit

    CALL config_windbreak_jsbach(nml_unit)

    module_status = IOR(module_status,2)

  END SUBROUTINE init_windbreak_jsbach

! Calculate burned area. Woody and grass types call the subroutine with different parameters 
  SUBROUTINE burned_frac_jsbach(nidx,ntiles,with_yasso,used_pft,is_present,                &
                                cover_type,cover_fract,veg_fract_correction,veg_ratio_max, &
                                rel_hum_air,fuel,burned_frac,                              &
                                litter_green_ag,litter_wood_ag,                            &
                                YCpool_acid_ag1,YCpool_water_ag1,                          &
                                YCpool_ethanol_ag1,YCpool_nonsoluble_ag1,                  &
                                YCpool_acid_ag2,YCpool_water_ag2,                          &
                                YCpool_ethanol_ag2,YCpool_nonsoluble_ag2)
  USE mo_utils, ONLY: average_tiles

! Input parameters
    INTEGER                 , INTENT(IN)  :: nidx
    INTEGER                 , INTENT(IN)  :: ntiles
    LOGICAL                 , INTENT(IN)  :: with_yasso
    LOGICAL , DIMENSION(  :), INTENT(IN)  :: used_pft
    LOGICAL , DIMENSION(:,:), INTENT(IN)  :: is_present
    INTEGER , DIMENSION(:,:), INTENT(IN)  :: cover_type
    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: cover_fract
    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: veg_fract_correction 
    REAL(dp), DIMENSION(:  ), INTENT(IN)  :: veg_ratio_max
    REAL(dp), DIMENSION(  :), INTENT(IN)  :: rel_hum_air
    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: litter_green_ag
    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: litter_wood_ag
! Optional input parameters    

    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: Ycpool_acid_ag1
    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: YCpool_water_ag1
    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: YCpool_ethanol_ag1
    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: YCpool_nonsoluble_ag1
    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: Ycpool_acid_ag2
    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: YCpool_water_ag2
    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: YCpool_ethanol_ag2
    REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: YCpool_nonsoluble_ag2

! Output parameters
    REAL(dp), DIMENSION(:)  , INTENT(OUT)   :: fuel
    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: burned_frac

! Locals
    INTEGER  :: itile, i, ct
    REAL(dp) :: delta_time_yr
    
    delta_time_yr = 1._dp / 365._dp
    IF (with_yasso) THEN
      CALL average_tiles( (YCpool_acid_ag1   (:,:)   + YCpool_water_ag1     (:,:)    &
                           + YCpool_ethanol_ag1(:,:) + YCpool_nonsoluble_ag1(:,:)    &
                           + YCpool_acid_ag2   (:,:) + YCpool_water_ag2     (:,:)    &
                           + YCpool_ethanol_ag2(:,:) + YCpool_nonsoluble_ag2(:,:)    &
                          ) * veg_fract_correction(:,:) *                &
                          SPREAD(veg_ratio_max(:),NCOPIES=ntiles,DIM=2), &
                          is_present(:,:),cover_fract(:,:), fuel(:))
    ELSE
      CALL average_tiles( (litter_green_ag(:,:) + litter_wood_ag(:,:) ) * veg_fract_correction(:,:) * &
                          SPREAD(veg_ratio_max(:),NCOPIES=ntiles,DIM=2), &
                          is_present(:,:),cover_fract(:,:), fuel(:))
    END IF


    DO itile = 1,ntiles
      DO i = 1,nidx
        ct = cover_type(i,itile)
        IF (used_pft(ct)) THEN
          burned_frac(i,itile) = fire_minimum(ct) * delta_time_yr
          IF (rel_hum_air(i) < rel_hum_threshold .AND. fuel(i) > litter_threshold) &
            burned_frac(i,itile) = (fire_minimum(ct) + ((rel_hum_threshold - rel_hum_air(i)) /  &
                                   (fire_tau(ct) * rel_hum_threshold))) * delta_time_yr
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE burned_frac_jsbach

! Calculate windbreak for woody types
  SUBROUTINE broken_woody_frac_jsbach(nidx,ntiles,used_pft,act_fpc,cover_type, &
                                      wind_prev_day,wind_clim,damaged_frac)

! Input parameters
    INTEGER                 , INTENT(IN)  :: nidx
    INTEGER                 , INTENT(IN)  :: ntiles
    LOGICAL , DIMENSION(  :), INTENT(IN)  :: used_pft
    INTEGER , DIMENSION(:,:), INTENT(IN)  :: cover_type
    REAL(dp), DIMENSION(:,:), INTENT(IN)  :: act_fpc
    REAL(dp), DIMENSION(  :), INTENT(IN)  :: wind_prev_day
    REAL(dp), DIMENSION(  :), INTENT(IN)  :: wind_clim

! Output parameters
    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: damaged_frac

! Locals
    INTEGER  :: itile, i
    REAL(dp) :: delta_time_yr

    delta_time_yr = 1._dp / 365._dp
    DO itile = 1,ntiles
      DO i = 1,nidx
        IF (used_pft(cover_type(i,itile)) .AND. (act_fpc(i,itile) > fract_small)) THEN
          IF (wind_prev_day(i) > wind_clim(i) * wnd_threshold)                                                       &
            damaged_frac(i,itile) = MIN(1._dp-fract_small/act_fpc(i,itile),                                          &
                                        delta_time_yr * wnd_damage_scale * wind_prev_day(i) ** 3._dp / wind_clim(i))
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE broken_woody_frac_jsbach

END MODULE mo_disturbance_jsbach
