!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cbalone_dynveg

  USE mo_kind,          ONLY: dp
  USE mo_land_surface,  ONLY: land_surface_type
  USE mo_exception,     ONLY: message

  IMPLICIT none

  PUBLIC :: init_offline_dynveg
  PUBLIC :: update_offline_climbuf

  PRIVATE



CONTAINS

  SUBROUTINE init_offline_dynveg(grid, lctlib, surface, dynveg_clim, climbuf, use_dynveg, use_disturbance)

    USE mo_dynveg,         ONLY: dynveg, dynveg_options, cover_fract_pot_to_cover_fract, &
                                 calc_veg_ratio_max, scale_fpc, act_fpc_min
    USE mo_cbalone_memory, ONLY: dynveg_clim_type, check_err
    USE mo_land_surface,   ONLY: fract_small
    USE mo_jsbach_lctlib,  ONLY: lctlib_type
    USE mo_jsbach_grid,    ONLY: grid_type
    USE mo_jsbach,         ONLY: debug_Cconservation
    USE mo_climbuf,        ONLY: climbuf_type
    USE mo_disturbance,    ONLY: dist_opts

    INCLUDE 'netcdf.inc'

    TYPE(grid_type),               INTENT(in)    :: grid
    TYPE(lctlib_type),             INTENT(in)    :: lctlib
    TYPE(land_surface_type),       INTENT(inout) :: surface
    TYPE(dynveg_clim_type),        INTENT(inout) :: dynveg_clim
    TYPE(climbuf_type),            INTENT(inout) :: climbuf
    logical,                       INTENT(in)    :: use_dynveg
    logical,                       INTENT(in)    :: use_disturbance

    REAL(dp), ALLOCATABLE :: hlp3d(:,:,:), hlp4d(:,:,:,:)
    INTEGER,  ALLOCATABLE :: ndimlen(:), nvartyp(:), nvardim(:)
    INTEGER,  ALLOCATABLE :: nvardimid(:,:), nvaratt(:)
    CHARACTER(LEN=50), ALLOCATABLE :: fdimnam(:), fvarnam(:)
    INTEGER :: iret, ndim, nvar, i, idim, itile, ncid
    INTEGER :: nattr, nunlimdim
    INTEGER :: nin(100)

    IF (use_dynveg) THEN
      ALLOCATE(dynveg%act_fpc(grid%nland,surface%ntiles))
      ALLOCATE(dynveg%pot_fpc(grid%nland,surface%ntiles))
      ALLOCATE(dynveg%bare_fpc(grid%nland))
      ALLOCATE(dynveg%desert_fpc(grid%nland))
      ALLOCATE(dynveg%max_green_bio(grid%nland,surface%ntiles))
      ALLOCATE(dynveg%sum_green_bio_memory(grid%nland))

      IF (debug_Cconservation) THEN
         ALLOCATE(dynveg%dynveg_testCconserv_1(grid%nland))
         ALLOCATE(dynveg%dynveg_testCconserv_2(grid%nland))
         dynveg%dynveg_testCconserv_1(:) = 0._dp
         dynveg%dynveg_testCconserv_2(:) = 0._dp
      END IF
    ENDIF

    IF (use_dynveg .OR. use_disturbance) THEN
       ALLOCATE(dynveg%bio_exist(grid%nland,surface%ntiles))

       ALLOCATE(dynveg_clim%ave_npp5(grid%nland,surface%ntiles,31))
       ALLOCATE(dynveg_clim%bio_exist(grid%nland,surface%ntiles,31))
       ALLOCATE(dynveg_clim%max_wind10(grid%nland,31))
       ALLOCATE(dynveg_clim%prev_day_max_wind10(grid%nland,31))
    ENDIF
    NULLIFY(dynveg_clim%prev_day_mean_wind10,dynveg_clim%rel_hum_air)
    NULLIFY(dynveg_clim%prev_day_temp_max,dynveg_clim%prev_day_temp_min,dynveg_clim%prev_day_precip_mean) 
    !NULLIFY(dynveg_clim%dew_point_temp)
    NULLIFY(dynveg_clim%prev_day_mean_vol_moist)
    IF (use_disturbance) THEN
       IF (dist_opts%fire_algorithm /= 3) then
          ALLOCATE(dynveg_clim%rel_hum_air(grid%nland,31))
       ELSE
          ALLOCATE(dynveg_clim%prev_day_mean_wind10   (grid%nland,31), &
                   dynveg_clim%prev_day_precip_mean   (grid%nland,31), &
                   dynveg_clim%prev_day_temp_min      (grid%nland,31), &
                   dynveg_clim%prev_day_temp_max      (grid%nland,31), &
!                   dynveg_clim%dew_point_temp         (grid%nland,31),&
                   dynveg_clim%prev_day_mean_vol_moist(grid%nland,31))
          dynveg_clim%prev_day_mean_wind10   (:,:)   = 0._dp
          dynveg_clim%prev_day_precip_mean   (:,:)   = 0._dp
          dynveg_clim%prev_day_mean_vol_moist(:,:)   = 0._dp
       ENDIF
    ENDIF
    IF (use_dynveg .AND. .NOT. ASSOCIATED(dynveg_clim%rel_hum_air)) &
      ALLOCATE(dynveg_clim%rel_hum_air(grid%nland,31))

    IF (use_dynveg .OR. use_disturbance) THEN
       dynveg_clim%bio_exist(:,:,:) = 1._dp

       climbuf%ave_npp5                => dynveg_clim%ave_npp5(:,:,1)
       IF (ASSOCIATED(dynveg_clim%rel_hum_air)) climbuf%rel_hum_air => dynveg_clim%rel_hum_air(:,1)
       climbuf%max_wind10              => dynveg_clim%max_wind10(:,1)
       climbuf%prev_day_max_wind10     => dynveg_clim%prev_day_max_wind10(:,1)
       IF (use_disturbance) THEN
         IF (ASSOCIATED(dynveg_clim%rel_hum_air)) THEN
           climbuf%rel_hum_air            => dynveg_clim%rel_hum_air(:,1)
         ELSE
           climbuf%rel_hum_air            => dynveg_clim%max_wind10(:,1) ! Data not used, but pointer must be assigned
         ENDIF
       ENDIF
    ENDIF
    IF (use_disturbance .AND. dist_opts%fire_algorithm == 3) THEN
      climbuf%prev_day_mean_wind10    => dynveg_clim%prev_day_mean_wind10(:,1)
      climbuf%prev_day_temp_max       => dynveg_clim%prev_day_temp_max(:,1)
      climbuf%prev_day_temp_min       => dynveg_clim%prev_day_temp_min(:,1)
      climbuf%prev_day_precip_mean    => dynveg_clim%prev_day_precip_mean(:,1)
!      climbuf%dew_point_temp          => dynveg_clim%dew_point_temp(:,1)
      climbuf%prev_day_mean_vol_moist => dynveg_clim%prev_day_mean_vol_moist(:,1)
    ENDIF

    IF (.NOT. use_disturbance) THEN
      ALLOCATE(dynveg%Dummy_null(grid%nland,surface%ntiles))
      dynveg%Dummy_null(:,:) = 0._dp
    ENDIF

    IF (use_dynveg .AND. dynveg_options%read_fpc) THEN    !! Read FPCs from file

       CALL message ('init_offline_dynveg','Initial FPCs read from file '//TRIM(dynveg_options%fpc_file_name))
       iret = nf_open(dynveg_options%fpc_file_name,nf_nowrite,ncid)
       CALL check_err(iret)

       ! check what is in the input-file
       iret = nf_inq(ncid,ndim,nvar,nattr,nunlimdim)
       CALL check_err(iret)

       ! get the dimension name and length
       ALLOCATE (ndimlen(ndim))
       ALLOCATE (fdimnam(ndim))
       DO i=1,ndim
          iret = nf_inq_dim(ncid,i,fdimnam(i),ndimlen(i))
          call check_err(iret)
       END DO

       ! get variable names, types and shapes
       ALLOCATE (fvarnam(nvar))
       ALLOCATE (nvartyp(nvar))
       ALLOCATE (nvardim(nvar))
       ALLOCATE (nvardimid(nvar,100))
       ALLOCATE (nvaratt(nvar))
       DO i=1,nvar
          iret = nf_inq_var(ncid,i,fvarnam(i),nvartyp(i),nvardim(i),nin,nvaratt(i))
          CALL check_err(iret)
          DO idim=1,nvardim(i)
             nvardimid(i,idim)=nin(idim)
          END DO
       END DO

       ! get the data
       DO i=1,nvar

          SELECT CASE (fvarnam(i))
          CASE ('act_fpc', 'max_green_bio', 'cover_fract', 'cover_fract_pot', 'pot_fpc')
             ALLOCATE(hlp4d(ndimlen(nvardimid(i,1)), ndimlen(nvardimid(i,2)), &
                  ndimlen(nvardimid(i,3)), 1))
             IF (ndimlen(nvardimid(i,1)) /= grid%nlon .OR. ndimlen(nvardimid(i,2)) /= grid%nlat &
                  .OR.  ndimlen(nvardimid(i,3)) /= surface%ntiles) THEN
                WRITE(*,*) 'init_offline_dynveg: dimensions of ', TRIM(fvarnam(i)), ' are wrong'
                STOP
             END IF
             IF (ndim==3) THEN
                iret = NF_GET_VARA_DOUBLE(ncid,i,(/1,1,1,1/), &
                     SHAPE(hlp4d), hlp4d(:,:,:,:))
             ELSE
                iret = NF_GET_VARA_DOUBLE(ncid,i,(/1,1,1,ndimlen(nvardimid(i,4))/), &
                     SHAPE(hlp4d), hlp4d(:,:,:,:))
             END IF
             CALL check_err(iret)

             SELECT CASE (fvarnam(i))
             CASE ('act_fpc')
                DO itile = 1,surface%ntiles
                   dynveg%act_fpc(:,itile) = PACK(hlp4d(:,:,itile,1),MASK=grid%mask)
                END DO
             CASE ('max_green_bio')
                DO itile = 1,surface%ntiles
                   dynveg%max_green_bio(:,itile) = PACK(hlp4d(:,:,itile,1),MASK=grid%mask)
                END DO
             CASE ('cover_fract')
                DO itile = 1,surface%ntiles
                   surface%cover_fract(:,itile) = PACK(hlp4d(:,:,itile,1),MASK=grid%mask)
                END DO
             CASE ('cover_fract_pot')
                DO itile = 1,surface%ntiles
                   surface%cover_fract_pot(:,itile) = PACK(hlp4d(:,:,itile,1),MASK=grid%mask)
                END DO
             CASE ('pot_fpc')
                DO itile = 1,surface%ntiles
                   dynveg%pot_fpc(:,itile) = PACK(hlp4d(:,:,itile,1),MASK=grid%mask)
                END DO
             END SELECT

             DEALLOCATE(hlp4d)
             CALL message('init_offline_dynveg','read '//TRIM(fvarnam(i)))
          CASE ('bare_fpc', 'desert_fpc', 'sum_green_bio_memory')
             ALLOCATE(hlp3d(ndimlen(nvardimid(i,1)), ndimlen(nvardimid(i,2)), 1))
             IF (ndimlen(nvardimid(i,1)) /= grid%nlon .OR. ndimlen(nvardimid(i,2)) /= grid%nlat) THEN
                WRITE(*,*) 'init_offline_dynveg: dimensions of ', TRIM(fvarnam(i)), ' are wrong'
                STOP
             END IF
             IF (ndim==3) THEN
                iret = NF_GET_VARA_DOUBLE(ncid,i,(/1,1,1/), &
                     SHAPE(hlp3d), hlp3d(:,:,:))
             ELSE
                iret = NF_GET_VARA_DOUBLE(ncid,i,(/1,1,ndimlen(nvardimid(i,3))/), &
                     SHAPE(hlp3d), hlp3d(:,:,:))
             END IF
             CALL check_err(iret)

             SELECT CASE (fvarnam(i))
             CASE ('bare_fpc')
                dynveg%bare_fpc(:) = PACK(hlp3d(:,:,1), MASK=grid%mask)
             CASE ('desert_fpc')
                dynveg%desert_fpc(:) = PACK(hlp3d(:,:,1), MASK=grid%mask)
                IF (dynveg_options%dynveg_feedback) THEN
                   surface%veg_ratio_max(:) = MAX(fract_small,1._dp - dynveg%desert_fpc(:))
                END IF
             CASE ('sum_green_bio_memory')
                dynveg%sum_green_bio_memory(:) = PACK(hlp3d(:,:,1), MASK=grid%mask)
             END SELECT

             DEALLOCATE(hlp3d)
             CALL message('init_offline_dynveg','read '//TRIM(fvarnam(i)))
          END SELECT
       END DO

       ! close the file
       iret = nf_close(ncid)
       CALL check_err(iret)

       ! deallocate
       DEALLOCATE(fdimnam,ndimlen)
       DEALLOCATE(fvarnam,nvartyp,nvardim,nvardimid,nvaratt)
    ELSE                  !! Initialize FPCs from echam cover fractions

       IF (use_dynveg) THEN
         dynveg%act_fpc(:,:) = surface%cover_fract_pot(:,:)
         dynveg%bare_fpc(:) = 0._dp
         dynveg%desert_fpc(:) = MIN(1._dp - 2._dp * EPSILON(1._dp),MAX(0._dp,1._dp - surface%veg_ratio_max(:)))
         dynveg%max_green_bio(:,:) = 0._dp
         dynveg%sum_green_bio_memory(:) = MIN(1._dp,MAX(0._dp,surface%veg_ratio_max(:)))

         WHERE (ANY(surface%is_glacier(:,:),DIM=2))
            surface%rock_fract(:) = 1._dp
         END WHERE
         DO itile=1,surface%ntiles
            WHERE (surface%rock_fract(:) >= 1._dp - EPSILON(1._dp))
               dynveg%act_fpc(:,itile) = 0._dp
            END WHERE
         END DO

         CALL scale_fpc(grid%nland, surface%ntiles, ANY(surface%is_glacier(:,:),DIM=2), lctlib%dynamic_pft(:), &
              surface%cover_type(:,:), dynveg%act_fpc(:,:), dynveg%bare_fpc(:))

         dynveg%pot_fpc(:,:) = MAX(act_fpc_min,dynveg%act_fpc(:,:))

         IF (dynveg_options%dynveg_feedback) THEN

            ! calculate cover fractions from FPCs to get consistent veg_ratio_max

            CALL calc_veg_ratio_max(dynveg%desert_fpc(:), surface%rock_fract(:), surface%is_glacier(:,:), &
                                    surface%veg_ratio_max(:))
            CALL cover_fract_pot_to_cover_fract(grid%nland, surface%ntiles, lctlib%nlct, surface%is_present(:,:), &
                 surface%is_glacier(:,:), surface%is_naturalveg(:,:), lctlib%woody_pft(:), lctlib%dynamic_pft(:), &
                 lctlib%pasture_pft(:), surface%is_crop(:,:), surface%cover_type(:,:), surface%cover_fract_pot(:,:), &
                 surface%cover_fract_pot(:,:), surface%cover_fract(:,:))

         END IF
       END IF
    END IF

  END SUBROUTINE init_offline_dynveg

!--------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE update_offline_climbuf(day, dynveg_clim, climbuf, dynveg)
!--------------------------------------------------------------------------------------------------------------------------
    USE mo_cbalone_memory, ONLY: dynveg_clim_type
    USE mo_climbuf,        ONLY: climbuf_type
    USE mo_dynveg,         ONLY: dynveg_type
    USE mo_disturbance,    ONLY: dist_opts

    INTEGER,                       INTENT(in)    :: day
    TYPE(dynveg_clim_type),        INTENT(inout) :: dynveg_clim
    TYPE(climbuf_type),            INTENT(inout) :: climbuf
    TYPE(dynveg_type),             INTENT(inout) :: dynveg

    climbuf%ave_npp5                => dynveg_clim%ave_npp5(:,:,day)
    IF (dist_opts%fire_algorithm /= 3) THEN
      climbuf%rel_hum_air => dynveg_clim%rel_hum_air(:,day)
    ELSE
      climbuf%rel_hum_air => dynveg_clim%max_wind10(:,day) ! Data not used, but pointer needs to be assigned to something
    ENDIF
    climbuf%prev_day_max_wind10     => dynveg_clim%prev_day_max_wind10(:,day)
    climbuf%max_wind10              => dynveg_clim%max_wind10(:,day)
    IF (dist_opts%fire_algorithm == 3) THEN
      climbuf%prev_day_mean_wind10    => dynveg_clim%prev_day_mean_wind10(:,day)
      climbuf%prev_day_temp_min       => dynveg_clim%prev_day_temp_min(:,day)
      climbuf%prev_day_temp_max       => dynveg_clim%prev_day_temp_max(:,day)
!      climbuf%dew_point_temp          => dynveg_clim%dew_point_temp(:,day)
      climbuf%prev_day_precip_mean    => dynveg_clim%prev_day_precip_mean(:,day)
      climbuf%prev_day_mean_vol_moist => dynveg_clim%prev_day_mean_vol_moist(:,day)
     ENDIF

    dynveg%bio_exist(:,:) = dynveg_clim%bio_exist(:,:,day)

  END SUBROUTINE update_offline_climbuf

END MODULE mo_cbalone_dynveg
