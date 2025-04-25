!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7_nucl_diag.f90
!!
!! \brief
!! mo_ham_nucl_diag provides routines for sampling the cloud-free volume as
!! function of temperature, RH, H2SO4 concentration, H2SO4 condensation sink,
!! and ion pair production rate.
!!
!! \author Jan Kazil (MPI-Met)
!!
!! \responsible_coder
!! Jan Kazil, Jan.Kazil@noaa.gov
!!
!! \revision_history
!!   -# J. Kazil (MPI-Met) - original code (2008-02-21)
!!   -# L. Kornblueh (MPI-Met) (2008-09)
!!
!! \limitations
!! None
!!
!! \details
!! mo_ham_nucl_diag provides routines for sampling the cloud-free volume as
!! function of temperature, RH, H2SO4 concentration, H2SO4 condensation sink,
!! and ion pair production rate. These are the variables that matter for
!! sulfate aerosol nucleation. The volume as function of these variables is
!! saved in a netCDF archive. The volume sampling is memory intensive.
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

MODULE mo_ham_m7_nucl_diag 

  USE mo_kind,            ONLY: dp
  USE mo_filename,        ONLY: out_datapath
  USE mo_io_units,        ONLY: nout,nerr
  USE mo_filename,        ONLY: standard_output_file
  USE mo_time_control,    ONLY: current_date
  USE mo_time_conversion, ONLY: TC_get,TC_convert,time_native
  USE mo_mpi,             ONLY: p_parallel,p_parallel_io,p_io,p_gather,p_nprocs
  USE mo_exception,       ONLY: message, finish
  USE mo_netcdf
!!  USE mo_netcdf,          ONLY: nf_check

  IMPLICIT NONE

!!   INCLUDE 'netcdf.inc'

  CHARACTER(512):: file

  !
  ! Variable names and units
  !

  CHARACTER(*), PRIVATE, PARAMETER :: time_name = 'time'
  CHARACTER(*), PRIVATE, PARAMETER :: time_units = 'day as %Y%m%d.%f'

  CHARACTER(*), PRIVATE, PARAMETER :: temp_name = 'temperature'
  CHARACTER(*), PRIVATE, PARAMETER :: temp_units = 'K'
  CHARACTER(*), PRIVATE, PARAMETER :: temp_long_name = 'temperature'

  CHARACTER(*), PRIVATE, PARAMETER :: rh_name = 'rh'
  CHARACTER(*), PRIVATE, PARAMETER :: rh_units = '%'
  CHARACTER(*), PRIVATE, PARAMETER :: rh_long_name = 'relative humidity'

  CHARACTER(*), PRIVATE, PARAMETER :: h2so4_name = 'h2so4'
  CHARACTER(*), PRIVATE, PARAMETER :: h2so4_units = 'cm-3'
  CHARACTER(*), PRIVATE, PARAMETER :: h2so4_long_name = 'H2SO4(g) concentration'

  CHARACTER(*), PRIVATE, PARAMETER :: ipr_name = 'ipr'
  CHARACTER(*), PRIVATE, PARAMETER :: ipr_units = 'cm-3 s-1'
  CHARACTER(*), PRIVATE, PARAMETER :: ipr_long_name = 'ion pair production rate'

  CHARACTER(*), PRIVATE, PARAMETER :: cs_name = 'cs'
  CHARACTER(*), PRIVATE, PARAMETER :: cs_units = 's-1'
  CHARACTER(*), PRIVATE, PARAMETER :: cs_long_name = 'H2SO4 condensation sink'

  CHARACTER(*), PRIVATE, PARAMETER :: vol_name = 'vol_cf'
  CHARACTER(*), PRIVATE, PARAMETER :: vol_units = 'm3'
  CHARACTER(*), PRIVATE, PARAMETER :: vol_long_name = &
       TRIM(vol_name)//'('// &
       TRIM(cs_name)//'_i,'// &
       TRIM(ipr_name)//'_i,'// &
       TRIM(h2so4_name)//'_i,'// &
       TRIM(rh_name)//'_i,'// &
       TRIM(temp_name)//'_i)'// &
       ' is the accumulated, cloud-free volume with conditions in the intervals '// &
       '['//TRIM(temp_name)//'('//TRIM(temp_name)//'_i),'//TRIM(temp_name)//'('//TRIM(temp_name)//'_i+1)], '// &
       '['//TRIM(rh_name)//'('//TRIM(rh_name)//'_i),'//TRIM(rh_name)//'('//TRIM(rh_name)//'_i+1)], ' // &
       '['//TRIM(h2so4_name)//'('//TRIM(h2so4_name)//'_i),'//TRIM(h2so4_name)//'('//TRIM(h2so4_name)//'_i+1)], ' // &
       '['//TRIM(ipr_name)//'('//TRIM(ipr_name)//'_i),'//TRIM(ipr_name)//'('//TRIM(ipr_name)//'_i+1)], ' // &
       '['//TRIM(cs_name)//'('//TRIM(cs_name)//'_i),'//TRIM(cs_name)//'('//TRIM(cs_name)//'_i+1)]'

  !  
  ! Grid of variables that matter for sulfate aerosol nucleation
  !

  ! Resolution

  INTEGER, PRIVATE, PARAMETER :: temp_n = 41
  INTEGER, PRIVATE, PARAMETER :: rh_n = 41
  INTEGER, PRIVATE, PARAMETER :: h2so4_n = 41
  INTEGER, PRIVATE, PARAMETER :: ipr_n = 21
  INTEGER, PRIVATE, PARAMETER :: cs_n = 40

  ! Min/max values

  REAL(dp), PRIVATE :: temp_min,temp_max, &
       rh_min,rh_max, &
       h2so4_min,h2so4_max, &
       ipr_min,ipr_max, &
       cs_min,cs_max

  ! Grid points

  REAL(dp), PRIVATE, ALLOCATABLE :: temp_grid(:), &
       rh_grid(:), &
       h2so4_grid(:), &
       ipr_grid(:), &
       cs_grid(:)

  REAL(dp), PRIVATE, ALLOCATABLE :: vol(:,:,:,:,:)

CONTAINS

  SUBROUTINE ham_nucl_diag_initialize()

    ! *ham_nucl_diag_initialize* generates a multidimensional grid spanning the
    ! space of the variables that are relevant for sulfate aerosol nucleation
    ! and allocates memory for a corresponding array. The elements of this array
    ! will hold the cloud-free volume in the course of this run as function of
    ! the variables temperature, RH, H2SO4 concentration, H2SO4 condensation
    ! sink, and ion pair production rate. In parallel mode, the array is
    ! broadcast to all processes.
    !
    ! Local variables:

    REAL(dp), ALLOCATABLE :: ztemp_grid(:), &
         zrh_grid(:), &
         zh2so4_grid(:), &
         zipr_grid(:), &
         zcs_grid(:)

    REAL(dp):: zkappa

    INTEGER :: itemp, &
         irh, &
         ih2so4, &
         iipr, &
         ics

    ! Set up the grids:

    ALLOCATE(temp_grid(temp_n))
    ALLOCATE(rh_grid(rh_n))
    ALLOCATE(h2so4_grid(h2so4_n))
    ALLOCATE(ipr_grid(ipr_n))
    ALLOCATE(cs_grid(cs_n))

    ALLOCATE(vol(temp_n,rh_n,h2so4_n,ipr_n,cs_n))

    itemp = temp_n - 1
    irh = rh_n - 1
    ih2so4 = h2so4_n - 1
    iipr = ipr_n - 1
    ics = cs_n - 1

    ALLOCATE(ztemp_grid(itemp))
    ALLOCATE(zrh_grid(irh))
    ALLOCATE(zh2so4_grid(ih2so4))
    ALLOCATE(zipr_grid(iipr))
    ALLOCATE(zcs_grid(ics))

    vol = 0.0_dp

    ! Temperature:

    temp_min = 180.0_dp ! K
    temp_max = 320.0_dp ! K

    zkappa = -1.0_dp/3.0_dp

    CALL arched_nodes(temp_min,temp_max,itemp,ztemp_grid,zkappa)

    temp_grid(1)         = 0.0_dp
    temp_grid(2:temp_n) = ztemp_grid(1:itemp)

    ! RH

    rh_min =   1.0_dp ! % over water
    rh_max = 101.0_dp ! % over water

    zkappa = 2.0_dp

    CALL hyper_geometric_nodes_1d(rh_min,rh_max,irh,zrh_grid,zkappa)

    rh_grid(1)      = 0.0_dp
    rh_grid(2:rh_n) = zrh_grid(1:irh)

    ! H2SO4 concentration

    h2so4_min = 1.0E5_dp ! cm-3 s-1
    h2so4_max = 5.0E9_dp ! cm-3 s-1

    zkappa = 1.0_dp

    CALL hyper_geometric_nodes_1d(h2so4_min,h2so4_max,ih2so4,zh2so4_grid,zkappa)

    h2so4_grid(1)         = 0.0_dp
    h2so4_grid(2:h2so4_n) = zh2so4_grid(1:ih2so4)

    ! Ionization rate

    ipr_min =  1.0_dp ! cm-3 s-1
    ipr_max = 55.0_dp ! cm-3 s-1

    zkappa = 1.0_dp

    CALL hyper_geometric_nodes_1d(ipr_min,ipr_max,iipr,zipr_grid,zkappa)

    ipr_grid(1)       = 0.0_dp
    ipr_grid(2:ipr_n) = zipr_grid(1:iipr)

    ! H2SO4 condensation sink of preexisting aerosol particles

    cs_min = 1.75E-5_dp ! s-1
    cs_max = 1.00E-1_dp ! s-1

    zkappa = 1.0_dp

    CALL hyper_geometric_nodes_1d(cs_min,cs_max,cs_n,cs_grid,zkappa)

    cs_grid(1) = 0.0_dp

    DEALLOCATE(ztemp_grid)
    DEALLOCATE(zrh_grid)
    DEALLOCATE(zh2so4_grid)
    DEALLOCATE(zipr_grid)
    DEALLOCATE(zcs_grid)

  END SUBROUTINE ham_nucl_diag_initialize

  !=============================================================================

  SUBROUTINE create_nucl_statistics_file()

    ! *create_nucl_statistics_file* creates a netCDF archive for saving data and
    ! unless the archive already exists, in which case nothing is done.
    !
    ! Local variables:
    !

    LOGICAL :: ll_exists

    INTEGER :: i_nc_id

    INTEGER :: i_time_dim_id,&
         i_temp_dim_id,&
         i_rh_dim_id,&
         i_h2so4_dim_id,&
         i_ipr_dim_id,&
         i_cs_dim_id

    INTEGER :: i_time_id,&
         i_temp_id,&
         i_rh_id,&
         i_h2so4_id,&
         i_ipr_id,&
         i_cs_id, &
         i_vol_id

    INTEGER :: i_length
    INTEGER :: i_vector(6)

    !
    ! Create the netCDF archive. If it already exists then do nothing.
    !

    INQUIRE(FILE=file,EXIST=ll_exists)

    IF (ll_exists) RETURN

    CALL nf_check(nf_create(file,NF_NOCLOBBER,i_nc_id),fname=TRIM(file))

    !
    ! Define the dimensions on which the data are defined:
    !

    ! Time:
    CALL nf_check(nf_def_dim(i_nc_id,time_name,NF_UNLIMITED,i_time_dim_id))

    ! Temperature:
    CALL nf_check(nf_def_dim(i_nc_id,temp_name,temp_n,i_temp_dim_id))

    ! RH:
    CALL nf_check(nf_def_dim(i_nc_id,rh_name,rh_n,i_rh_dim_id))

    ! Gas phase H2SO4 concentration:
    CALL nf_check(nf_def_dim(i_nc_id,h2so4_name,h2so4_n,i_h2so4_dim_id))

    ! Ion pair production rate:
    CALL nf_check(nf_def_dim(i_nc_id,ipr_name,ipr_n,i_ipr_dim_id))

    ! H2SO4 condensation sink:
    CALL nf_check(nf_def_dim(i_nc_id,cs_name,cs_n,i_cs_dim_id))

    !
    ! Define the variables holding the model input / output:
    !

    ! Time:

    CALL nf_check(nf_def_var(i_nc_id,time_name,NF_DOUBLE,1,i_time_dim_id,i_time_id))
    i_length = len(TRIM(time_units))
    CALL nf_check(nf_put_att_text(i_nc_id,i_time_id,'units',i_length,time_units))

    ! Temperature:

    CALL nf_check(nf_def_var(i_nc_id,temp_name,NF_DOUBLE,1,i_temp_dim_id,i_temp_id))

    i_length = len(TRIM(temp_units))
    CALL nf_check(nf_put_att_text(i_nc_id,i_temp_id,'units',i_length,temp_units))

    i_length = len(TRIM(temp_long_name))
    CALL nf_check(nf_put_att_text(i_nc_id,i_temp_id,'long_name',i_length,temp_long_name))

    ! Relative humidity:

    CALL nf_check(nf_def_var(i_nc_id,rh_name,NF_DOUBLE,1,i_rh_dim_id,i_rh_id))

    i_length = len(TRIM(rh_units))
    CALL nf_check(nf_put_att_text(i_nc_id,i_rh_id,'units',i_length,rh_units))

    i_length = len(TRIM(rh_long_name))
    CALL nf_check(nf_put_att_text(i_nc_id,i_rh_id,'long_name',i_length,rh_long_name))

    ! H2SO4 concentration:

    CALL nf_check(nf_def_var(i_nc_id,h2so4_name,NF_DOUBLE,1,i_h2so4_dim_id,i_h2so4_id))

    i_length = len(TRIM(h2so4_units))
    CALL nf_check(nf_put_att_text(i_nc_id,i_h2so4_id,'units',i_length,h2so4_units))

    i_length = len(TRIM(h2so4_long_name))
    CALL nf_check(nf_put_att_text(i_nc_id,i_h2so4_id,'long_name',i_length,h2so4_long_name))

    ! Ionization:

    CALL nf_check(nf_def_var(i_nc_id,ipr_name,NF_DOUBLE,1,i_ipr_dim_id,i_ipr_id))

    i_length = len(TRIM(ipr_units))
    CALL nf_check(nf_put_att_text(i_nc_id,i_ipr_id,'units',i_length,ipr_units))

    i_length = len(TRIM(ipr_long_name))
    CALL nf_check(nf_put_att_text(i_nc_id,i_ipr_id,'long_name',i_length,ipr_long_name))

    ! H2SO4 condensation sink:

    CALL nf_check(nf_def_var(i_nc_id,cs_name,NF_DOUBLE,1,i_cs_dim_id,i_cs_id))

    i_length = len(TRIM(cs_units))
    CALL nf_check(nf_put_att_text(i_nc_id,i_cs_id,'units',i_length,cs_units))

    i_length = len(TRIM(cs_long_name))
    CALL nf_check(nf_put_att_text(i_nc_id,i_cs_id,'long_name',i_length,cs_long_name))

    ! Volume:

    i_vector(1) = i_temp_dim_id
    i_vector(2) = i_rh_dim_id
    i_vector(3) = i_h2so4_dim_id
    i_vector(4) = i_ipr_dim_id
    i_vector(5) = i_cs_dim_id
    i_vector(6) = i_time_dim_id

    CALL nf_check(nf_def_var(i_nc_id,vol_name,NF_FLOAT,6,i_vector,i_vol_id))

    i_length = len(TRIM(vol_units))
    CALL nf_check(nf_put_att_text(i_nc_id,i_vol_id,'units',i_length,vol_units))

    i_length = len(TRIM(vol_long_name))
    CALL nf_check(nf_put_att_text(i_nc_id,i_vol_id,'long_name',i_length,vol_long_name))

    ! Switch from ''define mode'' into ''data mode'' :

    CALL nf_check(NF_ENDDEF(i_nc_id))

    !
    ! Fill the archive with the ambient conditions:
    !

    ! Temperature:

    CALL nf_check(nf_put_var_double(i_nc_id,i_temp_id,temp_grid))

    ! RH:

    CALL nf_check(nf_put_var_double(i_nc_id,i_rh_id,rh_grid))

    ! H2SO4 concentration:

    CALL nf_check(nf_put_var_double(i_nc_id,i_h2so4_id,h2so4_grid))

    ! Ionization:

    CALL nf_check(nf_put_var_double(i_nc_id,i_ipr_id,ipr_grid))

    ! condensation sink:

    CALL nf_check(nf_put_var_double(i_nc_id,i_cs_id,cs_grid))

    CALL nf_check(nf_close(i_nc_id))

  END SUBROUTINE create_nucl_statistics_file

  !=============================================================================

  SUBROUTINE sample_nucl_statistics(kproma,kbdim,klev,ptemp,psatrat,ph2so4,pipr,pcs,paclc,pvol)

    ! *sample_nucl_statistics* identifies the inteval on the temperature, RH,
    ! H2SO4 concentration, H2SO4 condensation sink, and ion pair production rate
    ! grid which the current cloud-free volume belongs to, and increments the
    ! corresponding accumulated volume by the volume at the current location.
    
    INTEGER, INTENT(in)  :: kproma,kbdim,klev

    REAL(dp), INTENT(in) :: ptemp(kbdim,klev), & ! Temperature (K)
         psatrat(kbdim,klev), & ! Saturation ratio (1)
         ph2so4(kbdim,klev), & ! H2SO4 concentration (cm-3)
         pipr(kbdim,klev), & ! Ion pair production rate (cm-3 s-1)
         pcs(kbdim,klev), & ! H2SO4 condensation sink (s-1)
         paclc(kbdim,klev), & ! cloud fraction in a grid box
         pvol(kbdim,klev) ! Grid box volume (m3)

    !
    ! Local variables:
    !

    INTEGER :: ik,jk,jl
    INTEGER :: i_temp_0,i_rh_0,i_h2so4_0,i_ipr_0,i_cs_0 ! Array indices
    INTEGER :: i_temp_1,i_rh_1,i_h2so4_1,i_ipr_1,i_cs_1 ! Array indices

    REAL :: zrh(kbdim,klev)

    zrh = 100.0_dp*psatrat ! RH (%)

    DO jk = 1, kproma
      DO jl = 1, klev

        !
        ! Identify the intervals in the grids arrays where the given values are located:
        !

        ! Temperature

        i_temp_0 = 1
        i_temp_1 = temp_n

        IF (ptemp(jk,jl).gt.temp_grid(temp_n)) THEN

          i_temp_0 = temp_n
          i_temp_1 = 1

        ELSE

          DO WHILE (i_temp_1-i_temp_0.gt.1)

            ik = (i_temp_1+i_temp_0)/2

            IF (temp_grid(ik).gt.ptemp(jk,jl)) THEN
              i_temp_1 = ik
            ELSE
              i_temp_0 = ik
            ENDIF

          ENDDO

        ENDIF

        ! RH

        i_rh_0 = 1
        i_rh_1 = rh_n

        IF (zrh(jk,jl).gt.rh_grid(rh_n)) THEN

          i_rh_0 = rh_n
          i_rh_1 = 1

        ELSE

          DO WHILE (i_rh_1-i_rh_0.gt.1)

            ik = (i_rh_1+i_rh_0)/2

            IF (rh_grid(ik).gt.zrh(jk,jl)) THEN
              i_rh_1 = ik
            ELSE
              i_rh_0 = ik
            ENDIF

          ENDDO

        ENDIF

        ! H2SO4(g) concentration

        i_h2so4_0 = 1
        i_h2so4_1 = h2so4_n

        IF (ph2so4(jk,jl).gt.h2so4_grid(h2so4_n)) THEN

          i_h2so4_0 = h2so4_n
          i_h2so4_1 = 1

        ELSE

          DO WHILE (i_h2so4_1-i_h2so4_0.gt.1)

            ik = (i_h2so4_1+i_h2so4_0)/2

            IF (h2so4_grid(ik).gt.ph2so4(jk,jl)) THEN
              i_h2so4_1 = ik
            ELSE
              i_h2so4_0 = ik
            ENDIF

          ENDDO

        ENDIF

        ! Ion pair production rate

        i_ipr_0 = 1
        i_ipr_1 = ipr_n

        IF (pipr(jk,jl).gt.ipr_grid(ipr_n)) THEN

          i_ipr_0 = ipr_n
          i_ipr_1 = 1

        ELSE

          DO WHILE (i_ipr_1-i_ipr_0.gt.1)

            ik = (i_ipr_1+i_ipr_0)/2

            IF (ipr_grid(ik).gt.pipr(jk,jl)) THEN
              i_ipr_1 = ik
            ELSE
              i_ipr_0 = ik
            ENDIF

          ENDDO

        ENDIF

        ! Condensation sink

        i_cs_0 = 1
        i_cs_1 = cs_n

        IF (pcs(jk,jl).gt.cs_grid(cs_n)) THEN

          i_cs_0 = cs_n
          i_cs_1 = 1

        ELSE

          DO WHILE (i_cs_1-i_cs_0.gt.1)

            ik = (i_cs_1+i_cs_0)/2

            IF (cs_grid(ik).gt.pcs(jk,jl)) THEN
              i_cs_1 = ik
            ELSE
              i_cs_0 = ik
            ENDIF

          ENDDO

        ENDIF

        ! Add up the volume of the cloud-free portion of the model gridbox in the
        ! corresponding location of the space spanned by the temperature, RH,
        ! H2SO4 concentration, ionization rate, and condensation sink:

        vol(i_temp_0,i_rh_0,i_h2so4_0,i_ipr_0,i_cs_0) = &
             vol(i_temp_0,i_rh_0,i_h2so4_0,i_ipr_0,i_cs_0) + (1.0_dp-paclc(jk,jl))*pvol(jk,jl)

      ENDDO
    ENDDO

  END SUBROUTINE sample_nucl_statistics

  !=============================================================================

  SUBROUTINE save_nucl_statistics()

    ! *save_nucl_statistics* opens a netCDF archive and adds the array
    ! holding the volume distribution as function of temperature, RH, H2SO4
    ! concentration, H2SO4 condensation sink, and ion pair production rate.
    !
    ! Local variables:
    !

    INTEGER :: i_nc_id,i_vol_id,i_time_dim_id,i_time_id

    TYPE (time_native) :: date
    INTEGER :: iyear,imonth,iday,ihour,iminute,isecond

    INTEGER :: i_time
    INTEGER :: istart(6),icount(6)

    REAL(dp) :: ztime

    ! Generate the output file name:

    file = TRIM(standard_output_file)//'_nucstat.nc'

    ! Generate the current time:

    CALL TC_convert(current_date,date)
    CALL TC_get(date,iyear,imonth,iday,ihour,iminute,isecond)

    ztime = 10000.0_dp*REAL(iyear,dp) &
         + 100.0_dp*REAL(imonth,dp) &
         + REAL(iday,dp) &
         + REAL(ihour,dp)/24.0_dp &
         + REAL(iminute,dp)/24.0_dp/60.0_dp &
         + REAL(isecond,dp)/24.0_dp/60.0_dp/60.0_dp ! day as %Y%m%d.%f

    ! Create the netCDF output file, if it does not yet exist:
    CALL create_nucl_statistics_file()

    ! Open the netCDF archive for reading/writing:
    CALL nf_check(nf_open(file,NF_WRITE,i_nc_id))

    ! Get the ID of the time dimension:
    CALL nf_check(nf_inq_unlimdim(i_nc_id,i_time_dim_id))

    ! Get the current length of the time dimension:
    CALL nf_check(nf_inq_dimlen(i_nc_id,i_time_dim_id,i_time))

    ! Get the ID of the time variable:
    CALL nf_check(nf_inq_varid(i_nc_id,time_name,i_time_id))

    ! Get the ID of the nucleation volume variable:
    CALL nf_check(nf_inq_varid(i_nc_id,vol_name,i_vol_id))

    i_time = i_time + 1

    ! Save the data:

    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = 1
    istart(5) = 1
    istart(6) = i_time

    icount(1) = temp_n
    icount(2) = rh_n
    icount(3) = h2so4_n
    icount(4) = ipr_n
    icount(5) = cs_n
    icount(6) = 1

    CALL nf_check(nf_put_vara_double(i_nc_id,i_vol_id,istart,icount,vol))

    ! Save the time:

    CALL nf_check(nf_put_var1_double(i_nc_id,i_time_id,i_time,ztime))

    ! Close the netCDF archive:
    CALL nf_check(nf_close(i_nc_id))

  END SUBROUTINE save_nucl_statistics

  !=============================================================================

  SUBROUTINE ham_nucl_diag_cleanup()

    ! *ham_nucl_diag_cleanup* gathers the arrays holding the volume
    ! distribution from the different processes, sums them up, saves the
    ! resulting array in a netCDF file, and deallocates all corresponding
    ! arrays.
    !
    ! Local variables:
    !

    INTEGER :: ji

    ! Making recvbuf allocatable will not work when using the NAG compiler;
    ! this compiler will complain that the allocatable array recvbuf is not
    ! allocated when passing it to p_gather. p_gather does, however, not need an
    ! allocated recvbuf array unless the CALL is issued in the receiving
    ! process. When using the NAG compiler, more memory will be used.

#if defined (NAG)
    REAL(dp) :: recvbuf(temp_n,rh_n,h2so4_n,ipr_n,cs_n,p_nprocs)
    recvbuf = 0.0_dp
#else
    REAL(dp), ALLOCATABLE :: recvbuf(:,:,:,:,:,:)
#endif

    IF (p_parallel) THEN ! In parallel mode

      IF (p_parallel_io) THEN ! We are on the I/O processor

#if defined (NAG)

#else
        ALLOCATE(recvbuf(temp_n,rh_n,h2so4_n,ipr_n,cs_n,p_nprocs))
#endif

      END IF

      CALL p_gather(vol,recvbuf,p_io)

      IF (p_parallel_io) THEN ! We are on the I/O processor

        vol(:,:,:,:,:) = recvbuf(:,:,:,:,:,1)

        DO ji = 2, p_nprocs
          vol(:,:,:,:,:) = vol(:,:,:,:,:) + recvbuf(:,:,:,:,:,ji)
        ENDDO

        CALL save_nucl_statistics()

#if defined (NAG)

#else
        DEALLOCATE(recvbuf)
#endif

      END IF

    ELSE ! In serial mode

      CALL save_nucl_statistics()

    END IF

    DEALLOCATE(temp_grid)
    DEALLOCATE(rh_grid)
    DEALLOCATE(h2so4_grid)
    DEALLOCATE(ipr_grid)
    DEALLOCATE(cs_grid)

    DEALLOCATE(vol)

  END SUBROUTINE ham_nucl_diag_cleanup

  !=============================================================================

  SUBROUTINE hyper_geometric_nodes_1d(p_x_0,p_x_1,k_x,p_t,pkappa)

    ! *hyper_geometric_nodes_1d* calculates n hypergeometric nodes t covering
    ! the interval [p_x_0,p_x_1]. p_x_0 and p_x_1 must be both > 0 or both < 0.
    ! For pkappa > 1, the nodes will be denser in the upper portion of the
    ! interval compared with a geometric node progression, and less dense in the
    ! lower portion. The opposite is the case for pkappa < 1. For pkappa = 1,
    ! one obtains a geometric node progression.

    INTEGER,  INTENT(in) :: k_x
    REAL(dp), INTENT(in) :: p_x_0,p_x_1
    REAL(dp), INTENT(in) :: pkappa

    REAL(dp), INTENT(out) :: p_t(k_x)

    !
    ! Local variables:
    !

    REAL(dp) :: zfactor,zexponent

    CHARACTER(512) :: routine,error_message

    ! Loop variables:

    INTEGER :: ji

    IF (p_x_1*p_x_0.le.0.0_dp) THEN

      routine = 'hyper_geometric_nodes_1d'
      error_message = 'Bounds of the interval on which the (hyper-) geometric nodes are requested contains zero.'
      CALL message(routine,error_message)
      error_message = 'Run terminated.'
      CALL finish(routine,error_message)

    END IF

    zexponent = 0.0_dp

    DO ji = 2, k_x
      zexponent = zexponent + 1.0_dp/pkappa &
           + (pkappa - 1.0_dp/pkappa)*REAL(ji-k_x,dp)/REAL(2-k_x,dp)
    ENDDO

    zfactor = (p_x_1/p_x_0)**(1.0_dp/zexponent)

    p_t(1) = p_x_0

    DO ji = 2, k_x - 1
      zexponent =  + 1.0_dp/pkappa &
           + (pkappa - 1.0_dp/pkappa)*REAL(ji-k_x,dp)/REAL(2-k_x,dp)
      p_t(ji) = p_t(ji-1)*zfactor**zexponent
    ENDDO

    p_t(k_x) = p_x_1

  END SUBROUTINE hyper_geometric_nodes_1d

  !=============================================================================

  SUBROUTINE equidistant_nodes(p_x_0,p_x_1,k_x,p_t)

    ! *equidistant_nodes* calculates n nodes t covering the interval
    ! [p_x_0,p_x_1] with equal spacing.

    INTEGER,  INTENT(in) :: k_x
    REAL(dp), INTENT(in) :: p_x_0,p_x_1

    REAL(dp), INTENT(out) :: p_t(k_x)

    ! Loop variables:

    INTEGER  :: ji

    REAL(dp) :: zdelta

    zdelta = (p_x_1-p_x_0)/REAL(k_x-1,dp)

    p_t(1) = p_x_0

    DO ji = 2, k_x - 1
      p_t(ji) = p_t(ji-1) + zdelta
    ENDDO

    p_t(k_x) = p_x_1

  END SUBROUTINE equidistant_nodes

  !=============================================================================

  SUBROUTINE arched_nodes(p_x_0,p_x_1,k_x,p_t,pkappa)

    ! *arched_nodes* calculates n nodes p_t covering the interval [p_x_0,p_x_1].
    ! pkappa must be in ]-1,1[. For pkappa > 1 one obtains increasing node
    ! distances on the lower half and decreasing node distances on the upper
    ! half of the interval [p_x_0,p_x_1]. For pkappa < 0 the opposite is the case.
    ! Hence in the first case the node distance in the center of the interval
    ! will be larger compared to the node distance at the interval edges, and
    ! vice versa in the latter case.
    !
    ! Status: tested and functional

    INTEGER,  INTENT(in) :: k_x
    REAL(dp), INTENT(in) :: p_x_0,p_x_1
    REAL(dp), INTENT(in) :: pkappa

    REAL(dp), INTENT(out)  :: p_t(k_x)

    !
    ! Local variables:
    !

    REAL(dp) :: zdx
    REAL(dp) :: zl_array(k_x)
    REAL(dp) :: zr_array(k_x)
    REAL(dp) :: zx_l,zx_r
    REAL(dp) :: zl_pkappa,zr_pkappa

    INTEGER :: il_len,ir_len

    CHARACTER(512) :: routine,error_message

    ! Loop variables:

    INTEGER :: ji

    IF (k_x.lt.3) then

      routine = 'arched_nodes'
      error_message = 'The array p_t contains less than 3 elements, which is too few.'
      CALL message(routine,error_message)
      error_message = 'Run terminated.'
      CALL finish(routine,error_message)

    END IF

    zdx = (p_x_1-p_x_0)/REAL(k_x-1,dp)

    zl_pkappa =  pkappa
    zr_pkappa = -pkappa

    IF (mod(k_x,2).eq.0) then

      il_len = k_x/2
      ir_len = k_x/2
      zx_l   = 0.5_dp*(p_x_0+p_x_1-(1.0_dp+pkappa)*zdx)
      zx_r   = 0.5_dp*(p_x_0+p_x_1+(1.0_dp+pkappa)*zdx)

      CALL slanted_nodes_a(p_x_0,zx_l,il_len,zl_array,zl_pkappa)
      CALL slanted_nodes_a(zx_r,p_x_1,ir_len,zr_array,zr_pkappa)

      DO ji = 1, il_len
        p_t(ji) = zl_array(ji)
      END DO

      DO ji = 1, ir_len
        p_t(il_len+ji) = zr_array(ji)
      END DO

    ELSE

      il_len = (k_x+1)/2
      ir_len = (k_x+1)/2
      zx_l   = 0.5_dp*(p_x_0+p_x_1)
      zx_r   = 0.5_dp*(p_x_0+p_x_1)

      CALL slanted_nodes_a(p_x_0,zx_l,il_len,zl_array,zl_pkappa)
      CALL slanted_nodes_a(zx_r,p_x_1,ir_len,zr_array,zr_pkappa)

      DO ji = 1, il_len
        p_t(ji) = zl_array(ji)
      END DO

      DO ji = 1, ir_len - 1
        p_t(il_len+ji) = zr_array(ji+1)
      END DO

    END IF

  END SUBROUTINE arched_nodes

  !=============================================================================

  SUBROUTINE slanted_nodes_a(p_x_0,p_x_1,k_x,p_t,pkappa)

    ! *slanted_nodes_a* calculates n nodes p_t covering the interval [p_x_0,p_x_1].
    ! pkappa must be in ]-1,1[. For pkappa = 0 one obtains equidistant nodes.
    ! For pkappa > 0 one obtains increasing and for pkappa < 0 decreasing node
    ! distances on the interval [p_x_0,p_x_1].
    !
    ! Status: tested and functional

    INTEGER,  INTENT(in) :: k_x
    REAL(dp), INTENT(in) :: p_x_0,p_x_1
    REAL(dp), INTENT(in) :: pkappa

    REAL(dp), INTENT(out) :: p_t(k_x)

    !
    ! Local variables:
    !

    REAL(dp) :: zdx

    ! Loop variables:

    INTEGER :: ji

    zdx = (p_x_1-p_x_0)/REAL(k_x-1,dp)

    p_t(1) = p_x_0

    DO ji = 2, k_x - 1
      p_t(ji) = p_t(ji-1) &
           + (1.0_dp + pkappa*REAL(2*(ji-1)-k_x,dp)/REAL(k_x-2,dp))*zdx
    END DO

    p_t(k_x) = p_x_1

  END SUBROUTINE slanted_nodes_a
  
END MODULE mo_ham_m7_nucl_diag 
