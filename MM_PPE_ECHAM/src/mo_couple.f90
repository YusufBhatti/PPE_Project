!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_couple
!
! !DESCRIPTION:
!
!   Contains additional code used for coupling via OASIS3 or OASIS3-MCT.
!   Note: coupling initialisation and retrieval of the local communicator for 
!   model internal parallelisation are treated in mo_mpi module.
!  
!   coupling strategies : 
!     compile flag cpl_mpiom    if used with MPIOM;
!     namelist flag lcouple_co2 if used with MPIOM, HAMOCC and interactive CO2
!
!   exchanged fields (see SBR collect):
!     over water:   pawsol  downwelling solar radiation (used for solar 
!                           penetration)
!                   pawsta  ustar**3 (C-grid HOPE only)
!                   pawhea  net heat flux
!                   pawfre  liquid freshwater flux
!                   pawust  zonal wind stress
!                   pawvst  meridional wind stress
!
!     over sea ice: paicon  conductive heat flux
!                   paiqre  residual heat flux (used for surface melt)
!                   paifre  solid freshwater flux
!                   paiust  zonal wind stress
!                   paivst  meridional wind stress
!
!     cpl_mpiom   defined: 
!                   aicon,aiqre not multiplied with sea ice concentration
!                   awhea       not multiplied with open water fraction
!                   awust,awvst not multiplied with open water fraction
!     always     multiplied with sea ice concentration: aifre
!     always not  "       "    "   "       "          : awsol,awsta,aiust,
!                                                       aivst
!                  weighted average of both           : awfre
!
! !USES:

  USE mo_kind,            ONLY: wp ! working precision (echam6)
  USE mo_control,         ONLY: nlon,ngl, lfractional_mask, lcouple_co2
  USE mo_time_control,    ONLY: delta_time,get_time_step,lresume,l_putocean, &
                                lbreak, stop_date, next_date, &
                                write_date, l_putrerun
  USE mo_memory_g3b,      ONLY: awsol, &
                                awhea,awust,awvst,awsta, &
                                aicon,aiqre,aifre,aiust,aivst, &
                                seaice,tsw,siced,sni,alake,slf,slm, &
                                ocu,ocv, apmebco, rain, resi, pmeb
  USE mo_hydrology,       ONLY: awfre, disch
  USE mo_io_units,        ONLY: nerr
  USE mo_co2,             ONLY: co2trans, co2ocean, co2atmos, co2flux_cpl
  USE mo_physical_constants,       ONLY: rhoh2o
  USE mo_mpi,             ONLY: p_bcast, p_io, p_parallel_io, p_parallel, p_all_comm, &
                                p_pe
  USE mo_decomposition,   ONLY: gd => global_decomposition, ld => local_decomposition

  USE mo_tr_gather,       ONLY: gather_field
  USE mo_tr_scatter,      ONLY: scatter_field     
  USE mo_exception,       ONLY: message_text, message, em_warn, em_error, finish

#if defined (__prism)
  USE mo_couple_wrap,     ONLY : oasis_init_comp, oasis_get_localcomm, &
                                 oasis_def_partition, oasis_def_var, oasis_enddef, &
                                 oasis_start_grids_writing, oasis_write_grid, &
                                 oasis_write_corner, oasis_write_mask, oasis_write_area, &
                                 oasis_terminate_grids_writing, oasis_get, oasis_put, &
                                 oasis_terminate, oasis_abort, &
                                 oasis_ok, oasis_Recvd, oasis_FromRest, &
                                 oasis_RecvOut, oasis_FromRestOut, oasis_Sent, oasis_ToRest, &
                                 oasis_SentOut, oasis_ToRestOut, oasis_Input, oasis_LocTrans, &
                                 oasis_In, oasis_Out, oasis_Output, oasis_Real, &
                                 ip_realwp_p, &
                                 clim_strategy, clim_serial, clim_length, clim_offset, &
                                 clim_orange, clim_segments, clim_sizex, clim_sizey, clim_ldx, &
                                 namsrcfld, namdstfld, namfldseq, namruntim, namrstfil, &
                                 namflddti, nnamcpl, namfldlag
#ifdef __oa3mct
  USE mo_couple_wrap,     ONLY : prism_nmodels, prism_modnam
#else
  USE mo_couple_wrap,     ONLY : ig_inidate, lg_ncdfrst, oasis_put_restart
#endif

  USE mo_gaussgrid,       ONLY: philon, philat, gridarea, gl_budw, gridarea
  USE mo_hydrology,       ONLY: diag_water_budget
  USE mo_filename,        ONLY: find_next_free_unit
  USE mo_mpi,             ONLY: prism_model_number, prism_model_name
#endif

  IMPLICIT NONE    

  PRIVATE

#if defined (__prism)
  INTEGER, PARAMETER :: pwp = ip_realwp_p   ! working precision (psmile)

  INTEGER, ALLOCATABLE :: paral(:) ! parallel strategy
  INTEGER :: nsegments             ! no. of segments
  INTEGER :: parsize               ! size of partition description array
  INTEGER :: nbtotproc             ! total no. of procs
  INTEGER :: nbcplproc             ! no. of procs involved in data exchange
  INTEGER :: info                  ! info codes
  INTEGER :: ierror                ! error codes
  INTEGER :: part_id               ! partition ID
  INTEGER :: var_nodims(2)             ! 1: rank of the coupling field array
                                       ! 2: number of bundles (1 for oasis3)
  INTEGER, ALLOCATABLE :: var_shape(:) ! min and max indices for each dimension
                                       ! of the coupling field array
  INTEGER :: debug_unit

#if defined __cpl_mpiom
  INTEGER :: nflda2o     ! no. of fields passed from atmosphere to ocean
  INTEGER :: nfldo2a     ! no. of fields passed from ocean to atmosphere
  PUBLIC  :: nflda2o, nfldo2a
#endif

  CHARACTER (len=8), ALLOCATABLE :: clstr8a2o(:) ! Port names of exchange fields sent
  CHARACTER (len=8), ALLOCATABLE :: clstr8o2a(:) ! Port names of exchange fields received

  INTEGER*4, ALLOCATABLE :: portout_id(:) ! Port IDs of exchange fields sent
  INTEGER*4, ALLOCATABLE :: portin_id(:)  ! Port IDs of exchange fields received

  REAL(pwp) :: couple_a2o_time = 0.0_pwp  ! time passed since last couple_put_a2o
#else
  REAL(wp)  :: couple_a2o_time = 0.0_wp   ! time passed since last couple_put_a2o
#endif

#if defined (__prism)
  REAL(pwp), POINTER :: exfld(:,:)       ! buffer for receiving exchange fields
#else
  REAL(wp), POINTER  :: exfld(:,:)       ! buffer for receiving exchange fields
#endif
  REAL(wp), POINTER  :: tmp_loco2a(:,:)  ! temporary buffer for local exchange fields

  INTEGER :: ngpblks                     ! number of blocks
  INTEGER :: nbdim, nproma               ! size of blocks (nproma)
  INTEGER :: nbdimz                      ! size of last block (npromz)

  INTEGER :: nmseq                       ! Maximum sequential number of exchange fields

  INTEGER :: run_date_secs  = 0          ! time (sec) passed since start of run

  INTEGER :: couple_resum = 0            ! model step when coupled run is resumed

  REAL(wp), POINTER :: gl_slf  (:,:)     ! global fractional land sea mask
  REAL(wp), POINTER :: gl_slm  (:,:)     ! global [1,0] land sea mask
  REAL(wp), POINTER :: gl_alake(:,:)     ! global fractional lake mask

  CHARACTER(LEN=20) :: decomp_strategy  ! decomposition strategy, possible values:
                                        !   serial (= none)
                                        !   orange

  LOGICAL :: lcouple_parallel           ! true if all PEs communicate with oasis
  LOGICAL :: p_parallel_oasis           ! PEs communicating with oasis
  LOGICAL :: ldebugcpl                  ! activate coupling debug output

  LOGICAL, ALLOCATABLE :: locean_cell(:,:)


  PUBLIC :: couple_init, couple_end, couple_get_o2a, couple_put_a2o, couple_calendar, &
            couple_a2o_time, lcouple_parallel, ldebugcpl
!
!
! !DESCRIPTION:
!
! - contains all code related to coupling with oasis3 or oasis3-mct
! - Note: this interface to the psmile library uses variables from
!   the psmile modules mod_comprism_proto (coupler oasis3) or 
!   mod_oasis_namcouple (coupler oasis3-mct)
!   
! !REVISION HISTORY:
! 03.05.14 Stephanie Legutke, MPI/M&D
!       - created
! 03.05.14 Veronika Gayler,   MPI/M&D
!       -
! 03.07.08 Luis Kornblueh,    MPI
!       - minor modifications
! 03.08.04 Johann Jungclaus MPI
!       - include call for water budget correction
! 10.07.15 Rene Redler, MPI
!       - cleaning up
!       - elimination of uncessary gather call during initialisation and put_o2a
!       - elimination of global 2d array for put_o2a
!       - elimination of redundant use of arrays
!       - pfield changed to pointer to avoid copying
! 10.08.25 Rene Redler, MPI
!       - replace writing to kout and kerr by call message and call message_text
!       - removed kout and kerr
!       - renamed __synout to __cpl_verbose
!       - additional output into fort.99 with __cpl_debug
! 12.04.03 Irina Fast, DKRZ
!       - allow for usage of oasis3-mct
!       - allow for parallel communication with oasis3
!       - replace CPP keys __cpl_verbose and __cpl_debug with logical switch
!         ldebugcpl
!
! EOP
!-----------------------------------------------------------------------
! 

CONTAINS

!-------------------------------------------------------------------------------
!>
!! Initialize coupling
!! This routine should not be called in uncoupled mode !
!!
!! 03.05.15  S. Legutke - created
  SUBROUTINE couple_init
    INTEGER :: k, is, jfld    ! loop indices
    INTEGER :: local_extent, global_offset, i1, i2, d1, d2, iost
    CHARACTER(LEN=20) :: debug_outfile, coupler

#if defined (__prism)

    IF (ldebugcpl) THEN
      CALL message('','Start couple_init')
      CALL message('','*****************')
      FLUSH(nerr)
    ENDIF

    CALL MPI_Comm_Size(p_all_comm, nbtotproc, ierror)

#ifdef __oa3mct
    coupler = 'OASIS3-MCT'
    lcouple_parallel = .TRUE.
#else
    coupler = 'OASIS3'
#endif

    IF (lcouple_parallel) THEN
      nbcplproc = nbtotproc
    ELSE
      nbcplproc = 1
    ENDIF

    nmseq = MAXVAL (namfldseq)

    p_parallel_oasis = (lcouple_parallel .OR. p_pe == p_io)
    IF (lcouple_parallel) THEN
      decomp_strategy = 'orange'
    ELSE
      decomp_strategy = 'serial'
    ENDIF

    IF (ldebugcpl) THEN
      debug_unit = find_next_free_unit (51,100)
      IF (p_parallel_io) THEN
        debug_outfile = 'echam_cpl_debug'
      ELSE
        WRITE(debug_outfile,'(A,"_",I4.4)') 'echam_cpl_debug',p_pe
      ENDIF
      IF (p_parallel_oasis) THEN
        iost = 0
        OPEN(debug_unit,FILE = debug_outfile, &
          STATUS='unknown',FORM ='formatted',IOSTAT = iost)
        IF (iost .NE. 0) &
          CALL finish ('mo_couple: couple_init','could not open output file '//debug_outfile)
      ENDIF
    ENDIF

    ! get control variables for selecting features on global arrays
    ! gl_slf, gl_slm and gl_alake are not needed in most cases and
    ! could be shifted into the local routines.

    IF (p_parallel_io) THEN
      ALLOCATE (gl_slf  (nlon,ngl))
      ALLOCATE (gl_slm  (nlon,ngl))
      ALLOCATE (gl_alake(nlon,ngl))
    ELSE
      gl_slf   => NULL()
      gl_slm   => NULL()
      gl_alake => NULL()
    END IF

    ngpblks = ld% ngpblks ! number of rows
    nbdim   = ld% nproma
    nbdimz  = ld% npromz


    ! Handling exfld and tmp_loco2a this way seems to give better
    ! performance compared to doing a local declaration and allocate
    ! deallocate within the subroutine.
    ALLOCATE(tmp_loco2a(nbdim,ngpblks))

    IF (lcouple_parallel) THEN
      !d1 = nbdim
      !d2 = ngpblks
      d1 = ld% nglon
      d2 = ld% nglat
    ELSE
      d1 = nlon
      d2 = ngl
    ENDIF

    IF (p_parallel_oasis) THEN
      ALLOCATE (exfld(d1,d2))
    ELSE
      exfld    => NULL()
    ENDIF

    ! Ocean cells
    ALLOCATE (locean_cell(nbdim, ngpblks))
    locean_cell = ( slf < 1.0_wp .AND. alake == 0.0_wp )
    IF (nbdimz < nbdim) locean_cell(nbdimz+1:,ngpblks) = .FALSE.

    CALL gather_field(gl_slf, slf)
    CALL gather_field(gl_slm, slm)
    CALL gather_field(gl_alake, alake)

    !
    !-- Get model step when coupled run is started.
    !  
    couple_resum = get_time_step()

    WRITE (message_text,'(a,i12)') 'Coupled run is started at model step ', couple_resum

    IF (ldebugcpl) THEN
      !
      !-- Check land/lake/sea masks. 
      !  
      ! CALL check_alake ! FIXME: check the routine
      !
      !-- Print land/lake/sea masks.
      !  
      ! CALL print_masks ! FIXME: check the routine
    ENDIF

    !
    !-- Write grids file for oasis
    !
    CALL grids_writing

    !
    !-- Set exchange-field locators.
    !
#if defined __cpl_mpiom
    IF (lcouple_co2) THEN
      nflda2o = 13 ! no. of fields passed from atmosphere to ocean
      nfldo2a =  8 ! no. of fields passed from ocean to atmosphere
    ELSE
      nflda2o = 11 ! no. of fields passed from atmosphere to ocean
      nfldo2a =  6 ! no. of fields passed from ocean to atmosphere
    END IF
    !
    ! Set locator for fields passed from the ocean to the atmosphere.
    !
    ALLOCATE(clstr8o2a(nfldo2a)) 
    ALLOCATE(portin_id(nfldo2a))  ! Port IDs of exchange fields received
    clstr8o2a(1:nfldo2a) = namdstfld(1:nfldo2a)
    !
    ! Set locator for fields passed from the atmosphere to the ocean.
    !
    ALLOCATE(clstr8a2o(nflda2o))
    ALLOCATE(portout_id(nflda2o)) ! Port IDs of exchange fields sent
    clstr8a2o(1:nflda2o) = namsrcfld(nfldo2a+1:nfldo2a+nflda2o)
#endif

    !> Describe decomposition strategy for field exchange
    SELECT CASE (TRIM(decomp_strategy))
      CASE('serial') ! none
        nsegments = 1
        parsize   = 3
        ALLOCATE(paral(parsize))
        paral ( clim_strategy ) = clim_serial
        paral ( clim_length   ) = nlon*ngl
        paral ( clim_offset   ) = 0
        var_nodims(1) = 1  ! rank of coupling field array
        var_nodims(2) = 1  ! always 1 for OASIS3
        ALLOCATE(var_shape(2*var_nodims(1)))
        var_shape(1)  = 1  ! min and max index for each dimension of the field array
        var_shape(2)  = paral (clim_length)

      CASE('orange')
        nsegments = size(ld%glat)
        parsize   = 2 + 2 * nsegments
        ALLOCATE(paral(parsize))
        paral ( clim_strategy ) = clim_orange
        paral ( clim_segments ) = nsegments
        DO is=1,nsegments
          i1 = clim_segments + 2*is - 1
          i2 = i1 + 1
          local_extent  = ld%nglon
          global_offset = (ld%glat(is) - 1) * nlon + ld%glon(is)
          paral(i1) = global_offset
          paral(i2) = local_extent
        ENDDO 
        var_nodims(1) = 2  ! rank of coupling field array
        var_nodims(2) = 1  ! always 1 for OASIS3
        ALLOCATE(var_shape(2*var_nodims(1)))
        var_shape(1)  = 1  ! min and max index for each dimension of the field array
        var_shape(2)  = ld%nglon
        var_shape(3)  = 1
        var_shape(4)  = ld%nglat
      CASE DEFAULT
        CALL oasis_abort(0,'couple_init', &
            'unhandled decomposition startegy for model '//prism_model_name)
    END SELECT

    !
    !-- Associate a part_id to the decomposition of the process.
    !  
    IF (p_parallel_oasis) THEN
       ierror   = oasis_ok
       CALL oasis_def_partition(part_id, paral, ierror)
       IF (ierror /= oasis_ok) &
          CALL oasis_abort (0,'couple_init', &
              'oasis_def_partition failed for model '//prism_model_name)

    !     
    !-- Print coupling parameters received from oasis3 or oasis3-mct
    !  
    IF (p_parallel_io) CALL print_couplparam

    !
    !-- Define ports of incoming fields.
    !  

       DO jfld = 1,nfldo2a

          ierror = oasis_ok

          CALL oasis_def_var(portin_id(jfld), clstr8o2a(jfld), &
                             part_id, var_nodims, oasis_In,    &
                             var_shape, oasis_Real, ierror )

          IF (ierror /= oasis_ok) THEN
            WRITE (message_text,'(a,a)')  'Problem with import port ',clstr8o2a(jfld)
            CALL message('',message_text, level=em_warn)
            WRITE (message_text,'(a,i4)') 'Port ID returned is : ', portin_id(jfld)
            CALL message('',message_text, level=em_warn)
            WRITE (message_text,'(a,i4)') ' =======   Error code number = ', ierror
            CALL message('',message_text, level=em_warn)
          ELSE
            WRITE (message_text,'(a,a)')   ' couple_init: Import port defined   ', clstr8o2a(jfld)
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)')  ' couple_init: Port ID returned is:  ', portin_id(jfld)
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)')  ' -     with port status             ', oasis_In
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)')  ' -     with port data type          ', oasis_Real
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)')  ' -     with decomposition  paral(1) ', paral(1)
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)') '  -                         paral(2) ', paral(2)
            CALL message('',message_text) 
            WRITE (message_text,'(a,i6)') '  -                         paral(3) ', paral(3)
            CALL message('',message_text)
          ENDIF 

       END DO

    !
    !--   Defined ports of outgoing fields.
    !  

       DO jfld = 1,nflda2o

          ierror = oasis_ok
          CALL oasis_def_var(portout_id(jfld), clstr8a2o(jfld), &
                             part_id, var_nodims, oasis_Out,    &
                             var_shape, oasis_Real, ierror )
          IF (ierror /= oasis_ok) THEN
            WRITE (message_text,'(a,a)')  'Problem with export port ',clstr8a2o(jfld)
            CALL message('',message_text, level=em_warn)
            WRITE (message_text,'(a,i4)') 'Port ID returned is : ', portout_id(jfld)
            CALL message('',message_text, level=em_warn)
            WRITE (message_text,'(a,i4)') ' =======   Error code number = ', ierror
            CALL message('',message_text, level=em_warn)
          ELSE

            WRITE (message_text,'(a,a)')  ' couple_init: Export port defined   ', clstr8a2o(jfld)
            CALL message('',message_text)
            WRITE (message_text,'(a,i4)') ' couple_init: Port ID returned is:  ', portout_id(jfld)
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)') ' -     with port status             ', oasis_Out
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)') ' -     with port data type          ', oasis_Real
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)') ' -     with decomposition  paral(1) ', paral(1)
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)') ' -                         paral(2) ', paral(2)
            CALL message('',message_text)
            WRITE (message_text,'(a,i6)') ' -                         paral(3) ', paral(3)
            CALL message('',message_text)
          ENDIF 

       END DO

       ierror = oasis_ok
       CALL oasis_enddef(ierror)
       IF (ierror /= oasis_ok) THEN
         CALL message('',' Problem with oasis_enddef', level=em_warn)
         WRITE (message_text,'(a,i6)')' =======   Error code number = ', ierror
         CALL message('',message_text)
       ENDIF 

    ENDIF

    !
    !-- Initilize exchange fields.
    !

    CALL ini_a2o

    !
    !-- Check coupling-control parameters
    !   against those received from oasis.
    !  

    CALL chck_par(lresume,nlon,ngl)

    IF (ldebugcpl) THEN
      CALL message('','End of couple_init')
      CALL message('','******************')
      FLUSH(nerr)
    ENDIF

  CONTAINS

    SUBROUTINE print_couplparam
    !     
    !-- Print coupling parameters
    !  
    CALL message ('',' Couple_init:')
    CALL message ('',' -----------------------------------------------------')
    
    WRITE (message_text,'(a,t20,a)') 'Coupler', '= '//TRIM(coupler)
    CALL message('',message_text)

    WRITE (message_text,'(a,t20,a)') 'Model name', '= '//TRIM(prism_model_name)
    CALL message('',message_text)
    
    WRITE (message_text,'(a,t20,a,f12.5)') 'Model time step', '=', delta_time
    CALL message('',message_text)
    
    CALL message ('',' -----------------------------------------------------')
    CALL message ('',' Parameters received from coupler:')

    WRITE (message_text,'(a,i0)')  '- Total no. of processors without I/O servers ', nbtotproc
    CALL message('',message_text)
    
    WRITE (message_text,'(a,t40,i0)')  '- No. of procs for data exchange', nbcplproc
    CALL message('',message_text)
    
    WRITE (message_text,'(a,t40,i0)')  '- Model number', prism_model_number
    CALL message('',message_text)
    
    WRITE (message_text,'(a,t40,i0)')  '- Depomposition ID', part_id
    CALL message('',message_text)

    WRITE (message_text,'(a,t40,i0)') '- Total time of the simulation', namruntim
    CALL message('',message_text)

    WRITE (message_text,'(a,t40,i0)')  '- Number of fields exchanged', nnamcpl
    CALL message('',message_text)
    
#ifndef __oa3mct
    WRITE (message_text,'(a,t40,6i4)') '- Initial date', ig_inidate
    CALL message('',message_text)
#endif
    
    DO jfld = 1, nfldo2a + nflda2o
      CALL message('','')
      IF (jfld <= nfldo2a) THEN 
         WRITE (message_text,'(a,t40,a)') '- Name of exchanged field (in)', TRIM(namdstfld(jfld))
      ELSE
         WRITE (message_text,'(a,t40,a)') '- Name of exchanged field (out)', TRIM(namsrcfld(jfld))
      ENDIF
      CALL message('',message_text)
      
      WRITE (message_text,'(a,t40,i0)') '-   lag of field', namfldlag(jfld)
      CALL message('',message_text)

      WRITE (message_text,'(a,t40,i0)') '-   coupling period of field', namflddti(jfld)
      CALL message('',message_text)
      
      WRITE (message_text,'(a,t40,i0)') '-   sequential index of field', namfldseq(jfld)
      CALL message('',message_text)
      
      WRITE (message_text,'(a,t40,a)') '-   name of restart file', TRIM(namrstfil(jfld))
      CALL message('',message_text)
    ENDDO
#ifndef __oa3mct
      WRITE (message_text,'(a,t40,l4)') '- Restart files in netcdf', lg_ncdfrst
      CALL message('',message_text)
#endif
    
    IF(nmseq /= 1) THEN
      CALL message ('',' The models run sequentially !')
      WRITE (message_text,'(a,i4)') ' No. of sequential exchange fields: ', nmseq
      CALL message('',message_text)
    ELSE
      CALL message('','  All fields have sequential no. 1 !')
      CALL message('','  The models run concurrently!'//new_line('x'))
    ENDIF

    END SUBROUTINE print_couplparam
#endif /* __prism */

  END SUBROUTINE couple_init

!-------------------------------------------------------------------------------
!>
!! Receives ocean data from the coupling system; called in stepon.
!! 03.05.15  S. Legutke - created
  SUBROUTINE couple_get_o2a

#if defined (__prism)

    INTEGER   :: istep ! model step since initialisation
    INTEGER   :: jfld  ! loop index

    IF (ldebugcpl) THEN
      istep = get_time_step()
      CALL message('','Start of couple_get_o2a')
      CALL message('','***********************')
      WRITE (message_text,'(a,i12)') ' model step since init. : ',istep
      CALL message('',message_text)
      WRITE (message_text,'(a,i12)') ' model run step         : ',istep-couple_resum
      CALL message('',message_text)
      WRITE (message_text,'(a,i12)') ' run time (sec) passed  : ',run_date_secs
      CALL message('',message_text)
      FLUSH(nerr)
    ENDIF

    DO jfld = 1,nfldo2a

      info = oasis_ok
      IF ( p_parallel_oasis ) THEN
        IF (ldebugcpl) THEN
          WRITE (message_text,'(a,a)') ' oasis_get field with ', clstr8o2a(jfld) 
          CALL message('',message_text)
          FLUSH(nerr)
        ENDIF
        CALL oasis_get(portin_id(jfld), run_date_secs, exfld, info)
      ENDIF

      IF (.NOT. lcouple_parallel) CALL p_bcast(info, p_io)
      CALL digest_get_Id (info,clstr8o2a(jfld),run_date_secs)

      IF ( info == oasis_Recvd        .OR. &
           info == oasis_FromRest     .OR. &
           info == oasis_RecvOut      .OR. &
           info == oasis_FromRestOut         ) THEN
        CALL put_o2a(clstr8o2a(jfld),exfld)
      ENDIF

    ENDDO

    !
    !-- Accumulate time passed since last averaging.
    !
    couple_a2o_time = couple_a2o_time + delta_time

    IF (ldebugcpl) THEN
      CALL message('','End of couple_get_o2a')
      CALL message('','*********************')
      FLUSH(nerr)
    ENDIF

#endif

  END SUBROUTINE couple_get_o2a

!-------------------------------------------------------------------------------
!>
!! Sends exchange fields (for ocean) to oasis
!!
!! 03.05.15  S. Legutke - created
  SUBROUTINE couple_put_a2o
 
#if defined (__prism)
    INTEGER :: istep ! model step since initialisation
    INTEGER :: jfld  ! loop index

    IF (ldebugcpl) THEN
      istep = get_time_step()
      CALL message('','Start of couple_put_a2o')
      CALL message('','***********************')
      WRITE (message_text,'(a,i18)') ' Run time (sec) passed                   : ',run_date_secs
      CALL message('',message_text)
      WRITE (message_text,'(a,e13.6)') ' Accum. time (sec) since last couple_put : ',couple_a2o_time
      CALL message('',message_text)
      WRITE (message_text,'(a,i12)') ' -                model run step         : ',istep-couple_resum
      CALL message('',message_text)
      WRITE (message_text,'(a,i12)') ' -                model step since init. : ',istep
      CALL message('',message_text)
      FLUSH(nerr)
    ENDIF
    !
    !-- Export exchange fields
    !  
    IF (l_putocean) THEN

      CALL water_budget_corr

      DO jfld = 1,nflda2o

        CALL get_a2o(clstr8a2o(jfld),exfld)

        IF (p_parallel_oasis) THEN
          IF (ldebugcpl) THEN
            WRITE (message_text,'(a,i4)') ' Exporting field no. ', jfld
            CALL message('',message_text)
            WRITE (message_text,'(a,a)')  ' -   with field name ', clstr8a2o(jfld)
            CALL message('',message_text)
          ENDIF

          info = oasis_ok
          CALL oasis_put(portout_id(jfld), run_date_secs, exfld, info)
          CALL digest_put_Id (info,clstr8a2o(jfld),run_date_secs)
        END IF

      ENDDO

    ELSE

      IF (ldebugcpl) &
         WRITE (message_text,'(a,i12)') ' No call of oasis_put at date (sec) ', run_date_secs
    END IF

    IF (ldebugcpl) THEN
      CALL message('','End of couple_put_a2o')
      CALL message('','*********************')
    END IF

#endif

  END SUBROUTINE couple_put_a2o
  
!-------------------------------------------------------------------------------
!>
!! Terminates the coupled run.
!! 03.05.15  S. Legutke - created
!!
  SUBROUTINE couple_end
 
#if defined (__prism)
 
    REAL(pwp) ::  zrcouple

    IF (ldebugcpl) THEN
      CALL message('','Start of couple_end')
      CALL message('','*******************')
      FLUSH(nerr)
    ENDIF

    IF (.NOT.l_putocean) THEN
      CALL message ('','This run does not stop at the end of coupled time step.', level=em_warn)
      WRITE (message_text,'(a,e13.6)')' Time (sec) passed since last put : ', couple_a2o_time
      CALL message('',message_text)
      WRITE (message_text,'(a,i12)')' Date in oasis_put                : ', namruntim
      CALL message('',message_text)
      zrcouple = 1.0_pwp/couple_a2o_time
      awust(:,:) = awust(:,:)*zrcouple
      awvst(:,:) = awvst(:,:)*zrcouple
      aiust(:,:) = aiust(:,:)*zrcouple
      aivst(:,:) = aivst(:,:)*zrcouple
      aifre(:,:) = aifre(:,:)*zrcouple/rhoh2o
      awfre(:,:) = awfre(:,:)*zrcouple/rhoh2o
      aiqre(:,:) = aiqre(:,:)*zrcouple
      aicon(:,:) = aicon(:,:)*zrcouple
      awhea(:,:) = awhea(:,:)*zrcouple
      awsol(:,:) = awsol(:,:)*zrcouple
      awsta(:,:) = awsta(:,:)*zrcouple
      IF (lcouple_co2) THEN
        co2atmos(:,:)    = co2atmos(:,:)*zrcouple
        co2flux_cpl(:,:) = co2flux_cpl(:,:)*zrcouple
       END IF
      CALL couple_put_a2o
    ENDIF
    !
    !-- Write restart file.
    !
    IF ( nmseq > 1) THEN
       CALL couple_restart_a2o
    END IF
    !     
    !-- Deallocate memory.
    !
    ierror = oasis_ok
    CALL oasis_terminate(ierror)

    IF (ierror /= oasis_ok) THEN
      WRITE (message_text,'(a,i8)') 'couple_end: rank = ', p_pe
      CALL message('', message_text)
      CALL message('', 'pb oasis_terminate ')
      CALL oasis_abort (0, 'couple_end', &
          'oasis_terminate failed for model '//prism_model_name)
    END IF

    IF (ldebugcpl) CLOSE(debug_unit)

    IF (ldebugcpl) THEN
      CALL message('','End of couple_end')
      CALL message('','*****************')
      FLUSH(nerr)
    ENDIF

#endif

  END SUBROUTINE couple_end

!-------------------------------------------------------------------------------
!> 
!! Writes restart file with exchange fields at end of coupled run
!! 03.08.19  S. Legutke - created
  SUBROUTINE couple_restart_a2o
 
#if defined (__prism)
    INTEGER :: jfld ! loop index
 
    IF (ldebugcpl) THEN
      CALL message('',' Start of couple_restart_a2o')
      CALL message('','****************************')
      FLUSH(nerr)
    ENDIF

    IF (p_parallel_oasis) THEN
      DO jfld = 1,nflda2o
        info = oasis_ok
#ifdef __oa3mct
        CALL message('',' ERROR: oasis_put_restart is not implemented in OASIS3-MCT yet!')
        CALL oasis_abort (0,'couple_restart_a2o', &
                            'oasis_put_restart is not avilable yet!')
#else
        CALL oasis_put_restart(portout_id(jfld), run_date_secs, info)
        CALL digest_put_Id(info,clstr8a2o(jfld), run_date_secs)
#endif
      ENDDO
    END IF

    IF (ldebugcpl) THEN
      CALL message('',' End of couple_restart_a2o')
      CALL message('','**************************')
      FLUSH(nerr)
    ENDIF

#endif

  END SUBROUTINE couple_restart_a2o

!-------------------------------------------------------------------------------
!>
!! Checks cell partitioning into lake, land, ocean surfaces.
!! 03.05.15  S. Legutke - created
!! TODO: check this routine!

  SUBROUTINE check_alake
 
#if defined (__prism)

    INTEGER  :: jl, jrow
    REAL(wp) :: cellarea

    IF (p_parallel_io) THEN

    ! lakes are only allowed on cells without ocean.
    ! therefore check whether this is fullfilled [else
    ! enforce by setting lake part of cell to zero if ocean part > 0]
    ! [Also: check whether slf+alake=1 in non-ocean cells.]

    DO jl=1,nlon
      DO jrow=1,ngl

        !  Print cells with lake
        IF ( gl_alake(jl,jrow) > 0.0_wp ) THEN
          WRITE (message_text,'(a,2i4)') 'lake in cell ', jl,jrow
          CALL message('',message_text)
          WRITE (message_text,'(a,2i4)') '  alake, slf ', gl_alake(jl,jrow), gl_slf(jl,jrow)
          CALL message('',message_text)
          !  Modify if inconsistent
          IF ( gl_slf(jl,jrow)+gl_alake(jl,jrow) > 1._wp ) THEN
            CALL message('','Inconsistent land/lake partioning at lon/lat !')
            cellarea=gl_slf(jl,jrow)+gl_alake(jl,jrow)
            gl_slf(jl,jrow)=gl_slf(jl,jrow)/cellarea
            gl_alake(jl,jrow)=gl_alake(jl,jrow)/cellarea
            WRITE (message_text,'(a,2i4)') '  alake, slf changed to :', gl_alake(jl,jrow), gl_slf(jl,jrow)
            CALL message('',message_text)
          END IF
        END IF
      END DO
    END DO

    ENDIF ! p_parallel_io
      
#endif
  END SUBROUTINE check_alake

!-------------------------------------------------------------------------------
!>
!! Prints window of SST (array tsw) around min/max.
!! 03.05.15  S. Legutke - created
!! TODO: check this routine!

  SUBROUTINE print_sst
    
#if defined (__prism)

    REAL(wp), POINTER :: gl_tsw  (:,:)

    REAL(wp) :: tswgmax, tswgmin, tswmax, tswmin
    INTEGER  :: igmax, igmin, jgmax, jgmin, imax, imin
    INTEGER  :: jx, jy, jpdx, jpdy

    IF (p_parallel_io) THEN
      ALLOCATE (gl_tsw  (nlon,ngl))
    ELSE
      gl_tsw   => NULL()
    END IF

    CALL gather_field(gl_tsw, tsw)

    jpdy =ngl
    jpdx =nlon
    tswgmax = -999.0_wp
    tswgmin = +999.0_wp
    igmax = 1
    igmin = 1
    jgmax = 1
    jgmin = 1
    DO jy = 1, jpdy
      tswmax = -999.0_wp
      tswmin = +999.0_wp
      imax = 1
      imin = 1
      DO jx = 1,jpdx
        IF(gl_tsw(jx,jy) > tswmax) THEN 
          IF(gl_slf(jx,jy) < 1.0_wp .AND. gl_alake(jx,jy) == 0.0_wp) THEN
            !                  IF(gl_tsw(jx,jy) > 1.) THEN
            tswmax = gl_tsw(jx,jy)
            imax = jx
          ENDIF
        ENDIF
        IF(gl_tsw(jx,jy)  <  tswmin) THEN 
          IF(gl_slf(jx,jy) < 1.0_wp .AND. gl_alake(jx,jy) == 0.0_wp) THEN
            !                  IF(gl_tsw(jx,jy) > 1.) THEN
            tswmin = gl_tsw(jx,jy)
            imin = jx
          ENDIF
        ENDIF
      ENDDO
      IF(tswmax > tswgmax) THEN 
        tswgmax = tswmax
        igmax = imax
        jgmax = jy
      ENDIF
      IF(tswmin < tswgmin) THEN 
        tswgmin = tswmin
        igmin = imin
        jgmin = jy
      ENDIF

      WRITE (message_text,'(a,i5,a)') 'Max/min of SST at j= ', jy,' : '
      CALL message ('',message_text)
      CALL message ('', ' ------------------------------------------------')
      WRITE (message_text,'(2i5,2f12.7)') imax,imin,gl_tsw(imax,jy),gl_tsw(imin,jy)
      CALL message ('',message_text)

    ENDDO

    WRITE (message_text,'(a,2i5,2f12.7)') ' SST near global max: ', igmax,jgmax,gl_tsw(igmax,jgmax)
    CALL message ('',message_text)
    CALL message ('', ' ------------------------------------------------')

    DO jy = MAX(jgmax-5,1),MIN(jgmax+5,jpdy)
       WRITE (message_text,'(11(1x,f7.2,1x))') (gl_tsw(jx,jy),jx=MAX(igmax-5,1),MIN(igmax+5,jpdx))
       CALL message('',message_text)
    ENDDO

    WRITE (message_text,'(a,2f7.4)') ' slf there (1/0=land/sea only): ',gl_slf(igmax,jgmax)
    CALL message ('',message_text)
    CALL message ('', ' ------------------------------------------------')

    DO jy = MAX(jgmax-5,1),MIN(jgmax+5,jpdy)
       WRITE (message_text,'(11(1x,f7.4,1x))') (gl_slf(jx,jy),jx=MAX(igmax-5,1),MIN(igmax+5,jpdx))
       CALL message('',message_text)
    ENDDO

    WRITE (message_text,'(a,2i6,f7.4)') ' SST near global min: ',igmin,jgmin, gl_tsw(igmin,jgmin)
    CALL message ('',message_text)
    CALL message ('', ' ------------------------------------------------')

    DO jy = MAX(jgmin-5,1),MIN(jgmin+5,jpdy)
      WRITE (message_text,'(11(1x,f7.2,1x))') (gl_tsw(jx,jy),jx=MAX(igmin-5,1),MIN(igmin+5,jpdx))
      CALL message('',message_text)
    ENDDO

    WRITE (message_text,'(a,f7.2)') ' slf there (1/0=land/sea only): ', gl_slf(igmin,jgmin)
    CALL message('',message_text)
    CALL message ('', ' ------------------------------------------------')

    DO jy = MAX(jgmin-5,1),MIN(jgmin+5,jpdy)
      WRITE (message_text,'(11(1x,f7.4,1x))') (gl_slf(jx,jy),jx=MAX(igmin-5,1),MIN(igmin+5,jpdx))
      CALL message('',message_text)
    ENDDO

    tswgmax = -999.0_wp
    tswgmin = +999.0_wp
    igmax = 1
    igmin = 1
    jgmax = 1
    jgmin = 1
    DO jy = 1, jpdy
      tswmax = -999.0_wp
      tswmin = +999.0_wp
      imax = 1
      imin = 1
      DO jx = 1,jpdx
        IF(gl_tsw(jx,jy) > tswmax) THEN 
          IF(gl_slf(jx,jy) < 1.0_wp .AND. gl_alake(jx,jy)==0.0_wp) THEN
            tswmax = gl_tsw(jx,jy)
            imax = jx
          ENDIF
        ENDIF
        IF(gl_tsw(jx,jy) < tswmin) THEN 
          IF(gl_slf(jx,jy) < 1.0_wp .AND. gl_alake(jx,jy) == 0.0_wp) THEN
            tswmin = gl_tsw(jx,jy)
            imin = jx
          ENDIF
        ENDIF
      ENDDO
      IF(tswmax > tswgmax) THEN 
        tswgmax = tswmax
        igmax = imax
        jgmax = jy
      ENDIF
      IF(tswmin < tswgmin) THEN 
        tswgmin = tswmin
        igmin = imin
        jgmin = jy
      ENDIF
      WRITE (message_text,'(a,i5,a)') 'Max/min of SST at j= ',jy,' : '
      CALL message('',message_text)
      WRITE (message_text,'(2i5,2f12.7)') imax,imin,gl_tsw(imax,jy),gl_tsw(imin,jy)
      CALL message('',message_text)
    ENDDO

#endif

  END SUBROUTINE print_sst

!-------------------------------------------------------------------------------
!>
!! Prints land/sea mask characteristics.
!! 03.05.15  S. Legutke - created
!! TODO: Output needs to be reformatted as it does not look nice this way.
!!       (R.Redler)

  SUBROUTINE print_masks

    REAL(wp)    :: slfgmin, slfgmax, slfmin, slfmax, slfsum, alakesum

    INTEGER :: igmin, igmax, imin, imax
    INTEGER :: jgmin, jgmax
    INTEGER :: jx, jy

    CHARACTER(len=3) :: cmask(nlon,ngl)

#if defined (__prism)

    !-- Land mask characteristics

    slfgmin = 2.0_wp
    slfgmax =-1.0_wp
    DO jy = 1, ngl
      imax = 0
      imin = 0
      slfsum = 0.0_wp
      alakesum = 0.0_wp
      slfmin = 2.0_wp
      slfmax =-1.0_wp
      DO jx = 1,nlon
        IF(gl_slf(jx,jy) < 1.0_wp .AND. gl_alake(jx,jy) == 0.0_wp) THEN
          IF(gl_slf(jx,jy) > slfmax) THEN
            slfmax = gl_slf(jx,jy)
            imax = jx
          ENDIF
        ENDIF
        IF(gl_slf(jx,jy) > 0.0_wp .AND. gl_alake(jx,jy) == 0.0_wp) THEN
          IF(gl_slf(jx,jy) < slfmin) THEN
            slfmin = gl_slf(jx,jy)
            imin = jx
          ENDIF
        ENDIF
        slfsum   = slfsum   + gl_slf(jx,jy)
        alakesum = alakesum + gl_alake(jx,jy)
        IF(slfmax > slfgmax) THEN 
          slfgmax = slfmax
          igmax = imax
          jgmax = jy
        ENDIF
        IF(slfmin  <  slfgmin) THEN 
          slfgmin = slfmin
          igmin = imin
          jgmin = jy
        ENDIF
      ENDDO
    ENDDO

    WRITE (message_text,'(a,3i4)') &
    'Global min of fractional land area on cells with dry part but without lakes: ', slfgmin, igmin, jgmin
    CALL message('',message_text)
    WRITE (message_text,'(a,3i4)') &
    'Global max of fractional land area on cells with wet part but without lakes: ', slfgmax, igmax, jgmax
    CALL message('',message_text)
       
    !       Print land mask

    DO jy = 1, ngl
      DO jx = 1,nlon
        IF     (gl_slf(jx,jy) == 0.0_wp) THEN
          cmask(jx,jy) = '  .'
        ELSEIF (gl_slf(jx,jy) == 1.0_wp) THEN
          cmask(jx,jy) = ' **'
        ELSEIF (gl_slf(jx,jy)<1.0_wp .AND. gl_slf(jx,jy) > 0.0_wp) THEN
          WRITE(cmask(jx,jy),'(1X,I2.2)') INT(gl_slf(jx,jy)*100.0_wp)
        ELSE
          WRITE (message_text,'(a,2i4,f7.4)') ' Invalid fractional land area on ', jx, jy, gl_slf(jx,jy)
          CALL message('',message_text)
          CALL oasis_abort (0,'print_masks', &
              'pb fractional land mask in model '//prism_model_name)
        ENDIF
      ENDDO
    ENDDO

    CALL message('',' Partial land mask (%, jx=1,64) :')
    CALL message('','   .: no land / **: no water / >%')

    WRITE (message_text,'(64(1X,I2))') (MOD(jx,10),jx=1,nlon/2)
    CALL message('',message_text)
    DO jy = 1, ngl
      WRITE (message_text,'(64a3)') (cmask(jx,jy),jx=1,nlon/2)
      CALL message('',message_text)
    ENDDO

    CALL message('',' Partial sea/land mask (%, jx=65,128) :')
    WRITE (message_text,'(64(1X,I2))') (MOD(jx,10),jx=nlon/2+1,nlon)
    CALL message('',message_text)

    DO jy = 1, ngl
      WRITE (message_text,'(64a3)') (cmask(jx,jy),jx=nlon/2+1,nlon)
      CALL message('',message_text)
    ENDDO

    !       Print lake mask

    DO jy = 1, ngl
      DO jx = 1,nlon
        cmask(jx,jy) = '  .'
        IF(gl_slf(jx,jy) > 0.0_wp) cmask(jx,jy) = ' **'
        IF(gl_alake(jx,jy) > 0.0_wp) THEN
          WRITE(cmask(jx,jy),'(1X,I2)') INT(gl_alake(jx,jy)*100.0_wp+1.0_wp)
        ENDIF
      ENDDO
    ENDDO

    CALL message('', ' Lake mask (%, jx=1,64) :')
    CALL message('', '   .: no land / **: no lake / <%')
    WRITE (message_text,'(64(1X,I2))') (MOD(jx,10),jx=1,nlon/2)
    CALL message('',message_text)

    DO jy = 1, ngl
      WRITE (message_text,'(64a3)') (cmask(jx,jy),jx=1,nlon/2)
      CALL message('',message_text)
    ENDDO

    CALL message('', ' Lake mask (%, jx=65,128) :')

    WRITE (message_text,'(64(1X,I2))') (MOD(jx,10),jx=nlon/2+1,nlon)
    CALL message('',message_text)

    DO jy = 1, ngl
      WRITE (message_text,'(64a3)') (cmask(jx,jy),jx=nlon/2+1,nlon)
      CALL message('',message_text)
    ENDDO

    !       Print ocean mask

    DO jy = 1, ngl
      DO jx = 1,nlon
        cmask(jx,jy) = '  .'
        IF(gl_slf(jx,jy) < 1.0_wp .AND. gl_alake(jx,jy) == 0.0_wp) THEN
          WRITE(cmask(jx,jy),'(1X,I2)') INT((1.-gl_slf(jx,jy))*100.0_wp)
        ENDIF
      ENDDO
    ENDDO

    CALL message('', ' Ocean mask (%, jx=1,64) :')
    CALL message('', '   .: no ocean / **: only / >%')

    WRITE (message_text,'(64(1X,I2))') (MOD(jx,10),jx=1,nlon/2)

    DO jy = 1, ngl
      WRITE (message_text,'(64a3)') (cmask(jx,jy),jx=1,nlon/2)
      CALL message('',message_text)
    ENDDO

    CALL message('', ' Ocean mask (%, jx=65,128) :')
    WRITE (message_text,'(64(1X,I2))') (MOD(jx,10),jx=nlon/2+1,nlon)
    CALL message('',message_text)

    DO jy = 1, ngl
      WRITE (message_text,'(64a3)') (cmask(jx,jy),jx=nlon/2+1,nlon)
      CALL message('',message_text)
    ENDDO

#endif

  END SUBROUTINE print_masks

!-------------------------------------------------------------------------------
!>
!! Checks coupled model control parameter against those of the calling model
!! ECHAM6.
!!
!! 03.05.15  S. Legutke - created
  SUBROUTINE chck_par(lmres,klon,kgl)

    USE mo_time_conversion, ONLY: print_date, OPERATOR(==), OPERATOR(<)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: klon, kgl ! dimensions of model grid
    LOGICAL, INTENT(in) :: lmres     ! logical switch for restart mode

#if defined (__prism)

    INTEGER             :: nfldsize  ! fields size of exchange fields

    CHARACTER(len=32)   :: date_text

    nfldsize = klon*kgl
    
    !     
    !-- Check whether month of restart file is ok.
    !   ------------------------------------------
    
    IF (stop_date < next_date .OR. stop_date == next_date) THEN
      CALL message('','stop date <= next date!')
      CALL print_date(next_date, mess=date_text)
      WRITE (message_text,'(a,a)') 'next date = ', TRIM(date_text)
      CALL message('',message_text)
      CALL print_date(stop_date, mess=date_text)
      WRITE (message_text,'(a,a)') 'stop date = ', TRIM(date_text)
      CALL message('',message_text)
      CALL message('', 'STOP in chck_par')
      CALL oasis_abort (0,'chck_par','pb with dates in model '//prism_model_name)
    ENDIF

    !     
    !-- Check restart mode.
    !   -------------------

    IF (lmres) THEN
      CALL message('',' Model is restarted!')
    ELSE
      CALL message('',' Model is initialized!')
    ENDIF

#endif

  END SUBROUTINE chck_par

!-------------------------------------------------------------------------------
!>
!! Write grids, masks and areas for oasis.
!!
!! July 7, 2003  V. Gayler - created
  SUBROUTINE grids_writing

  INTEGER, PARAMETER   :: nc = 4            ! number of corners per grid cell
  INTEGER              :: gwrite            ! flag to state whether grids writing
                                            ! is needed or not. (1 / 0)
  INTEGER              :: i, j              ! looping indices
  CHARACTER*4          :: grdacr            ! grid acronym (as used in namcouple)
  CHARACTER*4          :: grdacr_2          ! grid acronym (as used in namcouple)
  REAL(wp)             :: lon(nlon,ngl)     ! 2dim array of longitudes
  REAL(wp)             :: lat(nlon,ngl)     ! 2dim array of latitudes 
  REAL(wp)             :: clon(nlon,ngl,nc) ! 3dim array of corner longitudes
  REAL(wp)             :: clat(nlon,ngl,nc) ! 3dim array of corner latitudes 
  INTEGER              :: mask(nlon,ngl)    ! land sea mask (0 for all cells
                                            !    with wet area fraction > 0.)
  INTEGER              :: mask_slm(nlon,ngl)! land sea mask (corresponding to
                                            !    the echam land sea mask slm)
  REAL(wp)             :: area(nlon,ngl)    ! grid cell areas corresponding to slf
  REAL(wp)             :: area_slm(nlon,ngl)!           corresponding to slm
  REAL(wp)             :: slf_sn(nlon,ngl)  ! land sea mask in oasis order (S->N)

#if defined (__prism)
  !> Write grids out only on p_parallel_io
  IF (p_parallel_io) THEN

    IF (ldebugcpl) THEN
      CALL message('',' Start of grids_writing')
      CALL message('','***********************')
    ENDIF

    CALL oasis_start_grids_writing(gwrite)

    IF (gwrite == 1) THEN  ! currently always 1 for oasis3-mct
!
!-- create 2d arrays of longitudes and latitudes
!
      DO i = 1, nlon
        lat(i,:) = philat(:)
      ENDDO
      DO j = 1, ngl
        lon(:,j) = philon(:)
      ENDDO
      WHERE (lon(:,:) < 0._wp)
        lon = lon + 360._wp
      END WHERE
      WHERE (lon(:,:) >= 360._wp)
        lon = lon - 360._wp
      END WHERE
!
!-- create 3d arrays of grid cell corner longitudes and latitudes
!   (grid cell corners must be written in counterclockwise sense)
!      The writing of corners is optional. If they are missing in the grids
!      file they will be calculated by scrip. In this case it is better to
!      calculate the corners by the model. (Scrip-corner will not reach the
!      poles.) 

      DO i = 1, nlon-1
        clon(i,:,1) = 0.5_wp * (lon(i+1,1) + lon(i,1))
        clon(i,:,4) = 0.5_wp * (lon(i+1,1) + lon(i,1))
      ENDDO

      DO i = 2, nlon
        clon(i,:,2) = 0.5_wp * (lon(i-1,1) + lon(i,1))
        clon(i,:,3) = 0.5_wp * (lon(i-1,1) + lon(i,1))
      ENDDO

      clon(nlon,:,1) = 0.5_wp * (lon(1,1) + lon(nlon,1) + 360._wp)
      clon(nlon,:,4) = 0.5_wp * (lon(1,1) + lon(nlon,1) + 360._wp)
      clon(   1,:,2) = 0.5_wp * (lon(nlon,1) + lon(1,1) - 360._wp)
      clon(   1,:,3) = 0.5_wp * (lon(nlon,1) + lon(1,1) - 360._wp)

      IF (clon(nlon,1,1) >= 360._wp) clon(nlon,:,1) = clon(nlon,:,1) - 360._wp
      IF (clon(nlon,1,4) >= 360._wp) clon(nlon,:,4) = clon(nlon,:,4) - 360._wp
      IF (clon(nlon,1,2) < 0._wp)    clon(nlon,:,2) = clon(nlon,:,2) + 360._wp
      IF (clon(nlon,1,3) < 0._wp)    clon(nlon,:,3) = clon(nlon,:,3) + 360._wp

      DO j = 1, ngl-1
        clat(:,j,3) = 0.5_wp * (lat(1,j+1) + lat(1,j))
        clat(:,j,4) = 0.5_wp * (lat(1,j+1) + lat(1,j))
      ENDDO

      DO j = 2, ngl
        clat(:,j,1) = 0.5_wp * (lat(1,j-1) + lat(1,j))
        clat(:,j,2) = 0.5_wp * (lat(1,j-1) + lat(1,j))
      ENDDO

      clat(:,  1,1) =  90._wp
      clat(:,  1,2) =  90._wp
      clat(:,ngl,3) = -90._wp
      clat(:,ngl,4) = -90._wp
!
!-- create 2d land ocean mask
!
      slf_sn(:,:) = gl_slf(:,:) + gl_alake(:,:)
!
!-- create 2d integer sea land mask
!
      mask(:,:) = 0
      WHERE (slf_sn(:,:) > 0.9999999_wp)
        mask(:,:) = 1
      END WHERE

      mask_slm(:,:) = INT(gl_slm(:,:))
      WHERE (gl_alake(:,:) > 0.0001_wp)
        mask_slm(:,:) = 1
      END WHERE
!
!-- create 2d array of of grid cell area
!     
!     The area array is used for CONSERV analysis of OASIS. 
!     In ECHAM6, fluxes over water and fluxes over land are calculated
!     separately. But the final flux going into the lowest grid cell is either
!     the flux calculated over water or the flux calculated over land, 
!     depending on the 1/0-land-sea-mask (slm); (compare ioinitial).
!     In this configuration CONSERV analysis is used for arrays going
!     FROM THE ATMOSPHERE TO THE OCEAN. Here all grid cells with a wet 
!     area fraction greater 0 are considered.
!
!     If echam uses a fractional land sea mask (not only for coupling)
!     (lfractional_mask) the 'atml' arrays are no longer needed.

      DO j = 1, ngl
        area(:,j) = gridarea(j)
      ENDDO
      area_slm(:,:) = area(:,:)

      WHERE (slf_sn(:,:) > 0.0_wp .AND. mask(:,:) == 0)
        area(:,:) = (1.0_wp-slf_sn(:,:)) * area(:,:)
      END WHERE

      WHERE (mask(:,:) == 1)
        area(:,:) = 0.0_wp
      END WHERE

      WHERE (mask_slm(:,:) == 1)
        area_slm(:,:) = 0.0_wp
      END WHERE

      grdacr='atmo'
      grdacr_2='atml'   ! the mask on this grid corresponds to slm (0/1 mask)
      IF (ldebugcpl) THEN
        WRITE (message_text,'(a,a,a)') grdacr, ', ', grdacr_2
        CALL message('oasis_write_grid',message_text)
      ENDIF
      CALL oasis_write_grid (grdacr, nlon, ngl, lon(:,:), lat(:,:))
      CALL oasis_write_grid (grdacr_2, nlon, ngl, lon(:,:), lat(:,:))

      IF (ldebugcpl) THEN
        WRITE (message_text,'(a)') grdacr
        CALL message('oasis_write_corner',message_text)
      ENDIF
      CALL oasis_write_corner (grdacr, nlon, ngl, nc, clon(:,:,:), clat(:,:,:))
      CALL oasis_write_corner (grdacr_2, nlon, ngl, nc, clon(:,:,:), clat(:,:,:))
 
      IF (ldebugcpl) THEN
        WRITE (message_text,'(a)') grdacr
        CALL message('oasis_write_mask',message_text)
      ENDIF
      CALL oasis_write_mask (grdacr, nlon, ngl, mask(:,:))
      CALL oasis_write_mask (grdacr_2, nlon, ngl, mask_slm(:,:))

      IF (ldebugcpl) THEN
        WRITE (message_text,'(a)') grdacr
        CALL message('oasis_write_area',message_text)
      ENDIF
      CALL oasis_write_area (grdacr, nlon, ngl, area(:,:))
      CALL oasis_write_area (grdacr_2, nlon, ngl, area_slm(:,:))

      IF (ldebugcpl) CALL message('','oasis_terminate_grids_writing')
      CALL oasis_terminate_grids_writing

    ELSE
      IF (ldebugcpl) CALL message('','  grids files are existing, no writing needed')
    ENDIF  ! gwrite

    IF (ldebugcpl) THEN
      CALL message('',' End of grids_writing')
      CALL message('','*********************')
    ENDIF

  ENDIF ! p_parallel_io

#endif /*__prism*/

  END SUBROUTINE grids_writing

!-------------------------------------------------------------------------------
!>
!! Initializes outgoing exchange fields
!! 03.05.15  S. Legutke - created
  SUBROUTINE ini_a2o

#if defined (__prism)
    IF (ldebugcpl) THEN
      CALL message('',' Start of ini_a2o')
      CALL message('','*****************')
      FLUSH(nerr)
    ENDIF

    awhea(:,:) = 0.0_wp
    awsol(:,:) = 0.0_wp
    awfre(:,:) = 0.0_wp
    awust(:,:) = 0.0_wp
    awvst(:,:) = 0.0_wp
    awsta(:,:) = 0.0_wp
    aicon(:,:) = 0.0_wp
    aiqre(:,:) = 0.0_wp
    aifre(:,:) = 0.0_wp
    aiust(:,:) = 0.0_wp
    aivst(:,:) = 0.0_wp
    IF (lcouple_co2) THEN
      co2atmos(:,:)    = 0.0_wp
      co2flux_cpl(:,:) = 0.0_wp
    END IF
    
    IF (ldebugcpl) THEN
      CALL message('',' End of ini_a2o')
      CALL message('','***************')
      FLUSH(nerr)
    ENDIF
#endif
    
  END SUBROUTINE ini_a2o


!-------------------------------------------------------------------------------
!>
!! Move model field to exchange array.
!! 03.05.15  S. Legutke - created

  SUBROUTINE get_a2o(clfield,pfield)

  IMPLICIT NONE

  CHARACTER(len=8), INTENT(in) :: clfield

#if defined (__prism)
  REAL(pwp), POINTER        :: pfield(:,:)
#else
  REAL(wp), POINTER         :: pfield(:,:)
#endif

  INTEGER :: j

#if defined (__prism)

  IF (ldebugcpl) THEN
    CALL message('',' Start of get_a2o')
    CALL message('','*****************')
    FLUSH(nerr)
  ENDIF

  SELECT CASE (TRIM(clfield))
    !> zonal wind stress over water
    CASE ('TXWATMOU')
      IF (lcouple_parallel) THEN
        !pfield = awust
        pfield = reshape(awust,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, awust)
      ENDIF

    !> meridional wind stress over water
    CASE ('TYWATMOU')
      IF (lcouple_parallel) THEN
        !pfield = awvst
        pfield = reshape(awvst,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, awvst)
      ENDIF

    !> zonal wind stress over ice
    CASE ('TXIATMOU')
      IF (lcouple_parallel) THEN
        !pfield = aiust
        pfield = reshape(aiust,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, aiust)
      ENDIF
 
    !> meridional wind stress over ice
    CASE ('TYIATMOU')
      IF (lcouple_parallel) THEN
        !pfield = aivst
        pfield = reshape(aivst,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, aivst)
      ENDIF

    !> snow flux on ice
    CASE ('FRIATMOS')
      IF (lcouple_parallel) THEN
        !pfield = aifre
        pfield = reshape(aifre,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, aifre)
      ENDIF

    !> water flux into the ocean
    CASE ('FRWATMOS')
      IF (lcouple_parallel) THEN
        !pfield = awfre
        pfield = reshape(awfre,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, awfre)
      ENDIF

    !> residual heat flux (sea-ice topmelt heat flux)
    CASE ('RHIATMOS')
      IF (lcouple_parallel) THEN
        !pfield = aiqre
        pfield = reshape(aiqre,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, aiqre)
      ENDIF

    !> heat flux in sea ice
    CASE ('CHIATMOS')
      IF (lcouple_parallel) THEN
        !pfield = aicon
        pfield = reshape(aicon,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, aicon)
      ENDIF

    !> heat flux into the ocean
    CASE ('NHWATMOS')
      IF (lcouple_parallel) THEN
        !pfield = awhea
        pfield = reshape(awhea,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, awhea)
      ENDIF

    !> solar radiation into the ocean
    CASE ('SHWATMOS')
      IF (lcouple_parallel) THEN
        !pfield = awsol
        pfield = reshape(awsol,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, awsol)
      ENDIF

#if defined __cpl_mpiom
    !> 10m wind speed
    CASE ('WSVATMOS')
      IF (lcouple_parallel) THEN
        !pfield = awsta
        pfield = reshape(awsta,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, awsta)
      ENDIF
#endif

    !> CO2 concentration (used for diagnostics in hamocc)
    CASE ('CO2CONAT')
      IF (lcouple_parallel) THEN
        !pfield = co2atmos
        pfield = reshape(co2atmos,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, co2atmos)
      ENDIF

    !> CO2 flux
    CASE ('CO2FLXAT')
      IF (lcouple_parallel) THEN
        !pfield = co2flux_cpl
        pfield = reshape(co2flux_cpl,(/ld%nglon,ld%nglat/))
      ELSE
        CALL gather_field(pfield, co2flux_cpl)
      ENDIF

    CASE DEFAULT
      CALL oasis_abort (0,'get_a2o','Invalid locator string '//TRIM(clfield)// &
                        &' for model '//prism_model_name)
    END SELECT

    IF (ldebugcpl) THEN
      IF (lcouple_parallel) THEN
        WRITE(debug_unit,*) ' ',TRIM(clfield),' : ',MINVAL(pfield), MAXVAL(pfield)
      ELSEIF (p_parallel_io) THEN
        WRITE(debug_unit,*) ' ',TRIM(clfield),' : ',(pfield(nlon/2,j),j=1,ngl,14)
      ENDIF
      FLUSH(debug_unit)
      CALL message('',' End of get_a2o')
      CALL message('','***************')
      FLUSH(nerr)
    END IF

#endif

  END SUBROUTINE get_a2o
  
!-------------------------------------------------------------------------------
!>
!! Moves the data received from OASIS to the respective 
!!  fields defined on the model grid.
!!  Since seaice/siced/sni are also used for lake surfaces,
!!  use only those values where slf<1 and alake=0 (criterion for
!!  ocean surface in grid cell).
!!
!! 03.05.15  S. Legutke - created
  SUBROUTINE put_o2a(clfield,pfield)

    IMPLICIT NONE

    CHARACTER(LEN=8), INTENT(IN) :: clfield           ! exchange field symbolic name

#if defined (__prism)
    REAL(pwp), POINTER, INTENT(IN) :: pfield(:,:)       ! array with received field
#else
    ! gl_* is only declared as wp, so we can make the truncation already here
    REAL(wp), POINTER, INTENT(IN)  :: pfield(:,:)       ! array with received field
#endif

    INTEGER :: npad
    INTEGER :: i, j, jrow, jl   ! loop indices

#if defined (__prism)

    IF (ldebugcpl) THEN
      CALL message('',' Start of put_o2a')
      CALL message('','************************')
      FLUSH(nerr)
    ENDIF
 
    ! start receiving and processing
    IF (lcouple_parallel) THEN
       !tmp_loco2a = pfield   ! Pointer not possible due to different precision (wp vs. pwp)
       !npad = ld% nproma - ld% npromz + 1
       npad = size(pfield,1)
       tmp_loco2a = reshape(pfield,(/nbdim,ngpblks/),(/(-999.0_pwp, i=1,npad)/))
    ELSE
       CALL scatter_field (pfield, tmp_loco2a)
    ENDIF

    SELECT CASE (TRIM(clfield))
      !> sea surface temperature
      CASE ('SSTATMOS')
        WHERE (locean_cell)
          tsw = tmp_loco2a
        ENDWHERE

        DO jrow = 1, ngpblks
          nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
          DO jl = 1, nproma
              IF( tsw(jl,jrow) < 270.0_wp .AND. locean_cell(jl,jrow)) THEN
                WRITE (message_text,'(a,3f7.4,2i4)') 'slf,alake,tsw:', &
                   slf(jl,jrow), alake(jl,jrow), tsw(jl,jrow), jl, jrow
                CALL message('',message_text,level=em_warn)
              ENDIF
          ENDDO
        ENDDO
        
       ! DO jrow = 1, ngpblks
       !   nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
       !   DO jl = 1, nproma
       !     IF (locean_cell(jl,jrow)) THEN
       !       tsw(jl,jrow) = tmp_loco2a(jl,jrow)
       !       IF( tsw(jl,jrow) < 270.0_wp) THEN
       !         WRITE (message_text,'(a,3f7.4,2i4)') 'slf,alake,tsw:', &
       !            slf(jl,jrow), alake(jl,jrow), tsw(jl,jrow), jl, jrow
       !         CALL message('',message_text,level=em_warn)
       !       ENDIF
       !     ENDIF
       !   ENDDO
       ! ENDDO

      !> sea ice thickness
      CASE ('SITATMOS')
        WHERE (locean_cell)
          siced = tmp_loco2a
        ENDWHERE
        !DO jrow = 1, ngpblks
        !  nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
        !  WHERE (locean_cell(1:nproma,jrow))
        !    siced(1:nproma,jrow) = tmp_loco2a(1:nproma,jrow)
        !  END WHERE
        !ENDDO

      !> sea ice area fraction
      CASE ('SICATMOS')
        WHERE (locean_cell)
          seaice = tmp_loco2a
        ENDWHERE
        !DO jrow = 1, ngpblks
        !  nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
        !  WHERE (locean_cell(1:nproma,jrow))
        !    seaice(1:nproma,jrow) = tmp_loco2a(1:nproma,jrow)
        !  END WHERE
        !ENDDO

      !> snow thickness over sea ice
      CASE ('SNTATMOS')
        WHERE (locean_cell)
          sni = tmp_loco2a
        ENDWHERE
        !DO jrow = 1, ngpblks
        !  nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
        !  WHERE (locean_cell(1:nproma,jrow))
        !    sni(1:nproma,jrow) = tmp_loco2a(1:nproma,jrow)
        !  END WHERE
        !ENDDO

      !> ocean surface velocity (E-W component)
      CASE ('OCUATMOS')
        WHERE (locean_cell)
          ocu = tmp_loco2a
        ENDWHERE
        !DO jrow = 1, ngpblks
        !  nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
        !  WHERE (locean_cell(1:nproma,jrow))
        !    ocu(1:nproma,jrow) = tmp_loco2a(1:nproma,jrow)
        !  END WHERE
        !ENDDO

      !> ocean surface velocity (N-S component)
      CASE ('OCVATMOS')
        WHERE (locean_cell)
          ocv = tmp_loco2a
        ENDWHERE
        !DO jrow = 1, ngpblks
        !  nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
        !  WHERE (locean_cell(1:nproma,jrow))
        !    ocv(1:nproma,jrow) = tmp_loco2a(1:nproma,jrow)
        !  END WHERE
        !ENDDO

      !> CO2 transfer coefficient
      CASE ('CO2TRAAT')
        WHERE (locean_cell)
          co2trans = tmp_loco2a * 1.e6_wp 
        ENDWHERE
        !DO jrow = 1, ngpblks
        !  nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
        !  WHERE (locean_cell(1:nproma,jrow))
        !    co2trans(1:nproma,jrow) = tmp_loco2a(1:nproma,jrow) * 1.e6_wp 
        !  END WHERE
        !ENDDO

      !> pCO2 in uppermost ocean layer
      CASE ('CO2ATMOS')
        WHERE (locean_cell)
          co2ocean = tmp_loco2a
        ENDWHERE
        !DO jrow = 1, ngpblks
        !  nproma = MERGE(ld% npromz, ld% nproma, jrow == ngpblks)
        !  WHERE (locean_cell(1:nproma,jrow))
        !    co2ocean(1:nproma,jrow) = tmp_loco2a(1:nproma,jrow)
        !  END WHERE
        !ENDDO

      CASE DEFAULT
        CALL oasis_abort(0,'put_o2a','Invalid locator string '// &
                         &TRIM(clfield)//' for model '//prism_model_name)
      END SELECT

      IF (ldebugcpl) THEN
        IF (lcouple_parallel) THEN
          WRITE(debug_unit,*) ' ',TRIM(clfield),' : ',MINVAL(pfield,MASK=locean_cell), MAXVAL(pfield,MASK=locean_cell)
        ELSEIF (p_parallel_io) THEN
          WRITE(debug_unit,*) ' ',TRIM(clfield),' : ',(pfield(nlon/2,j),j=1,ngl,14)
        ENDIF
        FLUSH(debug_unit)
        CALL message('',' End of put_o2a')
        CALL message('','***************')
        FLUSH(nerr)
      END IF

#endif

  END SUBROUTINE put_o2a

!-------------------------------------------------------------------------------
!>
!! Time control of coupling aspects
!! 03.05.22  S. Legutke - created
  SUBROUTINE couple_calendar(pdt)
    REAL(wp), INTENT(in)   :: pdt   ! model time step length (secs)

#if defined (__prism)

    IF (l_putocean) THEN

      !
      !-- Reset exchange fields and accumulation interval.
      !  

      IF ( .NOT. lbreak .OR. .NOT. l_putrerun ) THEN
        CALL ini_a2o
        couple_a2o_time = 0.0_pwp
      END IF

    END IF

#endif

    !> Accumulation of time passed in this run.
    run_date_secs = run_date_secs + INT(pdt)

    IF (ldebugcpl) THEN
      CALL message('couple_calendar',' Calendar is updated.')
      FLUSH(nerr)
    ENDIF

  END SUBROUTINE couple_calendar


!-------------------------------------------------------------------------------
!>
!! Print info after oasis_get. Abort in case of error.
!! 03.05.15  S. Legutke - created

  SUBROUTINE digest_get_Id(kinfo,clfield,kdate,kcount)

    INTEGER,            INTENT(IN)  :: kinfo   ! info ID passed from psmile
    INTEGER,            INTENT(IN)  :: kdate   ! date (seconds)
    INTEGER, OPTIONAL,  INTENT(IN)  :: kcount  ! exchanges field count
    CHARACTER(LEN=8),   INTENT(IN)  :: clfield ! fields symbolic name

#if defined (__prism)
    INTEGER  :: icount = -1

    IF (ldebugcpl) THEN
      IF (PRESENT(kcount)) icount = kcount
      WRITE (message_text,'(a,i12,1x,a)') ' At date (seconds) ', kdate, clfield
      CALL message('',message_text)

      IF (icount /= -1 ) THEN
         WRITE (message_text,'(a,i14,a)') ' (field no. ',kcount,')'
         CALL message('',message_text)
      ENDIF
      FLUSH(nerr)
    ENDIF

    IF ( kinfo /= oasis_ok          .AND. &
         kinfo /= oasis_FromRest    .AND. &
         kinfo /= oasis_Input       .AND. &
         kinfo /= oasis_RecvOut     .AND. &
         kinfo /= oasis_FromRestOut .AND. &
         kinfo /= oasis_Recvd            )  THEN
      WRITE (message_text,'(a,a)') ' Problem with port = ', clfield
      CALL message('',message_text)
      WRITE (message_text,'(a,i12)') ' Error code number = ', kinfo
      CALL message('',message_text)
      WRITE (message_text,'(a,i12)') ' Seconds passed    = ', kdate
      CALL message('',message_text)
      CALL oasis_abort (0,'couple_get_o2a:digest_get_Id', &
          'oasis_get failed for model '//prism_model_name)
    ENDIF

    IF(ldebugcpl) THEN
      SELECT CASE (kinfo)
        CASE (oasis_Recvd)
          CALL message('', ' was received from another model')

        CASE (oasis_ok)
          CALL message('', ' was not received; no error.')

        CASE (oasis_FromRest)
          CALL message('', ' was received from restart file.')
         
        CASE (oasis_Input)
          CALL message('', ' was received from input file.')

        CASE (oasis_RecvOut)
          CALL message('', ' was received from input file and other model.')

        CASE (oasis_FromRestOut)
          CALL message('', ' was received from input file and written to an output file.')
      END SELECT
      FLUSH(nerr)
    ENDIF

#endif

  END SUBROUTINE digest_get_Id

!-------------------------------------------------------------------------------
! BOP
!
! !IROUTINE:  digest_put_Id
!
! !INTERFACE:

  SUBROUTINE digest_put_Id(kinfo,clfield,kdate,kcount)

  IMPLICIT NONE

  INTEGER,            INTENT(in)  :: kinfo   ! info ID passed from psmile
  INTEGER,            INTENT(in)  :: kdate   ! date (seconds)
  INTEGER, OPTIONAL,  INTENT(in)  :: kcount  ! exchanges field count

  CHARACTER(len=8), INTENT(in)    :: clfield ! fields symbolic name

!
! !DESCRIPTION:
!
! - Print info after oasis_put. Abort if error.
!
! !REVISION HISTORY:
! 03.05.15  S. Legutke - created
!
! EOP
!-----------------------------------------------------------------------

      INTEGER              :: icount = -1

#if defined (__prism)

      IF (ldebugcpl) THEN
        IF (PRESENT(kcount)) THEN
          icount = kcount
        ENDIF

        WRITE (message_text,'(a,i12,1x,a)') ' At date (seconds) ', kdate, clfield
        CALL message('',message_text)
        IF (icount /= -1 ) THEN
          WRITE (message_text,'(a,i12,a)') ' (field no. ',kcount,')'
          CALL message('',message_text)
          FLUSH(nerr)
        ENDIF
      ENDIF

      IF ( kinfo /= oasis_ok        .AND. &
           kinfo /= oasis_LocTrans  .AND. &
           kinfo /= oasis_ToRest    .AND. &
           kinfo /= oasis_Output    .AND. &
           kinfo /= oasis_SentOut   .AND. &
           kinfo /= oasis_ToRestOut .AND. &
           kinfo /= oasis_Sent            )  THEN
         WRITE (message_text,'(a,a)')   ' Problem with port = ',clfield
         CALL message('',message_text)
         WRITE (message_text,'(a,i12)') ' Error code number = ',kinfo
         CALL message('',message_text)
         WRITE (message_text,'(a,i12)') ' Seconds passed    = ',kdate
         CALL message('',message_text)
         CALL oasis_abort (0,'couple_put_a2o:digest_put_Id', &
             'oasis_put failed for model '//prism_model_name)
      ENDIF

      IF (ldebugcpl) THEN
        SELECT CASE (kinfo)
          CASE (oasis_Sent)
            CALL message('', ' was sent to another model')
 
          CASE (oasis_ok)
            CALL message('', ' was not sent; no error.')
 
          CASE (oasis_LocTrans)
            CALL message('', ' was used in local transformation.')
 
          CASE (oasis_ToRest)
            CALL message('', ' was written to a restart file.')
 
          CASE (oasis_Output)
            CALL message('', ' was output to a file.')
 
          CASE (oasis_SentOut)
            CALL message('', ' was sent to another model and output to a file.')

          CASE (oasis_ToRestOut)
            CALL message('', ' was sent to another model and written to a restart file.')
        END SELECT
        FLUSH(nerr)
      ENDIF

#endif

  END SUBROUTINE digest_put_Id
!----------------------------------------------------------------------

  SUBROUTINE water_budget_corr

    ! Two freshwater flux corrections are needed:
    !   - P-E correction: the echam water transport is not completely
    !     conservative. The budget change is diagnosed in physc (pmeb).
    !   - Land sea mask mismatch: jsbach uses a 1/0 land sea mask
    !     (SLM), whereas the fluxes to the ocean are calculated fractionally.
    ! Both corrections are scaled with the global rain and given to the ice-
    ! free ocean.
    !
    !  Author:
    !  U. Schlese, MPI, July 2003
    !  E. Roeckner, U. Schlese, MPI, February 2007
    !  T. Raddatz, MPI, November 2007
    !  V. Gayler,  MPI, December 2011

    IMPLICIT NONE

#if defined (__prism)

    ! local arrays

    REAL(wp), POINTER  :: zpmeb(:,:)   ! p-e correction, accumulated [kg/m**2]
    REAL(wp), POINTER  :: zrain(:,:)   ! total rain, accumulated     [kg/m**2]
    REAL(wp), POINTER  :: zawfre(:,:)  ! liquid freshwater flux [m/s]
    REAL(wp), POINTER  :: zaifre(:,:)  ! solid freshwater flux [m/s]
    REAL(wp), POINTER  :: zdisch(:,:)  ! river discharge [m/s]
    REAL(wp), POINTER  :: zresi(:,:)   ! freshwater flux due to land-sea-mask missmatch [m/s]
    REAL(wp), POINTER  :: zslf(:,:)    ! fractional land cover
    REAL(wp), POINTER  :: zslm(:,:)    ! echam land-sea mask (1. for land and 0. for ocean)
    REAL(wp), POINTER  :: zalake(:,:)  ! fractional lake mask
    REAL(wp), POINTER  :: zseaice(:,:) ! fraction of water covered by sea ice
    REAL(wp), POINTER  :: zwater(:,:)  ! fraction of grid area covered by ice-free ocean
    REAL(wp), POINTER  :: zocean(:,:)  ! fraction of grid area covered with ocean  
    REAL(wp), POINTER  :: olm(:,:)     ! ocean land mask [0,1]
    REAL(wp), POINTER  :: olf(:,:)     ! fractional ocean land mask

    REAL(wp) :: zoceanz(ngl)
    REAL(wp) :: zpmebz(ngl)
    REAL(wp) :: zrainz(ngl)
    REAL(wp) :: zresiz(ngl)

    ! Local scalars

    REAL(wp) :: zocean_glob
    REAL(wp) :: zpmeb_glob
    REAL(wp) :: zrain_glob
    REAL(wp) :: zresi_glob

    REAL(wp) :: rdiv

    INTEGER :: jlat

    IF (p_pe == p_io) THEN
      ALLOCATE (zpmeb(nlon,ngl))
      ALLOCATE (zrain(nlon,ngl))
      ALLOCATE (zawfre(nlon,ngl))
      ALLOCATE (zaifre(nlon,ngl))
      ALLOCATE (zdisch(nlon,ngl))
      ALLOCATE (zresi(nlon,ngl))
      ALLOCATE (zslf(nlon,ngl))
      ALLOCATE (zslm(nlon,ngl))
      ALLOCATE (zalake(nlon,ngl))
      ALLOCATE (zseaice(nlon,ngl))
      ALLOCATE (zwater(nlon,ngl))
      ALLOCATE (zocean(nlon,ngl))
      ALLOCATE (olm(nlon,ngl))
      ALLOCATE (olf(nlon,ngl))
    END IF

    CALL gather_field (zpmeb,          apmebco)
    CALL gather_field (zrain,          rain)
    CALL gather_field (zawfre,         awfre)
    CALL gather_field (zaifre,         aifre)
    CALL gather_field (zdisch,         disch)
    CALL gather_field (zslf,           slf)
    CALL gather_field (zslm,           slm)
    CALL gather_field (zalake,         alake)
    CALL gather_field (zseaice,        seaice)

!$OMP PARALLEL      
    IF (p_pe == p_io) THEN

    ! calculation of ocean masks and of the ice-free ocean fractional area

!$OMP DO PRIVATE(jlat)
      DO jlat = 1,ngl
        olf(:,jlat) = zslf(:,jlat) + zalake(:,jlat)                        ! fract. ocean land mask
        zocean(:,jlat) =  1._wp - olf(:,jlat)                              ! fract. ocean area
        zwater(:,jlat) = (1._wp - olf(:,jlat)) * (1._wp - zseaice(:,jlat)) ! ice-free ocean
        IF(.NOT. lfractional_mask)                              &
          olm(:,jlat) = MERGE(zslm(:,jlat),1._wp,zalake(:,jlat) < 0.5_wp)    ! 1/0 ocean land mask
      ENDDO
!$OMP END DO

      ! total rain and p-e correction [kg/m2] averaged over coupling time step [m/s]

      rdiv = 1.0_wp/(rhoh2o*couple_a2o_time)
!$OMP DO PRIVATE(jlat)
      DO jlat = 1,ngl
        zrain(:,jlat) = zrain(:,jlat)*rdiv
        zpmeb(:,jlat) = zpmeb(:,jlat)*rdiv
      ENDDO
!$OMP END DO

      ! global sums of rain and p-e correction [m/s]

!$OMP DO PRIVATE(jlat)
      DO jlat = 1,ngl
        zrainz(jlat) = SUM(zrain(:,jlat) * zwater(:,jlat)) * gl_budw(jlat)
        zpmebz(jlat) = SUM(zpmeb(:,jlat)) * gl_budw(jlat)
        zoceanz(jlat) = SUM(zocean(:,jlat)) * gl_budw(jlat)
        IF(lfractional_mask) THEN
          zresiz(jlat) = 0._wp
        ELSE
          zresiz(jlat) = SUM((zawfre(:,jlat)+zaifre(:,jlat)) * (olf(:,jlat)-olm(:,jlat))) * gl_budw(jlat)
        ENDIF
      ENDDO
!$OMP END DO

      zocean_glob = SUM(zoceanz(:))
      zrain_glob = SUM(zrainz(:))
      zpmeb_glob = SUM(zpmebz(:))
      zresi_glob = SUM(zresiz(:))

      IF (diag_water_budget) THEN
        WRITE (message_text,*) 'P-E correction:  ', zpmeb_glob * SUM(SPREAD(gridarea,1,nlon)) , ' m3/s'
        CALL message('water_budget_corr', message_text)
        WRITE (message_text,*) 'slm/slf residuum: ', zresi_glob * SUM(SPREAD(gridarea,1,nlon)) , ' m3/s'
        CALL message('water_budget_corr', message_text)
        IF (lfractional_mask) THEN
          WRITE (message_text,*) 'slm/slf residuum contribution by aifre: ', 0._wp, ' m3/s'
        ELSE
          WRITE (message_text,*) 'slm/slf residuum contribution by aifre: ', &
               SUM(zaifre(:,:) * (olf(:,:)-olm(:,:)) * SPREAD(gridarea,1,nlon)), ' m3/s'
        ENDIF
        CALL message('water_budget_corr', message_text)
      END IF

      ! Application of the corrections:
      !   add p-e correction and residual flux scaled with ocean area to
      !   the ocean freshwater flux [m/s]

      rdiv = 1.0_wp/zocean_glob

      ! Note: awfre is multiplied with the area of the wet grid cell fraction
      ! in oasis (CONSERV). The correction factor takes this into account:
      ! in contrast to above, zwater is calculated here without (fractional)
      ! land mask.  

      zwater(:,:) = 1._wp - zseaice(:,:)
      zresi(:,:) = zresi_glob * rdiv  ! diagnostic
      zpmeb(:,:) = zpmeb_glob * rdiv  ! diagnostic

      zawfre(:,:) = zawfre(:,:) + &
           (zpmeb_glob + zresi_glob) * rdiv

      IF (diag_water_budget) THEN
         ! awfre corrected =  awfre of mo_hydrology (precip-evap on ocean)
         !                  + zpmeb_glob (P-E correction)
         !                  + slm/slf residuum contribution by aifre
         WRITE (message_text,*) 'awfre corrected:  ', &
              SUM(zawfre(:,:) * (1._wp-olf(:,:)) * SPREAD(gridarea(:),1,nlon)), ' m3/s'
         CALL message('water_budget_corr', message_text)
      END IF

      ! Add river discharge to the freshwater
      ! Arrays defined on 1/0-mask are devided by 1-olf to account for later
      ! multiplication with the oasis areas file intended for fractional grid cells

      WHERE (olf(:,:) < 1._wp)
         zawfre(:,:) =   zawfre(:,:) +  zdisch(:,:) / (1._wp-olf(:,:))
      ELSEWHERE
         zawfre = 0._wp
         zaifre = 0._wp
      END WHERE

      IF (diag_water_budget) THEN
        WRITE (message_text,*) 'river discharge:  ', SUM(zdisch(:,:) * SPREAD(gridarea(:),1,nlon)), ' m3/s'
        CALL message('water_budget_corr', message_text)
        WRITE (message_text,*) 'water to ocean:   ', &
             SUM(zawfre(:,:) * (1._wp-olf(:,:)) * SPREAD(gridarea(:),1,nlon)), ' m3/s'
        CALL message('water_budget_corr', message_text)
      END IF

    ENDIF
!$OMP END PARALLEL

    CALL scatter_field (zawfre, awfre)
    CALL scatter_field (zaifre, aifre)
    CALL scatter_field (zresi, resi)
    CALL scatter_field (zpmeb, pmeb)

    IF (p_pe == p_io) THEN
      DEALLOCATE (zpmeb)
      DEALLOCATE (zawfre)
      DEALLOCATE (zaifre)
      DEALLOCATE (zdisch)
      DEALLOCATE (zresi)
      DEALLOCATE (zrain)
      DEALLOCATE (zslf)
      DEALLOCATE (zslm)
      DEALLOCATE (zalake)
      DEALLOCATE (zseaice)
      DEALLOCATE (zwater)
      DEALLOCATE (zocean)
      DEALLOCATE (olm)
      DEALLOCATE (olf)
    END IF
#endif

  END SUBROUTINE water_budget_corr

END MODULE mo_couple
