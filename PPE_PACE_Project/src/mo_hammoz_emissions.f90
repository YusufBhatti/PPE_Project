!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_emissions.f90
!!
!! \brief
!! Module for diagnostics of emissions in ECHAM submodels
!!
!! \author M. Schultz   (FZ Juelich)
!! \author S. Schroeder (FZ Juelich)
!!
!! \responsible_coder
!! M. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich), S. Schroeder (FZ Juelich) - original code (2010-02-01) ]
!!   -# A. Laakso (FMI) - Emissions for SALSA
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

MODULE mo_hammoz_emissions

  USE mo_kind,             ONLY: dp
  USE mo_memory_base,      ONLY: t_stream
  USE mo_submodel_diag,    ONLY: t_diag_list
  USE mo_species,          ONLY: nmaxspec
  USE mo_emi_matrix,       ONLY: maxsectors, EM_NONE, EM_SURFACE, EM_VOLUME, EM_LEVEL50M, EM_FIRE


  IMPLICIT NONE

  PRIVATE

  ! public variables  (see declaration below)

  ! subprograms
  PUBLIC                       :: init_emissions
  PUBLIC                       :: init_emi_stream
  PUBLIC                       :: emi_interface

  ! emi_stream
  INTEGER, PARAMETER           :: nemivars=1 
  CHARACTER(LEN=32)            :: emivars(1:nemivars)= &
                                (/'emi              '/)   ! emissions mass flux as 2D field
  ! Note: special submodels may add other variables to emi stream (e.g. megan or bgc_dust)
  LOGICAL            :: lhas_volcc, lhas_volce     ! sectors for continuous and explosive volcanoes defined
  LOGICAL            :: lhas_oceani                ! interactive ocean emissions active (DMS)
  LOGICAL            :: lhas_biogenic                ! interactive vegetation emissions active (organics)

  ! variable pointers and diagnostic lists
  TYPE (t_diag_list), PUBLIC   :: emi_diag      ! emissions mass flux diagnostics (total)
  TYPE (t_diag_list), PUBLIC   :: emi_diag_detail(maxsectors)   ! diagnostics per sector

  INTEGER               :: nsectors
  INTEGER               :: emimod(maxsectors,nmaxspec)     ! mode of application
  INTEGER               :: ibc_emis(maxsectors,nmaxspec)   ! index to boundary condition
  REAL(dp)              :: emfactor(maxsectors,nmaxspec)   ! scaling factor
  REAL(dp), ALLOCATABLE :: dgeo(:,:)                       ! delta(geopotential height)
  REAL(dp), ALLOCATABLE :: height(:,:)                     ! grid box height
  REAL(dp), ALLOCATABLE :: dens(:,:)                       ! air density

  CONTAINS

  !! ---------------------------------------------------------------------------
  !! subroutine to initialize emissions: read emi_spec.dat and set-up arrays

  SUBROUTINE init_emissions

  USE mo_emi_matrix,               ONLY: em_read, em_get_sector_info, em_get_bc_from_matrix
  USE mo_boundary_condition,       ONLY: bc_nml, bc_define, BC_BOTTOM
  USE mo_external_field_processor, ONLY: EF_MODULE, EF_INACTIVE
  USE mo_species,                  ONLY: speclist, query_species, itrprog
  USE mo_exception,                ONLY: message, message_text, em_info, em_warn, em_error, em_param
  USE mo_submodel,                 ONLY: lham, lemissions, lbioemi_stdalone
  USE mo_hammoz_emi_volcano,       ONLY: init_emi_volcano
  USE mo_hammoz_emi_ocean,         ONLY: init_emi_ocean
  USE mo_hammoz_emi_biogenic,      ONLY: init_emi_biogenic, lbioemi_dyn, start_biogenic_emissions
  USE mo_ham_m7_emissions,         ONLY: ham_m7_init_emissions
  USE mo_ham_species,              ONLY: id_so2, id_so4
  USE mo_ham,                      ONLY: ibc_dust, ibc_seasalt, &
                                         nham_subm,             &
                                         HAM_BULK,              &
                                         HAM_M7,                &
                                         HAM_SALSA 
  USE mo_tracdef,                  ONLY: trlist, ntrac
  USE mo_mpi,                      ONLY: p_parallel_io
  !++alaak
  USE mo_ham_salsa_emissions,      ONLY: ham_salsa_init_emissions
  !--alaak
  USE mo_ham_salsactl,             ONLY:  in2b, fn2b, in2a,fn2a,fn1a
  USE mo_decomposition,            ONLY: ldc=>local_decomposition

  TYPE(bc_nml)       :: bc_struc
  INTEGER            :: i, j, jspec, itrtype, nvars, idims, emtype, idum, jtrac, ii
  INTEGER            :: nproma, nlev
  REAL(dp)           :: zfactor, zeps
  CHARACTER(LEN=64)  :: secname, shortname
  CHARACTER(LEN=3)   :: cspec

  ! -- Initialisation
  zeps=EPSILON(1.0_dp)
  emimod(:,:) = EM_NONE
  emfactor(:,:) = 0._dp
  lhas_volcc = .FALSE.
  lhas_volce = .FALSE.
  lhas_oceani = .FALSE.
  lhas_biogenic = .FALSE.
  nlev    = ldc% nlev
  nproma  = ldc% nproma
  ALLOCATE(dgeo  (nproma,nlev)); dgeo   = 0.0_dp
  ALLOCATE(height(nproma,nlev)); height = 0.0_dp
  ALLOCATE(dens  (nproma,nlev)); dens   = 0.0_dp

  ! -- read emission matrix from emi_spec.dat file
  CALL em_read(nsectors)

  lemissions = .FALSE.
  DO i = 1, nsectors
    CALL em_get_sector_info(i, secname, nvars)
    DO j=1,nvars
      CALL em_get_bc_from_matrix(i, j, secname, bc_struc, emtype, zfactor, shortname)
      ! determine geometry of emissions
      ! emi_matrix only knows about option "surface" (corresponding to BC_BOTTOM)! 
      ! There's no correspondance to BC_TOP.
      IF (bc_struc%ef_type   == EF_MODULE) bc_struc%ef_actual_unit = 'kg m-2 s-1'
      IF (bc_struc%bc_domain == BC_BOTTOM) THEN
        idims = 2
      ELSE
        idims = 3
      END IF

      ! get species index
      CALL query_species(shortname=TRIM(shortname)//"_moz", index=jspec, itrtype=itrtype)
      IF (jspec <= 0) CALL query_species(shortname=shortname, index=jspec, itrtype=itrtype)
      IF (jspec > 0) THEN
        IF (itrtype == itrprog) THEN
          lemissions = .TRUE.
          write(cspec,'(i3)') jspec
          CALL message('','',level=em_param)
          CALL message('init_emissions', '  jspec='//cspec//': varname='//TRIM(shortname)//   &
                       ', sector='//TRIM(secname), level=em_info)
          !! use emi_matrix information to define a boundary condition
          ibc_emis(i,jspec) = bc_define(TRIM(secname)//' emissions of '//shortname, bc_struc, &
                                        idims, .TRUE.)
          ! mark dust/seasalt bcs for latter use in mo_ham_m7_emissions
          IF (TRIM(secname) == 'DUST') ibc_dust = ibc_emis(i,jspec)
          IF (TRIM(secname) == 'SEASALT') ibc_seasalt = ibc_emis(i,jspec)
          !! set emission handler
          emimod(i, jspec) = emtype
          emfactor(i, jspec) = zfactor    ! set scaling factor for emissions
          IF ((abs((zfactor - 1.0_dp)) > zeps) .AND. (p_parallel_io)) THEN
            write (message_text,'(a,e9.2)') 'emissions will be scaled with a factor of ', zfactor
            CALL message('init_emissions', message_text, level=em_info)
          ENDIF
          !! set emission flag in speclist
          speclist(jspec)%lemis = .TRUE.
          DO jtrac = 1,ntrac
            IF (trlist%ti(jtrac)%spid == jspec) trlist%ti(jtrac)%nemis=1 
          ENDDO

          !! if SO2 is found: pretend that also SO4 has been read from the emission matrix
          !! do NOT define a boundary condition like that for SO2 (if you do, the bc scheme expects
          !! another input file for SO4)

          IF ((jspec == id_so2) .AND. (id_so4 > 0)) THEN
            speclist(id_so4)%lemis = .true.
            emimod(i,id_so4) = emtype
            emfactor(i,id_so4) = zfactor
            DO jtrac = 1, ntrac
              IF (trlist%ti(jtrac)%spid == id_so4) trlist%ti(jtrac)%nemis=1 
            ENDDO
          ENDIF

          !! if "DUST by module" or "SEASALT by module" is found: another bc is needed!

          IF (((TRIM(secname) == 'DUST') .or. (TRIM(secname) == 'SEASALT')) .and. &
               (bc_struc%ef_type == EF_MODULE)) THEN

             IF (lham) THEN
                SELECT CASE(nham_subm)
                    CASE(HAM_BULK)
                        !do nothing
                    CASE(HAM_M7)
    
                        !! the second bc will get the next bc-index (therefore there is no need to store the index)
            
                        idum = bc_define(TRIM(secname)//' emissions of '//TRIM(shortname)//'2', bc_struc, &
                                         idims, .TRUE.)
                        IF (TRIM(secname) == 'SEASALT') THEN
                          idum = bc_define(TRIM(secname)//' emissions of '//TRIM(shortname)//'3', bc_struc, &
                                           idims, .TRUE.)
                          idum = bc_define(TRIM(secname)//' emissions of '//TRIM(shortname)//'4', bc_struc, &
                                           idims, .TRUE.)
                        ENDIF
    
                    CASE(HAM_SALSA)
    
                        !++alaak
                        DO ii=2,fn2b-fn2a !in2b has been defined already
                        idum = bc_define(TRIM(secname)//' emissions of '//TRIM(shortname)//CHAR(ii+48), bc_struc, &
                                         idims, .TRUE.)    
                        ENDDO
                        !--alaak
            
                        IF (TRIM(secname) == 'SEASALT') THEN
                           !For seasalt emissions BC is defined
                           !for number:
                           DO ii=1,fn2a-fn1a
                              idum = bc_define(TRIM(secname)//' emissions of '//TRIM(shortname)//'n'//CHAR(ii+48), &
                                               bc_struc, idims, .TRUE.) 
                           ENDDO 
                        ENDIF
            
                END SELECT
             ENDIF

          ENDIF

          !! set flags for special emissions
          IF ((TRIM(secname) == 'VOLCC') .AND. (bc_struc%ef_type /= EF_INACTIVE)) lhas_volcc = .TRUE.
          IF ((TRIM(secname) == 'VOLCE') .AND. (bc_struc%ef_type /= EF_INACTIVE)) lhas_volce = .TRUE.
          IF ((TRIM(secname) == 'OCEANI') .AND. (bc_struc%ef_type /= EF_INACTIVE)) lhas_oceani = .TRUE.
          IF ((TRIM(secname) == 'BIOGENIC') .AND. (bc_struc%ef_type /= EF_INACTIVE)) THEN
             lhas_biogenic = .TRUE.
             IF (bc_struc%ef_type == EF_MODULE) lbioemi_dyn = .TRUE. !SF replacement for former megan switch
                                                                      ! See IssueID #153
          ENDIF
        ELSE        ! itrtype /= itrprog --> no emissions
          WRITE(message_text,'(3a)') 'No emissions for species ',TRIM(shortname),' because itrtype/=itrprog.'
          CALL message('init_emissions',message_text,level=em_warn)
        END IF
!>>SF add more warning info here in case the current species is not defined
      ELSE
          WRITE(message_text,'(6a)') 'Species ',TRIM(shortname), &
                                ' (requested in the emissions matrix by sector ', &
                                TRIM(secname),') is not defined by any submodel! ', &
                                'You may check and revise your namelists settings!'   
          CALL message('init_emissions',TRIM(message_text),level=em_warn)
!<<SF
      END IF
    END DO
  END DO

  ! -- Initialize submodel specific emissions
  IF (lhas_volcc)  CALL init_emi_volcano(1)
  IF (lhas_volce)  CALL init_emi_volcano(2)
  IF (lhas_oceani) CALL init_emi_ocean
!++mgs 2015-02-24: bug fix
  IF (lhas_biogenic .AND. lbioemi_dyn) THEN
     ! start bioemi module here only if not standalone. In standalone mode it is started from init_subm
     ! in mo_submodel_interface
     IF (.NOT. lbioemi_stdalone) CALL start_biogenic_emissions
     CALL init_emi_biogenic
!--mgs
  ELSE
     IF (lbioemi_stdalone) THEN
        ! !!minor baustelle!! -- the following error should not occur, because lbioemi_standalone cannot
        ! run with either HAM or MOZ.
        WRITE(message_text,'(a)') 'Biogenic NMVOC emissions as standalone submodel is enabled '// &
                                  'but no BIOGENIC sector was found in the emissions matrix'
        CALL message('init_emissions',TRIM(message_text),level=em_error)
     ENDIF
  ENDIF

  IF (lham) THEN
     !++alaak
     !select emissions
     SELECT CASE (nham_subm)
         CASE(HAM_BULK)
             !CALL ham_bulk_init_emissions(nsectors)
         CASE(HAM_M7)
             CALL ham_m7_init_emissions(nsectors)
         CASE(HAM_SALSA)
             CALL ham_salsa_init_emissions(nsectors)
     END SELECT
     !--alaak
  ENDIF
  
  IF (.NOT. lemissions) &
    CALL message('init_emissions', 'lemissions=.TRUE., but no emissions defined! Now set to false.', &
                 level=em_warn)

  END SUBROUTINE init_emissions

  !! ---------------------------------------------------------------------------
  !! subroutine to define the emi diagnostic stream
  SUBROUTINE init_emi_stream

  USE mo_string_utls,         ONLY: st1_in_st2_proof
  USE mo_util_string,         ONLY: tolower
  USE mo_exception,           ONLY: finish
  USE mo_memory_base,         ONLY: new_stream, &
                                    default_stream_setting, &
                                    add_stream_reference, &
                                    AUTO
  USE mo_tracer,              ONLY: validate_traclist
  USE mo_tracdef,             ONLY: ln, ntrac, trlist, GAS, AEROSOL
  USE mo_species,             ONLY: nspec, speclist, ITRPROG
  USE mo_ham,                 ONLY: nclass
  USE mo_submodel,            ONLY: lham
  USE mo_submodel_streams,    ONLY: emi_lpost, emi_lpost_sector, emi_tinterval, &
                                    eminam, emi_gastrac, emi_keytype
  USE mo_submodel_diag,       ONLY: new_diag_list, new_diag, new_diag_element,  &
                                    BYTRACER, BYSPECIES
  USE mo_emi_matrix,          ONLY: em_get_sector_info
  USE mo_hammoz_emi_ocean,    ONLY: init_emi_ocean_stream
  USE mo_hammoz_emi_biogenic, ONLY: init_emi_biogenic_stream, lbioemi_dyn
  USE mo_ham_dust,            ONLY: bgc_dust_init_diag

   ! local variables
  INTEGER, PARAMETER             :: ndefault = 1
  CHARACTER(LEN=32)              :: defnam(1:ndefault)   = &   ! default output
                              (/ 'emi             ' /)         ! total emission mass flux

  CHARACTER(len=ln)              :: defaultgas(8)        = &   ! default gas-phase tracers for diagnostics
                              (/ 'SO2     ',               &
                                 'H2SO4   ',               &
                                 'DMS     ',               &
                                 'NO      ',               &
                                 'NO2     ',               &
                                 'CO      ',               &
                                 'C5H8    ',               &
                                 'CH4     '           /)
  LOGICAL                        :: tracflag(ntrac), specflag(nspec)
  CHARACTER(LEN=ln)              :: tracname(ntrac), specname(nspec)
  CHARACTER(LEN=32)              :: cunit 
  CHARACTER(LEN=64)              :: sectemp, secname, cdiagname
  CHARACTER(LEN=256)             :: cdiaglongname
  TYPE (t_stream), POINTER       :: emi_stream
  INTEGER                        :: ierr, jt, jsec, jspec, nvars
  LOGICAL                        :: lpost, ldiagdetail(maxsectors)

  ! default values and namelist read are done in init_submodel_streams !
  ldiagdetail(:) = .TRUE.      ! submodels may switch off detailed diagnostics for specific sectors

  !-- handle ALL and DEFAULT options for emi output variables
  IF (TRIM(tolower(eminam(1))) == 'all')     eminam(1:nemivars) = emivars(:)
  IF (TRIM(tolower(eminam(1))) == 'default') eminam(1:ndefault) = defnam(:)

  !-- check that all variable names from namelist are valid
  IF (.NOT. st1_in_st2_proof( eminam, emivars, ierr=ierr) ) THEN
    IF (ierr > 0) CALL finish ( 'ini_emi_stream', 'variable '// &
                                eminam(ierr)//' does not exist in emi stream' )
  END IF

  !-- find out which gas-phase tracers shall be included in diagnostics
  CALL validate_traclist(emi_gastrac, defaultgas, nphase=GAS,              &
                         ltran=.true., lemis=.true.)

  !-- define the flags and names for the diagnostic lists. We need one set of flags and
  !   names for each key_type (BYTRACER, BYSPECIES, BYMODE)
  !   gas-phase tracers will always be defined BYTRACER, for aerosol tracers one of the
  !   following lists will be empty.
  !   Note: vddep uses BYTRACER or BYMODE, ddep uses BYTRACER or BYSPECIES
  tracflag(:) = .FALSE.
  DO jt = 1,ntrac
    tracname(jt) = trlist%ti(jt)%fullname
    IF (trlist%ti(jt)%nphase == GAS) THEN
      tracflag(jt) = st1_in_st2_proof(trlist%ti(jt)%fullname, emi_gastrac)
    ELSE     ! aerosol tracer
      IF (emi_keytype == BYTRACER .AND. nclass > 0) THEN
        tracflag(jt) = trlist%ti(jt)%nemis > 0
      END IF
    END IF
  END DO
  specflag(:) = .FALSE.
  DO jt = 1,nspec
    specname(jt) = speclist(jt)%shortname
    IF (emi_keytype == BYSPECIES .AND.                        &
        IAND(speclist(jt)%nphase, AEROSOL) /= 0 .AND.         &
        speclist(jt)%itrtype == ITRPROG .AND.                 &
        nclass > 0) THEN
      specflag(jt) = speclist(jt)%lemis
    END IF
  END DO

  !-- open new stream
  CALL new_stream (emi_stream,'emi',lpost=emi_lpost,lrerun=.FALSE., &
                   interval=emi_tinterval)
  CALL default_stream_setting (emi_stream, lrerun = .FALSE., &
                   contnorest = .TRUE., table = 199, &
                   laccu = .false., code = AUTO)
   
  !-- add standard ECHAM variables
  IF (emi_lpost) THEN
    CALL add_stream_reference (emi_stream, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (emi_stream, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (emi_stream, 'aps'     ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (emi_stream, 'gboxarea','geoloc',lpost=.TRUE.)
  END IF

  !-- add emission mass flux diagnostics (accumulated)
  CALL default_stream_setting (emi_stream, lrerun=.FALSE., laccu=.TRUE.)

  lpost = st1_in_st2_proof( 'emi', eminam) .AND. emi_lpost
  CALL new_diag_list (emi_diag, emi_stream, diagname='emi', tsubmname='',    &
                      longname='accumulated emission mass flux',             &
                      units='kg m-2 s-1', ndims=2,                           &
                      nmaxkey=(/ntrac, nspec, 0, 0, 0 /), lpost=lpost )
  ! add diagnostic elements only when output is activated
  IF (lpost .AND. (ANY(tracflag) .OR. ANY(specflag))) THEN
    CALL new_diag(emi_diag, ntrac, tracflag, tracname, BYTRACER)
    CALL new_diag(emi_diag, nspec, specflag, specname, BYSPECIES)
  END IF

  !-- add submodel specific diagnostics
! sschr: this should now be obsolete
! IF (lham) CALL ham_m7_init_emi_stream(nsectors, emi_stream, ldiagdetail(1:nsectors))
! ustar_acrit should be added to the diagnostic stream
  IF (lham) CALL bgc_dust_init_diag(emi_stream)

  !-- add special output for interactive ocean emissions
  IF (lhas_oceani) CALL init_emi_ocean_stream(nsectors, emi_stream, ldiagdetail(1:nsectors))

  !-- add special output for interactive vegetation emissions
  IF (lbioemi_dyn) CALL init_emi_biogenic_stream !SF replaced the condition here: this call is relevant
                                                 ! ONLY in case of dynamic vegetation emissions

  !-- add detailed emission mass flux diagnostics (per sector)
  ! for now only BYSPECIES
  IF (emi_lpost_sector) THEN 
    CALL default_stream_setting (emi_stream, lrerun = .FALSE., &
                                 table = 199, &
                                 laccu = .true., code = AUTO)
  
    DO jsec = 1, nsectors
      CALL em_get_sector_info(jsec, secname, nvars)
      IF (nvars > 0) THEN
        sectemp=tolower(secname) !HK #529 circumvent Cray bug
        CALL new_diag_list (emi_diag_detail(jsec), emi_stream, diagname='emi_'//TRIM(sectemp), &
                            longname='accumulated '//TRIM(sectemp)//' emission mass flux',     &
                            units='kg m-2 s-1',                                                         &
                            tsubmname='', ndims=2, nmaxkey=(/0, nspec, 0, 0, 0 /) )
        DO jspec=1,nspec
          IF (ldiagdetail(jsec) .AND. emimod(jsec, jspec) /= EM_NONE) THEN
            cdiagname = 'emi_'//TRIM(speclist(jspec)%shortname)//'_'//TRIM(tolower(secname))
            cdiaglongname = 'accumulated emission mass flux of '//TRIM(speclist(jspec)%shortname) &
                            //' due to '//TRIM(tolower(secname))
            cunit = 'kg m-2 s-1'
            CALL new_diag_element(emi_diag_detail(jsec), cdiagname,                 &
                                  BYSPECIES, jspec,                                 &
                                  longname=cdiaglongname,                           &
                                  units=cunit                  )
          END IF
        END DO
      END IF
    END DO
  END IF

  END SUBROUTINE init_emi_stream

!gf #161  SUBROUTINE get_emi_field(ibc_emis,emimod,kproma,kbdim,krow,klev,lo3d,pbc2d,pbc3d)
  SUBROUTINE get_emi_field(ibc_emis,emimod,kproma,kbdim,krow,klev,ihpbl,lfirst,lo3d,pbc2d,pbc3d)
  USE mo_boundary_condition,    ONLY: bc_apply, bc_modify, bc_query, &
                                      BC_VERTICAL_WEIGHTED_INTERPOLATION
  USE mo_hammoz_emi_fire,       ONLY: distribute_emi_fire
  USE mo_exception,             ONLY: finish, message_text


  INTEGER,  INTENT(in)    :: ibc_emis                 ! index of boundary_condition to store to emission field
  INTEGER,  INTENT(in)    :: emimod                   ! mode of emission
  INTEGER,  INTENT(in)    :: kproma                   ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                    ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: krow                     ! geographic block number
  INTEGER,  INTENT(in)    :: klev                     ! number of levels
  INTEGER,  INTENT(in)    :: ihpbl(kbdim)             ! level of PBL top
  LOGICAL,  INTENT(in)    :: lfirst                   ! flag to indicate whether it is the first timestep
  LOGICAL,  INTENT(inout) :: lo3d                     ! flag to indicate whether the emissions are applied in 3D
  REAL(dp), INTENT(inout) :: pbc2d(kbdim)             ! for emission field
  REAL(dp), INTENT(inout) :: pbc3d(kbdim, klev)       ! for emission field
  CHARACTER*30            :: strunit
  CHARACTER(LEN=128)      :: bc_name

  CALL bc_query(ibc_emis, name=bc_name, ef_actual_unit=strunit)
  SELECT CASE (emimod)
  CASE (EM_SURFACE)
    IF ((index(strunit,'/m2') > 0) .OR. (index(strunit,'m-2') > 0) &
         .OR. (index(strunit,'/m**2') > 0) .OR. (index(strunit,'m**-2') > 0)) THEN
      CALL bc_apply (ibc_emis,kproma,krow,pbc2d)
      ! don't convert unit [kg m-2 s-1] --> it is needed like this
      lo3d = .FALSE.
    ELSE
      WRITE(message_text,'(3a)') 'Wrong unit for ',TRIM(bc_name), ' of type EM_SURFACE: '//strunit
      CALL finish('get_emi_field', message_text)
    ENDIF
  CASE (EM_VOLUME)
    IF ((index(strunit,'/m3') > 0) .OR. (index(strunit,'m-3') > 0) &
         .OR. (index(strunit,'/m**3') > 0) .OR. (index(strunit,'m**-3') > 0)) THEN
      IF (lfirst) THEN
        CALL bc_modify(ibc_emis, bc_vertint=BC_VERTICAL_WEIGHTED_INTERPOLATION)
      ENDIF
      CALL bc_apply (ibc_emis,kproma,krow,pbc3d)
!     convert to kg kg-1 s-1
      pbc3d(1:kproma,:) = pbc3d(1:kproma,:) / dens(1:kproma,:)
    ELSEIF ((index(strunit,'/m2') > 0) .OR. (index(strunit,'m-2') > 0) &
         .OR. (index(strunit,'/m**2') > 0) .OR. (index(strunit,'m**-2') > 0)) THEN
      CALL bc_apply (ibc_emis,kproma,krow,pbc3d)
!     convert to kg kg-1 s-1
      pbc3d(1:kproma,:) = pbc3d(1:kproma,:) / dens(1:kproma,:) / height(1:kproma,:)
    ELSE
      WRITE(message_text,'(2a)') 'Wrong unit for emissions of type EM_VOLUME: '//strunit
      CALL finish('get_emi_field', message_text)
    ENDIF
  CASE (EM_LEVEL50M)
    IF ((index(strunit,'/m2') > 0) .OR. (index(strunit,'m-2') > 0) &
         .OR. (index(strunit,'/m**2') > 0) .OR. (index(strunit,'m**-2') > 0)) THEN
      CALL bc_apply (ibc_emis,kproma,krow,pbc2d)
      ! 50m means: put emissions in second lowest model level
      pbc3d(1:kproma,klev-1) = pbc2d(1:kproma)
!     convert to kg kg-1 s-1
      pbc3d(1:kproma,:) = pbc3d(1:kproma,:) / dens(1:kproma,:) / height(1:kproma,:)
      pbc2d(:) = 0._dp
    ELSE
      WRITE(message_text,'(2a)') 'Wrong unit for emissions of type EM_SURFACE: '//strunit
      CALL finish('get_emi_field', message_text)
    ENDIF
  CASE (EM_FIRE)
    IF ((index(strunit,'/m2') > 0) .OR. (index(strunit,'m-2') > 0) &
         .OR. (index(strunit,'/m**2') > 0) .OR. (index(strunit,'m**-2') > 0)) THEN
      CALL bc_apply (ibc_emis,kproma,krow,pbc2d)
!gf #161   CALL distribute_emi_fire(kproma, kbdim, klev, pbc2d, pbc3d)
      CALL distribute_emi_fire(kproma, kbdim, klev, krow, ihpbl, pbc2d, pbc3d)
      pbc2d(:) = 0._dp
!     convert to kg kg-1 s-1
      pbc3d(1:kproma,:) = pbc3d(1:kproma,:) / dens(1:kproma,:) / height(1:kproma,:)
    ELSE
      WRITE(message_text,'(2a)') 'Wrong unit for emissions of type EM_SURFACE: '//strunit
      CALL finish('get_emi_field', message_text)
    ENDIF
  END SELECT
  END SUBROUTINE get_emi_field

  !@brief: the emi_interface routine provides the interface for all submodel emissions
  ! 
  SUBROUTINE emi_interface(kproma, kbdim, ktrac, klev, klevp1, krow,    &
!gf #161                   paphp1, pgeom1, loland, ptm1, pxtems, pxtte)
                           paphp1, pgeom1, loland, ptm1, ihpbl, pxtems, pxtte)

  USE mo_physical_constants,    ONLY: grav
  USE mo_species,               ONLY: nspec, nmaxtrspec, spec_idt, spec_ntrac
  USE mo_time_control,          ONLY: delta_time
  USE mo_submodel,              ONLY: lham, lmoz
  USE mo_submodel_diag,         ONLY: get_diag_pointer
  USE mo_ham,                   ONLY: npist
  USE mo_ham_species,           ONLY: id_so2, id_so4
  USE mo_ham_m7_emissions,      ONLY: ham_m7_emissions
  USE mo_hammoz_emi_volcano,    ONLY: read_emi_volcano
  USE mo_hammoz_emi_ocean,      ONLY: oceani_emissions     ! interactive ocean emissions
  USE mo_hammoz_emi_biogenic,   ONLY: lbioemi_dyn, calc_biogenic_emissions  ! interactive terrestrial 
                                                                            ! vegetation emissions
  !++alaak
  USE mo_ham_salsa_emissions,   ONLY: ham_salsa_emissions
  USE mo_ham,                   ONLY: nham_subm,HAM_BULK,HAM_M7,HAM_SALSA
  !--alaak
  USE mo_vphysc,                ONLY: vphysc
  USE mo_geoloc,                ONLY: gboxarea
  USE mo_physical_constants,    ONLY: grav
! arguments

  INTEGER,  INTENT(in)    :: kproma                   ! geographic block number of locations
  INTEGER,  INTENT(in)    :: kbdim                    ! geographic block maximum number of locations
  INTEGER,  INTENT(in)    :: ktrac                    ! number of tracers
  INTEGER,  INTENT(in)    :: klev                     ! number of levels
  INTEGER,  INTENT(in)    :: klevp1                   ! number of levels + 1
  INTEGER,  INTENT(in)    :: krow                     ! geographic block number
  INTEGER,  INTENT(in)    :: ihpbl(kbdim)             ! level of PBL top
  REAL(dp), INTENT(in)    :: paphp1(kbdim,klev+1)     ! interface pressure levels at t+1
  REAL(dp), INTENT(in)    :: pgeom1(kbdim,klev)       ! geopotential height
  LOGICAL,  INTENT(in)    :: loland(kbdim)            ! land mask
  REAL(dp), INTENT(in)    :: ptm1(kbdim,klev)         ! temperature
  REAL(dp), INTENT(inout) :: pxtems(kbdim,ktrac)      ! surface emissions
  REAL(dp), INTENT(inout) :: pxtte (kbdim,klev, ktrac)! tracer tendencies

  INTEGER           :: jl, jk, jt, jsec, jspec, jtrac, ierr, i, ispec, ibc_extra(nmaxtrspec)
  LOGICAL           :: lo3d
  LOGICAL, SAVE     :: lfirst = .TRUE.
  REAL(dp)          :: zbc2d(kbdim), zbc3d(kbdim, klev)  ! for boundary conditions
  REAL(dp), POINTER :: fld2d(:,:)                        ! for diagnostics
  REAL(dp)          :: zfac, zfactor(nmaxtrspec)         ! weighting factor for aerosol modes
  REAL(dp)          :: zemidiag(kbdim)                   ! 2D field for diagnostics

! should probably be in argument list and handled by ECHAM similar to pxtems
! REAL(dp), INTENT(inout) :: zxtems3d(kbdim,klev,ktrac)
  REAL(dp) :: zxtems3d(kbdim,klev)

  !-- read extra files at first time step where necessary 
  IF (lfirst) THEN
    ! read continuous volcano emissions
    IF (lhas_volcc) CALL read_emi_volcano(kproma, kbdim, klev, 1)
    ! read explosive volcano emissions
    IF (lhas_volce) CALL read_emi_volcano(kproma, kbdim, klev, 2)
  END IF

  !-- Conversion factor for 2D mass flux to 3D mass mixing ratio change
  DO jk=1,klev
    DO jl=1,kproma
  !-- model levels are topdown (highest level near surface)
      dgeo(jl,jk) = vphysc%geohm1(jl,jk,krow)-vphysc%geohm1(jl,jk+1,krow)
      height(jl,jk) = dgeo(jl,jk) /grav
      dens(jl,jk) = (paphp1(jl,jk+1)-paphp1(jl,jk)) / dgeo(jl,jk)
    END DO
  END DO

  ! try to keep statements here at minimum
  ! store mass fluxes from interactive ocean emissions in boundary condition
  IF (lhas_oceani) CALL oceani_emissions(kproma, kbdim, krow, npist)
  IF (lbioemi_dyn) CALL calc_biogenic_emissions(kproma, kbdim, krow, loland, ptm1(:,klev))
  
!-- Sector loop: Prepare species-independent stuff that may be sector-dependent
  DO jsec = 1, nsectors

!!!      CALL get_emi_time_factor(jsec, ...) !weekly cycle f.ex.

    !-- Species loop: obtain 2D or 3D emission flux fields
    !   and distribute vertically where needed
    DO jspec = 1, nspec

      IF (emimod(jsec,jspec) == EM_NONE) CYCLE     ! do nothing

      !-- Calculate emissions of type EF_MODULE (e.g. seasalt, dust, megan)
      !   and store them to the bc scheme

      !-- Also distribute species emissions among tracers
      !   This is trivial for gas-phase tracers but needs some attention for aerosols
      !   where species emissions are distributed among aerosol modes or bins
      !   Therefore we provide an explicit interface to various HAM emissions here

      ! set default weighting factors as equal distribution across all tracers
      ! (attention: This is NO equal distribution for HAM_M7 as ham_m7_init_emissions
      ! added NUM-fields to spec_ntrac and spec_idt; zfactor for HAM_M7 emissions
      ! are set via ham_m7_emissions!)
      zfac = 1._dp/spec_ntrac(jspec)
      zfactor(:) = zfac
      ibc_extra(:) = 0

      ! Handle special cases
      IF (lham) THEN
         !++alaak
         !select emissions
         SELECT CASE (nham_subm)
             CASE(HAM_BULK)
                 !CALL ham_bulk_init_emissions(nsectors)
             CASE(HAM_M7)
                 CALL ham_m7_emissions(kproma, kbdim, klev, krow, ktrac, jsec, jspec, &
                                       ibc_extra, zfactor)
             CASE(HAM_SALSA)
                 CALL ham_salsa_emissions(kproma, kbdim, klev, krow, ktrac, jsec, jspec, &
                                          ibc_extra, zfactor)
     
         END SELECT
         !--alaak
      ENDIF

      IF (lmoz) THEN 
        !! some settings for lightning NOx (??)
      END IF

      zbc2d(:) = 0._dp
      zbc3d(:,:) = 0._dp
      lo3d = .TRUE.        ! process as 3D emissions (default = true, because we don't know)

      ispec = jspec
      IF (ispec == id_so4) ispec = id_so2
!gf #161  CALL get_emi_field(ibc_emis(jsec,ispec),emimod(jsec,ispec),kproma,kbdim,krow,klev,lo3d,zbc2d,zbc3d)
      CALL get_emi_field(ibc_emis(jsec,ispec),emimod(jsec,ispec),kproma,kbdim,krow,klev,ihpbl,lfirst,lo3d,zbc2d,zbc3d)
      ! we should now have the total species emissions
      ! in either zbc2d or zbc3d. Next task is to distribute them among all tracers belonging
      ! to this species. Effectively this means copying zbc2d in pxtems and zbc3d in zxtems3d.
      ! zbc2d unit: [kg m-2 s-1]
      ! zbc3d unit: [kg kg-1]

      ! Default handling of emissions
      DO jtrac = 1, spec_ntrac(jspec)
!gf #161   IF (ibc_extra(jtrac) /= 0) CALL get_emi_field(ibc_extra(jtrac),emimod(jsec,ispec),kproma,kbdim,krow,klev,lo3d,zbc2d,zbc3d)
        IF (ibc_extra(jtrac) /= 0) &
           CALL get_emi_field(ibc_extra(jtrac),emimod(jsec,ispec),kproma,kbdim,krow,klev,ihpbl,lfirst,lo3d,zbc2d,zbc3d)
        jt = spec_idt(jspec, jtrac)
        IF (zfactor(jtrac) > 0._dp) THEN
          ! initialize zemidiag: will be zero if zbc2d is zero
          zemidiag(1:kproma) = zfactor(jtrac)*zbc2d(1:kproma)*emfactor(jsec, jspec)
          pxtems(1:kproma, jt) = pxtems(1:kproma, jt) + zemidiag(1:kproma)
          IF (lo3d) THEN
            zxtems3d(1:kproma,:) = zfactor(jtrac) * zbc3d(1:kproma,:) * emfactor(jsec, jspec)
            pxtte(1:kproma,:,jt) = pxtte(1:kproma,:,jt) + zxtems3d(1:kproma,:)
            DO i = 1, klev
              zemidiag(1:kproma) = zemidiag(1:kproma) + zxtems3d(1:kproma,i) * dens(1:kproma,i) * height(1:kproma,i)
            ENDDO
          END IF
          ! total mass flux
          CALL get_diag_pointer(emi_diag, fld2d, jt, ierr=ierr)
          IF (ierr == 0) fld2d(1:kproma,krow)=fld2d(1:kproma,krow)          &
                            +zemidiag(1:kproma)*delta_time
          ! detailed
          CALL get_diag_pointer(emi_diag_detail(jsec), fld2d, jt, ierr=ierr)
          IF (ierr == 0) fld2d(1:kproma,krow)=fld2d(1:kproma,krow)          &
                            +zemidiag(1:kproma)*delta_time
        END IF
      END DO   ! tracer loop

    END DO     ! species loop

  END DO       ! sector loop

  lfirst = .FALSE.

  END SUBROUTINE emi_interface

END MODULE mo_hammoz_emissions
