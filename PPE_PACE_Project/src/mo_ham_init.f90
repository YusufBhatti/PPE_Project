!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_init.f90
!!
!! \brief
!! Initialisation routines for HAM aerosol model
!!
!! \author Martin Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (2009-09-22)
!!   -# D. O'Donnell, ETH-Z  - added SOA support (2010-04-22)
!!
!! \limitations
!! None
!!
!! \details
!! This module contains the general HAM initialisation routine which is called 
!! from submodel_initialise.
!! It first calls setham to initialize ham module variables and read the hamctl namelist.
!! Then other namelists are read and the M7 module is started. Next, the HAM species are
!! defined and the aerosol modes are populated.
!! Adopted from former init_ham routine. Contributors to
!! the code from which the present routines were derived include:
!! P. Stier, D. O'Donnell, K. Zhang
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

MODULE mo_ham_init

  IMPLICIT NONE

  PRIVATE

  SAVE

  PUBLIC  :: start_ham,            &
             ham_define_tracer,    &
             ham_initialize,       &
             ham_init_memory,      &
             ham_free_memory

! Variable declarations


  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! start_ham: high-level initialisation routine; interface for HAM aerosol module
!! initialisation including species definition
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! init_submodels
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE start_ham
  
    USE mo_ham,             ONLY : setham,          &
                                   naerosol,        &
                                   nsoa,            &
                                   nham_subm,       &
                                   sizeclass,       &
                                   nclass,          &
                                   HAM_BULK,        &
                                   HAM_M7,          &
                                   HAM_SALSA,       &
                                   subm_naerospec
    USE mo_ham_dust,        ONLY : setdust
    USE mo_ham_m7ctl,       ONLY : sethamM7, m7_initialize
    USE mo_ham_salsactl,    ONLY : setham_salsa  
    USE mo_ham_salsa_init,  ONLY : salsa_initialize
    USE mo_ham_subm_species,ONLY : map_ham_subm_species
    USE mo_ham_species,     ONLY : ham_species
    USE mo_ham_soa,         ONLY : soa_species, start_soa_aerosols
    ! >> thk: volatility basis set (VBS)
    USE mo_ham_vbsctl, ONLY : setham_vbs, laqsoa
    USE mo_ham_vbs,    ONLY : vbs_species
    ! << thk
  
    !-- local variables
    LOGICAL    :: lsoainclass(nclass)    ! copy from sizeclass(:)%lsoainclass for code structure reasons
    INTEGER    :: jm
  
    ! -- 1. set default values and read hamctl namelist

    CALL setham

    ! -- read ham_dustctl namelist
    CALL setdust

    SELECT CASE(nham_subm)

       ! Initialization for bulk microphysics

    CASE(HAM_BULK)

       ! CALL sethambulk
       
       ! -- initialize bulk scheme
       ! CALL bulk_initialize
       
    CASE(HAM_M7)
       
       
       ! -- initialize M7 scheme
       CALL sethamM7
       CALL m7_initialize
       
    CASE(HAM_SALSA)
       
       ! -- initialize SALSA scheme
       CALL setham_salsa
       CALL salsa_initialize
       
    END SELECT
  
    ! -- 3. define aerosol species
    CALL ham_species

    SELECT CASE (nsoa)
       CASE(1)
          CALL soa_species
       CASE(2)
          CALL setham_vbs
          ! this had to be moved here from mo_ham_species,
          ! because it caused a circular dependence
          CALL vbs_species
    END SELECT

    SELECT CASE (nham_subm)
       
    CASE(HAM_BULK)
       
       ! CALL map_bulk_species
       
    CASE(HAM_M7)
       
       ! -- map general species list onto M7 condensed list
       CALL map_ham_subm_species
  
       ! -- 4. generate mode x species matrix
       CALL ham_define_modes(nclass, naerosol(nham_subm), subm_naerospec)
       IF (nsoa == 1) THEN
          DO jm=1,nclass
             lsoainclass(jm) = sizeclass(jm)%lsoainclass
          END DO
          CALL start_soa_aerosols(nclass, lsoainclass)
       END IF

    CASE(HAM_SALSA)

       ! CALL map_salsa_species
       CALL map_ham_subm_species

       CALL ham_define_bins(nclass, naerosol(nham_subm), subm_naerospec)
    END SELECT
  
  
  END SUBROUTINE start_ham

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! ham_define_tracer: create ECHAM tracers based on HAM species
!!
!! @author see module info
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info
!!
!! @par This subroutine is called by
!! init_submodels
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Notes
!! ### ToDo: complete allocation of burden diagnostics
!!
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ham_define_tracer
 
    ! 
    ! uses local grib table 131 (tracer)
    ! codes >= 160 : gas-phase species
    ! codes >= 201 : aerosol species
    !
   
    USE mo_exception,         ONLY: message, message_text, em_error
    USE mo_species,           ONLY: nspec, speclist
    USE mo_tracdef,           ONLY: OK, RESTART, CONSTANT, INITIAL,      &
                                    ON, OFF, itrprog, itrdiag, itrpresc, &
                                    itrnone, & !gf #57
                                    GAS, GAS_OR_AEROSOL, AEROSOLNUMBER, AEROSOLMASS, &
                                    SOLUBLE, INSOLUBLE, t_flag, ln, ll
    USE mo_tracer,            ONLY: new_tracer, new_diag_burden
    USE mo_submodel_diag,     ONLY: BYTRACER, BYSPECIES
    USE mo_physical_constants, ONLY: rhoh2o
    USE mo_advection,         ONLY: iadvec
    USE mo_ham,               ONLY: naerocomp, aerocomp, aerowater,      &
                                    naerosol,nham_subm,&
                                    ndrydep, nwetdep, nsoa,              &
                                    burden_keytype, sizeclass,nclass
    USE mo_ham,               ONLY: nsol
    USE mo_ham_m7_trac,       ONLY: idt_cdnc_ham, idt_icnc_ham
    USE mo_ham_soa,           ONLY: soaprop, set_soa_tracer_attr
    
    IMPLICIT NONE

    INTEGER :: jn, jclass, igribc, idt, ierr, jspid
    INTEGER :: init, iwrite, iburdenid, idrydep, iwetdep, itran, iconv, ivdiff, iint, isol, isedi

    CHARACTER(LEN=ll) :: trac_lname
    CHARACTER(LEN=ln) :: trac_sname
    CHARACTER(LEN=ln) :: trac_fullname
    CHARACTER(LEN=ln) :: str1
    CHARACTER(LEN=ln) :: csubname

    !>>dod fix for missing gas burdens
    CHARACTER(LEN=ln) :: burden_name
    !<<dod

    !---executable procedure
    
    !-------- 
    ! gases 
    !-------- 

    igribc=160      ! starting grib code for HAM gas species

    DO jn = 1,nspec

      ! look for HAM species
      IF (TRIM(speclist(jn)%tsubmname) /= 'HAM' .OR. IAND(speclist(jn)%nphase, GAS) == 0) CYCLE

      IF (speclist(jn)%itrtype == itrnone) CYCLE !gf #57

      ! units 
      SELECT CASE(TRIM(ADJUSTL(speclist(jn)%units)) )
      CASE('kg kg-1') 
         str1 = ' mass mixing ratio - gas phase'
      CASE('mole mole-1')
         str1 = ' volume mixing ratio'
      CASE DEFAULT
         str1=''
      END SELECT

      trac_sname = TRIM(ADJUSTL(speclist(jn)%shortname))
      trac_lname = TRIM(ADJUSTL(speclist(jn)%longname))//str1

      idrydep = OFF
      IF (speclist(jn)%ldrydep) idrydep = ndrydep   ! set drydep scheme to hamctl default
      iwetdep = OFF
      IF (speclist(jn)%lwetdep) iwetdep = nwetdep   ! set wetdep scheme to hamctl default      

      !>>dod SOA has unusual tracer handling, let the SOA module set the tracer attributes
      IF (nsoa == 1 .AND. ( (ANY(soaprop(:)%spid_tot == jn)) .OR. &
                            (ANY(soaprop(:)%spid_soa == jn)) ) ) THEN
         CALL set_soa_tracer_attr(jn, init, itran, iwrite, iconv, ivdiff, iint, iwetdep, idrydep)
      ELSE
        SELECT CASE(speclist(jn)%itrtype)
       
        CASE(itrprog)
           init     = RESTART+CONSTANT+INITIAL    ! Note: File(s) containing initial values for tracers must be in MMR!!!
           iwrite   = ON 
           itran    = iadvec  ! ++mgs ### was ON
           iconv    = ON
           ivdiff   = ON
           iint     = ON

        CASE(itrdiag)
           init     = -1
           iwrite   = ON 
           itran    = OFF
           iconv    = ON
           ivdiff   = ON
           iint     = OFF

        CASE(itrpresc)
           init     = -1
           iwrite   = OFF
           itran    = OFF
           iconv    = OFF
           ivdiff   = OFF
           iint     = OFF
           idrydep  = OFF
           iwetdep  = OFF

        END SELECT
      END IF
          

      csubname = ''
      !>>dod fix for missing gas burdens
      IF (speclist(jn)%nphase == GAS_OR_AEROSOL) THEN              ! if species exists in aerosol and gas phase
         csubname = 'gas'  
         burden_name=TRIM(trac_sname)//'_gas'
      ELSE
         burden_name=TRIM(trac_sname)
      END IF
      !<<dod

      IF (speclist(jn)%lburden .AND.   &
          (speclist(jn)%itrtype==itrprog .OR. speclist(jn)%itrtype==itrdiag) ) THEN
        iburdenid = new_diag_burden(burden_name, itype=1, lclobber=.false.)
      ELSE
        iburdenid = 0 
      END IF

      CALL new_tracer(trac_sname,  'HAM', ierr=ierr,                    &
                      spid          = jn,                               &
                      subname       = csubname,                         &
                      longname      = trac_lname,                       &
                      units         = speclist(jn)%units,               &
                      nphase        = GAS,                              &
                      ninit         = init,                             &
                      nwrite        = iwrite,                           &
                      code          = igribc,                           &
                      table         = 131,                              &   
                      burdenid      = iburdenid,                        &
                      ntran         = itran,                            &
                      nconv         = iconv,                            &
               !!     nconvmassfix  = ON,                               &
                      nvdiff        = ivdiff,                           &
                      nint          = iint,                             &
                      ndrydep       = idrydep,                          &
                      nwetdep       = iwetdep,                          &
                      moleweight    = speclist(jn)%moleweight,          &
                      idx           = idt                                  )

      IF (ierr == OK) THEN
        speclist(jn)%idt = idt
      ELSE
        WRITE(message_text,'(a,i0)') 'new_tracer '//TRIM(trac_sname)//' returned error code ',ierr
        CALL message('ham_define_tracer', message_text, level=em_error)
      END IF

      igribc = igribc + 1

    END DO
    

    !---------- 
    ! aerosols  
    !----------
    
    ! grib codes 201 .. 255 aerosol model

    !--- 2) Allocate aerosol masses according to the modes to obtain succeding tracer identifiers:

    igribc = 201

    DO jn=1,naerocomp

      ! sedimentation type 
      jclass = aerocomp(jn)%iclass
      jspid  = aerocomp(jn)%spid

      IF (sizeclass(jclass)%lsed) THEN
        isedi = ON
      ELSE
        isedi = OFF
      END IF

      ! solubility 
      IF (sizeclass(jclass)%lsoluble) THEN
        isol = SOLUBLE
      ELSE
        isol = INSOLUBLE
      END IF

      trac_sname = TRIM(ADJUSTL(aerocomp(jn)%species%shortname))
      trac_lname = TRIM(ADJUSTL(aerocomp(jn)%species%longname))&
                    &//'mass mixing ratio - mode '//sizeclass(jclass)%shortname    
      trac_fullname = TRIM(trac_sname) // '_' // TRIM(sizeclass(jclass)%shortname)

      idrydep = OFF
      IF (aerocomp(jn)%species%ldrydep) idrydep = ndrydep   ! set drydep scheme to hamctl default
      iwetdep = OFF
      IF (aerocomp(jn)%species%lwetdep) iwetdep = nwetdep   ! set wetdep scheme to hamctl default      

      !>>dod SOA has unusual tracer handling, let the SOA module set the tracer attributes
      !>>dod bugfix (redmine #59 and #60)
      IF (nsoa == 1 .AND.  (ANY(soaprop(:)%spid_soa == jspid) ) )  THEN
         CALL set_soa_tracer_attr(jspid, init, itran, iwrite, iconv, ivdiff, iint, iwetdep, idrydep)
      !<<dod   
      ELSE

        SELECT CASE(aerocomp(jn)%species%itrtype)
       
        CASE(itrprog)
          init     = RESTART+CONSTANT+INITIAL    ! Note: File(s) containing initial values for tracers must be in MMR!!!
          itran    = iadvec  ! ++mgs ### was ON
          iconv    = ON
          ivdiff   = ON
          iint     = ON

        CASE(itrdiag)
          init     = -1
          itran    = OFF
          iconv    = ON
          ivdiff   = ON
          iint     = OFF

        CASE(itrpresc)
          init     = -1
          itran    = OFF
          iconv    = OFF
          ivdiff   = OFF
          iint     = OFF
          idrydep  = OFF
          iwetdep  = OFF

        END SELECT
      END IF

!### not sure about this logic: please test!!
      iburdenid = 0 
      IF (aerocomp(jn)%species%lburden         &
          .AND. (aerocomp(jn)%species%itrtype==itrprog   &
                 .OR. aerocomp(jn)%species%itrtype==itrdiag) ) THEN
        IF (burden_keytype == BYTRACER) THEN       ! burden per species and per mode (i.e. per tracer)
          iburdenid = new_diag_burden(trac_fullname, itype=1, lclobber=.false.)
        ELSE IF (burden_keytype == BYSPECIES) THEN  ! burden per species
          iburdenid = new_diag_burden(aerocomp(jn)%species%shortname, itype=1, lclobber=.true.)
        END IF
      END IF

      CALL new_tracer(TRIM(ADJUSTL(trac_sname)),  'HAM', ierr=ierr,             &
                      spid        = aerocomp(jn)%spid,                          &
                      subname     = sizeclass(jclass)%shortname,                     & 
                      units       = aerocomp(jn)%species%units,                 &
                      nphase      = AEROSOLMASS,                                &
                      mode        = jclass,                                       &
                      ninit       = init,                                       &
                      nwrite      = ON,                                         &
                      code        = igribc,                                     &
                      table       = 131,                                        &
                      burdenid    = iburdenid,                                  &
                      ntran       = itran,                                      &
                      nvdiff      = ivdiff,                                     &
                      nconv       = iconv,                                      &
               !!     nconvmassfix  = ON,                                       &
                      nint        = iint,                                       &
                      nsoluble    = isol,                                       &
                      ndrydep     = idrydep,                                    &
                      nwetdep     = iwetdep,                                    &
                      nsedi       = isedi,                                      &
                      moleweight  = aerocomp(jn)%species%moleweight,            &
                      longname    = TRIM(ADJUSTL(trac_lname)),                  & 
                      idx         = idt                                           )

      IF (ierr == OK) THEN
        aerocomp(jn)%idt = idt
      ELSE
        WRITE(message_text,'(a,i0)') 'new_tracer '//TRIM(trac_sname)//' returned error code ',ierr
        CALL message('ham_define_tracer', message_text, level=em_error)
      END IF
    
      igribc = igribc+1

    END DO

    ! Note: additional tracer for transport of SOA species are defined in the SOA module  (dod)

    !--- 3) Special tracers for aerosol numbers and aerosol water

    DO jn=1,nclass
    
      IF (sizeclass(jn)%lsed) THEN
        isedi = ON
      ELSE
        isedi = OFF
      END IF

      IF (sizeclass(jn)%lsoluble) THEN
        isol = SOLUBLE
      ELSE
        isol = INSOLUBLE
      END IF

      idrydep = ndrydep    ! aerosol number is prognostic tracer: use hamctl flags
      iwetdep = nwetdep    ! to choose drydep and wetdep scheme as for progn. mass tracers

      !>>gf fix for missing aerosol number burden (#218 Redmine)
      IF (burden_keytype == BYTRACER) THEN       ! burden per species and per mode (i.e. per tracer)
         iburdenid = new_diag_burden('NUM_'//sizeclass(jn)%shortname, itype=1, lclobber=.false.)
      ELSE IF (burden_keytype == BYSPECIES) THEN  ! burden per species
         iburdenid = -1 !SF #299 this diag make no sense. Just set iburdenid to a negative value here
                        !        so that no burden is computed in this case
      END IF
      !<<gf

      CALL new_tracer('NUM',   'HAM', ierr=ierr,      &
                      subname=sizeclass(jn)%shortname,&
                      units='1 kg-1',                 &
                      table=131,                      &
                      code=igribc,                    &
                      nwrite=ON,                      &
                      burdenid=iburdenid,             &
                      nphase=AEROSOLNUMBER,           &
                      mode=jn,                        &
                      ndrydep=idrydep,                &
                      nwetdep=iwetdep,                &
                      nsoluble=isol,                  &
                      nsedi=isedi,                    &
                      nconv=ON,                       &
               !!     nconvmassfix=ON,                &
                      nvdiff=ON,                      &
                      ntran=iadvec,                   &
                      longname='number mixing ratio - aerosol mode ' &
                               //TRIM(ADJUSTL(sizeclass(jn)%classname)), &
                      idx=idt                            )

      IF (ierr == OK) THEN
        sizeclass(jn)%idt_no = idt
      ELSE
        WRITE(message_text,'(a,i0)') 'new_tracer '//TRIM(sizeclass(jn)%shortname)//' returned error code ',ierr
        CALL message('ham_define_tracer', message_text, level=em_error)
      END IF
    
      igribc = igribc + 1

    END DO

    DO jn=1,naerosol(nham_subm)

      jclass = aerowater(jn)%iclass

      !>>dod fix for missing aerosol water burden
      IF (burden_keytype == BYTRACER) THEN       ! burden per species and per mode (i.e. per tracer)
         iburdenid = new_diag_burden('WAT_'//sizeclass(jn)%shortname, itype=1, lclobber=.false.)
      ELSE IF (burden_keytype == BYSPECIES) THEN  ! burden per species
         iburdenid = new_diag_burden('WAT', itype=1, lclobber=.true.)
      END IF
      !<<dod

      CALL new_tracer('WAT',   'HAM', ierr=ierr,           &
                      subname=sizeclass(jn)%shortname,        & 
                      spid=aerowater(jn)%spid,             &
                      units=aerowater(jn)%species%units,   &
                      mode=sizeclass(jn)%self,                &
                      table=131,                           &
                      code=igribc,                         &
                      nwrite=ON,                           &
                      burdenid=iburdenid,                  &
                      ndrydep=OFF,                         &
                      nwetdep=OFF,                         &
                      nsedi=OFF,                           &
                      nphase=AEROSOLMASS,                  &
                      nconv=OFF,                           &
               !!     nconvmassfix=OFF,                    &
                      nvdiff=OFF,                          &
                      ntran=OFF,                           &
                      nint=ON,                             &
                      myflag=(/t_flag('density',rhoh2o)/), &
                      longname='Aerosol water mass mixing ratio - mode ' &
                                   //TRIM(ADJUSTL(sizeclass(jn)%classname)),     &
                      idx=idt                                 )

      IF (ierr == OK) THEN
        aerowater(jn)%idt = idt
      ELSE
        WRITE(message_text,'(a,i0)') 'new_tracer '//TRIM(sizeclass(jn)%shortname)//' returned error code ',ierr
        CALL message('ham_define_tracer', message_text, level=em_error)
      END IF
    
      igribc = igribc+1
  
    END DO

    !--- 5) Cloud properties: 
    
    iburdenid = 0 
    IF (burden_keytype > 0) THEN
      iburdenid = new_diag_burden('CDNC', itype=1, lclobber=.false.)
    END IF
    CALL new_tracer('CDNC', 'HAM', ierr=ierr,      &
                     idx=idt_cdnc_ham,             &
                     units='1 kg-1',               &
                     table=131,                    &
                     code=131,                     &
                     nwrite=ON,                    &
                     burdenid=iburdenid,           &
                     nconv=OFF,                    &
              !!     nconvmassfix=OFF,             &
                     nvdiff=ON,                    &
                     nint=ON,                      &
                     longname='cloud droplet number concentration')

    IF (ierr /= OK) THEN
      WRITE(message_text,'(a,i0)') 'new_tracer CDNC returned error code ',ierr
      CALL message('ham_define_tracer', message_text, level=em_error)
    END IF
    
    iburdenid = 0 
    IF (burden_keytype > 0) THEN
      iburdenid = new_diag_burden('ICNC', itype=1, lclobber=.false.)
    END IF
    CALL new_tracer('ICNC', 'HAM', ierr=ierr,      &
                     idx=idt_icnc_ham,             &
                     units='1 kg-1',               &
                     table=131,                    &
                     code=132,                     &
                     nwrite=ON,                    &
                     burdenid=iburdenid,           &
                     nconv=OFF,                    &
              !!     nconvmassfix=OFF,             &
                     nvdiff=ON,                    &
                     nint=ON,                      &
                     longname='ice crystal number concentration')

    IF (ierr /= OK) THEN
      WRITE(message_text,'(a,i0)') 'new_tracer ICNC returned error code ',ierr
      CALL message('ham_define_tracer', message_text, level=em_error)
    END IF

  END SUBROUTINE ham_define_tracer

  ! -- utility function to allow definition of tracer with same name for gas and aero species
  LOGICAL FUNCTION lsubname_needed(i)

    USE mo_species,      ONLY: speclist, naerospec, aero_idx

    IMPLICIT NONE

    !---function interface
    INTEGER, INTENT(IN) :: i

    !---function variables
    INTEGER  :: j

    !---executable procedure

    lsubname_needed = .FALSE.

    DO j=1,naerospec
       IF (speclist(aero_idx(j))%shortname == speclist(i)%shortname) lsubname_needed = .TRUE.
    END DO

  END FUNCTION lsubname_needed


  SUBROUTINE ham_define_modes(nclass, nsol, naerospec)

    USE mo_exception,               ONLY: finish
    USE mo_species,                 ONLY: speclist
    USE mo_ham_species,             ONLY: id_so4, &
                                          id_bc,  &
                                          id_oc,  &
                                          id_ss,  &
                                          id_du,  &
                                          id_wat
    USE mo_ham,                     ONLY: aerocomp, aerowater, naerocomp,        &
                                          new_aerocomp
    USE mo_ham_m7ctl,               ONLY: inucs, iaits, iaccs, icoas,          &
                                          iaiti, iacci, icoai,                 &
                                          iso4ns, iso4ks, iso4as, iso4cs,      &
                                          ibcks, ibcas, ibccs, ibcki,          &
                                          iocks, iocas, ioccs, iocki,          &
                                          issas, isscs,                        &
                                          iduas, iducs, iduai, iduci


    INTEGER, INTENT(in)          :: nclass    ! number of aerosol modes
    INTEGER, INTENT(in)          :: nsol      ! number of soluble modes
    INTEGER, INTENT(in)          :: naerospec ! number of aerosol species defined

    INTEGER           :: jn

    IF (ALLOCATED(aerocomp)) CALL finish('ham_define_modes', 'aerocomp already allocated!')    ! ###
    ALLOCATE(aerocomp(nclass * naerospec))
    ALLOCATE(aerowater(nsol))

    !---1) Aerosol compounds in the basic model (SO4, BC, OC, SS, DU) in applicable modes------
    naerocomp = 0
    
    !---sulphate in soluble modes    
    iso4ns = new_aerocomp(inucs, id_so4)
    iso4ks = new_aerocomp(iaits, id_so4)
    iso4as = new_aerocomp(iaccs, id_so4)
    iso4cs = new_aerocomp(icoas, id_so4)

    IF (.NOT. ALLOCATED(speclist(id_so4)%iaerocomp)) ALLOCATE(speclist(id_so4)%iaerocomp(nclass))
    speclist(id_so4)%iaerocomp(:)     = 0
    speclist(id_so4)%iaerocomp(inucs) = iso4ns     ! ### replaces former im7table construct
    speclist(id_so4)%iaerocomp(iaits) = iso4ks
    speclist(id_so4)%iaerocomp(iaccs) = iso4as
    speclist(id_so4)%iaerocomp(icoas) = iso4cs
    
    !---black carbon in aitken, accumulation and coarse soluble and aitken insoluble modes
    ibcks = new_aerocomp(iaits, id_bc)
    ibcas = new_aerocomp(iaccs, id_bc)
    ibccs = new_aerocomp(icoas, id_bc)
    ibcki = new_aerocomp(iaiti, id_bc)
    
    IF (.NOT. ALLOCATED(speclist(id_bc)%iaerocomp)) ALLOCATE(speclist(id_bc)%iaerocomp(nclass))
    speclist(id_bc)%iaerocomp(:)     = 0
    speclist(id_bc)%iaerocomp(iaits) = ibcks
    speclist(id_bc)%iaerocomp(iaccs) = ibcas
    speclist(id_bc)%iaerocomp(icoas) = ibccs
    speclist(id_bc)%iaerocomp(iaiti) = ibcki
    !---organic carbon in aitken, accumulation and coarse soluble and aitken insoluble modes
    iocks = new_aerocomp(iaits, id_oc)
    iocas = new_aerocomp(iaccs, id_oc)
    ioccs = new_aerocomp(icoas, id_oc)
    iocki = new_aerocomp(iaiti, id_oc)

    IF (.NOT. ALLOCATED(speclist(id_oc)%iaerocomp)) ALLOCATE(speclist(id_oc)%iaerocomp(nclass))
    speclist(id_oc)%iaerocomp(:)     = 0
    speclist(id_oc)%iaerocomp(iaits) = iocks
    speclist(id_oc)%iaerocomp(iaccs) = iocas
    speclist(id_oc)%iaerocomp(icoas) = ioccs
    speclist(id_oc)%iaerocomp(iaiti) = iocki

    !---sea salt in accumulation and coarse soluble modes
    issas = new_aerocomp(iaccs, id_ss)
    isscs = new_aerocomp(icoas, id_ss)

    IF (.NOT. ALLOCATED(speclist(id_ss)%iaerocomp)) ALLOCATE(speclist(id_ss)%iaerocomp(nclass))
    speclist(id_ss)%iaerocomp(:)     = 0
    speclist(id_ss)%iaerocomp(iaccs) = issas
    speclist(id_ss)%iaerocomp(icoas) = isscs

    !---dust in accumulation and coarse modes
    iduas = new_aerocomp(iaccs, id_du)
    iducs = new_aerocomp(icoas, id_du)
    iduai = new_aerocomp(iacci, id_du)
    iduci = new_aerocomp(icoai, id_du)

    IF (.NOT. ALLOCATED(speclist(id_du)%iaerocomp)) ALLOCATE(speclist(id_du)%iaerocomp(nclass))
    speclist(id_du)%iaerocomp(:)     = 0
    speclist(id_du)%iaerocomp(iaccs) = iduas
    speclist(id_du)%iaerocomp(icoas) = iducs
    speclist(id_du)%iaerocomp(iacci) = iduai
    speclist(id_du)%iaerocomp(icoai) = iduci

    !---2) Aerosol water ----------------------------------------------------------------------
!ham_ps: introduce aerosol water as species, but not component into aerocomp (no new_aerocomp call)
!        (this will need revisions)
    IF (.NOT. ALLOCATED(speclist(id_wat)%iaerocomp)) ALLOCATE(speclist(id_wat)%iaerocomp(nclass))
    speclist(id_wat)%iaerocomp(:)     = -1
!!$    speclist(id_wat)%iaerocomp(inucs) = 0
!!$    speclist(id_wat)%iaerocomp(iaits) = 0
!!$    speclist(id_wat)%iaerocomp(iaccs) = 0
!!$    speclist(id_wat)%iaerocomp(icoas) = 0
    DO jn = 1,nsol
       aerowater(jn)%iclass = jn
       aerowater(jn)%species => speclist(id_wat)
       aerowater(jn)%spid = id_wat
       ! note: %aero_idx and %idt are left undefined and should not be used.
    END DO
    
  END SUBROUTINE ham_define_modes 

  SUBROUTINE ham_define_bins(nclass, nsol, naerospec)
    ! --A. Laakso (FMI) 2013-05

    USE mo_exception,               ONLY: finish
    USE mo_species,                 ONLY: speclist
    USE mo_ham_species,             ONLY: id_so4, &
                                          id_bc,  &
                                          id_oc,  &
                                          id_ss,  &
                                          id_du,  &
                                          id_wat
    USE mo_ham,                     ONLY: aerocomp, aerowater, naerocomp,        &
                                          new_aerocomp
    USE mo_ham_salsactl,            ONLY: iso4b, ibcb,iocb,issb, idub,&
                                          in1a, in2a, in2b,  &
                                          fn1a, fn2a, fn2b

    ! >> thk: VBS
    USE mo_ham,                     ONLY:&
         sizeclass,                      &
         nsoa                            

    USE mo_ham_vbsctl,         ONLY:&
         vbs_ngroup,                &
         vbs_nvocs,                 &
         vbs_voc_prec,              &
         vbs_set,                   &
         nclass_vbs,                &                      
         laqsoa,                    &
         aqsoa_ngroup,              &
         aqsoa_set,                 &
         nclass_aqsoa
    ! << thk

    INTEGER, INTENT(in)          :: nclass    ! number of aerosol bins
    INTEGER, INTENT(in)          :: nsol      ! number of soluble bins
    INTEGER, INTENT(in)          :: naerospec ! number of aerosol species defined

    INTEGER           :: jn,i

    ! << thk: VBS
    INTEGER :: jv, jg, jg2, jc, spid
    ! >> thk

    iso4b(:) = 0
    ibcb(:) = 0
    iocb(:) = 0
    issb(:) = 0
    idub(:) = 0

    IF (ALLOCATED(aerocomp)) CALL finish('ham_define_bins', 'aerocomp already allocated!')    ! ###
    ALLOCATE(aerocomp(nclass * naerospec))
    ALLOCATE(aerowater(nsol))

    !---1) Aerosol compounds in the basic model (SO4, BC, OC, SS, DU) in applicable modes------
    naerocomp = 0
    
    !---sulphate  
!! write(0,*) '### new_ham_aerocomp: id_so4 = ',id_so4
    DO i = in1a,fn2b
        iso4b(i) = new_aerocomp(i, id_so4)
     END DO

    IF (.NOT. ALLOCATED(speclist(id_so4)%iaerocomp)) ALLOCATE(speclist(id_so4)%iaerocomp(nclass))
    speclist(id_so4)%iaerocomp(:)     = 0

    DO i = in1a,fn2b
        speclist(id_so4)%iaerocomp(i) = iso4b(i)    
    END DO

    !--oc
    DO i = in1a,fn2b
        iocb(i) = new_aerocomp(i, id_oc)
     END DO
    IF (.NOT. ALLOCATED(speclist(id_oc)%iaerocomp)) ALLOCATE(speclist(id_oc)%iaerocomp(nclass))

    speclist(id_oc)%iaerocomp(:)     = 0
    DO i = in1a,fn2b

        speclist(id_oc)%iaerocomp(i) = iocb(i)
    END DO

    !bc
    DO i = in2a,fn2b
        ibcb(i) = new_aerocomp(i, id_bc)
    END DO    
    IF (.NOT. ALLOCATED(speclist(id_bc)%iaerocomp)) ALLOCATE(speclist(id_bc)%iaerocomp(nclass))
    speclist(id_bc)%iaerocomp(:)     = 0
    DO i = in2a,fn2b
        speclist(id_bc)%iaerocomp(i) = ibcb(i)    
    END DO

    !dust
    DO i = in2a,fn2b
        idub(i) = new_aerocomp(i, id_du)
    END DO

    IF (.NOT. ALLOCATED(speclist(id_du)%iaerocomp)) ALLOCATE(speclist(id_du)%iaerocomp(nclass))
    speclist(id_du)%iaerocomp(:)     = 0

    DO i = in2a,fn2b  
        speclist(id_du)%iaerocomp(i) = idub(i)
    END DO

    !seasalt, only soluble
    DO i = in2a,fn2a
       issb(i) = new_aerocomp(i, id_ss)
    END DO
    IF (.NOT. ALLOCATED(speclist(id_ss)%iaerocomp)) ALLOCATE(speclist(id_ss)%iaerocomp(nclass))
    speclist(id_ss)%iaerocomp(:)     = 0
    DO i = in2a,fn2a
       speclist(id_ss)%iaerocomp(i) = issb(i)    
    END DO

    ! >> thk: VBS
    IF (nsoa == 2) THEN
       ! creating the per-class aerosol tracers
       DO jg = 1,vbs_ngroup
          spid = vbs_set(jg)%spid

          ! allocating memory for aerosol tracer indices in vbs_set
          IF (.NOT. ALLOCATED(vbs_set(jg)%idx)) THEN
             ALLOCATE(vbs_set(jg)%idx(nclass))
          END IF
          vbs_set(jg)%idx(:) = 0

          IF (vbs_set(jg)%lcreateaero) THEN

             ! allocating memory for aerosol tracer indices in speclist
             IF (.NOT. ALLOCATED(speclist(spid)%iaerocomp)) THEN
                ALLOCATE(speclist(spid)%iaerocomp(nclass))
             END IF
             speclist(spid)%iaerocomp(:) = 0

             ! creating tracers and saving id's
             DO i = 1, nclass
                IF (sizeclass(i)%lsoainclass) THEN
                   vbs_set(jg)%idx(i) = new_aerocomp(i, spid)
                   speclist(spid)%iaerocomp(i) = vbs_set(jg)%idx(i)
                END IF
             END DO
          ELSE

             ! getting the tracer ids from the given spid:
             DO i = 1, nclass
                IF (sizeclass(i)%lsoainclass) THEN
                   vbs_set(jg)%idx(i) = speclist(vbs_set(jg)%spid_aero)%iaerocomp(i)
                END IF
             END DO
          END IF
       END DO
       ! counting all modes/bins that include VBS soa:
       nclass_vbs = 0
       DO jc = 1,nclass
          IF (sizeclass(jc)%lsoainclass) nclass_vbs = nclass_vbs+1
       END DO

       IF (laqsoa) THEN
          ! creating the per-class aerosol tracers
          DO jg2 = 1,aqsoa_ngroup
             spid = aqsoa_set(jg2)%spid


             ! allocating memory for aerosol tracer indices in aqsoa_set
             IF (.NOT. ALLOCATED(aqsoa_set(jg2)%idx)) THEN
                ALLOCATE(aqsoa_set(jg2)%idx(nclass))
             END IF
             aqsoa_set(jg2)%idx(:) = 0

             IF (aqsoa_set(jg2)%lcreateaero) THEN

                ! allocating memory for aerosol tracer indices in speclist
                IF (.NOT. ALLOCATED(speclist(spid)%iaerocomp)) THEN
                   ALLOCATE(speclist(spid)%iaerocomp(nclass))
                END IF
                speclist(spid)%iaerocomp(:) = 0

                ! creating tracers and saving id's
                DO i = 1, nclass
                   IF (sizeclass(i)%lsoainclass) THEN
                      aqsoa_set(jg2)%idx(i) = new_aerocomp(i, spid)
                      speclist(spid)%iaerocomp(i) = aqsoa_set(jg2)%idx(i)
                   END IF
                END DO
             ELSE

                ! getting the tracer ids from the given spid:
                DO i = 1, nclass
                   IF (sizeclass(i)%lsoainclass) THEN
                      aqsoa_set(jg2)%idx(i) = speclist(aqsoa_set(jg2)%spid_aero)%iaerocomp(i)
                   END IF
                END DO
             END IF
          END DO
          ! counting all modes/bins that include wet soa:
          nclass_aqsoa = 0
          DO jc = 1,nclass
             IF (sizeclass(jc)%lsoainclass) nclass_aqsoa = nclass_aqsoa+1
          END DO

       END IF ! laqsoa
    END IF ! nsoa == 2
    ! << thk

    !---2) Aerosol water ----------------------------------------------------------------------
!ham_ps: introduce aerosol water as species, but not component into aerocomp (no new_aerocomp call)
!        (this will need revisions) no sepa kiva
    IF (.NOT. ALLOCATED(speclist(id_wat)%iaerocomp)) ALLOCATE(speclist(id_wat)%iaerocomp(nclass))
    speclist(id_wat)%iaerocomp(:)     = -1
    DO jn = 1,nsol
       aerowater(jn)%iclass = jn
       aerowater(jn)%species => speclist(id_wat)
       aerowater(jn)%spid = id_wat
       ! note: %aero_idx and %idt are left undefined and should not be used.
    END DO
    
  END SUBROUTINE ham_define_bins 


  SUBROUTINE ham_initialize

    ! Purpose:
    ! ---------
    ! Initializes constants and parameters used in the HAM aerosol model.
    ! Performs consistency checks.
    !
    ! Author:
    ! ---------
    ! Philip Stier, MPI                           03/2003
    ! Martin Schultz, FZJ                         09/2009 - renamed to ham_
    !                                                     - cleanup
    !
    ! Interface:
    ! ---------
    ! *ham_initialize*  is called from *init_subm* in mo_submodel_interface
    !                    needs to be called after initialization of the 
    !                    submodel as it may make use of parameters in mo_ham_m7ctl,
    !

    USE mo_exception,       ONLY: message, em_error
    USE mo_ham,             ONLY: naerorad,nclass, nham_subm,       &
                                  HAM_BULK,                         &
                                  HAM_M7,                           &
                                  HAM_SALSA
    USE mo_ham_m7ctl,       ONLY: nwater
    USE mo_ham_rad,         ONLY: ham_rad_initialize
    USE mo_ham_kappa,       ONLY: start_kappa
    USE mo_ham_m7_trac,     ONLY: ham_M7_set_idt
    USE mo_ham_salsa_trac,  ONLY: ham_salsa_set_idt
    USE mo_activ,           ONLY: nfrzmod
    !>>dod (redmine #44) import of seasalt emission schemes from HAM2
    USE mo_ham_m7_emi_seasalt, ONLY: start_emi_seasalt
    !<<dod

!>>SF #390 (for security)
    USE mo_param_switches, ONLY: lcdnc_progn
    USE mo_ham, ONLY: nwetdep
!<<SF #390

    IMPLICIT NONE

    !--- set tracer ids needed for example in wetdep scheme
    !    this must be done *after* definition of MOZ tracers in order to account for HAMMOZ coupling
    SELECT CASE(nham_subm)
        CASE(HAM_BULK)
        CASE(HAM_M7) 
            CALL ham_m7_set_idt
        CASE(HAM_SALSA)
            CALL ham_salsa_set_idt
    END SELECT
    
    
    !--- initialize ham_rad, kappa and seasalt emission bin schemes
    IF (naerorad>0) CALL ham_rad_initialize(nclass)

    IF (nwater == 1 .AND. nham_subm == HAM_M7) CALL start_kappa
    !>>dod (redmine #44) import of seasalt emission schemes from HAM2
    CALL start_emi_seasalt
    !<<dod

    !--- Set the number of freezing modes:
    nfrzmod = 1  !SF WARNING: no other value is possible for now (see mo_ham_freezing.f90)

!>>SF security
    IF (nfrzmod /= 1) THEN
       call message('ham_initialize','nfrzmod must be equal to 1!',level=em_error)
    ENDIF
!<<SF security

!>>SF #390
   IF (.NOT. lcdnc_progn .AND. nwetdep == 3) THEN
      call message('ham_initialize','nwetdep = 3 is not possible when lcdnc_progn is false!',level=em_error)
   ENDIF
!<<SF #390

  END SUBROUTINE ham_initialize

  ! --- initialisation of memory (output and diagnostic streams etc.)

  SUBROUTINE ham_init_memory

  USE mo_ham,              ONLY: nclass,naerorad, nsoa, lgcr, nham_subm, HAM_M7
  USE mo_ham_streams,      ONLY: new_stream_ham, new_stream_ham_rad
!!++mgs:removed - belongs in mo_submodel_interface      USE mo_aero_activ,     ONLY: construct_stream_activ
!!++mgs:removed - belongs in mo_submodel_interface      USE mo_conv,           ONLY: construct_stream_conv
  USE mo_ham_soa,          ONLY: construct_soa_streams
  USE mo_ham_m7ctl,        ONLY: nsnucl,     &
                                 nonucl,     &
                                 lnucl_stat
  USE mo_ham_gcrion,       ONLY: gcr_initialize
  USE mo_ham_m7_nucl,      ONLY: ham_nucl_initialize
  USE mo_ham_m7_nucl_diag, ONLY: ham_nucl_diag_initialize

  CALL new_stream_ham (nclass)     !>>dod<<
!!mgs!!  CALL construct_stream_input
  IF (naerorad>0)         CALL new_stream_ham_rad            !Philip
  IF (lgcr)               CALL gcr_initialize                !Jan
  IF (nham_subm == HAM_M7) THEN
     IF (nsnucl+nonucl.gt.0) CALL ham_nucl_initialize           !Jan
     IF (lnucl_stat)         CALL ham_nucl_diag_initialize      !Jan
  ENDIF
  IF (nsoa == 1)               CALL construct_soa_streams         !Declan

  END SUBROUTINE ham_init_memory


  ! --- release of memory

  SUBROUTINE ham_free_memory

    USE mo_ham,              ONLY: naerorad, nham_subm, HAM_M7, lgcr
    USE mo_ham_rad,          ONLY: ham_rad_mem_cleanup
    USE mo_ham_dust,         ONLY: bgc_dust_cleanup
    USE mo_ham_m7ctl,        ONLY: nsnucl, nonucl, lnucl_stat, nwater
    USE mo_ham_gcrion,       ONLY: gcr_cleanup
    USE mo_ham_m7_nucl,      ONLY: ham_nucl_cleanup
    USE mo_ham_m7_nucl_diag, ONLY: ham_nucl_diag_cleanup
    USE mo_ham_kappa,        ONLY: term_kappa

    IMPLICIT NONE

    IF (nham_subm == HAM_M7) THEN
       IF (nwater == 1) CALL term_kappa
       IF (nsnucl+nonucl.gt.0) CALL ham_nucl_cleanup
       IF (lnucl_stat)         CALL ham_nucl_diag_cleanup
    ENDIF
    IF (naerorad>0)         CALL ham_rad_mem_cleanup
    CALL bgc_dust_cleanup
    IF (lgcr)               CALL gcr_cleanup

  END SUBROUTINE ham_free_memory

END MODULE mo_ham_init






