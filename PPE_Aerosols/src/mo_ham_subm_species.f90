!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!! \filename 
!! mo_ham_subm_species.f90
!!
!! \brief
!! Species mapping from ECHAM species list to condensed gas-phase and aerosol lists in M7
!!
!! \author Martin G. Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# Martin G. Schultz (FZ Juelich) - original code (2009-10)
!!   -# Harri Kokkola (FMI) - Implementation of SALSA aerosol microphysics model (2014)
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

MODULE mo_ham_subm_species

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: map_ham_subm_species

  PUBLIC :: isubm_so2, isubm_so4g, isubm_oc, isubm_ocnv, isubm_so4, isubm_wat


  !-- species indices for microphysics routines
  INTEGER :: isubm_so2, isubm_so4g, isubm_oc, isubm_ocnv, isubm_so4, isubm_wat
  
  ! >> thk: volatility basis set (VBS)
  INTEGER :: isubm_bc, isubm_ss, isubm_du
  ! << thk

  CONTAINS

! ---------------------------------------------------------------------------
!  map_ham_subm_species: construct condensed species lists for M7 gas phase and 
!  aerosol species from the general speclist
! ---------------------------------------------------------------------------

SUBROUTINE map_ham_subm_species

  USE mo_species,         ONLY: nspec, speclist
  USE mo_ham,             ONLY: nsoa, nsoaspec, nham_subm, HAM_SALSA
  USE mo_ham_species,     ONLY: id_so2, id_so4g, id_ocnv,                   &
                                id_so4, id_bc, id_oc, id_ss, id_du, id_wat
  USE mo_ham_soa,         ONLY: soaprop
  USE mo_ham,             ONLY: immr2molec, ivmr2molec, immr2ug, &
                                subm_ngasspec, subm_gasspec, subm_naerospec, subm_aerospec, &
                                subm_naerospec_nowat, subm_aerospec_nowat
  USE mo_exception,       ONLY: message, message_text, em_info, em_param
  USE mo_util_string,     ONLY: separator
                                
  USE mo_ham_vbsctl,      ONLY: vbs_nvocs, vbs_voc_prec,  vbs_ngroup, &
                                vbs_set, laqsoa,aqsoa_ngroup, aqsoa_set

  INTEGER         :: jt, jm, jb, jv

  
  DO jt = 1, nspec
    IF (jt == id_so2)   CALL new_subm_gasspec(id_so2,  immr2molec, isubm_so2)
    IF (jt == id_so4g)  CALL new_subm_gasspec(id_so4g, immr2molec, isubm_so4g)
    IF(nham_subm == HAM_SALSA) THEN
       IF (jt == id_ocnv)  CALL new_subm_gasspec(id_ocnv, immr2molec, isubm_ocnv) 
    END IF

    IF (jt == id_so4)   CALL new_subm_aerospec(id_so4, immr2molec, isubm_so4)
    IF (jt == id_bc)    CALL new_subm_aerospec(id_bc,  immr2ug,    isubm_bc)
    IF (jt == id_oc)    CALL new_subm_aerospec(id_oc,  immr2ug,    isubm_oc)
    IF (jt == id_ss)    CALL new_subm_aerospec(id_ss,  immr2ug,    isubm_ss)
    IF (jt == id_du)    CALL new_subm_aerospec(id_du,  immr2ug,    isubm_du)
    IF (jt == id_wat)   CALL new_subm_aerospec(id_wat, immr2ug,    isubm_wat)

    SELECT CASE (nsoa)
       CASE(1) 

          DO jm = 1,nsoaspec
             IF (jt == soaprop(jm)%spid_soa) CALL new_subm_aerospec(soaprop(jm)%spid_soa, immr2ug)
          END DO

       CASE(2) !thk: VBS

          DO jb = 1, vbs_ngroup
             IF (jt == vbs_set(jb)%spid) THEN
                IF (vbs_set(jb)%lcreateaero) THEN
                   CALL new_subm_gasspec(jt, immr2molec, vbs_set(jb)%id_gasspec)
                   CALL new_subm_aerospec(jt, immr2ug,   vbs_set(jb)%id_aerospec)
                   vbs_set(jb)%id_vols = subm_naerospec_nowat
                ELSE
                   DO jm = 1,subm_ngasspec
                      IF (subm_gasspec(jm) == jt) vbs_set(jb)%id_gasspec = jm
                   END DO

!>>SF special case
                   IF (vbs_set(jb)%spid_aero /= 0) THEN
                      vbs_set(jb)%id_vols = subm_naerospec_nowat
                   ENDIF
!<<SF
                END IF
             END IF
          END DO
          DO jv = 1, vbs_nvocs
             if (jt == vbs_voc_prec(jv)%spid) THEN
                CALL new_subm_gasspec(jt, immr2molec, vbs_voc_prec(jv)%id_gasspec)
             END IF
          END DO

          IF (laqsoa) THEN
             DO jb = 1, aqsoa_ngroup
                IF (jt == aqsoa_set(jb)%spid) THEN
                   IF (aqsoa_set(jb)%lcreateaero) THEN
                      CALL new_subm_gasspec(jt, immr2molec, aqsoa_set(jb)%id_gasspec)
                      CALL new_subm_aerospec(jt, immr2ug,   aqsoa_set(jb)%id_aerospec)
                      aqsoa_set(jb)%id_aqsoa = subm_naerospec_nowat
                   ELSE
                      DO jm = 1,subm_ngasspec
                         IF (subm_gasspec(jm) == jt) aqsoa_set(jb)%id_gasspec = jm
                      END DO
                   END IF
                END IF
             END DO
          END IF !laqsoa

    END SELECT

  END DO

  !-- print status
  CALL message('',separator)
  CALL message('', 'Aerosol microphysics species lists', level=em_info)
  CALL message('','',level=em_param)
  CALL message('','Aerosol microphysics gas species id   Species id   Name',level=em_param)

  DO jt=1,subm_ngasspec
     WRITE(message_text,'(i0,t21,i0,t34,a)') jt, subm_gasspec(jt), speclist(subm_gasspec(jt))%longname
     CALL message('',message_text,level=em_param)
  END DO

  CALL message('','',level=em_param)
  CALL message('','Aerosol microphysics species id   Species id   Name',level=em_param)

  DO jt=1,subm_naerospec
     WRITE(message_text,'(i0,t21,i0,t34,a)') jt, subm_aerospec(jt), speclist(subm_aerospec(jt))%longname
     CALL message('',message_text,level=em_param)  
  END DO

  CALL message('','',level=em_param)
  CALL message('',separator)

END SUBROUTINE map_ham_subm_species
 
! ---------------------------------------------------------------------------
!  new_subm_gasspec: add a species id of a gas species to the list of species to
!  be considered in M7 processes.
! ---------------------------------------------------------------------------

SUBROUTINE new_subm_gasspec(nspid, nunitconv, idlocal)

    USE mo_exception,           ONLY : finish
    USE mo_species,             ONLY : nmaxspec
    USE mo_ham,                 ONLY : subm_ngasspec, subm_gasspec, subm_gasunitconv

    INTEGER, INTENT(in)             :: nspid, nunitconv
    INTEGER, INTENT(out), OPTIONAL  :: idlocal              ! local species id

    !### note: minimial error checking, because use of this routine is practically hardwired in HAM

    subm_ngasspec = subm_ngasspec + 1
    IF (subm_ngasspec > nmaxspec) CALL finish('new_subm_gasspec',   &
                                         'Number of gas species for HAM (subm_ngasspec) exceeds nmaxspec!')

    subm_gasspec(subm_ngasspec)     = nspid
    subm_gasunitconv(subm_ngasspec) = nunitconv

    IF (PRESENT(idlocal)) idlocal = subm_ngasspec

END SUBROUTINE new_subm_gasspec


! ---------------------------------------------------------------------------
!  new_subm_aerospec: add a species id of a aero species to the list of species to
!  be considered in aerosol processes.
! ---------------------------------------------------------------------------

  SUBROUTINE new_subm_aerospec(nspid, nunitconv, idlocal)
  
    USE mo_exception,   ONLY: finish
    USE mo_species,     ONLY: nmaxspec
    USE mo_ham,         ONLY: subm_naerospec, subm_aerospec, subm_aerounitconv, subm_aero_idx, &
                              subm_naerospec_nowat, subm_aerospec_nowat !SF for convenience
    USE mo_ham_species, ONLY: id_wat        
  
    INTEGER, INTENT(in)    :: nspid, nunitconv
    INTEGER, INTENT(out), OPTIONAL  :: idlocal              ! local species id for M7 routines
    !### note: minimial error checking, because use of this routine is practically hardwired in HAM
  
    subm_naerospec = subm_naerospec + 1 
    IF (subm_naerospec > nmaxspec) CALL finish('new_subm_aerospec',        &
                   'Number of aerosol species exceeds nmaxspec!')
  
    subm_aerospec(subm_naerospec)     = nspid 
    subm_aerounitconv(subm_naerospec) = nunitconv
    ! reverse mapping
    subm_aero_idx(nspid) = subm_naerospec
  
    IF (PRESENT(idlocal)) idlocal = subm_naerospec

    !>>SF for convenience only: introduce a special mapping for all but water aero species
    !     (mostly used in SALSA for now, but could be more widely used in future)
    IF (nspid /= id_wat) THEN 
       subm_naerospec_nowat = subm_naerospec_nowat + 1
       subm_aerospec_nowat(subm_naerospec_nowat) = nspid 
    ENDIF
    !<<SF

    END SUBROUTINE new_subm_aerospec

 

END MODULE mo_ham_subm_species
