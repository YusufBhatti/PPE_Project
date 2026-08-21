!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_hammoz_aerocom_MMPPEmon.f90
!!
!! \brief
!! Module for AeroCom Multi Model Perturbed Physics Ensemble (MMPPE) monthly diagnostics
!!
!! \author D. Neubauer (ETH Zurich)
!!  adapted from AerChemMIP diagnostics from S. Ferrachat
!!
!! D. Neubauer, david.neubauer@env.ethz.ch
!!
!! \revision_history
!!   -# D. Neubauer (ETH Zurich) - original code (2019-03-12)
!!   -# A. Arifi - monthly version of the MMPPE diagnostics (2026-08-21)
!!
!! \limitations
!! None
!!
!! \details
!! Monthly counterpart of the MMPPE stream (mo_hammoz_aerocom_MMPPE): the
!! diagnostics are meant to be accumulated over the output interval
!! (laccu = .TRUE.), so that the post-processed fields are monthly means rather
!! than instantaneous snapshots.
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
MODULE mo_hammoz_aerocom_MMPPEmon

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_MMPPEmon_stream
  PUBLIC :: update_MMPPEmon_diags

  TYPE (t_stream), PUBLIC, POINTER :: acmmppemon

  REAL(dp), PUBLIC, POINTER :: nmrcdnc(:,:,:)
  REAL(dp), PUBLIC, POINTER :: nmricnc(:,:,:)
  REAL(dp), PUBLIC, POINTER :: loadbc(:,:)
  REAL(dp), PUBLIC, POINTER :: loaddu(:,:)
  REAL(dp), PUBLIC, POINTER :: loadoc(:,:)
  REAL(dp), PUBLIC, POINTER :: loadso4(:,:)
  REAL(dp), PUBLIC, POINTER :: loadss(:,:)
  REAL(dp), PUBLIC, POINTER :: loadwat(:,:)

  !-- indices of the burden species in the local zload array (see update_MMPPEmon_diags)
  INTEGER, PARAMETER :: nload = 6
  INTEGER, PARAMETER :: iload_bc  = 1, &
                        iload_du  = 2, &
                        iload_oc  = 3, &
                        iload_so4 = 4, &
                        iload_ss  = 5, &
                        iload_wat = 6

  CONTAINS

  !------------------------------------------------

  SUBROUTINE construct_MMPPEmon_stream

    USE mo_hammoz_aerocom_data, ONLY: AEROCOM_FILETYPE
    USE mo_memory_base,       ONLY: new_stream, add_stream_element, AUTO,  &
                                    default_stream_setting, add_stream_reference
    USE mo_time_event,        ONLY: io_time_event

    TYPE(io_time_event) :: put_interval

    !-- set output interval
    put_interval%counter      = 1
    put_interval%unit         = 'months'
    put_interval%adjustment   = 'last'
    put_interval%offset       = 0

    !-- Create new stream:
    CALL new_stream (acmmppemon ,'acmmppemon', &
                     filetype = AEROCOM_FILETYPE, &
                     lrerun = .TRUE., &
                     interval = put_interval, &
                     lpost = .TRUE.)

    !-- Add standard fields for post-processing:
    CALL default_stream_setting (acmmppemon, &
                                 table = 199, &
                                 code = AUTO )

    !-- Basic vars (may sometimes be necessary)
    CALL add_stream_reference (acmmppemon, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (acmmppemon, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (acmmppemon, 'aps'     ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (acmmppemon, 'gboxarea','geoloc',lpost=.TRUE.)

    !-- Diagnostics table'
    CALL add_stream_element (acmmppemon, 'nmrcdnc', nmrcdnc, &
        longname = 'cloud droplet number mixing ratio', &
        units = '1 kg-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acmmppemon, 'nmricnc', nmricnc, &
        longname = 'ice crystal number mixing ratio', &
        units = '1 kg-1', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acmmppemon, 'loadbc', loadbc, &
        longname = 'Column black carbon mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acmmppemon, 'loaddu', loaddu, &
        longname = 'Column dust mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acmmppemon, 'loadoc', loadoc, &
        longname = 'Column organic carbon mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acmmppemon, 'loadso4', loadso4, &
        longname = 'Column sulphate mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acmmppemon, 'loadss', loadss, &
        longname = 'Column seasalt mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

    CALL add_stream_element (acmmppemon, 'loadwat', loadwat, &
        longname = 'Column aerosol water mass load', &
        units = 'kg m-2', &
        laccu = .TRUE., &
        lpost = .TRUE., &
        lrerun = .TRUE. )

  END SUBROUTINE construct_MMPPEmon_stream

  SUBROUTINE update_MMPPEmon_diags(kproma, kbdim, klev, krow)

    USE mo_time_control,   ONLY: delta_time
    USE mo_memory_g1a,     ONLY: xtm1
    USE mo_ham_m7_trac,    ONLY: idt_cdnc_ham, idt_icnc_ham
    USE mo_ham,            ONLY: naerocomp, aerocomp, aerowater, nsol
    USE mo_ham_species,    ONLY: id_bc, id_du, id_oc, id_so4, id_ss
    USE mo_vphysc,         ONLY: vphysc
    USE mo_physical_constants, ONLY: grav

    INTEGER, INTENT(in) :: kproma, kbdim, klev, krow

    INTEGER  :: jn, jk, jt, ispec, iload
    REAL(dp) :: zdpg(kbdim,klev)    !< delta p over g [kg m-2]
    REAL(dp) :: zload(kbdim,nload)  !< instantaneous burdens [kg m-2]

    ! NUMBER CONCENTATION (km-3): CDNC and IDNC

    !-- nmrcdnc, nmricnc (CDNC and ICNC tracers, see ham_define_tracer)
    nmrcdnc(1:kproma,:,krow) = nmrcdnc(1:kproma,:,krow) &
         + xtm1(1:kproma,:,idt_cdnc_ham,krow) * delta_time

    nmricnc(1:kproma,:,krow) = nmricnc(1:kproma,:,krow) &
         + xtm1(1:kproma,:,idt_icnc_ham,krow) * delta_time

    ! BURDEN (kg m-2): AEROSOL SPECIES

    !-- loadbc, loaddu, loadoc, loadso4, loadss, loadwat
    !   Column burdens as in xt_burden (mo_tracer_processes): the mass mixing ratio
    !   [kg kg-1] is weighted with dp/g and summed over all levels to give [kg m-2].
    !   The burdens are multiplied by dt here and divided by the output interval by
    !   the stream management (laccu), which gives a mean over the output interval.
    !   Difference to xt_burden: everything is taken at t-dt, i.e. the tracer is
    !   xtm1 and dp/g comes from the t-dt pressures held in vphysc, whereas
    !   xt_burden evaluates both at t+dt (xtm1+xtte*time_step_len and papp1/paphp1,
    !   see luse_p1_vars in mo_submodel_interface). This is consistent with the
    !   other AeroCom diagnostics, which all accumulate xtm1.

    !--- 1) Calculate auxiliary variable dp/g :

    !--- Uppermost level:
    zdpg(1:kproma,1) = 2._dp*(vphysc%aphm1(1:kproma,2,krow)-vphysc%apm1(1:kproma,1,krow))/grav
    !--- Other levels:
    zdpg(1:kproma,2:klev) = (vphysc%aphm1(1:kproma,3:klev+1,krow)   &
                            -vphysc%aphm1(1:kproma,2:klev,krow))/grav

    !--- 2) Sum the column integral over all tracers of a given species:

    zload(1:kproma,:) = 0._dp

    DO jn = 1, naerocomp
       ispec = aerocomp(jn)%spid
       jt    = aerocomp(jn)%idt

       IF (ispec == id_bc) THEN
          iload = iload_bc
       ELSE IF (ispec == id_du) THEN
          iload = iload_du
       ELSE IF (ispec == id_oc) THEN
          iload = iload_oc
       ELSE IF (ispec == id_so4) THEN
          iload = iload_so4
       ELSE IF (ispec == id_ss) THEN
          iload = iload_ss
       ELSE
          CYCLE    ! nothing to be done for this species
       END IF

       DO jk = 1, klev
          zload(1:kproma,iload) = zload(1:kproma,iload) &
                                + xtm1(1:kproma,jk,jt,krow)*zdpg(1:kproma,jk)
       END DO
    END DO

    !--- aerosol water is not part of aerocomp, but held in aerowater:
    DO jn = 1, nsol
       jt = aerowater(jn)%idt

       DO jk = 1, klev
          zload(1:kproma,iload_wat) = zload(1:kproma,iload_wat) &
                                    + xtm1(1:kproma,jk,jt,krow)*zdpg(1:kproma,jk)
       END DO
    END DO

    !--- 3) Accumulate:

    loadbc(1:kproma,krow)  = loadbc(1:kproma,krow)  + zload(1:kproma,iload_bc) *delta_time
    loaddu(1:kproma,krow)  = loaddu(1:kproma,krow)  + zload(1:kproma,iload_du) *delta_time
    loadoc(1:kproma,krow)  = loadoc(1:kproma,krow)  + zload(1:kproma,iload_oc) *delta_time
    loadso4(1:kproma,krow) = loadso4(1:kproma,krow) + zload(1:kproma,iload_so4)*delta_time
    loadss(1:kproma,krow)  = loadss(1:kproma,krow)  + zload(1:kproma,iload_ss) *delta_time
    loadwat(1:kproma,krow) = loadwat(1:kproma,krow) + zload(1:kproma,iload_wat)*delta_time

  END SUBROUTINE update_MMPPEmon_diags

END MODULE mo_hammoz_aerocom_MMPPEmon
