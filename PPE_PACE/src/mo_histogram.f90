MODULE mo_histogram

  USE mo_kind,                    ONLY: dp
  USE mo_memory_base,             ONLY: new_stream, add_stream_element, &
                                        default_stream_setting, AUTO, t_stream, &
                                        add_stream_reference
  USE mo_time_event,              ONLY: io_time_event
  USE mo_time_control,            ONLY: p_bcast_event
  USE mo_linked_list,             ONLY: HYBRID
  USE mo_netcdf,                  ONLY: HYBRID_H
  USE mo_p3_fields,               ONLY: lfull_diag
  USE mo_cloud_utils,             ONLY: effective_ice_radius, effective_liquid_radius

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: simple_histogram, double_histogram, init_histos, histogram_diagnostic
  PUBLIC :: add_multiple_stream_elements, t_globptr, t_surfptr, t_hist, LINSCALE, LOGSCALE
  PUBLIC :: frac_spec, mxt_spec, mmr_spec, prcp_spec, ten_spec, dp_spec
  PUBLIC :: cosp_ct_lf, cosp_ct_lf_time, cosp_ct_cc, cosp_ct_cc_time

  INTEGER, PARAMETER :: LINSCALE = 1 ! bin data on linearly scaled bins
  INTEGER, PARAMETER :: LOGSCALE = 2 ! bin data on logarithmically scaled bins

  TYPE t_globptr
    REAL(dp), POINTER :: v(:,:,:)
  END TYPE t_globptr

  TYPE t_surfptr
    REAL(dp), POINTER :: v(:,:)
  END TYPE t_surfptr

  TYPE :: t_hist
     INTEGER  :: nbin
     REAL(dp) :: vmin
     REAL(dp) :: vmax
     INTEGER  :: SCALE
  END TYPE t_hist

  ! HISTOGRAMS (x, y, z) being: x: number of bins, y minimum value, z maximum value
  TYPE(t_hist) :: frac_spec = t_hist(10, 0, 1, LINSCALE) !< supercooled liquid fraction
  TYPE(t_hist) :: ten_spec = t_hist(10, 0, 10, LINSCALE) !< histogram [1, ..., 10]
  TYPE(t_hist) :: dp_spec = t_hist(10, 0, 100000, LINSCALE) !< pressure difference (Pa)
  TYPE(t_globptr), POINTER :: clf_cc(:), ccc(:)
  TYPE(t_hist) :: mxt_spec = t_hist(16, 228.15_dp, 276.15_dp, LINSCALE) !< mixed-phase temperature (K)
  TYPE(t_hist) :: mmr_spec = t_hist(10, 1e-6_dp, 1e-2_dp, LOGSCALE) !< ice cloud content (kg/kg)
  TYPE(t_hist) :: prcp_spec = t_hist(10, 0.001_dp, 100._dp, LOGSCALE) !< precipitation (mm/h)
  TYPE(t_hist) :: wp_spec = t_hist(20, 0, 200, LINSCALE)              !< water path (kg/m-2)
  TYPE(t_hist) :: radl_spec = t_hist(20, 0, 20, LINSCALE)             !< liquid radius (um)
  TYPE(t_hist) :: radi_spec = t_hist(40, 0, 100, LINSCALE)            !< ice radius (um)
  TYPE(t_globptr), POINTER :: cxi(:), cxl(:)
  TYPE(t_surfptr), POINTER :: ctlf_cc(:), cttoplf_cc(:), cprcp_tot(:), cprcp_rain(:), cprcp_snow(:), &
                              cprcp_phase(:), cprcp_topphase(:), ct_lf(:), ct_cc(:), &
                              ctlf_hetf(:), cthetf_cci(:), ct_hetf(:), ct_cci(:), ctlf_cci(:), &
                              cthetf_lf(:), cttophetf_cci(:), cttop_hetf(:), cttop_cci(:), &
                              ctlf_cc_night(:), ctlf_cc_day(:), ct_lf_day(:), ct_lf_night(:), &
                              ct_cc_day(:), ct_cc_night(:), clwp_ccl(:), ciwp_cci(:), &
                              crleff_ccl_top(:), crieff_cci_top(:), csosif_cc(:), csocmf_cc(:), &
                              ctoliqf_cci(:), coliqf_cci(:), chetf_cci(:), ct_oliqf(:)
  
  ! quick fix to output cosp temperature scale
  TYPE(t_surfptr), POINTER :: cosp_ct_lf(:), cosp_ct_cc(:), cosp_ct_lf_time(:), cosp_ct_cc_time(:)

  REAL(dp), PUBLIC, POINTER :: cld_time(:,:,:)

  INTERFACE add_multiple_stream_elements
     MODULE PROCEDURE add_multiple_stream_elements_2d
     MODULE PROCEDURE add_multiple_stream_elements_3d
  END INTERFACE add_multiple_stream_elements

  INTERFACE linear_index
     MODULE PROCEDURE linear_index_3d
     MODULE PROCEDURE linear_index_2d
  END INTERFACE linear_index

  INTERFACE log_index
     MODULE PROCEDURE log_index_3d
     MODULE PROCEDURE log_index_2d
  END INTERFACE log_index

  INTERFACE simple_histogram
     MODULE PROCEDURE simple_histogram_3d
     MODULE PROCEDURE simple_histogram_2d
  END INTERFACE simple_histogram

  INTERFACE double_histogram
     MODULE PROCEDURE double_histogram_3d
     MODULE PROCEDURE double_histogram_2d
  END INTERFACE double_histogram

  CONTAINS

  SUBROUTINE add_multiple_stream_elements_3d(&
                !-- IN
                stream, basename, nelements, &
                !-- INOUT
                ptr_array, &
                !-- OPTIONAL
                leveltype, lpost, laccu, &
                reset, contnorest)

    USE mo_exception,             ONLY: finish

    TYPE(t_stream), INTENT(inout)           :: stream
    CHARACTER(*), INTENT(in)                :: basename
    INTEGER, INTENT(in)                     :: nelements
    INTEGER, INTENT(in), OPTIONAL           :: leveltype
    LOGICAL, INTENT(in), OPTIONAL           :: lpost
    LOGICAL, INTENT(in), OPTIONAL           :: laccu
    REAL(dp), INTENT(in), OPTIONAL          :: reset
    LOGICAL, INTENT(in), OPTIONAL           :: contnorest

    TYPE(t_globptr), POINTER, INTENT(inout) :: ptr_array(:)

    CHARACTER(LEN=32) :: ichar
    INTEGER           :: i, j, ilength, nlength

    ! allocate nelements 3-d variables on full levels (layer centres) and store pointers
    ALLOCATE(ptr_array(nelements))
    WRITE(ichar,*) nelements
    nlength=LEN_TRIM(ADJUSTL(ichar))
    IF ( nlength > 27 ) THEN
       CALL finish('mo_histogram, add_multiple_stream_elements_3d', 'too many 3d variables')
    END IF
    DO i=1,nelements
      WRITE(ichar,*) i
      ilength=LEN_TRIM(ADJUSTL(ichar))
      DO j=ilength+1,nlength
         ichar='0'//TRIM(ADJUSTL(ichar))
      END DO
      CALL add_stream_element (stream,basename//TRIM(ADJUSTL(ichar)),ptr_array(i)%v, &
                               leveltype=leveltype, lpost=lpost, laccu=laccu, reset=reset, &
                               contnorest=contnorest)
    END DO


  END SUBROUTINE add_multiple_stream_elements_3d

  SUBROUTINE add_multiple_stream_elements_2d(&
                !-- IN
                stream, basename, nelements, &
                !-- INOUT
                ptr_array, &
                !-- OPTIONAL
                leveltype, lpost, laccu, &
                reset, contnorest)

    USE mo_exception,             ONLY: finish

    TYPE(t_stream), INTENT(inout)           :: stream
    CHARACTER(*), INTENT(in)                :: basename
    INTEGER, INTENT(in)                     :: nelements
    INTEGER, INTENT(in), OPTIONAL           :: leveltype
    LOGICAL, INTENT(in), OPTIONAL           :: lpost
    LOGICAL, INTENT(in), OPTIONAL           :: laccu
    REAL(dp), INTENT(in), OPTIONAL          :: reset
    LOGICAL, INTENT(in), OPTIONAL           :: contnorest

    TYPE(t_surfptr), POINTER, INTENT(inout) :: ptr_array(:)

    CHARACTER(LEN=32) :: ichar
    INTEGER           :: i, j, ilength, nlength

    ! allocate nelements 3-d variables on full levels (layer centres) and store pointers
    ALLOCATE(ptr_array(nelements))
    WRITE(ichar,*) nelements
    nlength=LEN_TRIM(ADJUSTL(ichar))
    IF ( nlength > 27 ) THEN
       CALL finish('mo_histogram, add_multiple_stream_elements_2d', 'too many 2d variables')
    END IF
    DO i=1,nelements
      WRITE(ichar,*) i
      ilength=LEN_TRIM(ADJUSTL(ichar))
      DO j=ilength+1,nlength
         ichar='0'//TRIM(ADJUSTL(ichar))
      END DO
      CALL add_stream_element (stream,basename//TRIM(ADJUSTL(ichar)),ptr_array(i)%v, &
                               leveltype=leveltype, lpost=lpost, laccu=laccu, reset=reset, &
                               contnorest=contnorest)
    END DO


  END SUBROUTINE add_multiple_stream_elements_2d

  SUBROUTINE init_histos

    TYPE (t_stream), POINTER      :: hists

    CALL new_stream (hists,'hists',lrerun=.false.)

    !--- Add standard fields for post-processing:

    CALL add_stream_reference (hists, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (hists, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (hists, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (hists, 'gboxarea','geoloc',lpost=.TRUE.)
    CALL add_stream_reference (hists, 'slm','g3b',lpost=.TRUE.)


    CALL default_stream_setting (hists, lpost     = .TRUE. , &
                                        lrerun    = .FALSE. , &
                                        leveltype = HYBRID , &
                                        table     = 199,     &
                                        laccu     = .TRUE.,  &
                                        ldivi     = .TRUE.,  &
                                        code      = AUTO     )

    ! essential histograms
    CALL add_multiple_stream_elements(hists, 'clf_cc', frac_spec%nbin, clf_cc)
    CALL add_multiple_stream_elements(hists, 'coliqf_cci', frac_spec%nbin, coliqf_cci)
    CALL add_multiple_stream_elements(hists, 'chetf_cci', frac_spec%nbin, chetf_cci)
    CALL add_multiple_stream_elements(hists, 'csosif_cc', frac_spec%nbin, csosif_cc)
    CALL add_multiple_stream_elements(hists, 'csocmf_cc', ten_spec%nbin, csocmf_cc)
    CALL add_multiple_stream_elements(hists, 'ccc', frac_spec%nbin, ccc)
    CALL add_multiple_stream_elements(hists, 'ct_lf', mxt_spec%nbin, ct_lf)
    CALL add_multiple_stream_elements(hists, 'ct_lf_night', mxt_spec%nbin, ct_lf_night)
    CALL add_multiple_stream_elements(hists, 'ct_cc', mxt_spec%nbin, ct_cc)
    CALL add_multiple_stream_elements(hists, 'ct_cci', mxt_spec%nbin, ct_cci)
    CALL add_multiple_stream_elements(hists, 'ct_cc_night', mxt_spec%nbin, ct_cc_night)
    CALL add_multiple_stream_elements(hists, 'ctlf_cc', mxt_spec%nbin*frac_spec%nbin, ctlf_cc)
    CALL add_multiple_stream_elements(hists, 'ctlf_cc_night', mxt_spec%nbin*frac_spec%nbin, ctlf_cc_night)
    CALL add_multiple_stream_elements(hists, 'ct_hetf', mxt_spec%nbin, ct_hetf)
    CALL add_multiple_stream_elements(hists, 'ct_oliqf', mxt_spec%nbin, ct_oliqf)
    CALL add_multiple_stream_elements(hists, 'clwp_ccl', wp_spec%nbin, clwp_ccl)
    CALL add_multiple_stream_elements(hists, 'ciwp_cci', wp_spec%nbin, ciwp_cci)
    CALL add_multiple_stream_elements(hists, 'crleff_ccl_top', radl_spec%nbin, crleff_ccl_top)
    CALL add_multiple_stream_elements(hists, 'crieff_cci_top', radi_spec%nbin, crieff_cci_top)

    CALL add_multiple_stream_elements(hists, 'cosp_ct_lf', 40, cosp_ct_lf)
    CALL add_multiple_stream_elements(hists, 'cosp_ct_cc', 40, cosp_ct_cc)
    CALL add_multiple_stream_elements(hists, 'cosp_ct_lf_time', 40, cosp_ct_lf_time)
    CALL add_multiple_stream_elements(hists, 'cosp_ct_cc_time', 40, cosp_ct_cc_time)

    CALL add_stream_element(hists, 'cld_time', cld_time, &
                            longname='Time of cloud occurence', units='s')
    
    ! auxilliary histograms ONLY OUTPUT IF lfull_diag == .TRUE.
    CALL default_stream_setting (hists, lpost=lfull_diag)
    CALL add_multiple_stream_elements(hists, 'ct_lf_day', mxt_spec%nbin, ct_lf_day)
    CALL add_multiple_stream_elements(hists, 'cttop_hetf', mxt_spec%nbin, cttop_hetf)
    CALL add_multiple_stream_elements(hists, 'cxi', mmr_spec%nbin, cxi)
    CALL add_multiple_stream_elements(hists, 'cxl', mmr_spec%nbin, cxl)
    CALL add_multiple_stream_elements(hists, 'ct_cc_day', mxt_spec%nbin, ct_cc_day)
    CALL add_multiple_stream_elements(hists, 'cttop_cci', mxt_spec%nbin, cttop_cci)
    CALL add_multiple_stream_elements(hists, 'ctlf_cc_day', mxt_spec%nbin*frac_spec%nbin, ctlf_cc_day)
    CALL add_multiple_stream_elements(hists, 'ctlf_cci', mxt_spec%nbin*frac_spec%nbin, ctlf_cci)
    CALL add_multiple_stream_elements(hists, 'cttoplf_cc', mxt_spec%nbin*frac_spec%nbin, cttoplf_cc)
    CALL add_multiple_stream_elements(hists, 'ctlf_hetf', mxt_spec%nbin*frac_spec%nbin, ctlf_hetf)
    CALL add_multiple_stream_elements(hists, 'cprcp_tot', prcp_spec%nbin, cprcp_tot)
    CALL add_multiple_stream_elements(hists, 'cprcp_rain', prcp_spec%nbin, cprcp_rain)
    CALL add_multiple_stream_elements(hists, 'cprcp_snow', prcp_spec%nbin, cprcp_snow)
    CALL add_multiple_stream_elements(hists, 'cprcp_phase', frac_spec%nbin, cprcp_phase)
    CALL add_multiple_stream_elements(hists, 'cprcp_topphase', frac_spec%nbin, cprcp_topphase)
    CALL add_multiple_stream_elements(hists, 'cthetf_cci', mxt_spec%nbin*frac_spec%nbin, cthetf_cci)
    CALL add_multiple_stream_elements(hists, 'ctoliqf_cci', mxt_spec%nbin*frac_spec%nbin, ctoliqf_cci)
    CALL add_multiple_stream_elements(hists, 'cttophetf_cci', mxt_spec%nbin*frac_spec%nbin, cttophetf_cci)
    CALL add_multiple_stream_elements(hists, 'cthetf_lf', mxt_spec%nbin*frac_spec%nbin, cthetf_lf)


  END SUBROUTINE init_histos
                

  ! maps each grid-point onto a bin specified by
  ! h_spec (number of bins, minimum value, maximum value)
  ! and counts number of occurences
  SUBROUTINE simple_histogram_3d(&
                !-- IN
                kproma, kbdim, klev, &
                h_spec, &
                data, &
                !-- OUT
                hist)

    INTEGER, INTENT(in) :: kproma, kbdim, klev            ! block parameter
    TYPE(t_hist), INTENT(in) :: h_spec                    ! histogram parameter, number, range(min,max)
    REAL(dp), INTENT(in), DIMENSION(kbdim,klev)   :: data ! data to be binned

    REAL(dp), INTENT(out), DIMENSION(h_spec%nbin,kbdim,klev) :: hist ! array containing the counts

    INTEGER, DIMENSION(kbdim,klev) :: indices           ! each data-element mapped to bins
    LOGICAL, DIMENSION(kbdim,klev) :: mask              ! mask for each bin
    LOGICAL, DIMENSION(kbdim,klev) :: bnds
    INTEGER :: i

    ! initialize 
    hist = 0
    bnds = .TRUE.

    ! map data to bin indices
    SELECT CASE(h_spec%SCALE)
       CASE(LINSCALE)
          indices = linear_index(kproma, kbdim, klev, data, h_spec)
       CASE(LOGSCALE)
          indices = log_index(kproma, kbdim, klev, data, h_spec)
       CASE DEFAULT
          indices = linear_index(kproma, kbdim, klev, data, h_spec)
    END SELECT

    ! check histogram bounds
    bnds(1:kproma,:) = data(1:kproma,:) < h_spec%vmin .OR. data(1:kproma,:) > h_spec%vmax    
    indices(1:kproma,:) = MERGE(-1, indices(1:kproma,:), bnds(1:kproma,:))

    DO i=1,h_spec%nbin
       ! sort the bins into the corresponding array
       mask(1:kproma,:) = indices(1:kproma,:) == i-1
       hist(i,1:kproma,:) = MERGE(1._dp, 0._dp, mask(1:kproma,:))
    END DO

  END SUBROUTINE simple_histogram_3d
  SUBROUTINE simple_histogram_2d(&
                !-- IN
                kproma, kbdim, &
                h_spec, &
                data, &
                !-- OUT
                hist)

    INTEGER, INTENT(in) :: kproma, kbdim                  ! block parameter
    TYPE(t_hist), INTENT(in) :: h_spec                    ! histogram parameter, number, range(min,max)
    REAL(dp), INTENT(in), DIMENSION(kbdim) :: data ! data to be binned

    REAL(dp), INTENT(out), DIMENSION(h_spec%nbin,kbdim) :: hist ! array containing the counts

    INTEGER, DIMENSION(kbdim) :: indices           ! each data-element mapped to bins
    LOGICAL, DIMENSION(kbdim) :: mask              ! mask for each bin
    LOGICAL, DIMENSION(kbdim) :: bnds
    INTEGER :: i

    ! initialize 
    hist = 0
    bnds = .TRUE.

    ! map data to bin indices
    SELECT CASE(h_spec%SCALE)
       CASE(LINSCALE)
          indices = linear_index(kproma, kbdim, data, h_spec)
       CASE(LOGSCALE)
          indices = log_index(kproma, kbdim, data, h_spec)
       CASE DEFAULT
          indices = linear_index(kproma, kbdim, data, h_spec)
    END SELECT

    ! check histogram bounds
    bnds(1:kproma) = data(1:kproma) < h_spec%vmin .OR. data(1:kproma) > h_spec%vmax    
    indices(1:kproma) = MERGE(-1, indices(1:kproma), bnds(1:kproma))

    DO i=1,h_spec%nbin
       ! sort the bins into the corresponding array
       mask(1:kproma) = indices(1:kproma) == i-1
       hist(i,1:kproma) = MERGE(1._dp, 0._dp, mask(1:kproma))
    END DO

  END SUBROUTINE simple_histogram_2d

  ! creates a Histogram for the grid nbin1 x nbin2
  ! for corresponding data1 and data2
  SUBROUTINE double_histogram_3d(&
             !--IN
             kproma, kbdim, klev, &
             h_spec1, h_spec2, &
             data1, data2, &
             !--OUT
             hist)

    INTEGER, INTENT(in) :: kproma, kbdim, klev            ! block parameter
    TYPE(t_hist), INTENT(in) :: h_spec1, h_spec2          ! histogram parameter, number, range(min,max)
    REAL(dp), INTENT(in), DIMENSION(kbdim,klev) :: data1, data2 ! data to be binned

    REAL(dp), INTENT(out), DIMENSION(h_spec1%nbin*h_spec2%nbin,kbdim,klev) :: hist ! array containing the counts

    INTEGER, DIMENSION(kbdim,klev) :: indices, indices1, indices2 ! each data-element mapped to bins
    LOGICAL, DIMENSION(kbdim,klev) :: mask               ! mask for each bin
    LOGICAL, DIMENSION(kbdim,klev) :: bnds
    INTEGER :: i

    ! initialize 
    hist = 0
    bnds = .TRUE.

    ! map data to bin indices
    SELECT CASE(h_spec1%SCALE)
       CASE(LINSCALE)
          indices1 = linear_index(kproma, kbdim, klev, data1, h_spec1)
       CASE(LOGSCALE)
          indices1 = log_index(kproma, kbdim, klev, data1, h_spec1)
       CASE DEFAULT
          indices1 = linear_index(kproma, kbdim, klev, data1, h_spec1)
    END SELECT
    SELECT CASE(h_spec2%SCALE)
       CASE(LINSCALE)
          indices2 = linear_index(kproma, kbdim, klev, data2, h_spec2)
       CASE(LOGSCALE)
          indices2 = log_index(kproma, kbdim, klev, data2, h_spec2)
       CASE DEFAULT
          indices2 = linear_index(kproma, kbdim, klev, data2, h_spec2)
    END SELECT

    ! check histogram bounds
    bnds(1:kproma,:) = data1(1:kproma,:) < h_spec1%vmin .OR. data1(1:kproma,:) > h_spec1%vmax .OR. &
                       data2(1:kproma,:) < h_spec2%vmin .OR. data2(1:kproma,:) > h_spec2%vmax

    ! make one dimensional index array.
    indices(1:kproma,:) = indices1(1:kproma,:) + indices2(1:kproma,:)*h_spec1%nbin
    indices(1:kproma,:) = MERGE(-1, indices(1:kproma,:), bnds(1:kproma,:))

    ! assign bins and sum up vertically
    DO i=1,h_spec1%nbin*h_spec2%nbin
       ! sort the bins into the corresponding array
       mask(1:kproma,:) = indices(1:kproma,:) == i-1
       hist(i,1:kproma,:) = MERGE(1._dp, 0._dp, mask(1:kproma,:))
    END DO
END SUBROUTINE double_histogram_3d
SUBROUTINE double_histogram_2d(&
             !--IN
             kproma, kbdim, &
             h_spec1, h_spec2, &
             data1, data2, &
             !--OUT
             hist)

    INTEGER, INTENT(in) :: kproma, kbdim            ! block parameter
    TYPE(t_hist), INTENT(in) :: h_spec1, h_spec2          ! histogram parameter, number, range(min,max)
    REAL(dp), INTENT(in), DIMENSION(kbdim) :: data1, data2 ! data to be binned

    REAL(dp), INTENT(out), DIMENSION(h_spec1%nbin*h_spec2%nbin,kbdim) :: hist ! array containing the counts

    INTEGER, DIMENSION(kbdim) :: indices, indices1, indices2 ! each data-element mapped to bins
    LOGICAL, DIMENSION(kbdim) :: mask               ! mask for each bin
    LOGICAL, DIMENSION(kbdim) :: bnds
    INTEGER :: i

    ! initialize 
    hist = 0
    bnds = .TRUE.

    ! map data to bin indices
    SELECT CASE(h_spec1%SCALE)
       CASE(LINSCALE)
          indices1 = linear_index(kproma, kbdim, data1, h_spec1)
       CASE(LOGSCALE)
          indices1 = log_index(kproma, kbdim, data1, h_spec1)
       CASE DEFAULT
          indices1 = linear_index(kproma, kbdim, data1, h_spec1)
    END SELECT
    SELECT CASE(h_spec2%SCALE)
       CASE(LINSCALE)
          indices2 = linear_index(kproma, kbdim, data2, h_spec2)
       CASE(LOGSCALE)
          indices2 = log_index(kproma, kbdim, data2, h_spec2)
       CASE DEFAULT
          indices2 = linear_index(kproma, kbdim, data2, h_spec2)
    END SELECT

    ! check histogram bounds
    bnds(1:kproma) = data1(1:kproma) < h_spec1%vmin .OR. data1(1:kproma) > h_spec1%vmax .OR. &
                     data2(1:kproma) < h_spec2%vmin .OR. data2(1:kproma) > h_spec2%vmax

    ! make one dimensional index array.
    indices(1:kproma) = indices1(1:kproma) + indices2(1:kproma)*h_spec1%nbin
    indices(1:kproma) = MERGE(-1, indices(1:kproma), bnds(1:kproma))

    ! assign bins and sum up vertically
    DO i=1,h_spec1%nbin*h_spec2%nbin
       ! sort the bins into the corresponding array
       mask(1:kproma) = indices(1:kproma) == i-1
       hist(i,1:kproma) = MERGE(1._dp, 0._dp, mask(1:kproma))
    END DO
END SUBROUTINE double_histogram_2d

! maps data onto the discrete grid (bins) defined by h_spec
FUNCTION linear_index_3d(kproma, kbdim, klev, data, h_spec) RESULT(indices)

  INTEGER :: kproma, kbdim, klev
  REAL(dp), DIMENSION(kbdim,klev) :: data
  TYPE(t_hist) :: h_spec

  INTEGER, DIMENSION(kbdim,klev) :: indices
  REAL(dp) :: hw

  ! initialize
  indices = -1
  
  ! calculate bin width
  hw = (h_spec%vmax-h_spec%vmin)/h_spec%nbin

  ! map data to bins. Note that this way the case data == vmax is lost.
  indices(1:kproma,:) = INT((data(1:kproma,:)-h_spec%vmin)/hw)
  ! make the interval [vmin, vmax) inclusive --> [vmin, vmax]
  indices(1:kproma,:) = MERGE(h_spec%nbin-1, indices(1:kproma,:), data(1:kproma,:) == h_spec%vmax)

END FUNCTION linear_index_3d
FUNCTION linear_index_2d(kproma, kbdim, data, h_spec) RESULT(indices)

  INTEGER :: kproma, kbdim
  REAL(dp), DIMENSION(kbdim) :: data
  TYPE(t_hist) :: h_spec

  INTEGER, DIMENSION(kbdim) :: indices
  REAL(dp) :: hw

  ! initialize
  indices = -1
  
  ! calculate bin width
  hw = (h_spec%vmax-h_spec%vmin)/h_spec%nbin

  ! map data to bins. Note that this way the case data == vmax is lost.
  indices(1:kproma) = INT((data(1:kproma)-h_spec%vmin)/hw)
  ! make the interval [vmin, vmax) inclusive --> [vmin, vmax]
  indices(1:kproma) = MERGE(h_spec%nbin-1, indices(1:kproma), data(1:kproma) == h_spec%vmax)

END FUNCTION linear_index_2d

FUNCTION log_index_3d(kproma, kbdim, klev, data, h_spec) RESULT(indices)

  INTEGER :: kproma, kbdim, klev
  REAL(dp), DIMENSION(kbdim,klev) :: data
  TYPE(t_hist) :: h_spec

  INTEGER, DIMENSION(kbdim,klev) :: indices
  REAL(dp) :: hw

  ! initialize
  indices = -1
  
  ! calculate bin width in log-space
  hw = (LOG10(h_spec%vmax/h_spec%vmin))/h_spec%nbin

  ! map data to bins. Note that this way the case data == vmax is lost.
  indices(1:kproma,:) = INT(LOG10(MAX(data(1:kproma,:)/h_spec%vmin, 1._dp))/hw)
  ! make the interval [vmin, vmax) inclusive --> [vmin, vmax]
  indices(1:kproma,:) = MERGE(h_spec%nbin-1, indices(1:kproma,:), data(1:kproma,:) == h_spec%vmax)

END FUNCTION log_index_3d
FUNCTION log_index_2d(kproma, kbdim, data, h_spec) RESULT(indices)

  INTEGER :: kproma, kbdim
  REAL(dp), DIMENSION(kbdim) :: data
  TYPE(t_hist) :: h_spec

  INTEGER, DIMENSION(kbdim) :: indices
  REAL(dp) :: hw

  ! initialize
  indices = 0
  
  ! calculate bin width in log-space
  hw = (LOG10(h_spec%vmax/h_spec%vmin))/h_spec%nbin

  ! map data to bins. Note that this way the case data == vmax is lost.
  indices(1:kproma) = INT(LOG10(MAX(data(1:kproma)/h_spec%vmin, 1._dp))/hw)
  ! make the interval [vmin, vmax) inclusive --> [vmin, vmax]
  indices(1:kproma) = MERGE(h_spec%nbin-1, indices(1:kproma), data(1:kproma) == h_spec%vmax)

END FUNCTION log_index_2d

SUBROUTINE histogram_diagnostic(&
              !--IN
              kproma, kbdim, klev, klevp1, ktrac, ktdia, krow, &
              paphm1, papm1, ptvm1, pgeo, paclc, &
              pt, pq, pxl, pxi, pxt, pssfl, prsfl)

  USE mo_cloud_utils, ONLY: clc_min, cloud_type_helper, get_cloud_bounds, get_util_var, &
                            get_precip_fraction, eps
  USE mo_time_control, ONLY: delta_time
  USE mo_echam_cloud_params, ONLY: cqtmin
  USE mo_p3_fields, ONLY: idt_qihet, idt_qirim, idt_birim, idt_qsrc, idt_qprc, idt_qioliq
  USE mo_activ, ONLY: idt_icnc, idt_cdnc
  USE mo_geoloc, ONLY: rdayl_x
  USE mo_physical_constants, ONLY : rd

  INTEGER, INTENT(IN)  :: kproma, kbdim, klev, klevp1, ktrac, ktdia, krow
  REAL(dp), INTENT(IN) :: paphm1(kbdim,klevp1), papm1(kbdim,klev), ptvm1(kbdim,klev)
  REAL(dp), INTENT(IN) :: paclc(kbdim,klev) !< cloud cover (1)
  REAL(dp), INTENT(IN) :: pgeo(kbdim,klev)  !< geopotential height at full levels
  REAL(dp), INTENT(IN) :: pt(kbdim,klev)    !< temperature (K)
  REAL(dp), INTENT(IN) :: pq(kbdim,klev)    !< humidity mmr (1)
  REAL(dp), INTENT(IN) :: pxl(kbdim,klev)   !< liquid water mmr (1)
  REAL(dp), INTENT(IN) :: pxi(kbdim,klev)   !< ice water mmr (1)
  REAL(dp), INTENT(IN) :: pxt(kbdim,klev, ktrac)   !< tracers
  
  REAL(dp), INTENT(IN) :: pssfl(kbdim)      !< snow surface flux (kg m-2 s-1)
  REAL(dp), INTENT(IN) :: prsfl(kbdim)      !< rain surface flux (kg m-2 s-1)

  REAL(dp) :: zlf(kbdim,klev)               !< liquid fraction (1)
  REAL(dp) :: zhetf(kbdim,klev)             !< heterogeneously formed ice fraction (1)
  REAL(dp) :: zoliqf(kbdim,klev)            !< heterogeneously formed ice fraction (1)
  REAL(dp) :: zsosif(kbdim,klev)            !< source sink ratio (1)
  REAL(dp) :: zsocmf(kbdim,klev)            !< source to cloud mass ratio (1)
  REAL(dp) :: zrieff(kbdim,klev)            !< effective ice radius (um)
  REAL(dp) :: zrleff(kbdim,klev)            !< effective liquid radius (um)
  LOGICAL  :: ll_cc(kbdim,klev)             !< cloud flag
  REAL(dp) :: zxlb(kbdim,klev)              !< in-cloud liquid water mmr (1)
  REAL(dp) :: zxib(kbdim,klev)              !< in-cloud liquid water mmr (1)
  REAL(dp) :: zlwp(kbdim)                   !< column integrated liquid water content (g m-2)
  REAL(dp) :: ziwp(kbdim)                   !< column integrated ice water content (g m-2)
  REAL(dp) :: zaclci(kbdim,klev)            !< ice cloud cover (precip) (1)
  REAL(dp) :: zdtime                        !< time step length (s)

! ----------------------------GRIDBOX DESCRIPTION---------------------------------
  REAL(dp) :: zgeoh(kbdim,klevp1)  !< Geopotential height at half levels
  REAL(dp) :: zdz(kbdim,klev)      !< layer thickness [m]
  REAL(dp) :: zdp(kbdim,klev)      !< pressure difference of the layer [Pa]
  REAL(dp) :: zdpg(kbdim,klev)     !< delta p over g [kg/m2]
  REAL(dp) :: zaaa(kbdim,klev)     !< Air density correction needed for the ice crystal fall velocity
  REAL(dp) :: zviscos(kbdim,klev)  !< Dynamic viscosity of water in air
  REAL(dp) :: zrho(kbdim,klev)     !< air density

! --------------------------------HISTOGRAMS------------------------------------
  REAL(dp) :: zclf(frac_spec%nbin,kbdim,klev) !< supercooled liquid fraction histogram
  REAL(dp) :: zchetf(frac_spec%nbin,kbdim,klev) !< heterogeneous origin histogram
  REAL(dp) :: zcoliqf(frac_spec%nbin,kbdim,klev) !< liquid origin ice histogram
  REAL(dp) :: zcsosif(frac_spec%nbin,kbdim,klev) !< source sind ratio histogram
  REAL(dp) :: zcsocmf(ten_spec%nbin,kbdim,klev) !< source sind ratio histogram
  REAL(dp) :: zccc(frac_spec%nbin,kbdim,klev) !< cloud fraction histogram
  REAL(dp) :: zct(mxt_spec%nbin,kbdim,klev)  !< temperature histogram of all grid-boxes
  REAL(dp) :: zcttop(mxt_spec%nbin,kbdim,klev)  !< cloud top temp. of all grid-boxes
  REAL(dp) :: zcxi(mmr_spec%nbin,kbdim,klev)  !< xi histogram
  REAL(dp) :: zcxl(mmr_spec%nbin,kbdim,klev)  !< xl histogram
  REAL(dp) :: zctlf(mxt_spec%nbin*frac_spec%nbin,kbdim,klev)  !< supercooled liquid fraction avg by temperature
  REAL(dp) :: zcttoplf(mxt_spec%nbin*frac_spec%nbin,kbdim,klev)  !< supercooled liquid fraction avg by temperature
  REAL(dp) :: zcthetf(mxt_spec%nbin*frac_spec%nbin,kbdim,klev) !< het. formed ice fraction by temperature
  REAL(dp) :: zctoliqf(mxt_spec%nbin*frac_spec%nbin,kbdim,klev) !< het. formed ice fraction by temperature
  REAL(dp) :: zcttophetf(mxt_spec%nbin*frac_spec%nbin,kbdim,klev) !< het. formed ice fraction by cloud top temperature
  REAL(dp) :: zcprcp_tot(prcp_spec%nbin,kbdim)  !< precipitation histogram (tot)
  REAL(dp) :: zcprcp_rain(prcp_spec%nbin,kbdim)  !< precipitation histogram (rain)
  REAL(dp) :: zcprcp_snow(prcp_spec%nbin,kbdim)  !< precipitation histogram (snow)
  REAL(dp) :: zclwp(wp_spec%nbin,kbdim)          !< integrated liquid water histogram
  REAL(dp) :: zciwp(wp_spec%nbin,kbdim)          !< integrated ice water histogram
  REAL(dp) :: zcrleff(radl_spec%nbin,kbdim,klev)  !< ice effective radius histogram
  REAL(dp) :: zcrieff(radi_spec%nbin,kbdim,klev)  !< ice effective radius histogram
  
  REAL(dp), DIMENSION(kbdim,klev) :: ztmp1_2d, ztmp2_2d, ztmp3_2d, ztmp4_2d
  REAL(dp), DIMENSION(kbdim)      :: ztmp1, ztmp2
  LOGICAL, DIMENSION(kbdim,klev)  :: ll1_2d, ll2_2d, ll3_2d
  REAL(dp), DIMENSION(kbdim,klev) :: zcc_flag, zcci_flag, zccl_flag, zccil_flag, zcctop_flag
  INTEGER :: jk

! --------------------------------CLOUD TYPES-----------------------------------
  REAL(dp) :: zdp_cld(kbdim,klev) ! cloud pressure difference (Pa)
  REAL(dp) :: zt_top(kbdim,klev)  ! cloud top temperature (K)
  REAL(dp) :: zp_top(kbdim,klev)  ! cloud top pressure (Pa)
  INTEGER  :: zalltop(kbdim,klev) ! cloud top index
  INTEGER  :: zallbas(kbdim,klev) ! cloud base index
  INTEGER  :: itop(kbdim,klev),         & !< flag for cloud tops
              ibas(kbdim,klev),         & !< flag for cloud bases
              icl_minusbas(kbdim,klev), & !< flag for all cloud levels excepted their base
              icl_minustop(kbdim,klev), & !<  flag for all cloud levels excepted their top (useless for now)
              iclnb(kbdim)                !< number of clouds per column

!-------------------------------------------------------------------------------------
!                                      CLOUD TYPES
!-------------------------------------------------------------------------------------
ll1_2d(1:kproma,:) = (pxl(1:kproma,:)+pxi(1:kproma,:)) > cqtmin
ztmp1_2d(1:kproma,:) = MERGE(paclc(1:kproma,:), 0._dp, ll1_2d(1:kproma,:))

CALL get_cloud_bounds( &
        !-- IN
        kproma, kbdim, ktdia, klev, ztmp1_2d, &
        !- OUT
        itop, ibas, icl_minustop, icl_minusbas, iclnb)

CALL cloud_type_helper(&
        !--IN
        kproma, kbdim, ktdia, klev, klevp1, &
        ztmp1_2d, itop, ibas, pt, paphm1, &
        !--OUT
        zdp_cld, zt_top, zp_top, zalltop, zallbas)

!-------------------------------------------------------------------------------------
!                                  UTILITY VARIABLES
!-------------------------------------------------------------------------------------

!--- Get several utility variables:
  CALL get_util_var( &
          !-- IN
          kproma, kbdim, ktdia, klev, klevp1, &
          paphm1(:,:), pgeo(:,:), papm1(:,:), pt(:,:), &
          !-- OUT
          zgeoh(:,:), zdp(:,:), zdpg(:,:), &
          zdz(:,:), zaaa(:,:), zviscos(:,:) )
  zrho(1:kproma,:)     = papm1(1:kproma,:)/(rd*ptvm1(1:kproma,:))

  CALL get_precip_fraction(kproma, kbdim, ktdia, klev, paclc, zaclci)

! liquid fraction
ztmp1_2d(1:kproma,:) = pxi(1:kproma,:) + pxl(1:kproma,:)
zlf(1:kproma,:) = MERGE(pxl(1:kproma,:)/ztmp1_2d(1:kproma,:), 0._dp, &
                        ztmp1_2d(1:kproma,:) > eps)
zlf(1:kproma,:) = MAX(MIN(zlf(1:kproma,:), 1._dp), 0._dp)

! fraction of ice that formed through frozen liquid water
zoliqf(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qioliq)/pxi(1:kproma,:), 0._dp, &
                           pxi(1:kproma,:) > eps)
zoliqf(1:kproma,:) = MAX(MIN(1._dp, zoliqf(1:kproma,:)), 0._dp)

! heterogeneously formed ice fraction
zhetf(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qihet)/pxi(1:kproma,:), 0._dp, &
                          pxi(1:kproma,:) > eps)
zhetf(1:kproma,:) = MAX(MIN(zhetf(1:kproma,:), 1._dp), 0._dp)

! fraction of source and sink terms
ll1_2d(1:kproma,:) = pxt(1:kproma,:,idt_qsrc) > eps
zsosif(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qprc)/pxt(1:kproma,:,idt_qsrc), 0._dp, &
                           pxt(1:kproma,:,idt_qsrc) > eps)
zsosif(1:kproma,:) = MAX(0._dp, MIN(1._dp, zsosif(1:kproma,:)))

! fraction of source and cloud mass
ztmp1_2d(1:kproma,:) = pxi(1:kproma,:) + pxl(1:kproma,:)
zsocmf(1:kproma,:) = MERGE(pxt(1:kproma,:,idt_qsrc)/ztmp1_2d(1:kproma,:), 0._dp, &
                           ztmp1_2d(1:kproma,:) > eps)
zsocmf(1:kproma,:) = MAX(ten_spec%vmin, MIN(ten_spec%vmax, zsocmf(1:kproma,:)))

ztmp1_2d(:,:) = pxt(:,:,idt_qirim)
ztmp2_2d(:,:) = pxt(:,:,idt_icnc)
ztmp3_2d(:,:) = pxt(:,:,idt_birim)
ztmp4_2d(:,:) = pxt(:,:,idt_cdnc)
zrleff(1:kproma,:) = effective_liquid_radius(kbdim, kproma, klev, ztmp4_2d, pxl)
DO jk=1,klev
   zrieff(:,jk) = 1.e6_dp*effective_ice_radius(kbdim, kproma, pxi(:,jk), ztmp1_2d(:,jk), &
                                                         ztmp2_2d(:,jk), ztmp3_2d(:,jk))
END DO


!-------------------------------------------------------------------------------------
!                                      HISTOGRAMS
!-------------------------------------------------------------------------------------

ll_cc(1:kproma,:) = paclc(1:kproma,:) > clc_min
zdtime = delta_time

! compute in-cloud quantities
ztmp1_2d(1:kproma,:) = MERGE(1._dp/paclc(1:kproma,:), 0._dp, ll_cc(1:kproma,:))
zxib(1:kproma,:) = pxi(1:kproma,:)*ztmp1_2d(1:kproma,:)
zxlb(1:kproma,:) = pxi(1:kproma,:)*ztmp1_2d(1:kproma,:)

! compute column integrated quantities
zlwp(1:kproma) = SUM(pxl(1:kproma,:)*zdpg(1:kproma,:), 2)*1.e3_dp
ziwp(1:kproma) = SUM(pxi(1:kproma,:)*zdpg(1:kproma,:), 2)*1.e3_dp

ztmp2_2d(1:kproma,:) = pxl(1:kproma,:) + pxi(1:kproma,:)
ll1_2d(1:kproma,:) = ztmp2_2d(1:kproma,:) > cqtmin .AND. ll_cc(1:kproma,:)
ll2_2d(1:kproma,:) = pxi(1:kproma,:) > cqtmin .AND. ll_cc(1:kproma,:)
ll3_2d(1:kproma,:) = pxl(1:kproma,:) > cqtmin .AND. ll_cc(1:kproma,:)

! Used to average histograms over cloudy regions.
! E.g. the supercooled liquid fraction is computed in-cloud.
zcc_flag(1:kproma,:) = MERGE(1._dp, 0._dp, ll1_2d(1:kproma,:))
zcci_flag(1:kproma,:) = MERGE(1._dp, 0._dp, ll2_2d(1:kproma,:))
zccl_flag(1:kproma,:) = MERGE(1._dp, 0._dp, ll3_2d(1:kproma,:))
zccil_flag(1:kproma,:) = MERGE(1._dp, 0._dp, ll2_2d(1:kproma,:) .AND. ll3_2d(1:kproma,:))

cld_time(1:kproma,:,krow) = cld_time(1:kproma,:,krow) + zdtime*zcc_flag(1:kproma,:)

! compute cloud top flag
DO jk=1,klev
   zcctop_flag(1:kproma,jk) = MERGE(1._dp, 0._dp, itop(1:kproma,jk) == jk)
END DO

! construct histograms
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, frac_spec, zlf, &
          !-- OUT
          zclf)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, frac_spec, zoliqf, &
          !-- OUT
          zcoliqf)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, frac_spec, zhetf, &
          !-- OUT
          zchetf)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, frac_spec, zsosif, &
          !-- OUT
          zcsosif)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, ten_spec, zsocmf, &
          !-- OUT
          zcsocmf)
  
  CALL double_histogram(&
          !--IN
          kproma, kbdim, klev, mxt_spec, frac_spec, &
          pt, zlf, &
          !--OUT
          zctlf)
  CALL double_histogram(&
          !--IN
          kproma, kbdim, klev, mxt_spec, frac_spec, &
          zt_top, zlf, &
          !--OUT
          zcttoplf)
  CALL double_histogram(&
          !--IN
          kproma, kbdim, klev, mxt_spec, frac_spec, &
          pt, zhetf, &
          !--OUT
          zcthetf)
  CALL double_histogram(&
          !--IN
          kproma, kbdim, klev, mxt_spec, frac_spec, &
          pt, zoliqf, &
          !--OUT
          zctoliqf)
  CALL double_histogram(&
          !--IN
          kproma, kbdim, klev, mxt_spec, frac_spec, &
          zt_top, zhetf, &
          !--OUT
          zcttophetf)

  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, frac_spec, paclc, &
          !-- OUT
          zccc)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, mxt_spec, pt, &
          !-- OUT
          zct)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, mxt_spec, zt_top, &
          !-- OUT
          zcttop)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, mmr_spec, zxib, &
          !-- OUT
          zcxi)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, mmr_spec, zxlb, &
          !-- OUT
          zcxl)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, wp_spec, zlwp, &
          !-- OUT
          zclwp)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, wp_spec, ziwp, &
          !-- OUT
          zciwp)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, radl_spec, zrleff, &
          !-- OUT
          zcrleff)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, klev, radi_spec, zrieff, &
          !-- OUT
          zcrieff)

  ztmp1(:) = 0._dp
  ztmp1(1:kproma) = prsfl(1:kproma)+pssfl(1:kproma)
  ztmp2(:) = 0._dp
  ztmp2(1:kproma) = MERGE(1._dp/zaclci(1:kproma,klev), 0._dp, zaclci(1:kproma,klev) > clc_min)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, prcp_spec, ztmp1*3600._dp*ztmp2, &
          !-- OUT
          zcprcp_tot)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, prcp_spec, pssfl*3600._dp*ztmp2, &
          !-- OUT
          zcprcp_snow)
  CALL simple_histogram(&
          !-- IN
          kproma, kbdim, prcp_spec, prsfl*3600._dp*ztmp2, &
          !-- OUT
          zcprcp_rain)

! sum up histograms
  ! liquid fraction histogram
  DO jk=1,frac_spec%nbin
     ! liquid fraction
     clf_cc(jk)%v(1:kproma,:,krow) = clf_cc(jk)%v(1:kproma,:,krow) &
                                   + zdtime*zcc_flag(1:kproma,:)*zclf(jk,1:kproma,:)*zdpg(1:kproma,:)

     ! cloud fraction
     ccc(jk)%v(1:kproma,:,krow) = ccc(jk)%v(1:kproma,:,krow) + zdtime*zccc(jk,1:kproma,:)*zdpg(1:kproma,:)

     ! add up souce sink ratio
     ztmp1(1:kproma) = SUM(zcc_flag(1:kproma,:)*zcsosif(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     csosif_cc(jk)%v(1:kproma,krow) = csosif_cc(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! liquid origin fraction
     ztmp1(1:kproma) = SUM(zcci_flag(1:kproma,:)*zcoliqf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     coliqf_cci(jk)%v(1:kproma,krow) = coliqf_cci(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! heterogeneous origin fraction
     ztmp1(1:kproma) = SUM(zcci_flag(1:kproma,:)*zchetf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     chetf_cci(jk)%v(1:kproma,krow) = chetf_cci(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO

  ! ten histogram
  DO jk=1,ten_spec%nbin
     ! add up the source to cloud mass ratios
     ztmp1(1:kproma) = SUM(zcc_flag(1:kproma,:)*zcsocmf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     csocmf_cc(jk)%v(1:kproma,krow) = csocmf_cc(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO

  ! ice histogram
  DO jk=1,mmr_spec%nbin
     ! add up the occurences per water content (ice / liquid)
     cxi(jk)%v(1:kproma,:,krow) = cxi(jk)%v(1:kproma,:,krow) + zdtime*zcxi(jk,1:kproma,:)*zdpg(1:kproma,:)
     cxl(jk)%v(1:kproma,:,krow) = cxl(jk)%v(1:kproma,:,krow) + zdtime*zcxl(jk,1:kproma,:)*zdpg(1:kproma,:)
  END DO

  ! precip histogram
  DO jk=1,prcp_spec%nbin
     cprcp_tot(jk)%v(1:kproma,krow) = cprcp_tot(jk)%v(1:kproma,krow) + zdtime*zcprcp_tot(jk,1:kproma)
     cprcp_snow(jk)%v(1:kproma,krow) = cprcp_snow(jk)%v(1:kproma,krow) + zdtime*zcprcp_snow(jk,1:kproma)
     cprcp_rain(jk)%v(1:kproma,krow) = cprcp_rain(jk)%v(1:kproma,krow) + zdtime*zcprcp_rain(jk,1:kproma)
  END DO

  ! temperature histogram
  DO jk=1,mxt_spec%nbin
     ! count all clouds
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zcc_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ct_cc(jk)%v(1:kproma,krow) = ct_cc(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! count all clouds during day
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zcc_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ztmp1(1:kproma) = ztmp1(1:kproma)*rdayl_x(1:kproma,krow)
     ct_cc_day(jk)%v(1:kproma,krow) = ct_cc_night(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! count all clouds during night
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zcc_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ztmp1(1:kproma) = ztmp1(1:kproma)*(1._dp - rdayl_x(1:kproma,krow))
     ct_cc_night(jk)%v(1:kproma,krow) = ct_cc_night(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! count cloud containing ice
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zcci_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ct_cci(jk)%v(1:kproma,krow) = ct_cci(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
     
     ! clouds weighted by liquid fraction
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zlf(1:kproma,:)*zcc_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ct_lf(jk)%v(1:kproma,krow) = ct_lf(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! clouds weighted by liquid fraction during day
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zlf(1:kproma,:)*zcc_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ztmp1(1:kproma) = ztmp1(1:kproma)*rdayl_x(1:kproma,krow)
     ct_lf_day(jk)%v(1:kproma,krow) = ct_lf_day(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! clouds weighted by liquid fraction during night
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zlf(1:kproma,:)*zcc_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ztmp1(1:kproma) = ztmp1(1:kproma)*(1._dp - rdayl_x(1:kproma,krow))
     ct_lf_night(jk)%v(1:kproma,krow) = ct_lf_night(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! clouds weighted by het. formed ice fraction
     ! multiply by zcci_flag to mask ice outside clouds (aclc == 0)
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zhetf(1:kproma,:)*zcci_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ct_hetf(jk)%v(1:kproma,krow) = ct_hetf(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! clouds weighted by liquid origin ice fraction
     ! multiply by zcci_flag to mask ice outside clouds (aclc == 0)
     ztmp1(1:kproma) = SUM(zct(jk,1:kproma,:)*zoliqf(1:kproma,:)*zcci_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     ct_oliqf(jk)%v(1:kproma,krow) = ct_oliqf(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO

  ! cloud top temperature histogram
  DO jk=1,mxt_spec%nbin
     ! count all clouds containing ice
     ztmp1(1:kproma) = SUM(zcttop(jk,1:kproma,:)*zcci_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     cttop_cci(jk)%v(1:kproma,krow) = cttop_cci(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! clouds weighted by het. formed ice fraction
     ! multiply by zcci_flag to mask ice outside clouds (aclc == 0)
     ztmp1(1:kproma) = SUM(zcttop(jk,1:kproma,:)*zhetf(1:kproma,:)*zcci_flag(1:kproma,:)*zdpg(1:kproma,:), 2)
     cttop_hetf(jk)%v(1:kproma,krow) = cttop_hetf(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO
     
  ! temperature histogram - liquid fraction
  DO jk=1,frac_spec%nbin*mxt_spec%nbin
     ! count all clouds
     ztmp1(1:kproma) = SUM(zcc_flag(1:kproma,:)*zctlf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     ctlf_cc(jk)%v(1:kproma,krow) = ctlf_cc(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! count all clouds during day
     ztmp1(1:kproma) = SUM(zcc_flag(1:kproma,:)*zctlf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     ztmp1(1:kproma) = ztmp1(1:kproma)*rdayl_x(1:kproma,krow)
     ctlf_cc_day(jk)%v(1:kproma,krow) = ctlf_cc_day(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! count all clouds during night
     ztmp1(1:kproma) = SUM(zcc_flag(1:kproma,:)*zctlf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     ztmp1(1:kproma) = ztmp1(1:kproma)*(1._dp - rdayl_x(1:kproma,krow))
     ctlf_cc_night(jk)%v(1:kproma,krow) = ctlf_cc_night(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! count clouds containing ice
     ztmp1(1:kproma) = SUM(zcci_flag(1:kproma,:)*zctlf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     ctlf_cci(jk)%v(1:kproma,krow) = ctlf_cci(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
     
     ! clouds weighted by het. formed ice fraction
     ! multiply by zcci_flag to mask ice outside clouds
     ztmp1(1:kproma) = SUM(zcci_flag(1:kproma,:)*zhetf(1:kproma,:)*zctlf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     ctlf_hetf(jk)%v(1:kproma,krow) = ctlf_hetf(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO

  ! cloud top temperature - liquid fraction histogram
  DO jk=1,frac_spec%nbin*mxt_spec%nbin
     ! cout all clouds
     ztmp1(1:kproma) = SUM(zcc_flag(1:kproma,:)*zcttoplf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     cttoplf_cc(jk)%v(1:kproma,krow) = cttoplf_cc(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO

  ! heterogeneously formed ice fraction - temperature histogram
  DO jk=1,frac_spec%nbin*mxt_spec%nbin
     ! count all clouds containing ice
     ztmp1(1:kproma) = SUM(zcci_flag(1:kproma,:)*zcthetf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     cthetf_cci(jk)%v(1:kproma,krow) = cthetf_cci(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! count all clouds containing ice
     ztmp1(1:kproma) = SUM(zcci_flag(1:kproma,:)*zctoliqf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     ctoliqf_cci(jk)%v(1:kproma,krow) = ctoliqf_cci(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)

     ! count weighted by liquid fraction
     ! multiply by zcci_flag to mask ice outside clouds
     ztmp1(1:kproma) = SUM(zcci_flag(1:kproma,:)*zcthetf(jk,1:kproma,:)*zlf(1:kproma,:)*zdpg(1:kproma,:), 2)
     cthetf_lf(jk)%v(1:kproma,krow) = cthetf_lf(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO

  ! heterogeneously formed ice fraction - cloud top temperature histogram
  DO jk=1,frac_spec%nbin*mxt_spec%nbin
     ! count all clouds containing ice
     ztmp1(1:kproma) = SUM(zcci_flag(1:kproma,:)*zcttophetf(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     cttophetf_cci(jk)%v(1:kproma,krow) = cttophetf_cci(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO

  ! water path histogram
  DO jk=1,wp_spec%nbin
     ! liquid water, only count if liquid water is present in column
     ztmp1(1:kproma) = MERGE(1._dp, 0._dp, SUM(zccl_flag(1:kproma,:), 2) > 0.5_dp)
     clwp_ccl(jk)%v(1:kproma,krow) = clwp_ccl(jk)%v(1:kproma,krow) + zdtime*zclwp(jk,1:kproma)*ztmp1(1:kproma)

     ! ice water, only count if ice water is present in column
     ztmp1(1:kproma) = MERGE(1._dp, 0._dp, SUM(zcci_flag(1:kproma,:), 2) > 0.5_dp)
     ciwp_cci(jk)%v(1:kproma,krow) = ciwp_cci(jk)%v(1:kproma,krow) + zdtime*zciwp(jk,1:kproma)*ztmp1(1:kproma)
  END DO

  ! liquid radius histogram
  DO jk=1,radl_spec%nbin
     ! effective droplet radius
     ztmp1(1:kproma) = SUM(zcctop_flag(1:kproma,:)*zccl_flag(1:kproma,:)*zcrleff(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     crleff_ccl_top(jk)%v(1:kproma,krow) = crleff_ccl_top(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO

  ! ice radius histogram
  DO jk=1,radi_spec%nbin
     ! effective crystal radius
     ztmp1(1:kproma) = SUM(zcctop_flag(1:kproma,:)*zcci_flag(1:kproma,:)*zcrieff(jk,1:kproma,:)*zdpg(1:kproma,:), 2)
     crieff_cci_top(jk)%v(1:kproma,krow) = crieff_cci_top(jk)%v(1:kproma,krow) + zdtime*ztmp1(1:kproma)
  END DO
  
END SUBROUTINE histogram_diagnostic

END MODULE mo_histogram

    


