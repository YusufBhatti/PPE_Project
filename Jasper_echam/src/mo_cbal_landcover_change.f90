!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
module mo_cbal_landcover_change

  !! Christian Reick,    2006-08-15
  !! Christian Reick,    2009-06-03
  !! Christian Reick,    2009-10-++
  !! Stiig Wilkenskjeld, 2012-05-11
  !! Daniel Goll,        May 2013    yasso implementation

  use mo_kind,          only: dp
  use mo_mpi,           only: p_parallel_io
  use mo_exception,     only: finish, message, message_text
  use mo_linked_list,   only: t_stream
  use mo_jsbach,        only: new_year, new_day, debug, debug_Cconservation
  use mo_cbal_parameters, only: frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos, &
                                frac_mobile_2_atmos, frac_harvest_2_atmos

  implicit none

  type landcover_change_type
     ! Flux variables [mol m-2(grid box)s-1]
     real(dp),pointer :: LCC_flux_box_C2atmos(:)            !! Carbon flux to atmosphere from landcover change
     real(dp),pointer :: LCC_flux_box_C2litterGreenPools(:) !! Carbon flux from green and reserve pool to fast soil pool by lcc
     real(dp),pointer :: LCC_flux_box_C2litterWoodPool(:)   !! Carbon flux from wood pool to woody litter pool by lcc
     real(dp),pointer :: LCC_flux_box_N2atmos(:)            !! Nitrogen released to atmosphere by landcover changes
     real(dp),pointer :: LCC_flux_box_N2litterGreenPools(:) !! Nitrogen released from green and reserve pool to fast soil pool by
                                                            !!    landcover change (summed through whole output interval)
     real(dp),pointer :: LCC_flux_box_N2litterWoodPool(:)   !! Nitrogen released from wood pool to woody litter pool by landcover
                                                            !!    change (summed through whole output interval)
     real(dp),pointer :: LCC_flux_box_N2SMINNpool(:)        !! surplus nitrogen from wood stubbing that because of higher 
                                                            !!    C/N-ratio of woody litter as compared to living wood, is fastly 
                                                            !!    freed into the soil mineral N pool.
     ! Variables for conservation checks [mol m-2(grid box)]
     real(dp),pointer :: LCC_testCconserv(:)                !! For testing carbon conservation (should be zero if conserved)
     real(dp),pointer :: LCC_testNconserv(:)                !! For testing nitrogen conservation (should be zero if conserved)
     ! Various variables
     real(dp),pointer :: LCC_coverFract_target(:,:)         !! Fraction of vegetated part of a gridbox that due to landcover change
                                                            !!    should be reached at the end of the year. (nland x ntiles)
     real(dp),pointer :: C2atmos(:)                         !! Carbon flux to atmosphere from land cover change

     ! Variables below applies only to use of the anthropogenic 
     ! C-pools in accordance with the grand slam protocol (Houghton et al. 1983)

     ! Instant LCC induced fluxes to anthropogenic pools [mol m-2(veg) s-1]
     real(dp),pointer :: C_2_onSite(:)                      !! Carbon flux from LCC to onsite
     real(dp),pointer :: C_2_paper(:)                       !! Carbon flux from LCC to paper
     real(dp),pointer :: C_2_construction(:)                !! Carbon flux from LCC to construction
     real(dp),pointer :: C_2_paper_harv(:)                  !! Carbon flux from LCC to paper from harvest
     real(dp),pointer :: C_2_construction_harv(:)           !! Carbon flux from LCC to construction from harvest

     ! Instant fluxes from anthropogenic pools to atmosphere [mol m-2(veg) s-1]
     real(dp),pointer :: C_onSite_2_atmos(:)                !! Carbon flux from onSite to atmosphere
     real(dp),pointer :: C_paper_2_atmos(:)                 !! Carbon flux from paper to atmosphere
     real(dp),pointer :: C_construction_2_atmos(:)          !! Carbon flux from construction to atmosphere
     real(dp),pointer :: C_paper_harv_2_atmos(:)            !! Carbon flux from paper to atmosphere from harvest
     real(dp),pointer :: C_construction_harv_2_atmos(:)     !! Carbon flux from construction to atmosphere from harvest

     ! LCC induced fluxes to anthropogenic pools summed for temporal averages [mol m-2(grid box) s-1]
     real(dp),pointer :: boxC_2_onSite(:)                  !! Carbon flux from LCC to onsite green
     real(dp),pointer :: boxC_2_paper(:)                   !! Carbon flux from LCC to paper
     real(dp),pointer :: boxC_2_construction(:)            !! Carbon flux from LCC to construction
     real(dp),pointer :: boxC_2_paper_harv(:)              !! Carbon flux from LCC to paper from harvest
     real(dp),pointer :: boxC_2_construction_harv(:)       !! Carbon flux from LCC to construction from harvest

     ! Fluxes from anthropogenic pools to atmosphere summed for temporal averages [mol m-2(grid box) s-1]
     real(dp),pointer :: boxC_onSite_2_atmos(:)            !! Carbon flux from onSite to atmosphere
     real(dp),pointer :: boxC_paper_2_atmos(:)             !! Carbon flux from paper to atmosphere
     real(dp),pointer :: boxC_construction_2_atmos(:)      !! Carbon flux from construction to atmosphere
     real(dp),pointer :: boxC_paper_harv_2_atmos(:)        !! Carbon flux from paper to atmosphere from harvest
     real(dp),pointer :: boxC_construction_harv_2_atmos(:) !! Carbon flux from construction to atmosphere from harvest

  end type landcover_change_type

  type landuse_transitions_type
     !! Elements of the (reduced, 3x3) transition matrix of the New Hampshire Harmonized Landuse Protocol (NHHLP)
     !! .. as read in from file. Here only the landcover types crop (CROP), pasture (PAST), and natural vegetation (NATL) appear.
     real(dp),pointer :: TransMtrx_CROP_2_PAST(:) 
     real(dp),pointer :: TransMtrx_PAST_2_CROP(:)
     real(dp),pointer :: TransMtrx_NATL_2_PAST(:)
     real(dp),pointer :: TransMtrx_PAST_2_NATL(:)
     real(dp),pointer :: TransMtrx_CROP_2_NATL(:)
     real(dp),pointer :: TransMtrx_NATL_2_CROP(:)

     !! Additional elements of the (extended, 4x4) transition matrix of the New Hampshire protocol computed by the 
     !! .. JSBACH-implementation. The extra landcover types are: grasslands (GRAS) and Natural non-grassland vegetation (FRST)
     real(dp),pointer :: TransMtrx_FRST_2_PAST(:)
     real(dp),pointer :: TransMtrx_PAST_2_FRST(:)
     real(dp),pointer :: TransMtrx_GRAS_2_PAST(:)
     real(dp),pointer :: TransMtrx_PAST_2_GRAS(:)
     real(dp),pointer :: TransMtrx_FRST_2_CROP(:)
     real(dp),pointer :: TransMtrx_CROP_2_FRST(:)
     real(dp),pointer :: TransMtrx_GRAS_2_CROP(:)
     real(dp),pointer :: TransMtrx_CROP_2_GRAS(:)


     !! Cover fractions as parts of the whole grid box according to the NHHLP reached at the end of the last year for ..
     real(dp),pointer :: Grass_coverFract_lastYear(:)        !! grasslands (GRAS)
     real(dp),pointer :: NatWood_coverFract_lastYear(:)      !! natural non-grassland vegetation (FRST)
     real(dp),pointer :: Pasture_coverFract_lastYear(:)      !! pastures (PAST)
     real(dp),pointer :: Crop_coverFract_lastYear(:)         !! croplands (CROP)

     !! Test-arrays: They show wether the area tranfered by the tile-transitions is equal to the transitions between the 4 extended
     !! .. NHHLP transition types. If everything is correct, these arrays should be zero.  
     
     real(dp),pointer :: Test_NATL_2_PAST(:) 
     real(dp),pointer :: Test_PAST_2_NATL(:)
     real(dp),pointer :: Test_CROP_2_NATL(:)
     real(dp),pointer :: Test_NATL_2_CROP(:)

     real(dp),pointer :: Test_CROP_2_PAST(:) 
     real(dp),pointer :: Test_PAST_2_CROP(:)

     real(dp),pointer :: Test_FRST_2_PAST(:)
     real(dp),pointer :: Test_PAST_2_FRST(:)
     real(dp),pointer :: Test_GRAS_2_PAST(:)
     real(dp),pointer :: Test_PAST_2_GRAS(:)
     real(dp),pointer :: Test_FRST_2_CROP(:)
     real(dp),pointer :: Test_CROP_2_FRST(:)
     real(dp),pointer :: Test_GRAS_2_CROP(:)
     real(dp),pointer :: Test_CROP_2_GRAS(:)


     !! Diagnostics for inconsistency with New Hampshire protocol: cover fraction that could not be converted y 
     !! .. because of non-availability of potential GRAS or FRST area. Note that these arrays are
     !! .. are computed only once a year, i.e. these values indicate the ignored transition are
     !! .. per year. 

     real(dp),pointer :: CROP_2_NATL_ignored(:) 
     real(dp),pointer :: PAST_2_NATL_ignored(:)

     !! Harvest

     real(dp),pointer :: Box_harvest(:)             !! Prescribed harvest in grid box [mol(C)/m^2(gridbox)/s]]
     real(dp),pointer :: Box_flux_harvest(:)        !! Flux of total harvest in grid box [mol(C)/m^2(gridbox)/s]
     real(dp),pointer :: Box_flux_harvest_2atmos(:) !! Flux of harvest in grid box emitted directly as CO2 into 
                                                    !! .. the atmosphere [mol(C)/m^2(gridbox)/s]

     real(dp),pointer :: C2atmos_LUtrans(:)         !! Carbon flux to atmosphere from land use change
     real(dp),pointer :: C2atmos_harvest(:)         !! Carbon flux to atmosphere from harvest
     real(dp),pointer :: N2atmos_LUtrans(:)
     real(dp),pointer :: N2atmos_harvest(:)
  end type landuse_transitions_type

  TYPE(landcover_change_type)   , SAVE :: landcover_change
  TYPE(landuse_transitions_type), SAVE :: landuse_transitions

  public :: landcover_change
  public :: landcover_change_type
  public :: landuse_transitions
  public :: landuse_transitions_type
  public :: init_landcover_change
  public :: read_landcover_fractions
  public :: read_landuse_transitions
  public :: read_harvest
  public :: do_landcover_change
  public :: do_landuse_transitions

  private 

  REAL(dp), parameter :: sec_per_day            = 86400._dp    ! seconds per day
  character(len=*),parameter :: coverFractVarName = "cover_fract" !! name of variable for landcover fractions in netcdf input file

  !! --- private variables ---------------------------------------------------------

  TYPE(t_stream), POINTER, SAVE     :: LCC_stream, LCC_Nstream   !! Memory streams for model state

CONTAINS


  !! --- init_landcover_change() ----------------------------------------------------

  SUBROUTINE init_landcover_change(g_nland, l_nland, ntiles, isRestart, UseLanduseTransitions, with_nitrogen, &
                                   lcc_scheme, lcc, fileformat, fileztype, stream, nstream)
    USE mo_linked_list,  ONLY: LAND, TILES
    USE mo_output,       ONLY: veg_table
    USE mo_netCDF,       ONLY: max_dim_name
    USE mo_memory_base,  ONLY: new_stream,default_stream_setting, add =>add_stream_element
    USE mo_jsbach,       ONLY: missing_value

    INTEGER, INTENT(in)               :: g_nland, l_nland
    INTEGER, INTENT(in)               :: ntiles                !! number of tiles
    LOGICAL, INTENT(in)               :: isRestart             !! restart flag
    LOGICAL, INTENT(in)               :: UseLanduseTransitions !! landcover forcing by prescribing landuse transitions?
    LOGICAL, INTENT(in)               :: with_nitrogen         !! run jsbach with nitrogen?
    INTEGER, INTENT(in)               :: lcc_scheme            !! pool schema for lcc
    TYPE (landcover_change_type), INTENT(out) :: lcc
    INTEGER, INTENT(in)               :: fileformat            !! output file format
    INTEGER, INTENT(in)               :: fileztype             !! output file compression
    TYPE(t_stream), POINTER, OPTIONAL :: stream, Nstream       !! streams to which local streams can optionally be associated

    ! --- local variables

    INTEGER                     :: dim1p(1), dim1(1)
    INTEGER                     :: dim3p(2), dim3(2)
    CHARACTER(LEN=max_dim_name) :: dim1n(1), dim3n(2) 

    !! --- printout parameters 

    CALL message("init_landcover_change()","=== Parameters Landcover Change ============")
    WRITE (message_text,*)                 " ac_wood_2_atmos=",frac_wood_2_atmos
    CALL message("init_landcover_change()", message_text)
    WRITE (message_text,*)                 " frac_green_2_atmos=",frac_green_2_atmos
    CALL message("init_landcover_change()", message_text)
    WRITE (message_text,*)                 " frac_reserve_2_atmos=",frac_reserve_2_atmos
    CALL message("init_landcover_change()", message_text)
    WRITE (message_text,*)                 " frac_mobile_2_atmos=",frac_mobile_2_atmos
    CALL message("init_landcover_change()", message_text)
    WRITE (message_text,*)                 " frac_harvest_2_atmos=",frac_harvest_2_atmos
    CALL message("init_landcover_change()", message_text)
    CALL message("init_landcover_change()","============================================")

    IF (PRESENT(stream)) THEN 
       IF (.NOT. ASSOCIATED(stream)) THEN
          ! Add new stream
          CALL new_stream(stream, 'landCoverChange', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(stream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       END IF
       LCC_Stream => stream
    ELSE
       ! Add new stream
       CALL new_stream(LCC_Stream, 'landCoverChange', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(LCC_Stream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    END IF

    IF (PRESENT(Nstream)) THEN
       IF (.NOT. ASSOCIATED(Nstream)) THEN
          ! Add new stream
          CALL new_stream(Nstream, 'lccNitro', filetype=fileformat, ztype=fileztype)
          ! Set default stream options
          CALL default_stream_setting(Nstream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
       END IF
       LCC_NStream => Nstream
    ELSE
       ! Add new stream
       CALL new_stream(LCC_NStream, 'lccNitro', filetype=fileformat, ztype=fileztype)
       ! Set default stream options
       CALL default_stream_setting(LCC_NStream, table=veg_table, repr=LAND, lpost=.TRUE., lrerun=.TRUE., leveltype=TILES)
    END IF

    !! --- add stream elements

    dim1p = (/ l_nland /)
    dim1  = (/ g_nland /)
    dim1n = (/ 'landpoint' /)

    dim3p = (/ l_nland, ntiles /)
    dim3  = (/ g_nland, ntiles /)
    dim3n(1) = 'landpoint'
    dim3n(2) = 'tiles'

    ! code numbers of this routine range from 80-84, 89-130, 224, 226-228, 230-232, 234-236 and 238-246 using the GRIB veg_table
    ! code numbers of this routine range from 85-88 and 90 and use the GRIB nitrogen_table
    CALL add(LCC_Stream, 'LCC_coverFract_target', landcover_change%LCC_coverFract_target,                     &
             longname="Vegetation cover fractions that should be reached at the end of the year as a result of landcover change.", &
             units='', lmiss=.TRUE., missval=missing_value,                                                   &
             ldims=dim3p, gdims=dim3, dimnames=dim3n, code=81, lrerun=.TRUE., contnorest=.TRUE., lpost=.true.)

    CALL add(LCC_Stream, 'LCC_flux_box_C2atmos', landcover_change%LCC_flux_box_C2atmos,                       &
             longname="Carbon flux to atmosphere from land cover/land use change",                            &
             units='mol(C) m-2(grid box) s-1',  lmiss=.TRUE., missval=missing_value,                          &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=82, lpost=.TRUE., lrerun=.FALSE.)
    CALL add(LCC_Stream, 'LCC_flux_box_C2litterGreenPools', landcover_change%LCC_flux_box_C2litterGreenPools, &
             longname="Flux of green and reserve carbon relocated by landcover change to the two green litter pools", &
             units='mol(C) m-2(grid box) s-1',  lmiss=.TRUE., missval=missing_value,                          &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=83, lpost=.TRUE., lrerun=.FALSE.)
    CALL add(LCC_Stream, 'LCC_flux_box_C2litterWoodPool', landcover_change%LCC_flux_box_C2litterWoodPool,     &
             longname="Wood carbon relocated by landcover change to woody litter pool",                       &
             units='mol(C) m-2(grid box) s-1',  lmiss=.TRUE., missval=missing_value,                          &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=84, lpost=.TRUE., lrerun=.FALSE.)

    IF (debug_Cconservation) THEN
       CALL add(LCC_stream, 'LCC_testCconserv', landcover_change%LCC_testCconserv,                            &
             longname="For testing carbon conservation (should be zero if conserved)",                        &
             units='mol(C) m-2(grid box)',                                                                    &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=89, lrerun=.TRUE., contnorest=.TRUE., lpost=.TRUE.)
    END IF

    IF (with_nitrogen) THEN
       CALL add(LCC_NStream, 'LCC_flux_box_N2atmos', landcover_change%LCC_flux_box_N2atmos,                   &
             longname="Nitrogen emitted to atmosphere due to land cover change",                              &
             units='mol(N) m-2(vegetated area) s-1',                                                          &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=85, lpost=.TRUE., lrerun=.FALSE.)
       CALL add(LCC_NStream, 'LCC_flux_box_N2litterGreenPools', landcover_change%LCC_flux_box_N2litterGreenPools, &
             longname="Nitrogen relocated from green and reserve pools to the green litter pools due to LCC", &
             units='mol(N) m-2(vegetated area) s-1',                                                          &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=86, lpost=.TRUE., lrerun=.FALSE.)
       CALL add(LCC_NStream, 'LCC_flux_box_N2litterWoodPool', landcover_change%LCC_flux_box_N2litterWoodPool, &
             longname="Wood nitrogen relocated by landcover change to woody litter pool",                     &
             units='mol(N) m-2(vegetated area) s-1',                                                              &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=87, lpost=.TRUE., lrerun=.FALSE.)
       CALL add(LCC_NStream, 'LCC_flux_box_N2SMINNpool', landcover_change%LCC_flux_box_N2SMINNpool,           &
             longname="Wood nitrogen relocated by landcover change to soil mineral N pool",                   &
             units='mol(N) m-2(vegetated area) s-1',                                                              &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=88, lpost=.TRUE., lrerun=.FALSE.)
       CALL add(LCC_NStream, 'LCC_testNconserv', landcover_change%LCC_testNconserv,                           &
             longname="For testing nitrogen conservation (should be zero if conserved)",                      &
             units='mol(C) m-2(grid box)',                                                                    &
             ldims=dim1p, gdims=dim1, dimnames=dim1n, code=90, lrerun=.TRUE., contnorest=.TRUE., lpost=.TRUE.)
    END IF

    IF (UseLanduseTransitions) THEN

       CALL add(LCC_stream, 'TransMtrx_CROP_2_PAST', landuse_transitions%TransMtrx_CROP_2_PAST,                  &
                longname="Landuse transition matrix element for conversion of croplands to pastures",            &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=91, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_PAST_2_CROP', landuse_transitions%TransMtrx_PAST_2_CROP,                  &
                longname="Landuse transition matrix element for conversion of pastures to croplands",            &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=92, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_NATL_2_PAST', landuse_transitions%TransMtrx_NATL_2_PAST,                  &
                longname="Landuse transition matrix element for conversion of natural vegetation to pastures",   &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=93, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_PAST_2_NATL', landuse_transitions%TransMtrx_PAST_2_NATL,                  &
                longname="Landuse transition matrix element for conversion of pastures to natural vegetation",   &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=94, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_CROP_2_NATL', landuse_transitions%TransMtrx_CROP_2_NATL,                  &
                longname="Landuse transition matrix element for conversion of croplands to natural vegetation",  &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=95, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_NATL_2_CROP', landuse_transitions%TransMtrx_NATL_2_CROP,                  &
                longname="Landuse transition matrix element for conversion of natural vegetation to croplands",  &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=96, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_FRST_2_PAST', landuse_transitions%TransMtrx_FRST_2_PAST,                  &
                longname="Landuse transition matrix element for conversion of natural non-grassland vegetation to pastures",&
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=97, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_PAST_2_FRST', landuse_transitions%TransMtrx_PAST_2_FRST,                  &
                longname="Landuse transition matrix element for conversion of pastures to natural non-grassland vegetation",&
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=98, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_GRAS_2_PAST', landuse_transitions%TransMtrx_GRAS_2_PAST,                  &
                longname="Landuse transition matrix element for conversion of grasslands to pastures",           &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=99, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_PAST_2_GRAS', landuse_transitions%TransMtrx_PAST_2_GRAS,                  &
                longname="Landuse transition matrix element for conversion of pastures to grasslands",           &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=100, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_FRST_2_CROP', landuse_transitions%TransMtrx_FRST_2_CROP,                  &
                longname="Landuse transition matrix element for conversion of natural non-grassland vegetation to croplands",&
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=101, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_CROP_2_FRST', landuse_transitions%TransMtrx_CROP_2_FRST,                  &
                longname="Landuse transition matrix element for conversion of croplands to natural non-grassland vegetation",&
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=102, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_GRAS_2_CROP', landuse_transitions%TransMtrx_GRAS_2_CROP,                  &
                longname="Landuse transition matrix element for conversion of grasslands to croplands",          &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=103, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'TransMtrx_CROP_2_GRAS', landuse_transitions%TransMtrx_CROP_2_GRAS,                  &
                longname="Landuse transition matrix element for conversion of croplands to grasslands",          &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=104, lrerun=.TRUE., lpost=debug)

       CALL add(LCC_stream, 'Grass_coverFract_lastYear', landuse_transitions%Grass_coverFract_lastYear,          &
                longname="Cover fraction for grasslands at last time step of last year",                         &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=105, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'NatWood_coverFract_lastYear', landuse_transitions%NatWood_coverFract_lastYear,      &
                longname="Cover fraction for natural non-grassland vegetation at last time step of last year",   &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=106, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'Pasture_coverFract_lastYear', landuse_transitions%Pasture_coverFract_lastYear,      &
                longname="Cover fraction for pastures at last time step of last year",                           &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=107, lrerun=.TRUE., lpost=debug)
       CALL add(LCC_stream, 'Crop_coverFract_lastYear', landuse_transitions%Crop_coverFract_lastYear,            &
                longname="Cover fraction for croplands at last time step of last year",                          &
                units='-',                                                                                       &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=108, lrerun=.TRUE., lpost=debug)

       IF (debug_Cconservation) THEN
          CALL add(LCC_stream, 'Test_NATL_2_PAST', landuse_transitions%Test_NATL_2_PAST,                 &
               longname="Test of conversion from natural vegetation to pastures (should be zero)",       &
               units='-',                                                                                &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=109, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_PAST_2_NATL', landuse_transitions%Test_PAST_2_NATL,                 &
               longname="Test of conversion from pastures to natural vegetation (should be zero)",       &
               units='-',                                                                                &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=110, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_CROP_2_NATL', landuse_transitions%Test_CROP_2_NATL,                 &
               longname="Test of conversion from croplands to natural vegetation (should be zero)",      &
               units='-',                                                                                &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=111, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_NATL_2_CROP', landuse_transitions%Test_NATL_2_CROP,                 &
               longname="Test of conversion from natural vegetation to croplands (should be zero)",      &
               units='-',                                                                                &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=112, lrerun=.FALSE., lpost=.TRUE.)


          CALL add(LCC_stream, 'Test_CROP_2_PAST', landuse_transitions%Test_CROP_2_PAST,                 &
               longname="Test of conversion from croplands to pasture (should be zero)",                 &
               units='-',                                                                                &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=113, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_PAST_2_CROP', landuse_transitions%Test_PAST_2_CROP,                 &
               longname="Test of conversion from  pasture to croplands (should be zero)",                &
               units='-',                                                                                &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=114, lrerun=.FALSE., lpost=.TRUE.)


          CALL add(LCC_stream, 'Test_FRST_2_PAST', landuse_transitions%Test_FRST_2_PAST,                         &
               longname="Test of conversion from natural non-grassland vegetation to pastures (should be zero)", &
               units='-',                                                                                        &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=115, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_PAST_2_FRST', landuse_transitions%Test_PAST_2_FRST,                         &
               longname="Test of conversion from pastures to natural non-grassland vegetation (should be zero)", &
               units='-',                                                                                        &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=116, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_GRAS_2_PAST', landuse_transitions%Test_GRAS_2_PAST,                         &
               longname="Test of conversion from grasslands to pastures (should be zero)",                       &
               units='-',                                                                                        &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=117, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_PAST_2_GRAS', landuse_transitions%Test_PAST_2_GRAS,                         &
               longname="Test of conversion from pastures to grasslands (should be zero)",                       &
               units='-',                                                                                        &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=118, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_GRAS_2_CROP', landuse_transitions%Test_GRAS_2_CROP,                         &
               longname="Test of conversion from grasslands to croplands (should be zero)",                      &
               units='-',                                                                                        &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=119, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_CROP_2_GRAS', landuse_transitions%Test_CROP_2_GRAS,                         &
               longname="Test of conversion from croplands to grasslands (should be zero)",                      &
               units='-',                                                                                        &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=120, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_FRST_2_CROP', landuse_transitions%Test_FRST_2_CROP,          &
               longname="Test of conversion from natural non-grassland vegetation to croplands (should be zero)",&
               units='-',                                                                                        &
               ldims=dim1p, gdims=dim1, dimnames=dim1n, code=121, lrerun=.FALSE., lpost=.TRUE.)
          CALL add(LCC_stream, 'Test_CROP_2_FRST', landuse_transitions%Test_CROP_2_FRST,                         &
               longname="Test of conversion from croplands to natural non-grassland vegetation (should be zero)",&
               units='-',                                                                                        &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=122, lrerun=.FALSE., lpost=.TRUE.)
       END IF ! debug_Cconservation

       CALL add(LCC_stream, 'CROP_2_NATL_ignored', landuse_transitions%CROP_2_NATL_ignored,                       &
                longname="Cover fraction that could not be transfered from croplands to natural vegetation",      &
                units='-', lmiss=.TRUE., missval=missing_value,                                                   &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=123, lrerun=.TRUE., lpost=.TRUE.)
       CALL add(LCC_stream, 'PAST_2_NATL_ignored', landuse_transitions%PAST_2_NATL_ignored,                       &
                longname="Cover fraction that could not be transfered from pastures to natural vegetation",       &
                units='-', lmiss=.TRUE., missval=missing_value,                                                   &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=124, lrerun=.TRUE., lpost=.TRUE.)

       CALL add(LCC_stream, 'Box_harvest', landuse_transitions%Box_harvest,                                       &
                longname="Prescribed harvest from natural vegetation for full gridbox",                           &
                units='mol(C) m-2(gridbox) s-1', lmiss=.TRUE., missval=missing_value,                             &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=125, lrerun=.TRUE., lpost=.TRUE.)
       CALL add(LCC_stream, 'Box_flux_harvest', landuse_transitions%Box_flux_harvest,                             &
                longname="Total harvest flux from natural vegetation for full gridbox",                           &
                units='mol(C) m-2(gridbox) s-1', lmiss=.TRUE., missval=missing_value,                             &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=126, lrerun=.TRUE., lpost=.TRUE., laccu=.true.)
       CALL add(LCC_stream, 'Box_flux_harvest_2atmos', landuse_transitions%Box_flux_harvest_2atmos,               &
                longname="Harvest flux from natural vegetation emitted from gridbox as CO2 to atmosphere",        &
                units='mol(C) m-2(gridbox) s-1', lmiss=.TRUE., missval=missing_value,                             &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=80, lrerun=.TRUE., lpost=.TRUE., laccu=.true.)

       CALL add(LCC_stream, 'C2atmos_LUtrans', landuse_transitions%C2atmos_LUtrans,                               &
                longname="Carbon flux to atmosphere from land use change per time step",                          &
                units='mol(C) m-2(grid box)', CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,             &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=127, lpost=.FALSE., laccu=.FALSE., lrerun=.TRUE.)
       CALL add(LCC_stream, 'C2atmos_harvest', landuse_transitions%C2atmos_harvest,                               &
                longname="Carbon flux to atmosphere from harvest per time step",                                  &
                units='mol(C) m-2(grid box)', CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,             &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=128, lpost=.FALSE., laccu=.FALSE., lrerun=.TRUE.)
       IF (with_nitrogen) THEN
         CALL add(LCC_Nstream, 'N2atmos_LUtrans', landuse_transitions%N2atmos_LUtrans,                            &
                longname="Nitrogen flux to atmosphere from land use change per time step",                        &
                units='mol(N) m-2(grid box)', lmiss=.TRUE., missval=missing_value,                                &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=129, lpost=.FALSE., laccu=.FALSE., lrerun=.TRUE.)
         CALL add(LCC_Nstream, 'N2atmos_harvest', landuse_transitions%N2atmos_harvest,                            &
                longname="Nitrogen flux to atmosphere from harvest per time step",                                &
                units='mol(M) m-2(grid box)', lmiss=.TRUE., missval=missing_value,                                &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=130, lpost=.FALSE., laccu=.FALSE., lrerun=.TRUE.)	
       END IF
    ELSE

       CALL add(LCC_stream, 'C2atmos_LCC', landcover_change%C2atmos,                                              &
                longname="Carbon flux to atmosphere from land cover change per time step",                        &
                units='mol(C) m-2(grid box)', CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,             &
                ldims=dim1p, gdims=dim1, dimnames=dim1n, code=127, lpost=.FALSE., laccu=.FALSE., lrerun=.TRUE.)

    ENDIF ! UseLanduseTransitions

    IF (lcc_scheme==2) THEN
       CALL add(LCC_stream,'C_flux_2_onSite_LCC', landcover_change%C_2_onSite,units='mol(C) m-2(veg) s-1', &
                longname='Carbon flux from living plants to ground pool from land use change',code=224,          &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,          &
                lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.)
       CALL add(LCC_stream,'C_flux_2_paper_LCC', landcover_change%C_2_paper,units='mol(C) m-2(veg) s-1',                           &
                longname='Wood carbon flux from living plants to short/intermediate term anthropogenic pool from land use change', &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                      &
                lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.,code=226)
       CALL add(LCC_stream,'C_flux_2_construction_LCC', landcover_change%C_2_construction,units='mol(C) m-2(veg) s-1',             &
                longname='Wood carbon flux from living plants to long term anthropogenic pool from land use change',code=227,      &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                      &
                lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.)
       ! Fluxes from living pools to anthropogenic lcc cpools on grid box area averaged over output period 
       CALL add(LCC_stream,'boxC_flux_2_onSite_LCC',landcover_change%boxC_2_onSite,units='mol(C) m-2(grid box) s-1', &
                longname='Carbon flux from living plants to ground pool from land use change',code=228,                      &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                      &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
       CALL add(LCC_stream,'boxC_flux_2_paper_LCC', landcover_change%boxC_2_paper,units='mol(C) m-2(grid box) s-1',              &
                longname='Wood carbon flux from living plants to short/intermediate term anthropogenic pool from land use change', &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                      &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.,code=230)
       CALL add(LCC_stream,'boxC_flux_2_construction_LCC', landcover_change%boxC_2_construction,units='mol(C) m-2(grid box) s-1',&
                longname='Wood carbon flux from living plants to long term anthropogenic pool from land use change',code=231,    &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                    &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
       ! Fluxes from anthropogenic lcc cpools to atmosphere on vegetated area
       CALL add(LCC_stream,'C_flux_onSite_2_atmos_LCC', landcover_change%C_onSite_2_atmos,units='mol(C) m-2(veg) s-1', &
                longname='Carbon flux ground pool from land use change to atmosphere',code=232,                              &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                      &
                lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.)
       CALL add(LCC_stream,'C_flux_paper_2_atmos_LCC', landcover_change%C_paper_2_atmos,units='mol(C) m-2(veg) s-1',code=234,      &
                longname='Carbon flux from short/intermediate term anthropogenic pool from land use change to atmos',              &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                      &
                lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.)
       CALL add(LCC_stream,'C_flux_construction_2_atmos_LCC',landcover_change%C_construction_2_atmos,units='mol(C) m-2(veg) s-1',  &
                longname='Carbon flux from long term anthropogenic pool from land use change to atmos',code=235,                   &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                      &
                lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.)
       ! Fluxes from anthropogenic lcc cpools to atmosphere on grid box area averaged over output period 
       CALL add(LCC_stream,'boxC_flux_onSite_2_atmos_LCC',landcover_change%boxC_onSite_2_atmos,   &
                units='mol(C) m-2(grid box) s-1',                                                             &
                longname='Carbon flux from ground pool from land use change to atmos',code=236,         &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n, &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
       CALL add(LCC_stream,'boxC_flux_paper_2_atmos_LCC', landcover_change%boxC_paper_2_atmos,units='mol(C) m-2(grid box) s-1',  &
                longname='Carbon flux from short/intermediate term anthropogenic pool from land use change to atmos',            &
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                    &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.,code=238)
       CALL add(LCC_stream,'boxC_flux_construction_2_atmos_LCC', landcover_change%boxC_construction_2_atmos,    &
                units='mol(C) m-2(grid box) s-1',                                                               &
                longname='Carbon flux from long term anthropogenic pool from land use change to atmos',code=239,&
                CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,   &
                lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)

       IF (UseLanduseTransitions) THEN
          ! Fluxes from living pools to anthropogenic harvest cpools on vegetated area
          CALL add(LCC_stream,'C_flux_2_paper_harvest', landcover_change%C_2_paper_harv,units='mol(C) m-2(veg) s-1',longname =  &
                   'Wood carbon flux from living plants to short/intermediate term anthropogenic pool from harvest',            &
                   CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                &
                   lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.,code=240)
          CALL add(LCC_stream,'C_flux_2_construction_harvest',landcover_change%C_2_construction_harv,units='mol(C) m-2(veg) s-1'&
                   ,longname='Wood carbon flux from living plants to long term anthropogenic pool from harvest',code=241,       &
                   CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                &
                   lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.)
          ! Fluxes from living pools to anthropogenic harvest cpools on grid box area averaged over output period 
          CALL add(LCC_stream,'boxC_flux_2_paper_harvest', landcover_change%boxC_2_paper_harv,units='mol(C) m-2(grid box) s-1', &
                   longname='Wood carbon flux from living plants to short/intermediate term anthropogenic pool from harvest',   &
                   CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                &
                   lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.,code=242)
          CALL add(LCC_stream,'boxC_flux_2_construction_harvest', landcover_change%boxC_2_construction_harv,     &
                   units='mol(C) m-2(grid box) s-1',code=243,                                                    &
                   longname='Wood carbon flux from living plants to long term anthropogenic pool from harvest',  &
                   CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n, &
                   lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
          ! Fluxes from anthropogenic harvest cpools on vegetated area to atmosphere
          CALL add(LCC_stream,'C_flux_paper_2_atmos_harvest', landcover_change%C_paper_harv_2_atmos,units='mol(C) m-2(veg) s-1',&
                   longname = 'Carbon flux from short/intermediate term anthropogenic pool from harvest to atmos',              &
                   CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                &
                   lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.,code=244)
          CALL add(LCC_stream,'C_flux_construction_harvest_2_atmos',landcover_change%C_construction_harv_2_atmos, code=245,     &
                   units='mol(C) m-2(veg) s-1', longname='Carbon flux from long term anthropogenic pool from harvest to atmos', &
                   CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                &
                   lpost=.TRUE.,laccu=.FALSE.,lrerun=.FALSE.)
          ! Fluxes from anthropogenic harvest cpools on grid box area averaged over output period 
          CALL add(LCC_stream,'boxC_flux_paper_2_atmos_harvest', landcover_change%boxC_paper_harv_2_atmos,code=246,             &
                   units='mol(C) m-2(grid box) s-1',                                                                            &
                   longname='Carbon flux short/intermediate term anthropogenic pool from harvest to atmos',                     &
                   CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,                &
                   lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
          CALL add(LCC_stream,'boxC_flux_construction_2_atmos_harvest', landcover_change%boxC_construction_harv_2_atmos,   &
                   units='mol(C) m-2(grid box) s-1',code=247,                                                              &
                   longname='Carbon flux from long term anthropogenic pool from harvest to atmos',                         &
                   CONTNOREST=.TRUE., lmiss=.TRUE., missval=missing_value,ldims=dim1p,gdims=dim1,dimnames=dim1n,           &
                   lpost=.TRUE.,laccu=.TRUE.,lrerun=.FALSE.)
       ENDIF 

    ENDIF

    landcover_change%LCC_flux_box_C2atmos           = 0.0_dp
    landcover_change%LCC_flux_box_C2litterGreenPools= 0.0_dp
    landcover_change%LCC_flux_box_C2litterWoodPool  = 0.0_dp

    IF (debug_Cconservation) THEN
       landcover_change%LCC_testCconserv = 0.0_dp
    END IF

    IF (with_nitrogen) THEN
       landcover_change%LCC_flux_box_N2atmos            = 0.0_dp
       landcover_change%LCC_flux_box_N2litterGreenPools = 0.0_dp
       landcover_change%LCC_flux_box_N2litterWoodPool   = 0.0_dp
       landcover_change%LCC_flux_box_N2SMINNpool        = 0.0_dp
       landcover_change%LCC_testNconserv                = 0.0_dp
    END IF

    IF(UseLanduseTransitions) THEN
       landuse_transitions%Box_flux_harvest        = 0.0_dp
       landuse_transitions%Box_flux_harvest_2atmos = 0.0_dp
       landuse_transitions%C2atmos_LUtrans         = 0.0_dp
       landuse_transitions%C2atmos_harvest         = 0.0_dp
       IF (with_nitrogen) THEN
          landuse_transitions%N2atmos_LUtrans      = 0.0_dp
          landuse_transitions%N2atmos_harvest      = 0.0_dp
       END IF
    ELSE
       landcover_change%C2atmos = 0.0_dp
    END IF

    IF (lcc_scheme==2) THEN
       landcover_change%boxC_2_onSite            (:) = 0._dp
       landcover_change%boxC_2_paper             (:) = 0._dp
       landcover_change%boxC_2_construction      (:) = 0._dp
       landcover_change%boxC_onSite_2_atmos      (:) = 0._dp
       landcover_change%boxC_paper_2_atmos       (:) = 0._dp
       landcover_change%boxC_construction_2_atmos(:) = 0._dp
       IF (UseLanduseTransitions) THEN
          landcover_change%boxC_2_paper_harv             (:) = 0._dp
          landcover_change%boxC_2_construction_harv      (:) = 0._dp
          landcover_change%boxC_paper_harv_2_atmos       (:) = 0._dp
          landcover_change%boxC_construction_harv_2_atmos(:) = 0._dp
       ENDIF
    ENDIF

    IF(isRestart ) THEN
       landcover_change%LCC_coverFract_target = 0.0_dp
    END IF

    lcc = landcover_change ! Make landcover_change globally available

    !! ======= CONSISTENCY CHECKS ==========================
    !!
    !! Here the consistency of entries in the Lct-lib with assumptions underlying the implementation of the New Hampshire Harmonized
    !! Protocol (see subroutine do_landcover_transitions()) are checked.

    !! Exactly two grass types: C3+C4
    !! Exactly two pasture types: C3+C4

    !! CHR: CHECKS NEED TO BE IMPLEMENTED!!

  END SUBROUTINE init_landcover_change

  ! Accumulate anthropogenic pools and fluxes to/from those 
  ! while converting from m-2(veg) to m-2(grid box) and fluxes from d-1 to s-1
  ! Time averaged vars. are declared with laccu=.true. in the add_stream_element calls and therefore needs an extra factor 86400
  SUBROUTINE CumulateAnthroFluxes(lcc, veg_ratio_max)
    TYPE(landcover_change_type), INTENT(inout) :: lcc                !! all pools to cumulate
    REAL(dp),                    INTENT(in)    :: veg_ratio_max(:)   !! maximum fraction of grid cell that can be covered by

      ! Fluxes into pools
      lcc%boxC_2_onSite      (:) = lcc%boxC_2_onSite      (:) + lcc%C_2_onSite      (:) * veg_ratio_max(:)
      lcc%boxC_2_paper       (:) = lcc%boxC_2_paper       (:) + lcc%C_2_paper       (:) * veg_ratio_max(:)
      lcc%boxC_2_construction(:) = lcc%boxC_2_construction(:) + lcc%C_2_construction(:) * veg_ratio_max(:)

      ! Fluxes from pools to atmosphere
      lcc%boxC_onSite_2_atmos      (:)=lcc%boxC_onSite_2_atmos      (:)+lcc%C_onSite_2_atmos      (:)*veg_ratio_max(:)
      lcc%boxC_paper_2_atmos       (:)=lcc%boxC_paper_2_atmos       (:)+lcc%C_paper_2_atmos       (:)*veg_ratio_max(:)
      lcc%boxC_construction_2_atmos(:)=lcc%boxC_construction_2_atmos(:)+lcc%C_construction_2_atmos(:)*veg_ratio_max(:)

      IF (ASSOCIATED(lcc%boxC_2_paper_harv)) THEN

        ! Harvest fluxes into pools
        lcc%boxC_2_paper_harv       (:) = lcc%boxC_2_paper_harv       (:) + lcc%C_2_paper_harv       (:) * veg_ratio_max(:)
        lcc%boxC_2_construction_harv(:) = lcc%boxC_2_construction_harv(:) + lcc%C_2_construction_harv(:) * veg_ratio_max(:)

        ! Fluxes from harvest pools to atmosphere
        lcc%boxC_paper_harv_2_atmos       (:) = &
        lcc%boxC_paper_harv_2_atmos       (:) + lcc%C_paper_harv_2_atmos       (:) * veg_ratio_max(:)
        lcc%boxC_construction_harv_2_atmos(:) = &
        lcc%boxC_construction_harv_2_atmos(:) + lcc%C_construction_harv_2_atmos(:) * veg_ratio_max(:)
      ENDIF

  END SUBROUTINE CumulateAnthroFluxes

  !! --- do_landcover_change() ----------------------------------------------------
  !!
  !! This routine:
  !!    1. Replaces the old map of land cover fractions by a new map that is read from file
  !!    2. Performs the necessary changes in the carbon pools that go along
  !!       with landcover changes
  !!
  !! NOTE: This routine assumes that the restart files have already been
  !! read, because they contain the old landcover distribution that is needed
  !! to derive the changes in the carbon pools
  !!

  SUBROUTINE do_landcover_change(ntiles, lctlib, surface, cbalance, nbalance, with_nitrogen,   &
                                 with_yasso, LeafLit_coef, WoodLit_coef,                       &
                                 veg_ratio_max, veg_fract_correction, landcover_fract_current, &
                                 CO2_emission, nitrogen_2_atmos, lcc, lcc_scheme)
    use mo_jsbach_lctlib,    only: lctlib_type
    use mo_jsbach_grid,      only: kstart, kend, nidx
    use mo_cbal_bethy,       only: cbalance_type, nbalance_type
    use mo_cbal_cpools,      only: relocate_CarbonAndNitrogen
    USE mo_time_control,     ONLY: current_date, get_year_day, lstart
    use mo_jsbach_constants, only: molarMassCO2_kg
    use mo_time_conversion,  only: year_len
    use mo_time_control,     only: l_trigfiles, lstart
    USE mo_land_surface,     ONLY: land_surface_type

    integer,                  intent(in)    :: ntiles                       !! number of tiles of the land points
    type(lctlib_type),        intent(in)    :: lctlib                       !! PTF-specific parameters
    type(land_surface_type),  intent(in)    :: surface                      !! surface parameters
    type(cbalance_type),      intent(inout) :: cbalance                     !! the carbon pools that will be changed in this call
    type(nbalance_type),      intent(inout) :: nbalance                     !! the nitrogen pools that will be changed in this call
    LOGICAL,                  INTENT(in)    :: with_nitrogen                !! If .true. Nitrogen will also be cycled
    LOGICAL,                  INTENT(in)    :: with_yasso                   !! If .true. Yasso mode will be used
    real(dp),                 intent(in)    :: LeafLit_coef(:,:,:)          !! Factor to seperate non woody litter fall to Yasso
                                                                            !!     litter pool
    real(dp),                 intent(in)    :: WoodLit_coef(:,:,:)          !! Factor to seperate woody litter fall to Yasso litter
                                                                            !!     pools
    real(dp),                 intent(in)    :: veg_ratio_max(:)             !! maximum fraction of grid cell that can be covered by
                                                                            !!    vegetation
    real(dp),                 intent(in)    :: veg_fract_correction(:,:)    !! Corr. factor for cover fractions 1-exp(-LAI_max/2)
    real(dp),                 intent(inout) :: landcover_fract_current(:,:) !! On call old cover fractions, on return the new ones
    real(dp),                 intent(out)   :: CO2_emission(:)              !! The CO2 emission from landcover change in 
                                                                            !! kg(CO2)/(m^2 s) in this timestep
    type(landcover_change_type), optional, intent(inout) :: lcc             !! Access to anthropogenic pools

    real(dp),optional,intent(out)         :: Nitrogen_2_atmos(:)            !! Amount of nitrogen directly emitted to atmosphere 
    integer, optional,intent(in)          :: lcc_scheme                     !! Which lcc scheme is used 

    !! locals

    integer             :: dayNo_in_current_year
    real(dp)            :: C2litterWoodPool_ag(1:nidx)
    real(dp)            :: C2litterWoodPool_bg(1:nidx)
    real(dp)            :: C2litterGreenPools(1:nidx)
    real(dp)            :: N2atmos(1:nidx)
    real(dp)            :: N2litterWoodPool_ag(1:nidx)
    real(dp)            :: N2litterWoodPool_bg(1:nidx)
    real(dp)            :: N2litterGreenPools(1:nidx)
    real(dp)            :: N2SMINNpool(1:nidx)
    real(dp)            :: landcover_fract_new(1:nidx,1:ntiles)
    integer             :: noOfDays_in_current_year

    !! -- Restart summations after a restart

    if(l_trigfiles) then  !! Restart summations with begin of a new output interval
       IF (with_nitrogen) THEN
          landcover_change%LCC_flux_box_N2atmos(kstart:kend) = 0.0_dp
          lcc%LCC_flux_box_N2litterGreenPools(kstart:kend) = 0.0_dp
          lcc%LCC_flux_box_N2litterWoodPool(kstart:kend) = 0.0_dp
          lcc%LCC_flux_box_N2SMINNpool(kstart:kend) = 0.0_dp
       END IF
    endif


    !! --- Check for exiting the routine without landcover change
    IF (.NOT. new_day .AND. .NOT. lstart) THEN !! We are not at the first time step of a day ..
       !! ... so we use the same fluxes during the whole day ...
       CO2_emission(:) = lcc%C2atmos(kstart:kend) * veg_ratio_max(:) * molarMassCO2_kg/86400._dp
       IF (with_nitrogen) Nitrogen_2_atmos(:) = 0.0_dp
       RETURN                                           !! .. and exit the routine
    END IF

    !! --- Initialize Carbon and Nitrogen conservation tests

    IF (debug_Cconservation) THEN
       IF (lcc_scheme==2) THEN
          lcc%LCC_testCconserv(kstart:kend) = lcc%LCC_testCconserv(kstart:kend) &
               + (  cbalance%Cpool_onSite      (kstart:kend)  &
                  + cbalance%Cpool_paper       (kstart:kend)  &
                  + cbalance%Cpool_construction(kstart:kend)  &
                 ) * veg_ratio_max(1:nidx)
       ENDIF
       IF (.NOT. with_yasso) THEN
          lcc%LCC_testCconserv(kstart:kend) = veg_ratio_max(1:nidx)            &
               * SUM(  landcover_fract_current(1:nidx,1:ntiles)                &
                     * veg_fract_correction   (1:nidx,1:ntiles)                &
                     * (  cbalance%Cpool_green          (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_woods          (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_reserve        (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles) &
                        + cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles) &
                        + cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_slow           (kstart:kend,1:ntiles) &
                       ), DIM=2 )
       ELSE !! with yasso
          lcc%LCC_testCconserv(kstart:kend) = veg_ratio_max(1:nidx)            &
               * SUM(  landcover_fract_current          (     1:nidx,1:ntiles) &
                     * veg_fract_correction             (     1:nidx,1:ntiles) &
                     * (  cbalance%Cpool_green          (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_woods          (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_reserve        (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_acid_ag1      (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_acid_bg1      (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_water_ag1     (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_water_bg1     (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_ethanol_ag1   (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_ethanol_bg1   (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_nonsoluble_ag1(kstart:kend,1:ntiles) &
                        + cbalance%YCpool_nonsoluble_bg1(kstart:kend,1:ntiles) &
                        + cbalance%YCpool_humus_1       (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_acid_ag2      (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_acid_bg2      (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_water_ag2     (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_water_bg2     (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_ethanol_ag2   (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_ethanol_bg2   (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_nonsoluble_ag2(kstart:kend,1:ntiles) &
                        + cbalance%YCpool_nonsoluble_bg2(kstart:kend,1:ntiles) &
                        + cbalance%YCpool_humus_2       (kstart:kend,1:ntiles) &
                       ), DIM=2 )
       END IF

       IF(with_nitrogen) THEN
          lcc%LCC_testNconserv(kstart:kend) = veg_ratio_max(1:nidx)        &
           * SUM(  landcover_fract_current(1:nidx,1:ntiles)                &
                 * veg_fract_correction   (1:nidx,1:ntiles)                &
                 * (  nbalance%Npool_green          (kstart:kend,1:ntiles) &
                    + nbalance%Npool_woods          (kstart:kend,1:ntiles) &
                    + nbalance%Npool_mobile         (kstart:kend,1:ntiles) &
                    + nbalance%Npool_litter_green_ag(kstart:kend,1:ntiles) &
                    + nbalance%Npool_litter_green_bg(kstart:kend,1:ntiles) &
                    + nbalance%Npool_litter_wood_ag (kstart:kend,1:ntiles) &
                    + nbalance%Npool_litter_wood_bg (kstart:kend,1:ntiles) &
                    + nbalance%Npool_slow           (kstart:kend,1:ntiles) &
                    + nbalance%SMINN_pool           (kstart:kend,1:ntiles) &
                   ) , DIM=2 )
       ENDIF
    END IF

    !! --- determine new cover fractions (we are at the first time step of a day)

    noOfDays_in_current_year = year_len(current_date)
    dayNo_in_current_year = floor(get_year_day(current_date))
    landcover_fract_new(1:nidx,:) =                                                                 &
                 landcover_fract_current(:,:)                                                       &
                 + ( lcc%LCC_coverFract_target(kstart:kend,:) - landcover_fract_current(1:nidx,:) ) &
                 / real(1 + noOfDays_in_current_year - dayNo_in_current_year,dp)

    !! --- Do the changes in the carbon pools and compute amount of carbon to be released to atmosphere

    IF (with_nitrogen) THEN
       CALL relocate_CarbonAndNitrogen (lctlib, surface, landcover_fract_current(1:nidx,1:ntiles), &
                         landcover_fract_new           (     1:nidx,1:ntiles),                 &
                         veg_fract_correction          (     1:nidx,1:ntiles),                 &
                         cbalance%Cpool_green          (kstart:kend,1:ntiles),                 &
                         cbalance%Cpool_woods          (kstart:kend,1:ntiles),                 &
                         cbalance%Cpool_reserve        (kstart:kend,1:ntiles),                 &
                         Cpool_litter_green_ag        = cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles), &
                         Cpool_litter_green_bg        = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles), &
                         Cpool_litter_wood_ag         = cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles), &
                         Cpool_litter_wood_bg         = cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles), &
                         Cpool_slow                   = cbalance%Cpool_slow           (kstart:kend,1:ntiles), &
                         C_2_atmos                    = lcc%C2atmos        (kstart:kend),      &
                         C_2_litterGreenPools         = C2litterGreenPools (     1:nidx),      &
                         C_2_litterWoodPool_ag        = C2litterWoodPool_ag(     1:nidx),      &
                         C_2_litterWoodPool_bg        = C2litterWoodPool_bg(     1:nidx),      &
                         Npool_green                  = nbalance%Npool_green          (kstart:kend,1:ntiles), &
                         Npool_woods                  = nbalance%Npool_woods          (kstart:kend,1:ntiles), &
                         Npool_mobile                 = nbalance%Npool_mobile         (kstart:kend,1:ntiles), &
                         Npool_litter_green_ag        = nbalance%Npool_litter_green_ag(kstart:kend,1:ntiles), &
                         Npool_litter_green_bg        = nbalance%Npool_litter_green_bg(kstart:kend,1:ntiles), &
                         Npool_litter_wood_ag         = nbalance%Npool_litter_wood_ag (kstart:kend,1:ntiles), &
                         Npool_litter_wood_bg         = nbalance%Npool_litter_wood_bg (kstart:kend,1:ntiles), &
                         Npool_slow                   = nbalance%Npool_slow           (kstart:kend,1:ntiles), &
                         SMINN_pool                   = nbalance%SMINN_pool           (kstart:kend,1:ntiles), &
                         Nitrogen_2_atmos             = N2atmos            (     1:nidx),      &
                         Nitrogen_2_litterGreenPools  = N2litterGreenPools (     1:nidx),      &
                         Nitrogen_2_litterWoodPool_ag = N2litterWoodPool_ag(     1:nidx),      &
                         Nitrogen_2_litterWoodPool_bg = N2litterWoodPool_bg(     1:nidx),      &
                         Nitrogen_2_SMINNpool         = N2SMINNpool        (     1:nidx),      &
                         lcc_scheme                   = lcc_scheme )
                        

    ELSE !! Without nitrogen
       SELECT CASE (lcc_scheme)
       CASE(1) ! Standard JSBACH scheme
          IF (.NOT. with_yasso) THEN
             CALL relocate_CarbonAndNitrogen(lctlib,surface,              &
                  landcover_fract_current       (     1:nidx,1:ntiles),     &
                  landcover_fract_new           (     1:nidx,1:ntiles),     &
                  veg_fract_correction          (     1:nidx,1:ntiles),     &
                  cbalance%Cpool_green          (kstart:kend,1:ntiles),     &
                  cbalance%Cpool_woods          (kstart:kend,1:ntiles),     &
                  cbalance%Cpool_reserve        (kstart:kend,1:ntiles),     &
                  Cpool_litter_green_ag        = cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles), &
                  Cpool_litter_green_bg        = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles), &
                  Cpool_litter_wood_ag         = cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles), &
                  Cpool_litter_wood_bg         = cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles), &
                  Cpool_slow                   = cbalance%Cpool_slow           (kstart:kend,1:ntiles), &
                  C_2_atmos                  = lcc%C2atmos(kstart:kend),    &
                  C_2_litterGreenPools       = C2litterGreenPools (1:nidx), &
                  C_2_litterWoodPool_ag      = C2litterWoodPool_ag(1:nidx), &
                  C_2_litterWoodPool_bg      = C2litterWoodPool_bg(1:nidx), &
                  lcc_scheme                 = lcc_scheme)
          ELSE   !! with yasso
             CALL relocate_CarbonAndNitrogen(lctlib,surface,              &
                  landcover_fract_current       (     1:nidx,1:ntiles),     &
                  landcover_fract_new           (     1:nidx,1:ntiles),     &
                  veg_fract_correction          (     1:nidx,1:ntiles),     &
                  cbalance%Cpool_green          (kstart:kend,1:ntiles),     &
                  cbalance%Cpool_woods          (kstart:kend,1:ntiles),     &
                  cbalance%Cpool_reserve        (kstart:kend,1:ntiles),     &
                  YCpool_acid_ag1            = cbalance%YCpool_acid_ag1      (kstart:kend,1:ntiles), &
                  YCpool_water_ag1           = cbalance%YCpool_water_ag1     (kstart:kend,1:ntiles), &
                  YCpool_ethanol_ag1         = cbalance%YCpool_ethanol_ag1   (kstart:kend,1:ntiles), &
                  YCpool_nonsoluble_ag1      = cbalance%YCpool_nonsoluble_ag1(kstart:kend,1:ntiles), &
                  YCpool_acid_bg1            = cbalance%YCpool_acid_bg1      (kstart:kend,1:ntiles), &
                  YCpool_water_bg1           = cbalance%YCpool_water_bg1     (kstart:kend,1:ntiles), &
                  YCpool_ethanol_bg1         = cbalance%YCpool_ethanol_bg1   (kstart:kend,1:ntiles), &
                  YCpool_nonsoluble_bg1      = cbalance%YCpool_nonsoluble_bg1(kstart:kend,1:ntiles), &
                  YCpool_humus_1             = cbalance%YCpool_humus_1       (kstart:kend,1:ntiles), &
                  YCpool_acid_ag2            = cbalance%YCpool_acid_ag2      (kstart:kend,1:ntiles), &
                  YCpool_water_ag2           = cbalance%YCpool_water_ag2     (kstart:kend,1:ntiles), &
                  YCpool_ethanol_ag2         = cbalance%YCpool_ethanol_ag2   (kstart:kend,1:ntiles), &
                  YCpool_nonsoluble_ag2      = cbalance%YCpool_nonsoluble_ag2(kstart:kend,1:ntiles), &
                  YCpool_acid_bg2            = cbalance%YCpool_acid_bg2      (kstart:kend,1:ntiles), &
                  YCpool_water_bg2           = cbalance%YCpool_water_bg2     (kstart:kend,1:ntiles), &
                  YCpool_ethanol_bg2         = cbalance%YCpool_ethanol_bg2   (kstart:kend,1:ntiles), &
                  YCpool_nonsoluble_bg2      = cbalance%YCpool_nonsoluble_bg2(kstart:kend,1:ntiles), &
                  YCpool_humus_2             = cbalance%YCpool_humus_2       (kstart:kend,1:ntiles), &
                  C_2_atmos                  = lcc%C2atmos(kstart:kend),    &
                  C_2_litterGreenPools       = C2litterGreenPools (1:nidx), &
                  C_2_litterWoodPool_ag      = C2litterWoodPool_ag(1:nidx), &
                  C_2_litterWoodPool_bg      = C2litterWoodPool_bg(1:nidx), &
                  lcc_scheme                 = lcc_scheme,                  &
                  LeafLit_coef               = LeafLit_coef             (kstart:kend,1:ntiles,1:5), &
                  WoodLit_coef               = WoodLit_coef             (kstart:kend,1:ntiles,1:5))
                
          END IF
       CASE(2) ! LCC scheme Grand Slam after Houghton et al. (lcc_scheme==2) 
          IF (.NOT. with_yasso) THEN
             CALL relocate_CarbonAndNitrogen(lctlib,surface,            &
                  landcover_fract_current       (     1:nidx,1:ntiles),   &
                  landcover_fract_new           (     1:nidx,1:ntiles),   &
                  veg_fract_correction          (     1:nidx,1:ntiles),   &
                  cbalance%Cpool_green          (kstart:kend,1:ntiles),   &
                  cbalance%Cpool_woods          (kstart:kend,1:ntiles),   &
                  cbalance%Cpool_reserve        (kstart:kend,1:ntiles),   &
                  Cpool_litter_green_ag        = cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles), &
                  Cpool_litter_green_bg        = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles), &
                  Cpool_litter_wood_ag         = cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles), &
                  Cpool_litter_wood_bg         = cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles), &
                  Cpool_slow                   = cbalance%Cpool_slow           (kstart:kend,1:ntiles), &
                  C_2_atmos                     = lcc%C2atmos                (kstart:kend), &
                  C_2_litterGreenPools          = C2litterGreenPools         (     1:nidx), &
                  C_2_litterWoodPool_ag         = C2litterWoodPool_ag        (     1:nidx), &
                  C_2_litterWoodPool_bg         = C2litterWoodPool_bg        (     1:nidx), &                        
                  Cpool_onSite                  = cbalance%Cpool_onSite      (kstart:kend), &
                  Cpool_paper                   = cbalance%Cpool_paper       (kstart:kend), &
                  Cpool_construction            = cbalance%Cpool_construction(kstart:kend), &
                  C_onSite_2_atmos              = lcc%C_onSite_2_atmos       (kstart:kend), &
                  C_paper_2_atmos               = lcc%C_paper_2_atmos        (kstart:kend), &
                  C_construction_2_atmos        = lcc%C_construction_2_atmos (kstart:kend), &                         
                  C_2_onSite                    = lcc%C_2_onSite             (kstart:kend), &
                  C_2_paper                     = lcc%C_2_paper              (kstart:kend), &
                  C_2_construction              = lcc%C_2_construction       (kstart:kend), &
                  lcc_scheme                    = lcc_scheme )
          ELSE !! with yasso
             CALL relocate_CarbonAndNitrogen(lctlib,surface,            &
                  landcover_fract_current       (     1:nidx,1:ntiles),   &
                  landcover_fract_new           (     1:nidx,1:ntiles),   &
                  veg_fract_correction          (     1:nidx,1:ntiles),   &
                  cbalance%Cpool_green          (kstart:kend,1:ntiles),   &
                  cbalance%Cpool_woods          (kstart:kend,1:ntiles),   &
                  cbalance%Cpool_reserve        (kstart:kend,1:ntiles),   &
                  YCpool_acid_ag1            = cbalance%YCpool_acid_ag1      (kstart:kend,1:ntiles), &
                  YCpool_water_ag1           = cbalance%YCpool_water_ag1     (kstart:kend,1:ntiles), &
                  YCpool_ethanol_ag1         = cbalance%YCpool_ethanol_ag1   (kstart:kend,1:ntiles), &
                  YCpool_nonsoluble_ag1      = cbalance%YCpool_nonsoluble_ag1(kstart:kend,1:ntiles), &
                  YCpool_acid_bg1            = cbalance%YCpool_acid_bg1      (kstart:kend,1:ntiles), &
                  YCpool_water_bg1           = cbalance%YCpool_water_bg1     (kstart:kend,1:ntiles), &
                  YCpool_ethanol_bg1         = cbalance%YCpool_ethanol_bg1   (kstart:kend,1:ntiles), &
                  YCpool_nonsoluble_bg1      = cbalance%YCpool_nonsoluble_bg1(kstart:kend,1:ntiles), &
                  YCpool_humus_1             = cbalance%YCpool_humus_1      (kstart:kend,1:ntiles), &
                  YCpool_acid_ag2            = cbalance%YCpool_acid_ag2      (kstart:kend,1:ntiles), &
                  YCpool_water_ag2           = cbalance%YCpool_water_ag2     (kstart:kend,1:ntiles), &
                  YCpool_ethanol_ag2         = cbalance%YCpool_ethanol_ag2   (kstart:kend,1:ntiles), &
                  YCpool_nonsoluble_ag2      = cbalance%YCpool_nonsoluble_ag2(kstart:kend,1:ntiles), &
                  YCpool_acid_bg2            = cbalance%YCpool_acid_bg2      (kstart:kend,1:ntiles), &
                  YCpool_water_bg2           = cbalance%YCpool_water_bg2     (kstart:kend,1:ntiles), &
                  YCpool_ethanol_bg2         = cbalance%YCpool_ethanol_bg2   (kstart:kend,1:ntiles), &
                  YCpool_nonsoluble_bg2      = cbalance%YCpool_nonsoluble_bg2(kstart:kend,1:ntiles), &
                  YCpool_humus_2             = cbalance%YCpool_humus_2      (kstart:kend,1:ntiles), &
                  C_2_atmos                     = lcc%C2atmos                (kstart:kend), &
                  C_2_litterGreenPools          = C2litterGreenPools         (     1:nidx), &
                  C_2_litterWoodPool_ag         = C2litterWoodPool_ag        (     1:nidx), &
                  C_2_litterWoodPool_bg         = C2litterWoodPool_bg        (     1:nidx), &                        
                  Cpool_onSite                  = cbalance%Cpool_onSite      (kstart:kend), &
                  Cpool_paper                   = cbalance%Cpool_paper       (kstart:kend), &
                  Cpool_construction            = cbalance%Cpool_construction(kstart:kend), &
                  C_onSite_2_atmos              = lcc%C_onSite_2_atmos       (kstart:kend), &
                  C_paper_2_atmos               = lcc%C_paper_2_atmos        (kstart:kend), &
                  C_construction_2_atmos        = lcc%C_construction_2_atmos (kstart:kend), &                         
                  C_2_onSite                    = lcc%C_2_onSite             (kstart:kend), &
                  C_2_paper                     = lcc%C_2_paper              (kstart:kend), &
                  C_2_construction              = lcc%C_2_construction       (kstart:kend), &
                  lcc_scheme                    = lcc_scheme )
          END IF
       END SELECT

    END IF

    !! --- Perform carbon and nitrogen conservation test

    IF (debug_Cconservation) THEN
       IF (.NOT. with_yasso) THEN
          lcc%LCC_testCconserv(kstart:kend) = lcc%LCC_testCconserv(kstart:kend)  &
              - (  lcc%C2atmos(kstart:kend)                                      &
                 + SUM(  landcover_fract_new              (     1:nidx,1:ntiles) &
                       * veg_fract_correction             (     1:nidx,1:ntiles) &
                       * (  cbalance%Cpool_green          (kstart:kend,1:ntiles) &
                          + cbalance%Cpool_woods          (kstart:kend,1:ntiles) &
                          + cbalance%Cpool_reserve        (kstart:kend,1:ntiles) &
                          + cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles) &
                          + cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles) &
                          + cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles) &
                          + cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles) &
                          + cbalance%Cpool_slow           (kstart:kend,1:ntiles) &
                         ), DIM=2 )                                              &
                ) * veg_ratio_max(1:nidx)
       ELSE
          lcc%LCC_testCconserv(kstart:kend) = lcc%LCC_testCconserv(kstart:kend)  &
              - (  lcc%C2atmos(kstart:kend)                                      &
                 + SUM(  landcover_fract_new              (     1:nidx,1:ntiles) &
                       * veg_fract_correction             (     1:nidx,1:ntiles) &
                       * (  cbalance%Cpool_green          (kstart:kend,1:ntiles) &
                          + cbalance%Cpool_woods          (kstart:kend,1:ntiles) &
                          + cbalance%Cpool_reserve        (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_acid_ag1      (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_water_ag1     (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_ethanol_ag1   (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_nonsoluble_ag1(kstart:kend,1:ntiles) &
                          + cbalance%YCpool_acid_bg1      (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_water_bg1     (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_ethanol_bg1   (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_nonsoluble_bg1(kstart:kend,1:ntiles) &
                          + cbalance%YCpool_humus_1       (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_acid_ag2      (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_water_ag2     (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_ethanol_ag2   (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_nonsoluble_ag2(kstart:kend,1:ntiles) &
                          + cbalance%YCpool_acid_bg2      (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_water_bg2     (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_ethanol_bg2   (kstart:kend,1:ntiles) &
                          + cbalance%YCpool_nonsoluble_bg2(kstart:kend,1:ntiles) &
                          + cbalance%YCpool_humus_2       (kstart:kend,1:ntiles) &
                         ), DIM=2 )                                              &
                ) * veg_ratio_max(1:nidx)
       END IF

       IF (lcc_scheme==2) THEN
          lcc%LCC_testCconserv(kstart:kend) =                       &
          lcc%LCC_testCconserv(kstart:kend) - veg_ratio_max(1:nidx) &
               * (  cbalance%Cpool_onSite      (kstart:kend)        &
                  + cbalance%Cpool_paper       (kstart:kend)        &
                  + cbalance%Cpool_construction(kstart:kend)        &             
                 )
       ENDIF

       IF(with_nitrogen) then
          lcc%LCC_testNconserv(kstart:kend) = lcc%LCC_testNconserv(kstart:kend)   &
               - (  N2atmos(1:nidx)                                               &
                  + SUM(  landcover_fract_current(1:nidx,1:ntiles)                &
                        * veg_fract_correction   (1:nidx,1:ntiles)                &
                        * (  nbalance%Npool_green          (kstart:kend,1:ntiles) &
                           + nbalance%Npool_woods          (kstart:kend,1:ntiles) &
                           + nbalance%Npool_mobile         (kstart:kend,1:ntiles) &
                           + nbalance%Npool_litter_green_ag(kstart:kend,1:ntiles) &
                           + nbalance%Npool_litter_green_bg(kstart:kend,1:ntiles) &
                           + nbalance%Npool_litter_wood_ag (kstart:kend,1:ntiles) &
                           + nbalance%Npool_litter_wood_bg (kstart:kend,1:ntiles) &
                           + nbalance%Npool_slow           (kstart:kend,1:ntiles) &
                           + nbalance%SMINN_pool           (kstart:kend,1:ntiles) &
                          ) , DIM=2                                               &
                       )                                                          &
                 ) * veg_ratio_max(1:nidx) 
       END IF
    ENDIF

    IF (lcc_scheme==2) CALL CumulateAnthroFluxes(lcc, veg_ratio_max)

    IF (lcc_scheme==1) THEN  ! Standard JSBACH LCC 
       !! Calculate landcover change fluxes (conversion from [mol m-2 day-1] to [mol m-2 s-1])
       
       lcc%LCC_flux_box_C2atmos           (kstart:kend) = &
       lcc%LCC_flux_box_C2atmos           (kstart:kend) + lcc%C2atmos       (kstart:kend) * veg_ratio_max(1:nidx)/86400._dp
       lcc%LCC_flux_box_C2litterGreenPools(kstart:kend) = &
       lcc%LCC_flux_box_C2litterGreenPools(kstart:kend) + C2litterGreenPools(     1:nidx) * veg_ratio_max(1:nidx)/86400._dp
       lcc%LCC_flux_box_C2litterWoodPool  (kstart:kend) = &
       lcc%LCC_flux_box_C2litterWoodPool  (kstart:kend) +   (C2litterWoodPool_ag(1:nidx)+C2litterWoodPool_bg(1:nidx)) &
                                                          * veg_ratio_max(1:nidx)/86400._dp
       IF (with_nitrogen) THEN
          lcc%LCC_flux_box_N2atmos           (kstart:kend) = &
          lcc%LCC_flux_box_N2atmos           (kstart:kend) + N2atmos           (1:nidx) * veg_ratio_max(1:nidx)/86400._dp
          lcc%LCC_flux_box_N2litterGreenPools(kstart:kend) = &
          lcc%LCC_flux_box_N2litterGreenPools(kstart:kend) + N2litterGreenPools(1:nidx) * veg_ratio_max(1:nidx)/86400._dp
          lcc%LCC_flux_box_N2SMINNpool       (kstart:kend) = &
          lcc%LCC_flux_box_N2SMINNpool       (kstart:kend) + N2SMINNpool       (1:nidx) * veg_ratio_max(1:nidx)/86400._dp
          lcc%LCC_flux_box_N2litterWoodPool  (kstart:kend) = &
          lcc%LCC_flux_box_N2litterWoodPool  (kstart:kend) +   (N2litterWoodPool_ag(1:nidx)+N2litterWoodPool_bg(1:nidx)) &
                                                             * veg_ratio_max(1:nidx)/86400._dp
       END IF
    ENDIF

 !! Compute CO2 and nitrogen flux from landcover change: 
    !! .. i.e. conversion from [mol(C)/m^2(vegetated area)]during whole day to [kg(CO2)/m^2(grid box) s]
    CO2_emission(:) = lcc%C2atmos(kstart:kend) * veg_ratio_max(:) * molarMassCO2_kg/86400._dp

    IF (with_nitrogen) Nitrogen_2_atmos(:) = N2atmos(1:nidx)

    !! replace the old by the new landcover fractions

    landcover_fract_current(1:nidx,:) = landcover_fract_new(1:nidx,:)
    
  end subroutine do_landcover_change

  !! --- read_landcover_fractions() ----------------------------------------------------
  !!
  !! This routine reads (when necessary) the new landcover fractions into LCC_coverFract_target(:,:) 
  !!
  SUBROUTINE read_landcover_fractions(current_year, grid, domain, surface)

    USE mo_land_surface,     ONLY: scale_cover_fract, land_surface_type
    use mo_jsbach_grid,      only: grid_type,domain_type
    use mo_netcdf,           only: FILE_INFO,IO_inq_dimid,IO_inq_dimlen,IO_inq_varid,io_get_var_double
    USE mo_tr_scatter,       ONLY: scatter_field
    use mo_io,               only: IO_open, IO_READ, IO_close
    USE mo_temp,             ONLY: zreal2d, zreal3d, zreal2d_ptr

    integer,           intent(in)   :: current_year
    type(grid_type),   intent(in)   :: grid
    type(domain_type), intent(in)   :: domain
    TYPE(land_surface_type), INTENT(in)  :: surface

    !! --- locals

    character(len=1024)         :: filename
    integer                     :: IO_file_id, IO_var_id, IO_dim_id
    type(FILE_INFO)             :: IO_file
    integer                     :: znlon, znlat,zntiles
    integer                     :: i,status

    !! Generic name of land surface file with new landcover
    IF (current_year < 10000) THEN
       WRITE(filename,'(a,I4.4,a)') "cover_fract.", current_year, ".nc"
    ELSE
       CALL finish('read_landcover_fractions','Only years between 0 and 9999 supported currently')
    END IF

    if (p_parallel_io) then
       ! Open ini file
       call message('read_landcover_fractions','Reading new land surface fields from '//trim(filename))
       IO_file%opened = .false.
       call IO_open(trim(filename), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Check resolution
       call IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       call IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)
! Select one of the two lines below dependent on your data
       call IO_inq_dimid  (IO_file_id, 'ntiles', IO_dim_id)
!       call IO_inq_dimid  (IO_file_id, 'tiles', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, zntiles)

       if (znlon /= grid%nlon .or. znlat /= grid%nlat) then
          write (message_text,*) 'Unexpected grid resolution: ', znlon, 'x', znlat
          call finish('read_landcover_fractions()', message_text)
       endif
       if (zntiles /= surface%ntiles) then
          write (message_text,*) 'Unexpected number of tiles: ', zntiles
          call finish('read_landcover_fractions()', message_text)
       endif

       !! allocate temporary memory

       allocate(zreal3d(grid%nlon,grid%nlat,surface%ntiles),STAT=status)
       if(status .ne. 0) call finish('read_landcover_fractions()','Allocation failure (1)')

       !! read cover fractions

       call IO_inq_varid(IO_file_id, trim(coverFractVarName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal3d)

       call IO_close(IO_file)
    endif !! end parallel_io

    !! Bring cover fractions to the other processors

    allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_landcover_fractions()','Allocation failure (3)')
    nullify(zreal2d_ptr)
    do i=1,surface%ntiles
       if (p_parallel_io) zreal2d_ptr => zreal3d(:,:,i)
       call scatter_field(zreal2d_ptr, zreal2d)
       landcover_change%LCC_coverFract_target(:,i) = pack(zreal2d, MASK=domain%mask)
    enddo
    CALL scale_cover_fract(domain%nland, surface%ntiles, surface%is_present(:,:), surface%is_glacier(:,:), &
         surface%is_naturalveg(:,:), landcover_change%LCC_coverFract_target(:,:))
    
    !! Free temporary memory
    if (p_parallel_io) deallocate(zreal3d)
    deallocate(zreal2d)

  end subroutine read_landcover_fractions

  !! --- do_landuse_transitions() ----------------------------------------------------
  !!
  !! This routine computes every day the new landuse transitions for all PFTs from the 
  !! (reduced) landuse transitions of the New Hamphire Harmonized protocol (NHHLP). The term "reduced" 
  !! means here that instead of the "primary" and "secondary" lands of the original harmonized 
  !! protocol, here "primary" and "secondary" lands are collectively handled as "natural" lands. 
  !! See the JSBACH-documentation for more details.
  !!
  SUBROUTINE do_landuse_transitions(lctlib, surface, ntiles, is_naturalVeg, is_C4vegetation, is_grass, is_crop, is_pasture,      &
                                   is_vegetation, is_glacier, veg_ratio_max, veg_fract_correction, landcover_fract_pot,          &
                                   landcover_fract_current, cbalance, CO2_emission_LUtrans, CO2_emission_harvest,                &
                                   lcc, with_yasso, LeafLit_coef, WoodLit_coef,nbalance,with_nitrogen, lcc_scheme)

    USE mo_jsbach_grid,      ONLY: kstart, kend, nidx
    USE mo_jsbach_lctlib,    ONLY: lctlib_type
    USE mo_time_control,     ONLY: current_date, lstart
    USE mo_time_conversion,  ONLY: year_len
    USE mo_land_surface,     ONLY: fract_small
    USE mo_cbal_bethy,       ONLY: cbalance_type, nbalance_type
    USE mo_cbal_cpools,      ONLY: C_relocation_from_LUtransitions, C_relocation_from_harvest
    USE mo_jsbach_constants, ONLY: molarMassCO2_kg
    USE mo_land_surface,     ONLY: land_surface_type

    TYPE(land_surface_type), INTENT(in)  :: surface
    TYPE(lctlib_type),       INTENT(in)  :: lctlib  !! PTF-specific parameters
    INTEGER,  INTENT(in)    :: ntiles               !! number of tiles of the land points
    LOGICAL,  INTENT(in)    :: is_naturalVeg(:,:)   !! mask indicating natural vegetation (nidx x ntiles)
    LOGICAL,  INTENT(in)    :: is_C4vegetation(:,:) !! mask indicating photosynthetic pathway (true: C4; false: C3) (nidx x ntiles)
    LOGICAL,  INTENT(in)    :: is_grass(:,:)        !! mask indicating grassland (nidx x ntiles)
    LOGICAL,  INTENT(in)    :: is_crop(:,:)         !! mask indicating cropland (nidx x ntiles)
    LOGICAL,  INTENT(in)    :: is_pasture(:,:)      !! mask indicating pasture (nidx x ntiles)
    LOGICAL,  INTENT(in)    :: is_vegetation(:,:)   !! mask indicating vegetated tiles (nidx x ntiles)
    LOGICAL,  INTENT(in)    :: is_glacier(:,:)      !! mask indicating glacier tiles (nidx x ntiles)
    REAL(dp), INTENT(in)    :: veg_ratio_max(:)     !! maximum fraction of grid cell that can be covered by vegetation
    REAL(dp), INTENT(in)    :: veg_fract_correction(:,:)        !! Correction factor for cover fractions 1-exp(-LAI_max/2)
    REAL(dp), INTENT(in)    :: landcover_fract_pot(:,:)         !! Potential landcover fractions, i.e. potential maximum extent of  
                                                                !! .. natural PFTs in absence in absence of agricultural types. 
                                                                !! .. Note that these are fractions relative to the extent of
                                                                !! .. vegetation in a grid cell (sum over tiles is 1)
    REAL(dp), INTENT(inout) :: landcover_fract_current(:,:)     !! On call the old landcover fractions, on return the new ones
    TYPE(cbalance_type),         INTENT(inout) :: cbalance      !! The carbon pools to be changed along with the LU-transitions
    REAL(dp), INTENT(out)   :: CO2_emission_LUtrans(:)          !! The CO2 emission from land use transitions in kg(CO2)/(m^2 s)
    REAL(dp), INTENT(out)   :: CO2_emission_harvest(:)          !! The CO2 emission from harvest in kg(CO2)/(m^2 s)
    TYPE(landcover_change_type), INTENT(inout) :: lcc
    LOGICAL, INTENT(in)           :: with_yasso                 !! YASSO pools                      
    REAL(dp), INTENT(in)          :: LeafLit_coef(:,:,:)        !! Factor to seperate non woody litter into yasso pools[ ]          
    REAL(dp), INTENT(in)          :: WoodLit_coef(:,:,:)        !! Factor to seperate woody litter into yasso pools    [ ]
    TYPE(nbalance_type), INTENT(inout) :: nbalance
    LOGICAL, INTENT(in)     :: with_nitrogen
    INTEGER, OPTIONAL, INTENT(in) :: lcc_scheme

    !! locals

    !! Abreviations to increase verbosity of code related to the New Hampshire Harmonized Protocol of landuse transitions (NHHLP)

    INTEGER, PARAMETER :: CROP = 1 !! verbose index for crop types of landcover
    INTEGER, PARAMETER :: PAST = 2 !! verbose index for pasture types of landcover
    INTEGER, PARAMETER :: GRAS = 3 !! verbose index for natural graslands
    INTEGER, PARAMETER :: FRST = 4 !! verbose index for natural non-grasland types, summarized here as FoReST, although this class
                                   !! .. contains also shrub or tundra types.

    REAL(dp), DIMENSION(1:nidx) :: cf0_CROP, cf0_PAST, cf0_GRAS, cf0_FRST !! NHHLP-cover fractions in old year for crops, pastures,
                                                                          !! .. grass, and natural non-grasland (FRST) cover
    REAL(dp), DIMENSION(1:nidx) :: cf1_NATL                            !! NHHLP-cover fractions to be reached at end of new year for
                                                                       !! .. natural vegetation cover
    REAL(dp), DIMENSION(1:nidx) :: cf_pot_GRAS, cf_pot_FRST            !! Potential NHHLP-cover fractions in new year
    REAL(dp), DIMENSION(1:nidx) :: d_NATL_2_PAST, d_GRAS_2_PAST        !! Fractions of land cover converted from type A to type B

    REAL(dp) :: DailyTransMtrx(1:nidx,1:4,1:4)  !! Transition matrices of (extended) NHHLP, i.e. natural vegetation is splitted 
                                                !! .. into grasslands (GRAS) and non-grasslands (FRST).
                                                !! .. Daily_TransMtrx(n,i,j) is the fraction of area of cover type j that in 
                                                !! ..  grid box n is converted into cover type i.
 
    real(dp) :: Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles)   !! Matrices of the daily landuse transitions between PFTs at the 
                                                           !! .. different tiles as derived from the Daily_TransMtrx(). 
                                                           !! .. Tile_TransMtrx(n,i,j) is the fraction of area of PFT j that in 
                                                           !! .. grid box n is converted into cover type i.
    REAL(dp) :: landcover_fract_new(1:nidx,1:ntiles)       !! helper array for new cover fractions

    INTEGER  :: noOfDays_in_current_year
    REAL(dp) :: convFactor(1:nidx)

    REAL(dp) :: C2litterWoodPool_ag(1:nidx)
    REAL(dp) :: C2litterWoodPool_bg(1:nidx)
    REAL(dp) :: C2litterGreenPools(1:nidx)

    REAL(dp) :: C2litter_harvest(1:nidx)                   !! Carbon from harvest released to the below ground litter pools
                                                           !! .. (more realsitically this carbon should be put into anthropogenic
                                                           !! .. pools, but such pools do not exist in JSBACH, therefore this 
                                                           !! .. carbon is put into litter pools of similar lifetime)
    !! helpers

    REAL(dp) :: rhlp(1:nidx),ahlp(1:nidx,1:ntiles)
    INTEGER  :: n, i, j, opt

    !! === GO ========================================

    !! Set model running mode

    !! --- Check for exiting the routine without landcover change
    IF (.NOT. new_day .AND. .NOT. lstart) THEN !! We are not at the first time step of a day ..
       !! ... so we use the same fluxes during the whole day ...
       CO2_emission_LUtrans(:) = landuse_transitions%C2atmos_LUtrans(kstart:kend) * veg_ratio_max(:) * molarMassCO2_kg / 86400._dp
       CO2_emission_harvest(:) = landuse_transitions%C2atmos_harvest(kstart:kend) * veg_ratio_max(:) * molarMassCO2_kg / 86400._dp
       RETURN                                                !! .. and exit the routine
    END IF

    !! A fundamental check
    IF(ANY( (landcover_fract_current(:,:)-landcover_fract_pot(:,:)) > ntiles * epsilon(1._dp)    &
                                     .AND. is_naturalVeg(:,:) )) THEN !! something went wrong (wrong initial cover fractions?)
       ahlp(1:nidx,1:ntiles) = 0.0_dp
       WHERE( is_naturalVeg(:,:)) ahlp(1:nidx,1:ntiles)= landcover_fract_current(:,:)-landcover_fract_pot(:,:)
       
       WRITE (message_text, *) ' Actual cover_fract larger than potential cover_fract at (', &
            MAXLOC(ahlp(:,:)),') difference is ',MAXVAL(ahlp(:,:))
       CALL finish('do_landuse_transitions',message_text)
    END IF

    !! Assure cover fractions of natural lands to be not larger than potential vegetation cover up to numerical accuracy
    WHERE (is_naturalVeg(1:nidx,1:ntiles) .AND. &
           landcover_fract_current(1:nidx,1:ntiles) > landcover_fract_pot(1:nidx,1:ntiles) ) 
       landcover_fract_current(1:nidx,1:ntiles) = landcover_fract_pot(1:nidx,1:ntiles)
    END WHERE

    !! === Update of landcover transitions
 
    !! Determine NHHLP cover fractions (GRAS, FRST, PAST, CROP) at last time step
    
    cf0_CROP(1:nidx) = & !! Fraction of grid cell covered by crops
         SUM(landcover_fract_current(1:nidx,1:ntiles), &
                                           DIM=2,MASK=is_crop(1:nidx,1:ntiles))

    cf0_PAST(1:nidx) = & !! Fraction of grid cell covered by pastures
         SUM(landcover_fract_current(1:nidx,1:ntiles), &
                                           DIM=2,MASK=is_pasture(1:nidx,1:ntiles))

    cf0_GRAS(1:nidx) = &  !! Fraction of grid cell covered by natural grasslands
         SUM(landcover_fract_current(1:nidx,1:ntiles), &
                                           DIM=2,MASK=is_grass(1:nidx,1:ntiles))

    cf_pot_GRAS(1:nidx) = &  !! Fraction of grid cell covered potentially by natural grasslands
         SUM(landcover_fract_pot(1:nidx,1:ntiles), &
                                           DIM=2,MASK=is_grass(1:nidx,1:ntiles))    

    cf0_FRST(1:nidx) = &  !! Fraction of grid cell covered by natural non-grassland types
         SUM(landcover_fract_current(1:nidx,1:ntiles), &
                                          DIM=2,MASK=is_naturalVeg(1:nidx,1:ntiles) .AND. .NOT. is_grass(1:nidx,1:ntiles))

    cf_pot_FRST(1:nidx) = &  !! Fraction of grid cell covered potentially by natural non-grassland types
         SUM(landcover_fract_pot(1:nidx,1:ntiles), &
                                          DIM=2,MASK=is_naturalVeg(1:nidx,1:ntiles) .AND. .NOT. is_grass(1:nidx,1:ntiles))

    IF (new_year) THEN

       !! Save cover fractions that were reached at the end of the last year

       landuse_transitions%Grass_coverFract_lastYear(kstart:kend)       = cf0_GRAS(1:nidx)
       landuse_transitions%NatWood_coverFract_lastYear(kstart:kend)     = cf0_FRST(1:nidx)
       landuse_transitions%Pasture_coverFract_lastYear(kstart:kend)     = cf0_PAST(1:nidx)
       landuse_transitions%Crop_coverFract_lastYear(kstart:kend)        = cf0_CROP(1:nidx)


       !! Check whether extent of natural lands exceeds potential extent -- in such a case modify original transition matrix
       !! .. elements of NHHLP to match potential extent

       cf1_NATL(1:nidx) = cf0_FRST(1:nidx) + cf0_GRAS(1:nidx)&
            + landuse_transitions%TransMtrx_CROP_2_NATL(kstart:kend) * cf0_CROP(1:nidx) &
            + landuse_transitions%TransMtrx_PAST_2_NATL(kstart:kend) * cf0_PAST(1:nidx) &
            - (  landuse_transitions%TransMtrx_NATL_2_CROP(kstart:kend) &
               + landuse_transitions%TransMtrx_NATL_2_PAST(kstart:kend) &
              ) * (cf0_GRAS(1:nidx) + cf0_FRST(1:nidx))

       !! CHR --- code for checking is missing

       !!
       !! === Use elements of 3x3 transition matrix to determine the missing elements of the 4x4 transition matrix
       !!

       !! 1. EXPANSION OF AGRICULTURAL LANDS =====================================================================

       !! 1.1 Conversion of natural lands to pasture - establishment of pastures on grasslands has priority  

       !! Determine cover fraction of grid box converted from natural lands (GRAS+FRST) to pasture (PAST):
       d_NATL_2_PAST(1:nidx) = &
            landuse_transitions%TransMtrx_NATL_2_PAST(kstart:kend) * ( cf0_GRAS(1:nidx) + cf0_FRST(1:nidx) )
       !! Compute matrix element for the transition from grasslands to pastures (but: prevent TransMtrx>1 and divisions by zero)
       landuse_transitions%TransMtrx_GRAS_2_PAST(kstart:kend) = &
            MIN( 1._dp , d_NATL_2_PAST(1:nidx)/MAX(cf0_GRAS(1:nidx),fract_small) )
       !! Compute conversion fraction GRAS->PAST
       d_GRAS_2_PAST(1:nidx) = landuse_transitions%TransMtrx_GRAS_2_PAST(kstart:kend) * cf0_GRAS(1:nidx)
       !! Compute matrix element for the transition from non-grasslands (FRST) to pasture by distributing the rest of d_NATL_2_PAST:
       landuse_transitions%TransMtrx_FRST_2_PAST(kstart:kend) = &
            MIN( 1.0_dp , (d_NATL_2_PAST(1:nidx) - d_GRAS_2_PAST(1:nidx))/MAX(cf0_FRST(1:nidx),fract_small) )
       landuse_transitions%TransMtrx_FRST_2_PAST(kstart:kend) = &
            MAX(0._dp,landuse_transitions%TransMtrx_FRST_2_PAST(kstart:kend))

       !! 1.2 Conversion of natural lands to crops

       rhlp(1:nidx) =  landuse_transitions%TransMtrx_NATL_2_CROP(kstart:kend) * ( cf0_GRAS(1:nidx) + cf0_FRST(1:nidx) ) &
                       / MAX( cf0_GRAS(1:nidx) + cf0_FRST(1:nidx) - d_NATL_2_PAST(1:nidx) , fract_small)
       !! Compute matrix element for the transition from GRAS to CROP (but: prevent TransMtrx>1 and divisions by zero)
       landuse_transitions%TransMtrx_GRAS_2_CROP(kstart:kend) = &
             rhlp(1:nidx) * (1._dp - landuse_transitions%TransMtrx_GRAS_2_PAST(kstart:kend))
       !! Compute matrix element for the transition from non-grasslands (FRST) to PAST by distributing the rest of d_NATL_2_CROP:
       landuse_transitions%TransMtrx_FRST_2_CROP(kstart:kend) =  &
             rhlp(1:nidx) * (1._dp - landuse_transitions%TransMtrx_FRST_2_PAST(kstart:kend))


       !! 2. ABANDONING OF AGRICULTURAL LANDS (operate in reverse order as compared to 1.) =========================


       !! 2.1 Back-conversion of croplands to natural lands

       rhlp(1:nidx) = landuse_transitions%TransMtrx_CROP_2_NATL(kstart:kend) &
                    / MAX( cf_pot_GRAS(1:nidx) + cf_pot_FRST(1:nidx) - cf0_GRAS(1:nidx) - cf0_FRST(1:nidx), fract_small)
       rhlp(1:nidx) = MIN(1._dp/max(cf0_CROP(1:nidx),fract_small),rhlp(1:nidx)) 
       !! Compute matrix element for the transition from CROP to FRST (but: prevent TransMtrx>1 and divisions by zero)
       landuse_transitions%TransMtrx_CROP_2_GRAS(kstart:kend) = rhlp(1:nidx) * ( cf_pot_GRAS(1:nidx) - cf0_GRAS(1:nidx) )
       !! Compute matrix element for the transition from CROP to GRAS (prevent TransMtrx<0)
       landuse_transitions%TransMtrx_CROP_2_FRST(kstart:kend) = rhlp(1:nidx) * ( cf_pot_FRST(1:nidx) - cf0_FRST(1:nidx) )
            
       !! 2.2 Back-conversion of pastures to natural lands

       !! Compute matrix element for the transition from PAST to FRST (prevent divisions by zero)
       landuse_transitions%TransMtrx_PAST_2_FRST(kstart:kend) = &
            MIN( landuse_transitions%TransMtrx_PAST_2_NATL(kstart:kend) , &
                 MAX(0._dp, cf_pot_FRST(1:nidx) - cf0_FRST(1:nidx) &
                          - landuse_transitions%TransMtrx_CROP_2_FRST(kstart:kend) * cf0_CROP(1:nidx) ) /      &
                                                                         max(cf0_PAST(1:nidx),fract_small) &
               )
       landuse_transitions%TransMtrx_PAST_2_FRST(kstart:kend) = &
            MAX(0._dp,landuse_transitions%TransMtrx_PAST_2_FRST(kstart:kend))
       !! Compute matrix element for the transition from PAST to GRAS (prevent TransMtrx<0)
       landuse_transitions%TransMtrx_PAST_2_GRAS(kstart:kend) =                                                        &
            MAX( 0._dp ,                                                                                               &
                 MIN( landuse_transitions%TransMtrx_PAST_2_NATL(kstart:kend)                                           &
                                            - landuse_transitions%TransMtrx_PAST_2_FRST(kstart:kend),                  &
                      ( cf_pot_GRAS(1:nidx) - cf0_GRAS(1:nidx)                                                         &
                                            - landuse_transitions%TransMtrx_CROP_2_GRAS(kstart:kend) *cf0_CROP(1:nidx) &
                      ) / MAX(cf0_PAST(1:nidx),fract_small) )                                                      &
               )
       !! Diagnose inconsistency with New Hampshire protocol: Compute cover fraction that could not be converted  
       !! .. because of non-availability of potential GRAS or FRST area. Note: this is computed only once a year!

       landuse_transitions%CROP_2_NATL_ignored(kstart:kend) =                       &
                       ( landuse_transitions%TransMtrx_CROP_2_NATL(kstart:kend)     &
                         - landuse_transitions%TransMtrx_CROP_2_FRST(kstart:kend)   &
                           - landuse_transitions%TransMtrx_CROP_2_GRAS(kstart:kend) &
                       ) * cf0_CROP(1:nidx)

       landuse_transitions%PAST_2_NATL_ignored(kstart:kend) =                   &
                       ( landuse_transitions%TransMtrx_PAST_2_NATL(kstart:kend)     &
                         - landuse_transitions%TransMtrx_PAST_2_FRST(kstart:kend)   &
                           - landuse_transitions%TransMtrx_PAST_2_GRAS(kstart:kend) &
                       ) * cf0_PAST(1:nidx)

    END IF

    IF (debug_Cconservation) THEN

       !! Check consistency of 4x4 with original 3x3 transition matrix

       landuse_transitions%Test_NATL_2_PAST(kstart:kend) = 0.0_dp
       landuse_transitions%Test_NATL_2_CROP(kstart:kend) = 0.0_dp
       landuse_transitions%Test_PAST_2_NATL(kstart:kend) = 0.0_dp
       landuse_transitions%Test_CROP_2_NATL(kstart:kend) = 0.0_dp

       landuse_transitions%Test_NATL_2_PAST(kstart:kend) =                            &
              landuse_transitions%TransMtrx_GRAS_2_PAST(kstart:kend)                  &
                 * landuse_transitions%Grass_coverFract_lastYear(kstart:kend)         &
            + landuse_transitions%TransMtrx_FRST_2_PAST(kstart:kend)                  &
                 * landuse_transitions%NatWood_coverFract_lastYear(kstart:kend)       &
            - landuse_transitions%TransMtrx_NATL_2_PAST(kstart:kend)                  &
                 * (   landuse_transitions%Grass_coverFract_lastYear(kstart:kend)     &
                     + landuse_transitions%NatWood_coverFract_lastYear(kstart:kend) )

       landuse_transitions%Test_NATL_2_CROP(kstart:kend) =                            &
              landuse_transitions%TransMtrx_GRAS_2_CROP(kstart:kend)                  &
                 * landuse_transitions%Grass_coverFract_lastYear(kstart:kend)         &
            + landuse_transitions%TransMtrx_FRST_2_CROP(kstart:kend)                  &
                 * landuse_transitions%NatWood_coverFract_lastYear(kstart:kend)       &
            - landuse_transitions%TransMtrx_NATL_2_CROP(kstart:kend)                  &
                 * (   landuse_transitions%NatWood_coverFract_lastYear(kstart:kend)   &
                     + landuse_transitions%Grass_coverFract_lastYear(kstart:kend)  )

       landuse_transitions%Test_PAST_2_NATL(kstart:kend) =                            &
            ( landuse_transitions%TransMtrx_PAST_2_FRST(kstart:kend)                  &
            + landuse_transitions%TransMtrx_PAST_2_GRAS(kstart:kend)                  &
            - landuse_transitions%TransMtrx_PAST_2_NATL(kstart:kend)                  &
            ) * landuse_transitions%Pasture_coverFract_lastYear(kstart:kend)          &
            + landuse_transitions%PAST_2_NATL_ignored(kstart:kend)  !! <-- correct for ignored fraction

       landuse_transitions%Test_CROP_2_NATL(kstart:kend) =                            &
            ( landuse_transitions%TransMtrx_CROP_2_FRST(kstart:kend)                  &
            + landuse_transitions%TransMtrx_CROP_2_GRAS(kstart:kend)                  &
            - landuse_transitions%TransMtrx_CROP_2_NATL(kstart:kend)                  &
            ) * landuse_transitions%Crop_coverFract_lastYear(kstart:kend)             &
            + landuse_transitions%CROP_2_NATL_ignored(kstart:kend)  !! <-- correct for ignored fraction

    END IF

    !!
    !! === Compute the 4x4 elements of the landuse transition matrix for daily values
    !!

    noOfDays_in_current_year = year_len(current_date)

    !! GRAS -> OTHER
    convFactor(1:nidx) = landuse_transitions%Grass_coverFract_lastYear(kstart:kend) &
                         / max(cf0_GRAS(1:nidx),fract_small)/noOfDays_in_current_year
    convFactor(1:nidx) = MIN(1._dp,convFactor(1:nidx))

    DailyTransMtrx(1:nidx,PAST,GRAS) =  landuse_transitions%TransMtrx_GRAS_2_PAST(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,CROP,GRAS) =  landuse_transitions%TransMtrx_GRAS_2_CROP(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,FRST,GRAS) =  0.0_dp  
    DailyTransMtrx(1:nidx,GRAS,GRAS) = &
         MAX(0._dp,1._dp - DailyTransMtrx(1:nidx,PAST,GRAS) -DailyTransMtrx(1:nidx,CROP,GRAS))

    !! FRST -> OTHER
    convFactor(1:nidx) = landuse_transitions%NatWood_coverFract_lastYear(kstart:kend) &
                         / max(cf0_FRST(1:nidx),fract_small)/noOfDays_in_current_year
    convFactor(1:nidx) = MIN(1._dp,convFactor(1:nidx))

    DailyTransMtrx(1:nidx,PAST,FRST) =  landuse_transitions%TransMtrx_FRST_2_PAST(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,CROP,FRST) =  landuse_transitions%TransMtrx_FRST_2_CROP(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,GRAS,FRST) =  0.0_dp  
    DailyTransMtrx(1:nidx,FRST,FRST) = &
         MAX(0._dp,1._dp - DailyTransMtrx(1:nidx,PAST,FRST) - DailyTransMtrx(1:nidx,CROP,FRST))

    !! CROP -> OTHER
    convFactor(1:nidx) = landuse_transitions%Crop_coverFract_lastYear(kstart:kend) &
                         / max(cf0_CROP(1:nidx),fract_small)/noOfDays_in_current_year
    convFactor(1:nidx) = MIN(1._dp,convFactor(1:nidx))

    DailyTransMtrx(1:nidx,PAST,CROP) =  landuse_transitions%TransMtrx_CROP_2_PAST(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,GRAS,CROP) =  landuse_transitions%TransMtrx_CROP_2_GRAS(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,FRST,CROP) =  landuse_transitions%TransMtrx_CROP_2_FRST(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,CROP,CROP) = &
         MAX(0._dp,1._dp - DailyTransMtrx(1:nidx,PAST,CROP) - DailyTransMtrx(1:nidx,GRAS,CROP) - DailyTransMtrx(1:nidx,FRST,CROP))

    !! PAST -> OTHER
    convFactor(1:nidx) = landuse_transitions%Pasture_coverFract_lastYear(kstart:kend) &
                         / max(cf0_PAST(1:nidx),fract_small)/noOfDays_in_current_year
    convFactor(1:nidx) = MIN(1._dp,convFactor(1:nidx))

    DailyTransMtrx(1:nidx,GRAS,PAST) =  landuse_transitions%TransMtrx_PAST_2_GRAS(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,CROP,PAST) =  landuse_transitions%TransMtrx_PAST_2_CROP(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,FRST,PAST) =  landuse_transitions%TransMtrx_PAST_2_FRST(kstart:kend) * convFactor(1:nidx)
    DailyTransMtrx(1:nidx,PAST,PAST) = &
         MAX(0._dp,1._dp - DailyTransMtrx(1:nidx,GRAS,PAST) - DailyTransMtrx(1:nidx,CROP,PAST) - DailyTransMtrx(1:nidx,FRST,PAST))


    !!
    !! === Derive from 4x4 transition matrix DailyTransMtrx() the ntiles x ntiles transition matrix Tile_TransMtrx()
    !!

    !! Initialize  Tile_TransMtrx() to zero so that only non-zero transitions have to be considered below
    !! (Glacier points are implicitely handled: For these points the matrix is getting a unit matrix.)
    Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles) = 0._dp

    !! 0. CONVERSIONS BETWEEN AGRICULTURAL LANDS

    !! CROP (i) -> PAST (j)
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles, is_Crop(n,i) .AND. is_Pasture(n,j))
       Tile_TransMtrx(n,j,i)=&
            DailyTransMtrx(n,PAST,CROP) * landcover_fract_current(n,j) &
                                        / MAX(cf0_PAST(n),fract_small)
    END FORALL

    !! PAST (i) -> CROP (j)
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles,is_Pasture(n,i) .AND. is_Crop(n,j) )
       Tile_TransMtrx(n,j,i)=&
            DailyTransMtrx(n,CROP,PAST) * landcover_fract_current(n,j) &
                                        / MAX(cf0_CROP(n),fract_small)
    END FORALL

    !! 1. EXPANSION OF AGRICULTURAL LANDS

    !! C4 GRAS (i) -> C4 PAST (j)
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles, &
                    is_Grass(n,i) .AND. is_C4vegetation(n,i) .AND. is_Pasture(n,j) .AND. is_C4vegetation(n,j))
       Tile_TransMtrx(n,j,i)=DailyTransMtrx(n,PAST,GRAS)
    END FORALL
    !! C3 GRAS (i) -> C3 PAST (j)
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles, &
                    is_Grass(n,i) .AND. (.NOT. is_C4vegetation(n,i)) .AND. is_Pasture(n,j) .AND. (.NOT. is_C4vegetation(n,j)))
       Tile_TransMtrx(n,j,i)=DailyTransMtrx(n,PAST,GRAS)
    END FORALL

    !! GRAS (i) ->  CROP (j)
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles, is_Grass(n,i) .AND. is_Crop(n,j))
       Tile_TransMtrx(n,j,i)=&
            DailyTransMtrx(n,CROP,GRAS) * landcover_fract_current(n,j) &
                                        / MAX(cf0_CROP(n),fract_small)
    END FORALL
    !! FRST (i) ->  PAST (j)
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles,  ( is_NaturalVeg(n,i) .AND. .NOT. is_Grass(n,i) ) .AND. is_pasture(n,j) )
       Tile_TransMtrx(n,j,i)=&
            DailyTransMtrx(n,PAST,FRST) * landcover_fract_current(n,j) &
                                        / MAX(cf0_PAST(n),fract_small)
    END FORALL
    !! FRST (i) ->  CROP (j)
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles, ( is_NaturalVeg(n,i) .AND. .NOT. is_Grass(n,i) ) .AND. is_crop(n,j)  )
       Tile_TransMtrx(n,j,i)=&
            DailyTransMtrx(n,CROP,FRST) * landcover_fract_current(n,j) &
                                        / MAX(cf0_CROP(n),fract_small)
    END FORALL

    !! 2. SHRINKAGE OF AGRICULTURAL LANDS

    !! C4 PAST (i) -> C4 GRASS (j): convert only C4 pasture to C4 grass
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles,&
                          is_Pasture(n,i) .AND. is_C4vegetation(n,i) .AND. is_Grass(n,j) .AND. is_C4vegetation(n,j))
       Tile_TransMtrx(n,j,i)= DailyTransMtrx(n,GRAS,PAST)
       Tile_TransMtrx(n,j,i)= &  !! Assure that back-conversion does not surmount potential vegetation are
            MIN( Tile_TransMtrx(n,j,i), &
                     ( MAX(0.0_dp,landcover_fract_pot(n,j)-landcover_fract_current(n,j)) )/ &
                                                                      MAX(landcover_fract_current(n,i),fract_small))
    END FORALL
    !! C3 PAST (i) -> C3 GRASS (j): convert only C3 pasture to C3 grass
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles,&
                 is_Pasture(n,i) .AND. (.NOT. is_C4vegetation(n,i)) .AND. is_Grass(n,j) .AND. (.NOT. is_C4vegetation(n,j)))
       Tile_TransMtrx(n,j,i)= DailyTransMtrx(n,GRAS,PAST)
       Tile_TransMtrx(n,j,i)= &  !! Assure that back-conversion does not surmount potential vegetation are
            MIN( Tile_TransMtrx(n,j,i),&
                 ( landcover_fract_pot(n,j)-landcover_fract_current(n,j) )/ MAX(landcover_fract_current(n,i),fract_small) )
    END FORALL

    !! CROP (i) -> GRAS (j): convert CROP to GRASS proportional to potential extent of GRASS left
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles, is_crop(n,i) .AND. is_Grass(n,j))
       Tile_TransMtrx(n,j,i)= &
            DailyTransMtrx(n,GRAS,CROP) * (landcover_fract_pot(n,j) - landcover_fract_current(n,j)) &
                                        / MAX(cf_pot_GRAS(n) - cf0_GRAS(n),fract_small)
       Tile_TransMtrx(n,j,i)= &  !! Assure that back-conversion does not surmount potential vegetation are
            MIN( Tile_TransMtrx(n,j,i),&
                 ( landcover_fract_pot(n,j)-landcover_fract_current(n,j) )/ MAX(landcover_fract_current(n,i),fract_small) )
    END FORALL

    !! PAST (i) -> FRST (j): convert PAST to  FRST proportional to potential extent of FRST left
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles, is_pasture(n,i) .AND.  ( is_NaturalVeg(n,j) .AND. .NOT. is_Grass(n,j) ))
       Tile_TransMtrx(n,j,i)= &
            DailyTransMtrx(n,FRST,PAST) * (landcover_fract_pot(n,j)-landcover_fract_current(n,j))&
                                        / MAX(cf_pot_FRST(n) - cf0_FRST(n),fract_small)
       Tile_TransMtrx(n,j,i)= &  !! Assure that back-conversion does not surmount potential vegetation are
            MIN( Tile_TransMtrx(n,j,i),&
                 ( landcover_fract_pot(n,j)-landcover_fract_current(n,j) )/ MAX(landcover_fract_current(n,i),fract_small) )
    END FORALL

    !! CROP (i) -> FRST (j): convert CROP to  FRST proportional to potential extent of FRST left
    FORALL(n=1:nidx,i=1:ntiles,j=1:ntiles, is_crop(n,i) .AND.  ( is_NaturalVeg(n,j) .AND. .NOT. is_Grass(n,j) ))
       Tile_TransMtrx(n,j,i)= &
            DailyTransMtrx(n,FRST,CROP) * (landcover_fract_pot(n,j)-landcover_fract_current(n,j))&
                                        / MAX(cf_pot_FRST(n) - cf0_FRST(n),fract_small)
       Tile_TransMtrx(n,j,i)= &  !! Assure that back-conversion does not surmount potential vegetation are
            MIN( Tile_TransMtrx(n,j,i),&
                 ( landcover_fract_pot(n,j)-landcover_fract_current(n,j) )/ MAX(landcover_fract_current(n,i),fract_small) )
    END FORALL

    !! Compute diagonal elements of tile transition matrix (because the matrix was set to zero above, the
    !! .. matrix is getting a unit matrix at glacier points)

    do i=1,ntiles
       Tile_TransMtrx(1:nidx,i,i)= 0.0_dp !! this is to exclude diagonal elements from the summation below
       Tile_TransMtrx(1:nidx,i,i)= 1.0_dp - SUM(Tile_TransMtrx(1:nidx,1:ntiles,i),DIM=2)
    end do

    !! Do some checks for the matrix

    IF(ANY( Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles) > 1._dp+10._dp*EPSILON(1.0_dp))) THEN
       WRITE (message_text, *) ' Tile_TransMtrx(',MAXLOC(Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles)),')=',&
            MAXVAL(Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles)),' is larger than 1!!'
          CALL message ('do_landuse_transitions', message_text)
    END IF

    IF(ANY( Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles) < -10._dp*EPSILON(1.0_dp))) THEN
       WRITE (message_text, *) ' Tile_TransMtrx(',MINLOC(Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles)),')=',&
            MINVAL(Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles)),' is smaller than 0!!'
          CALL message ('do_landuse_transitions', message_text)
    END IF

    !! Assure Matrixelements to be in [0,1] (remove possible numerical inaccuracies in above computations)
    Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles) = MAX(0.0_dp,MIN(1.0_dp,Tile_TransMtrx(1:nidx,1:ntiles,1:ntiles)))


    IF (debug_Cconservation) THEN

       !! Check consistency with 4x4 transition matrix
       !!
       !! Note that the test is not applicable (i) in cases where in the computation of the tile-transition matrices
       !! .. the division by cover fractions (or difference between potential and actual cover fraction) is
       !! .. cut off to prevent division by zero. In these cases the test arrays computed below are set to
       !! .. zero to prevent non-zero values to be misinterpreted as a numerical inconsistency.
       !! .. The test is also not applicable when in the computation of the tile-transition matrices in back-conversion
       !! .. cases the transition is reduced so that potential vegetation extent is not exceeded. 
       !! .. THIS LATTER CASE IS NOT CORRECTED FOR IN THE TEST, SO THAT A NON_ZERO VALUE MEANS ALL CASES FAILURE!!

       !! CROP->PAST
       landuse_transitions%Test_CROP_2_PAST(kstart:kend) = 0._dp
       landuse_transitions%Test_CROP_2_PAST(kstart:kend) = - DailyTransMtrx(1:nidx,PAST,CROP)  * cf0_CROP(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(is_Crop(n,i) .AND. is_Pasture(n,j)) then
                   landuse_transitions%Test_CROP_2_PAST(kstart-1+n) = &
                        landuse_transitions%Test_CROP_2_PAST(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do
       !! Set Testarray to zero where test is not applicable
       WHERE(cf0_PAST(1:nidx) < fract_small) landuse_transitions%Test_CROP_2_PAST(kstart:kend)=0._dp

       !! PAST->CROP
       landuse_transitions%Test_PAST_2_CROP(kstart:kend) = 0._dp
       landuse_transitions%Test_PAST_2_CROP(kstart:kend) = - DailyTransMtrx(1:nidx,CROP,PAST)  * cf0_PAST(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(is_Pasture(n,i)  .AND. is_Crop(n,j)) then
                   landuse_transitions%Test_PAST_2_CROP(kstart-1+n) = &
                        landuse_transitions%Test_PAST_2_CROP(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do
       !! Set Testarray to zero where test is not applicable
       WHERE(cf0_CROP(1:nidx) < fract_small) landuse_transitions%Test_PAST_2_CROP(kstart:kend)=0._dp

       !! FRST->PAST
       landuse_transitions%Test_FRST_2_PAST(kstart:kend) = 0._dp
       landuse_transitions%Test_FRST_2_PAST(kstart:kend) = - DailyTransMtrx(1:nidx,PAST,FRST) * cf0_FRST(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(( is_NaturalVeg(n,i) .AND. .NOT. is_Grass(n,i) ) .AND. is_pasture(n,j)) then
                   landuse_transitions%Test_FRST_2_PAST(kstart-1+n) = &
                        landuse_transitions%Test_FRST_2_PAST(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do
       !! Set Testarray to zero where test is not applicable
       WHERE(cf0_PAST(1:nidx) < fract_small) landuse_transitions%Test_FRST_2_PAST(kstart:kend)=0._dp

       !! PAST->FRST
       landuse_transitions%Test_PAST_2_FRST(kstart:kend) = 0._dp
       landuse_transitions%Test_PAST_2_FRST(kstart:kend) = - DailyTransMtrx(1:nidx,FRST,PAST) * cf0_PAST(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(is_pasture(n,i) .AND. ( is_NaturalVeg(n,j) .AND. .NOT. is_Grass(n,j) )) then
                   landuse_transitions%Test_PAST_2_FRST(kstart-1+n) = &
                        landuse_transitions%Test_PAST_2_FRST(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do
       !! Set Testarray to zero where test is not applicable
       WHERE(cf_pot_FRST(1:nidx) - cf0_FRST(1:nidx) < fract_small) landuse_transitions%Test_PAST_2_FRST(kstart:kend)=0._dp

       !! GRAS->PAST

       landuse_transitions%Test_GRAS_2_PAST(kstart:kend) = 0._dp
       landuse_transitions%Test_GRAS_2_PAST(kstart:kend) = - DailyTransMtrx(1:nidx,PAST,GRAS) * cf0_GRAS(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(is_Grass(n,i) .AND. is_pasture(n,j)) then
                   landuse_transitions%Test_GRAS_2_PAST(kstart-1+n) = &
                        landuse_transitions%Test_GRAS_2_PAST(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do

       !! PAST->GRAS
       landuse_transitions%Test_PAST_2_GRAS(kstart:kend) = 0._dp
       landuse_transitions%Test_PAST_2_GRAS(kstart:kend) = - DailyTransMtrx(1:nidx,GRAS,PAST) * cf0_PAST(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(is_Pasture(n,i) .AND. is_Grass(n,j) ) then
                   landuse_transitions%Test_PAST_2_GRAS(kstart-1+n) = &
                        landuse_transitions%Test_PAST_2_GRAS(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do

       !! FRST->CROP
       landuse_transitions%Test_FRST_2_CROP(kstart:kend) = 0._dp
       landuse_transitions%Test_FRST_2_CROP(kstart:kend) = - DailyTransMtrx(1:nidx,CROP,FRST) * cf0_FRST(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(( is_NaturalVeg(n,i) .AND. .NOT. is_Grass(n,i) ) .AND. is_Crop(n,j)) then
                   landuse_transitions%Test_FRST_2_CROP(kstart-1+n) = &
                        landuse_transitions%Test_FRST_2_CROP(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do
       !! Set Testarray to zero where test is not applicable
       WHERE(cf0_CROP(1:nidx) < fract_small) landuse_transitions%Test_FRST_2_CROP(kstart:kend)=0._dp

       !! CROP->FRST
       landuse_transitions%Test_CROP_2_FRST(kstart:kend) = 0._dp
       landuse_transitions%Test_CROP_2_FRST(kstart:kend) = - DailyTransMtrx(1:nidx,FRST,CROP) * cf0_CROP(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(is_Crop(n,i) .AND. ( is_NaturalVeg(n,j) .AND. .NOT. is_Grass(n,j) )) then
                   landuse_transitions%Test_CROP_2_FRST(kstart-1+n) = &
                        landuse_transitions%Test_CROP_2_FRST(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do
       WHERE(cf_pot_FRST(1:nidx) - cf0_FRST(1:nidx) < fract_small) landuse_transitions%Test_CROP_2_FRST(kstart:kend) = 0._dp

       !! GRAS->CROP
       landuse_transitions%Test_GRAS_2_CROP(kstart:kend) = 0._dp
       landuse_transitions%Test_GRAS_2_CROP(kstart:kend) = - DailyTransMtrx(1:nidx,CROP,GRAS) * cf0_GRAS(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(is_Grass(n,i) .AND. is_Crop(n,j)) then
                   landuse_transitions%Test_GRAS_2_CROP(kstart-1+n) = &
                        landuse_transitions%Test_GRAS_2_CROP(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do
       WHERE(cf0_CROP(1:nidx) < fract_small) landuse_transitions%Test_GRAS_2_CROP(kstart:kend) = 0._dp

       !! CROP->GRAS
       landuse_transitions%Test_CROP_2_GRAS(kstart:kend) = 0._dp
       landuse_transitions%Test_CROP_2_GRAS(kstart:kend) = - DailyTransMtrx(1:nidx,GRAS,CROP) * cf0_CROP(1:nidx)
       do n=1,nidx
          do i=1,ntiles
             do j=1,ntiles
                if(is_Crop(n,i) .AND. is_Grass(n,j) ) then
                   landuse_transitions%Test_CROP_2_GRAS(kstart-1+n) = &
                        landuse_transitions%Test_CROP_2_GRAS(kstart-1+n) &
                        + Tile_TransMtrx(n,j,i) * landcover_fract_current(n,i)
                end if
             end do
          end do
       end do
       WHERE(cf_pot_GRAS(1:nidx) - cf0_GRAS(1:nidx) < fract_small) landuse_transitions%Test_CROP_2_GRAS(kstart:kend) = 0._dp

    END IF

    !!
    !! --- Update daily cover fractions ----------------------------------------
    !!

    DO n=1,nidx !! Go through all gridboxes

       !! Update cover fractions (this works also for glacier points because there the matrix is a unit matrix)

       landcover_fract_new(n,1:ntiles) = &
            MATMUL( Tile_TransMtrx(n,1:ntiles,1:ntiles),landcover_fract_current(n,1:ntiles) )

       IF(ANY(is_glacier(n,1:ntiles))) CYCLE !! nothing to do for glacier points

       !! Iif necessary modify Tile_TransMtrx() such that landcover_fract remains larger than fract_small)

       opt=0
       DO !! It needs maximally ntiles steps of modification because then Tile_TransMtrx(n,:,:) is the unit matrix
          opt=opt+1 !! Count the number of modification steps

          IF(opt>ntiles) THEN !! Emergency exit: In this case something went wrong (e.g. wrong initial cover fractions)
             WRITE (message_text, *) ' After opt=',opt,' modification steps: landcover_fract_new(',&
                  n,',:) either still smaller than fract_small or larger than fract_pot for a natural vegetation tile' 
             CALL finish('do_landuse_transitions', message_text)
          END IF

          IF( .NOT. ANY( landcover_fract_new(n,1:ntiles) < fract_small .AND. is_vegetation(n,1:ntiles) ) &
              .AND. &
              .NOT. ANY( (landcover_fract_new(n,1:ntiles)-landcover_fract_pot(n,1:ntiles))>10._dp*epsilon(1.0_dp) &
                         .AND. is_naturalVeg(n,1:ntiles)) ) EXIT !! exit opt-loop

          DO i=1,ntiles
             IF(  landcover_fract_new(n,i) < fract_small                                                &
                  .OR. ( landcover_fract_new(n,i) > landcover_fract_pot(n,i) .AND. is_naturalVeg(n,i) ) &
                ) THEN !! exclude i-th tile from landcover dynamics
                Tile_TransMtrx(n,i,1:ntiles)=0._dp
                Tile_TransMtrx(n,1:ntiles,i)=0._dp
             END IF
          END DO

          !! Recompute diagonal elements of tile transition matrix

          DO i=1,ntiles
             Tile_TransMtrx(n,i,i)= 0.0_dp !! this is to exclude summation of diagonal elements below
             Tile_TransMtrx(n,i,i)= MAX( 0._dp, 1.0_dp - SUM(Tile_TransMtrx(n,1:ntiles,i)) )
          END DO

          !! Assure Matrixelements to be in [0,1] (remove possible numerical inaccuracies in above computations)
          Tile_TransMtrx(n,1:ntiles,1:ntiles) = MAX(0.0_dp,MIN(1.0_dp,Tile_TransMtrx(n,1:ntiles,1:ntiles)))

          !! Recompute new landcover fractions

          landcover_fract_new(n,1:ntiles) = &
               MATMUL( Tile_TransMtrx(n,1:ntiles,1:ntiles),landcover_fract_current(n,1:ntiles) )

       END DO !! opt-loop

    END DO !! n-loop

    !! Correct for possible numerical errors so that natural vegetation cover remains smaller than potential cover

    WHERE (is_naturalVeg(1:nidx,1:ntiles) .AND. &
           landcover_fract_new(1:nidx,1:ntiles) > landcover_fract_pot(1:nidx,1:ntiles) ) 
       landcover_fract_new(1:nidx,1:ntiles) = landcover_fract_pot(1:nidx,1:ntiles)
    END WHERE

    !! Do some checks


    DO n=1,nidx
       IF (ANY(landcover_fract_new(n,1:ntiles) < fract_small .AND. is_vegetation(n,1:ntiles))) THEN
             WRITE (message_text, *) ' landcover_fract_new(',n,',',MINLOC(landcover_fract_new(n,1:ntiles)),')=',&
                  MINVAL(landcover_fract_new(n,1:ntiles)),' is smaller than fract_small' 
          CALL finish ('do_landuse_transitions', message_text)
       END IF
    END DO

    IF (debug_Cconservation) THEN

       !! --- Initialize Carbon conservation test
       IF (.NOT. with_yasso) THEN
          lcc%LCC_testCconserv(kstart:kend) = veg_ratio_max(1:nidx)         &
            * SUM(  landcover_fract_current(1:nidx,1:ntiles)                &
                  * veg_fract_correction   (1:nidx,1:ntiles)                &
                  * (  cbalance%Cpool_green          (kstart:kend,1:ntiles) &
                     + cbalance%Cpool_woods          (kstart:kend,1:ntiles) &
                     + cbalance%Cpool_reserve        (kstart:kend,1:ntiles) &
                     + cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles) &
                     + cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles) &
                     + cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles) &
                     + cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles) &
                     + cbalance%Cpool_slow           (kstart:kend,1:ntiles) &
                    ), DIM=2 )
       ELSE
          lcc%LCC_testCconserv(kstart:kend) = veg_ratio_max(1:nidx)          &
            * SUM(  landcover_fract_current(1:nidx,1:ntiles)                 &
                  * veg_fract_correction   (1:nidx,1:ntiles)                 &
                  * (  cbalance%Cpool_green           (kstart:kend,1:ntiles) &
                     + cbalance%Cpool_woods           (kstart:kend,1:ntiles) &
                     + cbalance%Cpool_reserve         (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_acid_ag1       (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_acid_bg1       (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_water_ag1      (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_water_bg1      (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_ethanol_ag1    (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_ethanol_bg1    (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_nonsoluble_ag1 (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_nonsoluble_bg1 (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_humus_1        (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_acid_ag2       (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_acid_bg2       (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_water_ag2      (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_water_bg2      (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_ethanol_ag2    (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_ethanol_bg2    (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_nonsoluble_ag2 (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_nonsoluble_bg2 (kstart:kend,1:ntiles) &
                     + cbalance%YCpool_humus_2        (kstart:kend,1:ntiles) &
                    ), DIM=2 )
       END IF
       IF (lcc_scheme==2) THEN 
          lcc%LCC_testCconserv(kstart:kend) =                &
          lcc%LCC_testCconserv(kstart:kend)                  &
               + (  cbalance%Cpool_onSite              (kstart:kend)  &
                  + cbalance%Cpool_paper               (kstart:kend)  &
                  + cbalance%Cpool_construction        (kstart:kend)  &
                  + cbalance%Cpool_paper_harvest       (kstart:kend)  &
                  + cbalance%Cpool_construction_harvest(kstart:kend)  &
                 ) * veg_ratio_max(1:nidx)
       ENDIF

    END IF

    !! === Perform the carbon relocations associated with the landuse transitions
    SELECT CASE (lcc_scheme)
    CASE (1)
       IF (with_nitrogen) THEN
          CALL C_relocation_from_LUtransitions(lctlib, surface, nidx, ntiles, is_vegetation,                  &
                                               landcover_fract_current       (     1:nidx,1:ntiles),          &
                                               landcover_fract_new           (     1:nidx,1:ntiles),          &
                                               Tile_TransMtrx                (     1:nidx,1:ntiles,1:ntiles), &
                                               veg_fract_correction          (     1:nidx,1:ntiles),          &
                                               frac_wood_2_atmos,                                             &
                                               frac_green_2_atmos,                                            &
                                               frac_reserve_2_atmos,                                          &
                                               landuse_transitions%C2atmos_LUtrans(kstart:kend),              &
                                               cbalance%Cpool_green          (kstart:kend,1:ntiles),          &
                                               cbalance%Cpool_woods          (kstart:kend,1:ntiles),          &
                                               cbalance%Cpool_reserve        (kstart:kend,1:ntiles),          &
                     Cpool_litter_green_ag   = cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles),          &
                     Cpool_litter_green_bg   = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles),          &
                     Cpool_litter_wood_ag    = cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles),          &
                     Cpool_litter_wood_bg    = cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles),          &
                     Cpool_slow              = cbalance%Cpool_slow           (kstart:kend,1:ntiles),          &
                     C_2_litterGreenPools    = C2litterGreenPools            (     1:nidx),                   &
                     C_2_litterWoodPool_ag   = C2litterWoodPool_ag           (     1:nidx),                   &
                     C_2_litterWoodPool_bg   = C2litterWoodPool_bg           (     1:nidx),                   &
                     Npool_green             = nbalance%Npool_green(kstart:kend,1:ntiles),                    &
                     Npool_woods             = nbalance%Npool_woods(kstart:kend,1:ntiles),                    &
                     Npool_mobile            = nbalance%Npool_mobile(kstart:kend,1:ntiles),                   &
                     Npool_litter_green_ag   = nbalance%Npool_litter_green_ag(kstart:kend,1:ntiles),          &
                     Npool_litter_green_bg   = nbalance%Npool_litter_green_bg(kstart:kend,1:ntiles),          &
                     Npool_litter_wood_ag    = nbalance%Npool_litter_wood_ag(kstart:kend,1:ntiles),           &
                     Npool_litter_wood_bg    = nbalance%Npool_litter_wood_bg(kstart:kend,1:ntiles),           &
                     Npool_slow              = nbalance%Npool_slow(kstart:kend,1:ntiles),                     &
                     SMINN_pool              = nbalance%SMINN_pool(kstart:kend,1:ntiles),                     & 
                     Nitrogen_2_atmos        = landuse_transitions%N2atmos_LUtrans(kstart:kend),              &
                     frac_mobile_2_atmos     = frac_mobile_2_atmos,                                           &
                     lcc_scheme              = lcc_scheme)
       ELSE
          IF (.NOT. with_yasso) THEN
             CALL C_relocation_from_LUtransitions(lctlib, surface, nidx, ntiles, is_vegetation,               &
                                               landcover_fract_current       (     1:nidx,1:ntiles),          &
                                               landcover_fract_new           (     1:nidx,1:ntiles),          &
                                               Tile_TransMtrx                (     1:nidx,1:ntiles,1:ntiles), &
                                               veg_fract_correction          (     1:nidx,1:ntiles),          &
                                               frac_wood_2_atmos,                                             &
                                               frac_green_2_atmos,                                            &
                                               frac_reserve_2_atmos,                                          &
                                               landuse_transitions%C2atmos_LUtrans(kstart:kend),              &
                                               cbalance%Cpool_green          (kstart:kend,1:ntiles),          &
                                               cbalance%Cpool_woods          (kstart:kend,1:ntiles),          &
                                               cbalance%Cpool_reserve        (kstart:kend,1:ntiles),          &
                     Cpool_litter_green_ag   = cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles),          &
                     Cpool_litter_green_bg   = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles),          &
                     Cpool_litter_wood_ag    = cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles),          &
                     Cpool_litter_wood_bg    = cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles),          &
                     Cpool_slow              = cbalance%Cpool_slow           (kstart:kend,1:ntiles),          &
                     C_2_litterGreenPools    = C2litterGreenPools            (     1:nidx),                   &
                     C_2_litterWoodPool_ag   = C2litterWoodPool_ag           (     1:nidx),                   &
                     C_2_litterWoodPool_bg   = C2litterWoodPool_bg           (     1:nidx),                   &
                     lcc_scheme              = lcc_scheme)
          ELSE  !! with yasso
             CALL C_relocation_from_LUtransitions(lctlib,surface,nidx,ntiles,is_vegetation,                    &
                                               landcover_fract_current       (     1:nidx,1:ntiles),           &
                                               landcover_fract_new           (     1:nidx,1:ntiles),           &
                                               Tile_TransMtrx                (     1:nidx,1:ntiles,1:ntiles),  &
                                               veg_fract_correction          (     1:nidx,1:ntiles),           &
                                               frac_wood_2_atmos,                                              &
                                               frac_green_2_atmos,                                             &
                                               frac_reserve_2_atmos,                                           &
                                               landuse_transitions%C2atmos_LUtrans(kstart:kend),               &
                                               cbalance%Cpool_green          (kstart:kend,1:ntiles),           &
                                               cbalance%Cpool_woods          (kstart:kend,1:ntiles),           &
                                               cbalance%Cpool_reserve        (kstart:kend,1:ntiles),           &
                     ! Yasso pools & parameters                          
                     YCpool_acid_ag1         = cbalance%YCpool_acid_ag1      (kstart:kend,1:ntiles),           &
                     YCpool_water_ag1        = cbalance%YCpool_water_ag1     (kstart:kend,1:ntiles),           &
                     YCpool_ethanol_ag1      = cbalance%YCpool_ethanol_ag1   (kstart:kend,1:ntiles),           &
                     YCpool_nonsoluble_ag1   = cbalance%YCpool_nonsoluble_ag1(kstart:kend,1:ntiles),           &
                     YCpool_acid_bg1         = cbalance%YCpool_acid_bg1      (kstart:kend,1:ntiles),           &
                     YCpool_water_bg1        = cbalance%YCpool_water_bg1     (kstart:kend,1:ntiles),           &
                     YCpool_ethanol_bg1      = cbalance%YCpool_ethanol_bg1   (kstart:kend,1:ntiles),           &
                     YCpool_nonsoluble_bg1   = cbalance%YCpool_nonsoluble_bg1(kstart:kend,1:ntiles),           &
                     YCpool_humus_1          = cbalance%YCpool_humus_1       (kstart:kend,1:ntiles),           &
                     YCpool_acid_ag2         = cbalance%YCpool_acid_ag2      (kstart:kend,1:ntiles),           &
                     YCpool_water_ag2        = cbalance%YCpool_water_ag2     (kstart:kend,1:ntiles),           &
                     YCpool_ethanol_ag2      = cbalance%YCpool_ethanol_ag2   (kstart:kend,1:ntiles),           &
                     YCpool_nonsoluble_ag2   = cbalance%YCpool_nonsoluble_ag2(kstart:kend,1:ntiles),           &
                     YCpool_acid_bg2         = cbalance%YCpool_acid_bg2      (kstart:kend,1:ntiles),           &
                     YCpool_water_bg2        = cbalance%YCpool_water_bg2     (kstart:kend,1:ntiles),           &
                     YCpool_ethanol_bg2      = cbalance%YCpool_ethanol_bg2   (kstart:kend,1:ntiles),           &
                     YCpool_nonsoluble_bg2   = cbalance%YCpool_nonsoluble_bg2(kstart:kend,1:ntiles),           &
                     YCpool_humus_2          = cbalance%YCpool_humus_2       (kstart:kend,1:ntiles),           &
                     LeafLitcoef             = LeafLit_coef                  (     1:nidx,1:ntiles,1:5),       &
                     WoodLitcoef             = WoodLit_coef                  (     1:nidx,1:ntiles,1:5),       &
                     !
                     C_2_litterGreenPools    = C2litterGreenPools            (     1:nidx),                    &
                     C_2_litterWoodPool_ag   = C2litterWoodPool_ag           (     1:nidx),                    &
                     C_2_litterWoodPool_bg   = C2litterWoodPool_bg           (     1:nidx),                    &
                     lcc_scheme              = lcc_scheme)
          END IF ! yasso
       END IF ! nitrogen
    CASE (2)
       IF (.NOT. with_yasso) THEN
          CALL C_relocation_from_LUtransitions(lctlib,surface,nidx,ntiles,is_vegetation,&
                         landcover_fract_current       (     1:nidx,1:ntiles),          &
                         landcover_fract_new           (     1:nidx,1:ntiles),          &
                         Tile_TransMtrx                (     1:nidx,1:ntiles,1:ntiles), &
                         veg_fract_correction          (     1:nidx,1:ntiles),          &
                         frac_wood_2_atmos,                                             &
                         frac_green_2_atmos,                                            &
                         frac_reserve_2_atmos,                                          &
                         landuse_transitions%C2atmos_LUtrans(kstart:kend),              &
                         cbalance%Cpool_green          (kstart:kend,1:ntiles),          &
                         cbalance%Cpool_woods          (kstart:kend,1:ntiles),          &
                         cbalance%Cpool_reserve        (kstart:kend,1:ntiles),          &
                         Cpool_litter_green_ag   = cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles),   &
                         Cpool_litter_green_bg   = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles),   &
                         Cpool_litter_wood_ag    = cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles),   &
                         Cpool_litter_wood_bg    = cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles),   &
                         Cpool_slow              = cbalance%Cpool_slow           (kstart:kend,1:ntiles),   &
                         Cpool_onSite             = cbalance%Cpool_onSite        (kstart:kend), &
                         Cpool_paper              = cbalance%Cpool_paper         (kstart:kend), &
                         Cpool_construction       = cbalance%Cpool_construction  (kstart:kend), &
                         C_onSite_2_atmos         = lcc%C_onSite_2_atmos         (kstart:kend), &
                         C_paper_2_atmos          = lcc%C_paper_2_atmos          (kstart:kend), &
                         C_construction_2_atmos   = lcc%C_construction_2_atmos   (kstart:kend), &                         
                         C_2_onSite               = lcc%C_2_onSite               (kstart:kend), &
                         C_2_paper                = lcc%C_2_paper                (kstart:kend), &
                         C_2_construction         = lcc%C_2_construction         (kstart:kend), &
                         lcc_scheme               = lcc_scheme)
       ELSE  !! with yasso
          CALL C_relocation_from_LUtransitions(lctlib,surface,nidx,ntiles,is_vegetation,&
                         landcover_fract_current       (     1:nidx,1:ntiles),          &
                         landcover_fract_new           (     1:nidx,1:ntiles),          &
                         Tile_TransMtrx                (     1:nidx,1:ntiles,1:ntiles), &
                         veg_fract_correction          (     1:nidx,1:ntiles),          &
                         frac_wood_2_atmos,                                             &
                         frac_green_2_atmos,                                            &
                         frac_reserve_2_atmos,                                          &
                         landuse_transitions%C2atmos_LUtrans(kstart:kend),              &
                         cbalance%Cpool_green          (kstart:kend,1:ntiles),          &
                         cbalance%Cpool_woods          (kstart:kend,1:ntiles),          &
                         cbalance%Cpool_reserve        (kstart:kend,1:ntiles),          &
                         YCpool_acid_ag1         = cbalance%YCpool_acid_ag1      (kstart:kend,1:ntiles),  &
                         YCpool_water_ag1        = cbalance%YCpool_water_ag1     (kstart:kend,1:ntiles),  &
                         YCpool_ethanol_ag1      = cbalance%YCpool_ethanol_ag1   (kstart:kend,1:ntiles),  &
                         YCpool_nonsoluble_ag1   = cbalance%YCpool_nonsoluble_ag1(kstart:kend,1:ntiles),  &
                         YCpool_acid_bg1         = cbalance%YCpool_acid_bg1      (kstart:kend,1:ntiles),  &
                         YCpool_water_bg1        = cbalance%YCpool_water_bg1     (kstart:kend,1:ntiles),  &
                         YCpool_ethanol_bg1      = cbalance%YCpool_ethanol_bg1   (kstart:kend,1:ntiles),  &
                         YCpool_nonsoluble_bg1   = cbalance%YCpool_nonsoluble_bg1(kstart:kend,1:ntiles),  &
                         YCpool_humus_1          = cbalance%YCpool_humus_1       (kstart:kend,1:ntiles),  &
                         YCpool_acid_ag2         = cbalance%YCpool_acid_ag2      (kstart:kend,1:ntiles),  &
                         YCpool_water_ag2        = cbalance%YCpool_water_ag2     (kstart:kend,1:ntiles),  &
                         YCpool_ethanol_ag2      = cbalance%YCpool_ethanol_ag2   (kstart:kend,1:ntiles),  &
                         YCpool_nonsoluble_ag2   = cbalance%YCpool_nonsoluble_ag2(kstart:kend,1:ntiles),  &
                         YCpool_acid_bg2         = cbalance%YCpool_acid_bg2      (kstart:kend,1:ntiles),  &
                         YCpool_water_bg2        = cbalance%YCpool_water_bg2     (kstart:kend,1:ntiles),  &
                         YCpool_ethanol_bg2      = cbalance%YCpool_ethanol_bg2   (kstart:kend,1:ntiles),  &
                         YCpool_nonsoluble_bg2   = cbalance%YCpool_nonsoluble_bg2(kstart:kend,1:ntiles),  &
                         YCpool_humus_2          = cbalance%YCpool_humus_2       (kstart:kend,1:ntiles),  &
                         LeafLitcoef             = LeafLit_coef                 (kstart:kend,1:ntiles,1:5), &
                         WoodLitcoef             = WoodLit_coef                 (kstart:kend,1:ntiles,1:5), &
                         Cpool_onSite            = cbalance%Cpool_onSite        (kstart:kend),              &
                         Cpool_paper             = cbalance%Cpool_paper         (kstart:kend),              &
                         Cpool_construction      = cbalance%Cpool_construction  (kstart:kend),              &
                         C_onSite_2_atmos        = lcc%C_onSite_2_atmos         (kstart:kend),              &
                         C_paper_2_atmos         = lcc%C_paper_2_atmos          (kstart:kend),              &
                         C_construction_2_atmos  = lcc%C_construction_2_atmos   (kstart:kend),              &
                         C_2_onSite              = lcc%C_2_onSite               (kstart:kend),              &
                         C_2_paper               = lcc%C_2_paper                (kstart:kend),              &
                         C_2_construction        = lcc%C_2_construction         (kstart:kend),              &
                         lcc_scheme              = lcc_scheme)
       END IF
    END SELECT

    !!
    !! Do harvest
    !!

    !! compute harvest per day relative to vegetated area
    rhlp(1:nidx) = landuse_transitions%Box_harvest(kstart:kend) * sec_per_day &
                     / MAX(veg_ratio_max(1:nidx),fract_small)

    IF (lcc_scheme==1) THEN

       IF (with_nitrogen) THEN

          CALL C_relocation_from_harvest(nidx, ntiles, lcc_scheme, lctlib,surface,                     &
                         landcover_fract_new(1:nidx,1:ntiles),                                         &
                         veg_fract_correction(1:nidx,1:ntiles),                                        &
                         is_naturalVeg(1:nidx,1:ntiles),                                               &
                         cbalance%Cpool_green(kstart:kend,1:ntiles),                                   &
                         cbalance%Cpool_woods(kstart:kend,1:ntiles),                                   &
                         cbalance%Cpool_reserve(kstart:kend,1:ntiles),                                 &
                         rhlp(1:nidx),                                                                 &
                         C2litter_harvest(1:nidx),                                                     &
                         landuse_transitions%C2atmos_harvest(kstart:kend),                             &
                         frac_harvest_2_atmos  = frac_harvest_2_atmos,                                 &
                         Cpool_litter_green_bg = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles), &
                         Cpool_litter_wood_bg  = cbalance%Cpool_litter_wood_bg( kstart:kend,1:ntiles), &
                         Npool_green           = nbalance%Npool_green(kstart:kend,1:ntiles),           &
                         Npool_woods           = nbalance%Npool_woods(kstart:kend,1:ntiles),           &
                         Npool_mobile          = nbalance%Npool_mobile(kstart:kend,1:ntiles),          &
                         Npool_litter_green_bg = nbalance%Npool_litter_green_bg(kstart:kend,1:ntiles), &
                         Npool_litter_wood_bg  = nbalance%Npool_litter_wood_bg(kstart:kend,1:ntiles),  &
                         SMINN_pool            = nbalance%SMINN_pool(kstart:kend,1:ntiles),            &
                         N2atmos               = landuse_transitions%N2atmos_harvest(kstart:kend),     &
                         frac_mobile_2_atmos   = frac_mobile_2_atmos )
       ELSE
          IF (.NOT. with_yasso) THEN
             CALL C_relocation_from_harvest(nidx,ntiles,lcc_scheme,lctlib,surface,  &
                         landcover_fract_new(1:nidx,1:ntiles),                      &
                         veg_fract_correction(1:nidx,1:ntiles),                     &
                         is_naturalVeg(1:nidx,1:ntiles),                            &
                         cbalance%Cpool_green(kstart:kend,1:ntiles),                &
                         cbalance%Cpool_woods(kstart:kend,1:ntiles),                &
                         cbalance%Cpool_reserve(kstart:kend,1:ntiles),              &
                         rhlp                               (     1:nidx),          &
                         C2litter_harvest                   (     1:nidx),          &
                         landuse_transitions%C2atmos_harvest(kstart:kend),          &
                         Cpool_litter_green_bg = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles), &
                         Cpool_litter_wood_bg  = cbalance%Cpool_litter_wood_bg( kstart:kend,1:ntiles), &
                         frac_harvest_2_atmos  = frac_harvest_2_atmos               &
                        )
          ELSE  !! with yasso
             CALL C_relocation_from_harvest(nidx,ntiles,lcc_scheme,lctlib,surface,                       &  
                                            landcover_fract_new                (1:nidx,     1:ntiles),   &
                                            veg_fract_correction               (1:nidx,     1:ntiles),   &
                                            is_naturalVeg                      (1:nidx,     1:ntiles),   &
                                            cbalance%Cpool_green               (kstart:kend,1:ntiles),   &
                                            cbalance%Cpool_woods               (kstart:kend,1:ntiles),   &
                                            cbalance%Cpool_reserve             (kstart:kend,1:ntiles),   &
                                            rhlp                               (     1:nidx),            &
                                            C2litter_harvest                   (     1:nidx),            &
                                            landuse_transitions%C2atmos_harvest(kstart:kend),            &
                    YCpool_acid_bg1       = cbalance%YCpool_acid_bg1           (kstart:kend,1:ntiles),   &
                    YCpool_water_bg1      = cbalance%YCpool_water_bg1          (kstart:kend,1:ntiles),   &
                    YCpool_ethanol_bg1    = cbalance%YCpool_ethanol_bg1        (kstart:kend,1:ntiles),   & 
                    YCpool_nonsoluble_bg1 = cbalance%YCpool_nonsoluble_bg1     (kstart:kend,1:ntiles),   &
                    YCpool_humus_1        = cbalance%YCpool_humus_1            (kstart:kend,1:ntiles),   &
                    YCpool_acid_bg2       = cbalance%YCpool_acid_bg2           (kstart:kend,1:ntiles),   &
                    YCpool_water_bg2      = cbalance%YCpool_water_bg2          (kstart:kend,1:ntiles),   &
                    YCpool_ethanol_bg2    = cbalance%YCpool_ethanol_bg2        (kstart:kend,1:ntiles),   & 
                    YCpool_nonsoluble_bg2 = cbalance%YCpool_nonsoluble_bg2     (kstart:kend,1:ntiles),   &
                    YCpool_humus_2        = cbalance%YCpool_humus_2            (kstart:kend,1:ntiles),   &
                    LeafLit_coef          = LeafLit_coef                       (     1:nidx,1:ntiles,:), &
                    WoodLit_coef          = WoodLit_coef                       (     1:nidx,1:ntiles,:), &
                    frac_harvest_2_atmos  = frac_harvest_2_atmos                                         &
                    )
          END IF  ! yasso
       END IF  ! with_nitrogen
    ELSE ! lcc_scheme==2
       IF (.NOT. with_yasso) THEN
          CALL C_relocation_from_harvest(nidx,ntiles,lcc_scheme,lctlib,surface,     &
                         landcover_fract_new                (     1:nidx,1:ntiles), &
                         veg_fract_correction               (     1:nidx,1:ntiles), &
                         is_naturalVeg                      (     1:nidx,1:ntiles), &
                         cbalance%Cpool_green               (kstart:kend,1:ntiles), &
                         cbalance%Cpool_woods               (kstart:kend,1:ntiles), &
                         cbalance%Cpool_reserve             (kstart:kend,1:ntiles), &
                         rhlp                               (     1:nidx),          &
                         C2litter_harvest                   (     1:nidx),          &
                         landuse_transitions%C2atmos_harvest(kstart:kend),          &
                         Cpool_litter_green_bg   = cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles), &
                         Cpool_litter_wood_bg    = cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles), &
                         Cpool_paper                 = cbalance%Cpool_paper_harvest       (kstart:kend), &
                         Cpool_construction          = cbalance%Cpool_construction_harvest(kstart:kend), &
                         Carbon_paper_2_atmos        = lcc%C_paper_harv_2_atmos           (kstart:kend), &
                         Carbon_construction_2_atmos = lcc%C_construction_harv_2_atmos    (kstart:kend), &
                         Carbon_2_paper              = lcc%C_2_paper_harv                 (kstart:kend), &
                         Carbon_2_construction       = lcc%C_2_construction_harv          (kstart:kend)  &
                        )
       ELSE    ! YASSO               
          CALL C_relocation_from_harvest(nidx,ntiles,lcc_scheme,lctlib,surface,     &
                         landcover_fract_new                (     1:nidx,1:ntiles), &
                         veg_fract_correction               (     1:nidx,1:ntiles), &
                         is_naturalVeg                      (     1:nidx,1:ntiles), &
                         cbalance%Cpool_green               (kstart:kend,1:ntiles), &
                         cbalance%Cpool_woods               (kstart:kend,1:ntiles), &
                         cbalance%Cpool_reserve             (kstart:kend,1:ntiles), &
                         rhlp                               (     1:nidx),          &
                         C2litter_harvest                   (     1:nidx),          &
                         landuse_transitions%C2atmos_harvest(kstart:kend),          &
                         YCpool_acid_bg1       = cbalance%YCpool_acid_bg1           (kstart:kend,1:ntiles),   &
                         YCpool_water_bg1      = cbalance%YCpool_water_bg1          (kstart:kend,1:ntiles),   &
                         YCpool_ethanol_bg1    = cbalance%YCpool_ethanol_bg1        (kstart:kend,1:ntiles),   & 
                         YCpool_nonsoluble_bg1 = cbalance%YCpool_nonsoluble_bg1     (kstart:kend,1:ntiles),   &
                         YCpool_humus_1        = cbalance%YCpool_humus_1            (kstart:kend,1:ntiles),   &
                         YCpool_acid_bg2       = cbalance%YCpool_acid_bg2           (kstart:kend,1:ntiles),   &
                         YCpool_water_bg2      = cbalance%YCpool_water_bg2          (kstart:kend,1:ntiles),   &
                         YCpool_ethanol_bg2    = cbalance%YCpool_ethanol_bg2        (kstart:kend,1:ntiles),   & 
                         YCpool_nonsoluble_bg2 = cbalance%YCpool_nonsoluble_bg2     (kstart:kend,1:ntiles),   &
                         YCpool_humus_2        = cbalance%YCpool_humus_2            (kstart:kend,1:ntiles),   &
                         LeafLit_coef          = LeafLit_coef                       (     1:nidx,1:ntiles,:), &
                         WoodLit_coef          = WoodLit_coef                       (     1:nidx,1:ntiles,:), &
                         Cpool_paper                 = cbalance%Cpool_paper_harvest       (kstart:kend),      &
                         Cpool_construction          = cbalance%Cpool_construction_harvest(kstart:kend),      &
                         Carbon_paper_2_atmos        = lcc%C_paper_harv_2_atmos           (kstart:kend),      &
                         Carbon_construction_2_atmos = lcc%C_construction_harv_2_atmos    (kstart:kend),      &
                         Carbon_2_paper              = lcc%C_2_paper_harv                 (kstart:kend),      &
                         Carbon_2_construction       = lcc%C_2_construction_harv          (kstart:kend)       &
                        )
       END IF  ! .NOT.with_YASSO   
    ENDIF

    IF (debug_Cconservation) THEN
       !! --- Finish carbon conservation test
       IF (.NOT. with_yasso) THEN
          lcc%LCC_testCconserv(kstart:kend) = lcc%LCC_testCconserv(kstart:kend)   &
            - (  landuse_transitions%C2atmos_LUtrans(kstart:kend)              &
               + landuse_transitions%C2atmos_harvest(kstart:kend)              &
               + SUM(  landcover_fract_new              (     1:nidx,1:ntiles) &
                     * veg_fract_correction             (     1:nidx,1:ntiles) &
                     * (  cbalance%Cpool_green          (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_woods          (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_reserve        (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_litter_green_ag(kstart:kend,1:ntiles) &
                        + cbalance%Cpool_litter_green_bg(kstart:kend,1:ntiles) &
                        + cbalance%Cpool_litter_wood_ag (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_litter_wood_bg (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_slow           (kstart:kend,1:ntiles) &
                       ), DIM=2 )                                              &
              ) * veg_ratio_max(1:nidx)
       ELSE ! with_yasso
          lcc%LCC_testCconserv(kstart:kend) = lcc%LCC_testCconserv(kstart:kend)   &
            - (  landuse_transitions%C2atmos_LUtrans(kstart:kend)              &
               + landuse_transitions%C2atmos_harvest(kstart:kend)              &
               + SUM(  landcover_fract_new              (     1:nidx,1:ntiles) &
                     * veg_fract_correction             (     1:nidx,1:ntiles) &
                     * (  cbalance%Cpool_green          (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_woods          (kstart:kend,1:ntiles) &
                        + cbalance%Cpool_reserve        (kstart:kend,1:ntiles) &
                        + cbalance%YCpool_acid_ag1      (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_acid_bg1      (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_water_ag1     (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_water_bg1     (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_ethanol_ag1   (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_ethanol_bg1   (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_nonsoluble_ag1(kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_nonsoluble_bg1(kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_humus_1       (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_acid_ag2      (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_acid_bg2      (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_water_ag2     (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_water_bg2     (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_ethanol_ag2   (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_ethanol_bg2   (kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_nonsoluble_ag2(kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_nonsoluble_bg2(kstart:kend,1:ntiles) & 
                        + cbalance%YCpool_humus_2       (kstart:kend,1:ntiles) & 
                       ), DIM=2 )                                              &
               ) * veg_ratio_max(1:nidx)
       END IF
       IF (lcc_scheme==2) THEN
          lcc%LCC_testCconserv(kstart:kend) =                        &
          lcc%LCC_testCconserv(kstart:kend) - veg_ratio_max(1:nidx)  &
               * (  cbalance%Cpool_onSite              (kstart:kend) &
                  + cbalance%Cpool_paper               (kstart:kend) &
                  + cbalance%Cpool_construction        (kstart:kend) &
                  + cbalance%Cpool_paper_harvest       (kstart:kend) &
                  + cbalance%Cpool_construction_harvest(kstart:kend) &
                 )
       ENDIF
    ENDIF

    !! Calculate and cumulate landcover change fluxes (conversion from [mol m-2(veg) day-1] to [mol m-2(grid box) s-1])
    IF (lcc_scheme==2) CALL CumulateAnthroFluxes(lcc, veg_ratio_max)

    IF (lcc_scheme==1) THEN

       lcc%LCC_flux_box_C2litterGreenPools(kstart:kend) = &
       lcc%LCC_flux_box_C2litterGreenPools(kstart:kend) + C2litterGreenPools(1:nidx) * veg_ratio_max(1:nidx)/86400._dp

       lcc%LCC_flux_box_C2litterWoodPool  (kstart:kend) = &
       lcc%LCC_flux_box_C2litterWoodPool  (kstart:kend) +   (C2litterWoodPool_ag(1:nidx)+C2litterWoodPool_bg(1:nidx)) &
                                                          * veg_ratio_max(1:nidx)/86400._dp

       landuse_transitions%Box_flux_harvest(kstart:kend) = &
       landuse_transitions%Box_flux_harvest(kstart:kend) + &
                   (landuse_transitions%C2atmos_harvest(kstart:kend) + C2litter_harvest(1:nidx)) * veg_ratio_max(:)
    ELSE ! lcc_scheme==2
       landuse_transitions%Box_flux_harvest(kstart:kend) = &
       landuse_transitions%Box_flux_harvest(kstart:kend) + veg_ratio_max(:) &
         * (  landuse_transitions%C2atmos_harvest(kstart:kend) + C2litter_harvest(1:nidx) &
            + lcc%C_2_paper_harv(1:nidx) + lcc%C_2_construction_harv(1:nidx) )
    ENDIF

    lcc%LCC_flux_box_C2atmos(kstart:kend) = &
    lcc%LCC_flux_box_C2atmos(kstart:kend) + landuse_transitions%C2atmos_LUtrans(kstart:kend) * veg_ratio_max(1:nidx)/86400._dp

    landuse_transitions%Box_flux_harvest_2atmos(kstart:kend) = landuse_transitions%Box_flux_harvest_2atmos(kstart:kend) &
          + landuse_transitions%C2atmos_harvest(kstart:kend)*veg_ratio_max(:)

    !! Compute fluxes from landcover transitions: 
    !! .. i.e. conversion from [mol(C)/m^2(vegetated area)] during whole day to [kg(CO2)/m^2(grid box) s]
    CO2_emission_LUtrans(:) = landuse_transitions%C2atmos_LUtrans(kstart:kend) * veg_ratio_max(:) * molarMassCO2_kg / 86400._dp
    CO2_emission_harvest(:) = landuse_transitions%C2atmos_harvest(kstart:kend) * veg_ratio_max(:) * molarMassCO2_kg / 86400._dp

    !! Return new cover fractions
    landcover_fract_current(:,:) = landcover_fract_new(:,:)

  END SUBROUTINE do_landuse_transitions


  !! --- read_landuse_transitions() ----------------------------------------------------
  !!
  !! This routine reads in the landuse transitions of the (reduced) New Hampshire Harmonized Landuse Protocol 
  !! for a new year from an external file. The transitions are stored in landuse_transitions%NHHLP_TransMtrx().

  subroutine  read_landuse_transitions(current_year,grid,domain)
    use mo_jsbach_grid,      only: grid_type,domain_type
    use mo_exception,        only: finish, message
    use mo_netcdf,           only: FILE_INFO,IO_inq_dimid,IO_inq_dimlen,IO_inq_varid,io_get_var_double
    USE mo_tr_scatter,       ONLY: scatter_field
    use mo_io,               only: IO_open, IO_READ, IO_close

    INTEGER,           intent(in)   :: current_year
    type(grid_type),   intent(in)   :: grid
    type(domain_type), intent(in)   :: domain

    !! --- locals

    integer                     :: IO_file_id, IO_var_id, IO_dim_id
    type(FILE_INFO)             :: IO_file
    integer                     :: znlon, znlat
    integer                     :: status
    character(len=64)           :: varName 
    character(len=1024)         :: filename
    REAL(dp),ALLOCATABLE,TARGET :: zreal2d_global(:,:)
    REAL(dp),ALLOCATABLE        :: zreal2d(:,:)
    REAL(dp), POINTER           :: zreal2d_ptr(:,:) => NULL()

    IF (current_year < 10000) THEN
       WRITE (filename,'(a,I4.4,a)') "landuseTransitions.", current_year, ".nc"
    ELSE
       CALL finish('read_landuse_transitions','Only years between 0 and 9999 supported currently')
    END IF

    if (p_parallel_io) then
       ! Open input file
       call message('read_landuse_transitions()','Reading new landuse transition fields from '//trim(filename))
       IO_file%opened = .false.
       call IO_open(trim(filename), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Check resolution
       call IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       call IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)

       if (znlon /= grid%nlon .or. znlat /= grid%nlat) then
          write(message_text,*) 'Unexpected grid resolution: ', znlon, 'x', znlat
          call finish('read_landuse_transitions()', message_text)
       endif

       !! allocate temporary memory

       allocate(zreal2d_global(grid%nlon,grid%nlat),STAT=status)
       if(status .ne. 0) call finish('read_landuse_transitions()','Allocation failure (1)')

    endif !! end parallel_io


    !! === Read in LU-transition Natural => Crops

    varName="natural2crop"

    if (p_parallel_io) then

       zreal2d_global(:,:) = 0.0_dp
       call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal2d_global)

    endif !! end parallel_io

    !! Bring landuse transitions to the other processors

    allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_landuse_transitions()','Allocation failure NATL=>CROPS')
    nullify(zreal2d_ptr)
    if (p_parallel_io) zreal2d_ptr => zreal2d_global(:,:)
    call scatter_field(zreal2d_ptr, zreal2d)
    landuse_transitions%TransMtrx_NATL_2_CROP(:) = pack(zreal2d, MASK=domain%mask)
    deallocate(zreal2d)

    !! Do range check
    if(       ANY(landuse_transitions%TransMtrx_NATL_2_CROP(:) < 0.0_dp) &
         .or. ANY(landuse_transitions%TransMtrx_NATL_2_CROP(:) > 1.0_dp) ) then
       call finish('read_landuse_transitions()',&
            'NHHLP-Transition matrix element '//trim(varName)//' outside allowed range [0,1].')
    end if

    !! === Read in LU-transition Crops => Natural

    varName="crop2natural"

    if (p_parallel_io) then

       zreal2d_global(:,:) = 0.0_dp
       call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal2d_global)

    endif !! end parallel_io

    !! Bring landuse transitions to the other processors

    allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_landuse_transitions()','Allocation failure CROPS=>NATL')
    nullify(zreal2d_ptr)
    if (p_parallel_io) zreal2d_ptr => zreal2d_global(:,:)
    call scatter_field(zreal2d_ptr, zreal2d)
    landuse_transitions%TransMtrx_CROP_2_NATL(:) = pack(zreal2d, MASK=domain%mask)
    deallocate(zreal2d)

    !! Do range check
    if(       ANY(landuse_transitions%TransMtrx_CROP_2_NATL(:) < 0.0_dp) &
         .or. ANY(landuse_transitions%TransMtrx_CROP_2_NATL(:) > 1.0_dp) ) then
       call finish('read_landuse_transitions()',&
            'NHHLP-Transition matrix element '//trim(varName)//' outside allowed range [0,1].')
    end if

    !! === Read in LU-transition Natural => Pastures

    varName="natural2pasture"

    if (p_parallel_io) then

       zreal2d_global(:,:) = 0.0_dp
       call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal2d_global)

    endif !! end parallel_io

    !! Bring landuse transitions to the other processors

    allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_landuse_transitions()','Allocation failure NATL=>PAST')
    nullify(zreal2d_ptr)
    if (p_parallel_io) zreal2d_ptr => zreal2d_global(:,:)
    call scatter_field(zreal2d_ptr, zreal2d)
    landuse_transitions%TransMtrx_NATL_2_PAST(:) = pack(zreal2d, MASK=domain%mask)
    deallocate(zreal2d)
    
    !! Do range check
    if(       ANY(landuse_transitions%TransMtrx_NATL_2_PAST(:) < 0.0_dp) &
         .or. ANY(landuse_transitions%TransMtrx_NATL_2_PAST(:) > 1.0_dp) ) then
       call finish('read_landuse_transitions()',&
            'NHHLP-Transition matrix element '//trim(varName)//' outside allowed range [0,1].')
    end if

    !! === Read in LU-transition Pastures ==> Natural 

    varName="pasture2natural"

    if (p_parallel_io) then

       zreal2d_global(:,:) = 0.0_dp
       call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal2d_global)

    endif !! end parallel_io

    !! Bring landuse transitions to the other processors

    allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_landuse_transitions()','Allocation failure PAST=>NATL')
    nullify(zreal2d_ptr)
    if (p_parallel_io) zreal2d_ptr => zreal2d_global(:,:)
    call scatter_field(zreal2d_ptr, zreal2d)
    landuse_transitions%TransMtrx_PAST_2_NATL(:) = pack(zreal2d, MASK=domain%mask)
    deallocate(zreal2d)
    
    !! Do range check
    if(       ANY(landuse_transitions%TransMtrx_PAST_2_NATL(:) < 0.0_dp) &
         .or. ANY(landuse_transitions%TransMtrx_PAST_2_NATL(:) > 1.0_dp) ) then
       call finish('read_landuse_transitions()',&
            'NHHLP-Transition matrix element '//trim(varName)//' outside allowed range [0,1].')
    end if

    !! === Read in LU-transition Crops ==> Pastures 

    varName="crop2pasture"

    if (p_parallel_io) then

       zreal2d_global(:,:) = 0.0_dp
       call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal2d_global)

    endif !! end parallel_io

    !! Bring landuse transitions to the other processors

    allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_landuse_transitions()','Allocation failure CROP=>PAST')
    nullify(zreal2d_ptr)
    if (p_parallel_io) zreal2d_ptr => zreal2d_global(:,:)
    call scatter_field(zreal2d_ptr, zreal2d)
    landuse_transitions%TransMtrx_CROP_2_PAST(:) = pack(zreal2d, MASK=domain%mask)
    deallocate(zreal2d)

    !! Do range check
    if(       ANY(landuse_transitions%TransMtrx_CROP_2_PAST(:) < 0.0_dp) &
         .or. ANY(landuse_transitions%TransMtrx_CROP_2_PAST(:) > 1.0_dp) ) then
       call finish('read_landuse_transitions()',&
            'NHHLP-Transition matrix element '//trim(varName)//' outside allowed range [0,1].')
    end if

    !! === Read in LU-transition Pastures => Crops

    varName="pasture2crop"

    if (p_parallel_io) then
       
       zreal2d_global(:,:) = 0.0_dp
       call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal2d_global)

    endif !! end parallel_io

    !! Bring landuse transitions to the other processors

    allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_landuse_transitions()','Allocation failure PAST=>CROP')
    nullify(zreal2d_ptr)
    if (p_parallel_io) zreal2d_ptr => zreal2d_global(:,:)
    call scatter_field(zreal2d_ptr, zreal2d)
    landuse_transitions%TransMtrx_PAST_2_CROP(:) = pack(zreal2d, MASK=domain%mask)
    deallocate(zreal2d)

    !! Do range check
    if(       ANY(landuse_transitions%TransMtrx_PAST_2_CROP(:) < 0.0_dp) &
         .or. ANY(landuse_transitions%TransMtrx_PAST_2_CROP(:) > 1.0_dp) ) then
       call finish('read_landuse_transitions()',&
            'NHHLP-Transition matrix element '//trim(varName)//' outside allowed range [0,1].')
    end if

    !! Free temporary memory
    
    if (p_parallel_io) deallocate(zreal2d_global)

    !! Close input file
    
    if (p_parallel_io) call IO_close(IO_file)

!!$       call message('read_landuse_transitions()','Reading of NHHLP-Transition matrix elements finished.') 

    !! Do sum checks for transition elements

    IF( ANY(landuse_transitions%TransMtrx_PAST_2_CROP(:) + landuse_transitions%TransMtrx_PAST_2_NATL(:) > 1.0_dp) ) THEN
       WRITE (message_text, *) ' Sum of transitions from pasture has value ', &
            MAXVAL(landuse_transitions%TransMtrx_PAST_2_CROP(:) + landuse_transitions%TransMtrx_PAST_2_NATL(:)), &
            ' larger 1 at grid cell ', &
            MAXLOC(landuse_transitions%TransMtrx_PAST_2_CROP(:) + landuse_transitions%TransMtrx_PAST_2_NATL(:))
       CALL message('read_landuse_transitions()', message_text)
    END IF
       
    IF( MAXVAL(landuse_transitions%TransMtrx_CROP_2_PAST(:) + landuse_transitions%TransMtrx_CROP_2_NATL(:)) > 1.0_dp ) THEN
       WRITE (message_text, *) ' Sum of transitions from crops has value ', &
            MAXVAL(landuse_transitions%TransMtrx_CROP_2_PAST(:) + landuse_transitions%TransMtrx_CROP_2_NATL(:)), &
            ' larger 1 at grid cell ', &
            MAXLOC(landuse_transitions%TransMtrx_CROP_2_PAST(:) + landuse_transitions%TransMtrx_CROP_2_NATL(:))
       CALL message('read_landuse_transitions()', message_text)
    END IF

    IF( MAXVAL(landuse_transitions%TransMtrx_NATL_2_PAST(:) + landuse_transitions%TransMtrx_NATL_2_CROP(:)) > 1.0_dp ) THEN
       WRITE (message_text, *) ' Sum of transitions from natural has value ', &
            MAXVAL(landuse_transitions%TransMtrx_NATL_2_PAST(:) + landuse_transitions%TransMtrx_NATL_2_CROP(:)), &
            ' larger 1 at grid cell ', &
            MAXLOC(landuse_transitions%TransMtrx_NATL_2_PAST(:) + landuse_transitions%TransMtrx_NATL_2_CROP(:))
       CALL message('read_landuse_transitions()', message_text)
    END IF

  end subroutine read_landuse_transitions

  !! --- read_harvest() ----------------------------------------------------
  !!
  !! This routine reads in the harvest data of the (reduced) New Hampshire Harmonized Landuse Protocol 
  !! for a new year from an external file. The transitions are stored in landuse_transitions%harvest().

  SUBROUTINE  read_harvest(current_year,grid,domain)
    use mo_jsbach_grid,      only: grid_type,domain_type
    use mo_exception,        only: finish, message
    use mo_netcdf,           only: FILE_INFO,IO_inq_dimid,IO_inq_dimlen,IO_inq_varid,io_get_var_double
    USE mo_tr_scatter,       ONLY: scatter_field
    use mo_io,               only: IO_open, IO_READ, IO_close

    INTEGER, intent(in)             :: current_year
    type(grid_type),   intent(in)   :: grid
    type(domain_type), intent(in)   :: domain

    !! --- parameters

    real(dp),parameter :: conversionFactorHarvest= 1._dp !! Used to convert harvest data from the unit in the input file
                                                         !! to mol(C)/m^2(gridbox)/s

    !! --- locals

    integer                     :: IO_file_id, IO_var_id, IO_dim_id
    type(FILE_INFO)             :: IO_file
    integer                     :: znlon, znlat
    integer                     :: status
    character(len=64)           :: varName 
    character(len=1024)         :: filename
    REAL(dp),ALLOCATABLE,TARGET :: zreal2d_global(:,:)
    REAL(dp),ALLOCATABLE        :: zreal2d(:,:)
    REAL(dp), POINTER           :: zreal2d_ptr(:,:) => NULL()

    IF (current_year < 10000) THEN
       WRITE (filename,'(a,I4.4,a)') "landuseHarvest.", current_year, ".nc"
    ELSE
       CALL finish('read_harvest','Only years between 0 and 9999 supported currently')
    END IF

    if (p_parallel_io) then
       ! Open input file
       call message('read_harvest()','Reading new harvest data from '//trim(filename))
       IO_file%opened = .false.
       call IO_open(trim(filename), IO_file, IO_READ)
       IO_file_id = IO_file%file_id

       ! Check resolution
       call IO_inq_dimid  (IO_file_id, 'lat', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlat)
       call IO_inq_dimid  (IO_file_id, 'lon', IO_dim_id)
       call IO_inq_dimlen (IO_file_id, IO_dim_id, znlon)

       if (znlon /= grid%nlon .or. znlat /= grid%nlat) then
          write(message_text,*) 'Unexpected grid resolution: ', znlon, 'x', znlat
          call finish('read_harvest()', message_text)
       endif

       !! allocate temporary memory

       allocate(zreal2d_global(grid%nlon,grid%nlat),STAT=status)
       if(status .ne. 0) call finish('read_harvest()','Allocation failure (1)')

    endif !! end parallel_io

    !! === Read in harvest data

    varName="harvest"

    if (p_parallel_io) then

       zreal2d_global(:,:) = 0.0_dp
       call IO_inq_varid(IO_file_id, trim(varName), IO_var_id)
       call IO_get_var_double(IO_file_id, IO_var_id, zreal2d_global)

       !! Do range check
       if( ANY(zreal2d_global < 0.0_dp) ) then
          call finish('read_harvest()',&
               'Found negative harvest in '//trim(varName))
       end if

    endif !! end parallel_io

    !! Bring harvest data to the other processors

    allocate(zreal2d(domain%ndim,domain%nblocks),STAT=status)
    if(status .ne. 0) call finish('read_harvest()','Allocation failure')
    nullify(zreal2d_ptr)
    if (p_parallel_io) zreal2d_ptr => zreal2d_global(:,:)
    call scatter_field(zreal2d_ptr, zreal2d)
    !! Save harvest into stream but convert also to mol(C)/m^2(gridbox)/s
    landuse_transitions%Box_harvest(:) = conversionFactorHarvest*pack(zreal2d, MASK=domain%mask)

    deallocate(zreal2d)

    !! Free temporary memory

    if (p_parallel_io) deallocate(zreal2d_global)

    !! Close input file

    if (p_parallel_io) call IO_close(IO_file)

  END subroutine read_harvest
  
end module mo_cbal_landcover_change
