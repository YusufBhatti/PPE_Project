!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_moz_lightning.f90
!!
!! \brief
!! Calculate lightning flashes and NOx emissions from lightning for chemistry
!!
!! \author Martin Schultz (FZ Juelich)
!!
!! \responsible_coder
!! Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# M. Schultz (FZ Juelich) - original code (XXXX-XX-XX)
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

MODULE mo_moz_lightning

  USE mo_kind,           ONLY: dp
  USE mo_time_event,     ONLY: io_time_event


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: start_lightning
  PUBLIC :: init_lightning_stream
  PUBLIC :: moz_lightning

  ! -- public diagnostic pointers
  REAL(dp), POINTER, PUBLIC  :: dpavupdr(:,:) => NULL()         ! cloud average updraft velocity
  REAL(dp), POINTER, PUBLIC  :: dpavupdr_mean(:,:) => NULL()    ! mean cloud average updraft velocity
  REAL(dp), POINTER, PUBLIC  :: dpff(:,:) => NULL()             ! flash frequency
  REAL(dp), POINTER, PUBLIC  :: dpff_mean(:,:) => NULL()        ! mean flash frequency
  REAL(dp), POINTER, PUBLIC  :: dpemino(:,:) => NULL()          ! NO emission flux (column integral)
  REAL(dp), POINTER, PUBLIC  :: dpemino_accu(:,:) => NULL()     ! accumulated NO emission flux
  REAL(dp), POINTER, PUBLIC  :: dpcloudbot(:,:) => NULL()       ! cloud bottom height
  REAL(dp), POINTER, PUBLIC  :: dpcloudtop(:,:) => NULL()       ! cloud top height
  REAL(dp), POINTER, PUBLIC  :: dpcloudthick(:,:) => NULL()     ! cold cloud thickness
  REAL(dp), POINTER, PUBLIC  :: dpcloudfreezbot(:,:) => NULL()  ! altitude of lowest cloud level below freezing
  REAL(dp), POINTER, PUBLIC  :: dpcgfrac(:,:) => NULL()         ! fraction of cloud-ground lightning
  REAL(dp), POINTER, PUBLIC  :: dpupdr(:,:,:) => NULL()         ! 3D updraft velocity
  REAL(dp), POINTER, PUBLIC  :: dpemino3d(:,:,:) => NULL()      ! 3D NO emissions

  ! -- private module variables
  REAL(dp)             :: ffscale       ! resolution dependent scaling parameter for flash frequency
  LOGICAL              :: ldetail       ! generate detailed diagnostic output
  TYPE(io_time_event), SAVE  :: tinterval     ! stream output interval


CONTAINS

!
! ---
!

SUBROUTINE start_lightning

  USE mo_time_control,          ONLY : putdata

! set resolution dependent parameters, read namelist, initialise diag

  tinterval = putdata              ! set defaut lightning stream output interval to general ECHAM value
  ldetail = .TRUE.

END SUBROUTINE start_lightning

!
! ---
!

SUBROUTINE init_lightning_stream

  USE mo_linked_list,         ONLY: t_stream, SURFACE, HYBRID
  USE mo_filename,            ONLY: trac_filetype
  USE mo_memory_base,         ONLY: AUTO, new_stream, default_stream_setting,     &
                                    add_stream_reference, add_stream_element

  TYPE (t_stream), POINTER   :: stream_lghtng
  LOGICAL                    :: lpost

  ! Create stream and define basic settings
  CALL new_stream(stream_lghtng,'lghtng', filetype=trac_filetype,interval=tinterval,lrerun=.FALSE.)
  CALL default_stream_setting (stream_lghtng,              &
                               lrerun    = .FALSE.,      &
                               laccu     =  .TRUE. ,     &
                               lpost     =  .TRUE.,      &
                               leveltype =  SURFACE,     &
                               table     =  199,         &
                               code      =  AUTO         )

  CALL add_stream_reference (stream_lghtng, 'geosp'   ,'g3b'   )
  CALL add_stream_reference (stream_lghtng, 'aps'     ,'g3b'   )
  CALL add_stream_reference (stream_lghtng, 'gboxarea','geoloc')

  ! Add cumulative diagnostics (will always be output)
  CALL add_stream_element (stream_lghtng, 'cloudavupdr_mean',  dpavupdr_mean,     &
                           longname='cloud averaged mean updraft velocity',       &
                           units='m s-1')
  CALL add_stream_element (stream_lghtng, 'flashfreq_mean',  dpff_mean,           &
                           longname='mean flash frequency',                       &
                           units='s-1')
  CALL add_stream_element (stream_lghtng, 'emi_NO_lghtng_accu',  dpemino_accu,    &
                           longname='accumulated lightning NO emissions',         &
                           units='kg(N) output_interval-1')

  ! Detailed diagnostics (no averaging or accumulation)
  CALL default_stream_setting (stream_lghtng, laccu=.FALSE., lpost=.TRUE.)
  ! dpemino needed for global diagnostics
  CALL add_stream_element (stream_lghtng, 'emi_NO_lghtng',  dpemino,              &
                           longname='lightning NO emissions',                     &
                           units='kg(N) s-1', lpost=ldetail)
  IF (ldetail) THEN
    ! 2D fields
    CALL add_stream_element (stream_lghtng, 'cloudavupdr',  dpavupdr,             &
                             longname='cloud averaged updraft velocity',          &
                             units='m s-1')
    CALL add_stream_element (stream_lghtng, 'flashfreq',  dpff,                   &
                             longname='flash frequency',                          &
                             units='s-1')
    CALL add_stream_element (stream_lghtng, 'cloudbot',  dpcloudbot,              &
                             longname='height of cloud bottom',                   &
                             units='m')
    CALL add_stream_element (stream_lghtng, 'cloudtop',  dpcloudtop,              &
                             longname='height of cloud top',                      &
                             units='m')
    CALL add_stream_element (stream_lghtng, 'cloudthick',  dpcloudthick,          &
                             longname='cold cloud thickness',                     &
                             units='m')
    CALL add_stream_element (stream_lghtng, 'cloudfreezbot',  dpcloudfreezbot,    &
                             longname='height of lowest cloud level below freezing temperature', &
                             units='m')
    CALL add_stream_element (stream_lghtng, 'cgfrac',  dpcgfrac,                  &
                             longname='fraction of cloud-to-ground lightning strikes',           &
                             units='1')
    ! 3D fields
    CALL default_stream_setting (stream_lghtng, leveltype = HYBRID, lpost=.TRUE.)
    CALL add_stream_element (stream_lghtng, 'cloudupdr',  dpupdr,                 &
                             longname='updraft velocity',                         &
                             units='m s-1')
    CALL add_stream_element (stream_lghtng, 'emi_NO_lghtng_3d',  dpemino3d,       &
                             longname='updraft velocity',                         &
                             units='kg(N) s-1')
  END IF

END SUBROUTINE init_lightning_stream

!
! ---
!

SUBROUTINE moz_lightning (kproma, kbdim,  klev,    krow,             &
                          ktype,  kcbot,  kctop,                     &
                          pfrl,   pten,   ptenh,  paphp1,            &
                          pmfu,   pnoems_3d                 )

!! Calculates lightning flash density and NOx emissions from lightning.
!! Subroutine moz_lightning is called from cuflx_subm in mo_submodel_interface

  USE mo_kind,                ONLY: dp
  USE mo_physical_constants,  ONLY: rd
  USE mo_vphysc,              ONLY: vphysc     ! grheightm1, grmassm1
  USE mo_time_control,        ONLY: time_step_len   
  USE mo_geoloc,              ONLY: philat_2d
  USE mo_exception,           ONLY: message, message_text, em_info

  INTEGER, INTENT(in)      :: kproma                       ! column index
  INTEGER, INTENT(in)      :: kbdim                        ! block size (= max(kproma))
  INTEGER, INTENT(in)      :: klev                         ! number of levels
  INTEGER, INTENT(in)      :: krow                         ! block number
  INTEGER, INTENT(in)      :: ktype(kbdim)                 ! cloud type
  INTEGER, INTENT(in)      :: kcbot(kbdim)                 ! model level for cloud bottom
  INTEGER, INTENT(in)      :: kctop(kbdim)                 ! model level for cloud top
  REAL(dp), INTENT(in)     :: pfrl(kbdim)                  ! land fraction
  REAL(dp), INTENT(in)     :: pten(kbdim, klev)            ! ???
  REAL(dp), INTENT(in)     :: ptenh(kbdim, klev)           ! ???
  REAL(dp), INTENT(in)     :: paphp1(kbdim, klev)          ! half level pressures
  REAL(dp), INTENT(in)     :: pmfu(kbdim, klev)            ! convective updraft mass flux
  REAL(dp), INTENT(inout)  :: pnoems_3d(kbdim, klev)       ! NO emission flux [kg(NO)/s]

  REAL(dp), PARAMETER      :: ffconst = 1.54e-5_dp         ! parameter for flash frequency according to Grewe
  REAL(dp), PARAMETER      :: ffexp   = 4.9_dp             ! parameter for flash frequency according to Grewe
  REAL(dp), PARAMETER      :: mw_N    = 14.00674_dp        ! molecular weight of N
  REAL(dp), PARAMETER      :: mw_NO   = 30.00614_dp        ! molecular weight of NO

  INTEGER      :: jl, jk 
  INTEGER      :: ncfreezbot(kbdim)                        ! bottom level of ice cloud (T <= 273.15 K)
  REAL(dp)     :: zupdr(kbdim, klev)                       ! updraft velocity [m/s]
  REAL(dp)     :: zavupdr(kbdim)                           ! average updraft velocity in cloud
  REAL(dp)     :: zcbot(kbdim), zctop(kbdim)               ! cloud bottom/topi level index as real value
  REAL(dp)     :: zcfreezbot(kbdim)                        ! bottom level of ice cloud (T <= 273.15 K)
  REAL(dp)     :: zrhoa                                    ! air density
  REAL(dp)     :: zcbh(kbdim), zcth(kbdim), zcdh(kbdim)    ! cloud bottom/top/cold thickness in m
  REAL(dp)     :: zcdhkm                                   ! zzcdh in km
  REAL(dp)     :: zff(kbdim)                               ! flash frequency [s-1]
  REAL(dp)     :: zfraccg(kbdim)                           ! fraction of cloud-ground flashes
  REAL(dp)     :: pnopcgf                                  ! NO production [kg(N)] per cloud-ground flash
  REAL(dp)     :: pnopccf                                  ! NO production [kg(N)] per cloud-cloud flash
  REAL(dp)     :: zpnocg(kbdim), zpnocc(kbdim), zpno(kbdim) ! NO produced per second
  REAL(dp)     :: znoems(klev)                             ! NO emission flux [kg(N)/s] in one column

  LOGICAL, SAVE     :: lfirst = .TRUE.

  ! ### move to lght_init (if needed)!!
  !!ffscale = 0.2_dp    ! estimate for T31  #########   should yield ~50 flashes/second globally
  ! should be higher according to Martin Schultz
  !ffscale = 0.35_dp    ! estimate for T63  #########   should yield ~50 flashes/second globally
  ffscale = 0.7_dp    ! estimate for T63  #########   should yield ~50 flashes/second globally
  !++ Scarlet: pnopcgf should be 16.8 according to EMAC
  !pnopcgf = 1._dp    ! estimate for T31  #########   should yield ~5 TgN/year globally
  pnopcgf = 8.4_dp ! Reduced to have half NOx emissions still yielding between 2 and 5 TgN/year globally 
  pnopccf = 0.1_dp * pnopcgf    ! assume that cloud-cloud flashes only produce 10% as much NO
  IF (lfirst) THEN  
    WRITE(message_text,*) "Flash frequency scale factor (ffscale)", ffscale
    CALL message('moz_lightning', message_text, level=em_info)
    WRITE(message_text,*) "NO production per cloud to ground flash (pnopcgf)", pnopcgf
    CALL message('moz_lightning', message_text, level=em_info)
    WRITE(message_text,*) "NO production per cloud to cloud flash (pnopccf)", pnopccf
    CALL message('moz_lightning', message_text, level=em_info)
  END IF
  lfirst = .FALSE.


  ! Initialisation
  ncfreezbot(1:kproma) = klev
  zupdr(1:kproma,:) =    0._dp
  zavupdr(1:kproma) =    0._dp
  zcbh(1:kproma) =       0._dp
  zcth(1:kproma) =       0._dp
  zcdh(1:kproma) =       0._dp
  zcfreezbot(1:kproma) = 0._dp
  zff(1:kproma) =        0._dp
  zfraccg(1:kproma) =    0._dp
  zpnocg(1:kproma) =     0._dp
  zpnocc(1:kproma) =     0._dp
  zpno(1:kproma) =       0._dp
  znoems(:) =            0._dp

  ! part 1 (formerly in cuflx)
  WHERE (ktype(1:kproma) == 1)
    zcbot(1:kproma) = REAL(kcbot(1:kproma),dp)
    zctop(1:kproma) = REAL(kctop(1:kproma),dp)
  ELSEWHERE
    zcbot(1:kproma) = 0._dp
    zctop(1:kproma) = 0._dp
  END WHERE

  DO jl=1,kproma
    IF (ktype(jl) == 1) THEN
      ! Determine cold cloud thickness
      ncfreezbot(jl) = kctop(jl)
      DO jk=kctop(jl), kcbot(jl)
        IF (pten(jl,jk) <= 273.15_dp ) ncfreezbot(jl) = jk
      END DO
      ! Calculate updraft velocity
      DO jk=kctop(jl), kcbot(jl)
        zrhoa = paphp1(jl,jk)/(ptenh(jl,jk)*rd)
        zupdr(jl,jk) = pmfu(jl,jk)/zrhoa
      END DO
    END IF
  END DO

  ! part 2
  ! Determine cloud bottom and top and cold cloud thickness in m
  DO jl=1, kproma
    IF (kctop(jl) < kcbot(jl) .AND. ncfreezbot(jl) < klev) THEN
      DO jk=klev, kcbot(jl), -1
        zcbh(jl) = zcbh(jl) + vphysc%grheightm1(jl, jk, krow)
      END DO 
      DO jk=klev, kctop(jl), -1
        zcth(jl) = zcth(jl) + vphysc%grheightm1(jl, jk, krow)
      END DO 
      DO jk=kctop(jl), ncfreezbot(jl)
        zcdh(jl) = zcdh(jl) + vphysc%grheightm1(jl, jk, krow)
      END DO 
      IF (ldetail) THEN
        DO jk=klev, ncfreezbot(jl), -1
          zcfreezbot(jl) = zcfreezbot(jl) + vphysc%grheightm1(jl, jk, krow)
        END DO
      END IF
      ! Compute average updraft velocity weighted by grid cell height
      DO jk=kctop(jl), kcbot(jl)
        zavupdr(jl) = zavupdr(jl) + zupdr(jl,jk)*vphysc%grheightm1(jl, jk, krow)
      END DO
      zavupdr(jl) = zavupdr(jl) / ( zcth(jl) - zcbh(jl) )
    END IF
  END DO

  ! part 3
  ! Calculate flash frequency according to Grewe et al., 2001 
  ! ffscale : resolution dependent scaling parameter
  ! ffconst : scale factor according to Grewe paper (actually also resolution dependent: needs tuning!)
  ! ffexp   : exponent according to Grewe paper
  ! ### need to re-adjust parameters - resolution dependent!! ###

  ! ++ Scarlet: bug fix, the unit conversion is missing here. This flashes/min,
  ! but it should be flashes/sec
  ! zff(1:kproma) = ffscale * ffconst * (zavupdr(1:kproma) * sqrt( zcth(1:kproma)-zcbh(1:kproma) ))**ffexp
  zff(1:kproma) = (ffscale * ffconst * (zavupdr(1:kproma) * sqrt( zcth(1:kproma)-zcbh(1:kproma) ))**ffexp)/60._dp
  ! -- Scarlet

  ! part 4
  ! Determine the fraction of cloud-ground lightning strikes (Price and Rind, 1995)
  ! and compute NO production
  DO jl=1, kproma
     IF (zcdh(jl) > 0._dp) THEN
        IF (zcdh(jl) < 5500._dp) THEN
           zfraccg(jl) = 1._dp
        ELSE IF (zcdh(jl) < 14000._dp) THEN
           zcdhkm = zcdh(jl) * 1.e-3_dp
           zfraccg(jl) = 0.021_dp * zcdhkm - 0.648_dp
           zfraccg(jl) = zfraccg(jl) * zcdhkm + 7.49_dp
           zfraccg(jl) = zfraccg(jl) * zcdhkm - 36.54_dp
           zfraccg(jl) = zfraccg(jl) * zcdhkm + 64.09_dp
           zfraccg(jl) = 1._dp/zfraccg(jl)
        ELSE
           zfraccg(jl) = 0.02_dp
        END IF
        ! evaluate NO production only if at least 1 flash per time step
        IF (zff(jl)*time_step_len >= 1._dp) THEN
           zpnocg(jl) = pnopcgf * zfraccg(jl) * zff(jl)
           zpnocc(jl) = pnopccf * (1._dp - zfraccg(jl)) * zff(jl)
           zpno(jl) = zpnocg(jl) + zpnocc(jl)
        END IF
     END IF
  END DO

  ! part 5
  ! Vertical distribution of NO
  DO jl=1, kproma
     CALL zprofile(klev, kcbot(jl), kctop(jl),                       &
                   philat_2d(jl, krow),                              &
                   pfrl(jl),                                         &
                   vphysc%grheightm1(jl, :, krow),                   &
                   zpno(jl),                                         &
                   znoems(:)                            )
     ! convert into mass mixing ratio
     pnoems_3d(jl,:) = mw_NO/mw_N * znoems(:) / vphysc%grmassm1(jl, :, krow)
  END DO

  ! diagnostics
  dpavupdr_mean(1:kproma, krow) = dpavupdr_mean(1:kproma, krow)      &
                                  + zavupdr(1:kproma) * time_step_len
  dpff_mean(1:kproma, krow)     = dpff_mean(1:kproma, krow)          &
                                  + zff(1:kproma) * time_step_len
  dpemino_accu(1:kproma, krow)  = dpemino_accu(1:kproma, krow)       &
                                  + zpno(1:kproma) * time_step_len      ! ### CONVERSION!! ###
                                  ! * gridcell_area(?)
  IF (ldetail) THEN
    dpavupdr(1:kproma, krow)        = zavupdr(1:kproma)
    dpff(1:kproma, krow)            = zff(1:kproma)
    dpemino(1:kproma, krow)         = zpno(1:kproma)      ! ### UNITS ? !! ###
    dpcloudbot(1:kproma, krow)      = zcbot(1:kproma) 
    dpcloudtop(1:kproma, krow)      = zctop(1:kproma)
    dpcloudthick(1:kproma, krow)    = zcdh(1:kproma)
    dpcloudfreezbot(1:kproma, krow) = zcfreezbot(1:kproma)
    dpcgfrac(1:kproma, krow)        = zfraccg(1:kproma)
    dpupdr(1:kproma,:, krow)        = zupdr(1:kproma,:)
    dpemino3d(1:kproma,:, krow)     = pnoems_3d(1:kproma,:)
  END IF

END SUBROUTINE moz_lightning


SUBROUTINE zprofile(klev, kcbot, kctop, plat, pfrl, pheight, ppno, pnoems)

  ! Estimate the vertical distribution of lightning generated NO 
  ! Pickering et al., JGR 103(D23), 31203-31216 (1998)
  USE mo_kind,        ONLY: dp

  INTEGER, INTENT(in)      :: klev
  INTEGER, INTENT(in)      :: kcbot            ! bottom level of cloud
  INTEGER, INTENT(in)      :: kctop            ! top level of cloud
  REAL(dp), INTENT(in)     :: plat             ! latitude
  REAL(dp), INTENT(in)     :: pfrl             ! land fraction of surface grid cell
  REAL(dp), INTENT(in)     :: pheight(klev)    ! model level heights in m
  REAL(dp), INTENT(in)     :: ppno             ! NO production in column
  REAL(dp), INTENT(inout)  :: pnoems(klev)     ! vertically resolved NO production

  INTEGER            :: jk, nreg
  INTEGER            :: kl_low, kl_high
  REAL(dp)           :: zh_low, zh_high
  REAL(dp)           :: zpp(16,3)              ! C-profile shape values
  REAL(dp)           :: zz(klev)               ! level altitude in m
  REAL(dp)           :: zqp(klev)              ! interpolate
  REAL(dp)           :: zqs                    ! interpolate norm
  REAL(dp)           :: zd

  ! Empirical profiles from Pickering study [1 value per km]
  ! midlatitude continental
  zpp(:,1) = (/20.1_dp,  2.3_dp,  0.8_dp,  1.5_dp,  3.4_dp,  5.3_dp,  3.6_dp,  3.8_dp, &
                5.4_dp,  6.6_dp,  8.3_dp,  9.6_dp, 12.8_dp, 10.0_dp,  6.2_dp,  0.3_dp  /)
  ! tropical marine
  zpp(:,2) = (/ 5.8_dp,  2.9_dp,  2.6_dp,  2.4_dp,  2.2_dp,  2.1_dp,  2.3_dp,  6.1_dp, &
               16.5_dp, 14.1_dp, 13.7_dp, 12.8_dp, 12.5_dp,  2.8_dp,  0.9_dp,  0.3_dp/)
  ! tropical continental
  zpp(:,3) = (/ 8.2_dp,  1.9_dp,  2.1_dp,  1.6_dp,  1.1_dp,  1.6_dp,  3.0_dp,  5.8_dp, &
                7.6_dp,  9.6_dp, 10.5_dp, 12.3_dp, 11.8_dp, 12.5_dp,  8.1_dp,  2.3_dp/)

  ! initialize pnoems
  pnoems(:) = 0._dp

  ! return immediately if no NO is produced in this column
  IF (abs(ppno) < 1.e-15_dp) RETURN

  ! calculate heights in meters
  zz(klev) = pheight(klev) * 0.5_dp
  DO jk = klev-1, kctop, -1
     zz(jk) = zz(jk+1) + (pheight(jk+1) + pheight(jk)) * 0.5_dp
  ENDDO
  ! adjust to interval [0.5,15.5] km
  zz(klev) = 0.5_dp
  IF (kctop < klev) THEN
     zd = 1.0_dp / (zz(kctop) - zz(klev))
     DO jk = kctop, klev-1
        zz(jk) = (zz(jk) - zz(klev)) * zd * 15.0_dp + 0.5_dp
     END DO
  END IF
  ! select vertical profile according to region
  IF ( abs(plat) > 30._dp ) THEN
     nreg = 1
  ELSE
     IF ( pfrl < 0.5_dp ) THEN   !land points
        nreg = 2
     ELSE
        nreg = 3
     END IF
  END IF

  ! linear interpolation
  zqs = 0._dp
  DO jk = klev, kctop+1, -1
     ! calculate indices 
     kl_low  = INT(zz(jk) - 0.4999999999999999_dp) + 1
     kl_high = kl_low + 1
     zh_low  = REAL(kl_low,dp) - 0.5_dp
     zh_high = zh_low + 1.0_dp
     ! interpolation formula, attention: zh(kl_high)-zh(kl_low)=1
     zqp(jk) = zpp(kl_high, nreg) * (zz(jk) - zh_low  ) + &
                 zpp(kl_low,  nreg) * (zh_high  - zz(jk))
     zqs = zqs + zqp(jk)
  ENDDO
  zqp(kctop) = zpp(16, nreg)
  zqs = zqs + zqp(kctop)
  zqs = 1.0_dp / zqs
  DO jk = kctop, klev
     pnoems(jk) = ppno * zqp(jk) * zqs
  END DO

END SUBROUTINE zprofile


END MODULE mo_moz_lightning

