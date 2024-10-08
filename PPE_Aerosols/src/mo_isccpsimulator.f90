!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_isccpsimulator
 !
 ! *mo_isccpsimulator*
 !
 ! Purpose:
 ! --------
 ! Model output for comparisons with ISCCP satellite data
 ! 
 ! The ISCCP simulator is written by Steve Klein (stephen.klein@noaa.gov)
 ! and Mark Webb (mark.webb@metoffice.com)
 ! See http://gcss-dime.giss.nasa.gov/simulator.html
 ! 
 ! Version v3.4.1 implemented by Johannes Quaas, 20/01/2005
 !
 !
  USE mo_kind,           ONLY: dp, wp
  USE mo_linked_list,    ONLY: t_stream,  HYBRID, SURFACE, BELOWSUR, &
                               ABOVESUR10, NETCDF
  USE mo_memory_base,    ONLY: new_stream, add_stream_element,     &
                               default_stream_setting, add_stream_reference, &
                               AUTO
  USE mo_random_numbers, ONLY: random_uniform => get_random

  IMPLICIT NONE

  TYPE, PUBLIC :: vmem2d
     REAL(dp), POINTER :: ptr(:,:)
  END TYPE vmem2d

  INTEGER, PARAMETER :: ncol=40

 ! Variables
  REAL(dp), SAVE :: tautab(0:255)      
  INTEGER, SAVE :: invtau(-20:45000)     
  
  TYPE (vmem2d), PRIVATE, POINTER :: isccp_cldtypes(:)   ! ISCCP cloud fraction ofr the 49 ISCCP 
                                                         !  cloud types
  REAL(dp), PRIVATE, POINTER :: isccp_cldfra(:, :)   ! ISCCP total cloud fraction
  REAL(dp), PRIVATE, POINTER :: isccp_cldtau(:, :)   ! ISCCP mean cloud optical thickness
  REAL(dp), PRIVATE, POINTER :: isccp_cldptop(:, :)  ! ISCCP mean cloud top pressure
  REAL(dp), PRIVATE, POINTER :: isccp_cldfreq(:, :)  ! Frequency of cloud occurrence
  REAL(dp), PUBLIC,  POINTER :: isccp_sunlit(:, :)   ! Gridbox sunlit (1.) or not (0.)
  REAL(dp), PUBLIC,  POINTER :: isccp_cldtau3d(:, :, :)  ! Cloud optical thickness
  REAL(dp), PUBLIC,  POINTER :: isccp_cldemi3d(:, :, :)  ! Cloud emissivity @ 10.5 µm
  REAL(dp), PUBLIC,  POINTER :: isccp_f3d(:, :, :) 
  

  TYPE (t_stream), PUBLIC, POINTER :: isccp

CONTAINS
  
!--------------------------------------------------------------------------------------------------

  SUBROUTINE isccp_initialize

    USE mo_mpi,         ONLY: p_parallel, p_parallel_io, p_bcast, p_io
    USE mo_filename,    ONLY: find_next_free_unit
    USE mo_submodel,    ONLY: loisccp
    USE mo_exception,   ONLY: finish

    IMPLICIT NONE
    
    !--- Local variables:

    INTEGER :: iunit
    
    CALL finish('isccp_initialize','introduce new namelist with satoverpasstime')

    !--- 1) If to use the ISCCP simulator, read in the data
    IF (loisccp .AND. p_parallel_io) THEN

       iunit = find_next_free_unit (20, 99)

       OPEN(iunit, FILE='tautab.formatted', FORM='FORMATTED')
       READ(iunit,'(f30.20)') tautab
       CLOSE(iunit)
            
       OPEN(iunit, FILE='invtau.formatted', FORM='FORMATTED')
       READ(iunit,'(i10)') invtau
       CLOSE(iunit)

    ENDIF

    IF (p_parallel) THEN
       CALL p_bcast (tautab, p_io)
       CALL p_bcast (invtau, p_io)
    END IF ! p_parallel

  END SUBROUTINE isccp_initialize

!--------------------------------------------------------------------------------------------------

  SUBROUTINE construct_stream_isccp

    IMPLICIT NONE
    INTEGER i
    CHARACTER*40 :: name
    !
    ! Local
    ! 


    !--- 1) Construct the isccp stream: -----------------------------------------------------------

    CALL new_stream (isccp,'isccp',filetype=NETCDF)

    !--- Add standard fields for post-processing:

    CALL add_stream_reference (isccp, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (isccp, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (isccp, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (isccp, 'gboxarea','geoloc',lpost=.TRUE.)

    !--- 2) Add stream elements: ------------------------------------------------------------------

    CALL default_stream_setting (isccp,  units     = '',          &
                                         lrerun    = .TRUE. ,     &
                                         lpost     = .TRUE. ,     &
                                         laccu     = .FALSE. ,    &
                                         leveltype = SURFACE ,    &
                                         contnorest = .TRUE. )

    ! 2d mean quantities
    CALL add_stream_element (isccp, 'sunlit', isccp_sunlit,  &
                             longname='Sunlit fraction')
    CALL add_stream_element (isccp, 'cldfra', isccp_cldfra,  &
                             longname='Total cloud fraction')
    CALL add_stream_element (isccp, 'cldtau', isccp_cldtau,  &
                             longname='Mean cloud optical thickness')
    CALL add_stream_element (isccp, 'cldptop', isccp_cldptop,  &
                             longname='Mean cloud top pressure', units='Pa')
    CALL add_stream_element (isccp, 'cldfreq', isccp_cldfreq,  &
                             longname='Frequency of occurrence of clouds')

    CALL add_stream_element (isccp, 'tau3d', isccp_cldtau3d, &
                             leveltype = HYBRID, lpost=.FALSE., &
                             longname='cloud optical thickness')

    CALL add_stream_element (isccp, 'emi3d', isccp_cldemi3d, &
                             leveltype = HYBRID, lpost=.FALSE., &
                             longname='cloud optical thickness')

    CALL add_stream_element (isccp, 'f3d', isccp_f3d, &
                             leveltype = HYBRID, lpost=.FALSE., &
                             longname='cloud fraction')

    ALLOCATE(isccp_cldtypes(49))
    DO i=1,49
       IF ( i.LT.10 ) THEN
          WRITE(name, '(a8,i1)') 'cldtype0',i
       ELSE
          WRITE(name, '(a7,i2)') 'cldtype',i
       ENDIF
       CALL add_stream_element (isccp, name, isccp_cldtypes(i)%ptr,  &
            longname='Fraction covered by ISCCP cloud type')
    ENDDO

  END SUBROUTINE construct_stream_isccp

!--------------------------------------------------------------------------------------------------

  SUBROUTINE call_isccpsimulator(&
  kproma, kbdim, klev, jrow,     &
  p, ph, qm, at, skt)
    !
    USE mo_radiation_parameters,  ONLY: emsfc_lw => cemiss
    IMPLICIT NONE
    !
    !
    ! Input
    !
    INTEGER                        :: kproma, kbdim, klev, jrow ! Dimensions
    REAL(dp), DIMENSION(kbdim, klev)   :: p                  ! Pressure @ full levels [Pa]
    REAL(dp), DIMENSION(kbdim, klev+1) :: ph                 ! Pressure @ layer interfaces [Pa]
    REAL(dp), DIMENSION(kbdim, klev)   :: qm                 ! Water vapour mixing ratio [kg kg-1]
    REAL(dp), DIMENSION(kbdim, klev)   :: at                 ! Air temperature [K]
    REAL(dp), DIMENSION(kproma)        :: skt                ! Surface skin temperature [K] 
    !
    ! Input to ISCCP-Simulator to be set here
    INTEGER :: debug, debugcol, top_height, overlap
    INTEGER, DIMENSION(kproma)    :: sunlit
    REAL(dp), DIMENSION(kproma, klev) :: dem_s               !  10.5 micron longwave emissivity of
                                                             !   stratiform
    REAL(dp), DIMENSION(kproma, klev) :: dem_c               !  10.5 micron longwave emissivity of
                                                             !   convective

    ! Local
    REAL(dp), DIMENSION(kbdim, klev)   :: fc                 ! Convective cloud cover [0-1]
    REAL(dp), DIMENSION(kbdim, klev)   :: tauc               ! In-cloud optical thickness 
                                                             !  (convective clouds)
    REAL(dp), DIMENSION(kbdim, klev)   :: taus               ! In-cloud optical thickness 
                                                             !  (stratiform clouds)
    REAL(dp), DIMENSION(kproma, 7, 7)  :: fq_isccp
    REAL(dp), POINTER :: ct(:,:)
    REAL(dp), DIMENSION(kproma) :: pcldfra
    REAL(dp), DIMENSION(kproma) :: pcldtau
    REAL(dp), DIMENSION(kproma) :: pcldptop

    INTEGER i,itau,ipres

    ! Avoid invalid cloud fraction
    WHERE  ( isccp_cldtau3d(1:kproma,1:klev,jrow).LE.0._dp ) 
       isccp_f3d(1:kproma,1:klev,jrow) = 0._dp
    ENDWHERE

    ! Set general input to ISCCP-Simulator
    debug=0
    debugcol=0
    top_height=1 ! adjustment of cloud top height
    overlap=3    ! random-maximum
    
    WHERE ( isccp_sunlit(1:kproma,jrow) > 0._dp )
       sunlit(1:kproma) = 1
    ELSEWHERE
       sunlit(1:kproma) = 0
    ENDWHERE
    
    fc(1:kproma,1:klev)    = 0._dp  ! convective clouds = 0
    tauc(1:kproma,1:klev)  = 0._dp  ! convective clouds = 0
    dem_c(1:kproma,1:klev) = 0._dp  ! convective clouds = 0
    
    
    ! In-cloud optical thickness and emissivity
    WHERE ( isccp_f3d(1:kproma,1:klev,jrow) > 0._dp ) 
       taus(1:kproma,1:klev) = isccp_cldtau3d(1:kproma,1:klev,jrow)/isccp_f3d(1:kproma,1:klev,jrow)
       dem_s(1:kproma,1:klev) = isccp_cldemi3d(1:kproma,1:klev,jrow)
    ELSEWHERE
       taus(1:kproma,1:klev) = 0._dp
       dem_s(1:kproma,1:klev) = 0._dp
    ENDWHERE

    CALL isccp_cloud_types(&
         debug, debugcol, &
         kproma, &
         sunlit, &
         klev, ncol, &
         p(1:kproma,:), ph(1:kproma,:), &
         qm(1:kproma,:), &
         isccp_f3d(1:kproma,:,jrow), fc(1:kproma,:), &
         taus(1:kproma,:), tauc(1:kproma,:), &
         top_height, overlap, &
         tautab, invtau, &
         skt(1:kproma),&
         emsfc_lw,&
         at(1:kproma,:),&
         dem_s,&
         dem_c,&
         fq_isccp, &
         pcldfra, pcldptop, pcldtau)
    
    WHERE ( pcldfra(1:kproma) > 0._dp .AND. pcldtau(1:kproma) > 0._dp ) 
       isccp_cldfreq(1:kproma,jrow) = 1._dp
    ELSEWHERE
       isccp_cldfreq(1:kproma,jrow) = 0._dp
    END WHERE

    i=0
    DO itau=1,7
       DO ipres=1,7
          i=i+1
          ct => isccp_cldtypes(i)%ptr
          ct(1:kproma,jrow) = fq_isccp(1:kproma,itau,ipres) 
       ENDDO
    ENDDO

    isccp_cldtau(1:kproma,jrow)  = pcldtau(1:kproma)  * isccp_cldfreq(1:kproma,jrow)
    isccp_cldptop(1:kproma,jrow) = pcldptop(1:kproma) * isccp_cldfreq(1:kproma,jrow)
    isccp_cldfra(1:kproma,jrow)  = pcldfra(1:kproma)  

  END SUBROUTINE call_isccpsimulator

!--------------------------------------------------------------------------------------------------

      SUBROUTINE isccp_cloud_types(&
          debug,&
          debugcol,&
          npoints,&
          sunlit,&
          nlev,&
          ncol,&
          pfull,&
          phalf,&
          qv,&
          cc,&
          conv,&
          dtau_s,&
          dtau_c,&
          top_height,&
          overlap,&
          tautab,&
          invtau,&
          skt,&
          emsfc_lw,&
          at,&
          dem_s,&
          dem_c,&
          fq_isccp,&
          totalcldarea,&
          meanptop,&
          meantaucld &!,&
!          boxtau,&
!          boxptop&
     )
        
!$Id: isccp_cloud_types.f,v 3.4 2003/05/06 08:45:28 hadmw Exp $

! Copyright Steve Klein and Mark Webb 2002 - all rights reserved.
!
! This code is available without charge with the following conditions:
!
!  1. The code is available for scientific purposes and is not for 
!     commercial use.
!  2. Any improvements you make to the code should be made available 
!     to the to the authors for incorporation into a future release.
!  3. The code should not be used in any way that brings the authors 
!     or their employers into disrepute.

  USE mo_exception      ,ONLY: finish, message, message_text

      IMPLICIT NONE

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER ncolprint
      
!     -----
!     Input 
!     -----

      INTEGER npoints                   !  number of model points in the horizontal
      INTEGER nlev                      !  number of model levels in column
      INTEGER ncol                      !  number of subcolumns

      INTEGER sunlit(npoints)           !  1 for day points, 0 for night time


      REAL(dp) ::  pfull(npoints,nlev)          !  pressure of full model levels (Pascals)
                                        !  pfull(npoints,1)    is    top level of model
                                        !  pfull(npoints,nlev) is bottom level of model

      REAL(dp) ::  phalf(npoints,nlev+1)        !  pressure of half model levels (Pascals)
                                        !  phalf(npoints,1)    is    top       of model
                                        !  phalf(npoints,nlev+1) is the surface pressure

      REAL(dp) ::  qv(npoints,nlev)             !  water vapor specific humidity (kg vapor/ kg air)
                                        !         on full model levels

      REAL(dp) ::  cc(npoints,nlev)             !  input cloud cover in each model level (fraction) 
                                        !  NOTE:  This is the HORIZONTAL area of each
                                        !         grid box covered by clouds

      REAL(dp) ::  conv(npoints,nlev)   !  input convective cloud cover in each model level
                                        !   (fraction) 
                                        !  NOTE:  This is the HORIZONTAL area of each
                                        !         grid box covered by convective clouds

      REAL(dp) ::  dtau_s(npoints,nlev)         !  mean 0.67 micron optical depth of stratiform
                                        !  clouds in each model level
                                        !  NOTE:  this the cloud optical depth of only the
                                        !         cloudy part of the grid box, it is not weighted
                                        !         with the 0 cloud optical depth of the clear
                                        !         part of the grid box

      REAL(dp) ::  dtau_c(npoints,nlev)         !  mean 0.67 micron optical depth of convective
                                        !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.

      INTEGER overlap                   !  overlap type
                                        !  1=max
                                        !  2=rand
                                        !  3=max/rand

      INTEGER top_height                !  1 = adjust top height using both a computed
                                        !  infrared brightness temperature and the visible
                                        !  optical depth to adjust cloud top pressure. Note
                                        !  that this calculation is most appropriate to compare
                                        !  to ISCCP data during sunlit hours.
                                        !  2 = do not adjust top height, that is cloud top
                                        !  pressure is the actual cloud top pressure
                                        !  in the model
                                        !  3 = adjust top height using only the computed
                                        !  infrared brightness temperature. Note that this
                                        !  calculation is most appropriate to compare to ISCCP
                                        !  IR only algortihm (i.e. you can compare to nighttime
                                        !  ISCCP data with this option)

      REAL(dp) ::  tautab(0:255)                !  ISCCP table for converting count value to 
                                        !  optical thickness

      INTEGER invtau(-20:45000)         !  ISCCP table for converting optical thickness 
                                        !  to count value
!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL(dp) ::  skt(npoints)                 !  skin Temperature (K)
      REAL(dp) ::  emsfc_lw                     !  10.5 micron emissivity of surface (fraction)
      REAL(dp) ::  at(npoints,nlev)     !  temperature in each model level (K)
      REAL(dp) ::  dem_s(npoints,nlev)  !  10.5 micron longwave emissivity of stratiform
                                        !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.
      REAL(dp) ::  dem_c(npoints,nlev)  !  10.5 micron longwave emissivity of convective
                                        !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.
!     ------
!     Output
!     ------

      REAL(dp) ::  fq_isccp(npoints,7,7)        !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types

      REAL(dp) ::  totalcldarea(npoints)        !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  This should
                                        !  equal the sum over all entries of fq_isccp
        
        
      ! The following three means are averages over the cloudy areas only.  If no
      ! clouds are in grid box all three quantities should equal zero.  
                                        
      REAL(dp) ::  meanptop(npoints)            !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
                                        
      REAL(dp) ::  meantaucld(npoints)          !  mean optical thickness 
                                        !  linear averaging in albedo performed.
      
!      REAL(dp) ::  boxtau(npoints,ncol)         !  optical thickness in each column
      
!      REAL(dp) ::  boxptop(npoints,ncol)    !  cloud top pressure (mb) in each column
                                        
                                                                                                                        
!
!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------

      REAL(dp) ::  frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                                        ! Equivalent of BOX in original version, but
                                        ! indexed by column then row, rather than
                                        ! by row then column

      REAL(dp) ::  tca(npoints,0:nlev) ! total cloud cover in each model level (fraction)
                                        ! with extra layer of zeroes on top
                                        ! in this version this just contains the values input
                                        ! from cc but with an extra level
      REAL(dp) ::  cca(npoints,nlev)    ! convective cloud cover in each model level (fraction)
                                        ! from conv 

      REAL(dp) ::  threshold(npoints,ncol)   ! pointer to position in gridbox
      REAL(dp) ::  maxocc(npoints,ncol)      ! Flag for max overlapped conv cld
      REAL(dp) ::  maxosc(npoints,ncol)      ! Flag for max overlapped strat cld

      REAL(dp) ::  boxpos(npoints,ncol)      ! ordered pointer to position in gridbox

      REAL(dp) ::  threshold_min(npoints,ncol) ! minimum value to define range in with new 
                                               !  threshold is chosen

      REAL(dp) ::  dem(npoints,ncol),bb(npoints)     !  working variables for 10.5 micron longwave 
                                        !  emissivity in part of
                                        !  gridbox under consideration

      REAL(dp) ::  ran(npoints)                 ! vector of random numbers
      REAL(dp) ::  ptrop(npoints)
      REAL(dp) ::  attrop(npoints)
      REAL(dp) ::  attropmin (npoints)
      REAL(dp) ::  atmax(npoints)
      REAL(dp) ::  atmin(npoints)
      REAL(dp) ::  btcmin(npoints)
      REAL(dp) ::  transmax(npoints)

      INTEGER i,j,ilev,ibox,itrop(npoints)
      INTEGER ipres(npoints)
      INTEGER itau(npoints),ilev2
      INTEGER acc(nlev,ncol)
      INTEGER match(npoints,nlev-1)
      INTEGER nmatch(npoints)
      INTEGER levmatch(npoints,ncol)
      
      !variables needed for water vapor continuum absorption
      REAL(dp) :: fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
      REAL(dp) :: taumin(npoints)
      REAL(dp) :: dem_wv(npoints,nlev), wtmair, wtmh20, Navo, grav, pstd, t0
      REAL(dp) :: press(npoints), dpress(npoints), atmden(npoints)
      REAL(dp) :: rvh20(npoints), wk(npoints), rhoave(npoints)
      REAL(dp) :: rh20s(npoints), rfrgn(npoints)
      REAL(dp) :: tmpexp(npoints),tauwv(npoints)
      
      character*1 cchar(6),cchar_realtops(6)
      integer icycle
      REAL(dp) ::  tau(npoints,ncol)
      LOGICAL box_cloudy(npoints,ncol)
      REAL(dp) ::  tb(npoints,ncol)
      REAL(dp) ::  ptop(npoints,ncol)
      REAL(dp) ::  emcld(npoints,ncol)
      REAL(dp) ::  fluxtop(npoints,ncol)
      REAL(dp) ::  trans_layers_above(npoints,ncol)
      REAL(dp) :: isccp_taumin,fluxtopinit(npoints),tauir(npoints)
      REAL(dp) :: meanalbedocld(npoints) 
      REAL(dp) ::  albedocld(npoints,ncol)
      REAL(dp) :: boxarea
      integer debug       ! set to non-zero value to print out inputs
                          ! with step debug
      integer debugcol    ! set to non-zero value to print out column
                          ! decomposition with step debugcol

      integer index1(npoints),num1,jj
      real rec2p13,tauchk

      character*10 ftn09
      
      DATA isccp_taumin / 0.3_dp /
      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

      tauchk = -1._dp*log(0.9999999_dp)
      rec2p13=1._dp/2.13_dp

      ncolprint=0

      if ( debug.ne.0 ) then
          j=1
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write(6,'(a10)') 'debug='
          write(6,'(8I10)') debug
          write(6,'(a10)') 'debugcol='
          write(6,'(8I10)') debugcol
          write(6,'(a10)') 'npoints='
          write(6,'(8I10)') npoints
          write(6,'(a10)') 'nlev='
          write(6,'(8I10)') nlev
          write(6,'(a10)') 'ncol='
          write(6,'(8I10)') ncol
          write(6,'(a10)') 'top_height='
          write(6,'(8I10)') top_height
          write(6,'(a10)') 'overlap='
          write(6,'(8I10)') overlap
          write(6,'(a10)') 'emsfc_lw='
          write(6,'(8f10.2)') emsfc_lw
          write(6,'(a10)') 'tautab='
          write(6,'(8f10.2)') tautab
          write(6,'(a10)') 'invtau(1:100)='
          write(6,'(8i10)') (invtau(i),i=1,100)
          do j=1,npoints,debug
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write(6,'(a10)') 'sunlit='
          write(6,'(8I10)') sunlit(j)
          write(6,'(a10)') 'pfull='
          write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
          write(6,'(a10)') 'phalf='
          write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
          write(6,'(a10)') 'qv='
          write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
          write(6,'(a10)') 'cc='
          write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
          write(6,'(a10)') 'conv='
          write(6,'(8f10.2)') (conv(j,i),i=1,nlev)
          write(6,'(a10)') 'dtau_s='
          write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
          write(6,'(a10)') 'dtau_c='
          write(6,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
          write(6,'(a10)') 'skt='
          write(6,'(8f10.2)') skt(j)
          write(6,'(a10)') 'at='
          write(6,'(8f10.2)') (at(j,i),i=1,nlev)
          write(6,'(a10)') 'dem_s='
          write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
          write(6,'(a10)') 'dem_c='
          write(6,'(8f10.2)') (dem_c(j,i),i=1,nlev)
          enddo
      endif

!     ---------------------------------------------------!

!     assign 2d tca array using 1d input array cc

      do j=1,npoints
        tca(j,0)=0._dp
      enddo
  
      do ilev=1,nlev
        do j=1,npoints
          tca(j,ilev)=cc(j,ilev)
        enddo
      enddo

!     assign 2d cca array using 1d input array conv

      do ilev=1,nlev
!jq        do ibox=1,ncol
          do j=1,npoints
             cca(j,ilev)=conv(j,ilev)
          enddo
!jq        enddo
      enddo

      if (ncolprint.ne.0) then
        do j=1,npoints,1000
        write(6,'(a10)') 'j='
        write(6,'(8I10)') j
        write (6,'(a)') 'tca_pp_rev:'
        write (6,'(8f5.2)') &
        ((tca(j,ilev)),&
           ilev=1,nlev)

        write (6,'(a)') 'cca_pp_rev:'
        write (6,'(8f5.2)') &
        ((cca(j,ilev),ibox=1,ncolprint),ilev=1,nlev)
        enddo
      endif

      if (top_height .eq. 1 .or. top_height .eq. 3) then 

      do j=1,npoints 
          ptrop(j)=5000._dp
          atmin(j) = 400._dp
          attropmin(j) = 400._dp
          atmax(j) = 0._dp
          attrop(j) = 120._dp
          itrop(j) = 1
      enddo 

      do 12 ilev=1,nlev
        do j=1,npoints 
         if (pfull(j,ilev) .lt. 40000._dp .and.&
               pfull(j,ilev) .gt.  5000._dp .and.&
               at(j,ilev) .lt. attropmin(j)) then
                ptrop(j) = pfull(j,ilev)
                attropmin(j) = at(j,ilev)
                attrop(j) = attropmin(j)
                itrop(j)=ilev
           end if
           if (at(j,ilev) .gt. atmax(j)) atmax(j)=at(j,ilev)
           if (at(j,ilev) .lt. atmin(j)) atmin(j)=at(j,ilev)
        enddo
12    continue

      end if

!     -----------------------------------------------------!

!     ---------------------------------------------------!

      do 13 ilev=1,nlev
        num1=0
        do j=1,npoints
            if (cc(j,ilev) .lt. 0._dp .or. cc(j,ilev) .gt. 1._dp) then
                num1=num1+1
                index1(num1)=j
            end if
        enddo
       IF (num1 > 0) THEN
         WRITE(message_text,'(a)') ' error = cloud fraction less than zero'
         CALL message('',message_text)
         WRITE(message_text,'(a)') ' or '
         CALL message('',message_text)
         WRITE(message_text,'(a)') ' error = cloud fraction greater than 1'
         CALL message('',message_text)
         DO jj=1,num1   
               j=index1(jj)
               WRITE(message_text,'(a,i4,a,i4,a,e30.19)') &
                'value at point ',j,',level ',ilev,' is ',cc(j,ilev)
         ENDDO
         CALL finish('isccp_cloud_types',' Wrong cloud fraction')
       ENDIF

        num1=0
        do j=1,npoints
            if (conv(j,ilev) .lt. 0._dp .or. conv(j,ilev) .gt. 1._dp) then
                num1=num1+1
                index1(num1)=j
            end if
        enddo
       IF (num1 > 0) THEN
         WRITE(message_text,'(a)') ' error = convective cloud fraction less than zero'
         CALL message('',message_text)
         WRITE(message_text,'(a)') ' or '
         CALL message('',message_text)
         WRITE(message_text,'(a)') ' error = convective cloud fraction greater than 1'
         CALL message('',message_text)
         DO jj=1,num1   
               j=index1(jj)
               WRITE(message_text,'(a,i4,a,i4,a,e30.19)') &
                'value at point ',j,',level ',ilev,' is ',conv(j,ilev)
         ENDDO
         CALL finish('isccp_cloud_types',' Wrong convective cloud fraction')
       ENDIF

        num1=0
        do j=1,npoints
            if (dtau_s(j,ilev) .lt. 0._dp) then
                num1=num1+1
                index1(num1)=j
            end if
        enddo
       IF (num1 > 0) THEN
         WRITE(message_text,'(a)') ' error = stratiform cloud opt. depth less than zero'
         CALL message('',message_text)
         DO jj=1,num1   
               j=index1(jj)
               WRITE(message_text,'(a,i4,a,i4,a,e30.19)') &
                'value at point ',j,',level ',ilev,' is ',dtau_s(j,ilev)
         ENDDO  
         CALL finish('isccp_cloud_types',' Wrong stratiform cloud opt. depth')
       ENDIF

        num1=0
        do j=1,npoints
            if (dtau_c(j,ilev) .lt. 0._dp) then
                num1=num1+1
                index1(num1)=j
            end if
        enddo
       IF (num1 > 0) THEN
         WRITE(message_text,'(a)') ' error = convective cloud opt. depth less than zero'
         CALL message('',message_text)
         DO jj=1,num1
               j=index1(jj)
               WRITE(message_text,'(a,i4,a,i4,a,e30.19)') &
                'value at point ',j,',level ',ilev,' is ',dtau_c(j,ilev)
         ENDDO
         CALL finish('isccp_cloud_types',' Wrong convective cloud opt. depth')
       ENDIF

        num1=0
        do j=1,npoints
            if (dem_s(j,ilev) .lt. 0._dp .or. dem_s(j,ilev) .gt. 1._dp) then
                num1=num1+1
                index1(num1)=j
            end if
        enddo
       IF (num1 > 0) THEN
         WRITE(message_text,'(a)') ' error = stratiform cloud emissivity less than zero'
         CALL message('',message_text)
         WRITE(message_text,'(a)') ' or '
         CALL message('',message_text)
         WRITE(message_text,'(a)') ' error = stratiform cloud emissivity greater than 1'
         CALL message('',message_text)
         DO jj=1,num1   
               j=index1(jj)
               WRITE(message_text,'(a,i4,a,i4,a,e30.19)') &
                'value at point ',j,',level ',ilev,' is ',dem_s(j,ilev)
         ENDDO
         CALL finish('isccp_cloud_types',' Wrong stratiform cloud emissivity')
       ENDIF

        num1=0
        do j=1,npoints
            if (dem_c(j,ilev) .lt. 0._dp .or. dem_c(j,ilev) .gt. 1._dp) then
                num1=num1+1
                index1(num1)=j
            end if
        enddo
       IF (num1 > 0) THEN
         WRITE(message_text,'(a)') ' error = convective cloud emissivity less than zero'
         CALL message('',message_text)
         WRITE(message_text,'(a)') ' or '
         CALL message('',message_text)
         WRITE(message_text,'(a)') ' error = convective cloud emissivity greater than 1'
         CALL message('',message_text)
         DO jj=1,num1   
               j=index1(jj)
               WRITE(message_text,'(a,i4,a,i4,a,e30.19)') &
                'value at point ',j,',level ',ilev,' is ',dem_c(j,ilev)
         ENDDO
         CALL finish('isccp_cloud_types',' Wrong convective cloud emissivity')
       ENDIF
13    continue


      do ibox=1,ncol
        do j=1,npoints 
          boxpos(j,ibox)=(ibox-.5_dp)/ncol
        enddo
      enddo

!     ---------------------------------------------------!
!     Initialise working variables
!     ---------------------------------------------------!

!     Initialised frac_out to zero

      do ibox=1,ncol
        do ilev=1,nlev
          do j=1,npoints
            frac_out(j,ibox,ilev)=0.0_dp
          enddo
        enddo
      enddo

      if (ncolprint.ne.0) then
        write (6,'(a)') 'frac_out_pp_rev:'
          do j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(8f5.2)') & 
          ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)

          enddo
        write (6,'(a)') 'ncol:'
        write (6,'(I3)') ncol
      endif
      if (ncolprint.ne.0) then
        write (6,'(a)') 'last_frac_pp:'
          do j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(8f5.2)') (tca(j,0))
          enddo
      endif

!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
!     frac_out is the array that contains the information 
!     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
!     convective cloud
      
      !loop over vertical levels
      DO 200 ilev = 1,nlev
                                  
!     Initialise threshold

        IF (ilev.eq.1) then
            ! If max overlap 
            IF (overlap.eq.1) then
              ! select pixels spread evenly
              ! across the gridbox
              DO ibox=1,ncol
                do j=1,npoints
                  threshold(j,ibox)=boxpos(j,ibox)
                enddo
              enddo
            ELSE
              DO ibox=1,ncol
                call random_uniform(ran)
                ! select random pixels from the non-convective
                ! part the gridbox ( some will be converted into
                ! convective pixels below )
                do j=1,npoints
                  threshold(j,ibox)= &
                 cca(j,ilev)+(1-cca(j,ilev))*ran(j)
                enddo
              enddo
            ENDIF
            IF (ncolprint.ne.0) then
              write (6,'(a)') 'threshold_nsf2:'
                do j=1,npoints,1000
                write(6,'(a10)') 'j='
                write(6,'(8I10)') j
                write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
                enddo
            ENDIF
        ENDIF

        IF (ncolprint.ne.0) then
            write (6,'(a)') 'ilev:'
            write (6,'(I2)') ilev
        ENDIF

        DO ibox=1,ncol

          ! All versions
          do j=1,npoints
            if (boxpos(j,ibox).le.cca(j,ilev)) then
              maxocc(j,ibox) = 1
            else
              maxocc(j,ibox) = 0
            end if
          enddo

          ! Max overlap
          if (overlap.eq.1) then 
            do j=1,npoints
              threshold_min(j,ibox)=cca(j,ilev)
              maxosc(j,ibox)=1
            enddo
          endif

          ! Random overlap
          if (overlap.eq.2) then 
            do j=1,npoints
              threshold_min(j,ibox)=cca(j,ilev)
              maxosc(j,ibox)=0
            enddo
          endif

          ! Max/Random overlap
          if (overlap.eq.3) then 
            do j=1,npoints
              threshold_min(j,ibox)=max(cca(j,ilev), &
               min(tca(j,ilev-1),tca(j,ilev)))
              if (threshold(j,ibox) &
               .lt.min(tca(j,ilev-1),tca(j,ilev)) &
               .and.(threshold(j,ibox).gt.cca(j,ilev))) then
                   maxosc(j,ibox)= 1
              else
                   maxosc(j,ibox)= 0
              end if
            enddo
          endif
    
          ! Reset threshold 

          call random_uniform(ran)
           
          do j=1,npoints
            threshold(j,ibox)= &
              !if max overlapped conv cloud
             maxocc(j,ibox) * ( &                                       
                 boxpos(j,ibox)  &                                             
             ) +                  &                                    
              !else
             (1-maxocc(j,ibox)) * ( &                                   
                  !if max overlapped strat cloud
                 (maxosc(j,ibox)) * ( &                                 
                      !threshold=boxpos 
                     threshold(j,ibox) &                                       
                 ) +                    &                              
                  !else
                 (1-maxosc(j,ibox)) * (  &                             
                      !threshold_min=random[thrmin,1]
                     threshold_min(j,ibox)+&
                       (1-threshold_min(j,ibox))*ran(j)  &
                ) &
             )
          enddo

        ENDDO ! ibox

!          Fill frac_out with 1's where tca is greater than the threshold

           DO ibox=1,ncol
             do j=1,npoints 
               if (tca(j,ilev).gt.threshold(j,ibox)) then
               frac_out(j,ibox,ilev)=1._dp
               else
               frac_out(j,ibox,ilev)=0._dp
               end if               
             enddo
           ENDDO

!          Code to partition boxes into startiform and convective parts
!          goes here

           DO ibox=1,ncol
             do j=1,npoints 
                if (threshold(j,ibox).le.cca(j,ilev)) then
                    ! = 2 IF threshold le cca(j)
                    frac_out(j,ibox,ilev) = 2 
                else
                    ! = the same IF NOT threshold le cca(j) 
                    frac_out(j,ibox,ilev) = frac_out(j,ibox,ilev)
                end if
             enddo
           ENDDO

!         Set last_frac to tca at this level, so as to be tca 
!         from last level next time round

          if (ncolprint.ne.0) then

            do j=1,npoints ,1000
            write(6,'(a10)') 'j='
            write(6,'(8I10)') j
            write (6,'(a)') 'last_frac:'
            write (6,'(8f5.2)') (tca(j,ilev-1))
    
            write (6,'(a)') 'cca:'
            write (6,'(8f5.2)') (cca(j,ilev),ibox=1,ncolprint)
    
            write (6,'(a)') 'max_overlap_cc:'
            write (6,'(8f5.2)') (maxocc(j,ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'max_overlap_sc:'
            write (6,'(8f5.2)') (maxosc(j,ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'threshold_min_nsf2:'
            write (6,'(8f5.2)') (threshold_min(j,ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'threshold_nsf2:'
            write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'frac_out_pp_rev:'
            write (6,'(8f5.2)') &
            ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
            enddo
          endif

200   CONTINUE    !loop over nlev

!
!     ---------------------------------------------------!

      
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !initialize tau and albedocld to zero
      do 15 ibox=1,ncol
        do j=1,npoints 
            tau(j,ibox)=0._dp
            albedocld(j,ibox)=0._dp
!            boxtau(j,ibox)=0._dp
!            boxptop(j,ibox)=0._dp
            box_cloudy(j,ibox)=.false.
        enddo
15    continue

      !compute total cloud optical depth for each column     
      do ilev=1,nlev
            !increment tau for each of the boxes
            do ibox=1,ncol
              do j=1,npoints 
                 if (frac_out(j,ibox,ilev).eq.1) then
                        tau(j,ibox)=tau(j,ibox) &
                          + dtau_s(j,ilev)
                 endif
                 if (frac_out(j,ibox,ilev).eq.2) then
                        tau(j,ibox)=tau(j,ibox) &
                          + dtau_c(j,ilev)
                 end if
              enddo
            enddo ! ibox
      enddo ! ilev
          if (ncolprint.ne.0) then

              do j=1,npoints ,1000
                write(6,'(a10)') 'j='
                write(6,'(8I10)') j
                write(6,'(i2,1X,8(f7.2,1X))') &
               ilev, &
               (tau(j,ibox),ibox=1,ncolprint)
              enddo
          endif 
!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
!             sky versions of these quantities.

      if (top_height .eq. 1 .or. top_height .eq. 3) then


        !----------------------------------------------------------------------
        !    
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1, 
        !or 10.47 microns 
        wtmair = 28.9644_dp
        wtmh20 = 18.01534_dp
        Navo = 6.023E+23_dp
        grav = 9.806650E+02_dp
        pstd = 1.013250E+06_dp
        t0 = 296._dp
        if (ncolprint .ne. 0) &
              write(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
        do 125 ilev=1,nlev
          do j=1,npoints 
               !press and dpress are dyne/cm2 = Pascals *10
               press(j) = pfull(j,ilev)*10.
               dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10
               !atmden = g/cm2 = kg/m2 / 10 
               atmden(j) = dpress(j)/grav
               rvh20(j) = qv(j,ilev)*wtmair/wtmh20
               wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
               rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
               rh20s(j) = rvh20(j)*rhoave(j)
               rfrgn(j) = rhoave(j)-rh20s(j)
               tmpexp(j) = exp(-0.02_dp*(at(j,ilev)-t0))
               tauwv(j) = wk(j)*1.e-20_dp*( &
                (0.0224697_dp*rh20s(j)*tmpexp(j)) + &
                     (3.41817e-7_dp*rfrgn(j)) )*0.98_dp
               dem_wv(j,ilev) = 1._dp - exp( -1._dp * tauwv(j))
          enddo
               if (ncolprint .ne. 0) then
               do j=1,npoints ,1000
               write(6,'(a10)') 'j='
               write(6,'(8I10)') j
               write(6,'(i2,1X,3(f8.3,3X))') ilev,&
                qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))/(grav/100._dp),&
                tauwv(j),dem_wv(j,ilev)
               enddo
               endif
125     continue

        !initialize variables
        do j=1,npoints 
          fluxtop_clrsky(j) = 0._dp
          trans_layers_above_clrsky(j)=1._dp
        enddo

        do ilev=1,nlev
          do j=1,npoints 
 
            ! Black body emission at temperature of the layer

                bb(j)=1 / ( exp(1307.27_dp/at(j,ilev)) - 1._dp )
                !bb(j)= 5.67e-8*at(j,ilev)**4

                ! increase TOA flux by flux emitted from layer
                ! times total transmittance in layers above

                fluxtop_clrsky(j) = fluxtop_clrsky(j) &
                 + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j) 
            
                ! update trans_layers_above with transmissivity
                ! from this layer for next time around loop

                trans_layers_above_clrsky(j)= &
                 trans_layers_above_clrsky(j)*(1._dp-dem_wv(j,ilev))
                   

          enddo   
            if (ncolprint.ne.0) then
             do j=1,npoints ,1000
              write(6,'(a10)') 'j='
              write(6,'(8I10)') j
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
    
              write (6,'(a)') &
             'emiss_layer,100.*bb(j),100.*f,total_trans:'
              write (6,'(4(f7.2,1X))') dem_wv(j,ilev),100._dp*bb(j),&
                  100._dp*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
             enddo   
            endif

        enddo   !loop over level
        
        do j=1,npoints 
          !add in surface emission
          bb(j)=1/( exp(1307.27_dp/skt(j)) - 1. )
          !bb(j)=5.67e-8*skt(j)**4

          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j) & 
          * trans_layers_above_clrsky(j)
        enddo

        if (ncolprint.ne.0) then
        do j=1,npoints ,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'

          write (6,'(a)') 'emsfc,100.*bb(j),100.*f,total_trans:'
          write (6,'(4(f7.2,1X))') emsfc_lw,100._dp*bb(j), &
           100._dp*fluxtop_clrsky(j), &
            trans_layers_above_clrsky(j)
        enddo
        endif
    

        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------



        if (ncolprint.ne.0) then

        do j=1,npoints ,1000
            write(6,'(a10)') 'j='
            write(6,'(8I10)') j
            write (6,'(a)') 'ts:'
            write (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)
    
            write (6,'(a)') 'ta_rev:'
            write (6,'(8f7.2)') &
            ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)

        enddo
        endif 
        !loop over columns 
        do ibox=1,ncol
          do j=1,npoints
            fluxtop(j,ibox)=0.
            trans_layers_above(j,ibox)=1.
          enddo
        enddo

        do ilev=1,nlev
              do j=1,npoints 
                ! Black body emission at temperature of the layer

                bb(j)=1 / ( exp(1307.27_dp/at(j,ilev)) - 1. )
                !bb(j)= 5.67e-8*at(j,ilev)**4
              enddo

            do ibox=1,ncol
              do j=1,npoints 

                ! emissivity for point in this layer
                if (frac_out(j,ibox,ilev).eq.1) then
                dem(j,ibox)= 1._dp - & 
                  ( (1._dp - dem_wv(j,ilev)) * (1._dp -  dem_s(j,ilev)) )
                else if (frac_out(j,ibox,ilev).eq.2) then
                dem(j,ibox)= 1._dp - &
                  ( (1._dp - dem_wv(j,ilev)) * (1._dp -  dem_c(j,ilev)) )
                else
                dem(j,ibox)=  dem_wv(j,ilev)
                end if
                

                ! increase TOA flux by flux emitted from layer
                ! times total transmittance in layers above

                fluxtop(j,ibox) = fluxtop(j,ibox) &
                 + dem(j,ibox) * bb(j) &
                 * trans_layers_above(j,ibox) 
            
                ! update trans_layers_above with transmissivity
                ! from this layer for next time around loop

                trans_layers_above(j,ibox)= &
                 trans_layers_above(j,ibox)*(1._dp-dem(j,ibox))

              enddo ! j
            enddo ! ibox

            if (ncolprint.ne.0) then
              do j=1,npoints,1000
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
    
              write(6,'(a10)') 'j='
              write(6,'(8I10)') j
              write (6,'(a)') 'emiss_layer:'
              write (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
        
              write (6,'(a)') '100.*bb(j):'
              write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
        
              write (6,'(a)') '100.*f:'
              write (6,'(8f7.2)') &
              (100.*fluxtop(j,ibox),ibox=1,ncolprint)
        
              write (6,'(a)') 'total_trans:'
              write (6,'(8f7.2)') &
               (trans_layers_above(j,ibox),ibox=1,ncolprint)
              enddo
          endif

        enddo ! ilev


          do j=1,npoints 
            !add in surface emission
            bb(j)=1/( exp(1307.27_dp/skt(j)) - 1._dp )
            !bb(j)=5.67e-8*skt(j)**4
          end do

        do ibox=1,ncol
          do j=1,npoints 

            !add in surface emission

            fluxtop(j,ibox) = fluxtop(j,ibox)  &
              + emsfc_lw * bb(j) & 
              * trans_layers_above(j,ibox) 
            
          end do
        end do

        if (ncolprint.ne.0) then

          do j=1,npoints ,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'

          write (6,'(a)') 'emiss_layer:'
          write (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') '100.*bb(j):'
          write (6,'(8f7.2)') (100._dp*bb(j),ibox=1,ncolprint)
    
          write (6,'(a)') '100.*f:'
          write (6,'(8f7.2)') (100._dp*fluxtop(j,ibox),ibox=1,ncolprint)
          end do
        endif
    
        !now that you have the top of atmosphere radiance account
        !for ISCCP procedures to determine cloud top temperature

        !account for partially transmitting cloud recompute flux 
        !ISCCP would see assuming a single layer cloud
        !note choice here of 2.13, as it is primarily ice
        !clouds which have partial emissivity and need the 
        !adjustment performed in this section
        !
        !If it turns out that the cloud brightness temperature
        !is greater than 260K, then the liquid cloud conversion
        !factor of 2.56 is used.
        !
        !Note that this is discussed on pages 85-87 of 
        !the ISCCP D level documentation (Rossow et al. 1996)
           
          do j=1,npoints  
            !compute minimum brightness temperature and optical depth
            btcmin(j) = 1._dp /  ( exp(1307.27_dp/(attrop(j)-5._dp)) - 1._dp ) 
          enddo 
        do ibox=1,ncol
          do j=1,npoints  
            transmax(j) = (fluxtop(j,ibox)-btcmin(j)) &
                     /(fluxtop_clrsky(j)-btcmin(j))
            !note that the initial setting of tauir(j) is needed so that
            !tauir(j) has a realistic value should the next if block be
            !bypassed
            tauir(j) = tau(j,ibox) * rec2p13
            taumin(j) = -1._dp * log(max(min(transmax(j),0.9999999_dp),0.001_dp))

          enddo 

          if (top_height .eq. 1) then
            do j=1,npoints  
              if (transmax(j) .gt. 0.001_dp .and. &
               transmax(j) .le. 0.9999999_dp) then
                fluxtopinit(j) = fluxtop(j,ibox)
                tauir(j) = tau(j,ibox) *rec2p13
              endif
            enddo
            do icycle=1,2
              do j=1,npoints  
                if (tau(j,ibox) .gt. (tauchk            )) then 
                if (transmax(j) .gt. 0.001_dp .and. & 
                 transmax(j) .le. 0.9999999_dp) then
                  emcld(j,ibox) = 1._dp - exp(-1._dp * tauir(j)  )
                  fluxtop(j,ibox) = fluxtopinit(j) -   &
                   ((1._dp-emcld(j,ibox))*fluxtop_clrsky(j))
                  fluxtop(j,ibox)=max(1.E-06_dp, &
                   (fluxtop(j,ibox)/emcld(j,ibox)))
                  tb(j,ibox)= 1307.27_dp &
                   / (log(1._dp + (1._dp/fluxtop(j,ibox))))
                  if (tb(j,ibox) .gt. 260._dp) then
                    tauir(j) = tau(j,ibox) / 2.56_dp
                  end if                         
                end if
                end if
              enddo
            enddo
                
          endif
        
          do j=1,npoints
            if (tau(j,ibox) .gt. (tauchk            )) then 
                !cloudy box
                tb(j,ibox)= 1307.27_dp/ (log(1._dp + (1._dp/fluxtop(j,ibox))))
                if (top_height.eq.1.and.tauir(j).lt.taumin(j)) then
                         tb(j,ibox) = attrop(j) - 5. 
                         tau(j,ibox) = 2.13_dp*taumin(j)
                end if
            else
                !clear sky brightness temperature
                tb(j,ibox) = 1307.27_dp/(log(1._dp+(1._dp/fluxtop_clrsky(j))))
            end if
          enddo ! j
        enddo ! ibox

        if (ncolprint.ne.0) then

          do j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j

          write (6,'(a)') 'attrop:'
          write (6,'(8f7.2)') (attrop(j))
    
          write (6,'(a)') 'btcmin:'
          write (6,'(8f7.2)') (btcmin(j))
    
          write (6,'(a)') 'fluxtop_clrsky*100:'
          write (6,'(8f7.2)') & 
           (100.*fluxtop_clrsky(j))

          write (6,'(a)') '100.*f_adj:'
          write (6,'(8f7.2)') (100._dp*fluxtop(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'transmax:'
          write (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'tau:'
          write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'emcld:'
          write (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)') &
          (trans_layers_above(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_emiss:'
          write (6,'(8f7.2)') &
          (1.0-trans_layers_above(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)') &
          (trans_layers_above(j,ibox),ibox=1,ncolprint)
    
          write (6,'(a)') 'ppout:'
          write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
          enddo ! j
        endif

      end if

!     ---------------------------------------------------!

!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!

      !compute cloud top pressure
      do 30 ibox=1,ncol
        !segregate according to optical thickness
        if (top_height .eq. 1 .or. top_height .eq. 3) then  
          !find level whose temperature
          !most closely matches brightness temperature
          do j=1,npoints 
            nmatch(j)=0
          enddo
          do 29 ilev=1,nlev-1
            !cdir nodep
            do j=1,npoints 
              if ((at(j,ilev)   .ge. tb(j,ibox) .and. &
               at(j,ilev+1) .lt. tb(j,ibox)) .or. &
               (at(j,ilev) .le. tb(j,ibox) .and. &
               at(j,ilev+1) .gt. tb(j,ibox))) then 
   
                nmatch(j)=nmatch(j)+1
                if(abs(at(j,ilev)-tb(j,ibox)) .lt. &
                 abs(at(j,ilev+1)-tb(j,ibox))) then
                  match(j,nmatch(j))=ilev
                else
                  match(j,nmatch(j))=ilev+1
                end if
              end if                        
            enddo
29        continue

          do j=1,npoints 
            if (nmatch(j) .ge. 1) then
              ptop(j,ibox)=pfull(j,match(j,nmatch(j)))
              levmatch(j,ibox)=match(j,nmatch(j))   
            else
              if (tb(j,ibox) .lt. atmin(j)) then
                ptop(j,ibox)=ptrop(j)
                levmatch(j,ibox)=itrop(j)
              end if
              if (tb(j,ibox) .gt. atmax(j)) then
                ptop(j,ibox)=pfull(j,nlev)
                levmatch(j,ibox)=nlev
              end if                                
            end if
          enddo ! j

        else ! if (top_height .eq. 1 .or. top_height .eq. 3) 
 
          do j=1,npoints     
            ptop(j,ibox)=0._dp
          enddo
          do ilev=1,nlev
            do j=1,npoints     
              if ((ptop(j,ibox) .eq. 0._dp ) &
                .and.(frac_out(j,ibox,ilev) .ne. 0)) then
                ptop(j,ibox)=pfull(j,ilev)
                levmatch(j,ibox)=ilev
              end if
            end do
          end do
        end if                            
          
        do j=1,npoints
          if (tau(j,ibox) .le. (tauchk            )) then
            ptop(j,ibox)=0._dp
            levmatch(j,ibox)=0      
          endif 
        enddo

30    continue
              
!
!
!     ---------------------------------------------------!


!     
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined, 
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy 
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.  
!

      !compute isccp frequencies

      !reset frequencies
      do 37 ilev=1,7
      do 38 ilev2=1,7
        do j=1,npoints ! 
             fq_isccp(j,ilev,ilev2)=0._dp
        enddo
38    continue
37    continue

      !reset variables need for averaging cloud properties
      do j=1,npoints 
        totalcldarea(j) = 0._dp
        meanalbedocld(j) = 0._dp
        meanptop(j) = 0._dp
        meantaucld(j) = 0._dp
      enddo ! j

      boxarea = 1./real(ncol)
     
      do 39 ibox=1,ncol
        do j=1,npoints 

          if (tau(j,ibox) .gt. (tauchk            ) &
           .and. ptop(j,ibox) .gt. 0._dp) then
              box_cloudy(j,ibox)=.true.
          endif

          if (box_cloudy(j,ibox)) then

              ! totalcldarea always diagnosed day or night
              totalcldarea(j) = totalcldarea(j) + boxarea

              if (sunlit(j).eq.1) then

                ! tau diagnostics only with sunlight

!                boxtau(j,ibox) = tau(j,ibox)

                !convert optical thickness to albedo
                albedocld(j,ibox) &
                 =real(invtau(min(nint(100._dp*tau(j,ibox)),45000)))
            
                !contribute to averaging
                meanalbedocld(j) = meanalbedocld(j)  &
                 +albedocld(j,ibox)*boxarea

             endif

          endif

          if (sunlit(j).eq.1 .or. top_height .eq. 3) then 

            !convert ptop to millibars
            ptop(j,ibox)=ptop(j,ibox) / 100._dp
            
            !save for output cloud top pressure and optical thickness
!            boxptop(j,ibox) = ptop(j,ibox)
    
            if (box_cloudy(j,ibox)) then
            
              meanptop(j) = meanptop(j) + ptop(j,ibox)*boxarea

              !reset itau(j), ipres(j)
              itau(j) = 0
              ipres(j) = 0

              !determine optical depth category
              if (tau(j,ibox) .lt. isccp_taumin) then
                  itau(j)=1
              else if (tau(j,ibox) .ge. isccp_taumin &
               .and. tau(j,ibox) .lt. 1.3_dp) then
                itau(j)=2
              else if (tau(j,ibox) .ge. 1.3_dp & 
               .and. tau(j,ibox) .lt. 3.6_dp) then
                itau(j)=3
              else if (tau(j,ibox) .ge. 3.6_dp &
               .and. tau(j,ibox) .lt. 9.4_dp) then
                  itau(j)=4
              else if (tau(j,ibox) .ge. 9.4_dp &
               .and. tau(j,ibox) .lt. 23._dp) then
                  itau(j)=5
              else if (tau(j,ibox) .ge. 23._dp &
               .and. tau(j,ibox) .lt. 60._dp) then
                  itau(j)=6
              else if (tau(j,ibox) .ge. 60._dp) then
                  itau(j)=7
              end if

              !determine cloud top pressure category
              if (    ptop(j,ibox) .gt. 0._dp  &
               .and.ptop(j,ibox) .lt. 180._dp) then
                  ipres(j)=1
              else if(ptop(j,ibox) .ge. 180._dp &
               .and.ptop(j,ibox) .lt. 310._dp) then
                  ipres(j)=2
              else if(ptop(j,ibox) .ge. 310._dp &
               .and.ptop(j,ibox) .lt. 440._dp) then
                  ipres(j)=3
              else if(ptop(j,ibox) .ge. 440._dp &
               .and.ptop(j,ibox) .lt. 560._dp) then
                  ipres(j)=4
              else if(ptop(j,ibox) .ge. 560._dp &
               .and.ptop(j,ibox) .lt. 680._dp) then
                  ipres(j)=5
              else if(ptop(j,ibox) .ge. 680._dp &
               .and.ptop(j,ibox) .lt. 800._dp) then
                  ipres(j)=6
              else if(ptop(j,ibox) .ge. 800._dp) then
                  ipres(j)=7
              end if 

              !update frequencies
              if(ipres(j) .gt. 0.and.itau(j) .gt. 0) then
              fq_isccp(j,itau(j),ipres(j))= &
               fq_isccp(j,itau(j),ipres(j))+ boxarea
              end if

            end if

          end if
                       
        enddo ! j
39    continue
      
      !compute mean cloud properties
      do j=1,npoints 
        if (totalcldarea(j) .gt. 0._dp) then
          meanptop(j) = meanptop(j) / totalcldarea(j)
          if (sunlit(j).eq.1) then
            meanalbedocld(j) = meanalbedocld(j) / totalcldarea(j)
            meantaucld(j) = tautab(min(255,max(1,nint(meanalbedocld(j))))) 
          end if
        end if
      enddo ! j
!
!     ---------------------------------------------------!

!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
      if (debugcol.ne.0) then
!     
         do j=1,npoints,debugcol

            !produce character output
            do ilev=1,nlev
              do ibox=1,ncol
                   acc(ilev,ibox)=0
              enddo
            enddo

            do ilev=1,nlev
              do ibox=1,ncol
                   acc(ilev,ibox)=frac_out(j,ibox,ilev)*2
                   if (levmatch(j,ibox) .eq. ilev) &
                      acc(ilev,ibox)=acc(ilev,ibox)+1
              enddo
            enddo

             !print test

          write(ftn09,11) j
11        format('ftn09.',i4.4)
          open(9, FILE=ftn09, FORM='FORMATTED')

             write(9,'(a1)') ' '
                    write(9,'(10i5)') &
                       (ilev,ilev=5,nlev,5)
             write(9,'(a1)') ' '
             
             do ibox=1,ncol
               write(9,'(40(a1),1x,40(a1))') &
                (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev) &
                ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev) 
             end do
             close(9)

             if (ncolprint.ne.0) then
               write(6,'(a1)') ' '
                    write(6,'(a2,1X,5(a7,1X),a50)') & 
                       'ilev',&
                       'pfull','at',&
                       'cc*100','dem_s','dtau_s',&
                       'cchar'

               write (6,'(a)') 'skt(j):'
               write (6,'(8f7.2)') skt(j)
                                      
               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
              
               write (6,'(a)') 'tau:'
               write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
    
               write (6,'(a)') 'tb:'
               write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
    
               write (6,'(a)') 'ptop:'
               write (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
             endif 
    
        enddo
       
      end if 

      return
      END SUBROUTINE isccp_cloud_types 

END MODULE mo_isccpsimulator
