MODULE mo_p3_fields

  USE mo_kind,    ONLY: dp
  USE mo_control, ONLY : lcolumn
  USE mo_submodel, ONLY : print_value
  USE mo_exception, ONLY: message

  IMPLICIT NONE

  PUBLIC :: access_lookup_table, &
            calculate_lookup_table_indices, &
            calculate_lookup_table_indices_1d, &
            access_my_lookup_table, &
            calculate_my_lookup_table_indices, &
            get_ni_limits

  INTERFACE get_ni_limits
     MODULE PROCEDURE get_ni_limits_1d
     MODULE PROCEDURE get_ni_limits_2d
  END INTERFACE get_ni_limits

  INTERFACE calculate_my_lookup_table_indices
     MODULE PROCEDURE calculate_my_lookup_table_indices_1d
     MODULE PROCEDURE calculate_my_lookup_table_indices_2d
  END INTERFACE calculate_my_lookup_table_indices

  INTERFACE access_my_lookup_table
     MODULE PROCEDURE access_my_lookup_table_1d
     MODULE PROCEDURE access_my_lookup_table_2d
  END INTERFACE access_my_lookup_table

  LOGICAL                                  :: p3isalloc = .FALSE.
  INTEGER                                  :: n_iceCat  = 1

  INTEGER :: idt_qirim, idt_birim, idt_qihet, idt_qiliq, idt_qioliq, &
             idt_qsrc, idt_qprc, idt_nihet, idt_nihom, idt_ninuc, idt_nidet

! -----------------------------LOOKUP-TABLE UTIL--------------------------------

  ! ice microphysics lookup table array dimensions
  integer, parameter :: isize     = 20
  integer, parameter :: jsize     = 20
  integer, parameter :: densize   =  5
  integer, parameter :: rimsize   =  4
  integer, parameter :: rcollsize = 30

  ! number of ice microphysical quantities used ookup table
  integer, parameter :: tabsize   = 12

  ! number of ice-rain collection microphysical quantities used from lookup table
  integer, parameter :: colltabsize = 2

  REAL(dp) :: itab(densize,rimsize,isize,jsize,tabsize)                    !< ice lookup-table
  REAL(dp) :: itabcoll(densize,rimsize,isize,jsize,rcollsize,colltabsize)  !< ice lookup-table

! -----------------------------MY LOOKUP-TABLE UTIL--------------------------------

  ! ice microphysics lookup table array dimensions
  ! normalized q parameters (qi/nitot)
  INTEGER,     parameter :: qnsize  = 50
  REAL(dp),    parameter :: qnmin   = 4e-16
  REAL(dp),    parameter :: qnmax   = 4e-4

  ! rime fraction parameters
  INTEGER,     parameter :: frsize  = 4
  REAL(dp),    parameter :: frmin   = 0
  REAL(dp),    parameter :: frmax   = 1

  ! rime density parameters
  INTEGER,     parameter :: rhosize =  5
  REAL(dp),    parameter :: rhomin  = 100
  REAL(dp),    parameter :: rhomax  = 900

  ! number of ice microphysical quantities used ookup table
  INTEGER, parameter :: mtabsize = 12

  REAL(dp) :: mytab(frsize,rhosize,qnsize,mtabsize) !< ice lookup-table


! -----------------------------HARDCODED NAMELIST-------------------------------
! CONTROL
  INTEGER :: iprog                    = 1       !< 1: P3, 2: 2 cat 2 moment, 3: old
  LOGICAL :: lpiggy                   = .FALSE. !< calculate 2m in case of iprog==1
  LOGICAL :: l2moment                 = .FALSE.  !< only use m and N, in lookup calculations
  LOGICAL :: ldisableupdraftcond      = .FALSE.  !< disable the condition on updraft for ice formation
  LOGICAL :: lconstnc                 = .FALSE.  !< set cloud droplet number to constant
  LOGICAL :: lfull_diag               = .FALSE. !< skip storing massive diagnostic streams
  LOGICAL :: lspeedrun                = .FALSE. !< skip secondary diagnostic tracers tracers
! ACTIVATION
  LOGICAL :: lprescribeaerosols       = .FALSE.  !< prescribe aerosols using mo_preaero
  LOGICAL :: lprescribeaerosols_every = .FALSE. !< prescribe aerosols every timestep
  LOGICAL :: lprescribe_numbers       = .FALSE.  !< prescribe numbers (-> nxi(l)mult)
  LOGICAL :: lprescribe_numbers_every = .FALSE. !< ... every timestep
! SEDIMENTATION
  LOGICAL :: lsnow_sed                = .TRUE.  !< calculate snow sedimentation
  LOGICAL :: lfalling_ice             = .TRUE.  !< calculate falling ice in iprog==1,2
  LOGICAL :: lrain_sed                = .TRUE.  !< calculate rain
  REAL(dp):: vt_rdc                   = 1._dp   !< factor by which fall speeds are reduced
! CLOUD COVER
  REAL(dp):: ccclpwr                  = 1._dp   !< power in weight calculation for cloud cover.
                                                !< >1 : more weight on liquid phase
                                                !< <1 : more weight on ice phase
  REAL(dp):: ccsupci                  = 1._dp   !< Factor by which the koop supersat is multiplied
                                                !< to delay cover == 1
  LOGICAL :: lallow_si                = .FALSE. !< Allow ice to grow down to ice saturation

! PROCESS RATES
  LOGICAL :: lact                     = .TRUE.  !< switch for nucleation processes bf. sed.
  LOGICAL :: lctrl_frz_below_238K     = .TRUE. !< Contorols freezing below 238K 
  LOGICAL :: lctrl_het_mxphase_frz    = .TRUE.  !< Controls heterogeneous freezing
  LOGICAL :: levap_rain               = .TRUE.
  LOGICAL :: lmelting                 = .TRUE.
  LOGICAL :: lriming                  = .TRUE.
  LOGICAL :: lself_collection         = .TRUE.
  LOGICAL :: ladjustment              = .TRUE.  !< saturation adjustment
  INTEGER :: inumberadjustment        = 1       !< 0: no adjustment
                                                !< 1: adjust radius
                                                !< 2: adjust fixed number
  LOGICAL :: lcover_adjustment        = .TRUE.
  LOGICAL :: lmass_transport          = .TRUE.
  LOGICAL :: lsubsat_cnd_dep          = .TRUE.
  LOGICAL :: lcirrus                  = .TRUE.
  LOGICAL :: lwbf                     = .TRUE.
  LOGICAL :: lconvice                 = .TRUE.
  INTEGER :: isublimation             = 1
  INTEGER :: ihomfrz                  = 2

! PRESCRIBE (only if prescribed or const nc)
  REAL(dp):: nconstnc                 = 100e6_dp !< N_c if set to constant
  REAL(dp):: nximult                  = 1e4   !< zicnc/zxi if prescribed
  REAL(dp):: nxlmult                  = 1e14   !< zcdnc/zxl if prescribed
  REAL(dp):: csubw                    = 10._dp !< in cm/s
  LOGICAL :: ldiabheat                = .TRUE. !< setting used for numerics test.
                                               !< do not update t and q

! MASS MINIMUM FOR CLIPPING
  REAL(dp):: xismall                  = 1.e-12_dp   !< xi below this will be sublimated
  REAL(dp):: xlsmall                  = 1.e-12_dp   !< xl below this will be evaporated

! TUNING
  REAL(dp):: rcmax                    = 25e-6  !< maximum mean liquid radius
  REAL(dp):: adjsupsat                = 0.001  !< assumed supersaturation in growth equation for number adj.
  REAL(dp):: ccnislf                  = 80    !< enhancment factor for self-collection of ice
  REAL(dp):: ccqccol                  = 1.3    !< enhancment factor for riming
  REAL(dp):: ccftau                   = 1  !< delay factor for snow autoconversion
  REAL(dp):: ccrsnow                  = 1e-4 !< minimum ice radius for snow conversion

! SUBSTEPS
  INTEGER :: nmicro                   = 0     !< number of microphysics substeps (0 = dynamic)
  INTEGER :: nmicro_max               = 0    !< maximal number of micro substeps (0 = unlimited, -1 = all levels)
  INTEGER :: nsedi                    = 0     !< number of sedimentation substeps (0 = dynamic)
  INTEGER :: iintscheme               = 0     !< which integration scheme to use (0: explicit euler
                                              !                                   1: implicit euler)
  REAL(dp):: microprct                = 0.2   ! fraction of the time-step that does not meet CFL criterion

  CONTAINS

  SUBROUTINE p3_init

    REAL(dp) :: dum
    integer :: i,j,k,ii,jj
    
    CALL message('','initializing p3-lookuptable')
    !------------------------------------------------------------------------------------------!
    ! read in ice microphysics table
    open(unit=10,file='../p3data/lookup_ice_prog_de_v13.dat',status='old')

    do jj = 1,densize
       do ii = 1,rimsize
          do i = 1,isize
             do k = 1,jsize
                read(10,*) dum,dum,dum,dum,itab(jj,ii,i,k,1),itab(jj,ii,i,k,2),             &
                     itab(jj,ii,i,k,3),itab(jj,ii,i,k,4),itab(jj,ii,i,k,5),                 &
                     itab(jj,ii,i,k,6),itab(jj,ii,i,k,7),itab(jj,ii,i,k,8),dum,             &
                     itab(jj,ii,i,k,9),itab(jj,ii,i,k,10),itab(jj,ii,i,k,11),               &
                     itab(jj,ii,i,k,12)
             enddo
          enddo
    ! read in table for ice-rain collection
          do i = 1,isize
             do k = 1,jsize
                do j = 1,rcollsize
                   read(10,*) dum,dum,dum,dum,dum,itabcoll(jj,ii,i,k,j,1),                  &
                        itabcoll(jj,ii,i,k,j,2),dum
                   itabcoll(jj,ii,i,k,j,1) = dlog10(itabcoll(jj,ii,i,k,j,1))
                   itabcoll(jj,ii,i,k,j,2) = dlog10(itabcoll(jj,ii,i,k,j,2))
                enddo
             enddo
          enddo
       enddo
    enddo

    close(unit=10)
    CALL message('','p3-lookup table is initialized!')


  END SUBROUTINE p3_init

  SUBROUTINE my_p3_init

    INTEGER :: ifr, irhor, iqnorm

    CALL message('','initializing my p3-lookuptable')
    OPEN(unit=10,file='../p3data/my_p3_table.dat',status='old')

    DO ifr=1,frsize
       DO irhor=1,rhosize
          DO iqnorm=1,qnsize
             READ(10,*) mytab(ifr,irhor,iqnorm,1),mytab(ifr,irhor,iqnorm,2),mytab(ifr,irhor,iqnorm,3), &
                        mytab(ifr,irhor,iqnorm,4),mytab(ifr,irhor,iqnorm,5),mytab(ifr,irhor,iqnorm,6), &
                        mytab(ifr,irhor,iqnorm,7),mytab(ifr,irhor,iqnorm,8),mytab(ifr,irhor,iqnorm,9), &
                        mytab(ifr,irhor,iqnorm,10),mytab(ifr,irhor,iqnorm,11),mytab(ifr,irhor,iqnorm,12)
          END DO
       END DO
    END DO

    CLOSE(unit=10)
    CALL message('','my p3-lookup table is initialized!')

  END SUBROUTINE my_p3_init

  SUBROUTINE access_lookup_table(dum1,dumi,dum2,dumk,dum4,dumii,dum5,dumjj,index,proc)

   IMPLICIT NONE

   INTEGER,  INTENT(in)  :: index
   REAL(dp), INTENT(in)  :: dum1,dum2,dum4,dum5
   INTEGER,  INTENT(in)  :: dumjj,dumii,dumi,dumk
   REAL(dp), INTENT(out) :: proc

   REAL(dp) :: dproc1,dproc2,iproc2,iproc1,gproc1,tmp1,tmp2
   LOGICAL  :: ll1, ll2

  ! first interpolate for current rimed fraction index
    dproc1 = itab(dumjj,dumii,dumi,dumk,index)          &
                     +(dum1-real(dumi))                                      &
                     *(itab(dumjj,dumii,dumi+1,dumk,index)     &
                     -itab(dumjj,dumii,dumi,dumk,index))

    dproc2 = itab(dumjj,dumii,dumi,dumk+1,index)        &
                     +(dum1-REAL(dumi))                                      &
                     *(itab(dumjj,dumii,dumi+1,dumk+1,index)   &
                     -itab(dumjj,dumii,dumi,dumk+1,index))

    iproc1 = dproc1+(dum2-REAL(dumk))                          &
                     *(dproc2-dproc1)

  ! linearly interpolate to get process rates for rimed fraction index + 1
    dproc1 = itab(dumjj,dumii+1,dumi,dumk,index)        &
                     +(dum1-REAL(dumi))                                      &
                     *(itab(dumjj,dumii+1,dumi+1,dumk,index)   &
                     -itab(dumjj,dumii+1,dumi,dumk,index))

    dproc2 = itab(dumjj,dumii+1,dumi,dumk+1,index)      &
                     +(dum1-REAL(dumi))                                      &
                     *(itab(dumjj,dumii+1,dumi+1,dumk+1,index) &
                     -itab(dumjj,dumii+1,dumi,dumk+1,index))

    gproc1 = dproc1+(dum2-REAL(dumk))                          &
         *(dproc2-dproc1)
    tmp1   = iproc1+(dum4-REAL(dumii))                         &
                     *(gproc1-iproc1)

  ! get value at density index + 1

  ! first interpolate for current rimed fraction index
    dproc1 = itab(dumjj+1,dumii,dumi,dumk,index)        &
                     +(dum1-REAL(dumi))                                      &
                     *(itab(dumjj+1,dumii,dumi+1,dumk,index)   &
                     -itab(dumjj+1,dumii,dumi,dumk,index))

    dproc2 = itab(dumjj+1,dumii,dumi,dumk+1,index)      &
                     +(dum1-REAL(dumi))                                      &
                     *(itab(dumjj+1,dumii,dumi+1,dumk+1,index) &
                     -itab(dumjj+1,dumii,dumi,dumk+1,index))

    iproc1 = dproc1+(dum2-REAL(dumk))                          &
                     *(dproc2-dproc1)

  ! linearly interpolate to get process rates for rimed fraction index + 1
    dproc1 = itab(dumjj+1,dumii+1,dumi,dumk,index)      &
                     +(dum1-REAL(dumi))                                      &
                     *(itab(dumjj+1,dumii+1,dumi+1,dumk,index) &
                     -itab(dumjj+1,dumii+1,dumi,dumk,index))

    dproc2 = itab(dumjj+1,dumii+1,dumi,dumk+1,index)    &
                     +(dum1-REAL(dumi))                                      &
                     *(itab(dumjj+1,dumii+1,dumi+1,dumk+1,index)&
                     -itab(dumjj+1,dumii+1,dumi,dumk+1,index))

    gproc1 = dproc1+(dum2-REAL(dumk))                          &
                     *(dproc2-dproc1)
    tmp2   = iproc1+(dum4-REAL(dumii))                         &
                     *(gproc1-iproc1)

  ! get final process rate
    proc   = tmp1+(dum5-REAL(dumjj))*(tmp2-tmp1)

  END SUBROUTINE access_lookup_table

  SUBROUTINE calculate_lookup_table_indices(&
                           ! IN
                           kbdim,klev,kproma, &
                           rhop,qifrac,qitot,nitot,      &
                           ! OUT
                           dum1, dum2, dum4, dum5, dumi, dumk, dumii, dumjj)

    IMPLICIT NONE
    integer, INTENT(in)  :: kbdim, klev, kproma
    REAL(dp), INTENT(in), DIMENSION(kbdim,klev) :: rhop, qitot, qifrac, nitot
    REAL(dp), INTENT(out), DIMENSION(kbdim,klev) :: dum1,dum2,dum4,dum5
    integer, INTENT(out),  DIMENSION(kbdim,klev) :: dumjj,dumii,dumi,dumk
    REAL(dp)             :: ztmp1(kbdim,klev), ztmp2(kbdim,klev)
    logical, DIMENSION(kbdim,klev) :: ll1, ll2
    ! find indices in 4D ice lookup table
    !------------------------------------------------------------------------------------------!

    ! find index for qi (total ice mass mixing ratio)
    dum1(1:kproma,:) = (LOG10(MAX(qitot(1:kproma,:),1.e-12_dp))+16._dp)*1.41328_dp
    dumi(1:kproma,:) = INT(dum1(1:kproma,:))

    ! set limits to make sure the calculated index doesn't exceed range of lookup table

    ll1(1:kproma,:) = (dum1(1:kproma,:) .GE. DBLE(isize))
    ll2(1:kproma,:) = (dum1(1:kproma,:) .LE. 1._dp)
    dum1(1:kproma,:) = MERGE(DBLE(isize), dum1(1:kproma,:), ll1(1:kproma,:))
    dum1(1:kproma,:) = MERGE(1._dp, dum1(1:kproma,:), ll2(1:kproma,:))
    ll1(1:kproma,:) = (dumi(1:kproma,:) .GE. isize-1)
    ll2(1:kproma,:) = (dumi(1:kproma,:) .LE. 1)
    dumi(1:kproma,:) = MERGE(isize-1, dumi(1:kproma,:), ll1(1:kproma,:))
    dumi(1:kproma,:) = MERGE(1, dumi(1:kproma,:), ll2(1:kproma,:))

    ! find index for Ni (ice number mixing ratio)
    dum2(1:kproma,:) = (log10(MAX(nitot(1:kproma,:),1._dp))+10._dp)*1.10731
    dumk(1:kproma,:) = INT(dum2(1:kproma,:))

    ! set limits to make sure the calculated index doesn't exceed range of lookup table

    ll1(1:kproma,:) = (dum2(1:kproma,:) .GE. DBLE(jsize))
    ll2(1:kproma,:) = (dum2(1:kproma,:) .LE. 1._dp)
    dum2(1:kproma,:) = MERGE(DBLE(jsize), dum2(1:kproma,:), ll1(1:kproma,:))
    dum2(1:kproma,:) = MERGE(1._dp, dum2(1:kproma,:), ll2(1:kproma,:))
    ll1(1:kproma,:) = (dumk(1:kproma,:) .GE. jsize-1)
    ll2(1:kproma,:) = (dumk(1:kproma,:) .LE. 1)
    dumk(1:kproma,:) = MERGE(jsize-1, dumk(1:kproma,:), ll1(1:kproma,:))
    dumk(1:kproma,:) = MERGE(1, dumk(1:kproma,:), ll2(1:kproma,:))

    ! find index for rime mass fraction
    IF(l2moment) THEN
       dum4(1:kproma,:) = 1._dp
       dumii(1:kproma,:) = 1
    ELSE
       dum4(1:kproma,:)  = qifrac(1:kproma,:)*3._dp + 1._dp
       dumii(1:kproma,:) = INT(dum4(1:kproma,:))
    ENDIF ! l2moment

    ! set limits
    ll1(1:kproma,:) = (dum4(1:kproma,:) .GE. DBLE(rimsize))
    ll2(1:kproma,:) = (dum4(1:kproma,:) .LE. 1._dp)
    dum4(1:kproma,:) = MERGE(DBLE(rimsize), dum4(1:kproma,:), ll1(1:kproma,:))
    dum4(1:kproma,:) = MERGE(1._dp, dum4(1:kproma,:), ll2(1:kproma,:))
    ll1(1:kproma,:) = (dumii(1:kproma,:) .GE. rimsize-1)
    ll2(1:kproma,:) = (dumii(1:kproma,:) .LE. 1)
    dumii(1:kproma,:) = MERGE(rimsize-1, dumii(1:kproma,:), ll1(1:kproma,:))
    dumii(1:kproma,:) = MERGE(1, dumii(1:kproma,:), ll2(1:kproma,:))

    ! find index for bulk rime density
    ! account for uneven spacing in lookup table for density
    IF(l2moment) THEN
       dum5(1:kproma,:) = 1._dp
       dumjj(1:kproma,:) = 1
    ELSE
       ll1(1:kproma,:) = (rhop(1:kproma,:) .LE. 650._dp)
       ztmp1(1:kproma,:) = (rhop(1:kproma,:)-50._dp)*0.005_dp + 1._dp
       ztmp2(1:kproma,:) = (rhop(1:kproma,:)-650._dp)*0.004_dp+ 4._dp
       dum5(1:kproma,:) = MERGE(ztmp1(1:kproma,:), ztmp2(1:kproma,:), ll1(1:kproma,:))
       dumjj(1:kproma,:) = INT(dum5(1:kproma,:))
    ENDIF ! l2moment

  ! set limits

    ll1(1:kproma,:) = (dum5(1:kproma,:) .GE. DBLE(densize))
    ll2(1:kproma,:) = (dum5(1:kproma,:) .LE. 1._dp)
    dum5(1:kproma,:) = MERGE(DBLE(densize), dum5(1:kproma,:), ll1(1:kproma,:))
    dum5(1:kproma,:) = MERGE(1._dp, dum5(1:kproma,:), ll2(1:kproma,:))
    ll1(1:kproma,:) = (dumjj(1:kproma,:) .GE. densize-1)
    ll2(1:kproma,:) = (dumjj(1:kproma,:) .LE. 1)
    dumjj(1:kproma,:) = MERGE(densize-1, dumjj(1:kproma,:), ll1(1:kproma,:))
    dumjj(1:kproma,:) = MERGE(1, dumjj(1:kproma,:), ll2(1:kproma,:))

  END SUBROUTINE calculate_lookup_table_indices

  SUBROUTINE calculate_lookup_table_indices_1d(&
                           ! IN
                           kbdim,kproma, &
                           rhop,qifrac,qitot,nitot, &
                           ! OUT
                           dum1, dum2, dum4, dum5, dumi, dumk, dumii, dumjj)

    IMPLICIT NONE
    integer, INTENT(in)  :: kbdim, kproma
    REAL(dp), INTENT(in) :: rhop(kbdim), qifrac(kbdim), qitot(kbdim), nitot(kbdim)
    REAL(dp), INTENT(out), DIMENSION(kbdim) :: dum1,dum2,dum4,dum5
    integer, INTENT(out),  DIMENSION(kbdim) :: dumjj,dumii,dumi,dumk
    REAL(dp)             :: ztmp1(kbdim), ztmp2(kbdim)
    logical, DIMENSION(kbdim) :: ll1, ll2
  ! find indices in 4D ice lookup table
  !------------------------------------------------------------------------------------------!

  ! find index for qi (total ice mass mixing ratio)
    dum1(1:kproma) = (LOG10(MAX(qitot(1:kproma),1.e-12_dp))+16._dp)*1.41328_dp
    dumi(1:kproma) = INT(dum1(1:kproma))

  ! set limits to make sure the calculated index doesn't exceed range of lookup table

    ll1(1:kproma) = (dum1(1:kproma) .GE. DBLE(isize))
    ll2(1:kproma) = (dum1(1:kproma) .LE. 1._dp)
    dum1(1:kproma) = MERGE(DBLE(isize), dum1(1:kproma), ll1(1:kproma))
    dum1(1:kproma) = MERGE(1._dp, dum1(1:kproma), ll2(1:kproma))
    ll1(1:kproma) = (dumi(1:kproma) .GE. isize-1)
    ll2(1:kproma) = (dumi(1:kproma) .LE. 1)
    dumi(1:kproma) = MERGE(isize-1, dumi(1:kproma), ll1(1:kproma))
    dumi(1:kproma) = MERGE(1, dumi(1:kproma), ll2(1:kproma))

  ! find index for Ni (ice number mixing ratio)
    dum2(1:kproma) = (log10(MAX(nitot(1:kproma),1._dp))+10._dp)*1.10731
    dumk(1:kproma) = INT(dum2(1:kproma))

  ! set limits to make sure the calculated index doesn't exceed range of lookup table

    ll1(1:kproma) = (dum2(1:kproma) .GE. DBLE(jsize))
    ll2(1:kproma) = (dum2(1:kproma) .LE. 1._dp)
    dum2(1:kproma) = MERGE(DBLE(jsize), dum2(1:kproma), ll1(1:kproma))
    dum2(1:kproma) = MERGE(1._dp, dum2(1:kproma), ll2(1:kproma))
    ll1(1:kproma) = (dumk(1:kproma) .GE. jsize-1)
    ll2(1:kproma) = (dumk(1:kproma) .LE. 1)
    dumk(1:kproma) = MERGE(jsize-1, dumk(1:kproma), ll1(1:kproma))
    dumk(1:kproma) = MERGE(1, dumk(1:kproma), ll2(1:kproma))

    ! find index for rime mass fraction
    IF(l2moment) THEN
       dum4(1:kproma) = 1._dp
       dumii(1:kproma) = 1
    ELSE
       dum4(1:kproma)  = qifrac(1:kproma)*3._dp + 1._dp
       dumii(1:kproma) = INT(dum4(1:kproma))
    ENDIF ! l2moment

  ! set limits

    ll1(1:kproma) = (dum4(1:kproma) .GE. DBLE(rimsize))
    ll2(1:kproma) = (dum4(1:kproma) .LE. 1._dp)
    dum4(1:kproma) = MERGE(DBLE(rimsize), dum4(1:kproma), ll1(1:kproma))
    dum4(1:kproma) = MERGE(1._dp, dum4(1:kproma), ll2(1:kproma))
    ll1(1:kproma) = (dumii(1:kproma) .GE. rimsize-1)
    ll2(1:kproma) = (dumii(1:kproma) .LE. 1)
    dumii(1:kproma) = MERGE(rimsize-1, dumii(1:kproma), ll1(1:kproma))
    dumii(1:kproma) = MERGE(1, dumii(1:kproma), ll2(1:kproma))

    ! find index for bulk rime density
    ! account for uneven spacing in lookup table for density
    IF(l2moment) THEN
       dum5(1:kproma) = 1._dp
       dumjj(1:kproma) = 1
    ELSE
       ll1(1:kproma) = (rhop(1:kproma) .LE. 650._dp)
       ztmp1(1:kproma) = (rhop(1:kproma)-50._dp)*0.005_dp + 1._dp
       ztmp2(1:kproma) = (rhop(1:kproma)-650._dp)*0.004_dp+ 4._dp
       dum5(1:kproma) = MERGE(ztmp1(1:kproma), ztmp2(1:kproma), ll1(1:kproma))
       dumjj(1:kproma) = INT(dum5(1:kproma))
    ENDIF ! l2moment

  ! set limits

    ll1(1:kproma) = (dum5(1:kproma) .GE. DBLE(densize))
    ll2(1:kproma) = (dum5(1:kproma) .LE. 1._dp)
    dum5(1:kproma) = MERGE(DBLE(densize), dum5(1:kproma), ll1(1:kproma))
    dum5(1:kproma) = MERGE(1._dp, dum5(1:kproma), ll2(1:kproma))
    ll1(1:kproma) = (dumjj(1:kproma) .GE. densize-1)
    ll2(1:kproma) = (dumjj(1:kproma) .LE. 1)
    dumjj(1:kproma) = MERGE(densize-1, dumjj(1:kproma), ll1(1:kproma))
    dumjj(1:kproma) = MERGE(1, dumjj(1:kproma), ll2(1:kproma))

  END SUBROUTINE calculate_lookup_table_indices_1d

  SUBROUTINE calculate_my_lookup_table_indices_2d(&
                !--IN
                kproma, kbdim, klev, &
                qnorm, fr, rhor, &
                !--OUT
                iqnorm, irhor, ifr)

    INTEGER, INTENT(IN)   :: kproma, kbdim, klev
    REAL(dp), INTENT(IN)  :: qnorm(kbdim,klev)  !< qi/ni       mass per particle
    REAL(dp), INTENT(IN)  :: fr(kbdim,klev)     !< qirim/qi    rime fraction
    REAL(dp), INTENT(IN)  :: rhor(kbdim,klev)   !< qirim/birim rime density

    REAL(dp), INTENT(OUT) :: iqnorm(kbdim,klev) !< exact position in iqnorm, irhor, ifr
    REAL(dp), INTENT(OUT) :: irhor(kbdim,klev)  !< space (used later to interpolate to
    REAL(dp), INTENT(OUT) :: ifr(kbdim,klev)    !< discrete indices in table

    ! log spaced indices in qnorm table
    iqnorm(1:kproma,:) = (LOG10(MAX(qnmin,qnorm(1:kproma,:)))-LOG10(qnmin))*qnsize &
                         /(LOG10(qnmax)-LOG10(qnmin))
    ! linearly spaced indices in rhor and fr table
    irhor(1:kproma,:)  = (rhor(1:kproma,:)-rhomin)*rhosize/(rhomax-rhomin)
    ifr(1:kproma,:)    = (fr(1:kproma,:)-frmin)*frsize/(frmax-frmin)

  END SUBROUTINE calculate_my_lookup_table_indices_2d

  SUBROUTINE calculate_my_lookup_table_indices_1d(&
                !--IN
                kproma, kbdim, &
                qnorm, fr, rhor, &
                !--OUT
                iqnorm, irhor, ifr)

    INTEGER, INTENT(IN)   :: kproma, kbdim
    REAL(dp), INTENT(IN)  :: qnorm(kbdim)  !< qi/ni       mass per particle
    REAL(dp), INTENT(IN)  :: fr(kbdim)     !< qirim/qi    rime fraction
    REAL(dp), INTENT(IN)  :: rhor(kbdim)   !< qirim/birim rime density

    REAL(dp), INTENT(OUT) :: iqnorm(kbdim) !< exact position in iqnorm, irhor, ifr
    REAL(dp), INTENT(OUT) :: irhor(kbdim)  !< space (used later to interpolate to
    REAL(dp), INTENT(OUT) :: ifr(kbdim)    !< discrete indices in table

    ! log spaced indices in qnorm table
    iqnorm(1:kproma) = (LOG10(MAX(qnmin,qnorm(1:kproma)))-LOG10(qnmin))*qnsize &
                         /(LOG10(qnmax)-LOG10(qnmin))
    ! linearly spaced indices in rhor and fr table
    irhor(1:kproma)  = (rhor(1:kproma)-rhomin)*rhosize/(rhomax-rhomin)
    ifr(1:kproma)    = (fr(1:kproma)-frmin)*frsize/(frmax-frmin)

  END SUBROUTINE calculate_my_lookup_table_indices_1d

  SUBROUTINE access_my_lookup_table_2d(&
                !--IN
                kproma, kbdim, klev, &
                iqnorm, ifr, irhor, idx, &
                !--OUT
                proc)

    INTEGER, INTENT(IN)      :: kproma, kbdim, klev
    REAL(dp), INTENT(INOUT)  :: iqnorm(kbdim,klev) !< exact position in iqnorm, irhor, ifr
    REAL(dp), INTENT(INOUT)  :: irhor(kbdim,klev)  !< space (used later to interpolate to
    REAL(dp), INTENT(INOUT)  :: ifr(kbdim,klev)    !< discrete indices in table
    INTEGER, INTENT(IN)      :: idx                !< index of the value from table

    REAL(dp), INTENT(OUT) :: proc(kbdim,klev)   !< interpolated value from table

    INTEGER, DIMENSION(kbdim,klev) :: iqnorm0, iqnorm1, ifr0, ifr1, irhor0, irhor1

    REAL(dp), DIMENSION(kbdim,klev) :: xd, yd, zd
    REAL(dp), DIMENSION(kbdim,klev) :: v000, v100, v010, v001, v110, v011, v101, v111
    REAL(dp), DIMENSION(kbdim,klev) :: c00, c10, c01, c11
    REAL(dp), DIMENSION(kbdim,klev) :: c0, c1
    REAL(dp), DIMENSION(kbdim,klev) :: c
    INTEGER :: jl, jk

    ! apply table bounds (TODO include check BEFORE call to this routine and adjust nitot/qirim/birim)
    iqnorm(1:kproma,:) = MIN(MAX(1._dp, iqnorm(1:kproma,:)), dble(qnsize))
    irhor(1:kproma,:)  = MIN(MAX(1._dp, irhor(1:kproma,:)), dble(rhosize))
    ifr(1:kproma,:)    = MIN(MAX(1._dp, ifr(1:kproma,:)), dble(frsize))

    ! calculate integers above and below the decimals
    iqnorm1(1:kproma,:) = MIN(MAX(1, CEILING(iqnorm(1:kproma,:))), qnsize) 
    iqnorm0(1:kproma,:) = MIN(MAX(1, CEILING(iqnorm(1:kproma,:))-1), qnsize)
    irhor1(1:kproma,:) = MIN(MAX(1, CEILING(irhor(1:kproma,:))), rhosize) 
    irhor0(1:kproma,:) = MIN(MAX(1, CEILING(irhor(1:kproma,:))-1), rhosize)
    ifr1(1:kproma,:) = MIN(MAX(1, CEILING(ifr(1:kproma,:))), frsize) 
    ifr0(1:kproma,:) = MIN(MAX(1, CEILING(ifr(1:kproma,:))-1), frsize)

    ! get distance to next table index
    xd(1:kproma,:) = iqnorm(1:kproma,:) - iqnorm0(1:kproma,:)
    yd(1:kproma,:) = irhor(1:kproma,:) - irhor0(1:kproma,:)
    zd(1:kproma,:) = ifr(1:kproma,:) - ifr0(1:kproma,:)

    ! get all corners of the cube in 3 dimensions (qnorm, rhor, fr)
    DO jk=1,klev
       DO jl=1,kproma
          v000(jl,jk) = mytab(ifr0(jl,jk), irhor0(jl,jk), iqnorm0(jl,jk), idx)
          v100(jl,jk) = mytab(ifr0(jl,jk), irhor0(jl,jk), iqnorm1(jl,jk), idx)
          v010(jl,jk) = mytab(ifr1(jl,jk), irhor0(jl,jk), iqnorm0(jl,jk), idx)
          v001(jl,jk) = mytab(ifr0(jl,jk), irhor1(jl,jk), iqnorm0(jl,jk), idx)
          v110(jl,jk) = mytab(ifr1(jl,jk), irhor0(jl,jk), iqnorm1(jl,jk), idx)
          v011(jl,jk) = mytab(ifr1(jl,jk), irhor1(jl,jk), iqnorm0(jl,jk), idx)
          v101(jl,jk) = mytab(ifr0(jl,jk), irhor1(jl,jk), iqnorm1(jl,jk), idx)
          v111(jl,jk) = mytab(ifr1(jl,jk), irhor1(jl,jk), iqnorm1(jl,jk), idx)
       END DO ! jl
    END DO ! jk

    ! do a trilinear interpolation (Note that this could be refined as qnorm is logarithmically
    ! spaced in the table.

    ! interpolate in x dimension
    c00(1:kproma,:)  = v000(1:kproma,:)*(1-xd(1:kproma,:)) + v100(1:kproma,:)*xd(1:kproma,:)
    c10(1:kproma,:)  = v010(1:kproma,:)*(1-xd(1:kproma,:)) + v110(1:kproma,:)*xd(1:kproma,:)
    c01(1:kproma,:)  = v001(1:kproma,:)*(1-xd(1:kproma,:)) + v101(1:kproma,:)*xd(1:kproma,:)
    c11(1:kproma,:)  = v011(1:kproma,:)*(1-xd(1:kproma,:)) + v111(1:kproma,:)*xd(1:kproma,:)

    ! interpolate in y dimension
    c0(1:kproma,:)   = c00(1:kproma,:)*(1-yd(1:kproma,:)) + c10(1:kproma,:)*yd(1:kproma,:)
    c1(1:kproma,:)   = c01(1:kproma,:)*(1-yd(1:kproma,:)) + c11(1:kproma,:)*yd(1:kproma,:)

    ! interpolate in z dimension
    c(1:kproma,:)    = c0(1:kproma,:)*(1-zd(1:kproma,:)) + c1(1:kproma,:)*zd(1:kproma,:)

    proc(1:kproma,:) = c(1:kproma,:)

  END SUBROUTINE access_my_lookup_table_2d

  SUBROUTINE access_my_lookup_table_1d(&
                !--IN
                kproma, kbdim, &
                iqnorm, ifr, irhor, idx, &
                !--OUT
                proc)

    INTEGER, INTENT(IN)      :: kproma, kbdim
    REAL(dp), INTENT(INOUT)  :: iqnorm(kbdim) !< exact position in iqnorm, irhor, ifr
    REAL(dp), INTENT(INOUT)  :: irhor(kbdim)  !< space (used later to interpolate to
    REAL(dp), INTENT(INOUT)  :: ifr(kbdim)    !< discrete indices in table
    INTEGER, INTENT(IN)      :: idx                !< index of the value from table

    REAL(dp), INTENT(OUT) :: proc(kbdim)   !< interpolated value from table

    INTEGER, DIMENSION(kbdim) :: iqnorm0, iqnorm1, ifr0, ifr1, irhor0, irhor1

    REAL(dp), DIMENSION(kbdim) :: xd, yd, zd
    REAL(dp), DIMENSION(kbdim) :: v000, v100, v010, v001, v110, v011, v101, v111
    REAL(dp), DIMENSION(kbdim) :: c00, c10, c01, c11
    REAL(dp), DIMENSION(kbdim) :: c0, c1
    REAL(dp), DIMENSION(kbdim) :: c
    INTEGER :: jl

    ! apply table bounds (TODO include check BEFORE call to this routine and adjust nitot/qirim/birim)
    iqnorm(1:kproma) = MIN(MAX(1._dp, iqnorm(1:kproma)), dble(qnsize))
    irhor(1:kproma)  = MIN(MAX(1._dp, irhor(1:kproma)), dble(rhosize))
    ifr(1:kproma)    = MIN(MAX(1._dp, ifr(1:kproma)), dble(frsize))

    ! calculate integers above and below the decimals
    iqnorm1(1:kproma) = MIN(MAX(1, CEILING(iqnorm(1:kproma))), qnsize) 
    iqnorm0(1:kproma) = MIN(MAX(1, CEILING(iqnorm(1:kproma))-1), qnsize)
    irhor1(1:kproma) = MIN(MAX(1, CEILING(irhor(1:kproma))), rhosize) 
    irhor0(1:kproma) = MIN(MAX(1, CEILING(irhor(1:kproma))-1), rhosize)
    ifr1(1:kproma) = MIN(MAX(1, CEILING(ifr(1:kproma))), frsize) 
    ifr0(1:kproma) = MIN(MAX(1, CEILING(ifr(1:kproma))-1), frsize)

    ! get distance to next table index
    xd(1:kproma) = iqnorm(1:kproma) - iqnorm0(1:kproma)
    yd(1:kproma) = irhor(1:kproma) - irhor0(1:kproma)
    zd(1:kproma) = ifr(1:kproma) - ifr0(1:kproma)

    ! get all corners of the cube in 3 dimensions (qnorm, rhor, fr)
    DO jl=1,kproma
       v000(jl) = mytab(ifr0(jl), irhor0(jl), iqnorm0(jl), idx)
       v100(jl) = mytab(ifr0(jl), irhor0(jl), iqnorm1(jl), idx)
       v010(jl) = mytab(ifr1(jl), irhor0(jl), iqnorm0(jl), idx)
       v001(jl) = mytab(ifr0(jl), irhor1(jl), iqnorm0(jl), idx)
       v110(jl) = mytab(ifr1(jl), irhor0(jl), iqnorm1(jl), idx)
       v011(jl) = mytab(ifr1(jl), irhor1(jl), iqnorm0(jl), idx)
       v101(jl) = mytab(ifr0(jl), irhor1(jl), iqnorm1(jl), idx)
       v111(jl) = mytab(ifr1(jl), irhor1(jl), iqnorm1(jl), idx)
    END DO ! jl

    ! do a trilinear interpolation (Note that this could be refined as qnorm is logarithmically
    ! spaced in the table.

    ! interpolate in x dimension
    c00(1:kproma)  = v000(1:kproma)*(1-xd(1:kproma)) + v100(1:kproma)*xd(1:kproma)
    c10(1:kproma)  = v010(1:kproma)*(1-xd(1:kproma)) + v110(1:kproma)*xd(1:kproma)
    c01(1:kproma)  = v001(1:kproma)*(1-xd(1:kproma)) + v101(1:kproma)*xd(1:kproma)
    c11(1:kproma)  = v011(1:kproma)*(1-xd(1:kproma)) + v111(1:kproma)*xd(1:kproma)

    ! interpolate in y dimension
    c0(1:kproma)   = c00(1:kproma)*(1-yd(1:kproma)) + c10(1:kproma)*yd(1:kproma)
    c1(1:kproma)   = c01(1:kproma)*(1-yd(1:kproma)) + c11(1:kproma)*yd(1:kproma)

    ! interpolate in z dimension
    c(1:kproma)    = c0(1:kproma)*(1-zd(1:kproma)) + c1(1:kproma)*zd(1:kproma)

    proc(1:kproma) = c(1:kproma)

  END SUBROUTINE access_my_lookup_table_1d
       
  ! lookup table is limited by the fraction qi/ni, therefore
  ! the maximum ni depends on qi. This function returns these ni limits
  SUBROUTINE get_ni_limits_2d(&
                !--IN
                kproma, kbdim, klev, &
                qi, nitot, &
                !--OUT
                nimin, nimax)

    INTEGER, INTENT(IN) :: kproma, kbdim, klev
    REAL(dp), INTENT(IN), DIMENSION(kbdim,klev) :: qi, nitot
  
    REAL(dp), INTENT(OUT), DIMENSION(kbdim,klev) :: nimin, nimax

    nimin(1:kproma,:) = qi(1:kproma,:)/qnmax
    nimax(1:kproma,:) = qi(1:kproma,:)/qnmin

  END SUBROUTINE get_ni_limits_2d

  SUBROUTINE get_ni_limits_1d(&
                !--IN
                kproma, kbdim, &
                qi, nitot, &
                !--OUT
                nimin, nimax)

    INTEGER, INTENT(IN) :: kproma, kbdim
    REAL(dp), INTENT(IN), DIMENSION(kbdim) :: qi, nitot
  
    REAL(dp), INTENT(OUT), DIMENSION(kbdim) :: nimin, nimax

    nimin(1:kproma) = qi(1:kproma)/qnmax
    nimax(1:kproma) = qi(1:kproma)/qnmin

  END SUBROUTINE get_ni_limits_1d



END MODULE mo_p3_fields
