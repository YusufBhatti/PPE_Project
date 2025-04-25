! ECHAM post processing
! Generates tables of averages
!
! Replaced gauaw by current echam version to allow to process horizontal resolutions > T213
! Implemented dynamic memory management to allow for single executable
! Lines per page is now calculated automatically
! fglac initialization error marked and fixed, see below

#ifdef HAVE_CONFIG_INC
#include <config.inc>
#else
#define QPLOT_TABLE_LINES 53
#endif

PROGRAM momitt

  USE mo_tables_diff
  USE mo_util_string

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  REAL(dp) , PARAMETER :: m_pi = 3.14159265358979323846D0
  REAL(dp), ALLOCATABLE :: phi(:), pw(:)

  INTEGER, PARAMETER :: npp = 15

  REAL, ALLOCATABLE :: fice(:,:),     fwater(:,:), fland(:,:)
  REAL, ALLOCATABLE :: fglac(:,:),    fprec(:,:),  xfac(:,:)
  REAL, ALLOCATABLE :: fcode91(:,:),  fcode92(:,:)
  REAL, ALLOCATABLE :: fcode93(:,:),  fcode94(:,:)
  REAL, ALLOCATABLE :: fcode95(:,:),  fcode96(:,:)
  REAL, ALLOCATABLE :: fcode110(:,:), fcode111(:,:)
  REAL, ALLOCATABLE :: fcode112(:,:), fcode113(:,:)
  REAL, ALLOCATABLE :: fcode114(:,:), fcode115(:,:)
  REAL, ALLOCATABLE :: fcode119(:,:), fcode120(:,:)
  REAL, ALLOCATABLE :: fcode121(:,:), fcode210(:,:)
  REAL, ALLOCATABLE :: fcode140(:,:), fcode144(:,:)
  REAL, ALLOCATABLE :: fcode160(:,:)
  REAL, ALLOCATABLE :: fcode146(:,:), fcode147(:,:)
  REAL, ALLOCATABLE :: fcode176(:,:), fcode177(:,:)
  REAL, ALLOCATABLE :: fcode178(:,:), fcode179(:,:)
  REAL, ALLOCATABLE :: fcode182(:,:), fcode187(:,:)
  REAL, ALLOCATABLE :: fcode188(:,:)
  REAL, ALLOCATABLE :: fcode218(:,:), fcode228(:,:)
  REAL, ALLOCATABLE :: fcode206(:,:), fcode229(:,:)
  REAL, ALLOCATABLE :: fcode221(:,:), fcode222(:,:)

  REAL, ALLOCATABLE :: field(:,:), dat(:,:,:), zm(:,:), xzm(:)
  REAL, ALLOCATABLE ::  zm1(:,:)

  REAL :: nh(kmax), sh(kmax), gm(kmax), faktor(kmax), addi(kmax)
  REAL :: nh1(kmax), sh1(kmax), gm1(kmax)

  REAL :: xnh, xsh, xgm, sumpnh, sumpsh, sumpgm 

  INTEGER :: year1, year2, ilauf, iunit, landcode
  INTEGER :: ihead(8), icode(kmax), imean

  INTEGER :: nlon, nlat, nlath

  INTEGER :: i, j, k, kk, n, ic, ida, kend, ja, jh, jo, ilayer207

  LOGICAL :: lland, lice, laccu, lglac, lprec
  LOGICAL :: lcode91,lcode92,lcode93,lcode94,lcode95,lcode96 
  LOGICAL :: lcode110,lcode111,lcode112,lcode113,lcode114,lcode115,lcode119 
  LOGICAL :: lcode120,lcode121,lcode140,lcode144,lcode146,lcode147,lcode160
  LOGICAL :: lcode176,lcode177,lcode178,lcode179,lcode182,lcode187,lcode188
  LOGICAL :: lcode191,lcode192,lcode206,lcode210
  LOGICAL :: lcode218,lcode221,lcode222,lcode228,lcode229


  CHARACTER(len=82) :: titel

  CHARACTER(len=80) :: average

  CHARACTER(len=49) :: trange
  CHARACTER(len=60) :: stringo(kmax)

  CHARACTER(len=40) :: half
  CHARACTER(len=20) :: expnam

  CHARACTER(len= 3) :: cmean 

  CHARACTER(len=14) :: month(17) = (/ &
       '  January     ','  February    ','  March       ', &
       '  April       ','  May         ','  June        ', &
       '  July        ','  August      ','  September   ', &
       '  October     ','  November    ','  December    ', &
       's DEC/JAN/FEB ','s MAR/APR/MAY ','s JUN/JUL/AUG ', &
       's SEP/OCT/NOV ','s JAN...DEC   ' /)

  NAMELIST / exper / expnam, average, year1, year2, landcode

  !     -----------------------------------------------------
  !     Read the experimentname from namelist experim

  READ (*,exper)
  write (6,*) 'namelist2 ',expnam, average, year1, year2, landcode
  
  cmean(1:3) = tolower(average)

average_periode: SELECT CASE (cmean)
  CASE ('jan') ! Januray
    imean =  1
  CASE ('feb') ! February
    imean =  2
  CASE ('mar') ! March
    imean =  3
  CASE ('apr') ! April
    imean =  4
  CASE ('may') ! May
    imean =  5
  CASE ('jun') ! June
    imean =  6
  CASE ('jul') ! July
    imean =  7
  CASE ('aug') ! August
    imean =  8
  CASE ('sep') ! September
    imean =  9
  CASE ('oct') ! October
    imean = 10
  CASE ('nov') ! November
    imean = 11
  CASE ('dec') ! December
    imean = 12
  CASE ('djf') ! dec/jan/feb
    imean = 13
  CASE ('mam') ! mar/apr/may
    imean = 14
  CASE ('jja') ! jun/jul/au
    imean = 15
  CASE ('son') ! sep/oct/nov
    imean = 16
  CASE ('ann') ! a whole year 
    imean = 17
  CASE DEFAULT
    WRITE (0,*) 'Unsupported average peride selected ...'
    STOP 'abort ...'
  END SELECT average_periode

  CALL   check_files (kmax,nlon,nlat)
  write(0,*) 'kmax,lon, lat: ',kmax,nlon,nlat

  ! Initilaize tables

  CALL init_tables
  

  OPEN(8,file='tabelle',form='FORMATTED')
  OPEN(9,file='tablecode',form='FORMATTED')
 
  OPEN (unit=10,file='BBOTa.srv',form='UNFORMATTED')
 
  OPEN (unit=20,file='BBOTb.srv',form='UNFORMATTED')
 


  ALLOCATE (phi(nlat), pw(nlat))

  CALL gauaw(phi,pw,nlat)
  DO n = 1, nlat
     phi(n) = ASIN(phi(n))/(2.0D0*m_pi)*360.0D0
  ENDDO

  ALLOCATE (fice    (nlon,nlat), fwater  (nlon,nlat), fland(nlon,nlat))
  ALLOCATE (fglac   (nlon,nlat), fprec   (nlon,nlat), xfac (nlon,nlat))
  ALLOCATE (fcode91 (nlon,nlat), fcode92 (nlon,nlat))
  ALLOCATE (fcode93 (nlon,nlat), fcode94 (nlon,nlat))
  ALLOCATE (fcode95 (nlon,nlat), fcode96 (nlon,nlat))
  ALLOCATE (fcode110(nlon,nlat), fcode111(nlon,nlat))
  ALLOCATE (fcode112(nlon,nlat), fcode113(nlon,nlat))
  ALLOCATE (fcode114(nlon,nlat), fcode115(nlon,nlat))
  ALLOCATE (fcode119(nlon,nlat), fcode120(nlon,nlat))
  ALLOCATE (fcode121(nlon,nlat), fcode210(nlon,nlat))
  ALLOCATE (fcode140(nlon,nlat), fcode144(nlon,nlat))
  ALLOCATE (fcode160(nlon,nlat))
  ALLOCATE (fcode146(nlon,nlat), fcode147(nlon,nlat))
  ALLOCATE (fcode176(nlon,nlat), fcode177(nlon,nlat))
  ALLOCATE (fcode178(nlon,nlat), fcode179(nlon,nlat))
  ALLOCATE (fcode182(nlon,nlat), fcode187(nlon,nlat))
  ALLOCATE (fcode188(nlon,nlat))
  ALLOCATE (fcode218(nlon,nlat), fcode228(nlon,nlat))
  ALLOCATE (fcode206(nlon,nlat), fcode229(nlon,nlat))
  ALLOCATE (fcode221(nlon,nlat), fcode222(nlon,nlat))
      
  ALLOCATE (field(nlon,nlat), dat(nlon,nlat,kmax), zm(nlat,kmax), xzm(nlat))
  ALLOCATE (zm1(nlat,kmax))

  ! fglac is  processed for couple
    fglac(:,:) = 0.0   


  iunit=10

  DO  ilauf = 1,2
  !      Get ready for zonal mean calculation
  k=0
  ilayer207=0
  lland=.FALSE.
  lice=.FALSE.
  lcode91=.FALSE.
  lcode92=.FALSE.
  lcode93=.FALSE.
  lcode94=.FALSE.
  lcode95=.FALSE.
  lcode96=.FALSE.
  lcode110=.FALSE.
  lcode111=.FALSE.
  lcode112=.FALSE.
  lcode113=.FALSE.
  lcode114=.FALSE.
  lcode115=.FALSE.
  lcode119=.FALSE.
  lcode120=.FALSE.
  lcode121=.FALSE.
  lcode140=.FALSE.
  lcode144=.FALSE.
  lcode146=.FALSE.
  lcode147=.FALSE.
  lcode160=.FALSE.
  lcode176=.FALSE.
  lcode177=.FALSE.
  lcode178=.FALSE.
  lcode179=.FALSE.
  lcode182=.FALSE.
  lcode187=.FALSE.
  lcode188=.FALSE.
  lcode191=.FALSE.
  lcode192=.FALSE.
  lcode206=.FALSE.
  lcode210=.FALSE.
  lcode218=.FALSE.
  lcode221=.FALSE.
  lcode222=.FALSE.
  lcode229=.FALSE.
  lglac=.FALSE.
  lprec=.FALSE.

  DO ic = 1, kmax

    READ(iunit,END=100) ihead

    ida        = ihead(1)

!    WRITE (0,'(a,i4,a,i4,a,i4)') 'header - code = ',ihead(1), &
!      ' longitudes = ', ihead(5), &
!      ' latitudes = ', ihead(6)


    READ(iunit,END=100) ((field(i,j),i=1,nlon),j=1,nlat)


    k = k+1

    

    icode(k)   = ihead(1)
    faktor(k)  = factor(ida)
    addi(k)    = offset(ida)
    stringo(k) = string(ida)

    dat(:,:,k) = field(:,:)

    IF (ida == 207) THEN
      ilayer207 = ilayer207 + 1
      IF (ilayer207 == 2) stringo(k)(51:51)='2'
      IF (ilayer207 == 3) stringo(k)(51:51)='3'
      IF (ilayer207 == 4) stringo(k)(51:51)='4'
      IF (ilayer207 == 5) stringo(k)(51:51)='5'
    ENDIF

    IF (ida == 97) THEN
      lice=.TRUE.
      WHERE (field(:,:) < 0.001)
        fice(:,:)  = 0.0
        field(:,:) = 0.0
      ELSEWHERE
        fice(:,:) = field(:,:)
      END WHERE
      dat(:,:,k) = field(:,:)
    ENDIF
    IF (ida == landcode) THEN
      lland=.TRUE.
      fland(:,:)  = field(:,:)
      fwater(:,:) = 1.0-(fice(:,:)+fland(:,:))
      dat(:,:,k)  = field(:,:)
    ENDIF

    IF (ida == 91)  fcode91(:,:)  = field(:,:)
    IF (ida == 91)  lcode91 = .TRUE.
    IF (ida == 92)  fcode92(:,:)  = field(:,:)
    IF (ida == 92)  lcode92 = .TRUE.
    IF (ida == 93)  fcode93(:,:)  = field(:,:)
    IF (ida == 93)  lcode93 = .TRUE.
    IF (ida == 94)  fcode94(:,:)  = field(:,:)
    IF (ida == 94)  lcode94 = .TRUE.
    IF (ida == 95)  fcode95(:,:)  = field(:,:)
    IF (ida == 95)  lcode95 = .TRUE.
    IF (ida == 96)  fcode96(:,:)  = field(:,:)
    IF (ida == 96)  lcode96 = .TRUE.
    IF (ida == 110) fcode110(:,:) = field(:,:)
    IF (ida == 110) lcode110 = .TRUE.
    IF (ida == 111) fcode111(:,:) = field(:,:)
    IF (ida == 111) lcode111 = .TRUE.
    IF (ida == 112) fcode112(:,:) = field(:,:)
    IF (ida == 112) lcode112 = .TRUE.
    IF (ida == 113) fcode113(:,:) = field(:,:)
    IF (ida == 113) lcode113 = .TRUE.
    IF (ida == 114) fcode114(:,:) = field(:,:)
    IF (ida == 114) lcode114 = .TRUE.
    IF (ida == 115) fcode115(:,:) = field(:,:)
    IF (ida == 115) lcode115 = .TRUE.
    IF (ida == 119) fcode119(:,:) = field(:,:)
    IF (ida == 119) lcode119 = .TRUE.
    IF (ida == 120) fcode120(:,:) = field(:,:)
    IF (ida == 120) lcode120 = .TRUE.
    IF (ida == 121) fcode121(:,:) = field(:,:)
    IF (ida == 121) lcode121 = .TRUE.
    IF (ida == 140) fcode140(:,:) = field(:,:)
    IF (ida == 140) lcode140 = .TRUE.
    IF (ida == 144) fcode144(:,:) = field(:,:)
    IF (ida == 144) lcode144 = .TRUE.
    IF (ida == 146) fcode146(:,:) = field(:,:)
    IF (ida == 146) lcode146 = .TRUE.
    IF (ida == 147) fcode147(:,:) = field(:,:)
    IF (ida == 147) lcode147 = .TRUE.
    IF (ida == 160) fcode160(:,:) = field(:,:)
    IF (ida == 160) lcode160 = .TRUE.
    IF (ida == 176) fcode176(:,:) = field(:,:)
    IF (ida == 176) lcode176 = .TRUE.
    IF (ida == 177) fcode177(:,:) = field(:,:)
    IF (ida == 177) lcode177 = .TRUE.
    IF (ida == 178) fcode178(:,:) = field(:,:)
    IF (ida == 178) lcode178 = .TRUE.
    IF (ida == 179) fcode179(:,:) = field(:,:)
    IF (ida == 179) lcode179 = .TRUE.
    IF (ida == 182) fcode182(:,:) = field(:,:)
    IF (ida == 182) lcode182 = .TRUE.
    IF (ida == 187) fcode187(:,:) = field(:,:)
    IF (ida == 187) lcode187 = .TRUE.
    IF (ida == 188) fcode188(:,:) = field(:,:)
    IF (ida == 188) lcode188 = .TRUE.
    IF (ida == 191) lcode191 = .TRUE.
    IF (ida == 192) lcode192 = .TRUE.
    IF (ida == 206) fcode206(:,:) = field(:,:)
    IF (ida == 206) lcode206 = .TRUE.


    IF (ida == 210) fcode210(:,:) = field(:,:)
    IF (ida == 210) lcode210 = .TRUE.
    IF (ida == 218) fcode218(:,:) = field(:,:)
    IF (ida == 218) lcode218 = .TRUE.
    IF (ida == 221) fcode221(:,:) = field(:,:)
    IF (ida == 221) lcode221 = .TRUE.
    IF (ida == 222) fcode222(:,:) = field(:,:)
    IF (ida == 222) lcode222 = .TRUE.
    IF (ida == 229) fcode228(:,:) = field(:,:)
    IF (ida == 229) lcode228 = .TRUE.
    IF (ida == 229) fcode229(:,:) = field(:,:)
    IF (ida == 229) lcode229 = .TRUE.

    IF (ida == 232) THEN
      WHERE (fland > 0.0)
        dat(:,:,k) = field(:,:)
      ELSEWHERE
        dat(:,:,k) = 0.0
      END WHERE
      fglac(:,:) = dat(:,:,k)
      lglac = .TRUE.
    END IF

    IF (ida == 260) fprec(:,:) = field(:,:)
    IF (ida == 260) lprec = .TRUE.


  ENDDO

100 CONTINUE

   IF (.NOT. lland) THEN
     WRITE (0,*) 'land-sea-mask missing, you need code ',landcode
     STOP 'abort ...'
   ENDIF
   IF (.NOT. lice) THEN
     WRITE (0,*) 'ice-mask missing, you need code 97'
     STOP 'abort ...'
   ENDIF
   IF (.NOT. lglac) THEN
     WRITE (0,*) 'glacier-mask missing, you need code 232'
     STOP 'abort ...'
   ENDIF
 

    kend = k
    write(6,*) 'kend ',kend

!--------------------------------------------------------

!    added codes : code191 and code192

     IF (.NOT. lcode191 .AND. lcode178 .AND. lcode187) THEN
      k = k+1
      kend = k
      ida        = 191
      icode(k)   = ida
      faktor(k)  = factor(ida)
      addi(k)    = offset(ida)
      stringo(k) = string(ida)
      dat(:,:,k) = fcode178(:,:)-fcode187(:,:)
     ENDIF

     IF (.NOT. lcode192 .AND. lcode179 .AND. lcode188 ) THEN
      k = k+1
      kend = k
      ida        = 192
      icode(k)   = ida
      faktor(k)  = factor(ida)
      addi(k)    = offset(ida)
      stringo(k) = string(ida)
      dat(:,:,k) = fcode179(:,:)-fcode188(:,:)
     ENDIF


  ! extra-codes fuer fwater (code1-code14)
  write(6,*) 'extra-codes fuer fwater (code1-code14)'

  IF (lprec.AND.lcode91.AND.lcode92.AND.lcode93.AND.lcode94 &
      .AND.lcode95.AND.lcode96.AND.lcode110.AND.lcode111.AND.lcode112 &
      .AND.lcode113.AND.lcode114.AND.lcode115.AND.lcode119.AND.lcode120 &
      .AND.lcode121.AND.lcode140.AND.lcode146.AND.lcode147.AND.lcode160 &
      .AND.lcode176.AND.lcode177.AND.lcode178.AND.lcode179.AND.lcode182 &
      .AND.lcode228.AND.lcode144.AND.lcode218 &
      .AND.lcode221.AND.lcode222.AND.lcode229) THEN
  dat(:,:,k+1) = fprec(:,:)+fcode182(:,:)
  dat(:,:,k+2) = fwater(:,:)
  dat(:,:,k+3) = fprec(:,:)*fice(:,:)+fcode113(:,:)
  dat(:,:,k+4) = fprec(:,:)*fwater(:,:)+fcode114(:,:)
  dat(:,:,k+5) = fprec(:,:)*fland(:,:)+fcode115(:,:)       &
       -fcode160(:,:)-fcode221(:,:)-fcode222(:,:)
  dat(:,:,k+6) = fcode178(:,:)+fcode179(:,:)
  dat(:,:,k+7) = fcode146(:,:)+fcode147(:,:)+fcode176(:,:)+fcode177(:,:)
  dat(:,:,k+8) = fcode140(:,:)/(fcode229(:,:)+0.001)*(1.-fglac(:,:))
  dat(:,:,k+9) = fcode121(:,:)+fcode112(:,:)               &
       +fcode96(:,:)+fcode93(:,:)                         &
       -(fcode218(:,:)+fcode228(:,:))*3.337E05
  dat(:,:,k+10)= fcode120(:,:)+fcode111(:,:)               &
       +fcode95(:,:)+fcode92(:,:)-fcode144(:,:)*fwater(:,:)*3.337E05
  dat(:,:,k+11)= fcode119(:,:)+fcode110(:,:)               &
       +fcode94(:,:)+fcode91(:,:)
  dat(:,:,k+12) = fprec(:,:)
  dat(:,:,k+13) = fprec(:,:)
  dat(:,:,k+14) = fprec(:,:)

  icode  (k+1)   =  1        !
  icode  (k+2)   =  2        !
  icode  (k+3)   =  3        !
  icode  (k+4)   =  4        !
  icode  (k+5)   =  5        !
  icode  (k+6)   =  6        !
  icode  (k+7)   =  7        !
  icode  (k+8)   =  8        !
  icode  (k+9)   =  9        !
  icode  (k+10)  = 10       !
  icode  (k+11)  = 11       !
  icode  (k+12)  = 12       ! total precipitation - ICE
  icode  (k+13)  = 13       ! total precipitation - WATER
  icode  (k+14)  = 14       ! total precipitation - LAND

  faktor (k+1)   = factor(266)
  faktor (k+2)   = factor(267)
  faktor (k+3)   = factor(268)
  faktor (k+4)   = factor(269)
  faktor (k+5)   = factor(270)
  faktor (k+6)   = factor(261)
  faktor (k+7)   = factor(263)
  faktor (k+8)   = factor(265)
  faktor (k+9)   = factor(271)
  faktor (k+10)  = factor(272)
  faktor (k+11)  = factor(273)
  faktor (k+12)  = factor(274)
  faktor (k+13)  = factor(275)
  faktor (k+14)  = factor(276)

  addi   (k+1)   = offset(266)
  addi   (k+2)   = offset(267)
  addi   (k+3)   = offset(268)
  addi   (k+4)   = offset(269)
  addi   (k+5)   = offset(270)
  addi   (k+6)   = offset(261)
  addi   (k+7)   = offset(263)
  addi   (k+8)   = offset(265)
  addi   (k+9)   = offset(271)
  addi   (k+10)  = offset(272)
  addi   (k+11)  = offset(273)
  addi   (k+12)  = offset(274)
  addi   (k+13)  = offset(275)
  addi   (k+14)  = offset(276)

  stringo(k+1)   = string(266)
  stringo(k+2)   = string(267)
  stringo(k+3)   = string(268)
  stringo(k+4)   = string(269)
  stringo(k+5)   = string(270)
  stringo(k+6)   = string(261)
  stringo(k+7)   = string(263)
  stringo(k+8)   = string(265)
  stringo(k+9)   = string(271)
  stringo(k+10)  = string(272)
  stringo(k+11)  = string(273)
  stringo(k+12)  = string(274)
  stringo(k+13)  = string(275)
  stringo(k+14)  = string(276)

  kend = k+14

  ELSE
     write(6,*) 'would you compute extra-codes fuer fwater (code1-code14)'
     write(6,*) 'you need the following codes: 260,232,91,92,93,94,95,96,'
     write(6,*) '    110,111,112,113,114,115,119,120,121,140,146,147,160,'
     write(6,*) '    176,177,178,179,182,208,209,218,221,222,228,229'
  ENDIF ! extra-codes fuer fwater (code1-code14)
!---------------------------------------------------------------



  WRITE (0,*) 'number of datasets found: ', kend



  !      Zonal and global mean calculation

  zm(:,:) = 0.0

  DO k=1,kend

    WRITE (0,'(a,i4)') 'process ', icode(k)

    laccu = .FALSE.

    nh(k) = 0.0
    xnh   = 0.0
    sh(k) = 0.0
    xsh   = 0.0
    gm(k) = 0.0
    xgm   = 0.0
    sumpnh = 0.0000000001
    sumpsh = 0.0000000001
    sumpgm = 0.0000000001

    IF (stringo(k)(60:60) == 'A') laccu = .TRUE.

    IF (stringo(k)(59:59) == 'L') THEN
      xfac(:,:) = fland(:,:)
    ELSE IF (stringo(k)(59:59) == 'W') THEN
      xfac(:,:) = fwater(:,:)
    ELSE IF (stringo(k)(59:59) == 'I') THEN
      xfac(:,:) = fice(:,:)
    ELSE IF (icode(k) == 213) THEN
      xfac(:,:) = fland(:,:)*fglac(:,:)
    ELSE IF (icode(k) == 210 .OR. icode(k) == 219) THEN
      xfac(:,:) = 1.0-fland(:,:)
    ELSE IF (icode(k) == 214) THEN
      xfac(:,:) = fcode210(:,:)*(1.0-fland(:,:))
    ELSE
      xfac(:,:) = 1.0
    ENDIF

    xzm(:)  = 0.0
    DO j = 1, nlat
      IF (laccu) THEN
        zm(j,k) = zm(j,k) + SUM(dat(:,j,k))
      ELSE
        zm(j,k) = zm(j,k) + SUM(dat(:,j,k)*xfac(:,j))
      END IF
      xzm(j) = xzm(j) + SUM(xfac(:,j))

      IF (phi(j) > 0.0) THEN
        nh(k) = nh(k)+zm(j,k)*pw(j)
        xnh = xnh+xzm(j)
        sumpnh = sumpnh+pw(j)*xzm(j)
      ELSE
        sh(k) = sh(k)+zm(j,k)*pw(j)
        xsh = xsh+xzm(j)
        sumpsh = sumpsh+pw(j)*xzm(j)
      ENDIF
      IF (xzm(j) > 0.0) THEN
        zm(j,k) = zm(j,k)/xzm(j)
        if(ilauf.eq.2) zm(j,k)=zm1(j,k)-zm(j,k)
      ELSE
        zm(j,k) = -1000000000
      ENDIF
    ENDDO

    !       mean values

    gm(k) = nh(k)+sh(k)
    xgm = xnh+xsh
    sumpgm = sumpnh+sumpsh
    IF (xnh > 0.0) THEN
      nh(k) = nh(k)/sumpnh
      if(ilauf.eq.2) nh(k)=nh1(k)-nh(k)
    ELSE
      nh(k) = -1000000000
    ENDIF
    IF (xsh > 0.0) THEN
      sh(k) = sh(k)/sumpsh
      if(ilauf.eq.2) sh(k)=sh1(k)-sh(k)
    ELSE
      sh(k) = -1000000000
    ENDIF
    IF (xgm > 0.0) THEN
      gm(k) = gm(k)/sumpgm
      if(ilauf.eq.2) gm(k)=gm1(k)-gm(k)
    ELSE
      gm(k) = -1000000000
    ENDIF

  END DO
!
!  save values of first dataset on zm1, nh1,sh1 and gm1
!
      IF(ilauf.eq.1) then
         do k=1,kend
           do j=1,nlat
             zm1(j,k)=zm(j,k)
           enddo
           sh1(k)=sh(k)
           nh1(k)=nh(k)
           gm1(k)=gm(k)
         enddo
         iunit=20
      ENDIF
  END DO !ENDDO ilauf

  !      and now print the tables

  nlath = nlat/2
  trange = ' Mean of Month               Years      to     '

  WRITE(trange(15:28),'(a14)') month(imean)
  WRITE(trange(37:49),'(i4,a4,i4)')  year1,' to ',year2

  half = 'northern hemisphere                     '

!!  jo = nlat/(nlat/30)
!!  IF (jo > nlath) jo = nlath

  IF (nlat.eq.32) THEN
     jo = nlat
  ELSEIF (nlat.eq.48.or.nlat.eq.64.or.nlat.eq.96) THEN
     jo = nlat/2
  ELSEIF (nlat.eq.128) THEN
     jo = 32
  ELSEIF (nlat.eq.160.or.nlat.eq.240.or.nlat.eq.320) THEN
     jo = 40
  ELSEIF (nlat.eq.480) THEN
     jo = 48
  ELSE
     jo=50
  ENDIF

  ja = 1
  jh = jo
  DO WHILE (jh <= nlat)

    titel(:) = ' '
    WRITE(titel( 1:20),'(a20)') expnam
    WRITE(titel(34:  ),'(a49)') trange

    CALL output(titel, zm, faktor, addi, nh, sh, gm, phi, icode, stringo, &
         kend, kmax, ja, jh, nlat, npp, half)

    ja = ja+jo
    jh = jh+jo

    IF (jh <= nlath) THEN
      half = 'northern hemisphere                     '
    ELSE
      half = 'southern hemisphere                     '
    ENDIF
  ENDDO

  WRITE(9,'(a1,a82)') ' ', titel
  WRITE(9,*)
  WRITE(9,*)
  WRITE(9,*)' code   name          unit   description'
  WRITE(9,*)
  DO kk=1,kend
    WRITE(9,'(1x,1i5,3x,1a60)') icode(kk), stringo(kk)
  ENDDO
END PROGRAM momitt

SUBROUTINE output(title,zonal,faktor,addi,north,south,global,breite, &
     icode, cstring, k, kmax, ja, je, nlat, npp, half)

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  !    Ausgabe von Zonal- und Globalmittelwerten in Tabellenform

  INTEGER, INTENT(in) :: nlat, kmax, k, ja, je, npp

  REAL, INTENT(in) :: zonal(nlat,kmax), global(kmax), faktor(kmax), addi(kmax)
  REAL(dp), INTENT(in) :: breite(nlat)
  REAL, INTENT(in) :: north(kmax), south(kmax)

  INTEGER, INTENT(in) :: icode(kmax)

  CHARACTER(len=*), INTENT(in) :: cstring(*), title, half

  INTEGER :: ka, ke, nlath, kk, j, kzah

  ka = 1
  ke = npp

  DO

    WRITE(8,'(a1,a82)') ' ', title
!    WRITE(8,*) half
    WRITE(8,*)
    WRITE(8,*)

    ke = MIN(ke,k)
    nlath=nlat/2

    ! Print zonal mean fields

    ! formats are fixed for 15 codes per single page !!!!

    WRITE(8,'(1x," code:  |",15(i5,"   |"))') (icode(kk),kk=ka,ke)
    WRITE(8,'(1x,"--------",15a9,"|")') (('+--------'),kk=ka,ke)
    kzah=0
    DO j=ja,je
      WRITE(8,'(1x,16(1f8.3,"|"))') breite(j),(zonal(j,kk)*faktor(kk)+addi(kk),kk=ka,ke)
      kzah=kzah+1
    ENDDO

!new    WRITE(8,'(1x,"--------",15a9,"|")') (('+--------'),kk=ka,ke)
!new     IF (je <= nlath) THEN
!new       WRITE(8,'(1x," north: |",15(1f8.3,"|"))') (north(kk)*faktor(kk)+addi(kk),kk=ka,ke)
!new     ELSE
!new       WRITE(8,'(1x," south: |",15(1f8.3,"|"))') (south(kk)*faktor(kk)+addi(kk),kk=ka,ke)
!new     ENDIF
!new     WRITE(8,'(1x,"global: |",15(1f8.3,"|"))') (global(kk)*faktor(kk)+addi(kk),kk=ka,ke)
!new 

    WRITE(8,'(1x,"--------",15a9,"|")') (('+--------'),kk=ka,ke)
    IF (je <= nlath) THEN
      WRITE(8,'(1x," north: |",15(1f8.3,"|"))') (north(kk)*faktor(kk)+addi(kk),kk=ka,ke)
      kzah=kzah+1
    ELSEIF (nlat.eq.(je-1+ja)) THEN
      WRITE(8,'(1x," south: |",15(1f8.3,"|"))') (south(kk)*faktor(kk)+addi(kk),kk=ka,ke)
      WRITE(8,'(1x," north: |",15(1f8.3,"|"))') (north(kk)*faktor(kk)+addi(kk),kk=ka,ke)
      kzah=kzah+2
    ELSE
      WRITE(8,'(1x," south: |",15(1f8.3,"|"))') (south(kk)*faktor(kk)+addi(kk),kk=ka,ke)
      kzah=kzah+1
    ENDIF
    WRITE(8,'(1x,"global: |",15(1f8.3,"|"))') (global(kk)*faktor(kk)+addi(kk),kk=ka,ke)


      DO kk=kzah+1,QPLOT_TABLE_LINES
        WRITE(8,'(1x)')
      ENDDO


    IF ( ke == k ) EXIT
    ka=ke+1
    ke=ke+npp

  ENDDO

END SUBROUTINE output

SUBROUTINE gauaw (pa, pw, nlat)

  ! Compute abscissas and weights for gaussian integration.

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  REAL(dp), PARAMETER :: api = 3.14159265358979323846_dp

  !  Scalar arguments
  INTEGER :: nlat

  !  Array arguments
  REAL(dp) :: pa(nlat), pw(nlat)
  ! pa  - array, length at least k, to receive abscissas.
  ! pw  - array, length at least k, to receive weights.

  !  Local scalars:
  REAL(dp), PARAMETER :: eps = EPSILON(0.0_dp)
  INTEGER, PARAMETER :: itemax = 20

  INTEGER :: iter, ins2, isym, jn, jgl
  REAL(dp) :: za, zw, z, zan
  REAL(dp) :: zk, zkm1, zkm2, zx, zxn, zldn, zmod

  !  Intrinsic functions
  INTRINSIC ABS, COS, MOD, TAN

  !  Executable statements

  ins2 = nlat/2+MOD(nlat,2)

  ! Find first approximation of the roots of the
  ! Legendre polynomial of degree nlat

  DO jgl = 1, ins2
     z = REAL(4*jgl-1,dp)*api/REAL(4*nlat+2,dp)
     pa(jgl) = COS(z+1.0_dp/(TAN(z)*REAL(8*nlat**2,dp)))
  END DO

  ! Computes roots and weights
  ! Perform the Newton loop
  ! Find 0 of Legendre polynomial with Newton loop

  DO jgl = 1, ins2

     za = pa(jgl)

     DO iter = 1, itemax+1
        zk = 0.0_dp

        ! Newton iteration step

        zkm2 = 1.0_dp
        zkm1 = za
        zx = za
        DO jn = 2, nlat
           zk = (REAL(2*jn-1,dp)*zx*zkm1-REAL(jn-1,dp)*zkm2)/REAL(jn,dp)
           zkm2 = zkm1
           zkm1 = zk
        END DO
        zkm1 = zkm2
        zldn = (REAL(nlat,dp)*(zkm1-zx*zk))/(1.0_dp-zx*zx)
        zmod = -zk/zldn
        zxn = zx+zmod
        zan = zxn

        ! computes weight

        zkm2 = 1.0_dp
        zkm1 = zxn
        zx = zxn
        DO jn = 2,nlat
           zk = (REAL(2*jn-1,dp)*zx*zkm1-REAL(jn-1,dp)*zkm2)/REAL(jn,dp)
           zkm2 = zkm1
           zkm1 = zk
        END DO
        zkm1 = zkm2
        zw = (1.0_dp-zx*zx)/(REAL(nlat*nlat,dp)*zkm1*zkm1)
        za = zan
        IF (ABS(zmod) <= eps) EXIT
     END DO

     pa(jgl) = zan
     pw(jgl) = 2*zw

  ENDDO

!DIR$ IVDEP
!OCL NOVREC

  DO jgl = 1, nlat/2
     isym = nlat-jgl+1
     pa(isym) = -pa(jgl)
     pw(isym) = pw(jgl)
  ENDDO

END SUBROUTINE gauaw


SUBROUTINE check_files (kmax,nlon,nlat)
!----------------------------------!

! check the allowed codes
! read the input file BOT.srv and write a new file BBOT.srv
! with allowed codes

  INTEGER, PARAMETER :: nnlat = 960
  INTEGER, PARAMETER :: nnlon = 480
  INTEGER :: kmax,nlon,nlat
  INTEGER :: iheada(8),iheadb(8)
  INTEGER :: i, ii, jj
  REAL    :: dummy(nnlon,nnlat)

 
  OPEN (unit=10,file='BOTa.srv',form='UNFORMATTED')
  OPEN (unit=11,file='BBOTa.srv',form='UNFORMATTED')
 
  OPEN (unit=20,file='BOTb.srv',form='UNFORMATTED')
  OPEN (unit=21,file='BBOTb.srv',form='UNFORMATTED')
 
  READ(10,END=100) iheada

  nlon       = iheada(5)
  nlat       = iheada(6)
  READ(20,END=100) iheadb

  IF (iheada(5) == iheadb(5) .and. iheada(6) == iheadb(6)) THEN

    nlon       = iheada(5)
    nlat       = iheada(6)
    IF (nlat > nnlat) THEN
      write(6,*) 'number of latitude ',nlat,' is to big!!'
      STOP
    ENDIF
    IF (nlon > nnlon) THEN
      write(6,*) 'number of longitude ',nlon,' is to big!!'
      STOP
    ENDIF
  ELSE
    WRITE (0,*) 'File A longitudes = ',iheada(5),' latitudes = ', iheada(6)
    WRITE (0,*) 'File B longitudes = ',iheadb(5),' latitudes = ', iheadb(6)
    WRITE (0,*) '-----the resolution grid is not the same grid------------'
    STOP
  ENDIF
  WRITE (0,'(a,i4,a,i4)') 'header - longitudes = ', iheada(5), &
      ' latitudes = ', iheada(6)
  REWIND 10
  REWIND 20

  k=0
  DO i = 1, kmax 
    READ(10,END=100) iheada
    READ(20,END=100) iheadb

    IF (iheada(1) .NE. iheadb(1)) THEN
      WRITE (0,*) 'File A code = ',iheada(1)
      WRITE (0,*) 'File B code = ',iheadb(1)
      WRITE (0,*) '-----DIFFERENT codes------------'
      STOP
    END IF
    IF (iheada(1) == 4 ) THEN
        WRITE (0,*) 'code: ', iheada(1),'  Updated to code: ', iheada(1)+256
        iheada(1) = iheada(1)+256 
        iheadb(1) = iheadb(1)+256 
    END IF
    IF (iheada(1) <  34) WRITE (0,*) 'WARNING: Code',iheada(1),' <  34 not supported!!'
    IF (iheada(1) > 260) WRITE (0,*) 'WARNING: Code',iheada(1),' > 260 not supported!!'

    IF (iheada(1) > 33 .AND. iheada(1) < 261) THEN
        k=k+1
        READ(10,END=100) ((dummy(ii,jj),ii=1,nlon),jj=1,nlat)
        WRITE(11) iheada
       
        WRITE(11) ((dummy(ii,jj),ii=1,nlon),jj=1,nlat)
        READ(20,END=100) ((dummy(ii,jj),ii=1,nlon),jj=1,nlat)
        WRITE(21) iheadb
       
        WRITE(21) ((dummy(ii,jj),ii=1,nlon),jj=1,nlat)
    ENDIF

  ENDDO
100 CONTINUE
    write(6,*) 'k,kmax :  ',k,kmax
  
  CLOSE (10)
  CLOSE (11)
  CLOSE (20)
  CLOSE (21)

END SUBROUTINE check_files
