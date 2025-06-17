!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_fastj
!!
!! \brief
!!  This module provides an implementation of the fast-JX code of Bian and Prather (2002)
!!
!! \author Luca Pozzoli (EPFL)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - LP: original version based on O. Wild's code (2003-05)
!!  - LP: upgraded to fastj-X code (2008)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!
!! \limitations
!!  not yet interfaced with HAMMOZ and working
!!
!! \details
!!
!! \bibliographic_references
!!  - Bian HS; Prather MJ, J. Atmos. Chem., 41 (3), 281-296, doi: 10.1023/A: 1014980619462
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TODO:
!! * use some information from ECHAM to determine current solar flux (solar min, solar max, etc.)
!!   Variable zsolf. This would entail to remove use of caldayn which is currently used to estimate
!!   Excentricity of Earth's orbit.
!!

MODULE mo_moz_fastj
  !-----------------------------------------------------------------------
  ! Module containing variables in fastctl namelist and input of fastj data
  !-----------------------------------------------------------------------
  USE mo_kind,           ONLY: dp
  USE mo_math_constants, ONLY: pi
  USE mo_exception,      ONLY: finish
!++ost: not available anymore, change write(nout, to write(*, throughout the code
! USE mo_doctor,         ONLY: nout
!--ost

  IMPLICIT NONE

  PRIVATE

  PUBLIC     :: id_fjx
  PUBLIC     :: setfjx, fjx_interface, jmap, jvn_
!++lp temporary nwv for aerosol   
  PUBLIC     :: nwv
!--lp

  INTEGER :: id_fjx  ! submodel ID
  INTEGER :: nwv = 5 ! temporary  nwv: number of additional wavelengths
!++mgs 20130304 : for debugging purposes
  INTEGER :: error_count = 0
!--mgs
!!  LOGICAL :: lcloud = .FALSE.! T=cloud optical properties are included in the calculation  
!!  LOGICAL :: laero  = .FALSE.! T=aerosols optical properties are included in the calculation
  CHARACTER(LEN=54) :: ratj,spec,scat 

  ! for coupling with HAM aerosols
  REAL(dp), PUBLIC, ALLOCATABLE :: odaer_fj(:,:,:,:)
  REAL(dp), PUBLIC, ALLOCATABLE :: ssaer_fj(:,:,:,:)
  REAL(dp), PUBLIC, ALLOCATABLE :: pp_fj(:,:,:,:,:)

  CHARACTER (len=256)         :: cerr1, cerr2, cerr3

  INTEGER, PARAMETER :: iph = 8           ! unit number to read the Fast-J.2 input files 
  REAL(dp),    PARAMETER :: rad = 6375.e5_dp     ! Radius of earth (cm)     
  REAL(dp),    PARAMETER :: zzht = 5.e5_dp       ! Effective scale height above top of atmosphere (cm)
  REAL(dp),    PARAMETER :: dtaumax = 1._dp    ! Maximum opt.depth above which a new level should be inserted
  REAL(dp),    PARAMETER :: dtausub = 1.0_dp   ! Number of additional levels to add at top of cloud 
  REAL(dp),    PARAMETER :: dsubdiv =10._dp    ! No. of opt.depths at top of cloud requiring subdivision
  REAL(dp),    PARAMETER :: szamax = 98.0_dp   ! Solar zenith angle cut-off, above which to skip calculation

  REAL(dp), PARAMETER :: pi180=pi/180._dp  

!!  INTEGER, PARAMETER :: JVN_= 35 !Number of photolysis reactions in test_chem_Js.dat
  INTEGER, PARAMETER :: JVN_= 58 !Number of photolysis reactions in test_noalias_chem_Js.dat 
  INTEGER, PARAMETER :: W_  = 18 !Number of wavelength bins
  INTEGER, PARAMETER :: X_  = 64 !Number of X-section data sets
  INTEGER, PARAMETER :: A_  = 40 !Number of aerosol/cloud Mie sets 
  INTEGER, PARAMETER :: MX  =  4 !Number of aerosol/cloud types provided by GCM/CTM
  INTEGER, PARAMETER :: N_  =501 !Number of levels in Mie scattering arrays 
  INTEGER, PARAMETER :: M_  =  4 !Number of Gauss points used, must = 4 in fast_JX (no option) 
 
  CHARACTER(len=20) :: TITLAA(A_), TITLUM(33)
  CHARACTER(len=78) :: TITLE0
  CHARACTER(len=7)  :: TITLEJ(X_), TITLEJ2, TITLEJ3
  CHARACTER(len=7)  :: JLABEL(JVN_)    ! fast-J label
  CHARACTER(len=9)  :: JMOZLABEL(JVN_) ! MOZART name
  CHARACTER(len=60) :: JLONLABEL(JVN_) ! longname label (for output in photo stream)
  INTEGER           :: JMAP(JVN_)      ! mapping of fastjx rates onto MOZART indices

  REAL(dp) :: ATAU, ATAU0

  INTEGER :: JIND(JVN_)
  INTEGER :: NJVAL,NRATJ,NW1,NW2,NAA,JTAUMX
  REAL(dp)    :: TANHT
  REAL(dp)    :: WBIN(W_+1),WL(W_),FL(W_),QO2(W_,3),QO3(W_,3),Q1D(W_,3),QQQ(W_,2,X_), &
                 QRAYL(W_+1),TQQ(3,X_)
  REAL(dp)    :: WAA(5,A_),QAA(5,A_),PAA(8,5,A_)
  REAL(dp)    :: RAA(A_),SAA(5,A_),DAA(A_)
  REAL(dp)    :: JFACTA(JVN_)


   contains
!-----------------------------------------------------------------------
      subroutine inphot_fjx
!-----------------------------------------------------------------------
!  Routine to initialise photolysis rate data, called directly from the
!  cinit routine in ASAD. Currently use it to read the JPL spectral data
!  and standard O3 and T profiles and to set the appropriate reaction index.
!-----------------------------------------------------------------------

     implicit none

! Read in fast-J X-sections (spectral data) <<<<<<<<<<<<<< new fast-JX
     call RD_XXX(IPH,spec)

! Read in aerosol/cloud scattering data <<<<<<<<<<<<<<<<<< new fast-JX
     call RD_MIE(IPH,scat)

! Read in labels of photolysis rates required	>>>>> keyed to users chem code
!   this is a tranfer map from the J's automatically calculated in fast-JX
!   onto the names and order in the users chemistry code
     call RD_JS(IPH,ratj)
     
      return
      end subroutine inphot_fjx

!-----------------------------------------------------------------------
     subroutine RD_MIE(NJ1,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NJ1      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     NK       Number of wavelengths at which functions supplied (set as 4)
!     WAA      Wavelengths for the NK supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NJ1
      character(*), intent(in) ::  NAMFIL
      integer  I, J, K

      open (NJ1,FILE=NAMFIL,status='old',form='formatted')

      read (NJ1,'(i2,a78)') NAA,TITLE0
        if (NAA .gt. A_) then
         WRITE(cerr1,*) NAA
         WRITE(cerr2,*) A_
         CALL finish('mo_fastj, RD_MIE:', &
         'NAA too high='//TRIM(cerr1)// &
         ' but A_='//TRIM(cerr2))
        endif
      read (NJ1,'(5x,i5,2f10.5)') JTAUMX,ATAU,ATAU0
      write(*,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX
      read (NJ1,*)
      do J = 1,NAA
          read (NJ1,'(3x,a20,32x,f5.3,15x,f5.3)') &
                TITLAA(J),RAA(J),DAA(J)
          write (*,'(3x,a20,32x,f5.3,15x,f5.3)')  &
                TITLAA(J),RAA(J),DAA(J)
        do K = 1,5     ! ver 6.0 extend to 5 ref wavelengths for mie-scat data
          read (NJ1,'(f4.0,f7.4,f7.4,7f6.3)') &
            WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          write (*,'(f4.0,f7.4,f7.4,7f6.3)')  &
            WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1.d0
        enddo
      enddo
 
      close(NJ1)

      write(*,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):' &
                  ,(WAA(K,1),K=1,5)
      write(*,*) TITLE0
      do J=1,NAA
      write(*,'(i3,1x,a8,7f8.3)') &
                   J,TITLAA(J),RAA(J),DAA(J),(QAA(K,J),K=1,5)
      enddo

      return
      end subroutine rd_mie

!-----------------------------------------------------------------------
     subroutine RD_XXX(NJ1,NAMFIL)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh parameters, 
!      T-dependent X-sections. 

!  >>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<

!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (j2_spec.dat) >> j2 for fast-J2
!     NJ1      Channel number for reading data file
!
!     NJVAL    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NJ1
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, IW, NQQQ, NWWW

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data J-ver8.3------------------
!	  note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects read in
!	  for 2005a data, NJVAL = 62 (including a spare XXXX) and 
!	       NQQQ = 64 so that 4 wavelength datasets read in for acetone
!	  note NQQQ is not used outside this subroutine!

      open (NJ1,FILE=NAMFIL,status='old',form='formatted')
      read (NJ1,100) TITLE0
      read (NJ1,101) NJVAL,NQQQ, NWWW,NW1,NW2
      if (NJVAL.gt.X_ .or. NQQQ.gt.X_) then
         WRITE(cerr1,*) NJVAL
         WRITE(cerr2,*) NQQQ
         WRITE(cerr3,*) X_
         CALL finish('mo_fastj, RD_XXX:', &
         'NJVALS or NQQQ too high='//TRIM(cerr1)// &
         ',NQQQ='//TRIM(cerr2)//' but X_='//TRIM(cerr2))
      endif
      write(*,'(1X,A)') TITLE0
!----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      read (NJ1,102) (WL(IW),IW=1,NWWW)
      read (NJ1,102) (FL(IW),IW=1,NWWW)
      read (NJ1,102) (QRAYL(IW),IW=1,NWWW)

!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      read (NJ1,103) TITLEJ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

      read (NJ1,103) TITLEJ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

      read (NJ1,103) TITLEJ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NJ1,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NJ1,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

      do J = 1,3
        write(*,200) TITLEJ(J),(TQQ(I,J),I=1,3)
      enddo

!---Read remaining species:  X-sections at 2 T_s
      do J = 4,NQQQ
        read (NJ1,103) TITLEJ(J),TQQ(1,J),(QQQ(IW,1,J),IW=1,NWWW)
        read (NJ1,103) TITLEJ2,  TQQ(2,J),(QQQ(IW,2,J),IW=1,NWWW)
        write(*,200) TITLEJ(J),(TQQ(I,J),I=1,2)
      enddo

!  Reset the titles for NJVAL-1 & NJVAL to be the two acetone J_s
!   61: C3H6O  = Acet-a     (CH3CO + CH3) 
!   62: Q2-Ac  = Acet-b     (CH3 + CO + CH3)

      TITLEJ(NJVAL-1) = 'Acet-a'
      TITLEJ(NJVAL)   = 'Acet-b'
      
      close(NJ1)
      
  100 format(a)
  101 format(10x,5i5)
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
  103 format(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(' Number of x-sections supplied to Fast-JX: ',i3,/, &
             ' Maximum number allowed (X_) only set to: ',i3,    &
             ' - increase in cmn_jv.f')

      return
      end subroutine rd_xxx

!-----------------------------------------------------------------------
      subroutine RD_JS(NJ1,NAMFIL)
!-----------------------------------------------------------------------
!  Reread the chem_Js.dat file to map photolysis rate to reaction
!  Read in quantum yield 'jfacta' and fastj2 label 'jlabel'
!-----------------------------------------------------------------------
!
!     jfacta	Quantum yield (or multiplication factor) for photolysis
!     jlabel	Reference label identifying appropriate J-value to use
!     ipr	Photolysis reaction counter - should total 'JVN_'
!
!-----------------------------------------------------------------------
      implicit none
!
      integer, intent(in) ::  NJ1
      character(*), intent(in) ::  NAMFIL

      integer  IPR, I, J, K
      character*120 CLINE
!
! Reread the chem_Js.dat file to map photolysis rate to reaction
!                     Read in quantum yield jfacta and fastj2 label jlabel
      IPR = 0
      open (NJ1,file=NAMFIL,status='old',form='formatted')
 10   read (NJ1,'(A)',err=20)  CLINE
      if (IPR .eq. JVN_) goto 20

      if (CLINE(2:5).eq.'9999') then
        go to 20
      elseif (CLINE(1:1).eq.'#') then
        go to 10
      elseif (CLINE(5:5).eq.'$') then
        go to 10
      else
        IPR = IPR+1
        read (CLINE(101:105),'(F5.1)') JFACTA(IPR)
        read (CLINE(108:114),'(A7)')   JLABEL(IPR)
        read (CLINE( 3:11),'(A9)')     JMOZLABEL(IPR)     !++mgs - new
        read (CLINE(19:78),'(A60)')    JLONLABEL(IPR)     !++mgs - new
        JFACTA(IPR) = JFACTA(IPR)/100.d0
        go to 10
      endif
 20   close(NJ1)

      NRATJ = IPR

!-----------------------------------------------------------------------
!  compare Xsections titles with J-values listed in chem code (jratd.dat)
!  map J-values needed for chemistry (chem_Js.dat) onto the fast-JX rates
!  >>>>>>>>>>>>>>>>current code revised to JPL-02 ver 8.5 (5/05)<<<<<<<<<
!	   >>>this must now follow the read in of Xsects, etc<<<
!-----------------------------------------------------------------------

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do J = 1,JVN_
        JIND(J) = 0
      enddo
      do J = 1,NJVAL
      do K = 1,NRATJ
        if (JLABEL(K) .eq. TITLEJ(J)) JIND(K)=J
      enddo
      enddo

      write(*,'(a,i4,a)') ' Photochemistry Scheme with ',IPR,' J-values'
      do K=1,NRATJ
        J = JIND(K)
        if (J.eq.0) then
         write(*,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K), &
              ' has no mapping onto onto fast-JX'
        else
         write(*,'(i5,a9,f6.2,a,i4,a9)') K,JLABEL(K),JFACTA(K), &
              ' mapped onto fast-JX:',J,TITLEJ(J)
        endif
      enddo  

      return
      end subroutine rd_js

!--------------------------------------------------------------------------------------------
  SUBROUTINE setfjx(lng_indexer)

    ! *setfjx* namelist control of fastjx module: control flags and
    !          definition of output
    ! 
    ! Authors:
    ! --------
    ! Luca Pozzoli, EPFL                       01/2004
    !
    ! *setfjx* is called from *init_submodels_1* in *mo_submodel_interface*
    !

    USE mo_mpi,         ONLY: p_parallel, p_parallel_io, p_bcast, p_io
    USE mo_namelist,    ONLY: open_nml, position_nml, POSITIONED
    USE mo_decomposition, ONLY: ldc => local_decomposition 
    USE mo_exception,   ONLY: finish, message     !++mgs
    USE mo_string_utls, ONLY: st1_in_st2_proof, st1_in_st2_idx
    USE mo_submodel,    ONLY: lmoz, lham
!++mgs 20130228
    USE mo_moz_mods,    ONLY: phtcnt, pht_alias_lst, rxt_tag_lst
!--mgs
    USE mo_moz_photo,   ONLY: nphotoreac, clonglabels, jdefined
!--mgs
!!    USE mo_radiation_parameters,   ONLY: io3
!!    USE mo_moz,         ONLY: lfastj, idt_o3

    IMPLICIT NONE

!!    INCLUDE 'moz_fastjctl.inc'

    !--- Dummy arguments

    integer, intent(inout) :: lng_indexer(:)        ! FJX to reaction_index mapping

    !--- Local variables:

    INTEGER :: i, ierr, ivar, iphotoextra, inml, iunit
!!mgs!! QUICK FIX
    INTEGER, SAVE :: Nwv_opt = 1

    !--- 0) set default values
    ! input files
    ratj = 'chem_Js.dat'
    spec = 'FJX_spec.dat'
    scat = 'FJX_scat.dat'
 
    JMAP(:) = -1

    !--- 1) Read namelist:
    
!!    IF (p_parallel_io) THEN
!!       inml = open_nml('namelist.echam')
!!       iunit = position_nml ('MOZ_FASTJCTL', inml, status=ierr)
!!       SELECT CASE (ierr)
!!       CASE (POSITIONED)
!!          READ (iunit, moz_fastjctl)
!!       END SELECT

    !--- check if it makes sense to run fastj

!++mgs: changed from finish to message - turn off fastj automatically for convenience
!!     IF (lfastj .AND. .NOT. lmoz) CALL finish('setfjx',  &
!!          'Fastjx activated with MOZ turned off. Maybe possible in the future, but not now!')
!!       IF (lfastj .AND. .NOT. lmoz) THEN
!!          CALL message('setfjx', 'Fastjx was activated with MOZ turned off. Switching off fast J now!' )
!!          lfastj = .false.
!!       END IF
!--mgs
!!    ENDIF
    
    !--- 2) Broadcast over processors:

!!    IF (p_parallel) THEN
!!     CALL p_bcast (lcloud,        p_io)
!!     CALL p_bcast (laero,         p_io)
!!     CALL p_bcast (ratj,          p_io)
!!     CALL p_bcast (spec,          p_io)
!!     CALL p_bcast (scat,          p_io)     
!!    END IF
   
      ! Read Fast-J files
      write(*,*) '--- FJX: reading fastj input files ...'
      CALL inphot_fjx

      ! test JMOZLABELS for validity
!      IF (.NOT. st1_in_st2_proof( JMOZLABEL(1:JVN_), cphotovars(1:nphotovars) )) THEN
!         CALL finish ( 'setfjx', 'MOZ label '// &
!                       JMOZLABEL(ierr)// ' in chem_Js.dat not defined by MOZART preprocessor! (set_sim_dat)' )
!      END IF
      ! map j names, set long names for photo stream output and mark rates by fastj
      DO ivar = 1,JVN_
         IF (.NOT. st1_in_st2_proof( JMOZLABEL(ivar), rxt_tag_lst(1:phtcnt) )) THEN
            CALL finish ( 'setfjx', 'MOZ label '// &
                         JMOZLABEL(ivar)// ' in chem_Js.dat not defined by MOZART preprocessor! (set_sim_dat)' )
         END IF
         ierr = st1_in_st2_idx( JMOZLABEL(ivar), rxt_tag_lst(1:phtcnt) )
         IF ( ierr /= 0 ) THEN
           clonglabels(ierr) = JLONLABEL(ivar)
           JMAP(ivar) = ierr
           jdefined(ierr) = .true.
         END IF
      END DO
      ! now do the reverse mapping
      ! first reactions explicitly calculated by fastj
      DO ivar = 1,JVN_
        lng_indexer(JMAP(ivar)) = ivar
      END DO
      ! now add reactions with pht_alias
      DO ivar = 1,phtcnt
        IF ( pht_alias_lst(ivar,2) /= ' ' ) THEN
          ierr = st1_in_st2_idx( pht_alias_lst(ivar,2), rxt_tag_lst(1:phtcnt) )
          IF ( ierr /= 0 ) THEN
              DO i = 1,JVN_
                 IF (JMAP(i) == ierr) THEN
                    lng_indexer(ivar) = i
                    cycle
                 END IF
              END DO
          END IF
          IF( pht_alias_lst(ivar,2) == 'userdefined' ) THEN
            lng_indexer(ivar) = -1
          END IF
        END IF
      END DO
!!    ! safety check
!!    IF ( any(lng_indexer(:) > JVN_) ) THEN
!!       write(0,*) "**** lng_indexer with values > JVN_!", lng_indexer
!!       CALL finish ("mo_moz_fastj", "lng_indexer")
!!    END IF
      IF (p_parallel_io) THEN
        write(*,*) 'Info: FJX: JMAP = ', JMAP
        write(*,*) 'Info: FJX: lng_indexer = ',lng_indexer
      END IF

!!      IF (p_parallel_io) THEN
!!        write(*,*) '----------------------------------------------------------------------'
!!        write(*,*) '--- FJX                 =',lfastj
!!        write(*,*) '--- FJX clouds          =',lcloud
!!        write(*,*) '--- FJX aerosols        =',laero
!!        write(*,*) '--- FJX index file      =',ratj
!!        write(*,*) '--- FJX spectral file   =',spec
!!        write(*,*) '--- FJX scattering data =',scat
!!        write(*,*) '--- index   moz-index   MOZ label     FJX label     long name'
!!        DO ivar = 1,JVN_
!!           write(*,'(A4,i5," : ",i8,4x,A9,3x,A8,3x,A60)') '--- ', ivar, JMAP(ivar), &
!!                      JMOZLABEL(ivar), JLABEL(ivar), JLONLABEL(ivar)
!!        END DO
!--mgs
!        IF (io3 == 1)   write(*,*) '--- FJX use chemical ozone field'
!        IF (io3 == 1 .AND. idt_o3 <= 0) THEN
!           write(*,*) '--- *** io3==1 but idt_o3 == 0 !! ***'
!           write(0,*) 'setfjx: --- *** io3==1 but idt_o3 == 0 !! ***'
!        END IF

!++lp move out (hammoz module)
!      !-- allocate variables for HAM-fastJ coupling
!      IF (lham) THEN
!         IF ( .NOT. ALLOCATED(odaer_fj)) ALLOCATE(odaer_fj(ldc% nglon,ldc% nlev,Nwv_opt,nclass))
!         IF ( .NOT. ALLOCATED(ssaer_fj)) ALLOCATE(ssaer_fj(ldc% nglon,ldc% nlev,Nwv_opt,nclass))
!         IF ( .NOT. ALLOCATED(pp_fj))    ALLOCATE(pp_fj(ldc% nglon,ldc% nlev,Nwv_opt,nclass,8)) 
!         ! 8 = number of expansion coefficients; to be parameterized
!      END IF

  END SUBROUTINE setfjx

!!baustelle!! double-check use of pmid and pfull: mid=39 levels, full=40 levels!

  SUBROUTINE fjx_interface(plonl,klev,krow,  & !ECHAM indices
                           caldayn,           & !calendar day
                           palbedo,           & !surface albedo
                           pmid,              & !half level pressure
                           pfull,             & !full level pressure
                           pdel,              & !layer thickness
                           tfld,              & !full level temperature
                           tsurf,             & !surface temperature
                           pcwat,             & !totel water content (water+ice)
                           pcice,             & !ice content
                           cldfr,             & !cloud fraction
                           pcdnc,             & !cdnc field
                           po3,               & !ozone vmr
                           pjrates            ) !photolysis values
  !
  !**** *fjx_interface* serves as interface between Fast-J.X and ECHAM5
  !
  !      Authors:
  !      --------
  !      Luca Pozzoli     Feb. 2004

  !**    Interface:
  !      ----------
  !      *fjx_interface*  is called from   *chemdr*
  !
  !      Method:
  !      -------
  !      
  !      This subroutine computes various quantities, partly the same as in radiation.
  !      As a final step, the order of all input variables needed for FJX is reversed
  !      (bottom->top).
  !

  USE mo_math_constants,   ONLY: pi
  USE mo_physical_constants, ONLY: grav, rd, rhoh2o
  USE mo_echam_cloud_params, ONLY: ceffmin, ceffmax
!++ost: amu0 is not needed in the subroutine, delete?
!++ost: was amu0 before, attention: krow!
  USE mo_geoloc,      ONLY: amu0_x                !cosine of the solar zenith angle
!--ost
!--ost
  USE mo_debugs
!!  USE mo_mpi, only : p_pe      ! ### DEBUG ##

IMPLICIT NONE
!!++lp  INTEGER, PARAMETER :: nwv = 1     !!mgs!! QUICK FIX
! --- input arrays -------------------------------
  INTEGER, INTENT(in)     :: plonl, klev, krow

  REAL(dp), INTENT(in)    :: caldayn,                &
                             palbedo(plonl),        &
                             pmid(plonl,klev),     &
                             pfull(plonl,klev+1),       &
                             pdel(plonl,klev),      &
                             tfld(plonl,klev),       &
                             tsurf(plonl),          &
                             pcwat(plonl,klev),     &
                             pcice(plonl,klev),     &
                             cldfr(plonl,klev),     &
                             pcdnc(plonl,klev),     &
                             po3(plonl,klev)       

  REAL(dp), INTENT(out)   :: pjrates(plonl,klev,jvn_)
!-Local variables---------------------------------------------
  REAL(dp)    :: zht(plonl,klev+1)
  REAL(dp)    :: zcwat(plonl,klev)    ! liquid water only
  REAL(dp)    :: zcice(plonl,klev)    ! ice
  REAL(dp)    :: zfac(plonl,klev)     ! pmid/tmid

  REAL(dp)    :: zwgkg(plonl,klev),  &
                 zlwgkg(plonl,klev), &
                 zlwp(plonl,klev),   &
                 zlwc(plonl,klev),   &
                 ziwp(plonl,klev),   &
                 ziwc(plonl,klev),   &
                 zradip(plonl,klev), &
                 zradlp(plonl,klev), &
                 zho3(plonl,klev+1), &
                 zkap(plonl),        &
                 zrex,                &
                 zref,                &
                 zsolf

  REAL(dp)    :: cldfr_fj(plonl,klev),     &
                 od_fj(plonl,klev+1,nwv),  &
!!                 od_fj(plonl,klev+1,8),  &        !!mgs!! QUICK FIX
                 od600_fj(plonl,klev+1),   &
                 ssa_fj(plonl,klev+1,nwv), &
                 sleg_fj(plonl,klev+1,nwv,8)

  REAL(dp)    :: jrates_fj(plonl,klev,jvn_)

  INTEGER     :: jk,jl,i,j

  LOGICAL     :: locldlyr(plonl,klev)  ! indicate presence of cloud

! write(0,*) "##### starting fjx_interface, p_pe = ", p_pe
!-- Initialize
   od600_fj(:,:)=0._dp
   od_fj(:,:,:)=0._dp
   ssa_fj(:,:,:)=0._dp
   sleg_fj(:,:,:,:)=0._dp
   jrates_fj(:,:,:) = 0._dp
!++lp convert O3 vmr to mmr
!!   po3(:,:)=po3(:,:)*(47.9_dp/28.97_dp)

!++mgs: changed from pfull to pmid (appeared to be an error in argument list ??
!! zfac(:,:)   = pfull(:,:) / tfld(:,:) / rd
   zfac(1:plonl,1:klev)   = pmid(1:plonl,1:klev) / tfld(1:plonl,1:klev) / rd
!--mgs

!-- calculate LWP, IWP, and the effecitved radius of cloud droplets and ice crystals
!   needed to calculate the water and ice cloud optical depth.
!   Calculated as in rad_int

   ! cloud liquid water
   zcwat(1:plonl,1:klev) = MAX(pcwat(1:plonl,1:klev)-pcice(1:plonl,1:klev), 0._dp)
   zcice(1:plonl,1:klev) = MAX(pcice(1:plonl,1:klev), 0._dp)

   ! Clear/cloudy flag
   WHERE (cldfr(1:plonl,1:klev)>EPSILON(1._dp))
      locldlyr(1:plonl,1:klev)=.true.
   ELSEWHERE
      locldlyr(1:plonl,1:klev)=.false.
   END WHERE

   ! Secure cloud fraction
   cldfr_fj(1:plonl,1:klev)=MAX(cldfr(1:plonl,1:klev),EPSILON(1._dp))

   ! Specific ice water content, g/kg
   WHERE (locldlyr(1:plonl,1:klev))
      zwgkg(1:plonl,1:klev)=zcice(1:plonl,1:klev)*1000._dp/cldfr_fj(1:plonl,1:klev)
   ELSEWHERE
      zwgkg(1:plonl,1:klev)=0._dp
   END WHERE

   ! Ice water content per volume g/m3
   ziwc(1:plonl,1:klev)=zwgkg(1:plonl,1:klev)*zfac(1:plonl,1:klev)

   ! Ice water path, g/m2
   ziwp(1:plonl,1:klev)=zwgkg(1:plonl,1:klev)*pdel(1:plonl,1:klev)/grav

   ! Specific liquid water content, g/kg
   WHERE (locldlyr(1:plonl,1:klev))
    zlwgkg(1:plonl,1:klev)=zcwat(1:plonl,1:klev)*1000._dp/cldfr_fj(1:plonl,1:klev)
   ELSEWHERE
    zlwgkg(1:plonl,1:klev)=0._dp
   END WHERE

   ! Liquid water content per volume, g/m3
   zlwc(1:plonl,1:klev)=zlwgkg(1:plonl,1:klev)*zfac(1:plonl,1:klev)

   ! Liquid water path, g/m2
   zlwp(1:plonl,1:klev)=zlwgkg(1:plonl,1:klev)*pdel(1:plonl,1:klev)/grav

   ! Boucher/Lohmann (1995) and Moss et al. (1995)
   ! Effective radii for ice particles [um]:

    !look at mo_cloud_optics, changed ceffmin
    zradip(1:plonl,1:klev)=MAX(ceffmin,MIN(ceffmax,83.8_dp*ziwc(1:plonl,1:klev)**0.216_dp))

    ! - liquid
!>>>mgs: the following is erroneous!
!!  WHERE (loland.AND.(.NOT.loglac))
!!     zkap=1.143_dp
!!  ELSEWHERE
!!     zkap=1.077_dp
!!  END WHERE
    zkap(1:plonl) = 1.100_dp !look at mo_cloud_optics
!<<<
    zrex=1._dp/3._dp
    zref=1.e6_dp*(3.e-9_dp/(4._dp*pi*rhoh2o))**zrex
!++lp commented out, needed for cloud optical depth calc.
    DO jk = 1, klev
      DO jl = 1, plonl
        !IF (ncd_activ>0) THEN
        !  !changed zkap(jl) to take dispersion effect (Peng&Lohmann, GRL, 2003) into account
        !  zkap(jl)=0.00045_dp*pcdnc(jl,jk)+1.18_dp
        !ENDIF
        zradlp(jl,jk) = zref*zkap(jl)*(zlwc(jl,jk)/pcdnc(jl,jk))**zrex
      END DO
    END DO
    zradlp(1:plonl,1:klev)=MAX(4._dp,MIN(24._dp,zradlp(1:plonl,1:klev)))
!--lp

! write(0,*) "##### fjx_interface: th interpolation, p_pe = ", p_pe
   ! temperature at half levels
   DO jk=2,klev
!++mgs: I don't understand this formula. Should be linear interpolation, no?
! changed to   Th2 = Tm1 + (Tm2-Tm1)/(pm2-pm1)*(ph2-pm1)
       zht(:,jk) = tfld(:,jk-1) + ( tfld(:,jk) - tfld(:,jk-1) )   &
                                  / ( pmid(:,jk) - pmid(:,jk-1) )   &
                                  * ( pfull(:,jk) - pmid(:,jk-1) )
!!!    zht(:,jk)   = ( tfld(:,jk-1)*pfull(:,jk-1)*(pfull(:,jk)-pmid(:,jk)  )   &
!!!                   + tfld(:,jk)  *pfull(:,jk)  *(pmid(:,jk)-pfull(:,jk-1)) ) &
!!!                  /(             pmid(:,jk)  *(pfull(:,jk)-pfull(:,jk-1)) )
!--mgs
   END DO
   zht(:,klev+1)=tsurf(:) !tsurf is temp2 passed from physc but the values are not correct
!++lp pfull(:,1)=0 replace with 0.1
!!   zht(:,1)=tfld(:,1)-pfull(:,1)*(tfld(:,1)-zht(:,2))/(pfull(:,1)-pmid(:,2))
   zht(:,1)=tfld(:,1)-0.1_dp*(tfld(:,1)-zht(:,2))/(0.1_dp-pmid(:,2))
!!
   !solar flux factor
   zsolf  = 1._dp-(0.034_dp*cos(caldayn-186._dp)*2._dp*pi/365._dp)

! calculate the ozone profile at half level
   DO jk = 2,klev
!++mgs: same as above for zht...
      zho3(:,jk) = po3(:,jk-1) + ( po3(:,jk) - po3(:,jk-1) )   &
                                 / ( pmid(:,jk) - pmid(:,jk-1) )   &
                                 * ( pfull(:,jk) - pmid(:,jk-1) )
!!!  zho3(:,jk) = (po3(:,jk-1)*pfull(:,jk-1)*(pfull(:,jk)-pmid(:,jk))   &
!!!              + po3(:,jk)*pfull(:,jk)*(pmid(:,jk)-pfull(:,jk-1)))/   &
!!!                (pmid(:,jk)*(pfull(:,jk)-pfull(:,jk-1)))
!--mgs
   ENDDO

!++lp pfull(:,1)=0 replace with 0.1
!!   zho3(:,1) =po3(:,1)+(pmid(:,1)-pfull(:,1))*(po3(:,1)-zho3(:,2))/(pfull(:,1)-pmid(:,2))
   zho3(:,1) =po3(:,1)+(pmid(:,1)-0.1_dp)*(po3(:,1)-zho3(:,2))/(0.1_dp-pmid(:,2))
!!
!++mgs: the following statement caused a range error. Confusion pfull/pmid??
! probably safest to set surface cocnentration to lowest mid-level value anyhow (???)
!! zho3(:,klev+1) =po3(:,klev)+(pfull(:,klev)-pmid(:,klev+1))*(po3(:,klev)-zho3(:,klev))/ &
!!                 (pmid(:,klev)-pfull(:,klev))
   zho3(:,klev+1) =po3(:,klev)
!--mgs
   !!zho3(:,:) = 1.e-30_dp
   zho3(:,:) = MAX(zho3(:,:),1.e-30_dp)
   !
   ! ===============================================
   ! reverse levs for fastj   !!
   ! ===============================================
   !

   !water cloud optical depth
   CALL OPTICDW(plonl,klev,od_fj(:,1:klev,:),ssa_fj(:,1:klev,:), &
                sleg_fj(:,1:klev,:,:),    &
                vrev2(zlwp),vrev2(zradlp),vrev2(cldfr) )
    
   !ice cloud optical depth
   CALL OPTICDI(plonl,klev,od_fj(:,1:klev,:),ssa_fj(:,1:klev,:), &
                sleg_fj(:,1:klev,:,:),    &
                vrev2(ziwp),vrev2(zradip),vrev2(cldfr) )

    !use OD of clouds (not aerosols) at 600 nm to determine added layers
    od600_fj(1:plonl,1:klev)=od_fj(1:plonl,1:klev,4)

    !add aerosol optical depth
!++lp let couple aerosols later
!!    CALL OPTICDA(plonl,klev,od_fj(:,1:klev,:),ssa_fj(:,1:klev,:),sleg_fj(:,1:klev,:,:))

    do J=1,nwv
     do jl=1,klev
      do jk=1,plonl
       if (od_fj(jk,jl,J) .gt. 0._dp) then
        do i=1,8
         sleg_fj(jk,jl,j,i)=sleg_fj(jk,jl,j,i)/od_fj(jk,jl,j)
        enddo
       endif
      enddo
     enddo
    enddo


!- call Fast-J.X--------------------------------------------

! write(0,*) "##### fjx_interface: debug output, p_pe = ", p_pe
   !zdf1(:,krow)=tsurf(:) 
   !zdf2(:,krow)=amu0_x(:,krow) !checked OK
   !zdf3(:,krow)=palbedo !checked OK
   !ddf1(:,:,krow)=vrev2(pfull(:,2:klev+1)) !checked OK
   !ddf2(:,:,krow)=vrev2(pfull(:,1:klev)) !checked OK
   !ddf3(:,:,krow)=vrev2(tfld) !checked OK
   !ddf4(:,:,krow)=vrev2(zht(:,2:klev+1)) !checked OK
   !ddf5(:,:,krow)=vrev2(zht(:,1:klev))   !checked OK
   !ddf6(:,:,krow)=vrev2(zho3(:,2:klev+1))
   !ddf7(:,:,krow)=vrev2(zho3(:,1:klev))
   !ddf8(:,:,krow)=vrev2(pmid(:,1:klev)) !checked OK
  CALL fjx(plonl,klev,krow,   &
           amu0_x(:,krow),&      !++ost: was amu0 before, added krow!, not needed in fjx, delete?
           palbedo,       & 
           vrev2(pfull),  &      !++mgs: changed from pmid !
           vrev2(tfld),   &
           vrev2(zht),    &
           od_fj,         & 
           od600_fj,      & 
           ssa_fj,        & 
           sleg_fj,       & 
           vrev2(zho3),   &
           zsolf,          &
           jrates_fj        )
 
  !-- reverse photorates
  DO J=1,jvn_
     pjrates(:,:,J) = vrev2(jrates_fj(:,:,J))
  END DO

  RETURN 

  CONTAINS

      FUNCTION vrev2(p)

        IMPLICIT NONE

        ! inverts array p(plonl,klev) in the second index

        REAL(dp), INTENT(in), DIMENSION(:,:)     :: p
        REAL(dp), DIMENSION(SIZE(p,1),SIZE(p,2)) :: vrev2

        INTEGER :: klev, jk

        klev=SIZE(p,2)

        DO jk=1,klev
           vrev2(:,jk)=p(:,klev+1-jk)
        END DO

      END FUNCTION vrev2

  END SUBROUTINE fjx_interface


  SUBROUTINE FJX (plonl,klev,krow,  &
                u0,            &                   !++ost, amu0 -> u0 not needed in subroutine, delete?
                alb_s,         &
                h_pres,        &
                t_fld,h_temp,  &
                od, 	       &
                od600,	       &
                ssalb,	       &
                sleg,	       &
                ho3,           &
		solf,          &
                zpj            )

!c  >>>>>>>>>>>>>>>>current code revised to JX ver 6.2 (6/08)<<<<<<<<<<<<
!c
!c version 6.2 corrects a long-standing problem at SZA > 89 degrees.
!c   In prior versions the ray-tracing of the path (and air-mass functions)
!c   back to the sun was done at the edges of the CTM layers (it was developed
!c   for the grid-point J-value code at Harvard/GISS/UCI).  This left the 
!c   interpolation to the mid-layer (needed for J's) open.  The prior method
!c   gave irregular fluctuations in the direct solar beam at mid-layer for
!c   large SZA > 88.  This is now corrected with exact ray-tracing from
!c   the mid-pt of each CTM layer.  For small SZA, there is no effective 
!c   difference, for large SZA, results could be erratic.
!c
!c   6.2 fix should be easy if you have migrated to v6.1, else some minor
!c   caution may be needed:
!c      replace sub SPHERE with SPHERE2, AMF2 report factors for mid and egdes.
!c      replace sub OPMIE with new OPMIE, this uses the new AMF2 correctly.
!c      replace sub PHOTOJ with new PHOTOJ, this just hands off AMF2 from 
!c            SPHERE2 to OPMIE.
!c
!c version 6.1 adds
!c      6.1b simplifies calling sequences feeds solar factor, albedo, to PHOTOJ
!c         and read LAT, LNG directly.  No substantive changes.
!c      new read-in of scat data for clouds/aerosols to allow for UMich data
!c      This has required substantial rewrite of some of the core subroutines:
!c         OPMIE is now called for each wavelength and without aersol/cloud data
!c              all subs below OPMIE are unchanged
!c         OPTICD & OPTICM are new subs to convert path (g/m2) to OD and phase fn
!c              D is std UCI scat data (re-ordered for clouds 1st)
!c              M is U Michigan data tables for aerosols, includes Rel Hum effect
!c         PHOTOJ now assembles the aerosol data (better for CTM implementation)
!c      This version can reproduce earlier versions exactly, but the test input 
!c         is changed from OD and NDX to PATH (g/m2) and NDX.
!c version 6.0 adds
!c      new 200-nm scattering data so that stratospheric aerosols can be done!
!c version 5.7
!c     adds the new flux diagnostics (including heating rates)
!c        accurate fluxes for spherical atmos and SZA > 90 !
!c     recommend geometric delta-tau factor from 1.18 to 1.12 for more accurate
!c        heating rates (but more layers!)
!c     tuned and corrected to be almost flux conserving (1.e-5), except
!c        deep clouds, where diffusive flux is created (1.e-4)
!c     still needs to return to the original 1970-code for the block-tri solution
!c        after extensive profiling with F95 and 'modern' versions
!c        it was found that they are much more expensive!!!
!c     corrects typo in JAC(2000) fast-J paper on I+ (reflected from l.b.):
!c        I+(lb) = refl/(1+refl) * (4*Integ[j(lb)*mu*dmu] + mu0*Fdirect(lb))
!c version 5.6 adds
!c      clean up problems with thick clouds does correct solar attenuation
!c	  into cloud sub-layers and into the mid-point of the CTM level
!c	New calculated upward and downward FLUXES at each wavelength at TOP/BOT
!c	Correct deposition of solar flux in each CTM layer (spherical)
!c	  awaits new diagnostics of the h's for heating rates.
!c	back to old matrix solver (UCI blocksolver and matinv-4)
!c version 5.5 adds
!c	new code for generating and solving the block tri-diagonal scattering
!c	     problem.  Uses single call to GEM and general 4x4 block-tri solver.
!c version 5.3c adds
!c	calculates reflected UV-vis solar energy (relative to 1.0)
!c	new solar spectrum (J-O2 increases in strat by 10%, J-NO by 15+%)
!c
!c version 5.3b changes include:
!c	new data files for specral Xsection and mie-scattering.
!c	add sub-layers (JXTRA) to thick cloud/aerosol layers,
!c	     sets up log-spaced sub-layers of increasing thickness ATAU
!c	correction 'b' does massive clean up of the linking code,
!c	     now the only subroutine that has access to CTM arrays is PHOTOJ
!c	     Also, the access to the cmn_JVdat.f is 'read-only' after init.
!c	     This should enable safe openMP/MPI coding.
!c
!c common files and what they mean:
!c   parm_CTM.f  dimensions & params for code (CTM and fast-JX)
!c   parm_MIE.f  dimensions for mie code variables.
!c   cmn_metdat.f  CTM 3-D arrays, time of day, grid,  etc.
!c   cmn_JVdat.f   Xsects, Mie, etc., (initialized and then read-only)
!c
!c<<<<<<<<<<<<<<<<<<<<<begin CTM-specific subroutines<<<<<<<<<<<<<<<<<<<<
!c subroutines:
!c
!c     SET_ATM(GMTAU)
!c	     set ups atmosphere (p,T,O3,airmass, etc) for time GMTAU
!c		COMMON BLOCKS: cmn_metdat.f
!c
!c     SET_AER(GMTAU)
!c		set ups aerosols for time GMTAU = DUMMY
!c			    
!c     SET_CLD(GMTAU)
!c	     set ups clouds for time GMTAU = DUMMY
!c
!c     INPHOT:  Init. photolysis rate, called once by CHMSET
!c		COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f
!c		Input files: ECT42_grid.dat
!c
!c     RD_JS(NJ1,NAMFIL):  Read labels of photo. rates, called once by INPHOT.
!c		COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f
!c		Input files: chem_Js.dat
!c
!c     RD_PROF(NJ2,NAMFIL):  Read T & O3 climatology, called once by INPHOT.
!c		COMMON BLOCKS: cmn_metdat.f
!c		Input files: atmos_std.dat
!c
!c     SET_CLD0(TINIT,ODINIT,ODINDX,ALBEDO)
!c		Initialize cloud and surface properties, called by MAIN.
!c		COMMON BLOCKS: cmn_metdat.f
!c
!c     SET_AER0:  Iniitalize (climatology) aerosol OD and types (3 arrays)
!c		       called by MAIN, CHMSET
!c		COMMON BLOCKS: cmn_metdat.f
!c
!c     SET_ATM0:  Initialize climatologies for T & O3, set up atmospheric profiles
!c		COMMON BLOCKS: cmn_metdat.f
!c
!c<<<<<<<<<<<<<<<<<<<<<begin CTM-fastJX linking subroutines<<<<<<<<<<<<<<
!c  
!c     PHOTOJ(UTIME,IDAY,ILNG,JLAT, SOLF,SZA,U0,FREFL,ZPJ)
!c		Gateway to fast-JX, Update the photolysis rates
!c		COMMON BLOCKS: cmn_metdat.f, cmn_JVdat.f
!c
!c<<<<<<<<<<<<<<<<<<<<<begin core fast-J subroutines<<<<<<<<<<<<<<<<<<<<<
!c  N.B. all these need access to cmn_JVdat.f, but do NOT write into it.
!c	     also have no need to access cmn_metdat.f
!c
!c     OPTICD (OPTD,SSALB,SLEG, PATH,DENS,L)
!c	  UCI aerosol/CLOUD data sets, calculate scattering properties
!c		COMMON BLOCKS: cmn_JVdat.f
!c
!c     OPTICM (OPTD,SSALB,SLEG, PATH,RELH,L)
!c	   U Michigan aerosol data sets, generates fast-JX data formats
!c		COMMON BLOCKS: cmn_JVdat.f
!c
!c     JRATET(PPJ,TTJ,FFF, VALJL):  Calculate J-value, called by PTOTOJ.
!c		COMMON BLOCKS: cmn_JVdat.f
!c
!c     JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,DTAUX,POMEGAX,JXTRA)
!c		print out atmosphere used in J-value calc.
!c		>>>will be superseded by CTM routines
!c		COMMON BLOCKS: cmn_JVdat.f
!c
!c     RD_XXX(NJ1,NAMFIL):  Read wavelength bins, solar fluxes, Rayleigh
!c	       parameters, TEM-dependent X-sections, called once by INPHOT.
!c		COMMON BLOCKS: cmn_JVdat.f
!c		Input files: FJX_spec.dat
!c
!c     RD_MIE(NJ1,NAMFIL):  Set aerosols/cloud scattering, called once by INPHOT
!c		COMMON BLOCKS: cmn_JVdat.f
!c		Input files: FJX_scat.dat
!c
!c     RD_UM (NJ1,NAMFIL):  UMich aerosol optical data, called once by INPHOT
!c		COMMON BLOCKS: cmn_JVdat.f
!c		Input files: FJX_UMaer.dat
!c
!c     FUNCTION FLINT(TINT,T1,T2,T3,F1,F2,F3)
!c
!c     SOLARZ(GMTIME,NDAY,YGRDJ,XGRDI, SZA,COSSZA,SOLFX)
!c		calc SZA and Solar Flux factor for given lat/lon/UT
!c
!c     SPHERE2(U0,RAD,ZHL,ZZHT,AMF2,L1_):  
!c		calculate spherical geometry, air-mass factors (v 6.2)
!c
!c     EXTRAL(DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0, JXTRA)
!c		add sub-layers (JXTRA) to thick cloud/aerosol layers
!c
!c     OPMIE (KW, DTAUX,POMEGAX,U0,RFLECT,AMF,JXTRA,
!c    &       AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)
!c		calculate mean intensity (actinic) at each CTM levels
!c		calculate fluxes and deposition (heating rates)
!c		COMMON BLOCKS: cmn_JVdat.f
!c
!c<<<<<<<<<<<<<<<<<<<<<<<begin core scattering subroutines<<<<<<<<<<<<<<<
!c
!c	MIESCT (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
!c	      include 'parm_MIE.f' = dimension parameters
!c
!c	BLKSLV (FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,M,N,MFIT,ND)
!c		PARAMETER FILE: parm_MIE.f
!c
!c	GEN (POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,B,CC,AA,A,H,C1
!c	      ,M,N,MFIT,ND,ID)
!c		PARAMETER FILE: parm_MIE.f
!c
!c	LEGND0 (X,PL,N)
!c
!c	MATIN4 (A)
!c
!c	GAUSSP (N,XPT,XWT)
!c
!c-----------------------------------------------------------------------

  USE mo_kind,               ONLY: dp
  USE mo_physical_constants, ONLY: avo, grav, amd
  USE mo_debugs 
!! USE mo_geoloc,   ONLY: philat_2d,philon_2d
        
IMPLICIT NONE
!! INTEGER, PARAMETER :: nwv = 1     !!mgs!! QUICK FIX
! dummy arguments
      INTEGER, intent(in)     :: plonl
      INTEGER, intent(in)     :: klev,krow  
!++ost: amu0 -> u0 is not needed in the subroutine, delete?
      REAL(dp), intent(in)    :: u0(plonl)
!--ost
      REAL(dp), intent(in)    :: alb_s(plonl) 
      REAL(dp), intent(in)    :: h_pres(plonl,klev+1)
      REAL(dp), intent(in)    :: t_fld(plonl,klev)
      REAL(dp), intent(in)    :: h_temp(plonl,klev+1)
      REAL(dp), intent(in)    :: ho3(plonl,klev+1)
      REAL(dp), intent(in)    :: od(plonl,klev+1,nwv)
      REAL(dp), intent(in)    :: od600(plonl,klev+1)
      REAL(dp), intent(in)    :: ssalb(plonl,klev+1,nwv)
      REAL(dp), intent(in)    :: sleg(plonl,klev+1,nwv,8)
      REAL(dp), intent(in)    :: SOLF
      REAL(dp), intent(out) :: zpj(plonl,klev,jvn_)

! local variables
      CHARACTER*3 silng,sjlat
      INTEGER :: L,RATIO(w_)
      
      INTEGER  :: JL,JK,klevp2,klevp1
      INTEGER  :: KMIE, K, I, J
      REAL(dp) :: MASFAC, SCALEH, WAVE, &
                  TTT, XQO2, XQO3, ODABS, ODRAY!!, FLINT
      REAL(dp), dimension(plonl,klev,w_) :: fff, avgf, fjflx
      REAL(dp), dimension(plonl)         :: frefi, frefl, frefs, sza, &
                                             rflect
      REAL(dp), dimension(plonl,w_)      :: fjtop,fjbot,fsbot,flxd0, &
                                             flxup,dirup,flxdn,dirdn, &
					     fabot,fxbot,ffx0   
      REAL(dp), dimension(plonl,klev+2)  :: ppj, ttj, ddj, zzj, zhl  
      REAL(dp), dimension(plonl,klev+1,w_)   :: dtaux,flxd,ffx,flxj 
      REAL(dp), dimension(plonl,klev+1,w_,8) :: pomegax
      REAL(dp), dimension(plonl,2*klev+3,2*klev+3) :: amf2
      integer,  dimension(plonl,2*klev+3)          :: jxtra 
      REAL(dp), dimension(plonl,w_,8)    :: ffxnet          
      real(dp), dimension(plonl,klev,njval) ::  valjl
!C-----------------------------------------------------------------------      
      klevp2=klev+2
      klevp1=klev+1
!C-----------------------------------------------------------------------      
      ppj(:,:)=0._dp
      ttj(:,:)=0._dp
      zzj(:,:)=0._dp
      ddj(:,:)=0._dp
      zhl(:,:)=0._dp
      ZPJ(:,:,:) = 0._dp
      FFF(:,:,:) = 0._dp
!++mgs 20130304 : initialize FREFI - this seems to be a real FASTJ bug!
      FREFI(:)   = 0._dp
!--mgs
      FREFL(:)   = 0._dp
      FREFS(:)   = 0._dp
!C-----------------------------------------------------------------------
!     SOLARZ skipped as SZA is calculated by the model
!c---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
!c                        or         99.                  80 km
!c skip the calculation if all vector is in "dark"
      DO JK=1,plonl
       SZA(JK)=ACOS(U0(JK))/PI180 
      ENDDO
!!    zdf1(:,krow)=sza(:)
       
      DO JK=1,plonl
       if (SZA(JK) .le. SZAMAX) goto 198
      ENDDO
      write(*,*) 'all lat band skipped' 
      goto 199
 198  CONTINUE            

!c---load the amtospheric column data
      ppj(:,1:klevp1)=h_pres(:,1:klevp1)/100. !hPa
! last two levels
      ppj(:,klevp1)=0.01_dp !hPa
! one level above 
      ppj(:,klevp2)=0.005_dp  !hPA

      ttj(:,1:klevp1)=h_temp(:,1:klevp1) 
!c  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
     ! masfac = 100._dp*6.022e+23_dp/(28.97_dp*9.8_dp*10._dp)
      masfac = 100._dp*avo/(amd*grav*10._dp)

      do jl=1,klevp1
       do jk=1,plonl
        ddj(jk,jl)=(ppj(jk,jl)-ppj(jk,jl+1))*masfac
       enddo
      enddo 

      do jl=1,klevp1
       do jk=1,plonl
        !from VMR to molec/cm^2
        zzj(jk,jl) = ho3(jk,jl)*ddj(jk,jl)
       enddo
      enddo  		     
	
!c---calculate spherical weighting functions (AMF: Air Mass Factor)
      zhl(:,1)=1._dp !surface 
      do jl=1,klev
       do jk=1,plonl
        scaleh = 1.3806e-19_dp*masfac*h_temp(jk,jl)
        zhl(jk,jl+1) = zhl(jk,jl)-(log(ppj(jk,jl+1)/ppj(jk,jl))*scaleh)
       enddo 
      enddo
      zhl(:,klevp2) = zhl(:,klevp1) + zzht !zzht is a constant

!C-----------------------------------------------------------------------
!++mgs: deactivated following routine because it crashes !!baustelle!!
   call SPHERE2_V(plonl,U0,RAD,ZHL,ZZHT,AMF2,KLEVP1)
!!!      AMF2(:,:,:) = 1._dp
!--mgs
!C-----------------------------------------------------------------------
!c---Now given the aerosol+cloud OD/layer in visible (600 nm) can calculate
!c        how to add additonal levels at top of clouds (now uses log spacing)
!c---Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add 
!C       additonal levels at top of clouds (now uses log spacing)
      call EXTRAL_V(plonl,OD600,KLEVP1,2*KLEV+2,N_,JTAUMX,ATAU,ATAU0,JXTRA)
!C-----------------------------------------------------------------------

!c---set surface reflectance
      RFLECT(1:plonl) = ALB_S(1:plonl)
      DO JK=1,plonl
        RFLECT(JK) = max(0._dp,min(1._dp,RFLECT(JK)))
      ENDDO

!C---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
      do K = NW1,NW2
        WAVE = WL(K)
!C---Pick nearest Mie wavelength to get scattering properites------------
                                KMIE=1  ! use 200 nm prop for <255 nm
        if( WAVE .gt. 255._dp ) KMIE=2  ! use 300 nm prop for 255-355 nm
        if( WAVE .gt. 355._dp ) KMIE=3  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500._dp ) KMIE=4
        if( WAVE .gt. 800._dp ) KMIE=5 

!c---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
!c---values at KLEVP1=KLEV+1 are a pseudo/climatol layer above the top CTM layer (L_)
	do JL = 1,KLEVP1
	 DO JK=1,plonl
          TTT  = TTJ(JK,JL)
	  XQO3 = FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2) &
        		     ,QO3(K,1),QO3(K,2),QO3(K,3))

	  XQO2 = FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1) &
        		     ,QO2(K,1),QO2(K,2),QO2(K,3))

          ODABS = XQO3*ZZJ(JK,JL) + XQO2*DDJ(JK,JL)*0.20948_dp
          ODRAY = DDJ(JK,JL)*QRAYL(K)

          DTAUX(JK,JL,K) = OD(JK,JL,KMIE) + ODABS + ODRAY

          do I=1,8
           POMEGAX(JK,JL,K,I) = SLEG(JK,JL,KMIE,I)*OD(JK,JL,KMIE)
          enddo
           POMEGAX(JK,JL,K,1) = POMEGAX(JK,JL,K,1) + 1.0_dp*ODRAY
           POMEGAX(JK,JL,K,3) = POMEGAX(JK,JL,K,3) + 0.5_dp*ODRAY
          do I=1,8
           POMEGAX(JK,JL,K,I) = POMEGAX(JK,JL,K,I)/DTAUX(JK,JL,K)
          enddo
	 ENDDO
	enddo
       ENDDO !K

!C-----------------------------------------------------------------------
!#debug#
!!write(0,*) '#debug fastJ: shape(DTAUX) = ',shape(DTAUX)
!!write(0,*) '#debug fastJ: shape(POMEGAX) = ',shape(POMEGAX)
!!write(0,*) '#debug fastJ: shape(U0) = ',shape(U0)
!!write(0,*) '#debug fastJ: shape(RFLECT) = ',shape(RFLECT)
!!write(0,*) '#debug fastJ: shape(AMF2) = ',shape(AMF2)
!!write(0,*) '#debug fastJ: shape(JXTRA) = ',shape(JXTRA)
!!write(0,*) '#debug fastJ: shape(AVGF) = ',shape(AVGF)
!!write(0,*) '#debug fastJ: shape(FJTOP) = ',shape(FJTOP)
!!write(0,*) '#debug fastJ: shape(FJBOT) = ',shape(FJBOT)
!!write(0,*) '#debug fastJ: shape(FSBOT) = ',shape(FSBOT)
!!write(0,*) '#debug fastJ: shape(FJFLX) = ',shape(FJFLX)
!!write(0,*) '#debug fastJ: shape(FLXD) = ',shape(FLXD)
!!write(0,*) '#debug fastJ: shape(FLXD0) = ',shape(FLXD0)

    call OPMIE_V2 (plonl,KLEV,                        &
                   DTAUX,POMEGAX,U0,RFLECT,AMF2,JXTRA, &
                   AVGF,FJTOP,FJBOT,                   &
                   FSBOT,FJFLX,FLXD,FLXD0)
!C-----------------------------------------------------------------------
!c----direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
!c----     also at bottom (DN), does not include diffuse reflected flux.

!C SEEMS TO BE SOME DIAGNOSTIC OUTPUT DONT CARE FOR THE MOMENT OR ANYWAY CALC IN OPMIE
	DO K=NW1,NW2
	 FLXUP(:,K) =  FJTOP(:,K)
	 DIRUP(:,K) = -FLXD0(:,K)
	 FLXDN(:,K) = -FJBOT(:,K)
	 DIRDN(:,K) = -FSBOT(:,K)
	ENDDO 
!C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	DO K=NW1,NW2
	 do JL = 1,KLEV
	  DO JK=1,plonl
	   FFF(JK,JL,K) = FFF(JK,JL,K) + SOLF*FL(K)*AVGF(JK,JL,K) !SOLF IS CONSTANT
	  ENDDO   
	 enddo
	ENDDO 

!++mgs: !!baustelle!! deactivated the following block because it led to code crash	
      DO K=NW1,NW2
	 DO JK=1,plonl 
	  FREFI(JK) = FREFI(JK) + SOLF*FL(K)*FLXD0(JK,K)/WL(K) 
	  FREFL(JK) = FREFL(JK) + SOLF*FL(K)*FJTOP(JK,K)/WL(K) 
	  FREFS(JK) = FREFS(JK) + SOLF*FL(K)/WL(K)
	 ENDDO
	ENDDO
!--mgs

!c---for each wavelength calculate the flux budget/heating rates:
!c  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L)]
!c	       but for spherical atmosphere!
!c  FJFLX(L) = diffuse flux across top of layer L

!c---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
!c---     need special fix at top and bottom: 
!c---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
	DO K=NW1,NW2
	 DO JK=1,plonl 
	  FABOT(JK,K) = (1.d0-RFLECT(JK))*(FJBOT(JK,K)+FSBOT(JK,K))
	  FXBOT(JK,K) = -FJBOT(JK,K) + RFLECT(JK)*(FJBOT(JK,K)+FSBOT(JK,K))
	  FLXJ(JK,1,K) = FJFLX(JK,1,K) - FXBOT(JK,K)
	 ENDDO
	ENDDO  
	
	DO K=NW1,NW2
	 do JL=2,KLEV
	  DO JK=1,plonl 
	   FLXJ(JK,JL,K) = FJFLX(JK,JL,K) - FJFLX(JK,JL-1,K)
	  enddo
	 ENDDO
	ENDDO  
	
	DO K=NW1,NW2
	 DO JK=1,plonl 
	  FLXJ(JK,KLEV+1,K) = FJTOP(JK,K) - FJFLX(JK,KLEV,K)
	 ENDDO
	ENDDO
	  
!c---calculate net flux deposited in each CTM layer (direct & diffuse):
	FFX0(:,:) = 0._dp
	DO K=NW1,NW2
	 do JL=1,KLEVP1
	  DO JK=1,plonl 
	   FFX(JK,JL,K) = FLXD(JK,JL,K) - FLXJ(JK,JL,K)
	   FFX0(JK,K) = FFX0(JK,K) + FFX(JK,JL,K)
	  ENDDO
	enddo
       ENDDO 

!c  NB: the radiation level ABOVE the top CTM level is included in these budgets
!c      these are the flux budget/heating terms for the column:
!c  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
!c  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface) 
!c  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
!c  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos  
!c  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos 
!c  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
!c       these are surface fluxes to compare direct vs. diffuse:
!c  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
!c  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)

      DO K=NW1,NW2
       FFXNET(:,K,1) = FLXD0(:,K)
       FFXNET(:,K,2) = FSBOT(:,K)
       FFXNET(:,K,3) = FLXD0(:,K)+FSBOT(:,K)
       FFXNET(:,K,4) = FJTOP(:,K)
       FFXNET(:,K,5) = FFX0(:,K)
       FFXNET(:,K,6) = FABOT(:,K)
       FFXNET(:,K,7) = FSBOT(:,K)
       FFXNET(:,K,8) = FJBOT(:,K)
      ENDDO		 

!++mgs: !!baustelle!! changed the following statements to avoid crashing
      FREFL(:) = FREFL(:)/FREFS(:)	!calculate reflected flux (energy weighted)
      FREFI(:) = FREFI(:)/FREFS(:)
!!      WHERE( FREFS(:) > 1.e-4_dp )
!!         FREFL(:) = FREFL(:)/FREFS(:)	!calculate reflected flux (energy weighted)
!!         FREFI(:) = FREFI(:)/FREFS(:)
!!      END WHERE
!--mgs

!c---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)
!C-----------------------------------------------------------------------
      call JRATET_V(plonl,KLEV,PPJ,TTJ,FFF,VALJL)
!C-----------------------------------------------------------------------
!c---map the J-values from fast-JX onto ASAD ones (use JIND & JFACTA)
	do J = 1,NRATJ
	 if (JIND(J).gt.0) then
	   ZPJ(:,:,J) = VALJL(:,:,JIND(J))*JFACTA(J)
	 else
	   ZPJ(:,:,J) = 0._dp
	 endif
	enddo

!---------------------------------------------------------------------------
!C--WRITE OUTPUT DATA
!!!      DO I=1,plonl

!!!      write(silng,'(i3.3)') I 
!!!      write(sjlat,'(i3.3)') KROW 
!!!      OPEN(2,file='FJX_ref_csk_LNG'//TRIM(silng)//'_LAT'//TRIM(sjlat)//'.dat')

!!!      WRITE(2,*)'fast-JX-(6.1)----PHOTOJ internal print: Atmosphere---'
!WRITE DTAUX AND POMEGAX AT LAST WAVELENGHT WHICH SHOULD BE AR 600 NM
!THE SAME WAS DOES DONE IN NOT VECTORIZED VERSION

!!!      call JP_ATM(KLEV,PPJ(I,:),TTJ(I,:),DDJ(I,:),ZZJ(I,:),ZHL(I,:), &
!!!                 ZZHT,DTAUX(I,:,NW2),POMEGAX(I,:,NW2,:),JXTRA(I,:))

!C---PRINT SUMMARY of mean intensity, flux, heating rates:
!!!      WRITE(2,*)
!!!      WRITE(2,*)'fast-JX(6.1)----PHOTOJ internal print: Mean Intens---'
!!!      WRITE(2,'(a,5f10.4)') ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/', &
!!!                            RFLECT(I),SZA(I),U0(I),FREFI(I),FREFL(I)

!!!      WRITE(2,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
!!!      WRITE(2,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
!!!      WRITE(2,'(a)') ' ---  100000=Fsolar   MEAN INTENSITY per wvl bin'
!!!      do L = KLEV,1,-1
!!!       do K=NW1,NW2
!!!        RATIO(K) = (1.d5*FFF(I,L,K)/FL(K))
!!!       enddo
!!!        WRITE(2,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
!!!      enddo
      
!!!      WRITE(2,*)
!!!      WRITE(2,*)'fast-JX(6.1)----PHOTOJ internal print: Net Fluxes----'
!!!      WRITE(2,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
!!!      WRITE(2,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
!c      WRITE(2,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), K=NW2,NW1,-1)
!c      WRITE(2,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), K=NW2,NW1,-1)
!!!      WRITE(2,*) ' ---NET FLUXES--- '
!!!      WRITE(2,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(I,K,3),K=NW2,NW1,-1)
!!!      WRITE(2,'(a11,18f8.4)') ' dif outtop',(FFXNET(I,K,4),K=NW2,NW1,-1)
!!!      WRITE(2,'(a11,18f8.4)') ' abs in atm',(FFXNET(I,K,5),K=NW2,NW1,-1)
!!!      WRITE(2,'(a11,18f8.4)') ' abs at srf',(FFXNET(I,K,6),K=NW2,NW1,-1)
!!!      WRITE(2,*) ' ---SRF FLUXES--- '
!!!      WRITE(2,'(a11,18f8.4)') ' srf direct',(FFXNET(I,K,7),K=NW2,NW1,-1)
!!!      WRITE(2,'(a11,18f8.4)') ' srf diffus',(FFXNET(I,K,8),K=NW2,NW1,-1)
!!!      WRITE(2,'(2a)') '  ---NET ABS per layer:       10000=Fsolar', &
!!!                      '  [NB: values <0 = numerical error w/clouds or SZA>90, colm OK]'
!!!      do L = KLEV,1,-1
!!!      do K=NW1,NW2
!!!        RATIO(K) = 1.d5*FFX(I,L,K)
!!!       enddo
!!!        WRITE(2,'(i9,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
!!!      enddo

!c---PRINT j VALUES
!!!      WRITE(2,'(a)')
!!!      WRITE(2,'(a)') ' fast-JX (6.1)----J-values----'
!!!	write(2,'(a,f7.2,a,2f7.2,a,f8.5)') &
!!!           '  SZA=',SZA(I), &
!!!           '  LAT x LONG=',philat_2d(I,KROW),philon_2d(I,KROW),' SolarFlx=',SOLF
!!!	write(2,'(1x,a,64(a7,2x))') 'L=  ',(JLABEL(K), K=1,JVN_)
!!!       do L=KLEV,1,-1
!!!	write(2,'(i3,1p, 64e9.2)') L,(ZPJ(I,L,K),K=1,JVN_)
!!!       enddo
       
!!!      CLOSE(2) 
!!!      ENDDO
!------------------------------------------------------------------------------

 199  CONTINUE
      RETURN
         
END SUBROUTINE FJX       

!c-----------------------------------------------------------------------
      subroutine SPHERE2_V(plonl,GMU,RAD,ZHL,ZZHT,AMF2,L1_)
!C-----------------------------------------------------------------------
!c----new v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
!c  This new AMF2 does each of the half-layers of the CTM separately,
!c     whereas the original, based on the pratmo code did the whole layers
!c     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
!c  Since fast-JX is meant to calculate the intensity at the mid-layer, the
!c     solar beam at low sun (interpolated between layer edges) was incorrect.
!c  This new model does make some approximations of the geometry of the layers:
!c     the CTM layer is split evenly in mass (good) and in height (approx).
!c
!c  Calculation of spherical geometry; derive tangent heights, slant path
!c  lengths and air mass factor for each layer. Not called when
!c  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!c  beam (where tangent height is below altitude J-value desired at).
!C-----------------------------------------------------------------------
!c in:
!c     GMU     = MU0 = cos(solar zenith angle)
!c     RAD     radius of Earth mean sea level (cm)
!c     ZHL(L)  height (cm) of the bottome edge of CTM level L
!c     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1)
!c     L1_     dimension of CTM = levels +1 (L+1 = above-CTM level)
!c out:
!c     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
!C-----------------------------------------------------------------------
      USE mo_kind,   ONLY: dp

      implicit none

      INTEGER,INTENT(IN)  ::   plonl
      integer, intent(in) ::   L1_
      real(dp), intent(in)  ::   GMU(plonl),RAD,ZHL(plonl,L1_+1),ZZHT
      real(dp), intent(out) ::   AMF2(plonl,2*L1_+1,2*L1_+1)

!c     RZ      Distance from centre of Earth to each point (cm)
!c     RQ      Square of radius ratios
!c     SHADHT  Shadow height for the current SZA
!c     XL      Slant path between points

      integer  I, J, K, II, L2
      real(dp)   XMU1(plonl),XMU2(plonl,2*L1_+1),XL(plonl,2*L1_+1)
      REAL(dp)   DIFF(plonl,2*L1_+1),SHADHT(plonl),RZ(plonl,L1_+1)
      real(dp)   RZ2(plonl,2*L1_+1),RQ2(plonl,2*L1_+1)
      INTEGER  JK
!c
!c--- must have top-of-atmos (NOT top-of-CTM) defined
!c      ZHL(L1_+1) = ZHL(L1_) + ZZHT

      RZ(1:plonl,1) = RAD + ZHL(1:plonl,1)
      do II = 2,L1_+1
        RZ(1:plonl,II)   = RAD + ZHL(1:plonl,II)
      enddo

!c---calculate heights for edges of split CTM-layers
      L2 = 2*L1_
      do II = 2,L2,2
        I = II/2
        RZ2(1:plonl,II-1) = RZ(1:plonl,I)
        RZ2(1:plonl,II) = 0.5_dp*(RZ(1:plonl,I)+RZ(1:plonl,I+1))
      enddo
      RZ2(1:plonl,L2+1) = RZ(1:plonl,L1_+1)
      do II = 1,L2
        RQ2(1:plonl,II) = (RZ2(1:plonl,II)/RZ2(1:plonl,II+1))**2 
      enddo

!c---shadow height for SZA > 90
      DO JK=1,plonl
       if (GMU(JK) .lt. 0.0_dp)  then
        SHADHT(JK) = RZ2(JK,1)/dsqrt(1.0_dp-GMU(JK)**2)
       else
        SHADHT(JK) = 0._dp
       endif
      ENDDO  
!c---up from the surface calculating the slant paths between each level
!c---  and the level above, and deriving the appropriate Air Mass Factor
      AMF2(:,:,:) = 0._dp

      XMU1(:) = abs(GMU(:))
      do J = 1,2*L1_+1
        do I = J,2*L1_
         DO JK=1,plonl
!c  Air Mass Factors all zero if below the tangent height
          if (RZ2(JK,J) .lt. SHADHT(JK)) goto 16
!c  Ascend from layer J calculating AMF2s
!++mgs 20130304: the following line caused an arithmetic exception (sqrt (<0)) - now defused
!!ori!!    XMU2(JK,I) = dsqrt(1.0_dp - RQ2(JK,I)*(1.0_dp-XMU1(JK)**2))
           if (RQ2(JK,I)*(1.0_dp-XMU1(JK)**2) > 1._dp) then
              if (error_count < 50) then
                 write(0,*) "ERROR in fast-j|sphere_v2: negative root: jk, i, XMU1(JK), ", &
                   "RQ2(JK,I), (1.0_dp - RQ2(JK,I)*(1.0_dp-XMU1(JK)**2)) = ", &
                   jk, i, XMU1(JK), RQ2(JK,I), (1.0_dp - RQ2(JK,I)*(1.0_dp-XMU1(JK)**2))
              end if
              error_count = error_count + 1
              XMU2(JK,I) = 0._dp
           else
              XMU2(JK,I) = dsqrt(1.0_dp - RQ2(JK,I)*(1.0_dp-XMU1(JK)**2))
           end if
!--mgs
           XL(JK,I)   = RZ2(JK,I+1)*XMU2(JK,I) - RZ2(JK,I)*XMU1(JK)
!++mgs !!baustelle!! code crashes in following line !
           AMF2(JK,I,J) = XL(JK,I) / (RZ2(JK,I+1)-RZ2(JK,I))
!#debug
!           IF ( ABS((RZ2(JK,I+1)-RZ2(JK,I))) < 1.e-4_dp ) THEN
!              write(0,*) '#debug# fast-j, sphere_v2: jk, i, RZ2(JK,I), RZ2(JK,I+1) = ', &
!                              jk, i, RZ2(JK,I), RZ2(JK,I+1)
!              AMF2(JK,I,J) = 1.e4_dp * XL(JK,I)
!           ELSE
!              AMF2(JK,I,J) = XL(JK,I) / (RZ2(JK,I+1)-RZ2(JK,I))
!           END IF
!--mgs   (end #debug#) ---------------
           XMU1(JK)     = XMU2(JK,I)
  16     ENDDO
	enddo
!c--fix above top-of-atmos (L=L1_+1), must set DTAU(L1_+1)=0
        AMF2(:,2*L1_+1,J) = 1._dp

!c  Twilight case - Emergent Beam, calc air mass factors below layer
        XMU1(:)       = abs(GMU(:))
!c  Descend from layer J 
         do II = J-1,1,-1
          DO JK=1,plonl
           if (GMU(JK) .ge. 0.0_dp) goto 18
            DIFF(JK,II) = RZ2(JK,II+1)*sqrt(1.0_dp-XMU1(JK)**2) &
      	                  -RZ2(JK,II)
            if (II.eq.1)  DIFF(JK,II) = max(DIFF(JK,II),0._dp)   ! filter
!c  Tangent height below current level - beam passes through twice
            if (DIFF(JK,II) .lt. 0.0_dp)  then
             XMU2(JK,II)= sqrt(1.0_dp - (1.0_dp-XMU1(JK)**2)/RQ2(JK,II))
             XL(JK,II)  = abs(RZ2(JK,II+1)*XMU1(JK)-RZ2(JK,II) &
     	                  *XMU2(JK,II))
             AMF2(JK,II,J) = 2._dp*XL(JK,II)/(RZ2(JK,II+1)-RZ2(JK,II))
             XMU1(JK)      = XMU2(JK,II)
!c  Lowest level intersected by emergent beam
            else
             XL(JK,II)     = RZ2(JK,II+1)*XMU1(JK)*2.0_dp
             AMF2(JK,II,J) = XL(JK,II)/(RZ2(JK,II+1)-RZ2(JK,II))
            goto 18
          endif
   18	 ENDDO 
        enddo
      
      ENDDO !J
        
      return
      end subroutine sphere2_v
!c-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
      subroutine EXTRAL_V(plonl,DTAUX,L1X,L2X,NX, &
                          JTAUMX,ATAU,ATAU0,JXTRA)
!C-----------------------------------------------------------------------
!c
!c    new version 6.1, add sub-layers (JXTRA) to thick cloud/aerosol layers
!c    this version sets up log-spaced sub-layers of increasing thickness ATAU
!c
!c     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
!c        This can be just cloud or cloud+aerosol, it is used only to set
!c        the number in levels to insert in each layer L
!c        Set for log-spacing of tau levels, increasing top-down.
!c
!c     N.B. the TTAU, etc calculated here are NOT used elsewhere

!c---The log-spacing parameters have been tested for convergence and chosen
!c---  to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
!c---  use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100 
!c---  ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
!C-----------------------------------------------------------------------
 
      USE mo_kind, ONLY: dp 

      implicit none
      
      INTEGER, INTENT(IN) ::  plonl
      integer, intent(in) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real(dp),intent(in) ::  DTAUX(plonl,L1X)      !cloud+3aerosol OD in each layer
      real(dp),intent(in) ::  ATAU,ATAU0
      integer, intent(out)::  JXTRA(plonl,L2X+1)    !number of sub-layers to be added

      integer JTOTL(plonl),I,L,L2 
      real(dp)  TTAU(plonl,L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1
      INTEGER JK

!C---Reinitialize arrays
      TTAU(:,:)  = 0._dp
      JXTRA(:,:) = 0

!c---combine these edge- and mid-layer points into grid of size:
!c---              L2X+1 = 2*L1X+1 = 2*L_+3
!c---calculate column optical depths above each level, TTAU(1:L2X+1)
!c---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
!c
!c---Divide thick layers to achieve better accuracy in the scattering code
!c---In the original fast-J, equal sub-layers were chosen, this is wasteful
!c---and this new code (ver 5.3) uses log-scale:  
!c---        Each succesive layer (down) increase thickness by ATAU > 1
!c---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
!c---        4 sub-layers with ODs = 1 - 2 - 4 - 8
!c---The key parameters are:
!c---        ATAU = factor increase from one layer to the next
!c---        ATAUMN = the smallest OD layer desired
!c---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
!c---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1._dp
      ATAULN = log(ATAU)
      TTAU(:,L2X+1)  = 0.0_dp
      do L2 = L2X,1,-1
       L         = (L2+1)/2
       DO JK=1,plonl
        DTAUJ     = 0.5_dp * DTAUX(JK,L)
        TTAU(JK,L2)  = TTAU(JK,L2+1) + DTAUJ
!c---Now compute the number of log-spaced sub-layers to be added in
!c---   the interval TTAU(L2) > TTAU(L2+1)
!c---The objective is to have successive TAU-layers increasing by factor ATAU >1
!c---the number of sub-layers + 1
        if (TTAU(JK,L2) .lt. ATAU0) then
          JXTRA(JK,L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(JK,L2+1))
          ATAUN1 = log(TTAU(JK,L2)/ATAUM) / ATAULN
          JXTRA(JK,L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5_dp)))
        endif
       ENDDO 	
      enddo

!c---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL(:)    = L2X + 2
      do L2 = L2X,1,-1
       DO JK=1,plonl
        JTOTL(JK)  = JTOTL(JK) + JXTRA(JK,L2)
        if (JTOTL(JK) .gt. NX/2)  then
          write(*,'(A,3I5)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
          do L = L2,1,-1
            JXTRA(JK,L) = 0
          enddo
        endif
       ENDDO
      enddo

      return
      end subroutine extral_v
!c-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
      subroutine OPMIE_V2(plonl,KLEV,                                  &
                          DTAUX,POMEGAX,U0,RFLECT,AMF2,                 & 
                          JXTRA,FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0)
!C-----------------------------------------------------------------------
!C  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
!c in:    
!c     KW = wavelength bin # (1:18)  - used only for some diagnostics !!
!c     DTAUX(1:L1_) = optical depth of each lauer
!c     POMEGAX(1:8,1:L1_) = scattering phase fn (multiplied by s-s abledo)
!c     U0  = cos (SZA)
!c     RFLECT = Lambertian albedo of surface
!c     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
!c        AMF2 now does both edges and middle of CTM layers
!c     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
!c out:
!c     FJACT(1:L_) = mean actinic flux(diff+direct) at std CTM levels(mid-layer)
!c  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
!c     FJTOP = diffuse flux out top-of-atmosphere (TAU=0, above top model lauer)
!c     FJBOT = diffuse flux onto surface (<0 by definition)
!c     FSBOT = direct/solar flux onto surface  (<0 by definition)
!c     FJFLX(1:L_) = diffuse flux across top of model layer L
!C        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
!c     FLXD(1:L_+1) = solar flux deposited in layer L (includes layer above CTM)
!c        this should take into account sphericity, and is not just = mu0
!c     FLXD0 = sum of solar flux deposited in atmos
!c        does NOT include flux on lower surface, does NOT mean absorbed!
!C-----------------------------------------------------------------------
!c
!c     DTAU     Local optical depth of each CTM level
!c     TTAU     Optical depth of air vertically above each point (to top of atm)
!c     FTAU2     Attenuation of solar beam
!c     POMEGAJ  Scattering phase function
!c
!c---new Ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the 
!c   factor increase from sub-layer to sub-layer 
!C-----------------------------------------------------------------------
      USE mo_kind,       ONLY: dp
      USE mo_exception,  ONLY: finish
      
      implicit none

      integer,  intent(in) ::  plonl,KLEV
      real(dp), intent(in) ::  DTAUX(plonl,KLEV+1,W_)
      real(dp), intent(in) ::  POMEGAX(plonl,KLEV+1,W_,8)
      real(dp), intent(in) ::  AMF2(plonl,2*KLEV+3,2*KLEV+3)
      real(dp), intent(in) ::  U0(plonl),RFLECT(plonl)
      integer,  intent(in) ::  JXTRA(plonl,2*KLEV+3)
      real(dp), intent(out)::  FJACT(plonl,KLEV,W_),FJTOP(plonl,W_)
      real(dp), intent(out)::  FJBOT(plonl,W_),FSBOT(plonl,W_)
      real(dp), intent(out)::  FJFLX(plonl,KLEV,W_),FLXD(plonl,KLEV+1,W_)
      real(dp), intent(out)::  FLXD0(plonl,W_)

      integer JTOTL,I,II,J,K,L,LL,IX,JK,L2,L2L,L22,LZ,LZZ
      integer LZ0,LZMID
      real(dp)  SUMT,SUMJ

      real(dp)  DTAU(plonl,KLEV+2,W_)
      real(dp)  POMEGAB(2*M_),FBTM,FTOP
      real(dp)  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,TAUBTM,TAUTOP
!c--- variables used in mie code-----------------------------------------
      real(dp)   FJFLX0
      integer  MFIT
!c--- vectorized --------------------------------------------------------
      REAL(dp), DIMENSION(plonl,W_)           :: ZFLUX,FJT,FJB 
      REAL(dp), DIMENSION(plonl,N_,W_)        :: FZ,ZTAU,FJ
      REAL(dp), DIMENSION(plonl,2*M_,N_,W_)   :: POMEGA
!!      REAL(dp), ALLOCATABLE  :: POMEGA(:,:,:,:),ZTAU(:,:,:),FZ(:,:,:),FJ(:,:,:)  

      REAL(dp), DIMENSION(plonl,2*M_,2*KLEV+3,W_):: POMEGAJ
      REAL(dp), DIMENSION(plonl,2*KLEV+3,W_)     :: TTAU,FTAU2
      REAL(dp), DIMENSION(plonl,KLEV,W_)         :: FJACT2
      REAL(dp), DIMENSION(plonl,2*KLEV+2,W_)     :: FLXD2

      INTEGER,DIMENSION(plonl,2*KLEV+3)     :: JADDLV,JADDTO,L2LEV
      INTEGER,DIMENSION(plonl,KLEV)         :: JNDLEV
      INTEGER,DIMENSION(plonl,KLEV+1)       :: JNELEV      
      INTEGER,DIMENSION(plonl)              :: NDZ,LZ1

      INTEGER :: L1_, L2_, L_  

!C---Set level indices
     L_ =KLEV
     L1_=KLEV+1
     L2_=2*L1_
            
!C---Reinitialize arrays
      ZTAU(:,:,:)     = 0._dp
      FZ(:,:,:)       = 0._dp
      POMEGA(:,:,:,:) = 0._dp

!C---Set up optical depth DTAU(L)
      do L = 1,L1_
       DTAU(1:plonl,L,NW1:NW2) = DTAUX(1:plonl,L,NW1:NW2)
      enddo
      DTAU(1:plonl,L1_+1,NW1:NW2) = 0._dp
!c---Define the total scattering phase fn for each CTM layer L=1:L_+1
!c---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
!C---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
      MFIT = 2*M_
      do L = 1,L1_
        do I = 1,MFIT
         POMEGAJ(1:plonl,I,L,NW1:NW2) = POMEGAX(1:plonl,L,NW1:NW2,I)
        enddo
      enddo

!c---version 6.2 fix to do mid-layers correctly at large SZA (mp, 6/2008)

!C---Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
!c---      at the middle & edges of the CTM layers L=1:2*L1_+1
!c---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
!c---  note that DTAU(L1_) is optical depth in the FULL CTM layer just above
      

      FTAU2(:,:,:) = 0._dp
      FTAU2(:,L2_+1,:) = 1.0_dp

      DO K=NW1,NW2
       do LL = 1,2*L1_+1
        L = (LL+1)/2
        DO JK =1,plonl 
         if (AMF2(JK,LL,LL) .gt. 0.0_dp) then
          XLTAU = 0.0_dp
          do II = 1,2*L1_+1
           I = (II+1)/2
           XLTAU = XLTAU + 0.5_dp*DTAU(JK,I,K)*AMF2(JK,II,LL)
          enddo
          if (XLTAU .lt. 76._dp) then   ! zero out flux at 1e-33
           FTAU2(JK,LL,K) = exp(-XLTAU)
          endif
         endif
        ENDDO 	
       enddo
      ENDDO 

!c---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
!c---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
      FLXD2(1:plonl,:,NW1:NW2) = 0._dp
      DO K=NW1,NW2
       do LL = 1,2*L1_
        DO JK =1,plonl 
         if (AMF2(JK,LL,LL) .gt. 0._dp) then 
          FLXD2(JK,LL,K) = (FTAU2(JK,LL+1,K) -  &
      	                    FTAU2(JK,LL,K))/AMF2(JK,LL,LL)
         endif
        ENDDO 
       enddo
      ENDDO
      
      DO K=NW1,NW2
       DO JK =1,plonl 
        if (AMF2(JK,1,1) .gt. 0._dp) then 
         FSBOT(JK,K) = FTAU2(JK,1,K)/AMF2(JK,1,1)
        else
         FSBOT(JK,K) = 0._dp
        endif
       ENDDO
      ENDDO
       
      DO K=NW1,NW2
       do LL = 2,2*L1_,2
       L=LL/2
        DO JK =1,plonl 
         FLXD(JK,L,K) = FLXD2(JK,LL,K)+FLXD2(JK,LL-1,K)
        ENDDO   
       enddo
      ENDDO
      
!c---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
!c---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
      FLXD0(1:plonl,NW1:NW2) = 0._dp
      DO K=NW1,NW2
       DO JK =1,plonl 
        if (AMF2(JK,2*L1_,2*L1_) .gt. 0._dp) then
         do L=1,L1_
          FLXD0(JK,K) = FLXD0(JK,K) + FLXD(JK,L,K)
         enddo
        endif
       ENDDO
      ENDDO 
!C------------------------------------------------------------------------
!c  Take optical properties on CTM layers and convert to a photolysis
!c  level grid corresponding to layer centres and boundaries. This is
!c  required so that J-values can be calculated for the centre of CTM
!c  layers; the index of these layers is kept in the JNDLEV array.
!C------------------------------------------------------------------------
!c
!c---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
!c---    points (1:L_) plus 1 for the mid point of added top layer.

!c---combine these edge- and mid-layer points into grid of size:
!c---              L2_+1 = 2*L1_+1 = 2*L_+3
!c---calculate column optical depths above each level, TTAU(1:L2_+1)
!c---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD
      TTAU(1:plonl,L2_+1,NW1:NW2) = 0.0_dp
      DO K=NW1,NW2
       do L2 = L2_,1,-1
        DO JK=1,plonl
         L          = (L2+1)/2
         DTAUJ      = 0.5_dp * DTAU(JK,L,K)
         TTAU(JK,L2,K)   = TTAU(JK,L2+1,K) + DTAUJ
        ENDDO
       enddo
      ENDDO 
!c
!c----solar flux incident on lower boundary & Lambertian reflect factor:
      DO K=NW1,NW2
       DO JK=1,plonl
        if (FSBOT(JK,K) .gt. 0._dp) then
!c---        FSBOT = U0*FTAU(1)
!c---        ZFLUX = U0*FTAU(1)*RFLECT/(1.d0+RFLECT)
         ZFLUX(JK,K) = FSBOT(JK,K)*RFLECT(JK)/(1._dp+RFLECT(JK))
        else
         ZFLUX(JK,K) = 0._dp
        endif
       ENDDO
      ENDDO  

!c  Calculate scattering properties, level centres then level boundaries
!c ***be careful of order, we are overwriting/shifting the 'POMEGAJ' upward in index***
      DO K=NW1,NW2
       do L2 = L2_,2,-2
        DO JK=1,plonl
         L   = L2/2
         do I = 1,MFIT
          POMEGAJ(JK,I,L2,K) = POMEGAJ(JK,I,L,K)
         enddo
        ENDDO 
       enddo
      ENDDO
!c---lower boundary value is set (POMEGAJ(I,1), but set upper:
      DO K=NW1,NW2
       do I = 1,MFIT
        DO JK=1,plonl
          POMEGAJ(JK,I,L2_+1,K) = POMEGAJ(JK,I,L2_,K)
        ENDDO
       enddo
      ENDDO 
!c---now have POMEGAJ filled at even points from L2=3:L2_-1
!c---use inverse interpolation for correct tau-weighted values at edges
      DO K=NW1,NW2
       do L2 = 3,L2_-1,2
        DO JK=1,plonl
         TAUDN = TTAU(JK,L2-1,K)-TTAU(JK,L2,K)
         TAUUP = TTAU(JK,L2,K)-TTAU(JK,L2+1,K)
         do I = 1,MFIT
          POMEGAJ(JK,I,L2,K) = (POMEGAJ(JK,I,L2-1,K)*TAUDN + & 
                   POMEGAJ(JK,I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
         enddo
        ENDDO
       enddo
      ENDDO 
!C---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
!c---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface
!c
!C------------------------------------------------------------------------
!c  Calculate cumulative total and define levels we want J-values at.
!c  Sum upwards for levels, and then downwards for Mie code readjustments.
!c
!c     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
!c           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
!c     JADDLV(L2)  Number of new levels actually added at each wavelength
!c            where JADDLV = 0 when there is effectively no FTAU2 
!c     JADDTO(L2)   Total number of new levels to add to and above level (L2)
!c     JNDLEV(L) = L2 index that maps on CTM mid-layer L
!c
!C------------------------------------------------------------------------

!c---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
!c---    JADDLV is taken from JXTRA, which is based on visible OD.
!c---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
!c---note that JADDLV and JADDTO will change with wavelength and solar zenith

!c--now try to insert additional levels for thick clouds, ONLY IF FTAU2 > 1.e-8
!c-- this will cut off additional levels where the solar beam is negligible.

!c---new v5.6-----keep all wavelengths the same for now
!c      do L2 = 1,L2_,1
!c        if (FTAU2(L2+1) .lt. 1.d-30) then
!c          JADDLV(L2) = 0
!c        else
!c          JADDLV(L2) = JXTRA(L2)
!c        endif
!c      enddo
      do L2 = 1,L2_,1
       DO JK=1,plonl      
        JADDLV(JK,L2) = JXTRA(JK,L2)
       ENDDO
      enddo
      JADDTO(:,L2_+1) = 0
      do L2 = L2_,1,-1
       DO JK=1,plonl      
        JADDTO(JK,L2) = JADDTO(JK,L2+1) + JADDLV(JK,L2)
       ENDDO
      enddo

!c---expanded grid now included CTM edge and mid layers plus expanded 
!c---    grid to allow for finer delta-tau at tops of clouds.
!c---    DIM of new grid = L2_ + JADDTO(1) + 1

!c---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
!c     in absence of JADDLV, L2LEV(L2) = L2
      L2LEV(:,1)  = 1
      do L2 = 2,L2_+1
       DO JK=1,plonl      
        L2LEV(JK,L2) = L2LEV(JK,L2-1) + 1 + JADDLV(JK,L2-1)
       ENDDO
      enddo

!c---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
!c---JNELEV(L=1:L_) = L2-index for top of layer L
      do L = 1,L_
       DO JK=1,plonl      
        JNDLEV(JK,L) = L2LEV(JK,2*L)
        JNELEV(JK,L) = L2LEV(JK,2*L+1)
       ENDDO 
      enddo
      JNELEV(:,L_+1) = 0  !need to set this to top-of-atmosphere
!C---------------------SET UP FOR MIE CODE-------------------------------
!c
!c  Transpose the ascending TTAU grid to a descending ZTAU grid.
!c  Double the resolution - TTAU points become the odd points on the
!c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
!c  Odd point added at top of grid for unattenuated beam   (Z='inf')
!c
!c  The following mapping holds for JADDLV=0
!c        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
!c        Top:       TTAU(L2_)  ==> ZTAU(3)
!c        Infinity:     0.0     ==> ZTAU(1)
!c        index: 2*(L2_+1-L2)+1 ==> LZ
!c
!c  Mie scattering code only used from surface to level L2_
!C------------------------------------------------------------------------
!c
!C------------------------------------------------------------------------
!c  Insert new levels, working downwards from the top of the atmosphere
!c  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
!c  to be incremented linearly, and the flux fz to be attenuated top-down 
!c    (avoiding problems where lower level fluxes are zero).
!C------------------------------------------------------------------------
!c
!c  Ascend through atmosphere transposing grid and adding extra points
!c  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
!c  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
!c    because we need to insert the intermediate layers (even LZ) for the 
!c    asymmetric scattering code.


!c  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
!c    order, expanded, doubled-level scatter grid. 
!c    Note that we need to deal with the expansion by JADD levels (L2L).
!c      These JADDLV levels are skipped and need to be interpolated later.
!c    Note that only odd LZ levels are filled, 

      NDZ(:) = 2*L2_ + 2*JADDTO(:,1) + 1
!!fix the size of N_
!!      N_ = MAXVAL(NDZ(1:plonl))+2      
!!      ALLOCATE(POMEGA(plonl,2*M_,N_,W_))
!!      ALLOCATE(ZTAU(plonl,N_,W_))
!!      ALLOCATE(FZ(plonl,N_,W_))
!!      ALLOCATE(FJ(plonl,N_,W_))      
!!initialize
!!      POMEGA(:,:,:,:)=0._dp
!!      ZTAU(:,:,:)=0._dp
!!      FZ(:,:,:)=0._dp
      
!c   Note that the successive sub-layers have the ratio in OD of ATAU
!c      ATAUA = (ATAU - 1.d0)/ATAU     ! this is the limit for L22=>inf

      DO K=NW1,NW2
       do L2 = 1,L2_+1          ! L2 = index of CTM edge- and mid-layers
        DO JK=1,plonl      
         L2L = L2LEV(JK,L2)        ! L2L = index for L2 in expanded scale(JADD)
         LZ  = NDZ(JK) + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
         ZTAU(JK,LZ,K) = TTAU(JK,L2,K)
         FZ(JK,LZ,K)   = FTAU2(JK,L2,K)
         do I=1,MFIT
          POMEGA(JK,I,LZ,K) = POMEGAJ(JK,I,L2,K)
         enddo
        ENDDO
       enddo
      ENDDO
      
!c   Now go thru the pairs of L2 levels to see if we need JADD levels
      DO K=NW1,NW2
      do L2 = 1,L2_             ! L2 = index of CTM edge- and mid-layers
       DO JK=1,plonl      
        L2L = L2LEV(JK,L2)         ! L2L = index for L2 in expanded scale(JADD)
        LZ  = NDZ(JK) + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
        L22 = L2LEV(JK,L2+1) - L2LEV(JK,L2) - 1   ! L22 = 0 if no added levels
        if (L22 .gt. 0) then
         TAUBTM = TTAU(JK,L2,K)
         TAUTOP = TTAU(JK,L2+1,K)
         FBTM   = FTAU2(JK,L2,K)
         FTOP   = FTAU2(JK,L2+1,K)
         do I = 1,MFIT
          POMEGAB(I) = POMEGAJ(JK,I,L2,K)
         enddo

!c---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
!c---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM

         ATAUZ = exp(-log(TAUBTM/max(TAUTOP,ATAU0))/float(L22+1))

         do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt arrays
          ZTAU(JK,LZZ,K) = TAUBTM * ATAUZ

          ATAUA=(TAUBTM-ZTAU(JK,LZZ,K))/(TAUBTM-TAUTOP) !fraction from TAUBTM=>TAUTOP

!c---version 5.6 fix   (mp, 3/2007)
!c---solar flux at interp-levels: use exp(TAU/U0) if U0>0, else scale by TAU
!c---solar flux at interp-levels: use exp(TAU/U0) if U0>0.02 (89 deg), 
!c---   else scale by TAU

          write(*,*) ''
          if (U0(JK) .gt. 0.02_dp) then             
           !not working with -Chopt, now -Cvopt is used 
           FZ(JK,LZZ,K) = FTOP*exp((TAUTOP-ZTAU(JK,LZZ,K))/U0(JK))
          else
           if (FBTM .lt. 1.e-32_dp) then
             FZ(JK,LZZ,K) = 0._dp
           else    
             FZ(JK,LZZ,K) = FBTM * (FTOP/FBTM)**ATAUA
           endif
          endif

!          if (U0(JK) .gt. 0.02_dp) then 
!             FZ(JK,LZZ,K) = FTOP*exp((TAUTOP-ZTAU(JK,LZZ,K))/U0(JK))
!          endif
!          if (U0(JK) .le. 0.02_dp .and. FBTM .lt. 1.e-32_dp) then
!             FZ(JK,LZZ,K) = 0._dp
!          endif    
!          if (U0(JK) .le. 0.02_dp .and. FBTM .ge. 1.e-32_dp) then         
!             FZ(JK,LZZ,K) = FBTM * (FTOP/FBTM)**ATAUA
!          endif

          do I = 1,MFIT
           POMEGA(JK,I,LZZ,K) = POMEGAB(I) + &
                    ATAUA*(POMEGAJ(JK,I,L2+1,K)-POMEGAB(I))
          enddo
          TAUBTM = ZTAU(JK,LZZ,K)
          FBTM   = FZ(JK,LZZ,K)
          do I = 1,MFIT
           POMEGAB(I) = POMEGA(JK,I,LZZ,K)
          enddo
         enddo
        endif
       ENDDO !plonl
      enddo !L2_
      ENDDO !K      

!c   Now fill in the even points with simple interpolation in scatter arrays:
      DO K=NW1,NW2
       DO JK=1,plonl
        do LZ = 2,NDZ(JK)-1,2
         ZTAU(JK,LZ,K) = 0.5_dp*(ZTAU(JK,LZ-1,K)+ZTAU(JK,LZ+1,K))
         FZ(JK,LZ,K)   = sqrt(FZ(JK,LZ-1,K)*FZ(JK,LZ+1,K))
         do I=1,MFIT
          POMEGA(JK,I,LZ,K) = 0.5_dp*(POMEGA(JK,I,LZ-1,K)+ &
       	                             POMEGA(JK,I,LZ+1,K))
         enddo
        enddo
       ENDDO 
      ENDDO        
!c---PRINT diagnostics
!c----now check integral of FZ over each CTM layer to ensure that it equals FLXD
!c      if (KW.eq.18) then
!c       WRITE(2,'(A,3I6)') 'OPMIE levels: L,TAU,Fs,pomega',KW, ND,N_
!c      do L=1,ND
!c       WRITE(2,'(i5,f15.5,1p,e15.5,0p,10f10.5)') L,
!c     & ZTAU(L),FZ(L),(POMEGA(I,L),I=1,3) 
!c      enddo
!c
!c      WRITE(2,'(a,2i6)') '  net flux in each CTM layer',KW,L_
!c      WRITE(2,'(2i5,a8,1p,e12.5)')(L,JNELEV(L),' flxd ',FLXD(L),L=1,L1_)
!c
!c      endif

       DO JK=1,plonl
       if(NDZ(JK) .gt. N_) then
        WRITE(cerr1,*) NDZ(JK)
        WRITE(cerr2,*) N_
        CALL finish('fjx, OPMIE:', &
         'overflow of scatter arrays: NDZ='//TRIM(cerr1)// &
         ', N_='//TRIM(cerr2))
       endif
       ENDDO 
!c--------------------------------------------------------------------------------
       call MIESCT_V(plonl,FJ,FJT,FJB,POMEGA, &
                     FZ,ZTAU,ZFLUX,RFLECT, &
                     U0,MFIT,NDZ)
!c--------------------------------------------------------------------------------
        

!c---Move mean intensity from scatter array FJ(LZ=1:ND) 
!c---              to CTM mid-level array FJACT(L=1:L_)

!c---mean intensity:  4*<I> + solar at mid-layer
      DO K=NW1,NW2
      do L = 1,L_
       DO JK=1,plonl 
        L2L = JNDLEV(JK,L)
        LZ  = NDZ(JK)+2 - 2*L2L
        FJACT(JK,L,K) = 4._dp*FJ(JK,LZ,K) + FZ(JK,LZ,K)
       ENDDO
      enddo
      ENDDO
!c---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
!c---      average (tau-wtd) the h's just above and below the L-edge
      DO K=NW1,NW2
      do L = 1,L_
       DO JK=1,plonl 
        L2L = JNELEV(JK,L)
        LZ  = NDZ(JK)+2 - 2*L2L
!c---       FJFLX(L) = 2.0d0*(FJ(LZ-1) + FJ(LZ+1))
         FJFLX0 = (ZTAU(JK,LZ+1,K)-ZTAU(JK,LZ,K))/(ZTAU(JK,LZ+1,K) &
       	          -ZTAU(JK,LZ-1,K))
         FJFLX(JK,L,K) = 4._dp*(FJ(JK,LZ-1,K)*FJFLX0 + &
        	       FJ(JK,LZ+1,K)*(1._dp-FJFLX0))
       ENDDO
      enddo
      ENDDO  
!c---NB if one needs the mean intensity throughout layer L (instead of mid-pt)
!c---   then average (tau-weighted) the odd-points from: NELEV(L-1) to NELEV(L)
!c---NB This is NOT now used.  Errors in cloudy layers are 0-10%  (low sun)
!c---   the results are stored locally in FJACT2 (actinic)
      DO K=NW1,NW2
      LZ1(:)=NDZ(:)
      do L = 1,L_
       DO JK=1,plonl
        LZMID = NDZ(JK)+2-2*JNDLEV(JK,L)
        LZ0 = NDZ(JK)+2-2*JNELEV(JK,L)
        SUMT = 0._dp
        SUMJ = 0._dp
	do L2 = LZ0,LZ1(JK)-2,2
         SUMT = SUMT + ZTAU(JK,L2+2,K)-ZTAU(JK,L2,K)
         SUMJ = SUMJ + (ZTAU(JK,L2+2,K)-ZTAU(JK,L2,K))*(FJ(JK,L2,K) &
      	             +FJ(JK,L2+2,K))*0.5_dp
        enddo
        FJACT2(JK,L,K) = 4._dp*SUMJ/SUMT + FZ(JK,LZMID,K)
!c-----LZ are indices to the 1:ND array     LZ1 > LZMID > LZ0
!c       WRITE(2,'(4i5,3f10.3)') L, LZ0,LZMID,LZ1,
!c     &      ZTAU(LZMID),FJACT(L),FJACT2(L)
        LZ1(JK) = LZ0
       ENDDO
      enddo
      ENDDO
!c---diffuse fluxes reflected at top, incident at bottom 
      FJTOP(:,:) = FJT(:,:)
      FJBOT(:,:) = FJB(:,:)

!++lp
!!      DEALLOCATE(POMEGA,ZTAU,FZ,FJ)
!--lp
    
      return
      end subroutine opmie_v2
!C-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
      subroutine MIESCT_V(plonl,FJ,FJT,FJB, &
                          POMEGA,FZ,ZTAU,ZFLUX,ZREFL,ZU0,MFIT,ND)
!C-----------------------------------------------------------------------
!C   This is an adaption of the Prather radiative transfer code, (mjp, 10/95)
!C     Prather, 1974, Astrophys. J. 192, 787-792.
!C         Sol_n of inhomogeneous Rayleigh scattering atmosphere. 
!C         (original Rayleigh w/ polarization)
!C     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!C         Raman scattering in the atmospheres of the major planets.
!C         (first use of anisotropic code)
!C     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!C         Chemistry of a polluted cloudy boundary layer,
!C         (documentation of extension to anisotropic scattering)
!C
!C    takes atmospheric structure and source terms from std J-code
!C    ALSO limited to 4 Gauss points, only calculates mean field!
!C
!C   mean rad. field ONLY (M=1)
!C-----------------------------------------------------------------------
      USE mo_kind,  ONLY: dp

      implicit none

!c--- expect parameters M_, N_ in parm_MIE.f------------------------------
!c
      integer, intent(in) :: MFIT, plonl, ND(plonl)
      real(dp), intent(in)  :: POMEGA(plonl,2*M_,N_,W_),FZ(plonl,N_,W_),&
                             ZTAU(plonl,N_,W_),ZREFL(plonl),&
                             ZU0(plonl),ZFLUX(plonl,W_)
      real(dp), intent(out) :: FJ(plonl,N_,W_),FJT(plonl,W_),&
                             FJB(plonl,W_)

      real(dp) :: WT(M_),EMU(M_),PM(plonl,M_,2*M_),PM0(plonl,2*M_),CMEQ1
      integer I, ID, IM, M, N
      INTEGER JK,K
      integer count1,count2
!C-----------------------------------------------------------------------
!C---fix scattering to 4 Gauss pts = 8-stream
      call GAUSSP (N,EMU,WT)
!c---calc in OPMIE:  ZFLUX = (ZU0*FZ(ND)*ZREFL)/(1.0d0+ZREFL)
      M = 1
      DO JK=1,plonl
       do I = 1,N
        call LEGND0 (EMU(I),PM0(JK,:),MFIT)
        do IM = M,MFIT
          PM(JK,I,IM) = PM0(JK,IM)
        enddo
       enddo
      ENDDO

      CMEQ1 = 0.25_dp
      DO JK=1,plonl
       call LEGND0 (-ZU0(JK),PM0(JK,:),MFIT)
       do IM=M,MFIT
        PM0(JK,IM) = CMEQ1*PM0(JK,IM)
       enddo
      ENDDO
      
      call system_clock(count=count1)
      DO K=NW1,NW2
       DO JK=1,plonl        
!C-----------------------------------------------------------------------
      call BLKSLV(FJ(JK,:,K),POMEGA(JK,:,:,K),FZ(JK,:,K),ZTAU(JK,:,K),&
		  ZFLUX(JK,K),ZREFL(JK),WT,EMU,PM(JK,:,:),PM0(JK,:),&
		  FJT(JK,K),FJB(JK,K),M,N,MFIT,ND(JK))
!C-----------------------------------------------------------------------
       ENDDO
!!      call BLKSLV_V(FJ(:,:,K),POMEGA(:,:,:,K),FZ(:,:,K),ZTAU(:,:,K),&
!!                  ZFLUX(:,K),ZREFL(:),WT,EMU,PM(:,:,:),PM0(:,:),&
!!                  FJT(:,K),FJB(:,K),M,N,MFIT,ND(:))
      ENDDO 
      call system_clock(count=count2)
      write(*,*) 'count blkslv',count2-count1
!C
!c      do ID=1,ND,2             !!!! this is now done in MIESCT()
!c        FJ(ID) = 4.0d0*FJ(ID) + FZ(ID)
!c      enddo

      return
      end subroutine miesct_v
!C-----------------------------------------------------------------------
!c-----------------------------------------------------------------------
      subroutine JRATET_V(plonl,KLEV,PPJ,TTJ,FFF,VALJL)
!c-----------------------------------------------------------------------
!c in:
!c        PPJ(L1_+1) = pressure profile at edges
!c        TTJ(L1_) = = temperatures at mid-level
!c        FFF(K=1:NW, L=1:JVL_) = mean actinic flux 
!c out:
!c        VALJL(JVL_,NJVAL)  JVL_ = no of levels
!C-----------------------------------------------------------------------
      USE mo_kind,   ONLY: dp

      implicit none

      INTEGER,  INTENT(IN)  ::  plonl,KLEV
      real(dp), intent(in)  ::  PPJ(plonl,KLEV+2),TTJ(plonl,KLEV+1)
      real(dp), intent(in)  ::  FFF(plonl,KLEV,W_)
      real(dp), intent(out) ::  VALJL(plonl,KLEV,NJVAL)

      !!real(dp)  FLINT             ! external function for X-sections
      real(dp)  VALJ(plonl,KLEV,NJVAL)
      real(dp)  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT(plonl,KLEV)
      real(dp)  TT(plonl,KLEV),PP(plonl,KLEV),DD(plonl,KLEV)
      REAL(dp)  TT200,TFACA(plonl,KLEV),TFAC0(plonl,KLEV)
      REAL(dp)  TFAC1(plonl,KLEV),TFAC2(plonl,KLEV),QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV
      INTEGER JK,JL

      do JL = 1,KLEV    
       DO JK=1,plonl
!c---need temperature and density (for some quantum yields):
!c---in this case the Pressures PPJ are defined at the boundaries,
!c---                Temperatures in the middle of each layer
        TT(JK,JL)  = TTJ(JK,JL)
        PP(JK,JL)  = (PPJ(JK,JL)+PPJ(JK,JL+1))*0.5_dp
        if (JL.eq.1) PP(JK,JL) = PPJ(JK,1)
!++mgs: !!baustelle!! TT becomes zero here ?? interface problem! 
!#debug# check for zero...
!          IF (TT(JK,JL) > 100._dp) THEN
!            DD(JK,JL) = 7.24e18_dp*PP(JK,JL)/TT(JK,JL)
!          ELSE
!            write(0,*) '#debug# fastj: TT<100.  JK, JL = ',jk,jl
!            DD(JK, JL) = 7.24e18_dp
!          END IF
!--- original statement:
       DD(JK,JL) = 7.24e18_dp*PP(JK,JL)/TT(JK,JL)
!--mgs
       ENDDO
      enddo    	

      VALJ(:,:,:) = 0._dp
        
      do K = NW1,NW2                    ! Using model 'T's here
       do JL = 1,KLEV   
        DO JK=1,plonl
           QO3TOT = FLINT(TT(JK,JL),TQQ(1,2),TQQ(2,2),TQQ(3,2) &
                             ,QO3(K,1),QO3(K,2),QO3(K,3))
           QO2TOT = FLINT(TT(JK,JL),TQQ(1,1),TQQ(2,1),TQQ(3,1) &
                             ,QO2(K,1),QO2(K,2),QO2(K,3))
           QO31DY = FLINT(TT(JK,JL),TQQ(1,3),TQQ(2,3),TQQ(3,3) &
                             ,Q1D(K,1),Q1D(K,2),Q1D(K,3))
           QO31D  = QO31DY*QO3TOT
          VALJ(JK,JL,1) = VALJ(JK,JL,1) + QO2TOT*FFF(JK,JL,K)
          VALJ(JK,JL,2) = VALJ(JK,JL,2) + QO3TOT*FFF(JK,JL,K)
          VALJ(JK,JL,3) = VALJ(JK,JL,3) + QO31D*FFF(JK,JL,K)
	ENDDO
       enddo	  
      enddo

      do J = 4,NJVAL
       if (TQQ(2,J) .gt. TQQ(1,J)) then
        DO JL = 1,KLEV 
         DO JK=1,plonl
	   TFACT(JK,JL) = max(0._dp,min(1._dp,(TT(JK,JL)-TQQ(1,J))/ &
                                            (TQQ(2,J)-TQQ(1,J))))
	 ENDDO
	ENDDO   
       else
        TFACT(:,:) = 0._dp
       endif

       do K = NW1,NW2
        DO JL = 1,KLEV 
         DO JK=1,plonl       
           QQQT    = QQQ(K,1,J) + (QQQ(K,2,J) - QQQ(K,1,J))*TFACT(JK,JL)
           VALJ(JK,JL,J) = VALJ(JK,JL,J) + QQQT*FFF(JK,JL,K)
         ENDDO
	ENDDO
       enddo


!c #52 Methylvinyl ketone   'MeVK  '     q(M) = 1/(1 + 1.67e-19*[M])
       if (TITLEJ(J).eq.'MeVK  ') then
        DO JL = 1,KLEV 
         DO JK=1,plonl       
         VALJ(JK,JL,J) = VALJ(JK,JL,J)/(1.0_dp + 1.67e-19_dp*DD(JK,JL))
	 ENDDO
	ENDDO 
       endif
!c #55 Methylethyl ketone   MEKeto     q(M) = 1/(1 + 2.0*[M/2.5e19])
       if (TITLEJ(J).eq.'MEKeto') then
        DO JL = 1,KLEV 
         DO JK=1,plonl       
         VALJ(JK,JL,J) = VALJ(JK,JL,J)/(1.0_dp + 0.80E-19_dp*DD(JK,JL))
	 ENDDO
	ENDDO 
       endif
!c #57 Methyl glyoxal       MGlyxl     q(M) = 1/(1 + 4.15*[M/2.5E19])
       if (TITLEJ(J).eq.'MGlyxl') then
        DO JL = 1,KLEV 
         DO JK=1,plonl       
         VALJ(JK,JL,J) = VALJ(JK,JL,J)/(1.0_dp + 1.66e-19_dp*DD(JK,JL))
	 ENDDO
	ENDDO 
       endif

      enddo !njval

      if (TITLEJ(NJVAL-1).eq.'Acet-a') then
!c--------------J-ref v8.3 includes Blitz ACETONE q-yields--------------
!c---Acetone is a special case:   (as per Blitz et al GRL, 2004)
!c---	  61 = NJVAL-1 = J1(acetone-a) ==> CH3CO + CH3
!c---	  62 = NJVAL   = J2(acetone-b) ==> CH3 + CO + CH3
        VALJ(:,:,NJVAL-1) = 0._dp
        VALJ(:,:,NJVAL)   = 0._dp
!c---IV=NJVAL-1 = Xsect (total abs) for Acetone - pre-calc Temp interp factors
	 IV    = NJVAL-1
	 DO JL = 1,KLEV 
	  DO JK=1,plonl	      
	  TFACA(JK,JL) = (TT(JK,JL)-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
	   TFACA(JK,JL) = max(0._dp, min(1._dp, TFACA(JK,JL)))
	 ENDDO
	ENDDO  
!c---IV=NJVAL = Q2 for Acetone=>(2), specifically designed for quadratic interp.
!c---	   but force to Q2=0 by 210K
	 IV    = NJVAL
	 DO JL = 1,KLEV 
	  DO JK=1,plonl	      
	   TFAC0(JK,JL) = ( (TT(JK,JL)-TQQ(1,IV))/&
			   (TQQ(2,IV)-TQQ(1,IV)) )**2
	   if (TT(JK,JL) .lt. TQQ(1,IV)) then
	    TFAC0(JK,JL) = (TT(JK,JL) - 210._dp)/(TQQ(1,IV)-210._dp)
	   endif
	   TFAC0(JK,JL) = max(0._dp, min(1._dp, TFAC0(JK,JL)))
	 ENDDO
	ENDDO  
!c---IV=NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
	 IV    = NJVAL+1
	 DO JL = 1,KLEV 
	  DO JK=1,plonl	      
	   TT200 = min(300._dp, max(200._dp, TT(JK,JL)))
	   TFAC1(JK,JL) = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
	  ENDDO
	ENDDO 
!c---IV=NJVAL+2 = Q1B for Acetone => (1)
	 IV    = NJVAL+2
	 DO JL = 1,KLEV 
	  DO JK=1,plonl	      
	   TT200 = min(300._dp, max(200._dp, TT(JK,JL)))
	   TFAC2(JK,JL) = (TT200-TQQ(1,IV))/(TQQ(2,IV)-TQQ(1,IV))
	  ENDDO
	ENDDO  
!c---now integrate over wavelengths
	do K = NW1,NW2
	 DO JL = 1,KLEV 
	  DO JK=1,plonl	      
!c---NJVAL-1 = Xsect (total abs) for Acetone
	  IV   = NJVAL-1
	  QQQA = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFACA(JK,JL)
!c---NJVAL   = Q2 for Acetone=>(2), specifically designed for quadratic interp.
	  IV   = NJVAL
	  QQ2  = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC0(JK,JL)
	  if (TT(JK,JL) .lt. TQQ(1,IV)) then
	    QQ2 = QQQ(K,1,IV)*TFAC0(JK,JL)
	  endif
!c---NJVAL+1 = Q1A for Acetone => (1), allow full range of T = 200K-300K
	  IV   = NJVAL+1
	  QQ1A = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC1(JK,JL)
!c---NJVAL+2 = Q1B for Acetone => (1)	! scaled to [M]=2.5e19
	  IV   = NJVAL+2
	  QQ1B = QQQ(K,1,IV) + (QQQ(K,2,IV)-QQQ(K,1,IV))*TFAC2(JK,JL)
	  QQ1B = QQ1B*4.e-20_dp
!c---J(61)
	  VALJ(JK,JL,NJVAL-1) = VALJ(JK,JL,NJVAL-1) &
	       + FFF(JK,JL,K)*QQQA*(1._dp-QQ2)/(QQ1A + QQ1B*DD(JK,JL))
!c---J(62)
	  VALJ(JK,JL,NJVAL) = VALJ(JK,JL,NJVAL) + FFF(JK,JL,K)*QQQA*QQ2
	 ENDDO 
	 ENDDO
	enddo	 !K
!c-----------end v-8.3 includes Blitz ACETONE q-yields--------------
     endif

!c----Load array of J-values in native order, need to be indexed/scaled
!c    by ASAD-related code later: ZPJ(L,JJ) = VALJL(L,JIND(JJ))*JFACTA(JJ)
      VALJL(:,:,:) = VALJ(:,:,:)

      return
      end subroutine jratet_v
!c-----------------------------------------------------------------------
!C-----------------------------------------------------------------------
      subroutine BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
                       ,FJTOP,FJBOT,  M,N,MFIT,ND)
!C-----------------------------------------------------------------------
!C  Sets up and solves the block tri-diagonal system:  
!C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!c  This goes back to the old, dumb, fast version 5.3
!C-----------------------------------------------------------------------
      USE mo_kind,   ONLY: dp

      implicit none

!c--- expect parameters M_, N_ in parm_MIE.f------------------------------

      integer, intent(in) ::  M, N, MFIT, ND
      real(dp), intent(in)  ::  POMEGA(2*M_,N_),FZ(N_),ZTAU(N_) &
                             ,WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_) &
                             ,ZFLUX,ZREFL
      real(dp), intent(out) ::  FJ(N_),FJTOP,FJBOT

      real(dp), dimension(M_)    :: A, C1, H
      real(dp), dimension(M_,M_) :: B, AA, CC
      real(dp)                      DD(M_,M_,N_), RR(M_,N_)
      real(dp)  SUMM, FIPLUS
      integer I, J, K, ID

!c
!C-----------UPPER BOUNDARY ID=1
      call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
              ,B,CC,AA,A,H,C1,M,N,MFIT,ND,1)
      call MATIN4 (B)
      do I = 1,N
         RR(I,1) = 0.0_dp
        do J = 1,N
          SUMM = 0.0_dp
         do K = 1,N
          SUMM = SUMM - B(I,K)*CC(K,J)
         enddo
         DD(I,J,1) = SUMM
         RR(I,1) = RR(I,1) + B(I,J)*H(J)
        enddo
      enddo
!C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      do ID = 2,ND-1
        call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
                ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ID)
        do I = 1,N
          do J = 1,N
          B(I,J) = B(I,J) + A(I)*DD(I,J,ID-1)
          enddo
          H(I) = H(I) - A(I)*RR(I,ID-1)
        enddo
        call MATIN4 (B)
        do I = 1,N
          RR(I,ID) = 0.0_dp
          do J = 1,N
          RR(I,ID) = RR(I,ID) + B(I,J)*H(J)
          DD(I,J,ID) = - B(I,J)*C1(J)
          enddo
        enddo
      enddo
!C---------FINAL DEPTH POINT: ND
      call GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
              ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ND)
      do I = 1,N
        do J = 1,N
          SUMM = 0.0_dp
          do K = 1,N
          SUMM = SUMM + AA(I,K)*DD(K,J,ND-1)
          enddo
        B(I,J) = B(I,J) + SUMM
        H(I) = H(I) - AA(I,J)*RR(J,ND-1)
        enddo
      enddo
      call MATIN4 (B)
      do I = 1,N
        RR(I,ND) = 0.0_dp
        do J = 1,N
        RR(I,ND) = RR(I,ND) + B(I,J)*H(J)
        enddo
      enddo
!C-----------BACK SOLUTION
      do ID = ND-1,1,-1
       do I = 1,N
        do J = 1,N
         RR(I,ID) = RR(I,ID) + DD(I,J,ID)*RR(J,ID+1)
        enddo
       enddo
      enddo

!C----------MEAN J & H
      do ID = 1,ND,2
        FJ(ID) = 0.0_dp
       do I = 1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)
       enddo
      enddo
      do ID = 2,ND,2
        FJ(ID) = 0.0_dp
       do I = 1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)*EMU(I)
       enddo
      enddo

!c---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
!c---FJBOT = scaled diffuse flux onto surface:  
        FJTOP = 0.0_dp
        FJBOT = 0.0_dp
       do I = 1,N
        FJTOP = FJTOP + RR(I,1)*WT(I)*EMU(I)
        FJBOT = FJBOT + RR(I,ND)*WT(I)*EMU(I)
       enddo
        FJTOP = 4._dp*FJTOP
!c---         ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
        FIPLUS = 4._dp*ZREFL*FJBOT/(1.0_dp + ZREFL) + ZFLUX
        FJBOT = 4._dp*FJBOT - FIPLUS

      return
      end subroutine blkslv
!----------------------------------------------------------------
!C-----------------------------------------------------------------------
      subroutine GEN(POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0 &
                    ,B,CC,AA,A,H,C1,M,N,MFIT,ND,ID)
!C-----------------------------------------------------------------------
!C  Generates coefficient matrices for the block tri-diagonal system:
!C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!C-----------------------------------------------------------------------
      USE mo_kind,   ONLY: dp

      implicit none

!c--- expect parameters M_, N_ in parm_MIE.f------------------------------
      integer, intent(in) ::  M, N, MFIT, ND, ID
      real(dp), intent(in)  ::  POMEGA(2*M_,N_),FZ(N_),ZTAU(N_) &
                             ,WT(M_),EMU(M_),PM(M_,2*M_),PM0(2*M_) &
                             ,ZFLUX,ZREFL
      real(dp), intent(out) ::  B(M_,M_),AA(M_,M_),CC(M_,M_),A(M_),C1(M_) &
                             ,H(M_)

      integer ID0, ID1, IM, I, J, K, MSTART
      real(dp)  SUM0, SUM1, SUM2, SUM3
      real(dp)  DELTAU, D1, D2, SURFAC

      real(dp)  S(M_,M_), W(M_,M_), U1(M_,M_), V1(M_)
!C---------------------------------------------
      if (ID.eq.1 .or. ID.eq.ND) then
!C---------calculate generic 2nd-order terms for boundaries
       ID0 = ID
       ID1 = ID+1
       if (ID .ge. ND) ID1 = ID-1

       do I = 1,N
          SUM0 = 0.0_dp
          SUM1 = 0.0_dp
          SUM2 = 0.0_dp
          SUM3 = 0.0_dp
        do IM = M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        enddo
        do IM = M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM0(IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM0(IM)
        enddo
         H(I) = 0.5_dp*(SUM0*FZ(ID0) + SUM2*FZ(ID1))
         A(I) = 0.5_dp*(SUM1*FZ(ID0) + SUM3*FZ(ID1))
        do J = 1,I
          SUM0 = 0.0_dp
          SUM1 = 0.0_dp
          SUM2 = 0.0_dp
          SUM3 = 0.0_dp
         do IM = M,MFIT,2
          SUM0 = SUM0 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM2 = SUM2 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         enddo
         do IM = M+1,MFIT,2
          SUM1 = SUM1 + POMEGA(IM,ID0)*PM(I,IM)*PM(J,IM)
          SUM3 = SUM3 + POMEGA(IM,ID1)*PM(I,IM)*PM(J,IM)
         enddo
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         W(I,J) = - SUM1*WT(J)
         W(J,I) = - SUM1*WT(I)
         U1(I,J) = - SUM3*WT(J)
         U1(J,I) = - SUM3*WT(I)
          SUM0 = 0.5_dp*(SUM0 + SUM2)
         B(I,J) = - SUM0*WT(J)
         B(J,I) = - SUM0*WT(I)
        enddo
         S(I,I) = S(I,I) + 1.0_dp
         W(I,I) = W(I,I) + 1.0_dp
         U1(I,I) = U1(I,I) + 1.0_dp
         B(I,I) = B(I,I) + 1.0_dp
       enddo

       do I = 1,N
         SUM0 = 0.0_dp
        do J = 1,N
         SUM0 = SUM0 + S(I,J)*A(J)/EMU(J)
        enddo
        C1(I) = SUM0
       enddo
       do I = 1,N
        do J = 1,N
          SUM0 = 0.0_dp
          SUM2 = 0.0_dp
         do K = 1,N
          SUM0 = SUM0 + S(J,K)*W(K,I)/EMU(K)
          SUM2 = SUM2 + S(J,K)*U1(K,I)/EMU(K)
         enddo
         A(J) = SUM0
         V1(J) = SUM2
        enddo
        do J = 1,N
         W(J,I) = A(J)
         U1(J,I) = V1(J)
        enddo
       enddo
       if (ID .eq. 1) then
!C-------------upper boundary, 2nd-order, C-matrix is full (CC)
        DELTAU = ZTAU(2) - ZTAU(1)
        D2 = 0.25_dp*DELTAU
        do I = 1,N
          D1 = EMU(I)/DELTAU
          do J = 1,N
           B(I,J) = B(I,J) + D2*W(I,J)
           CC(I,J) = D2*U1(I,J)
          enddo
          B(I,I) = B(I,I) + D1
          CC(I,I) = CC(I,I) - D1
!C         H(I) = H(I) + 2.0_dp*D2*C1(I) + D1*SISOTP
          H(I) = H(I) + 2.0_dp*D2*C1(I)
          A(I) = 0.0_dp
        enddo
       else
!C-------------lower boundary, 2nd-order, A-matrix is full (AA)
        DELTAU = ZTAU(ND) - ZTAU(ND-1)
        D2 = 0.25_dp*DELTAU
        SURFAC = 4.0_dp*ZREFL/(1.0_dp + ZREFL)
        do I = 1,N
          D1 = EMU(I)/DELTAU
          H(I) = H(I) - 2.0_dp*D2*C1(I)
           SUM0 = 0.0_dp
          do J = 1,N
           SUM0 = SUM0 + W(I,J)
          enddo
           SUM0 = D1 + D2*SUM0
           SUM1 = SURFAC*SUM0
          do J = 1,N
           B(I,J) = B(I,J) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
          enddo
          B(I,I) = B(I,I) + D1
          H(I) = H(I) + SUM0*ZFLUX
          do J = 1,N
           AA(I,J) = - D2*U1(I,J)
          enddo
           AA(I,I) = AA(I,I) + D1
           C1(I) = 0.0_dp
        enddo
       endif
!C------------intermediate points:  can be even or odd, A & C diagonal
      else
        DELTAU = ZTAU(ID+1) - ZTAU(ID-1)
        MSTART = M + MOD(ID+1,2)
        do I = 1,N
          A(I) = EMU(I)/DELTAU
          C1(I) = -A(I)
           SUM0 = 0.0_dp
          do IM = MSTART,MFIT,2
           SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM0(IM)
          enddo
          H(I) = SUM0*FZ(ID)
          do J=1,I
            SUM0 = 0.0_dp
           do IM = MSTART,MFIT,2
            SUM0 = SUM0 + POMEGA(IM,ID)*PM(I,IM)*PM(J,IM)
           enddo
            B(I,J) =  - SUM0*WT(J)
            B(J,I) =  - SUM0*WT(I)
          enddo
          B(I,I) = B(I,I) + 1.0_dp
        enddo
      endif
      return
      end subroutine gen
!C-----------------------------------------------------------------------
      subroutine MATIN4 (A)
!C-----------------------------------------------------------------------
!C  invert 4x4 matrix A(4,4) in place with L-U decomposition (mjp, old...)
      USE mo_kind, ONLY: dp
 
      implicit none

      real(dp), intent(inout)  ::  A(4,4)
!C---SETUP L AND U
      A(2,1) = A(2,1)/A(1,1)
      A(2,2) = A(2,2)-A(2,1)*A(1,2)
      A(2,3) = A(2,3)-A(2,1)*A(1,3)
      A(2,4) = A(2,4)-A(2,1)*A(1,4)
      A(3,1) = A(3,1)/A(1,1)
      A(3,2) = (A(3,2)-A(3,1)*A(1,2))/A(2,2)
      A(3,3) = A(3,3)-A(3,1)*A(1,3)-A(3,2)*A(2,3)
      A(3,4) = A(3,4)-A(3,1)*A(1,4)-A(3,2)*A(2,4)
      A(4,1) = A(4,1)/A(1,1)
      A(4,2) = (A(4,2)-A(4,1)*A(1,2))/A(2,2)
      A(4,3) = (A(4,3)-A(4,1)*A(1,3)-A(4,2)*A(2,3))/A(3,3)
      A(4,4) = A(4,4)-A(4,1)*A(1,4)-A(4,2)*A(2,4)-A(4,3)*A(3,4)
!C---INVERT L
      A(4,3) = -A(4,3)
      A(4,2) = -A(4,2)-A(4,3)*A(3,2)
      A(4,1) = -A(4,1)-A(4,2)*A(2,1)-A(4,3)*A(3,1)
      A(3,2) = -A(3,2)
      A(3,1) = -A(3,1)-A(3,2)*A(2,1)
      A(2,1) = -A(2,1)
!C---INVERT U
      A(4,4) = 1._dp/A(4,4)
      A(3,4) = -A(3,4)*A(4,4)/A(3,3)
      A(3,3) = 1._dp/A(3,3)
      A(2,4) = -(A(2,3)*A(3,4)+A(2,4)*A(4,4))/A(2,2)
      A(2,3) = -A(2,3)*A(3,3)/A(2,2)
      A(2,2) = 1._dp/A(2,2)
      A(1,4) = -(A(1,2)*A(2,4)+A(1,3)*A(3,4)+A(1,4)*A(4,4))/A(1,1)
      A(1,3) = -(A(1,2)*A(2,3)+A(1,3)*A(3,3))/A(1,1)
      A(1,2) = -A(1,2)*A(2,2)/A(1,1)
      A(1,1) = 1._dp/A(1,1)
!C---MULTIPLY (U-INVERSE)*(L-INVERSE)
      A(1,1) = A(1,1)+A(1,2)*A(2,1)+A(1,3)*A(3,1)+A(1,4)*A(4,1)
      A(1,2) = A(1,2)+A(1,3)*A(3,2)+A(1,4)*A(4,2)
      A(1,3) = A(1,3)+A(1,4)*A(4,3)
      A(2,1) = A(2,2)*A(2,1)+A(2,3)*A(3,1)+A(2,4)*A(4,1)
      A(2,2) = A(2,2)+A(2,3)*A(3,2)+A(2,4)*A(4,2)
      A(2,3) = A(2,3)+A(2,4)*A(4,3)
      A(3,1) = A(3,3)*A(3,1)+A(3,4)*A(4,1)
      A(3,2) = A(3,3)*A(3,2)+A(3,4)*A(4,2)
      A(3,3) = A(3,3)+A(3,4)*A(4,3)
      A(4,1) = A(4,4)*A(4,1)
      A(4,2) = A(4,4)*A(4,2)
      A(4,3) = A(4,4)*A(4,3)

      return
      end subroutine matin4
!C----------------------------------------------------------------------

      FUNCTION XSECO3(K,TTT)
!-----------------------------------------------------------------------
!  Cross-sections for O3 for all processes interpolated across 3 temps
! ----------------------------------------------------------------------
!     include 'cmn_h.f'
!     include 'jv_cmn.h'
      USE mo_kind,     ONLY: dp

      integer :: k
      REAL(dp) :: ttt, xseco3!!, flint
      XSECO3  = &
       FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))
      return
      end function xseco3


      FUNCTION XSEC1D(K,TTT)
!-----------------------------------------------------------------------
!  Quantum yields for O3 --> O2 + O(1D) interpolated across 3 temps
!-----------------------------------------------------------------------
!     include 'cmn_h.f'
!     include 'jv_cmn.h'
      USE mo_kind,     ONLY: dp

      integer :: k
      REAL(dp) :: ttt, xsec1d!!, flint
      XSEC1D = FLINT(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))
      return
      end function xsec1d


      FUNCTION XSECO2(K,TTT)
!-----------------------------------------------------------------------
!  Cross-sections for O2 interpolated across 3 temps; No S_R Bands yet!
!-----------------------------------------------------------------------
!      include 'cmn_h.f'
!      include 'jv_cmn.h'
      USE mo_kind,     ONLY: dp

      integer :: k
      REAL(dp) :: ttt, xseco2!!, flint
      XSECO2 = &
       FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))
      return
      end function xseco2


      FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
!-----------------------------------------------------------------------
!  Three-point linear interpolation function
!-----------------------------------------------------------------------
      USE mo_kind,     ONLY: dp

      REAL(dp) :: TINT,T1,T2,T3,F1,F2,F3,flint

      IF (TINT .LE. T2)  THEN
        IF (TINT .LE. T1)  THEN
          FLINT  = F1
        ELSE
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT .GE. T3)  THEN
          FLINT  = F3
        ELSE
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      return
      end function flint

!-----------------------------------------------------------------------
      SUBROUTINE LEGND0(X,PL,N)
!---Calculates ORDINARY LEGENDRE fns of X (real) 
!---   from P[0] = PL(1) = 1,  P[1] = X, .... P[N-1] = PL(N)
      USE mo_kind,  ONLY: dp   

      IMPLICIT NONE

      INTEGER:: N,I
      REAL(dp) :: X,PL(N),DEN
!---Always does PL(2) = P[1]
        PL(1) = 1._dp
        PL(2) = X
        DO I=3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2._dp-1._dp/DEN) - PL(I-2)*(1._dp-1._dp/DEN)
        ENDDO
      RETURN
      END SUBROUTINE LEGND0
!-----------------------------------------------------------------------
      subroutine GAUSSP (N,XPT,XWT)
!C-----------------------------------------------------------------------
!C  Loads in pre-set Gauss points for 4 angles from 0 to +1 in cos(theta)=mu
      USE mo_kind,  ONLY: dp   
 
      implicit none

      integer, intent(out) :: N
      real(dp), intent(out)  ::  XPT(*),XWT(*)
      real(dp)   GPT4(4),GWT4(4)
      integer  I
      data GPT4/.06943184420297_dp,.33000947820757_dp,.66999052179243_dp,&
                .93056815579703_dp/
      data GWT4/.17392742256873_dp,.32607257743127_dp,.32607257743127_dp,&
                .17392742256873_dp/
      N = 4
      do I = 1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      enddo
      return
      end subroutine gaussp
!C-----------------------------------------------------------------------

!c------------------------------------------------------------------------------
      subroutine OPTICDW (plonl,KLEV,OPTD,SSALB,SLEG,PATH,REFF,CCOVER)
!c------------------------------------------------------------------------------
!c---for the UCI aerosol/CLOUD data sets, 
!c---calculate fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!c---UCI aersols optical data  v-6.1:
!c
!c   04 W_H01 (H1/Deirm.)GAMMA:r-mode=0.1/alfa=2 r=.25 s=.20 n=1.335 ref=0.25
!c   05 W_H04 (H1/Deirm.)GAMMA:r-mode=0.4/alfa=2 r=1.0 s=.20 n=1.335 ref=1.00
!c   06 W_H40 (H1/Deir)GAMMA:r-m=4.0/alf=2 n=1.335     reff=10.00____rho=1.000
!c   07 W_C02 (C1/Deirm.)GAMMA:r-mode=2.0/alfa=6 r=3.0 s=.111n=1.335 ref=3.0   
!c   08 W_C04 (C1/Deirm.)GAMMA:r-mode=4.0/alfa=6 r=6.0 s=.111n=1.335 ref=6.0   
!c   09 W_C08 (C1/Deirm.)GAMMA:r-mode=8.0/alfa=6 r=12. s=.111n=1.335 ref=12.   
!c   10 W_C13 (C1/Deirm.)GAMMA:r-mode=13.3alfa=6 r=3.0 s=.111n=1.335 ref=20.   
!c   11 W_L06 (W/Lacis)GAM:r-mode=5.5 alfa=11/3 r=10. s=.15 n=1.335 ref=10.   
      
      USE mo_kind,     ONLY: dp

      implicit none  

!!      INTEGER, PARAMETER :: nwv = 1     !!mgs!! QUICK FIX

      integer,  intent(in) ::    plonl,KLEV
      real(dp), intent(inout)::  OPTD(plonl,klev,nwv)  ! optical depth of layer
      real(dp), intent(inout)::  SSALB(plonl,klev,nwv) ! single-scattering albedo
      real(dp), intent(inout)::  SLEG(plonl,klev,nwv,8)! scatt phase fn (Leg coeffs)
      real(dp), intent(in)::     PATH(plonl,klev)      ! path (g/m2) of cloud
      real(dp), intent(in)::     REFF(plonl,klev)      ! water cloud effective radius
      real(dp), intent(in)::     CCOVER(plonl,klev)      ! cloud cover
        
      integer :: I,J,L, JL, JK
      real(dp):: XTINCT, RHO
      integer :: IDXCLD(plonl,klev)
      real(dp):: CCFACT(plonl,klev)

      !use fixed L for all clouds for the moment
!      L=9 !c   09 W_C08 (C1/Deirm.)GAMMA:r-mode=8.0/alfa=6 r=12. s=.111n=1.335 ref=12.
      IDXCLD(:,:)=0  
      WHERE(reff(1:plonl,1:klev) <= 8._dp)
       IDXCLD(1:plonl,1:klev)=8       
      ENDWHERE
      WHERE(reff(1:plonl,1:klev) > 8._dp .AND. reff(1:plonl,1:klev) <= 16._dp)
       IDXCLD(1:plonl,1:klev)=9       
      ENDWHERE
      WHERE(reff(1:plonl,1:klev) > 16._dp)
       IDXCLD(1:plonl,1:klev)=10       
      ENDWHERE
         
      RHO = 1._dp !liquid water density

      do J=1,nwv
!c---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
!       XTINCT = 0.75_dp*QAA(J,L)/(REFF*RHO)
       do jl=1,klev
        do jk=1,plonl 
         
	 !select a cloud droplet distribution (C1/Deirm.) depending on the eff. radius
	 L=IDXCLD(jk,jl)

         !approx of the max/random overlap scheme (Feng,2004)    
         ccfact(jk,jl)=CCOVER(jk,jl)**1.5_dp
	 
!!baustelle!! code crashes here because J(=nwv) exceeds 5 ! 
         XTINCT = 0.75_dp*QAA(J,L)/(REFF(jk,jl)*RHO)

         OPTD(jk,jl,J)   = OPTD(jk,jl,J) + PATH(jk,jl)*XTINCT*CCFACT(jk,jl)
 
         SSALB(jk,jl,J)  = SSALB(jk,jl,J)+ SAA(J,L)*PATH(jk,jl)*XTINCT*CCFACT(jk,jl)

         do I=1,8
          SLEG(jk,jl,J,I) = SLEG(jk,jl,J,I)+PAA(I,J,L)*SAA(J,L)*PATH(jk,jl)*XTINCT*CCFACT(jk,jl)
         enddo
	 
	enddo
       enddo                     
      enddo

      return
      end subroutine opticdw

!--------------------------------------------------------------------------------
!c------------------------------------------------------------------------------
      subroutine OPTICDI (plonl,KLEV,OPTD,SSALB,SLEG,PATH,REFF,CCOVER)
!c------------------------------------------------------------------------------
!c---for the UCI aerosol/CLOUD data sets, 
!c---calculate fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
!c---UCI aersols optical data  v-6.1:
!c   12 Ice-H = hexagonal ice cloud (Mishchenko)                                     
!c   13 Ice-I = irregular ice cloud (Mishchenko)                                     
      
      USE mo_kind,     ONLY: dp

      implicit none  
!!   INTEGER, PARAMETER :: nwv = 1     !!mgs!! QUICK FIX

      integer,  intent(in) ::    plonl,KLEV
      real(dp), intent(inout)::  OPTD(plonl,klev,nwv)  ! optical depth of layer
      real(dp), intent(inout)::  SSALB(plonl,klev,nwv) ! single-scattering albedo
      real(dp), intent(inout)::  SLEG(plonl,klev,nwv,8)! scatt phase fn (Leg coeffs)
      real(dp), intent(in)::     PATH(plonl,klev)      ! path (g/m2) of cloud
      real(dp), intent(in)::     REFF(plonl,klev)      ! water cloud effective radius
      real(dp), intent(in)::     CCOVER(plonl,klev)      ! cloud cover
        
      integer :: I,J,L, JL, JK
      real(dp):: XTINCT, RHO
      real(dp):: CCFACT(plonl,klev)

      !use fixed L for all ice clouds
      L=13 !c   13 Ice-I = irregular ice cloud (Mishchenko)
         
      RHO = DAA(L) !ice crystal density

      do J=1,nwv
!c---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
!       XTINCT = 0.75_dp*QAA(J,L)/(REFF*RHO)
       do jl=1,klev
        do jk=1,plonl 
         
         !approx of the max/random overlap scheme (Feng,2004)    
         ccfact(jk,jl)=CCOVER(jk,jl)**1.5_dp
	  
         XTINCT = 0.75_dp*QAA(J,L)/(REFF(jk,jl)*RHO)

         OPTD(jk,jl,J)   = OPTD(jk,jl,J) + PATH(jk,jl)*XTINCT*CCFACT(jk,jl)
 
         SSALB(jk,jl,J)  = SSALB(jk,jl,J)+ SAA(J,L)*PATH(jk,jl)*XTINCT*CCFACT(jk,jl)

         do I=1,8
          SLEG(jk,jl,J,I) = SLEG(jk,jl,J,I)+PAA(I,J,L)*SAA(J,L)*PATH(jk,jl)*XTINCT*CCFACT(jk,jl)
         enddo
	 
	enddo
       enddo                     
      enddo

      return
      end subroutine opticdi

!--------------------------------------------------------------------------------
!c------------------------------------------------------------------------------
      subroutine OPTICDA (plonl,KLEV,OPTD,SSALB,SLEG)
!c------------------------------------------------------------------------------
! Add aerosol optical properties to the input fields of FASTJX
! Loop over the 7 modes of HAM
! Note: reversing levels of HAM fields on the fly
      
    USE mo_submodel,    ONLY: lham
!++mgs 20130304 ### !!! Must remove all explicit references to HAM !!
    USE mo_ham,         ONLY: nclass
!--mgs

      implicit none  

      integer,  intent(in) ::    plonl,KLEV
      real(dp), intent(inout)::  OPTD(plonl,klev,nwv)  ! optical depth of layer
      real(dp), intent(inout)::  SSALB(plonl,klev,nwv) ! single-scattering albedo
      real(dp), intent(inout)::  SLEG(plonl,klev,nwv,8)! scatt phase fn (Leg coeffs)
        
      integer :: I,J,JL,JK,JCLASS

      !!baustelle!! For now only compute aerosol OD if HAM is active. Later extend to use of Tanre climatology ?!
      IF ( lham ) THEN
         do jclass=1,nclass
          do J=1,nwv
           do jl=1,klev
            do jk=1,plonl 
!! write(0,*) '#debug#  OPTICDA in mo_fastj: crazy message...'
            
             OPTD(jk,jl,J)    = OPTD(jk,jl,J) + ODAER_FJ(jk,klev-jl+1,J,jclass)
    
             SSALB(jk,jl,J)   = SSALB(jk,jl,J)+ SSAER_FJ(jk,klev-jl+1,J,jclass)  &
                                               *ODAER_FJ(jk,klev-jl+1,J,jclass)
   
             do I=1,8
              SLEG(jk,jl,J,I) = SLEG(jk,jl,J,I)    &
                          + PP_FJ(jk,klev-jl+1,J,jclass,I)  &
                           *SSAER_FJ(jk,klev-jl+1,J,jclass)*ODAER_FJ(jk,klev-jl+1,J,jclass)
             enddo
	    
	    enddo
	   enddo
          enddo                     
         enddo
      END IF
   
      return
      end subroutine opticda

!c-----------------------------------------------------------------------
      subroutine BLKSLV_V(plonl, &
                          FJ,POMEGA,FZ,ZTAU,ZFLUX,ZREFL,WT,EMU,PM,PM0,&
                          FJTOP,FJBOT,  M,N,MFIT,ND)
!C-----------------------------------------------------------------------
!C  Sets up and solves the block tri-diagonal system:  
!C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!c  This goes back to the old, dumb, fast version 5.3
!C-----------------------------------------------------------------------
      USE mo_kind,  ONLY : dp

      implicit none
!c--- expect parameters M_, N_ in parm_MIE.f------------------------------

      INTEGER,  INTENT(IN) :: plonl
      integer,  intent(in) :: M, N, MFIT, ND(plonl)
      real(dp), intent(in) :: POMEGA(plonl,2*M_,N_,W_),FZ(plonl,N_,W_),&
                            ZTAU(plonl,N_,W_),&
                            WT(M_),EMU(M_),PM(plonl,M_,2*M_),&
                            PM0(plonl,2*M_),&
                            ZFLUX(plonl,W_),ZREFL(plonl)
      real(dp), intent(out):: FJ(plonl,N_,W_),&
                            FJTOP(plonl,W_),FJBOT(plonl,W_)

      real(dp), dimension(plonl,M_,W_)    :: A, C1, H
      real(dp), dimension(plonl,M_,M_,W_) :: B, AA, CC
      real(dp)  DD(plonl,M_,M_,N_,W_), RR(plonl,M_,N_,W_)
      real(dp)  SUMM, FIPLUS
      integer I, J, K, ID
      INTEGER JK,NW

!C-----------UPPER BOUNDARY ID=1
      DO NW=NW1,NW2
      DO JK=1,plonl
      call GEN(POMEGA(JK,:,:,NW),FZ(JK,:,NW),ZTAU(JK,:,NW),ZFLUX(JK,NW),&
              ZREFL(JK),WT,EMU,PM(JK,:,:),PM0(JK,:),&
              B(JK,:,:,NW),CC(JK,:,:,NW),AA(JK,:,:,NW),&
              A(JK,:,NW),H(JK,:,NW),C1(JK,:,NW),M,N,MFIT,ND(JK),1)
      ENDDO
      ENDDO

      DO NW=NW1,NW2
       DO JK=1,plonl
        call MATIN4 (B(JK,:,:,NW))
       ENDDO
      ENDDO
      do I = 1,N
        RR(:,I,1,:) = 0.0_dp
        do J = 1,N
         DO NW=NW1,NW2
	  DO JK=1,plonl
  	   SUMM = 0.0_dp
           do K = 1,N
            SUMM = SUMM - B(JK,I,K,NW)*CC(JK,K,J,NW)
           enddo
           DD(JK,I,J,1,NW) = SUMM
           RR(JK,I,1,NW) = RR(JK,I,1,NW) + B(JK,I,J,NW)*H(JK,J,NW)
          ENDDO
	 ENDDO 
	enddo
      enddo
!C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      DO NW=NW1,NW2
       DO JK=1,plonl
       do ID = 2,ND(JK)-1
        call GEN(POMEGA(JK,:,:,NW),FZ(JK,:,NW),ZTAU(JK,:,NW),&
            ZFLUX(JK,NW),ZREFL(JK),WT,EMU,PM(JK,:,:),PM0(JK,:),&
            B(JK,:,:,NW),CC(JK,:,:,NW),AA(JK,:,:,NW),&
            A(JK,:,NW),H(JK,:,NW),C1(JK,:,NW),M,N,MFIT,ND(JK),ID)

        do I = 1,N
          do J = 1,N
          B(JK,I,J,NW) = B(JK,I,J,NW) + A(JK,I,NW)*DD(JK,I,J,ID-1,NW)
          enddo
          H(JK,I,NW) = H(JK,I,NW) - A(JK,I,NW)*RR(JK,I,ID-1,NW)
        enddo
        call MATIN4 (B(JK,:,:,NW))
        do I = 1,N
          RR(JK,I,ID,NW) = 0.0_dp
          do J = 1,N
          RR(JK,I,ID,NW) = RR(JK,I,ID,NW) + B(JK,I,J,NW)*H(JK,J,NW)
          DD(JK,I,J,ID,NW) = - B(JK,I,J,NW)*C1(JK,J,NW)
          enddo
        enddo
       enddo
       ENDDO
      ENDDO
!C---------FINAL DEPTH POINT: ND
      DO NW=NW1,NW2
       DO JK=1,plonl
        call GEN(POMEGA(JK,:,:,NW),FZ(JK,:,NW),ZTAU(JK,:,NW),&
           ZFLUX(JK,NW),ZREFL(JK),WT,EMU,PM(JK,:,:),PM0(JK,:),&
           B(JK,:,:,NW),CC(JK,:,:,NW),AA(JK,:,:,NW),&
           A(JK,:,NW),H(JK,:,NW),C1(JK,:,NW),M,N,MFIT,ND(JK),ND(JK))
       ENDDO
      ENDDO 
      do I = 1,N
        do J = 1,N
         DO NW=NW1,NW2
	  DO JK=1,plonl
          SUMM = 0.0_dp
          do K = 1,N
          SUMM = SUMM + AA(JK,I,K,NW)*DD(JK,K,J,ND(JK)-1,NW)
          enddo
          B(JK,I,J,NW) = B(JK,I,J,NW) + SUMM
          H(JK,I,NW) = H(JK,I,NW) - AA(JK,I,J,NW)*RR(JK,J,ND(JK)-1,NW)
	  ENDDO
	 ENDDO 
        enddo
      enddo
      DO NW=NW1,NW2
       DO JK=1,plonl
        call MATIN4 (B(JK,:,:,NW))
       ENDDO
      ENDDO
      DO NW=NW1,NW2
       DO JK=1,plonl
        do I = 1,N
        RR(JK,I,ND(JK),NW) = 0.0_dp
        do J = 1,N
           RR(JK,I,ND(JK),NW) = RR(JK,I,ND(JK),NW) + &
     	                        B(JK,I,J,NW)*H(JK,J,NW)
	enddo
        enddo
       ENDDO
      ENDDO 
!C-----------BACK SOLUTION
      DO NW=NW1,NW2
       DO JK=1,plonl
        do ID = ND(JK)-1,1,-1
         do I = 1,N
          do J = 1,N
           RR(JK,I,ID,NW) = RR(JK,I,ID,NW) + &
                           DD(JK,I,J,ID,NW)*RR(JK,J,ID+1,NW)
          enddo
         enddo
        enddo
       ENDDO
      ENDDO 	

!C----------MEAN J & H
      DO NW=NW1,NW2
       DO JK=1,plonl
        do ID = 1,ND(JK),2
         FJ(JK,ID,NW) = 0.0_dp
         do I = 1,N
          FJ(JK,ID,NW) = FJ(JK,ID,NW) + RR(JK,I,ID,NW)*WT(I)
         enddo
        enddo
        do ID = 2,ND(JK),2
         FJ(JK,ID,NW) = 0.0_dp
         do I = 1,N
          FJ(JK,ID,NW) = FJ(JK,ID,NW) + RR(JK,I,ID,NW)*WT(I)*EMU(I)
         enddo
        enddo
       ENDDO
      ENDDO 	

!c---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
!c---FJBOT = scaled diffuse flux onto surface:  
       FJTOP(:,:) = 0.0_dp
       FJBOT(:,:) = 0.0_dp
       do I = 1,N
        DO NW=NW1,NW2
         DO JK=1,plonl
          FJTOP(JK,NW) = FJTOP(JK,NW) + RR(JK,I,1,NW)*WT(I)*EMU(I)
          FJBOT(JK,NW) = FJBOT(JK,NW) + RR(JK,I,ND(JK),NW)*WT(I)*EMU(I)
         ENDDO
	ENDDO 
       enddo
       DO NW=NW1,NW2
        DO JK=1,plonl
         FJTOP(JK,NW) = 4._dp*FJTOP(JK,NW)
!c---         ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
         FIPLUS = 4._dp*ZREFL(JK)*FJBOT(JK,NW)/(1.0_dp + ZREFL(JK)) &
                  + ZFLUX(JK,NW)
         FJBOT(JK,NW) = 4._dp*FJBOT(JK,NW) - FIPLUS
	ENDDO 
       ENDDO

      return
      end subroutine blkslv_v
!C-----------------------------------------------------------------------
      subroutine JP_ATM(KLEV,PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,DTAUX,POMEGAX,JXTRA)
!c-----------------------------------------------------------------------
      use mo_kind, only : dp

      implicit none
!c-----------------------------------------------------------------------
!c--------key amtospheric data needed to solve plane-parallel J----------
      integer :: KLEV
      real(dp), dimension(KLEV+2) :: TTJ,DDJ,ZZJ,ZHL
      real(dp), dimension(KLEV+2) :: PPJ 
      integer,dimension(2*KLEV+3) :: JXTRA
      real(dp)                 :: ZZHT
      real(dp)                 :: DTAUX(KLEV+1),POMEGAX(KLEV+1,8)
!c-----------------------------------------------------------------------
      integer  I,J,K,L
      real(dp) ::  COLO2,COLO3,ZKM,DELZ,ZTOP

      write(2,'(4a)') '   L z(km)     p      T   ', &
      '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb',&
      '  g(cos) CTM lyr=>'

          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZHL(KLEV+1) + ZZHT

        do L = KLEV+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0  
          COLO3 = COLO3 + ZZJ(L)
          DELZ = ZTOP-ZHL(L)
          ZTOP = ZHL(L)
          ZKM = ZHL(L)*1.d-5

      write(2,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)')& 
           L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,ZZJ(L)/DELZ,&
           COLO2,COLO3,DTAUX(L),POMEGAX(L,1),POMEGAX(L,2)/3.d0,&
           JXTRA(L+L),JXTRA(L+L-1)


        enddo

      RETURN
  END SUBROUTINE jp_atm

!--------------------------------------------------------------------------------
END MODULE mo_moz_fastj

