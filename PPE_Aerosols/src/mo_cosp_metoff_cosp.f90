!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
! (c) COPYRIGHT British Crown / Met Office 2008  
! Please refer to Met_Office_license.txt for details. 
MODULE mo_cosp_metoff_cosp
  USE mo_cosp_metoff_cosp_simulator
  USE mo_cosp_constants
  USE mo_kind,           ONLY: wp
  USE mo_random_numbers, ONLY: random_uniform=>get_random

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP ----------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE COSP( kbdim, klev, Ncolumns,                       & 
                  Lisccp_sim,  Llidar_sim, Llidar_cfad,  Lstats,            &
!gbx
                  use_reff, lidar_ice_type,                    &
                  isccp_top_heigh, isccp_top_height_direction, &
                  isccp_overlap, isccp_emsfc_lw,               &
                  ph, p, T, Reff, tca,                         &             
                  psfc, mr_hydro, sh, dem_s, dtau_s,           &
                  sunlit, skt, zlev, zlev_half, land,          &
!sglidar
                  beta_mol, beta_tot, refl,  tau_tot,          &
!stlidar
                  srbval, cfad_sr, lidarcld,                   &
                  cldlayer, parasolrefl,                       &
!isccp
                  fq_isccp,  totalcldarea,                     &
                  meanptop,  meantaucld,                       &
                  meantb, meantbclr, boxtau,  boxptop,         &
                  meanalbedocld                                )
  ! Arguments


  LOGICAL, INTENT(IN) ::  Lisccp_sim
  LOGICAL, INTENT(IN) ::  Llidar_sim
  LOGICAL, INTENT(IN) ::  Llidar_cfad
  LOGICAL, INTENT(IN) ::  Lstats

!gbx
  LOGICAL, INTENT(IN) :: use_reff

  INTEGER, INTENT(IN) :: lidar_ice_type

  INTEGER, INTENT(IN) :: isccp_top_heigh
  INTEGER, INTENT(IN) :: isccp_top_height_direction
  INTEGER, INTENT(IN) :: isccp_overlap !  overlap type in SCOPS: 1=max, 2=rand, 3=max/rand

  INTEGER, INTENT(IN) :: kbdim         ! Number of gridpoints
  INTEGER, INTENT(IN) :: klev          ! Number of levels
  INTEGER, INTENT(IN) :: Ncolumns      ! Number of columns

  REAL(wp), INTENT(IN) :: isccp_emsfc_lw 

  ! check intent statements !
  REAL(wp), INTENT(INOUT) :: p(kbdim,klev) 
  REAL(wp), INTENT(INOUT) :: ph(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: T(kbdim,klev) 
  REAL(wp), INTENT(INOUT) :: Reff(kbdim,klev,N_hydro)
  REAL(wp), INTENT(INOUT) :: tca(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: psfc(kbdim)
  REAL(wp), INTENT(INOUT) :: mr_hydro(kbdim,klev,N_hydro)
  REAL(wp), INTENT(INOUT) :: sh(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: dem_s(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: dtau_s(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: sunlit(kbdim)
  REAL(wp), INTENT(INOUT) :: skt(kbdim)
  REAL(wp), INTENT(INOUT) :: zlev(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: zlev_half(kbdim,klev)
  REAL(wp), INTENT(INOUT) :: land(kbdim)

!sglidar
  REAL(wp) :: beta_mol(kbdim,klev)
  REAL(wp) :: beta_tot(kbdim,Ncolumns,klev)
  REAL(wp) :: refl(kbdim,Ncolumns,PARASOL_NREFL)
  REAL(wp) :: tau_tot(kbdim,Ncolumns,klev)

! isccp
  REAL(wp) :: fq_isccp(kbdim,7,7), totalcldarea(kbdim)
  REAL(wp) :: meanptop(kbdim), meantaucld(kbdim)
  REAL(wp) :: meantb(kbdim), meantbclr(kbdim)
  REAL(wp) :: boxtau(kbdim,Ncolumns), boxptop(kbdim,Ncolumns)
  REAL(wp) :: meanalbedocld(kbdim)
    
! stlidar lidar stats
  REAL(wp) :: srbval(SR_BINS)
  REAL(wp) :: cfad_sr(kbdim,SR_BINS,Nlr) 
  REAL(wp) :: lidarcld(kbdim,Nlr)
  REAL(wp) :: cldlayer(kbdim,LIDAR_NCAT)
  REAL(wp) :: parasolrefl(kbdim,PARASOL_NREFL)

! Local variables 

! sghydro
  REAL(wp) :: sg_mr_hydro(kbdim,Ncolumns,klev,N_hydro)

! subgrid sgx
  REAL(wp) :: frac_out(kbdim,Ncolumns,klev)

  INTEGER :: i,j,k

  LOGICAL :: reff_zero

  REAL(wp),DIMENSION(kbdim,klev) :: column_frac_out ! Array with one column of frac_out
  INTEGER :: scops_debug = 0    !  set to non-zero value to print out inputs for debugging in SCOPS

  INTEGER,DIMENSION(:),allocatable,save :: & ! Dimensions nPoints
                  seed    !  It is recommended that the seed is set to a different value for each model
                          !  gridbox it is called on, as it is possible that the choice of the same 
                          !  seed value every time may introduce some statistical bias in the results, 
                          !  particularly for low values of NCOL.
!$OMP THREADPRIVATE(seed)

  REAL(wp),DIMENSION(Kbdim,Klev) :: tca_scops ! Cloud cover in each model level 
                                                ! (HORIZONTAL gridbox fraction) of total cloud.
                                                ! Levels are from TOA to SURFACE. (kbdim, 0:nLev)
  REAL(wp),DIMENSION(Kbdim,Klev) :: frac_ls     ! Cloud fraction in each model level
                                                ! Levels are from SURFACE to TOA
  REAL(wp) :: maxp,minp
  
    frac_out  = 0.0_wp

  !++++++++++ Dimensions +++++++++++
  
  reff_zero=.true.
  IF (any(Reff > 1.e-8_wp)) THEN
     reff_zero=.false.
      ! reff_zero == .false.
      !     and gbx%use_reff == .true.   Reff use in radar and lidar
      !     and reff_zero    == .false.  Reff use in lidar and set to 0 for radar
  ENDIF
  IF ((.not. use_reff) .and. (reff_zero)) THEN ! No Reff in radar. Default in lidar
        Reff = DEFAULT_LIDAR_REFF
        print *, '---------- COSP WARNING ------------'
        print *, ''
        print *, 'Using default Reff in lidar simulations'
        print *, ''
        print *, '----------------------------------'
  ENDIF

  !++++++++++ Subgrid sampling ++++++++++
  ! Allocate arrays before calling SCOPS
  allocate(seed(Kbdim))
  ! Cloud fractions for SCOPS from TOA to SFC
  ! tca_scops(:,0) = 0.0_wp ! Zeros on top   CNam: Must keep separate from below (unlike v1.0 release)
  tca_scops(:,1:Klev) = tca(:,Klev:1:-1)

  frac_ls=0.0_wp
  
   ! We base the seed in the decimal part of the surface pressure.
  seed = int(psfc) ! This is to avoid division by zero when Kbdim = 1
  minp = minval(psfc)
  maxp = maxval(psfc)
  IF (Kbdim .gt. 1) seed=int((psfc-minp)/(maxp-minp)*1000000._wp)+1
  

  ! SCOPS is called
  ! strat and conv arrays are passed with levels from TOA to SURFACE.  
  !call scops(Kbdim,Klev,Ncolumns,seed,tca_scops,isccp_overlap,frac_out,scops_debug)
  call scops(Kbdim,Klev,Ncolumns,tca_scops,isccp_overlap,frac_out,scops_debug)
  
  
  DO j=1,Kbdim,1
   DO k=1,Klev,1
    DO i=1,Ncolumns,1
     IF (frac_out (j,i,Klev+1-k) .eq. 1._wp) frac_ls(j,k)=frac_ls(j,k)+1._wp
    ENDDO  !i
    frac_ls(j,k)=frac_ls(j,k)/REAL(Ncolumns,wp)
   ENDDO  !k
  ENDDO  !j
   
  ! Levels from SURFACE to TOA for output.
  ! This can be done within a loop (unvectorized) over kbdim to save memory,
  ! because frac_out in vector mode can be very big
  DO j=1,Kbdim
   frac_out(j,:,1:Klev)  = frac_out(j,:,Klev:1:-1)
  ENDDO

  ! Deallocate arrays that will no longer be used
  deallocate(seed)

  
  !++++++++++ Simulator ++++++++++
    sg_mr_hydro = 0.0_wp
   

  ! Populate the subgrid arrays
  DO k=1,Ncolumns
      !--------- Mixing ratios for clouds and Reff for Clouds -------
      column_frac_out(:,:) = frac_out(:,k,:)
      where (column_frac_out == 1._wp)     !+++++++++++ LS clouds ++++++++
          sg_mr_hydro(:,k,:,I_LSCLIQ) = mr_hydro(:,:,I_LSCLIQ)
          sg_mr_hydro(:,k,:,I_LSCICE) = mr_hydro(:,:,I_LSCICE)
      endwhere 
  ENDDO

  ! convert the mixing ratio  from gridbox mean to the fraction-based values
   DO j=1,kbdim
    DO k=1,Klev
     !--------- Clouds -------
     IF (frac_ls(j,k) .ne. 0._wp) THEN
     sg_mr_hydro(j,:,k,I_LSCLIQ) = sg_mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
     sg_mr_hydro(j,:,k,I_LSCICE) = sg_mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
     ENDIF
    ENDDO !k
   ENDDO !j
  

  !++++++++++ Simulator ++++++++++
  call cosp_simulator(kbdim,Ncolumns,klev,      &
       Lisccp_sim,  Llidar_sim, Llidar_cfad, Lstats,        &
!gbx 
       lidar_ice_type,                              &
       isccp_top_heigh, isccp_top_height_direction, &
       isccp_overlap, isccp_emsfc_lw,               &
       ph, p, T, Reff, tca,                         &             
       sh, dem_s, dtau_s,                           &
       sunlit, skt, zlev, zlev_half, land,          &
!sgx
      frac_out,                                     &
      sg_mr_hydro,                                  & 
!sglidar
       beta_mol, beta_tot, refl, tau_tot,          &
!stlidar
      srbval, cfad_sr, lidarcld,                   &
      cldlayer, parasolrefl,                       &
!isccp
      fq_isccp,  totalcldarea,                     &
      meanptop,  meantaucld,                       &              
      meantb, meantbclr, boxtau,  boxptop, &
      meanalbedocld )

 
END SUBROUTINE COSP


!CNam: *** INSERTED BECAUSE WOULD NOT READ FROM /SRC/

!SUBROUTINE scops(npoints,nlev,ncol,seed,cc,  &
SUBROUTINE scops(npoints,nlev,ncol,cc,  &
                             overlap,frac_out,ncolprint)

! *****************************COPYRIGHT*******************************
! (c) COPYRIGHT Steve Klein and Mark Webb 2004, All Rights Reserved.
! Steve Klein klein21@mail.llnl.gov
! Mark Webb mark.webb@metoffice.gov.uk markwebb@mail.com
! ISCCP SIMULATOR icarus-scops version 3.5

! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.

! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
! See also http://www.gnu.org/copyleft/lesser.txt
!
! The Met Office hereby disclaims all copyright interest in
! the ISCCP Simulator ( a library which converts model clouds into
! direct equivalents to the satellite observations of the ISCCP)
! written by Steve Klein and Mark Webb
!
! Catherine Senior, August 2004
! Manager Climate Sensitivy Group
! Met Office   Hadley Centre for Climate Prediction and Research
! FitzRoy Road  Exeter  EX1 3PB  United Kingdom
! *****************************COPYRIGHT*******************************

      USE mo_kind,   ONLY: wp

      IMPLICIT NONE

      INTEGER :: npoints       !  number of model points in the horizontal
      INTEGER :: nlev          !  number of model levels in column
      INTEGER :: ncol          !  number of subcolumns


      INTEGER :: overlap       !  overlap type
                               !  1=max
                               !  2=rand
                               !  3=max/rand

      REAL(wp) :: cc(npoints,nlev)
                  !  input cloud cover in each model level (fraction)
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by clouds


      INTEGER :: j,ilev,ibox,ncolprint,ilev2

      REAL(wp) ::  frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column


      !INTEGER :: seed(npoints)
      !  seed values for marsaglia  random number generator
      !  It is recommended that the seed is set
      !  to a different value for each model
      !  gridbox it is called on, as it is
      !  possible that the choice of the same
      !  seed value every time may introduce some
      !  statistical bias in the results, particularly
      !  for low values of NCOL.

      REAL(wp) ::  tca(npoints,0:nlev)  ! total cloud cover in each model level (fraction)
                                        ! with extra layer of zeroes on top
                                        ! in this version this just contains the values input
                                        ! from cc but with an extra level

      REAL(wp) ::  threshold(npoints,ncol)     ! pointer to position in gridbox
      REAL(wp) ::  maxosc(npoints,ncol)        ! Flag for max overlapped strat cld

      REAL(wp) ::  boxpos(npoints,ncol)        ! ordered pointer to position in gridbox

      REAL(wp) ::  threshold_min(npoints,ncol) ! minimum value to define range in with new threshold
                                               ! is chosen

      REAL(wp) :: ran(npoints)                 ! vector of random numbers

!!$      INTEGER  :: irand,i2_16,overflow_32  ! variables for RNG
!!$      INTEGER, PARAMETER :: huge32=2147483647

 !!$      i2_16=65536

      DO ibox=1,ncol
        DO j=1,npoints 
        boxpos(j,ibox)=(REAL(ibox,wp)-.5_wp)/REAL(ncol,wp)
        ENDDO
      ENDDO

!     ---------------------------------------------------!
!     Initialise working variables
!     ---------------------------------------------------!

!     Initialised frac_out to zero

      DO ilev=1,nlev
        DO ibox=1,ncol
          DO j=1,npoints
          frac_out(j,ibox,ilev)=0.0_wp
          ENDDO
        ENDDO
      ENDDO

!     assign 2d tca array using 1d input array cc

      DO j=1,npoints
        tca(j,0)=0._wp
      ENDDO

      DO ilev=1,nlev
        DO j=1,npoints
          tca(j,ilev)=cc(j,ilev)
        ENDDO
      ENDDO

      IF (ncolprint.ne.0) THEN
        write (6,'(a)') 'frac_out_pp_rev:'
          DO j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(8f5.2)') &
           ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)

          ENDDO
        write (6,'(a)') 'ncol:'
        write (6,'(I3)') ncol
      ENDIF
      IF (ncolprint.ne.0) THEN
        write (6,'(a)') 'last_frac_pp:'
          DO j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(8f5.2)') (tca(j,0))
          ENDDO
      ENDIF

!     ---------------------------------------------------!
!     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, KLEV
!     frac_out is the array that contains the information 
!     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
!     convective cloud
      
      !loop over vertical levels
      DO 200 ilev = 1,nlev
                                  
!     Initialise threshold

        IF (ilev.eq.1) THEN
          ! If max overlap 
          IF (overlap.eq.1) THEN
            ! select pixels spread evenly
            ! across the gridbox
              DO ibox=1,ncol
                DO j=1,npoints
                  threshold(j,ibox)=boxpos(j,ibox)
                ENDDO
              ENDDO
          ELSE
              DO ibox=1,ncol
!**************************************************************8
!***            include 'congvec.f'
! --- Generate Random Number: Originally in congvec.f
!!$      do irand = 1, npoints
!!$          ! Marsaglia CONG algorithm
!!$          seed(irand) = 69069*seed(irand)+1234567
!!$          ! mod 32 bit overflow
!!$          seed(irand) = mod(seed(irand),2**30)   
!!$          ran(irand) = seed(irand)*0.931322574615479E-09
!!$      end do
!!$
!!$
!!$      ! convert to range 0-1 (32 bit only)
!!$      overflow_32=i2_16*i2_16
!!$      if ( overflow_32 .le. huge32 ) then
!!$          do irand = 1, npoints
!!$              ran(irand)=ran(irand)+1
!!$              ran(irand)=(ran(irand))-int(ran(irand))
!!$         end do
!!$      end if

  CALL random_uniform(ran)

! ---
!***********************************************8


                ! select random pixels from the non-convective
                ! part the gridbox ( some will be converted into
                ! convective pixels below )
                DO j=1,npoints
                  threshold(j,ibox)= ran(j)
                ENDDO
              ENDDO
            ENDIF
            IF (ncolprint.ne.0) THEN
              write (6,'(a)') 'threshold_nsf2:'
                DO j=1,npoints,1000
                write(6,'(a10)') 'j='
                write(6,'(8I10)') j
                write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
                ENDDO
            ENDIF
        ENDIF

        IF (ncolprint.ne.0) THEN
            write (6,'(a)') 'ilev:'
            write (6,'(I2)') ilev
        ENDIF

        DO ibox=1,ncol


          ! Max overlap
          IF (overlap.eq.1) THEN 
            DO j=1,npoints
              threshold_min(j,ibox)=0._wp
              maxosc(j,ibox)=1._wp
            ENDDO
          ENDIF

          ! Random overlap
          IF (overlap.eq.2) THEN 
            DO j=1,npoints
              threshold_min(j,ibox)=0._wp
              maxosc(j,ibox)=0._wp
            ENDDO
          ENDIF

          ! Max/Random overlap
          IF (overlap.eq.3) THEN 
            DO j=1,npoints
              threshold_min(j,ibox)=max(0._wp, &
                min(tca(j,ilev-1),tca(j,ilev)))
              IF (threshold(j,ibox) &
               .lt.min(tca(j,ilev-1),tca(j,ilev)) &
               .and.(threshold(j,ibox).gt.0._wp)) THEN
                   maxosc(j,ibox)= 1._wp
              else
                   maxosc(j,ibox)= 0._wp
              ENDIF
            ENDDO
          ENDIF
    
          ! Reset threshold 
!********************************************************

! ***       include 'congvec.f'
! --- Generate Random Number: Originally in congvec.f
!!$      do irand = 1, npoints
!!$          ! Marsaglia CONG algorithm
!!$          seed(irand) = 69069*seed(irand)+1234567
!!$          ! mod 32 bit overflow
!!$          seed(irand) = mod(seed(irand),2**30)   
!!$          ran(irand) = seed(irand)*0.931322574615479E-09_wp
!!$      end do
!!$
!!$
!!$      ! convert to range 0-1 (32 bit only)
!!$      overflow_32=i2_16*i2_16
!!$      if ( overflow_32 .le. huge32 ) then
!!$          do irand = 1, npoints
!!$              ran(irand)=ran(irand)+1
!!$              ran(irand)=(ran(irand))-int(ran(irand))
!!$         end do
!!$      end if

  CALL random_uniform(ran)
! ---
!********************************************************
          DO j=1,npoints
            threshold(j,ibox)= &
                ( &                                   
                  !if max overlapped strat cloud
                  (maxosc(j,ibox)) * ( &                                
                      !threshold=boxpos
                      threshold(j,ibox) &
                  ) + &                                                 
                  !else
                  (1._wp-maxosc(j,ibox)) * ( &                              
                      !threshold_min=random[thrmin,1]
                      threshold_min(j,ibox)+ &
                        (1._wp-threshold_min(j,ibox))*ran(j) &  
                 ) &
              )
          ENDDO

        ENDDO ! ibox

!          Fill frac_out with 1's where tca is greater than the threshold

           DO ibox=1,ncol
             DO j=1,npoints 
               IF (tca(j,ilev).gt.threshold(j,ibox)) THEN
               frac_out(j,ibox,ilev)=1._wp
               else
               frac_out(j,ibox,ilev)=0._wp
               ENDIF               
             ENDDO
           ENDDO



!         Set last_frac to tca at this level, so as to be tca 
!         from last level next time round

          IF (ncolprint.ne.0) THEN

            DO j=1,npoints ,1000
            write(6,'(a10)') 'j='
            write(6,'(8I10)') j
            write (6,'(a)') 'last_frac:'
            write (6,'(8f5.2)') (tca(j,ilev-1))
    
            write (6,'(a)') 'max_overlap_sc:'
            write (6,'(8f5.2)') (maxosc(j,ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'threshold_min_nsf2:'
            write (6,'(8f5.2)') (threshold_min(j,ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'threshold_nsf2:'
            write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
    
            write (6,'(a)') 'frac_out_pp_rev:'
            write (6,'(8f5.2)') &
             ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
          ENDDO
          ENDIF

200   CONTINUE    !loop over nlev

      END SUBROUTINE scops
!***

END MODULE mo_cosp_metoff_cosp
