!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!------------------------------------------------------------------------------------
! Authors: Sandrine Bony and Helene Chepfer (LMD/IPSL, CNRS, UPMC, France).
!------------------------------------------------------------------------------------
MODULE mo_cosp_lmd_ipsl_stats

  USE mo_cosp_llnl_stats
  USE mo_kind,        ONLY: wp
  USE mo_memory_g3b,  ONLY: slm

  IMPLICIT NONE

CONTAINS

      SUBROUTINE diag_lidar(npoints,ncol,llm,max_bin,nrefl & ! parasol
                  ,pnorm,pmol,refl,land,pplay,undef,ok_lidar_cfad & !parasol,
                  ,cfad2,srbval,ncat,lidarcld,cldlayer,parasolrefl) ! parasol

! c-----------------------------------------------------------------------------------
! Lidar outputs :
! 
! Diagnose cloud fraction (3D cloud fraction + low/middle/high/total cloud fraction
! from the lidar signals (ATB and molecular ATB) computed from model outputs
!      +
! Compute CFADs of lidar scattering ratio SR and of depolarization index
! 
! Authors: Sandrine Bony and Helene Chepfer (LMD/IPSL, CNRS, UPMC, France).
!
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - change of the cloud detection threshold S_cld from 3 to 5, for better
! with both day and night observations. The optical thinest clouds are missed.
! - remove of the detection of the first fully attenuated layer encountered from above.
! December 2008, A. Bodas-Salcedo:
! - Dimensions of pmol reduced to (npoints,llm)
!
! Version 1.0 (June 2007)
! Version 1.1 (May 2008)
! Version 1.2 (June 2008)
! Version 2.0 (October 2008)
! Version 2.1 (December 2008)
! c------------------------------------------------------------------------------------
! c inputs :
      INTEGER :: npoints
      INTEGER :: ncol
      INTEGER :: llm
      INTEGER :: max_bin               ! nb of bins for SR CFADs
      INTEGER :: ncat                  ! nb of cloud layer types (low,mid,high,total)
      INTEGER :: nrefl                 ! nb of solar zenith angles for parasol reflectances ! parasol

      REAL(wp) :: undef                    ! undefined value
      REAL(wp) :: pnorm(npoints,ncol,llm)  ! lidar ATB 
      REAL(wp) :: pmol(npoints,llm)   ! molecular ATB
      REAL(wp) :: land(npoints)            ! Land-Sea mask [0:Ocean 1:Land]
      REAL(wp) :: pplay(npoints,llm)       ! pressure on model levels (Pa)
      LOGICAL ::  ok_lidar_cfad         ! true if lidar CFAD diagnostics need to be computed
      REAL(wp) :: refl(npoints,ncol,nrefl) ! subgrid parasol reflectance ! parasol

! c outputs :
      REAL(wp) :: lidarcld(npoints,llm)     ! 3D "lidar" cloud fraction 
      REAL(wp) :: cldlayer(npoints,ncat)    ! "lidar" cloud fraction (low, mid, high, total)
      REAL(wp) :: cfad2(npoints,max_bin,llm) ! CFADs of SR  
      REAL(wp) :: srbval(max_bin)           ! SR bins in CFADs  
      REAL(wp) :: parasolrefl(npoints,nrefl)! grid-averaged parasol reflectance ! parasol


! c threshold for cloud detection :
      REAL(wp), PARAMETER :: S_clr = 1.2_wp 
      REAL(wp), PARAMETER :: S_cld = 5.0_wp !Threshold for cloud detection (Dec 2008)
!      parameter (S_cld = 3.0_wp)
      REAL(wp), PARAMETER :: S_att = 0.01_wp

! c local variables :
      INTEGER :: ic,k
      REAL(wp) :: x3d(npoints,ncol,llm)
      REAL(wp) :: x3d_c(npoints,llm),pnorm_c(npoints,llm)
      REAL(wp) :: xmax
!
! c -------------------------------------------------------
! c 0- Initializations
! c -------------------------------------------------------
!
! Parasol reflectance algorithm is not valid over land. Write
! a warning if there is no land. Landmask [0 - Ocean, 1 - Land] 
!CNam:      IF ( MAXVAL(land(:)) .EQ. 0.0_wp) THEN
!      IF ( MAXVAL(land(:)) .GT. 0.0_wp) THEN
!          WRITE (*,*) 'WARNING. PARASOL reflectance is not valid over land' &
!            & ,' and there is only land'
!      END IF

      xmax=undef-1.0_wp

! c -------------------------------------------------------
! c 1- Lidar scattering ratio :
! c -------------------------------------------------------
!
!       where ((pnorm.lt.xmax) .and. (pmol.lt.xmax) .and. (pmol.gt. 0.0 )) 
!          x3d = pnorm/pmol
!       elsewhere
!           x3d = undef
!       end where
! A.B-S: pmol reduced to 2D (npoints,llm) (Dec 08)
      DO ic = 1, ncol
        pnorm_c = pnorm(:,ic,:)
        WHERE ((pnorm_c.lt.xmax) .and. (pmol.lt.xmax) .and. (pmol.gt. 0.0_wp )) 
            x3d_c = pnorm_c/pmol
        ELSEWHERE
            x3d_c = undef
        ENDWHERE
        x3d(:,ic,:) = x3d_c
      ENDDO

! c -------------------------------------------------------
! c 2- Diagnose cloud fractions (3D, low, middle, high, total)
! c from subgrid-scale lidar scattering ratios :
! c -------------------------------------------------------

      CALL cosp_cldfrac(npoints,ncol,llm,ncat,  &
              x3d,pplay,S_att,S_cld,undef,lidarcld, &
              cldlayer)

! c -------------------------------------------------------
! c 3- CFADs 
! c -------------------------------------------------------
      IF (ok_lidar_cfad) THEN
!
! c CFADs of subgrid-scale lidar scattering ratios :
! c -------------------------------------------------------
      CALL cosp_cfad_sr(npoints,ncol,llm,max_bin,undef, &
                 x3d,S_att,S_clr,xmax,cfad2,srbval)

      ENDIF   ! ok_lidar_cfad
! c -------------------------------------------------------
! c -------------------------------------------------------
! c 4- Compute grid-box averaged Parasol reflectances
! c -------------------------------------------------------
      parasolrefl(:,:) = 0.0_wp

      DO k = 1, nrefl
       DO ic = 1, ncol
         parasolrefl(:,k) = parasolrefl(:,k) + refl(:,ic,k)
       ENDDO
      ENDDO

      DO k = 1, nrefl
        parasolrefl(:,k) = parasolrefl(:,k) / REAL(ncol,wp)
! if land=1 -> parasolrefl=undef
! if land=0 -> parasolrefl=parasolrefl
        parasolrefl(:,k) = parasolrefl(:,k) * MAX(1.0_wp-land(:),0.0_wp) &
                           + (1.0_wp - MAX(1.0_wp-land(:),0.0_wp))*undef 
      ENDDO

      END SUBROUTINE diag_lidar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- FUNCTION COSP_CFAD_SR ------------------------
! Author: Sandrine Bony (LMD/IPSL, CNRS, Paris)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cosp_cfad_sr(Npoints,Ncolumns,Nlevels,Nbins,undef, &
                      x,S_att,S_clr,xmax,cfad,srbval)

      IMPLICIT NONE

!--- Input arguments
! Npoints: Number of horizontal points
! Ncolumns: Number of subcolumns
! Nlevels: Number of levels
! Nbins: Number of x axis bins
! xmax: maximum value allowed for x
! S_att: Threshold for full attenuation
! S_clr: Threshold for clear-sky layer
!
!--- Input-Outout arguments
! x: variable to process (Npoints,Ncolumns,Nlevels), modified where saturation occurs
!
! -- Output arguments
! srbval : values of the histogram bins
! cfad: 2D histogram on each horizontal point

! Input arguments
      INTEGER :: Npoints,Ncolumns,Nlevels,Nbins
      REAL(wp) :: xmax,S_att,S_clr,undef 
! Input-outout arguments
      REAL(wp) :: x(Npoints,Ncolumns,Nlevels)
! Output :
      REAL(wp) :: cfad(Npoints,Nbins,Nlevels)
      REAL(wp) :: srbval(Nbins)
! Local variables
      INTEGER :: i, j, k, ib
      REAL(wp) :: srbval_ext(0:Nbins)

! c -------------------------------------------------------
! c 0- Initializations
! c -------------------------------------------------------
      IF ( Nbins .lt. 6) return

      srbval(1) =  S_att
      srbval(2) =  S_clr
      srbval(3) =  3.0_wp
      srbval(4) =  5.0_wp
      srbval(5) =  7.0_wp
      srbval(6) = 10.0_wp
      DO i = 7, MIN(10,Nbins)
       srbval(i) = srbval(i-1) + 5.0_wp
      ENDDO
      DO i = 11, MIN(13,Nbins)
       srbval(i) = srbval(i-1) + 10.0_wp
      ENDDO
      srbval(MIN(14,Nbins)) = 80.0_wp
      srbval(Nbins) = xmax
      cfad(:,:,:) = 0.0_wp

! c -------------------------------------------------------
! c c- Compute CFAD
! c -------------------------------------------------------

!cms++ replace from v1.3:
      srbval_ext(1:Nbins) = srbval
      srbval_ext(0) = -1.0_wp

      DO j = 1, Nlevels
         DO ib = 1, Nbins
            DO k = 1, Ncolumns
               DO i = 1, Npoints
                  IF (x(i,k,j) /= undef) THEN
                     IF ((x(i,k,j).gt.srbval_ext(ib-1)).and.(x(i,k,j).le.srbval_ext(ib))) &
                          cfad(i,ib,j) = cfad(i,ib,j) + 1.0_wp
                  ELSE 
                     cfad(i,ib,j) = undef
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      WHERE (cfad .ne. undef)  cfad = cfad / REAL(Ncolumns,wp)

!!$
!!$        do j = Nlevels, 1, -1 
!!$          do k = 1, Ncolumns
!!$              where ( x(:,k,j).le.srbval(1) ) &
!!$                        cfad(:,1,j) = cfad(:,1,j) + 1.0_wp
!!$          enddo  !k
!!$        enddo  !j
!!$
!!$      do ib = 2, Nbins
!!$        do j = Nlevels, 1, -1 
!!$          do k = 1, Ncolumns
!!$              where ( ( x(:,k,j).gt.srbval(ib-1) ) .and. ( x(:,k,j).le.srbval(ib) ) ) &
!!$                        cfad(:,ib,j) = cfad(:,ib,j) + 1.0_wp
!!$          enddo  !k
!!$        enddo  !j
!!$      enddo  !ib
!!$
!!$      cfad(:,:,:) = cfad(:,:,:) / REAL(Ncolumns,wp)
!cms --
! c -------------------------------------------------------

      END SUBROUTINE cosp_cfad_sr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- SUBROUTINE COSP_CLDFRAC -------------------
! c Purpose: Cloud fraction diagnosed from lidar measurements 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cosp_cldfrac(Npoints,Ncolumns,Nlevels,Ncat, &
                  x,pplay,S_att,S_cld,undef,lidarcld, &
                  cldlayer)
      IMPLICIT NONE
! Input arguments
      INTEGER :: Npoints,Ncolumns,Nlevels,Ncat
      REAL(wp) :: x(Npoints,Ncolumns,Nlevels)
      REAL(wp) :: pplay(Npoints,Nlevels)
      REAL(wp) :: S_att,S_cld
      REAL(wp) :: undef
! Output :
      REAL(wp) :: lidarcld(Npoints,Nlevels) ! 3D cloud fraction
      REAL(wp) :: cldlayer(Npoints,Ncat)    ! low, middle, high, total cloud fractions
! Local variables
      INTEGER :: ip, k, iz, ic
      REAL(wp) :: p1
      REAL(wp) :: cldy(Npoints,Ncolumns,Nlevels)
      REAL(wp) :: srok(Npoints,Ncolumns,Nlevels)
      REAL(wp) :: cldlay(Npoints,Ncolumns,Ncat)
      REAL(wp) :: nsublay(Npoints,Ncolumns,Ncat), nsublayer(Npoints,Ncat)
      REAL(wp) :: nsub(Npoints,Nlevels)

! ---------------------------------------------------------------
! 1- initialization 
! ---------------------------------------------------------------

      IF ( Ncat .ne. 4 ) THEN
         print *,'Error in lmd_ipsl_stats.cosp_cldfrac, Ncat must be 4, not',Ncat
         stop
      ENDIF

      lidarcld = 0.0_wp
      nsub = 0.0_wp
      cldlay = 0.0_wp
      nsublay = 0.0_wp

! ---------------------------------------------------------------
! 2- Cloud detection
! ---------------------------------------------------------------
      DO k = 1, Nlevels

! cloud detection at subgrid-scale:
         WHERE ( (x(:,:,k).gt.S_cld) .and. (x(:,:,k).ne. undef) )
           cldy(:,:,k)=1.0_wp
         ELSEWHERE
           cldy(:,:,k)=0.0_wp
         ENDWHERE

! number of usefull sub-columns:
         WHERE ( (x(:,:,k).gt.S_att) .and. (x(:,:,k).ne. undef)  ) 
           srok(:,:,k)=1.0_wp
         ELSEWHERE
           srok(:,:,k)=0.0_wp
         ENDWHERE

      ENDDO ! k

! ---------------------------------------------------------------
! 3- grid-box 3D cloud fraction and layered cloud fractions (ISCCP pressure
! categories) :
! ---------------------------------------------------------------

      DO k = Nlevels, 1, -1
       DO ic = 1, Ncolumns
        DO ip = 1, Npoints

          iz=1
          p1 = pplay(ip,k)
!jq          if ( p1.gt.0. .and. p1.lt.(440.*100.)) then ! high clouds          
          IF ( p1.ge.0._wp .and. p1.lt.(440._wp*100._wp)) THEN ! high clouds
            iz=3
          ELSEIF (p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) THEN  ! mid clouds
            iz=2
         ENDIF

         cldlay(ip,ic,iz) = MAX(cldlay(ip,ic,iz),cldy(ip,ic,k))
         cldlay(ip,ic,4) = MAX(cldlay(ip,ic,4),cldy(ip,ic,k))
         lidarcld(ip,k)=lidarcld(ip,k) + cldy(ip,ic,k)

         nsublay(ip,ic,iz) = MAX(nsublay(ip,ic,iz),srok(ip,ic,k))
         nsublay(ip,ic,4) = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
         nsub(ip,k)=nsub(ip,k) + srok(ip,ic,k)

        ENDDO
       ENDDO
      ENDDO

! -- grid-box 3D cloud fraction

      WHERE ( nsub(:,:).gt.0.0_wp )
         lidarcld(:,:) = lidarcld(:,:)/nsub(:,:)
      ELSEWHERE
         lidarcld(:,:) = undef
      ENDWHERE

! -- layered cloud fractions

      cldlayer = 0.0_wp
      nsublayer = 0.0_wp

      DO iz = 1, Ncat
       DO ic = 1, Ncolumns
          cldlayer(:,iz)=cldlayer(:,iz) + cldlay(:,ic,iz)    
          nsublayer(:,iz)=nsublayer(:,iz) + nsublay(:,ic,iz) 
       ENDDO
      ENDDO

      WHERE ( nsublayer(:,:).gt.0.0_wp )
         cldlayer(:,:) = cldlayer(:,:)/nsublayer(:,:)
      ELSEWHERE
         cldlayer(:,:) = undef
      ENDWHERE

      END SUBROUTINE cosp_cldfrac
! ---------------------------------------------------------------
	  
END MODULE mo_cosp_lmd_ipsl_stats
