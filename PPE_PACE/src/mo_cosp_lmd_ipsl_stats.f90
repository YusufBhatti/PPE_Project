!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
! Copyright (c) 2009, Centre National de la Recherche Scientifique
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 15:08:38 +0100 (Mi, 13. Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/actsim/lmd_ipsl_stats.F90 $
!
! Redistribution and use in source and binary forms, with or without modification, are permitted
! provided that the following conditions are met:
!
!     * Redistributions of source code must retain the above copyright notice, this list
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials
!       provided with the distribution.
!     * Neither the name of the LMD/IPSL/CNRS/UPMC nor the names of its
!       contributors may be used to endorse or promote products derived from this software without
!       specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!------------------------------------------------------------------------------------
! Authors: Sandrine Bony and Helene Chepfer (LMD/IPSL, CNRS, UPMC, France).
!------------------------------------------------------------------------------------
MODULE mo_cosp_lmd_ipsl_stats

  USE mo_cosp_llnl_stats
  USE mo_kind,        ONLY: wp
  USE mo_memory_g3b,  ONLY: slm

  IMPLICIT NONE

CONTAINS

      SUBROUTINE diag_lidar(npoints,ncol,llm,max_bin,nrefl &
                  ,tmp,pnorm,pnorm_perp,pmol,refl,land,pplay,undef,ok_lidar_cfad &
                  ,cfad2,srbval,ncat,lidarcld,lidarcldphase,cldlayer,cldlayerphase &
                  ,lidarcldtmp,parasolrefl)

! -----------------------------------------------------------------------------------
! Lidar outputs :
!
! Diagnose cloud fraction (3D cloud fraction + low/middle/high/total cloud fraction)
! and phase cloud fraction (3D, low/mid/high/total and 3D temperature)
! from the lidar signals (ATB, ATBperp and molecular ATB) computed from model outputs
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
! August 2009, A. Bodas-Salcedo:
! - Warning message regarding PARASOL being valid only over ocean deleted.
! February 2010, A. Bodas-Salcedo:
! - Undef passed into cosp_cfad_sr
! June 2010, T. Yokohata, T. Nishimura and K. Ogochi
! Optimisation of COSP_CFAD_SR
!
! January 2013, G. Cesana, H. Chepfer:
! - Add the perpendicular component of the backscattered signal (pnorm_perp) in the arguments
! - Add the temperature (tmp) in the arguments
! - Add the 3D Phase cloud fraction (lidarcldphase) in the arguments
! - Add the Phase low mid high cloud fraction (cldlayerphase) in the arguments
! - Add the 3D Phase cloud fraction as a function of temperature (lidarcldtmp) in the arguments
! - Modification of the phase diagnosis within the COSP_CLDFRAC routine to integrate the phase
!   diagnosis (3D, low/mid/high, 3D temperature)
! Reference: Cesana G. and H. Chepfer (2013): Evaluation of the cloud water phase
! in a climate model using CALIPSO-GOCCP, J. Geophys. Res., doi: 10.1002/jgrd.50376
!
! ------------------------------------------------------------------------------------

! c inputs :
      INTEGER :: npoints
      INTEGER :: ncol
      INTEGER :: llm   !CNam: note orig llm was Nlr vgrid 40, now Nlevels
      INTEGER :: max_bin               ! nb of bins for SR CFADs
      INTEGER :: ncat                  ! nb of cloud layer types (low,mid,high,total)
      INTEGER :: nrefl                 ! nb of solar zenith angles for parasol reflectances ! parasol

      REAL(wp) :: undef                    ! undefined value
      REAL(wp) :: pnorm(npoints,ncol,llm)  ! lidar ATB 
      REAL(wp) :: pmol(npoints,llm)        ! molecular ATB
      REAL(wp) :: land(npoints)            ! Land-Sea mask [0:Ocean 1:Land]
      REAL(wp) :: pplay(npoints,llm)       ! pressure on model levels (Pa)
      LOGICAL ::  ok_lidar_cfad            ! true if lidar CFAD diagnostics need to be computed
      REAL(wp) :: refl(npoints,ncol,nrefl) ! subgrid parasol reflectance ! parasol
      REAL(wp) :: tmp(npoints,llm)         ! temp at each levels
      REAL(wp) :: pnorm_perp(npoints,ncol,llm)  ! lidar perpendicular ATB

! c outputs :
      REAL(wp) :: lidarcld(npoints,llm)     ! 3D "lidar" cloud fraction 
      REAL(wp) :: sub(npoints,llm)          ! 3D "lidar" indice
      REAL(wp) :: cldlayer(npoints,ncat)    ! "lidar" cloud_layer fraction (low, mid, high, total)

      REAL(wp) :: cfad2(npoints,max_bin,llm)! CFADs of SR  
      REAL(wp) :: srbval(max_bin)           ! SR bins in CFADs  
      REAL(wp) :: parasolrefl(npoints,nrefl)! grid-averaged parasol reflectance ! parasol


! c threshold for cloud detection :
      REAL(wp), PARAMETER :: S_clr = 1.2_wp 
      REAL(wp), PARAMETER :: S_cld = 5.0_wp !Threshold for cloud detection (Dec 2008)
      REAL(wp), PARAMETER :: S_att = 0.01_wp

! c local variables :
      INTEGER :: ic,k,i,j
      REAL(wp) :: x3d(npoints,ncol,llm)
      REAL(wp) :: x3d_c(npoints,llm),pnorm_c(npoints,llm)
      REAL(wp) :: xmax

! Output variables
      INTEGER,PARAMETER :: nphase = 6 ! nb of cloud layer phase types (ice,liquid,undefined,false ice,false liquid,Percent of ice)
      REAL(wp) :: lidarcldphase(npoints,llm,nphase)   ! 3D "lidar" phase cloud fraction
      REAL(wp) :: lidarcldtmp(npoints,llm,5)          ! 3D "lidar" phase cloud fraction as a function of temp
      REAL(wp) :: cldlayerphase(npoints,ncat,nphase)  ! "lidar" phase low mid high cloud fraction 

! SR detection threshold
      REAL(wp),PARAMETER  ::  S_cld_att = 30._wp ! New threshold for undefine cloud phase detection	



! c -------------------------------------------------------
! c 0- Initializations
! c -------------------------------------------------------
!
      xmax=undef-1.0_wp

! c -------------------------------------------------------
! c 1- Lidar scattering ratio :
! c -------------------------------------------------------
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

      CALL cosp_cldfrac(npoints,ncol,llm,ncat,nphase, &
              tmp,x3d,pnorm,pnorm_perp,pplay,S_att,S_cld,S_cld_att,undef,lidarcld, &
              cldlayer,lidarcldphase,sub,cldlayerphase,lidarcldtmp)

! c -------------------------------------------------------
! c 3- CFADs 
! c -------------------------------------------------------
      IF (ok_lidar_cfad) THEN

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
        parasolrefl(:,k) = parasolrefl(:,k) / REAL(ncol,wp) !cms
! if land=1 -> parasolrefl=undef
! if land=0 -> parasolrefl=parasolrefl
        parasolrefl(:,k) = parasolrefl(:,k) * MAX(1.0_wp-land(:),0.0_wp) &
                           + (1.0_wp - MAX(1.0_wp-land(:),0.0_wp))*undef 
      ENDDO

      RETURN
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
!cms--
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

      WHERE (cfad .ne. undef)  cfad = cfad / REAL(Ncolumns,wp)  !cms

! c -------------------------------------------------------
      RETURN
      END SUBROUTINE cosp_cfad_sr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- SUBROUTINE COSP_CLDFRAC -------------------
! c Purpose: Cloud fraction diagnosed from lidar measurements 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cosp_cldfrac( Npoints,Ncolumns,Nlevels,Ncat,NPhase, &
                  tmp,x,ATB,ATBperp,pplay,S_att,S_cld,S_cld_att,undef,lidarcld, &
                  cldlayer,lidarcldphase,nsub,cldlayerphase,lidarcldtemp )

      IMPLICIT NONE
! Input arguments
      INTEGER :: Npoints, Ncolumns, Nlevels, Ncat
      REAL(wp) :: x(Npoints,Ncolumns,Nlevels)
      REAL(wp) :: tmp(Npoints,Nlevels)		! temperature
      REAL(wp) :: ATB(Npoints,Ncolumns,Nlevels) ! 3D Attenuated backscatter
      REAL(wp) :: ATBperp(Npoints,Ncolumns,Nlevels) ! 3D perpendicular attenuated backscatter
      REAL(wp) :: pplay(Npoints,Nlevels)
      REAL(wp) :: S_att,S_cld
      REAL(wp) :: undef

! Local variables
      INTEGER :: ip, k, iz, ic, ncol, nlev, i, itemp  ! loop indice

      INTEGER :: nphase ! nb of cloud layer phase types 
                                      ! (ice,liquid,undefined,false ice,false liquid,Percent of ice)
      INTEGER,PARAMETER  ::  Ntemp = 40 ! indice of the temperature vector !CNam: Ntemp not model Nlevels. See cosp_constants.
      REAL(wp) :: S_cld_att ! New threshold for undefine cloud phase detection (SR=30)
      INTEGER  :: toplvlsat  ! level of the first cloud with SR>30
      REAL(wp) :: alpha50, beta50, gamma50, delta50, epsilon50, zeta50 ! Polynomial Coef of the phase
                                                              ! discrimination line   

      REAL(wp) :: tmpi(Npoints,Ncolumns,Nlevels)	! temperature of ice cld
      REAL(wp) :: tmpl(Npoints,Ncolumns,Nlevels)	! temperature of liquid cld
      REAL(wp) :: tmpu(Npoints,Ncolumns,Nlevels)	! temperature of undef cld

      REAL(wp) :: checktemp, ATBperp_tmp ! temporary variable
      REAL(wp) :: checkcldlayerphase, checkcldlayerphase2 ! temporary variable
      REAL(wp) :: sumlidarcldtemp(Npoints,Ntemp) ! temporary variable

      REAL(wp) :: cldlayphase(Npoints,Ncolumns,Ncat,Nphase) ! subgrided low mid high phase cloud fraction
      REAL(wp) :: cldlayerphasetmp(Npoints,Ncat) ! temporary variable
      REAL(wp) :: cldlayerphasesum(Npoints,Ncat) ! temporary variable
      REAL(wp) :: lidarcldtempind(Npoints,Ntemp) ! 3D Temperature indice
      REAL(wp) :: lidarcldphasetmp(Npoints,Nlevels)  ! 3D sum of ice and liquid cloud occurences

      REAL(wp) :: p1
      REAL(wp) :: cldy(Npoints,Ncolumns,Nlevels)
      REAL(wp) :: srok(Npoints,Ncolumns,Nlevels)
      REAL(wp) :: cldlay(Npoints,Ncolumns,Ncat)
      REAL(wp) :: nsublay(Npoints,Ncolumns,Ncat), nsublayer(Npoints,Ncat)
      REAL(wp) :: nsub(Npoints,Nlevels)

!CNam: commented out SYS_SX related stuff.
!#ifdef SYS_SX
!      REAL(wp) cldlay1(Npoints,Ncolumns)
!      REAL(wp) cldlay2(Npoints,Ncolumns)
!      REAL(wp) cldlay3(Npoints,Ncolumns)
!      REAL(wp) nsublay1(Npoints,Ncolumns)
!      REAL(wp) nsublay2(Npoints,Ncolumns)
!      REAL(wp) nsublay3(Npoints,Ncolumns)
!#endif

! Output :
      REAL(wp) :: lidarcldtemp(Npoints,Ntemp,5) ! 3D cloud fraction
      REAL(wp) :: tempmod(Ntemp+1)     ! temperature bins
      REAL(wp) :: lidarcldphase(Npoints,Nlevels,Nphase)    ! 3D cloud phase fraction
      REAL(wp) :: cldlayerphase(Npoints,Ncat,Nphase) ! low, middle, high, total cloud fractions for ice liquid and undefine phase
      REAL(wp) :: lidarcld(Npoints,Nlevels) ! 3D cloud fraction
      REAL(wp) :: cldlayer(Npoints,Ncat)    ! low, middle, high, total cloud fractions


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

      ATBperp_tmp = 0._wp
      lidarcldphase(:,:,:) = 0._wp
      cldlayphase(:,:,:,:) = 0._wp
      cldlayerphase(:,:,:) = 0._wp
      tmpi(:,:,:) = 0._wp
      tmpl(:,:,:) = 0._wp
      tmpu(:,:,:) = 0._wp
      cldlayerphasesum(:,:) = 0._wp
      lidarcldtemp(:,:,:) = 0._wp
      lidarcldtempind(:,:) = 0._wp
      sumlidarcldtemp(:,:) = 0._wp
      toplvlsat=0
      lidarcldphasetmp(:,:) = 0._wp
      cldy(:,:,:) = 0._wp

! temperature bins
      tempmod=(/-273.15_wp,-90._wp,-87._wp,-84._wp,-81._wp,-78._wp,-75._wp,-72._wp,-69._wp,-66._wp,-63._wp,-60._wp,-57._wp, &
                -54._wp,-51._wp,-48._wp,-45._wp,-42._wp,-39._wp,-36._wp,-33._wp,-30._wp,-27._wp,-24._wp,-21._wp,-18._wp,  &
                -15._wp,-12._wp,-9._wp,-6._wp,-3._wp,0._wp,3._wp,6._wp,9._wp,12._wp,15._wp,18._wp,21._wp,24._wp,200._wp /)
	
! convert C to K
      tempmod=tempmod+273.15_wp

! Polynomial coefficient of the phase discrimination line used to separate liquid from ice
! (Cesana and Chepfer, JGR, 2013)
! ATBperp = ATB^5*alpha50 + ATB^4*beta50 + ATB^3*gamma50 + ATB^2*delta50 + ATB*epsilon50 + zeta50
      alpha50   = 9.0322e+15_wp
      beta50    = -2.1358e+12_wp
      gamma50   = 173.3963e06_wp
      delta50   = -3.9514e03_wp
      epsilon50 = 0.2559_wp
      zeta50    = -9.4776e-07_wp

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
      lidarcld = 0.0_wp
      nsub = 0.0_wp

!CNam: commented out SYS_SX related stuff.
!#ifdef SYS_SX
!! Use cldlay[1-3] and nsublay[1-3] to avoid bank-conflicts.
!      cldlay1 = 0.0_wp
!      cldlay2 = 0.0_wp
!      cldlay3 = 0.0_wp
!      cldlay(:,:,4) = 0.0_wp ! Ncat == 4
!      nsublay1 = 0.0_wp
!      nsublay2 = 0.0_wp
!      nsublay3 = 0.0_wp
!      nsublay(:,:,4) = 0.0_wp
!
!      do k = Nlevels, 1, -1
!       do ic = 1, Ncolumns
!        do ip = 1, Npoints
!
!        if(srok(ip,ic,k).gt.0_wp.)then
!           ! Computation of the cloud fraction as a function of the temperature
!           ! instead of height, for ice,liquid and all clouds
!           do itemp=1,Ntemp
!             if( (tmp(ip,k).ge.tempmod(itemp)).and.(tmp(ip,k).lt.tempmod(itemp+1)) )then
!               lidarcldtempind(ip,itemp)=lidarcldtempind(ip,itemp)+1._wp
!             endif
!           enddo
!         endif
!
!         if (cldy(ip,ic,k).eq.1._wp) then
!           do itemp=1,Ntemp
!             if( (tmp(ip,k).ge.tempmod(itemp)).and.(tmp(ip,k).lt.tempmod(itemp+1)) )then
!               lidarcldtemp(ip,itemp,1)=lidarcldtemp(ip,itemp,1)+1._wp
!             endif
!           enddo
!         endif
!
!         p1 = pplay(ip,k)
!
!         if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then ! high clouds !CNam: Orig p1.gt.0._wp
!            cldlay3(ip,ic) = MAX(cldlay3(ip,ic), cldy(ip,ic,k))
!            nsublay3(ip,ic) = MAX(nsublay3(ip,ic), srok(ip,ic,k))
!         else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then  ! mid clouds
!            cldlay2(ip,ic) = MAX(cldlay2(ip,ic), cldy(ip,ic,k))
!            nsublay2(ip,ic) = MAX(nsublay2(ip,ic), srok(ip,ic,k))
!         else
!            cldlay1(ip,ic) = MAX(cldlay1(ip,ic), cldy(ip,ic,k))
!            nsublay1(ip,ic) = MAX(nsublay1(ip,ic), srok(ip,ic,k))
!         endif
!
!         cldlay(ip,ic,4) = MAX(cldlay(ip,ic,4), cldy(ip,ic,k))
!         lidarcld(ip,k)=lidarcld(ip,k) + cldy(ip,ic,k)
!         nsublay(ip,ic,4) = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
!         nsub(ip,k)=nsub(ip,k) + srok(ip,ic,k)
!        enddo
!       enddo
!      enddo
!      cldlay(:,:,1) = cldlay1
!      cldlay(:,:,2) = cldlay2
!      cldlay(:,:,3) = cldlay3
!      nsublay(:,:,1) = nsublay1
!      nsublay(:,:,2) = nsublay2
!      nsublay(:,:,3) = nsublay3
! #else
      cldlay = 0.0_wp
      nsublay = 0.0_wp
      do k = Nlevels, 1, -1
       do ic = 1, Ncolumns
        do ip = 1, Npoints

          ! Computation of the cloud fraction as a function of the temperature
          ! instead of height, for ice,liquid and all clouds
        if(srok(ip,ic,k).gt. 0._wp)then
         do itemp=1,Ntemp
           if( (tmp(ip,k).ge.tempmod(itemp)).and.(tmp(ip,k).lt.tempmod(itemp+1)) )then
             lidarcldtempind(ip,itemp)=lidarcldtempind(ip,itemp)+1._wp
           endif
         enddo
        endif

        if(cldy(ip,ic,k).eq.1._wp)then
         do itemp=1,Ntemp
           if( (tmp(ip,k).ge.tempmod(itemp)).and.(tmp(ip,k).lt.tempmod(itemp+1)) )then
             lidarcldtemp(ip,itemp,1)=lidarcldtemp(ip,itemp,1)+1._wp
           endif
         enddo
         endif


          iz=1
          p1 = pplay(ip,k)         
          IF ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) THEN ! high clouds !CNam: Orig p1.gt.0._wp
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
!#endif ! SYS_SX


! -- grid-box 3D cloud fraction

      WHERE ( nsub(:,:).gt. 0.0_wp )
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

      WHERE ( nsublayer(:,:).gt. 0.0_wp )
         cldlayer(:,:) = cldlayer(:,:)/nsublayer(:,:)
      ELSEWHERE
         cldlayer(:,:) = undef
      ENDWHERE


! ---------------------------------------------------------------
! 4- grid-box 3D cloud Phase :
! ---------------------------------------------------------------
! ---------------------------------------------------------------
! 4.1 - For Cloudy pixels with 8.16km < z < 19.2km
! ---------------------------------------------------------------
  do ncol=1,Ncolumns
  do i=1,Npoints
      do nlev=18,1,-1             !CNam Orig: Nlevels,18,-1  ! from 19.2km until 8.16km
         p1 = pplay(i,nlev)

! Avoid zero values
	if( (cldy(i,ncol,nlev).eq. 1._wp) .and. (ATBperp(i,ncol,nlev).gt. 0._wp) )then
! Computation of the ATBperp along the phase discrimination line
           ATBperp_tmp = (ATB(i,ncol,nlev)**5)*alpha50 + (ATB(i,ncol,nlev)**4)*beta50 + &
                         (ATB(i,ncol,nlev)**3)*gamma50 + (ATB(i,ncol,nlev)**2)*delta50 + &
                          ATB(i,ncol,nlev)*epsilon50 + zeta50

!____________________________________________________________________________________________________
!
!4.1.a Ice: ATBperp above the phase discrimination line
!____________________________________________________________________________________________________
!
           if( (ATBperp(i,ncol,nlev)-ATBperp_tmp).ge. 0._wp )then   ! Ice clouds
             ! ICE with temperature above 273,15°K = Liquid (false ice)
            if(tmp(i,nlev).gt. 273.15_wp)then                ! Temperature above 273,15 K
              ! Liquid: False ice corrected by the temperature to Liquid
               lidarcldphase(i,nlev,2)=lidarcldphase(i,nlev,2)+1._wp   ! false ice detection ==> added to Liquid
               tmpl(i,ncol,nlev)=tmp(i,nlev)
               lidarcldphase(i,nlev,5)=lidarcldphase(i,nlev,5)+1._wp   ! keep the information "temperature criterium used"
                                                    ! to classify the phase cloud
         	   cldlayphase(i,ncol,4,2) = 1._wp                       ! tot cloud
                if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
        	   cldlayphase(i,ncol,3,2) = 1._wp
         	else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,2) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,2) = 1._wp
                endif
         	   cldlayphase(i,ncol,4,5) = 1._wp                       ! tot cloud
         	if ( p1.ge.0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
        	   cldlayphase(i,ncol,3,5) = 1._wp
         	else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,5) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,5) = 1._wp
                endif

             else
             ! ICE with temperature below 273,15°K
              lidarcldphase(i,nlev,1)=lidarcldphase(i,nlev,1)+1._wp
              tmpi(i,ncol,nlev)=tmp(i,nlev)
         	   cldlayphase(i,ncol,4,1) = 1._wp                       ! tot cloud
         	if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
        	   cldlayphase(i,ncol,3,1) = 1._wp
         	else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,1) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,1) = 1._wp
                endif

              endif

!____________________________________________________________________________________________________
!
! 4.1.b Liquid: ATBperp below the phase discrimination line
!____________________________________________________________________________________________________
!
             else                                        ! Liquid clouds
              ! Liquid with temperature above 231,15°K
            if(tmp(i,nlev).gt. 231.15_wp)then 
               lidarcldphase(i,nlev,2)=lidarcldphase(i,nlev,2)+1._wp
               tmpl(i,ncol,nlev)=tmp(i,nlev)
         	   cldlayphase(i,ncol,4,2) = 1._wp                       ! tot cloud
         	if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
         	   cldlayphase(i,ncol,3,2) = 1._wp  
         	else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,2) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,2) = 1._wp
	 	endif

             else
             ! Liquid with temperature below 231,15°K = Ice (false liquid)
               tmpi(i,ncol,nlev)=tmp(i,nlev)
               lidarcldphase(i,nlev,1)=lidarcldphase(i,nlev,1)+1._wp   ! false liquid detection ==> added to ice
               lidarcldphase(i,nlev,4)=lidarcldphase(i,nlev,4)+1._wp   ! keep the information "temperature criterium used"
                                                    ! to classify the phase cloud
         	   cldlayphase(i,ncol,4,4) = 1._wp                       ! tot cloud
        	if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
         	   cldlayphase(i,ncol,3,4) = 1._wp  
         	else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,4) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,4) = 1._wp
	 	endif
         	   cldlayphase(i,ncol,4,1) = 1._wp                       ! tot cloud
        	if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
         	   cldlayphase(i,ncol,3,1) = 1._wp  
         	else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,1) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,1) = 1._wp
	 	endif

             endif

            endif  ! end of discrimination condition 
	 endif  ! end of cloud condition
      enddo ! end of altitude loop



! ---------------------------------------------------------------
! 4.2 - For Cloudy pixels with 0km < z < 8.16km
! ---------------------------------------------------------------

      toplvlsat=0
      do nlev=Nlevels,17,-1    !CNam: Orig 17,1,-1  ! from 8.16km until 0km

	if( (cldy(i,ncol,nlev).eq. 1._wp) .and. (ATBperp(i,ncol,nlev).gt. 0._wp) )then
! Phase discrimination line : ATBperp = ATB^5*alpha50 + ATB^4*beta50 + ATB^3*gamma50 + ATB^2*delta50 
!                                  + ATB*epsilon50 + zeta50
! Computation of the ATBperp of the phase discrimination line
           ATBperp_tmp = (ATB(i,ncol,nlev)**5)*alpha50 + (ATB(i,ncol,nlev)**4)*beta50 + &
                         (ATB(i,ncol,nlev)**3)*gamma50 + (ATB(i,ncol,nlev)**2)*delta50 + &
                          ATB(i,ncol,nlev)*epsilon50 + zeta50
!____________________________________________________________________________________________________
!
! 4.2.a Ice: ATBperp above the phase discrimination line
!____________________________________________________________________________________________________
!
            ! ICE with temperature above 273,15°K = Liquid (false ice)
          if( (ATBperp(i,ncol,nlev)-ATBperp_tmp).ge. 0._wp )then   ! Ice clouds
            if(tmp(i,nlev).gt. 273.15_wp)then 
               lidarcldphase(i,nlev,2)=lidarcldphase(i,nlev,2)+1._wp  ! false ice ==> liq
               tmpl(i,ncol,nlev)=tmp(i,nlev)
               lidarcldphase(i,nlev,5)=lidarcldphase(i,nlev,5)+1._wp

         	   cldlayphase(i,ncol,4,2) = 1._wp                       ! tot cloud
                if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then      ! high cloud !CNam: Orig p1.gt.0._wp
        	   cldlayphase(i,ncol,3,2) = 1._wp
         	elseif (p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,2) = 1._wp 
	 	else                                                       ! low cloud
         	   cldlayphase(i,ncol,1,2) = 1._wp
                endif

         	   cldlayphase(i,ncol,4,5) = 1._wp                       ! tot cloud
        	if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
        	   cldlayphase(i,ncol,3,5) = 1._wp
         	elseif (p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,5) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,5) = 1._wp
                endif

             else
              ! ICE with temperature below 273,15°K
              lidarcldphase(i,nlev,1)=lidarcldphase(i,nlev,1)+1._wp
              tmpi(i,ncol,nlev)=tmp(i,nlev)

          	   cldlayphase(i,ncol,4,1) = 1._wp                       ! tot cloud
        	if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
        	   cldlayphase(i,ncol,3,1) = 1._wp
         	else if (p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,1) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,1) = 1._wp
                endif

              endif

!____________________________________________________________________________________________________
!
! 4.2.b Liquid: ATBperp below the phase discrimination line
!____________________________________________________________________________________________________
!
          else  
             ! Liquid with temperature above 231,15°K
            if(tmp(i,nlev).gt. 231.15_wp)then 
               lidarcldphase(i,nlev,2)=lidarcldphase(i,nlev,2)+1._wp
               tmpl(i,ncol,nlev)=tmp(i,nlev)

         	   cldlayphase(i,ncol,4,2) = 1._wp                       ! tot cloud
         	if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
         	   cldlayphase(i,ncol,3,2) = 1._wp  
         	else if (p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,2) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,2) = 1._wp
	 	endif

             else
             ! Liquid with temperature below 231,15°K = Ice (false liquid)
               tmpi(i,ncol,nlev)=tmp(i,nlev)
               lidarcldphase(i,nlev,1)=lidarcldphase(i,nlev,1)+1._wp  ! false liq ==> ice
               lidarcldphase(i,nlev,4)=lidarcldphase(i,nlev,4)+1._wp  ! false liq ==> ice

         	   cldlayphase(i,ncol,4,4) = 1._wp                       ! tot cloud
         	if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
         	   cldlayphase(i,ncol,3,4) = 1._wp  
         	else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,4) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,4) = 1._wp
	 	endif

         	   cldlayphase(i,ncol,4,1) = 1._wp                       ! tot cloud
        	if ( p1.ge.0._wp .and. p1.lt.(440._wp*100._wp)) then    ! high cloud !CNam: Orig p1.gt.0._wp
         	   cldlayphase(i,ncol,3,1) = 1._wp  
         	else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid cloud
         	   cldlayphase(i,ncol,2,1) = 1._wp
	 	else                                                    ! low cloud
         	   cldlayphase(i,ncol,1,1) = 1._wp
	 	endif

             endif
           endif  ! end of discrimination condition 

       	    toplvlsat=0

           ! Find the level of the highest cloud with SR>30
	    if(x(i,ncol,nlev).gt.S_cld_att)then	 ! SR > 30.
      		toplvlsat=nlev-1
       		goto 99 
    	    endif

	endif  ! end of cloud condition
       enddo  ! end of altitude loop

99 continue

!____________________________________________________________________________________________________
!
! Undefined phase: For a cloud located below another cloud with SR>30 
! see Cesana and Chepfer 2013 Sect.III.2
!____________________________________________________________________________________________________
!
  if(toplvlsat.ne.0)then     	
      do nlev=toplvlsat,1,-1
         p1 = pplay(i,nlev)
	if(cldy(i,ncol,nlev).eq. 1._wp)then
             lidarcldphase(i,nlev,3)=lidarcldphase(i,nlev,3)+1._wp
             tmpu(i,ncol,nlev)=tmp(i,nlev)
             cldlayphase(i,ncol,4,3) = 1._wp                        ! tot cloud
          if ( p1.ge. 0._wp .and. p1.lt.(440._wp*100._wp)) then     ! high cloud !CNam: Orig p1.gt.0._wp
             cldlayphase(i,ncol,3,3) = 1._wp
          else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then  ! mid cloud
             cldlayphase(i,ncol,2,3) = 1._wp
	  else                                                     ! low cloud
             cldlayphase(i,ncol,1,3) = 1._wp
	  endif

        endif	
      enddo
  endif
     
      toplvlsat=0

  enddo
  enddo



!____________________________________________________________________________________________________
!
! Computation of final cloud phase diagnosis
!____________________________________________________________________________________________________
!

! Compute the Ice percentage in cloud = ice/(ice+liq) as a function
! of the occurrences
  lidarcldphasetmp(:,:)=lidarcldphase(:,:,1)+lidarcldphase(:,:,2);
  WHERE (lidarcldphasetmp(:,:).gt. 0._wp)
     lidarcldphase(:,:,6)=lidarcldphase(:,:,1)/lidarcldphasetmp(:,:)
  ELSEWHERE
     lidarcldphase(:,:,6) = undef
  ENDWHERE

! Compute Phase 3D Cloud Fraction
     WHERE ( nsub(:,:).gt.0.0_wp )
       lidarcldphase(:,:,1)=lidarcldphase(:,:,1)/nsub(:,:)
       lidarcldphase(:,:,2)=lidarcldphase(:,:,2)/nsub(:,:)
       lidarcldphase(:,:,3)=lidarcldphase(:,:,3)/nsub(:,:)
       lidarcldphase(:,:,4)=lidarcldphase(:,:,4)/nsub(:,:)
       lidarcldphase(:,:,5)=lidarcldphase(:,:,5)/nsub(:,:)
     ELSEWHERE
       lidarcldphase(:,:,1) = undef
       lidarcldphase(:,:,2) = undef
       lidarcldphase(:,:,3) = undef
       lidarcldphase(:,:,4) = undef
       lidarcldphase(:,:,5) = undef
     ENDWHERE


! Compute Phase low mid high cloud fractions
    do iz = 1, Ncat
       do i=1,Nphase-3
       do ic = 1, Ncolumns
          cldlayerphase(:,iz,i)=cldlayerphase(:,iz,i) + cldlayphase(:,ic,iz,i)
          cldlayerphasesum(:,iz)=cldlayerphasesum(:,iz)+cldlayphase(:,ic,iz,i)
       enddo
      enddo
    enddo

    do iz = 1, Ncat
       do i=4,5
       do ic = 1, Ncolumns
          cldlayerphase(:,iz,i)=cldlayerphase(:,iz,i) + cldlayphase(:,ic,iz,i)          
       enddo
       enddo
    enddo
    
! Compute the Ice percentage in cloud = ice/(ice+liq)
    cldlayerphasetmp(:,:)=cldlayerphase(:,:,1)+cldlayerphase(:,:,2)
    WHERE (cldlayerphasetmp(:,:).gt. 0._wp)
       cldlayerphase(:,:,6)=cldlayerphase(:,:,1)/cldlayerphasetmp(:,:)
    ELSEWHERE
       cldlayerphase(:,:,6) = undef
    ENDWHERE

    do i=1,Nphase-1
      WHERE ( cldlayerphasesum(:,:).gt. 0.0_wp )
         cldlayerphase(:,:,i) = (cldlayerphase(:,:,i)/cldlayerphasesum(:,:)) * cldlayer(:,:) 
      ENDWHERE
    enddo


    do i=1,Npoints
       do iz=1,Ncat
          checkcldlayerphase=0._wp
          checkcldlayerphase2=0._wp

          if (cldlayerphasesum(i,iz).gt. 0.0_wp )then
             do ic=1,Nphase-3
                checkcldlayerphase=checkcldlayerphase+cldlayerphase(i,iz,ic)  
             enddo
             checkcldlayerphase2=cldlayer(i,iz)-checkcldlayerphase
             if( (checkcldlayerphase2.gt. 0.01_wp).or.(checkcldlayerphase2.lt.-0.01_wp) ) print *, checkcldlayerphase,cldlayer(i,iz)

          endif

       enddo
    enddo

    do i=1,Nphase-1
      WHERE ( nsublayer(:,:).eq. 0.0_wp )
         cldlayerphase(:,:,i) = undef
      ENDWHERE
   enddo

!____________________________________________________________________________________________________

 
! Compute Phase 3D as a function of temperature
 do nlev=1,Nlevels
   do ncol=1,Ncolumns     
     do i=1,Npoints
       do itemp=1,Ntemp
          if(tmpi(i,ncol,nlev).gt. 0._wp)then
             if( (tmpi(i,ncol,nlev).ge.tempmod(itemp)).and.(tmpi(i,ncol,nlev).lt.tempmod(itemp+1)) )then
               lidarcldtemp(i,itemp,2)=lidarcldtemp(i,itemp,2)+1._wp
             endif
          elseif(tmpl(i,ncol,nlev).gt. 0._wp)then
             if( (tmpl(i,ncol,nlev).ge.tempmod(itemp)).and.(tmpl(i,ncol,nlev).lt.tempmod(itemp+1)) )then
               lidarcldtemp(i,itemp,3)=lidarcldtemp(i,itemp,3)+1._wp
             endif
          elseif(tmpu(i,ncol,nlev).gt. 0._wp)then
             if( (tmpu(i,ncol,nlev).ge.tempmod(itemp)).and.(tmpu(i,ncol,nlev).lt.tempmod(itemp+1)) )then
       	lidarcldtemp(i,itemp,4)=lidarcldtemp(i,itemp,4)+1._wp
             endif
         endif
       enddo
     enddo
   enddo
 enddo

! Check temperature cloud fraction
 do i=1,Npoints
    do itemp=1,Ntemp
       checktemp=lidarcldtemp(i,itemp,2)+lidarcldtemp(i,itemp,3)+lidarcldtemp(i,itemp,4)

       if(checktemp.NE.lidarcldtemp(i,itemp,1))then
         print *, "in lmd_ispl_stats.f90: checktemp.NE.lidarcldtemp ",i,itemp
       endif

    enddo
 enddo

! Compute the Ice percentage in cloud = ice/(ice+liq)
!>>DN add undef to ice fraction computation in denomiator
!  sumlidarcldtemp=sum(lidarcldtemp(:,:,2:3),3)
!  sumlidarcldtemp(:,:)=lidarcldtemp(:,:,2)+lidarcldtemp(:,:,3)
 sumlidarcldtemp(:,:)=lidarcldtemp(:,:,2)+lidarcldtemp(:,:,3)+lidarcldtemp(:,:,4)
!<<DN

 WHERE(sumlidarcldtemp(:,:)>0._wp)
!>>DN add undef to ice fraction computation in nomiator
!   lidarcldtemp(:,:,5)=lidarcldtemp(:,:,2)/sumlidarcldtemp(:,:)
   lidarcldtemp(:,:,5)=(lidarcldtemp(:,:,2)+lidarcldtemp(:,:,4))/sumlidarcldtemp(:,:)
!<<DN
 ELSEWHERE
   lidarcldtemp(:,:,5)=undef
 ENDWHERE

 do i=1,4
   WHERE(lidarcldtempind(:,:).gt.0._wp)
      lidarcldtemp(:,:,i) = lidarcldtemp(:,:,i)/lidarcldtempind(:,:)
   ELSEWHERE
      lidarcldtemp(:,:,i) = undef
   ENDWHERE
 enddo


  RETURN
  END SUBROUTINE cosp_cldfrac
! ---------------------------------------------------------------
	  
END MODULE mo_cosp_lmd_ipsl_stats
