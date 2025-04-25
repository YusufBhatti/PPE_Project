! (c) 2008, Lawrence Livermore National Security Limited Liability Corporation.
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 15:08:38 +0100 (Mi, 13. Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/llnl/llnl_stats.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list 
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Lawrence Livermore National Security Limited Liability Corporation 
!       nor the names of its contributors may be used to endorse or promote products derived from 
!       this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
!
! Jan 2013 - G. Cesana        - Added betaperp_tot and temp_tot arguments 
!

MODULE mo_cosp_llnl_stats

  USE mo_kind, ONLY: wp
  USE mo_cosp_constants

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- FUNCTION COSP_CFAD ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION cosp_cfad(Npoints,Ncolumns,Nlevels,Nbins,x,xmin,xmax,bmin,bwidth)

   ! Input arguments
   INTEGER, INTENT(in) :: Npoints,Ncolumns,Nlevels,Nbins
   REAL(wp), INTENT(in) :: x(Npoints,Ncolumns,Nlevels)
   REAL(wp), INTENT(in) :: xmin,xmax 
   REAL(wp), INTENT(in) :: bmin,bwidth

   REAL(wp) :: cosp_cfad(Npoints,Nbins,Nlevels)
   REAL(wp) :: dbbval_ext(0:Nbins)
   ! Local variables
   INTEGER :: i, j, k
   INTEGER :: ibin
   
   !--- Input arguments
   ! Npoints: Number of horizontal points
   ! Ncolumns: Number of subcolumns
   ! Nlevels: Number of levels
   ! Nbins: Number of x axis bins
   ! x: variable to process (Npoints,Ncolumns,Nlevels)
   ! xmin: minimum value allowed for x
   ! xmax: maximum value allowed for x
   ! bmin: mimumum value of first bin
   ! bwidth: bin width
   !
   ! Output: 2D histogram on each horizontal point (Npoints,Nbins,Nlevels)

   cosp_cfad(:,:,:) = 0.0_wp
   ! bwidth intervals in the range [bmin,bmax=bmin+Nbins*hwidth]
   ! Valid x values smaller than bmin and larger than bmax are set 
   ! into the smallest bin and largest bin, respectively.

   DO j = 1, Nlevels, 1
      DO k = 1, Ncolumns, 1
         DO i = 1, Npoints, 1 
            IF ((x(i,k,j) >= xmin) .and. (x(i,k,j) <= xmax)) THEN !CNam switched order
                  ibin = ceiling((x(i,k,j) - bmin)/bwidth)
                  IF (ibin > Nbins) ibin = Nbins
                  IF (ibin < 1)     ibin = 1
                  cosp_cfad(i,ibin,j) = cosp_cfad(i,ibin,j) + 1.0_wp 
            ELSEIF (x(i,k,j) == R_GROUND) THEN
                  cosp_cfad(i,:,j) = R_UNDEF
            ENDIF ! 
         ENDDO  !i
      ENDDO  !k
   ENDDO  !j

   !CNam: was   cosp_cfad = cosp_cfad / REAL(Ncolumns,wp)
   where ((cosp_cfad /= R_UNDEF).and.(cosp_cfad /= 0.0_wp)) cosp_cfad = cosp_cfad / REAL(Ncolumns,wp)

END FUNCTION cosp_cfad


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------- SUBROUTINE COSP_LIDAR_ONLY_CLOUD -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE cosp_lidar_only_cloud(Npoints,Ncolumns,Nlevels,temp_tot,beta_tot, &
                   betaperp_tot,beta_mol,Ze_tot,lidar_only_freq_cloud,tcc)

   USE mo_kind,           ONLY: wp

   implicit none

   ! Input arguments
   INTEGER,intent(in) :: Npoints,Ncolumns,Nlevels
   REAL(wp), dimension(Npoints,Nlevels),intent(in) :: beta_mol   ! Molecular backscatter
   REAL(wp), dimension(Npoints,Ncolumns,Nlevels),intent(in) :: beta_tot   ! Total backscattered signal
   REAL(wp), dimension(Npoints,Ncolumns,Nlevels),intent(in) :: temp_tot   ! Total backscattered signal
   REAL(wp), dimension(Npoints,Ncolumns,Nlevels),intent(in) :: betaperp_tot   ! perpendicular Total backscattered signal
   REAL(wp), dimension(Npoints,Ncolumns,Nlevels),intent(in) :: Ze_tot     ! Radar reflectivity
   ! Output arguments
   REAL(wp), dimension(Npoints,Nlevels),intent(out) :: lidar_only_freq_cloud
   REAL(wp), dimension(Npoints),intent(out) :: tcc

   ! local variables
   REAL(wp) :: sc_ratio
   REAL(wp) :: s_cld, s_att
   parameter (S_cld = 5.0_wp)
   parameter (s_att = 0.01_wp)
   INTEGER :: flag_sat !first saturated level encountered from top
   INTEGER :: flag_cld !cloudy column
   INTEGER :: pr,i,j

   lidar_only_freq_cloud = 0.0_wp
   tcc = 0.0_wp
   do pr=1,Npoints
     do i=1,Ncolumns
       flag_sat = 0
       flag_cld = 0
       do j=Nlevels,1,-1 !top->surf
        sc_ratio = beta_tot(pr,i,j)/beta_mol(pr,j)
        if ((sc_ratio .le. s_att) .and. (flag_sat .eq. 0)) flag_sat = j
        if (Ze_tot(pr,i,j) .lt. -30._wp) then  !radar can't detect cloud
         if ( (sc_ratio .gt. s_cld) .or. (flag_sat .eq. j) ) then  !lidar sense cloud
            lidar_only_freq_cloud(pr,j)=lidar_only_freq_cloud(pr,j)+1._wp !top->surf
            flag_cld=1
         endif
        else  !radar sense cloud (z%Ze_tot(pr,i,j) .ge. -30.)
           flag_cld=1
        endif
       enddo !levels
       if (flag_cld .eq. 1) tcc(pr)=tcc(pr)+1._wp
     enddo !columns
   enddo !points
   lidar_only_freq_cloud=lidar_only_freq_cloud/Ncolumns
   tcc=tcc/Ncolumns

END SUBROUTINE COSP_LIDAR_ONLY_CLOUD

END MODULE mo_cosp_llnl_stats
