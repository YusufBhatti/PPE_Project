! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 15:08:38 +0100 (Wed, 13 Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/cosp_simulator.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
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
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Jan 2013 - G. Cesana - Add new variables linked to the lidar cloud phase 
!

#include "cosp_defs.h" 
MODULE mo_cosp_v1p4_cosp_simulator
  USE mo_cosp_constants, ONLY: I_RADAR, I_LIDAR, I_ISCCP, I_MISR, I_MODIS, &
                                I_RTTOV, I_STATS, tsim
  USE mo_cosp_types
  USE mo_cosp_radar
  USE mo_cosp_lidar
  USE mo_cosp_isccp
  USE mo_cosp_modis
  USE mo_cosp_misr
!#ifdef RTTOV
!  USE mo_cosp_rttov
!#endif
  USE mo_cosp_stats
  USE mo_kind, only:wp

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_SIMULATOR ------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!CNam: #ifdef RTTOV
! SUBROUTINE COSP_SIMULATOR(kbdim,Ncolumns,klev,gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar) 
! #else
SUBROUTINE COSP_SIMULATOR(kbdim,Ncolumns,klev,gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
! #endif

  ! Arguments
  INTEGER, INTENT(IN) :: kbdim,Ncolumns,klev
  type(cosp_gridbox),intent(inout) :: gbx      ! Grid-box inputs
  type(cosp_subgrid),intent(in) :: sgx      ! Subgrid inputs
  type(cosp_sghydro),intent(in) :: sghydro  ! Subgrid info for hydrometeors
  type(cosp_config),intent(in)  :: cfg      ! Configuration options
  type(cosp_vgrid),intent(in)   :: vgrid    ! Information on vertical grid of stats 
  type(cosp_sgradar),intent(inout) :: sgradar ! Output from radar simulator
  type(cosp_sglidar),intent(inout) :: sglidar ! Output from lidar simulator
  type(cosp_isccpstats),intent(inout)   :: isccp   ! Output from ISCCP simulator
  type(cosp_misrstats),intent(inout)    :: misr    ! Output from MISR simulator
  type(cosp_modisstats),intent(inout)   :: modis   ! Output from MODIS simulator
!CNam: #ifdef RTTOV
!  type(cosp_rttov),intent(inout)    :: rttov    ! Output from RTTOV
!#endif
  type(cosp_radarstats),intent(inout) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats),intent(inout) :: stlidar ! Summary statistics from lidar simulator
  ! Local variables
  integer :: i,j,k
  logical :: inconsistent

  inconsistent=.false.


  !+++++++++ Radar model ++++++++++
  if (cfg%Lradar_sim) then
    call cosp_radar(gbx,sgx,sghydro,sgradar)
  endif

  !+++++++++ Lidar model ++++++++++
  if (cfg%Llidar_sim) then
    call cosp_lidar(gbx,sgx,sghydro,sglidar)
  endif

  !+++++++++ ISCCP simulator ++++++++++
  if (cfg%Lisccp_sim) then
    call cosp_isccp(gbx,sgx,isccp)
  endif

  !+++++++++ MISR simulator ++++++++++
  if (cfg%Lmisr_sim) then
    call cosp_misr(gbx,sgx,misr)
  endif

  !+++++++++ MODIS simulator ++++++++++
  if (cfg%Lmodis_sim) then
    call cosp_modis(gbx,sgx,sghydro,isccp, modis)
  endif

  !+++++++++ RTTOV ++++++++++ 
!CNam: #ifdef RTTOV
!  if (cfg%Lrttov_sim) then 
!    call cosp_rttov_simulator(gbx,rttov)
!  endif
!#endif

  !+++++++++++ Summary statistics +++++++++++
  if (cfg%Lstats) then
    call cosp_stats(gbx,sgx,cfg,sgradar,sglidar,vgrid,stradar,stlidar)
  endif

  !+++++++++++ Change of units after computation of statistics +++++++++++ 
  ! This avoids using UDUNITS in CMOR
  ! CNam: removed unit conversion & extra logicals.

END SUBROUTINE cosp_simulator

END MODULE mo_cosp_v1p4_cosp_simulator
