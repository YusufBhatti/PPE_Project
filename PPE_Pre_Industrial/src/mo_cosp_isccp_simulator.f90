!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
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

MODULE mo_cosp_isccp_simulator

  USE mo_cosp_constants
  USE mo_kind, only:wp

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------- SUBROUTINE COSP_ISCCP_SIMULATOR -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE cosp_isccp_simulator(                          &
             kbdim,  klev,  Ncolumns,                     &   
!gbx
             isccp_top_height, isccp_top_height_direction, &
             isccp_overlap, sh, tca, dtau_s, T, dem_s,    &
             p, ph, sunlit, skt, isccp_emsfc_lw,          &
!sgx
             frac_out,                                    &
!isccp
             fq_isccp,  totalcldarea,                     &
             meanptop,  meantaucld,                       &
             meantb, meantbclr, boxtau,  boxptop,         &
             meanalbedocld )

  ! Arguments
   INTEGER, INTENT(IN) :: kbdim, klev, Ncolumns

!gbx
   INTEGER, INTENT(IN) :: isccp_top_height
   INTEGER, INTENT(IN) :: isccp_top_height_direction
   INTEGER, INTENT(IN) :: isccp_overlap

   REAL(wp), INTENT(INOUT) :: p(kbdim,klev) 
   REAL(wp), INTENT(INOUT) :: ph(kbdim,klev)
   REAL(wp), INTENT(IN)    :: isccp_emsfc_lw 
   REAL(wp), INTENT(INOUT) :: sh(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: tca(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: dem_s(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: dtau_s(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: sunlit(kbdim)
   REAL(wp), INTENT(INOUT) :: T(kbdim,klev) 
   REAL(wp), INTENT(INOUT) :: skt(kbdim)

! isccp
   REAL(wp), INTENT(INOUT) :: fq_isccp(kbdim,7,7), totalcldarea(kbdim)
   REAL(wp), INTENT(INOUT) :: meanptop(kbdim), meantaucld(kbdim)
   REAL(wp), INTENT(INOUT) :: meantb(kbdim), meantbclr(kbdim)
   REAL(wp), INTENT(INOUT) :: boxtau(kbdim,Ncolumns), boxptop(kbdim,Ncolumns)
   REAL(wp), INTENT(INOUT) :: meanalbedocld(kbdim)
  
! subgrid sgx
   REAL(wp), INTENT(IN) :: frac_out(kbdim,Ncolumns,klev)

  ! Local variables 
  REAL(wp) :: pfull(kbdim, klev)
  REAL(wp) :: phalf(kbdim, klev + 1)
  REAL(wp) :: qv(kbdim, klev)
  REAL(wp) :: cc(kbdim, klev)
  REAL(wp) :: dtau_s_l(kbdim, klev)
  REAL(wp) :: at(kbdim, klev)
  REAL(wp) :: dem_s_l(kbdim, klev)
  REAL(wp) :: frac_out_l(kbdim, Ncolumns, klev)
  INTEGER :: sunlit_l(kbdim)
  
  ! Flip inputs. Levels from TOA to surface
  pfull  = p(:,klev:1:-1) 
  phalf(:,1)         = 0.0_wp ! Top level
  phalf(:,2:klev+1) = ph(:,klev:1:-1)
  qv     = sh(:,klev:1:-1) 
  cc     = 0.999999_wp*tca(:,klev:1:-1) 
  dtau_s_l = dtau_s(:,klev:1:-1) 
  at     = T(:,klev:1:-1) 
  dem_s_l  = dem_s(:,klev:1:-1) 
  frac_out_l(1:Kbdim,:,1:klev) = frac_out(1:Kbdim,:,klev:1:-1)
  sunlit_l = int(sunlit)

  CALL icarus(0,0,kbdim,sunlit_l,klev,ncolumns, &
            pfull,phalf,qv,cc,dtau_s_l, &
            isccp_top_height,isccp_top_height_direction, &
            isccp_overlap,frac_out_l, &
            skt,isccp_emsfc_lw,at,dem_s_l,fq_isccp,totalcldarea, &
            meanptop,meantaucld,meanalbedocld, &
            meantb,meantbclr,boxtau,boxptop)

  ! Flip outputs. Levels from surface to TOA
  ! --- (kbdim,tau=7,pressure=7)
  fq_isccp(:,:,:) = fq_isccp(:,:,7:1:-1)
      
  ! Check if there is any value slightly greater than 1
  WHERE ((totalcldarea > 1.0_wp-1.e-5_wp) .and. (totalcldarea < 1.0_wp+1.e-5_wp))
    totalcldarea = 1.0_wp
  ENDWHERE
              
END SUBROUTINE cosp_isccp_simulator

END MODULE mo_cosp_isccp_simulator
