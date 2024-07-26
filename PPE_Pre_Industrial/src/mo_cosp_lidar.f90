!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
! (c) COPYRIGHT British Crown / Met Office 2008  
! Please refer to Met_Office_license.txt for details.
!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Oct 2008 - S. Bony          - Instructions "Call for large-scale cloud" removed  -> sgx%frac_out is used instead.
!                               Call lidar_simulator changed (lsca, gbx%cca and depol removed; 
!                               frac_out changed in sgx%frac_out)
!
! 
MODULE mo_cosp_lidar

  USE mo_cosp_constants
  USE mo_kind,          ONLY: wp

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_LIDAR ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE cosp_lidar( kbdim,Ncolumns,klev,            &
!gbx
                      lidar_ice_type, Reff,            &
                      ph, mr_hydro, p, T,         &
!sgx, sglidar                       
!                      frac_out,                        &
                      beta_mol, beta_tot, tau_tot, refl )
  
  ! Arguments
   INTEGER, INTENT(IN) :: kbdim,Ncolumns,klev
 
!gbx
   INTEGER, INTENT(IN) :: lidar_ice_type
   REAL(wp), INTENT(INOUT) :: Reff(kbdim,klev,N_hydro)
   REAL(wp), INTENT(INOUT) :: p(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: ph(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: T(kbdim,klev) 

   REAL(wp), INTENT(IN) :: mr_hydro(kbdim,Ncolumns,klev,N_hydro) !sghydro?

!sglidar
   REAL(wp) :: beta_tot(kbdim,Ncolumns,klev)
   REAL(wp) :: tau_tot(kbdim,Ncolumns,klev)
   REAL(wp) :: refl(kbdim,Ncolumns,PARASOL_NREFL)
   REAL(wp) :: beta_mol(kbdim,klev)

! subgrid sgx
!   REAL(wp) :: frac_out(kbdim,Ncolumns,klev)

  ! Local variables 
  INTEGER :: i
  REAL(wp) :: presf(kbdim, klev + 1)
  REAL(wp), DIMENSION(kbdim, klev) :: mr_ll,mr_li
  REAL(wp), DIMENSION(kbdim, klev) :: beta_tot_l,tau_tot_l
  REAL(wp), DIMENSION(kbdim, PARASOL_NREFL)  :: refle_l
 
  presf(:,1:klev) = ph
  presf(:,klev + 1) = 0.0_wp
 
   beta_tot_l=0._wp
   tau_tot_l=0._wp
   refle_l=0._wp


  DO i=1,Ncolumns
      ! Temporary arrays for simulator call
      mr_ll(:,:) = mr_hydro(:,i,:,I_LSCLIQ)
      mr_li(:,:) = mr_hydro(:,i,:,I_LSCICE)
      CALL cosp_lidar_simulator(kbdim, klev, 2 ,PARASOL_NREFL  &
                 , p, presf, T & 
                 , mr_ll, mr_li &
                 , Reff(:,:,I_LSCLIQ), Reff(:,:,I_LSCICE) &
!!                 , frac_out, lidar_ice_type, beta_mol, beta_tot_l, tau_tot_l  &
                 , lidar_ice_type, beta_mol, beta_tot_l, tau_tot_l  &
                 , refle_l ) ! reflectance
      
      beta_tot(:,i,:) = beta_tot_l(:,:)
      tau_tot(:,i,:)  = tau_tot_l(:,:)
      refl(:,i,:)     = refle_l(:,:)
  ENDDO
  
END SUBROUTINE cosp_lidar

END MODULE mo_cosp_lidar
