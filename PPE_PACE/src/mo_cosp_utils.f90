! (c) COPYRIGHT British Crown / Met Office 2008 
! Please refer to Met_Office_license.txt for details. 
!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version

MODULE mo_cosp_utils
  USE mo_cosp_constants
  USE mo_kind,		ONLY: wp
  IMPLICIT NONE

  INTERFACE Z_TO_DBZ
    MODULE PROCEDURE Z_TO_DBZ_2D,Z_TO_DBZ_3D,Z_TO_DBZ_4D
  END INTERFACE

  INTERFACE COSP_CHECK_INPUT
    MODULE PROCEDURE COSP_CHECK_INPUT_1D,COSP_CHECK_INPUT_2D,COSP_CHECK_INPUT_3D
  END INTERFACE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_PRECIP_MXRATIO --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_PRECIP_MXRATIO(Npoints,Nlevels,Ncolumns,p,T,prec_frac,prec_type, &
                          n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,gamma3,gamma4, &
                          flux,mxratio,reff)

    ! Input arguments, (IN)
    integer,intent(in) :: Npoints, Nlevels, Ncolumns
    real(wp),intent(in),dimension(Npoints,Nlevels) :: p, T, flux
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: prec_frac
    real(wp),intent(in) :: n_ax, n_bx, alpha_x, c_x, d_x, g_x, a_x, b_x, gamma1, gamma2, gamma3, gamma4, prec_type
    ! Input arguments,  (OUT)
    real(wp),intent(out),dimension(Npoints,Ncolumns,Nlevels) :: mxratio
    real(wp),intent(inout),dimension(Npoints,Ncolumns,Nlevels) :: reff
    ! Local variables
    integer :: i, j, k
    real(wp) :: sigma, one_over_xip1, xi, rho0, rho, lambda_x, gamma_4_3_2, delta
    
    mxratio = 0.0_wp

    if (n_ax >= 0.0_wp) then ! N_ax is used to control which hydrometeors need to be computed
        xi      = d_x/(alpha_x + b_x - n_bx + 1.0_wp)
        rho0    = 1.29_wp
        sigma   = (gamma2/(gamma1*c_x))*(n_ax*a_x*gamma2)**xi
        one_over_xip1 = 1.0_wp/(xi + 1.0_wp)
        gamma_4_3_2 = 0.5_wp*gamma4/gamma3
        delta = (alpha_x + b_x + d_x - n_bx + 1.0_wp)
        
        do k=1,Nlevels
            do j=1,Ncolumns
                do i=1,Npoints
                    if ((prec_frac(i,j,k)==prec_type).or.(prec_frac(i,j,k)==3._wp)) then
                        rho = p(i,k)/(287.05_wp*T(i,k))
                        mxratio(i,j,k)=(flux(i,k)*((rho/rho0)**g_x)*sigma)**one_over_xip1
                        mxratio(i,j,k)=mxratio(i,j,k)/rho
                        ! Compute effective radius
                        if ((reff(i,j,k) <= 0.0_wp).and.(flux(i,k) /= 0.0_wp)) then
                           lambda_x = (a_x*c_x*((rho0/rho)**g_x)*n_ax*gamma1/flux(i,k))**(1._wp/delta)
                           reff(i,j,k) = gamma_4_3_2/lambda_x
                        endif
                    endif
                enddo
            enddo
        enddo
    endif
END SUBROUTINE COSP_PRECIP_MXRATIO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE ZERO_INT -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELEMENTAL SUBROUTINE ZERO_INT(x, y01, y02, y03, y04, y05, y06, y07, y08, y09, y10,  &
                                 y11, y12, y13, y14, y15, y16, y17, y18, y19, y20,  &
                                 y21, y22, y23, y24, y25, y26, y27, y28, y29, y30)

  integer, intent(inout) :: x
  integer, intent(inout), optional :: y01, y02, y03, y04, y05, y06, y07, y08, y09, y10,  &
                                    y11, y12, y13, y14, y15, y16, y17, y18, y19, y20,  &
                                    y21, y22, y23, y24, y25, y26, y27, y28, y29,y30
  x = 0
  if (present(y01)) y01 = 0
  if (present(y02)) y02 = 0
  if (present(y03)) y03 = 0
  if (present(y04)) y04 = 0
  if (present(y05)) y05 = 0
  if (present(y06)) y06 = 0
  if (present(y07)) y07 = 0
  if (present(y08)) y08 = 0
  if (present(y09)) y09 = 0
  if (present(y10)) y10 = 0
  if (present(y11)) y11 = 0
  if (present(y12)) y12 = 0
  if (present(y13)) y13 = 0
  if (present(y14)) y14 = 0
  if (present(y15)) y15 = 0
  if (present(y16)) y16 = 0
  if (present(y17)) y17 = 0
  if (present(y18)) y18 = 0
  if (present(y19)) y19 = 0
  if (present(y20)) y20 = 0
  if (present(y21)) y21 = 0
  if (present(y22)) y22 = 0
  if (present(y23)) y23 = 0
  if (present(y24)) y24 = 0
  if (present(y25)) y25 = 0
  if (present(y26)) y26 = 0
  if (present(y27)) y27 = 0
  if (present(y28)) y28 = 0
  if (present(y29)) y29 = 0
  if (present(y30)) y30 = 0
END SUBROUTINE  ZERO_INT

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE ZERO_REAL ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELEMENTAL SUBROUTINE ZERO_REAL(x, y01, y02, y03, y04, y05, y06, y07, y08, y09, y10,  &
                                 y11, y12, y13, y14, y15, y16, y17, y18, y19, y20,  &
                                 y21, y22, y23, y24, y25, y26, y27, y28, y29, y30)

  real(wp), intent(inout) :: x
  real(wp), intent(inout), optional :: y01, y02, y03, y04, y05, y06, y07, y08, y09, y10,  &
                                 y11, y12, y13, y14, y15, y16, y17, y18, y19, y20,  &
                                 y21, y22, y23, y24, y25, y26, y27, y28, y29, y30
  x = 0.0_wp
  if (present(y01)) y01 = 0.0_wp
  if (present(y02)) y02 = 0.0_wp
  if (present(y03)) y03 = 0.0_wp
  if (present(y04)) y04 = 0.0_wp
  if (present(y05)) y05 = 0.0_wp
  if (present(y06)) y06 = 0.0_wp
  if (present(y07)) y07 = 0.0_wp
  if (present(y08)) y08 = 0.0_wp
  if (present(y09)) y09 = 0.0_wp
  if (present(y10)) y10 = 0.0_wp
  if (present(y11)) y11 = 0.0_wp
  if (present(y12)) y12 = 0.0_wp
  if (present(y13)) y13 = 0.0_wp
  if (present(y14)) y14 = 0.0_wp
  if (present(y15)) y15 = 0.0_wp
  if (present(y16)) y16 = 0.0_wp
  if (present(y17)) y17 = 0.0_wp
  if (present(y18)) y18 = 0.0_wp
  if (present(y19)) y19 = 0.0_wp
  if (present(y20)) y20 = 0.0_wp
  if (present(y21)) y21 = 0.0_wp
  if (present(y22)) y22 = 0.0_wp
  if (present(y23)) y23 = 0.0_wp
  if (present(y24)) y24 = 0.0_wp
  if (present(y25)) y25 = 0.0_wp
  if (present(y26)) y26 = 0.0_wp
  if (present(y27)) y27 = 0.0_wp
  if (present(y28)) y28 = 0.0_wp
  if (present(y29)) y29 = 0.0_wp
  if (present(y30)) y30 = 0.0_wp
END SUBROUTINE  ZERO_REAL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_2D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_2D(mdi,z)
    real(wp),intent(in) :: mdi
    real(wp),dimension(:,:),intent(inout) :: z
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]
    
    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18_wp*z
    where (z > 1.0e-6_wp) ! Limit to -60 dBZ
      z = 10.0_wp*log10(z)
    elsewhere
      z = mdi
    end where  
  END SUBROUTINE Z_TO_DBZ_2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_3D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_3D(mdi,z)
    real(wp),intent(in) :: mdi
    real(wp),dimension(:,:,:),intent(inout) :: z
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]
    
    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18_wp*z
    where (z > 1.0e-6_wp) ! Limit to -60 dBZ
      z = 10.0_wp*log10(z)
    elsewhere
      z = mdi
    end where  
  END SUBROUTINE Z_TO_DBZ_3D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE Z_TO_DBZ_4D --------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE Z_TO_DBZ_4D(mdi,z)
    real(wp),intent(in) :: mdi
    real(wp),dimension(:,:,:,:),intent(inout) :: z
    integer :: i, N
    ! Reflectivity Z:
    ! Input in [m3]
    ! Output in dBZ, with Z in [mm6 m-3]
    
!     N = size(z,4)
!     do i=1,N
!        call z_to_dbz_3d(mdi,z(:,:,:,i))
!     enddo
        
    ! 1.e18 to convert from [m3] to [mm6 m-3]
    z = 1.e18_wp*z
    where (z > 1.0e-6_wp) ! Limit to -60 dBZ
      z = 10.0_wp*log10(z)
    elsewhere
      z = mdi
    end where  
  END SUBROUTINE Z_TO_DBZ_4D


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_1D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_1D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real(wp),intent(inout) :: x(:)
    real(wp),intent(in),optional :: min_val, max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_1D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_1D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_2D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_2D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real(wp),intent(inout) :: x(:,:)
    real(wp),intent(in),optional :: min_val, max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_2D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_3D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_3D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real(wp),intent(inout) :: x(:,:,:)
    real(wp),intent(in),optional :: min_val, max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_3D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_3D


END MODULE mo_cosp_utils
