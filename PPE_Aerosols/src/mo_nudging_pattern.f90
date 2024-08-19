!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_nudging_pattern
!BOP

  ! !MODULE: mo_nudging_pattern (layer 2)

  ! !DESCRIPTION: 

  ! \newcommand{\CLIM}{\overline{\Psi_{cl}}^{\tau}(x)}
  ! \newcommand{\PT}{\Phi_{T}(x)}
  ! \newcommand{\NORM}[1]{\frac{N[#1,\PT]}{N[\PT,\PT]}}

  ! This module contains all procedures for pattern nudging. The pattern nudging
  ! is switched on with {\it LNDGPAT=.TRUE.} and was developed by
  ! Martin Widmann and Ingo Kirchner.
  !
  ! In pattern nudging mode a single pattern can be assimilated into
  ! the model based on projections. (An extension to nudge several
  ! patterns is planned.)
  ! 
  ! Every field can be expanded using an orthogonal basis system
  ! $\Phi_i$. Note that normality of the basis vectors is not required.
  ! 
  ! \begin{equation}
  ! \Psi(x,t) = \CLIM{} + \sum^{\infty}_{i=1}\alpha_i(t)\Phi_i(x) 
  ! \end{equation}
  ! 
  ! Pattern nudging forces the model towards a {\sl target state}
  ! 
  ! \begin{equation}
  ! \Psi_{T}(x,t) = \CLIM{} + \alpha_{T}(t)\PT
  ! \end{equation}
  ! 
  ! While conventional nudging methods would prescribe $\Psi_{T}(x,t)$ and
  ! nudge the model towards this state at each point (thus nudging the
  ! model towards the climatology where $\PT = 0$, which in some
  ! applications are exactly the locations where no information is
  ! available), pattern nudging considers every state
  ! 
  ! \begin{equation}
  ! \Psi(x,t) = \CLIM{} + \alpha_{T}(t)\PT 
  !                     + \sum^{\infty}_{i=2}\alpha_i(t)\Phi_i(x)
  ! \end{equation}
  ! 
  ! to be consistent with $\Psi_{T}(x,t)$.
  ! 
  ! With $\Phi_1(x) = \PT$, any model state can be expanded as
  ! 
  ! \begin{equation}
  ! \Psi_{mod}(x,t) = \CLIM{} + \alpha_{mod}(t)\PT 
  !                           + \sum^{\infty}_{i=2}\alpha_i(t)\Phi_i(x)
  ! \end{equation}
  ! 
  ! The goal is to nudge the model such that $\alpha_{mod}(t)$ is close to
  ! $\alpha_T(t)$. In order to do this, one needs to calculate
  ! 
  ! \begin{equation}
  ! \alpha_{mod}(t) = \NORM{\Delta\Psi(x,t)}
  ! \end{equation}
  ! 
  ! with
  ! 
  ! \begin{equation}
  ! \Delta\Psi(x,t) = \Psi_{mod}(x,t) - \CLIM
  ! \end{equation}
  ! 
  ! and the norm or scalar product $N$ of two spatial fields 
  ! $X(x),Y(x)$ defined as
  ! 
  ! \begin{equation}
  ! N[X(x),Y(x)] = \frac{1}{A}\int{X(x)Y(x)} \, dA
  ! \end{equation}
  ! 
  ! In the discrete version used for calculations on the model grid this
  ! becomes 
  ! 
  ! \begin{eqnarray}
  ! N[X(x),Y(x)] &=& \frac{1}{A}\sum_{i,j}X(x_{ij})Y(x_{ij}) \Delta a_{ij} \\
  !              &=& \sum_{i,j} \frac{\Delta a_{ij}}{A} X(x_{ij})Y(x_{ij}) \\
  !              &=& \sum_{i,j} \frac{w_i}{nlon}X(x_{ij})Y(x_{ij}) 
  ! \end{eqnarray}
  ! 
  ! where $w_i$ are the Gaussian weigts and the number of longitudes
  ! 
  ! \begin{equation}
  ! \frac{\Delta a_{ij}}{A} = \frac{w_i}{nlon}
  ! \end{equation}
  ! 
  ! Equivalent the norm can be calculated in the spectral space using
  ! 
  ! \begin{equation}
  ! N[X_s, Y_s] = \sum_{i=1}^{nsp}\kappa_i\cdot X_i \cdot Y_i
  ! \end{equation}
  ! 
  ! with $\kappa_i=1$ for zonal wavenumber 0 coefficients and $\kappa_i=2$
  ! for all other wavenumbers
  ! 
  ! The  relaxation term or additional forcing term $R$ is given by
  ! 
  ! \begin{equation}
  ! R = G\,\left( \alpha_{T} - \alpha_{mod}(t)\right)\PT
  ! \end{equation}
  ! 
  ! and correction of the prognostic field ($\Psi_{mod,old}(x,t)$) with
  ! 
  ! \begin{equation}
  ! \Psi_{mod,new}(x,t) = \Psi_{mod,old}(x,t) + 2\Delta{}t\cdot R
  ! \end{equation}
  ! 
  ! Since a different $\PT$ is used for every month, a temporal
  ! interpolation between two pattern at time 1 and time 2 is needed at
  ! each timestep/
  ! 
  ! \begin{equation}
  ! R_{t_1,t_2} = G\left( \alpha_{T}(t_1,t_2) + \NORM{\CLIM} 
  !                         - \NORM{\Psi_{mod}(x,t)} \right)_{\Phi_{T}(t_1,t_2)}
  ! \end{equation}
  ! 
  ! The calculation of the relaxation term is implemented in two steps:
  ! 
  ! {\bf STEP 1}, after reading new pattern
  ! 
  ! \begin{equation}
  ! \alpha_{t_i} = \alpha_{T}(t_i) + \left(\NORM{\CLIM}\right)_{\Phi_{T}(t_i)}
  ! \end{equation}
  ! 
  ! {\bf STEP 2}, at each model time step
  ! 
  ! \begin{equation}
  ! \beta_1(t) = \alpha_{t_1} - \left(\NORM{\Psi_{mod}(x,t)}\right)_{\Phi_{T}(t_1)}
  ! \end{equation}
  ! 
  ! \begin{equation}
  ! \beta_2(t) = \alpha_{t_2} - \left(\NORM{\Psi_{mod}(x,t)}\right)_{\Phi_{T}(t_2)}
  ! \end{equation}
  ! 
  ! linear interpolation between
  ! 
  ! \begin{equation}
  ! R_1 = \beta_1(t)\Phi_{T}(t_1)
  ! \end{equation}
  ! 
  ! and
  ! 
  ! \begin{equation}
  ! R_2 = \beta_2(t)\Phi_{T}(t_2)
  ! \end{equation}
  ! 
  ! {\bf Changes in the data format}
  ! 
  ! In the pattern nudging mode additional to the default data field an
  ! second record with the climatological pattern and as a third record
  ! the level dependend time evolution coefficient is read in.

  ! !REVISION HISTORY: 
  ! Ingo Kirchner, MPI Hamburg, April-2001
  ! Ingo Kirchner, MPI Hamburg, Aug-2002 revision

!BOX
  IMPLICIT NONE

  PRIVATE
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:

  PUBLIC  :: Nudg_Init_Alpha    ! Step 1
  PUBLIC  :: Nudg_Update_Alpha  ! Step 2
!EOP

  CHARACTER(len=256) :: mess

CONTAINS

  !======================================================================
!BOP
  ! !IROUTINE:  Nudg_Update_Alpha
  ! !INTERFACE:

  SUBROUTINE Nudg_Update_Alpha

    ! !DESCRIPTION: 
    ! calculates b\_pat\_* = a\_pat\_* + NORM(model,obs)


    ! !USES:

    !* general modules
    USE mo_kind,       ONLY: dp
    USE mo_memory_sp,  ONLY: sd, svo, stp
    USE mo_control,    ONLY: nnp1, nsp, nlevp1
    USE mo_exception,  ONLY: message

    !* nudging modules
    USE mo_nudging_constants, ONLY: lnudgdbx, &
         linp_vor, ilev_v_min, ilev_v_max, &
         linp_div, ilev_d_min, ilev_d_max, &
         linp_tem, ilev_t_min, ilev_t_max, &
         linp_lnp
    USE mo_nudging_buffer, ONLY: &
         a_pat_vor, b_pat_vor, a_nrm_vor, &
         a_pat_div, b_pat_div, a_nrm_div, &
         a_pat_tem, b_pat_tem, a_nrm_tem, &
         a_pat_lnp, b_pat_lnp, a_nrm_lnp, &
         sdobs1, svobs1, stobs1, &
         sdobs2, svobs2, stobs2
!EOP

    INTEGER       :: i
    REAL(kind=dp) :: aa(2,nsp), bb(2,nsp), nnab, nnbb

    b_pat_vor(:,:) = 0.0_dp
    IF (linp_vor) THEN

       DO i=ilev_v_min,ilev_v_max
          aa(:,:) = svo   (i,:,:)
          bb(:,:) = svobs1(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          nnbb = a_nrm_vor(i,2)
          b_pat_vor(i,1) = a_pat_vor(i,2) - nnab/nnbb

          bb(:,:) = svobs2(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          nnbb = a_nrm_vor(i,3)
          b_pat_vor(i,2) = a_pat_vor(i,3) - nnab/nnbb
       END DO

       IF (lnudgdbx) THEN
          WRITE(mess,'(a,40f12.4)') 'BETA(vor)1= ',b_pat_vor(:,1)
          CALL message('Nudg_Update_Alpha',mess)
          WRITE(mess,'(a,40f12.4)') 'BETA(vor)2= ',b_pat_vor(:,2)
          CALL message('Nudg_Update_Alpha',mess)
       END IF

    END IF
    
    b_pat_div(:,:) = 0.0_dp
    IF (linp_div) THEN

       DO i=ilev_d_min,ilev_d_max
          aa(:,:) = sd    (i,:,:)
          bb(:,:) = sdobs1(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          nnbb = a_nrm_div(i,2)
          b_pat_div(i,1) = a_pat_div(i,2) - nnab/nnbb

          bb(:,:) = sdobs2(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          nnbb = a_nrm_div(i,3)
          b_pat_div(i,2) = a_pat_div(i,3) - nnab/nnbb
       END DO

       IF (lnudgdbx) THEN
          WRITE(mess,'(a,40f12.4)') 'BETA(div)1= ',b_pat_div(:,1)
          CALL message('Nudg_Update_Alpha',mess)
          WRITE(mess,'(a,40f12.4)') 'BETA(div)2= ',b_pat_div(:,2)
          CALL message('Nudg_Update_Alpha',mess)
       END IF

    END IF

    b_pat_tem(:,:) = 0.0_dp
    IF (linp_tem) THEN

       DO i=ilev_t_min,ilev_t_max
          aa(:,:) = stp   (i,:,:)
          bb(:,:) = stobs1(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          nnbb = a_nrm_tem(i,2)
          b_pat_tem(i,1) = a_pat_tem(i,2) - nnab/nnbb

          bb(:,:) = stobs2(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          nnbb = a_nrm_tem(i,3)
          b_pat_tem(i,2) = a_pat_tem(i,3) - nnab/nnbb
       END DO

       IF (lnudgdbx) THEN
          WRITE(mess,'(a,40f12.4)') 'BETA(tem)1= ',b_pat_tem(:,1)
          CALL message('Nudg_Update_Alpha',mess)
          WRITE(mess,'(a,40f12.4)') 'BETA(tem)2= ',b_pat_tem(:,2)
          CALL message('Nudg_Update_Alpha',mess)
       END IF

    END IF

    b_pat_lnp(:) = 0.0_dp
    IF (linp_lnp) THEN
       aa(:,:) = stp   (nlevp1,:,:)
       bb(:,:) = stobs1(nlevp1,:,:)
       CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
       nnbb = a_nrm_lnp(2)
       b_pat_lnp(1) = a_pat_lnp(2) - nnab/nnbb

       bb(:,:) = stobs2(nlevp1,:,:)
       CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
       nnbb = a_nrm_lnp(3)
       b_pat_lnp(2) = a_pat_lnp(3) - nnab/nnbb

       IF (lnudgdbx) THEN
          WRITE(mess,'(a,40f12.4)') 'BETA(lnps)1,2= ',b_pat_lnp(:)
          CALL message('Nudg_Update_Alpha',mess)
       END IF

    END IF

  END SUBROUTINE Nudg_Update_Alpha

  !======================================================================
!BOP
  ! !IROUTINE:  Nudg_Init_Alpha
  ! !INTERFACE:

  SUBROUTINE Nudg_Init_Alpha

    ! !DESCRIPTION: 
    ! calculates b\_pat\_* = a\_pat\_* + NORM(model,obs)

    ! !USES:

    !* general modules
    USE mo_kind,     ONLY: dp
    USE mo_control,  ONLY: nnp1, nsp, nlevp1

    !* nudging modules
    USE mo_nudging_constants, ONLY: &
         linp_vor, ilev_v_min, ilev_v_max, &
         linp_div, ilev_d_min, ilev_d_max, &
         linp_tem, ilev_t_min, ilev_t_max, &
         linp_lnp
    USE mo_nudging_buffer, ONLY: &
         a_pat_vor, a_nrm_vor, &
         a_pat_div, a_nrm_div, &
         a_pat_tem, a_nrm_tem, &
         a_pat_lnp, a_nrm_lnp, &
         sdobs0, svobs0, stobs0, &    ! block 0 contains climatology
         sdobs3, svobs3, stobs3       ! block 3 contains pattern
!EOP

    INTEGER       :: i
    REAL(kind=dp) :: aa(2,nsp), bb(2,nsp), nnab, nnbb
    
    IF (linp_vor) THEN

       DO i=ilev_v_min,ilev_v_max
          aa(:,:) = svobs3(i,:,:)
          bb(:,:) = svobs3(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnbb)
          a_nrm_vor(i,4) = nnbb

          aa(:,:) = svobs0(i,:,:)
          bb(:,:) = svobs3(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          a_pat_vor(i,4) = a_pat_vor(i,4) + nnab/nnbb
       END DO

    END IF

    IF (linp_div) THEN

       DO i=ilev_d_min,ilev_d_max
          aa(:,:) = sdobs3(i,:,:)
          bb(:,:) = sdobs3(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnbb)
          a_nrm_div(i,4) = nnbb

          aa(:,:) = sdobs0(i,:,:)
          bb(:,:) = sdobs3(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          a_pat_div(i,4) = a_pat_div(i,4) + nnab/nnbb
       END DO

    END IF

    IF (linp_tem) THEN

       DO i=ilev_t_min,ilev_t_max
          aa(:,:) = stobs3(i,:,:)
          bb(:,:) = stobs3(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnbb)
          a_nrm_tem(i,4) = nnbb

          aa(:,:) = stobs0(i,:,:)
          bb(:,:) = stobs3(i,:,:)
          CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
          a_pat_tem(i,4) = a_pat_tem(i,4) + nnab/nnbb
       END DO

    END IF

    IF (linp_lnp) THEN

       aa(:,:) = stobs3(nlevp1,:,:)
       bb(:,:) = stobs3(nlevp1,:,:)
       CALL Nudg_Norm(aa,bb,nnp1,nsp,nnbb)
       a_nrm_lnp(4) = nnbb

       aa(:,:) = stobs0(nlevp1,:,:)
       bb(:,:) = stobs3(nlevp1,:,:)
       CALL Nudg_Norm(aa,bb,nnp1,nsp,nnab)
       a_pat_lnp(4) = a_pat_lnp(4) + nnab/nnbb

    END IF

  END SUBROUTINE Nudg_Init_Alpha

  !======================================================================
!BOP
  ! !IROUTINE:  Nudg_Norm
  ! !INTERFACE:

  SUBROUTINE Nudg_Norm(a,b,no_w1,no_sp,nab)

    ! !DESCRIPTION: 
    ! calculates the norm using spectral coefficients

    ! !USES:
    USE mo_kind,    ONLY: dp

    ! !INPUT PARAMETERS: 
    REAL(kind=dp), INTENT(IN)  :: a(:,:) ! first spectral field
    REAL(kind=dp), INTENT(IN)  :: b(:,:) ! second spectral field
    INTEGER, INTENT(in)        :: no_w1  ! no. of coeffs for first wavenumber
    INTEGER, INTENT(in)        :: no_sp  ! total number of spec.-coeffs

    ! !OUTPUT PARAMETERS:
    REAL(kind=dp), INTENT(out) :: nab    ! norm  in spectralspace
!EOP

    INTEGER            :: i
    
    nab = 0.0_dp

    DO i=1,no_w1
       nab = nab + a(1,i)*b(1,i) + a(2,i)*b(2,i)
    END DO

    DO i=no_w1+1,no_sp
       nab = nab + 2.0_dp*( a(1,i)*b(1,i) + a(2,i)*b(2,i) )
    END DO
      
  END SUBROUTINE Nudg_Norm

  !======================================================================

END MODULE mo_nudging_pattern
