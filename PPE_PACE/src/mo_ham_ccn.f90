!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! mo_ham_ccn.f90
!!
!! \brief
!! mo_ham_ccn contains calculations relative to Cloud Condensation Nuclei diagnostics
!!
!! \author Zak Kipling (Univ. Oxford)
!!
!! \responsible_coder
!! Duncan Watson-Parris, duncan.watson-parris@physics.ox.ac.uk
!!
!! \revision_history
!!   -# Z. Kipling (Univ. Oxford) - original version - (2014)
!!   -# D. Watson-Parris (Univ. Oxford) - (2017-01)
!!
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_ccn

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream
  USE mo_ham,           ONLY: nclass, nccndiag
  USE mo_submodel_diag, ONLY: vmem3d,vmem2d
  USE mo_ham_ccnctl,    ONLY: nsat, zsat

  IMPLICIT NONE

  PUBLIC ham_ccn
 
  PRIVATE

CONTAINS

  SUBROUTINE ham_ccn(kproma, kbdim, klev,  klevp1, krow,  ktrac, &
                     pxtm1,  ptvm1, papm1, paphm1                )

    ! *ham_ccn* calculates the number Cloud Condensation Nuclei
    !           CCN at fixed supersaturation x
    !
    ! Author:
    ! -------
    ! Philip Stier, University of Oxford  2009
    !
    ! Method:
    ! -------
    ! The calculation of the CCNx can be reduced to 3 tasks:
    ! 
    ! I)   Use prescribed supersaturation
    ! II)  Calculate the corresponding radius of activation
    !      for each mode
    ! III) Calculate the number of particles that are larger
    !      then the radius of activation for each mode.
    ! 
    ! III) Calculation of the number of activated particles:
    !      See the routine aero_activ_tail below.
    !
    ! References:
    ! -----------
    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.

    USE mo_ham,                ONLY: sizeclass
    USE mo_ham_tools,          ONLY: ham_m7_logtail
    USE mo_ham_streams,        ONLY: a, b, rdry, &
                                     ccn_2d, ccn_burden, ccn_3d, &
                                     cn_2d, cn_burden, cn_3d
    USE mo_physical_constants, ONLY: grav, rd


    IMPLICIT NONE

    !--- Arguments:

    INTEGER, INTENT(in) :: kproma, kbdim, klev, klevp1, krow, ktrac

 
    REAL(dp), INTENT(in) :: ptvm1(kbdim,klev),   & ! virtual temperature
                            papm1(kbdim,klev),   & ! pressure 
                            paphm1(kbdim,klevp1)   ! half-level pressure 

    REAL(dp), INTENT(in) :: pxtm1(kbdim,klev,ktrac)

    !--- Local variables:

    INTEGER :: jclass, jt, jk

    REAL(dp):: zeps, ztmp

    REAL(dp):: zrho(kbdim,klev),           & ! air density
               zdpg(kbdim,klev)              ! auxiliary variable dp/g

    REAL(dp):: zra(kbdim,klev,nclass),     & ! radius of activation, i.e. radius of smallest particle with 
                                             ! critical supersaturation larger than prescribed supersaturation [m] 
               zfracn(kbdim,klev,nclass),  & ! fraction of activated aerosol numbers for each mode - stratiform
               zn(kbdim,klev,nclass)         ! aerosol number concentration for each mode [m-3]

    REAL(dp):: zccn(kbdim,klev),           & ! CCN at supersaturations zsat(jsat) [m-3]
               zcn(kbdim,klev)               ! CN(r>zrmin) at ambient conditions [m-3]

    REAL(dp) :: zrmin(kbdim,klev,nclass)     ! cut-off radius for CN calculation (typical measurement limit 3nm)

    LOGICAL, PARAMETER :: ll_numb = .TRUE.  ! switch between number/mass in logtail calculation

    INTEGER :: jsat

    INTEGER :: itoplev                        ! highest level considered: ilev=klev for 2D surface diagnostics
                                              !                           ilev=1 for full 3D diagnostics

    !--- 0) Initializations:

    zrmin(1:kproma,:,:)=3.0E-9_dp             ! 3nm cut-off radius for CN calculation

    zeps=EPSILON(1._dp)


    IF (nccndiag==1 .OR. nccndiag==3) THEN
       itoplev=klev
    ELSE
       itoplev=1

       IF (nccndiag==5 .OR. nccndiag==6) THEN
          !--- Calculate zdpg:
          zdpg(1:kproma,1)=2._dp*(paphm1(1:kproma,2)-papm1(1:kproma,1))/grav
          zdpg(1:kproma,2:klev)=(paphm1(1:kproma,3:klev+1)-paphm1(1:kproma,2:klev))/grav
       END IF
    END IF

    !--- Number per unit volume for each mode:

    zrho(1:kproma,itoplev:klev) = papm1(1:kproma,itoplev:klev)/(rd*ptvm1(1:kproma,itoplev:klev))

    DO jclass=1, nclass
       jt = sizeclass(jclass)%idt_no
       zn(1:kproma,itoplev:klev,jclass)=pxtm1(1:kproma,itoplev:klev,jt)*zrho(1:kproma,itoplev:klev)
    END DO

    !--- 1) Calculate properties for each aerosol mode:
    !--- 1.1) Calculate the auxiliary parameters A and B of the Koehler equation:
    !       Now done in ham_activ_koehler_ab once so that they can be used in
    !       convective and stratiform activation 


    !--- 2) Calculate the radius of activation [m], i.e. radius of smallest particle 
    !       with critical supersaturation larger than prescribed supersaturations [%]:
    !       AR et al. (2000) (Eq.4)

    DO jsat=1, nsat

       zcn(1:kproma,itoplev:klev)  = 0._dp
       zccn(1:kproma,itoplev:klev) = 0._dp

       ztmp=(2._dp/zsat(jsat))**(2._dp/3._dp)

       DO jclass=1, nclass
       
          IF (sizeclass(jclass)%lactivation) THEN

              zra(1:kproma,:,jclass) = 1.0_dp ! Set to large value of 1.0 m to set CCN count to zero

              WHERE(b(jclass)%ptr(1:kproma,itoplev:klev,krow)>zeps)
                 zra(1:kproma,itoplev:klev,jclass)                                          &
                   = a(jclass)%ptr(1:kproma,itoplev:klev,krow)                              &
                       / (3._dp * b(jclass)%ptr(1:kproma,itoplev:klev,krow)**(1._dp/3._dp)) &
                       * ztmp

              END WHERE

             !--- 3) Calculate the fractional number of each mode
             !       larger than the mode critical radius:

             CALL ham_m7_logtail(kproma,  kbdim,   klev,   krow,   jclass, &
                                 ll_numb, rdry(jclass)%ptr(:,:,krow),      &
                                 zra(:,:,jclass), zfracn(:,:,jclass)       )

             !--- 4) Sum up the total number of CCN [m-3]:

             zccn(1:kproma,itoplev:klev) = zccn(1:kproma,itoplev:klev)              &
                                             + zfracn(1:kproma,itoplev:klev,jclass) &
                                                 * zn(1:kproma,itoplev:klev,jclass)
         END IF
       END DO ! jclass

       !--- Store diagnostics in aero stream:

       IF (nccndiag==1 .OR. nccndiag==3 .OR. nccndiag==5) THEN

          ccn_2d(jsat)%ptr(1:kproma,krow)=zccn(1:kproma,klev)

       ELSE IF (nccndiag==2 .OR. nccndiag==4 .OR. nccndiag==6) THEN

          ccn_3d(jsat)%ptr(1:kproma,itoplev:klev,krow)=zccn(1:kproma,itoplev:klev)

       END IF

       IF (nccndiag==5 .OR. nccndiag==6) THEN

          DO jk=itoplev, klev
             ccn_burden(jsat)%ptr(1:kproma,krow)=ccn_burden(jsat)%ptr(1:kproma,krow) + &
                                                 zccn(1:kproma,jk)/zrho(1:kproma,jk)*zdpg(1:kproma,jk)
          END DO ! jk

       END IF

    END DO ! jsat

    IF (nccndiag > 2) THEN 

       DO jclass=1, nclass
          !--- 3) Calculate the fractional number of each mode
          !       larger than the minimum radius zrmin

          CALL ham_m7_logtail(kproma,  kbdim,   klev,   krow,   jclass, &
                              ll_numb, rdry(jclass)%ptr(:,:,krow),      &
                              zrmin(:,:,jclass), zfracn(:,:,jclass)     )

          !--- 4) Sum up the total number of CCN [m-3]:

          zcn(1:kproma,itoplev:klev) = zcn(1:kproma,itoplev:klev)               &
                                         + zfracn(1:kproma,itoplev:klev,jclass) &
                                             * zn(1:kproma,itoplev:klev,jclass)
       END DO ! jclass

       !--- Store diagnostics in aero stream:

       IF (nccndiag==3 .OR. nccndiag==5) THEN

          cn_2d(1:kproma,krow)=zcn(1:kproma,klev)

       ELSE IF (nccndiag==4 .OR. nccndiag==6) THEN

          cn_3d(1:kproma,itoplev:klev,krow)=zcn(1:kproma,itoplev:klev)

       END IF

       !--- 4) Sum up the CN over levels to burdens:

       IF (nccndiag==5 .OR. nccndiag==6) THEN

          DO jk=itoplev, klev
             cn_burden(1:kproma,krow)=cn_burden(1:kproma,krow) + &
                                      zcn(1:kproma,jk)/zrho(1:kproma,jk)*zdpg(1:kproma,jk)
          END DO ! jk

       END IF

    END IF

  END SUBROUTINE ham_ccn

END MODULE mo_ham_ccn
