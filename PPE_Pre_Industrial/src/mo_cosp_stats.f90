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
! Jul 2008 - A. Bodas-Salcedo - Added capability of producing outputs in standard grid
! Oct 2008 - J.-L. Dufresne   - Bug fixed. Assignment of Npoints,Nlevels,Nhydro,Ncolumns in COSP_STATS
! Oct 2008 - H. Chepfer       - Added PARASOL reflectance arguments
!
! 

MODULE mo_cosp_stats

  USE mo_kind,  ONLY: wp
  USE mo_cosp_constants
  USE mo_cosp_llnl_stats
  USE mo_cosp_lmd_ipsl_stats

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_STATS ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE cosp_stats(kbdim,Ncolumns,klev,           &
                      Llidar_sim, Llidar_cfad,       &
!gbx
                      ph, zlev, zlev_half, land,     &
!sglidar
                      beta_mol, beta_tot, refl,      &
                      srbval, cfad_sr, lidarcld,     &
                      cldlayer, parasolrefl          )
   ! Input arguments
   INTEGER, INTENT(IN) :: kbdim,Ncolumns,klev

   LOGICAL, INTENT(IN ) :: Llidar_sim, Llidar_cfad

!gbx
   REAL(wp), INTENT(INOUT) :: ph(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: zlev(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: zlev_half(kbdim,klev)
   REAL(wp), INTENT(INOUT) :: land(kbdim)

!sglidar
   REAL(wp) :: beta_mol(kbdim,klev) 
   REAL(wp) :: beta_tot(kbdim,Ncolumns,klev)
   REAL(wp) :: refl(kbdim,Ncolumns,PARASOL_NREFL)

! stlidar lidar stats
   REAL(wp), INTENT(INOUT) :: srbval(SR_BINS)
   REAL(wp), INTENT(INOUT) :: cfad_sr(kbdim,SR_BINS,Nlr) 
   REAL(wp), INTENT(INOUT) :: lidarcld(kbdim,Nlr)
   REAL(wp), INTENT(INOUT) :: cldlayer(kbdim,LIDAR_NCAT)
   REAL(wp), INTENT(INOUT) :: parasolrefl(kbdim,PARASOL_NREFL)

   ! Local variables 
   REAL(wp) :: betatot_out(kbdim,Ncolumns,Nlr)
   REAL(wp) :: betamol_in(kbdim,1,klev)
   REAL(wp) :: betamol_out(kbdim,1,Nlr)
   REAL(wp) :: betamol_c(kbdim,Nlr)
   REAL(wp) :: ph_in(kbdim,1,klev)
   REAL(wp) :: ph_out(kbdim,1,Nlr)
   REAL(wp) :: ph_c(kbdim,Nlr)

! vgrid
   REAL(wp) :: zl(Nlr), zu(Nlr)
   REAL(wp) :: zstep
   
   INTEGER :: i

   IF (use_vgrid) THEN ! Statistics in a different vertical grid

       ! cloudsat grid
        zstep = 480.0_wp
      DO i=1,Nlr
         zl(i) = REAL((i-1),wp)*zstep
         zu(i) =  REAL(i,wp)*zstep
      ENDDO

        betatot_out  = 0.0_wp
        betamol_out= 0.0_wp
        betamol_c  = 0.0_wp
        ph_in(:,1,:)  = ph(:,:)
        ph_out  = 0.0_wp
        ph_c    = 0.0_wp

        !++++++++++++ Lidar CFAD ++++++++++++++++
        IF (Llidar_sim) THEN
            betamol_in(:,1,:) = beta_mol(:,:)
            CALL cosp_change_vertical_grid(kbdim,1,klev,zlev,zlev_half,betamol_in, &
                                           Nlr,zl,zu,betamol_out)
            CALL cosp_change_vertical_grid(kbdim,Ncolumns,klev,zlev,zlev_half,beta_tot, &
                                           Nlr,zl,zu,betatot_out)
            CALL cosp_change_vertical_grid(kbdim,1,klev,zlev,zlev_half,ph_in, &
                                           Nlr,zl,zu,ph_out)
            ph_c(:,:) = ph_out(:,1,:)
            betamol_c(:,:) = betamol_out(:,1,:)
            ! Stats from lidar_stat_summary
            CALL diag_lidar(kbdim,Ncolumns,Nlr,SR_BINS,PARASOL_NREFL &
                            ,betatot_out,betamol_c,refl,land,ph_c &
                            ,LIDAR_UNDEF,Llidar_cfad &
                            ,cfad_sr,srbval &
                            ,LIDAR_NCAT,lidarcld,cldlayer,parasolrefl)
        ENDIF
    
   ELSE ! Statistics in model levels
    !++++++++++++ Lidar CFAD ++++++++++++++++
    ! Stats from lidar_stat_summary
    IF (Llidar_sim) CALL diag_lidar(kbdim,Ncolumns,Nlr,SR_BINS,PARASOL_NREFL &
                         ,beta_tot,beta_mol,refl,land,ph &
                         ,LIDAR_UNDEF,Llidar_sim & 
                         ,cfad_sr,srbval &
                         ,LIDAR_NCAT,lidarcld,cldlayer,parasolrefl)
    ENDIF
 ! Replace undef
   WHERE (cfad_sr   == LIDAR_UNDEF) cfad_sr   = R_UNDEF 
   WHERE (lidarcld  == LIDAR_UNDEF) lidarcld  = R_UNDEF 
   WHERE (cldlayer  == LIDAR_UNDEF) cldlayer  = R_UNDEF 
   WHERE (parasolrefl == LIDAR_UNDEF) parasolrefl = R_UNDEF 
   
END SUBROUTINE cosp_stats

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------- SUBROUTINE COSP_CHANGE_VERTICAL_GRID ----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE cosp_change_vertical_grid(Npoints,Ncolumns,Nlevels,zfull,zhalf,y,M,zl,zu,r,log_units)

   IMPLICIT NONE

   ! Input arguments
   INTEGER, INTENT(in) :: Npoints  !# of grid points
   INTEGER, INTENT(in) :: Nlevels  !# of levels
   INTEGER, INTENT(in) :: Ncolumns !# of columns
   REAL(wp), DIMENSION(Npoints,Nlevels), INTENT(in) :: zfull ! Height at model levels [m] (Bottom of model layer)
   REAL(wp), DIMENSION(Npoints,Nlevels), INTENT(in) :: zhalf ! Height at half model levels [m] (Bottom of model layer)
   REAL(wp), DIMENSION(Npoints,Ncolumns,Nlevels), INTENT(in) :: y     ! Variable to be changed to a different grid
   INTEGER, INTENT(in) :: M  !# levels in the new grid
   REAL(wp), DIMENSION(M), INTENT(in) :: zl ! Lower boundary of new levels  [m]
   REAL(wp), DIMENSION(M), INTENT(in) :: zu ! Upper boundary of new levels  [m]
   LOGICAL, OPTIONAL, INTENT(in) :: log_units ! log units, need to convert to linear units
   ! Output
   REAL(wp), DIMENSION(Npoints,Ncolumns,M), INTENT(out) :: r ! Variable on new grid

   ! Local variables
   INTEGER :: i,j,k
   LOGICAL :: lunits

   REAL(wp) :: ws
   REAL(wp), DIMENSION(Nlevels) :: xl,xu ! Lower and upper boundaries of model grid
   REAL(wp), DIMENSION(M) :: dz          ! Layer depth
   REAL(wp), DIMENSION(Nlevels,M) :: w   ! Weights to do the mean at each point
   REAL(wp), DIMENSION(Ncolumns,Nlevels) :: yp  ! Variable to be changed to a different grid.
                                           ! Local copy at a particular point.
                                           ! This allows for change of units.
   lunits=.false.
   IF (present(log_units)) lunits=log_units

   r = R_UNDEF
   DO i=1,Npoints
     ! Vertical grid at that point
     xl = zhalf(i,:)
     xu(1:Nlevels-1) = xl(2:Nlevels)
     xu(Nlevels) = zfull(i,Nlevels) +  zfull(i,Nlevels) - zhalf(i,Nlevels) ! Top level symmetric
     dz = zu - zl
     yp = y(i,:,:) ! Temporary variable to regrid
     ! Find weights
     w = 0.0_wp
     DO k=1,M
       DO j=1,Nlevels
         IF ((xl(j) < zl(k)).and.(xu(j) > zl(k)).and.(xu(j) <= zu(k))) THEN
           !xl(j)-----------------xu(j)
           !      zl(k)------------------------------zu(k)
           w(j,k) = xu(j) - zl(k)
         ELSEIF ((xl(j) >= zl(k)).and.(xu(j) <= zu(k))) THEN
           !           xl(j)-----------------xu(j)
           !      zl(k)------------------------------zu(k)
           w(j,k) = xu(j) - xl(j)
         ELSEIF ((xl(j) >= zl(k)).and.(xl(j) < zu(k)).and.(xu(j) >= zu(k))) THEN
           !                           xl(j)-----------------xu(j)
           !      zl(k)------------------------------zu(k)
           w(j,k) = zu(k) - xl(j)
         ELSEIF ((xl(j) <= zl(k)).and.(xu(j) >= zu(k))) THEN
           !  xl(j)---------------------------xu(j)
           !        zl(k)--------------zu(k)
           w(j,k) = dz(k)
         ENDIF
       ENDDO
     ENDDO
     ! Check for dBZ and change if necessary
     IF (lunits) THEN
        WHERE (yp /= R_UNDEF)
          yp = 10.0_wp**(yp/10.0_wp)
        ELSEWHERE
          yp = 0.0_wp
        ENDWHERE
     ENDIF
     ! Do the weighted mean
     DO j=1,Ncolumns
       DO k=1,M
          IF (zu(k) <= zhalf(i,1)) THEN ! Level below model bottom level
             r(i,j,k) = R_GROUND
          ELSE
            ws = sum(w(:,k))
            IF ((ws > 0.0_wp).and.(r(i,j,k) /= R_GROUND)) r(i,j,k) = sum(w(:,k)*yp(j,:))/ws
            ! Check for dBZ and change if necessary
            IF ((lunits).and.(r(i,j,k) /= R_GROUND)) THEN
                IF (r(i,j,k) <= 0.0_wp) THEN
                    r(i,j,k) = R_UNDEF
                ELSE
                    r(i,j,k) = 10.0_wp*log10(r(i,j,k))
                ENDIF
            ENDIF
          ENDIF
       ENDDO
     ENDDO
   ENDDO

END SUBROUTINE cosp_change_vertical_grid

END MODULE mo_cosp_stats
