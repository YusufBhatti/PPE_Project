!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_cosp_llnl_stats

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- FUNCTION COSP_CFAD ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION cosp_cfad(Npoints,Ncolumns,Nlevels,Nbins,x,xmin,xmax,bmin,bwidth)

   INTEGER, INTENT(in) :: Npoints,Ncolumns,Nlevels,Nbins
   real(dp), DIMENSION(Npoints,Nbins,Nlevels) :: cosp_cfad
   ! Input arguments
   REAL(dp), DIMENSION(Npoints,Ncolumns,Nlevels), INTENT(in) :: x
   REAL(dp), INTENT(in) :: xmin,xmax 
   REAL(dp), INTENT(in) :: bmin,bwidth
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

   cosp_cfad = 0.0_dp
   ! bwidth intervals in the range [bmin,bmax=bmin+Nbins*hwidth]
   ! Valid x values smaller than bmin and larger than bmax are set 
   ! into the smallest bin and largest bin, respectively.

   DO j = 1, Nlevels, 1
      DO k = 1, Ncolumns, 1
         DO i = 1, Npoints, 1 
            IF ((x(i,k,j) >= xmin) .and. (x(i,k,j) <= xmax)) THEN 
               ibin = ceiling((x(i,k,j) - bmin)/bwidth)
               IF (ibin > Nbins) ibin = Nbins
               IF (ibin < 1)     ibin = 1
               cosp_cfad(i,ibin,j) = cosp_cfad(i,ibin,j) + 1.0_dp 
            ENDIF
         ENDDO  !i
      ENDDO  !k
   ENDDO  !j

   cosp_cfad = cosp_cfad / REAL(Ncolumns,dp)

END FUNCTION cosp_cfad

END MODULE mo_cosp_llnl_stats
