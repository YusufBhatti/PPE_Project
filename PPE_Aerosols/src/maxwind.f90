#ifdef __xlC__
@PROCESS STRICT
#endif
!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
SUBROUTINE maxwind

  ! Description:
  !
  ! Computes maximum winds for horizontal diffusion
  ! and diagnostics.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, April 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, December 2000, time control
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_semi_impl,     ONLY: nulev, uvmax, vcheck, vmax 
  USE mo_time_control,  ONLY: get_time_step
  USE mo_exception,     ONLY: message, message_text

  IMPLICIT NONE

  ! Local arrays
  INTEGER :: nmaxlev(1)   

  !  Local scalars: 
  INTEGER ::  istep

  !  Executable statements 

  !-- 1. Maximum level winds for diffusion

!  vmax(:) = maxval_latit (vmaxz(:,:))

  istep = get_time_step()

  !-- 2. Maximum wind for diagnostics

  nmaxlev = MAXLOC(vmax)                    ! vertical index of max
  nulev   = nmaxlev(1)                      ! for postatd
  uvmax   = vmax(nulev)                     ! max value (absolute)

  !-- 3. Check for high windspeeds

  IF (uvmax > vcheck) THEN
    WRITE (message_text,'(a,f4.0,a)') ' WARNING! high wind speed: ',uvmax,' m/s'
    CALL message ('',message_text)
    WRITE (message_text,'(a,i4,a,i8)') ' Level: ', nulev, ' NSTEP= ', istep
    CALL message ('',message_text)
  END IF

END SUBROUTINE maxwind

SUBROUTINE zonal_winds_fs

  ! Description:
  !
  ! Computes zonal and latitude maximum winds for horizontal diffusion and diagnostics.
  !
  ! Authors:
  !
  ! M. Puetz, IBM, November 2010, all computations in FS space
  ! 
  ! for more details see file AUTHORS
  !
  
  USE mo_kind,          ONLY: dp
  USE mo_buffer_fft,    ONLY: fvmax,fulz,fu,fv
  USE mo_decomposition, ONLY: dc => local_decomposition
  USE mo_gaussgrid,     ONLY: gl_rsqcst
  USE mo_global_op,     ONLY: maxval_latit
  
  IMPLICIT NONE

  REAL(dp) :: vmaxz(dc%nflev,dc%nflat)
  INTEGER ::  nlon, nlat, nlev, jlat, jglat, k

  !-- Compute maximum !u!+!v!                                                                                                        
  
  nlat = dc%nflat
  nlev = dc%nflev
  nlon = dc%nlon

!$OMP PARALLEL DO PRIVATE(jlat,k)
  DO jlat = 1, nlat
    jglat = dc% glat(jlat)          ! global continuous north -> south
    ! use FS mapping to global latitude index
    ! jglat = dc%glat_fs(jlat)
    DO k = 1,nlev
      ! zonal max/min over longitudes
      fulz(k,jlat)   = 0.5_dp * (MAXVAL(fu(1:nlon,k,jlat)) + MINVAL(fu(1:nlon,k,jlat)))
      vmaxz(k,jlat) = gl_rsqcst(jglat) * SQRT(MAXVAL(fu(1:nlon,k,jlat)**2 + fv(1:nlon,k,jlat)**2))
    END DO
  END DO
!$OMP END PARALLEL DO

  !- Maximum level winds for diffusion

  fvmax(:) = maxval_latit (vmaxz(:,:))
  
  ! NOTE: fvmax and fulz will be transposed to vmax and ulz in GP space within fourier_to_gridspace()

END SUBROUTINE zonal_winds_fs

SUBROUTINE zonal_winds_gp
  
  ! Description:
  !
  ! Computes zonal and latitude maximum winds for horizontal diffusion and diagnostics.
  !
  ! Authors:
  !
  ! M. Puetz, IBM, November 2010, all computations in FS space
  ! 
  ! for more details see file AUTHORS
  !
  
  USE mo_kind,          ONLY: dp
  USE mo_scan_buffer,   ONLY: ul, ulz, u, v
  USE mo_semi_impl,     ONLY: vmax
  USE mo_decomposition, ONLY: dc => local_decomposition
  USE mo_transpose,     ONLY: reorder
  USE mo_gaussgrid,     ONLY: gl_rsqcst
  USE mo_global_op,     ONLY: maxval_latit,maxval_zonal,minval_zonal
  
  IMPLICIT NONE
  
  REAL(dp) :: zru(dc%nglon,dc%nlev,dc%nglat)
  REAL(dp) :: zrv(dc%nglon,dc%nlev,dc%nglat)
  REAL(dp) :: vmaxz(dc%nlev,dc%nglat)
  REAL(dp) :: zrcst(dc%nglat)
  INTEGER ::  nglat, nlev, jrow, jglat

  nlev  = dc%nlev
  nglat = dc%nglat

  !-- Compute maximum |u|+|v|

  CALL reorder(zru, u)
  CALL reorder(zrv, v)
  
  DO jrow = 1, nglat
    jglat = dc% glat(jrow)          ! global continuous north -> south
    zrcst(jrow) = gl_rsqcst(jglat)
  END DO

  ! NOTE: MPI collective operations -> do not OMP parallelize
  
  ulz(:,:) = (maxval_zonal(zru(:,:,:))  + &
              minval_zonal(zru(:,:,:))) * 0.5_dp

  ul(:,:) = ulz(:,:)

  vmaxz(:,:) = SQRT(maxval_zonal(zru(:,:,:)*zru(:,:,:)                 &
                               + zrv(:,:,:)*zrv(:,:,:)))               &
                 * SPREAD(zrcst,dim=1,ncopies=nlev)

  !-- Maximum level winds for diffusion

  vmax(:) = maxval_latit (vmaxz(:,:))
  
END SUBROUTINE zonal_winds_gp
