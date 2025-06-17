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
SUBROUTINE sccd

  ! Description:
  !
  ! This subroutine computes the final value of divergence
  !
  ! Method:
  !
  ! *sccd* is called from *stepon*
  !
  ! Results:
  ! The implicit equation is inverted with the help of
  ! temporary array *zd*
  !
  ! Reference:
  ! See echam3 manual  eq. 2.5.44
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, August 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! T. Diehl, DKRZ, July 1999, parallel version 
  ! S. Borowski, HLRS/NEC, July 2006, improved vector length
  !
  ! for more details see file AUTHORS
  !
  
  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_memory_sp,     ONLY: sd
  USE mo_tmp_buffer,    ONLY: cn
  USE mo_control,       ONLY: nlev
  
  IMPLICIT NONE
  
  !  Local scalars: 
  INTEGER :: ic, is, snsp
  
  !  Local arrays: 
  REAL(dp) :: zd(nlev,2)
  INTEGER :: np1(lc%snsp)
  
  !  Executable statements 
  
  snsp = lc%snsp
  np1 = lc%np1
  
  !-- 1. Invert divergence equation

!$OMP PARALLEL PRIVATE(is,ic,zd)
!$OMP DO
!cdir novector
  DO  is = 1, snsp
    
    ic = np1(is)
        
!lk
!lk      zd(jl,1) = DOT_PRODUCT(cn(1:nlev,jl,ic),sd(1:nlev,1,is))
!lk      zd(jl,2) = DOT_PRODUCT(cn(1:nlev,jl,ic),sd(1:nlev,2,is))
!lk    END DO

    CALL dgemv('t',nlev,nlev,1.0_dp,cn(1,1,ic),nlev,sd(1,1,is),1,0.0_dp,zd(1,1),1)
    CALL dgemv('t',nlev,nlev,1.0_dp,cn(1,1,ic),nlev,sd(1,2,is),1,0.0_dp,zd(1,2),1)
    
    ! store implicit part of divergence equation
    
    sd(:,:,is) = zd(:,:)
    
  END DO
!$OMP END DO
!$OMP END PARALLEL
  
END SUBROUTINE sccd
