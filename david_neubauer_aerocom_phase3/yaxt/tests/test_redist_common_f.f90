!>
!! @file test_redist_common_f.f90
!! @brief common routines for Fortran test of redist classes
!!
!! @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://doc.redmine.dkrz.de/yaxt/html/
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
MODULE test_redist_common
  USE mpi
  USE yaxt, ONLY: xt_idxlist, xt_int_kind, xt_idxvec_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_all2all_new
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: build_odd_selection_xmap
CONTAINS
  ! build xmap for destination list containing all odd elements of
  ! source list dimensioned 1 to src_slice_len
  FUNCTION build_odd_selection_xmap(src_slice_len) RESULT(xmap)
    INTEGER, INTENT(in) :: src_slice_len
    TYPE(xt_xmap) :: xmap
    INTEGER :: i, j, dst_slice_len
    INTEGER, PARAMETER :: dst_step = 2
    INTEGER(xt_int_kind), ALLOCATABLE :: index_list(:)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist

    dst_slice_len = (src_slice_len + dst_step - 1)/dst_step
    ALLOCATE(index_list(src_slice_len))
    DO i = 1, src_slice_len
      index_list(i) = INT(i, xt_int_kind)
    END DO
    src_idxlist = xt_idxvec_new(index_list)
    j = 1
    DO i = 1, src_slice_len, dst_step
      index_list(j) = INT(i, xt_int_kind)
      j = j + 1
    END DO
    dst_idxlist = xt_idxvec_new(index_list, dst_slice_len)
    DEALLOCATE(index_list)

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world);
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END FUNCTION build_odd_selection_xmap

END MODULE test_redist_common
