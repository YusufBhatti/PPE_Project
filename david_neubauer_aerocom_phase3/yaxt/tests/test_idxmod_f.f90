!>
!! @file test_idxmod.c
!!
!! @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
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
PROGRAM test_idxmod_f
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_idxlist, xt_idxempty_new, &
       xt_int_kind, xi => xt_int_kind, xt_stripe, &
       xt_modifier, xt_idxvec_new, xt_idxstripes_new, &
       xt_idxlist_delete, xt_idxmod_new
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count
  IMPLICIT NONE

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL test_idxvec_modifier
  CALL test_idxstripes_modifier
  CALL test_multimod

  CALL xt_finalize
  CALL finish_mpi

  IF (test_err_count() /= 0) CALL test_abort("non-zero error count", &
       __FILE__, &
       __LINE__)

CONTAINS
  SUBROUTINE test_idxvec_modifier
    INTEGER, PARAMETER :: g_src_num = 9, g_dst_num=9, patch_num=7
    INTEGER(xi) :: i
    INTEGER(xi), PARAMETER :: g_src_idx(g_src_num) = (/ (i, i=1,g_src_num) /), &
         g_dst_idx(g_dst_num) = (/ (i, i=g_dst_num,1,-1) /), &
         patch_idx(patch_num) = (/ 3_xi, 4_xi, 4_xi, 4_xi, 7_xi, 7_xi, 8_xi /)
    ! idx:{3,4,4,4,7,7,8} -> pos:{2,3,3,3,6,6,7} => idx:{7,6,6,6,3,3,2}
    INTEGER(xi), PARAMETER :: &
         ref_mpatch_idx(patch_num) = (/ 7_xi, 6_xi, 6_xi, 6_xi, 3_xi, 3_xi, &
         2_xi /)

    TYPE(xt_idxlist) :: g_src_idxlist,  g_dst_idxlist, patch_idxlist, &
         mpatch_idxlist
    TYPE(xt_modifier) :: modifier(1)

    g_src_idxlist = xt_idxvec_new(g_src_idx, g_src_num)

    g_dst_idxlist = xt_idxvec_new(g_dst_idx, g_dst_num)

    modifier(1) = xt_modifier(g_src_idxlist, g_dst_idxlist, 0)

    patch_idxlist = xt_idxvec_new(patch_idx, patch_num)

    mpatch_idxlist = xt_idxmod_new(patch_idxlist, modifier)

    CALL check_idxlist(mpatch_idxlist, ref_mpatch_idx)

    CALL xt_idxlist_delete(mpatch_idxlist)
    CALL xt_idxlist_delete(patch_idxlist)
    CALL xt_idxlist_delete(g_dst_idxlist)
    CALL xt_idxlist_delete(g_src_idxlist)
  END SUBROUTINE test_idxvec_modifier

  SUBROUTINE test_idxstripes_modifier
    INTEGER :: i
    INTEGER, PARAMETER :: patch_num = 8
    TYPE(xt_stripe), PARAMETER :: &
         g_src_stripe = xt_stripe(1, 1, 20), &
         g_dst_stripe = xt_stripe(100, -1, 20)
    TYPE(xt_idxlist) :: g_src_idxlist, g_dst_idxlist, patch_idxlist, &
         mpatch_idxlist
    TYPE(xt_modifier) :: modifier(1)
    INTEGER :: mstate(patch_num)
    INTEGER, PARAMETER :: ref_mstate(patch_num) &
         = (/ 1, IOR(2, 32), IOR(3, 32), IOR(4, 32), IOR(5, 32), 6, 7, 8 /)
    ! inter:{1,3,3,5} => extract_pos:{0,2,2,4} => subst_idx:{100,98,98,96},
    ! patch_pos:{1,2,3,4} = > mpatch:{0,100,98,98,96,50,100,150}
    INTEGER(xi), PARAMETER :: &
         patch_idx(patch_num) = (/   0_xi,   1_xi,   3_xi,   3_xi, &
         &                           5_xi,  50_xi, 100_xi, 150_xi /), &
         ref_mpatch_idx(patch_num) = (/   0, 100,  98,  98, &
         &                               96,  50, 100, 150 /)

    g_src_idxlist = xt_idxstripes_new(g_src_stripe)
    g_dst_idxlist = xt_idxstripes_new(g_dst_stripe)

    modifier(1) = xt_modifier(g_src_idxlist, g_dst_idxlist, 32)

    patch_idxlist = xt_idxvec_new(patch_idx)

    DO i = 1, 8
      mstate(i) = i
    END DO

    mpatch_idxlist = xt_idxmod_new(patch_idxlist, modifier, mstate)

    CALL check_idxlist(mpatch_idxlist, ref_mpatch_idx)

    ! check mstate
    IF (ANY(mstate(:) /= ref_mstate(:))) &
         CALL test_abort("mstate(:) /= ref_mstate(:)", &
         __FILE__, &
         __LINE__)
    CALL xt_idxlist_delete(mpatch_idxlist)
    CALL xt_idxlist_delete(patch_idxlist)
    CALL xt_idxlist_delete(g_dst_idxlist)
    CALL xt_idxlist_delete(g_src_idxlist)
  END SUBROUTINE test_idxstripes_modifier

  ! track modifier usage
  SUBROUTINE test_multimod
    INTEGER, PARAMETER :: g1_src_num = 9, g1_dst_num = 9, &
         g2_src_num = 5, g2_dst_num = 5, patch_num = 6
    INTEGER(xi) :: i
    INTEGER(xi), PARAMETER :: &
         g1_src_idx(g1_src_num) = (/ (i, i = 1_xi,g1_src_num) /), &
         g1_dst_idx(g1_dst_num) = (/ (i, i = g1_dst_num,1_xi,-1_xi) /), &
         g2_src_idx(g2_src_num) = (/ 1_xi, 2_xi, 8_xi, 9_xi, 10_xi /), &
         g2_dst_idx(g2_dst_num) = (/ 8_xi, 2_xi, 8_xi, 2_xi, 5_xi /), &
         patch_idx(patch_num) = (/ 6_xi, 7_xi, 25_xi, 8_xi, 9_xi, 10_xi /), &
         ! mod1: idx:{6,7,25,8,9,10} -> pos:{5,6,nil,7,8,nil} => idx:{4,3,25,2,1,10}
         ! mod2: idx:{4,3,25,2,1,10} -> pos:{nil,nil,nil,1,0,4}
         !                           => idx:{4,3,25,2,8,5}
         ref_mpatch_idx(patch_num) = (/ 4_xi, 3_xi, 25_xi, 2_xi, 8_xi, 5_xi /)

    TYPE(xt_idxlist) :: g1_src_idxlist, g1_dst_idxlist, g2_src_idxlist, &
         g2_dst_idxlist, patch_idxlist, mpatch_idxlist
    TYPE(xt_modifier) ::  modifier(2)
    INTEGER :: mstate(patch_num)
    INTEGER, PARAMETER :: ref_mstate(patch_num) &
         = (/ IOR(1, 0), IOR(1, 0), IOR(0, 0), IOR(1,2), IOR(1,2), IOR(0,2) /)

    g1_src_idxlist = xt_idxvec_new(g1_src_idx, g1_src_num)
    g1_dst_idxlist = xt_idxvec_new(g1_dst_idx)
    g2_src_idxlist = xt_idxvec_new(g2_src_idx)
    g2_dst_idxlist = xt_idxvec_new(g2_dst_idx, g2_dst_num)

    patch_idxlist = xt_idxvec_new(patch_idx)
    ! reset mstate
    mstate(:) = 0

    modifier(1) = xt_modifier(g1_src_idxlist, g1_dst_idxlist, 1)
    modifier(2) = xt_modifier(g2_src_idxlist, g2_dst_idxlist, 2)
    mpatch_idxlist = xt_idxmod_new(patch_idxlist, modifier, 2, mstate)
    CALL check_idxlist(mpatch_idxlist, ref_mpatch_idx)

    ! check mstate
    IF (ANY(mstate(:) /= ref_mstate(:))) &
         CALL test_abort("mstate(:) /= ref_mstate(:)", &
         __FILE__, &
         __LINE__)
    CALL xt_idxlist_delete(mpatch_idxlist)
    CALL xt_idxlist_delete(patch_idxlist)
    CALL xt_idxlist_delete(g2_dst_idxlist)
    CALL xt_idxlist_delete(g2_src_idxlist)
    CALL xt_idxlist_delete(g1_dst_idxlist)
    CALL xt_idxlist_delete(g1_src_idxlist)
  END SUBROUTINE test_multimod

END PROGRAM test_idxmod_f
