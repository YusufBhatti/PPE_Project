!>
!! @file xt_xmap_intersection_f.f90
!! @brief Fortran interface to yaxt xmap declarations
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
MODULE xt_xmap_intersection
  USE iso_c_binding, ONLY: c_int, c_loc, c_null_ptr, c_ptr
  USE xt_core, ONLY: xt_abort, xt_mpi_fint_kind, xt_slice_c_loc
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_f2c
  USE xt_xmap_abstract, ONLY: xt_xmap, xt_xmap_c2f
  IMPLICIT NONE
  PRIVATE

  TYPE, BIND(c), PUBLIC :: xt_com_list
    TYPE(xt_idxlist) :: list
    INTEGER(c_int) :: rank
  END TYPE xt_com_list

  PUBLIC :: xt_xmap_intersection_new, xt_xmap_intersection_ext_new

  INTERFACE
    FUNCTION xmi_new_f2c(num_src_intersections, src_com, &
         num_dst_intersections, dst_com, &
         src_idxlist, dst_idxlist, comm) RESULT(xmap) &
         BIND(c, name='xt_xmap_intersection_new_f2c')
      IMPORT :: c_int, c_ptr, xt_mpi_fint_kind
      INTEGER(c_int), VALUE, INTENT(in) :: num_src_intersections, &
           num_dst_intersections
      TYPE(c_ptr), VALUE, INTENT(in) :: src_com, dst_com, src_idxlist, &
           dst_idxlist
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
      TYPE(c_ptr) :: xmap
    END FUNCTION xmi_new_f2c

    FUNCTION xmi_ext_new_f2c(num_src_intersections, src_com, &
         num_dst_intersections, dst_com, &
         src_idxlist, dst_idxlist, comm) RESULT(xmap) &
         BIND(c, name='xt_xmap_intersection_ext_new_f2c')
      IMPORT :: c_int, c_ptr, xt_mpi_fint_kind
      INTEGER(c_int), VALUE, INTENT(in) :: num_src_intersections, &
           num_dst_intersections
      TYPE(c_ptr), VALUE, INTENT(in) :: src_com, dst_com, src_idxlist, &
           dst_idxlist
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
      TYPE(c_ptr) :: xmap
    END FUNCTION xmi_ext_new_f2c
  END INTERFACE

  INTERFACE xt_xmap_intersection_new
    MODULE PROCEDURE xmi_new_i_a_i_a
    MODULE PROCEDURE xmi_new_a_a
  END INTERFACE xt_xmap_intersection_new

  INTERFACE xt_xmap_intersection_ext_new
    MODULE PROCEDURE xmi_ext_new_i_a_i_a
    MODULE PROCEDURE xmi_ext_new_a_a
  END INTERFACE xt_xmap_intersection_ext_new

CONTAINS
  FUNCTION xmi_new_i_a_i_a(num_src_intersections, src_com, &
       num_dst_intersections, dst_com, &
       src_idxlist, dst_idxlist, comm) RESULT(xmap)
    INTEGER(c_int), VALUE, INTENT(in) :: num_src_intersections, &
         num_dst_intersections
    TYPE(xt_com_list), TARGET, INTENT(in) :: src_com(num_src_intersections), &
         dst_com(num_dst_intersections)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap
    TYPE(c_ptr) :: src_com_p, dst_com_p

    src_com_p = C_LOC(src_com)
    dst_com_p = C_LOC(dst_com)
    xmap = xt_xmap_c2f(xmi_new_f2c(&
         num_src_intersections, src_com_p, &
         num_dst_intersections, dst_com_p, &
         xt_idxlist_f2c(src_idxlist), xt_idxlist_f2c(dst_idxlist), comm))
  END FUNCTION xmi_new_i_a_i_a

  FUNCTION xmi_new_a_a(src_com, dst_com, src_idxlist, dst_idxlist, comm) &
       RESULT(xmap)
    TYPE(xt_com_list), TARGET, INTENT(in) :: src_com(:), dst_com(:)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    TYPE(xt_com_list), ALLOCATABLE, TARGET :: src_com_a(:), dst_com_a(:)
    TYPE(c_ptr) :: src_com_p, dst_com_p

    CALL com_p_arg(src_com, src_com_a, src_com_p)
    CALL com_p_arg(dst_com, dst_com_a, dst_com_p)

    xmap = xt_xmap_c2f(xmi_new_f2c(&
         INT(SIZE(src_com), c_int), src_com_p, &
         INT(SIZE(dst_com), c_int), dst_com_p, &
         xt_idxlist_f2c(src_idxlist), xt_idxlist_f2c(dst_idxlist), comm))
  END FUNCTION xmi_new_a_a

  FUNCTION xmi_ext_new_i_a_i_a(num_src_intersections, src_com, &
       num_dst_intersections, dst_com, &
       src_idxlist, dst_idxlist, comm) RESULT(xmap)
    INTEGER(c_int), VALUE, INTENT(in) :: num_src_intersections, &
         num_dst_intersections
    TYPE(xt_com_list), TARGET, INTENT(in) :: src_com(num_src_intersections), &
         dst_com(num_dst_intersections)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    xmap = xt_xmap_c2f(xmi_ext_new_f2c(&
         num_src_intersections, C_LOC(src_com), &
         num_dst_intersections, C_LOC(dst_com), &
         xt_idxlist_f2c(src_idxlist), xt_idxlist_f2c(dst_idxlist), comm))
  END FUNCTION xmi_ext_new_i_a_i_a

  FUNCTION xmi_ext_new_a_a(src_com, dst_com, src_idxlist, dst_idxlist, comm) &
       RESULT(xmap)
    TYPE(xt_com_list), TARGET, INTENT(in) :: src_com(:), dst_com(:)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    TYPE(xt_com_list), ALLOCATABLE, TARGET :: src_com_a(:), dst_com_a(:)
    TYPE(c_ptr) :: src_com_p, dst_com_p

    CALL com_p_arg(src_com, src_com_a, src_com_p)
    CALL com_p_arg(dst_com, dst_com_a, dst_com_p)

    xmap = xt_xmap_c2f(xmi_ext_new_f2c(&
         INT(SIZE(src_com), c_int), src_com_p, &
         INT(SIZE(dst_com), c_int), dst_com_p, &
         xt_idxlist_f2c(src_idxlist), xt_idxlist_f2c(dst_idxlist), comm))
  END FUNCTION xmi_ext_new_a_a

  SUBROUTINE com_p_arg(com, com_a, com_p)
    TYPE(xt_com_list), TARGET, INTENT(in) :: com(:)
    TYPE(xt_com_list), TARGET, ALLOCATABLE :: com_a(:)
    TYPE(c_ptr), INTENT(out) :: com_p

    INTEGER :: com_size
    INTERFACE
      FUNCTION xt_com_list_contiguous(com_a, com_b) RESULT(p) &
            BIND(c, name='xt_com_list_contiguous')
        IMPORT :: c_int, xt_com_list
        TYPE(xt_com_list), INTENT(in) :: com_a, com_b
        INTEGER(c_int) :: p
      END FUNCTION xt_com_list_contiguous
    END INTERFACE

    com_size = SIZE(com)
    IF (com_size > HUGE(1_c_int)) &
         CALL xt_abort('invalid size', &
         __FILE__, &
         __LINE__)
    IF (com_size > 0) THEN
      IF (com_size > 1) THEN
        IF (xt_com_list_contiguous(com(1), com(2)) /= 0) THEN
          CALL xt_slice_c_loc(com(1), com_p)
        ELSE
          ALLOCATE(com_a(com_size))
          com_a(:) = com(:)
          com_p = C_LOC(com_a)
        END IF
      ELSE
        CALL xt_slice_c_loc(com(1), com_p)
      END IF
    ELSE
      com_p = c_null_ptr
    END IF
  END SUBROUTINE com_p_arg

END MODULE xt_xmap_intersection
