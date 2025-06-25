!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_couple_wrap_externals
#ifdef __prism
  IMPLICIT NONE
  PRIVATE
!> Specify explicit interfaces of PRISM non-module subroutines
#ifndef __oa3mct
  INTERFACE
    SUBROUTINE prism_init_comp_proto(il_mynummod, cdnam, kinfo)
      USE mod_kinds_model
      CHARACTER*(*) cdnam
      INTEGER (kind=ip_intwp_p)     kinfo
      INTEGER (kind=ip_intwp_p)     il_mynummod
    END SUBROUTINE prism_init_comp_proto

    SUBROUTINE prism_get_localcomm_proto(il_local_comm, kinfo)
      USE mod_kinds_model
      INTEGER (kind=ip_intwp_p)     il_local_comm, kinfo
    END SUBROUTINE prism_get_localcomm_proto

    SUBROUTINE prism_def_var_proto(id_nports, cdport, id_part, &
          id_var_nodims, kinout, id_var_shape, ktype, kinfo)
      USE mod_kinds_model
      INTEGER (kind=ip_intwp_p)       kinout, ktype, kinfo, id_nports, &
                                      id_part
      INTEGER (kind=ip_intwp_p)       id_var_nodims(2), &
                                      id_var_shape(2*id_var_nodims(1))
      CHARACTER*(*) cdport
    END SUBROUTINE prism_def_var_proto

    SUBROUTINE prism_enddef_proto(kinfo)
      USE mod_kinds_model
      INTEGER (kind=ip_intwp_p) kinfo
    END SUBROUTINE prism_enddef_proto

    SUBROUTINE prism_abort_proto(id_compid, cd_routine, cd_message)
      INTEGER,          INTENT(in) :: id_compid
      CHARACTER(len=*), INTENT(in) :: cd_routine
      CHARACTER(len=*), INTENT(in) :: cd_message
    END SUBROUTINE prism_abort_proto

    SUBROUTINE prism_terminate_proto(kinfo)
      USE mod_kinds_model
      INTEGER (kind=ip_intwp_p) kinfo
    END SUBROUTINE prism_terminate_proto

    SUBROUTINE prism_get_freq(id_fldid, id_freq, id_info)
      USE mod_kinds_model
      INTEGER (kind=ip_intwp_p)       id_fldid, id_freq, id_info
    END SUBROUTINE prism_get_freq

    SUBROUTINE prism_put_restart_proto(id_port_id, kstep, kinfo)
      USE mod_kinds_model
      USE mod_prism_proto
      USE mod_comprism_proto
      INTEGER (kind=ip_intwp_p)      kstep, kinfo, id_port_id
    END SUBROUTINE prism_put_restart_proto
  END INTERFACE

  PUBLIC :: prism_init_comp_proto, prism_get_localcomm_proto, &
            prism_def_var_proto, prism_enddef_proto, &
            prism_abort_proto, prism_terminate_proto, prism_get_freq, &
            prism_put_restart_proto
#endif
#endif /*__prism*/
END MODULE mo_couple_wrap_externals

MODULE mo_couple_wrap
#ifdef __prism
#ifdef __oa3mct
  USE mod_oasis, ONLY : oasis_init_comp, oasis_get_localcomm, oasis_set_couplcomm, &
                        oasis_def_partition, oasis_def_var, oasis_enddef, &
                        oasis_start_grids_writing, oasis_write_grid, &
                        oasis_write_corner, oasis_write_mask, oasis_write_area, &
                        oasis_terminate_grids_writing, oasis_get, oasis_put, &
                        oasis_terminate, oasis_abort

  USE mod_oasis, ONLY : oasis_ok, oasis_Recvd, oasis_FromRest, &
                        oasis_RecvOut, oasis_FromRestOut, oasis_Sent, oasis_ToRest, &
                        oasis_SentOut, oasis_ToRestOut, oasis_Input, oasis_LocTrans, &
                        oasis_In, oasis_Out, oasis_Output, oasis_Real

  USE mod_oasis, ONLY : clim_strategy, clim_serial, clim_length, clim_offset, &
                        clim_orange, clim_segments, clim_sizex, clim_sizey, clim_ldx

  USE mod_oasis, ONLY : ip_realwp_p
  USE mod_oasis_namcouple, ONLY: &
                        namsrcfld, namdstfld, namfldseq, &
                        namruntim, namrstfil, namflddti, &
                        nnamcpl, namfldlag, prism_nmodels, prism_modnam
  USE mod_oasis_sys, ONLY: oasis_abort_noarg
#else
  USE mod_kinds_model, ONLY: &
                             ip_realwp_p

  USE mo_couple_wrap_externals, ONLY: &
      oasis_init_comp     => prism_init_comp_proto, &
      oasis_get_localcomm => prism_get_localcomm_proto, &
      oasis_def_var       => prism_def_var_proto, &
      oasis_enddef        => prism_enddef_proto, &
      oasis_terminate     => prism_terminate_proto, &
      oasis_abort         => prism_abort_proto, &
      oasis_put_restart   => prism_put_restart_proto

  USE mod_prism_proto, ONLY: &
      oasis_ok           => prism_Ok, &
      oasis_Recvd        => prism_Recvd, &
      oasis_FromRest     => prism_FromRest, &
      oasis_RecvOut      => prism_RecvOut, &
      oasis_FromRestOut  => prism_FromRestOut, &
      oasis_Sent         => prism_Sent, &
      oasis_ToRest       => prism_ToRest, &
      oasis_SentOut      => prism_SentOut, &
      oasis_ToRestOut    => prism_ToRestOut, &
      oasis_Input        => prism_Input, &
      oasis_LocTrans     => prism_LocTrans, &
      oasis_In           => prism_In, &
      oasis_Out          => prism_Out, &
      oasis_Output       => prism_Output, &
      oasis_Real         => prism_Real

  USE mod_prism_proto, ONLY: &
                        clim_strategy, clim_serial, clim_length, clim_offset, &
                        clim_orange, clim_segments, clim_sizex, clim_sizey, clim_ldx

  USE mod_comprism_proto, ONLY: &
      namsrcfld          => cg_cnaminp, &
      namdstfld          => cg_cnamout, &
      namfldseq          => ig_clim_seq, &
      !namfldseq         => ig_def_seq, &
      namruntim          => ig_ntime, &
      namrstfil          => cg_clim_rstfile, &
      !namrstfil         => cg_def_rstfile, &
      namflddti          => ig_def_freq, &
      nnamcpl            => ig_clim_nfield, &
      namfldlag          => ig_def_lag, &
      ig_inidate, lg_ncdfrst

  USE mod_prism_get_proto, ONLY: &
      oasis_get          => prism_get_proto

  USE mod_prism_put_proto, ONLY: &
      oasis_put          => prism_put_proto

  USE mod_prism_def_partition_proto, ONLY: &
      oasis_def_partition           => prism_def_partition_proto

  USE mod_prism_grids_writing, ONLY: &
      oasis_start_grids_writing     => prism_start_grids_writing, &
      oasis_write_grid              => prism_write_grid, &
      oasis_write_corner            => prism_write_corner, &
      oasis_write_mask              => prism_write_mask, &
      oasis_write_area              => prism_write_area, &
      oasis_terminate_grids_writing => prism_terminate_grids_writing
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: oasis_init_comp, oasis_get_localcomm, &
            oasis_def_partition, oasis_def_var, oasis_enddef, &
            oasis_start_grids_writing, oasis_write_grid, &
            oasis_write_corner, oasis_write_mask, oasis_write_area, &
            oasis_terminate_grids_writing, oasis_get, oasis_put, &
            oasis_terminate, oasis_abort
  PUBLIC :: oasis_ok, oasis_Recvd, oasis_FromRest, &
            oasis_RecvOut, oasis_FromRestOut, oasis_Sent, oasis_ToRest, &
            oasis_SentOut, oasis_ToRestOut, oasis_Input, oasis_LocTrans, &
            oasis_In, oasis_Out, oasis_Output, oasis_Real
  PUBLIC :: ip_realwp_p
  PUBLIC :: clim_strategy, clim_serial, clim_length, clim_offset, &
            clim_orange, clim_segments, clim_sizex, clim_sizey, clim_ldx
  PUBLIC :: namsrcfld, namdstfld, namfldseq, namruntim, namrstfil, &
            namflddti, nnamcpl, namfldlag
#ifdef __oa3mct
  PUBLIC :: prism_nmodels, prism_modnam, oasis_set_couplcomm
#else
  PUBLIC :: ig_inidate, lg_ncdfrst, oasis_put_restart
#endif

#endif /*__prism*/
END MODULE mo_couple_wrap

