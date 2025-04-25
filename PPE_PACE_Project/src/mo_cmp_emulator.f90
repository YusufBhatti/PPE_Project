MODULE mo_cmp_emulator
! Based on mo_active
! Author: UP
! Date: 07 2020

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream

  IMPLICIT NONE

  LOGICAL :: lemuphase_nic_cirrus ! Switch to phase the cirrus scheme in/out (between nic_cirrus = 1 and = 2)

  REAL(dp) :: eta_emu_nic_cirrus ! Phase cirrus scheme 1/2 in/out by this much, between 0 and 1


END MODULE mo_cmp_emulator
