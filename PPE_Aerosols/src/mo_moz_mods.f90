      module mo_moz_mods
      use mo_kind, only : dp
      use mo_exception, only : finish, message_text
!---------------------------------------------------------------------
! basic grid point resolution parameters
!---------------------------------------------------------------------
      implicit none
      save
      integer, parameter :: &
               pcnst = 242 +1, & ! number of advected constituents including cloud water
               pcnstm1 = 242 ! number of advected constituents excluding cloud water
!---------------------------------------------------------------------
! dimensions are set by ECHAM (moz_init)
!---------------------------------------------------------------------
      integer :: &
                plev, & ! number of vertical levels
                plevp, & ! plev plus 1
                plevm, & ! plev minus 1
                plon, & ! number of longitudes
                plat ! number of latitudes
      integer, parameter :: &
                pnats = 0 ! number of non-advected trace species
      integer :: nodes ! mpi task count
      integer :: plonl ! longitude tile dimension
      integer :: pplon ! longitude tile count
      integer :: plnplv ! plonl * plev
!--------------------------------------------------------------
! ... basic chemistry array parameters
!--------------------------------------------------------------
      integer, parameter :: hetcnt = 0, & ! number of heterogeneous processes
                            phtcnt = 142, & ! number of photo processes
                            rxntot = 733, & ! number of total reactions
                            gascnt = 591, & ! number of gas phase reactions
                            nfs = 2, & ! number of "fixed" species
                            relcnt = 0, & ! number of relationship species
                            grpcnt = 0, & ! number of group members
                            imp_nzcnt = 3482, & ! number of non-zero implicit matrix entries
                            rod_nzcnt = 0, & ! number of non-zero rodas matrix entries
                            extcnt = 0, & ! number of species with external forcing
                            clscnt1 = 0, & ! number of species in explicit class
                            clscnt2 = 0, & ! number of species in hov class
                            clscnt3 = 0, & ! number of species in ebi class
                            clscnt4 = 242, & ! number of species in implicit class
                            clscnt5 = 0, & ! number of species in rodas class
                            indexm = 1, & ! index of total atm density in invariant array
                            ncol_abs = 2, & ! number of column densities
                            indexh2o = 0, & ! index of water vapor density
                            clsze = 1 ! loop length for implicit chemistry
      integer :: ngrp = 0
      integer :: drydep_cnt = 0
      integer :: srfems_cnt = 0
      integer :: rxt_tag_cnt = 0
      integer :: fbc_cnt(2) = 0
      integer, allocatable :: grp_mem_cnt(:)
      integer, allocatable :: rxt_tag_map(:)
      real(dp) :: adv_mass(max(1,pcnstm1))
      real(dp) :: nadv_mass(max(1,grpcnt))
      real(dp), allocatable :: pht_alias_mult(:,:)
      character(len=32), allocatable :: rxt_tag_lst(:)
      character(len=16), allocatable :: pht_alias_lst(:,:)
      character(len=8), allocatable :: drydep_lst(:)
      character(len=8), allocatable :: srfems_lst(:)
      character(len=8), allocatable :: grp_lst(:)
      character(len=8), allocatable :: flbc_lst(:)
      character(len=8), allocatable :: fubc_lst(:)
      character(len=8) :: het_lst(max(1,hetcnt))
      character(len=8) :: extfrc_lst(max(1,extcnt))
      character(len=8) :: inv_lst(max(1,nfs))
      logical :: frc_from_dataset(max(1,extcnt))
      logical :: inv_from_dataset(max(1,nfs))
      type solver_class
         integer :: clscnt
         integer :: lin_rxt_cnt
         integer :: nln_rxt_cnt
         integer :: indprd_cnt
         integer :: iter_max
         integer :: cls_rxt_cnt(4)
         integer, pointer :: permute(:)
         integer, pointer :: diag_map(:)
         integer, pointer :: clsmap(:)
      end type solver_class
      type(solver_class) :: explicit, implicit, rodas
      contains
      subroutine chem_mods_inti
!--------------------------------------------------------------
! ... intialize the class derived type
!--------------------------------------------------------------
      implicit none
      integer :: astat
      explicit%clscnt = 0
      explicit%indprd_cnt = 0
      implicit%clscnt = 242
      implicit%lin_rxt_cnt = 175
      implicit%nln_rxt_cnt = 558
      implicit%indprd_cnt = 0
      implicit%iter_max = 11
      rodas%clscnt = 0
      rodas%lin_rxt_cnt = 0
      rodas%nln_rxt_cnt = 0
      rodas%indprd_cnt = 0
      if( explicit%clscnt > 0 ) then
         allocate( explicit%clsmap(explicit%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(message_text,*) ' failed to allocate explicit%clsmap ; error = ',astat
            call finish( 'chem_mods_inti:',message_text )
         end if
         explicit%clsmap(:) = 0
      end if
      if( implicit%clscnt > 0 ) then
         allocate( implicit%permute(implicit%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(message_text,*) ' failed to allocate implicit%permute ; error = ',astat
            call finish( 'chem_mods_inti:',message_text )
         end if
         implicit%permute(:) = 0
         allocate( implicit%diag_map(implicit%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(message_text,*) ' failed to allocate implicit%diag_map ; error = ',astat
            call finish( 'chem_mods_inti:',message_text )
         end if
         implicit%diag_map(:) = 0
         allocate( implicit%clsmap(implicit%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(message_text,*) ' failed to allocate implicit%clsmap ; error = ',astat
            call finish( 'chem_mods_inti:',message_text )
         end if
         implicit%clsmap(:) = 0
      end if
      if( rodas%clscnt > 0 ) then
         allocate( rodas%permute(rodas%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(message_text,*) ' failed to allocate rodas%permute ; error = ',astat
            call finish( 'chem_mods_inti:',message_text )
         end if
         rodas%permute(:) = 0
         allocate( rodas%diag_map(rodas%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(message_text,*) ' failed to allocate rodas%diag_map ; error = ',astat
            call finish( 'chem_mods_inti:',message_text )
         end if
         rodas%diag_map(:) = 0
         allocate( rodas%clsmap(rodas%clscnt),stat=astat )
         if( astat /= 0 ) then
            write(message_text,*) ' failed to allocate rodas%clsmap ; error = ',astat
            call finish( 'chem_mods_inti:',message_text )
         end if
         rodas%clsmap(:) = 0
      end if
      end subroutine chem_mods_inti
      end module mo_moz_mods
