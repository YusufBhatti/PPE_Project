!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!!  mo_moz_imp_sol
!!
!! \brief
!!  This module implements the implicit backward Euler chemistry solver of MOZART
!!
!! \author Stacy Walters (NCAR)
!!
!! \responsible_coder
!!  Martin Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!  - SW: original version (pre 2004)
!!
!! \limitations
!!  none
!!
!! \details
!!
!! \bibliographic_references
!!  none
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!!  Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!!  licencing agreement to be found at:
!!  https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!!  The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mo_moz_imp_sol

      use mo_kind, only : dp
      use mo_exception, only: finish, message_text
      use mo_moz_mods, only : clscnt4

      implicit none

      private
      public :: imp_slv_inti
      public :: imp_sol

      save

      real(dp), parameter    ::  rel_err      = 1.e-3
      real(dp), parameter    ::  high_rel_err = 1.e-4
!-----------------------------------------------------------------------      
!   	    newton-raphson iteration limit, spatial loop lengths
!-----------------------------------------------------------------------      
! if doing a kproma-test: set cut_limit=0
      integer, parameter :: cut_limit    = 5
!      integer, parameter :: vec_len      = 1

      integer :: ndx_o3, ndx_no, ndx_no2, ndx_no3, &
                 ndx_hno3, ndx_ho2no2, ndx_n2o5, &
                 ndx_oh, ndx_ho2
!     logical :: do_ox_pl = .true.

      real(dp), private                      ::   small
      real(dp), private                      ::   epsilon(clscnt4)
!     type(hst_pl), private, allocatable ::   imp_hst_prod(:)
!     type(hst_pl), private, allocatable ::   imp_hst_loss(:)
      logical, private, allocatable      ::   factor(:)

      contains

      subroutine imp_slv_inti
!-----------------------------------------------------------------------      
!	... initialize the implict solver
!-----------------------------------------------------------------------      

      use mo_moz_mods,    only : clscnt4, implicit
      use mo_moz_mods,    only : pcnstm1
      use mo_moz_util,    only : get_spc_ndx

      implicit none

!-----------------------------------------------------------------------      
!	... local variables
!-----------------------------------------------------------------------      
      integer :: astat
      integer :: il, iu, m
      integer :: wrk(21)
      real(dp)    :: eps(pcnstm1)

      allocate( factor(implicit%iter_max),stat=astat )
      if( astat /= 0 ) then
         write(message_text,'(a,i0)') 'failed to allocate factor array; error = ',astat
         CALL finish('imp_slv_inti', message_text)
      end if

      factor(:) = .true.
      eps(:)    = rel_err
! set stricter convergence criteria for key species
      ndx_o3 = get_spc_ndx('O3')
      if( ndx_o3 > 0 ) then
         eps(ndx_o3) = high_rel_err
      end if
      ndx_no = get_spc_ndx('NO')
      if( ndx_no > 0 ) then
         eps(ndx_no) = high_rel_err
      end if
      ndx_no2 = get_spc_ndx('NO2')
      if( ndx_no2 > 0 ) then
         eps(ndx_no2) = high_rel_err
      end if
      ndx_no3 = get_spc_ndx('NO3')
      if( ndx_no3 > 0 ) then
         eps(ndx_no3) = high_rel_err
      end if
      ndx_hno3 = get_spc_ndx('HNO3')
      if( ndx_hno3 > 0 ) then
         eps(ndx_hno3) = high_rel_err
      end if
      ndx_ho2no2 = get_spc_ndx('HO2NO2')
      if( ndx_ho2no2 > 0 ) then
         eps(ndx_ho2no2) = high_rel_err
      end if
      ndx_n2o5 = get_spc_ndx('N2O5')
      if( ndx_n2o5 > 0 ) then
         eps(ndx_n2o5) = high_rel_err
      end if
      ndx_oh = get_spc_ndx('OH')
      if( ndx_oh > 0 ) then
         eps(ndx_oh) = high_rel_err
      end if
      ndx_ho2 = get_spc_ndx('HO2')
      if( ndx_ho2 > 0 ) then
         eps(ndx_ho2) = high_rel_err
      end if
! enter eps values and respect permutation (clsmap)
      do m = 1,max(1,clscnt4)
         epsilon(m) = eps(implicit%clsmap(m))
      end do

!++mgs : removed do_ox_pl stuff
!++mgs : removed history stuff

      small = 1.e-40_dp 

      end subroutine imp_slv_inti

!++mgs : removed subroutine imp_slv_lt_inti (= history stuff)

      subroutine imp_sol( base_sol, reaction_rates, het_rates, extfrc, &
                          delt, lat, plonl, plnplv )
!-----------------------------------------------------------------------
!      	... imp_sol advances the volumetric mixing ratio
!           forward one time step via the fully implicit euler scheme.
!-----------------------------------------------------------------------

      use mo_moz_mods,       only : clscnt4, imp_nzcnt, clsze, rxntot, &
                                    hetcnt, extcnt, implicit
      use mo_moz,            only : tracnam 
!     use mo_control,        only : t_pdiag, pdiags
      use mo_moz_mods,       only : plev, pcnstm1
      use mo_moz_mat,        only : indprd
      use mo_moz_mat,        only : imp_prod_loss
      use mo_moz_mat,        only : imp_linmat
      use mo_moz_mat,        only : imp_nlnmat
      use mo_moz_mat,        only : imp_lu_fac
      use mo_moz_mat,        only : imp_lu_slv
      use mo_exception,      only: message, message_text, em_warn
!+++sschr: for testing
      use mo_mpi, only: p_parallel_io, p_pe
!---sschr: for testing
! to be put to diagnostic stream (at the moment in debug stream)
!     use mo_debugs, only: ddf1, ddf2, ddf3
!     if commented out: ddf1: Convergence of data (1.0: all data converged at this location,
!                                                  0.0: some data didn't converge at this location)
!                       ddf2: smallest dt at this location
!                       ddf3: number of iterations at this location
      implicit none

!-----------------------------------------------------------------------
!     	... dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::   lat                       ! lat index
      integer, intent(in) ::   plonl                     ! longitude tile dimension
      integer, intent(in) ::   plnplv                    ! plonl*plev
      real(dp), intent(in)    ::   delt                      ! time step (seconds)
! note: dims are identical in scalar and vec version
      real(dp), intent(in)    ::   reaction_rates(plnplv,rxntot)
      real(dp), intent(in)    ::   het_rates(plnplv,max(1,hetcnt))
      real(dp), intent(in)    ::   extfrc(plnplv,extcnt)
      real(dp), intent(inout) ::   base_sol(plnplv,pcnstm1)

!-----------------------------------------------------------------------
!     	... local variables
!-----------------------------------------------------------------------
      integer, parameter :: indx_debug = 311
      integer ::   nr_iter
      integer ::   indx, &
                   i, j, k, l, m, n
      integer ::  base_ndx
      integer ::  cls_ndx
!      integer ::  ofl, ofu
!      integer ::  glon
      integer ::  lev
!      integer ::  max_loc(1)
!      integer ::  max_ind(max(1,clscnt4))
      integer ::  cut_cnt
      integer ::  stp_con_cnt
      real(dp)    ::  interval_done 
      real(dp)    ::  dt
      real(dp)    ::  dti
      real(dp)    ::  wrk
      real(dp)    ::  max_delta(max(1,clscnt4))
      real(dp)    ::  sys_jac(max(1,imp_nzcnt))
      real(dp)    ::  lin_jac(max(1,imp_nzcnt))
      real(dp)    ::  solution(max(1,clscnt4))
      real(dp)    ::  forcing(max(1,clscnt4))
      real(dp)    ::  iter_invariant(max(1,clscnt4))
      real(dp)    ::  ind_prd(plnplv, max(1,clscnt4))
      real(dp)    ::  prod(max(1,clscnt4))
      real(dp)    ::  loss(max(1,clscnt4))
      real(dp)    ::  lrxt(max(1,rxntot))
      real(dp)    ::  lsol(max(1,pcnstm1))
      real(dp)    ::  lhet(max(1,hetcnt))

      real(dp)    ::  wrk_buff(plnplv)
!      real, allocatable ::  base_sol_sav(:,:)
      logical ::  convergence
      logical ::  dump_buff
      logical ::  fill_buff
      logical ::  frc_mask
      logical ::  iter_conv
!      logical ::  non_conv(plnplv)
!      logical ::  wrk_mask(plnplv)
!      logical ::  time_looping(plnplv)
      logical ::  converged(max(1,clscnt4))

!++mgs: removed history file stuff

      if( implicit%indprd_cnt > 0 .or. extcnt > 0 ) then
!-----------------------------------------------------------------------      
!        ... class independent forcing
!-----------------------------------------------------------------------      
         call indprd( 4, ind_prd, base_sol, extfrc, reaction_rates )
      else
         do m = 1,max(1,clscnt4)
            ind_prd(:,m) = 0._dp
         end do
      end if
! Note: zeroing done at different places in vector and cache version!

!-----------------------------------------------------------------------      
! spatial loop: lev and lon in scalar, 1 to plnplv in vector version
!-----------------------------------------------------------------------      
level_loop : &
      do lev = 1,plev
lon_tile_loop : &
         do i = 1,plonl
            indx = (lev - 1)*plonl + i
!-----------------------------------------------------------------------      
!       ... transfer from base to local work arrays
!-----------------------------------------------------------------------      
         do m = 1,rxntot
            lrxt(m) = reaction_rates(indx,m) 
         end do
         if( hetcnt > 0 ) then
            do m = 1,hetcnt
               lhet(m) = het_rates(indx,m) 
            end do
         end if
         dt            = delt
         cut_cnt       = 0
         stp_con_cnt   = 0
         interval_done = 0._dp

!-----------------------------------------------------------------------      
!        ... time step loop
!-----------------------------------------------------------------------      
time_step_loop : &
         do
            dti = 1._dp / dt
!-----------------------------------------------------------------------      
!        ... transfer from base to local work arrays
!-----------------------------------------------------------------------      
            do m = 1,pcnstm1
               lsol(m) = base_sol(indx,m) 
            end do
!-----------------------------------------------------------------------      
!        ... transfer from base to class array
!-----------------------------------------------------------------------      
            do k = 1,clscnt4
               j = implicit%clsmap(k)
               m = implicit%permute(k)
               solution(m) = lsol(j)
            end do

!-----------------------------------------------------------------------      
!       ... set the iteration invariant part of the function Ff(y)
!       ... (if there is "independent" production put it in the forcing)
!-----------------------------------------------------------------------      
            if( implicit%indprd_cnt > 0 .or. extcnt > 0 ) then
               do m = 1,clscnt4
                  iter_invariant(m) = dti * solution(m) + ind_prd(indx,m)
               end do
            else
               do m = 1,clscnt4
                  iter_invariant(m) = dti * solution(m)
               end do
            end if

!-----------------------------------------------------------------------      
!        ... the linear component
!-----------------------------------------------------------------------      
            do m = 1,imp_nzcnt
               lin_jac(m) = 0._dp
            end do
            call imp_linmat( lin_jac, lsol, lrxt, lhet )

!=======================================================================
!        the newton-raphson iteration for f(y) = 0
!=======================================================================
iter_loop : &
            do nr_iter = 1,implicit%iter_max
!-----------------------------------------------------------------------      
!       ... the non-linear component
!-----------------------------------------------------------------------      
               if( factor(nr_iter) ) then
                  do m = 1,imp_nzcnt
                     sys_jac(m) = 0._dp
                  end do
                  call imp_nlnmat( sys_jac, lsol, lrxt, lin_jac, dti )
!-----------------------------------------------------------------------      
!        ... factor the "system" matrix
!-----------------------------------------------------------------------      
                  call imp_lu_fac( sys_jac )
               end if      
!-----------------------------------------------------------------------      
!        ... form f(y)
!-----------------------------------------------------------------------      
               call imp_prod_loss( prod, loss, lsol, lrxt, lhet )
               do m = 1,clscnt4
                  forcing(m) = solution(m)*dti - (iter_invariant(m) + prod(m) - loss(m))
               end do
!-----------------------------------------------------------------------      
!        ... solve for the mixing ratio at t(n+1)
!-----------------------------------------------------------------------      
               call imp_lu_slv( sys_jac, forcing )
               do m = 1,clscnt4
                  solution(m) = solution(m) + forcing(m)
               end do
!-----------------------------------------------------------------------      
!        ... convergence measures
!-----------------------------------------------------------------------      
               if( nr_iter > 1 ) then
                  do k = 1,clscnt4
                     m = implicit%permute(k)
                     if( abs(solution(m)) > 1.e-40_dp ) then
                        max_delta(k) = abs( forcing(m)/solution(m) )
                     else
                        max_delta(k) = 0._dp
                     end if
                  end do
               end if

!-----------------------------------------------------------------------      
!   ... limit iterate
!-----------------------------------------------------------------------      
               where( solution(:) < 0._dp )
                  solution(:) = 0._dp
               endwhere
!-----------------------------------------------------------------------      
!    ... transfer latest solution back to work or base array
!-----------------------------------------------------------------------      
               do k = 1,clscnt4
                  j = implicit%clsmap(k)
                  m = implicit%permute(k)
                  lsol(j) = solution(m)
               end do

!-----------------------------------------------------------------------      
!    	... check for convergence
!-----------------------------------------------------------------------      
               if( nr_iter > 1 ) then
                  do k = 1,clscnt4
                     m = implicit%permute(k)
                     frc_mask = abs( forcing(m) ) > small
                     if( frc_mask ) then
                        converged(k) =  abs(forcing(m)) <= epsilon(k)*abs(solution(m))
                     else
                        converged(k) =  .true.
                     end if
                  end do
                  convergence = all( converged(:) )
                  if( convergence ) then
                     exit iter_loop
                  end if
               end if
            end do iter_loop

!-----------------------------------------------------------------------      
!    ... check for newton-raphson convergence
!-----------------------------------------------------------------------      
            if( .not. convergence ) then
!+++sschr
! The following lines show how to debug imp_sol
!             IF (indx == indx_debug) THEN
!               write(message_text,'(a,3i4,e11.3,i4)') 'imp_sol: iteration step failed to converge; p_pe, i, lev, dt, cut_cnt = ', &
!                         p_pe, i, lev, dt, cut_cnt   ! #debug# info
!               CALL message('moz_imp_sol', message_text, level=em_warn)
!             END IF
!             ddf1(i,lev,lat)=0.0_dp
!             IF (indx == indx_debug) THEN
!               write(message_text,'(a6,a5,a16)') " ######",'m','reaction_rates'
!               CALL message('moz_imp_sol', message_text, level=em_warn)
!               DO m = 1, rxntot
!                 write(message_text,'(a6,i5,e11.3)') " ######",m,reaction_rates(indx,m)
!                 CALL message('moz_imp_sol', message_text, level=em_warn)
!               END DO
!               write(message_text,'(a6,1x,'a4',1x,a8,1x,a5,1x,a1,1x,7a11)') " ######",'p_pe','tracnam', &
!                                         'm','c','max_delta', 'lsol', 'base_sol', 'prod', &
!                                         'loss', 'dt', 'int_done'
!               CALL message('moz_imp_sol', message_text, level=em_warn)
!               DO m = 1,clscnt4
!                 j = implicit%permute(m)
!                 write(message_text,'(a6,1x,i4,1x,a8,1x,i5,1x,L1,1x,7e11.3)') " ######",p_pe,tracnam(implicit%clsmap(m)), &
!                         m,converged(m),max_delta(m), lsol(m), base_sol(indx, m), prod(j), &
!                         loss(j), dt, interval_done
!                 CALL message('moz_imp_sol', message_text, level=em_warn)
!               END DO
!             ELSE
!               DO m = 1,clscnt4
!                 IF (.not. converged(m)) THEN
!                   j = implicit%permute(m)
!                   write(message_text,'(a6,1x,i4,1x,a8,1x,i5,1x,L1,1x,7e11.3,i5)') " ######",p_pe,tracnam(implicit%clsmap(m)), &
!                           m,converged(m),max_delta(m), lsol(m), base_sol(indx, m), prod(j), &
!                           loss(j), dt, interval_done, indx
!                   CALL message('moz_imp_sol', message_text, level=em_warn)
!                   END IF
!               END DO
!             END IF
!---sschr

!-----------------------------------------------------------------------      
!   ... non-convergence;
!       - reduce time step by factor of two and try again. 
!       - warning message if time step was reduced cut_limit times
!-----------------------------------------------------------------------      
               stp_con_cnt = 0
               if( cut_cnt < cut_limit ) then
                  cut_cnt = cut_cnt + 1
                  dt = .5_dp * dt
                  cycle      ! time step loop
               else
                  WRITE(message_text,'(a,4i5,e21.13)') ' imp_sol: failed to converge @ (lon,lat,lev,indx,dt) = ',i,lat,lev,indx,dt
                  CALL message('moz_imp_sol', message_text, level=em_warn, all_print=.true.)
                  write(message_text,'(a6,1x,a4,1x,a8,1x,a5,1x,a1,1x,7a11,1x,a5)') " ######",'p_pe','tracnam', &
                                            'm','c','max_delta', 'lsol', 'base_sol', 'prod', &
                                            'loss', 'dt', 'int_done','indx'
                  CALL message('moz_imp_sol', message_text, level=em_warn, all_print=.true.)
                  do m = 1,clscnt4
                    if (.not. converged(m)) THEN
                      j = implicit%permute(m)
                      write(message_text,'(a6,1x,i4,1x,a8,1x,i5,1x,L1,1x,7e11.3,i5)') " ######",p_pe,tracnam(implicit%clsmap(m)), &
                              m,converged(m),max_delta(m), lsol(m), base_sol(indx, m), prod(j), &
                              loss(j), dt, interval_done, indx
                      CALL message('moz_imp_sol', message_text, level=em_warn, all_print=.true.)
                      ENDIF
                  end do
               end if ! cut_cnt < cut_limit
            else
!             ddf1(i,lev,lat)=1.0_dp
            end if ! .not. convergence

!-----------------------------------------------------------------------      
!   ... check for interval done
!-----------------------------------------------------------------------      
            interval_done = interval_done + dt
!           ddf2(i,lev,lat)=dt
!           ddf3(i,lev,lat)=nr_iter
            if( abs( delt - interval_done ) <= .0001_dp ) then
               exit time_step_loop

!-----------------------------------------------------------------------      
!   ... transfer latest solution back to base array
!-----------------------------------------------------------------------      
            else
               if( convergence ) then
                  stp_con_cnt = stp_con_cnt + 1
               end if
               do m = 1,pcnstm1
                  base_sol(indx,m) = lsol(m)
               end do
               if( stp_con_cnt >= 2 ) then
                  dt = 2._dp*dt
                  stp_con_cnt = 0
               end if
               dt = min( dt,delt-interval_done )
!              if( pdiags%imp_slv ) then
!                 write(*,'('' imp_sol: new time step '',1p,e21.13)') dt
!              end if
            end if
!+++sschr
!           indx = (lev - 1)*plonl + i
!           if ((i == 1) .and. (lat == 3) .and. (lev == 9) .and. (indx == 65)) then
!             do m = 1,pcnstm1
!               write(*,'(1x,a8,1x,1pe10.3,1pe10.3)') tracnam(implicit%clsmap(m)), &
!                     max_delta(m), base_sol(indx, m) 
!             end do
!           end if
!---sschr

         end do time_step_loop

!-----------------------------------------------------------------------      
!    ... transfer latest solution back to base array
!-----------------------------------------------------------------------      
         do k = 1,clscnt4
            j = implicit%clsmap(k)
            m = implicit%permute(k)
            base_sol(indx,j) = solution(m)
         end do

!++mgs: removed history file stuff in scalar version
         end do lon_tile_loop
      end do level_loop


!++mgs: removed history file stuff

      end subroutine imp_sol

end module mo_moz_imp_sol

