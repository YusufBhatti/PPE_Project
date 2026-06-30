!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_timer

  USE mo_real_timer, ONLY: new_timer, timer_report, &
                           timer_start, timer_stop, &
                           timer_reset_all, timer_reset_top

  IMPLICIT NONE

  PRIVATE

  ! Time integration timer

  INTEGER, PUBLIC :: timer_total
  INTEGER, PUBLIC :: timer_loop
  INTEGER, PUBLIC :: timer_single_timestep

  ! Transport timer

  INTEGER, PUBLIC :: timer_slt
  INTEGER, PUBLIC :: timer_tpcore

  ! I/O timer

  INTEGER, PUBLIC :: timer_output
  INTEGER, PUBLIC :: timer_restart

  ! Non gpc calls

  INTEGER, PUBLIC :: timer_dyn
  INTEGER, PUBLIC :: timer_ewd
  INTEGER, PUBLIC :: timer_si1
  INTEGER, PUBLIC :: timer_si2
  INTEGER, PUBLIC :: timer_sym1
  INTEGER, PUBLIC :: timer_sym2
  INTEGER, PUBLIC :: timer_ffti
  INTEGER, PUBLIC :: timer_fftd
  INTEGER, PUBLIC :: timer_lti
  INTEGER, PUBLIC :: timer_ltd
  INTEGER, PUBLIC :: timer_tf1
  INTEGER, PUBLIC :: timer_tf2

  ! physics timer

  INTEGER, PUBLIC :: timer_radiation
  INTEGER, PUBLIC :: timer_rrtm
  INTEGER, PUBLIC :: timer_sw
  INTEGER, PUBLIC :: timer_cloud
  INTEGER, PUBLIC :: timer_cover
  INTEGER, PUBLIC :: timer_vdiff
  INTEGER, PUBLIC :: timer_radheat
  INTEGER, PUBLIC :: timer_gwdrag
  INTEGER, PUBLIC :: timer_cucall
  INTEGER, PUBLIC :: timer_bclistread

  ! submodel timers (general level)

  INTEGER, PUBLIC :: timer_physc_subm_1
  INTEGER, PUBLIC :: timer_physc_subm_2
  INTEGER, PUBLIC :: timer_physc_subm_3
  INTEGER, PUBLIC :: timer_physc_subm_4
  INTEGER, PUBLIC :: timer_cloud_subm
  INTEGER, PUBLIC :: timer_vdiff_subm
  INTEGER, PUBLIC :: timer_radiation_subm_1
  INTEGER, PUBLIC :: timer_radiation_subm_2
  INTEGER, PUBLIC :: timer_cuflx_subm
  INTEGER, PUBLIC :: timer_scan1_subm

  ! transpose timer

  INTEGER, PUBLIC :: timer_s2l
  INTEGER, PUBLIC :: timer_l2s
  INTEGER, PUBLIC :: timer_f2l
  INTEGER, PUBLIC :: timer_l2f
  INTEGER, PUBLIC :: timer_g2f
  INTEGER, PUBLIC :: timer_f2g

  INTEGER, PUBLIC :: timer_g2f_pack_gp_buf  
  INTEGER, PUBLIC :: timer_g2f_unpack_buf_fs
  
  INTEGER, PUBLIC :: timer_f2l_pack_fs_buf  
  INTEGER, PUBLIC :: timer_f2l_unpack_buf_ls
  
  INTEGER, PUBLIC :: timer_l2s_pack_ls_buf  
  INTEGER, PUBLIC :: timer_l2s_unpack_buf_sp
  
  INTEGER, PUBLIC :: timer_s2l_pack_sp_buf  
  INTEGER, PUBLIC :: timer_s2l_unpack_buf_ls
  
  INTEGER, PUBLIC :: timer_l2f_pack_ls_buf  
  INTEGER, PUBLIC :: timer_l2f_unpack_buf_fs
  
  INTEGER, PUBLIC :: timer_f2g_pack_fs_buf  
  INTEGER, PUBLIC :: timer_f2g_unpack_buf_gp

  INTEGER, PUBLIC :: timer_gather2
  INTEGER, PUBLIC :: timer_gather3

  ! stepon timer

  INTEGER, PUBLIC :: timer_time_set
  INTEGER, PUBLIC :: timer_time_reset
  INTEGER, PUBLIC :: timer_trigfiles
  INTEGER, PUBLIC :: timer_bcond
  INTEGER, PUBLIC :: timer_subm
  INTEGER, PUBLIC :: timer_jsbach
  INTEGER, PUBLIC :: timer_scan1
  INTEGER, PUBLIC :: timer_sccd
  INTEGER, PUBLIC :: timer_scctp
  INTEGER, PUBLIC :: timer_uspnge
  INTEGER, PUBLIC :: timer_hdiff
  INTEGER, PUBLIC :: timer_inv_legendre
  INTEGER, PUBLIC :: timer_prerad
  
  INTEGER, PUBLIC :: timer_couple_get
  INTEGER, PUBLIC :: timer_couple_put
  INTEGER, PUBLIC :: timer_couple_end
  INTEGER, PUBLIC :: timer_hd          ! hd model

  PUBLIC :: init_timer, print_timer, timer_start, timer_stop, cleanup_timer

CONTAINS

  SUBROUTINE init_timer

    timer_total             = new_timer('total')
    timer_loop              = new_timer('stepon: loop')
    timer_single_timestep   = new_timer('single_timestep')
 
    timer_radiation         = new_timer('physc:radiation')
    timer_rrtm              = new_timer('radiation:rrtm')
    timer_sw                = new_timer('radiation:sw')
    timer_cloud             = new_timer('physc:cloud')
    timer_cover             = new_timer('physc:cover')
    timer_vdiff             = new_timer('physc:vdiff')
    timer_radheat           = new_timer('physc:radheat')
    timer_gwdrag            = new_timer('physc:gwdrag')
    timer_cucall            = new_timer('physc:cucall')

    timer_bclistread        = new_timer('bc_list_read')

  ! submodel timers (general level)
    timer_physc_subm_1      = new_timer('physc_subm_1')
    timer_physc_subm_2      = new_timer('physc_subm_2')
    timer_physc_subm_3      = new_timer('physc_subm_3')
    timer_physc_subm_4      = new_timer('physc_subm_4')
    timer_cloud_subm        = new_timer('cloud_subm')
    timer_vdiff_subm        = new_timer('vdiff_subm')
    timer_radiation_subm_1  = new_timer('radiation_subm_1')
    timer_radiation_subm_2  = new_timer('radiation_subm_2')
    timer_cuflx_subm        = new_timer('cuflx_subm')
    timer_scan1_subm        = new_timer('scan1_subm')

  ! Transpose timer

    timer_slt               = new_timer('slt')
    timer_tpcore            = new_timer('tpcore')

    timer_dyn               = new_timer('scan1:dyn')
    timer_ewd               = new_timer('scan1:ewd')
    timer_si1               = new_timer('scan1:si1')
    timer_si2               = new_timer('scan1:si2')
    timer_sym1              = new_timer('scan1:sym1')
    timer_sym2              = new_timer('scan1:sym2')
    timer_ffti              = new_timer('scan1:ffti')
    timer_fftd              = new_timer('scan1:fftd')
    timer_lti               = new_timer('scan1:lti')
    timer_ltd               = new_timer('scan1:ltd')
    timer_tf1               = new_timer('scan1:tf1')
    timer_tf2               = new_timer('scan1:tf2')

    timer_s2l               = new_timer('spectral to legendre')    
    timer_l2s               = new_timer('legendre to spectral')
    timer_f2l               = new_timer('fourier to legendre')
    timer_l2f               = new_timer('legendre to fourier')
    timer_g2f               = new_timer('gridpoint to fourier')
    timer_f2g               = new_timer('fourier to gridpoint')

    timer_g2f_pack_gp_buf   = new_timer('g2f: pack gp send buffer')
    timer_g2f_unpack_buf_fs = new_timer('g2f: unpack fs recv buffer')
    timer_f2l_pack_fs_buf   = new_timer('f2l: pack fs send buffer')
    timer_f2l_unpack_buf_ls = new_timer('f2l: unpack ls recv buffer')
    timer_l2s_pack_ls_buf   = new_timer('l2s: pack ls send buffer')
    timer_l2s_unpack_buf_sp = new_timer('l2s: unpack sp recv buffer')
    timer_s2l_pack_sp_buf   = new_timer('s2l: pack sp send buffer')
    timer_s2l_unpack_buf_ls = new_timer('s2l: unpack ls recv buffer')
    timer_l2f_pack_ls_buf   = new_timer('l2f: pack ls send buffer') 
    timer_l2f_unpack_buf_fs = new_timer('l2f: unpack fs recv buffer')
    timer_f2g_pack_fs_buf   = new_timer('f2g: pack fs send buffer')
    timer_f2g_unpack_buf_gp = new_timer('f2g: unpack gp recv buffer')

    timer_gather2           = new_timer('gather 2d')
    timer_gather3           = new_timer('gather 3d')    

    timer_time_set          = new_timer('stepon: time manager set')    
    timer_time_reset        = new_timer('stepon: time manager reset')    
    timer_trigfiles         = new_timer('stepon: trigger files')    
    timer_bcond             = new_timer('stepon: boundary conditions')    
    timer_subm              = new_timer('stepon: sub models')    
    timer_jsbach            = new_timer('stepon: jsbach')    
    timer_scan1             = new_timer('stepon: scan1')    
    timer_sccd              = new_timer('stepon: sccd')    
    timer_scctp             = new_timer('stepon: scctp')    
    timer_uspnge            = new_timer('stepon: uspgne')    
    timer_hdiff             = new_timer('stepon: hdiff')    
    timer_inv_legendre      = new_timer('stepon: inverse legendre')    
    timer_prerad            = new_timer('stepon: prerad')    
    timer_output            = new_timer('stepon: output IO')
    timer_restart           = new_timer('stepon: restart IO')

    timer_couple_get        = new_timer('stepon: couple get') 
    timer_couple_put        = new_timer('stepon: couple put') 
    timer_couple_end        = new_timer('stepon: couple end') 
    timer_hd                = new_timer('stepon: hydrology model') 

  END SUBROUTINE init_timer

  SUBROUTINE print_timer(short)
    LOGICAL, OPTIONAL, INTENT(in) :: short

    CALL timer_report(short=short)

  END SUBROUTINE print_timer

  SUBROUTINE cleanup_timer
    
    CALL timer_reset_all
    CALL timer_reset_top

  END SUBROUTINE cleanup_timer

END MODULE mo_timer



