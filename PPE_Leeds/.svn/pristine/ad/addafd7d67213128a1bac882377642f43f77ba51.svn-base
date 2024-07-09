!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_streams

! Description:
!
! Initialisation and deinitialisation routine for the ECHAM memory buffer
!
! Method:
!
! Initialise:
! Call 'new_stream' for each memory buffer to create it.
! Afterwards call the specific construct routines to fill each buffer
!
! Authors: A. Rhodin, MPI, May 2001
!          L. Kornblueh, MPI, April 2003, added port test
! 
  USE mo_memory_base,        ONLY: new_stream, print_memory_use, print_sinfo, &
                                   set_stream, delete_streams,                &
                                   ostreams, nstreams
  USE mo_exception,          ONLY: message
  USE mo_sub_nml,            ONLY: set_stream_element_nml, set_stream_nml
  USE mo_memory_sp,          ONLY: sp  ,construct_sp  ,destruct_sp
  USE mo_memory_ls,          ONLY: ls  ,construct_ls  ,destruct_ls
  USE mo_memory_gl,          ONLY: gl  ,construct_gl  ,destruct_gl
  USE mo_memory_f,           ONLY: f   ,construct_f   ,destruct_f
  USE mo_memory_g1a,         ONLY: g1a ,construct_g1a ,destruct_g1a
  USE mo_memory_g1b,         ONLY: g1b ,construct_g1b ,destruct_g1b
  USE mo_memory_g2a,         ONLY: g2a ,construct_g2a ,destruct_g2a
  USE mo_memory_g2b,         ONLY: g2b ,construct_g2b ,destruct_g2b
  USE mo_memory_g3a,         ONLY:      construct_g3a
  USE mo_memory_g3b,         ONLY: g3b ,construct_g3b ,destruct_g3b
  USE mo_buffer_fft,         ONLY:      construct_fft ,destruct_fft
  USE mo_memory_cfdiag,      ONLY: locfdiag,construct_cfdiag,destruct_cfdiag
  USE mo_tracdef,            ONLY: ntrac
  USE mo_control,            ONLY: ngl, nhgl, nlev, nlon, nmp1, nnp1, nsp,   &
                                   ldebugmem, ltdiag, lcolumn
  USE mo_decomposition,      ONLY: lc => local_decomposition
  USE mo_mpi,                ONLY: p_parallel_io
  USE mo_filename,           ONLY: ftyp => out_filetype
  USE mo_port_test,          ONLY: init_port_test
  USE mo_surface_memory,     ONLY: surf, construct_surface, destruct_surface
  USE mo_submodel_interface, ONLY: init_subm_memory
  USE mo_cosp_simulator,     ONLY: locosp, construct_stream_cosp
  USE mo_cosp_offline,       ONLY: locospoffl, construct_stream_cospoffl, &
                                   destruct_stream_cospoffl
  USE mo_diag_tendency_new,  ONLY: init_tdiag
  USE mo_mvstream,           ONLY: init_mvstream, mvstream_update_cache
  USE mo_iocolumn,           ONLY: scm, construct_scm, destruct_scm
!++mgs
  USE mo_submodel,           ONLY: lburden
  USE mo_tracer_processes,   ONLY: xt_burden_init_mem 
!--mgs
  USE mo_radiation_forcing,  ONLY: construct_forcing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_memory
  PUBLIC :: free_memory
!==============================================================================
CONTAINS
!==============================================================================
  SUBROUTINE init_memory
    INTEGER :: i
    !--------------------------------
    ! declare standard output streams
    !--------------------------------
    CALL new_stream (sp   ,'sp'  ,ftyp,lrerun=.FALSE.,              post_suf='_echam')
    CALL new_stream (ls   ,'ls'  ,ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='_echam')
    CALL new_stream (f    ,'f'   ,ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='_echam')
    CALL new_stream (gl   ,'gl'  ,ftyp,lrerun=.FALSE.,              post_suf='_echam')
    CALL new_stream (g1a  ,'g1a' ,ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='_echam')
    CALL new_stream (g1b  ,'g1b' ,ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='_echam')
    CALL new_stream (g2a  ,'g2a' ,ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='_echam')
    CALL new_stream (g2b  ,'g2b' ,ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='_echam')
    CALL new_stream (g3b  ,'g3b' ,ftyp,lrerun=.FALSE.,              post_suf='_echam')
    CALL new_stream (surf ,'surf',ftyp,lrerun=.FALSE.,              post_suf='_surf')
    CALL new_stream (scm  ,'scm' ,ftyp,lrerun=.FALSE.,              post_suf='_echam')

    !-------------------------------------------------------------
    ! set restart flags and restart file suffixes (new convention)
    !-------------------------------------------------------------
    CALL set_stream (f   ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 31
    CALL set_stream (g1a ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 35
    CALL set_stream (g2a ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 36
    CALL set_stream (g3b ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 37
    CALL set_stream (gl  ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 32
    CALL set_stream (surf,lrerun=.TRUE.  ,rest_suf= '_surf' ) !
    !----------------------------------
    ! declare fields within each buffer
    !----------------------------------
    CALL construct_sp  (nlev,       lc%nsnm0, lc%snsp, &
                        nlev,       nnp1,     nsp)
    CALL construct_ls  (lc%nllev, lc%nllevp1, lc%nlnm0, lc%lnsp, &
                        nlev,                 nnp1,     nsp)
    CALL construct_f   (lc%nllev, lc%nllevp1, lc%nlm,  lc%nlat/2, &  
                        nlev,                 nmp1,    nhgl)
    CALL construct_g1b (lc%nproma,   lc%nlev, ntrac, lc%ngpblks, &
                        nlon,       nlev,    ntrac, ngl)
    CALL construct_g2a (lc%nproma,   lc%nlev, lc%ngpblks, &
                        nlon,       nlev,    ngl)
    CALL construct_g2b (lc%nproma,   lc%nlev, lc%ngpblks, &
                        nlon,       nlev,    ngl)
    CALL construct_g3b

    CALL construct_g3a

    CALL construct_gl  (lc%nproma,   lc%nlev, ntrac, lc%ngpblks, &
                        nlon,       nlev,    ntrac, ngl)
    CALL construct_g1a (lc%nproma,   lc%nlev, ntrac, lc%ngpblks, &
                        nlon,       nlev,    ntrac, ngl)
    CALL construct_fft (lc)
 
    CALL construct_surface

    IF ( locfdiag ) THEN
     CALL construct_cfdiag
    END IF

    IF (lcolumn) CALL construct_scm

    CALL construct_forcing

!++mgs
    IF (lburden) CALL xt_burden_init_mem
!--mgs

!!mgs&jsr!!
!!    IF (lanysubmodel) CALL init_subm_memory
    CALL init_subm_memory

  ! Initialize tendency diagnostic
    IF (ltdiag) CALL init_tdiag
    
    IF (locosp) CALL construct_stream_cosp

    IF (locospoffl) CALL construct_stream_cospoffl
!!mgs&jsr!!
    CALL init_mvstream

    CALL init_port_test

    CALL set_stream_element_nml
    CALL set_stream_nml

    CALL mvstream_update_cache

    IF (ldebugmem) THEN
      CALL message('',' ')
      CALL message('',' Global memory buffers:')
      DO i = 1, nstreams
        CALL print_sinfo(ostreams(i))
      ENDDO
    END IF

    IF (p_parallel_io) THEN
      CALL message('',' ')
      CALL message('',' Global memory buffers:')
      DO i = 1, nstreams
        CALL print_memory_use(ostreams(i))
      ENDDO
    END IF

  END SUBROUTINE init_memory
!------------------------------------------------------------------------------
  !
  ! deallocate all memory allocated in the ECHAM memory buffer
  !
  SUBROUTINE free_memory

  USE mo_submodel_interface, ONLY: free_subm_memory
!!mgs&jsr!!
  USE mo_submodel,           ONLY: lanysubmodel

    CALL destruct_sp
    CALL destruct_ls
    CALL destruct_f
    CALL destruct_gl
    CALL destruct_g1a
    CALL destruct_g1b
    CALL destruct_g2a
    CALL destruct_g2b
    CALL destruct_g3b
    CALL destruct_fft
    CALL destruct_surface
    IF(lcolumn) CALL destruct_scm
    IF ( locfdiag ) THEN
     CALL destruct_cfdiag
    END IF
    IF ( locospoffl ) THEN
     CALL  destruct_stream_cospoffl
    END IF

!!mgs&jsr!!
    IF (lanysubmodel) CALL free_subm_memory

    CALL delete_streams

  END SUBROUTINE free_memory
!------------------------------------------------------------------------------
END MODULE mo_memory_streams
