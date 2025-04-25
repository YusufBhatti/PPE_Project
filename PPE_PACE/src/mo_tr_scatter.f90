!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_tr_scatter

  ! Luis Kornblueh, MPIM, March 2010, initial version
  !
  !
  ! Experimental code for scattering data this code is not 
  ! intended for production purposes but to explore the scaling
  ! of some aspects of the MPI communication and the combination 
  ! with OpenMP

  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, message_text, finish
  USE mo_decomposition, ONLY: ldc => local_decomposition,  &
                              gdc => global_decomposition
  USE mo_mpi,           ONLY: p_nprocs, p_pe, p_io,        &
                              p_recv, p_send
  USE mo_transpose,     ONLY: reorder

  IMPLICIT NONE

  PRIVATE

  INTERFACE scatter_field
    MODULE PROCEDURE scatter_gp432
    MODULE PROCEDURE scatter_gp32
    MODULE PROCEDURE scatter_gp2
  END INTERFACE scatter_field

  INTERFACE scatter_spectral
    MODULE PROCEDURE scatter_sp4
    MODULE PROCEDURE scatter_sp3
    MODULE PROCEDURE scatter_sp0
  END INTERFACE scatter_spectral

  INTERFACE scatter_legendre
    MODULE PROCEDURE scatter_ls3
    MODULE PROCEDURE scatter_ls0
  END INTERFACE scatter_legendre

  INTERFACE scatter_fourier
    MODULE PROCEDURE scatter_sa42
  END INTERFACE scatter_fourier

  PUBLIC :: scatter_field
  PUBLIC :: scatter_spectral
  PUBLIC :: scatter_fourier

  INTEGER, PARAMETER :: tag_scatter_gp = 302 ! gridpoint
  INTEGER, PARAMETER :: tag_scatter_sp = 402 ! spectral
  INTEGER, PARAMETER :: tag_scatter_sa = 502 ! symmetric/anti-symmetric Fourier
  INTEGER, PARAMETER :: tag_scatter_ls = 602 ! spectral

CONTAINS

  SUBROUTINE scatter_gp432 (global, global_dims, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:,:,:)
    INTEGER,           INTENT(in)  :: global_dims(:)
    REAL(dp), TARGET,  INTENT(out) :: local(:,:,:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    REAL(dp), POINTER :: global3d(:,:,:)

    INTEGER :: i, nt
    INTEGER :: dimsize4

    LOGICAL :: lselect

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    dimsize4 = global_dims(4)

    NULLIFY(global3d)
    IF (dimsize4 == 1) THEN
      IF (p_pe == p_io .OR. lselect) global3d => global(:,:,:,1)
      CALL scatter_gp32 (global3d, global_dims, local(:,:,:,1), lselect)
    ELSE
      nt = global_dims(3)
      DO i = 1, nt
        IF (p_pe == p_io .OR. lselect) global3d => global(:,:,i,:)
        CALL scatter_gp32 (global3d, global_dims, local(:,:,i,:), lselect)
      ENDDO
    ENDIF
  END SUBROUTINE scatter_gp432
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_gp32 (global, global_dims, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:,:)
    INTEGER,           INTENT(in)  :: global_dims(:)
    REAL(dp), TARGET,  INTENT(out) :: local(:,:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    REAL(dp), POINTER :: global2d(:,:)

    INTEGER :: dimsize3

    LOGICAL :: lselect

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    dimsize3 = global_dims(3)

    NULLIFY(global2d)
    IF (dimsize3 == 1) THEN
      IF (p_pe == p_io .OR. lselect) global2d => global(:,:,1)
      CALL scatter_gp2 (global2d, local(:,:,1), lselect)
    ELSE
      CALL scatter_gp3 (global, local, lselect)
    ENDIF

  END SUBROUTINE scatter_gp32
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_gp2 (global, local, lall_read)
#ifdef STANDALONE
    USE mo_jsbach_transpose, ONLY: jsb_scatter_gp2
#endif
    REAL(dp), POINTER, INTENT(in)  :: global(:,:)
    REAL(dp), TARGET , INTENT(out) :: local(:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    REAL(dp), POINTER :: sendbuf(:,:)
    REAL(dp), POINTER :: recvbuf(:,:)

    INTEGER :: i, pe, nlon
    LOGICAL :: lreg, lselect

#ifdef STANDALONE
    CALL jsb_scatter_gp2(global, local(:,:))
    RETURN
#endif

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    nlon  = ldc%nlon
    lreg  = ldc%lreg

    IF (lreg) THEN
      recvbuf => local
    ELSE
      ALLOCATE (recvbuf(ldc%nglon,ldc%nglat))
    ENDIF

    IF (lselect) THEN
      recvbuf(:,:ldc%nglh(1)) = global(ldc%glons(1):ldc%glone(1), &
                                       ldc%glats(1):ldc%glate(1))

      IF (ldc%nglh(2) > 0) THEN
        IF (ldc%glone(2) > ldc%glons(2)) THEN
          recvbuf(:,ldc%nglat/2+1:) = global(ldc%glons(2):ldc%glone(2), &
                                             ldc%glats(2):ldc%glate(2))
        ELSE
          recvbuf(:nlon-ldc%glons(2)+1,ldc%nglat/2+1:) = &
               global(ldc%glons(2):nlon,ldc%glats(2):ldc%glate(2))
          recvbuf(ldc%nglon-ldc%glone(2)+1:,ldc%nglat/2+1:) = &
               global(1:ldc%glone(2),ldc%glats(2):ldc%glate(2))
        ENDIF
      ENDIF
    ELSE
      IF (p_pe == p_io) THEN
        DO i = 1, p_nprocs
          pe = gdc(i)%pe
          ALLOCATE (sendbuf(gdc(i)%nglon,gdc(i)%nglat))
          
          sendbuf(:,:gdc(i)%nglh(1)) = global(gdc(i)%glons(1):gdc(i)%glone(1), &
                                              gdc(i)%glats(1):gdc(i)%glate(1))
          
          IF (gdc(i)%nglh(2) > 0) THEN
            IF (gdc(i)%glone(2) > gdc(i)%glons(2)) THEN
              sendbuf(:,gdc(i)%nglat/2+1:) = global(gdc(i)%glons(2):gdc(i)%glone(2), &
                                                    gdc(i)%glats(2):gdc(i)%glate(2))
            ELSE
              sendbuf(:nlon-gdc(i)%glons(2)+1,gdc(i)%nglat/2+1:) = &
                   global(gdc(i)%glons(2):nlon,gdc(i)%glats(2):gdc(i)%glate(2))
              sendbuf(gdc(i)%nglon-gdc(i)%glone(2)+1:,gdc(i)%nglat/2+1:) = &
                   global(1:gdc(i)%glone(2),gdc(i)%glats(2):gdc(i)%glate(2))
            ENDIF
          ENDIF

          IF (pe /= p_pe) THEN
            CALL p_send( sendbuf, pe, tag_scatter_gp)
          ELSE
            recvbuf = sendbuf
          ENDIF
          DEALLOCATE (sendbuf)
        ENDDO
      ELSE
        CALL p_recv (recvbuf, p_io, tag_scatter_gp)
      ENDIF
    ENDIF

    IF (.NOT. lreg) THEN
      CALL reorder (local, recvbuf)
      DEALLOCATE (recvbuf)
    ENDIF

  END SUBROUTINE scatter_gp2
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_gp3 (global, local, lall_read)
#ifdef STANDALONE
    USE mo_jsbach_transpose, ONLY: jsb_scatter_gp3
#endif
    REAL(dp), POINTER, INTENT(in)  :: global(:,:,:)
    REAL(dp), TARGET , INTENT(out) :: local(:,:,:)
    LOGICAL,           INTENT(in)  :: lall_read

    REAL(dp), POINTER :: sendbuf(:,:,:)
    REAL(dp), POINTER :: recvbuf(:,:,:)

    INTEGER :: i, pe, nlon
    LOGICAL :: lreg

#ifdef STANDALONE
    CALL jsb_scatter_gp3(global, local(:,:,:))
    RETURN
#endif
    nlon = ldc%nlon
    lreg = ldc%lreg

    IF (lreg) THEN
      recvbuf => local
    ELSE
      ALLOCATE (recvbuf(ldc%nglon,SIZE(local,2),ldc%nglat))
    ENDIF

    IF (lall_read) THEN
      recvbuf(:,:,:ldc%nglh(1)) = global(ldc%glons(1):ldc%glone(1),:, &
                                         ldc%glats(1):ldc%glate(1))

      IF (ldc%nglh(2) > 0) THEN
        IF (ldc%glone(2) > ldc%glons(2)) THEN
          recvbuf(:,:,ldc%nglh(1)+1:) = global(ldc%glons(2):ldc%glone(2),:, &
                                               ldc%glats(2):ldc%glate(2))
        ELSE
          recvbuf(:nlon-ldc%glons(2)+1,:,ldc%nglh(1)+1:) = &
               global(ldc%glons(2):nlon,:,ldc%glats(2):ldc%glate(2))
          recvbuf(ldc%nglon-ldc%glone(2)+1:,:,ldc%nglh(1)+1:) = &
               global (1:ldc%glone(2),:,ldc%glats(2):ldc%glate(2))
        ENDIF
      ENDIF
    ELSE
      IF (p_pe == p_io) THEN
        DO i = 1, p_nprocs
          pe = gdc(i)%pe
          ALLOCATE (sendbuf (gdc(i)%nglon, SIZE(global,2), gdc(i)%nglat))
          
          sendbuf(:,:,:gdc(i)%nglh(1)) = global(gdc(i)%glons(1):gdc(i)%glone(1),:, &
                                                gdc(i)%glats(1):gdc(i)%glate(1))

          IF (gdc(i)%nglh(2) > 0) THEN
            IF (gdc(i)%glone(2) > gdc(i)%glons(2)) THEN
              sendbuf(:,:,gdc(i)%nglh(1)+1:) = global(gdc(i)%glons(2):gdc(i)%glone(2),:, &
                                                      gdc(i)%glats(2):gdc(i)%glate(2))
            ELSE
              sendbuf(:nlon-gdc(i)%glons(2)+1,:,gdc(i)%nglh(1)+1:) = &
                   global(gdc(i)%glons(2):nlon,:,gdc(i)%glats(2):gdc(i)%glate(2))
              sendbuf(gdc(i)%nglon-gdc(i)%glone(2)+1:,:,gdc(i)%nglh(1)+1:) = &
                   global (1:gdc(i)%glone(2),:,gdc(i)%glats(2):gdc(i)%glate(2))
            ENDIF
          ENDIF
          
          IF (pe /= p_pe) THEN
            CALL p_send( sendbuf, pe, tag_scatter_gp)
          ELSE
            recvbuf = sendbuf
          ENDIF
          DEALLOCATE (sendbuf)
        ENDDO
      ELSE
        CALL p_recv (recvbuf, p_io, tag_scatter_gp)
      ENDIF
    ENDIF

    IF (.NOT.lreg) THEN
      CALL reorder (local,recvbuf)
      DEALLOCATE (recvbuf)
    ENDIF

  END SUBROUTINE scatter_gp3
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_ls3 (global, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:,:)
    REAL(dp), TARGET , INTENT(out) :: local(:,:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    REAL(dp), POINTER :: sendbuf(:,:,:)

    INTEGER :: i, im, pe, mp1, llevs, lleve, lnsp, nlm
    INTEGER :: ke, nk, mp, np, mpgl

    LOGICAL :: lselect

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    IF (lselect) THEN
      IF (SIZE(local) > 0) THEN
        llevs = ldc%llevs
        lleve = ldc%lleve
        nlm   = ldc%nlm
        lnsp  = ldc%lnsp
        
        ke = MIN (lleve,SIZE(global,1))
        nk = ke-llevs+1
        IF (nk*lnsp < 1) RETURN

        mp = 0
        DO im = 1, nlm
          mp1  = ldc%lm(im)+1
          np   = ldc%nnp(mp1)
          mpgl = ldc%nmp(mp1)
          local(:,:,mp+1:mp+np) = global(llevs:ke,:,mpgl+1:mpgl+np)
          mp = mp+np
        ENDDO
      ENDIF
    ELSE
      IF (p_pe == p_io) THEN
        DO i = 1, p_nprocs
          pe    = gdc(i)%pe
          llevs = gdc(i)%llevs
          lleve = gdc(i)%lleve
          nlm   = gdc(i)%nlm
          lnsp  = gdc(i)%lnsp
          
          ke = MIN (lleve,SIZE(global,1))
          nk = ke-llevs+1
          IF (nk*lnsp < 1) CYCLE
          
          ALLOCATE (sendbuf (nk, 2, lnsp))
          
          mp = 0
          DO im = 1, nlm
            mp1  = gdc(i)%lm(im)+1
            np   = gdc(i)%nnp(mp1)
            mpgl = gdc(i)%nmp(mp1)
            sendbuf(:,:,mp+1:mp+np) = global(llevs:ke,:,mpgl+1:mpgl+np)
            mp = mp+np
          ENDDO
          
          IF (pe /= p_pe) THEN
            CALL p_send(sendbuf, pe, tag_scatter_ls)
          ELSE
            local = sendbuf
          ENDIF
          DEALLOCATE (sendbuf)
        ENDDO
      ELSE
        IF (SIZE(local) > 0) CALL p_recv(local, p_io, tag_scatter_ls)
      ENDIF
    ENDIF

  END SUBROUTINE scatter_ls3
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_ls0 (global, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:)
    REAL(dp), TARGET , INTENT(out) :: local(:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    REAL(dp), POINTER :: sendbuf(:,:)

    INTEGER :: i, pe, llevs, lleve, ke, nk, nlnm0

    LOGICAL :: lselect

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    IF (lselect) THEN
      IF (SIZE(local) > 0) THEN
        llevs  = ldc%llevs
        lleve  = ldc%lleve
        nlnm0  = ldc%nlnm0
        IF (nlnm0 == 0) RETURN
        ke = MIN(lleve,SIZE(global,1))
        nk = ke-llevs+1
        IF (nk < 1) RETURN
        local(:,:) = global(llevs:ke,:)
      ENDIF
    ELSE
      IF (p_pe == p_io) THEN
        DO i = 1, p_nprocs
          pe     = gdc(i)%pe
          llevs  = gdc(i)%llevs
          lleve  = gdc(i)%lleve
          nlnm0  = gdc(i)%nlnm0
          IF (nlnm0 == 0) CYCLE
          ke = MIN(lleve,SIZE(global,1))
          nk = ke-llevs+1
          IF (nk < 1) CYCLE
          
          ALLOCATE (sendbuf (nk, nlnm0))
          
          sendbuf(:,:) = global(llevs:ke,:)
          
          IF (pe /= p_pe) THEN
            CALL p_send(sendbuf, pe, tag_scatter_ls)
          ELSE
            local = sendbuf
          ENDIF
          DEALLOCATE (sendbuf)
        ENDDO
      ELSE
        IF (SIZE(local) > 0) CALL p_recv(local, p_io, tag_scatter_ls)
      ENDIF
    ENDIF

  END SUBROUTINE scatter_ls0
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_sp4 (global, global_dims, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:,:,:)
    INTEGER,           INTENT(in)  :: global_dims(:)
    REAL(dp),          INTENT(out) :: local(:,:,:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    INTEGER :: dimsize3, dimsize4

    REAL(dp), POINTER :: global3d(:,:,:)
    REAL(dp), POINTER :: global2d(:,:)
    
    LOGICAL :: lselect

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    dimsize3 = global_dims(3)
    dimsize4 = global_dims(4)

    IF (dimsize4 /= 1) &
         CALL finish('scatter_sp4','dimsize4 /= 1')

    NULLIFY (global2d)
    NULLIFY (global3d) 
    IF (dimsize3 == 1) THEN
      IF (p_pe == p_io .OR. lselect) global2d => global(:,:,1,1)
      CALL scatter_sp0 (global2d, local(:,:,1,1), lselect)
    ELSE
      IF (p_pe == p_io .OR. lselect) global3d => global(:,:,:,1)
      CALL scatter_sp3 (global3d, local(:,:,:,1), lselect)
    ENDIF
  END SUBROUTINE scatter_sp4
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_sp3 (global, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:,:)
    REAL(dp),          INTENT(out) :: local(:,:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    REAL(dp), POINTER :: sendbuf(:,:,:)

    INTEGER :: i, im, pe, nlev, mp1, snsp, nsm, mp, np, mpgl

    LOGICAL :: lselect

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    IF (lselect) THEN
      IF (SIZE(local) > 0) THEN
        nsm  = ldc%nsm
        snsp = ldc%snsp
        nlev = SIZE(global,1)
        IF (snsp < 1) RETURN
        mp = 0
        DO im = 1, nsm
          mp1  = ldc%sm(im)+1
          np   = ldc%snnp(im)
          mpgl = ldc%nmp(mp1)+ldc%snn0(im)
          local(:,:,mp+1:mp+np) = global(:,:,mpgl+1:mpgl+np)
          mp = mp+np
        ENDDO
        IF (mp /= snsp) THEN
          WRITE(message_text,*) 'scatter_sp: PE',pe,',mp/=snsp:',mp,snsp
          CALL message('', TRIM(message_text))
          CALL finish('scatter_sp','mp/=snsp')
        ENDIF
      ENDIF
    ELSE
      IF (p_pe == p_io) THEN
        DO i = 1, p_nprocs
          pe   = gdc(i)%pe
          nsm  = gdc(i)%nsm
          snsp = gdc(i)%snsp
          nlev = SIZE(global,1)
          IF (snsp < 1) CYCLE
          
          ALLOCATE (sendbuf (nlev, 2, snsp))
          
          mp = 0
          DO im = 1, nsm
            mp1  = gdc(i)%sm(im)+1
            np   = gdc(i)%snnp(im)
            mpgl = gdc(i)%nmp(mp1)+gdc(i)%snn0(im)
            sendbuf(:,:,mp+1:mp+np) = global(:,:,mpgl+1:mpgl+np)
            mp = mp+np
          ENDDO
          IF (mp /= snsp) THEN
            WRITE(message_text,*) 'scatter_sp: PE',pe,',mp/=snsp:',mp,snsp
            CALL message('', TRIM(message_text))
            CALL finish('scatter_sp','mp/=snsp')
          ENDIF
          
          IF (pe /= p_pe) THEN
            CALL p_send(sendbuf, pe, tag_scatter_sp)
          ELSE
            local = sendbuf
          ENDIF
          DEALLOCATE (sendbuf)
        ENDDO
      ELSE
        IF (SIZE(local) > 0) CALL p_recv (local(:,:,:), p_io, tag_scatter_sp)
      ENDIF
    ENDIF

  END SUBROUTINE scatter_sp3
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_sp0 (global, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:)
    REAL(dp),          INTENT(out) :: local(:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    REAL(dp), POINTER :: sendbuf (:,:)

    INTEGER :: i, pe, nsnm0, snn0, nlev

    LOGICAL :: lselect

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    IF (lselect) THEN
      IF (SIZE(local) > 0) THEN
        nlev   = SIZE(global,1)
        nsnm0  = ldc%nsnm0
        IF (nsnm0 == 0) RETURN
        snn0   = ldc%snn0(1)
        local(:,:) = global(:,1+snn0:nsnm0+snn0)
      ENDIF
    ELSE
      IF (p_pe == p_io) THEN
        DO i = 1, p_nprocs
          pe     = gdc(i)%pe
          nlev   = SIZE(global,1)
          nsnm0  = gdc(i)%nsnm0
          IF (nsnm0 == 0) CYCLE
          snn0   = gdc(i)%snn0(1)
          
          ALLOCATE (sendbuf (nlev, nsnm0))
          
          sendbuf(:,:) = global(:,1+snn0:nsnm0+snn0)
          
          IF (pe /= p_pe) THEN
            CALL p_send(sendbuf, pe, tag_scatter_sp)
          ELSE
            local = sendbuf
          ENDIF
          DEALLOCATE (sendbuf)
        ENDDO
      ELSE
        IF (SIZE(local) > 0) CALL p_recv (local, p_io, tag_scatter_sp)
      ENDIF
    ENDIF

  END SUBROUTINE scatter_sp0
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_sa42 (global, global_dims, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:,:,:)
    INTEGER,           INTENT(in)  :: global_dims(:)
    REAL(dp),          INTENT(out) :: local(:,:,:,:)
    LOGICAL, OPTIONAL, INTENT(in)  :: lall_read

    INTEGER :: dimsize, dimsize3, dimsize4

    REAL(dp), POINTER :: global2d(:,:)

    LOGICAL :: lselect

    IF (PRESENT(lall_read)) THEN
      lselect = lall_read
    ELSE
      lselect = .FALSE.
    ENDIF

    dimsize3 = global_dims(3)
    dimsize4 = global_dims(4)
    dimsize = dimsize3*dimsize4

    NULLIFY(global2d)
    IF (dimsize == 1) THEN
      IF (p_pe == p_io .OR. lselect) global2d => global(:,:,1,1)
      CALL scatter_sa2 (global2d, local(:,:,1,1), lselect)
    ELSE
      CALL scatter_sa4 (global, local, lselect)
    ENDIF

  END SUBROUTINE scatter_sa42
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_sa4 (global, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:,:,:)
    REAL(dp),          INTENT(out) :: local(:,:,:,:)
    LOGICAL,           INTENT(in)  :: lall_read

    REAL(dp), ALLOCATABLE :: sendbuf(:,:,:,:)

    INTEGER :: i, im, pe, mp1, llevs, lleve, nlm, nhgl, ke, nk

    IF (lall_read) THEN
      IF (SIZE(local,1) > 0) THEN
        llevs  = ldc%llevs
        lleve  = ldc%lleve
        nlm    = ldc%nlm
        nhgl   = ldc%nlat/2
        ke = MIN (lleve,SIZE(global,1))
        nk = ke-llevs+1
        IF (nk < 1) RETURN
        DO im = 1, nlm
          mp1 = ldc%lm(im)+1
          local(:,:,im,:) = global(llevs:ke,:,mp1,:)
        ENDDO
      ENDIF
    ELSE
      IF (p_pe == p_io) THEN
        DO i = 1, p_nprocs
          pe     = gdc(i)%pe
          llevs  = gdc(i)%llevs
          lleve  = gdc(i)%lleve
          nlm    = gdc(i)%nlm
          nhgl   = gdc(i)%nlat/2
          ke = MIN (lleve,SIZE(global,1))
          nk = ke-llevs+1
          IF (nk < 1) CYCLE
          
          ALLOCATE (sendbuf (nk,2,nlm,nhgl))
          
          DO im = 1, nlm
            mp1 = gdc(i)%lm(im)+1
            sendbuf(:,:,im,:) = global(llevs:ke,:,mp1,:)
          ENDDO
          
          IF (pe /= p_pe) THEN
            CALL p_send( sendbuf, pe, tag_scatter_sa)
          ELSE
            local = sendbuf
          ENDIF
          DEALLOCATE (sendbuf)
        ENDDO
      ELSE
        IF (SIZE(local,1) > 0) THEN
          CALL p_recv (local, p_io, tag_scatter_sa)
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE scatter_sa4
  !-----------------------------------------------------------------------------
  SUBROUTINE scatter_sa2 (global, local, lall_read)
    REAL(dp), POINTER, INTENT(in)  :: global(:,:)
    REAL(dp),          INTENT(out) :: local(:,:)
    LOGICAL,           INTENT(in)  :: lall_read

    REAL(dp), POINTER :: sendbuf(:,:)

    INTEGER :: i, pe, llevs, lleve, nlnm0, nhgl, ke, nk

    IF (lall_read) THEN
      nlnm0 = nlnm0_index(p_pe)
      IF (SIZE(local,1) > 0 .AND. nlnm0 > 0 ) THEN
          llevs  = ldc%llevs
          lleve  = ldc%lleve
          nlnm0  = ldc%nlnm0
          nhgl   = ldc%nlat/2
          ke = MIN (lleve,SIZE(global,1))
          nk = ke-llevs+1
          IF (nk < 1 .OR. nlnm0 < 1) RETURN
          local(:,:) = global(llevs:ke,:)
      ENDIF
    ELSE
      IF (p_pe == p_io) THEN
        DO i = 1, p_nprocs
          pe = gdc(i)%pe
          llevs  = gdc(i)%llevs
          lleve  = gdc(i)%lleve
          nlnm0  = gdc(i)%nlnm0
          nhgl   = gdc(i)%nlat/2
          ke = MIN (lleve,SIZE(global,1))
          nk = ke-llevs+1
          IF (nk < 1 .OR. nlnm0 < 1) CYCLE
          
          ALLOCATE (sendbuf (nk,nhgl))
          
          sendbuf(:,:) = global(llevs:ke,:)
          
          IF (pe /= p_pe) THEN
            CALL p_send( sendbuf, pe, tag_scatter_sa)
          ELSE
            local = sendbuf
          ENDIF
          DEALLOCATE (sendbuf)
        ENDDO
      ELSE
        nlnm0 = nlnm0_index(p_pe)
        IF (SIZE(local,1) > 0 .AND. nlnm0 > 0 ) THEN
          CALL p_recv (local, p_io, tag_scatter_sa)
        ENDIF
      ENDIF
    ENDIF

    CONTAINS

      INTEGER FUNCTION nlnm0_index(lpe)
        INTEGER, INTENT(in) :: lpe
        INTEGER :: i, indx
        indx = -1
        search_index: DO i = 1, SIZE(gdc)
          IF(gdc(i)%pe == lpe) THEN
            indx = i
            EXIT search_index
          ENDIF
        ENDDO search_index
        nlnm0_index = gdc(indx)%nlnm0
      END FUNCTION nlnm0_index

  END SUBROUTINE scatter_sa2

END MODULE mo_tr_scatter
