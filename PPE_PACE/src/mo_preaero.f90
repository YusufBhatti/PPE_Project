MODULE MO_PREAERO

  USE mo_ham_m7_trac,     ONLY: idt_nas, idt_ms4as, idt_mbcas, &
                		idt_ncs, idt_ms4cs, idt_mbccs, &
                		idt_mduas, idt_mducs, idt_mocas, &
                                idt_nks, idt_ms4ks, idt_mduai, &
                                idt_nai, idt_mbcks, idt_mduci, &
                                idt_nci
  USE mo_kind,            ONLY: dp
  USE mo_time_control,    ONLY: current_date, get_date_components

  IMPLICIT NONE

  public :: prescribe_aerosols_mpace

  contains

SUBROUTINE prescribe_aerosols_mpace(&
                       !--IN
                       kproma, kbdim, klev, ktrac, method,                                      &
                       !--INOUT
                       pxtm1, pxtte)
  INTEGER, INTENT(in)                                  :: kproma, kbdim, klev, ktrac, method
  REAL(dp), INTENT(inout), DIMENSION(kbdim,klev,ktrac) :: pxtm1, pxtte


  REAL(dp) :: zaermin, zaermassmin, zaermincs, zaermassmincs
  REAL(dp) :: fact
  INTEGER  :: ithismonth, ithisday, ithishour, ithissecond,                                     &
              ithisyear, ithisminute, ibltop
  INTEGER  :: jk
  !Baustelle
  !>>LI [Prescribe aerosols, Included to prescribe frist day aerosol concentration for the idealised SCM study]
  !  IF (lmpace) THEN      !eventuell Switch noch einfügen
  !     print*,lmpace, '<-------------lmpace'
  !Möglichkeit auswählen:
  !1: Aerosole am ersten Tag vorgeben, MPACE Konzentration
  !2: Aerosole jeden Tag, MPACE Konzentration
  !3: Aersole am ersten Tag vorgeben, MPACE*10
  !4: Aerosole jeden Tag, MPACE Konzentration*10
  !5: Aerosole jeden Tag, Zahlen in Anlehnung an den globalen Lauf
  !6: Keine Aerosole vorschreiben
  !7: Aerosole jeden Tag, MPACE Konzentration*100 -> Scenario 5
  !8: Aerosole jeden Tag, Zahlen in Anlehnung an den globalen Lauf, nur Mineral Dust
  !9: Aerosole jeden Tag, MPACE Konzentration, nur Mineral Dust
  !10: ...
  !11: Remos case
  select case(method)
  case(1)
     zaermin= 72.2E6       !Minimum aerosol number in m-3 (MPACE HHPC-6 mes.)
     zaermassmin= 5.E-12   !Minimum aerosol mass in kg m-3 for 50nm particles in acc.mode
     zaermincs= 1.8E6      !Minimum aerosol number  in m-3 for particles in coarse mode
     zaermassmincs= 6.5E-9 !Minimum aerosol mass in kg m-3 for 1.3µm particles in coarse mode

     !only for first day
     CALL get_date_components (current_date,  month=ithismonth, year=ithisyear,  &
          day=ithisday,  hour=ithishour,   minute=ithisminute, &
          second=ithissecond                                   )
     IF (ithismonth .EQ. 01 .AND. ithisday .EQ. 01) THEN
        ibltop = 30
        DO jk=klev,ibltop,-1
           pxtm1(:,jk,idt_nas) =zaermin
           pxtte(:,jk,idt_nas) = 0._dp
           pxtm1(:,jk,idt_ms4as)=0.7*zaermassmin
           pxtte(:,jk,idt_ms4as)=0._dp
           pxtm1(:,jk,idt_mbcas)=0.3*zaermassmin
           pxtte(:,jk,idt_mbcas)=0._dp

           pxtm1(:,jk,idt_ncs) = zaermincs
           pxtte(:,jk,idt_ncs) = 0._dp
           pxtm1(:,jk,idt_ms4cs) = 0.7*zaermassmincs
           pxtte(:,jk,idt_ms4cs) = 0._dp
           pxtm1(:,jk,idt_mbccs) = 0.3*zaermassmincs
           pxtte(:,jk,idt_mbccs) = 0._dp
        END DO
     END IF !first day

  case(2)
     zaermin= 72.2E6       !Minimum aerosol number in m-3 (MPACE HHPC-6 mes.)
     zaermassmin= 5.0E-12  !Minimum aerosol mass in kg m-3 for 50nm particles in acc.mode
     zaermincs= 1.8E6      !Minimum aerosol number  in m-3 for particles in coarse mode
     zaermassmincs= 6.5E-9 !Minimum aerosol mass in kg m-3 for 1.3µm particles in coarse mode

     !only for first day
     !     CALL get_date_components (current_date,  month=ithismonth, year=ithisyear,  &
     !          day=ithisday,  hour=ithishour,   minute=ithisminute, &
     !          second=ithissecond                                   )
     !     IF (ithismonth .EQ. 01 .AND. ithisday .EQ. 01) THEN
     ibltop = 30
     DO jk=klev,ibltop,-1
        pxtm1(:,jk,idt_nas) =zaermin
        pxtte(:,jk,idt_nas) = 0._dp
        pxtm1(:,jk,idt_ms4as)=0.7*zaermassmin
        pxtte(:,jk,idt_ms4as)=0._dp
        pxtm1(:,jk,idt_mbcas)=0.3*zaermassmin
        pxtte(:,jk,idt_mbcas)=0._dp

        pxtm1(:,jk,idt_ncs) = zaermincs
        pxtte(:,jk,idt_ncs) = 0._dp
        pxtm1(:,jk,idt_ms4cs) = 0.7*zaermassmincs
        pxtte(:,jk,idt_ms4cs) = 0._dp
        pxtm1(:,jk,idt_mbccs) = 0.3*zaermassmincs
        pxtte(:,jk,idt_mbccs) = 0._dp
     END DO
     !     END IF !first day

  case(3)
     zaermin= 722.E6 !72.2E6       !Minimum aerosol number in m-3 (MPACE HHPC-6 mes.)
     zaermassmin= 50.E-12 !5.E-12 !Minimum aerosol mass in kg m-3 for 50nm particles in acc.mode
     zaermincs= 18.E6 !1.8E6      !Minimum aerosol number  in m-3 for particles in coarse mode
     zaermassmincs= 65.E-9 !6.5E-9!Minimum aerosol mass in kg m-3 for 1.3µm particles in coarse mode

     !only for first day
     CALL get_date_components (current_date,  month=ithismonth, year=ithisyear,  &
          day=ithisday,  hour=ithishour,   minute=ithisminute, &
          second=ithissecond                                   )
     IF (ithismonth .EQ. 01 .AND. ithisday .EQ. 01) THEN
        ibltop = 30
        DO jk=klev,ibltop,-1
           pxtm1(:,jk,idt_nas) =zaermin
           pxtte(:,jk,idt_nas) = 0._dp
           pxtm1(:,jk,idt_ms4as)=0.7*zaermassmin
           pxtte(:,jk,idt_ms4as)=0._dp
           pxtm1(:,jk,idt_mbcas)=0.3*zaermassmin
           pxtte(:,jk,idt_mbcas)=0._dp

           pxtm1(:,jk,idt_ncs) = zaermincs
           pxtte(:,jk,idt_ncs) = 0._dp
           pxtm1(:,jk,idt_ms4cs) = 0.7*zaermassmincs
           pxtte(:,jk,idt_ms4cs) = 0._dp
           pxtm1(:,jk,idt_mbccs) = 0.3*zaermassmincs
           pxtte(:,jk,idt_mbccs) = 0._dp
        END DO
     END IF !first day

  case(4)
     zaermin= 722.E6 !72.2E6       !Minimum aerosol number in m-3 (MPACE HHPC-6 mes.)
     zaermassmin= 50.E-12 !5.E-12 !Minimum aerosol mass in kg m-3 for 50nm particles in acc.mode
     zaermincs= 18.E6 !1.8E6      !Minimum aerosol number  in m-3 for particles in coarse mode
     zaermassmincs= 65.E-9 !6.5E-9!Minimum aerosol mass in kg m-3 for 1.3µm particles in coarse mode

     !only for first day
     !     CALL get_date_components (current_date,  month=ithismonth, year=ithisyear,  &
     !          day=ithisday,  hour=ithishour,   minute=ithisminute, &
     !          second=ithissecond                                   )
     !     IF (ithismonth .EQ. 01 .AND. ithisday .EQ. 01) THEN
     ibltop = 30
     DO jk=klev,ibltop,-1
        pxtm1(:,jk,idt_nas) =zaermin
        pxtte(:,jk,idt_nas) = 0._dp
        pxtm1(:,jk,idt_ms4as)=0.7*zaermassmin
        pxtte(:,jk,idt_ms4as)=0._dp
        pxtm1(:,jk,idt_mbcas)=0.3*zaermassmin
        pxtte(:,jk,idt_mbcas)=0._dp

        pxtm1(:,jk,idt_ncs) = zaermincs
        pxtte(:,jk,idt_ncs) = 0._dp
        pxtm1(:,jk,idt_ms4cs) = 0.7*zaermassmincs
        pxtte(:,jk,idt_ms4cs) = 0._dp
        pxtm1(:,jk,idt_mbccs) = 0.3*zaermassmincs
        pxtte(:,jk,idt_mbccs) = 0._dp
     END DO
     !     END IF !first day

  case(5)
     zaermin= 4E7        !Minimum aerosol number in acc mode 
     zaermassmin= 3.7E-10 !4.5E-10  !Minimum aerosol mass in acc mode
     zaermincs= 3E4 !3.5E3       !Minimum aerosol number  in coarse mode
     zaermassmincs= 1E-11 !1.5E-11 !Minimum aerosol mass in coarse mode

     ibltop = 30
     DO jk=klev,ibltop,-1
        pxtm1(:,jk,idt_nas) =zaermin
        pxtte(:,jk,idt_nas) = 0._dp
        pxtm1(:,jk,idt_ms4as)=0.81*zaermassmin !0.9*zaermassmin
        pxtte(:,jk,idt_ms4as)=0._dp
        pxtm1(:,jk,idt_mbcas)=0.03*zaermassmin !0.07*zaermassmin
        pxtte(:,jk,idt_mbcas)=0._dp
        pxtm1(:,jk,idt_mocas)=0.14*zaermassmin 
        pxtte(:,jk,idt_mocas)=0._dp
        pxtm1(:,jk,idt_mduas)=0.02*zaermassmin !0.03*zaermassmin
        pxtte(:,jk,idt_mduas)=0._dp

        pxtm1(:,jk,idt_ncs) = zaermincs
        pxtte(:,jk,idt_ncs) = 0._dp
        pxtm1(:,jk,idt_ms4cs) = 0.04*zaermassmincs !0.05*zaermassmincs
        pxtte(:,jk,idt_ms4cs) = 0._dp
        !        pxtm1(:,jk,idt_mbccs) = 0.01*zaermassmincs
        !        pxtte(:,jk,idt_mbccs) = 0._dp
        pxtm1(:,jk,idt_mducs)=0.96*zaermassmincs
        pxtte(:,jk,idt_mducs)=0._dp
     END DO

  case(6)

  case(7)
     zaermin= 7220.E6 !72.2E6       !Minimum aerosol number in m-3 (MPACE HHPC-6 mes.)
     zaermassmin= 500.E-12 !5.E-12 !Minimum aerosol mass in kg m-3 for 50nm particles in acc.mode
     zaermincs= 180.E6 !1.8E6      !Minimum aerosol number  in m-3 for particles in coarse mode
     zaermassmincs= 650.E-9 !6.5E-9!Minimum aerosol mass in kg m-3 for 1.3µm particles in coarse mode

     !only for first day
     !     CALL get_date_components (current_date,  month=ithismonth, year=ithisyear,  &
     !          day=ithisday,  hour=ithishour,   minute=ithisminute, &
     !          second=ithissecond                                   )
     !     IF (ithismonth .EQ. 01 .AND. ithisday .EQ. 01) THEN
     ibltop = 30
     DO jk=klev,ibltop,-1
        pxtm1(:,jk,idt_nas) =zaermin
        pxtte(:,jk,idt_nas) = 0._dp
        pxtm1(:,jk,idt_ms4as)=0.7*zaermassmin
        pxtte(:,jk,idt_ms4as)=0._dp
        pxtm1(:,jk,idt_mbcas)=0.3*zaermassmin
        pxtte(:,jk,idt_mbcas)=0._dp

        pxtm1(:,jk,idt_ncs) = zaermincs
        pxtte(:,jk,idt_ncs) = 0._dp
        pxtm1(:,jk,idt_ms4cs) = 0.7*zaermassmincs
        pxtte(:,jk,idt_ms4cs) = 0._dp
        pxtm1(:,jk,idt_mbccs) = 0.3*zaermassmincs
        pxtte(:,jk,idt_mbccs) = 0._dp
     END DO
     !     END IF !first day
     !<<LI

     !>>CZ [only mineral dust for heterogeneous freezing from Niemand et. al. 2012]
  case(8)
     zaermin       = 4E7_dp     !Minimum aerosol number in acc mode 
     zaermassmin   = 3.7E-10_dp !Minimum aerosol mass in acc mode
     zaermincs     = 3E4_dp     !Minimum aerosol number  in coarse mode
     zaermassmincs = 1E-11_dp   !Minimum aerosol mass in coarse mode

     ibltop = 1
     DO jk=klev,ibltop,-1
        pxtm1(:,jk,idt_nas)   = zaermin
        pxtte(:,jk,idt_nas)   = 0._dp
        pxtm1(:,jk,idt_ms4as) = 0._dp*zaermassmin
        pxtte(:,jk,idt_ms4as) = 0._dp
        pxtm1(:,jk,idt_mbcas) = 0._dp*zaermassmin
        pxtte(:,jk,idt_mbcas) = 0._dp
        pxtm1(:,jk,idt_mocas) = 0._dp*zaermassmin 
        pxtte(:,jk,idt_mocas) = 0._dp
        pxtm1(:,jk,idt_mduas) = 1._dp*zaermassmin
        pxtte(:,jk,idt_mduas) = 0._dp
        pxtm1(:,jk,idt_ncs)   = zaermincs
        pxtte(:,jk,idt_ncs)   = 0._dp
        pxtm1(:,jk,idt_ms4cs) = 0._dp*zaermassmincs
        pxtte(:,jk,idt_ms4cs) = 0._dp
        pxtm1(:,jk,idt_mducs) = 1._dp*zaermassmincs
        pxtte(:,jk,idt_mducs) = 0._dp
     END DO
     !<<CZ

  case(9)
     zaermin= 72.2E6       !Minimum aerosol number in m-3 (MPACE HHPC-6 mes.)
     zaermassmin= 5.0E-12  !Minimum aerosol mass in kg m-3 for 50nm particles in acc.mode
     zaermincs= 1.8E6      !Minimum aerosol number  in m-3 for particles in coarse mode
     zaermassmincs= 6.5E-9 !Minimum aerosol mass in kg m-3 for 1.3µm particles in coarse mode

     ibltop = 1
     DO jk=klev,ibltop,-1     
        !                pxtm1(:,klev,idt_nas) =zaermin
        !                pxtte(:,klev,idt_nas) = 0._dp
        !                pxtm1(:,klev,idt_ms4as)=0.0*zaermassmin
        !                pxtte(:,klev,idt_ms4as)=0._dp
        !                pxtm1(:,klev,idt_mbcas)=0.0*zaermassmin
        !                pxtte(:,klev,idt_mbcas)=0._dp
        !                pxtm1(:,klev,idt_mduas) = 1._dp*zaermassmin
        !                pxtte(:,klev,idt_mduas) = 0._dp
        !
        !                pxtm1(:,klev,idt_ncs) = zaermincs
        !                pxtte(:,klev,idt_ncs) = 0._dp
        !                pxtm1(:,klev,idt_ms4cs) = 0.0*zaermassmincs
        !                pxtte(:,klev,idt_ms4cs) = 0._dp
        !                pxtm1(:,klev,idt_mbccs) = 0.0*zaermassmincs
        !                pxtte(:,klev,idt_mbccs) = 0._dp
        !                pxtm1(:,klev,idt_mducs) = 1._dp*zaermassmincs
        !                pxtte(:,klev,idt_mducs) = 0._dp

        pxtm1(:,jk,idt_nas) =zaermin
        pxtte(:,jk,idt_nas) = 0._dp
        pxtm1(:,jk,idt_ms4as)=0.0*zaermassmin
        pxtte(:,jk,idt_ms4as)=0._dp
        pxtm1(:,jk,idt_mbcas)=0.0*zaermassmin
        pxtte(:,jk,idt_mbcas)=0._dp
        pxtm1(:,jk,idt_mduas) = 1._dp*zaermassmin
        pxtte(:,jk,idt_mduas) = 0._dp

        pxtm1(:,jk,idt_ncs) = zaermincs
        pxtte(:,jk,idt_ncs) = 0._dp
        pxtm1(:,jk,idt_ms4cs) = 0.0*zaermassmincs
        pxtte(:,jk,idt_ms4cs) = 0._dp
        pxtm1(:,jk,idt_mbccs) = 0.0*zaermassmincs
        pxtte(:,jk,idt_mbccs) = 0._dp
        pxtm1(:,jk,idt_mducs) = 1._dp*zaermassmincs
        pxtte(:,jk,idt_mducs) = 0._dp
     END DO

  case(10)
     zaermin= 72.2E6       !Minimum aerosol number in m-3 (MPACE HHPC-6 mes.)
     zaermassmin= 5.0E-12  !Minimum aerosol mass in kg m-3 for 50nm particles in acc.mode
     zaermincs= 1.8E6      !Minimum aerosol number  in m-3 for particles in coarse mode
     zaermassmincs= 6.5E-9 !Minimum aerosol mass in kg m-3 for 1.3µm particles in coarse mode

     ibltop = 1
     DO jk=klev,ibltop,-1     
        pxtm1(:,jk,idt_nas) =zaermin
        pxtte(:,jk,idt_nas) = 0._dp
        pxtm1(:,jk,idt_ms4as)=0.2_dp*zaermassmin
        pxtte(:,jk,idt_ms4as)=0._dp
        pxtm1(:,jk,idt_mbcas)=0.0*zaermassmin
        pxtte(:,jk,idt_mbcas)=0._dp
        pxtm1(:,jk,idt_mduas) = 0.8_dp*zaermassmin
        pxtte(:,jk,idt_mduas) = 0._dp
     END DO

  case(11)
     zaermin= 72.2E6       !Minimum aerosol number in m-3 (MPACE HHPC-6 mes.)
     zaermassmin= 5.0E-12  !Minimum aerosol mass in kg m-3 for 50nm particles in acc.mode
     zaermincs= 1.8E6      !Minimum aerosol number  in m-3 for particles in coarse mode
     zaermassmincs= 6.5E-9 !Minimum aerosol mass in kg m-3 for 1.3µm particles in coarse mode

     pxtm1(:,:,idt_nas) =zaermin
     pxtte(:,:,idt_nas) = 0._dp
     pxtm1(:,:,idt_ms4as)=0.7*zaermassmin
     pxtte(:,:,idt_ms4as)=0._dp
     pxtm1(:,:,idt_mbcas)=0.3*zaermassmin
     pxtte(:,:,idt_mbcas)=0._dp

     pxtm1(:,:,idt_ncs) = zaermincs
     pxtte(:,:,idt_ncs) = 0._dp
     pxtm1(:,:,idt_ms4cs) = 0.7*zaermassmincs
     pxtte(:,:,idt_ms4cs) = 0._dp
     pxtm1(:,:,idt_mbccs) = 0.3*zaermassmincs
     pxtte(:,:,idt_mbccs) = 0._dp

  case(12)
     fact = 0.001

     ! pxtm1(1:kproma,:,idt_nks) = 1e8_dp
     ! ! pxtm1(1:kproma,:,idt_ms4ks) = 1e-10_dp

     ! pxtm1(1:kproma,:,idt_nas) = 1e6_dp*fact
     ! pxtm1(1:kproma,:,idt_mduas) = 1e-12_dp*fact

     ! pxtm1(1:kproma,:,idt_mduai) = 1e-11_dp*fact
     ! pxtm1(1:kproma,:,idt_nai) = 1e5_dp*fact

     ! pxtm1(1:kproma,:,idt_mbcks) = 1e-10_dp

     pxtm1(1:kproma,:,idt_nai) = 0._dp*fact
     pxtm1(1:kproma,:,idt_nci) = 1e5_dp*fact

     pxtm1(1:kproma,:,idt_mduai) = 0._dp*fact
     pxtm1(1:kproma,:,idt_mduci) = 1e-12_dp*fact

  case default

  end select

END SUBROUTINE prescribe_aerosols_mpace

END MODULE mo_preaero
