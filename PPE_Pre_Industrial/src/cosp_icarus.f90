!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
      SUBROUTINE icarus(                      &
!cms reduced number of continuation lines
     &     debug, debugcol, npoints, sunlit,  &
     &     nlev, ncol, pfull, phalf,qv,         & 
     &     cc, dtau_s, top_height,  &
     &     top_height_direction, overlap,   &
     &     frac_out, skt, emsfc_lw,   &
     &     at,  &
     &     dem_s,  &
     &     fq_isccp, &
     &     totalcldarea, &
     &     meanptop, &
     &     meantaucld, &
     &     meanalbedocld, &
     &     meantb, &
     &     meantbclr, &
     &     boxtau, &
     &     boxptop &
     &)

     USE mo_kind,  ONLY: wp

!$Id: icarus.f,v 4.1 2010/05/27 16:30:18 hadmw Exp $

! *****************************COPYRIGHT****************************
! (c) 2009, Lawrence Livermore National Security Limited Liability 
! Corporation.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the
! following conditions are met:
! 
!     * Redistributions of source code must retain the above 
!       copyright  notice, this list of conditions and the following 
!       disclaimer.
!     * Redistributions in binary form must reproduce the above 
!       copyright notice, this list of conditions and the following 
!       disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Lawrence Livermore National Security 
!       Limited Liability Corporation nor the names of its 
!       contributors may be used to endorse or promote products
!       derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
! 
! *****************************COPYRIGHT*******************************

      IMPLICIT NONE

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER :: ncolprint
      
!     -----
!     Input 
!     -----

      INTEGER :: npoints       !  number of model points in the horizontal
      INTEGER :: nlev          !  number of model levels in column
      INTEGER :: ncol          !  number of subcolumns

      INTEGER :: sunlit(npoints) !  1 for day points, 0 for night time

      REAL(wp) :: pfull(npoints,nlev)
                       !  pressure of full model levels (Pascals)
                  !  pfull(npoints,1) is top level of model
                  !  pfull(npoints,nlev) is bot of model

      REAL(wp) :: phalf(npoints,nlev+1)
                  !  pressure of half model levels (Pascals)
                  !  phalf(npoints,1) is top of model
                  !  phalf(npoints,nlev+1) is the surface pressure

      REAL(wp) :: qv(npoints,nlev)
                  !  water vapor specific humidity (kg vapor/ kg air)
                  !         on full model levels

      REAL(wp) :: cc(npoints,nlev)   
                  !  input cloud cover in each model level (fraction) 
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by clouds


      REAL(wp) :: dtau_s(npoints,nlev) 
                  !  mean 0.67 micron optical depth of stratiform
                !  clouds in each model level
                  !  NOTE:  this the cloud optical depth of only the
                  !  cloudy part of the grid box, it is not weighted
                  !  with the 0 cloud optical depth of the clear
                  !         part of the grid box


      INTEGER :: overlap                   !  overlap type
                              !  1=max
                              !  2=rand
                              !  3=max/rand

      INTEGER :: top_height                !  1 = adjust top height using both a computed
                                        !  infrared brightness temperature and the visible
                              !  optical depth to adjust cloud top pressure. Note
                              !  that this calculation is most appropriate to compare
                              !  to ISCCP data during sunlit hours.
                                        !  2 = do not adjust top height, that is cloud top
                                        !  pressure is the actual cloud top pressure
                                        !  in the model
                              !  3 = adjust top height using only the computed
                              !  infrared brightness temperature. Note that this
                              !  calculation is most appropriate to compare to ISCCP
                              !  IR only algortihm (i.e. you can compare to nighttime
                              !  ISCCP data with this option)

      INTEGER :: top_height_direction ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 !
                                 ! 1 = find the *lowest* altitude (highest pressure) level
                                 ! with interpolated temperature equal to the radiance
                                 ! determined cloud-top temperature
                                 !
                                 ! 2 = find the *highest* altitude (lowest pressure) level
                                 ! with interpolated temperature equal to the radiance 
                                 ! determined cloud-top temperature
                                 ! 
                                 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
                                 !                                 !
                                 ! 1 = old setting: matches all versions of 
                                 ! ISCCP simulator with versions numbers 3.5.1 and lower
                                 !
                                 ! 2 = default setting: for version numbers 4.0 and higher
!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL(wp) :: skt(npoints)                 !  skin Temperature (K)
      REAL(wp) :: emsfc_lw                     !  10.5 micron emissivity of surface (fraction)
      REAL(wp) :: at(npoints,nlev)             !  temperature in each model level (K)
      REAL(wp) :: dem_s(npoints,nlev)          !  10.5 micron longwave emissivity of stratiform
                                               !  clouds in each
                                               !  model level.  Same note applies as in dtau_s.

      REAL(wp) :: frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column



!     ------
!     Output
!     ------

      REAL(wp) :: fq_isccp(npoints,7,7)        !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types

      REAL(wp) :: totalcldarea(npoints)        !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  NOTE: This diagnostic
                                        ! does not count model clouds with tau < isccp_taumin
                              ! Thus this diagnostic does not equal the sum over all entries of 
                              !  fq_isccp. However, this diagnostic does equal the sum over entries
                              !   of fq_isccp with itau = 2:7 (omitting itau = 1)
      
      
     ! The following three means are averages only over the cloudy areas with tau > isccp_taumin.  
     ! If no clouds with tau > isccp_taumin are in grid box all three quantities should equal zero.      
                              
      REAL(wp) :: meanptop(npoints)            !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
                              
      REAL(wp) :: meantaucld(npoints)          !  mean optical thickness 
                                        !  linear averaging in albedo performed.
      
      REAL(wp) :: meanalbedocld(npoints)        ! mean cloud albedo
                                        ! linear averaging in albedo performed
                                        
      REAL(wp) :: meantb(npoints)              ! mean all-sky 10.5 micron brightness temperature
      
      REAL(wp) :: meantbclr(npoints)           ! mean clear-sky 10.5 micron brightness temperature
      
      REAL(wp) :: boxtau(npoints,ncol)         !  optical thickness in each column
      
      REAL(wp) :: boxptop(npoints,ncol)        !  cloud top pressure (mb) in each column
                              
                                                                                          
!
!     ------
!     Working variables added when program updated to mimic Mark Webb's PV-Wave code
!     ------

      REAL(wp) :: dem(npoints,ncol),bb(npoints)     !  working variables for 10.5 micron longwave 
                              !  emissivity in part of
                              !  gridbox under consideration

      REAL(wp) :: ptrop(npoints)
      REAL(wp) :: attrop(npoints)
      REAL(wp) :: attropmin (npoints)
      REAL(wp) :: atmax(npoints)
      REAL(wp) :: btcmin(npoints)
      REAL(wp) :: transmax(npoints)

      INTEGER :: i,j,ilev,ibox,itrop(npoints)
      INTEGER :: ipres(npoints)
      INTEGER :: itau(npoints),ilev2
      INTEGER :: acc(nlev,ncol)
      INTEGER :: match(npoints,nlev-1)
      INTEGER :: nmatch(npoints)
      INTEGER :: levmatch(npoints,ncol)
      
      !variables needed for water vapor continuum absorption
      REAL(wp) :: fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
      REAL(wp) :: taumin(npoints)
      REAL(wp) :: dem_wv(npoints,nlev), wtmair, wtmh20, Navo, grav, pstd, t0
      REAL(wp) :: press(npoints), dpress(npoints), atmden(npoints)
      REAL(wp) :: rvh20(npoints), wk(npoints), rhoave(npoints)
      REAL(wp) :: rh20s(npoints), rfrgn(npoints)
      REAL(wp) :: tmpexp(npoints),tauwv(npoints)
      
      CHARACTER*1 :: cchar(6),cchar_realtops(6)
      INTEGER :: icycle
      REAL(wp) :: tau(npoints,ncol)
      LOGICAL :: box_cloudy(npoints,ncol)
      REAL(wp) :: tb(npoints,ncol)
      REAL(wp) :: ptop(npoints,ncol)
      REAL(wp) :: emcld(npoints,ncol)
      REAL(wp) :: fluxtop(npoints,ncol)
      REAL(wp) :: trans_layers_above(npoints,ncol)
      REAL(wp) :: isccp_taumin,fluxtopinit(npoints),tauir(npoints)
      REAL(wp) :: albedocld(npoints,ncol)
      REAL(wp) :: boxarea
      INTEGER :: debug       ! set to non-zero value to print out inputs
                    ! with step debug
      INTEGER :: debugcol    ! set to non-zero value to print out column
                    ! decomposition with step debugcol
      INTEGER :: rangevec(npoints),rangeerror

      INTEGER :: k1,k2
      REAL(wp) :: rec2p13,tauchk,logp,logp1,logp2,atd
      REAL(wp) :: output_missing_value

      CHARACTER*10 :: ftn09
      
      DATA isccp_taumin / 0.3_wp /
!cms      DATA output_missing_value / -1.E+30 /
      DATA output_missing_value / 0._wp /
      DATA cchar / ' ','-','1','+','I','+'/
      DATA cchar_realtops / ' ',' ','1','1','I','I'/

!     ------ End duplicate definitions common to wrapper routine

      tauchk = -1._wp*log(0.9999999_wp)
      rec2p13=1._wp/2.13_wp

      ncolprint=0

      IF ( debug.ne.0 ) THEN
          j=1
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write(6,'(a10)') 'debug='
          write(6,'(8I10)') debug
          write(6,'(a10)') 'debugcol='
          write(6,'(8I10)') debugcol
          write(6,'(a10)') 'npoints='
          write(6,'(8I10)') npoints
          write(6,'(a10)') 'nlev='
          write(6,'(8I10)') nlev
          write(6,'(a10)') 'ncol='
          write(6,'(8I10)') ncol
          write(6,'(a11)') 'top_height='
          write(6,'(8I10)') top_height
          write(6,'(a21)') 'top_height_direction='
          write(6,'(8I10)') top_height_direction
          write(6,'(a10)') 'overlap='
          write(6,'(8I10)') overlap
          write(6,'(a10)') 'emsfc_lw='
          write(6,'(8f10.2)') emsfc_lw
        DO j=1,npoints,debug
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write(6,'(a10)') 'sunlit='
          write(6,'(8I10)') sunlit(j)
          write(6,'(a10)') 'pfull='
          write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
          write(6,'(a10)') 'phalf='
          write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
          write(6,'(a10)') 'qv='
          write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
          write(6,'(a10)') 'cc='
          write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
          write(6,'(a10)') 'dtau_s='
          write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
          write(6,'(a10)') 'skt='
          write(6,'(8f10.2)') skt(j)
          write(6,'(a10)') 'at='
          write(6,'(8f10.2)') (at(j,i),i=1,nlev)
          write(6,'(a10)') 'dem_s='
          write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
        ENDDO
      ENDIF

!     ---------------------------------------------------!

      IF (ncolprint.ne.0) THEN
      DO j=1,npoints,1000
        write(6,'(a10)') 'j='
        write(6,'(8I10)') j
      ENDDO
      ENDIF

      IF (top_height .eq. 1 .or. top_height .eq. 3) THEN 

      DO j=1,npoints 
          ptrop(j)=5000._wp
          attropmin(j) = 400._wp
          atmax(j) = 0._wp
          attrop(j) = 120._wp
          itrop(j) = 1
      ENDDO 

      DO 12 ilev=1,nlev
        DO j=1,npoints 
         IF (pfull(j,ilev) .lt. 40000._wp .and.   &
     &          pfull(j,ilev) .gt.  5000._wp .and. &
     &          at(j,ilev) .lt. attropmin(j)) THEN
                ptrop(j) = pfull(j,ilev)
                attropmin(j) = at(j,ilev)
                attrop(j) = attropmin(j)
                itrop(j)=ilev
           ENDIF
        ENDDO
12    CONTINUE

      DO 13 ilev=1,nlev
        DO j=1,npoints 
           IF (at(j,ilev) .gt. atmax(j) .and. &
     &              ilev  .ge. itrop(j)) atmax(j)=at(j,ilev)
        ENDDO
13    CONTINUE

      ENDIF


      IF (top_height .eq. 1 .or. top_height .eq. 3) THEN
          DO j=1,npoints
              meantb(j) = 0._wp
              meantbclr(j) = 0._wp 
          ENDDO
      ELSE
          DO j=1,npoints
              meantb(j) = output_missing_value
              meantbclr(j) = output_missing_value
          ENDDO
      ENDIF
      
!     -----------------------------------------------------!

!     ---------------------------------------------------!

      DO ilev=1,nlev
        DO j=1,npoints

          rangevec(j)=0

          IF (cc(j,ilev) .lt. 0._wp .or. cc(j,ilev) .gt. 1._wp) THEN
!           error = cloud fraction less than zero
!           error = cloud fraction greater than 1
            rangevec(j)=rangevec(j)+1
          ENDIF


          IF (dtau_s(j,ilev) .lt. 0._wp) THEN
!           ' error = stratiform cloud opt. depth less than zero'
            rangevec(j)=rangevec(j)+4
          ENDIF


          IF (dem_s(j,ilev) .lt. 0._wp .or. dem_s(j,ilev) .gt. 1._wp) THEN
!             ' error = stratiform cloud emissivity less than zero'
!             ' error = stratiform cloud emissivity greater than 1'
            rangevec(j)=rangevec(j)+16
          ENDIF
        ENDDO 

        rangeerror=0
        DO j=1,npoints
            rangeerror=rangeerror+rangevec(j)
        ENDDO

        IF (rangeerror.ne.0) THEN 
              write (6,*) 'Input variable out of range'
              write (6,*) 'rangevec:'
              write (6,*) rangevec
!              call flush(6)
              STOP
        ENDIF
      ENDDO
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !initialize tau and albedocld to zero
      DO 15 ibox=1,ncol
        DO j=1,npoints 
            tau(j,ibox)=0._wp
          albedocld(j,ibox)=0._wp
          boxtau(j,ibox)=output_missing_value
          boxptop(j,ibox)=output_missing_value
          box_cloudy(j,ibox)=.false.
        ENDDO
15    CONTINUE

      !compute total cloud optical depth for each column     
      DO ilev=1,nlev
            !increment tau for each of the boxes
            DO ibox=1,ncol
              DO j=1,npoints 
                 IF (frac_out(j,ibox,ilev).eq.1._wp) THEN
                        tau(j,ibox)=tau(j,ibox) &
     &                     + dtau_s(j,ilev)
                 ENDIF
              ENDDO
            ENDDO ! ibox
      ENDDO ! ilev
          IF (ncolprint.ne.0) THEN

              DO j=1,npoints ,1000
                write(6,'(a10)') 'j='
                write(6,'(8I10)') j
                write(6,'(i2,1X,8(f7.2,1X))')  &
     &          ilev,                          &
     &          (tau(j,ibox),ibox=1,ncolprint)
              ENDDO
          ENDIF 
!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done IF top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
!             sky versions of these quantities.

      IF (top_height .eq. 1 .or. top_height .eq. 3) THEN
        !----------------------------------------------------------------------
        !    
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1, 
        !or 10.47 microns 
        wtmair = 28.9644_wp
        wtmh20 = 18.01534_wp
        Navo = 6.023E+23_wp
        grav = 9.806650E+02_wp
        pstd = 1.013250E+06_wp
        t0 = 296._wp
        IF (ncolprint .ne. 0)  &
     &         write(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
        DO 125 ilev=1,nlev
          DO j=1,npoints 
               !press and dpress are dyne/cm2 = Pascals *10
               press(j) = pfull(j,ilev)*10._wp
               dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10._wp
               !atmden = g/cm2 = kg/m2 / 10 
               atmden(j) = dpress(j)/grav
               rvh20(j) = qv(j,ilev)*wtmair/wtmh20
               wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
               rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
               rh20s(j) = rvh20(j)*rhoave(j)
               rfrgn(j) = rhoave(j)-rh20s(j)
               tmpexp(j) = exp(-0.02_wp*(at(j,ilev)-t0))
               tauwv(j) = wk(j)*1.e-20_wp*(  &
     &           (0.0224697_wp*rh20s(j)*tmpexp(j)) +  &
     &                (3.41817e-7_wp*rfrgn(j)) )*0.98_wp
               dem_wv(j,ilev) = 1._wp - exp( -1._wp * tauwv(j))
          ENDDO
               IF (ncolprint .ne. 0) THEN
               DO j=1,npoints ,1000
               write(6,'(a10)') 'j='
               write(6,'(8I10)') j
               write(6,'(i2,1X,3(f8.3,3X))') ilev, &
     &           qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))/(grav/100._wp), &
     &           tauwv(j),dem_wv(j,ilev)
               ENDDO
             ENDIF
125     CONTINUE

        !initialize variables
        DO j=1,npoints 
          fluxtop_clrsky(j) = 0._wp
          trans_layers_above_clrsky(j)=1._wp
        ENDDO

        DO ilev=1,nlev
          DO j=1,npoints 
 
            ! Black body emission at temperature of the layer

              bb(j)=1._wp / ( exp(1307.27_wp/at(j,ilev)) - 1._wp )
              !bb(j)= 5.67e-8*at(j,ilev)**4

              ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above

                fluxtop_clrsky(j) = fluxtop_clrsky(j)  &
     &            + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j) 
            
                ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop

                trans_layers_above_clrsky(j)= &
     &            trans_layers_above_clrsky(j)*(1._wp-dem_wv(j,ilev))
          ENDDO   
            IF (ncolprint.ne.0) THEN
             DO j=1,npoints ,1000
              write(6,'(a10)') 'j='
              write(6,'(8I10)') j
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
              write (6,'(a)') &
     &        'emiss_layer,100._wp*bb(j),100._wp*f,total_trans:'
              write (6,'(4(f7.2,1X))') dem_wv(j,ilev),100._wp*bb(j), &
     &             100._wp*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
             ENDDO   
            ENDIF
        ENDDO   !loop over level
        
        DO j=1,npoints 
          !add in surface emission
          bb(j)=1._wp/( exp(1307.27_wp/skt(j)) - 1._wp )
          !bb(j)=5.67e-8*skt(j)**4

          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j)  &
     &     * trans_layers_above_clrsky(j)
     
          !clear sky brightness temperature
          meantbclr(j) = 1307.27_wp/(log(1._wp+(1._wp/fluxtop_clrsky(j))))
        ENDDO

        IF (ncolprint.ne.0) THEN
        DO j=1,npoints ,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'
          write (6,'(a)') 'emsfc,100._wp*bb(j),100._wp*f,total_trans:'
          write (6,'(5(f7.2,1X))') emsfc_lw,100._wp*bb(j), &
     &      100._wp*fluxtop_clrsky(j), &
     &       trans_layers_above_clrsky(j), meantbclr(j)
        ENDDO
      ENDIF

        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------

        IF (ncolprint.ne.0) THEN

        DO j=1,npoints ,1000
            write(6,'(a10)') 'j='
            write(6,'(8I10)') j
            write (6,'(a)') 'ts:'
            write (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)
            write (6,'(a)') 'ta_rev:'
            write (6,'(8f7.2)') &
     &       ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
        ENDDO
        ENDIF 
        !loop over columns 
        DO ibox=1,ncol
          DO j=1,npoints
            fluxtop(j,ibox)=0._wp
            trans_layers_above(j,ibox)=1._wp
          ENDDO
        ENDDO

        DO ilev=1,nlev
              DO j=1,npoints 
                ! Black body emission at temperature of the layer

              bb(j)=1._wp / ( exp(1307.27_wp/at(j,ilev)) - 1._wp )
              !bb(j)= 5.67e-8*at(j,ilev)**4
              ENDDO

            DO ibox=1,ncol
              DO j=1,npoints 

              ! emissivity for point in this layer
                IF (frac_out(j,ibox,ilev).eq.1._wp) THEN
                dem(j,ibox)= 1._wp - &
     &             ( (1._wp - dem_wv(j,ilev)) * (1._wp -  dem_s(j,ilev)) )
                ELSE
                dem(j,ibox)=  dem_wv(j,ilev)
                ENDIF

                ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above

                fluxtop(j,ibox) = fluxtop(j,ibox) &
     &            + dem(j,ibox) * bb(j) &
     &            * trans_layers_above(j,ibox) 
            
                ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop

                trans_layers_above(j,ibox)= &
     &            trans_layers_above(j,ibox)*(1._wp-dem(j,ibox))

              ENDDO ! j
            ENDDO ! ibox

            IF (ncolprint.ne.0) THEN
              DO j=1,npoints,1000
              write (6,'(a)') 'ilev:'
              write (6,'(I2)') ilev
              write(6,'(a10)') 'j='
              write(6,'(8I10)') j
              write (6,'(a)') 'emiss_layer:'
              write (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
              write (6,'(a)') '100.*bb(j):'
              write (6,'(8f7.2)') (100._wp*bb(j),ibox=1,ncolprint)
              write (6,'(a)') '100.*f:'
              write (6,'(8f7.2)')  &
     &         (100._wp*fluxtop(j,ibox),ibox=1,ncolprint)
              write (6,'(a)') 'total_trans:'
              write (6,'(8f7.2)') &
     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
            ENDDO
          ENDIF

        ENDDO ! ilev

          DO j=1,npoints 
            !add in surface emission
            bb(j)=1._wp/( exp(1307.27_wp/skt(j)) - 1._wp )
            !bb(j)=5.67e-8*skt(j)**4
          ENDDO

        DO ibox=1,ncol
          DO j=1,npoints 

            !add in surface emission

            fluxtop(j,ibox) = fluxtop(j,ibox) &
     &         + emsfc_lw * bb(j) &
     &         * trans_layers_above(j,ibox) 
            
          ENDDO
        ENDDO

        !calculate mean infrared brightness temperature
        DO ibox=1,ncol
          DO j=1,npoints 
            meantb(j) = meantb(j)+1307.27_wp/(log(1._wp+(1._wp/fluxtop(j,ibox))))
          ENDDO
        ENDDO
          DO j=1, npoints
            meantb(j) = meantb(j) / REAL(ncol,wp)
          ENDDO        

        IF (ncolprint.ne.0) THEN

          DO j=1,npoints ,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'id:'
          write (6,'(a)') 'surface'
          write (6,'(a)') 'emiss_layer:'
          write (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
          write (6,'(a)') '100.*bb(j):'
          write (6,'(8f7.2)') (100._wp*bb(j),ibox=1,ncolprint)
          write (6,'(a)') '100.*f:'
          write (6,'(8f7.2)') (100._wp*fluxtop(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'meantb(j):'
          write (6,'(8f7.2)') (meantb(j),ibox=1,ncolprint)
      
          ENDDO
      ENDIF
    
        !now that you have the top of atmosphere radiance account
        !for ISCCP procedures to determine cloud top temperature

        !account for partially transmitting cloud recompute flux 
        !ISCCP would see assuming a single layer cloud
        !note choice here of 2.13, as it is primarily ice
        !clouds which have partial emissivity and need the 
        !adjustment performed in this section
        !
      !If it turns out that the cloud brightness temperature
      !is greater than 260K, then the liquid cloud conversion
        !factor of 2.56 is used.
      !
        !Note that this is discussed on pages 85-87 of 
        !the ISCCP D level documentation (Rossow et al. 1996)
           
          DO j=1,npoints  
            !compute minimum brightness temperature and optical depth
            btcmin(j) = 1._wp /  ( exp(1307.27_wp/(attrop(j)-5._wp)) - 1._wp ) 
          ENDDO 
        DO ibox=1,ncol
          DO j=1,npoints  
            transmax(j) = (fluxtop(j,ibox)-btcmin(j)) &
     &                /(fluxtop_clrsky(j)-btcmin(j))
          !note that the initial setting of tauir(j) is needed so that
          !tauir(j) has a realistic value should the next if block be
          !bypassed
            tauir(j) = tau(j,ibox) * rec2p13
            taumin(j) = -1._wp * log(max(min(transmax(j),0.9999999_wp),0.001_wp))

          ENDDO 

          IF (top_height .eq. 1) THEN
            DO j=1,npoints  
              IF (transmax(j) .gt. 0.001_wp .and. &
     &          transmax(j) .le. 0.9999999_wp) THEN
                fluxtopinit(j) = fluxtop(j,ibox)
              tauir(j) = tau(j,ibox) *rec2p13
              ENDIF
            ENDDO
            DO icycle=1,2
              DO j=1,npoints  
                IF (tau(j,ibox) .gt. (tauchk            )) THEN 
                IF (transmax(j) .gt. 0.001_wp .and. &
     &            transmax(j) .le. 0.9999999_wp) THEN
                  emcld(j,ibox) = 1._wp - exp(-1._wp * tauir(j)  )
                  fluxtop(j,ibox) = fluxtopinit(j) -  &  
     &              ((1._wp-emcld(j,ibox))*fluxtop_clrsky(j))
                  fluxtop(j,ibox)=max(1.E-06_wp, &
     &              (fluxtop(j,ibox)/emcld(j,ibox)))
                  tb(j,ibox)= 1307.27_wp &
     &              / (log(1._wp + (1._wp/fluxtop(j,ibox))))
                  IF (tb(j,ibox) .gt. 260._wp) THEN
                  tauir(j) = tau(j,ibox) / 2.56_wp
                  ENDIF                   
                 ENDIF
                ENDIF
              ENDDO
            ENDDO
                
          ENDIF
        
          DO j=1,npoints
            IF (tau(j,ibox) .gt. (tauchk            )) THEN 
                !cloudy box 
                !NOTE: tb is the cloud-top temperature not infrared brightness temperature 
                !at this point in the code
                tb(j,ibox)= 1307.27_wp/ (log(1._wp + (1._wp/fluxtop(j,ibox))))
                IF (top_height.eq.1.and.tauir(j).lt.taumin(j)) THEN
                         tb(j,ibox) = attrop(j) - 5._wp
                   tau(j,ibox) = 2.13_wp*taumin(j)
                ENDIF
            ELSE
                !clear sky brightness temperature
                tb(j,ibox) = meantbclr(j)
            ENDIF
          ENDDO ! j
        ENDDO ! ibox

        IF (ncolprint.ne.0) THEN

          DO j=1,npoints,1000
          write(6,'(a10)') 'j='
          write(6,'(8I10)') j
          write (6,'(a)') 'attrop:'
          write (6,'(8f7.2)') (attrop(j))
          write (6,'(a)') 'btcmin:'
          write (6,'(8f7.2)') (btcmin(j))
          write (6,'(a)') 'fluxtop_clrsky*100:'
          write (6,'(8f7.2)') &
     &      (100._wp*fluxtop_clrsky(j))
          write (6,'(a)') '100.*f_adj:'
          write (6,'(8f7.2)') (100._wp*fluxtop(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'transmax:'
          write (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)
          write (6,'(a)') 'tau:'
          write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'emcld:'
          write (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)') &
     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'total_emiss:'
          write (6,'(8f7.2)') &
     &        (1.0_wp-trans_layers_above(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'total_trans:'
          write (6,'(8f7.2)') &
     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
          write (6,'(a)') 'ppout:'
          write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
          ENDDO ! j
      ENDIF

      ENDIF
!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!

      !compute cloud top pressure
      DO 30 ibox=1,ncol
        !segregate according to optical thickness
        IF (top_height .eq. 1 .or. top_height .eq. 3) THEN  
          !find level whose temperature
          !most closely matches brightness temperature
          DO j=1,npoints 
            nmatch(j)=0
          ENDDO
          DO 29 k1=1,nlev-1
            IF (top_height_direction .eq. 2) THEN
              ilev = nlev - k1 
            ELSE
              ilev = k1
            ENDIF
            !cdir nodep
            DO j=1,npoints 
             IF (ilev .ge. itrop(j)) THEN
              IF ((at(j,ilev)   .ge. tb(j,ibox) .and. &
     &          at(j,ilev+1) .le. tb(j,ibox)) .or. &
     &          (at(j,ilev) .le. tb(j,ibox) .and.  &
     &          at(j,ilev+1) .ge. tb(j,ibox))) THEN 
                nmatch(j)=nmatch(j)+1
                match(j,nmatch(j))=ilev
              ENDIF  
             ENDIF                         
            ENDDO
29        CONTINUE

          DO j=1,npoints 
            IF (nmatch(j) .ge. 1) THEN
              k1 = match(j,nmatch(j))
              k2 = k1 + 1
              logp1 = log(pfull(j,k1))
              logp2 = log(pfull(j,k2))
              atd = max(tauchk,abs(at(j,k2) - at(j,k1)))
              logp=logp1+(logp2-logp1)*abs(tb(j,ibox)-at(j,k1))/atd
              ptop(j,ibox) = exp(logp)
              IF(abs(pfull(j,k1)-ptop(j,ibox)) .lt. &
     &            abs(pfull(j,k2)-ptop(j,ibox))) THEN
                 levmatch(j,ibox)=k1
              ELSE
                 levmatch(j,ibox)=k2
              ENDIF   
            ELSE
              IF (tb(j,ibox) .le. attrop(j)) THEN
                ptop(j,ibox)=ptrop(j)
                levmatch(j,ibox)=itrop(j)
              ENDIF
              IF (tb(j,ibox) .ge. atmax(j)) THEN
                ptop(j,ibox)=pfull(j,nlev)
                levmatch(j,ibox)=nlev
              ENDIF                                
            ENDIF
          ENDDO ! j

        ELSE ! if (top_height .eq. 1 .or. top_height .eq. 3) 
 
          DO j=1,npoints     
            ptop(j,ibox)=0._wp
          ENDDO
          DO ilev=1,nlev
            DO j=1,npoints     
              IF ((ptop(j,ibox) .eq. 0._wp ) &
     &           .and.(frac_out(j,ibox,ilev) .ne. 0._wp)) THEN
                ptop(j,ibox)=phalf(j,ilev)
              levmatch(j,ibox)=ilev
              ENDIF
            ENDDO
          ENDDO
        ENDIF                            
          
        DO j=1,npoints
          IF (tau(j,ibox) .le. (tauchk            )) THEN
            ptop(j,ibox)=0._wp
            levmatch(j,ibox)=0      
          ENDIF 
        ENDDO

30    CONTINUE
              
!     
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined, 
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy 
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.  
!
      !compute isccp frequencies

      !reset frequencies
      DO 37 ilev=1,7
      DO 38 ilev2=1,7
        DO j=1,npoints ! 
             IF (sunlit(j).eq.1 .or. top_height .eq. 3) THEN 
                fq_isccp(j,ilev,ilev2)= 0._wp
             ELSE
                fq_isccp(j,ilev,ilev2)= output_missing_value
             ENDIF
        ENDDO
38    CONTINUE
37    CONTINUE

      !reset variables need for averaging cloud properties
      DO j=1,npoints 
        IF (sunlit(j).eq.1 .or. top_height .eq. 3) THEN 
             totalcldarea(j) = 0._wp
             meanalbedocld(j) = 0._wp
             meanptop(j) = 0._wp
             meantaucld(j) = 0._wp
        ELSE
             totalcldarea(j) = output_missing_value
             meanalbedocld(j) = output_missing_value
             meanptop(j) = output_missing_value
             meantaucld(j) = output_missing_value
        ENDIF
      ENDDO ! j

      boxarea = 1._wp/REAL(ncol,wp)
     
      DO 39 ibox=1,ncol
        DO j=1,npoints 

          IF (tau(j,ibox) .gt. (tauchk            ) &
     &      .and. ptop(j,ibox) .gt. 0._wp) THEN
              box_cloudy(j,ibox)=.true.
          ENDIF

          IF (box_cloudy(j,ibox)) THEN

              IF (sunlit(j).eq.1 .or. top_height .eq. 3) THEN

                boxtau(j,ibox) = tau(j,ibox)

                IF (tau(j,ibox) .ge. isccp_taumin) THEN

                   totalcldarea(j) = totalcldarea(j) + boxarea
                
                   !convert optical thickness to albedo
                   albedocld(j,ibox) &
     &                   = (tau(j,ibox)**0.895_wp)/((tau(j,ibox)**0.895_wp)+6.82_wp)
         
                   !contribute to averaging
                   meanalbedocld(j) = meanalbedocld(j) &
     &                                +albedocld(j,ibox)*boxarea

                ENDIF
            ENDIF
          ENDIF

          IF (sunlit(j).eq.1 .or. top_height .eq. 3) THEN 

           IF (box_cloudy(j,ibox)) THEN
          
              !convert ptop to millibars
              ptop(j,ibox)=ptop(j,ibox) / 100._wp
            
              !save for output cloud top pressure and optical thickness
              boxptop(j,ibox) = ptop(j,ibox)
    
              IF (tau(j,ibox) .ge. isccp_taumin) THEN
               meanptop(j) = meanptop(j) + ptop(j,ibox)*boxarea
              ENDIF

              !reset itau(j), ipres(j)
              itau(j) = 0
              ipres(j) = 0

              !determine optical depth category
              IF (tau(j,ibox) .lt. isccp_taumin) THEN
                  itau(j)=1
              ELSE IF (tau(j,ibox) .ge. isccp_taumin &
     &          .and. tau(j,ibox) .lt. 1.3_wp) THEN
                itau(j)=2
              ELSE IF (tau(j,ibox) .ge. 1.3_wp  &
     &          .and. tau(j,ibox) .lt. 3.6_wp) THEN
                itau(j)=3
              ELSE IF (tau(j,ibox) .ge. 3.6_wp &
     &          .and. tau(j,ibox) .lt. 9.4_wp) THEN
                  itau(j)=4
              ELSE IF (tau(j,ibox) .ge. 9.4_wp &
     &          .and. tau(j,ibox) .lt. 23._wp) THEN
                  itau(j)=5
              ELSE IF (tau(j,ibox) .ge. 23._wp &
     &          .and. tau(j,ibox) .lt. 60._wp) THEN
                  itau(j)=6
              ELSE IF (tau(j,ibox) .ge. 60._wp) THEN
                  itau(j)=7
              ENDIF

              !determine cloud top pressure category
              IF (    ptop(j,ibox) .gt. 0._wp  &
     &          .and.ptop(j,ibox) .lt. 180._wp) THEN
                  ipres(j)=1
              ELSE IF(ptop(j,ibox) .ge. 180._wp &
     &          .and.ptop(j,ibox) .lt. 310._wp) THEN
                  ipres(j)=2
              ELSE IF(ptop(j,ibox) .ge. 310._wp &
     &          .and.ptop(j,ibox) .lt. 440._wp) THEN
                  ipres(j)=3
              ELSE IF(ptop(j,ibox) .ge. 440._wp &
     &          .and.ptop(j,ibox) .lt. 560._wp) THEN
                  ipres(j)=4
              ELSE IF(ptop(j,ibox) .ge. 560._wp &
     &          .and.ptop(j,ibox) .lt. 680._wp) THEN
                  ipres(j)=5
              ELSE IF(ptop(j,ibox) .ge. 680._wp &
     &          .and.ptop(j,ibox) .lt. 800._wp) THEN
                  ipres(j)=6
              ELSE IF(ptop(j,ibox) .ge. 800._wp) THEN
                  ipres(j)=7
              ENDIF 

              !update frequencies
              IF(ipres(j) .gt. 0.and.itau(j) .gt. 0) THEN
              fq_isccp(j,itau(j),ipres(j))= &
     &          fq_isccp(j,itau(j),ipres(j))+ boxarea
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! j
39    CONTINUE
      
      !compute mean cloud properties
      DO j=1,npoints 
        IF (totalcldarea(j) .gt. 0._wp) THEN
          ! code above guarantees that totalcldarea > 0 
          ! only if sunlit .eq. 1 .or. top_height = 3 
          ! and applies only to clouds with tau > isccp_taumin
          meanptop(j) = meanptop(j) / totalcldarea(j)
          meanalbedocld(j) = meanalbedocld(j) / totalcldarea(j)
          meantaucld(j) = (6.82_wp/((1._wp/meanalbedocld(j))-1._wp))**(1._wp/0.895_wp)
        ELSE
          ! this code is necessary so that in the case that totalcldarea = 0.,
          ! that these variables, which are in-cloud averages, are set to missing
          ! note that totalcldarea will be 0. if all the clouds in the grid box have
          ! tau < isccp_taumin 
          meanptop(j) = output_missing_value
          meanalbedocld(j) = output_missing_value
          meantaucld(j) = output_missing_value
        ENDIF
      ENDDO ! j

!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
      IF (debugcol.ne.0) THEN
!     
         DO j=1,npoints,debugcol

            !produce character output
            DO ilev=1,nlev
              DO ibox=1,ncol
                   acc(ilev,ibox)=0
              ENDDO
            ENDDO

            DO ilev=1,nlev
              DO ibox=1,ncol
                   acc(ilev,ibox)=INT(frac_out(j,ibox,ilev)*2._wp)
                   IF (levmatch(j,ibox) .eq. ilev)  &
     &                 acc(ilev,ibox)=acc(ilev,ibox)+1
              ENDDO
            ENDDO

          write(ftn09,11) j
11        format('ftn09.',i4.4)
          open(9, FILE=ftn09, FORM='FORMATTED')
             write(9,'(a1)') ' '
             write(9,'(10i5)') &
     &                  (ilev,ilev=5,nlev,5)
             write(9,'(a1)') ' '
             DO ibox=1,ncol
               write(9,'(40(a1),1x,40(a1))') &
     &           (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev) &
     &           ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev) 
             ENDDO
             close(9)

             IF (ncolprint.ne.0) THEN
               write(6,'(a1)') ' '
                    write(6,'(a2,1X,5(a7,1X),a50)') &
     &                  'ilev',&
     &                  'pfull','at',&
     &                  'cc*100','dem_s','dtau_s',&
     &                  'cchar'

!               DO 4012 ilev=1,nlev
!                    write(6,'(60i2)') (box(i,ilev),i=1,ncolprint)
!                   write(6,'(i2,1X,5(f7.2,1X),50(a1))') 
!     &                  ilev,
!     &                  pfull(j,ilev)/100.,at(j,ilev),
!     &                  cc(j,ilev)*100.0,dem_s(j,ilev),dtau_s(j,ilev)
!     &                  ,(cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
!4012           continue
               write (6,'(a)') 'skt(j):'
               write (6,'(8f7.2)') skt(j)
               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
               write (6,'(a)') 'tau:'
               write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
               write (6,'(a)') 'tb:'
               write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
               write (6,'(a)') 'ptop:'
               write (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
             ENDIF 
        ENDDO
      ENDIF 

      END SUBROUTINE icarus 
