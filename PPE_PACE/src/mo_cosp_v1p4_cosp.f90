! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!#include "cosp_defs.h"
MODULE mo_cosp_v1p4_cosp
  USE mo_cosp_types
  USE mo_cosp_v1p4_cosp_simulator
  USE mo_cosp_llnl_precip
  USE mo_cosp_modis_simulator
  USE mo_cosp_utils
  USE mo_kind,           ONLY: wp

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP ---------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!CNam: #ifdef RTTOV
!SUBROUTINE COSP(overlap,kbdim, klev, Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
!#else
SUBROUTINE COSP(overlap,kbdim, klev, Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
!#endif
  ! Arguments
  INTEGER, INTENT(in) :: overlap !  overlap type in SCOPS: 1=max, 2=rand, 3=max/rand

  INTEGER, INTENT(IN) :: kbdim         ! Number of gridpoints
  INTEGER, INTENT(IN) :: klev          ! Number of levels
  INTEGER, INTENT(IN) :: Ncolumns      ! Number of columns

  type(cosp_config),INTENT(in) :: cfg   ! Configuration options
  type(cosp_vgrid),INTENT(in) :: vgrid   ! Information on vertical grid of stats
  type(cosp_gridbox),INTENT(inout) :: gbx
  type(cosp_subgrid),INTENT(inout) :: sgx   ! Subgrid info
  type(cosp_sgradar),INTENT(inout) :: sgradar ! Output from radar simulator
  type(cosp_sglidar),INTENT(inout) :: sglidar ! Output from lidar simulator
  type(cosp_isccpstats),INTENT(inout)   :: isccp   ! Output from ISCCP simulator
  type(cosp_misrstats),INTENT(inout)    :: misr    ! Output from MISR simulator
  type(cosp_modisstats),INTENT(inout)   :: modis   ! Output from MODIS simulator
!CNam: #ifdef RTTOV
!  type(cosp_rttov),INTENT(inout)   :: rttov   ! Output from RTTOV
!#endif
  type(cosp_radarstats),INTENT(inout) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats),INTENT(inout) :: stlidar ! Summary statistics from lidar simulator

  ! Local variables 
  INTEGER :: Npoints   ! Number of gridpoints
  INTEGER :: Nlevels   ! Number of levels
  INTEGER :: Nhydro    ! Number of hydrometeors
  INTEGER :: Niter     ! Number of calls to cosp_simulator
  INTEGER :: i_first,i_last ! First and last gridbox to be processed in each iteration
  INTEGER :: i,Ni
  INTEGER,DIMENSION(2) :: ix,iy
  logical :: reff_zero
  REAL(wp) :: maxp,minp
  INTEGER,DIMENSION(:),allocatable :: & ! DIMENSIONs nPoints
                  seed    !  It is recommended that the seed is set to a different value for each model
                          !  gridbox it is called on, as it is possible that the choice of the same 
                          !  seed value every time may introduce some statistical bias in the results, 
                          !  particularly for low values of NCOL.

!CNam: removed 'one iteration' 

!++++++++++ DIMENSIONs ++++++++++++
  Npoints  = gbx%Npoints
  Nlevels  = gbx%Nlevels
  Nhydro   = gbx%Nhydro

!++++++++++ Depth of model layers ++++++++++++
  do i=1,Nlevels-1
    gbx%dlev(:,i) = gbx%zlev_half(:,i+1) - gbx%zlev_half(:,i)
  enddo
  gbx%dlev(:,Nlevels) = 2.0*(gbx%zlev(:,Nlevels) - gbx%zlev_half(:,Nlevels))

!++++++++++ Apply sanity checks to inputs ++++++++++
  call cosp_check_input('longitude',gbx%longitude,min_val=0.0_wp,max_val=360.0_wp)
  call cosp_check_input('latitude',gbx%latitude,min_val=-90.0_wp,max_val=90.0_wp)
  call cosp_check_input('dlev',gbx%dlev,min_val=0.0_wp)
  call cosp_check_input('p',gbx%p,min_val=0.0_wp)
  call cosp_check_input('ph',gbx%ph,min_val=0.0_wp)
  call cosp_check_input('T',gbx%T,min_val=0.0_wp)
  call cosp_check_input('q',gbx%q,min_val=0.0_wp)
  call cosp_check_input('sh',gbx%sh,min_val=0.0_wp)
  call cosp_check_input('dtau_s',gbx%dtau_s,min_val=0.0_wp)
  call cosp_check_input('dtau_c',gbx%dtau_c,min_val=0.0_wp)
  call cosp_check_input('dem_s',gbx%dem_s,min_val=0.0_wp,max_val=1.0_wp)
  call cosp_check_input('dem_c',gbx%dem_c,min_val=0.0_wp,max_val=1.0_wp)
  ! Point information (Npoints)
  call cosp_check_input('land',gbx%land,min_val=0.0_wp,max_val=1.0_wp)
  call cosp_check_input('psfc',gbx%psfc,min_val=0.0_wp)
  call cosp_check_input('sunlit',gbx%sunlit,min_val=0.0_wp,max_val=1.0_wp)
  call cosp_check_input('skt',gbx%skt,min_val=0.0_wp)
  ! TOTAL and CONV cloud fraction for SCOPS
  call cosp_check_input('tca',gbx%tca,min_val=0.0_wp,max_val=1.0_wp)
  call cosp_check_input('cca',gbx%cca,min_val=0.0_wp,max_val=1.0_wp)
  ! Precipitation fluxes on model levels
  call cosp_check_input('rain_ls',gbx%rain_ls,min_val=0.0_wp)
!CNam:  call cosp_check_input('rain_cv',gbx%rain_cv,min_val=0.0_wp)   ECHAM convective rain can = 0.
  call cosp_check_input('snow_ls',gbx%snow_ls,min_val=0.0_wp)
!CNam:  call cosp_check_input('snow_cv',gbx%snow_cv,min_val=0.0_wp)   ECHAM convective snow can = 0.
  call cosp_check_input('grpl_ls',gbx%grpl_ls,min_val=0.0_wp)
  ! Hydrometeors concentration and distribution parameters
!CNam:  call cosp_check_input('mr_hydro',gbx%mr_hydro,min_val=0.0_wp) ECHAM mr_hydro can = 0.
  ! Effective radius [m]. (Npoints,Nlevels,Nhydro)
  call cosp_check_input('Reff',gbx%Reff,min_val=0.0_wp)
  reff_zero=.true.
  if (any(gbx%Reff > 1.e-8_wp)) then
     reff_zero=.false.
      ! reff_zero == .false.
      !     and gbx%use_reff == .true.   Reff use in radar and lidar
      !     and reff_zero    == .false.  Reff use in lidar and set to 0 for radar
  endif
  if ((.not. gbx%use_reff) .and. (reff_zero)) then ! No Reff in radar. Default in lidar
        gbx%Reff = DEFAULT_LIDAR_REFF
        print *, '---------- COSP WARNING ------------'
        print *, ''
        print *, 'Using default Reff in lidar simulations'
        print *, ''
        print *, '----------------------------------'
  endif
  
  ! Aerosols concentration and distribution parameters
  call cosp_check_input('conc_aero',gbx%conc_aero,min_val=0.0_wp)
  ! Checks for CRM mode !CNam: for ICON_LES
  if (Ncolumns == 1) then
     if (gbx%use_precipitation_fluxes) then
        print *, '---------- COSP ERROR ------------'
        print *, ''
        print *, 'Use of precipitation fluxes not supported in CRM mode (Ncolumns=1)'
        print *, ''
        print *, '----------------------------------'
        stop
     endif
     if ((maxval(gbx%dtau_c) > 0.0_wp).or.(maxval(gbx%dem_c) > 0.0_wp)) then
        print *, '---------- COSP ERROR ------------'
        print *, ''
        print *, ' dtau_c > 0.0 or dem_c > 0.0. In CRM mode (Ncolumns=1), '
        print *, ' the optical depth (emmisivity) of all clouds must be '
        print *, ' passed through dtau_s (dem_s)'
        print *, ''
        print *, '----------------------------------'
        stop
     endif
  endif

   ! We base the seed in the decimal part of the surface pressure.
   allocate(seed(kbdim))
   seed = int(gbx%psfc) ! This is to avoid division by zero when Npoints = 1   
      ! Roj Oct/2008 ... Note: seed value of 0 caused me some problems + I want to 
      ! randomize for each call to COSP even when Npoints ==1
   minp = minval(gbx%psfc)
   maxp = maxval(gbx%psfc)
   if (kbdim .gt. 1) seed=int((gbx%psfc-minp)/(maxp-minp)*100000) + 1
   ! Below it's how it was done in the original implementation of the ISCCP simulator. 
   ! The one above is better for offline data, when you may have packed data 
   ! that subsamples the decimal fraction of the surface pressure. 
!    if (Npoints .gt. 1) seed=(gbx%psfc-int(gbx%psfc))*1000000 

!CNam: Removed 'one iteration gbx%Npoints'  
!CNam: #ifdef RTTOV
!        call cosp_iter(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
!#else
        call cosp_iter(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
!#endif

   deallocate(seed)
    
END SUBROUTINE COSP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_ITER ----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!CNam: #ifdef RTTOV
! SUBROUTINE COSP_ITER(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
! #else
SUBROUTINE COSP_ITER(overlap,seed,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
! #endif

  ! Arguments
  INTEGER,INTENT(in) :: overlap !  overlap type in SCOPS: 1=max, 2=rand, 3=max/rand
  INTEGER,DIMENSION(:),INTENT(in) :: seed
  type(cosp_config),INTENT(in) :: cfg   ! Configuration options
  type(cosp_vgrid),INTENT(in) :: vgrid   ! Information on vertical grid of stats 
  type(cosp_gridbox),INTENT(inout) :: gbx
  type(cosp_subgrid),INTENT(inout) :: sgx   ! Subgrid info
  type(cosp_sgradar),INTENT(inout) :: sgradar ! Output from radar simulator
  type(cosp_sglidar),INTENT(inout) :: sglidar ! Output from lidar simulator
  type(cosp_isccpstats),INTENT(inout)   :: isccp   ! Output from ISCCP simulator
  type(cosp_misrstats),INTENT(inout)    :: misr    ! Output from MISR simulator
  type(cosp_modisstats),INTENT(inout)   :: modis   ! Output from MODIS simulator
!CNam: #ifdef RTTOV
!  type(cosp_rttov),INTENT(inout)   :: rttov   ! Output from RTTOV
!#endif
  type(cosp_radarstats),INTENT(inout) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats),INTENT(inout) :: stlidar ! Summary statistics from lidar simulator

  ! Local variables 
  INTEGER :: Npoints   ! Number of gridpoints
  INTEGER :: Ncolumns  ! Number of subcolumns
  INTEGER :: Nlevels   ! Number of levels
  INTEGER :: Nhydro    ! Number of hydrometeors
  INTEGER :: i,j,k
  INTEGER :: I_HYDRO 
  REAL(wp),DIMENSION(:,:),POINTER :: column_frac_out ! Array with one column of frac_out
  REAL(wp),DIMENSION(:,:),POINTER :: column_prec_out ! Array with one column of prec_frac
  INTEGER :: scops_debug=0    !  set to non-zero value to print out inputs for debugging in SCOPS
  REAL(wp),DIMENSION(:, :),allocatable :: cca_scops,ls_p_rate,cv_p_rate, &
                     tca_scops ! Cloud cover in each model level (HORIZONTAL gridbox fraction) of total cloud.
                               ! Levels are from TOA to SURFACE. (nPoints, nLev)
  REAL(wp),DIMENSION(:,:),allocatable :: frac_ls,prec_ls,frac_cv,prec_cv ! Cloud/Precipitation fraction in each model level
                                                                     ! Levels are from SURFACE to TOA
  REAL(wp),DIMENSION(:,:),allocatable :: rho ! (Npoints, Nlevels). Atmospheric density
  type(cosp_sghydro) :: sghydro   ! Subgrid info for hydrometeors en each iteration


  !++++++++++ DIMENSIONs ++++++++++++
  Npoints  = gbx%Npoints
  Ncolumns = gbx%Ncolumns
  Nlevels  = gbx%Nlevels
  Nhydro   = gbx%Nhydro
    
  !++++++++++ Climate/NWP mode ++++++++++  
  if (Ncolumns > 1) then
        !++++++++++ Subgrid sampling ++++++++++
        ! Allocate arrays before calling SCOPS
        allocate(frac_ls(Npoints,Nlevels),frac_cv(Npoints,Nlevels),prec_ls(Npoints,Nlevels),prec_cv(Npoints,Nlevels))
        allocate(tca_scops(Npoints,Nlevels),cca_scops(Npoints,Nlevels), &
                ls_p_rate(Npoints,Nlevels),cv_p_rate(Npoints,Nlevels))
        ! Initialize to zero
        frac_ls=0.0_wp
        prec_ls=0.0_wp
        frac_cv=0.0_wp
        prec_cv=0.0_wp
        ! Cloud fractions for SCOPS from TOA to SFC
        tca_scops = gbx%tca(:,Nlevels:1:-1)
        cca_scops = gbx%cca(:,Nlevels:1:-1)
        
        ! Call to SCOPS
        ! strat and conv arrays are passed with levels from TOA to SURFACE.
        call cosp_scops(Npoints,Nlevels,Ncolumns,tca_scops,overlap,sgx%frac_out,scops_debug) !CNam: no seed, cca_scops
        ! temporarily use prec_ls/cv to transfer information about precipitation flux into prec_scops
        if(gbx%use_precipitation_fluxes) then
            ls_p_rate(:,Nlevels:1:-1)=gbx%rain_ls(:,1:Nlevels)+gbx%snow_ls(:,1:Nlevels) !CNam:  +gbx%grpl_ls(:,1:Nlevels)
            cv_p_rate(:,Nlevels:1:-1)=gbx%rain_cv(:,1:Nlevels)+gbx%snow_cv(:,1:Nlevels)
        else
            ls_p_rate(:,Nlevels:1:-1)=gbx%mr_hydro(:,1:Nlevels,I_LSRAIN)+ &
                                      gbx%mr_hydro(:,1:Nlevels,I_LSSNOW)	!CNam: +gbx%mr_hydro(:,1:Nlevels,I_LSGRPL)
            cv_p_rate(:,Nlevels:1:-1)=gbx%mr_hydro(:,1:Nlevels,I_CVRAIN)+ &
                                      gbx%mr_hydro(:,1:Nlevels,I_CVSNOW)
        endif
        
        call prec_scops(Npoints,Nlevels,Ncolumns,ls_p_rate,cv_p_rate,sgx%frac_out,sgx%prec_frac)

        ! Precipitation fraction
        do j=1,Npoints,1	
        do k=1,Nlevels,1	
            do i=1,Ncolumns,1
                if (sgx%frac_out (j,i,Nlevels+1-k) == I_LSC) frac_ls(j,k)=frac_ls(j,k)+1._wp
                if (sgx%frac_out (j,i,Nlevels+1-k) == I_CVC) frac_cv(j,k)=frac_cv(j,k)+1._wp
                if (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 1) prec_ls(j,k)=prec_ls(j,k)+1._wp
                if (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 2) prec_cv(j,k)=prec_cv(j,k)+1._wp
                if (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 3) then
                    prec_cv(j,k)=prec_cv(j,k)+1._wp
                    prec_ls(j,k)=prec_ls(j,k)+1._wp
                endif
            enddo  !i
            frac_ls(j,k)=frac_ls(j,k)/Ncolumns
            frac_cv(j,k)=frac_cv(j,k)/Ncolumns
            prec_ls(j,k)=prec_ls(j,k)/Ncolumns
            prec_cv(j,k)=prec_cv(j,k)/Ncolumns
        enddo  !k
        enddo  !j
        
         ! Levels from SURFACE to TOA.
        if (Npoints*Ncolumns*Nlevels < 10000) then
            sgx%frac_out(1:Npoints,:,1:Nlevels)  = sgx%frac_out(1:Npoints,:,Nlevels:1:-1)
            sgx%prec_frac(1:Npoints,:,1:Nlevels) = sgx%prec_frac(1:Npoints,:,Nlevels:1:-1)
        else
            ! This is done within a loop (unvectorized) over nPoints to save memory
            do j=1,Npoints
                sgx%frac_out(j,:,1:Nlevels)  = sgx%frac_out(j,:,Nlevels:1:-1)
                sgx%prec_frac(j,:,1:Nlevels) = sgx%prec_frac(j,:,Nlevels:1:-1)
            enddo
        endif

       ! Deallocate arrays that will no longer be used
        deallocate(tca_scops,cca_scops,ls_p_rate,cv_p_rate)

        ! Populate the subgrid arrays
        call construct_cosp_sghydro(Npoints,Ncolumns,Nlevels,Nhydro,sghydro)

        do k=1,Ncolumns
            !--------- Mixing ratios for clouds and Reff for Clouds and precip -------
            column_frac_out => sgx%frac_out(:,k,:)
            where (column_frac_out == I_LSC)     !+++++++++++ LS clouds ++++++++
                sghydro%mr_hydro(:,k,:,I_LSCLIQ) = gbx%mr_hydro(:,:,I_LSCLIQ)
                sghydro%mr_hydro(:,k,:,I_LSCICE) = gbx%mr_hydro(:,:,I_LSCICE)

                sghydro%Reff(:,k,:,I_LSCLIQ)     = gbx%Reff(:,:,I_LSCLIQ)
                sghydro%Reff(:,k,:,I_LSCICE)     = gbx%Reff(:,:,I_LSCICE)

                sghydro%Np(:,k,:,I_LSCLIQ)     = gbx%Np(:,:,I_LSCLIQ)
                sghydro%Np(:,k,:,I_LSCICE)     = gbx%Np(:,:,I_LSCICE)

            elsewhere (column_frac_out == I_CVC) !+++++++++++ CONV clouds ++++++++
                sghydro%mr_hydro(:,k,:,I_CVCLIQ) = gbx%mr_hydro(:,:,I_CVCLIQ)
                sghydro%mr_hydro(:,k,:,I_CVCICE) = gbx%mr_hydro(:,:,I_CVCICE)

                sghydro%Reff(:,k,:,I_CVCLIQ)     = gbx%Reff(:,:,I_CVCLIQ)
                sghydro%Reff(:,k,:,I_CVCICE)     = gbx%Reff(:,:,I_CVCICE)

                sghydro%Np(:,k,:,I_CVCLIQ)     = gbx%Np(:,:,I_CVCLIQ)
                sghydro%Np(:,k,:,I_CVCICE)     = gbx%Np(:,:,I_CVCICE)

            end where
            column_prec_out => sgx%prec_frac(:,k,:)
            where ((column_prec_out == 1) .or. (column_prec_out == 3) )  !++++ LS precip ++++
                sghydro%Reff(:,k,:,I_LSRAIN) = gbx%Reff(:,:,I_LSRAIN)
                sghydro%Reff(:,k,:,I_LSSNOW) = gbx%Reff(:,:,I_LSSNOW)
!                sghydro%Reff(:,k,:,I_LSGRPL) = gbx%Reff(:,:,I_LSGRPL)

                sghydro%Np(:,k,:,I_LSRAIN)     = gbx%Np(:,:,I_LSRAIN)
                sghydro%Np(:,k,:,I_LSSNOW)     = gbx%Np(:,:,I_LSSNOW)
!                sghydro%Np(:,k,:,I_LSGRPL)     = gbx%Np(:,:,I_LSGRPL)
            elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3)) !++++ CONV precip ++++
                sghydro%Reff(:,k,:,I_CVRAIN) = gbx%Reff(:,:,I_CVRAIN)
                sghydro%Reff(:,k,:,I_CVSNOW) = gbx%Reff(:,:,I_CVSNOW)

                sghydro%Np(:,k,:,I_CVRAIN)     = gbx%Np(:,:,I_CVRAIN)
                sghydro%Np(:,k,:,I_CVSNOW)     = gbx%Np(:,:,I_CVSNOW)
            end where
            !--------- Precip -------
            if (.not. gbx%use_precipitation_fluxes) then
                where (column_frac_out == I_LSC)  !+++++++++++ LS Precipitation ++++++++
                    sghydro%mr_hydro(:,k,:,I_LSRAIN) = gbx%mr_hydro(:,:,I_LSRAIN)
                    sghydro%mr_hydro(:,k,:,I_LSSNOW) = gbx%mr_hydro(:,:,I_LSSNOW)
!                    sghydro%mr_hydro(:,k,:,I_LSGRPL) = gbx%mr_hydro(:,:,I_LSGRPL)
                elsewhere (column_frac_out == I_CVC) !+++++++++++ CONV Precipitation ++++++++
                    sghydro%mr_hydro(:,k,:,I_CVRAIN) = gbx%mr_hydro(:,:,I_CVRAIN)
                    sghydro%mr_hydro(:,k,:,I_CVSNOW) = gbx%mr_hydro(:,:,I_CVSNOW)
                end where 
            endif
        enddo
        ! convert the mixing ratio and precipitation flux from gridbox mean to the fraction-based values
        do k=1,Nlevels
            do j=1,Npoints
                !--------- Clouds -------
                if (frac_ls(j,k) .ne. 0._wp) then
                    sghydro%mr_hydro(j,:,k,I_LSCLIQ) = sghydro%mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                    sghydro%mr_hydro(j,:,k,I_LSCICE) = sghydro%mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
                endif
                if (frac_cv(j,k) .ne. 0._wp) then
                    sghydro%mr_hydro(j,:,k,I_CVCLIQ) = sghydro%mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                    sghydro%mr_hydro(j,:,k,I_CVCICE) = sghydro%mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
                endif
                !--------- Precip -------
                if (gbx%use_precipitation_fluxes) then
                    if (prec_ls(j,k) .ne. 0._wp) then
                        gbx%rain_ls(j,k) = gbx%rain_ls(j,k)/prec_ls(j,k)
                        gbx%snow_ls(j,k) = gbx%snow_ls(j,k)/prec_ls(j,k)
!CNam:                  gbx%grpl_ls(j,k) = gbx%grpl_ls(j,k)/prec_ls(j,k)
                    endif
                    if (prec_cv(j,k) .ne. 0._wp) then
                        gbx%rain_cv(j,k) = gbx%rain_cv(j,k)/prec_cv(j,k)
                        gbx%snow_cv(j,k) = gbx%snow_cv(j,k)/prec_cv(j,k)
                    endif
                else
                    if (prec_ls(j,k) .ne. 0._wp) then
                        sghydro%mr_hydro(j,:,k,I_LSRAIN) = sghydro%mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                        sghydro%mr_hydro(j,:,k,I_LSSNOW) = sghydro%mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
!CNam:                  sghydro%mr_hydro(j,:,k,I_LSGRPL) = sghydro%mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
                    endif
                    if (prec_cv(j,k) .ne. 0._wp) then
                        sghydro%mr_hydro(j,:,k,I_CVRAIN) = sghydro%mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                        sghydro%mr_hydro(j,:,k,I_CVSNOW) = sghydro%mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
                    endif
                endif  
            enddo !k
        enddo !j
        deallocate(frac_ls,prec_ls,frac_cv,prec_cv)
        
        if (gbx%use_precipitation_fluxes) then
        
#ifdef MMF_V3p5_TWO_MOMENT

        write(*,*) 'Precipitation Flux to Mixing Ratio conversion not (yet?) supported ', &
               'for MMF3.5 Two Moment Microphysics'
        stop
#else
            ! Density
            allocate(rho(Npoints,Nlevels))
            I_HYDRO = I_LSRAIN
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,1._wp, &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO), &
                    gamma_1(I_HYDRO),gamma_2(I_HYDRO),gamma_3(I_HYDRO),gamma_4(I_HYDRO), &
                    gbx%rain_ls,sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
            I_HYDRO = I_LSSNOW
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,1._wp, &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO), &
                    gamma_1(I_HYDRO),gamma_2(I_HYDRO),gamma_3(I_HYDRO),gamma_4(I_HYDRO), &
                    gbx%snow_ls,sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
            I_HYDRO = I_CVRAIN
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,2._wp, &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO), &
                    gamma_1(I_HYDRO),gamma_2(I_HYDRO),gamma_3(I_HYDRO),gamma_4(I_HYDRO), &
                    gbx%rain_cv,sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
            I_HYDRO = I_CVSNOW
            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,2._wp, &
                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO), &
                    gamma_1(I_HYDRO),gamma_2(I_HYDRO),gamma_3(I_HYDRO),gamma_4(I_HYDRO), &
                    gbx%snow_cv,sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
!CNam:       I_HYDRO = I_LSGRPL
!            call cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,sgx%prec_frac,1., &
!                    n_ax(I_HYDRO),n_bx(I_HYDRO),alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO), &
!                    g_x(I_HYDRO),a_x(I_HYDRO),b_x(I_HYDRO), &
!                    gamma_1(I_HYDRO),gamma_2(I_HYDRO),gamma_3(I_HYDRO),gamma_4(I_HYDRO), &
!                    gbx%grpl_ls,sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
            if(allocated(rho)) deallocate(rho)
#endif
        endif
   !++++++++++ CRM mode ++++++++++ !CNam: Note for ICON_LES
   else
      call construct_cosp_sghydro(Npoints,Ncolumns,Nlevels,Nhydro,sghydro)
      sghydro%mr_hydro(:,1,:,:) = gbx%mr_hydro
      sghydro%Reff(:,1,:,:) = gbx%Reff
      sghydro%Np(:,1,:,:) = gbx%Np      ! added by Roj with Quickbeam V3.0
      
      !--------- Clouds -------
      where ((gbx%dtau_s > 0.0_wp))
             sgx%frac_out(:,1,:) = 1  ! Subgrid cloud array. DIMENSIONs (Npoints,Ncolumns,Nlevels)
      endwhere
   endif ! Ncolumns > 1

   !++++++++++ Simulator ++++++++++
!CNam: #ifdef RTTOV
!    call cosp_simulator(kbdim,Ncolumns,klev,gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
!#else
    call cosp_simulator(kbdim,Ncolumns,klev,gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,misr,modis,stradar,stlidar) 
!#endif

    ! Deallocate subgrid arrays
    call free_cosp_sghydro(sghydro)

END SUBROUTINE COSP_ITER

END MODULE mo_cosp_v1p4_cosp
