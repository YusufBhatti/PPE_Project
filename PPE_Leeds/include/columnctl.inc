!Namelist for Single Column Model 
!
!Author: Suvarchal K. Cheedela  Feb 2010
!
!used in mo_column
!For detailed description of switches and flags 
!consult ECHAM Single Column model user guide.
!
NAMELIST /columnctl/ & 
!
!name      ! & description         & default
!
! (/a,b,c/) 
!         a= tendency option ( 0 no prescribed LS tendency , 1 prescribed LS tendency)
!         b= relaxation time in seconds
!         c= extrapolation option ( 0 Fixed period of forcing ,1 cycle forcings in time)     
nfor_t,    &!forcing temparature
nfor_ps,   &!forcing surface pressure
nfor_uv,   &!forcing horizontal winds
nfor_q,    &!forcing water vapor
nfor_xi,   &!forcing ice water 
nfor_xl,   &!forcing liquid water
!
!(/a,b/)
!a= read option ( 0 no external state prescribed , 1 state prescribed)
!b= interpolation option ( 0 Fixed period of forcing ,1 cycle forcings in time)

nfor_uvgeo, &!forcing geostrophic wind 
nfor_omega, &!forcing with pressure velocity
nfor_div,   &!forcing with divergence
nfor_ts,    &!forcing surface temparature
nfor_lhf,   &!forcing latent heat flux
nfor_shf,   &!forcing sensible heat flux
forcingfile,&!forcing file name
mld,        &!mixed layer depth
ml_input     !if true, the initial temperature of the mixed layer ocean is set to the value  
             !of the first entry of the surface temperature ts of the forcingfile 

