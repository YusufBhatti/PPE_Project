MODULE MO_COSP_LLNL_PRECIP
CONTAINS

      subroutine pf_to_mr (npoints,nlev,ncol,rain_ls,snow_ls,grpl_ls,rain_cv,snow_cv,&
            prec_frac,app1,ptm1,mx_rain_ls,mx_snow_ls,mx_grpl_ls,mx_rain_cv,mx_snow_cv)

      USE MO_KIND,		ONLY: wp

      implicit none

      INTEGER npoints       !  number of model points in the horizontal
      INTEGER nlev          !  number of model levels in column
      INTEGER ncol          !  number of subcolumns

      INTEGER i,j,ilev,ibox
      
      REAL(wp) rain_ls(npoints,nlev)
      REAL(wp) snow_ls(npoints,nlev) ! large-scale precipitation flux
      REAL(wp) grpl_ls(npoints,nlev)  
      REAL(wp) rain_cv(npoints,nlev)
      REAL(wp) snow_cv(npoints,nlev) ! convective precipitation flux

      REAL(wp) prec_frac(npoints,ncol,nlev) ! 0 -> clear sky
                                        ! 1 -> LS precipitation
                                        ! 2 -> CONV precipitation
					! 3 -> both
      REAL(wp) mx_rain_ls(npoints,ncol,nlev)
      REAL(wp) mx_snow_ls(npoints,ncol,nlev)
      REAL(wp) mx_grpl_ls(npoints,ncol,nlev)
      REAL(wp) mx_rain_cv(npoints,ncol,nlev)
      REAL(wp) mx_snow_cv(npoints,ncol,nlev)
      REAL(wp) app1(npoints,nlev),ptm1(npoints,nlev)
      REAL(wp) ar,as,ag,br,bs,bg,nr,ns,ng,rho0,rhor,rhos,rhog,rho
      REAL(wp) term1r,term1s,term1g,term2r,term2s,term2g,term3
      REAL(wp) term4r_ls,term4s_ls,term4g_ls,term4r_cv,term4s_cv
      REAL(wp) term1x2r,term1x2s,term1x2g,t123r,t123s,t123g

      ! method from Khairoutdinov and Randall (2003 JAS)		
					
      ar=842._wp
      as=4.84_wp
      ag=94.5_wp
      br=0.8_wp
      bs=0.25_wp
      bg=0.5_wp
      nr=8._wp*1000._wp*1000._wp
      ns=3._wp*1000._wp*1000._wp
      ng=4._wp*1000._wp*1000._wp
      rho0=1.29_wp
      rhor=1000._wp
      rhos=100._wp
      rhog=400._wp
      term1r=ar*17.8379_wp/6._wp
      term1s=as*8.28508_wp/6._wp
      term1g=ag*11.6317_wp/6._wp
      term2r=(3.14159265_wp*rhor*nr)**(-br/4._wp)
      term2s=(3.14159265_wp*rhos*ns)**(-bs/4._wp)
      term2g=(3.14159265_wp*rhog*ng)**(-bg/4._wp)

      term1x2r=term1r*term2r
      term1x2s=term1s*term2s
      term1x2g=term1g*term2g
      do ilev=1,nlev
        do j=1,npoints
            rho=app1(j,ilev)/(287.05_wp*ptm1(j,ilev))
            term3=(rho0/rho)**0.5_wp
            ! Term 4 of Eq. (A19).
            t123r=term1x2r*term3
            t123s=term1x2s*term3
            t123g=term1x2g*term3
            term4r_ls=rain_ls(j,ilev)/(t123r)
            term4s_ls=snow_ls(j,ilev)/(t123s)
            term4g_ls=grpl_ls(j,ilev)/(t123g)
            term4r_cv=rain_cv(j,ilev)/(t123r)
            term4s_cv=snow_cv(j,ilev)/(t123s)
            do ibox=1,ncol
                mx_rain_ls(j,ibox,ilev)=0._wp
                mx_snow_ls(j,ibox,ilev)=0._wp
                mx_grpl_ls(j,ibox,ilev)=0._wp
                mx_rain_cv(j,ibox,ilev)=0._wp
                mx_snow_cv(j,ibox,ilev)=0._wp
                if ((prec_frac(j,ibox,ilev) .eq. 1._wp) .or. &
                    (prec_frac(j,ibox,ilev) .eq. 3._wp)) then 
                    mx_rain_ls(j,ibox,ilev)=&
                      (term4r_ls**(1._wp/(1._wp+br/4._wp)))/rho
                    mx_snow_ls(j,ibox,ilev)=&
                      (term4s_ls**(1._wp/(1._wp+bs/4._wp)))/rho
                    mx_grpl_ls(j,ibox,ilev)=&
                      (term4g_ls**(1._wp/(1._wp+bg/4._wp)))/rho
                endif
                if ((prec_frac(j,ibox,ilev) .eq. 2._wp) .or. &
                    (prec_frac(j,ibox,ilev) .eq. 3._wp)) then 
                    mx_rain_cv(j,ibox,ilev)=&
                      (term4r_cv**(1._wp/(1._wp+br/4._wp)))/rho
                    mx_snow_cv(j,ibox,ilev)=&
                      (term4s_cv**(1._wp/(1._wp+bs/4._wp)))/rho
                endif
            enddo ! loop over ncol
        enddo ! loop over npoints
      enddo ! loop over nlev     

      end subroutine pf_to_mr 
! ---------------------------------------------------------------------
      subroutine prec_scops(npoints,nlev,ncol,ls_p_rate,cv_p_rate,&
                           frac_out,prec_frac)

      USE MO_KIND,		ONLY: wp
      implicit none
	
      INTEGER npoints       !  number of model points in the horizontal
      INTEGER nlev          !  number of model levels in column
      INTEGER ncol          !  number of subcolumns

      INTEGER i,j,ilev,ibox,cv_col
      
      REAL (wp) ls_p_rate(npoints,nlev)
      REAL (wp) cv_p_rate(npoints,nlev)

      REAL (wp) frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column
                              !TOA to SURFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL (wp) prec_frac(npoints,ncol,nlev) ! 0 -> clear sky
                                        ! 1 -> LS precipitation
                                        ! 2 -> CONV precipitation
					! 3 -> both
                                        !TOA to SURFACE!!!!!!!!!!!!!!!!!!
					
      INTEGER flag_ls, flag_cv
      INTEGER frac_out_ls(npoints,ncol),frac_out_cv(npoints,ncol) !flag variables for 
                       ! stratiform cloud and convective cloud in the vertical column

     cv_col = 0.05_wp*ncol
     if (cv_col .eq. 0) cv_col=1

      do ilev=1,nlev
      do ibox=1,ncol
        do j=1,npoints 
        prec_frac(j,ibox,ilev) = 0._wp
        enddo
      enddo
      enddo
      
      do j=1,npoints
       do ibox=1,ncol
       frac_out_ls(j,ibox)=0
       frac_out_cv(j,ibox)=0
       flag_ls=0
       flag_cv=0
        do ilev=1,nlev
	 if (frac_out(j,ibox,ilev) .eq. 1) then 
	  flag_ls=1
	 endif
	 if (frac_out(j,ibox,ilev) .eq. 2) then 
	  flag_cv=1
	 endif
	enddo !loop over nlev
	if (flag_ls .eq. 1) then
	 frac_out_ls(j,ibox)=1
	endif
	if (flag_cv .eq. 1) then
	 frac_out_cv(j,ibox)=1
	endif
       enddo  ! loop over ncol
      enddo ! loop over npoints

!      initialize the top layer      
       do j=1,npoints
        flag_ls=0
	flag_cv=0
	
        if (ls_p_rate(j,1) .gt. 0._wp) then 
         do ibox=1,ncol ! possibility ONE
          if (frac_out(j,ibox,1) .eq. 1) then 
             prec_frac(j,ibox,1) = 1
             flag_ls=1
	  endif
	 enddo ! loop over ncol
	 if (flag_ls .eq. 0) then ! possibility THREE
	  do ibox=1,ncol
	   if (frac_out(j,ibox,2) .eq. 1) then 
               prec_frac(j,ibox,1) = 1
               flag_ls=1
	   endif
	  enddo ! loop over ncol
	 endif
	 if (flag_ls .eq. 0) then ! possibility Four
	  do ibox=1,ncol
	   if (frac_out_ls(j,ibox) .eq. 1) then 
              prec_frac(j,ibox,1) = 1
	    flag_ls=1
	   endif
	  enddo ! loop over ncol
	 endif
	 if (flag_ls .eq. 0) then ! possibility Five
	  do ibox=1,ncol
!	    prec_frac(j,1:ncol,1) = 1
            prec_frac(j,ibox,1) = 1
	  enddo ! loop over ncol
       	 endif
	endif
       ! There is large scale precipitation
	 
        if (cv_p_rate(j,1) .gt. 0._wp) then 
         do ibox=1,ncol ! possibility ONE
          if (frac_out(j,ibox,1) .eq. 2) then 
           if (prec_frac(j,ibox,1) .eq. 0) then
               prec_frac(j,ibox,1) = 2
	   else
               prec_frac(j,ibox,1) = 3
	   endif
	   flag_cv=1
	  endif
	 enddo ! loop over ncol
	 if (flag_cv .eq. 0) then ! possibility THREE
	  do ibox=1,ncol
	   if (frac_out(j,ibox,2) .eq. 2) then 
            if (prec_frac(j,ibox,1) .eq. 0) then
	     prec_frac(j,ibox,1) = 2
	    else
	     prec_frac(j,ibox,1) = 3
	    endif
	    flag_cv=1
	   endif
	  enddo ! loop over ncol
	 endif
	 if (flag_cv .eq. 0) then ! possibility Four
	  do ibox=1,ncol
	   if (frac_out_cv(j,ibox) .eq. 1) then 
            if (prec_frac(j,ibox,1) .eq. 0) then
                prec_frac(j,ibox,1) = 2
	    else
                prec_frac(j,ibox,1) = 3
	    endif
	    flag_cv=1
	   endif
	  enddo ! loop over ncol
	 endif
	 if (flag_cv .eq. 0) then  ! possibility Five
	  do ibox=1,cv_col
            if (prec_frac(j,ibox,1) .eq. 0) then
                prec_frac(j,ibox,1) = 2
	    else
                prec_frac(j,ibox,1) = 3
	    endif 
	  enddo !loop over cv_col
       	 endif 
	endif 
       ! There is convective precipitation
	
       enddo ! loop over npoints
!      end of initializing the top layer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     working on the levels from top to surface
      do ilev=2,nlev
       do j=1,npoints
        flag_ls=0
	flag_cv=0
	
        if (ls_p_rate(j,ilev) .gt. 0._wp) then 
         do ibox=1,ncol ! possibility ONE&TWO
          if ((frac_out(j,ibox,ilev) .eq. 1) .or. &
             ((prec_frac(j,ibox,ilev-1) .eq. 1) .or. &
             (prec_frac(j,ibox,ilev-1) .eq. 3))) then 
              prec_frac(j,ibox,ilev) = 1
              flag_ls=1
          endif
	 enddo ! loop over ncol
	 if ((flag_ls .eq. 0) .and. (ilev .lt. nlev)) then ! possibility THREE
	  do ibox=1,ncol
	   if (frac_out(j,ibox,ilev+1) .eq. 1) then 
               prec_frac(j,ibox,ilev) = 1
               flag_ls=1
	   endif
	  enddo ! loop over ncol
	 endif
	 if (flag_ls .eq. 0) then ! possibility Four
	  do ibox=1,ncol
	   if (frac_out_ls(j,ibox) .eq. 1) then 
	    prec_frac(j,ibox,ilev) = 1
	    flag_ls=1
	   endif
	  enddo ! loop over ncol
	 endif
	 if (flag_ls .eq. 0) then ! possibility Five
	  do ibox=1,ncol
!	  prec_frac(j,1:ncol,ilev) = 1
	  prec_frac(j,ibox,ilev) = 1
	  enddo ! loop over ncol
       	 endif
	endif ! There is large scale precipitation
	
        if (cv_p_rate(j,ilev) .gt. 0._wp) then 
         do ibox=1,ncol ! possibility ONE&TWO
          if ((frac_out(j,ibox,ilev) .eq. 2) .or. &
             ((prec_frac(j,ibox,ilev-1) .eq. 2)  .or. &
             (prec_frac(j,ibox,ilev-1) .eq. 3))) then 
            if (prec_frac(j,ibox,ilev) .eq. 0) then
	     prec_frac(j,ibox,ilev) = 2
	    else
	     prec_frac(j,ibox,ilev) = 3
	    endif 
	   flag_cv=1
          endif
	 enddo ! loop over ncol
	 if ((flag_cv .eq. 0) .and. (ilev .lt. nlev)) then ! possibility THREE
	  do ibox=1,ncol
	   if (frac_out(j,ibox,ilev+1) .eq. 2) then 
            if (prec_frac(j,ibox,ilev) .eq. 0) then
	     prec_frac(j,ibox,ilev) = 2
	    else
	     prec_frac(j,ibox,ilev) = 3
	    endif
	    flag_cv=1
	   endif
	  enddo ! loop over ncol
	 endif
	 if (flag_cv .eq. 0) then ! possibility Four
	  do ibox=1,ncol
	   if (frac_out_cv(j,ibox) .eq. 1) then 
            if (prec_frac(j,ibox,ilev) .eq. 0) then
	     prec_frac(j,ibox,ilev) = 2
	    else
	     prec_frac(j,ibox,ilev) = 3
	    endif
	    flag_cv=1
	   endif
	  enddo ! loop over ncol
	 endif
	 if (flag_cv .eq. 0) then  ! possibility Five 
	  do ibox=1,cv_col
            if (prec_frac(j,ibox,ilev) .eq. 0) then
	     prec_frac(j,ibox,ilev) = 2
	    else
	     prec_frac(j,ibox,ilev) = 3
	    endif 
	  enddo !loop over cv_col 
       	 endif 
	endif ! There is convective precipitation

       enddo ! loop over npoints
      enddo ! loop over nlev

      end subroutine prec_scops


END MODULE MO_COSP_LLNL_PRECIP


!-----------------------
