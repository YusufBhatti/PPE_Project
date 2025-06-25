      module mo_moz_subs
      use mo_kind, only : dp
      private
      public :: setrxt, phtadj, adjrxt, rxt_mod, mak_grp_vmr, set_sim_dat
      contains
      subroutine setrxt( rate, temp, m, plonl )
      use mo_moz_mods, only : plev, plnplv
      use mo_moz_mods, only : rxntot
      use mo_moz_jpl, only : jpl
      implicit none
!-------------------------------------------------------
! ... Dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: plonl
      real(dp), intent(in) :: temp(plonl,plev), m(plonl,plev)
      real(dp), intent(inout) :: rate(plonl,plev,rxntot)
!-------------------------------------------------------
! ... Local variables
!-------------------------------------------------------
      real(dp) :: itemp(plonl,plev), exp_fac(plonl,plev)
      real(dp), dimension(plonl,plev) :: ko, kinf
      rate(:,:,152) = 7.200e-11_dp
      rate(:,:,153) = 1.600e-12_dp
      rate(:,:,154) = 6.900e-12_dp
      rate(:,:,158) = 1.800e-12_dp
      rate(:,:,166) = 1.800e-12_dp
      rate(:,:,167) = 5.000e-11_dp
      rate(:,:,179) = 1.000e-11_dp
      rate(:,:,180) = 2.200e-11_dp
      rate(:,:,181) = 3.500e-12_dp
      rate(:,:,205) = 4.000e-13_dp
      rate(:,:,219) = 1.000e-12_dp
      rate(:,:,222) = 2.000e-13_dp
      rate(:,:,224) = 7.600e-14_dp
      rate(:,:,233) = 4.000e-14_dp
      rate(:,:,237) = 1.000e-11_dp
      rate(:,:,238) = 1.000e-14_dp
      rate(:,:,241) = 2.500e-12_dp
      rate(:,:,242) = 1.000e-11_dp
      rate(:,:,246) = 1.000e-11_dp
      rate(:,:,247) = 1.000e-11_dp
      rate(:,:,249) = 4.000e-12_dp
      rate(:,:,250) = 2.730e-12_dp
      rate(:,:,251) = 6.190e-12_dp
      rate(:,:,252) = 1.000e-11_dp
      rate(:,:,253) = 1.000e-11_dp
      rate(:,:,256) = 4.000e-12_dp
      rate(:,:,257) = 1.230e-11_dp
      rate(:,:,258) = 1.580e-11_dp
      rate(:,:,268) = 1.000e-11_dp
      rate(:,:,271) = 2.500e-12_dp
      rate(:,:,273) = 8.300e-13_dp
      rate(:,:,274) = 1.000e-11_dp
      rate(:,:,276) = 3.000e-12_dp
      rate(:,:,282) = 1.000e-11_dp
      rate(:,:,285) = 2.500e-12_dp
      rate(:,:,287) = 1.000e-12_dp
      rate(:,:,288) = 1.000e-11_dp
      rate(:,:,289) = 7.000e-12_dp
      rate(:,:,290) = 1.300e-13_dp
      rate(:,:,291) = 5.400e-11_dp
      rate(:,:,296) = 1.000e-12_dp
      rate(:,:,297) = 1.000e-11_dp
      rate(:,:,301) = 1.000e-12_dp
      rate(:,:,302) = 1.000e-11_dp
      rate(:,:,306) = 1.000e-11_dp
      rate(:,:,307) = 1.000e-11_dp
      rate(:,:,310) = 4.000e-12_dp
      rate(:,:,312) = 9.200e-14_dp
      rate(:,:,313) = 1.000e-11_dp
      rate(:,:,316) = 2.500e-12_dp
      rate(:,:,319) = 1.800e-11_dp
      rate(:,:,320) = 1.800e-11_dp
      rate(:,:,322) = 3.200e-11_dp
      rate(:,:,323) = 1.510e-11_dp
      rate(:,:,324) = 1.870e-11_dp
      rate(:,:,327) = 1.000e-12_dp
      rate(:,:,328) = 1.000e-11_dp
      rate(:,:,332) = 2.500e-12_dp
      rate(:,:,333) = 5.600e-12_dp
      rate(:,:,334) = 5.000e-11_dp
      rate(:,:,335) = 4.500e-12_dp
      rate(:,:,336) = 2.450e-11_dp
      rate(:,:,338) = 1.000e-11_dp
      rate(:,:,339) = 1.000e-11_dp
      rate(:,:,342) = 4.000e-12_dp
      rate(:,:,343) = 1.000e-12_dp
      rate(:,:,347) = 1.000e-12_dp
      rate(:,:,348) = 1.000e-11_dp
      rate(:,:,352) = 1.000e-12_dp
      rate(:,:,353) = 1.000e-11_dp
      rate(:,:,354) = 3.500e-12_dp
      rate(:,:,359) = 1.000e-17_dp
      rate(:,:,364) = 1.000e-12_dp
      rate(:,:,365) = 1.000e-11_dp
      rate(:,:,367) = 5.e-12_dp
      rate(:,:,368) = 2.0e-12_dp
      rate(:,:,372) = 2.500e-12_dp
      rate(:,:,373) = 2.400e-12_dp
      rate(:,:,374) = 1.000e-11_dp
      rate(:,:,378) = 1.540e-10_dp
      rate(:,:,379) = 9.300e-11_dp
      rate(:,:,380) = 6.000e-11_dp
      rate(:,:,381) = 1.500e-11_dp
      rate(:,:,382) = 1.760e-11_dp
      rate(:,:,386) = 4.000e-12_dp
      rate(:,:,387) = 1.040e-11_dp
      rate(:,:,390) = 2.300e-12_dp
      rate(:,:,394) = 2.500e-12_dp
      rate(:,:,395) = 8.000e-13_dp
      rate(:,:,396) = 1.000e-11_dp
      rate(:,:,397) = 7.500e-11_dp
      rate(:,:,398) = 3.850e-11_dp
      rate(:,:,399) = 1.360e-11_dp
      rate(:,:,403) = 2.500e-12_dp
      rate(:,:,404) = 2.900e-12_dp
      rate(:,:,405) = 1.000e-11_dp
      rate(:,:,406) = 1.180e-10_dp
      rate(:,:,407) = 7.380e-11_dp
      rate(:,:,408) = 6.100e-11_dp
      rate(:,:,411) = 2.500e-12_dp
      rate(:,:,412) = 1.300e-12_dp
      rate(:,:,413) = 1.000e-11_dp
      rate(:,:,414) = 1.030e-10_dp
      rate(:,:,415) = 4.160e-11_dp
      rate(:,:,416) = 2.400e-17_dp
      rate(:,:,421) = 1.000e-12_dp
      rate(:,:,422) = 1.000e-11_dp
      rate(:,:,423) = 2.650e-11_dp
      rate(:,:,424) = 4.520e-11_dp
      rate(:,:,425) = 2.400e-17_dp
      rate(:,:,428) = 2.500e-12_dp
      rate(:,:,430) = 1.000e-12_dp
      rate(:,:,431) = 1.000e-11_dp
      rate(:,:,432) = 3.160e-11_dp
      rate(:,:,434) = 1.000e-11_dp
      rate(:,:,435) = 1.000e-11_dp
      rate(:,:,438) = 4.000e-12_dp
      rate(:,:,440) = 2.520e-11_dp
      rate(:,:,441) = 2.880e-11_dp
      rate(:,:,443) = 2.520e-11_dp
      rate(:,:,444) = 3.810e-11_dp
      rate(:,:,445) = 1.000e-12_dp
      rate(:,:,446) = 1.000e-11_dp
      rate(:,:,448) = 2.500e-12_dp
      rate(:,:,450) = 9.700e-12_dp
      rate(:,:,453) = 1.000e-11_dp
      rate(:,:,456) = 1.400e-11_dp
      rate(:,:,459) = 1.000e-12_dp
      rate(:,:,460) = 1.000e-11_dp
      rate(:,:,461) = 1.000e-12_dp
      rate(:,:,464) = 2.400e-12_dp
      rate(:,:,465) = 1.000e-12_dp
      rate(:,:,466) = 1.000e-11_dp
      rate(:,:,470) = 1.000e-12_dp
      rate(:,:,471) = 1.000e-11_dp
      rate(:,:,472) = 5.200e-11_dp
      rate(:,:,473) = 4.720e-11_dp
      rate(:,:,477) = 2.500e-12_dp
      rate(:,:,478) = 8.000e-13_dp
      rate(:,:,479) = 1.000e-11_dp
      rate(:,:,483) = 2.500e-12_dp
      rate(:,:,484) = 8.000e-13_dp
      rate(:,:,485) = 1.000e-11_dp
      rate(:,:,487) = 2.104e-11_dp
      rate(:,:,489) = 1.515e-11_dp
      rate(:,:,490) = 8.916e-12_dp
      rate(:,:,493) = 3.800e-12_dp
      rate(:,:,496) = 1.000e-12_dp
      rate(:,:,497) = 1.000e-11_dp
      rate(:,:,499) = 2.100e-12_dp
      rate(:,:,500) = 2.800e-13_dp
      rate(:,:,502) = 2.300e-12_dp
      rate(:,:,504) = 1.000e-12_dp
      rate(:,:,505) = 1.000e-11_dp
      rate(:,:,509) = 1.000e-12_dp
      rate(:,:,510) = 1.000e-11_dp
      rate(:,:,512) = 1.000e-10_dp
      rate(:,:,513) = 9.900e-11_dp
      rate(:,:,514) = 2.100e-12_dp
      rate(:,:,515) = 2.800e-13_dp
      rate(:,:,518) = 2.300e-12_dp
      rate(:,:,519) = 1.000e-12_dp
      rate(:,:,520) = 1.000e-11_dp
      rate(:,:,523) = 4.700e-11_dp
      rate(:,:,524) = 1.400e-11_dp
      rate(:,:,527) = 1.000e-12_dp
      rate(:,:,528) = 1.000e-11_dp
      rate(:,:,532) = 1.000e-12_dp
      rate(:,:,533) = 1.000e-11_dp
      rate(:,:,540) = 1.000e-12_dp
      rate(:,:,541) = 1.000e-11_dp
      rate(:,:,542) = 1.700e-11_dp
      rate(:,:,543) = 8.400e-11_dp
      rate(:,:,544) = 3.200e-11_dp
      rate(:,:,547) = 1.000e-12_dp
      rate(:,:,548) = 1.000e-11_dp
      rate(:,:,552) = 1.000e-12_dp
      rate(:,:,553) = 1.000e-11_dp
      rate(:,:,558) = 2.100e-10_dp
      rate(:,:,559) = 2.000e-10_dp
      rate(:,:,563) = 4.700e-16_dp
      rate(:,:,564) = 1.200e-14_dp
      rate(:,:,566) = 2.500e-12_dp
      rate(:,:,567) = 1.100e-11_dp
      rate(:,:,568) = 1.200e-11_dp
      rate(:,:,569) = 1.900e-11_dp
      rate(:,:,573) = 1.000e-11_dp
      rate(:,:,574) = 3.300e-11_dp
      rate(:,:,575) = 5.700e-11_dp
      rate(:,:,576) = 1.000e-12_dp
      rate(:,:,577) = 3.5e-12_dp
      rate(:,:,581) = 1.000e-11_dp
      rate(:,:,582) = 2.300e-11_dp
      rate(:,:,583) = 3.400e-11_dp
      rate(:,:,587) = 1.000e-11_dp
      rate(:,:,588) = 2.400e-12_dp
      rate(:,:,589) = 3.5e-12_dp
      rate(:,:,590) = 1.000e-11_dp
      rate(:,:,593) = 1.200e-10_dp
      rate(:,:,594) = 2.020e-10_dp
      rate(:,:,595) = 1.204e-10_dp
      rate(:,:,596) = 1.500e-10_dp
      rate(:,:,597) = 9.750e-11_dp
      rate(:,:,598) = 1.500e-11_dp
      rate(:,:,599) = 7.200e-11_dp
      rate(:,:,600) = 1.794e-10_dp
      rate(:,:,601) = 1.628e-10_dp
      rate(:,:,602) = 2.840e-10_dp
      rate(:,:,603) = 1.674e-10_dp
      rate(:,:,604) = 9.600e-11_dp
      rate(:,:,605) = 4.100e-11_dp
      rate(:,:,606) = 1.012e-10_dp
      rate(:,:,607) = 1.200e-10_dp
      rate(:,:,608) = 4.490e-10_dp
      rate(:,:,609) = 2.570e-10_dp
      rate(:,:,610) = 2.140e-11_dp
      rate(:,:,611) = 1.900e-10_dp
      rate(:,:,612) = 1.310e-10_dp
      rate(:,:,613) = 3.500e-11_dp
      rate(:,:,614) = 9.000e-12_dp
      rate(:,:,615) = 1.200e-10_dp
      rate(:,:,616) = 1.500e-10_dp
      rate(:,:,617) = 1.200e-10_dp
      rate(:,:,637) = 1.700e-13_dp
      rate(:,:,661) = 5.805e-11_dp
      rate(:,:,663) = 1.400e-11_dp
      rate(:,:,673) = 1.600e-10_dp
      rate(:,:,675) = 5.900e-11_dp
      rate(:,:,676) = 8.00e-11_dp
      rate(:,:,677) = 7.600e-11_dp
      rate(:,:,678) = 3.441e-11_dp
      rate(:,:,679) = 1.400e-10_dp
      rate(:,:,681) = 5.400e-11_dp
      rate(:,:,682) = 1.935e-10_dp
      rate(:,:,692) = 3.3e-10_dp
      rate(:,:,694) = 4.4e-13_dp
      rate(:,:,695) = 1.e-10_dp
      rate(:,:,697) = 3.e-13_dp
      rate(:,:,698) = 5.e-11_dp
      itemp(:,:) = 1._dp / temp(:,:)
      rate(:,:,144) = 8.000e-12_dp * exp( -2060._dp * itemp(:,:) )
      rate(:,:,146) = 2.150e-11_dp * exp( 110._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 55._dp * itemp(:,:) )
      rate(:,:,147) = 3.135e-11_dp * exp_fac(:,:)
      rate(:,:,148) = 1.650e-12_dp * exp_fac(:,:)
      rate(:,:,149) = 1.630e-10_dp * exp( 60._dp * itemp(:,:) )
      rate(:,:,151) = 1.400e-10_dp * exp( -470._dp * itemp(:,:) )
      rate(:,:,155) = 1.600e-11_dp * exp( -4570._dp * itemp(:,:) )
      rate(:,:,156) = 2.800e-12_dp * exp( -1800._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 180._dp * itemp(:,:) )
      rate(:,:,157) = 1.800e-11_dp * exp_fac(:,:)
      rate(:,:,234) = 4.200e-12_dp * exp_fac(:,:)
      rate(:,:,270) = 4.200e-12_dp * exp_fac(:,:)
      rate(:,:,293) = 4.032e-12_dp * exp_fac(:,:)
      rate(:,:,294) = 1.680e-13_dp * exp_fac(:,:)
      rate(:,:,299) = 4.200e-12_dp * exp_fac(:,:)
      rate(:,:,361) = 3.780e-12_dp * exp_fac(:,:)
      rate(:,:,362) = 0.420e-12_dp * exp_fac(:,:)
      rate(:,:,570) = 4.200e-12_dp * exp_fac(:,:)
      rate(:,:,578) = 4.200e-12_dp * exp_fac(:,:)
      rate(:,:,584) = 4.200e-12_dp * exp_fac(:,:)
      rate(:,:,160) = 1.700e-12_dp * exp( -940._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 200._dp * itemp(:,:) )
      rate(:,:,161) = 3.000e-11_dp * exp_fac(:,:)
      rate(:,:,204) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,225) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,269) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,275) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,283) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,298) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,366) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,455) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,486) = 7.600e-12_dp * exp_fac(:,:)
      rate(:,:,488) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,498) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,506) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,511) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,529) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,534) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,549) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,554) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,657) = 5.500e-12_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( 250._dp * itemp(:,:) )
      rate(:,:,162) = 4.800e-11_dp * exp_fac(:,:)
      rate(:,:,648) = 1.700e-11_dp * exp_fac(:,:)
      rate(:,:,165) = 1.400e-12_dp * exp( -2000._dp * itemp(:,:) )
      rate(:,:,168) = 1.500e-11_dp * exp( -3600._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 100._dp * itemp(:,:) )
      rate(:,:,169) = 2.100e-11_dp * exp_fac(:,:)
      rate(:,:,618) = 7.700e-11_dp * exp_fac(:,:)
      rate(:,:,170) = 5.800e-12_dp * exp( 220._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -1500._dp * itemp(:,:) )
      rate(:,:,172) = 3.000e-12_dp * exp_fac(:,:)
      rate(:,:,658) = 5.800e-12_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( 270._dp * itemp(:,:) )
      rate(:,:,173) = 3.300e-12_dp * exp_fac(:,:)
      rate(:,:,248) = 8.100e-12_dp * exp_fac(:,:)
      rate(:,:,255) = 8.100e-12_dp * exp_fac(:,:)
      rate(:,:,341) = 8.100e-12_dp * exp_fac(:,:)
      rate(:,:,437) = 8.100e-12_dp * exp_fac(:,:)
      rate(:,:,622) = 1.400e-11_dp * exp_fac(:,:)
      rate(:,:,626) = 7.400e-12_dp * exp_fac(:,:)
      rate(:,:,174) = 5.100e-12_dp * exp( 210._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -2450._dp * itemp(:,:) )
      rate(:,:,176) = 1.200e-13_dp * exp_fac(:,:)
      rate(:,:,631) = 3.000e-11_dp * exp_fac(:,:)
      rate(:,:,177) = 4.000e-10_dp * exp( -340._dp * itemp(:,:) )
      rate(:,:,182) = 1.500e-11_dp * exp( 170._dp * itemp(:,:) )
      rate(:,:,185) = 1.800e-11_dp * exp( -390._dp * itemp(:,:) )
      rate(:,:,187) = 1.300e-12_dp * exp( 380._dp * itemp(:,:) )
      rate(:,:,191) = 1.700E-12_dp * exp( -710._dp * itemp(:,:) )
      rate(:,:,194) = 2.450e-12_dp * exp( -1775._dp * itemp(:,:) )
      rate(:,:,195) = 2.900e-12_dp * exp( -345._dp * itemp(:,:) )
      rate(:,:,196) = 1.960e-12_dp * exp( 403._dp * itemp(:,:) )
      rate(:,:,197) = 3.800e-13_dp * exp( 730._dp * itemp(:,:) )
      rate(:,:,198) = 7.400e-13_dp * exp( -520._dp * itemp(:,:) )
      rate(:,:,199) = 2.330e-14_dp * exp( 678._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -1600._dp * itemp(:,:) )
      rate(:,:,200) = 3.400e-11_dp * exp_fac(:,:)
      rate(:,:,708) = 1.050e-12_dp * exp_fac(:,:)
      rate(:,:,709) = 1.250e-12_dp * exp_fac(:,:)
      rate(:,:,201) = 5.500e-12_dp * exp( 125._dp * itemp(:,:) )
      rate(:,:,202) = 9.700e-15_dp * exp( 625._dp * itemp(:,:) )
      rate(:,:,203) = 6.000e-13_dp * exp( -2058._dp * itemp(:,:) )
      rate(:,:,206) = 2.400e+12_dp * exp( -7000._dp * itemp(:,:) )
      rate(:,:,207) = 2.600e-12_dp * exp( 265._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 700._dp * itemp(:,:) )
      rate(:,:,208) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,221) = 7.400e-13_dp * exp_fac(:,:)
      rate(:,:,235) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,266) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,272) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,280) = 8.600e-13_dp * exp_fac(:,:)
      rate(:,:,295) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,300) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,363) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,454) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,495) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,503) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,508) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,516) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,525) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,530) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,546) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,550) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,571) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,579) = 7.500e-13_dp * exp_fac(:,:)
      rate(:,:,585) = 7.500e-13_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( -400._dp * itemp(:,:) )
      rate(:,:,209) = 1.200e-13_dp * exp_fac(:,:)
      rate(:,:,360) = 4.600e-14_dp * exp_fac(:,:)
      rate(:,:,214) = 9.100e-15_dp * exp( -2580._dp * itemp(:,:) )
      rate(:,:,215) = 7.660e-12_dp * exp( -1020._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 0._dp * itemp(:,:) )
      rate(:,:,216) = 3.350e-12_dp * exp_fac(:,:)
      rate(:,:,219) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,222) = 2.000e-13_dp * exp_fac(:,:)
      rate(:,:,224) = 7.600e-14_dp * exp_fac(:,:)
      rate(:,:,233) = 4.000e-14_dp * exp_fac(:,:)
      rate(:,:,237) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,238) = 1.000e-14_dp * exp_fac(:,:)
      rate(:,:,241) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,242) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,246) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,247) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,249) = 4.000e-12_dp * exp_fac(:,:)
      rate(:,:,250) = 2.730e-12_dp * exp_fac(:,:)
      rate(:,:,251) = 6.190e-12_dp * exp_fac(:,:)
      rate(:,:,252) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,253) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,256) = 4.000e-12_dp * exp_fac(:,:)
      rate(:,:,257) = 1.230e-11_dp * exp_fac(:,:)
      rate(:,:,258) = 1.580e-11_dp * exp_fac(:,:)
      rate(:,:,268) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,271) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,273) = 8.300e-13_dp * exp_fac(:,:)
      rate(:,:,274) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,276) = 3.000e-12_dp * exp_fac(:,:)
      rate(:,:,282) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,285) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,287) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,288) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,289) = 7.000e-12_dp * exp_fac(:,:)
      rate(:,:,290) = 1.300e-13_dp * exp_fac(:,:)
      rate(:,:,291) = 5.400e-11_dp * exp_fac(:,:)
      rate(:,:,296) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,297) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,301) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,302) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,306) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,307) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,310) = 4.000e-12_dp * exp_fac(:,:)
      rate(:,:,312) = 9.200e-14_dp * exp_fac(:,:)
      rate(:,:,313) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,316) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,319) = 1.800e-11_dp * exp_fac(:,:)
      rate(:,:,320) = 1.800e-11_dp * exp_fac(:,:)
      rate(:,:,322) = 3.200e-11_dp * exp_fac(:,:)
      rate(:,:,323) = 1.510e-11_dp * exp_fac(:,:)
      rate(:,:,324) = 1.870e-11_dp * exp_fac(:,:)
      rate(:,:,327) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,328) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,332) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,333) = 5.600e-12_dp * exp_fac(:,:)
      rate(:,:,334) = 5.000e-11_dp * exp_fac(:,:)
      rate(:,:,335) = 4.500e-12_dp * exp_fac(:,:)
      rate(:,:,336) = 2.450e-11_dp * exp_fac(:,:)
      rate(:,:,338) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,339) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,342) = 4.000e-12_dp * exp_fac(:,:)
      rate(:,:,343) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,347) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,348) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,352) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,353) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,354) = 3.500e-12_dp * exp_fac(:,:)
      rate(:,:,359) = 1.000e-17_dp * exp_fac(:,:)
      rate(:,:,364) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,365) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,367) = 5.e-12_dp * exp_fac(:,:)
      rate(:,:,368) = 2.0e-12_dp * exp_fac(:,:)
      rate(:,:,372) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,373) = 2.400e-12_dp * exp_fac(:,:)
      rate(:,:,374) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,378) = 1.540e-10_dp * exp_fac(:,:)
      rate(:,:,379) = 9.300e-11_dp * exp_fac(:,:)
      rate(:,:,380) = 6.000e-11_dp * exp_fac(:,:)
      rate(:,:,381) = 1.500e-11_dp * exp_fac(:,:)
      rate(:,:,382) = 1.760e-11_dp * exp_fac(:,:)
      rate(:,:,386) = 4.000e-12_dp * exp_fac(:,:)
      rate(:,:,387) = 1.040e-11_dp * exp_fac(:,:)
      rate(:,:,390) = 2.300e-12_dp * exp_fac(:,:)
      rate(:,:,394) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,395) = 8.000e-13_dp * exp_fac(:,:)
      rate(:,:,396) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,397) = 7.500e-11_dp * exp_fac(:,:)
      rate(:,:,398) = 3.850e-11_dp * exp_fac(:,:)
      rate(:,:,399) = 1.360e-11_dp * exp_fac(:,:)
      rate(:,:,403) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,404) = 2.900e-12_dp * exp_fac(:,:)
      rate(:,:,405) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,406) = 1.180e-10_dp * exp_fac(:,:)
      rate(:,:,407) = 7.380e-11_dp * exp_fac(:,:)
      rate(:,:,408) = 6.100e-11_dp * exp_fac(:,:)
      rate(:,:,411) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,412) = 1.300e-12_dp * exp_fac(:,:)
      rate(:,:,413) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,414) = 1.030e-10_dp * exp_fac(:,:)
      rate(:,:,415) = 4.160e-11_dp * exp_fac(:,:)
      rate(:,:,416) = 2.400e-17_dp * exp_fac(:,:)
      rate(:,:,421) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,422) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,423) = 2.650e-11_dp * exp_fac(:,:)
      rate(:,:,424) = 4.520e-11_dp * exp_fac(:,:)
      rate(:,:,425) = 2.400e-17_dp * exp_fac(:,:)
      rate(:,:,428) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,430) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,431) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,432) = 3.160e-11_dp * exp_fac(:,:)
      rate(:,:,434) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,435) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,438) = 4.000e-12_dp * exp_fac(:,:)
      rate(:,:,440) = 2.520e-11_dp * exp_fac(:,:)
      rate(:,:,441) = 2.880e-11_dp * exp_fac(:,:)
      rate(:,:,443) = 2.520e-11_dp * exp_fac(:,:)
      rate(:,:,444) = 3.810e-11_dp * exp_fac(:,:)
      rate(:,:,445) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,446) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,448) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,450) = 9.700e-12_dp * exp_fac(:,:)
      rate(:,:,453) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,456) = 1.400e-11_dp * exp_fac(:,:)
      rate(:,:,459) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,460) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,461) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,464) = 2.400e-12_dp * exp_fac(:,:)
      rate(:,:,465) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,466) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,470) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,471) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,472) = 5.200e-11_dp * exp_fac(:,:)
      rate(:,:,473) = 4.720e-11_dp * exp_fac(:,:)
      rate(:,:,477) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,478) = 8.000e-13_dp * exp_fac(:,:)
      rate(:,:,479) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,483) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,484) = 8.000e-13_dp * exp_fac(:,:)
      rate(:,:,485) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,487) = 2.104e-11_dp * exp_fac(:,:)
      rate(:,:,489) = 1.515e-11_dp * exp_fac(:,:)
      rate(:,:,490) = 8.916e-12_dp * exp_fac(:,:)
      rate(:,:,493) = 3.800e-12_dp * exp_fac(:,:)
      rate(:,:,496) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,497) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,499) = 2.100e-12_dp * exp_fac(:,:)
      rate(:,:,500) = 2.800e-13_dp * exp_fac(:,:)
      rate(:,:,502) = 2.300e-12_dp * exp_fac(:,:)
      rate(:,:,504) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,505) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,509) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,510) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,512) = 1.000e-10_dp * exp_fac(:,:)
      rate(:,:,513) = 9.900e-11_dp * exp_fac(:,:)
      rate(:,:,514) = 2.100e-12_dp * exp_fac(:,:)
      rate(:,:,515) = 2.800e-13_dp * exp_fac(:,:)
      rate(:,:,518) = 2.300e-12_dp * exp_fac(:,:)
      rate(:,:,519) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,520) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,523) = 4.700e-11_dp * exp_fac(:,:)
      rate(:,:,524) = 1.400e-11_dp * exp_fac(:,:)
      rate(:,:,527) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,528) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,532) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,533) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,540) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,541) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,542) = 1.700e-11_dp * exp_fac(:,:)
      rate(:,:,543) = 8.400e-11_dp * exp_fac(:,:)
      rate(:,:,544) = 3.200e-11_dp * exp_fac(:,:)
      rate(:,:,547) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,548) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,552) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,553) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,558) = 2.100e-10_dp * exp_fac(:,:)
      rate(:,:,559) = 2.000e-10_dp * exp_fac(:,:)
      rate(:,:,563) = 4.700e-16_dp * exp_fac(:,:)
      rate(:,:,564) = 1.200e-14_dp * exp_fac(:,:)
      rate(:,:,566) = 2.500e-12_dp * exp_fac(:,:)
      rate(:,:,567) = 1.100e-11_dp * exp_fac(:,:)
      rate(:,:,568) = 1.200e-11_dp * exp_fac(:,:)
      rate(:,:,569) = 1.900e-11_dp * exp_fac(:,:)
      rate(:,:,573) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,574) = 3.300e-11_dp * exp_fac(:,:)
      rate(:,:,575) = 5.700e-11_dp * exp_fac(:,:)
      rate(:,:,576) = 1.000e-12_dp * exp_fac(:,:)
      rate(:,:,577) = 3.5e-12_dp * exp_fac(:,:)
      rate(:,:,581) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,582) = 2.300e-11_dp * exp_fac(:,:)
      rate(:,:,583) = 3.400e-11_dp * exp_fac(:,:)
      rate(:,:,587) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,588) = 2.400e-12_dp * exp_fac(:,:)
      rate(:,:,589) = 3.5e-12_dp * exp_fac(:,:)
      rate(:,:,590) = 1.000e-11_dp * exp_fac(:,:)
      rate(:,:,593) = 1.200e-10_dp * exp_fac(:,:)
      rate(:,:,594) = 2.020e-10_dp * exp_fac(:,:)
      rate(:,:,595) = 1.204e-10_dp * exp_fac(:,:)
      rate(:,:,596) = 1.500e-10_dp * exp_fac(:,:)
      rate(:,:,597) = 9.750e-11_dp * exp_fac(:,:)
      rate(:,:,598) = 1.500e-11_dp * exp_fac(:,:)
      rate(:,:,599) = 7.200e-11_dp * exp_fac(:,:)
      rate(:,:,600) = 1.794e-10_dp * exp_fac(:,:)
      rate(:,:,601) = 1.628e-10_dp * exp_fac(:,:)
      rate(:,:,602) = 2.840e-10_dp * exp_fac(:,:)
      rate(:,:,603) = 1.674e-10_dp * exp_fac(:,:)
      rate(:,:,604) = 9.600e-11_dp * exp_fac(:,:)
      rate(:,:,605) = 4.100e-11_dp * exp_fac(:,:)
      rate(:,:,606) = 1.012e-10_dp * exp_fac(:,:)
      rate(:,:,607) = 1.200e-10_dp * exp_fac(:,:)
      rate(:,:,608) = 4.490e-10_dp * exp_fac(:,:)
      rate(:,:,609) = 2.570e-10_dp * exp_fac(:,:)
      rate(:,:,610) = 2.140e-11_dp * exp_fac(:,:)
      rate(:,:,611) = 1.900e-10_dp * exp_fac(:,:)
      rate(:,:,612) = 1.310e-10_dp * exp_fac(:,:)
      rate(:,:,613) = 3.500e-11_dp * exp_fac(:,:)
      rate(:,:,614) = 9.000e-12_dp * exp_fac(:,:)
      rate(:,:,615) = 1.200e-10_dp * exp_fac(:,:)
      rate(:,:,616) = 1.500e-10_dp * exp_fac(:,:)
      rate(:,:,617) = 1.200e-10_dp * exp_fac(:,:)
      rate(:,:,637) = 1.700e-13_dp * exp_fac(:,:)
      rate(:,:,661) = 5.805e-11_dp * exp_fac(:,:)
      rate(:,:,663) = 1.400e-11_dp * exp_fac(:,:)
      rate(:,:,673) = 1.600e-10_dp * exp_fac(:,:)
      rate(:,:,675) = 5.900e-11_dp * exp_fac(:,:)
      rate(:,:,676) = 8.00e-11_dp * exp_fac(:,:)
      rate(:,:,677) = 7.600e-11_dp * exp_fac(:,:)
      rate(:,:,678) = 3.441e-11_dp * exp_fac(:,:)
      rate(:,:,679) = 1.400e-10_dp * exp_fac(:,:)
      rate(:,:,681) = 5.400e-11_dp * exp_fac(:,:)
      rate(:,:,682) = 1.935e-10_dp * exp_fac(:,:)
      rate(:,:,692) = 3.3e-10_dp * exp_fac(:,:)
      rate(:,:,694) = 4.4e-13_dp * exp_fac(:,:)
      rate(:,:,695) = 1.e-10_dp * exp_fac(:,:)
      rate(:,:,697) = 3.e-13_dp * exp_fac(:,:)
      rate(:,:,698) = 5.e-11_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( 350._dp * itemp(:,:) )
      rate(:,:,217) = 4.630e-12_dp * exp_fac(:,:)
      rate(:,:,265) = 2.900e-12_dp * exp_fac(:,:)
      rate(:,:,218) = 1.400e-12_dp * exp( -1900._dp * itemp(:,:) )
      rate(:,:,220) = 2.620e-12_dp * exp( 373._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 500._dp * itemp(:,:) )
      rate(:,:,223) = 1.800e-12_dp * exp_fac(:,:)
      rate(:,:,229) = 2.000e-12_dp * exp_fac(:,:)
      rate(:,:,230) = 2.900e-12_dp * exp_fac(:,:)
      rate(:,:,281) = 7.100e-13_dp * exp_fac(:,:)
      rate(:,:,572) = 2.000e-12_dp * exp_fac(:,:)
      rate(:,:,580) = 2.000e-12_dp * exp_fac(:,:)
      rate(:,:,586) = 2.000e-12_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( 290._dp * itemp(:,:) )
      rate(:,:,226) = 7.500e-12_dp * exp_fac(:,:)
      rate(:,:,309) = 8.700e-12_dp * exp_fac(:,:)
      rate(:,:,345) = 7.500e-12_dp * exp_fac(:,:)
      rate(:,:,350) = 7.500e-12_dp * exp_fac(:,:)
      rate(:,:,385) = 7.500e-12_dp * exp_fac(:,:)
      rate(:,:,468) = 7.500e-12_dp * exp_fac(:,:)
      rate(:,:,538) = 7.500e-12_dp * exp_fac(:,:)
      rate(:,:,628) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,629) = 6.400e-12_dp * exp_fac(:,:)
      rate(:,:,654) = 4.100e-13_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( 980._dp * itemp(:,:) )
      rate(:,:,228) = 5.200e-13_dp * exp_fac(:,:)
      rate(:,:,384) = 5.200e-13_dp * exp_fac(:,:)
      rate(:,:,231) = 3.150e-14_dp * exp( 920._dp * itemp(:,:) )
      rate(:,:,236) = 4.000e-12_dp * exp( 1000._dp * itemp(:,:) )
      rate(:,:,239) = 1.600e+11_dp * exp( -4150._dp * itemp(:,:) )
      rate(:,:,240) = 3.100e-12_dp * exp( 340._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -1862._dp * itemp(:,:) )
      rate(:,:,243) = 1.440e-12_dp * exp_fac(:,:)
      rate(:,:,305) = 2.880e-12_dp * exp_fac(:,:)
      rate(:,:,337) = 5.760e-12_dp * exp_fac(:,:)
      rate(:,:,417) = 6.120e-12_dp * exp_fac(:,:)
      rate(:,:,426) = 6.120e-12_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( 1040._dp * itemp(:,:) )
      rate(:,:,245) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,254) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,308) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,340) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,346) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,349) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,436) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,458) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,462) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,467) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,539) = 4.300e-13_dp * exp_fac(:,:)
      rate(:,:,259) = 7.800e-13_dp * exp( -1050._dp * itemp(:,:) )
      rate(:,:,261) = 5.500e-15_dp * exp( -1880._dp * itemp(:,:) )
      rate(:,:,262) = 4.600e-13_dp * exp( -1156._dp * itemp(:,:) )
      rate(:,:,263) = 7.600e-12_dp * exp( -585._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -40._dp * itemp(:,:) )
      rate(:,:,267) = 3.750e-13_dp * exp_fac(:,:)
      rate(:,:,452) = 3.750e-13_dp * exp_fac(:,:)
      rate(:,:,277) = 8.400e-13_dp * exp( 830._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -1860._dp * itemp(:,:) )
      rate(:,:,278) = 1.400e-12_dp * exp_fac(:,:)
      rate(:,:,383) = 1.050e-11_dp * exp_fac(:,:)
      rate(:,:,279) = 2.900e-12_dp * exp( 300._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 360._dp * itemp(:,:) )
      rate(:,:,284) = 2.540e-12_dp * exp_fac(:,:)
      rate(:,:,286) = 1.320e-12_dp * exp_fac(:,:)
      rate(:,:,314) = 2.464e-12_dp * exp_fac(:,:)
      rate(:,:,315) = 7.620e-14_dp * exp_fac(:,:)
      rate(:,:,330) = 2.438e-12_dp * exp_fac(:,:)
      rate(:,:,331) = 1.020e-13_dp * exp_fac(:,:)
      rate(:,:,370) = 2.235e-12_dp * exp_fac(:,:)
      rate(:,:,371) = 0.305e-12_dp * exp_fac(:,:)
      rate(:,:,389) = 2.700e-12_dp * exp_fac(:,:)
      rate(:,:,392) = 2.235e-12_dp * exp_fac(:,:)
      rate(:,:,393) = 0.305e-12_dp * exp_fac(:,:)
      rate(:,:,401) = 2.235e-12_dp * exp_fac(:,:)
      rate(:,:,402) = 0.305e-12_dp * exp_fac(:,:)
      rate(:,:,410) = 2.540e-12_dp * exp_fac(:,:)
      rate(:,:,419) = 4.270e-12_dp * exp_fac(:,:)
      rate(:,:,420) = 3.302e-12_dp * exp_fac(:,:)
      rate(:,:,427) = 2.540e-12_dp * exp_fac(:,:)
      rate(:,:,447) = 2.540e-12_dp * exp_fac(:,:)
      rate(:,:,475) = 2.540e-12_dp * exp_fac(:,:)
      rate(:,:,476) = 0.305e-12_dp * exp_fac(:,:)
      rate(:,:,481) = 2.235e-12_dp * exp_fac(:,:)
      rate(:,:,482) = 0.305e-12_dp * exp_fac(:,:)
      rate(:,:,292) = 2.300e-12_dp * exp( -170._dp * itemp(:,:) )
      rate(:,:,303) = 1.860e-11_dp * exp( 175._dp * itemp(:,:) )
      rate(:,:,304) = 1.360e-15_dp * exp( -2112._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 1300._dp * itemp(:,:) )
      rate(:,:,317) = 1.820e-13_dp * exp_fac(:,:)
      rate(:,:,329) = 1.820e-13_dp * exp_fac(:,:)
      rate(:,:,369) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,388) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,391) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,400) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,409) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,418) = 1.930e-13_dp * exp_fac(:,:)
      rate(:,:,429) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,449) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,474) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,480) = 2.050e-13_dp * exp_fac(:,:)
      rate(:,:,318) = 2.900e+07_dp * exp( -5297._dp * itemp(:,:) )
      rate(:,:,325) = 4.130e-12_dp * exp( 452._dp * itemp(:,:) )
      rate(:,:,326) = 7.510e-16_dp * exp( -1521._dp * itemp(:,:) )
      rate(:,:,355) = 2.700e-11_dp * exp( 390._dp * itemp(:,:) )
      rate(:,:,356) = 7.860e-15_dp * exp( -1913._dp * itemp(:,:) )
      rate(:,:,357) = 3.030e-12_dp * exp( -446._dp * itemp(:,:) )
      rate(:,:,358) = 8.100e-12_dp * exp( 610._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -7700._dp * itemp(:,:) )
      rate(:,:,375) = 4.1e+08_dp * exp_fac(:,:)
      rate(:,:,376) = 4.1e+08_dp * exp_fac(:,:)
      rate(:,:,377) = 4.1e+08_dp * exp_fac(:,:)
      rate(:,:,433) = 4.1e+08_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( 365._dp * itemp(:,:) )
      rate(:,:,451) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,457) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,463) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,494) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,501) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,507) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,517) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,526) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,531) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,545) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,551) = 2.600e-12_dp * exp_fac(:,:)
      rate(:,:,491) = 2.300e-12_dp * exp( -193._dp * itemp(:,:) )
      rate(:,:,492) = 4.700e-13_dp * exp( 1220._dp * itemp(:,:) )
      rate(:,:,521) = 1.900e-12_dp * exp( 190._dp * itemp(:,:) )
      rate(:,:,522) = 1.700e-12_dp * exp( 352._dp * itemp(:,:) )
      rate(:,:,535) = 5.900e-12_dp * exp( 225._dp * itemp(:,:) )
      rate(:,:,555) = 1.200e-11_dp * exp( 440._dp * itemp(:,:) )
      rate(:,:,556) = 1.600e-11_dp * exp( 470._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 400._dp * itemp(:,:) )
      rate(:,:,557) = 4.200e-11_dp * exp_fac(:,:)
      rate(:,:,666) = 6.000e-12_dp * exp_fac(:,:)
      rate(:,:,560) = 6.300e-16_dp * exp( -580._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -1300._dp * itemp(:,:) )
      rate(:,:,561) = 1.700e-15_dp * exp_fac(:,:)
      rate(:,:,699) = 2.350e-12_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( -780._dp * itemp(:,:) )
      rate(:,:,562) = 3.000e-15_dp * exp_fac(:,:)
      rate(:,:,643) = 1.600e-11_dp * exp_fac(:,:)
      rate(:,:,565) = 1.200e-12_dp * exp( 490._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 20._dp * itemp(:,:) )
      rate(:,:,591) = 7.250e-11_dp * exp_fac(:,:)
      rate(:,:,592) = 4.630e-11_dp * exp_fac(:,:)
      rate(:,:,619) = 2.300e-11_dp * exp( -200._dp * itemp(:,:) )
      rate(:,:,620) = 3.050e-11_dp * exp( -2270._dp * itemp(:,:) )
      rate(:,:,621) = 1.100e-11_dp * exp( -980._dp * itemp(:,:) )
      rate(:,:,623) = 3.600e-11_dp * exp( -375._dp * itemp(:,:) )
      rate(:,:,625) = 2.800e-11_dp * exp( 85._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 230._dp * itemp(:,:) )
      rate(:,:,627) = 6.000e-13_dp * exp_fac(:,:)
      rate(:,:,647) = 1.900e-11_dp * exp_fac(:,:)
      rate(:,:,632) = 1.000e-12_dp * exp( -1590._dp * itemp(:,:) )
      rate(:,:,633) = 3.500e-13_dp * exp( -1370._dp * itemp(:,:) )
      rate(:,:,635) = 1.800e-12_dp * exp( -250._dp * itemp(:,:) )
      rate(:,:,636) = 1.000e-11_dp * exp( -3300._dp * itemp(:,:) )
      rate(:,:,638) = 3.400e-12_dp * exp( -130._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -500._dp * itemp(:,:) )
      rate(:,:,639) = 3.000e-12_dp * exp_fac(:,:)
      rate(:,:,664) = 1.400e-10_dp * exp_fac(:,:)
      exp_fac(:,:) = exp( -840._dp * itemp(:,:) )
      rate(:,:,640) = 3.600e-12_dp * exp_fac(:,:)
      rate(:,:,701) = 2.000e-12_dp * exp_fac(:,:)
      rate(:,:,641) = 1.200e-12_dp * exp( -330._dp * itemp(:,:) )
      rate(:,:,642) = 6.500e-12_dp * exp( 135._dp * itemp(:,:) )
      rate(:,:,644) = 4.800e-12_dp * exp( -310._dp * itemp(:,:) )
      rate(:,:,646) = 1.648e+11_dp * exp( -7399._dp * itemp(:,:) )
      rate(:,:,649) = 4.500e-12_dp * exp( 460._dp * itemp(:,:) )
      exp_fac(:,:) = exp( 260._dp * itemp(:,:) )
      rate(:,:,650) = 8.800e-12_dp * exp_fac(:,:)
      rate(:,:,653) = 2.300e-12_dp * exp_fac(:,:)
      rate(:,:,652) = 9.500e-13_dp * exp( 550._dp * itemp(:,:) )
      rate(:,:,655) = 2.400e-12_dp * exp( 40._dp * itemp(:,:) )
      rate(:,:,656) = 2.800e-14_dp * exp( 860._dp * itemp(:,:) )
      rate(:,:,659) = 1.200e-10_dp * exp( -430._dp * itemp(:,:) )
      rate(:,:,660) = 1.900e-11_dp * exp( 215._dp * itemp(:,:) )
      rate(:,:,662) = 2.100e-11_dp * exp( 240._dp * itemp(:,:) )
      rate(:,:,665) = 1.600e-10_dp * exp( -260._dp * itemp(:,:) )
      rate(:,:,667) = 8.100e-11_dp * exp( -30._dp * itemp(:,:) )
      rate(:,:,668) = 7.300e-12_dp * exp( -1280._dp * itemp(:,:) )
      rate(:,:,669) = 1.6e-11_dp * exp( -2140._dp * itemp(:,:) )
      rate(:,:,672) = 7.20e-11_dp * exp( -70._dp * itemp(:,:) )
      rate(:,:,674) = 7.100e-11_dp * exp( -75._dp * itemp(:,:) )
      rate(:,:,680) = 1.5e-11_dp * exp( -590._dp * itemp(:,:) )
      rate(:,:,683) = 3.800e-11_dp * exp( 16._dp * itemp(:,:) )
      rate(:,:,684) = 3.300e-12_dp * exp( -115._dp * itemp(:,:) )
      exp_fac(:,:) = exp( -800._dp * itemp(:,:) )
      rate(:,:,685) = 1.700e-11_dp * exp_fac(:,:)
      rate(:,:,703) = 6.300e-12_dp * exp_fac(:,:)
      rate(:,:,686) = 1.800e-11_dp * exp( -460._dp * itemp(:,:) )
      rate(:,:,687) = 2.420e-14_dp * exp( 1617._dp * itemp(:,:) )
      rate(:,:,689) = 1.13e-11_dp * exp( -253._dp * itemp(:,:) )
      rate(:,:,691) = 1.9e-13_dp * exp( 520._dp * itemp(:,:) )
      rate(:,:,693) = 9.e-11_dp * exp( -2386._dp * itemp(:,:) )
      rate(:,:,696) = 1.8e13_dp * exp( -8661._dp * itemp(:,:) )
      rate(:,:,700) = 1.400e-11_dp * exp( -1030._dp * itemp(:,:) )
      rate(:,:,702) = 1.350e-12_dp * exp( -600._dp * itemp(:,:) )
      rate(:,:,704) = 4.850e-12_dp * exp( -850._dp * itemp(:,:) )
      rate(:,:,705) = 2.170e-11_dp * exp( -1130._dp * itemp(:,:) )
      rate(:,:,706) = 2.400e-12_dp * exp( -1250._dp * itemp(:,:) )
      rate(:,:,707) = 1.640e-12_dp * exp( -1520._dp * itemp(:,:) )
      rate(:,:,710) = 1.300e-12_dp * exp( -1770._dp * itemp(:,:) )
      itemp(:,:) = 300._dp * itemp(:,:)
      ko(:,:) = 4.400e-32_dp * itemp(:,:)**1.3_dp
      kinf(:,:) = 7.500e-11_dp * itemp(:,:)**(-0.2_dp)
      call jpl( rate(1,1,150), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 6.900e-31_dp * itemp(:,:)**1._dp
      kinf(:,:) = 2.600e-11_dp
      call jpl( rate(1,1,159), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 9.000e-32_dp * itemp(:,:)**1.5_dp
      kinf(:,:) = 3.000e-11_dp
      call jpl( rate(1,1,171), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 2.500e-31_dp * itemp(:,:)**1.8_dp
      kinf(:,:) = 2.200e-11_dp * itemp(:,:)**0.7_dp
      call jpl( rate(1,1,175), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 1.800e-30_dp * itemp(:,:)**3._dp
      kinf(:,:) = 2.800e-11_dp
      call jpl( rate(1,1,178), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 7.000e-31_dp * itemp(:,:)**2.6_dp
      kinf(:,:) = 3.600e-11_dp * itemp(:,:)**0.1_dp
      call jpl( rate(1,1,184), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 2.000e-31_dp * itemp(:,:)**3.4_dp
      kinf(:,:) = 2.900e-12_dp * itemp(:,:)**1.1_dp
      call jpl( rate(1,1,186), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 2.000e-30_dp * itemp(:,:)**4.4_dp
      kinf(:,:) = 1.400e-12_dp * itemp(:,:)**0.7_dp
      call jpl( rate(1,1,189), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 5.900e-33_dp * itemp(:,:)**1.4_dp
      kinf(:,:) = 1.100e-12_dp * itemp(:,:)**(-1.3_dp)
      call jpl( rate(1,1,193), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 1.000e-30_dp * itemp(:,:)**4.8_dp
      kinf(:,:) = 7.200e-12_dp * itemp(:,:)**2.1_dp
      call jpl( rate(1,1,210), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 5.500e-30_dp
      kinf(:,:) = 8.300e-13_dp * itemp(:,:)**(-2._dp)
      call jpl( rate(1,1,212), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 1.000e-28_dp * itemp(:,:)**4.5_dp
      kinf(:,:) = 7.500e-12_dp * itemp(:,:)**0.85_dp
      call jpl( rate(1,1,213), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 2.700e-28_dp * itemp(:,:)**7.1_dp
      kinf(:,:) = 1.200e-11_dp * itemp(:,:)**0.9_dp
      call jpl( rate(1,1,227), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 2.700e-28_dp * itemp(:,:)**7.1_dp
      kinf(:,:) = 1.200e-11_dp * itemp(:,:)**0.9_dp
      call jpl( rate(1,1,244), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 8.000e-27_dp * itemp(:,:)**3.5_dp
      kinf(:,:) = 3.000e-11_dp
      call jpl( rate(1,1,260), m, 0.5_dp, ko, kinf, plnplv )
      ko(:,:) = 2.700e-28_dp * itemp(:,:)**7.1_dp
      kinf(:,:) = 1.200e-11_dp * itemp(:,:)**0.9_dp
      call jpl( rate(1,1,311), m, 0.3_dp, ko, kinf, plnplv )
      ko(:,:) = 2.700e-28_dp * itemp(:,:)**7.1_dp
      kinf(:,:) = 1.200e-11_dp * itemp(:,:)**0.9_dp
      call jpl( rate(1,1,344), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 2.700e-28_dp * itemp(:,:)**7.1_dp
      kinf(:,:) = 1.200e-11_dp * itemp(:,:)**0.9_dp
      call jpl( rate(1,1,351), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 2.700e-28_dp * itemp(:,:)**7.1_dp
      kinf(:,:) = 1.200e-11_dp * itemp(:,:)**0.9_dp
      call jpl( rate(1,1,439), m, 0.3_dp, ko, kinf, plnplv )
      ko(:,:) = 2.700e-28_dp * itemp(:,:)**7.1_dp
      kinf(:,:) = 1.200e-11_dp * itemp(:,:)**0.9_dp
      call jpl( rate(1,1,469), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 2.700e-28_dp * itemp(:,:)**7.1_dp
      kinf(:,:) = 1.200e-11_dp * itemp(:,:)**0.9_dp
      call jpl( rate(1,1,536), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 1.800e-31_dp * itemp(:,:)**3.4_dp
      kinf(:,:) = 1.500e-11_dp * itemp(:,:)**1.9_dp
      call jpl( rate(1,1,630), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 1.600e-32_dp * itemp(:,:)**4.5_dp
      kinf(:,:) = 3.000e-12_dp * itemp(:,:)**2._dp
      call jpl( rate(1,1,634), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 4.200e-31_dp * itemp(:,:)**2.4_dp
      kinf(:,:) = 2.700e-11_dp
      call jpl( rate(1,1,645), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 5.200e-31_dp * itemp(:,:)**3.2_dp
      kinf(:,:) = 6.900e-12_dp * itemp(:,:)**2.9_dp
      call jpl( rate(1,1,651), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 5.20e-30_dp * itemp(:,:)**2.4_dp
      kinf(:,:) = 2.2e-10_dp * itemp(:,:)**0.7_dp
      call jpl( rate(1,1,670), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 1.60e-29_dp * itemp(:,:)**3.3_dp
      kinf(:,:) = 3.1e-10_dp * itemp(:,:)
      call jpl( rate(1,1,671), m, 0.6_dp, ko, kinf, plnplv )
      ko(:,:) = 3.3e-31_dp * itemp(:,:)**4.3_dp
      kinf(:,:) = 1.6e-12_dp
      call jpl( rate(1,1,688), m, 0.6_dp, ko, kinf, plnplv )
      end subroutine setrxt
      subroutine adjrxt( rate, inv, m, plnplv )
      use mo_moz_mods, only : nfs, rxntot
      implicit none
!--------------------------------------------------------------------
! ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real(dp), intent(in) :: inv(plnplv,nfs)
      real(dp), intent(in) :: m(plnplv)
      real(dp), intent(inout) :: rate(plnplv,rxntot)
!--------------------------------------------------------------------
! ... Local variables
!--------------------------------------------------------------------
      real(dp) :: im(plnplv)
      rate(:,143) = rate(:,143) * inv(:, 1)
      rate(:,145) = rate(:,145) * inv(:, 1)
      rate(:,146) = rate(:,146) * inv(:, 2)
      rate(:,150) = rate(:,150) * inv(:, 1)
      rate(:,159) = rate(:,159) * inv(:, 1)
      rate(:,171) = rate(:,171) * inv(:, 1)
      rate(:,175) = rate(:,175) * inv(:, 1)
      rate(:,178) = rate(:,178) * inv(:, 1)
      rate(:,186) = rate(:,186) * inv(:, 1)
      rate(:,188) = rate(:,188) * inv(:, 1)
      rate(:,189) = rate(:,189) * inv(:, 1)
      rate(:,190) = rate(:,190) * inv(:, 1)
      rate(:,193) = rate(:,193) * inv(:, 1)
      rate(:,210) = rate(:,210) * inv(:, 1)
      rate(:,211) = rate(:,211) * inv(:, 1)
      rate(:,212) = rate(:,212) * inv(:, 1)
      rate(:,213) = rate(:,213) * inv(:, 1)
      rate(:,227) = rate(:,227) * inv(:, 1)
      rate(:,232) = rate(:,232) * inv(:, 1)
      rate(:,244) = rate(:,244) * inv(:, 1)
      rate(:,260) = rate(:,260) * inv(:, 1)
      rate(:,311) = rate(:,311) * inv(:, 1)
      rate(:,321) = rate(:,321) * inv(:, 1)
      rate(:,344) = rate(:,344) * inv(:, 1)
      rate(:,351) = rate(:,351) * inv(:, 1)
      rate(:,439) = rate(:,439) * inv(:, 1)
      rate(:,442) = rate(:,442) * inv(:, 1)
      rate(:,469) = rate(:,469) * inv(:, 1)
      rate(:,536) = rate(:,536) * inv(:, 1)
      rate(:,537) = rate(:,537) * inv(:, 1)
      rate(:,624) = rate(:,624) * inv(:, 1)
      rate(:,630) = rate(:,630) * inv(:, 1)
      rate(:,634) = rate(:,634) * inv(:, 1)
      rate(:,645) = rate(:,645) * inv(:, 1)
      rate(:,646) = rate(:,646) * inv(:, 1)
      rate(:,651) = rate(:,651) * inv(:, 1)
      rate(:,670) = rate(:,670) * inv(:, 1)
      rate(:,671) = rate(:,671) * inv(:, 1)
      rate(:,688) = rate(:,688) * inv(:, 1)
      rate(:,143) = rate(:,143) * m(:)
      rate(:,144) = rate(:,144) * m(:)
      rate(:,145) = rate(:,145) * m(:)
      rate(:,147) = rate(:,147) * m(:)
      rate(:,148) = rate(:,148) * m(:)
      rate(:,149) = rate(:,149) * m(:)
      rate(:,150) = rate(:,150) * m(:)
      rate(:,151) = rate(:,151) * m(:)
      rate(:,152) = rate(:,152) * m(:)
      rate(:,153) = rate(:,153) * m(:)
      rate(:,154) = rate(:,154) * m(:)
      rate(:,155) = rate(:,155) * m(:)
      rate(:,156) = rate(:,156) * m(:)
      rate(:,157) = rate(:,157) * m(:)
      rate(:,158) = rate(:,158) * m(:)
      rate(:,159) = rate(:,159) * m(:)
      rate(:,160) = rate(:,160) * m(:)
      rate(:,161) = rate(:,161) * m(:)
      rate(:,162) = rate(:,162) * m(:)
      rate(:,163) = rate(:,163) * m(:)
      rate(:,164) = rate(:,164) * m(:)
      rate(:,165) = rate(:,165) * m(:)
      rate(:,166) = rate(:,166) * m(:)
      rate(:,167) = rate(:,167) * m(:)
      rate(:,168) = rate(:,168) * m(:)
      rate(:,169) = rate(:,169) * m(:)
      rate(:,170) = rate(:,170) * m(:)
      rate(:,171) = rate(:,171) * m(:)
      rate(:,172) = rate(:,172) * m(:)
      rate(:,173) = rate(:,173) * m(:)
      rate(:,174) = rate(:,174) * m(:)
      rate(:,175) = rate(:,175) * m(:)
      rate(:,176) = rate(:,176) * m(:)
      rate(:,177) = rate(:,177) * m(:)
      rate(:,178) = rate(:,178) * m(:)
      rate(:,179) = rate(:,179) * m(:)
      rate(:,180) = rate(:,180) * m(:)
      rate(:,181) = rate(:,181) * m(:)
      rate(:,182) = rate(:,182) * m(:)
      rate(:,183) = rate(:,183) * m(:)
      rate(:,184) = rate(:,184) * m(:)
      rate(:,185) = rate(:,185) * m(:)
      rate(:,186) = rate(:,186) * m(:)
      rate(:,187) = rate(:,187) * m(:)
      rate(:,189) = rate(:,189) * m(:)
      rate(:,191) = rate(:,191) * m(:)
      rate(:,192) = rate(:,192) * m(:)
      rate(:,193) = rate(:,193) * m(:)
      rate(:,194) = rate(:,194) * m(:)
      rate(:,195) = rate(:,195) * m(:)
      rate(:,196) = rate(:,196) * m(:)
      rate(:,197) = rate(:,197) * m(:)
      rate(:,198) = rate(:,198) * m(:)
      rate(:,199) = rate(:,199) * m(:)
      rate(:,200) = rate(:,200) * m(:)
      rate(:,201) = rate(:,201) * m(:)
      rate(:,202) = rate(:,202) * m(:)
      rate(:,203) = rate(:,203) * m(:)
      rate(:,204) = rate(:,204) * m(:)
      rate(:,205) = rate(:,205) * m(:)
      rate(:,207) = rate(:,207) * m(:)
      rate(:,208) = rate(:,208) * m(:)
      rate(:,209) = rate(:,209) * m(:)
      rate(:,210) = rate(:,210) * m(:)
      rate(:,212) = rate(:,212) * m(:)
      rate(:,213) = rate(:,213) * m(:)
      rate(:,214) = rate(:,214) * m(:)
      rate(:,215) = rate(:,215) * m(:)
      rate(:,216) = rate(:,216) * m(:)
      rate(:,217) = rate(:,217) * m(:)
      rate(:,218) = rate(:,218) * m(:)
      rate(:,219) = rate(:,219) * m(:)
      rate(:,220) = rate(:,220) * m(:)
      rate(:,221) = rate(:,221) * m(:)
      rate(:,222) = rate(:,222) * m(:)
      rate(:,223) = rate(:,223) * m(:)
      rate(:,224) = rate(:,224) * m(:)
      rate(:,225) = rate(:,225) * m(:)
      rate(:,226) = rate(:,226) * m(:)
      rate(:,227) = rate(:,227) * m(:)
      rate(:,228) = rate(:,228) * m(:)
      rate(:,229) = rate(:,229) * m(:)
      rate(:,230) = rate(:,230) * m(:)
      rate(:,231) = rate(:,231) * m(:)
      rate(:,233) = rate(:,233) * m(:)
      rate(:,234) = rate(:,234) * m(:)
      rate(:,235) = rate(:,235) * m(:)
      rate(:,236) = rate(:,236) * m(:)
      rate(:,237) = rate(:,237) * m(:)
      rate(:,238) = rate(:,238) * m(:)
      rate(:,240) = rate(:,240) * m(:)
      rate(:,241) = rate(:,241) * m(:)
      rate(:,242) = rate(:,242) * m(:)
      rate(:,243) = rate(:,243) * m(:)
      rate(:,244) = rate(:,244) * m(:)
      rate(:,245) = rate(:,245) * m(:)
      rate(:,246) = rate(:,246) * m(:)
      rate(:,247) = rate(:,247) * m(:)
      rate(:,248) = rate(:,248) * m(:)
      rate(:,249) = rate(:,249) * m(:)
      rate(:,250) = rate(:,250) * m(:)
      rate(:,251) = rate(:,251) * m(:)
      rate(:,252) = rate(:,252) * m(:)
      rate(:,253) = rate(:,253) * m(:)
      rate(:,254) = rate(:,254) * m(:)
      rate(:,255) = rate(:,255) * m(:)
      rate(:,256) = rate(:,256) * m(:)
      rate(:,257) = rate(:,257) * m(:)
      rate(:,258) = rate(:,258) * m(:)
      rate(:,259) = rate(:,259) * m(:)
      rate(:,260) = rate(:,260) * m(:)
      rate(:,261) = rate(:,261) * m(:)
      rate(:,262) = rate(:,262) * m(:)
      rate(:,263) = rate(:,263) * m(:)
      rate(:,264) = rate(:,264) * m(:)
      rate(:,265) = rate(:,265) * m(:)
      rate(:,266) = rate(:,266) * m(:)
      rate(:,267) = rate(:,267) * m(:)
      rate(:,268) = rate(:,268) * m(:)
      rate(:,269) = rate(:,269) * m(:)
      rate(:,270) = rate(:,270) * m(:)
      rate(:,271) = rate(:,271) * m(:)
      rate(:,272) = rate(:,272) * m(:)
      rate(:,273) = rate(:,273) * m(:)
      rate(:,274) = rate(:,274) * m(:)
      rate(:,275) = rate(:,275) * m(:)
      rate(:,276) = rate(:,276) * m(:)
      rate(:,277) = rate(:,277) * m(:)
      rate(:,278) = rate(:,278) * m(:)
      rate(:,279) = rate(:,279) * m(:)
      rate(:,280) = rate(:,280) * m(:)
      rate(:,281) = rate(:,281) * m(:)
      rate(:,282) = rate(:,282) * m(:)
      rate(:,283) = rate(:,283) * m(:)
      rate(:,284) = rate(:,284) * m(:)
      rate(:,285) = rate(:,285) * m(:)
      rate(:,286) = rate(:,286) * m(:)
      rate(:,287) = rate(:,287) * m(:)
      rate(:,288) = rate(:,288) * m(:)
      rate(:,289) = rate(:,289) * m(:)
      rate(:,290) = rate(:,290) * m(:)
      rate(:,291) = rate(:,291) * m(:)
      rate(:,292) = rate(:,292) * m(:)
      rate(:,293) = rate(:,293) * m(:)
      rate(:,294) = rate(:,294) * m(:)
      rate(:,295) = rate(:,295) * m(:)
      rate(:,296) = rate(:,296) * m(:)
      rate(:,297) = rate(:,297) * m(:)
      rate(:,298) = rate(:,298) * m(:)
      rate(:,299) = rate(:,299) * m(:)
      rate(:,300) = rate(:,300) * m(:)
      rate(:,301) = rate(:,301) * m(:)
      rate(:,302) = rate(:,302) * m(:)
      rate(:,303) = rate(:,303) * m(:)
      rate(:,304) = rate(:,304) * m(:)
      rate(:,305) = rate(:,305) * m(:)
      rate(:,306) = rate(:,306) * m(:)
      rate(:,307) = rate(:,307) * m(:)
      rate(:,308) = rate(:,308) * m(:)
      rate(:,309) = rate(:,309) * m(:)
      rate(:,310) = rate(:,310) * m(:)
      rate(:,311) = rate(:,311) * m(:)
      rate(:,312) = rate(:,312) * m(:)
      rate(:,313) = rate(:,313) * m(:)
      rate(:,314) = rate(:,314) * m(:)
      rate(:,315) = rate(:,315) * m(:)
      rate(:,316) = rate(:,316) * m(:)
      rate(:,317) = rate(:,317) * m(:)
      rate(:,319) = rate(:,319) * m(:)
      rate(:,320) = rate(:,320) * m(:)
      rate(:,322) = rate(:,322) * m(:)
      rate(:,323) = rate(:,323) * m(:)
      rate(:,324) = rate(:,324) * m(:)
      rate(:,325) = rate(:,325) * m(:)
      rate(:,326) = rate(:,326) * m(:)
      rate(:,327) = rate(:,327) * m(:)
      rate(:,328) = rate(:,328) * m(:)
      rate(:,329) = rate(:,329) * m(:)
      rate(:,330) = rate(:,330) * m(:)
      rate(:,331) = rate(:,331) * m(:)
      rate(:,332) = rate(:,332) * m(:)
      rate(:,333) = rate(:,333) * m(:)
      rate(:,334) = rate(:,334) * m(:)
      rate(:,335) = rate(:,335) * m(:)
      rate(:,336) = rate(:,336) * m(:)
      rate(:,337) = rate(:,337) * m(:)
      rate(:,338) = rate(:,338) * m(:)
      rate(:,339) = rate(:,339) * m(:)
      rate(:,340) = rate(:,340) * m(:)
      rate(:,341) = rate(:,341) * m(:)
      rate(:,342) = rate(:,342) * m(:)
      rate(:,343) = rate(:,343) * m(:)
      rate(:,344) = rate(:,344) * m(:)
      rate(:,345) = rate(:,345) * m(:)
      rate(:,346) = rate(:,346) * m(:)
      rate(:,347) = rate(:,347) * m(:)
      rate(:,348) = rate(:,348) * m(:)
      rate(:,349) = rate(:,349) * m(:)
      rate(:,350) = rate(:,350) * m(:)
      rate(:,351) = rate(:,351) * m(:)
      rate(:,352) = rate(:,352) * m(:)
      rate(:,353) = rate(:,353) * m(:)
      rate(:,354) = rate(:,354) * m(:)
      rate(:,355) = rate(:,355) * m(:)
      rate(:,356) = rate(:,356) * m(:)
      rate(:,357) = rate(:,357) * m(:)
      rate(:,358) = rate(:,358) * m(:)
      rate(:,359) = rate(:,359) * m(:)
      rate(:,360) = rate(:,360) * m(:)
      rate(:,361) = rate(:,361) * m(:)
      rate(:,362) = rate(:,362) * m(:)
      rate(:,363) = rate(:,363) * m(:)
      rate(:,364) = rate(:,364) * m(:)
      rate(:,365) = rate(:,365) * m(:)
      rate(:,366) = rate(:,366) * m(:)
      rate(:,367) = rate(:,367) * m(:)
      rate(:,368) = rate(:,368) * m(:)
      rate(:,369) = rate(:,369) * m(:)
      rate(:,370) = rate(:,370) * m(:)
      rate(:,371) = rate(:,371) * m(:)
      rate(:,372) = rate(:,372) * m(:)
      rate(:,373) = rate(:,373) * m(:)
      rate(:,374) = rate(:,374) * m(:)
      rate(:,378) = rate(:,378) * m(:)
      rate(:,379) = rate(:,379) * m(:)
      rate(:,380) = rate(:,380) * m(:)
      rate(:,381) = rate(:,381) * m(:)
      rate(:,382) = rate(:,382) * m(:)
      rate(:,383) = rate(:,383) * m(:)
      rate(:,384) = rate(:,384) * m(:)
      rate(:,385) = rate(:,385) * m(:)
      rate(:,386) = rate(:,386) * m(:)
      rate(:,387) = rate(:,387) * m(:)
      rate(:,388) = rate(:,388) * m(:)
      rate(:,389) = rate(:,389) * m(:)
      rate(:,390) = rate(:,390) * m(:)
      rate(:,391) = rate(:,391) * m(:)
      rate(:,392) = rate(:,392) * m(:)
      rate(:,393) = rate(:,393) * m(:)
      rate(:,394) = rate(:,394) * m(:)
      rate(:,395) = rate(:,395) * m(:)
      rate(:,396) = rate(:,396) * m(:)
      rate(:,397) = rate(:,397) * m(:)
      rate(:,398) = rate(:,398) * m(:)
      rate(:,399) = rate(:,399) * m(:)
      rate(:,400) = rate(:,400) * m(:)
      rate(:,401) = rate(:,401) * m(:)
      rate(:,402) = rate(:,402) * m(:)
      rate(:,403) = rate(:,403) * m(:)
      rate(:,404) = rate(:,404) * m(:)
      rate(:,405) = rate(:,405) * m(:)
      rate(:,406) = rate(:,406) * m(:)
      rate(:,407) = rate(:,407) * m(:)
      rate(:,408) = rate(:,408) * m(:)
      rate(:,409) = rate(:,409) * m(:)
      rate(:,410) = rate(:,410) * m(:)
      rate(:,411) = rate(:,411) * m(:)
      rate(:,412) = rate(:,412) * m(:)
      rate(:,413) = rate(:,413) * m(:)
      rate(:,414) = rate(:,414) * m(:)
      rate(:,415) = rate(:,415) * m(:)
      rate(:,416) = rate(:,416) * m(:)
      rate(:,417) = rate(:,417) * m(:)
      rate(:,418) = rate(:,418) * m(:)
      rate(:,419) = rate(:,419) * m(:)
      rate(:,420) = rate(:,420) * m(:)
      rate(:,421) = rate(:,421) * m(:)
      rate(:,422) = rate(:,422) * m(:)
      rate(:,423) = rate(:,423) * m(:)
      rate(:,424) = rate(:,424) * m(:)
      rate(:,425) = rate(:,425) * m(:)
      rate(:,426) = rate(:,426) * m(:)
      rate(:,427) = rate(:,427) * m(:)
      rate(:,428) = rate(:,428) * m(:)
      rate(:,429) = rate(:,429) * m(:)
      rate(:,430) = rate(:,430) * m(:)
      rate(:,431) = rate(:,431) * m(:)
      rate(:,432) = rate(:,432) * m(:)
      rate(:,434) = rate(:,434) * m(:)
      rate(:,435) = rate(:,435) * m(:)
      rate(:,436) = rate(:,436) * m(:)
      rate(:,437) = rate(:,437) * m(:)
      rate(:,438) = rate(:,438) * m(:)
      rate(:,439) = rate(:,439) * m(:)
      rate(:,440) = rate(:,440) * m(:)
      rate(:,441) = rate(:,441) * m(:)
      rate(:,443) = rate(:,443) * m(:)
      rate(:,444) = rate(:,444) * m(:)
      rate(:,445) = rate(:,445) * m(:)
      rate(:,446) = rate(:,446) * m(:)
      rate(:,447) = rate(:,447) * m(:)
      rate(:,448) = rate(:,448) * m(:)
      rate(:,449) = rate(:,449) * m(:)
      rate(:,450) = rate(:,450) * m(:)
      rate(:,451) = rate(:,451) * m(:)
      rate(:,452) = rate(:,452) * m(:)
      rate(:,453) = rate(:,453) * m(:)
      rate(:,454) = rate(:,454) * m(:)
      rate(:,455) = rate(:,455) * m(:)
      rate(:,456) = rate(:,456) * m(:)
      rate(:,457) = rate(:,457) * m(:)
      rate(:,458) = rate(:,458) * m(:)
      rate(:,459) = rate(:,459) * m(:)
      rate(:,460) = rate(:,460) * m(:)
      rate(:,461) = rate(:,461) * m(:)
      rate(:,462) = rate(:,462) * m(:)
      rate(:,463) = rate(:,463) * m(:)
      rate(:,464) = rate(:,464) * m(:)
      rate(:,465) = rate(:,465) * m(:)
      rate(:,466) = rate(:,466) * m(:)
      rate(:,467) = rate(:,467) * m(:)
      rate(:,468) = rate(:,468) * m(:)
      rate(:,469) = rate(:,469) * m(:)
      rate(:,470) = rate(:,470) * m(:)
      rate(:,471) = rate(:,471) * m(:)
      rate(:,472) = rate(:,472) * m(:)
      rate(:,473) = rate(:,473) * m(:)
      rate(:,474) = rate(:,474) * m(:)
      rate(:,475) = rate(:,475) * m(:)
      rate(:,476) = rate(:,476) * m(:)
      rate(:,477) = rate(:,477) * m(:)
      rate(:,478) = rate(:,478) * m(:)
      rate(:,479) = rate(:,479) * m(:)
      rate(:,480) = rate(:,480) * m(:)
      rate(:,481) = rate(:,481) * m(:)
      rate(:,482) = rate(:,482) * m(:)
      rate(:,483) = rate(:,483) * m(:)
      rate(:,484) = rate(:,484) * m(:)
      rate(:,485) = rate(:,485) * m(:)
      rate(:,486) = rate(:,486) * m(:)
      rate(:,487) = rate(:,487) * m(:)
      rate(:,488) = rate(:,488) * m(:)
      rate(:,489) = rate(:,489) * m(:)
      rate(:,490) = rate(:,490) * m(:)
      rate(:,491) = rate(:,491) * m(:)
      rate(:,492) = rate(:,492) * m(:)
      rate(:,493) = rate(:,493) * m(:)
      rate(:,494) = rate(:,494) * m(:)
      rate(:,495) = rate(:,495) * m(:)
      rate(:,496) = rate(:,496) * m(:)
      rate(:,497) = rate(:,497) * m(:)
      rate(:,498) = rate(:,498) * m(:)
      rate(:,499) = rate(:,499) * m(:)
      rate(:,500) = rate(:,500) * m(:)
      rate(:,501) = rate(:,501) * m(:)
      rate(:,502) = rate(:,502) * m(:)
      rate(:,503) = rate(:,503) * m(:)
      rate(:,504) = rate(:,504) * m(:)
      rate(:,505) = rate(:,505) * m(:)
      rate(:,506) = rate(:,506) * m(:)
      rate(:,507) = rate(:,507) * m(:)
      rate(:,508) = rate(:,508) * m(:)
      rate(:,509) = rate(:,509) * m(:)
      rate(:,510) = rate(:,510) * m(:)
      rate(:,511) = rate(:,511) * m(:)
      rate(:,512) = rate(:,512) * m(:)
      rate(:,513) = rate(:,513) * m(:)
      rate(:,514) = rate(:,514) * m(:)
      rate(:,515) = rate(:,515) * m(:)
      rate(:,516) = rate(:,516) * m(:)
      rate(:,517) = rate(:,517) * m(:)
      rate(:,518) = rate(:,518) * m(:)
      rate(:,519) = rate(:,519) * m(:)
      rate(:,520) = rate(:,520) * m(:)
      rate(:,521) = rate(:,521) * m(:)
      rate(:,522) = rate(:,522) * m(:)
      rate(:,523) = rate(:,523) * m(:)
      rate(:,524) = rate(:,524) * m(:)
      rate(:,525) = rate(:,525) * m(:)
      rate(:,526) = rate(:,526) * m(:)
      rate(:,527) = rate(:,527) * m(:)
      rate(:,528) = rate(:,528) * m(:)
      rate(:,529) = rate(:,529) * m(:)
      rate(:,530) = rate(:,530) * m(:)
      rate(:,531) = rate(:,531) * m(:)
      rate(:,532) = rate(:,532) * m(:)
      rate(:,533) = rate(:,533) * m(:)
      rate(:,534) = rate(:,534) * m(:)
      rate(:,535) = rate(:,535) * m(:)
      rate(:,536) = rate(:,536) * m(:)
      rate(:,538) = rate(:,538) * m(:)
      rate(:,539) = rate(:,539) * m(:)
      rate(:,540) = rate(:,540) * m(:)
      rate(:,541) = rate(:,541) * m(:)
      rate(:,542) = rate(:,542) * m(:)
      rate(:,543) = rate(:,543) * m(:)
      rate(:,544) = rate(:,544) * m(:)
      rate(:,545) = rate(:,545) * m(:)
      rate(:,546) = rate(:,546) * m(:)
      rate(:,547) = rate(:,547) * m(:)
      rate(:,548) = rate(:,548) * m(:)
      rate(:,549) = rate(:,549) * m(:)
      rate(:,550) = rate(:,550) * m(:)
      rate(:,551) = rate(:,551) * m(:)
      rate(:,552) = rate(:,552) * m(:)
      rate(:,553) = rate(:,553) * m(:)
      rate(:,554) = rate(:,554) * m(:)
      rate(:,555) = rate(:,555) * m(:)
      rate(:,556) = rate(:,556) * m(:)
      rate(:,557) = rate(:,557) * m(:)
      rate(:,558) = rate(:,558) * m(:)
      rate(:,559) = rate(:,559) * m(:)
      rate(:,560) = rate(:,560) * m(:)
      rate(:,561) = rate(:,561) * m(:)
      rate(:,562) = rate(:,562) * m(:)
      rate(:,563) = rate(:,563) * m(:)
      rate(:,564) = rate(:,564) * m(:)
      rate(:,565) = rate(:,565) * m(:)
      rate(:,566) = rate(:,566) * m(:)
      rate(:,567) = rate(:,567) * m(:)
      rate(:,568) = rate(:,568) * m(:)
      rate(:,569) = rate(:,569) * m(:)
      rate(:,570) = rate(:,570) * m(:)
      rate(:,571) = rate(:,571) * m(:)
      rate(:,572) = rate(:,572) * m(:)
      rate(:,573) = rate(:,573) * m(:)
      rate(:,574) = rate(:,574) * m(:)
      rate(:,575) = rate(:,575) * m(:)
      rate(:,576) = rate(:,576) * m(:)
      rate(:,577) = rate(:,577) * m(:)
      rate(:,578) = rate(:,578) * m(:)
      rate(:,579) = rate(:,579) * m(:)
      rate(:,580) = rate(:,580) * m(:)
      rate(:,581) = rate(:,581) * m(:)
      rate(:,582) = rate(:,582) * m(:)
      rate(:,583) = rate(:,583) * m(:)
      rate(:,584) = rate(:,584) * m(:)
      rate(:,585) = rate(:,585) * m(:)
      rate(:,586) = rate(:,586) * m(:)
      rate(:,587) = rate(:,587) * m(:)
      rate(:,588) = rate(:,588) * m(:)
      rate(:,589) = rate(:,589) * m(:)
      rate(:,590) = rate(:,590) * m(:)
      rate(:,591) = rate(:,591) * m(:)
      rate(:,592) = rate(:,592) * m(:)
      rate(:,593) = rate(:,593) * m(:)
      rate(:,594) = rate(:,594) * m(:)
      rate(:,595) = rate(:,595) * m(:)
      rate(:,596) = rate(:,596) * m(:)
      rate(:,597) = rate(:,597) * m(:)
      rate(:,598) = rate(:,598) * m(:)
      rate(:,599) = rate(:,599) * m(:)
      rate(:,600) = rate(:,600) * m(:)
      rate(:,601) = rate(:,601) * m(:)
      rate(:,602) = rate(:,602) * m(:)
      rate(:,603) = rate(:,603) * m(:)
      rate(:,604) = rate(:,604) * m(:)
      rate(:,605) = rate(:,605) * m(:)
      rate(:,606) = rate(:,606) * m(:)
      rate(:,607) = rate(:,607) * m(:)
      rate(:,608) = rate(:,608) * m(:)
      rate(:,609) = rate(:,609) * m(:)
      rate(:,610) = rate(:,610) * m(:)
      rate(:,611) = rate(:,611) * m(:)
      rate(:,612) = rate(:,612) * m(:)
      rate(:,613) = rate(:,613) * m(:)
      rate(:,614) = rate(:,614) * m(:)
      rate(:,615) = rate(:,615) * m(:)
      rate(:,616) = rate(:,616) * m(:)
      rate(:,617) = rate(:,617) * m(:)
      rate(:,618) = rate(:,618) * m(:)
      rate(:,619) = rate(:,619) * m(:)
      rate(:,620) = rate(:,620) * m(:)
      rate(:,621) = rate(:,621) * m(:)
      rate(:,622) = rate(:,622) * m(:)
      rate(:,623) = rate(:,623) * m(:)
      rate(:,625) = rate(:,625) * m(:)
      rate(:,626) = rate(:,626) * m(:)
      rate(:,627) = rate(:,627) * m(:)
      rate(:,628) = rate(:,628) * m(:)
      rate(:,629) = rate(:,629) * m(:)
      rate(:,630) = rate(:,630) * m(:)
      rate(:,631) = rate(:,631) * m(:)
      rate(:,632) = rate(:,632) * m(:)
      rate(:,633) = rate(:,633) * m(:)
      rate(:,634) = rate(:,634) * m(:)
      rate(:,635) = rate(:,635) * m(:)
      rate(:,636) = rate(:,636) * m(:)
      rate(:,637) = rate(:,637) * m(:)
      rate(:,638) = rate(:,638) * m(:)
      rate(:,639) = rate(:,639) * m(:)
      rate(:,640) = rate(:,640) * m(:)
      rate(:,641) = rate(:,641) * m(:)
      rate(:,642) = rate(:,642) * m(:)
      rate(:,643) = rate(:,643) * m(:)
      rate(:,644) = rate(:,644) * m(:)
      rate(:,645) = rate(:,645) * m(:)
      rate(:,647) = rate(:,647) * m(:)
      rate(:,648) = rate(:,648) * m(:)
      rate(:,649) = rate(:,649) * m(:)
      rate(:,650) = rate(:,650) * m(:)
      rate(:,651) = rate(:,651) * m(:)
      rate(:,652) = rate(:,652) * m(:)
      rate(:,653) = rate(:,653) * m(:)
      rate(:,654) = rate(:,654) * m(:)
      rate(:,655) = rate(:,655) * m(:)
      rate(:,656) = rate(:,656) * m(:)
      rate(:,657) = rate(:,657) * m(:)
      rate(:,658) = rate(:,658) * m(:)
      rate(:,659) = rate(:,659) * m(:)
      rate(:,660) = rate(:,660) * m(:)
      rate(:,661) = rate(:,661) * m(:)
      rate(:,662) = rate(:,662) * m(:)
      rate(:,663) = rate(:,663) * m(:)
      rate(:,664) = rate(:,664) * m(:)
      rate(:,665) = rate(:,665) * m(:)
      rate(:,666) = rate(:,666) * m(:)
      rate(:,667) = rate(:,667) * m(:)
      rate(:,668) = rate(:,668) * m(:)
      rate(:,669) = rate(:,669) * m(:)
      rate(:,670) = rate(:,670) * m(:)
      rate(:,671) = rate(:,671) * m(:)
      rate(:,672) = rate(:,672) * m(:)
      rate(:,673) = rate(:,673) * m(:)
      rate(:,674) = rate(:,674) * m(:)
      rate(:,675) = rate(:,675) * m(:)
      rate(:,676) = rate(:,676) * m(:)
      rate(:,677) = rate(:,677) * m(:)
      rate(:,678) = rate(:,678) * m(:)
      rate(:,679) = rate(:,679) * m(:)
      rate(:,680) = rate(:,680) * m(:)
      rate(:,681) = rate(:,681) * m(:)
      rate(:,682) = rate(:,682) * m(:)
      rate(:,683) = rate(:,683) * m(:)
      rate(:,684) = rate(:,684) * m(:)
      rate(:,685) = rate(:,685) * m(:)
      rate(:,686) = rate(:,686) * m(:)
      rate(:,687) = rate(:,687) * m(:)
      rate(:,688) = rate(:,688) * m(:)
      rate(:,689) = rate(:,689) * m(:)
      rate(:,690) = rate(:,690) * m(:)
      rate(:,691) = rate(:,691) * m(:)
      rate(:,692) = rate(:,692) * m(:)
      rate(:,693) = rate(:,693) * m(:)
      rate(:,694) = rate(:,694) * m(:)
      rate(:,695) = rate(:,695) * m(:)
      rate(:,697) = rate(:,697) * m(:)
      rate(:,698) = rate(:,698) * m(:)
      rate(:,699) = rate(:,699) * m(:)
      rate(:,700) = rate(:,700) * m(:)
      rate(:,701) = rate(:,701) * m(:)
      rate(:,702) = rate(:,702) * m(:)
      rate(:,703) = rate(:,703) * m(:)
      rate(:,704) = rate(:,704) * m(:)
      rate(:,705) = rate(:,705) * m(:)
      rate(:,706) = rate(:,706) * m(:)
      rate(:,707) = rate(:,707) * m(:)
      rate(:,708) = rate(:,708) * m(:)
      rate(:,709) = rate(:,709) * m(:)
      rate(:,710) = rate(:,710) * m(:)
      rate(:,720) = rate(:,720) * m(:)
      rate(:,721) = rate(:,721) * m(:)
      rate(:,722) = rate(:,722) * m(:)
      rate(:,725) = rate(:,725) * m(:)
      rate(:,726) = rate(:,726) * m(:)
      rate(:,731) = rate(:,731) * m(:)
      rate(:,732) = rate(:,732) * m(:)
      rate(:,733) = rate(:,733) * m(:)
      end subroutine adjrxt
      subroutine phtadj( p_rate, inv, m, plnplv )
      use mo_moz_mods, only : nfs, phtcnt
      use mo_kind, only : dp
      implicit none
!--------------------------------------------------------------------
! ... Dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real(dp), intent(in) :: inv(plnplv,nfs)
      real(dp), intent(in) :: m(plnplv)
      real(dp), intent(inout) :: p_rate(plnplv,phtcnt)
!--------------------------------------------------------------------
! ... Local variables
!--------------------------------------------------------------------
      real(dp) :: im(plnplv)
      end subroutine phtadj
      subroutine rxt_mod( rate, het_rates, grp_ratios, plnplv )
      use mo_moz_mods, only : rxntot, hetcnt, grpcnt
      implicit none
!---------------------------------------------------------------------------
! ... Dummy arguments
!---------------------------------------------------------------------------
      integer, intent(in) :: plnplv
      real(dp), intent(inout) :: rate(plnplv,rxntot)
      real(dp), intent(inout) :: het_rates(plnplv,hetcnt)
      real(dp), intent(in) :: grp_ratios(plnplv,grpcnt)
      end subroutine rxt_mod
      subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, plonl )
      use mo_moz_mods, only : plev, pcnstm1
      use mo_moz_mods, only : grpcnt
      implicit none
!----------------------------------------------------------------------------
! ... Dummy arguments
!----------------------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in) :: vmr(plonl,plev,pcnstm1)
      real, intent(in) :: group_ratios(plonl,plev,grpcnt)
      real, intent(out) :: group_vmrs(plonl,plev,grpcnt)
!----------------------------------------------------------------------------
! ... Local variables
!----------------------------------------------------------------------------
      integer :: k
      end subroutine mak_grp_vmr
      subroutine set_sim_dat
      use mo_moz_mods, only : explicit, implicit, rodas
      use mo_moz_mods, only : pcnstm1, grpcnt, ngrp, grp_mem_cnt, nfs
      use mo_moz_mods, only : hetcnt, drydep_cnt, srfems_cnt, extcnt, rxt_tag_cnt, fbc_cnt
      use mo_moz_mods, only : nadv_mass, adv_mass
      use mo_moz_mods, only : drydep_lst, srfems_lst, het_lst, extfrc_lst, grp_lst, inv_lst
      use mo_moz_mods, only : flbc_lst, fubc_lst
      use mo_moz_mods, only : rxt_tag_map, rxt_tag_lst
      use mo_moz_mods, only : phtcnt, pht_alias_lst, pht_alias_mult
      use mo_moz_mods, only : inv_from_dataset
      use mo_moz_mods, only : frc_from_dataset
      use mo_moz, only : tracnam, natsnam
      use mo_exception, only : finish, message_text
      implicit none
!--------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------
      integer :: ios
      tracnam(:242) = (/ 'ACBZO2                  ','ALKNO3                  ','ALKO2                   ', &
                         'ALKOH                   ','ALKOOH                  ','APIN                    ', &
                         'BCARY                   ','BENZ                    ','BENZO2                  ', &
                         'BENZOOH                 ','BEPOMUC                 ','BIACETOH                ', &
                         'BIGALD1                 ','BIGALD2                 ','BIGALD3                 ', &
                         'BIGALD4                 ','BIGALKANE               ','BIGENE                  ', &
                         'BPIN                    ','BR                      ','BRCL                    ', &
                         'BRO                     ','BRONO2                  ','BZALD                   ', &
                         'BZOO                    ','BZOOH                   ','C2H2                    ', &
                         'C2H4                    ','C2H5O2                  ','C2H5OH                  ', &
                         'C2H5OOH                 ','C2H6                    ','C3H6                    ', &
                         'C3H7O2                  ','C3H7OOH                 ','C3H8                    ', &
                         'C59O2                   ','C59OOH                  ','C5H8                    ', &
                         'C6H5O                   ','C6H5O2                  ','C6H5OOH                 ', &
                         'CCL4                    ','CF2CLBR                 ','CF3BR                   ', &
                         'CFC11                   ','CFC113                  ','CFC114                  ', &
                         'CFC115                  ','CFC12                   ','CH2BR2                  ', &
                         'CH2O                    ','CH3BR                   ','CH3CCL3                 ', &
                         'CH3CHO                  ','CH3CL                   ','CH3CN                   ', &
                         'CH3CO3                  ','CH3COCH2O2              ','CH3COCH3                ', &
                         'CH3COCHO                ','CH3COOH                 ','CH3COOOH                ', &
                         'CH3O2                   ','CH3OH                   ','CH3OOH                  ', &
                         'CH4                     ','CHBR3                   ','CL                      ', &
                         'CL2                     ','CL2O2                   ','CLO                     ', &
                         'CLONO2                  ','CO                      ','CO2                     ', &
                         'CO2H3CHO                ','CO2H3CO3                ','CO2H3CO3H               ', &
                         'COF2                    ','COFCL                   ','CRESOL                  ', &
                         'DICARBO2                ','ENEO2                   ','EO                      ', &
                         'EO2                     ','EOOH                    ','F                       ', &
                         'GLYALD                  ','GLYOXAL                 ','H                       ', &
                         'H1202                   ','H2                      ','H2402                   ', &
                         'H2O                     ','H2O2                    ','HBR                     ', &
                         'HCFC141B                ','HCFC142B                ','HCFC22                  ', &
                         'HCL                     ','HCN                     ','HCOC5                   ', &
                         'HCOCO2H                 ','HCOCO3                  ','HCOCO3H                 ', &
                         'HCOOH                   ','HF                      ','HNO3                    ', &
                         'HO2                     ','HO2NO2                  ','HOBR                    ', &
                         'HOCH2CO2H               ','HOCH2CO3                ','HOCH2CO3H               ', &
                         'HOCH2OO                 ','HOCL                    ','HONO                    ', &
                         'HPALD                   ','HYAC                    ','IBUTALOH                ', &
                         'IBUTALOHO2              ','IBUTALOHOOH             ','IEC1O2                  ', &
                         'ISOPAOH                 ','ISOPBNO3                ','ISOPBO2                 ', &
                         'ISOPBOH                 ','ISOPBOOH                ','ISOPDNO3                ', &
                         'ISOPDO2                 ','ISOPDOH                 ','ISOPDOOH                ', &
                         'LC578O2                 ','LC578OOH                ','LC5PAN1719              ', &
                         'LHC4ACCHO               ','LHC4ACCO2H              ','LHC4ACCO3               ', &
                         'LHC4ACCO3H              ','LHMVKABO2               ','LHMVKABOOH              ', &
                         'LIECHO                  ','LIECO3                  ','LIECO3H                 ', &
                         'LIEPOX                  ','LIMON                   ','LISOPACNO3              ', &
                         'LISOPACO2               ','LISOPACOOH              ','LNISO3                  ', &
                         'LNISOOH                 ','MACO2H                  ','MACO3H                  ', &
                         'MACR                    ','MACRO2                  ','MACROH                  ', &
                         'MACROOH                 ','MALO2                   ','MBO                     ', &
                         'MBONO3O2                ','MBOO2                   ','MBOOOH                  ', &
                         'MCO3                    ','MDIALO2                 ','MEK                     ', &
                         'MEKO2                   ','MEKOOH                  ','MPAN                    ', &
                         'MVK                     ','MYRC                    ','N                       ', &
                         'N2O                     ','N2O5                    ','NC4CHO                  ', &
                         'NISOPO2                 ','NISOPOOH                ','NO                      ', &
                         'NO2                     ','NO3                     ','NOA                     ', &
                         'NTERPO2                 ','O                       ','O1D                     ', &
                         'O2                      ','O3                      ','OCLO                    ', &
                         'OH                      ','PAN                     ','PBZNIT                  ', &
                         'PHENO2                  ','PHENOL                  ','PHENOOH                 ', &
                         'PO2                     ','POOH                    ','PR2O2HNO3               ', &
                         'PRONO3BO2               ','ROOH                    ','SF6                     ', &
                         'TEPOMUC                 ','TERP2O2                 ','MEKNO3                  ', &
                         'TERP2OOH                ','TERPNO3                 ','TERPO2                  ', &
                         'TERPOOH                 ','TERPROD1                ','TERPROD2                ', &
                         'TOL                     ','TOLO2                   ','TOLOOH                  ', &
                         'XYL                     ','XYLENO2                 ','XYLENOOH                ', &
                         'XYLOL                   ','XYLOLO2                 ','XYLOLOOH                ', &
                         'BRONO                   ','BRNO2                   ','NTERPNO3                ', &
                         'BR2                     ','CH3O2NO2                ','ELVOC                   ', &
                         'PACALD                  ','LISOPOOHO2              ','LISOPNO3O2              ', &
                         'LISOPNO3NO3             ','LISOPNO3OOH             ','LISOPOOHOOH             ', &
                         'CATECHOL                ','CATEC1OOH               ','CATEC1O2                ', &
                         'CATEC1O                 ','MACRN                   ','MVKN                    ', &
                         'NH3                     ','SO2                     ','H2SO4                   ', &
                         'DMS                     ','DMSO                    ','CH3SO2                  ', &
                         'CH3SO3                  ','CH3SO3H                 ' /)
      adv_mass(:242) = (/ 137.112198_dp , 133.141342_dp , 103.135201_dp , 88.1432037_dp , 104.142601_dp , &
                          136.228394_dp , 204.342590_dp , 78.1103973_dp , 159.114807_dp , 160.122192_dp , &
                          126.108597_dp , 102.086594_dp , 84.0724030_dp , 98.0982056_dp , 98.0982056_dp , &
                          112.123993_dp , 72.1437988_dp , 56.1031990_dp , 136.228394_dp , 79.9039993_dp , &
                          115.356705_dp , 95.9033966_dp , 141.908936_dp , 106.120796_dp , 123.127594_dp , &
                          124.134995_dp , 26.0368004_dp , 28.0515995_dp , 61.0578003_dp , 46.0657997_dp , &
                          62.0652008_dp , 30.0663986_dp , 42.0773964_dp , 75.0836029_dp , 76.0909958_dp , &
                          44.0921974_dp , 149.118591_dp , 150.126007_dp , 68.1141968_dp , 93.1023941_dp , &
                          109.101791_dp , 110.109192_dp , 153.821808_dp , 165.364517_dp , 148.910202_dp , &
                          137.367508_dp , 187.375320_dp , 170.921021_dp , 154.466721_dp , 120.913208_dp , &
                          173.833801_dp , 30.0251999_dp , 94.9372025_dp , 133.402313_dp , 44.0509987_dp , &
                          50.4859009_dp , 41.0509415_dp , 75.0424042_dp , 89.0681915_dp , 58.0767975_dp , &
                          72.0614014_dp , 60.0503998_dp , 76.0498047_dp , 47.0319977_dp , 32.0400009_dp , &
                          48.0393982_dp , 16.0405998_dp , 252.730408_dp , 35.4527016_dp , 70.9054031_dp , &
                          102.904205_dp , 51.4521027_dp , 97.4576416_dp , 28.0103989_dp , 44.0098000_dp , &
                          102.086594_dp , 133.078003_dp , 134.085403_dp , 66.0072098_dp , 82.4615097_dp , &
                          108.135597_dp , 129.089600_dp , 105.108795_dp , 61.0578003_dp , 77.0571976_dp , &
                          78.0645981_dp , 18.9984035_dp , 60.0503998_dp , 58.0355988_dp , 1.00740004_dp , &
                          209.815811_dp , 2.01480007_dp , 259.823608_dp , 18.0142002_dp , 34.0135994_dp , &
                          80.9113998_dp , 116.948013_dp , 100.493713_dp , 86.4679108_dp , 36.4601021_dp , &
                          27.0251389_dp , 100.112999_dp , 74.0350037_dp , 89.0269928_dp , 90.0343933_dp , &
                          46.0245972_dp , 20.0058041_dp , 63.0123405_dp , 33.0061989_dp , 79.0117416_dp , &
                          96.9107971_dp , 76.0498047_dp , 91.0417938_dp , 92.0491943_dp , 63.0314026_dp , &
                          52.4595032_dp , 47.0129395_dp , 116.112396_dp , 74.0762024_dp , 88.1019974_dp , &
                          119.093399_dp , 120.100800_dp , 149.118591_dp , 102.127800_dp , 147.125946_dp , &
                          117.119797_dp , 102.127800_dp , 118.127197_dp , 147.125946_dp , 117.119797_dp , &
                          102.127800_dp , 118.127197_dp , 149.118591_dp , 150.126007_dp , 177.109940_dp , &
                          100.112999_dp , 116.112396_dp , 131.104401_dp , 132.111801_dp , 119.093399_dp , &
                          120.100800_dp , 116.112396_dp , 147.103806_dp , 148.111206_dp , 118.127197_dp , &
                          136.228394_dp , 147.125946_dp , 117.119797_dp , 118.127197_dp , 176.102539_dp , &
                          177.109940_dp , 86.0872040_dp , 102.086594_dp , 70.0877991_dp , 119.093399_dp , &
                          104.101395_dp , 120.100800_dp , 115.063797_dp , 86.1284027_dp , 180.132141_dp , &
                          135.134003_dp , 136.141403_dp , 101.079193_dp , 129.089600_dp , 72.1026001_dp , &
                          103.093994_dp , 104.101395_dp , 147.084747_dp , 70.0877991_dp , 136.228394_dp , &
                          14.0067396_dp , 44.0128784_dp , 108.010483_dp , 145.111145_dp , 162.117950_dp , &
                          163.125336_dp , 30.0061398_dp , 46.0055389_dp , 62.0049400_dp , 119.074341_dp , &
                          230.232147_dp , 15.9994001_dp , 15.9994001_dp , 31.9988003_dp , 47.9981995_dp , &
                          67.4514999_dp , 17.0067997_dp , 121.047943_dp , 183.117737_dp , 175.114197_dp , &
                          94.1097946_dp , 176.121597_dp , 91.0829926_dp , 92.0903931_dp , 137.088531_dp , &
                          136.081146_dp , 90.0755920_dp , 146.056427_dp , 140.134399_dp , 199.218597_dp , &
                          149.099548_dp , 200.225998_dp , 215.240143_dp , 185.234009_dp , 186.241394_dp , &
                          168.227203_dp , 126.149796_dp , 92.1362000_dp , 173.140594_dp , 174.147995_dp , &
                          106.161995_dp , 187.166397_dp , 188.173798_dp , 122.161400_dp , 203.165802_dp , &
                          204.173187_dp , 125.909538_dp , 125.909538_dp , 231.239532_dp , 159.807999_dp , &
                          93.0375443_dp , 264.223602_dp , 130.097000_dp , 167.132797_dp , 196.131531_dp , &
                          226.137695_dp , 197.138947_dp , 168.140198_dp , 110.109192_dp , 126.108597_dp , &
                          125.101196_dp , 109.101791_dp , 149.099548_dp , 149.099548_dp , 17.0289402_dp , &
                          64.0648041_dp , 98.0783997_dp , 62.1324005_dp , 78.1318054_dp , 79.0980072_dp , &
                          95.0974045_dp , 96.1048050_dp /)
      implicit%cls_rxt_cnt(:) = (/ 0, 175, 558, 0 /)
      implicit%clsmap(:242) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
                                   11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                                   21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
                                   31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
                                   41, 42, 43, 44, 45, 46, 47, 48, 49, 50, &
                                   51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
                                   61, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
                                   71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
                                   81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
                                   91, 92, 93, 94, 95, 96, 97, 98, 99, 100, &
                                  101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
                                  111, 112, 113, 114, 115, 116, 117, 118, 119, 120, &
                                  121, 122, 123, 124, 125, 126, 127, 128, 129, 130, &
                                  131, 132, 133, 134, 135, 136, 137, 138, 139, 140, &
                                  141, 142, 143, 144, 145, 146, 147, 148, 149, 150, &
                                  151, 152, 153, 154, 155, 156, 157, 158, 159, 160, &
                                  161, 162, 163, 164, 165, 166, 167, 168, 169, 170, &
                                  171, 172, 173, 174, 175, 176, 177, 178, 179, 180, &
                                  181, 182, 183, 184, 185, 186, 187, 188, 189, 190, &
                                  191, 192, 193, 194, 195, 196, 197, 198, 199, 200, &
                                  201, 202, 203, 204, 205, 206, 207, 208, 209, 210, &
                                  211, 212, 213, 214, 215, 216, 217, 218, 219, 220, &
                                  221, 222, 223, 224, 225, 226, 227, 228, 229, 230, &
                                  231, 232, 233, 234, 235, 236, 237, 238, 239, 240, &
                                  241, 242 /)
      implicit%permute(:242) = (/ 146, 99, 199, 61, 111, 149, 152, 23, 141, 28, &
                                    24, 75, 76, 62, 108, 79, 134, 5, 153, 232, &
                                    43, 225, 130, 50, 137, 26, 65, 125, 176, 58, &
                                    70, 55, 156, 167, 90, 56, 195, 107, 155, 105, &
                                   174, 27, 7, 15, 16, 8, 17, 9, 18, 10, &
                                    82, 228, 100, 11, 202, 68, 77, 224, 200, 182, &
                                   209, 113, 81, 229, 192, 121, 207, 74, 235, 25, &
                                     6, 212, 145, 201, 181, 198, 171, 89, 39, 41, &
                                    63, 165, 170, 101, 172, 19, 138, 219, 211, 230, &
                                    12, 188, 13, 226, 157, 139, 40, 42, 51, 237, &
                                    52, 69, 123, 222, 126, 151, 53, 227, 234, 83, &
                                   148, 115, 223, 84, 85, 135, 57, 132, 208, 144, &
                                   179, 86, 177, 29, 95, 220, 30, 116, 87, 214, &
                                    31, 117, 205, 128, 92, 210, 129, 221, 118, 215, &
                                   169, 162, 143, 102, 93, 154, 114, 216, 103, 186, &
                                   104, 124, 88, 206, 203, 180, 133, 173, 131, 194, &
                                   191, 119, 193, 183, 140, 187, 71, 72, 213, 150, &
                                   120, 46, 54, 204, 185, 59, 242, 239, 233, 161, &
                                   184, 236, 240, 241, 231, 20, 238, 96, 21, 158, &
                                    66, 47, 197, 122, 97, 189, 73, 1, 34, 190, &
                                    32, 127, 78, 178, 98, 196, 168, 33, 163, 106, &
                                    38, 164, 35, 67, 159, 36, 44, 22, 94, 60, &
                                    45, 80, 91, 217, 218, 109, 136, 142, 64, 48, &
                                   175, 110, 147, 160, 4, 37, 2, 166, 49, 112, &
                                    14, 3 /)
      implicit%diag_map(:242) = (/ 1, 2, 3, 4, 6, 9, 12, 16, 20, 24, &
                                     28, 32, 36, 41, 43, 48, 53, 58, 63, 66, &
                                     69, 72, 75, 81, 85, 87, 90, 93, 96, 99, &
                                    102, 105, 109, 116, 120, 123, 127, 131, 139, 142, &
                                    147, 152, 157, 160, 165, 170, 174, 178, 184, 188, &
                                    192, 198, 204, 207, 213, 219, 225, 230, 235, 240, &
                                    244, 249, 252, 260, 265, 274, 282, 290, 298, 304, &
                                    310, 316, 322, 328, 334, 337, 340, 349, 354, 359, &
                                    364, 371, 378, 385, 392, 399, 406, 413, 420, 427, &
                                    434, 441, 447, 453, 458, 466, 474, 482, 490, 498, &
                                    507, 512, 520, 529, 538, 544, 553, 560, 564, 573, &
                                    577, 588, 593, 598, 604, 611, 622, 632, 642, 653, &
                                    661, 669, 677, 683, 689, 700, 707, 718, 727, 737, &
                                    746, 761, 772, 778, 785, 793, 803, 814, 823, 830, &
                                    840, 852, 867, 877, 886, 897, 907, 916, 926, 947, &
                                    967, 974, 995,1017,1038,1060,1075,1083,1095,1106, &
                                   1115,1123,1139,1159,1174,1189,1205,1218,1228,1239, &
                                   1256,1270,1283,1299,1312,1325,1341,1356,1370,1381, &
                                   1388,1391,1400,1414,1431,1452,1472,1490,1505,1525, &
                                   1548,1565,1576,1598,1618,1638,1657,1676,1699,1720, &
                                   1737,1744,1763,1784,1811,1836,1860,1873,1883,1900, &
                                   1920,1937,1965,2002,2035,2068,2104,2132,2151,2184, &
                                   2217,2239,2257,2348,2375,2397,2414,2435,2546,2566, &
                                   2648,2678,2775,2921,2980,3009,3034,3247,3302,3346, &
                                   3370,3482 /)
      rxt_tag_cnt = 733
      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(message_text,*) 'failed to allocate rxt_tag_lst; error = ',ios
         CALL finish( 'set_sim_dat',message_text )
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(message_text,*) 'failed to allocate rxt_tag_map; error = ',ios
         CALL finish( 'set_sim_dat',message_text )
      end if
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo3_a                           ', 'jo3_b                           ', &
                                     'jh2o2                           ', 'jn2o                            ', &
                                     'jno                             ', 'jno2                            ', &
                                     'jno3_a                          ', 'jno3_b                          ', &
                                     'jhno3                           ', 'jhono                           ', &
                                     'jho2no2_a                       ', 'jho2no2_b                       ', &
                                     'jn2o5_a                         ', 'jn2o5_b                         ', &
                                     'jco2                            ', 'jch4_a                          ', &
                                     'jch4_b                          ', 'jch2o_a                         ', &
                                     'jch2o_b                         ', 'jch3ooh                         ', &
                                     'jch3o2no2_a                     ', 'jch3o2no2_b                     ', &
                                     'jch3cho                         ', 'jch3coooh                       ', &
                                     'jc2h5ooh                        ', 'jpan                            ', &
                                     'jeooh                           ', 'jglyoxal                        ', &
                                     'jglyald                         ', 'jhoch2co3h                      ', &
                                     'jhcoco2h                        ', 'jhcoco3h                        ', &
                                     'jch3coch3                       ', 'jc3h7ooh                        ', &
                                     'jpooh                           ', 'jhyac                           ', &
                                     'jch3cocho                       ', 'jrooh                           ', &
                                     'jpr2o2hno3                      ', 'jnoa                            ', &
                                     'jmek                            ', 'jmekooh                         ', &
                                     'jmekno3                         ', 'jmacr_a                         ', &
                                     'jmacr_b                         ', 'jmacrooh                        ', &
                                     'jmacroh                         ', 'jmpan                           ', &
                                     'jmaco3h                         ', 'jmvk                            ', &
                                     'jlhmvkabooh                     ', 'jmvkn                           ', &
                                     'jmacrn                          ', 'jco2h3cho                       ', &
                                     'jco2h3co3h                      ', 'jbiacetoh                       ', &
                                     'jalkooh                         ', 'jalkno3                         ', &
                                     'jlisopacooh                     ', 'jlisopacno3                     ', &
                                     'jhpald                          ', 'jpacald                         ', &
                                     'jliecho                         ', 'jlieco3h                        ', &
                                     'jisopbooh                       ', 'jisopbno3                       ', &
                                     'jisopdooh                       ', 'jisopdno3                       ', &
                                     'jnisopooh                       ', 'jnc4cho                         ', &
                                     'jlnisooh                        ', 'jlhc4accho                      ', &
                                     'jlc578ooh                       ', 'jlhc4acco3h                     ', &
                                     'jhcoc5                          ', 'jc59ooh                         ', &
                                     'jlisopoohooh                    ', 'jlisopno3ooh                    ', &
                                     'jlisopno3no3                    ', 'jmboooh                         ', &
                                     'jibutaloh                       ', 'jibutalohooh                    ', &
                                     'jbepomuc                        ', 'jbigald1                        ', &
                                     'jtolooh                         ', 'jtepomuc                        ', &
                                     'jcatec1ooh                      ', 'jbigald2                        ', &
                                     'jbigald3                        ', 'jbigald4                        ', &
                                     'jterpooh                        ', 'jterprod1                       ', &
                                     'jterp2ooh                       ', 'jterprod2                       ', &
                                     'jterpno3                        ', 'jnterpno3                       ', &
                                     'jelvoc                          ', 'jo2_a                           ', &
                                     'jo2_b                           ', 'jh2o_a                          ', &
                                     'jh2o_b                          ', 'jh2o_c                          ', &
                                     'jcl2                            ', 'jcl2o2                          ', &
                                     'jclo                            ', 'jhcl                            ', &
                                     'jhocl                           ', 'jclono2_a                       ', &
                                     'jclono2_b                       ', 'joclo                           ', &
                                     'jbro                            ', 'jhbr                            ', &
                                     'jhobr                           ', 'jbrono_a                        ', &
                                     'jbrono_b                        ', 'jbrno2                          ', &
                                     'jbrono2_a                       ', 'jbrono2_b                       ', &
                                     'jbr2                            ', 'jbrcl                           ', &
                                     'jhf                             ', 'jsf6                            ', &
                                     'jch3br                          ', 'jch2br2                         ', &
                                     'jchbr3                          ', 'jch3cl                          ', &
                                     'jch3ccl3                        ', 'jcf3br                          ', &
                                     'jcf2clbr                        ', 'jccl4                           ', &
                                     'jcfc11                          ', 'jcfc12                          ', &
                                     'jcfc113                         ', 'jcfc114                         ', &
                                     'jcfc115                         ', 'jhcfc22                         ', &
                                     'jhcfc141b                       ', 'jhcfc142b                       ', &
                                     'jh1202                          ', 'jh2402                          ', &
                                     'jcof2                           ', 'jcofcl                          ', &
                                     'O_O2                            ', 'O_O3                            ', &
                                     'O_O                             ', 'O1D_N2                          ', &
                                     'O1D_O2_a                        ', 'O1D_O2_b                        ', &
                                     'O1D_H2O                         ', 'H_O2                            ', &
                                     'H_O3                            ', 'H_HO2_a                         ', &
                                     'H_HO2_b                         ', 'H_HO2_c                         ', &
                                     'H2_O                            ', 'H2_OH                           ', &
                                     'OH_O                            ', 'OH_OH_a                         ', &
                                     'OH_OH_b                         ', 'OH_O3                           ', &
                                     'HO2_O                           ', 'HO2_OH                          ', &
                                     'HO2_O3                          ', 'HO2_HO2                         ', &
                                     'H2O2_O                          ', 'H2O2_OH                         ', &
                                     'N_OH                            ', 'N_O2                            ', &
                                     'N_NO                            ', 'N_NO2                           ', &
                                     'NO_O                            ', 'NO_O3                           ', &
                                     'NO_HO2                          ', 'NO2_O_a                         ', &
                                     'NO2_O_b                         ', 'NO2_O3                          ', &
                                     'NO2_H                           ', 'NO2_OH                          ', &
                                     'NO3_O                           ', 'NO3_OH                          ', &
                                     'NO3_HO2                         ', 'NO3_NO                          ', &
                                     'HNO3_OH                         ', 'NO_OH                           ', &
                                     'HONO_OH                         ', 'NO2_HO2                         ', &
                                     'HO2NO2_OH                       ', 'HO2NO2                          ', &
                                     'NO2_NO3                         ', 'N2O5                            ', &
                                     'NH3_OH                          ', 'CO_OH_a                         ', &
                                     'CO_OH_b                         ', 'CH4_OH                          ', &
                                     'CH3OH_OH                        ', 'CH3O2_NO                        ', &
                                     'CH3O2_HO2                       ', 'CH3O2_CH3O2_a                   ', &
                                     'CH3O2_CH3O2_b                   ', 'CH2O_O                          ', &
                                     'CH2O_OH                         ', 'CH2O_HO2                        ', &
                                     'CH2O_NO3                        ', 'CH3OOH_OH                       ', &
                                     'HCOOH_OH                        ', 'HOCH2OO                         ', &
                                     'HOCH2OO_NO                      ', 'HOCH2OO_HO2                     ', &
                                     'HCN_OH                          ', 'CH3O2_NO2                       ', &
                                     'CH3O2NO2                        ', 'C2H2_OH                         ', &
                                     'C2H4_OH                         ', 'C2H4_O3                         ', &
                                     'C2H6_OH                         ', 'C2H5OH_OH                       ', &
                                     'CH3CHO_OH                       ', 'CH3CHO_NO3                      ', &
                                     'CH3COOOH_OH                     ', 'C2H5O2_NO                       ', &
                                     'C2H5O2_HO2                      ', 'C2H5O2_CH3O2                    ', &
                                     'C2H5O2_CH3CO3                   ', 'C2H5O2_C2H5O2                   ', &
                                     'C2H5OOH_OH                      ', 'CH3CO3_NO                       ', &
                                     'CH3CO3_NO2                      ', 'CH3CO3_HO2                      ', &
                                     'CH3CO3_CH3O2                    ', 'CH3CO3_CH3CO3                   ', &
                                     'CH3COOH_OH                      ', 'PAN                             ', &
                                     'PAN_OH                          ', 'EO2_NO                          ', &
                                     'EO2_HO2                         ', 'EO2_CH3O2                       ', &
                                     'EO2_CH3CO3                      ', 'EO_O2                           ', &
                                     'EO                              ', 'GLYOXAL_OH                      ', &
                                     'GLYOXAL_NO3                     ', 'GLYALD_OH                       ', &
                                     'GLYALD_NO3                      ', 'HOCH2CO3_NO2                    ', &
                                     'HOCH2CO3_HO2                    ', 'HOCH2CO3_CH3O2                  ', &
                                     'HOCH2CO3_CH3CO3                 ', 'HOCH2CO3_NO                     ', &
                                     'HOCH2CO3_NO3                    ', 'HOCH2CO2H_OH                    ', &
                                     'HOCH2CO3H_OH                    ', 'HCOCO3_CH3O2                    ', &
                                     'HCOCO3_CH3CO3                   ', 'HCOCO3_HO2                      ', &
                                     'HCOCO3_NO                       ', 'HCOCO3_NO3                      ', &
                                     'HCOCO2H_OH                      ', 'HCOCO3H_OH                      ', &
                                     'CH3CN_OH                        ', 'C3H6_OH                         ', &
                                     'C3H6_O3                         ', 'C3H6_NO3                        ', &
                                     'C3H8_OH                         ', 'CH3COCH3_OH                     ', &
                                     'C3H7O2_NO                       ', 'C3H7O2_HO2                      ', &
                                     'C3H7O2_CH3O2                    ', 'C3H7O2_CH3CO3                   ', &
                                     'C3H7OOH_OH                      ', 'PO2_NO                          ', &
                                     'PO2_NO3                         ', 'PO2_HO2                         ', &
                                     'PO2_CH3O2                       ', 'PO2_CH3CO3                      ', &
                                     'POOH_OH                         ', 'HYAC_OH                         ', &
                                     'CH3COCHO_OH                     ', 'CH3COCHO_NO3                    ', &
                                     'CH3COCH2O2_NO                   ', 'CH3COCH2O2_HO2                  ', &
                                     'CH3COCH2O2_CH3O2                ', 'CH3COCH2O2_CH3CO3               ', &
                                     'ROOH_OH                         ', 'PRONO3BO2_NO                    ', &
                                     'PRONO3BO2_NO3                   ', 'PRONO3BO2_HO2                   ', &
                                     'PRONO3BO2_CH3O2                 ', 'PRONO3BO2_CH3CO3                ', &
                                     'PR2O2HNO3_OH                    ', 'NOA_OH                          ', &
                                     'BIGENE_OH                       ', 'MEK_OH                          ', &
                                     'MEKO2_NO                        ', 'MEKO2_NOb                       ', &
                                     'MEKO2_HO2                       ', 'MEKO2_CH3O2                     ', &
                                     'MEKO2_CH3CO3                    ', 'MEKOOH_OH                       ', &
                                     'ENEO2_NO                        ', 'ENEO2_HO2                       ', &
                                     'ENEO2_CH3O2                     ', 'ENEO2_CH3CO3                    ', &
                                     'MACR_OH                         ', 'MACR_O3                         ', &
                                     'MACR_NO3                        ', 'MCO3_CH3O2                      ', &
                                     'MCO3_CH3CO3                     ', 'MCO3_HO2                        ', &
                                     'MCO3_NO                         ', 'MCO3_NO3                        ', &
                                     'MCO3_NO2                        ', 'MACRO2_CH3O2                    ', &
                                     'MACRO2_CH3CO3                   ', 'MACRO2_NO                       ', &
                                     'MACRO2_NOb                      ', 'MACRO2_NO3                      ', &
                                     'MACRO2_HO2                      ', 'MACRO2                          ', &
                                     'MACROOH_OH                      ', 'MACROH_OH                       ', &
                                     'MPAN                            ', 'MPAN_OH                         ', &
                                     'MACO2H_OH                       ', 'MACO3H_OH                       ', &
                                     'MVK_OH                          ', 'MVK_O3                          ', &
                                     'LHMVKABO2_CH3O2                 ', 'LHMVKABO2_CH3CO3                ', &
                                     'LHMVKABO2_HO2                   ', 'LHMVKABO2_NO                    ', &
                                     'LHMVKABO2_NOb                   ', 'LHMVKABO2_NO3                   ', &
                                     'MVKN_OH                         ', 'MACRN_OH                        ', &
                                     'LHMVKABOOH_OH                   ', 'CO2H3CHO_OH                     ', &
                                     'CO2H3CHO_NO3                    ', 'CO2H3CO3_CH3O2                  ', &
                                     'CO2H3CO3_CH3CO3                 ', 'CO2H3CO3_HO2                    ', &
                                     'CO2H3CO3_NO                     ', 'CO2H3CO3_NO3                    ', &
                                     'CO2H3CO3H_OH                    ', 'MALO2_NO2                       ', &
                                     'MALO2_NO                        ', 'MALO2_HO2                       ', &
                                     'MALO2_CH3O2                     ', 'MALO2_CH3CO3                    ', &
                                     'MDIALO2_HO2                     ', 'MDIALO2_NO                      ', &
                                     'MDIALO2_NO2                     ', 'MDIALO2_CH3O2                   ', &
                                     'MDIALO2_CH3CO3                  ', 'BIGALKANE_OH                    ', &
                                     'C5H8_OH                         ', 'C5H8_O3                         ', &
                                     'C5H8_NO3                        ', 'MBO_OH                          ', &
                                     'MBO_O3                          ', 'MBO_NO3                         ', &
                                     'ALKO2_NO                        ', 'ALKO2_NOb                       ', &
                                     'ALKO2_HO2                       ', 'ALKO2_CH3O2                     ', &
                                     'ALKO2_CH3CO3                    ', 'ALKOOH_OH                       ', &
                                     'ALKOH_OH                        ', 'ALKNO3_OH                       ', &
                                     'LISOPACO2_HO2                   ', 'LISOPACO2_NO                    ', &
                                     'LISOPACO2_NOb                   ', 'LISOPACO2_NO3                   ', &
                                     'LISOPACO2_CH3O2                 ', 'LISOPACO2_CH3CO3                ', &
                                     'LISOPACO2                       ', 'ISOPBO2                         ', &
                                     'ISOPDO2                         ', 'LISOPACOOH_OH                   ', &
                                     'ISOPAOH_OH                      ', 'LISOPACNO3_OH                   ', &
                                     'LIEPOX_OH                       ', 'LIECHO_OH                       ', &
                                     'LIECHO_NO3                      ', 'LIECO3_HO2                      ', &
                                     'LIECO3_NO                       ', 'LIECO3_NO3                      ', &
                                     'LIECO3H_OH                      ', 'IEC1O2_HO2                      ', &
                                     'IEC1O2_NO                       ', 'IEC1O2_NO3                      ', &
                                     'ISOPBO2_HO2                     ', 'ISOPBO2_NO                      ', &
                                     'ISOPBO2_NOb                     ', 'ISOPBO2_NO3                     ', &
                                     'ISOPBO2_CH3O2                   ', 'ISOPBO2_CH3CO3                  ', &
                                     'ISOPBOOH_OH                     ', 'ISOPBOH_OH                      ', &
                                     'ISOPBNO3_OH                     ', 'ISOPDO2_HO2                     ', &
                                     'ISOPDO2_NO                      ', 'ISOPDO2_NOb                     ', &
                                     'ISOPDO2_NO3                     ', 'ISOPDO2_CH3O2                   ', &
                                     'ISOPDO2_CH3CO3                  ', 'ISOPDOOH_OH                     ', &
                                     'ISOPDOH_OH                      ', 'ISOPDNO3_OH                     ', &
                                     'NISOPO2_HO2                     ', 'NISOPO2_NO                      ', &
                                     'NISOPO2_NO3                     ', 'NISOPO2_CH3O2                   ', &
                                     'NISOPO2_CH3CO3                  ', 'NISOPOOH_OH                     ', &
                                     'NC4CHO_OH                       ', 'NC4CHO_O3                       ', &
                                     'NC4CHO_NO3                      ', 'LNISO3_HO2                      ', &
                                     'LNISO3_NO                       ', 'LNISO3_NO3                      ', &
                                     'LNISO3_CH3O2                    ', 'LNISO3_CH3CO3                   ', &
                                     'LNISOOH_OH                      ', 'LHC4ACCHO_OH                    ', &
                                     'LHC4ACCHO_O3                    ', 'LHC4ACCHO_NO3                   ', &
                                     'LC578O2_NO                      ', 'LC578O2_NO3                     ', &
                                     'LC578O2_HO2                     ', 'LC578O2_CH3O2                   ', &
                                     'LC578O2_CH3CO3                  ', 'LC578OOH_OH                     ', &
                                     'LHC4ACCO3                       ', 'LHC4ACCO3_CH3O2                 ', &
                                     'LHC4ACCO3_CH3CO3                ', 'LHC4ACCO3_HO2                   ', &
                                     'LHC4ACCO3_NO                    ', 'LHC4ACCO3_NO3                   ', &
                                     'LHC4ACCO3_NO2                   ', 'LHC4ACCO2H_OH                   ', &
                                     'LHC4ACCO3H_OH                   ', 'LC5PAN1719                      ', &
                                     'LC5PAN1719_OH                   ', 'HCOC5_OH                        ', &
                                     'C59O2_CH3O2                     ', 'C59O2_CH3CO3                    ', &
                                     'C59O2_NO                        ', 'C59O2_NO3                       ', &
                                     'C59O2_HO2                       ', 'C59OOH_OH                       ', &
                                     'MBOO2_NO                        ', 'MBOO2_CH3O2                     ', &
                                     'MBOO2_CH3CO3                    ', 'MBOO2_HO2                       ', &
                                     'MBOOOH_OH                       ', 'IBUTALOH_OH                     ', &
                                     'IBUTALOHO2_NO                   ', 'IBUTALOHO2_HO2                  ', &
                                     'IBUTALOHO2_CH3O2                ', 'IBUTALOHO2_CH3CO3               ', &
                                     'IBUTALOHOOH_OH                  ', 'MBONO3O2_HO2                    ', &
                                     'MBONO3O2_NO                     ', 'MBONO3O2_NO3                    ', &
                                     'MBONO3O2_CH3O2                  ', 'MBONO3O2_CH3CO3                 ', &
                                     'DICARBO2_HO2                    ', 'DICARBO2_NO                     ', &
                                     'DICARBO2_NO2                    ', 'DICARBO2_CH3O2                  ', &
                                     'DICARBO2_CH3CO3                 ', 'HPALD_OH                        ', &
                                     'PACALD_OH                       ', 'LISOPOOHO2_HO2                  ', &
                                     'LISOPOOHO2_NO                   ', 'LISOPOOHO2_NOb                  ', &
                                     'LISOPOOHO2_NO3                  ', 'LISOPOOHO2_CH3O2                ', &
                                     'LISOPOOHO2_CH3CO3               ', 'LISOPNO3O2_HO2                  ', &
                                     'LISOPNO3O2_NO                   ', 'LISOPNO3O2_NOb                  ', &
                                     'LISOPNO3O2_NO3                  ', 'LISOPNO3O2_CH3O2                ', &
                                     'LISOPNO3O2_CH3CO3               ', 'LISOPOOHOOH_OH_a                ', &
                                     'LISOPOOHOOH_OH_b                ', 'LISOPNO3OOH_OH_a                ', &
                                     'LISOPNO3OOH_OH_b                ', 'LISOPNO3NO3_OH                  ', &
                                     'BENZ_OH                         ', 'PHENOL_OH                       ', &
                                     'PHENOL_NO3                      ', 'PHENO2_NO                       ', &
                                     'PHENO2_HO2                      ', 'PHENO2_CH3O2                    ', &
                                     'PHENO2_CH3CO3                   ', 'PHENOOH_OH                      ', &
                                     'C6H5O_NO2                       ', 'C6H5O_O3                        ', &
                                     'C6H5O2_NO                       ', 'C6H5O2_NO3                      ', &
                                     'C6H5O2_HO2                      ', 'C6H5O2_CH3O2                    ', &
                                     'C6H5O2_CH3CO3                   ', 'C6H5OOH_OH                      ', &
                                     'BENZO2_NO                       ', 'BENZO2_HO2                      ', &
                                     'BENZO2_CH3O2                    ', 'BENZO2_CH3CO3                   ', &
                                     'BENZOOH_OH                      ', 'CATECHOL_OH                     ', &
                                     'CATECHOL_NO3                    ', 'CATEC1O_NO2                     ', &
                                     'CATEC1O_O3                      ', 'CATEC1O2_HO2                    ', &
                                     'CATEC1O2_NO                     ', 'CATEC1O2_NO3                    ', &
                                     'CATEC1O2_CH3O2                  ', 'CATEC1O2_CH3CO3                 ', &
                                     'CATEC1OOH_OH                    ', 'TOL_OH                          ', &
                                     'CRESOL_OH                       ', 'CRESOL_NO3                      ', &
                                     'TOLO2_HO2                       ', 'TOLO2_NO                        ', &
                                     'TOLO2_CH3O2                     ', 'TOLO2_CH3CO3                    ', &
                                     'TOLOOH_OH                       ', 'BZOO_HO2                        ', &
                                     'BZOO_NO                         ', 'BZOO_CH3O2                      ', &
                                     'BZOO_CH3CO3                     ', 'BZOOH_OH                        ', &
                                     'BZALD_OH                        ', 'ACBZO2_NO2                      ', &
                                     'PBZNIT                          ', 'ACBZO2_NO                       ', &
                                     'ACBZO2_HO2                      ', 'ACBZO2_CH3O2                    ', &
                                     'ACBZO2_CH3CO3                   ', 'XYL_OH                          ', &
                                     'XYLOL_OH                        ', 'XYLOL_NO3                       ', &
                                     'XYLOLO2_NO                      ', 'XYLOLO2_HO2                     ', &
                                     'XYLOLO2_CH3O2                   ', 'XYLOLO2_CH3CO3                  ', &
                                     'XYLOLOOH_OH                     ', 'XYLENO2_HO2                     ', &
                                     'XYLENO2_NO                      ', 'XYLENO2_CH3O2                   ', &
                                     'XYLENO2_CH3CO3                  ', 'XYLENOOH_OH                     ', &
                                     'APIN_OH                         ', 'BPIN_OH                         ', &
                                     'LIMON_OH                        ', 'MYRC_OH                         ', &
                                     'BCARY_OH                        ', 'APIN_O3                         ', &
                                     'BPIN_O3                         ', 'LIMON_O3                        ', &
                                     'MYRC_O3                         ', 'BCARY_O3                        ', &
                                     'APIN_NO3                        ', 'BPIN_NO3                        ', &
                                     'LIMON_NO3                       ', 'MYRC_NO3                        ', &
                                     'BCARY_NO3                       ', 'TERPO2_NO                       ', &
                                     'TERPO2_HO2                      ', 'TERPO2_CH3O2                    ', &
                                     'TERPO2_CH3CO3                   ', 'TERPOOH_OH                      ', &
                                     'TERPROD1_OH                     ', 'TERPROD1_NO3                    ', &
                                     'TERPNO3_OH                      ', 'TERP2O2_NO                      ', &
                                     'TERP2O2_HO2                     ', 'TERP2O2_CH3O2                   ', &
                                     'TERP2O2_CH3CO3                  ', 'TERP2OOH_OH                     ', &
                                     'TERPROD2_OH                     ', 'NTERPO2_NO                      ', &
                                     'NTERPO2_HO2                     ', 'NTERPO2_CH3O2                   ', &
                                     'NTERPO2_CH3CO3                  ', 'NTERPO2_NO3                     ', &
                                     'NTERPNO3_OH                     ', 'ELVOC_OH                        ', &
                                     'O1D_N2O_a                       ', 'O1D_N2O_b                       ', &
                                     'O1D_O3                          ', 'O1D_CFC11                       ', &
                                     'O1D_CFC12                       ', 'O1D_CFC113                      ', &
                                     'O1D_CFC114                      ', 'O1D_CFC115                      ', &
                                     'O1D_HCFC22                      ', 'O1D_HCFC141B                    ', &
                                     'O1D_HCFC142B                    ', 'O1D_CCL4                        ', &
                                     'O1D_CH3BR                       ', 'O1D_CF2CLBR                     ', &
                                     'O1D_CF3BR                       ', 'O1D_H1202                       ', &
                                     'O1D_H2402                       ', 'O1D_CHBR3                       ', &
                                     'O1D_CH2BR2                      ', 'O1D_COF2                        ', &
                                     'O1D_COFCL                       ', 'O1D_CH4_a                       ', &
                                     'O1D_CH4_b                       ', 'O1D_CH4_c                       ', &
                                     'O1D_H2                          ', 'O1D_HCL                         ', &
                                     'O1D_HBR                         ', 'O1D_HCN                         ', &
                                     'CL_O3                           ', 'CL_H2                           ', &
                                     'CL_H2O2                         ', 'CL_HO2_a                        ', &
                                     'CL_HO2_b                        ', 'CL2O2                           ', &
                                     'CLO_O                           ', 'CLO_OH_a                        ', &
                                     'CLO_OH_b                        ', 'CLO_HO2                         ', &
                                     'CLO_NO                          ', 'CLO_NO2                         ', &
                                     'CLO_CLO_a                       ', 'CLO_CLO_b                       ', &
                                     'CLO_CLO_c                       ', 'CLO_CLO_d                       ', &
                                     'HCL_OH                          ', 'HCL_O                           ', &
                                     'HOCL_O                          ', 'HOCL_CL                         ', &
                                     'HOCL_OH                         ', 'CLONO2_O                        ', &
                                     'CLONO2_OH                       ', 'CLONO2_CL                       ', &
                                     'BR_O3                           ', 'BR_HO2                          ', &
                                     'BR_NO2                          ', 'BRONO                           ', &
                                     'BRO_O                           ', 'BRO_OH                          ', &
                                     'BRO_HO2                         ', 'BRO_NO                          ', &
                                     'BRO_NO2                         ', 'BRO_CLO_a                       ', &
                                     'BRO_CLO_b                       ', 'BRO_CLO_c                       ', &
                                     'BRO_BRO_a                       ', 'BRO_BRO_b                       ', &
                                     'HBR_OH                          ', 'HBR_O                           ', &
                                     'HOBR_O                          ', 'BRONO2_O                        ', &
                                     'BRONO2_BR                       ', 'BR2_OH                          ', &
                                     'F_H2O                           ', 'F_H2                            ', &
                                     'F_CH4                           ', 'F_HNO3                          ', &
                                     'CL_CH2O                         ', 'CL_CH4                          ', &
                                     'CL_CH3CN                        ', 'CL_C2H2                         ', &
                                     'CL_C2H4                         ', 'CL_C2H6                         ', &
                                     'CL_CH3O2                        ', 'CL_CH3OH                        ', &
                                     'CL_CH3OOH                       ', 'CL_CH3CHO                       ', &
                                     'CL_GLYALD                       ', 'CL_GLYOXAL                      ', &
                                     'CL_C3H8                         ', 'CL_CH3COCH3                     ', &
                                     'CL_HYAC                         ', 'CL_BIGALKANE                    ', &
                                     'CL_MEK                          ', 'CLO_CH3O2                       ', &
                                     'BR_CH2O                         ', 'BR_CH3CHO                       ', &
                                     'BRO_CH3O2                       ', 'SO2_OH                          ', &
                                     'DMS_OH_a                        ', 'DMS_OH_b                        ', &
                                     'DMS_NO3                         ', 'DMS_CL                          ', &
                                     'DMS_BR                          ', 'DMS_BRO                         ', &
                                     'DMSO_OH                         ', 'CH3SO2_M                        ', &
                                     'CH3SO2_O3                       ', 'CH3SO3_NO2                      ', &
                                     'CH3BR_OH                        ', 'CH3BR_CL                        ', &
                                     'CH2BR2_OH                       ', 'CHBR3_OH                        ', &
                                     'CH2BR2_CL                       ', 'CHBR3_CL                        ', &
                                     'CH3CL_CL                        ', 'CH3CL_OH                        ', &
                                     'CH3CCL3_OH                      ', 'HCFC22_OH                       ', &
                                     'HCFC141B_OH                     ', 'HCFC142B_OH                     ', &
                                     'het_O3_tropo                    ', 'het_HO2_tropo                   ', &
                                     'het_NO3_tropo                   ', 'het_NO2_tropo                   ', &
                                     'het_HNO3_tropo                  ', 'het_N2O5_tropo                  ', &
                                     'het_N2O5_strat_sulfate          ', 'het_CLONO2_strat_sulfate        ', &
                                     'het_BRONO2_strat_sulfate        ', 'het_CLONO2_HCL_strat_sulfate    ', &
                                     'het_HOCL_HCL_strat_sulfate      ', 'het_HOBR_HCL_strat_sulfate      ', &
                                     'het_N2O5_strat_NAD              ', 'het_CLONO2_strat_NAD            ', &
                                     'het_CLONO2_HCL_strat_NAD        ', 'het_HOCL_HCL_strat_NAD          ', &
                                     'het_BRONO2_strat_NAD            ', 'het_N2O5_strat_ice              ', &
                                     'het_CLONO2_strat_ice            ', 'het_BRONO2_strat_ice            ', &
                                     'het_CLONO2_HCL_strat_ice        ', 'het_HOCL_HCL_strat_ice          ', &
                                     'het_HOBR_HCL_strat_ice          ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
                                       11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
                                       31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50, &
                                       51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
                                       61, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
                                       71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
                                       81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
                                       91, 92, 93, 94, 95, 96, 97, 98, 99, 100, &
                                      101, 102, 103, 104, 105, 106, 107, 108, 109, 110, &
                                      111, 112, 113, 114, 115, 116, 117, 118, 119, 120, &
                                      121, 122, 123, 124, 125, 126, 127, 128, 129, 130, &
                                      131, 132, 133, 134, 135, 136, 137, 138, 139, 140, &
                                      141, 142, 143, 144, 145, 146, 147, 148, 149, 150, &
                                      151, 152, 153, 154, 155, 156, 157, 158, 159, 160, &
                                      161, 162, 163, 164, 165, 166, 167, 168, 169, 170, &
                                      171, 172, 173, 174, 175, 176, 177, 178, 179, 180, &
                                      181, 182, 183, 184, 185, 186, 187, 188, 189, 190, &
                                      191, 192, 193, 194, 195, 196, 197, 198, 199, 200, &
                                      201, 202, 203, 204, 205, 206, 207, 208, 209, 210, &
                                      211, 212, 213, 214, 215, 216, 217, 218, 219, 220, &
                                      221, 222, 223, 224, 225, 226, 227, 228, 229, 230, &
                                      231, 232, 233, 234, 235, 236, 237, 238, 239, 240, &
                                      241, 242, 243, 244, 245, 246, 247, 248, 249, 250, &
                                      251, 252, 253, 254, 255, 256, 257, 258, 259, 260, &
                                      261, 262, 263, 264, 265, 266, 267, 268, 269, 270, &
                                      271, 272, 273, 274, 275, 276, 277, 278, 279, 280, &
                                      281, 282, 283, 284, 285, 286, 287, 288, 289, 290, &
                                      291, 292, 293, 294, 295, 296, 297, 298, 299, 300, &
                                      301, 302, 303, 304, 305, 306, 307, 308, 309, 310, &
                                      311, 312, 313, 314, 315, 316, 317, 318, 319, 320, &
                                      321, 322, 323, 324, 325, 326, 327, 328, 329, 330, &
                                      331, 332, 333, 334, 335, 336, 337, 338, 339, 340, &
                                      341, 342, 343, 344, 345, 346, 347, 348, 349, 350, &
                                      351, 352, 353, 354, 355, 356, 357, 358, 359, 360, &
                                      361, 362, 363, 364, 365, 366, 367, 368, 369, 370, &
                                      371, 372, 373, 374, 375, 376, 377, 378, 379, 380, &
                                      381, 382, 383, 384, 385, 386, 387, 388, 389, 390, &
                                      391, 392, 393, 394, 395, 396, 397, 398, 399, 400, &
                                      401, 402, 403, 404, 405, 406, 407, 408, 409, 410, &
                                      411, 412, 413, 414, 415, 416, 417, 418, 419, 420, &
                                      421, 422, 423, 424, 425, 426, 427, 428, 429, 430, &
                                      431, 432, 433, 434, 435, 436, 437, 438, 439, 440, &
                                      441, 442, 443, 444, 445, 446, 447, 448, 449, 450, &
                                      451, 452, 453, 454, 455, 456, 457, 458, 459, 460, &
                                      461, 462, 463, 464, 465, 466, 467, 468, 469, 470, &
                                      471, 472, 473, 474, 475, 476, 477, 478, 479, 480, &
                                      481, 482, 483, 484, 485, 486, 487, 488, 489, 490, &
                                      491, 492, 493, 494, 495, 496, 497, 498, 499, 500, &
                                      501, 502, 503, 504, 505, 506, 507, 508, 509, 510, &
                                      511, 512, 513, 514, 515, 516, 517, 518, 519, 520, &
                                      521, 522, 523, 524, 525, 526, 527, 528, 529, 530, &
                                      531, 532, 533, 534, 535, 536, 537, 538, 539, 540, &
                                      541, 542, 543, 544, 545, 546, 547, 548, 549, 550, &
                                      551, 552, 553, 554, 555, 556, 557, 558, 559, 560, &
                                      561, 562, 563, 564, 565, 566, 567, 568, 569, 570, &
                                      571, 572, 573, 574, 575, 576, 577, 578, 579, 580, &
                                      581, 582, 583, 584, 585, 586, 587, 588, 589, 590, &
                                      591, 592, 593, 594, 595, 596, 597, 598, 599, 600, &
                                      601, 602, 603, 604, 605, 606, 607, 608, 609, 610, &
                                      611, 612, 613, 614, 615, 616, 617, 618, 619, 620, &
                                      621, 622, 623, 624, 625, 626, 627, 628, 629, 630, &
                                      631, 632, 633, 634, 635, 636, 637, 638, 639, 640, &
                                      641, 642, 643, 644, 645, 646, 647, 648, 649, 650, &
                                      651, 652, 653, 654, 655, 656, 657, 658, 659, 660, &
                                      661, 662, 663, 664, 665, 666, 667, 668, 669, 670, &
                                      671, 672, 673, 674, 675, 676, 677, 678, 679, 680, &
                                      681, 682, 683, 684, 685, 686, 687, 688, 689, 690, &
                                      691, 692, 693, 694, 695, 696, 697, 698, 699, 700, &
                                      701, 702, 703, 704, 705, 706, 707, 708, 709, 710, &
                                      711, 712, 713, 714, 715, 716, 717, 718, 719, 720, &
                                      721, 722, 723, 724, 725, 726, 727, 728, 729, 730, &
                                      731, 732, 733 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(message_text,*) 'failed to allocate pht_alias_lst; error = ',ios
         CALL finish( 'set_sim_dat',message_text )
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(message_text,*) 'failed to allocate pht_alias_mult; error = ',ios
         CALL finish( 'set_sim_dat',message_text )
      end if
      pht_alias_lst(:,1) = (/ '                ', '                ', '                ', '                ', &
                              'userdefined     ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', 'userdefined     ', 'userdefined     ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ' /)
      pht_alias_mult(:,1) = (/ 1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', 'jno3_a          ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              'jho2no2_a       ', 'jho2no2_b       ', '                ', 'jh2o2           ', &
                              'jch3ooh         ', '                ', 'jch3ooh         ', '                ', &
                              '                ', 'jh2o2           ', 'jglyoxal        ', 'jch3ooh         ', &
                              'jacet           ', 'jch3ooh         ', 'jch3ooh         ', '                ', &
                              'jglyoxal        ', 'jch3ooh         ', 'jch3ooh         ', 'jglyoxal        ', &
                              'jacet           ', 'jch3ooh         ', 'jglyoxal        ', '                ', &
                              '                ', 'jch3ooh         ', 'jhyac           ', 'jpan            ', &
                              'jch3ooh         ', '                ', 'jch3ooh         ', 'jglyoxal        ', &
                              'jglyoxal        ', 'jhyac           ', 'jhyac           ', 'jglyoxal        ', &
                              'jch3ooh         ', 'jhno3           ', 'jch3ooh         ', 'jhno3           ', &
                              'jmacr_a         ', 'jmacr_a         ', 'jglyoxal        ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jpan            ', 'jch3ooh         ', 'jpan            ', &
                              'jch3ooh         ', 'jglyoxal        ', 'jch3ooh         ', 'jglyoxal        ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jglyoxal        ', 'jch3ooh         ', &
                              'jch3ooh         ', 'jch3ooh         ', 'jpan            ', 'jch3ooh         ', &
                              'jglyald         ', 'jh2o2           ', 'jmacr_a         ', 'jno2            ', &
                              'jch3ooh         ', 'jmacr_a         ', 'jch3ooh         ', 'jno2            ', &
                              'jno2            ', 'jno2            ', 'jch3ooh         ', 'jch3cho         ', &
                              'jch3ooh         ', 'jch3cho         ', 'jpan            ', 'jch3ooh         ', &
                              'jch3ooh         ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', 'jbro            ', 'jbro            ', 'jbrono2_b       ', &
                              '                ', '                ', 'jbro            ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', 'jcfcl3          ', 'jcf2cl2         ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ' /)
      pht_alias_mult(:,2) = (/ 1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , .28_dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 0.28_dp , &
                               1._dp , 2._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 0.33_dp , &
                               1._dp , 1._dp , 0.33_dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 0.5_dp , 3.5_dp , 2._dp , 2._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               500._dp , 1000._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 5.6_dp , &
                               1._dp , 1._dp , 2._dp , 2._dp , 1._dp , &
                               2._dp , 2._dp , 1._dp , 2._dp , 1._dp , &
                               1._dp , 0.28_dp , 4._dp , .14_dp , 1._dp , &
                               4._dp , 1._dp , .20_dp , .20_dp , .006_dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 3._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 0.5_dp , 0.5_dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp , 1._dp , 1._dp , 1._dp , &
                               1._dp , 1._dp /)
      inv_lst(:nfs) = (/ 'M                       ', 'N2                      ' /)
      end subroutine set_sim_dat
      end module mo_moz_subs
