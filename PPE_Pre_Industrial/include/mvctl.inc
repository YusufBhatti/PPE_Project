! mvctl.inc
!
!    -----------------------------------------------------------------------
NAMELIST /mvctl/ &
  putmean, &       ! time interval for writing mean values
  meannam, &       ! list of species for which mean values are calculated
  stddev           ! mean of squares list (integers -1, 0 ,1)
