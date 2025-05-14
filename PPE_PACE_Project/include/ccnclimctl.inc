!------------------------------------------------------------------------------------------
!
! ccnclimctl.inc
!
!------------------------------------------------------------------------------------------
!
! variable name   & description & default \\
!
!------------------------------------------------------------------------------------------
NAMELIST /ccnclimctl/ &
  fracdusol,     & !fraction of dust particules in all soluble modes
  fracduai,      & !fraction of dust particules in the insoluble accumulation mode
  fracduci,      & !fraction of dust particules in the insoluble coarse mode
  fracbcsol,     & !fraction of BC particules in all soluble modes
  fracbcinsol,   & !fraction of BC particules in all insoluble modes
  rwetki,        & !wet radius for the insoluble aitken mode
  rwetai,        & !wet radius for the insoluble accumulation mode
  rwetci,        & !wet radius for the insoluble coarse mode
  rwetas           !wet radius for the soluble accumulation mode
