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
  fracdusol_0d,     & !fraction of dust particules in all soluble modes (as scalar)
  fracduai_0d,      & !fraction of dust particules in the insoluble accumulation mode (as scalar)
  fracduci_0d,      & !fraction of dust particules in the insoluble coarse mode (as scalar)
  fracbcsol_0d,     & !fraction of BC particules in all soluble modes (as scalar)
  fracbcinsol_0d,   & !fraction of BC particules in all insoluble modes (as scalar)
  rwetki_0d,        & !wet radius for the insoluble aitken mode (as scalar)
  rwetai_0d,        & !wet radius for the insoluble accumulation mode (as scalar)
  rwetci_0d,        & !wet radius for the insoluble coarse mode (as scalar)
  rwetas_0d,        & !wet radius for the soluble accumulation mode (as scalar):
  l_variable_ss,    & !switch for taking into account variable supersat (-> variable updraft velocity)
  ninput_type_CCN,  & !paramter to describe the CCN input type
                      ! ninput_type_CCN = 2 (2D CCN input @1km and 8km)
                      ! ninput_type_CCN = 3 (3D CCN input) 
  ninput_type_IN      !paramter to describe the IN input type
                      ! ninput_type_IN = 0 (0D IN input, ie scalar)
                      ! ninput_type_IN = 3 (3D IN input)
