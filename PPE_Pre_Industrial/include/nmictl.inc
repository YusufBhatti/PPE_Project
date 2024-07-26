! nmictl.inc
! 
!     ------------------------------------------------------------------
NAMELIST /nmictl/ &
    ntpre,        &! INTEGER number of timesteps skiped before accumulation of tendencies
    ntdia,        &! INTEGER number of accumulation time steps for diabatic tendencies
    ntiter,       &! INTEGER number of iteration time steps
    pcut,         &! INTEGER cut off period for fast gravity modes in hours (nudging)
    pcutd,        &! INTEGER cut off period for filtering of diabatic tendencies (initialization)
    lnmi_cloud,   &! LOGICAL run NMI with cloud forcing
    dt_nmi_start   ! start date of initialization procedure
!     ------------------------------------------------------------------


