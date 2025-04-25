!*    NAMELIST *DISTURBANCE_CTL* To control JSBACH vegetation disturbance algorithms.
!
!
!   NAME                   TYPE       PURPOSE
!   ----                   ----       -------
!
!  *FRAC_WOOD_2_ATMOS*     REAL       Fraction of wood    pool to release to atmosphere by shrinking tiles by LCC
!  *FRAC_GREEN_2_ATMOS*    REAL       Fraction of green   pool to release to atmosphere by shrinking tiles by LCC
!  *FRAC_RESERVE_2_ATMOS*  REAL       Fraction of reserve pool to release to atmosphere by shrinking tiles by LCC 
!  *FRAC_MOBILE_2_ATMOS*   REAL       Fraction of mobile plant nitrogen to release to atmosphere by shrinking tiles by LCC
!  *FRAC_HARVEST_2_ATMOS*  REAL       Fraction of carbon harvested from natural pfts to release to atmosphere
!  *cn_green*              REAL       Ratio of carbon to nitrogen in green pool
!  *cn_woods*              REAL       Ratio of carbon to nitrogen in wood pool
!  *cn_litter_green*       REAL       Ratio of carbon to nitrogen in green litter pool
!  *cn_litter_wood*        REAL       Ratio of carbon to nitrogen in wood litter pool
!  *cn_slow*               REAL       Ratio of carbon to nitrogen in slow pool
!  *tau_onSite*            REAL       10% decay time of anthropogenic onSite       pool               [days]
!  *tau_paper*             REAL       10% decay time of anthropogenic paper        pool               [days]
!  *tau_construction*      REAL       10% decay time of anthropogenic construction pool               [days]


NAMELIST /CBAL_PARAMETERS_CTL/                                      &
    frac_wood_2_atmos, frac_green_2_atmos, frac_reserve_2_atmos,    &
    frac_mobile_2_atmos, frac_harvest_2_atmos,                      &
    cn_green, cn_woods, cn_litter_green, cn_litter_wood, cn_slow,   &
    tau_onSite, tau_paper, tau_construction

