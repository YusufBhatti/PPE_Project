!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_memory_g3a

  USE mo_kind,        ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g3a ! construct the g3a table

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: geospm(:,:)
  REAL(dp), POINTER, PUBLIC :: wsm(:,:)
  REAL(dp), POINTER, PUBLIC :: wlm(:,:)
  REAL(dp), POINTER, PUBLIC :: snm(:,:)
  REAL(dp), POINTER, PUBLIC :: slmm(:,:)
  REAL(dp), POINTER, PUBLIC :: az0m(:,:)
  REAL(dp), POINTER, PUBLIC :: albm(:,:)
  REAL(dp), POINTER, PUBLIC :: forestm(:,:)
  REAL(dp), POINTER, PUBLIC :: vgratm(:,:)
  REAL(dp), POINTER, PUBLIC :: vltm(:,:)
  REAL(dp), POINTER, PUBLIC :: wsmxm(:,:)
  REAL(dp), POINTER, PUBLIC :: faom(:,:)
  REAL(dp), POINTER, PUBLIC :: apsm(:,:)
  REAL(dp), POINTER, PUBLIC :: aprlm(:,:)
  REAL(dp), POINTER, PUBLIC :: aprcm(:,:)
  REAL(dp), POINTER, PUBLIC :: aprsm(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrgwm(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrgwm(:,:)
  REAL(dp), POINTER, PUBLIC :: vdisgwm(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcovm(:,:)
  REAL(dp), POINTER, PUBLIC :: temp2m(:,:)
  REAL(dp), POINTER, PUBLIC :: dew2m(:,:)
  REAL(dp), POINTER, PUBLIC :: wind10m(:,:)
  REAL(dp), POINTER, PUBLIC :: u10m(:,:)
  REAL(dp), POINTER, PUBLIC :: v10m(:,:)
  REAL(dp), POINTER, PUBLIC :: sradsm(:,:)
  REAL(dp), POINTER, PUBLIC :: tradsm(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0m(:,:)
  REAL(dp), POINTER, PUBLIC :: trad0m(:,:)
  REAL(dp), POINTER, PUBLIC :: vdism(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrm(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrm(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfsm(:,:)
  REAL(dp), POINTER, PUBLIC :: evapm(:,:)
  REAL(dp), POINTER, PUBLIC :: ahflm(:,:)
  REAL(dp), POINTER, PUBLIC :: emterm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsolm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: srad0um(:,:)
  REAL(dp), POINTER, PUBLIC :: sradsum(:,:)
  REAL(dp), POINTER, PUBLIC :: tradsum(:,:)
  REAL(dp), POINTER, PUBLIC :: albedom(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo_vism(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo_nirm(:,:)
  REAL(dp), POINTER, PUBLIC :: tsurfm(:,:)
  REAL(dp), POINTER, PUBLIC :: seaicem(:,:)
  REAL(dp), POINTER, PUBLIC :: sicedm(:,:)
  REAL(dp), POINTER, PUBLIC :: wind10wm(:,:)
  REAL(dp), POINTER, PUBLIC :: glacm(:,:)
  REAL(dp), POINTER, PUBLIC :: sswdifnirm(:,:)
  REAL(dp), POINTER, PUBLIC :: sswvism(:,:)
  REAL(dp), POINTER, PUBLIC :: sswdifvism(:,:)
  REAL(dp), POINTER, PUBLIC :: sswparm(:,:)
  REAL(dp), POINTER, PUBLIC :: sswdifparm(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: aclcacm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: t2maxm(:,:)
  REAL(dp), POINTER, PUBLIC :: t2minm(:,:)
  REAL(dp), POINTER, PUBLIC :: wimaxm(:,:)
  REAL(dp), POINTER, PUBLIC :: topmaxm(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcvm(:,:)
  REAL(dp), POINTER, PUBLIC :: qvim(:,:)
  REAL(dp), POINTER, PUBLIC :: xlvim(:,:)
  REAL(dp), POINTER, PUBLIC :: xivim(:,:)
  REAL(dp), POINTER, PUBLIC :: runlndm(:,:)
  REAL(dp), POINTER, PUBLIC :: rgcgnm(:,:)
  REAL(dp), POINTER, PUBLIC :: sodifm(:,:)
  REAL(dp), POINTER, PUBLIC :: srafsm(:,:)
  REAL(dp), POINTER, PUBLIC :: trafsm(:,:)
  REAL(dp), POINTER, PUBLIC :: sraf0m(:,:)
  REAL(dp), POINTER, PUBLIC :: traf0m(:,:)
  REAL(dp), POINTER, PUBLIC :: emtefm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsofm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: srad0dm(:,:)
  REAL(dp), POINTER, PUBLIC :: acdncm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: aprfluxm(:,:)     ! for middle atmosphere only

CONTAINS

  SUBROUTINE construct_g3a ! (lnlon, lnlev, lngl, nlon, nlev, ngl)

    USE mo_memory_g3b

    ! construct the g3a table
    !
    ! all information specific to this table is set in this subroutine

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    geospm    => geosp
    wsm       => ws
    wlm       => wl
    snm       => sn
    slmm      => slm
    az0m      => az0
!!$    albm      => alb
    forestm   => forest
    vgratm    => vgrat
    vltm      => vlt
    wsmxm     => wsmx
    faom      => fao
    apsm      => aps
    aprlm     => aprl
    aprcm     => aprc
    aprsm     => aprs
    ustrgwm   => ustrgw
    vstrgwm   => vstrgw
    vdisgwm   => vdisgw
    aclcovm   => aclcov
    temp2m    => temp2
    dew2m     => dew2
    wind10m   => wind10
    u10m      => u10
    v10m      => v10
    sradsm    => srads
    tradsm    => trads
    srad0m    => srad0
    trad0m    => trad0
    vdism     => vdis
    ustrm     => ustr
    vstrm     => vstr
    ahfsm     => ahfs
    evapm     => evap
    ahflm     => ahfl
    emterm    => emter
    trsolm    => trsol
    srad0um   => srad0u
    sradsum   => sradsu
    tradsum   => tradsu
    albedom   => albedo
    albedo_vism => albedo_vis
    albedo_nirm => albedo_nir
    tsurfm    => tsurf
    seaicem   => seaice
    sicedm    => siced
    wind10wm  => wind10w
    glacm     => glac
    sswdifnirm   => sswdifnir
    sswvism   => sswvis_frc
    sswdifvism   => sswdifvis
    sswparm   => sswpar
    sswdifparm   => sswdifpar
    aclcm     => aclc
    aclcacm   => aclcac
    t2maxm    => t2max
    t2minm    => t2min
    wimaxm    => wimax
    topmaxm   => topmax
    aclcvm    => aclcv
    qvim      => qvi
    xlvim     => xlvi
    xivim     => xivi
    runlndm   => runlnd
    rgcgnm    => rgcgn
    sodifm    => sodif
    srafsm    => srafs
    trafsm    => trafs
    sraf0m    => sraf0
    traf0m    => traf0
    emtefm    => emtef
    trsofm    => trsof
    srad0dm   => srad0d
    acdncm    => acdnc
    aprfluxm  => aprflux

  END SUBROUTINE construct_g3a

END MODULE mo_memory_g3a
