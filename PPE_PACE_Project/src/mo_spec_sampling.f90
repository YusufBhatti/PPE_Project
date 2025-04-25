!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
!! @brief Module containing various strategies for sampling the electromagnetic spectrum
!!
!! @remarks 
!!   This module describes various strategies for sampling the longwave and shortwave spectra
!!   for radiation calculations in dymamical models. The g-points within each spectrum are 
!!   grouped into teams of identical length. One or more teams are then chosen randomly and the fluxes
!!   from those g-points are used as a proxy for the full radiation calculation.
!!
!! @author Robert Pincus, U. Colorado, while visiting MPI-M, Hamburg (2011-07)
!!           Additions for spectral teams, 2012-03
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Written by Robert Pincus 
!!

MODULE mo_spec_sampling 
  USE mo_exception,   ONLY : em_warn, message_text, message
  USE mo_rrtm_params, ONLY : ngptsw, ngptlw
  USE mo_random_numbers, & 
                      ONLY: get_random
  USE mo_kind,        ONLY: wp
  IMPLICIT NONE
  PRIVATE 

  INTEGER, PRIVATE :: i
  
  !
  ! Team choices - Longwave
  !
  INTEGER, DIMENSION(ngptlw, 1), TARGET :: &
    & lw_teams_all = RESHAPE( (/ (i, i = 1, ngptlw) /), SHAPE = (/ ngptlw, 1 /) )
  INTEGER, DIMENSION(1, ngptlw), TARGET :: &
    & lw_teams_one = RESHAPE( (/ (i, i = 1, ngptlw) /), SHAPE = (/ 1, ngptlw /) )
  INTEGER, DIMENSION(4, 35), TARGET :: &
    & lw_4 = RESHAPE(               &
    &        (/   3, 134,  36, 133, &
    &             1, 122,  94, 106, &
    &            70,  22, 113,  60, &
    &             2,  86,  75, 125, &
    &            11,  95,  48,  21, &
    &            69,  10, 112, 128, &
    &             4,  76,  84, 126, &
    &            12, 108,  50, 129, &
    &            71,  51, 105,  44, &
    &            23,   9,  49,  83, &
    &            24,  67,  66,  82, &
    &            13, 136,  65,  34, &
    &            25,  88, 121,  80, &
    &            14, 114,  85,  79, &
    &            26,  68, 107, 127, &
    &            15,  38,  87,  99, &
    &            27,  52, 132,  90, &
    &            72,  93, 111,  64, &
    &            55,  20,  19,  63, &
    &             5, 123,  74,   7, &
    &            56,  35, 103,  30, &
    &            54, 120,  31, 102, &
    &            16, 130,  89,  62, &
    &            28, 138,  91,  32, &
    &            57,  37, 104, 101, &
    &            53,  47,  45,  33, &
    &           117,   8,  61,  78, &
    &            58,  96,  73,  18, &
    &            39, 140,  92,  42, &
    &            29, 137,  77,  59, &
    &            17, 139,  97, 100, &
    &            40, 131, 115, 135, &
    &           118,  81,   6, 109, &
    &            98,  46, 110,  43, &
    &            41, 124, 116, 119 /), SHAPE = (/ 4, 35 /) ) 
  INTEGER, DIMENSION(5, 28), TARGET :: &
    & lw_5 = RESHAPE(                  &
    &        (/   3,   9, 137, 132, 106, &
    &             1,  76,  87,  94, 121, &
    &            70,  52, 113,  43,  48, &
    &             2, 134, 123,  97,  39, &
    &            11, 130,  65,  21,  64, &
    &            69,  96, 112,  60,  33, &
    &             4,  86,  75, 125, 126, &
    &            12, 114,  50, 129, 118, &
    &            71,  22,  46, 127,  18, &
    &            23,  67,  36, 115,  32, &
    &            24,  37,  66, 105,  30, &
    &            13,  38,  84,  89, 135, &
    &            25,  51,  95,  80,  42, &
    &            14, 138,  85,  77,  63, &
    &            26,  10, 140,  99,  82, &
    &            15,  68, 107,  98,  81, &
    &            27,  88, 124, 100, 133, &
    &            72, 120, 103, 102,   7, &
    &            55,  20,  31, 119,  29, &
    &             5, 131,  74,  92,  59, &
    &            56,  35,   6, 109, 128, &
    &            54,  93,  45, 101,  62, &
    &            16, 139,  90,  73,  49, &
    &            28, 136,  91, 111,  40, &
    &            57, 108, 104,  44,  34, &
    &            53,  47,  19,  79,  17, &
    &           117,   8,  61,  78,  83, &
    &            58, 122, 110, 116,  41 /), SHAPE = (/ 5, 28 /) ) 
  INTEGER, DIMENSION(7, 20), TARGET :: &
    & lw_7 = RESHAPE(                  & 
    &        (/   3, 108, 137,  84, 132, 126,  54, &
    &             1,  38, 131, 107, 121,  90,  56, &
    &            70,  76, 113,  19,  43,   7,  36, &
    &             2, 136, 123,  94,  53, 133,  64, &
    &            11, 114,  52,  21, 102,  82,  16, &
    &            69,  10,  47, 103,  18,  62,  17, &
    &             4, 122, 124,  75,  73,  99,  35, &
    &            12,   9,  86, 129, 101,  93, 117, &
    &            71,  37,  46, 116, 119,  32,  29, &
    &            23,  88,   8,  45,  80,  41,  50, &
    &            24,  22,  51, 111,  79,  30,  33, &
    &            13,  68, 140,  89,  81,  59, 106, &
    &            25, 130,  66, 104,  60,  92,  49, &
    &            14, 134,  85,  97, 105,  58,  65, &
    &            26, 138,  95,  61,  91,  40,  48, &
    &            15,  96,  87,  77, 100, 128,  39, &
    &            27,  67, 125, 110,  78,  83,  57, &
    &            72, 120,   6,  44, 112,  42,  63, &
    &            55,  20,  31, 115, 109, 135,  28, &
    &             5, 139,  74,  98, 127, 118,  34 /), SHAPE = (/ 7, 20 /) ) 
  INTEGER, DIMENSION(10, 14), TARGET :: &
    & lw_10 = RESHAPE(                  &
    &         (/   3,  96,  76, 124,  95,  94,  90,  55, 121,  48, &
    &              1, 122, 140,  86, 107,  89,  72,  33,  27,  85, &
    &             70,  22,  20,  19, 116, 109, 128,  17,  63,  15, &
    &              2, 134, 123,  74, 125,  73,  54, 129, 117,  50, &
    &             11,  10,  52,   7, 110, 101,  82,  40, 106,  87, &
    &             69,  88, 112,  45, 115,  18,  93, 100,  21,  16, &
    &              4, 136, 137,  75,  53, 126,  99,  62,  58,  84, &
    &             12,  38, 138, 113, 119,  80,  59,  98,   8,  66, &
    &             71, 108,  47,  31,  44, 127,  30,  32,  39,  65, &
    &             23,  37,  35,   6, 111,  60,  92, 132,  29,  26, &
    &             24,   9, 114,  46, 103,  81,  41,  77,  34,   5, &
    &             13,  68, 131,  83, 102,  79, 118, 120,  56,  36, &
    &             25,  67,  51,  61,  78, 105,  42, 135,  57,  49, &
    &             14, 139, 130,  97, 104,  43,  91, 133,  28,  64 /), SHAPE = (/ 10, 14 /) ) 
  INTEGER, DIMENSION(14, 10), TARGET :: &
    & lw_14 = RESHAPE(                  & 
    &         (/   3,  88, 136,  76, 137,  74,  97,  90,  55, 132,  72, 106,  87,  13, &
    &              1, 134, 138, 123,  37,  89,  98,  92,  41, 129, 117,  94,  84,  24, &
    &             70, 108,  20,  19, 110, 115,  43,  82,  59,  83,  63,   5,  49, 107, &
    &              2,  10, 139, 140,  53,  73, 135,  99, 120,  58, 121,  65,  27,  51, &
    &             11, 122, 130,   7, 111,  60,  78,  80, 118,  21,  29,  48,  85,  25, &
    &             69,   9, 112,  31,   6, 109, 101,  81,  40, 126,  62,  16,  36,  95, &
    &              4,  96, 131, 124, 125,  54, 127, 100, 128,  39,  75,  28,  86,  52, &
    &             12,  68, 114, 113, 103, 116,  79,  30,  91,  56,  33,  35,  26,  66, &
    &             71,  38,  47,  45,  44, 104, 119,  18,  93,  17, 133,  64,  15,  22, &
    &             23,  67,   8,  61,  46, 102, 105,  42,  77,  32,  57,  34,  50,  14 /), SHAPE = (/14, 10/) )
  
  !
  ! Team choices - Shortwave
  !
  INTEGER, DIMENSION(ngptsw, 1), TARGET :: &
    & sw_teams_all = RESHAPE( (/ (i, i = 1, ngptsw) /), SHAPE = (/ ngptsw, 1 /) )
  INTEGER, DIMENSION(1, ngptsw), TARGET :: &
    & sw_teams_one = RESHAPE( (/ (i, i = 1, ngptsw) /), SHAPE = (/ 1, ngptsw /) )
  INTEGER, DIMENSION(7, 16), TARGET :: &
    & sw_7 = RESHAPE(                  &
    &        (/  57,  14,  65,  92, 110,  96,  39, &
    &            58,  25,  44,  64, 109,   4,  72, &
    &            77,   6,  33,  99, 106, 111,  90, &
    &            78,  41,  40,  74, 108, 112, 104, &
    &            79,  13,  53,  94, 107,  91,  49, &
    &            67,  26,  42,  73, 100,  95,  63, &
    &            68,  12,  43,  93, 105,  50,   7, &
    &            69,   5,  51,  97,  80,  31,  62, &
    &            76,  11,  54,  98,  48,  30, 103, &
    &            75,   9,  52,   8,  71,  55,   3, &
    &            59,  15,  66,  89,  19,  84, 101, &
    &            83,  10,  24,  87,  20,  61,  32, &
    &            60,  17,  28,  37,  45,  21,  56, &
    &            81,  18,  38,  36,  46,  29,   2, &
    &            70,  34,  27,  88,  85,   1, 102, &
    &            82,  16,  23,  35,  47,  22,  86 /), SHAPE = (/ 7, 16 /) ) 
  INTEGER, DIMENSION(8, 14), TARGET :: &
    & sw_8 = RESHAPE(                  &
    &        (/  57,  13,  52,  92, 100, 112,  90,  32, &
    &            58,  15,  41,  64, 108,  96,   4,  72, &
    &            77,   5,  40,  97,  98, 110,  95, 104, &
    &            78,  17,  33,  66, 105, 111,  39,   7, &
    &            79,  25,  44,  51,  99,  91,  63, 101, &
    &            67,  14,  54,  94, 107, 109,  49, 102, &
    &            68,  12,  53,  73, 106,  50,  80,  82, &
    &            69,  18,  34,  74,  23,  71,  55,  31, &
    &            76,   9,  42,  93,  48,  61,  28,  86, &
    &            75,   6,  65,   8,  88,  45,  62,   2, &
    &            59,  10,  43,  89,  85,  19,  21,   3, &
    &            83,  26,  24,  35,  47,   1,  22,  56, &
    &            60,  11,  27,  36,  46,  20,  30,  70, &
    &            81,  16,  38,  87,  37,  84,  29, 103 /), SHAPE = (/ 8, 14 /) ) 
  INTEGER, DIMENSION(16, 7), TARGET :: &
    & sw_16 = RESHAPE(                 &
    &         (/  57,  10,  15,  34,  65,  93,  98, 108, 111, 112,  24,  28,  30,  83,  39,   4, &
    &             58,   9,  26,  33,  53,  94, 100, 107, 110,   7,  23,  55,  31, 103,  59,  95, &
    &             77,  13,  16,  44,  54,  97,  99, 106,  38,  19,  84,  62,   2,  89,  75,  96, &
    &             78,   5,  25,  43,  52,  92, 105, 109,  48,  20,  45,  22,  70,   3,  32,   8, &
    &             79,   6,  18,  42,  73,  64,  63,  35,  36,  61,  21,  27,  81,  86,  72,  90, &
    &             67,  12,  14,  51,  40,  74,  49,  87,  71,  46,   1,  82,  91, 102, 104,  76, &
    &             68,  11,  17,  41,  66,  50,  37,  47,  88,  85,  80,  29,  60, 101,  56,  69 /), SHAPE = (/ 16, 7 /) ) 
  
  !
  ! Encapsulate the strategy 
  !
  TYPE spec_sampling_strategy
    PRIVATE
    INTEGER, DIMENSION(:, :), POINTER :: teams => NULL()
    INTEGER :: num_gpts_ts            ! How many g points at each time step
  END TYPE spec_sampling_strategy
  
  PUBLIC :: spec_sampling_strategy, &
          & set_spec_sampling_lw, set_spec_sampling_sw, get_num_gpoints, get_gpoint_set
CONTAINS
  ! -----------------------------------------------------------------------------------------------
  !>
  !! @brief Sets a spectral sampling strategy
  !! 
  !! @remarks: Choose a set of g-point teams to use. 
  !!   Two end-member choices: 
  !!   strategy = 1 : a single team comprising all g-points, i.e. broadband integration
  !!   strategy = 2 : ngpts teams of a single g-point each, i.e. a single randomly chosen g-point
  !!     This can be modified to choose m samples at each time step (with or without replacement, eventually) 
  !!   Other strategies must combine n teams of m gpoints each such that m * n = ngpts
  !!   strategy 1 (broadband) is the default
  !!
  !
  FUNCTION set_spec_sampling_lw(strategy, num_gpts_ts) RESULT (this) 
    INTEGER, INTENT(IN) :: strategy
    INTEGER, OPTIONAL, &
             INTENT(IN) :: num_gpts_ts
    TYPE(spec_sampling_strategy) :: this
	!
	! Default is broadband integration 
	!
    this%teams => lw_teams_all
    SELECT CASE(strategy) 
      CASE (1) 
        !
        ! All gpoints, i.e. full broadband integration
        !
        this%teams => lw_teams_all 
      CASE (2) 
        !
        ! Randomly chosen g-points ("MCSI: Monte Carlo spectral integration") 
        !   Can choose more than one g-point at a time; 
        !   Can choose to sample with or without replacement (without is slower) 
        !
        this%teams => lw_teams_one
      !
      ! Teams - only some values are possible 
      !
      CASE ( 4) 
        this%teams => lw_4
      CASE ( 5) 
        this%teams => lw_5
      CASE ( 7) 
        this%teams => lw_7
      CASE (10) 
        this%teams => lw_10
      CASE (14) 
        this%teams => lw_14
      CASE DEFAULT
		WRITE (message_text, '("Trying to set LW sampling strategy to unknown value ", i2)') &
			   strategy
		CALL message('',message_text, level = em_warn)
    END SELECT
    
    this%num_gpts_ts = SIZE(this%teams, DIM = 1)    
    !
    ! Sampling options for MCSI, ignored otherwise
    !
    IF (strategy == 2) THEN
      IF(PRESENT(num_gpts_ts)) this%num_gpts_ts = num_gpts_ts
    END IF 

  END FUNCTION set_spec_sampling_lw
  ! -----------------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------------
  !>
  !! @brief Sets a spectral sampling strategy
  !! 
  !! @remarks: Choose a set of g-point teams to use. 
  !!   Two end-member choices: 
  !!   strategy = 1 : a single team comprising all g-points, i.e. broadband integration
  !!   strategy = 2 : ngpts teams of a single g-point each, i.e. a single randomly chosen g-point
  !!     This can be modified to choose m samples at each time step (with or without replacement, eventually) 
  !!   Other strategies must combine n teams of m gpoints each such that m * n = ngpts
  !!   strategy 1 (broadband) is the default
  !!
  !
  FUNCTION set_spec_sampling_sw(strategy, num_gpts_ts) RESULT (this) 
    INTEGER, INTENT(IN) :: strategy
    INTEGER, OPTIONAL, &
             INTENT(IN) :: num_gpts_ts
    TYPE(spec_sampling_strategy) :: this
    
    
	!
	! Default is broadband integration 
	!
    this%teams => sw_teams_all
    SELECT CASE(strategy) 
      CASE (1) 
        !
        ! All gpoints, i.e. full broadband integration
        !
        this%teams => sw_teams_all 
      CASE (2) 
        !
        ! Randomly chosen g-points ("MCSI: Monte Carlo spectral integration") 
        !   Can choose more than one g-point at a time; 
        !   Can choose to sample with or without replacement (without is slower) 
        !
        this%teams => sw_teams_one
      !
      ! Teams - only some values are possible 
      !
      CASE ( 7) 
        this%teams => sw_7
      CASE ( 8) 
        this%teams => sw_8
      CASE (16) 
        this%teams => sw_16
      CASE DEFAULT
		WRITE (message_text, '("Trying to set sw sampling strategy to unknown value ", i2)') &
			   strategy
		CALL message('',message_text, level = em_warn)
    END SELECT
    
    this%num_gpts_ts = SIZE(this%teams, DIM = 1) 
    
    !
    ! Sampling options for MCSI, ignored otherwise
    !
    IF (strategy == 2) THEN
      IF(PRESENT(num_gpts_ts)) this%num_gpts_ts = num_gpts_ts
    END IF 

  END FUNCTION set_spec_sampling_sw
  ! -----------------------------------------------------------------------------------------------
  !>
  !! @brief Returns the number of g-points to compute at each time step 
  !! 
  INTEGER FUNCTION get_num_gpoints(strategy) 
    TYPE (spec_sampling_strategy), INTENT(IN) :: strategy
    
    get_num_gpoints = strategy%num_gpts_ts
  END FUNCTION get_num_gpoints
  ! -----------------------------------------------------------------------------------------------
  !>
  !! @brief Returns one set of g-points consistent with sampling strategy
  !! 
  FUNCTION get_gpoint_set(kproma, kbdim, strategy, seeds) 
    INTEGER,                      INTENT(IN)    :: kproma, kbdim
    TYPE(spec_sampling_strategy), INTENT(IN)    :: strategy
    INTEGER,                      INTENT(INOUT) :: seeds(:,:) ! dimensions kbdim, rng seed_size
    INTEGER, DIMENSION(kproma, strategy%num_gpts_ts) :: get_gpoint_set
    
    REAL(WP):: rn(kbdim)
    INTEGER :: team(kbdim) 
    INTEGER :: num_teams, num_gpts_team, it, jl
    ! --------
    
    num_teams      = SIZE(strategy%teams, 2)
    num_gpts_team  = SIZE(strategy%teams, 1)
    
    IF(num_teams == 1) THEN 
      !
      ! Broadband integration
      !
      get_gpoint_set(1:kproma,:)  = SPREAD(strategy%teams(:, 1), DIM = 1, NCOPIES = kproma) 
    ELSE IF(num_gpts_team > 1) THEN 
      !
      ! Mutiple g-points per team, including broadband integration
      !   Return just one team
      !
      CALL get_random(kproma, kbdim, seeds, rn) 
      team(1:kproma) = MIN(INT(rn(1:kproma) * num_teams) + 1, num_teams)
      DO jl = 1, kproma
        get_gpoint_set(jl, :) = strategy%teams(:,team(jl))
      END DO  
    ELSE
      !
      ! MCSI - return one or more individual points chosen randomly
      !   Need to add option for sampling without replacement  
      !
      DO it = 1, strategy%num_gpts_ts
        CALL get_random(kproma, kbdim, seeds, rn) 
        team(1:kproma) = MIN(INT(rn(1:kproma) * num_teams) + 1, num_teams)
        get_gpoint_set(1:kproma, it) = strategy%teams(1, team(1:kproma)) 
      END DO 
    END IF 
    
  END FUNCTION get_gpoint_set
  ! -----------------------------------------------------------------------------------------------
END MODULE mo_spec_sampling 
