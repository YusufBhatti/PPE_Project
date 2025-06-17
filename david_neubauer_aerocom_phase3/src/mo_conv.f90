MODULE mo_conv

  ! U. Lohmann, ETHZ, 2007-01-24, additional cloud output

  USE mo_kind,          ONLY: dp
  USE mo_linked_list,   ONLY: t_stream

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_stream_conv

  TYPE (t_stream), PUBLIC, POINTER :: conv

  REAL(dp), PUBLIC, POINTER :: na_cv(:,:,:)
  REAL(dp), PUBLIC, POINTER :: cdncact_cv(:,:,:)
  REAL(dp), PUBLIC, POINTER :: twc_conv(:,:,:)
  REAL(dp), PUBLIC, POINTER :: conv_time(:,:,:)

CONTAINS

  SUBROUTINE construct_stream_conv

    ! construct_stream_conv: allocates output streams
    !                        for the activation schemes
    !
    ! Author:
    ! -------
    ! Philip Stier, 2004
    !

    USE mo_memory_base,   ONLY: new_stream, add_stream_element, AUTO,  &
                                default_stream_setting, add_stream_reference
    USE mo_linked_list,   ONLY: HYBRID
    USE mo_filename,      ONLY: out_filetype
 
    IMPLICIT NONE

    !--- Create new stream:

    CALL new_stream (conv ,'conv', filetype=out_filetype)

    ! add standard fields for post-processing:
    
    CALL add_stream_reference (conv, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
    CALL add_stream_reference (conv, 'lsp'     ,'sp'    ,lpost=.TRUE.)
    CALL add_stream_reference (conv, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
    CALL add_stream_reference (conv, 'gboxarea','geoloc',lpost=.TRUE.)
    
    CALL default_stream_setting (conv,                &
                                 lpost     = .TRUE. , &
                                 laccu     = .TRUE. , &         !>>dod<<
                                 lrerun    = .TRUE. , &
                                 leveltype = HYBRID , &
                                 table     = 199,     &
                                 code      = AUTO     )

    ! cloud Properties:
    
    CALL add_stream_element (conv,   'TWC_CONV',    twc_conv,                 &
         longname='LWC+IWC from detr.+ cloud weighted',  units='kg m-3')
    
    CALL add_stream_element (conv,   'CONV_TIME', conv_time,                  &
         longname='acc. cloud occ. conv time fraction', units='1'     )
    
    CALL add_stream_element (conv,   'CDNCACT_CV', cdncact_cv, lpost=.false., laccu=.false.,  &
         longname='activated CDNC in conv', units='m-3'      )
    
    CALL add_stream_element (conv,   'NA_CV', na_cv,  lpost=.false., laccu=.false.,  &
         longname='aerosols for nucleation in conv', units='m-3'      )
    
  END SUBROUTINE construct_stream_conv
  
END MODULE mo_conv

