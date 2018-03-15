!$
!===================================================================================================
!
!   this module is for input file's key word and pre-defined input parameter
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    Set_section_keyword
!                               Set_card_keyword
!                               Set_option_keyword
!                               Is_keyword
!
!   Public type lists:          No
!
!===================================================================================================
module input_keyword
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none
    
    ! Note: parameter in this module is public by default
    public
    
    ! --------------------------------------------------------------------------
    ! key word for input parameter
    integer, parameter, private  :: N_SECTION        = 20                       ! max section number of input file
    integer, parameter, private  :: N_KEYWORD        = 50                       ! max keyword number per section
    integer, parameter, private  :: N_OPTION         = 10                       ! max number of option for input parameter selection
    
    ! for section name
    character(len=MAX_WORD_LEN)  :: INP_SECTION(N_SECTION)    = CHAR_SENTINEL

    ! --------------------------------------------------------------------------
    ! key word for per section
    character(len=MAX_WORD_LEN)  :: CONTROL_CARD(N_KEYWORD)      = CHAR_SENTINEL                    ! card key word for control section
    character(len=MAX_WORD_LEN)  :: METHOD_CARD(N_KEYWORD)       = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: LINK_CARD(N_KEYWORD)         = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: MATERIAL_CARD(N_KEYWORD)     = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: GEOMETRY_CARD(N_KEYWORD)     = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: TRANSIENT_CARD(N_KEYWORD)    = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: FEEDBACK_CARD(N_KEYWORD)     = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: PERTURBATION_CARD(N_KEYWORD) = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: PKMODEL_CARD(N_KEYWORD)      = CHAR_SENTINEL
    
    ! --------------------------------------------------------------------------
    ! use for character option input checking 
    character(len=MAX_WORD_LEN)  :: OPTION_ADJOINT_TYPE(N_OPTION)               = CHAR_SENTINEL     ! option to select
    character(len=MAX_WORD_LEN)  :: OPTION_KINETICS_TYPE(N_OPTION)              = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: OPTION_LINK_TYPE(N_OPTION)                  = CHAR_SENTINEL
    
    character(len=MAX_WORD_LEN)  :: OPTION_TRANSIENT_SCHEME(N_OPTION)           = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: OPTION_TRANSIENT_EXTRAPOLATION(N_OPTION)    = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: OPTION_FEEDBACK_PARAMETER(N_OPTION)         = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: OPTION_BOUNDARY_CONDITION(N_OPTION)         = CHAR_SENTINEL
    
    character(len=MAX_WORD_LEN)  :: OPTION_XSEC_TYPE(N_OPTION)                  = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: OPTION_PERTURB_TYPE(N_OPTION)               = CHAR_SENTINEL
    character(len=MAX_WORD_LEN)  :: OPTION_FEEDBACK_MODEL(N_OPTION)             = CHAR_SENTINEL
                
contains
    !$
    !===============================================================================================
    ! set section key word
    !===============================================================================================
    subroutine Set_section_keyword ()

        INP_SECTION(1:11) =  ['CASENAME:            ',    &                     ! for case name
                          &   'CONTROL:             ',    &                     ! for control logical
                          &   'METHOD:              ',    &                     ! for calcultaion scheme
                          &   'LINK:                ',    &                     ! for cross section link
                          &   'MATERIAL:            ',    &                     ! for cross section and kinetics parameter
                          &   'GEOMETRY:            ',    &                     ! for geometry
                          &   'TRANSIENT:           ',    &                     ! for transient process
                          &   'FEEDBACK:            ',    &                     ! for feedback
                          &   'PERTURBATION:        ',    &                     ! for perturbation
                          &   'PKMODEL:             ',    &                     ! for point-kinetics 
                          &   'END:                 ']                          ! for end of input
    
    end subroutine Set_section_keyword

    !$
    !===============================================================================================
    ! set card key word
    !===============================================================================================
    subroutine Set_card_keyword ()
    
        CONTROL_CARD(1:17) =     ['PK_MODEL             ',      &
                          &       'STEADY               ',      &
                          &       'PPM_SEARCH           ',      &
                          &       'ROD_SEARCH           ',      &
                          &       'ADJOINT              ',      &
                          &       'TRANSIENT            ',      &
                          &       'LINK                 ',      &
                          &       'FEEDBACK             ',      &
                          &       'PERTURBATION         ',      &
                          &       'DECAY_HEAT           ',      &
                          &       'ENERGY_GROUP         ',      &
                          &       'DELAYD_GROUP         ',      &
                          &       'DELAYED_GROUP        ',      &
                          &       'OUTPUT               ',      &
                          &       'THREAD               ',      &
                          &       'INITVAL              ',      &
                          &       'MISC                 '] 
                
        METHOD_CARD(1: 7) =      ['QUDARATURE           ',      &
                          &       'QUADRATURE           ',      &
                          &       'STEADY               ',      &
                          &       'TRANSIENT            ',      &
                          &       'ERROR_TYPE           ',      &
                          &       'ERROR_EIGEN          ',      &
                          &       'ERROR_FSP            '] 
                          
        LINK_CARD(1: 5) =        ['DATA                 ',      &
                          &       'PARAMETER            ',      &
                          &       'MAX_LIMIT            ',      &
                          &       'REFERENCE            ',      &
                          &       'MIN_LIMIT            '] 
                
        MATERIAL_CARD(1: 9) =    ['SCALE                ',      &
                          &       'IMPORT               ',      &
                          &       'TYPE                 ',      &
                          &       'MAT                  ',      &
                          &       'MAT_D                ',      &
                          &       'Q_EXTERNAL           ',      &
                          &       'DETINFO              ',      &
                          &       'DETXS                ',      &
                          &       'DETPOINT             '] 
                
        GEOMETRY_CARD(1:32) =    ['MESH                 ',      &               ! 
                          &       'TYPE                 ',      &               ! ansys-mesh
                          &       'ANSYS_FILE           ',      &
                          &       'XY_EXPAND            ',      &
                          &       'SCALE                ',      &
                          &       'BC_RADIAL            ',      &
                          &       'REC_XDIM             ',      &               ! rec-mesh
                          &       'REC_YDIM             ',      &
                          &       'REC_CONF             ',      &
                          &       'REC_BC               ',      &
                          &       'REC_MESH             ',      &
                          &       'HEX_DIM              ',      &               ! hex-mesh
                          &       'HEX_CONF             ',      &
                          &       'HEX_BC               ',      &
                          &       'HEX_MESH             ',      &
                          &       'BC_AXIAL             ',      &               ! configure
                          &       'LAYER                ',      &
                          &       'REF_AXIAL            ',      &
                          &       'HEIGHT               ',      &
                          &       'PLANE_MAT            ',      &
                          &       'ASSIGN_MAT           ',      &
                          &       'PLANE_SOURCE         ',      &
                          &       'ASSIGN_SOURCE        ',      &
                          &       'FA_TYPE              ',      &
                          &       'ASSIGN_FA            ',      &
                          &       'EXTQ_TYPE            ',      &
                          &       'ASSIGN_EXTQ          ',      &
                          &       'PRINT_MASK           ',      &
                          &       'CR_STEP              ',      &               ! CR bank
                          &       'CR_BANK              ',      &
                          &       'CR_CONFIGURE         ',      &
                          &       'CR_POSITION          ']
            
        TRANSIENT_CARD(1:18) =   ['STEP_TR              ',      &
                          &       'SECTION              ',      &
                          &       'INTERVAL             ',      &
                          &       'PERTURB_XSEC         ',      &
                          &       'SIGMA_S              ',      &
                          &       'SIGMA_T              ',      &
                          &       'SIGMA_F_NU           ',      &
                          &       'PERTURB_SOURCE       ',      &
                          &       'INTENSITY            ',      &
                          &       'PERTURB_CR           ',      &
                          &       'CR_MOVE              ',      &
                          &       'CR_TRIP              ',      &
                          &       'PERTURB_MAT          ',      &
                          &       'MAT                  ',      &
                          &       'PERTURB_FLOW         ',      &
                          &       'FLOW                 ',      &
                          &       'PERTURB_TM           ',      &
                          &       'TM                   ']
            
        FEEDBACK_CARD(1: 15) =   ['FDBK                 ',      &
                          &       'RELAXATION           ',      &
                          &       'EFFTF                ',      &
                          &       'DESIGN               ',      &
                          &       'FLOW                 ',      &
                          &       'SEARCH               ',      &
                          &       'GAMMA                ',      &
                          &       'FQ_LATTICE           ',      &
                          &       'CRITERIA             ',      &
                          &       'MESH_SIZE            ',      &
                          &       'COOLANT              ',      &
                          &       'GEOM_INFO            ',      &
                          &       'ASSEMBLY             ',      &
                          &       'PROPERTY_TYPE        ',      &
                          &       'TH_CONFIGURE         ']
        
        PERTURBATION_CARD(1: 7) =['PERTURB              ',      &
                          &       'PERTURB_INFO         ',      &
                          &       'PLANE_MAT            ',      &
                          &       'ASSIGN_MAT           ',      &
                          &       'MAT                  ',      &
                          &       'CR_WORTH             ',      &
                          &       'CR_POS               ']
        
        PKMODEL_CARD(1: 4)      =['REACTIVITY           ',      &
                          &       'COEFFICIENT          ',      &
                          &       'BLOCK                ',      &
                          &       'PKPARAMETER          ']
        
    end subroutine Set_card_keyword

    !$
    !===============================================================================================
    ! set option key word
    !===============================================================================================
    subroutine Set_option_keyword ()
    
        OPTION_ADJOINT_TYPE(1: 4) =       ['HOMOGENEOUS     ',      &
                              &            'ONE             ',      &
                              &            'FISSION_F_NU    ',      &
                              &            'FISSION_F_KAPPA ']
        OPTION_KINETICS_TYPE(1: 2) =      ['NORMAL          ',      &
                              &            'SRAC            ']
        OPTION_LINK_TYPE(1: 2) =          ['FITTING         ',      &
                              &            'INTERPOLATION   ']
                              
                                    
        OPTION_TRANSIENT_SCHEME(1: 4) =   ['THETA           ',      &
                              &            'SCM             ',      &
                              &            'PCQS            ',      &
                              &            'PK              ']
        OPTION_TRANSIENT_EXTRAPOLATION(1: 2) = ['LINEAR     ',      &
                              &                 'LOG        ']
        OPTION_FEEDBACK_PARAMETER(1: 5) = ['BU              ',      &
                              &            'CB              ',      &
                              &            'TF              ',      &
                              &            'TM              ',      &
                              &            'RHO_M           ']
        OPTION_BOUNDARY_CONDITION(1: 2) = ['VACCUM          ',      &
                              &            'REFLECT         ']
                              
        OPTION_XSEC_TYPE(1: 3) =          ['SIGMA_S         ',      &
                              &            'SIGMA_T         ',      &
                              &            'SIGMA_F_NU      ']
                                    
        OPTION_PERTURB_TYPE(1: 2) =       ['STEP            ',      &
                              &            'RAMP            ']
                                    
        OPTION_FEEDBACK_MODEL(1: 2) =     ['LAR             ',      &
                              &            'NEACRP          ']
                                    
    end subroutine Set_option_keyword
    
    !$
    !===============================================================================================
    ! require 'key' is in 'list' ?
    !===============================================================================================
    function Is_keyword(list, key)  result(is_true)
        
        character(MAX_WORD_LEN), intent(in)  :: list(:)
        character(MAX_WORD_LEN), intent(in)  :: key
        
        logical  :: is_true
        integer  :: i
        
        is_true = .FALSE.
        
        do i = 1, SIZE(list)
            if (TRIM(ADJUSTL(list(i))) == TRIM(ADJUSTL(key)))  then
                is_true = .TRUE.
                exit
            end if
        end do
        
    end function Is_keyword
                
end module input_keyword
