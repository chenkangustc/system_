!$
!===================================================================================================
!
!   type for global scope parameter
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          SteadyState
!                               TransientState
!
!===================================================================================================
module state_header
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    implicit none 
    private
    public  :: SteadyState, TransientState
        
    ! --------------------------------------------------------------------------
    ! global parameters for steady calculation
    
    ! type for steady state of the problem
    type, private  :: state_steady_tp
        integer, public           :: ng                                         ! number of energy group
        integer, public           :: point                                      ! number of radial nodal point
        integer, public           :: nodal                                      ! number of radial nodal
        integer, public           :: layer                                      ! number of axial layer
        integer, public           :: segment                                    ! number of boudary segment
        integer, public           :: zone                   = 1                 ! number of zone per layer
        integer, public           :: mat                    = 1                 ! number of material used
        integer, public           :: source                 = 0                 ! number of extenal source kind
        integer, public           :: sn                     = 4                 ! number of the polar division (N in term 'SN')
        integer, public           :: scat_order             = 0                 ! number of anisotropic scatter order
        integer, public           :: layer_top              = 0                 ! number of upper reflector layers
        integer, public           :: layer_bottom           = 0                 ! number of bottom reflector layers
        integer, public           :: upscatter              = 0                 ! number of energy group with upscatter
    end type state_steady_tp
    
    ! type for deduced steady parameter 
    type, private  :: deduced_state_steady_tp
        integer, public           :: scat_xs                                    ! number of scatter matrix
        integer, public           :: nodal_total                                ! number of total nodal for the core
        integer, public           :: direction                                  ! number of total direction
    end type deduced_state_steady_tp
    
    ! type for steady flag
    type, private  :: flag_steady_tp
        logical, public           :: is_eigen               = .TRUE.            ! is eigenvalue problem ?    
        logical, public           :: is_bevel_edge          = .FALSE.           ! is bevel edge ?
        logical, public           :: is_60degree            = .FALSE.           ! is 60 degree symmetry ?
        integer, public           :: n_theta                = 0                 ! degree for the symmetry edge (n*60, n=1,5)
        logical, public           :: is_square              = .FALSE.           ! is the core only consisit of square assembly ?
        logical, public           :: is_hexagonal           = .FALSE.           ! is the core only consisit of hexgonal assembly ?
        logical, public                     :: is_link      = .FALSE.           ! is use Link module ?
        character(len=MAX_WORD_LEN), public :: link_type    = 'fitting'         ! Link module type ?
        character(len=MAX_WORD_LEN), public :: case_title   = 'sample'          ! title name of the problem
        real(KREAL), public             :: rated_power  = 1.0                   ! rated power level, in watt
        real(KREAL), public             :: power_frac   = 1.0                   ! initial power percentage, in 1.0
        real(KREAL), public             :: power_level  = 1.0                   ! actual power level, in Watt
    end type flag_steady_tp
    
    ! type for steady method selection
    type, private  :: method_steady_tp
        logical, public           :: is_LW                  = .FALSE.           ! is LW extrapolate ?
        logical, public           :: is_upscatter_cycle     = .FALSE.           ! is perform upscatter cycle ?
        integer, public           :: n_upscatter_cycle      = 3                 ! number of upscatter cycle per out iteration
        logical, public           :: is_Ks                  = .FALSE.           ! is Ks method to FSP problem
    end type method_steady_tp
    
    ! type for output selection
    type, private  :: output_steady_tp
        logical, public             :: is_HDF5              = .FALSE.           ! is output HDF5 format distribution parameter
        logical, public             :: is_vtk               = .FALSE.           ! is output vtk file for visualization
        logical, public             :: is_log               = .TRUE.            ! is output iteration log 
    end type output_steady_tp
    
    ! type for feedback model
    type, private  :: feedback_steady_tp
        logical, public                      :: is_feedback = .FALSE.           ! is feedback or not ?
        logical, public                      :: is_nested   = .FALSE.           ! feedback perform one time or nested cycle ?
        logical, public                      :: is_inner    = .FALSE.           ! is inner feedback ?
        logical, public                      :: is_coupled  = .FALSE.           ! is coupled feedback ? (by real calculation)
        logical, public                      :: is_model    = .FALSE.           ! is model feedback, for benchamrk verification
        character(len=MAX_WORD_LEN), public  :: model_name  = 'lra'             ! feedback model used
    end type feedback_steady_tp
 
    ! type for steady miscellaneous selction
    ! @TODO use for extention
    type, private :: miscellaneous_steady_tp
        logical, public                      :: is_external  = .FALSE.          ! is xsec come from external file ? 
        logical, public                      :: is_import    = .FALSE.          ! is xsec come from burnup output ?
        character(len=MAX_WORD_LEN), public  :: file         = 'to_kinetics.h5' ! burnup output filename 
        integer, public                      :: burnstep     = 0                ! burn point in output file
        integer, public                      :: nthread      = 0                ! number of threads, '0' means closing it 
        logical, public                      :: is_angwise   = .TRUE.
    end type miscellaneous_steady_tp
    
    ! --------------------------------------------------------------------------
    ! global parameters for transient calculation
        
    ! type for state of the transient problem
    type, private  :: state_transient_tp
        integer, public             :: dg                   = 6                 ! number of delay neutron precusor group
        integer, public             :: hg                   = 6                 ! number of decay heat precursor group
    end type state_transient_tp
    
    ! type of the problem for transient
    type, private  :: flag_transient_tp
        logical, public                      :: is_perturb      = .FALSE.       ! perform perturbation calculation ?
        logical, public                      :: is_adjoint      = .FALSE.       ! perform adjoint calculation ?
        character(len=MAX_WORD_LEN), public  :: adjoint_type    = 'homogeneous' ! select adjoint type
        logical, public                      :: is_transient    = .FALSE.       ! perform transint calculation ?
        logical, public                      :: is_pkmodel      = .FALSE.       ! perform point kinetics analysis ? 
        character(len=MAX_WORD_LEN), public  :: kinetics_type   = 'normal'      ! kinetics parameter type
        logical, public                      :: is_CR_rod       = .FALSE.       ! is CR_rod in the core ?
        logical, public                      :: is_decusping    = .FALSE.       ! perform control rod decusping ?
        logical, public                      :: is_xesm         = .FALSE.       ! consider Xe/Sm in transient ?
        logical, public                      :: is_decay_heat   = .FALSE.       ! consider decay heat in transient ?
        logical, public                      :: is_boron_search = .FALSE.       ! perform critical boron search ?
    end type flag_transient_tp
    
    ! type for transient method selection
    type, private :: method_transient_tp
        character(len=MAX_WORD_LEN), public  :: scheme                  = 'theta'                   ! transient method, select THETA ?
        logical, public                      :: is_extrapolation        = .TRUE.                    ! perform time step extrapolation ?
        character(len=MAX_WORD_LEN), public  :: extrapolation_scheme    = 'linear'                  ! scheme of extrapolation for time flux
    end type method_transient_tp
    
    ! type for perturbation type and details
    type, private  :: perturbation_tp
        logical, public             :: is_xsec              = .FALSE.           ! is due to direct cross section change ?
        logical, public             :: is_source            = .FALSE.           ! is due external neutron source ?
        logical, public             :: is_CR_move           = .FALSE.           ! is due to control rod remove ?
        logical, public             :: is_mat               = .FALSE.           ! is due to material change ?
        logical, public             :: is_flow              = .FALSE.           ! is due to inlet flow change ?
        logical, public             :: is_Tm                = .FALSE.           ! is due to inlet Tm change ? 
    end type perturbation_tp
        
    ! type for transient miscellaneous selection
    ! @TODO use for extention
    type, private :: miscellaneous_transient_tp
        integer, public             :: nothing
    end type miscellaneous_transient_tp

    ! --------------------------------------------------------------------------
    ! global objects for steady
    type  :: SteadyState
        type(state_steady_tp), public             :: state                      ! problem state
        type(deduced_state_steady_tp), public     :: deduce                     ! deduced parameters
        type(flag_steady_tp), public              :: flag                       ! problem flag
        type(method_steady_tp), public            :: method                     ! method selection
        type(output_steady_tp), public            :: output                     ! output option
        type(feedback_steady_tp), public          :: feedback                   ! feedback option
        type(miscellaneous_steady_tp), public     :: misc                       ! miscellaneous option
    end type SteadyState                    
                                          
    ! global objects for transient        
    type  :: TransientState                 
        type(state_transient_tp), public          :: state                      ! problem state
        type(flag_transient_tp), public           :: flag                       ! problem flag
        type(method_transient_tp), public         :: method                     ! method selection
        type(perturbation_tp), public             :: perturb                    ! perturbation option
        type(miscellaneous_transient_tp), public  :: misc                       ! miscellaneous option
    end type TransientState
    
end module state_header
