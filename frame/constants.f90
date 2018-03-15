!$
!===================================================================================================
!
!   define some general constants for the whole program scope
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          No
!
!===================================================================================================
module constants

    implicit none
    public
    
    ! --------------------------------------------------------------------------
    ! kind number for integer, real and character 
    integer, parameter  :: KINT            = SELECTED_INT_KIND(r=17)            ! high precision integer kind
    integer, parameter  :: KREAL           = SELECTED_REAL_KIND(p=14)           ! normal percision real kind, (6/14)
    integer, parameter  :: KDOUBLE         = SELECTED_REAL_KIND(p=14)           ! high percision real kind, this should larger than KREAL
    integer, parameter  :: MAX_LINE_LEN  = 800                                  ! max length of a line, as group increasing
    integer, parameter  :: MAX_WORD_LEN  = 800                                  ! max length of a word
    integer, parameter  :: MAX_WORDS     = 800                                  ! max number of word per line
    
    real(KREAL), parameter  :: INT_DEFAULT_BYTE = KIND(0) / (1024.0 * 1024.0)   ! size in 'MB' for default integer kind
    real(KREAL), parameter  :: REAL_DEFAULT_BYTE = KIND(0.0) / (1024.0 * 1024.0)! size in 'MB' for default real kind
    real(KREAL), parameter  :: INT_BYTE = KIND(0_KINT) / (1024.0 * 1024.0)      ! size in 'MB' for self-define integer kind
    real(KREAL), parameter  :: REAL_BYTE = KIND(0.0_KREAL) / (1024.0 * 1024.0)  ! size in 'MB' for self-define real kind

    ! --------------------------------------------------------------------------
    ! comparison criteria
    real(KREAL), parameter  :: EPS_LOWER      = 1.0D-03                         ! loer precision criteria
    real(KREAL), parameter  :: EPS_LOW        = 1.0D-06                         ! low precision criteria
    real(KREAL), parameter  :: EPS_HIGH       = 1.0D-09                         ! high precision criteria
    real(KREAL), parameter  :: EPS_HIGHER     = 1.0D-12                         ! higher precision criteria
                                                                                
    real(KREAL), parameter  :: EPS_EQUAL      = 1.0D-12                         ! two real number equal limitation, relative
    real(KREAL), parameter  :: EPS_ZERO       = 1.0D-12                         ! limit to set as positive zero
                                                                                
    real(KREAL), parameter  :: EQUAL_MARGIN   = 1.0D-12                         ! two real number equal limitation, replace the old 'EQUATION'

    ! --------------------------------------------------------------------------    
    ! mathmatic constants
    real(KREAL), parameter  :: PI              =  3.1415926535898D0             ! pi
    real(KREAL), parameter  :: HALF            =  0.5                           ! half of something
    real(KREAL), parameter  :: SQRT3           =  SQRT(3.0)                     ! square root of real number 3
    integer, parameter      :: REVERSE_ORDER   = -1                             ! reversed do cycle order
    
    real(KREAL), parameter  :: SIN30           = SIN(30.0D0*PI/180.0D0)
    real(KREAL), parameter  :: COS30           = COS(30.0D0*PI/180.0D0)
    real(KREAL), parameter  :: SIN60           = SIN(60.0D0*PI/180.0D0)
    real(KREAL), parameter  :: COS60           = COS(60.0D0*PI/180.0D0)
    
    ! --------------------------------------------------------------------------    
    ! physics constants
    real(KREAL), parameter  :: AVOGADRO   = 0.602214129                         ! Avogadro's number (10^24/mol)
    real(KREAL), parameter  :: BOLTZMANN  = 8.6173324D-11                       ! Boltzmann constant (Mev/K)
    real(KREAL), parameter  :: CKELVIN    = 273.15                              ! temperature transfer between C and K

    ! --------------------------------------------------------------------------    
    ! zero, one, infinity and null
    integer, parameter            :: INT_ERROR    = -HUGE(0)
    integer, parameter            :: INT_INFINITY =  HUGE(0)                    ! positive infinity in integer
    integer, parameter            :: INT_ZERO     =  0
    integer, parameter            :: INT_ONE      =  1
    
    real(KREAL), parameter        :: REAL_ERROR   = -HUGE(REAL(0.0, KREAL))
    real(KREAL), parameter        :: REAL_INFINITY=  HUGE(REAL(0.0, KREAL))     ! positive infinity in real
    real(KREAL), parameter        :: REAL_ZERO    =  REAL(0.0, KREAL)
    real(KREAL), parameter        :: REAL_ONE     =  REAL(1.0, KREAL)

    logical, parameter            :: LOGI_YES     = .TRUE.                      ! logical true
    logical, parameter            :: LOGI_NO      = .FALSE.                     ! logical false

    character(len=MAX_WORD_LEN), parameter  :: CHAR_NULL           = ''         ! null character 
    character(len=MAX_WORD_LEN), parameter  :: CHAR_SPACE          = ' '        ! space character
    character(len=MAX_WORD_LEN), parameter  :: CHAR_SENTINEL       = '----'     ! sentinel character
    character(len=MAX_WORD_LEN), parameter  :: CHAR_MARK           = '==================================================================================='
    character(len=MAX_WORD_LEN), parameter  :: CHAR_SUBMARK        = '___________________________________________________________________________________'
    character(len=MAX_WORD_LEN), parameter  :: CHAR_SSUBMARK       = '_____________________________________________'

    ! --------------------------------------------------------------------------
    ! warning, error and information kind lists
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_FRAMEWORK = 'SELF FRAMEWORK:'              ! for self-implement general function framework 
    
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_SYSTEM    = 'OPERATION SYSTEM:'            ! for os system
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_FILE      = 'FILE OPERATION:'              ! for file operation: open, read, or something else
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_INPUT     = 'INPUT FORMAT:'                ! for input parameter format

    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_GEOMETRY  = 'GEOMETRY:'                    ! for geometry information
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_XSEC      = 'CROSS SECTION:'               ! for cross section information    
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_TIMESTEP  = 'TRANSIENT STEP:'              ! for transient time step division
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_CRBANK    = 'CONTROL ROD BANK:'            ! for control rod movement
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_PROPERTY  = 'THERMAL PROPERTY:'            ! for thermal property
    character(len=MAX_WORD_LEN), parameter  :: INFO_LIST_FEEDBACK  = 'THERMAL FEEDBACK:'            ! for feedback calculation
    
    ! --------------------------------------------------------------------------
    ! following is used in DAISY code system only
    ! --------------------------------------------------------------------------
    
    ! code information   
    integer, parameter  :: VERSION_MAJOR    = 1
    integer, parameter  :: VERSION_MINOR    = 0
    integer, parameter  :: VERSION_RELEASE  = 0
    character(len=MAX_WORD_LEN), parameter  :: COPYRIGHT  = 'Copyright (c) 2017 Institute of Modern Physics,CAS'
    character(len=MAX_WORD_LEN), parameter  :: LABORATORY = 'Group of Reactor multi-physics Coupling(RMPC)'
    character(len=MAX_WORD_LEN), parameter  :: CODENAME   = 'IMPC-transient'
	character(len=MAX_WORD_LEN), parameter  :: IMP        = 'The Ready,Motivating,Persistent and Collaborative Team'

    
    ! --------------------------------------------------------------------------
    ! file unit management
    type, private  ::  file_tp
        integer, public  :: CASENAME                                            ! input--case define (steady & transient)

        integer, public  :: MAIN                                                ! output--main output files
        integer, public  :: MEMORY                                              ! output--array size in memory
        integer, public  :: TIMELIST                                            ! output--integral parameter for transient process
        integer, public  :: DET                                                 ! output--detector response
        integer, public  :: REACTIVITY                                          ! output--reactivity & kinetics parameter
        integer, public  :: POINTKINETICS                                       ! output--point kinetics value
        integer, public  :: PT                                                  ! output--perturbation calculation for reactivity coefficient
        
        integer, public  :: TH_WARNING                                          ! output--warning information for thermal property
        integer, public  :: TH_RBFD                                             ! output--thermal information for the consequent RBFD analysis
        integer, public  :: TH_HOT                                              ! output--hot value for thermal output
        integer, public  :: TH_AVERAGE                                          ! output--average value for thermal output
    end type file_tp
    
    ! this is a public files unit object
    type(file_tp) :: FILES
    
    ! --------------------------------------------------------------------------
    ! directory
    character(len=MAX_WORD_LEN), parameter  :: DIR_XSEC             = 'xsec/'                       ! directory which contain xsec file from lattics code
    character(len=MAX_WORD_LEN), parameter  :: DIR_LINK             = 'link/'                       ! directory which contain file for lilac

    ! name of input file
    character(len=MAX_WORD_LEN), parameter  :: PREFIX_INPUT         = 'input'                       ! input prefix for file    
    character(len=MAX_WORD_LEN), parameter  :: PREFIX_TMP           = 'tmp'                         ! add to file name when preprocess
    character(len=MAX_WORD_LEN), parameter  :: PREFIX_OUTPUT        = 'output'                      ! output prefix for file
                                                                    
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_XML           = '.xml'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_H5            = '.h5'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_CASENAME      = '.case'
                                                                    
    ! name of output file
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_MAIN          = '.main'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_MEMORY        = '.memory'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_TIMELIST      = '.timelist'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_DET           = '.detector'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_REACTIVITY    = '.reactivity'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_POINTKINETICS = '.pointkinetics'
    character(len=MAX_WORD_LEN), parameter  :: EXTEND_PT            = '.perturbation'
    
    character(len=MAX_WORD_LEN), parameter  :: TH_WARNING           = '.thwarning'
    character(len=MAX_WORD_LEN), parameter  :: TH_HOT               = '.thhotchannel'
    character(len=MAX_WORD_LEN), parameter  :: TH_AVERAGE           = '.thavgchannel'
    character(len=MAX_WORD_LEN), parameter  :: TH_RBFD              = '.RBFD.bin'
                                                                        
end module constants
