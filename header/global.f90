!$
!===================================================================================================
!
!   module for global array parameters
!---------------------------------------------------------------------------------------------------
!   Public subroutine lists:    No
!
!   Public type lists:          No
!
!===================================================================================================
module global
    
    use constants
    use, intrinsic  :: ISO_FORTRAN_ENV
    
    use geometry_header,            only : Meshing, Geometry, Boundary, VTKMeshing
    use material_header,            only : CrossSection, CrossSectionInput, Material,       &
                                        &   KineticsParameter, KineticsParameterInput, ExternalSource
    use quadrature_header,          only : QuadratureSet
    use sweeping_header,            only : Sweeping

    use iteration_header,           only : IterationSource, IterationFlux, IterationCounter, AnisotropicScatterFlux, IterationCriterion
    use initvalue_header,           only : InitialValue
    use adjoint_header,             only : AdjointIteration
    use contain_header,             only : TimeListParameter, DistributionParameter, GroupsFlux
    
    use CRbank_header,              only : ControlRodBank
    use perturbation_header,        only : XsecPerturbation, SourcePerturbation, ControlRodPerturbation, MaterialPerturbation, &
                                        &   CasePerturbation, THPerturbation
    use worth_header,               only : CRWorth
    use transient_header,           only : DelayedNeutron, ShapeFunction, AmplitudeFunction, PointKineticsParameter
    use link_header,                only : LinkParameter
    use feedback_header,            only : FeedbackParameter
    use LRAmodel_header,            only : LRAmodel
    
    use detector_header,            only : Detector
    use precursor_solver_header,    only : PrecursorSolver
    use pksolver_header,            only : PKSolver, PKParameter
    use pkreactivity_header,        only : PKReactivity
                                        
    implicit none
    public
    
    ! --------------------------------------------------------------------------
    ! Note:
    !       dependece: geometry <-- material <-- quadrature <-- sweeping
    !       others may depend on all above
    ! --------------------------------------------------------------------------
    
    ! --------------------------------------------------------------------------
    ! geometry
    type(Meshing)                               :: mesh                         ! ansys generate mesh
    type(Geometry)                              :: geom                         ! scale of geometry, also coordinate
    type(Boundary)                              :: bound                        ! boundary condition
    type(VTKMeshing)                            :: mesh_vtk
    
    ! --------------------------------------------------------------------------
    ! material
    type(Material)                              :: mat_info                     ! material information
                                                
    type(CrossSectionInput)                     :: xsec_inp                     ! cross section for input
    type(CrossSection)                          :: xsec_iter                    ! cross section used in iteration
    type(CrossSection)                          :: xsec_init                    ! for initial time value
    type(CrossSection)                          :: xsec                         ! cross section for current status
                                                
    type(KineticsParameterInput)                :: param_inp                    ! kinetics parameters for input, also for initial time value
    type(KineticsParameter)                     :: param                        ! kinetics parameters for really used in code
                                                
    type(ExternalSource)                        :: Q_ext                        ! external source
    
    ! --------------------------------------------------------------------------
    ! quadrature sets
    type(QuadratureSet)                         :: quad                         ! angular qudarature set
    type(Sweeping)                              :: sweep                        ! meshing sweeeping
    
    ! --------------------------------------------------------------------------
    ! iteration parameter
    type(IterationSource)                       :: iter_q                       ! source in iteration
    type(IterationFlux)                         :: iter_flux                    ! flux in iteration
    type(IterationCounter)                      :: iter_count                   ! count for iteration
    type(InitialValue)                          :: iter_init                    ! iteration initial flux 
    
    type(IterationCriterion), target            :: criteria_eigen              ! criteration for eigen problem
    type(IterationCriterion), target            :: criteria_fsp                ! criteration for fixed source problem
    type(IterationCriterion), target            :: criteria_upscat_eigen       ! criteration for upscatter cycle of eigen problem
    type(IterationCriterion), target            :: criteria_upscat_fsp         ! criteration for upscatter cycle of fsp problem
                                                
    ! adjoint iteration                         
    type(AdjointIteration)                      :: iter_adjoint                 ! parameter for adjoint iteration
                                                
    ! anisotropic scatter flux                  
    type(AnisotropicScatterFlux)                :: flux_scat                    ! anisotropic scatter flux 
        
    ! --------------------------------------------------------------------------
    ! result container
    type(TimeListParameter)                     :: timelist                     ! timelist for state point
                                                
    type(GroupsFlux)                            :: flux_forward                 ! forward scalar and angular flux
    type(GroupsFlux)                            :: flux_adjoint                 ! adjoint scalar and angular flux
    
    type(DistributionParameter)                 :: dist_flux                    ! flux distribution    
    type(DistributionParameter)                 :: dist_power                   ! power DENSITY distribution
    type(DistributionParameter)                 :: dist_fission_rate            ! fission rate distribution
    type(DistributionParameter), allocatable    :: dist_dnps(:)                 ! precursor distribution
    
    ! --------------------------------------------------------------------------
    ! control rod information
    type(ControlRodBank)                        :: cr_bank                      ! control rod bank information

    ! --------------------------------------------------------------------------
    ! perturbation    
    type(SourcePerturbation)                    :: pert_q                       ! external source perturbation
    type(XsecPerturbation)                      :: pert_xsec                    ! cross section perturbation
    type(ControlRodPerturbation)                :: pert_cr                      ! control rod perturbation
    type(MaterialPerturbation)                  :: pert_mat                     ! material perturbation
    type(CasePerturbation)                      :: pert_case                    ! case perturbation
    type(THPerturbation)                        :: pert_th                      ! thermal-hydraulic perturbation
    
    type(CRWorth)                               :: worth_cr 
    
    ! --------------------------------------------------------------------------
    ! transient process
    type(ShapeFunction)                         :: shape_last                   ! last time point shape function
    type(ShapeFunction)                         :: shape_current                ! current time point shape function
    type(ShapeFunction)                         :: shape_predict                ! predict shape function for next time
        
    type(AmplitudeFunction)                     :: amplitude                    ! amplitude function
    type(PointKineticsParameter)                :: pk_parameter                 ! point kinetics parameter
    
    ! --------------------------------------------------------------------------
    ! perturbation related
    type(CrossSection)                          :: xsec_unpert                  ! unperturbated xsec
    type(IterationCounter)                      :: iter_count_unpert            ! unperturbated iteration information
    
    type(GroupsFlux)                            :: flux_forward_unpert          ! unperturbated forward flux
    type(GroupsFlux)                            :: flux_adjoint_pt              ! unperturbated adjoint flux by PT
    type(GroupsFlux)                            :: flux_adjoint_gpt             ! unperturbated adjoint flux by GPT
    
    ! --------------------------------------------------------------------------
    ! feedback parameter
    type(LinkParameter)                         :: self_link
    type(FeedbackParameter)                     :: self_fdbk
    type(LRAmodel)                              :: self_lra
    type(Detector)                              :: det 
    
    ! pk-model 
    type(PrecursorSolver)                       :: dnp_solver 
    type(PKParameter)                           :: pk_param
    type(PKSolver)                              :: pk_solver
    type(PKReactivity)                          :: pk_rho 
    
end module global
