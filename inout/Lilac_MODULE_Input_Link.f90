!--------------------------------------------------------------------------------------------------
!                           Copyright(c)2013 Xi'an Jiaotong University
!--------------------------------------------------------------------------------------------------
!
!  Developed by Nuclear Engineering Computational Physics (NECP) Laboratory.
!
!**************************************************************************************************
!  Version: v1.0                                                                    Date: 4/22/2013
!  Purpose:
!     This is a input reading & file preparation module of Link methods.
!
!  Authors list:
!     Name                       Unit                              Email
!     =========================  ================================  ================================
!     Shengnan Gao               XJTU NECP                         gaoshengnan1989@163.com
!     Yunzhao Li                 XJTU NECP                         Yunzhao@mail.xjtu.edu.cn
!
!  Description:
!     This module reads user input files, lattice code output files included.
!     Only dragon output files can be read for instant.
!     Two main subroutines in this module, InputLink_ReadCR & InputLink_ReadDragon,
!     the first one reads control parameters provided by users and the second one reads abstracted
!     dragon output files. 
!     InputLink_ReadCR also does calculations needed to prepare for fitting/Interpolation,
!     and prepares input files for Link.
!
!  Subroutine list:
!     SUBROUTINE InputLink_ReadCR(User_InputLink,INunit1,InputELALL,INunit2,INunit3,Examine,RCoef)
!     SUBROUTINE InputLink_ReadDragon(User_InputLink,FileNames,LengthFN,NumSInF,NFiles,INunitXSF)
!     SUBROUTINE ParaInfo_GetLength(User_ParaInfo,AuxiDim,Length)
!     SUBROUTINE ParaInfo_SetLength(User_ParaInfo,AuxiDim,Length)
!     SUBROUTINE ParaInfo_GetValues(User_ParaInfo,Paras,ParaLength)
!     SUBROUTINE ParaInfo_SetValues(User_ParaInfo,Paras,ParaLength)
!     SUBROUTINE ParaInfo_Define(User_ParaInfo,AuxiDim,Length)
!     SUBROUTINE ParaInfo_GetStorage(User_ParaInfo,NumINT,NumREAL,NumLOG,NumCHR)
!     SUBROUTINE ParaInfo_Print(User_ParaInfo,OUTunit)
!     SUBROUTINE ParaInfo_Read(User_ParaInfo,INunit,OUTunit)
!     SUBROUTINE ParaInfo_Void(User_ParaInfo)
!     SUBROUTINE TermInfo_SetTerm(User_TermInfo,Operator,Parentheses,TOption, &
!                                 IndexXvars,Option,PolyOrder,NumXvars)
!     SUBROUTINE TermInfo_Define(User_TermInfo,NumXvars)
!     SUBROUTINE TermInfo_GetStorage(User_TermInfo,NumINT,NumREAL,NumLOG,NumCHR)
!     SUBROUTINE TermInfo_Print(User_TermInfo,OUTunit)
!     SUBROUTINE TermInfo_Read(User_TermInfo,INunit,OUTunit)
!     SUBROUTINE TermInfo_Void(User_TermInfo)
!     SUBROUTINE InputLink_DefineCR(User_InputLink,NumTXvars,NumSubDs,NumTterms)
!     SUBROUTINE InputLink_DefineXS(User_InputLink,NumTXvars,NumStates,NumMacs,NumMicNs,LengthMics)
!     SUBROUTINE InputLink_GetStorage(User_TermInfo,NumINT,NumREAL,NumLOG,NumCHR)
!     SUBROUTINE InputLink_Print(User_InputLink,OUTunit)
!     SUBROUTINE InputLink_Read(User_InputLink,INunit,OUTunit)
!     SUBROUTINE InputLink_Void(User_InputLink)
!************************************************************************************************** 

MODULE Input_Link
! USE DFLIB
IMPLICIT NONE

!--------------------------------------------------------------------------------------------------
! To store principal information of each few-group parameter 
!--------------------------------------------------------------------------------------------------
TYPE,PUBLIC :: ParaInfo
   LOGICAL             :: Defined  = .FALSE.
   INTEGER             :: ParaType = 0     ! The type of the parameter considered
   INTEGER             :: AuxiDim  = 0     ! The length of IGstart and IGend
   INTEGER,ALLOCATABLE :: IGstart(:)       ! (AuxiDim) of no-zero begin group numbers
   INTEGER,ALLOCATABLE :: IGend(:)         ! (AuxiDim) of no-zero end group numbers
   INTEGER             :: Length   = 0     ! The length of ParaValues
   REAL(8),ALLOCATABLE :: ParaValues(:)    ! (Length) Values of the few-group parameter    
   CHARACTER(32)       :: ParaName         ! Name of the parameter, not used for calculation
END TYPE ParaInfo

!--------------------------------------------------------------------------------------------------
! Information of how parameters to treat are developped into terms
!--------------------------------------------------------------------------------------------------
TYPE,PUBLIC :: TermInfo
   LOGICAL             :: Defined  = .FALSE.
   INTEGER             :: Operator = 0     ! The calculation operator of the term &
                                           ! 1:+,2:-,3:X,4:/
   INTEGER             :: TOption  = 0     ! The variable to multiply when &
                                           ! the value of the term considered is obtained
   INTEGER             :: Parentheses(2)   ! Number of parentheses before/after the term
   INTEGER             :: NumXvars  = 0    ! The length of IndexXvars/Option/PolyOrder
   INTEGER,ALLOCATABLE :: IndexXvars(:)    ! (NumXvars) Index of X variable positions &
                                           ! in the total list
   INTEGER,ALLOCATABLE :: Option(:)        ! (NumXvars) X variable treatment options
   INTEGER,ALLOCATABLE :: PolyOrder(:)     ! (NumXvars) Fitting orders of X variables
END TYPE TermInfo

!--------------------------------------------------------------------------------------------------
! Three parts: control parameters,term information and few-group parameter information
!--------------------------------------------------------------------------------------------------
TYPE,PUBLIC :: InputLink_Type   
   LOGICAL             :: CRDefined = .FALSE.
   LOGICAL             :: XSDefined = .FALSE.
   INTEGER             :: NumTXvars         ! Number of X variables
   CHARACTER(32),ALLOCATABLE :: XvarNames(:)    ! (NumTXvars) Names of X variables

   INTEGER             :: NumSubDs          ! Number of sub-divisions  
   INTEGER,ALLOCATABLE :: SubDStates(:,:,:) ! (2,NumTXvars,NumSubDs) Start&End state number of &
                                            ! sub divisions
   INTEGER,ALLOCATABLE :: TypStates(:,:)    ! (NumTXvars,NumSubDs) Typical states of each X

   INTEGER,ALLOCATABLE :: NumTerms(:)       ! (NumSubDs) Numbers of terms
   INTEGER             :: NumTterms         ! Equals to the sum of NumTerms(:)
   TYPE(TermInfo),ALLOCATABLE :: Terms(:)       ! (NumTterms) Term information

   INTEGER             :: NumStates
   INTEGER,ALLOCATABLE :: StateSelected(:)  ! (NumStates) Whether or not the state is used 

   TYPE(ParaInfo),ALLOCATABLE :: XvarParas(:,:) ! (NumTXvars,NumStates) X variable values per state

   INTEGER             :: NumEGrps          ! Number of energy groups

   INTEGER             :: NumMacs           ! Number of macroscopic cross sections to be fitted
   TYPE(ParaInfo),ALLOCATABLE :: MacParas(:,:)  ! (NumMacs,NumStates) Few-group parameters

   INTEGER             :: NumMicNs          ! Number of nuclides for micro depletion
   INTEGER,ALLOCATABLE :: TypeNs(:)         ! (NumMicNs) Nuclide types (f/a/general)
   INTEGER,ALLOCATABLE :: NumMicParas(:)    ! (NumMicNs) Number of few-group parameters to be fitted
   INTEGER             :: LengthMics        ! Equals to the sum of NumMicParas(:)
   TYPE(ParaInfo),ALLOCATABLE :: MicParas(:,:)  ! (LengthMics,NumStates) Few-group parameters
END TYPE InputLink_Type


CONTAINS

   !---------------------------------------------------------------------------------------------- 
   ! Read information from user input file CR and write into Link input files
   !----------------------------------------------------------------------------------------------    
   SUBROUTINE InputLink_ReadCR(User_InputLink,INunitCR,InputELALL,INunitXSF,INunitSEG, &
                               Examine,RCoef)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)  :: INunitCR   ! Channel number for the control parameter file "CR"
      INTEGER,INTENT(IN)  :: InputELALL ! Channel number for the "GenaralInfo" , an input file    &
                                        ! for module Link  
      INTEGER,INTENT(IN)  :: INunitXSF  ! Channel number for lattice code output files, also used &
                                        ! for "UALLAXXXXXX", a file with all parameters used for  &
                                        ! fitting/interpolation
      INTEGER,INTENT(IN)  :: INunitSEG  ! Channel number for "INPAXXXXXXXSXXXXXXXX", input files  &
                                        ! for module Link, also used for "LEFTAXXXXXX", a file    &
                                        ! with all parameters not used for fitting/interpolation

      ! Passed out
      TYPE (InputLink_Type),INTENT(OUT) :: User_InputLink ! Storage of information read
      INTEGER,INTENT(OUT) :: Examine                      ! Read from user input files & passed out
      INTEGER,INTENT(OUT) :: RCoef                        ! Read from user input files & passed out

      ! Locals
      INTEGER :: res         ! Returned value of function MAKEDIRQQ in library DFLIB to show if   &
                             ! folders have been made successfully
      INTEGER :: NumFuels    ! Number of fuel (assembly) types
      INTEGER :: NumTXvars   ! Total number of X varibles
      INTEGER :: NumSubDs    ! Number of sub divisions for functionalization
      INTEGER :: NumTterms   ! Total number of terms 
      INTEGER :: NumStatesNU ! Number of states not used for functionalization
      INTEGER :: NumStates   ! Number of states in lattice code output file
      INTEGER :: LengthFN    ! Total number of lattice code output files
      INTEGER :: NumTParts   ! Total number of "parts" in which used state numbers are laid out
      INTEGER :: LTypiStates ! Length of the array containing all typical state numbers
      INTEGER :: TermInSub   ! In which sub division the term is
      INTEGER :: LengthIR    ! Number of X that don't exist in the term considered
      INTEGER :: YLength     ! Total number of Y varibales
      INTEGER :: SameY = 1   ! Whether all macro-scopic XS have the same AuxiDim, 1:Yes
      INTEGER :: SelectAll=1 ! Whether all states in lattice code input files are selected, 1:Yes
      INTEGER :: SameTerms=1 ! Whether in all divisions the terms are the same, 1:Yes

      ! The following variables are local ones used when reading information of terms.
      ! Their significations are exactly the same with those in type "TermInfo"
      INTEGER :: Operator,NumXvars,TOption
      INTEGER :: Parentheses(2)

      ! The following variables are local ones used when reading information of User_InputLink.
      ! Their significations are exactly the same with those in type "InputLink_Type"
      INTEGER :: NumEGrps,NumMacs,NumMicNs,LengthMics
      INTEGER :: ParaType,AuxiDim,Length

      ! The following variables are all used for "DO" blocks
      INTEGER :: i,j,k,m,n,jj,kk,mm,nn
      INTEGER :: IDFileName,IDStartF,IDBG,IDNUsed,IDMIC,IDSubDs,IDMicPara
      INTEGER :: IDIR,IDTerm,IDPart,IDXbar,ID,IDStateInv,IDPartInv,IDTS

      ! Whether XXX micro-scopic XS is found
      ! XXX stands for: A-absorption,F-fission,C-capture,NuF-nusigma fission,KF-kappasigma fission
      LOGICAL :: AFound,FFound,CFound
      LOGICAL :: NuFFound,KFFound

      ! The following variables are used to find out if a term has common X variables with the
      ! basic term. 
      ! If yes, typical states are found by searching for the ones which have equal X values.
      LOGICAL :: XEXIST,CommonX,INCLUDED,XTYPICAL,EquValue

      REAL(8) :: min,Temp
      REAL(8) :: FormatPrint = 0.0D0

      REAL(8),PARAMETER          :: NECP_RealZero  = 1.0E-8
      REAL(8),PARAMETER          :: NECP_RealEqual = 1.0E-8

      ! Names of folders and files
      CHARACTER(32)              :: ParaName
      CHARACTER(12)              :: Folder_Name1 = "ExamineFiles"    ! Folder of "UALLAXXXXXXXX"  &
                                                                     ! "LEFTAXXXXXXXX"
      CHARACTER(14)              :: Folder_Name2 = "LinkInputFiles"  ! Folder of "INPA...S..."
      CHARACTER(15)              :: Folder_Name3 = "LinkOutputFiles" ! Folder of output files for &
                                                                     ! module Link
      CHARACTER(11)              :: Folder_Name4 = "CompareRSLT"     ! Folder of files containing &
                                                                     ! comparaison  of            &
                                                                     ! fitting/interpolation and  &
                                                                     ! lattice code results
      ! Input files for module Link
      CHARACTER(36)              :: LINFileName  = "LinkInputFiles/INPA00000000S00000000"
      ! File with all states used for functionalization, for examining
      CHARACTER(26)              :: EXUFileName  = "ExamineFiles/UALLA00000000"
      ! File with all states not used for functionalization, for examining
      CHARACTER(26)              :: EXLFileName  = "ExamineFiles/LEFTA00000000"

      ! The following arrays are local ones used when reading information of User_InputLink.
      ! Their significations are exactly the same with those in type "InputLink_Type"/"TermInfo"
      INTEGER,ALLOCATABLE        :: IndexXvars(:),Option(:),PolyOrder(:),StateSelected(:)
      INTEGER,ALLOCATABLE        :: IGstart(:),IGend(:),TypeNs(:),NumParts(:),NumMicParas(:)
      INTEGER,ALLOCATABLE        :: SubDStates(:,:,:),TypStates(:,:),NumTerms(:),MicTags(:)

      INTEGER,ALLOCATABLE        :: NumFiles(:)      ! Number of lattice code files for each fuel
      INTEGER,ALLOCATABLE        :: NumSInF(:)       ! Number of states in each file
      INTEGER,ALLOCATABLE        :: TermStates(:,:)  ! (2,NumTParts), start & end state number  &
                                                     ! for the terms
      INTEGER,ALLOCATABLE        :: IDXDIF(:)        ! Number of different values of each X
      INTEGER,ALLOCATABLE        :: SelfPosition(:)  ! Specially for scatter vectors, positions &
                                                     ! of self-scatter parameter
      INTEGER,ALLOCATABLE        :: BaseTerm(:)      ! Basic term in each subdivision
      INTEGER,ALLOCATABLE        :: IndexRest(:)     ! Index of X not existing in the term 
      INTEGER,ALLOCATABLE        :: NStatesInTerm(:) ! Number of states in terms
      INTEGER,ALLOCATABLE        :: TypicalStates(:)
      REAL(8),ALLOCATABLE        :: Xav(:),STD(:)    ! Average value & 方差
      REAL(8),ALLOCATABLE        :: XDifValues(:,:)  ! Different values of X
      REAL(8),ALLOCATABLE        :: Xbar(:)          ! 变量代换后的X
      LOGICAL,ALLOCATABLE        :: IndexCommon(:)   ! Index of common X in typical states
      CHARACTER(32),ALLOCATABLE  :: FuelName(:),XvarNames(:)  
      CHARACTER(8),ALLOCATABLE   :: FileNames(:)
      TYPE(TermInfo),ALLOCATABLE :: Terms(:)
      TYPE(ParaInfo),ALLOCATABLE :: XvarParas(:,:),MacParas(:,:),MicParas(:,:)      
      TYPE(ParaInfo),ALLOCATABLE :: CopyMacParas(:,:),CopyMicParas(:,:)      

      ! Creat folders by calling functions of library DFLIB
!      res = MAKEDIRQQ(Folder_Name1)
!      res = MAKEDIRQQ(Folder_Name2)
!      res = MAKEDIRQQ(Folder_Name3)
!      res = MAKEDIRQQ(Folder_Name4)

      ! Open user files and read information
      open (99,file="testfile")
      OPEN (INunitCR,file="CR")
      OPEN (InputELALL,file="LinkInputFiles/GenaralInfo") ! A temporary file needed by Link module
      ! Read CR file in which control parameters are available & write into file "GenaralInfo"
      READ (INunitCR,*) NumFuels
      WRITE(InputELALL,*) NumFuels
      WRITE(InputELALL,*) "**************************"
      ALLOCATE (FuelName(NumFuels),NumFiles(NumFuels))
      READ (INunitCR,*) (FuelName(i),i=1,NumFuels)
      READ (INunitCR,*) (NumFiles(i),i=1,NumFuels)
      LengthFN = 0
      DO i=1,NumFuels
         LengthFN = LengthFN + NumFiles(i)
      END DO
      ALLOCATE (NumSInF(LengthFN),FileNames(LengthFN))
      READ (INunitCR,*) (FileNames(i),i=1,LengthFN)
      READ (INunitCR,*) (NumSInf(i),i=1,LengthFN)
      IDFileName = 0
      IDBG = 1
      IDStartF = 0
      READ (INunitCR,*)
      ! Block for reading information of each fuel provided in file "CR"
      DO i=1,NumFuels
         ! Read the first part of file "CR" : general information
         READ (INunitCR,*) NumTXvars 
         !write (*,*) NumTXvars 
         ALLOCATE (XvarNames(NumTXvars))
         READ (INunitCR,*) (XvarNames(j),j=1,NumTXvars)
         !write (*,*) (XvarNames(j),j=1,NumTXvars)
         READ (INunitCR,*) NumSubDs
         !write (*,*) NumSubDs
         WRITE(InputELALL,*) NumSubDs
         ALLOCATE (SubDStates(2,NumTXvars,NumSubDs))
         READ (INunitCR,*) (((SubDStates(j,m,n),j=1,2),m=1,NumTXvars),n=1,NumSubDs)
         !write (*,*) (((SubDStates(j,m,n),j=1,2),m=1,NumTXvars),n=1,NumSubDs)
         ALLOCATE (TypStates(NumTXvars,NumSubDs))
         READ (INunitCR,*) ((TypStates(m,n),m=1,NumTXvars),n=1,NumSubDs)
         !write (*,*) ((TypStates(m,n),m=1,NumTXvars),n=1,NumSubDs)
         READ (INunitCR,*)
         ALLOCATE (NumTerms(NumSubDs))
         READ (INunitCR,*) SameTerms
         IF(SameTerms.EQ.1) THEN
            READ (INunitCR,*) NumTerms(1)
            DO j=1,NumSubDs
               NumTerms(j) = NumTerms(1)
            END DO
         ELSE
            READ (INunitCR,*) (NumTerms(j),j=1,NumSubDs)
         END IF
         !write (*,*) (NumTerms(j),j=1,NumSubDs)
         NumTterms = 0
         DO j=1,NumSubDs
            NumTterms = NumTterms + NumTerms(j)
         END DO
         !write (*,*) "NumTterms",NumTterms
         ! Set control parameter values in User_InputLink
         CALL InputLink_DefineCR(User_InputLink,NumTXvars,NumSubDs,NumTterms)
         DO j=1,NumTXvars
            User_InputLink%XvarNames(j) = XvarNames(j)
         END DO
         DO m=1,NumSubDs
            DO n=1,NumTXvars
               DO j=1,2
                  User_InputLink%SubDStates(j,n,m) = SubDStates(j,n,m)
               END DO
               User_InputLink%TypStates(n,m) = TypStates(n,m)
            END DO
            User_InputLink%NumTerms(m) = NumTerms(m)
         END DO
         DEALLOCATE (XvarNames,SubDStates,TypStates,NumTerms)
         ! Read the second part of "CR" : term information
         IF(SameTerms.NE.1) THEN
            ! Read information of each term
            DO j=1,NumTterms
               READ (INunitCR,*) Operator,(Parentheses(m),m=1,2),TOption
               READ (INunitCR,*) NumXvars
               ALLOCATE (IndexXvars(NumXvars),Option(NumXvars),PolyOrder(NumXvars))
               READ (INunitCR,*) (IndexXvars(m),m=1,NumXvars), &
                                 (Option(m),m=1,NumXvars),     &
                                 (PolyOrder(m),m=1,NumXvars)
               CALL TermInfo_Define(User_InputLink%Terms(j),NumXvars)
               CALL TermInfo_SetTerm(User_InputLink%Terms(j),Operator,Parentheses,TOption, &
                                     IndexXvars,Option,PolyOrder,NumXvars)
               DEALLOCATE (IndexXvars,Option,PolyOrder)
            END DO
            READ (INunitCR,*)
         ELSE
            ! Since all subdivisions have same terms, information of only terms &
            ! of the first subdivision needs to be read
            DO j=1,User_InputLink%NumTerms(1)
               READ (INunitCR,*) Operator,(Parentheses(m),m=1,2),TOption
               READ (INunitCR,*) NumXvars
               ALLOCATE (IndexXvars(NumXvars),Option(NumXvars),PolyOrder(NumXvars))
               READ (INunitCR,*) (IndexXvars(m),m=1,NumXvars), &
                                 (Option(m),m=1,NumXvars),     &
                                 (PolyOrder(m),m=1,NumXvars)
               CALL TermInfo_Define(User_InputLink%Terms(j),NumXvars)
               CALL TermInfo_SetTerm(User_InputLink%Terms(j),Operator,Parentheses,TOption, &
                                     IndexXvars,Option,PolyOrder,NumXvars)
               DO m=2,NumSubDs
                  CALL TermInfo_Define(User_InputLink%Terms((m-1)*User_InputLink%NumTerms(1)+j),  &
                                       NumXvars)
                  CALL TermInfo_SetTerm(User_InputLink%Terms((m-1)*User_InputLink%NumTerms(1)+j), &
                                        Operator,Parentheses,TOption, &
                                        IndexXvars,Option,PolyOrder,NumXvars)
               END DO
               DEALLOCATE (IndexXvars,Option,PolyOrder)
            END DO
            READ (INunitCR,*) 
         END IF
         ! Read information about how few-group parameters are provided in input files
         ALLOCATE (NumParts(NumTterms))
         READ (INunitCR,*) (NumParts(j),j=1,NumTterms)
         NumTParts = 0
         DO j=1,NumTterms          
            NumTParts = NumTParts + NumParts(j)
         END DO
         ALLOCATE (TermStates(2,NumTParts))
         READ (INunitCR,*) ((TermStates(m,n),m=1,2),n=1,NumTParts)
         READ (INunitCR,*) 
         READ (INunitCR,*) NumStates,SelectAll
         !write (*,*) NumStates,SelectAll
         READ (INunitCR,*)
         ALLOCATE (StateSelected(NumStates))
         DO j=1,NumStates
            StateSelected(j) = -1
         END DO
         ! Get the length of typical case array for not typical cases
         ! Define for which term each state is used 
         IDPart = 1
         LTypiStates = 0
         DO j=1,NumTterms 
            DO m=IDPart,IDPart+NumParts(j)-1
               DO n=TermStates(1,m),TermStates(2,m)
                  StateSelected(n) = j
               END DO
               IF(User_InputLink%Terms(j)%Operator.NE.0) THEN
                  LTypiStates = LTypiStates + TermStates(2,m) - TermStates(1,m) + 1
               END IF
            END DO
            IDPart = IDPart + NumParts(j)
         END DO
         ! Read the last part of CR : information about the number of FG parameters & their types
         ! Necessary information for difining the length of arrays comes first
         ALLOCATE (TypicalStates(LTypiStates))
         write (99,*) "test",LTypiStates
         READ (INunitCR,*) NumEGrps,NumMacs,NumMicNs
         WRITE(InputELALL,"(3(2X,I4))") NumEGrps,NumMacs,NumMicNs
         ALLOCATE (TypeNs(NumMicNs),NumMicParas(NumMicNs),MicTags(NumMicNs))
         READ (INunitCR,*) (MicTags(j),j=1,NumMicNs)
         WRITE(InputELALL,"(5(2X,I6))") (MicTags(j),j=1,NumMicNs)
         DEALLOCATE (MicTags)
         READ (INunitCR,*) (TypeNs(j),j=1,NumMicNs)
         WRITE(InputELALL,"(5(2X,I4))") (TypeNs(j),j=1,NumMicNs)
         READ (INunitCR,*) (NumMicParas(j),j=1,NumMicNs)
         WRITE(InputELALL,"(5(2X,I4))") (NumMicParas(j),j=1,NumMicNs)
         LengthMics = 0
         DO j=1,NumMicNs
            LengthMics = LengthMics + NumMicParas(j)
         END DO
         ! X variable information
         ALLOCATE (XvarParas(NumTXvars,NumStates))
         DO j=1,NumStates
            DO m=1,NumTXvars
               AuxiDim = 1
               Length  = 1
               CALL ParaInfo_Define(XvarParas(m,j),AuxiDim,Length)
               XvarParas(m,j)%IGstart(1) = 1
               XvarParas(m,j)%IGend(1) = 1
               XvarParas(m,j)%ParaName = User_InputLink%XvarNames(m)
               XvarParas(m,j)%ParaType = 0
            END DO
         END DO
         ! Y variable information
         READ (INunitCR,*) SameY
         WRITE(InputELALL,*) SameY
         ! Two types of FG parameters mainly : macro XS related & micro XS related
         ALLOCATE (MacParas(NumMacs,NumStates),MicParas(LengthMics,NumStates))
         ALLOCATE (CopyMacParas(NumMacs,NumStates),CopyMicParas(LengthMics,NumStates))
         IF(SameY.EQ.1) THEN
            ! Macro XS related information
            DO j=1,NumMacs
               READ (INunitCR,*) AuxiDim
               ALLOCATE (IGstart(AuxiDim),IGend(AuxiDim))
               DO m=1,AuxiDim
                  READ (INunitCR,*) IGstart(m),IGend(m)
               END DO
               Length = 0
               DO m=1,AuxiDim
                  Length = Length + IGend(m) - IGstart(m) + 1
               END DO
               READ (INunitCR,*) ParaType,ParaName
               WRITE (InputELALL,"(3(2X,I4))") j,ParaType,Length
               ALLOCATE (SelfPosition(AuxiDim))
               IF(ParaType.EQ.108) THEN
                  READ (INunitCR,*) (SelfPosition(m),m=1,AuxiDim)
                  WRITE (InputELALL,"(2X,I4)") AuxiDim
                  DO m=1,AuxiDim
                     WRITE (InputELALL,"(3(2X,I4))") IGstart(m),IGend(m),SelfPosition(m)
                  END DO
               END IF
               DEALLOCATE (SelfPosition)
               DO m=1,NumStates
                  CALL ParaInfo_Define(MacParas(j,m),AuxiDim,Length)
                  CALL ParaInfo_Define(CopyMacParas(j,m),AuxiDim,Length)
                  DO n=1,AuxiDim
                     MacParas(j,m)%IGstart(n) = IGstart(n)
                     MacParas(j,m)%IGend(n) = IGend(n)
                     CopyMacParas(j,m)%IGstart(n) = IGstart(n)
                     CopyMacParas(j,m)%IGend(n) = IGend(n)
                  END DO
                  MacParas(j,m)%ParaName = ParaName
                  MacParas(j,m)%ParaType = ParaType
                  CopyMacParas(j,m)%ParaName = ParaName
                  CopyMacParas(j,m)%ParaType = ParaType
               END DO
               DEALLOCATE (IGstart,IGend)
            END DO
            WRITE (InputELALL,*) "******************************"
            ! Micro XS related information
            DO j=1,LengthMics
               READ (INunitCR,*) AuxiDim
               ALLOCATE (IGstart(AuxiDim),IGend(AuxiDim))
               DO m=1,AuxiDim
                  READ (INunitCR,*) IGstart(m),IGend(m)
               END DO
               Length = 0
               DO m=1,AuxiDim
                  Length = Length + IGend(m) - IGstart(m) + 1
               END DO
               READ (INunitCR,*) ParaType,ParaName
               WRITE (InputELALL,"(3(2X,I4))") j,ParaType,Length
               DO m=1,NumStates
                  CALL ParaInfo_Define(MicParas(j,m),AuxiDim,Length)
                  CALL ParaInfo_Define(CopyMicParas(j,m),AuxiDim,Length)
                  DO n=1,AuxiDim
                     MicParas(j,m)%IGstart(n) = IGstart(n)
                     MicParas(j,m)%IGend(n) = IGend(n)
                     CopyMicParas(j,m)%IGstart(n) = IGstart(n)
                     CopyMicParas(j,m)%IGend(n) = IGend(n)
                  END DO
                  MicParas(j,m)%ParaName = ParaName
                  MicParas(j,m)%ParaType = ParaType
                  CopyMicParas(j,m)%ParaName = ParaName
                  CopyMicParas(j,m)%ParaType = ParaType
               END DO
               DEALLOCATE (IGstart,IGend)
            END DO
         ELSE
            ! Macro XS related information
            DO m=1,NumStates
               DO n=1,NumMacs
                  READ (INunitCR,*) AuxiDim
                  ALLOCATE (IGstart(AuxiDim),IGend(AuxiDim))
                  DO k=1,AuxiDim
                     READ (INunitCR,*) IGstart(k),IGend(k)
                  END DO
                  Length = 0
                  DO k=1,AuxiDim
                     Length = Length + IGend(m) - IGstart(m) + 1
                  END DO
                  READ (INunitCR,*) ParaType,ParaName
                  WRITE (InputELALL,"(3(2X,I4))") n,ParaType,Length
                  CALL ParaInfo_Define(MacParas(m,n),AuxiDim,Length)
                  CALL ParaInfo_Define(CopyMacParas(m,n),AuxiDim,Length)
                  DO k=1,AuxiDim
                     MacParas(m,j)%IGstart(k) = IGstart(k)
                     MacParas(m,j)%IGend(k) = IGend(k)
                     CopyMacParas(m,j)%IGstart(k) = IGstart(k)
                     CopyMacParas(m,j)%IGend(k) = IGend(k)
                  END DO 
                  MacParas(m,j)%ParaName = ParaName
                  MacParas(m,j)%ParaType = ParaType
                  CopyMacParas(m,j)%ParaName = ParaName
                  CopyMacParas(m,j)%ParaType = ParaType
                  DEALLOCATE (IGstart,IGend)
               END DO
            END DO
            ! Micro XS related information
            DO m=1,NumStates
               DO n=1,LengthMics
                  READ (INunitCR,*) AuxiDim
                  ALLOCATE (IGstart(AuxiDim),IGend(AuxiDim))
                  DO k=1,AuxiDim
                     READ (INunitCR,*) IGstart(k),IGend(k)
                  END DO
                  Length = 0
                  DO k=1,AuxiDim
                     Length = Length + IGend(m) - IGstart(m) + 1
                  END DO
                  READ (INunitCR,*) ParaType,ParaName
                  WRITE (InputELALL,"(3(2X,I4))") n,ParaType,Length
                  CALL ParaInfo_Define(MacParas(m,n),AuxiDim,Length)
                  CALL ParaInfo_Define(CopyMacParas(m,n),AuxiDim,Length)
                  DO k=1,AuxiDim
                     MicParas(j,m)%IGstart(k) = IGstart(k)
                     MicParas(j,m)%IGend(k) = IGend(k)
                     CopyMicParas(j,m)%IGstart(k) = IGstart(k)
                     CopyMicParas(j,m)%IGend(k) = IGend(k)
                  END DO
                  MicParas(m,j)%ParaName = ParaName
                  MicParas(m,j)%ParaType = ParaType
                  CopyMicParas(m,j)%ParaName = ParaName
                  CopyMicParas(m,j)%ParaType = ParaType
                  DEALLOCATE (IGstart,IGend)
               END DO
            END DO
         END IF
         READ (INunitCR,*) 

         ! Read XS file
         ! Define the XS part of User_InputLink
         write(*,*) "IDStartF",IDStartF
         CALL InputLink_DefineXS(User_InputLink,NumTXvars,NumStates,NumMacs,NumMicNs,LengthMics)
         IDBG = 1
         ! Define ParaInfo type objects in User_InputLink
         DO j=1,NumFiles(i)
            IDFileName = IDStartF + 1
            DO jj=IDBG,IDBG+NumSInF(IDFileName)-1
               DO m=1,NumTXvars
                  CALL ParaInfo_Define(User_InputLink%XvarParas(m,jj),XvarParas(m,jj)%AuxiDim, &
                                       XvarParas(m,jj)%Length)
               END DO
               DO m=1,NumMacs
                  CALL ParaInfo_Define(User_InputLink%MacParas(m,jj),MacParas(m,jj)%AuxiDim,   &
                                       MacParas(m,jj)%Length)
               END DO
               DO m=1,LengthMics
                  CALL ParaInfo_Define(User_InputLink%MicParas(m,jj),MicParas(m,jj)%AuxiDim,   &
                                       MicParas(m,jj)%Length)
               END DO
            END DO
            IDBG = IDBG + NumSInF(IDFileName)
         END DO
         ! Read XS files & store information in User_InputLink
         CALL InputLink_ReadDragon(User_InputLink,FileNames,LengthFN,IDStartF,NumSInF, &
                                   NumFiles(i),INunitXSF)
         IDBG = 1
         ! Re-put X & Y values in local variables for later use
         DO j=1,NumFiles(i)
            IDFileName = IDStartF + 1
            !write (*,*) INunitXSF,IDFileName,FileNames(IDFileName)
            !WRITE (*,*) IDBG,IDBG+NumSInF(IDFileName)-1
            DO jj=IDBG,IDBG+NumSInF(IDFileName)-1
               DO m=1,NumTXvars
                  CALL ParaInfo_GetValues(User_InputLink%XvarParas(m,jj),XvarParas(m,jj)%ParaValues,&
                                          XvarParas(m,jj)%Length)
               END DO
               DO m=1,NumMacs
                  CALL ParaInfo_GetValues(User_InputLink%MacParas(m,jj),MacParas(m,jj)%ParaValues,  &
                                          MacParas(m,jj)%Length)
               END DO
               DO m=1,LengthMics
                  CALL ParaInfo_GetValues(User_InputLink%MicParas(m,jj),MicParas(m,jj)%ParaValues,  &
                                          MicParas(m,jj)%Length)
               END DO
            END DO
            IDBG = IDBG + NumSInF(IDFileName)
         END DO

         ! Micro-depletion preparation
         ! The idea is to see for each macro-scopic XS if there are N*(micro-scopic XS) to minus
         DO j=1,NumStates
            DO m=1,NumMacs
               IF(MacParas(m,j)%ParaType.EQ.101) THEN ! Type 101 : Total XS
                  IDMIC = 1 
                  ! For all nuclides, do:
                  DO n=1,NumMicNs
                     AFound = .FALSE. ! Whether absorption micro-scopic XS exists
                     FFound = .FALSE. ! Whether fission micro-scopic XS exists
                     CFound = .FALSE. ! Whether capture micro-scopic XS exists
                     ! If absorption micro-scopic XS exists, then only absorption micro-scopic XS &
                     ! is taken out, else capture & fission micro-scopic XS will be searched for  &
                     ! and taken out
                     DO k=IDMIC,IDMIC+NumMicParas(n)-2
                        IF(MicParas(k,j)%ParaType.EQ.207) THEN ! absorption micro-scopic XS
                           DO jj=1,MacParas(m,j)%Length
                              MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                              MicParas(k,j)%ParaValues(jj) * &
                              MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                           END DO
                           AFound = .TRUE.
                        END IF
                     END DO
                     IF(.NOT.AFound) THEN
                        DO k=IDMIC,IDMIC+NumMicParas(n)-2
                           IF(MicParas(k,j)%ParaType.EQ.203) THEN ! Fission micro-scopic XS
                              DO jj=1,MacParas(m,j)%Length
                                 MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                                 MicParas(k,j)%ParaValues(jj) * &
                                 MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                              END DO
                              FFound = .TRUE.
                           END IF
                           IF(MicParas(k,j)%ParaType.EQ.206) THEN ! Capture micro-scopic XS
                              DO jj=1,MacParas(m,j)%Length
                                 MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                                 MicParas(k,j)%ParaValues(jj) * &
                                 MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                              END DO
                              CFound = .TRUE.
                           END IF
                        END DO
                        IF((.NOT.FFound).AND.(TypeNs(n).EQ.2)) THEN ! TypeNs=2 : Fission nuclide
                           ! WRITE(*,*) "[INPUT_LINK] WARNING: No fission XS for nuclide",n 
                        END IF
                        IF(.NOT.CFound) THEN
                           ! WRITE(*,*) "[INPUT_LINK] WARNING: No capture XS for nuclide",n 
                        END IF
                     END IF
                     IDMIC = IDMIC + NumMicParas(n)
                  END DO                   
               ELSE IF(MacParas(m,j)%ParaType.EQ.102) THEN ! Type 102 : Transport XS
                  IDMIC = 1 
                  DO n=1,NumMicNs
                     AFound = .FALSE.
                     FFound = .FALSE.
                     CFound = .FALSE.
                     DO k=IDMIC,IDMIC+NumMicParas(n)-2
                        IF(MicParas(k,j)%ParaType.EQ.207) THEN
                           DO jj=1,MacParas(m,j)%Length
                              MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                              MicParas(k,j)%ParaValues(jj) * &
                              MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                           END DO
                           AFound = .TRUE.
                        END IF
                     END DO
                     IF(.NOT.AFound) THEN
                        DO k=IDMIC,IDMIC+NumMicParas(n)-2
                           IF(MicParas(k,j)%ParaType.EQ.203) THEN
                              DO jj=1,MacParas(m,j)%Length
                                 MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                                 MicParas(k,j)%ParaValues(jj) * &
                                 MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                              END DO
                              FFound = .TRUE.
                           END IF
                           IF(MicParas(k,j)%ParaType.EQ.206) THEN
                              DO jj=1,MacParas(m,j)%Length
                                 MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                                 MicParas(k,j)%ParaValues(jj) * &
                                 MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                              END DO
                              CFound = .TRUE.
                           END IF
                        END DO
                        IF((.NOT.FFound).AND.(TypeNs(n).EQ.2)) THEN
                           ! WRITE(*,*) "[INPUT_LINK] WARNING: No fission XS for nuclide",n 
                        END IF
                        IF(.NOT.CFound) THEN
                           ! WRITE(*,*) "[INPUT_LINK] WARNING: No capture XS for nuclide",n 
                        END IF
                     END IF
                     IDMIC = IDMIC + NumMicParas(n)
                  END DO
               ELSE IF(MacParas(m,j)%ParaType.EQ.107) THEN ! Type 107 : Absorption XS
                  IDMIC = 1 
                  DO n=1,NumMicNs
                     ! Fission nuclide
                     AFound = .FALSE.
                     FFound = .FALSE.
                     CFound = .FALSE.
                     DO k=IDMIC,IDMIC+NumMicParas(n)-2
                        IF(MicParas(k,j)%ParaType.EQ.207) THEN
                           DO jj=1,MacParas(m,j)%Length
                              MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                              MicParas(k,j)%ParaValues(jj) * &
                              MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                           END DO
                           AFound = .TRUE.
                        END IF
                     END DO
                     IF(.NOT.AFound) THEN
                        DO k=IDMIC,IDMIC+NumMicParas(n)-2
                           IF(MicParas(k,j)%ParaType.EQ.203) THEN
                              DO jj=1,MacParas(m,j)%Length
                                 MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                                 MicParas(k,j)%ParaValues(jj) * &
                                 MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                              END DO
                              FFound = .TRUE.
                           END IF
                           IF(MicParas(k,j)%ParaType.EQ.206) THEN
                              DO jj=1,MacParas(m,j)%Length
                                 MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                                 MicParas(k,j)%ParaValues(jj) * &
                                 MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                                 ! write (*,*) 'mic',k,jj,IDMIC+NumMicParas(n)-1
                              END DO
                              CFound = .TRUE.
                           END IF
                        END DO
                        IF((.NOT.FFound).AND.(TypeNs(n).EQ.2)) THEN
                           ! WRITE(*,*) "[INPUT_LINK] WARNING: No fission XS for nuclide",n 
                        END IF
                        IF(.NOT.CFound) THEN
                           ! WRITE(*,*) "[INPUT_LINK] WARNING: No capture XS for nuclide",n 
                        END IF
                     END IF
                     IDMIC = IDMIC + NumMicParas(n)
                  END DO                  
               ELSE IF(MacParas(m,j)%ParaType.EQ.103) THEN ! Type 103 : Fission XS
                  IDMIC = 1 
                  FFound = .FALSE.
                  DO n=1,NumMicNs
                     ! Fission nuclide
                     DO k=IDMIC,IDMIC+NumMicParas(n)-2
                        IF(MicParas(k,j)%ParaType.EQ.203) THEN
                           DO jj=1,MacParas(m,j)%Length
                              MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                              MicParas(k,j)%ParaValues(jj) * &
                              MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                           END DO
                           FFound = .TRUE.
                        END IF
                     END DO
                     IF((.NOT.FFound).AND.(TypeNs(n).EQ.2)) THEN
                        ! WRITE(*,*) "[INPUT_LINK] WARNING: No fission XS for nuclide",n 
                     END IF
                     IDMIC = IDMIC + NumMicParas(n)
                  END DO                  
               ELSE IF(MacParas(m,j)%ParaType.EQ.104) THEN ! Type 104 : NuSigF 
                  IDMIC = 1 
                  DO n=1,NumMicNs
                     NuFFound = .FALSE.
                     ! Fission nuclide
                     DO k=IDMIC,IDMIC+NumMicParas(n)-2
                        IF(MicParas(k,j)%ParaType.EQ.204) THEN
                           DO jj=1,MacParas(m,j)%Length
                              MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                              MicParas(k,j)%ParaValues(jj) * &
                              MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                           END DO
                           NuFFound = .TRUE.
                        END IF
                     END DO
                     IF((.NOT.NuFFound).AND.(TypeNs(n).EQ.2)) THEN
                        ! WRITE(*,*) "[INPUT_LINK] WARNING: No nu-fission XS for nuclide",n 
                     END IF
                     IDMIC = IDMIC + NumMicParas(n)
                  END DO                
               ELSE IF(MacParas(m,j)%ParaType.EQ.105) THEN ! Type 105 : KappaSigF 
                  IDMIC = 1 
                  DO n=1,NumMicNs
                     KFFound = .FALSE.
                     ! Fission nuclide
                     DO k=IDMIC,IDMIC+NumMicParas(n)-2
                        IF(MicParas(k,j)%ParaType.EQ.205) THEN
                           DO jj=1,MacParas(m,j)%Length
                              MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                              MicParas(k,j)%ParaValues(jj) * &
                              MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                           END DO
                           KFFound = .TRUE.
                        END IF
                     END DO
                     IF((.NOT.KFFound).AND.(TypeNs(n).EQ.2)) THEN
                        ! WRITE(*,*) "[INPUT_LINK] WARNING: No kappa-fission XS for nuclide",n 
                     END IF
                     IDMIC = IDMIC + NumMicParas(n)
                  END DO                               
               ELSE IF(MacParas(m,j)%ParaType.EQ.106) THEN ! Type 106 : Capture XS 
                  IDMIC = 1 
                  DO n=1,NumMicNs
                     ! Fission nuclide
                     CFound = .FALSE.
                     DO k=IDMIC,IDMIC+NumMicParas(n)-2
                        IF(MicParas(k,j)%ParaType.EQ.206) THEN
                           DO jj=1,MacParas(m,j)%Length
                              MacParas(m,j)%ParaValues(jj) = MacParas(m,j)%ParaValues(jj) - &
                              MicParas(k,j)%ParaValues(jj) * &
                              MicParas(IDMIC+NumMicParas(n)-1,j)%ParaValues(1)
                           END DO
                           CFound = .TRUE.
                        END IF
                     END DO
                     IF(.NOT.CFound) THEN
                        ! WRITE(*,*) "[INPUT_LINK] WARNING: No capture XS for nuclide",n 
                     END IF
                     IDMIC = IDMIC + NumMicParas(n)
                  END DO                               
               END IF
            END DO
         END DO

         ! Write into examine files
         WRITE (EXUFileName(19:26),"(I8.8)") i
         WRITE (EXLFileName(19:26),"(I8.8)") i
         OPEN (INunitXSF,file=EXUFileName)
         OPEN (INunitSEG,file=EXLFileName)
         NumStatesNU = 0
         DO j=1,NumStates
            IF(StateSelected(j).LT.0) THEN
               NumStatesNU = NumStatesNU + 1
            END IF
         END DO
         WRITE (INunitXSF,'(2(I8))') NumTXvars,NumStates - NumStatesNU
         WRITE (INunitXSF,*) "**************************************"
         WRITE (INunitSEG,'(2(I8))') NumTXvars,NumStatesNU
         WRITE (INunitSEG,*) "**************************************"
         DO j=1,NumStates
            ! States not used are written into LeftAXXXXXXXX file through INunitSEG
            IF(StateSelected(j).LT.0) THEN
               WRITE (INunitSEG,'(I8)') j
               WRITE (INunitSEG,100) (XvarParas(m,j)%ParaValues(1),m=1,NumTXvars)
               DO m=1,NumMacs
                  WRITE (INunitSEG,100) (MacParas(m,j)%ParaValues(n),n=1,MacParas(m,j)%Length)
               END DO
               DO m=1,LengthMics
                  WRITE (INunitSEG,100) (MicParas(m,j)%ParaValues(n),n=1,MicParas(m,j)%Length)
               END DO
            ! States used are written into UsedAXXXXXXXX file through INunitXSF
            ELSE
               WRITE (INunitXSF,'(I8)') j
               WRITE (INunitXSF,100) (XvarParas(m,j)%ParaValues(1),m=1,NumTXvars)
               DO m=1,NumMacs
                  WRITE (INunitXSF,100) (MacParas(m,j)%ParaValues(n),n=1,MacParas(m,j)%Length)
               END DO
               DO m=1,LengthMics
                  WRITE (INunitXSF,100) (MicParas(m,j)%ParaValues(n),n=1,MicParas(m,j)%Length)
               END DO
            END IF
         END DO
         CLOSE (INunitXSF)
         CLOSE (INunitSEG)
         ! 考虑paratype的位置换一下喔
         ! Segments         
         ALLOCATE (XDifValues(NumTXvars,NumStates),IDXDIF(NumTXvars))
         ! 注意，这里XDifValues(NumStates,NumTXvars)初值给负值
         DO m=1,NumTXvars
            DO n=1,NumStates
               XDifValues(m,n) = -1.0D0
            END DO
         END DO
         ! Stock all different X values in XDifValues
         DO m=1,NumTXvars
            IDXDIF(m) = 1
            DO n=1,NumStates
               XEXIST = .FALSE.
               DO k=1,IDXDIF(m)
                  !write (*,*) 'x',m,'value',XvarParas(m,n)%ParaValues(1),XDifValues(m,k),XEXIST
                  IF(ABS(XDifValues(m,k)-0.0).LT.NECP_RealZero) THEN
                     IF(ABS(XDifValues(m,k)-XvarParas(m,n)%ParaValues(1)).LT.1.0E-18) THEN
                        XEXIST = .TRUE.
                        EXIT
                     END IF
                  ELSE
                     IF(ABS((XDifValues(m,k)-XvarParas(m,n)%ParaValues(1))/XDifValues(m,k)).LT. &
                        1.0E-8) THEN
                        XEXIST = .TRUE.
                        EXIT
                     END IF
                  END IF
               END DO
               IF(.NOT.XEXIST) THEN
                  XDifValues(m,IDXDIF(m)) = XvarParas(m,n)%ParaValues(1)
                  IDXDIF(m) = IDXDIF(m) + 1
               END IF
            END DO
         END DO
         ! Range different X values
         DO m=1,NumTXvars
            DO k=1,IDXDIF(m)-1
               min = XDifValues(m,k)
               DO j=k+1,IDXDIF(m)-1
                  IF(XDifValues(m,j).LT.min) THEN
                     Temp = XDifValues(m,j)             
                     XDifValues(m,j) = XDifValues(m,k)             
                     XDifValues(m,k) = Temp             
                     min = XDifValues(m,k)
                  END IF
               END DO
            END DO
         END DO
         ! Put in terms
         DO jj=1,NumStates
            DO m=1,NumMacs
               DO n=1,MacParas(m,jj)%Length
                  CopyMacParas(m,jj)%ParaValues(n) = MacParas(m,jj)%ParaValues(n)
               END DO
            END DO
            DO m=1,LengthMics
               DO n=1,MicParas(m,jj)%Length
                  CopyMicParas(m,jj)%ParaValues(n) = MicParas(m,jj)%ParaValues(n)
               END DO
            END DO
         END DO
         ! Find for each subdivision which terms are typical ones
         ALLOCATE (BaseTerm(NumSubDs))
         IDSubDs = 1
         DO m=1,NumSubDs
            DO n=IDSubDs,IDSubDs+User_InputLink%NumTerms(m)-1
               IF(User_InputLink%Terms(n)%Operator.EQ.0) THEN
                  BaseTerm(m) = n 
               END IF
            END DO
            IDSubDs = IDSubDs + User_InputLink%NumTerms(m)
         END DO
         ! For used but not typical states, find their corresponding typical state numbers 
         IDPart = 1
         IDTS   = 1
         DO n=1,NumTterms
            IF(User_InputLink%Terms(n)%Operator.NE.0) THEN
               ! Find in which subdivision the term is
               IDSubDs = 1
               DO k=1,NumSubDs
                  IF((n.LE.IDSubDs+User_InputLink%NumTerms(k)-1).AND.(n.GE.IDSubDs)) THEN
                     TermInSub = k
                  END IF
                  IDSubDs = IDSubDs + User_InputLink%NumTerms(k)
               END DO
               ! Find out if the actual term has common X variables with the typical term
               ALLOCATE (IndexCommon(User_InputLink%Terms(n)%NumXvars))
               CommonX = .FALSE.
               DO j=1,User_InputLink%Terms(n)%NumXvars
                  IndexCommon(j) = .FALSE.
               END DO
               DO j=1,User_InputLink%Terms(n)%NumXvars
                  DO jj=1,User_InputLink%Terms(BaseTerm(TermInSub))%NumXvars
                     IF(User_InputLink%Terms(n)%IndexXvars(j).EQ.   &
                        User_InputLink%Terms(BaseTerm(TermInSub))%IndexXvars(jj)) THEN
                        CommonX = .TRUE.
                        IndexCommon(j) = .TRUE.
                        EXIT
                     END IF
                  END DO
               END DO
               ! Get the index of X variables which have not appeared in the actual term 
               LengthIR = NumTXvars-User_InputLink%Terms(n)%NumXvars
               ALLOCATE (IndexRest(LengthIR))
               IDIR = 1
               DO k=1,NumTXvars
                  INCLUDED = .FALSE.
                  DO j=1,User_InputLink%Terms(n)%NumXvars
                     IF(User_InputLink%Terms(n)%IndexXvars(j).EQ.k) THEN
                        INCLUDED = .TRUE.
                        EXIT
                     END IF
                  END DO
                  IF(.NOT.INCLUDED) THEN
                     IndexRest(IDIR) = k
                     IDIR = IDIR + 1
                  END IF
               END DO
               ! If the term has no common X variables with the typical term,  &
               ! to find typical states, one can simply compare X variables    &
               ! not contained in the actual term. If they get typical values, &
               ! the typical state is found.
               ! However, when the term get same X variables with the typical  &
               ! term, then the common X variables should have the same values. 
               IF(.NOT.CommonX) THEN
                  IDStateInv = 1
                  DO j=1,BaseTerm(TermInSub)-1
                     IDStateInv = IDStateInv + NumParts(j)
                  END DO
                  DO IDPartInv=IDPart,IDPart+NumParts(n)-1
                     DO m=TermStates(1,IDPartInv),TermStates(2,IDPartInv)
                        ! Find the index of typical states, and analyse one by one    & 
                        ! to see if it is the typical state for the state considered.
                        DO jj=IDStateInv,IDStateInv+NumParts(BaseTerm(TermInSub))-1
                           DO k=TermStates(1,jj),TermStates(2,jj)
                              XTYPICAL = .TRUE.
                              DO j=1,IDIR-1
                                 ID = User_InputLink%TypStates(IndexRest(j),TermInSub)
                                 IF(ABS(XDifValues(IndexRest(j),ID)-0.000).LT.NECP_RealZero) THEN
                                    IF(ABS(XvarParas(IndexRest(j),k)%ParaValues(1) -   &
                                       XDifValues(IndexRest(j),ID)).GT.1.0E-18) THEN
                                       XTYPICAL = .FALSE.
                                       EXIT
                                    END IF
                                 ELSE 
                                    IF(ABS((XvarParas(IndexRest(j),k)%ParaValues(1) -  &
                                       XDifValues(IndexRest(j),ID)) /                  &
                                       XDifValues(IndexRest(j),ID)).GT.1.0E-18) THEN
                                       XTYPICAL = .FALSE.
                                       EXIT
                                    END IF
                                 END IF
                                 ! write (*,*) XTYPICAL
                              END DO
                              IF(XTYPICAL) THEN
                                 TypicalStates(IDTS) = k
                              END IF
                           END DO                           
                        END DO
                        IDTS = IDTS + 1
                     END DO
                  END DO
               ELSE
                  IDStateInv = 1
                  DO j=1,BaseTerm(TermInSub)-1
                     IDStateInv = IDStateInv + NumParts(j)
                  END DO
                  DO IDPartInv=IDPart,IDPart+NumParts(n)-1
                     DO m=TermStates(1,IDPartInv),TermStates(2,IDPartInv)
                        ! Find the index of typical states, and analyse one by one     & 
                        ! to see if it is the typical state for the state considered.
                        DO jj=IDStateInv,IDStateInv+NumParts(BaseTerm(TermInSub))-1
                           DO k=TermStates(1,jj),TermStates(2,jj)
                              XTYPICAL = .TRUE.
                              DO j=1,IDIR-1
                                 ID = User_InputLink%TypStates(IndexRest(j),TermInSub)
                                 IF(ABS(XDifValues(IndexRest(j),ID)-0.000).LT.NECP_RealZero) THEN
                                    IF(ABS(XvarParas(IndexRest(j),k)%ParaValues(1) -   &
                                       XDifValues(IndexRest(j),ID)).GT.NECP_RealEqual) THEN
                                       XTYPICAL = .FALSE.
                                       EXIT
                                    END IF
                                 ELSE 
                                    IF(ABS((XvarParas(IndexRest(j),k)%ParaValues(1) -  &
                                       XDifValues(IndexRest(j),ID)) /                  &
                                       XDifValues(IndexRest(j),ID)).GT.1.0E-8) THEN
                                       XTYPICAL = .FALSE.
                                       EXIT
                                    END IF
                                 END IF
                              END DO
                              ! Now that the state has typical X values, one should consider &
                              ! if common X variables have the same values
                              IF(XTYPICAL) THEN
                                 EquValue = .TRUE.
                                 DO j=1,User_InputLink%Terms(StateSelected(m))%NumXvars
                                    IF(IndexCommon(j)) THEN
                                       ID = User_InputLink%Terms(StateSelected(m))%IndexXvars(j)
                                       IF(ABS(XvarParas(ID,k)%ParaValues(1)-0.000).LT. &
                                          NECP_RealZero) THEN
                                          IF(ABS(XvarParas(ID,k)%ParaValues(1)-   &
                                             XvarParas(ID,m)%ParaValues(1)).GT.(1.0E-18)) THEN
                                             EquValue = .FALSE.
                                          END IF
                                       ELSE
                                          IF(ABS((XvarParas(ID,k)%ParaValues(1)-  &
                                             XvarParas(ID,m)%ParaValues(1))/XvarParas(ID,k)%  &
                                             ParaValues(1)).GT.(1.0E-3)) THEN
                                             EquValue = .FALSE.
                                          END IF
                                       END IF
                                    END IF
                                 END DO
                              END IF
                              IF(XTYPICAL.AND.EquValue) THEN
                                 TypicalStates(IDTS) = k
                              END IF
                           END DO    
                        END DO
                        ! OPEN(1,FILE='TEMP')
                        ! CLOSE(1)
                        IDTS = IDTS + 1
                     END DO
                  END DO
               END IF
               DEALLOCATE (IndexCommon)
               DEALLOCATE (IndexRest)
            END IF
            IDPart = IDPart + NumParts(n)
            !pause
         END DO
         ! Get the total number of states in each term
         ALLOCATE (NStatesInTerm(NumTterms))
         DO jj=1,NumTterms
            NStatesInTerm(jj) = 0
         END DO
         IDTerm = 1
         IDPart = 1
         DO j=1,NumSubDs
            DO jj=IDTerm,IDTerm+User_InputLink%NumTerms(j)-1
               DO n=IDPart,IDPart+NumParts(jj)-1
                  NStatesInTerm(jj) = NStatesInTerm(jj) + TermStates(2,n) - TermStates(1,n) + 1
               END DO
               IDPart = IDPart + NumParts(jj)
            END DO
            IDTerm = IDTerm + User_InputLink%NumTerms(j)
         END DO
         ! Write into link input files
         IDTerm = 1
         IDPart = 1
         LINFileName = "LinkInputFiles/INPA00000000S00000000"
         WRITE (InputELALL,*) "*****************************"
         WRITE (InputELALL,"(I4)") User_InputLink%NumTXvars
         WRITE (InputELALL,*) "*****************************"
         DO j=1,NumSubDs
            WRITE (LINFileName(20:27),"(I8.8)") i
            WRITE (LINFileName(29:36),"(I8.8)") j
            OPEN  (INunitSEG,file=LINFileName)
            WRITE (INunitSEG,*) "************************************************************"
            WRITE (INunitSEG,"(A33,A14,I4)") FuelName(i),"Sub-Division",j
            WRITE (INunitSEG,*) "************************************************************"
            WRITE (INunitSEG,*) "Term Information"
            WRITE (INunitSEG,*) "************************************************************"
            YLength = 0 
            DO jj=1,NumMacs
               YLength = YLength + MacParas(jj,1)%Length
            END DO
            IDMicPara = 1
            DO m=1,NumMicNs
               DO kk=IDMicPara,IDMicPara+NumMicParas(m)-2
                  YLength = YLength + MicParas(kk,1)%Length
               END DO
               IDMicPara = IDMicPara + NumMicParas(m)
            END DO
            WRITE (INunitSEG,"(3(2X,I6))") User_InputLink%NumTerms(j),User_InputLink%NumTXvars, &
                                        YLength
            WRITE (INunitSEG,'(5(2X,A8))') (User_InputLink%XvarNames(k),k=1,NumTXvars)
            DO m=1,NumMacs
               WRITE (INunitSEG,'(5(2X,A8))') (MacParas(m,1)%ParaName,jj=1,MacParas(m,1)%Length)
            END DO
            IDMicPara = 1
            DO m=1,NumMicNs
               DO kk=IDMicPara,IDMicPara+NumMicParas(m)-2
                  WRITE (INunitSEG,'(5(2X,A8))') (MicParas(kk,1)%ParaName,jj=1,MicParas(kk,1)%Length)
               END DO
               IDMicPara = IDMicPara + NumMicParas(m)
            END DO
            DO jj=IDTerm,IDTerm+User_InputLink%NumTerms(j)-1
               WRITE (INunitSEG,'(3(2X,I6))') jj,User_InputLink%Terms(jj)%Parentheses(1), &
                                           User_InputLink%Terms(jj)%Operator,          &
                                           User_InputLink%Terms(jj)%TOption
            END DO
            DO m=1,User_InputLink%NumTXvars
               WRITE (INunitSEG,'(2(2X,ES15.8))') XDifValues(m,User_InputLink%SubDStates(1,m,j)), &
                                                  XDifValues(m,User_InputLink%SubDStates(2,m,j))
               WRITE (InputELALL,'(2(2X,ES15.8))') XDifValues(m,User_InputLink%SubDStates(1,m,j)),&
                                                   XDifValues(m,User_InputLink%SubDStates(2,m,j))
            END DO
            WRITE (InputELALL,*) "*************************" 
            WRITE (INunitSEG,*) "************************************************************"
            WRITE (INunitSEG,*) "In Each Term"
            WRITE (INunitSEG,*) "************************************************************"
            ! Write cross sections treated in each term in link input files
            DO jj=IDTerm,IDTerm+User_InputLink%NumTerms(j)-1
               ! X variable values
               WRITE (INunitSEG,"(2(I4))") User_InputLink%Terms(jj)%NumXvars,NStatesInTerm(jj)
               WRITE (INunitSEG,*) "************************************************************"
               WRITE (INunitSEG,"(3(I4))") (User_InputLink%Terms(jj)%Option(m),     &
                                            m=1,User_InputLink%Terms(jj)%NumXvars)
               WRITE (INunitSEG,"(3(I4))") (User_InputLink%Terms(jj)%PolyOrder(m),  &
                                            m=1,User_InputLink%Terms(jj)%NumXvars)
               WRITE (INunitSEG,"(3(I4))") (User_InputLink%Terms(jj)%IndexXvars(m), &
                                            m=1,User_InputLink%Terms(jj)%NumXvars)
               WRITE (INunitSEG,*) "************************************************************"
               ALLOCATE (Xbar(NStatesInTerm(jj)))
               ALLOCATE (Xav(User_InputLink%Terms(jj)%NumXvars))
               ALLOCATE (STD(User_InputLink%Terms(jj)%NumXvars))
               Xbar(:) = 0.0D0
               DO m=1,User_InputLink%Terms(jj)%NumXvars
                  Xav(m) = 0.0D0
                  STD(m) = 0.0D0
                  DO n=IDPart,IDPart+NumParts(jj)-1                     
                     DO k=TermStates(1,n),TermStates(2,n)
                        Xav(m) = Xav(m) + XvarParas(User_InputLink%Terms(jj)%IndexXvars(m),k)%  &
                              ParaValues(1)
                     END DO
                  END DO
                  Xav(m) = Xav(m) / NStatesInTerm(jj)
                  DO n=IDPart,IDPart+NumParts(jj)-1                     
                     DO k=TermStates(1,n),TermStates(2,n)
                        STD(m) = STD(m) + (XvarParas(User_InputLink%Terms(jj)%IndexXvars(m),k)% &
                                     ParaValues(1)-Xav(m))**2
                     END DO
                  END DO
                  STD(m) = SQRT(STD(m)/NStatesInTerm(jj))
                  IDXbar = 1
                  DO n=IDPart,IDPart+NumParts(jj)-1                     
                     DO k=TermStates(1,n),TermStates(2,n)
                        IF(ABS(STD(m)-0.0D0).GT.1.0E-18) THEN
                           Xbar(IDXbar) = (XvarParas(User_InputLink%Terms(jj)%  &
                           IndexXvars(m),k)%ParaValues(1)-Xav(m)) / STD(m)
                        ELSE
                           Xbar(IDXbar) = (XvarParas(User_InputLink%Terms(jj)%  &
                           IndexXvars(m),k)%ParaValues(1)-Xav(m))
                        END IF
                        IDXbar = IDXbar + 1
                     END DO
                  END DO
                  WRITE(INunitSEG,100) (Xbar(k),k=1,NStatesInTerm(jj))
                  WRITE(INunitSEG,*) "*****************************"
               END DO
               DO m=1,User_InputLink%Terms(jj)%NumXvars
                  WRITE(INunitSEG,'(2(2X,ES15.8))') Xav(m),STD(m)
               END DO
               WRITE(INunitSEG,*) "*****************************"
               WRITE(INunitSEG,*) "*****************************"
               DEALLOCATE (Xbar,Xav,STD)
               ! Cross sections:
               ! For terms whose operator equals zero, write parameters directly. 
               ! For terms whose operator equals one, their parameters should minus     &
               ! corresponding typical state values.
               ! For terms whose operator equals three, their parameters should divise  &
               ! corresponding typical state values. 
               ! No raisonable solutions for other operators at the moment.
               IF(User_InputLink%Terms(jj)%Operator.EQ.0) THEN
                  DO m=1,NumMacs
                     DO kk=1,MacParas(m,1)%Length
                        DO n=IDPart,IDPart+NumParts(jj)-1
                           WRITE(INunitSEG,100) (MacParas(m,k)%ParaValues(kk),     &
                                                 k=TermStates(1,n),TermStates(2,n))
                        END DO
                        WRITE(INunitSEG,*) "*****************************"
                     END DO
                  END DO
                  IDMicPara = 1
                  DO m=1,NumMicNs
                     DO kk=IDMicPara,IDMicPara+NumMicParas(m)-2
                        DO mm=1,MicParas(kk,1)%Length
                           DO n=IDPart,IDPart+NumParts(jj)-1
                              WRITE(INunitSEG,100) (MicParas(kk,k)%ParaValues(mm), &
                                                    k=TermStates(1,n),TermStates(2,n))
                           END DO
                           WRITE(INunitSEG,*) "*****************************"
                        END DO
                     END DO
                     IDMicPara = IDMicPara + NumMicParas(m)
                  END DO
               ELSE IF(User_InputLink%Terms(jj)%Operator.EQ.1) THEN
                  ! Write macroscopic cross sections
                  DO m=1,NumMacs
                     DO kk=1,MacParas(m,1)%Length
                        ! Get the index of typical state, IDTS
                        IDTS = 1
                        ID = 1
                        DO mm=1,jj-1
                           IF(User_InputLink%Terms(mm)%Operator.NE.0) THEN
                              DO nn=ID,ID+NumParts(mm)-1
                                 IDTS = IDTS + TermStates(2,nn) - TermStates(1,nn) + 1
                              END DO
                           END IF
                           ID = ID + NumParts(mm)
                        END DO
                        !if(m==1.and.kk==1) write (*,*) "seeIDTS,oper1",IDTS
                        DO n=IDPart,IDPart+NumParts(jj)-1
                           DO k=TermStates(1,n),TermStates(2,n)-1
                              WRITE(INunitSEG,200,ADVANCE="NO") MacParas(m,k)%ParaValues(kk) - &
                              MacParas(m,TypicalStates(IDTS))%ParaValues(kk)
                              IDTS = IDTS + 1
                           END DO
                           WRITE(INunitSEG,200) MacParas(m,k)%ParaValues(kk) -    &
                           MacParas(m,TypicalStates(IDTS))%ParaValues(kk)
                           IDTS = IDTS + 1
                        END DO
                        WRITE(INunitSEG,*) "*****************************"
                     END DO
                  END DO
                  ! Write microscopic cross sections
                  IDMicPara = 1
                  DO m=1,NumMicNs
                     DO kk=IDMicPara,IDMicPara+NumMicParas(m)-2
                        DO mm=1,MicParas(kk,1)%Length
                           ! Get the index of typical state, IDTS
                           IDTS = 1
                           ID = 1
                           DO k=1,jj-1
                              IF(User_InputLink%Terms(k)%Operator.NE.0) THEN
                                 DO nn=ID,ID+NumParts(k)-1
                                    IDTS = IDTS + TermStates(2,nn) - TermStates(1,nn) + 1
                                 END DO
                              END IF
                              ID = ID + NumParts(k)
                           END DO
                           DO n=IDPart,IDPart+NumParts(jj)-1
                              DO k=TermStates(1,n),TermStates(2,n)-1
                                 WRITE(INunitSEG,200,ADVANCE="NO") MicParas(kk,k)%ParaValues(mm) &
                                 - MicParas(kk,TypicalStates(IDTS))%ParaValues(mm)
                                 IDTS = IDTS + 1
                              END DO
                              WRITE(INunitSEG,200) MicParas(kk,k)%ParaValues(mm) -    &
                              MicParas(kk,TypicalStates(IDTS))%ParaValues(mm)
                              IDTS = IDTS + 1
                           END DO
                           WRITE(INunitSEG,*) "*****************************"
                        END DO
                     END DO
                     IDMicPara = IDMicPara + NumMicParas(m)
                  END DO
               ELSE IF(User_InputLink%Terms(jj)%Operator.EQ.3) THEN
                  FormatPrint = 1.0D0
                  ! Write macroscopic cross sections
                  DO m=1,NumMacs
                     DO kk=1,MacParas(m,1)%Length
                        ! Get the index of typical state, IDTS
                        IDTS = 1
                        ID = 1
                        DO mm=1,jj-1
                           IF(User_InputLink%Terms(mm)%Operator.NE.0) THEN
                              DO nn=ID,ID+NumParts(mm)-1
                                 IDTS = IDTS + TermStates(2,nn) - TermStates(1,nn) + 1
                              END DO
                           END IF
                           ID = ID + NumParts(mm)
                        END DO
                        DO n=IDPart,IDPart+NumParts(jj)-1
                           DO k=TermStates(1,n),TermStates(2,n)-1
                              IF(ABS(MacParas(m,TypicalStates(IDTS))%ParaValues(kk)-0.0D0) &
                                 .LT.1.0E-18) THEN
                                 WRITE(INunitSEG,200,ADVANCE="NO") FormatPrint
                              ELSE
                                 WRITE(INunitSEG,200,ADVANCE="NO") MacParas(m,k)%ParaValues(kk) / &
                                 MacParas(m,TypicalStates(IDTS))%ParaValues(kk)
                              END IF
                              IDTS = IDTS + 1
                           END DO
                           write (99,*) "test",jj,n,IDTS,TypicalStates(IDTS)
                           IF(ABS(MacParas(m,TypicalStates(IDTS))%ParaValues(kk)-0.0D0) &
                              .LT.1.0E-18) THEN ! Must be < to E-16 since KSigF is too small
                              WRITE(INunitSEG,200) FormatPrint
                           ELSE
                              WRITE(INunitSEG,200) MacParas(m,k)%ParaValues(kk) /    &
                              MacParas(m,TypicalStates(IDTS))%ParaValues(kk)
                           END IF
                           IDTS = IDTS + 1
                        END DO
                        WRITE(INunitSEG,*) "*****************************"
                     END DO
                  END DO
                  ! Write microscopic cross sections
                  IDMicPara = 1
                  DO m=1,NumMicNs
                     DO kk=IDMicPara,IDMicPara+NumMicParas(m)-2
                        DO mm=1,MicParas(kk,1)%Length
                           ! Get the index of typical state, IDTS
                           IDTS = 1
                           ID = 1
                           DO k=1,jj-1
                              IF(User_InputLink%Terms(k)%Operator.NE.0) THEN
                                 DO nn=ID,ID+NumParts(k)-1
                                    IDTS = IDTS + TermStates(2,nn) - TermStates(1,nn) + 1
                                 END DO
                              END IF
                              ID = ID + NumParts(k)
                           END DO
                           DO n=IDPart,IDPart+NumParts(jj)-1
                              DO k=TermStates(1,n),TermStates(2,n)-1
                                 IF(ABS(MicParas(kk,TypicalStates(IDTS))%ParaValues(mm)-0.0D0) &
                                    .LT.NECP_RealZero) THEN
                                    WRITE(INunitSEG,200,ADVANCE="NO") FormatPrint
                                 ELSE
                                    WRITE(INunitSEG,200,ADVANCE="NO") MicParas(kk,k)% &
                                    ParaValues(mm) / &
                                    MicParas(kk,TypicalStates(IDTS))%ParaValues(mm)
                                 END IF
                                 IDTS = IDTS + 1
                              END DO
                              IF(ABS(MicParas(kk,TypicalStates(IDTS))%ParaValues(mm)-0.0D0) &
                                 .LT.1.0E-18) THEN
                                 WRITE(INunitSEG,200) FormatPrint
                              ELSE
                                 WRITE(INunitSEG,200) MicParas(kk,k)%ParaValues(mm) /    &
                                 MicParas(kk,TypicalStates(IDTS))%ParaValues(mm)
                              END IF
                              IDTS = IDTS + 1
                           END DO
                           WRITE(INunitSEG,*) "*****************************"
                        END DO
                     END DO
                     IDMicPara = IDMicPara + NumMicParas(m)
                  END DO
               ELSE
                  WRITE (*,*) "[INPUT_LINK] ERROR : Sorry, no raisonable treatment for term ",jj                  
               END IF
               IDPart = IDPart + NumParts(jj)
            END DO
            CLOSE (INunitSEG)
            IDTerm = IDTerm + User_InputLink%NumTerms(j)
         END DO
         ! Deallocate all arrays in X/Mac/Mic parameters
         DO j=1,NumStates
            DO m=1,NumTXvars
               CALL ParaInfo_Void(XvarParas(m,j))
            END DO
            DO m=1,NumMacs
               CALL ParaInfo_Void(MacParas(m,j))
               CALL ParaInfo_Void(CopyMacParas(m,j))
            END DO
            DO m=1,LengthMics
               CALL ParaInfo_Void(MicParas(m,j))
               CALL ParaInfo_Void(CopyMicParas(m,j))
            END DO
         END DO
         DO j=1,NumTterms
            CALL TermInfo_Void(User_InputLink%Terms(j))
         END DO
         DEALLOCATE (TypicalStates)
         DEALLOCATE (XDifValues,IDXDIF)
         DEALLOCATE (TermStates,BaseTerm,NStatesInTerm)
         DEALLOCATE (StateSelected,TypeNs,NumMicParas,XvarParas,MacParas,MicParas,NumParts)
         DEALLOCATE (CopyMacParas,CopyMicParas)
         CALL InputLink_Void(User_InputLink)
         IDStartF = IDStartF + NumFiles(i)
      END DO
100   FORMAT (5(2X,ES15.8))  
200   FORMAT (2X,ES15.8)  
      DEALLOCATE (NumSInF,FileNames,FuelName,NumFiles)
      READ (INunitCR,*) Examine,RCoef
      CLOSE (INunitCR)
      CLOSE (InputELALL)
      close (99)
   END SUBROUTINE InputLink_ReadCR

   !-----------------------------------------------------------------------------------------------
   ! Read dragon output files
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE InputLink_ReadDragon(User_InputLink,FileNames,LengthFN,IDStartF,NumSInF,NFiles, &
                                   INunitXSF)
      ! Passed in
      INTEGER,INTENT(IN)      :: INunitXSF
      INTEGER,INTENT(IN)      :: NFiles
      INTEGER,INTENT(IN)      :: LengthFN
      INTEGER,INTENT(IN)      :: IDStartF
      CHARACTER(8),INTENT(IN) :: FileNames(LengthFN)
      INTEGER,INTENT(IN)      :: NumSInF(LengthFN)
      ! Passed out
      TYPE (InputLink_Type) :: User_InputLink
      ! Locals
      INTEGER :: j,jj,m,n,IDFileName,IDBG,ID
      INTEGER :: NumStatesFuel
      INTEGER :: XTLength,MacTLength,MicTLength
      INTEGER,ALLOCATABLE :: XLength(:,:)
      INTEGER,ALLOCATABLE :: MacLength(:,:)
      INTEGER,ALLOCATABLE :: MicLength(:,:)
      REAL(8),ALLOCATABLE :: XParaValues(:)
      REAL(8),ALLOCATABLE :: MacParaValues(:)
      REAL(8),ALLOCATABLE :: MicParaValues(:)
      REAL(8),ALLOCATABLE :: XParas(:),MacParas(:),MicParas(:)
 
      NumStatesFuel = 0
      DO j=1,NFiles
         IDFileName = IDStartF + 1
         !WRITE (*,*) IDBG,IDBG+NumSInF(IDFileName)-1
         NumStatesFuel = NumStatesFuel + NumSInF(IDFileName)
      END DO
      ALLOCATE (XLength(User_InputLink%NumTXvars,NumStatesFuel))
      ALLOCATE (MacLength(User_InputLink%NumMacs,NumStatesFuel))
      ALLOCATE (MicLength(User_InputLink%LengthMics,NumStatesFuel))

      CALL ParaInfo_GetLength(User_InputLink%XvarParas,User_InputLink%NumTXvars,NumStatesFuel,  &
                              XLength)
      CALL ParaInfo_GetLength(User_InputLink%MacParas,User_InputLink%NumMacs,NumStatesFuel,     &
                              MacLength)
      CALL ParaInfo_GetLength(User_InputLink%MicParas,User_InputLink%LengthMics,NumStatesFuel, &
                              MicLength)

      IDBG = 1
      DO j=1,NFiles
         IDFileName = IDStartF + 1
         !write (*,*) INunitXSF,IDFileName,FileNames(IDFileName)
         OPEN (INunitXSF,file=FileNames(IDFileName))
         !WRITE (*,*) IDBG,IDBG+NumSInF(IDFileName)-1
         DO jj=IDBG,IDBG+NumSInF(IDFileName)-1
            XTLength   = 0
            MacTLength = 0
            MicTLength = 0
            DO m=1,User_InputLink%NumTXvars
               XTLength   = XTLength + XLength(m,jj)
            END DO
            DO m=1,User_InputLink%NumMacs
               MacTLength = MacTLength + MacLength(m,jj)
            END DO
            DO m=1,User_InputLink%LengthMics
               MicTLength = MicTLength + MicLength(m,jj)
            END DO
            ALLOCATE (XParaValues(XTLength),MacParaValues(MacTLength),MicParaValues(MicTLength))
            READ (INunitXSF,*)  
            READ (INunitXSF,*) (XParaValues(n),n=1,XTLength)
            READ (INunitXSF,*) (MacParaValues(n),n=1,MacTLength)
            READ (INunitXSF,*) (MicParaValues(n),n=1,MicTLength)
            ID = 1
            DO m=1,User_InputLink%NumTXvars
               ALLOCATE (XParas(XLength(m,jj)))
               DO n=1,XLength(m,jj)
                  XParas(n) = XParaValues(ID+n-1)
               END DO
               CALL ParaInfo_SetValues(User_InputLink%XvarParas(m,jj),XParas,XLength(m,jj))
               DEALLOCATE (XParas)
               ID = ID + XLength(m,jj)
            END DO
            ID = 1
            DO m=1,User_InputLink%NumMacs
               ALLOCATE (MacParas(MacLength(m,jj)))
               DO n=1,MacLength(m,jj)
                  MacParas(n) = MacParaValues(ID+n-1)
               END DO
               CALL ParaInfo_SetValues(User_InputLink%MacParas(m,jj),MacParas,MacLength(m,jj))
               DEALLOCATE (MacParas)
               ID = ID + MacLength(m,jj)
            END DO
            ID = 1
            DO m=1,User_InputLink%LengthMics
               ALLOCATE (MicParas(MicLength(m,jj)))
               DO n=1,MicLength(m,jj)
                  MicParas(n) = MicParaValues(ID+n-1)
               END DO
               CALL ParaInfo_SetValues(User_InputLink%MicParas(m,jj),MicParas,MicLength(m,jj))
               DEALLOCATE (MicParas)
               ID = ID + MicLength(m,jj)
            END DO
            DEALLOCATE (XParaValues,MacParaValues,MicParaValues)
         END DO
         IDBG = IDBG + NumSInF(IDFileName)
         CLOSE (INunitXSF)
      END DO
      DEALLOCATE (XLength,MacLength,MicLength)
   END SUBROUTINE InputLink_ReadDragon

   !-----------------------------------------------------------------------------------------------
   ! Get AuxiDim & Length of ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_GetLength(User_ParaInfo,NumParas,NumStates,Length)
      ! Passed in 
      INTEGER,INTENT(IN)        :: NumParas,NumStates
      TYPE(ParaInfo),INTENT(IN) :: User_ParaInfo(NumParas,NumStates)
      ! Passed out
      INTEGER :: Length(NumParas,NumStates)    
      ! Locals
      INTEGER :: i,j

      DO i=1,NumParas
         DO j=1,NumStates
            Length(i,j)  = User_ParaInfo(i,j)%Length
         END DO
      END DO
   END SUBROUTINE ParaInfo_GetLength

   !-----------------------------------------------------------------------------------------------
   ! Set AuxiDim & Length of ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_SetLength(User_ParaInfo,AuxiDim,Length)
      ! Passed in 
      INTEGER,INTENT(IN) :: AuxiDim,Length     
      ! Passed out
      TYPE(ParaInfo)     :: User_ParaInfo
      User_ParaInfo%AuxiDim = AuxiDim 
      User_ParaInfo%Length  = Length 
   END SUBROUTINE ParaInfo_SetLength

   !-----------------------------------------------------------------------------------------------
   ! Set ParaValues of ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_SetValues(User_ParaInfo,Paras,ParaLength)
      ! Passed in
      INTEGER,INTENT(IN) :: ParaLength
      REAL(8),INTENT(IN) :: Paras(ParaLength)
      ! Passed out
      TYPE(ParaInfo)     :: User_ParaInfo
      ! Locals
      INTEGER :: I
      IF(User_ParaInfo%Length.NE.ParaLength) THEN
         WRITE (*,*) "[Link_Input] FATAL ERROR : Cannot set paravalues, given length incorrect"
         RETURN
      END IF
      DO I=1,ParaLength
         User_ParaInfo%ParaValues(I) = Paras(I) 
      END DO
   END SUBROUTINE ParaInfo_SetValues

   !-----------------------------------------------------------------------------------------------
   ! Set ParaValues of ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_GetValues(User_ParaInfo,Paras,ParaLength)
      ! Passed in
      TYPE(ParaInfo),INTENT(IN) :: User_ParaInfo
      INTEGER,INTENT(IN)        :: ParaLength
      ! Passed out
      REAL(8)                   :: Paras(ParaLength)
      ! Locals
      INTEGER :: I
      IF(User_ParaInfo%Length.NE.ParaLength) THEN
         WRITE (*,*) "[Link_Input] FATAL ERROR : Cannot get paravalues, given length incorrect"
         RETURN
      END IF
      DO I=1,ParaLength
         Paras(I) = User_ParaInfo%ParaValues(I) 
      END DO
   END SUBROUTINE ParaInfo_GetValues

   !-----------------------------------------------------------------------------------------------
   ! Define ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_Define(User_ParaInfo,AuxiDim,Length)
      IMPLICIT NONE
      ! Passed in
      INTEGER :: AuxiDim,Length
      ! Passed out
      TYPE (ParaInfo) :: User_ParaInfo
      IF((AuxiDim.LE.0).OR.(Length.LE.0)) THEN
         WRITE (*,*) "[Input_Link] FATAL ERROR : Cannot define User_ParaInfo, length incorrect"
         RETURN
      END IF
      IF(User_ParaInfo%Defined) CALL ParaInfo_Void(User_ParaInfo)
      User_ParaInfo%AuxiDim = AuxiDim
      User_ParaInfo%Length  = Length
      ALLOCATE (User_ParaInfo%IGstart(AuxiDim),User_ParaInfo%IGend(AuxiDim), &
                User_ParaInfo%ParaValues(Length))
      User_ParaInfo%Defined = .TRUE.
   END SUBROUTINE ParaInfo_Define

   !-----------------------------------------------------------------------------------------------
   ! Get storage of ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_GetStorage(User_ParaInfo,NumINT,NumREAL,NumLOG,NumCHR)
      IMPLICIT NONE
      ! Passed in
      TYPE(ParaInfo),INTENT(IN) :: User_ParaInfo
      ! Passed out
      INTEGER :: NumINT,NumREAL,NumLOG,NumCHR
      NumINT  = 3
      NumREAL = 0
      NumLOG  = 1
      NumCHR  = 0
      IF (User_ParaInfo%Defined) THEN
         NumREAL = NumREAL + User_ParaInfo%Length
         NumINT  = NumINT  + User_ParaInfo%AuxiDim*2
      END IF
   END SUBROUTINE ParaInfo_GetStorage

   !-----------------------------------------------------------------------------------------------
   ! Print objects of ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_Print(User_ParaInfo,OUTunit)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)           :: OUTunit
      TYPE(ParaInfo),INTENT(IN)    :: User_ParaInfo
      ! Locals       
      INTEGER                      :: I,J
      IF(User_ParaInfo%Defined) THEN
         WRITE(OUTunit,10)
         WRITE(OUTunit,20)
         WRITE(OUTunit,90)  User_ParaInfo%Defined
         WRITE(OUTunit,30)  User_ParaInfo%ParaType,User_ParaInfo%AuxiDim,User_ParaInfo%Length
         WRITE(OUTunit,40) 
         WRITE(OUTunit,50)  (User_ParaInfo%IGstart(i),i=1,User_ParaInfo%AuxiDim)
         WRITE(OUTunit,80)
         WRITE(OUTunit,50)  (User_ParaInfo%IGend(i),i=1,User_ParaInfo%AuxiDim)
         WRITE(OUTunit,100)
         WRITE(OUTunit,110) (User_ParaInfo%ParaValues(i),i=1,User_ParaInfo%Length)
         WRITE(OUTunit,10)         
      ELSE
         WRITE(*,60)
         WRITE(OUTunit,60)
         WRITE(*,70)
         WRITE(OUTunit,70)
      END IF
10    FORMAT("[Link_Input] ",90("-")) 
20    FORMAT("[Link_Input] Object User_ParaInfo")
30    FORMAT("[Link_Input] ParaType,AuxiDim,Length:",I,I,I)
40    FORMAT("[Link_Input] IGstart")
50    FORMAT("[Link_Input]",100(2X,I))
60    FORMAT("[Link_Input] FATAL ERROR !!!")
70    FORMAT("[Link_Input] Sorry, UNDEFINED object User_ParaInfo.")
80    FORMAT("[Link_Input] IGend")
90    FORMAT("[Link_Input] User_ParaInfo defined?",L)
100   FORMAT("[Link_Input] ParaValues")
110   FORMAT("[Link]",100(2X,1PE16.6))
   END SUBROUTINE ParaInfo_Print

   !-----------------------------------------------------------------------------------------------
   ! Read objects of ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_Read(User_ParaInfo,INunit,OUTunit)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)           :: INunit,OUTunit
      TYPE (ParaInfo),INTENT(IN)   :: User_ParaInfo
      ! Locals       
      INTEGER                      :: i,j
      INTEGER                      :: ParaType,AuxiDim,Length
      LOGICAL                      :: Defined = .FALSE.
      READ (INunit,*)
      READ (INunit,*)
      READ (INunit,90) Defined
      IF(Defined) THEN
         READ (INunit,30) ParaType,AuxiDim,Length
         CALL ParaInfo_Define(User_ParaInfo,AuxiDim,Length)
         READ (INunit,*)
         READ (INunit,50)  (User_ParaInfo%IGstart(i),i=1,User_ParaInfo%AuxiDim)
         READ (INunit,*)
         READ (INunit,50)  (User_ParaInfo%IGend(i),i=1,User_ParaInfo%AuxiDim)
         READ (INunit,*)
         READ (INunit,110) (User_ParaInfo%ParaValues(i),i=1,User_ParaInfo%Length)
         READ (INunit,*)         
      ELSE
         WRITE(*,130)
         WRITE(OUTunit,130)
         WRITE(*,140)
         WRITE(OUTunit,140)
         WRITE(*,150)
         WRITE(OUTunit,150)
      END IF
30    FORMAT("[Link_Input] ParaType,AuxiDim,Length:",I,I,I)
50    FORMAT("[Link_Input]",100(2X,I))
90    FORMAT("[Link] User_ParaInfo defined?",L)
110   FORMAT("[Link]",100(2X,1PE16.6))
130   FORMAT("[Link] Non-fatal error !")   
140   FORMAT("[Link] UNDEFINED object User_Link in file INunit !")   
150   FORMAT("[Link] Object User_Link stays UNDEFINED !")   
    END SUBROUTINE ParaInfo_Read

   !-----------------------------------------------------------------------------------------------
   ! Void ParaInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE ParaInfo_Void(User_ParaInfo)
      IMPLICIT NONE
      ! Passed out
      TYPE (ParaInfo) :: User_ParaInfo
      IF(User_ParaInfo%Defined) THEN
         DEALLOCATE (User_ParaInfo%IGstart,User_ParaInfo%IGend,User_ParaInfo%ParaValues)
         User_ParaInfo%AuxiDim = 0
         User_ParaInfo%Length  = 0
         User_ParaInfo%Defined = .FALSE.
      END IF
   END SUBROUTINE ParaInfo_Void

   !-----------------------------------------------------------------------------------------------
   ! Define TermInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE TermInfo_Define(User_TermInfo,NumXvars)
      IMPLICIT NONE
      ! Passed in 
      INTEGER :: NumXvars
      ! Passed out
      TYPE (TermInfo) :: User_TermInfo
      IF(NumXvars.LE.0) THEN
         WRITE (*,*) "[Input_Link] FATAL ERROR : Cannot define User_TermInfo, length incorrect"
         RETURN
      END IF
      IF(User_TermInfo%Defined) CALL TermInfo_Void(User_TermInfo)
      User_TermInfo%NumXvars = NumXvars
      ALLOCATE (User_TermInfo%IndexXvars(NumXvars),User_TermInfo%Option(NumXvars), &
                User_TermInfo%PolyOrder(NumXvars))
      User_TermInfo%Defined = .TRUE.
   END SUBROUTINE TermInfo_Define

   !-----------------------------------------------------------------------------------------------
   ! Get storage of TermInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE TermInfo_GetStorage(User_TermInfo,NumINT,NumREAL,NumLOG,NumCHR)
      IMPLICIT NONE
      ! Passed in
      TYPE(TermInfo),INTENT(IN) :: User_TermInfo
      ! Passed out
      INTEGER :: NumINT,NumREAL,NumLOG,NumCHR
      NumINT  = 3
      NumREAL = 0
      NumLOG  = 1
      NumCHR  = 0
      IF(User_TermInfo%Defined) THEN
         NumINT  = NumINT + User_TermInfo%NumXvars*3
      END IF
   END SUBROUTINE TermInfo_GetStorage

   !-----------------------------------------------------------------------------------------------
   ! Print objects of TermInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE TermInfo_Print(User_TermInfo,OUTunit)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)           :: OUTunit
      TYPE(TermInfo),INTENT(IN)    :: User_TermInfo
      ! Locals       
      INTEGER                      :: I,J
      IF(User_TermInfo%Defined) THEN
         WRITE(OUTunit,10)
         WRITE(OUTunit,20)
         WRITE(OUTunit,90)  User_TermInfo%Defined
         WRITE(OUTunit,30)  User_TermInfo%Operator,User_TermInfo%TOption,User_TermInfo%NumXvars
         WRITE(OUTunit,40) 
         WRITE(OUTunit,50)  (User_TermInfo%IndexXvars(i),i=1,User_TermInfo%NumXvars)
         WRITE(OUTunit,80)
         WRITE(OUTunit,50)  (User_TermInfo%Option(i),i=1,User_TermInfo%NumXvars)
         WRITE(OUTunit,100)
         WRITE(OUTunit,50)  (User_TermInfo%PolyOrder(i),i=1,User_TermInfo%NumXvars)
         WRITE(OUTunit,10)         
      ELSE
         WRITE(*,60)
         WRITE(OUTunit,60)
         WRITE(*,70)
         WRITE(OUTunit,70)
      END IF
10    FORMAT("[Link_Input] ",90("-")) 
20    FORMAT("[Link_Input] Object User_TermInfo")
30    FORMAT("[Link_Input] Operator,TOption,NumXvars:",I,I,I)
40    FORMAT("[Link_Input] IndexXvars")
50    FORMAT("[Link_Input]",100(2X,I))
60    FORMAT("[Link_Input] FATAL ERROR !!!")
70    FORMAT("[Link_Input] Sorry, UNDEFINED object User_TermInfo.")
80    FORMAT("[Link_Input] Option")
90    FORMAT("[Link_Input] User_ParaInfo defined?",L)
100   FORMAT("[Link_Input] PolyOrder")
110   FORMAT("[Link_Input]",100(2X,1PE16.6))
   END SUBROUTINE TermInfo_Print

   !-----------------------------------------------------------------------------------------------
   ! Read objects of TermInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE TermInfo_Read(User_TermInfo,INunit,OUTunit)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)           :: INunit,OUTunit
      TYPE(TermInfo),INTENT(IN)    :: User_TermInfo
      ! Locals       
      INTEGER                      :: i,j
      INTEGER                      :: Operator,TOption,NumXvars
      LOGICAL                      :: Defined = .FALSE.
      READ (INunit,*)
      READ (INunit,*)
      READ (INunit,90) Defined
      IF(Defined) THEN
         READ (INunit,30) Operator,TOption,NumXvars
         CALL TermInfo_Define(User_TermInfo,NumXvars)
         READ (INunit,*)
         READ (INunit,50)  (User_TermInfo%IndexXvars(i),i=1,User_TermInfo%NumXvars)
         READ (INunit,*)
         READ (INunit,50)  (User_TermInfo%Option(i),i=1,User_TermInfo%NumXvars)
         READ (INunit,*)
         READ (INunit,50)  (User_TermInfo%PolyOrder(i),i=1,User_TermInfo%NumXvars)
         READ (INunit,*)         
      ELSE
         WRITE(*,130)
         WRITE(OUTunit,130)
         WRITE(*,140)
         WRITE(OUTunit,140)
         WRITE(*,150)
         WRITE(OUTunit,150)
      END IF
30    FORMAT("[Link_Input] ParaType,AuxiDim,Length:",I,I,I)
50    FORMAT("[Link_Input]",100(2X,I))
90    FORMAT("[Link_Input] User_ParaInfo defined?",L)
130   FORMAT("[Link_Input] Non-fatal error !")   
140   FORMAT("[Link_Input] UNDEFINED object User_Link in file INunit !")   
150   FORMAT("[Link_Input] Object User_Link stays UNDEFINED !")   
   END SUBROUTINE TermInfo_Read

   !-----------------------------------------------------------------------------------------------
   ! Set viables of TermInfo type objects
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE TermInfo_SetTerm(User_TermInfo,Operator,Parentheses,TOption, &
                               IndexXvars,Option,PolyOrder,NumXvars)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN) :: Operator,TOption
      INTEGER,INTENT(IN) :: Parentheses(2)
      INTEGER,INTENT(IN) :: NumXvars
      INTEGER,INTENT(IN) :: IndexXvars(NumXvars),Option(NumXvars),PolyOrder(NumXvars)
      ! Passed out
      TYPE (TermInfo)    :: User_TermInfo
      ! Locals
      INTEGER            :: I
      IF(NumXvars.NE.User_TermInfo%NumXvars) THEN
         WRITE (*,*) '[Link_Input] FATAL ERROR : Cannot set term, NumXvars incorrect'
         RETURN
      END IF
      User_TermInfo%Operator = Operator
      User_TermInfo%TOption = TOption
      User_TermInfo%Parentheses(1) = Parentheses(1)
      User_TermInfo%Parentheses(2) = Parentheses(2)
      DO I=1,NumXvars
         User_TermInfo%IndexXvars(I) = IndexXvars(I)
         User_TermInfo%Option(I)     = Option(I)
         User_TermInfo%PolyOrder(I)  = PolyOrder(I)
      END DO
   END SUBROUTINE TermInfo_SetTerm

   !-----------------------------------------------------------------------------------------------
   ! Void TermInfo type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE TermInfo_Void(User_TermInfo)   
      IMPLICIT NONE
      ! Passed out
      TYPE (TermInfo) :: User_TermInfo
      IF(User_TermInfo%Defined) THEN
         DEALLOCATE (User_TermInfo%IndexXvars,User_TermInfo%Option,User_TermInfo%PolyOrder)
         User_TermInfo%NumXvars = 0
         User_TermInfo%Defined  = .FALSE.
      END IF
   END SUBROUTINE TermInfo_Void

   !-----------------------------------------------------------------------------------------------
   ! Define control parameters of InputLink type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE InputLink_DefineCR(User_InputLink,NumTXvars,NumSubDs,NumTterms)
      IMPLICIT NONE
      ! Passed in
      INTEGER :: NumTXvars,NumSubDs,NumTterms
      ! Passed out
      TYPE (InputLink_Type) :: User_InputLink
      IF((NumTXvars.LE.0).OR.(NumSubDs.LE.0).OR.(NumTterms.LE.0)) THEN
         WRITE (*,*) "[Input_Link] FATAL ERROR : Cannot define User_InputLink CR, length incorrect"
         RETURN
      END IF
      IF(User_InputLink%CRDefined) CALL InputLink_Void(User_InputLink) 
      IF((User_InputLink%XSDefined).AND.(User_InputLink%NumTXvars.NE.NumTXvars)) THEN
         WRITE (*,*) "[Input_Link] FATAL ERROR : Cannot define User_InputLink CR, &
                      NumTXvars incorrect"
         RETURN
      END IF
      User_InputLink%NumTXvars = NumTXvars
      User_InputLink%NumSubDs  = NumSubDs
      User_InputLink%NumTterms = NumTterms
      ALLOCATE (User_InputLink%XvarNames(NumTXvars),             &
                User_InputLink%SubDStates(2,NumTXvars,NumSubDs), &
                User_InputLink%TypStates(NumTXvars,NumSubDs),    &
                User_InputLink%NumTerms(NumSubDs),               &
                User_InputLink%Terms(NumTterms))
      User_InputLink%CRDefined = .TRUE.
   END SUBROUTINE InputLink_DefineCR

   !-----------------------------------------------------------------------------------------------
   ! Define XS parameters of InputLink type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE InputLink_DefineXS(User_InputLink,NumTXvars,NumStates,NumMacs,NumMicNs,LengthMics)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN) :: NumTXvars,NumStates,NumMacs,NumMicNs,LengthMics
      ! Passed out
      TYPE (InputLink_Type) :: User_InputLink
      ! Locals
      IF((NumTXvars.LE.0).OR.(NumStates.LE.0).OR.(NumMacs.LE.0) &
         .OR.(NumMicNs.LE.0).OR.(LengthMics.LE.0)) THEN
         WRITE (*,*) "[Input_Link] FATAL ERROR : Cannot define User_InputLink XS, length incorrect"
         RETURN
      END IF
      IF(User_InputLink%XSDefined) CALL InputLink_Void(User_InputLink) 
      IF((User_InputLink%CRDefined).AND.(User_InputLink%NumTXvars.NE.NumTXvars)) THEN
         WRITE (*,*) "[Input_Link] FATAL ERROR : Cannot define User_InputLink XS, &
                      NumTXvars incorrect"
         RETURN
      END IF
      User_InputLink%NumTXvars  = NumTXvars
      User_InputLink%NumStates  = NumStates
      User_InputLink%NumMacs    = NumMacs
      User_InputLink%NumMicNs   = NumMicNs
      User_InputLink%LengthMics = LengthMics
      ALLOCATE (User_InputLink%StateSelected(NumStates),       &
                User_InputLink%XvarParas(NumTXvars,NumStates), &
                User_InputLink%MacParas(NumMacs,NumStates),    &
                User_InputLink%TypeNs(NumMicNs),               &
                User_InputLink%NumMicParas(NumMicNs),          &
                User_InputLink%MicParas(LengthMics,NumStates))
      User_InputLink%XSDefined = .TRUE.
   END SUBROUTINE InputLink_DefineXS

   !-----------------------------------------------------------------------------------------------
   ! Get storage of InputLink type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE InputLink_GetStorage(User_InputLink,NumINT,NumREAL,NumLOG,NumCHR)
      IMPLICIT NONE
      ! Passed in
      TYPE(InputLink_Type),INTENT(IN) :: User_InputLink
      ! Passed out
      INTEGER :: NumINT,NumREAL,NumLOG,NumCHR
      ! Locals
      INTEGER :: I,J
      INTEGER :: NINT,NREAL,NLOG,NCHR
      NumCHR  = 0
      NumINT  = 8
      NumLOG  = 2
      NumREAL = 0
      NINT    = 0
      NREAL   = 0
      NLOG    = 0
      NCHR    = 0
      IF(User_InputLink%CRDefined) THEN
         NumINT  = NumINT + User_InputLink%NumTXvars +                  &
                   User_InputLink%NumTXvars*User_InputLink%NumSubDs*3 + &
                   User_InputLink%NumSubDs + User_InputLink%NumTterms
         DO I=1,User_InputLink%NumTterms
            CALL TermInfo_GetStorage(User_InputLink%Terms(I),NINT,NREAL,NLOG,NCHR)
            NumINT  = NumINT + NINT
            NumREAL = NumREAL + NREAL
            NumLOG  = NumLOG + NLOG
            NumCHR = NumCHR + NCHR
         END DO
      END IF
      IF(User_InputLink%XSDefined) THEN
         NumINT = NumINT + User_InputLink%NumStates + User_InputLink%NumMicNs*2
         DO I=1,User_InputLink%NumStates
            DO J=1,User_InputLink%NumTXvars
               CALL ParaInfo_GetStorage(User_InputLink%XvarParas(J,I),NINT,NREAL,NLOG,NCHR)
               NumINT  = NumINT + NINT
               NumREAL = NumREAL + NREAL
               NumLOG  = NumLOG + NLOG
               NumCHR = NumCHR + NCHR
            END DO
            DO J=1,User_InputLink%NumMacs
               CALL ParaInfo_GetStorage(User_InputLink%MacParas(J,I),NINT,NREAL,NLOG,NCHR)
               NumINT  = NumINT + NINT
               NumREAL = NumREAL + NREAL
               NumLOG  = NumLOG + NLOG
               NumCHR = NumCHR + NCHR
            END DO
            DO J=1,User_InputLink%LengthMics
               CALL ParaInfo_GetStorage(User_InputLink%MicParas(J,I),NINT,NREAL,NLOG,NCHR)
               NumINT  = NumINT + NINT
               NumREAL = NumREAL + NREAL
               NumLOG  = NumLOG + NLOG
               NumCHR = NumCHR + NCHR
            END DO
         END DO
      END IF
   END SUBROUTINE InputLink_GetStorage

   !-----------------------------------------------------------------------------------------------
   ! Print objects of InputLink type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE InputLink_Print(User_InputLink,OUTunit)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)              :: OUTunit
      TYPE(InputLink_Type),INTENT(IN) :: User_InputLink
      ! Locals       
      INTEGER :: i,j,k
      WRITE(OUTunit,10)
      WRITE(OUTunit,20)
      WRITE(OUTunit,90)  User_InputLink%CRDefined
      WRITE(OUTunit,150) User_InputLink%XSDefined
      IF((User_InputLink%CRDefined).OR.(User_InputLink%XSDefined)) THEN
         WRITE(OUTunit,160) User_InputLink%NumTXvars        
      END IF
      IF(User_InputLink%CRDefined) THEN
         WRITE(OUTunit,30)  User_InputLink%NumSubDs,User_InputLink%NumTterms
         WRITE(OUTunit,40) 
         WRITE(OUTunit,50)  (User_InputLink%XvarNames(i),i=1,User_InputLink%NumTXvars)
         WRITE(OUTunit,80)
         WRITE(OUTunit,50)  (((User_InputLink%SubDStates(i,j,k),i=1,2),       &
                               j=1,User_InputLink%NumTXvars),                 &
                               k=1,User_InputLink%NumSubDs)
         WRITE(OUTunit,100)
         WRITE(OUTunit,120) ((User_InputLink%TypStates(i,j),i=1,User_InputLink%NumTXvars), &
                              j=1,User_InputLink%NumSubDs)
         WRITE(OUTunit,130)
         WRITE(OUTunit,120) (User_InputLink%NumTerms(i),i=1,User_InputLink%NumSubDs)
         DO i=1,User_InputLink%NumTterms
            WRITE(OUTunit,130) i            
            CALL TermInfo_Print(User_InputLink%Terms(i),OUTunit)
         END DO
         WRITE(OUTunit,10)         
      END IF
      IF(User_InputLink%XSDefined) THEN
         WRITE(OUTunit,60)  User_InputLink%NumEGrps,User_InputLink%NumStates,User_InputLink%NumMacs
         WRITE(OUTunit,70)  User_InputLink%NumMicNs,User_InputLink%LengthMics
         WRITE(OUTunit,170)  
         WRITE(OUTunit,120) (User_InputLink%StateSelected(i),i=1,User_InputLink%NumStates)
         WRITE(OUTunit,180)
         WRITE(OUTunit,120) (User_InputLink%TypeNs(i),i=1,User_InputLink%NumMicNs)
         WRITE(OUTunit,190)
         WRITE(OUTunit,120) (User_InputLink%NumMicParas(i),i=1,User_InputLink%NumMicNs)
         DO i=1,User_InputLink%NumStates
            DO j=1,User_InputLink%NumTXvars
               CALL ParaInfo_Print(User_InputLink%XvarParas(j,i),OUTunit)
            END DO
            DO j=1,User_InputLink%NumMacs
               CALL ParaInfo_Print(User_InputLink%MacParas(j,i),OUTunit)
            END DO
            DO j=1,User_InputLink%LengthMics
               CALL ParaInfo_Print(User_InputLink%MicParas(j,i),OUTunit)
            END DO
         END DO
      END IF
10    FORMAT("[Link_Input] ",90("-")) 
20    FORMAT("[Link_Input] Object User_InputLink")
30    FORMAT("[Link_Input] NumSubDs,NumTterms:",I,I)
40    FORMAT("[Link_Input] XvarNames")
50    FORMAT("[Link_Input]",100(2X,A))
60    FORMAT("[Link_Input] NumEGrps,NumStates,NumMacs:",I,I,I)
70    FORMAT("[Link_Input] NumMicNs,LengthMics:",I,I)
80    FORMAT("[Link_Input] SubDStates")
90    FORMAT("[Link_Input] User_InputLink CR defined?",L)
100   FORMAT("[Link_Input] TypStates")
110   FORMAT("[Link_Input]",100(2X,1PE16.6))
120   FORMAT("[Link_Input]",100(2X,I))
130   FORMAT("[Link_Input] NumTerms")
140   FORMAT("[Link_Input] Term",I)
150   FORMAT("[Link_Input] User_InputLink XS defined?",L)
160   FORMAT("[Link_Input] NumTXvars",I)
170   FORMAT("[Link_Input] StateSelected")
180   FORMAT("[Link_Input] TypeNs")
190   FORMAT("[Link_Input] NumMicParas")
   END SUBROUTINE InputLink_Print

   !-----------------------------------------------------------------------------------------------
   ! Read objects of InputLink type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE InputLink_Read(User_InputLink,INunit,OUTunit)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)              :: INunit,OUTunit
      ! Passed out
      TYPE(InputLink_Type) :: User_InputLink
      ! Locals       
      LOGICAL :: CRDefined = .FALSE.
      LOGICAL :: XSDefined = .FALSE.
      INTEGER :: NumTXvars,NumSubDs,NumTterms
      INTEGER :: NumStates,NumMacs,NumMicNs,LengthMics,NumEGrps
      INTEGER :: i,j,k
      READ(INunit,*)
      READ(INunit,*)
      READ(INunit,*) CRDefined
      READ(INunit,*) XSDefined
      IF((User_InputLink%CRDefined).OR.(User_InputLink%XSDefined)) THEN
         READ(INunit,160) NumTXvars        
      ELSE 
         WRITE(*,130)
         WRITE(OUTunit,130)
         WRITE(*,140)
         WRITE(OUTunit,140)
         WRITE(*,170)
         WRITE(OUTunit,170)         
      END IF
      IF(User_InputLink%CRDefined) THEN
         READ(INunit,30)  NumSubDs,NumTterms
         CALL InputLink_DefineCR(User_InputLink,NumTXvars,NumSubDs,NumTterms)
         READ(INunit,*) 
         READ(INunit,50)  (User_InputLink%XvarNames(i),i=1,User_InputLink%NumTXvars)
         READ(INunit,*)
         READ(INunit,50)  (((User_InputLink%SubDStates(i,j,k),i=1,2),       &
                             j=1,User_InputLink%NumTXvars),                 &
                             k=1,User_InputLink%NumSubDs)
         READ(INunit,*)
         READ(INunit,120) ((User_InputLink%TypStates(i,j),i=1,User_InputLink%NumTXvars), &
                              j=1,User_InputLink%NumSubDs)
         READ(INunit,*)
         READ(INunit,120) (User_InputLink%NumTerms(i),i=1,User_InputLink%NumSubDs)
         DO i=1,User_InputLink%NumTterms
            READ(INunit,*)             
            CALL TermInfo_Read(User_InputLink%Terms(i),INunit,OUTunit)
         END DO
         READ(INunit,*)
      END IF
      IF(User_InputLink%XSDefined) THEN
         READ(INunit,60)  NumEGrps,NumStates,NumMacs
         READ(INunit,70)  NumMicNs,LengthMics
         CALL InputLink_DefineXS(User_InputLink,NumTXvars,NumStates,NumMacs,NumMicNs,LengthMics)
         User_InputLink%NumEGrps = NumEGrps
         READ(INunit,*)  
         READ(INunit,120) (User_InputLink%StateSelected(i),i=1,User_InputLink%NumStates)
         READ(INunit,*)
         READ(INunit,120) (User_InputLink%TypeNs(i),i=1,User_InputLink%NumMicNs)
         READ(INunit,*)
         READ(INunit,120) (User_InputLink%NumMicParas(i),i=1,User_InputLink%NumMicNs)
         DO i=1,User_InputLink%NumStates
            DO j=1,User_InputLink%NumTXvars
               CALL ParaInfo_Print(User_InputLink%XvarParas(j,i),OUTunit)
            END DO
            DO j=1,User_InputLink%NumMacs
               CALL ParaInfo_Print(User_InputLink%MacParas(j,i),OUTunit)
            END DO
            DO j=1,User_InputLink%LengthMics
               CALL ParaInfo_Print(User_InputLink%MicParas(j,i),OUTunit)
            END DO
         END DO
      END IF
30    FORMAT("[Link_Input] NumSubDs,NumTterms:",I,I)
50    FORMAT("[Link_Input]",100(2X,A))
60    FORMAT("[Link_Input] NumEGrps,NumStates,NumMacs:",I,I,I)
70    FORMAT("[Link_Input] NumMicNs,LengthMics:",I,I)
90    FORMAT("[Link_Input] User_InputLink CR defined?",L)
110   FORMAT("[Link_Input]",100(2X,1PE16.6))
120   FORMAT("[Link_Input]",100(2X,I))
150   FORMAT("[Link_Input] User_InputLink XS defined?",L)
160   FORMAT("[Link_Input] NumTXvars",I)
130   FORMAT("[Link_Input] Non-fatal error !")   
140   FORMAT("[Link_Input] UNDEFINED object User_InputLink in file INunit !")   
170   FORMAT("[Link_Input] Object User_InputLink stays UNDEFINED !")   
   END SUBROUTINE InputLink_Read

   !-----------------------------------------------------------------------------------------------
   ! Void InputLink type variables
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE InputLink_Void(User_InputLink)
      IMPLICIT NONE
      ! Passed out
      TYPE (InputLink_Type) :: User_InputLink
      ! Locals
      INTEGER :: I,J
      IF(User_InputLink%CRDefined) THEN
         DO I=1,User_InputLink%NumTterms
            CALL TermInfo_Void(User_InputLink%Terms(I))
         END DO
         DEALLOCATE (User_InputLink%XvarNames,User_InputLink%SubDStates, &
                     User_InputLink%TypStates,                           &
                     User_InputLink%NumTerms,User_InputLink%Terms)         
         User_InputLink%NumTXvars = 0
         User_InputLink%NumSubDs  = 0
         User_InputLink%NumTterms = 0         
         User_InputLink%CRDefined = .FALSE.
      END IF
      IF(User_InputLink%XSDefined) THEN
         DO I=1,User_InputLink%NumStates
            DO J=1,User_InputLink%NumTXvars
               CALL ParaInfo_Void(User_InputLink%XvarParas(J,I))
            END DO
            DO J=1,User_InputLink%NumMacs
               CALL ParaInfo_Void(User_InputLink%MacParas(J,I))
            END DO
            DO J=1,User_InputLink%LengthMics
               CALL ParaInfo_Void(User_InputLink%MicParas(J,I))
            END DO           
         END DO
         DEALLOCATE (User_InputLink%StateSelected,User_InputLink%XvarParas, &
                     User_InputLink%MacParas,User_InputLink%TypeNs,         &
                     User_InputLink%NumMicParas,User_InputLink%MicParas)
         User_InputLink%NumTXvars  = 0
         User_InputLink%NumStates  = 0
         User_InputLink%NumMacs    = 0
         User_InputLink%NumMicNs   = 0
         User_InputLink%LengthMics = 0
         User_InputLink%XSDefined = .FALSE.
      END IF
   END SUBROUTINE InputLink_Void

END MODULE Input_Link