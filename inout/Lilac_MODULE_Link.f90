!--------------------------------------------------------------------------------------------------
!                           Copyright(c)2013 Xi'an Jiaotong University
!--------------------------------------------------------------------------------------------------
!
!  Developed by Nuclear Engineering Computational Physics (NECP) Laboratory.
!  Only fitting function available for instant.
!
!**************************************************************************************************
!  Version: v1.0                                                                    Date: 8/23/2013
!  Purpose:
!     This is a fitting/interpolating module of Link methods.
!
!  Authors list:
!     Name                       Unit                              Email
!     =========================  ================================  ================================
!     Shengnan Gao               XJTU NECP                         gaoshengnan1989@163.com
!     Yunzhao Li                 XJTU NECP                         Yunzhao@mail.xjtu.edu.cn
!
!  Description:
!     Only least-square fitting for instant.
!     Following parameters values depending on user's choice:
!        Number of X variables, their highest orders & combinations
!        Number of Y variables & Number of their types
!     Functions:
!        Two main functions:
!           Link_Fit: Fitting & Providing coefficients 
!           Link_Interpolate: Getting fitted values at certain states
!        Link_Fit reads information from input files, based on which a linear system is &
!        obtained. The left side of this system is performed by calling function        &
!        YXdiscrete_AmatrixFromX, while the right side by calling function              &
!        YXdiscrete_bVectorFromY. These procedures need several other functions, namly  &
!        YXdiscrete_OrderMatrix, etc.
!        Link_Interpolate reads information from link output files and gives fitted     &
!        results as demanded. This is made possible mainly on two steps: getting the    &
!        terms' values by Link_TermValueFromC and obtaining the final results by        &
!        Link_YfromTerms.
!
!  Subroutine list:
!     SUBROUTINE Link_Functionalization(User_Link,INunit,OUTunit,RCoef)
!     SUBROUTINE Link_Interpolate(User_Link,Fitted,Xvars,NumXvars,Yvars,NumTYvars)
!     SUBROUTINE Link_TermValueFromC(TermValue,C,OrderMatrix,XVars,NumXvars,NumCoefs)
!     SUBROUTINE Link_Examine(User_Link,INunit,INunitG,OUTunit,INunit2,RCoef)
!     SUBROUTINE YXdiscrete_AmatrixFromX(User_YXdiscrete,Amatrix,Asize)
!     SUBROUTINE YXdiscrete_bVectorFromY(User_YXdiscrete,bVector,bSize)
!     SUBROUTINE YXdiscrete_OrderMatrix(User_YXdiscrete)
!     SUBROUTINE YXdiscrete_Define(User_YXdiscrete,NumXvars,NumStates,NumCoefs)
!     SUBROUTINE YXdiscrete_Print(User_YXdiscrete,OUTunit)
!     SUBROUTINE YXdiscrete_Read(User_YXdiscrete,INunit,OUTunit)
!     SUBROUTINE YXdiscrete_GetStorage(User_YXdiscrete,NumINT,NumREAL,NumLOG,NumCHR)
!     SUBROUTINE YXdiscrete_Void(User_YXdiscrete)
!     SUBROUTINE YXdiscrete_SetXOptions(User_YX,XOrder,NumXvars,Option)
!     SUBROUTINE YXdiscrete_SetX(User_YX,X,NumStates,NumXvars)
!     SUBROUTINE YXdiscrete_GetOrderMatrix(User_YX,OrderMatrix,NumCoefs,NumXvars)
!     SUBROUTINE YXdiscrete_SetOrderMatrix(User_YX,OrderMatrix,NumCoefs,NumXvars)
!     SUBROUTINE YXdiscrete_SetY(User_YX,Y,NumStates)
!     SUBROUTINE Link_Define(User_Link,NumXvars,NumStates,NumCoefs)
!     SUBROUTINE Link_Print(User_Link,OUTunit)
!     SUBROUTINE Link_Read(User_Link,INunit,OUTunit)
!     SUBROUTINE Link_GetStorage(User_Link,NumINT,NumREAL,NumLOG,NumCHR)
!     SUBROUTINE Link_Void(User_Link)
!************************************************************************************************** 

MODULE Link
USE LinearSystem
IMPLICIT NONE

!--------------------------------------------------------------------------------------------------
! Storage of values of all points available for fitting
!--------------------------------------------------------------------------------------------------
TYPE,PUBLIC :: YXdiscrete   
   INTEGER             :: NumStates = 0      ! Number of states
   INTEGER             :: NumXvars  = 0      ! Number of X variables
   INTEGER,ALLOCATABLE :: XOptions(:,:)      !(NumXvars,2) XOptions, fit/interpolate options 
   LOGICAL             :: XDefined  = .FALSE.
   REAL(8),ALLOCATABLE :: X(:,:)             !(NumStates,NumXvars) X variables 
   LOGICAL             :: YDefined  = .FALSE.
   REAL(8),ALLOCATABLE :: Y(:)               !(NumStates) Y values obtained on X points
   INTEGER             :: NumCoefs  = 0      ! Number of coefficients
   INTEGER,ALLOCATABLE :: OrderMatrix(:,:)   !(NumXvars,NumCoefs) Orders of basic functions
END TYPE YXdiscrete

!--------------------------------------------------------------------------------------------------
! Storage of all fitting information, including how the function is 
!--------------------------------------------------------------------------------------------------
TYPE,PUBLIC :: Link_Type
   LOGICAL                  :: Defined   = .FALSE.
   INTEGER                  :: NumTYvars = 0    ! Total number of Y variables
   CHARACTER(8),ALLOCATABLE :: NameYvars(:)     !(NumTYvars) Y variable names
   INTEGER                  :: NumTXvars = 0    ! Total number of X variables
   CHARACTER(8),ALLOCATABLE :: NameXvars(:)     !(NumTXvars) X variable names
   REAL(8),ALLOCATABLE      :: Xvars(:)         !(NumTXvars) X variables
   INTEGER                  :: NumTerms  = 0
   CHARACTER(8),ALLOCATABLE :: NameTerm(:)      !(NumTerms) Sub-XS term names
   REAL(8),ALLOCATABLE      :: Terms(:,:)       !(NumTYvars,NumTerms) Term values
   INTEGER,ALLOCATABLE      :: Parentheses(:)   !(NumTerms) Parentheses, whether or not    &
                                                ! a term has parentheses before/after it 
   INTEGER,ALLOCATABLE      :: Operator(:)      !(NumTerms) Term operators, caculation operator  &
                                                ! after each term
   TYPE (YXdiscrete)        :: YX               ! All points used for fitting
   TYPE (LinearSystem_Type) :: LSystem          
END TYPE Link_Type

CONTAINS

   !-----------------------------------------------------------------------------------------------
   ! Calculate fitting coefficients
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE Link_Functionalization(User_Link,INunit,OUTunit,RCoef)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)  :: INunit   ! Channel number for "GenaralInfo" & "INPAXXXXXXXXSXXXXXXXX"
      INTEGER,INTENT(IN)  :: OUTunit  ! Channel number for "OUTAXXXXXXXXSXXXXXXXX"
      INTEGER,INTENT(IN)  :: RCoef    ! Whether the coefficients should be reset
      ! Passed out
      TYPE(Link_Type),INTENT(OUT)    :: User_Link
      ! Locals
      INTEGER             :: NumFuels ! Number of fuels
      INTEGER             :: NumGrps  ! Number of groups
      INTEGER             :: NumMacs  ! Number of macro-scopic cross sections
      INTEGER             :: NumMics  ! Number of nuclids
      INTEGER             :: NumLines ! Used for formatted read/write
      INTEGER             :: NumTMics ! Total number of few-group parameters of nuclids

      ! The following variables are used for reading information from input files
      ! They have exactly the same significations as in types defined in this module and in module
      ! "Input_Link"
      INTEGER             :: NumTerms,NumTXvars,NumTYvars,NumXvars
      INTEGER             :: NPara,ParaType,ParaLength,AuxiDim
      INTEGER             :: NumStates,NumCoefs,Asize,Size,NLIT
      INTEGER,ALLOCATABLE :: Option(:),XOrder(:),IndexXvars(:),TOption(:),OrderMatrix(:,:)
      INTEGER,ALLOCATABLE :: NumMicParas(:),NumSubDs(:)
      REAL(8),ALLOCATABLE :: Amatrix(:,:),bVector(:),Solution(:),SubDStates(:,:),Y(:),X(:,:)

      INTEGER             :: i,j,k,m,n,L,IDFuels,IDSubDs
      REAL(8)             :: BNav,STD,Yav ! Average values & 方差
      REAL(8)             :: RES          ! The residua of G-S iterations 
      REAL(8)             :: EPS = 1.0E-8 ! Precision of G-S iterations 
      CHARACTER(36)       :: InputFile   = "LinkInputFiles/INPA00000000S00000000"
      CHARACTER(37)       :: OutputFile  = "LinkOutputFiles/OUTA00000000S00000000"
      CHARACTER(37)       :: OUTFileName = "LinkOutputFiles/RSTA00000000S00000000" ! 变量代换后的系数           
      REAL(8),PARAMETER   :: NECP_RealZero  = 1.0E-18
      REAL(8),PARAMETER   :: NECP_RealEqual = 1.0E-8

      ! Read in control parameters through INunit
      OPEN (INunit,file="LinkInputFiles/GenaralInfo")
      READ (INunit,*) NumFuels
      READ (INunit,*) 
      ALLOCATE (NumSubDs(NumFuels))
      DO i=1,NumFuels
         READ (INunit,*) NumSubDs(i)
         READ (INunit,*) NumGrps,NumMacs,NumMics
         IF(MOD(NumMics,5).EQ.0) THEN
            NumLines = NumMics / 5
         ELSE
            NumLines = NumMics / 5 + 1
         END IF
         DO j=1,NumLines
            READ (INunit,*) 
         END DO
         DO j=1,NumLines
            READ (INunit,*) 
         END DO
         ALLOCATE (NumMicParas(NumMics))
         READ (INunit,*) (NumMicParas(j),j=1,NumMics)
         NumTMics = 0
         DO j=1,NumMics
            NumTMics = NumTMics + NumMicParas(j)
         END DO
         DEALLOCATE (NumMicParas)
         READ (INunit,*) 
         DO j=1,NumMacs
            READ (INunit,*) NPara,ParaType,ParaLength
            IF(ParaType.EQ.108) THEN
               READ (INunit,*) AuxiDim
               DO m=1,AuxiDim
                  READ (INunit,*)
               END DO
            END IF
         END DO
         READ (INunit,*) 
         DO j=1,NumTMics
            READ (INunit,*) 
         END DO
         READ (INunit,*) 
         READ (INunit,*) NumTXvars
         READ (INunit,*) 
         DO j=1,NumSubDs(i)
            DO m=1,NumTXvars
               READ (INunit,*) 
            END DO
            READ (INunit,*) 
         END DO
      END DO
      CLOSE (INunit)

      ! Fit and write results into corresponding files
      OPEN (90,file='bVector')                                                                  
      DO IDFuels = 1,NumFuels
         DO IDSubDs = 1,NumSubDs(IDFuels)
            ! Information for each subdivision of each fuel is stocked in a separated file
            WRITE(InputFile(20:27),"(I8.8)")   IDFuels                                                  
            WRITE(OutputFile(21:28),"(I8.8)")  IDFuels                                                  
            WRITE(OUTFileName(21:28),"(I8.8)") IDFuels                                                  
            WRITE(InputFile(29:36),"(I8.8)")   IDSubDs                                                  
            WRITE(OutputFile(30:37),"(I8.8)")  IDSubDs                                                  
            WRITE(OUTFileName(30:37),"(I8.8)") IDSubDs                                                  
            ! Read in LinkInput file                                            
            OPEN (INunit,file=InputFile)                                                                
            OPEN (OUTunit,file=OutputFile)                                                               
            ! Read the header of LinkInput                                                              
            READ (INunit,*)                                                                             
            WRITE(OUTunit,*) '*************************************************************************'
            READ (INunit,*)                                                                             
            WRITE(OUTunit,*) 'NECP: Output File of Module Link'                                         
            READ (INunit,*)                                                                             
            WRITE(OUTunit,*) '*************************************************************************'
            READ (INunit,*)                                                                             
            READ (INunit,*)                                                                             
            ! Read in term information                                                                  
            READ (INunit,*) NumTerms,NumTXvars,NumTYvars                                                
            WRITE(OUTunit,'(3(2X,I6))') NumTerms,NumTXvars,NumTYvars                                       
            CALL Link_Define(User_Link,NumTYvars,NumTXvars,NumTerms)                                    
            READ (INunit,*) (User_Link%NameXvars(i),i=1,User_Link%NumTXvars)                            
            WRITE(OUTunit,'(5(2X,A8))') (User_Link%NameXvars(i),i=1,User_Link%NumTXvars)                   
            READ (INunit,*) (User_Link%NameYvars(i),i=1,User_Link%NumTYvars)                            
            WRITE(OUTunit,'(5(2X,A8))') (User_Link%NameYvars(i),i=1,User_Link%NumTYvars)                   
            ALLOCATE (TOption(NumTerms))                                                                
            DO i=1,NumTerms                                                                             
               READ (INunit,*)  User_Link%NameTerm(i),User_Link%Parentheses(i),User_Link%Operator(i)    
               READ (INunit,*)  TOption(i)                                                              
               WRITE(OUTunit,*) User_Link%NameTerm(i),User_Link%Parentheses(i),User_Link%Operator(i), & 
                                TOption(i)                                                              
            END DO                                                                                      
            DEALLOCATE (TOption)                                                                        
            ALLOCATE (SubDStates(2,User_Link%NumTXvars))                                                
            READ (INunit,*) ((SubDStates(i,j),i=1,2),j=1,User_Link%NumTXvars)                           
            DEALLOCATE (SubDStates)                                                                     
            ! Read in X and Y variables used for fitting                                                                        
            READ (INunit,*)                                                                             
            WRITE(OUTunit,*) '*************************************************************************'
            READ (INunit,*)                                                                             
            WRITE(OUTunit,*) 'X Y Variables'                                                            
            READ (INunit,*)                                                                             
            WRITE(OUTunit,*) '*************************************************************************'
            ! Block for getting information in each term
            DO i=1,NumTerms                                                                             
               ! General information comes first
               READ (INunit,*)  NumXvars,NumStates                                                      
               WRITE(OUTunit,*) NumXvars                                                                
               READ (INunit,*)                                                                          
               ALLOCATE (Option(NumXvars),XOrder(NumXvars),IndexXvars(NumXvars))                        
               ! Then information of fitting orders and X varibales
               READ (INunit,*)  (Option(j),j=1,NumXvars)                                                
               WRITE(OUTunit,*) (Option(j),j=1,NumXvars)                                                
               DO j=1,NumXvars                                                                          
                  IF (Option(j)/=1) THEN                                                                
                     WRITE (*,*) 'ERROR: Sorry, Least Square only, Option must be 1'                    
                  END IF                                                                                
               END DO                                                                                   
               READ (INunit,*)  (XOrder(j),j=1,NumXvars)                                                
               WRITE(OUTunit,*) (XOrder(j),j=1,NumXvars)                                                
               READ (INunit,*)  (IndexXvars(j),j=1,NumXvars)                                            
               WRITE(OUTunit,*) (IndexXvars(j),j=1,NumXvars)                                            
               WRITE(OUTunit,*) '**********************'                                                
               ! Determine the number of coefficients according to fitting orders
               NumCoefs = 1                                                                             
               DO j=1,NumXvars                                                                          
                  NumCoefs = NumCoefs*(XOrder(j)+1)                                                     
               END DO                                                                                   
               ! Set corresponding array values in User_Link
               CALL YXdiscrete_Define(User_Link%YX,NumXvars,NumStates,NumCoefs)                         
               CALL YXdiscrete_SetXOptions(User_Link%YX,XOrder,NumXvars,Option)                                                                           
               DEALLOCATE (Option,XOrder,IndexXvars)                                                    
               READ (INunit,*)                                                                          
               ! Get fitting points X&Y values
               ALLOCATE (X(NumStates,NumXvars))
               DO j=1,NumXvars                                                                          
                  READ (INunit,*) (X(k,j),k=1,NumStates)                                   
                  !READ (INunit,*) (User_Link%YX%X(k,j),k=1,NumStates)                                   
                  READ (INunit,*)                                                                       
               END DO 
               CALL YXdiscrete_SetX(User_Link%YX,X,NumStates,NumXvars)
               DEALLOCATE (X)
               DO j=1,NumXvars                                                                          
                  READ (INunit,*) BNav,STD                                                              
                  WRITE(OUTunit,'(2(2X,ES15.8))') BNav,STD                                              
               END DO                                                                                   
               WRITE(OUTunit,*) '************************'                                              
               ! Solve least square problem                                                             
               CALL YXdiscrete_OrderMatrix(User_Link%YX)                                                
               !DO j=1,NumXvars                                                                         
               !   WRITE(*,*) (User_Link%YX%OrderMatrix(j,k),k=1,NumCoefs)                              
               !END DO                                                                                  
               ! Get the matrix A on the left-hand side of the equation
               ALLOCATE (Amatrix(NumCoefs,NumCoefs))                                                    
               CALL YXdiscrete_AmatrixFromX(User_Link%YX,Amatrix,NumCoefs)                              
               DO m=1,NumCoefs                                                                         
                  DO n=1,NumCoefs-1                                                                    
                     WRITE(90,'(2X,ES12.5)',ADVANCE='NO') Amatrix(m,n)                                 
                  END DO                                                                               
                  WRITE(90,'(2X,ES12.5)') Amatrix(m,n)                                                 
               END DO                                                                                  
               CALL LinearSystem_Define(User_Link%LSystem,NumCoefs)                                     
               CALL LinearSystem_SetA(User_Link%LSystem,Amatrix,NumCoefs)                               
               DEALLOCATE (Amatrix)                                                                     
               !CALL YXdiscrete_XVoid(User_Link%YX)                                                     
               ! Order matrix stocks all powers of X variables and their combinations
               ALLOCATE (OrderMatrix(NumXvars,NumCoefs))
               WRITE(OUTunit,*) 'ORDER MATRIX:'
               CALL YXdiscrete_GetOrderMatrix(User_Link%YX,OrderMatrix,NumCoefs,NumXvars)                                                         
               WRITE(OUTunit,*) ((OrderMatrix(m,n),n=1,NumCoefs),m=1,NumXvars)             
               DEALLOCATE (OrderMatrix)
               WRITE(OUTunit,*) '**********************'                                                
               READ (INunit,*)                                                                          
               READ (INunit,*) 
               ALLOCATE (Y(NumStates))
               DO j=1,NumTYvars                                                                         
                  WRITE(OUTunit,*) j                                                                    
                  ! Get the b vector on the right-hand side of the equation
                  READ (INunit,*) (Y(k),k=1,NumStates)
                  !READ (INunit,*) (User_Link%YX%Y(k),k=1,NumStates)                                     
                  ! When the average value of Y is not zero, all Y are divised by Y average
                  Yav = 0.0D0
                  DO k=1,NumStates
                     Yav = Yav + Y(k)
                  END DO
                  Yav = Yav / NumStates
                  IF(ABS(Yav-0.0D0).GT.NECP_RealZero) THEN
                     DO k=1,NumStates
                        Y(k) = Y(k) / Yav
                     END DO
                  ELSE
                     Yav = 1.0D0
                  END IF
                  CALL YXdiscrete_SetY(User_Link%YX,Y,NumStates)
                  ALLOCATE (bVector(NumCoefs))                                                          
                  CALL YXdiscrete_bVectorFromY(User_Link%YX,bVector,NumCoefs)                           
                  WRITE(90,'(I4)') j                                                                   
                  DO k=1,NumCoefs-1                                                                    
                     WRITE(90,'(2X,ES15.8)') bVector(k)                                                
                  END DO                                                                               
                  WRITE(90,'(2X,ES15.8)') bVector(k)                                                   
                  ! Solve the equation
                  CALL LinearSystem_Setb(User_Link%LSystem,bVector,NumCoefs)                            
                  !CALL LinearSystem_SolverGS(User_Link%LSystem,RES,EPS,NLIT)
                  CALL AGSDL(User_Link%LSystem%Amatrix,User_Link%LSystem%bVector,NumCoefs, &            
                             User_Link%LSystem%xVector,EPS,L)  
                  WRITE(90,*) "RESULT",j,EPS,L
                  DEALLOCATE (bVector)                                                                  
                  ALLOCATE (Solution(NumCoefs))                                                         
                  CALL LinearSystem_GetSolution(User_Link%LSystem,Solution,NumCoefs)                    
                  DO k=1,NumCoefs
                     Solution(k) = Solution(k) * Yav
                  END DO
                  ! Write into LinkOutput file                                        
                  WRITE (OUTunit,200) (Solution(k),k=1,NumCoefs)                                        
                  WRITE (OUTunit,*) '*************'                                                     
                  READ (INunit,*)                                                                       
                  DEALLOCATE (Solution)                                                         
               END DO                                                                                   
               DEALLOCATE (Y)
               CALL YXdiscrete_Void(User_Link%YX)                                                       
            END DO                                                                                      
        200 FORMAT(5(2X,ES15.8))                                                                        
            CALL  LinearSystem_Void(User_Link%LSystem)                                                  
            CALL  Link_Void(User_Link)                                                                  
            CLOSE (OUTunit)                                                                             
            CLOSE (INunit)                                                                              
            ! Reset coefficients if needed
            IF(RCoef.EQ.1) THEN                                                                         
               CALL Link_ResetCoefs(INunit,OUTunit,OutputFile,OUTFileName)                               
            END IF                                                                                      
         END DO
      END DO
      CLOSE (90)                                                                                 
      DEALLOCATE (NumSubDs)
   END SUBROUTINE Link_Functionalization

   !-----------------------------------------------------------------------------------------------
   ! Set X options of User_YX
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_SetXOptions(User_YX,XOrder,NumXvars,Option)
      ! Passed in
      INTEGER,INTENT(IN) :: NumXvars
      INTEGER,INTENT(IN) :: XOrder(NumXvars),Option(NumXvars)
      ! Passed out
      TYPE(YXdiscrete)   :: User_YX
      ! Locals
      INTEGER            :: j     
      IF(.NOT.User_YX%XDefined) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set X options, YX not defined'
         RETURN
      END IF
      IF(User_YX%NumXvars.NE.NumXvars) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set X options, given NumXvars incorrect'
         RETURN
      END IF
      DO j=1,NumXvars
         User_YX%XOptions(j,1) = Option(j)
         User_YX%XOptions(j,2) = XOrder(j)
      END DO
   END SUBROUTINE YXdiscrete_SetXOptions

   !-----------------------------------------------------------------------------------------------
   ! Read X values through INunit
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_SetX(User_YX,X,NumStates,NumXvars)
      ! Passed in
      INTEGER,INTENT(IN) :: NumStates,NumXvars
      REAL(8),INTENT(IN) :: X(NumStates,NumXvars)
      ! Passed out
      TYPE(YXdiscrete)   :: User_YX 
      ! Locals
      INTEGER            :: i,j
      IF(.NOT.User_YX%XDefined) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set X, YX not defined'
         RETURN
      END IF
      IF(User_YX%NumXvars.NE.NumXvars) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set X, given NumXvars incorrect'
         RETURN
      END IF
      IF(User_YX%NumStates.NE.NumStates) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set X, given NumStates incorrect'
         RETURN
      END IF
      DO i=1,NumStates
         DO j=1,NumXvars
            User_YX%X(i,j) = X(i,j)
         END DO
      END DO
   END SUBROUTINE YXdiscrete_SetX

   !-----------------------------------------------------------------------------------------------
   ! Get the OrderMatrix of User_YX
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_GetOrderMatrix(User_YX,OrderMatrix,NumCoefs,NumXvars)
      ! Passed in
      TYPE(YXdiscrete),INTENT(IN) :: User_YX
      INTEGER,INTENT(IN)          :: NumCoefs,NumXvars
      ! Passed out
      INTEGER,INTENT(OUT)         :: OrderMatrix(NumXvars,NumCoefs)
      ! Locals
      INTEGER                     :: m,n
      IF(.NOT.User_YX%XDefined) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot get order matrix, YX not defined'
         RETURN
      END IF      
      IF(User_YX%NumXvars.NE.NumXvars) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot get order matrix, given NumXvars incorrect'
         RETURN
      END IF
      IF(User_YX%NumCoefs.NE.NumCoefs) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot get order matrix, given NumCoefs incorrect'
         RETURN
      END IF
      DO m=1,NumXvars
         DO n=1,NumCoefs
            OrderMatrix(m,n) = User_YX%OrderMatrix(m,n)
         END DO
      END DO
   END SUBROUTINE YXdiscrete_GetOrderMatrix

   !-----------------------------------------------------------------------------------------------
   ! Set the OrderMatrix of User_YX
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_SetOrderMatrix(User_YX,OrderMatrix,NumXinTerm,NumCoefs)
      ! Passed in
      INTEGER,INTENT(IN) :: NumXinTerm,NumCoefs
      INTEGER,INTENT(IN) :: OrderMatrix(NumXinTerm,NumCoefs)
      ! Passed out
      TYPE(YXdiscrete)   :: User_YX
      ! Locals      
      INTEGER            :: m,n
      IF(.NOT.User_YX%XDefined) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set order matrix, YX not defined'
         RETURN
      END IF      
      IF(User_YX%NumXvars.NE.NumXinTerm) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set order matrix, given NumXvars incorrect'
         RETURN
      END IF
      IF(User_YX%NumCoefs.NE.NumCoefs) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set order matrix, given NumCoefs incorrect'
         RETURN
      END IF
      DO m=1,NumXinTerm
         DO n=1,NumCoefs
            User_YX%OrderMatrix(m,n) = OrderMatrix(m,n) 
         END DO
      END DO
   END SUBROUTINE YXdiscrete_SetOrderMatrix

   !-----------------------------------------------------------------------------------------------
   ! Set Y of User_YX 
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_SetY(User_YX,Y,NumStates)
      ! Passed in
      INTEGER,INTENT(IN) :: NumStates
      REAL(8),INTENT(IN) :: Y(NumStates)
      ! Passed out
      TYPE(YXdiscrete)   :: User_YX
      ! Locals
      INTEGER            :: k
      IF(.NOT.User_YX%YDefined) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set Y, YX not defined'
         RETURN
      END IF      
      IF(User_YX%NumStates.NE.NumStates) THEN
         WRITE (*,*) '[Link] FATAL ERROR: Cannot set Y, given NumStates incorrect'
         RETURN
      END IF
      DO k=1,NumStates
         User_YX%Y(k) = Y(k) 
      END DO
   END SUBROUTINE YXdiscrete_SetY

   !-----------------------------------------------------------------------------------------------
   ! Re-organise output information to match SIMME input file format
   ! Format change only
   ! Control parameters which are at the top of link output files are replaced in left colonnes of 
   ! the file "OUTOrange" to fit in SIMME
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE Link_ORange(INunitG,OUTunit)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN)       :: INunitG,OUTunit ! Channel numbers
      ! Locals
      CHARACTER(37)            :: INFileName  = "LinkOutputFiles/RSTA00000000S00000000"
      CHARACTER(9)             :: OUTFileName = "OUTORange"
      INTEGER                  :: NumFuels,NumSubDs,NumEGrps,NumMacs,NumNuclds
      INTEGER                  :: NumTMics,AuxiDim,NumXInTerm,NumMacsBS,NumLines
      INTEGER                  :: NumTerms,NumTXvars,NumTYRead,NumCoefs,INTTERM
      INTEGER                  :: i,j,k,m,n,kk,IDMacTag,IDX,IDY,IDk
      LOGICAL                  :: SigSExist = .FALSE.
      INTEGER,ALLOCATABLE      :: TypeNs(:),NMics(:),NumMicParas(:),IGstart(:),IGend(:)
      INTEGER,ALLOCATABLE      :: MacParaType(:),MacLength(:),SelfPosition(:)
      INTEGER,ALLOCATABLE      :: MicParaType(:),MicLength(:),TOption(:),OrderMatrix(:,:)
      INTEGER,ALLOCATABLE      :: XOrders(:)
      INTEGER                  :: BU = 1
      INTEGER                  :: CB = 7
      INTEGER                  :: Pr = 3
      INTEGER                  :: Ur = 4
      INTEGER                  :: zero = 0
      CHARACTER(4)             :: XSNames(15)
      TYPE (Link_Type)         :: User_Link
      TYPE (LinearSystem_Type),ALLOCATABLE :: MacCoefs(:,:)

      ! Define cross section names 
      XSNames(1)  = "1.00"
      XSNames(2)  = "2.00"
      XSNames(3)  = "3.00"
      XSNames(4)  = "4.00"
      XSNames(5)  = "5.00"
      XSNames(6)  = "6.00"
      XSNames(7)  = "7.00"
      XSNames(8)  = "1.11"
      XSNames(9)  = "1.12"
      XSNames(10) = "1.13"
      XSNames(11) = "1.14"
      XSNames(12) = "2.21"
      XSNames(13) = "2.22"
      XSNames(14) = "2.23"
      XSNames(15) = "2.24"
      ! Read GenaralInfo file
      OPEN (INunitG,file="LinkInputFiles\GenaralInfo")
      READ (INunitG,*) NumFuels
      READ (INunitG,*) 
      READ (INunitG,*) NumSubDs
      READ (INunitG,*) NumEGrps,NumMacs,NumNuclds
      ALLOCATE (TypeNs(NumNuclds),NMics(NumNuclds))
      READ (INunitG,*) (NMics(j),j=1,NumNuclds)
      READ (INunitG,*) (TypeNs(j),j=1,NumNuclds)
      ALLOCATE (NumMicParas(NumNuclds))
      READ (INunitG,*) (NumMicParas(j),j=1,NumNuclds)
      READ (INunitG,*) 
      ALLOCATE (MacParaType(NumMacs),MacLength(NumMacs))
      DO j=1,NumMacs
         READ (INunitG,"(8X,I4,2X,I4)") MacParaType(j),MacLength(j) 
         IF(MacParaType(j).EQ.108) THEN
            SigSExist = .TRUE.
            READ (INunitG,*) AuxiDim
            ALLOCATE (IGstart(AuxiDim),IGend(AuxiDim),SelfPosition(AuxiDim))
            DO m=1,AuxiDim
               READ (INunitG,*) IGstart(m),IGend(m),SelfPosition(m)
            END DO
         END IF
      END DO
      READ (INunitG,*) 
      NumTMics = 0
      DO j=1,NumNuclds
         NumTMics = NumTMics + NumMicParas(j)
      END DO
      ALLOCATE (MicParaType(NumTMics),MicLength(NumTMics))
      DO j=1,NumTMics
         READ (INunitG,"(8X,I4,2X,I4)") MicParaType(j),MicLength(j)
      END DO
      CLOSE (INunitG)
      NumMacsBS = 0
      DO j=1,NumMacs
         IF(MacParaType(j).EQ.108) EXIT
         NumMacsBS = NumMacsBS + 1
      END DO
      OPEN (OUTunit,file=OUTFileName)       
      ! Get coefficients from files and range them in the order demanded by program MG-Cate 
      DO i=1,NumSubDs
         DO j=1,NumFuels
            ! Read coefficients
            WRITE (INFileName(21:28),"(I8.8)") j
            WRITE (INFileName(30:37),"(I8.8)") i
            OPEN (INunitG,file=INFileName)       
            READ (INunitG,*)
            READ (INunitG,*)
            READ (INunitG,*)
            READ (INunitG,*)  NumTerms,NumTXvars,NumTYRead
            CALL Link_Define(User_Link,NumTYRead,NumTXvars,NumTerms)
            READ (INunitG,*)  (User_Link%NameXvars(m),m=1,User_Link%NumTXvars)
            READ (INunitG,*)  (User_Link%NameYvars(m),m=1,User_Link%NumTYvars)
            ALLOCATE (TOption(NumTerms))
            DO m=1,NumTerms
               READ (INunitG,*) User_Link%NameTerm(m),User_Link%Parentheses(m), &
                                User_Link%Operator(m),TOption(m)
            END DO
            DEALLOCATE (TOption)
            CALL Link_Void(User_Link)
            READ (INunitG,*)
            READ (INunitG,*)
            READ (INunitG,*)
            ALLOCATE (MacCoefs(NumTYRead,NumTerms))
            ALLOCATE (XOrders(7))
            IDX = 1
            DO m=1,NumTerms
               READ (INunitG,*) NumXInTerm
               READ (INunitG,*)
               READ (INunitG,*) (XOrders(n),n=IDX,IDX+NumXInTerm-1)
               NumCoefs = 1
               DO n=IDX,IDX+NumXInTerm-1
                  NumCoefs = NumCoefs * (XOrders(n)+1)
               END DO
               IDX = IDX + NumXInTerm
               READ (INunitG,*)
               READ (INunitG,*)
               DO n=1,NumXInTerm
                  READ (INunitG,*)
               END DO
               ALLOCATE (OrderMatrix(NumXInTerm,NumCoefs))
               READ (INunitG,*) 
               READ (INunitG,*) 
               READ (INunitG,*) ((OrderMatrix(k,n),n=1,NumCoefs),k=1,NumXInTerm)
               READ (INunitG,*) 
               DEALLOCATE (OrderMatrix)
               DO n=1,NumTYRead
                  CALL LinearSystem_Define(MacCoefs(n,m),NumCoefs)
                  READ (INunitG,*)
                  READ (INunitG,*) (MacCoefs(n,m)%xVector(k),k=1,NumCoefs) 
                  READ (INunitG,*)
               END DO
               !write (*,*) (MacCoefs(100,2)%xVector(k),k=1,NumCoefs) 
            END DO
            ! Write coefficients in the order demanded
            DO n=1,NumEGrps
               ! All coefficients of macro-scopic XS before scatter vector
               WRITE (OUTunit,300,ADVANCE="NO") zero,n,j,i
               WRITE (OUTunit,"(4X,I4,I4)") SelfPosition(n),IGend(n)-IGstart(n)+1
               DO m=1,NumMacsBS
                  IDY = (m-1)*NumEGrps + n
                  IF(MOD(MacCoefs(IDY,1)%Size,4).NE.0) THEN
                     NumLines = MacCoefs(IDY,1)%Size/4 + 1
                  ELSE
                     NumLines = MacCoefs(IDY,1)%Size/4
                  END IF
                  INTTERM = 1
                  IDk = 1
                  WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                  WRITE (OUTunit,"(6X,A4,1X,6(I3))") XSNames(m),BU,CB,Pr,zero,   &
                                                     XOrders(2),XOrders(1)
                  DO k=1,NumLines
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     IF(k.NE.NumLines) THEN
                        WRITE (OUTunit,400) (MacCoefs(IDY,1)%xVector(kk),kk=IDk,IDk+3)
                     ELSE
                        WRITE (OUTunit,400) (MacCoefs(IDY,1)%xVector(kk), &
                                             kk=IDk,MacCoefs(IDY,1)%Size)
                     END IF
                     IDk = IDk + 4 
                  END DO
                  WRITE (OUTunit,"(5(I3))",ADVANCE="NO") INTTERM,n,j,i,XOrders(3)
                  WRITE (OUTunit,500) (MacCoefs(IDY,2)%xVector(k),k=1,MacCoefs(IDY,2)%Size)
                  IF(MOD(MacCoefs(IDY,3)%Size,4).NE.0) THEN
                     NumLines = MacCoefs(IDY,3)%Size/4 + 1
                  ELSE
                     NumLines = MacCoefs(IDY,3)%Size/4
                  END IF
                  INTTERM = 2
                  IDk = 1
                  WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                  WRITE (OUTunit,"(6X,A4,1X,6(I3))") XSNames(m),Ur,CB,zero,   &
                                                     XOrders(5),XOrders(4)
                  DO k=1,NumLines
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     IF(k.NE.NumLines) THEN
                        WRITE (OUTunit,400) (MacCoefs(IDY,3)%xVector(kk),kk=IDk,IDk+3)
                     ELSE
                        WRITE (OUTunit,400) (MacCoefs(IDY,3)%xVector(kk), &
                                             kk=IDk,MacCoefs(IDY,3)%Size)
                     END IF
                     IDk = IDk + 4 
                  END DO
                  IF(NumTerms.EQ.4) THEN
                     IF(MOD(MacCoefs(IDY,4)%Size,4).NE.0) THEN
                        NumLines = MacCoefs(IDY,4)%Size/4 + 1
                     ELSE
                        NumLines = MacCoefs(IDY,4)%Size/4
                     END IF
                     INTTERM = 3
                     IDk = 1
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     WRITE (OUTunit,"(5X,A1,A4,1X,6(I3))") "-",XSNames(m),Ur,CB,zero, &
                                                           XOrders(7),XOrders(6)
                     DO k=1,NumLines
                        WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                        IF(k.NE.NumLines) THEN
                           WRITE (OUTunit,400) (0.0D0-MacCoefs(IDY,4)%xVector(kk),kk=IDk,IDk+3)
                        ELSE
                           WRITE (OUTunit,400) (0.0D0-MacCoefs(IDY,4)%xVector(kk), &
                                                kk=IDk,MacCoefs(IDY,4)%Size)
                        END IF
                        IDk = IDk + 4 
                     END DO
                  END IF
               END DO
               ! Scatter vector
               DO m=1,IGend(n)-IGstart(n)+1
                  IDY = NumMacsBS*NumEGrps
                  DO k=1,n-1
                     IDY = IDY + (IGend(k)-IGstart(k)+1)
                  END DO
                  IDY = IDY + m
                  IF(MOD(MacCoefs(IDY,1)%Size,4).NE.0) THEN
                     NumLines = MacCoefs(IDY,1)%Size/4 + 1
                  ELSE
                     NumLines = MacCoefs(IDY,1)%Size/4
                  END IF
                  INTTERM = 1
                  IDk = 1
                  WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                  WRITE (OUTunit,"(6X,A4,1X,6(I3))") XSNames(7),BU,CB,Pr,zero,   &
                                                     XOrders(2),XOrders(1)
                  DO k=1,NumLines
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     IF(k.NE.NumLines) THEN
                        WRITE (OUTunit,400) (MacCoefs(IDY,1)%xVector(kk),kk=IDk,IDk+3)
                     ELSE
                        WRITE (OUTunit,400) (MacCoefs(IDY,1)%xVector(kk), &
                                             kk=IDk,MacCoefs(IDY,1)%Size)
                     END IF
                     IDk = IDk + 4 
                  END DO
                  WRITE (OUTunit,"(5(I3))",ADVANCE="NO") INTTERM,n,j,i,XOrders(3)
                  WRITE (OUTunit,500) (MacCoefs(IDY,2)%xVector(k),k=1,MacCoefs(IDY,2)%Size)
                  IF(MOD(MacCoefs(IDY,3)%Size,4).NE.0) THEN
                     NumLines = MacCoefs(IDY,3)%Size/4 + 1
                  ELSE
                     NumLines = MacCoefs(IDY,3)%Size/4
                  END IF
                  INTTERM = 2
                  IDk = 1
                  WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                  WRITE (OUTunit,"(6X,A4,1X,6(I3))") XSNames(7),Ur,CB,zero,   &
                                                     XOrders(5),XOrders(4)
                  DO k=1,NumLines
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     IF(k.NE.NumLines) THEN
                        WRITE (OUTunit,400) (MacCoefs(IDY,3)%xVector(kk),kk=IDk,IDk+3)
                     ELSE
                        WRITE (OUTunit,400) (MacCoefs(IDY,3)%xVector(kk), &
                                             kk=IDk,MacCoefs(IDY,3)%Size)
                     END IF
                     IDk = IDk + 4 
                  END DO
                  IF(NumTerms.EQ.4) THEN
                     IF(MOD(MacCoefs(IDY,4)%Size,4).NE.0) THEN
                        NumLines = MacCoefs(IDY,4)%Size/4 + 1
                     ELSE
                        NumLines = MacCoefs(IDY,4)%Size/4
                     END IF
                     INTTERM = 3
                     IDk = 1
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     WRITE (OUTunit,"(5X,A1,A4,1X,6(I3))") "-",XSNames(7),Ur,CB,zero, &
                                                           XOrders(7),XOrders(6)
                     DO k=1,NumLines
                        WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                        IF(k.NE.NumLines) THEN
                           WRITE (OUTunit,400) (0.0D0-MacCoefs(IDY,4)%xVector(kk),kk=IDk,IDk+3)
                        ELSE
                           WRITE (OUTunit,400) (0.0D0-MacCoefs(IDY,4)%xVector(kk), &
                                                kk=IDk,MacCoefs(IDY,4)%Size)
                        END IF
                        IDk = IDk + 4 
                     END DO
                  END IF
               END DO
               IF((NumMacs-NumMacsBS).GT.1) THEN
                  DO m=1,NumMacs-NumMacsBS-1
                     IDY = NumMacsBS*NumEGrps
                     DO k=1,NumEGrps
                        IDY = IDY + (IGend(k)-IGstart(k)+1)
                     END DO
                     IDY = IDY + (m-1)*NumEGrps + n
                     IF(MOD(MacCoefs(IDY,1)%Size,4).NE.0) THEN
                        NumLines = MacCoefs(IDY,1)%Size/4 + 1
                     ELSE
                        NumLines = MacCoefs(IDY,1)%Size/4
                     END IF
                     INTTERM = 1
                     IDk = 1
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     WRITE (OUTunit,"(6X,A4,1X,6(I3))") XSNames(7+m),BU,CB,Pr,zero,   &
                                                        XOrders(2),XOrders(1)
                     DO k=1,NumLines
                        WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                        IF(k.NE.NumLines) THEN
                           WRITE (OUTunit,400) (MacCoefs(IDY,1)%xVector(kk),kk=IDk,IDk+3)
                        ELSE
                           WRITE (OUTunit,400) (MacCoefs(IDY,1)%xVector(kk), &
                                                kk=IDk,MacCoefs(IDY,1)%Size)
                        END IF
                        IDk = IDk + 4 
                     END DO
                     WRITE (OUTunit,"(5(I3))",ADVANCE="NO") INTTERM,n,j,i,XOrders(3)
                     WRITE (OUTunit,500) (MacCoefs(IDY,2)%xVector(k),k=1,MacCoefs(IDY,2)%Size)
                     IF(MOD(MacCoefs(IDY,3)%Size,4).NE.0) THEN
                        NumLines = MacCoefs(IDY,3)%Size/4 + 1
                     ELSE
                        NumLines = MacCoefs(IDY,3)%Size/4
                     END IF
                     INTTERM = 2
                     IDk = 1
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     WRITE (OUTunit,"(6X,A4,1X,6(I3))") XSNames(7+m),Ur,CB,zero,   &
                                                        XOrders(5),XOrders(4)
                     DO k=1,NumLines
                        WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                        IF(k.NE.NumLines) THEN
                           WRITE (OUTunit,400) (MacCoefs(IDY,3)%xVector(kk),kk=IDk,IDk+3)
                        ELSE
                           WRITE (OUTunit,400) (MacCoefs(IDY,3)%xVector(kk), &
                                                kk=IDk,MacCoefs(IDY,3)%Size)
                        END IF
                        IDk = IDk + 4 
                     END DO
                     IF(NumTerms.EQ.4) THEN
                        IF(MOD(MacCoefs(IDY,4)%Size,4).NE.0) THEN
                           NumLines = MacCoefs(IDY,4)%Size/4 + 1
                        ELSE
                           NumLines = MacCoefs(IDY,4)%Size/4
                        END IF
                        INTTERM = 3
                        IDk = 1
                        WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                        WRITE (OUTunit,"(5X,A1,A4,1X,6(I3))") "-",XSNames(7+m),Ur,CB,zero, &
                                                              XOrders(7),XOrders(6)
                        DO k=1,NumLines
                           WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                           IF(k.NE.NumLines) THEN
                              WRITE (OUTunit,400) (0.0D0-MacCoefs(IDY,4)%xVector(kk),kk=IDk,IDk+3)
                           ELSE
                              WRITE (OUTunit,400) (0.0D0-MacCoefs(IDY,4)%xVector(kk), &
                                                   kk=IDk,MacCoefs(IDY,4)%Size)
                           END IF
                           IDk = IDk + 4 
                        END DO
                     END IF
                  END DO
               END IF
               DO m=1,NumNuclds
                  IDY = NumMacsBS*NumEGrps
                  DO k=1,NumEGrps
                     IDY = IDY + (IGend(k)-IGstart(k)+1)
                  END DO
                  IF((NumMacs-NumMacsBS).GT.1) IDY = IDY + (NumMacs-NumMacsBS-1)*NumEGrps
                  IDk = 1
                  DO k=1,m-1
                     DO kk=IDk,IDk+NumMicParas(k)-2
                        IDY = IDY + MicLength(kk)
                     END DO
                     IDk = IDk + NumMicParas(k)
                  END DO
                  IDY = IDY + n
                  IF(MOD(MacCoefs(IDY,1)%Size,4).NE.0) THEN
                     NumLines = MacCoefs(IDY,1)%Size/4 + 1
                  ELSE
                     NumLines = MacCoefs(IDY,1)%Size/4
                  END IF
                  INTTERM = 1
                  IDk = 1
                  WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                  WRITE (OUTunit,"(6X,A4,1X,6(I3))") XSNames(11+m),BU,CB,Pr,zero,   &
                                                     XOrders(2),XOrders(1)
                  DO k=1,NumLines
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     IF(k.NE.NumLines) THEN
                        WRITE (OUTunit,400) (MacCoefs(IDY,1)%xVector(kk),kk=IDk,IDk+3)
                     ELSE
                        WRITE (OUTunit,400) (MacCoefs(IDY,1)%xVector(kk), &
                                             kk=IDk,MacCoefs(IDY,1)%Size)
                     END IF
                     IDk = IDk + 4 
                  END DO
                  WRITE (OUTunit,"(5(I3))",ADVANCE="NO") INTTERM,n,j,i,XOrders(3)
                  WRITE (OUTunit,500) (MacCoefs(IDY,2)%xVector(k),k=1,MacCoefs(IDY,2)%Size)
                  IF(MOD(MacCoefs(IDY,3)%Size,4).NE.0) THEN
                     NumLines = MacCoefs(IDY,3)%Size/4 + 1
                  ELSE
                     NumLines = MacCoefs(IDY,3)%Size/4
                  END IF
                  INTTERM = 2
                  IDk = 1
                  WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                  WRITE (OUTunit,"(6X,A4,1X,6(I3))") XSNames(11+m),Ur,CB,zero,   &
                                                     XOrders(5),XOrders(4)
                  DO k=1,NumLines
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     IF(k.NE.NumLines) THEN
                        WRITE (OUTunit,400) (MacCoefs(IDY,3)%xVector(kk),kk=IDk,IDk+3)
                     ELSE
                        WRITE (OUTunit,400) (MacCoefs(IDY,3)%xVector(kk), &
                                             kk=IDk,MacCoefs(IDY,3)%Size)
                     END IF
                     IDk = IDk + 4 
                  END DO
                  IF(NumTerms.EQ.4) THEN
                     IF(MOD(MacCoefs(IDY,4)%Size,4).NE.0) THEN
                        NumLines = MacCoefs(IDY,4)%Size/4 + 1
                     ELSE
                        NumLines = MacCoefs(IDY,4)%Size/4
                     END IF
                     INTTERM = 3
                     IDk = 1
                     WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                     WRITE (OUTunit,"(5X,A1,A4,1X,6(I3))") "-",XSNames(11+m),Ur,CB,zero, &
                                                           XOrders(7),XOrders(6)
                     DO k=1,NumLines
                        WRITE (OUTunit,300,ADVANCE="NO") INTTERM,n,j,i
                        IF(k.NE.NumLines) THEN
                           WRITE (OUTunit,400) (0.0D0-MacCoefs(IDY,4)%xVector(kk),kk=IDk,IDk+3)
                        ELSE
                           WRITE (OUTunit,400) (0.0D0-MacCoefs(IDY,4)%xVector(kk), &
                                                kk=IDk,MacCoefs(IDY,4)%Size)
                        END IF
                        IDk = IDk + 4 
                     END DO
                  END IF
               END DO
            END DO
            DO m=1,NumTerms
               DO n=1,NumMacs
                  CALL LinearSystem_Void(MacCoefs(n,m))
               END DO
            END DO
            DEALLOCATE (MacCoefs)
            DEALLOCATE (XOrders)
            CLOSE (INunitG)
         END DO
      END DO
      CLOSE (OUTunit)
300   FORMAT (4(I3))    
400   FORMAT (4X,4(1X,ES12.5))
500   FORMAT (1X,4(1X,ES12.5))
      DEALLOCATE (TypeNs,NMics,NumMicParas)
      DEALLOCATE (MacParaType,MacLength)
      DEALLOCATE (MicParaType,MicLength)
      IF(SigSExist) THEN
         DEALLOCATE (IGstart,IGend,SelfPosition)
      END IF
   END SUBROUTINE Link_ORange

   !-----------------------------------------------------------------------------------------------
   ! Calculate Y value with term values
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE Link_YfromTerms(User_Link,Yvar,IndexY)
      IMPLICIT NONE
      ! Passed in
      TYPE (Link_Type) User_Link
      INTEGER,INTENT(IN) :: IndexY
      ! Passed out 
      REAL(8)            :: Yvar
      ! Locals
      INTEGER            :: i,j,k,IDPN,IDPP,NumTerms,NLIT
      INTEGER            :: Parentheses(User_Link%NumTerms)
      INTEGER            :: Operator(User_Link%NumTerms)
      REAL(8)            :: RSLT
      LOGICAL            :: NFOUND,PFOUND
      NumTerms = User_Link%NumTerms
      DO i=1,User_Link%NumTerms
         Parentheses(i) = User_Link%Parentheses(i)
         Operator(i) = User_Link%Operator(i)
      END DO
      RSLT  = 0.0
      NLIT  = 0
      IF(NumTerms.EQ.1) THEN
         RSLT = User_Link%Terms(IndexY,1)
      END IF
      DO WHILE (NumTerms.GT.1)
         i = 1 
         NFOUND = .FALSE.
         ! Find the first close-parenthese
         DO WHILE ((NFOUND==.FALSE.).AND.(i<=User_Link%NumTerms))
            IF (Parentheses(i)<0) THEN
               IDPN = i
               NFOUND = .TRUE.
            END IF
            i = i+1
         END DO
         i = IDPN-1
         PFOUND = .FALSE.
         ! Find the first open-parenthese before the first close-parenthese found
         DO WHILE ((PFOUND==.FALSE.).AND.(i>=1))
            IF (Parentheses(i)>0) THEN
               IDPP = i
               PFOUND = .TRUE.
            END IF 
            i = i-1
         END DO
         ! Get the calculation result between the actual parentheses
         !WRITE(*,*) 'ID',NLIT,IDPN,IDPP
         RSLT = User_Link%Terms(IndexY,IDPP)
         !WRITE(*,*) RSLT
         DO i=IDPP+1,IDPN
            IF (Operator(i)==1) THEN
               IF ((i+1)<=IDPN) THEN 
                  IF (Operator(i+1)==3) THEN 
                     User_Link%Terms(IndexY,i+1) = User_Link%Terms(IndexY,i+1)                  &
                     *User_Link%Terms(IndexY,i)
                     Operator(i+1)   = 1
                     User_Link%Terms(IndexY,i) = 0
                  END IF
               END IF
               IF ((i+1)<=IDPN) THEN
                  IF (Operator(i+1)==4) THEN
                     User_Link%Terms(IndexY,i+1) = User_Link%Terms(IndexY,i)                    &
                     /User_Link%Terms(IndexY,i+1)
                     Operator(i)     = 1
                     User_Link%Terms(IndexY,i) = 0
                  END IF
               END IF
               RSLT = RSLT+User_Link%Terms(IndexY,i)
            ELSE IF (Operator(i)==2) THEN
               IF ((i+1)<=IDPN) THEN                  
                  IF (Operator(i+1)==3) THEN                  
                     User_Link%Terms(IndexY,i+1) = User_Link%Terms(IndexY,i+1)                  &
                     *User_Link%Terms(IndexY,i)
                     Operator(i+1)   = 2
                     User_Link%Terms(IndexY,i) = 0
                  END IF
               END IF
               IF ((i+1)<=IDPN) THEN
                  IF (Operator(i+1)==4) THEN
                     User_Link%Terms(IndexY,i+1) = User_Link%Terms(IndexY,i)                    &
                     /User_Link%Terms(IndexY,i+1)
                     Operator(i)     = 2
                     User_Link%Terms(IndexY,i) = 0
                  END IF
               END IF
               RSLT = RSLT-User_Link%Terms(IndexY,i)
            ELSE IF (Operator(i)==3) THEN
               RSLT = RSLT*User_Link%Terms(IndexY,i)
            ELSE IF (Operator(i)==4) THEN
               RSLT = RSLT/User_Link%Terms(IndexY,i)
            END IF
            !WRITE(*,*) i,RSLT             
         END DO
         User_Link%Terms(IndexY,IDPP) = RSLT
         Parentheses(IDPP)  = Parentheses(IDPP)+Parentheses(IDPN)
         ! Advance all the values in the vector so that in the next iteration fewer calculations &
         ! will be done 
         DO i=IDPN+1,NumTerms
            User_Link%Terms(IndexY,i-(IDPN-IDPP)) = User_Link%Terms(IndexY,i)
            Parentheses(i-(IDPN-IDPP))  = Parentheses(i)
            Operator(i-(IDPN-IDPP))     = Operator(i)
         END DO
         DO i=NumTerms,-1,NumTerms-(IDPN-IDPP)+1
            User_Link%Terms(IndexY,i) = 0.0
            Parentheses(i)  = 0
            Operator(i)     = 0
         END DO
         NumTerms = NumTerms-(IDPN-IDPP)
         !WRITE(*,*)  NumTerms
         NLIT = NLIT+1
         !WRITE (*,*) (User_Link%Terms(IndexY,i),i=1,NumTerms)
         !WRITE (*,*) (User_Link%Parentheses(i),i=1,NumTerms)
      END DO
      Yvar = RSLT
      !WRITE(*,*)  Yvar
   END SUBROUTINE Link_YfromTerms
    
   !-----------------------------------------------------------------------------------------------
   ! Calculate term value at given X with coefficients read
   !-----------------------------------------------------------------------------------------------
    SUBROUTINE Link_YfromXC(XOrders,C,XVars,Yvalue,NumCoefs,NumXvars)
       IMPLICIT NONE
       ! Passed in 
       INTEGER,INTENT(IN)  :: NumXvars,NumCoefs
       INTEGER,DIMENSION(NumXvars),INTENT(IN) :: XOrders
       REAL(8),INTENT(IN)  :: C(NumCoefs),XVars(NumXvars)
       ! Passed out
       REAL(8),INTENT(OUT) :: Yvalue
       ! Locals
       REAL(8),ALLOCATABLE :: BasicFunction(:)
       INTEGER             :: i,j
!#MA The calculation part has been replaced according to the work of Yaodong HE.
!    This replacement mainly aims at the improvement of efficiency.
       INTEGER             :: k,m,N,KK,IDE,Orderi,ID,IDB,Multi,NLit
       REAL(8)             ::T1,T2,XX
       REAL(8)             ::XXMatrix(NumXvars,NumCoefs)
       INTEGER             :: Order1,Order2,Order3,Order4,L
       REAL(8)             ::YValue1,YValue2,YValue3,YValue4,X1,X2,X3,X4

        
        IF(NumXvars==4) THEN
            Order1=XOrders(1)
            Order2=XOrders(2)+1
            Order3=XOrders(3)+1
            Order4=XOrders(4)+1
            L=(Order1+1)*Order2*Order3*Order4
            X1=Xvars(1)
            X2=Xvars(2)
            X3=Xvars(3)
            X4=Xvars(4)
            YValue4=0
            DO I=1,Order4
                YValue3=0
                DO J=1,Order3
                    YValue2=0
                     DO K=1,Order2  
                        YValue1=C(L)
                        L=L-1
                         DO M=1,Order1
                             YValue1=YValue1*X1+C(L)
                             L=L-1
                         END DO
                         YValue2=YValue2*X2+YValue1                
                     END DO
                     YValue3=YValue3*X3+YValue2
                END DO
                YValue4=YValue4*X4+YValue3
            END DO
            YValue=YValue4
        ELSE IF(NumXvars==3) THEN
            Order1=XOrders(1)
            Order2=XOrders(2)+1
            Order3=XOrders(3)+1
            L=(Order1+1)*Order2*Order3
            YValue3=0
            X1=Xvars(1)
            X2=Xvars(2)
            X3=Xvars(3)
            DO J=1,Order3
                YValue2=0
                DO K=1,Order2  
                    YValue1=C(L)
                    L=L-1
                    DO M=1,Order1
                        YValue1=YValue1*X1+C(L)
                        L=L-1
                    END DO
                    YValue2=YValue2*X2+YValue1                
                END DO
                YValue3=YValue3*X3+YValue2
            END DO
            YValue=YValue3
        ELSE IF(NumXvars==2) THEN
            YValue2=0
            Order1=XOrders(1)
            Order2=XOrders(2)+1
            L=(Order1+1)*Order2
            X1=Xvars(1)
            X2=Xvars(2)
            DO K=1,Order2  
                YValue1=C(L)
                L=L-1
                DO M=1,Order1
                    YValue1=YValue1*X1+C(L)
                    L=L-1
                END DO
                YValue2=YValue2*X2+YValue1                
            END DO
            YValue=YValue2
        ELSE IF(NumXvars==1) THEN
            Order1=XOrders(1)
            L=Order1+1
            X1=Xvars(1)
            YValue1=C(L)
            L=L-1
            DO M=1,Order1
                YValue1=YValue1*X1+C(L)
                L=L-1
            END DO
            YValue=YValue1
        ELSE
            WRITE(*,*) "Supported number of x variables is up to 4."
            STOP
        END IF
       !YValue = 0
       !ALLOCATE (BasicFunction(NumCoefs))
       !DO i=1,NumCoefs
       !   BasicFunction(i)=XVars(1)**User_YXdiscrete%OrderMatrix(1,i)
       !   DO j=2,NumXvars
       !      BasicFunction(i)=BasicFunction(i)*(XVars(j)**User_YXdiscrete%OrderMatrix(j,i))
       !   END DO
       !   YValue = YValue+C(i)*BasicFunction(i)
       !END DO
       !DEALLOCATE (BasicFunction)
    END SUBROUTINE Link_YfromXC


   !-----------------------------------------------------------------------------------------------
   ! Determine the A matrix with X values used for fitting
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_AmatrixFromX(User_YXdiscrete,Amatrix,Asize)
      IMPLICIT NONE
      ! Passed in
      TYPE (YXdiscrete)  User_YXdiscrete
      INTEGER,INTENT(IN)  :: Asize
      ! Passed out
      REAL(8)             :: Amatrix(Asize,Asize)
      ! Locals
      INTEGER             :: i,j,k,JA
      REAL(8),ALLOCATABLE :: FuncVari(:),BaseFunc(:)
      ! Elements of A:
      ALLOCATE (FuncVari(User_YXdiscrete%NumXvars),BaseFunc(Asize))
      ! Inisialization
      DO i=1,Asize
         DO j=1,Asize
            Amatrix(i,j) = 0.0D0
         END DO
      END DO
      ! Obtain elements of A by using ordermatrix composing basic function terms
      DO k=1,User_YXdiscrete%NumStates
         DO i=1,User_YXdiscrete%NumXvars
            FuncVari(i) = User_YXdiscrete%X(k,i)
         END DO
         DO j=1,Asize
            BaseFunc(j) = FuncVari(1)**User_YXdiscrete%OrderMatrix(1,j)
            DO i=2,User_YXdiscrete%NumXvars
               BaseFunc(j) = BaseFunc(j)*(FuncVari(i)**User_YXdiscrete%OrderMatrix(i,j))
            END DO
         END DO
         DO j=1,Asize
            DO JA=j,Asize ! Calculs effected only for the up-side triangle since A symetrique
               Amatrix(j,JA) = Amatrix(j,JA)+BaseFunc(j)*BaseFunc(JA)
            END DO
         END DO
      END DO
      ! Set the left-side triangle of A by using its symetry
      DO i=2,Asize
         DO j=1,i-1
            Amatrix(i,j) = Amatrix(j,i)
         END DO
      END DO
      ! Amatrix obtained, User_YXdiscrete%X should be voided
      ! CALL YXdiscrete_XVoid(User_YXdiscrete)
      DEALLOCATE (FuncVari,BaseFunc)       
   END SUBROUTINE YXdiscrete_AmatrixFromX

   !-----------------------------------------------------------------------------------------------
   ! Get b Vector from Y
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_bVectorFromY(User_YXdiscrete,bVector,bSize)
      IMPLICIT NONE
      ! Passed in
      TYPE(YXdiscrete),INTENT(IN) :: User_YXdiscrete
      ! Passed out       
      INTEGER             :: bSize
      REAL(8)             :: bVector(bSize)
      ! Locals
      INTEGER             :: i,j,k
      REAL(8),ALLOCATABLE :: FuncVari(:),BaseFunc(:)
      ! Elements of b:
      bSize = User_YXdiscrete%NumCoefs 
      ALLOCATE (FuncVari(User_YXdiscrete%NumXvars),BaseFunc(bSize))
      ! Inisialization
      DO i=1,bSize
         bVector(i) = 0.0D0
      END DO
      ! Actually, b(j) = sum(basicfunction(k)*y(k)),k=1,NumStates
      DO k=1,User_YXdiscrete%NumStates
         ! Get basic functions
         DO i=1,User_YXdiscrete%NumXvars
            FuncVari(i) = User_YXdiscrete%X(k,i)
         END DO
         DO j=1,bSize
            BaseFunc(j) = FuncVari(1)**User_YXdiscrete%OrderMatrix(1,j)
            DO i=2,User_YXdiscrete%NumXvars
               BaseFunc(j) = BaseFunc(j)*(FuncVari(i)**User_YXdiscrete%OrderMatrix(i,j))
            END DO
         END DO
         ! Get b elements
         DO j=1,bSize
            bVector(j)  = bVector(j)+BaseFunc(j)*User_YXdiscrete%Y(k)
         END DO
      END DO       
      DEALLOCATE (FuncVari,BaseFunc)
   END SUBROUTINE YXdiscrete_bVectorFromY

   !-----------------------------------------------------------------------------------------------
   ! Get OrderMatrix from User_YXdiscrete
   ! OrderMatrix contains all possible combinations of X variables powers
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_OrderMatrix(User_YXdiscrete)
      IMPLICIT NONE
      ! Passed in
      TYPE (YXdiscrete)  User_YXdiscrete
      ! Locals
      INTEGER         :: i,j,k,m,IDE,Orderi,ID,IDB,Multi,NLit
      IDE = User_YXdiscrete%NumCoefs/(1+User_YXdiscrete%XOptions(1,2))
      DO i=1,IDE
         User_YXdiscrete%OrderMatrix(1,(User_YXdiscrete%XOptions(1,2)+1)*(i-1)+1) = 0
         DO j=1,User_YXdiscrete%XOptions(1,2)
            User_YXdiscrete%OrderMatrix(1,(User_YXdiscrete%XOptions(1,2)+1)*(i-1)+1+j) = j
         END DO
      END DO
      DO i=2,User_YXdiscrete%NumXvars
         Multi = 1+User_YXdiscrete%XOptions(1,2)
         DO k=2,i
            Multi = Multi*(User_YXdiscrete%XOptions(k,2)+1)
         END DO
         NLit = User_YXdiscrete%NumCoefs/Multi
         !WRITE(*,*) Multi,NLit
         DO k=1,NLit
            IDB=Multi*(k-1)
            DO m=1,User_YXdiscrete%XOptions(i,2)+1 !进入几个数字都等于m-1的循环
               ID=(m-1)*Multi/(User_YXdiscrete%XOptions(i,2)+1)+1
               DO j=ID,ID+Multi/(User_YXdiscrete%XOptions(i,2)+1)-1
                  User_YXdiscrete%OrderMatrix(i,IDB+j) = m-1
               END DO
            END DO
         END DO 
      END DO
   END SUBROUTINE YXdiscrete_OrderMatrix

   !-----------------------------------------------------------------------------------------------
   ! Reset fitting coefficients for facilitating interpolation
   ! Since original coefficients obtained demande interpolation with X treated, this function
   ! develops the equation and gets coefficients reset that can be used directly 
   ! without X treatment
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE Link_ResetCoefs(INunit,OUTunit,INFileName,OUTFileName)
      IMPLICIT NONE
      ! Passed in
      INTEGER             :: INunit,OUTunit
      CHARACTER(37)       :: INFileName
      CHARACTER(37)       :: OUTFileName
      ! Locals
      TYPE (Link_Type)    User_Link
      INTEGER,ALLOCATABLE :: Option(:),XOrder(:),IndexXvars(:),TOption(:),OrderEC(:)
      REAL(8),ALLOCATABLE :: AV(:),STD(:)
      REAL(8),ALLOCATABLE :: DEVCoef(:,:),CoefReset(:)
      REAL(8)             :: Multi,DIV
      INTEGER             :: SUMUP,IDDEV
      INTEGER             :: NumTerms,NumTXvars,NumTYRead,NOfY
      INTEGER             :: NumXinTerm,NumCoefs,NumStates,LengthDEV
      INTEGER             :: i,j,k,m,n,mm
      INTEGER             :: IDBG
      LOGICAL             :: SMLER = .FALSE.
      REAL(8),PARAMETER   :: NECP_RealZero  = 1.0E-18
      REAL(8),PARAMETER   :: NECP_RealEqual = 1.0E-8

      OPEN (INunit,file=INFileName)       
      OPEN (OUTunit,file=OUTFileName)       
      ! Read in term information
      READ (INunit,*)
      READ (INunit,*)
      READ (INunit,*)
      WRITE(OUTunit,*) '**************************************************************************'
      WRITE(OUTunit,*) 'NECP: Output File of Module Link, Coefficients reset'
      WRITE(OUTunit,*) '**************************************************************************'
      READ (INunit,*)  NumTerms,NumTXvars,NumTYRead
      WRITE(OUTunit,'(3(2X,I4))') NumTerms,NumTXvars,NumTYRead
      CALL Link_Define(User_Link,NumTYRead,NumTXvars,NumTerms)
      READ (INunit,*)  (User_Link%NameXvars(i),i=1,User_Link%NumTXvars)
      WRITE(OUTunit,'(5(2X,A8))') (User_Link%NameXvars(i),i=1,User_Link%NumTXvars)
      READ (INunit,*)  (User_Link%NameYvars(i),i=1,User_Link%NumTYvars)
      WRITE(OUTunit,'(5(2X,A8))') (User_Link%NameYvars(i),i=1,User_Link%NumTYvars)
      ALLOCATE (TOption(NumTerms))
      DO i=1,NumTerms
         READ (INunit,*)  User_Link%NameTerm(i),User_Link%Parentheses(i),User_Link%Operator(i), &
                          TOption(i)
         WRITE(OUTunit,*) User_Link%NameTerm(i),User_Link%Parentheses(i),User_Link%Operator(i), &
                          TOption(i)
      END DO
      DEALLOCATE (TOption)
      READ (INunit,*)
      WRITE(OUTunit,*) '*************************************************************************'
      READ (INunit,*)
      WRITE(OUTunit,*) 'X Y Variables'
      READ (INunit,*)
      WRITE(OUTunit,*) '*************************************************************************'
      DO i=1,User_Link%NumTerms
         ! Read in information of X and Y variables
         READ (INunit,*) NumXinTerm
         WRITE(OUTunit,'(I4)') NumXinTerm
         ALLOCATE(Option(NumXinTerm),XOrder(NumXinTerm),IndexXvars(NumXinTerm))
         READ (INunit,*)  (Option(j),j=1,NumXinTerm)
         WRITE(OUTunit,*) (Option(j),j=1,NumXinTerm)
         READ (INunit,*)  (XOrder(j),j=1,NumXinTerm)
         WRITE(OUTunit,*) (XOrder(j),j=1,NumXinTerm)
         READ (INunit,*)  (IndexXvars(j),j=1,NumXinTerm)
         WRITE(OUTunit,*) (IndexXvars(j),j=1,NumXinTerm)
         READ (INunit,*)
         WRITE(OUTunit,*) '**********************************************************************'
         ALLOCATE (AV(NumXinTerm),STD(NumXinTerm))
         DO j=1,NumXinTerm
            READ (INunit,*)  AV(j),STD(j)
            WRITE(OUTunit,*)'0.0','   ','1.0'
         END DO
         NumCoefs = 1
         LengthDEV = 0
         DO j=1,NumXinTerm
            NumCoefs = NumCoefs * (XOrder(j)+1)
            LengthDEV = LengthDEV + (XOrder(j)+1)
         END DO
         NumStates = 1
         CALL YXdiscrete_Define(User_Link%YX,NumXinTerm,NumStates,NumCoefs)
         CALL LinearSystem_Define(User_Link%LSystem,NumCoefs)
         DO j=1,NumXinTerm
            User_Link%YX%XOptions(j,1) = Option(j) 
            User_Link%YX%XOptions(j,2) = XOrder(j) 
         END DO
         DEALLOCATE (Option,IndexXvars)
         READ (INunit,*)
         WRITE(OUTunit,*) '**********************************************************************'
         READ (INunit,*)
         WRITE(OUTunit,*) 'ORDER MATRIX'
         READ (INunit,*) ((User_Link%YX%OrderMatrix(k,j),j=1,NumCoefs),k=1,NumXinTerm)
         READ (INunit,*)
         DO k=1,NumXinTerm
            DO j=1,NumCoefs-1            
               WRITE (OUTunit,'(I4)',ADVANCE='NO') User_Link%YX%OrderMatrix(k,j)
            END DO
            WRITE (OUTunit,'(I4)') User_Link%YX%OrderMatrix(k,j)
         END DO
         WRITE(OUTunit,*) '***************'
         ALLOCATE (DEVCoef(NumCoefs,LengthDEV),CoefReset(NumCoefs),OrderEC(NumXinTerm))
         ! Develop each element according to each coefficient, results stocked in DEVCoef
         DO m=1,NumCoefs
            DO n=1,1+XOrder(1)
               IF(User_Link%YX%OrderMatrix(1,m).GE.(n-1)) THEN
                  CALL GetCmn(User_Link%YX%OrderMatrix(1,m),(n-1),DEVCoef(m,n))
                  DEVCoef(m,n) = DEVCoef(m,n) * ((-AV(1))**(User_Link%YX%OrderMatrix(1,m)-n+1)) 
               ELSE
                  DEVCoef(m,n) = 0.0D0
               END IF
               IF(ABS(STD(1)-0.0).GT.NECP_RealEqual) THEN
                  DEVCoef(m,n) = DEVCoef(m,n) / (STD(1)**User_Link%YX%OrderMatrix(1,m))
               END IF
            END DO
            DO k=2,NumXinTerm
               IDBG = 0
               DO n=1,k-1
                  IDBG = 1 + XOrder(n) + IDBG
               END DO
               DO n=IDBG+1,IDBG+1+XOrder(k)
                  IF(User_Link%YX%OrderMatrix(k,m).GE.(n-IDBG-1)) THEN
                     CALL GetCmn(User_Link%YX%OrderMatrix(k,m),(n-IDBG-1),DEVCoef(m,n))
                     DEVCoef(m,n) = DEVCoef(m,n) * ((-AV(k))**                               &
                                    (User_Link%YX%OrderMatrix(k,m)-(n-IDBG-1)))
                  ELSE
                     DEVCoef(m,n) = 0.0D0
                  END IF
                  IF(ABS(STD(k)-0.0).GT.NECP_RealEqual) THEN
                     DEVCoef(m,n) = DEVCoef(m,n) / (STD(k)**User_Link%YX%OrderMatrix(k,m))
                  END IF
               END DO
            END DO
         END DO
         !DO m=1,NumCoefs
         !   DO n=1,LengthDEV
         !      write(*,*) 'DEV',i,m,n,DEVCoef(m,n)
         !   END DO
         !END DO
         ! Read in coefficients 
         DO j=1,NumTYRead
            !READ (INunit,*)
            !WRITE(OUTunit,*) '*************'
            READ (INunit,*)  NOfY
            WRITE(OUTunit,*) NOfY
            READ (INunit,*) (User_Link%LSystem%xVector(k),k=1,NumCoefs)
            READ (INunit,*)
            DO n=1,NumCoefs
               CoefReset(n) = 0.0D0
            END DO
            DO n=1,NumCoefs
               DO m=1,NumXinTerm
                  ! Get orders of the variables for CoefReset(n)
                  OrderEC(m) = User_Link%YX%OrderMatrix(m,n)
               END DO
               ! Skip cases whose last X values are surely smaller then in the actual case
               IDBG = NumCoefs/(1+XOrder(NumXinTerm)) * OrderEC(NumXinTerm) + 1
               ! Find cases in which searched X-power may exist
               DO m=IDBG,NumCoefs
                  SMLER = .FALSE.
                  DO k=1,NumXinTerm
                     IF(User_Link%YX%OrderMatrix(k,m).LT.OrderEC(k)) SMLER = .TRUE.
                     EXIT
                  END DO
                  IF(.NOT.SMLER) THEN
                     Multi = DEVCoef(m,OrderEC(1)+1)
                     !DIV   = STD(1)**OrderEC(1)
                     IDDEV = 2 + User_Link%YX%XOptions(1,2)
                     DO k=2,NumXinTerm
                        Multi = Multi * DEVCoef(m,IDDEV+OrderEC(k))
                        !DIV = DIV * (STD(k)**OrderEC(k))
                        IDDEV = IDDEV + User_Link%YX%XOptions(k,2)
                     END DO
                     CoefReset(n) = CoefReset(n) + User_Link%LSystem%xVector(m) * Multi !/ DIV
                  END IF
               END DO
            END DO
            WRITE (OUTunit,'(5(2X,ES15.8))') (CoefReset(k),k=1,NumCoefs)
            WRITE (OUTunit,*) '*****************'
         END DO
         DEALLOCATE(AV,STD,XOrder,CoefReset,DEVCoef,OrderEC)
         CALL YXdiscrete_Void(User_Link%YX)
         CALL LinearSystem_Void(User_Link%LSystem)
      END DO
      CALL Link_Void(User_Link)
   CLOSE (INunit)
   CLOSE (OUTunit)
   END SUBROUTINE Link_ResetCoefs

   !-----------------------------------------------------------------------------------------------
   ! Calculate C(m,n)
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE GetCmn(PBelow,PTop,CValue)
      IMPLICIT NONE
      ! Passed in
      INTEGER :: PBelow,PTop
      ! Passed out
      REAL(8) :: CValue
      ! Locals
      INTEGER :: i,j,k
      INTEGER :: MultiB,MultiT,MultiM
      IF(PTop.GT.PBelow) THEN
         WRITE(*,*) '[Link] FATAL ERROR: Below < Top!'
         RETURN
      END IF
      IF(PBelow.EQ.PTop) THEN
         CValue = 1.0D0
      ELSE
         MultiB = 1
         MultiT = 1
         MultiM = 1
         DO i=2,PBelow
            MultiB = MultiB * i
         END DO
         DO i=2,PTop
            MultiT = MultiT * i
         END DO
         DO i=2,PBelow-PTop
            MultiM = MultiM * i
         END DO
         CValue = MultiB / MultiT / MultiM
      END IF
   END SUBROUTINE

   !-----------------------------------------------------------------------------------------------
   ! Define object User_YXdiscrete
   !-----------------------------------------------------------------------------------------------
    SUBROUTINE YXdiscrete_Define(User_YXdiscrete,NumXvars,NumStates,NumCoefs)
       IMPLICIT NONE
       ! Passed in
       INTEGER,INTENT(IN) :: NumXvars,NumStates,NumCoefs
       TYPE (YXdiscrete) User_YXdiscrete
       IF((NumXvars.LE.0).OR.(NumStates.LE.0).OR.(NumCoefs.LE.0)) THEN
          WRITE (*,*) '[Link] FATAL ERROR: Cannot define YX, length incorrect'
          RETURN
       END IF
       IF (User_YXdiscrete%XDefined.OR.User_YXdiscrete%YDefined) THEN
          CALL YXdiscrete_Void(User_YXdiscrete)
       END IF
       User_YXdiscrete%NumXvars  = NumXvars
       User_YXdiscrete%NumStates = NumStates
       User_YXdiscrete%NumCoefs  = NumCoefs
       ALLOCATE (User_YXdiscrete%X(NumStates,NumXvars),User_YXdiscrete%XOptions(NumXvars,2),  &
                User_YXdiscrete%Y(NumStates),User_YXdiscrete%OrderMatrix(NumXvars,NumCoefs)) 
       User_YXdiscrete%YDefined  = .TRUE.
       User_YXdiscrete%XDefined  = .TRUE.
    END SUBROUTINE YXdiscrete_Define

   !-----------------------------------------------------------------------------------------------
   ! Define object User_Link
   !-----------------------------------------------------------------------------------------------
    SUBROUTINE Link_Define(User_Link,NumTYvars,NumTXvars,NumTerms)
       IMPLICIT NONE
       ! Passed in
       INTEGER,INTENT(IN) :: NumTYvars,NumTXvars,NumTerms
       TYPE (Link_Type) User_Link
       IF((NumTXvars.LE.0).OR.(NumTYvars.LE.0).OR.(NumTerms.LE.0)) THEN
          WRITE (*,*) '[Link] FATAL ERROR: Cannot define User_Link, length incorrect'
          RETURN
       END IF
       IF(User_Link%Defined) CALL Link_Void(User_Link)
       User_Link%NumTerms   = NumTerms
       User_Link%NumTYvars  = NumTYvars
       User_Link%NumTXvars  = NumTXvars
       ALLOCATE (User_Link%NameYvars(NumTYvars),User_Link%NameXvars(NumTXvars),               &
                User_Link%Xvars(NumTXvars),User_Link%NameTerm(NumTerms),                      &
                User_Link%Terms(NumTYvars,NumTerms),User_Link%Parentheses(NumTerms),          &
                User_Link%Operator(NumTerms))
       User_Link%Defined = .TRUE.
    END SUBROUTINE Link_Define

   !-----------------------------------------------------------------------------------------------
   ! Print User_YXdiscrete through OUTunit
   !-----------------------------------------------------------------------------------------------
    SUBROUTINE YXdiscrete_Print(User_YXdiscrete,OUTunit)
       IMPLICIT NONE
       ! Passed in
       INTEGER,INTENT(IN)           :: OUTunit
       TYPE (YXdiscrete),INTENT(IN) :: User_YXdiscrete
       ! Locals
       INTEGER i,j
       IF(User_YXdiscrete%YDefined.AND.User_YXdiscrete%XDefined) THEN
         WRITE(OUTunit,10)
         WRITE(OUTunit,20)
         WRITE(OUTunit,110) User_YXdiscrete%XDefined
         WRITE(OUTunit,120) User_YXdiscrete%YDefined
         WRITE(OUTunit,30) User_YXdiscrete%NumStates,User_YXdiscrete%NumCoefs,                &
                           User_YXdiscrete%NumXvars
         WRITE(OUTunit,40) 
         WRITE(OUTunit,50) ((User_YXdiscrete%OrderMatrix(i,j),j=1,User_YXdiscrete%NumCoefs),  &
                            i=1,User_YXdiscrete%NumStates)
         WRITE(OUTunit,80)
         WRITE(OUTunit,50) ((User_YXdiscrete%XOptions(i,j),j=1,2),i=1,User_YXdiscrete%NumXvars)
         WRITE(OUTunit,50) ((User_YXdiscrete%X(i,j),j=1,User_YXdiscrete%NumXvars),            &
                           i=1,User_YXdiscrete%NumStates)             
         WRITE(OUTunit,90)
         WRITE(OUTunit,50) (User_YXdiscrete%Y(i),i=1,User_YXdiscrete%NumStates)              
         WRITE(OUTunit,10)
      ELSE
         WRITE(*,60)
         WRITE(OUTunit,60)
         WRITE(*,70)
         WRITE(OUTunit,70)
      END IF
10    FORMAT("[Link:YXdiscrete] ",83("-")) 
20    FORMAT("[Link:YXdiscrete] Object User_YXdiscrete")
30    FORMAT("[Link:YXdiscrete] NumStates,NumCoefs,NumXvars:",I,I,I)
40    FORMAT("[Link:YXdiscrete] OrderMatrix")
50    FORMAT("[Link:YXdiscrete]",100(2X,1PE16.6))
60    FORMAT("[Link:YXdiscrete] FATAL ERROR !!!")
70    FORMAT("[Link:YXdiscrete] Sorry, UNDEFINED object User_YXdiscrete.")
80    FORMAT("[Link:YXdiscrete] XOptions,X")
90    FORMAT("[Link:YXdiscrete] Y")
110   FORMAT("[Link:YXdiscrete] X defined?",L)
120   FORMAT("[Link:YXdiscrete] Y defined?",L)
    END SUBROUTINE YXdiscrete_Print

   !-----------------------------------------------------------------------------------------------
   ! Read back object User_YXdiscrete through OUTunit
   !-----------------------------------------------------------------------------------------------
    SUBROUTINE YXdiscrete_Read(User_YXdiscrete,INunit,OUTunit)
       IMPLICIT NONE
       ! Passed in
       INTEGER,INTENT(IN)           :: INunit,OUTunit
       TYPE (YXdiscrete),INTENT(IN) :: User_YXdiscrete
       ! Locals
       INTEGER :: i,j
       INTEGER :: NumStates,NumCoefs,NumXvars
       LOGICAL :: XDefined = .FALSE.
       LOGICAL :: YDefined = .FALSE.
       READ (INunit,*)
       READ (INunit,*)
       READ (INunit,110) XDefined
       READ (INunit,120) YDefined
       IF(YDefined.AND.XDefined) THEN
         READ (INunit,30) NumStates,NumCoefs,NumXvars
         CALL YXdiscrete_Define(User_YXdiscrete,NumXvars,NumStates,NumCoefs)
         READ (INunit,*)
         READ (INunit,50) ((User_YXdiscrete%OrderMatrix(i,j),j=1,User_YXdiscrete%NumCoefs),  &
                            i=1,User_YXdiscrete%NumStates)
         READ (INunit,*)
         READ (INunit,50) ((User_YXdiscrete%XOptions(i,j),j=1,2),i=1,User_YXdiscrete%NumXvars)
         READ (INunit,50) ((User_YXdiscrete%X(i,j),j=1,User_YXdiscrete%NumXvars),            &
                           i=1,User_YXdiscrete%NumStates)             
         READ (INunit,*)
         READ (INunit,50) (User_YXdiscrete%Y(i),i=1,User_YXdiscrete%NumStates)              
         READ (INunit,*)
      ELSE
         WRITE(*,130)
         WRITE(OUTunit,130)
         WRITE(*,140)
         WRITE(OUTunit,140)
         WRITE(*,150)
         WRITE(OUTunit,150)
      END IF
30    FORMAT("[Link:YXdiscrete] NumStates,NumCoefs,NumXvars:",I,I,I)
50    FORMAT("[Link:YXdiscrete]",100(2X,1PE16.6))
80    FORMAT("[Link:YXdiscrete] XOptions,X")
110   FORMAT("[Link:YXdiscrete] X defined?",L)
120   FORMAT("[Link:YXdiscrete] Y defined?",L)
130   FORMAT("[Link:YXdiscrete] Non-fatal error !")   
140   FORMAT("[Link:YXdiscrete] UNDEFINED object User_YXdiscrete in file INunit !")   
150   FORMAT("[Link:YXdiscrete] Object User_YXdiscrete stays UNDEFINED !")   
    END SUBROUTINE YXdiscrete_Read

   !-----------------------------------------------------------------------------------------------
   ! Print User_Link through OUTunit
   !-----------------------------------------------------------------------------------------------
    SUBROUTINE Link_Print(User_Link,OUTunit)
       IMPLICIT NONE
       ! Passed in
       INTEGER,INTENT(IN)           :: OUTunit
       TYPE (Link_Type),INTENT(IN)  :: User_Link
       ! Locals       
       INTEGER                      :: i,j
       IF(User_Link%Defined) THEN
         WRITE(OUTunit,10)
         WRITE(OUTunit,20)
         WRITE(OUTunit,160) User_Link%Defined
         WRITE(OUTunit,30)  User_Link%NumTYvars,User_Link%NumTXvars,User_Link%NumTerms
         WRITE(OUTunit,40) 
         WRITE(OUTunit,50)  (User_Link%NameYvars(i),i=1,User_Link%NumTYvars)
         WRITE(OUTunit,60)
         WRITE(OUTunit,50)  (User_Link%NameXvars(i),i=1,User_Link%NumTXvars)
         WRITE(OUTunit,70) 
         WRITE(OUTunit,80)  (User_Link%Xvars(i),i=1,User_Link%NumTXvars)
         WRITE(OUTunit,90) 
         WRITE(OUTunit,50)  (User_Link%NameTerm(i),i=1,User_Link%NumTerms)
         WRITE(OUTunit,100) 
         WRITE(OUTunit,80)  ((User_Link%Terms(i,j),i=1,User_Link%NumTYvars),j=1,User_Link%NumTerms)
         WRITE(OUTunit,110) 
         WRITE(OUTunit,120) (User_Link%Parentheses(i),i=1,User_Link%NumTerms)
         WRITE(OUTunit,130) 
         WRITE(OUTunit,120) (User_Link%Parentheses(i),i=1,User_Link%NumTerms)
         CALL YXdiscrete_Print(User_Link%YX,OUTunit)
         CALL LinearSystem_Print(User_Link%LSystem,OUTunit)
         WRITE(OUTunit,10)         
      ELSE
         WRITE(*,140)
         WRITE(OUTunit,140)
         WRITE(*,150)
         WRITE(OUTunit,150)
      END IF
10    FORMAT("[Link] ",94("-")) 
20    FORMAT("[Link] Object User_Link")
30    FORMAT("[Link] NumTYvars,NumTXvars,NumTerms:",I,I,I)
40    FORMAT("[Link] NameYvars")
50    FORMAT("[Link]",100(2X,A))
60    FORMAT("[Link] NameXvars")
70    FORMAT("[Link] Xvars")
80    FORMAT("[Link]",100(2X,1PE16.6))
90    FORMAT("[Link] NameTerm")
100   FORMAT("[Link] Terms")
110   FORMAT("[Link] Parentheses")
120   FORMAT("[Link]",100(2X,I))
130   FORMAT("[Link] Operator")
140   FORMAT("[Link] FATAL ERROR !!!")
150   FORMAT("[Link] Sorry, UNDEFINED object User_Link.")
160   FORMAT("[Link] User_Link defined?",L)
    END SUBROUTINE Link_Print

   !-----------------------------------------------------------------------------------------------
   ! Read User_Link through INunit
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE Link_Read(User_Link,INunit,OUTunit)
      IMPLICIT NONE
      ! Passed in
      INTEGER,INTENT(IN) :: INunit,OUTunit
      ! Passed out
      TYPE (Link_Type)   :: User_Link
      ! Locals       
      INTEGER            :: i,j
      LOGICAL            :: Defined = .FALSE.
      INTEGER            :: NumTYvars,NumTXvars,NumTerms
      READ(INunit,*)
      READ(INunit,*)
      READ(INunit,160) Defined
      IF(Defined) THEN
         READ(INunit,30)  NumTYvars,NumTXvars,NumTerms
         READ(INunit,*) 
         CALL Link_Define(User_Link,NumTYvars,NumTXvars,NumTerms)
         READ(INunit,50)  (User_Link%NameYvars(i),i=1,User_Link%NumTYvars)
         READ(INunit,*)
         READ(INunit,50)  (User_Link%NameXvars(i),i=1,User_Link%NumTXvars)
         READ(INunit,*) 
         READ(INunit,80)  (User_Link%Xvars(i),i=1,User_Link%NumTXvars)
         READ(INunit,*) 
         READ(INunit,50)  (User_Link%NameTerm(i),i=1,User_Link%NumTerms)
         READ(INunit,*) 
         READ(INunit,80)  ((User_Link%Terms(i,j),i=1,User_Link%NumTYvars),j=1,User_Link%NumTerms)
         READ(INunit,*) 
         READ(INunit,120) (User_Link%Parentheses(i),i=1,User_Link%NumTerms)
         READ(INunit,*) 
         READ(INunit,120) (User_Link%Parentheses(i),i=1,User_Link%NumTerms)
         CALL YXdiscrete_Read(User_Link%YX,INunit,OUTunit)
         CALL LinearSystem_Read(User_Link%LSystem,INunit,OUTunit)
         READ(INunit,*) 
      ELSE
         WRITE(*,140)
         WRITE(OUTunit,140)
         WRITE(*,150)
         WRITE(OUTunit,150)
      END IF
10    FORMAT("[Link] ",94("-")) 
20    FORMAT("[Link] Object User_Link")
30    FORMAT("[Link] NumTYvars,NumTXvars,NumTerms:",I,I,I)
40    FORMAT("[Link] NameYvars")
50    FORMAT("[Link]",100(2X,A))
60    FORMAT("[Link] NameXvars")
70    FORMAT("[Link] Xvars")
80    FORMAT("[Link]",100(2X,1PE16.6))
90    FORMAT("[Link] NameTerm")
100   FORMAT("[Link] Terms")
110   FORMAT("[Link] Parentheses")
120   FORMAT("[Link]",100(2X,I))
130   FORMAT("[Link] Operator")
140   FORMAT("[Link] FATAL ERROR !!!")
150   FORMAT("[Link] Sorry, UNDEFINED object User_Link.")
160   FORMAT("[Link] User_Link defined?",L)
   END SUBROUTINE Link_Read

   !-----------------------------------------------------------------------------------------------
   ! Provide numbers of INTEGER, REAL, LOGICAL and character variables in object User_YXdiscrete
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE YXdiscrete_GetStorage(User_YXdiscrete,NumINT,NumREAL,NumLOG,NumCHR)
      IMPLICIT NONE
      ! Passed in 
      TYPE(YXdiscrete),INTENT(IN) :: User_YXdiscrete
      ! Passed out
      INTEGER                     :: NumINT,NumREAL,NumLOG,NumCHR
      NumINT  = 2
      NumREAL = 0
      NumLOG  = 2
      NumCHR  = 0
      IF (User_YXdiscrete%XDefined.AND.User_YXdiscrete%YDefined) THEN
         NumREAL = NumREAL+User_YXdiscrete%NumStates*User_YXdiscrete%NumXvars                &
                   + User_YXdiscrete%NumStates
         NumINT  = NumINT+User_YXdiscrete%NumXvars*2                                         &
                   + User_YXdiscrete%NumXvars*User_YXdiscrete%NumCoefs
      ELSE IF((.NOT.User_YXdiscrete%XDefined).AND.User_YXdiscrete%YDefined) THEN
         NumREAL = NumREAL+User_YXdiscrete%NumStates
         NumINT  = NumINT+User_YXdiscrete%NumXvars*2                                         &
                   + User_YXdiscrete%NumXvars*User_YXdiscrete%NumCoefs          
      END IF
   END SUBROUTINE YXdiscrete_GetStorage

   !-----------------------------------------------------------------------------------------------
   ! Provide numbers of INTEGER, REAL, LOGICAL and character variables in object User_LinearLink
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE Link_GetStorage(User_Link,NumINT,NumREAL,NumLOG,NumCHR)
   IMPLICIT NONE
      ! Passed in
      TYPE(Link_Type),INTENT(IN) :: User_Link
      ! Passed out
      INTEGER                    :: NumINT,NumREAL,NumLOG,NumCHR
      ! Locals
      INTEGER                    :: NINT,NREAL,NLOG,NCHR
      NumCHR  = 0
      NumINT  = 3
      NumLOG  = 1
      NumREAL = 0
      NINT    = 0
      NREAL   = 0
      NLOG    = 0
      NCHR    = 0
      IF (User_Link%YX%XDefined.OR.User_Link%YX%YDefined) THEN
         !CALL linearSystem_GetStorage(User_Link%YX,NumINT,NumREAL,NumLOG,NumCHR)
         NumINT  = NumINT+3+User_Link%NumTerms
         NumLOG  = NumLOG+1+User_Link%NumTerms*2
         NumREAL = NumREAL+User_Link%NumTerms+User_Link%NumTXvars
         NumCHR  = NumCHR+User_Link%NumTYvars+User_Link%NumTerms
         CALL YXdiscrete_GetStorage(User_Link%YX,NINT,NREAL,NLOG,NCHR)
         NumINT  = NumINT + NINT
         NumREAL = NumREAL + NREAL
         NumLOG  = NumLOG + NLOG
         NumCHR = NumCHR + NCHR
      END IF
   END SUBROUTINE Link_GetStorage


   !-----------------------------------------------------------------------------------------------
   ! Void object User_YXdiscrete except for X
   !-----------------------------------------------------------------------------------------------       
   SUBROUTINE YXdiscrete_Void(User_YXdiscrete) !OrderMatrix 
      IMPLICIT NONE
      ! Passed in
      TYPE (YXdiscrete) User_YXdiscrete
      IF (User_YXdiscrete%XDefined) THEN
         DEALLOCATE (User_YXdiscrete%X)
         User_YXdiscrete%XDefined  = .FALSE.
      END IF
      IF (User_YXdiscrete%YDefined) THEN
         DEALLOCATE (User_YXdiscrete%Y) 
         User_YXdiscrete%YDefined  = .FALSE.
      END IF
      DEALLOCATE (User_YXdiscrete%XOptions,User_YXdiscrete%OrderMatrix)
      User_YXdiscrete%NumXvars  = 0
      User_YXdiscrete%NumStates = 0
      User_YXdiscrete%NumCoefs  = 0
   END SUBROUTINE YXdiscrete_Void

   !-----------------------------------------------------------------------------------------------
   ! Void X of object User_YXdiscrete
   !-----------------------------------------------------------------------------------------------       
   SUBROUTINE YXdiscrete_XVoid(User_YXdiscrete)
      IMPLICIT NONE
      ! To deallocate arrays of User_YXdiscrete, a Link_Paramtr-type variable
      ! Passed in
      TYPE (YXdiscrete) User_YXdiscrete
      IF (User_YXdiscrete%XDefined) THEN
         ! User_YXdiscrete%NumXvars = 0
         DEALLOCATE (User_YXdiscrete%X)
         User_YXdiscrete%XDefined = .FALSE.
      END IF
   END SUBROUTINE YXdiscrete_XVoid

   !-----------------------------------------------------------------------------------------------
   ! Void object User_Link
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE Link_Void(User_Link)
      IMPLICIT NONE
      ! Passed in
      TYPE (Link_Type) User_Link
      IF (User_Link%Defined) THEN
         User_Link%NumTerms   = 0
         User_Link%NumTYvars  = 0
         User_Link%NumTXvars  = 0
         DEALLOCATE (User_Link%NameYvars,User_Link%NameXvars,User_Link%Xvars,                &
                    User_Link%NameTerm,User_Link%Terms,User_Link%Parentheses,                &
                    User_Link%Operator)
         User_Link%Defined = .FALSE.
      END IF
   END SUBROUTINE Link_Void
   
!-----------------------------------------------------------------------------------------------
! Guass-Sidel method for Ax=b from book of Xu Shiliang
!-----------------------------------------------------------------------------------------------
    SUBROUTINE AGSDL (A, B, N, X, EPS, L)
        INTEGER       :: N, L
        REAL(8), INTENT(IN OUT) :: A(N,N), B(N), X(N)
        REAL(8), INTENT(IN) :: EPS
        real(8)    :: T, S, P, Q
        INTEGER    :: I, J 
        
        DO 5 I=1,N
          IF (ABS(A(I,I))+1.0.EQ.1.0) THEN
            L=0
            WRITE(*,100)
            RETURN
          END IF
5       CONTINUE
100     FORMAT(1X,'  FAIL')
        
        L=100
        DO 10 I=1,N
10      X(I)=0.0
20      P=0.0
        L=L-1
        DO 50 I=1,N
          T=X(I)
          S=0.0
          DO 30 J=1,N
            IF (J.NE.I) S=S+A(I,J)*X(J)
30        CONTINUE
          X(I)=(B(I)-S)/A(I,I)
          Q=ABS(X(I)-T)/(1+ABS(X(I)))
          IF (Q.GT.P) P=Q
50      CONTINUE
        IF ((P.GE.EPS).AND.(L.NE.0)) GOTO 20
        IF (L.EQ.0) THEN
         WRITE(*,100)
        END IF
        RETURN
    
    END SUBROUTINE AGSDL
       

END MODULE Link
