!--------------------------------------------------------------------------------------------------
!                           Copyright(c)2013 Xi'an Jiaotong University
!--------------------------------------------------------------------------------------------------
!
!  Developed by Nuclear Engineering Computational Physics (NECP) Laboratory.
!
!**************************************************************************************************
!  Version: v1.0                                                                    Date: 4/22/2013
!  Purpose:
!     This is a LinearSystem module for Link methods.
!
!  Authors list:
!     Name                       Unit                              Email
!     =========================  ================================  ================================
!     Shengnan Gao               XJTU NECP                         gaoshengnan1989@163.com
!     Yunzhao Li                 XJTU NECP                         Yunzhao@mail.xjtu.edu.cn
!
!  Subroutines list:
!     SUBROUTINE LinearSystem_Define(User_LinearSystem,Size)
!     SUBROUTINE LinearSystem_GetStorage(User_LinearSystem,NumINT,NumREAL,NumLOG,NumChar)
!     SUBROUTINE LinearSystem_Read(User_LinearSystem,INunit,OUTunit)
!     SUBROUTINE LinearSystem_Print(User_LinearSystem,OUTunit)
!     SUBROUTINE LinearSystem_Void(User_LinearSystem)
!     SUBROUTINE LinearSystem_SetA(User_LinearSystem,Matrix,Size)
!     SUBROUTINE LinearSystem_Setb(User_LinearSystem,Vector,Size)
!     SUBROUTINE LinearSystem_GetSize(User_LinearSystem,Size)
!     SUBROUTINE LinearSystem_GetAandb(User_LinearSystem,Matrix,Vector,Size)
!     SUBROUTINE LinearSystem_GetSolution(User_LinearSystem,Solution,Size)
!     SUBROUTINE LinearSystem_SolverJB(User_LinearSystem)
!     SUBROUTINE LinearSystem_SolverGS(User_LinearSystem)
!************************************************************************************************** 

MODULE LinearSystem
IMPLICIT NONE

TYPE,PUBLIC :: LinearSystem_Type
   LOGICAL             :: Defined  = .FALSE.
   INTEGER             :: Size     = 0           !Dimension of C 
   REAL(8),ALLOCATABLE :: Amatrix(:,:)           !(Size,Size) A, fitting matrix
   REAL(8),ALLOCATABLE :: bVector(:)             !(Size) b, fitting vector
   REAL(8),ALLOCATABLE :: xVector(:)             !(Size) x, coefficient vector to be calculated, A*x=b
END TYPE LinearSystem_Type

CONTAINS

   !-----------------------------------------------------------------------------------------------
   ! Define object User_LinearSystem
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_Define(User_LinearSystem,Size)
      IMPLICIT NONE
      !Passed in
      INTEGER Size
      TYPE(LinearSystem_Type) User_LinearSystem
      IF(Size.LE.0) THEN
         WRITE(*,*) '[LinearSystem] ERROR: Size must be positive when defining User_LinearSystem'
         RETURN
      ELSE IF(User_LinearSystem%Defined) THEN
         CALL LinearSystem_Void(User_LinearSystem)
      END IF
      User_LinearSystem%Size = Size
      ALLOCATE (User_LinearSystem%Amatrix(Size,Size),User_LinearSystem%bVector(Size),        &
               User_LinearSystem%xVector(Size))
      User_LinearSystem%Defined = .TRUE.
   END SUBROUTINE LinearSystem_Define

   !-----------------------------------------------------------------------------------------------
   ! Provide numbers of INTEGER, REAL, LOGICAL and character variables in object User_LinearSystem
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_GetStorage(User_LinearSystem,NumINT,NumREAL,NumLOG,NumChar)
      IMPLICIT NONE
      !Passed in
      TYPE (LinearSystem_Type) User_LinearSystem
      !Passed out
      INTEGER NumINT,NumREAL,NumLOG,NumChar
      NumChar = 0
      IF (User_LinearSystem%Defined) THEN
         NumINT  = 1
         NumREAL = User_LinearSystem%Size*User_LinearSystem%Size+User_LinearSystem%Size      &
                  +User_LinearSystem%Size
         NumLOG  = 1
      ELSE
         NumINT  = 0
         NumREAL = 0
         NumLOG  = 0
      END IF
   END SUBROUTINE LinearSystem_GetStorage

   !-----------------------------------------------------------------------------------------------
   ! Read back object User_LinearSystem through INunit
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_Read(User_LinearSystem,INunit,OUTunit)
      IMPLICIT NONE
      ! Passed in
      TYPE (LinearSystem_Type),INTENT (INOUT) :: User_LinearSystem
      INTEGER,INTENT(IN)                      :: INunit,OUTunit 
      ! Locals
      INTEGER                                 :: Size,i,j
   IF(User_LinearSystem%Defined) THEN
     WRITE(*,*)"[LinearSystem] Non-fatal error !"
     WRITE(*,*)"[LinearSystem] Object User_LinearSystem is going to be overlapped !"
     WRITE(OUTunit,*)"[LinearSystem] Non-fatal error !"
     WRITE(OUTunit,*)"[LinearSystem] Object User_LinearSystem is going to be overlapped !"
   END IF
      READ(INunit,*) 
      READ(INunit,*) 
      READ(INunit,20) Size 
      CALL LinearSystem_Define(User_LinearSystem,Size)
      READ(INunit,*)((User_LinearSystem%Amatrix(i,j),j=1,Size),i=1,Size)
      READ(INunit,*)(User_LinearSystem%bVector(i),i=1,Size)
      READ(INunit,*)(User_LinearSystem%xVector(i),i=1,Size)
      READ(INunit,*) 
20    FORMAT("[LinearSystem] Matrix/Vector size:",I)       
   END SUBROUTINE LinearSystem_Read

   !-----------------------------------------------------------------------------------------------
   ! Print User_LinearSystem through OUTunit
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_Print(User_LinearSystem,OUTunit)
      IMPLICIT NONE
      ! Passed in
      TYPE(LinearSystem_Type),INTENT (IN) :: User_LinearSystem
      INTEGER,INTENT(IN)                  :: OUTunit
      ! Locals
      INTEGER             :: Size,i,j
      REAL(8),ALLOCATABLE :: Amatrix(:,:),bVector(:),xVector(:)
      IF(User_LinearSystem%Defined) THEN
         CALL LinearSystem_GetSize(User_LinearSystem,Size)
         ALLOCATE(Amatrix(Size,Size),bVector(Size),xVector(Size))
         WRITE(OUTunit,10)
         WRITE(OUTunit,20)
         WRITE(OUTunit,30) User_LinearSystem%Size
         WRITE(OUTunit,40) 
         Amatrix = 0.0D0
         bVector = 0.0D0
         xVector = 0.0D0
         CALL LinearSystem_GetAandb(User_LinearSystem,Amatrix,bVector,Size)
         WRITE(OUTunit,50)((Amatrix(i,j),j=1,Size),i=1,Size)
         WRITE(OUTunit,50)(bVector(i),i=1,Size)
         WRITE(OUTunit,50)(xVector(i),i=1,Size)
         WRITE(OUTunit,10)
         DEALLOCATE(Amatrix,bVector,xVector)
      ELSE
         WRITE(*,60)
         WRITE(OUTunit,60)
         WRITE(*,70)
         WRITE(OUTunit,70)
      END IF
10    FORMAT("[LinearSystem] ",88("-")) 
20    FORMAT("[LinearSystem] Object User_LinearSystem")
30    FORMAT("[LinearSystem] Matrix/Vector size:",I)
40    FORMAT("[LinearSystem]       Amatrix    bVector   xVector")
50    FORMAT("[LinearSystem]",100(2X,1PE16.6))
60    FORMAT("[LinearSystem] FATAL ERROR !!!")
70    FORMAT("[LinearSystem] Sorry, UNDEFINED object User_LinearSystem.")
   END SUBROUTINE LinearSystem_Print

   !-----------------------------------------------------------------------------------------------
   ! Void object User_LinearSystem
   !-----------------------------------------------------------------------------------------------         
   SUBROUTINE LinearSystem_Void(User_LinearSystem)
   IMPLICIT NONE
      !Passed in
      TYPE (LinearSystem_Type) User_LinearSystem
      IF (User_LinearSystem%Defined) THEN
         User_LinearSystem%Size = 1
         DEALLOCATE (User_LinearSystem%Amatrix,User_LinearSystem%bVector,User_LinearSystem   &
                    %xVector)
         User_LinearSystem%Defined = .FALSE.
      END IF
   END SUBROUTINE LinearSystem_Void

   !-----------------------------------------------------------------------------------------------
   ! Set the array values in object User_LinearSystem%Amatix according to array Matrix(*,*)
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_SetA(User_LinearSystem,Matrix,Size)
      IMPLICIT NONE
      ! Passed in
      INTEGER :: Size
      REAL(8) :: Matrix(Size,Size)   
      ! Passed out
      TYPE (LinearSystem_Type) User_LinearSystem
      ! Locals
      INTEGER :: i,j
      IF(.NOT.User_LinearSystem%Defined) CALL LinearSystem_Define(User_LinearSystem,Size)
      DO i=1,Size
         DO j=1,Size
         User_LinearSystem%Amatrix(i,j) = Matrix(i,j)
         END DO
      END DO
   END SUBROUTINE LinearSystem_SetA

   !-----------------------------------------------------------------------------------------------
   ! Set the array values in object User_LinearSystem%bVector according to array Vector(*)
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_Setb(User_LinearSystem,Vector,Size)
      IMPLICIT NONE
      ! Passed in
      INTEGER :: Size
      REAL(8) :: Vector(Size)     
      ! Passed out
      TYPE (LinearSystem_Type) User_LinearSystem
      ! Locals
      INTEGER :: i
      IF(.NOT.User_LinearSystem%Defined) CALL LinearSystem_Define(User_LinearSystem,Size)
      DO i=1,Size
      User_LinearSystem%bVector(i) = Vector(i)
      END DO
   END SUBROUTINE LinearSystem_Setb

   !-----------------------------------------------------------------------------------------------
   ! Get the size of arrays from object User_LinearSystem
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_GetSize(User_LinearSystem,Size)
      IMPLICIT NONE
      ! Passed in
      TYPE (LinearSystem_Type) User_LinearSystem
      ! Passed out
      INTEGER Size
      IF(.NOT.User_LinearSystem%Defined) THEN
         WRITE(*,*) '[LinearSystem] ERROR: cannot GetSize, object not difined!'
         RETURN
      END IF
      Size = User_LinearSystem%Size
   END SUBROUTINE LinearSystem_GetSize

   !-----------------------------------------------------------------------------------------------
   ! Get Amatrix and bVector of object User_LinearSystem
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_GetAandb(User_LinearSystem,Matrix,Vector,Size)
      IMPLICIT NONE
      ! Passed in
      TYPE (LinearSystem_Type) User_LinearSystem
      ! Passed out
      INTEGER             :: Size
      REAL(8)             :: Matrix(Size,Size),Vector(Size)
      ! Locals
      INTEGER             :: i,j
      IF(.NOT.User_LinearSystem%Defined) THEN
         WRITE(*,*) '[LinearSystem] ERROR: cannot GetAandb, object not difined!'
         RETURN
      END IF
      IF(User_LinearSystem%Size==Size) THEN
         DO i=1,Size
            DO j=1,Size
               Matrix(i,j) = User_LinearSystem%Amatrix(i,j)
            END DO
            Vector(i) = User_LinearSystem%bVector(i)
         END DO
      ELSE
         WRITE(*,*) 'FAULT: Size must be equal to User_LinearSystem%Size'
      END IF
   END SUBROUTINE LinearSystem_GetAandb

   !-----------------------------------------------------------------------------------------------
   ! Get the solution of the linear system from object User_LinearSystem
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_GetSolution(User_LinearSystem,Solution,Size)
      IMPLICIT NONE
      ! Passed in 
      TYPE (LinearSystem_Type) User_LinearSystem
      INTEGER             :: Size  
      ! Passed out
      REAL(8)             :: Solution(Size)    
      ! Local
      INTEGER             :: i
      IF(.NOT.User_LinearSystem%Defined) THEN
         WRITE(*,*) '[LinearSystem] ERROR: cannot GetSolution, object not difined!'
         RETURN
      END IF
      IF(User_LinearSystem%Size/=Size) THEN
         WRITE(*,*) '[LinearSystem] ERROR: Size must be equal to User_LinearSystem%Size'
         RETURN
      END IF
      DO i=1,Size
         Solution(i) = User_LinearSystem%xVector(i)
      END DO
   END SUBROUTINE LinearSystem_GetSolution

   !-----------------------------------------------------------------------------------------------
   ! Solve the linear problem by Jacobi litteration
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_SolverJB(User_LinearSystem)
      IMPLICIT NONE
      ! Passed in 
      TYPE (LinearSystem_Type) User_LinearSystem
      ! Locals
      REAL(8),ALLOCATABLE :: X(:),Xplus1(:)
      INTEGER             :: i,j,NLIT,Size
      REAL(8)             :: RES,STD,SUMUP
      Size = User_LinearSystem%Size
      ALLOCATE (X(Size),Xplus1(Size))
      DO i=1,Size
        X(i) = 0
        Xplus1(i) = 0
      END DO
      RES  = 1.0
      STD  = 0.001
      NLIT = 0 
      DO WHILE (NLIT<10000.AND.RES>=STD)
         RES = 0.0
         DO i=1,Size
            SUMUP = 0.0 
            DO j=1,i-1
               SUMUP = SUMUP+User_LinearSystem%Amatrix(i,j)*X(j)
            END DO
            DO j=i+1,Size
               SUMUP = SUMUP+User_LinearSystem%Amatrix(i,j)*X(j)
            END DO
            Xplus1(i) = (User_LinearSystem%bVector(i)-SUMUP)/User_LinearSystem%Amatrix(i,i)
            RES  = RES+(Xplus1(i)-X(i))*(Xplus1(i)-X(i))
         END DO
         RES = RES**(0.5)
         DO i=1,Size
            X(i) = Xplus1(i)
         END DO
         NLIT = NLIT+1
      END DO
      ! WRITE(*,*) NLIT,RES
      DO i=1,Size
         User_LinearSystem%xVector(i) = Xplus1(i)
      END DO
      ! WRITE(*,*) User_LinearSystem%xVector
      DEALLOCATE (X,Xplus1)
   END SUBROUTINE LinearSystem_SolverJB

   !-----------------------------------------------------------------------------------------------
   ! Solve the linear problem by Gauss_Sied litteration
   !-----------------------------------------------------------------------------------------------
   SUBROUTINE LinearSystem_SolverGS(User_LinearSystem,RES,STD,NLIT)
      IMPLICIT NONE
      ! Passed in 
      REAL(8)             :: STD
      TYPE (LinearSystem_Type) User_LinearSystem
      ! Passed out
      REAL(8)             :: RES
      INTEGER             :: NLIT
      ! Locals
      REAL(8),ALLOCATABLE :: X(:),Xplus1(:)
      INTEGER             :: i,j,Size
      REAL(8)             :: SUMUP
      Size = User_LinearSystem%Size
      ALLOCATE (X(Size),Xplus1(Size))
      DO i=1,Size
        X(i) = 0
        Xplus1(i) = 0
      END DO
      RES  = 1.0
      NLIT = 0 
      DO WHILE (NLIT<10000.AND.RES>=STD)
         RES = 0.0
         DO i=1,Size
            SUMUP = 0.0 
            DO j=1,i-1
               SUMUP = SUMUP+User_LinearSystem%Amatrix(i,j)*Xplus1(j)
            END DO
            DO j=i+1,Size
               SUMUP = SUMUP+User_LinearSystem%Amatrix(i,j)*X(j)
            END DO
            IF(User_LinearSystem%Amatrix(i,i)==0) THEN
               WRITE(*,*) '[LinearSystem] SolverGS: NO ANSWER! A(i,i)=0'
            END IF
            Xplus1(i) = (User_LinearSystem%bVector(i)-SUMUP)/User_LinearSystem%Amatrix(i,i)
            RES = RES+(Xplus1(i)-X(i))*(Xplus1(i)-X(i))
         END DO
         RES  = RES**(0.5)
         DO i=1,Size
            X(i) = Xplus1(i)
         END DO
         NLIT = NLIT+1
      END DO
      ! WRITE(*,*) NLIT,RES
      DO i=1,Size
         User_LinearSystem%xVector(i) = Xplus1(i)
      END DO
      ! WRITE(*,*) User_LinearSystem%xVector
      DEALLOCATE (X,Xplus1)
   END SUBROUTINE LinearSystem_SolverGS

   !SUBROUTINE LinearSystem_SolverGMRES(User_LinearSystem,Solution,Size)
   !IMPLICIT NONE
   !END SUBROUTINE LinearSystem_SolverGMRES

END MODULE LinearSystem
