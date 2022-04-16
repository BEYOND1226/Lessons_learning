C====================================================================
C             Program : Grid
C             Author  : Luzuo
C             Date    : 2021.6
C====================================================================
      
C====================================================================
C             Module Main -- luzuo
C====================================================================
      MODULE MAIN
      IMPLICIT NONE
      Real , parameter :: PI = 3.1415926
      
      ! INPUT
      INTEGER ::  NX ,NY, NFAR , CMAX 
      REAL(KIND=8) ::  ERRMAX,ETA_STEP,XI_STEP 
      
      INTEGER :: C
      REAL(KIND=8) :: ALPHA,BETA,GAMMA , X_ETA ,X_XI ,Y_XI,Y_ETA,ERR
      REAL(KIND=8),ALLOCATABLE :: X(:,:) , Y(:,:)
      
      REAL(KIND=8) :: V_FAR  , PHI_XI , PHI_ETA , JAC
      REAL(KIND=8),ALLOCATABLE :: PHI(:,:) ,U(:,:) ,V(:,:),CP(:,:)
      
      
      END ! MODULE MAIN
      
C====================================================================
C             Program Main -- luzuo
C====================================================================
      USE MAIN
      
      CALL INPUT
      CALL INITIALIZATION
      CALL GRIDGENERATION
      CALL SOLVER
      CALL OUTPUT

      END