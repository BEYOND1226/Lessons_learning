C====================================================================
C             Program : Laval
C             Author  : Luzuo
C             Date    : 2021.6
C====================================================================
      
C====================================================================
C             Module Main -- luzuo
C====================================================================
      MODULE MAIN
      IMPLICIT NONE
      
      ! INPUT
      INTEGER ::  NX  ,CMAX
      REAL(KIND=8) :: P_0 , RHO_0 , P_E , NCFL , NCX, ERRMAX ,M1
      
      REAL(KIND=8) :: X_STEP  , AX_A ,GAMMA , A , ERR
      REAL(KIND=8),ALLOCATABLE :: RHO(:) , P(:) , V(:) , X(:) ,T_STEP(:)
      REAL(KIND=8),ALLOCATABLE :: U(:,:) , F(:,:) , H(:,:) , S(:,:)
      
      REAL(KIND=8),ALLOCATABLE ::DU_DT(:,:),DU_DT_PRE(:,:),DU_DT_AV(:,:)
      REAL(KIND=8),ALLOCATABLE :: P_PRE(:) , RHO_PRE(:) , V_PRE(:) 
      REAL(KIND=8),ALLOCATABLE :: U_PRE(:,:) , F_PRE(:,:)
      REAL(KIND=8),ALLOCATABLE :: H_PRE(:,:) , S_PRE(:,:)
      REAL(KIND=8),ALLOCATABLE :: M(:),FLUX(:)
      
      END ! MODULE MAIN
      
C====================================================================
C             Program Main -- luzuo
C====================================================================
      USE MAIN
      
      CALL INPUT
      CALL INITIALIZATION
      CALL MACCORMAK
      CALL OUTPUT

      END