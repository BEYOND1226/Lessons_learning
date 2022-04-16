C===================================================================
C     FUNCTION: 计算流体力学课程编程题2,求解Burgers方程
C     AUTHOR  ：luzuo
C     DATE    : 2021.5.4
C===================================================================
      
      
C===================================================================
C     MODULE MAIN --LUZUO
C===================================================================
      MODULE MAIN
      IMPLICIT NONE
      
      INTEGER :: NX , NT , NTECHNIQUE
      
      REAL(KIND = 8) :: T
      REAL(KIND = 8) :: DELTA_X , DELTA_T
      REAL(KIND = 8) :: MU
      
      REAL(KIND = 8),ALLOCATABLE :: U(:) , F(:) , X(:) 
      
      REAL(KIND = 8),ALLOCATABLE :: DU_DT(:) , DU_DX(:) 
      REAL(KIND = 8),ALLOCATABLE :: DF_DT(:) , DF_DX(:)
      
      ! FOR LAX_WENDROFF
      REAL(KIND = 8) :: DU2_DT2 , DU2_DXDT , DU3_DX2DT , DF2_DXDT
      REAL(KIND = 8),ALLOCATABLE :: U_NEXT(:)
      
      ! FOR MACCORMACK
      REAL(KIND = 8),ALLOCATABLE :: U_PRE(:) , F_PRE(:) 
      REAL(KIND = 8),ALLOCATABLE :: DU_DT_AV(:) , DU_DT_FIX(:)
      
      END ! MODULE MAIN
      
C===================================================================
C     MAIN --LUZUO
C===================================================================
      USE MAIN
      WRITE(*,*) 'BEGIN'
      CALL INPUT
      WRITE(*,*) 'END INPUT'
      CALL INITIALIZATION
      WRITE(*,*) 'END INITIALIZATION'

      IF( NTECHNIQUE == 0 )THEN
          CALL LAX_WENDROFF
          WRITE(*,*) 'END LAX_WENDROFF'
      ELSE
          CALL MACCORMAK
          WRITE(*,*) 'END MACCORMAK'
      END IF ! IF( NTECHNIQUE = 0 )
      
      CALL OUTPUT
      
      END ! MAIN
C===================================================================
C     Initialization --LUZUO
C===================================================================
      SUBROUTINE INITIALIZATION
      USE MAIN
      
      NX = 2.0 / DELTA_X + 1
      NT = 2.0 / DELTA_T
      WRITE(*,*) 'NX=' ,NX
      WRITE(*,*) 'NT=' ,NT
      ALLOCATE( U(NX) , F(NX) , X(NX) )
      ALLOCATE( DU_DT(NX) , DF_DT(NX) , DU_DX(NX) , DF_DX(NX) )
      ALLOCATE( DU_DT_AV(NX) , DU_DT_FIX(NX) , U_PRE(NX),F_PRE(NX) )
      ALLOCATE( U_NEXT(NX) )
      
      T = 0.D0
      DO I = 1 , NX
          IF( I == 1 )THEN
              X(I) = -1.D0
          ELSE IF( I == NX )THEN
              X(I) = 1.D0
          ELSE
              X(I) = X(I-1) + DELTA_X
          END IF ! IF( I == 1 )
          U(I) = X(I) * (- 0.5 )
          F(I) = 0.5 * U(I)**2
      END DO ! DO I = 1 , NX
      CALL OUTPUT
      END ! SUBROUTINE INITIALIZATION