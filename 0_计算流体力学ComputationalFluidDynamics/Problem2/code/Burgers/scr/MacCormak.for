C===================================================================
C     MACCORMAK --LUZUO
C===================================================================
      SUBROUTINE MACCORMAK
      USE MAIN
      INTEGER I 

      DO I = 1 , NT
          !WRITE(*,*) i
          CALL PREDICTOR
          CALL CORRECTOR
          T = T + DELTA_T
          WRITE(*,*) 'T=',t
          CALL OUTPUT
          !IF( T == 0.8 .OR. T == 1.6 .OR. T == 2.0 ) THEN
          !    CALL OUTPUT
          !END IF ! IF( T == 0.8 .OR. T == 1.6 .OR. T == 2.0 ) THEN
          
      END DO ! DO I = 1 , NT
      
      END ! SUBROUTINE MACCORMACK
      
C===================================================================
C     Predictor step --LUZUO
C===================================================================
      SUBROUTINE PREDICTOR
      USE MAIN
      INTEGER I 
      WRITE(*,*) 'BEGIN PRE'
      DO I = 1 , NX
          IF( I == 1 .OR. I == NX) THEN
              U_PRE(I) = U(I)
              F_PRE(I) = F(I)
          ELSE
              IF( I == NX - 1 ) THEN
                  ! 边界反射条件
                  DU_DT(I) = MU *( U(I+1)-2*U(I)+U(I-1) )/( DELTA_X**2 )
     1            - ( F(I+1)-F(I) ) / DELTA_X
              ELSE
                  DU_DT(I) = MU *( U(I+1)-2*U(I)+U(I-1) )/( DELTA_X**2 )
     1            - ( F(I+1)-F(I) ) / DELTA_X
              END IF ! IF( I == NX -1 )THEN
              !WRITE(*,*) ' F(I+1)-F(I) = ', ( F(I+1)-F(I))
              !WRITE(*,*) 'DU_DT(I)=',DU_DT(I)
              U_PRE(I) = U(I) + DU_DT(I) * DELTA_T
              F_PRE(I) = 0.5 * U_PRE(I) ** 2
              
              !WRITE(*,*) 'U_PRE=',U_PRE(I)
              !WRITE(*,*) 'F_PRE=',F_PRE(I)
          END IF ! IF( I == 1 .OR. I == NX) THEN
      END DO ! DO I = 1 , NX
      WRITE(*,*) 'END PRE'
      END ! SUBROUTINE PREDICTOR
      
C===================================================================
C     Corrector step --LUZUO
C===================================================================
      SUBROUTINE CORRECTOR
      USE MAIN
      INTEGER I 
      WRITE(*,*) 'BEGIN CORR'
      DO I = 2 , NX -1
          IF( I == 2 ) THEN
              ! 边界反射条件
              DU_DT_FIX(I) = MU*( U_PRE(I+1)-2*U_PRE(I)+U_PRE(I-1) )
     1        /( DELTA_X**2 ) - ( F_PRE(I)-F_PRE(I-1) ) / DELTA_X
              
          ELSE
              DU_DT_FIX(I) = MU*( U_PRE(I+1)-2*U_PRE(I)+U_PRE(I-1) )
     1        /( DELTA_X**2 ) - ( F_PRE(I)-F_PRE(I-1) ) / DELTA_X
          END IF ! IF( I == NX -1 )THEN
c          DF_DT_FIX(I) = U_PRE(I) * ( U_PRE(I) - U_PRE(I-1) ) /
c     1    DELTA_X
          
          DU_DT_AV(I) = 0.5 * ( DU_DT(I) + DU_DT_FIX(I) )
          !DF_DT_AV(I) = 0.5 * ( DF_DT(I) + DF_DT_FIX(I) )
          
          U(I) = U(I) + DU_DT_AV(I) * DELTA_T
          F(I) =0.5 * U(I) ** 2
          !WRITE(*,*) 'U=',U(I)
          !WRITE(*,*) 'F=',F(I)
      END DO ! DO I = 2 , NX -1
      WRITE(*,*) 'END CORR'
      END ! SUBROUTINE CORRECTOR