C====================================================================
C             SUBROUTINE MACCORMAK -- luzuo
C====================================================================
      SUBROUTINE MACCORMAK
      USE MAIN
      INTEGER C
      
      WRITE(*,*) 'BEGIN MACCORMAK'
      C = 0
      ERR = 1000.D0
      DO WHILE ( ERR > ERRMAX .and. C < CMAX )
          ERR = 0.D0
          CALL PREDICTOR
          CALL CORRECTOR
          C = C + 1
          !WRITE(*,*) 'ERR=',ERR
          !WRITE(*,*) 'C=',C
      END DO ! DO WHILE ( ERR > ERRMAX )
      WRITE(*,*) 'END MACCORMAK'
      WRITE(*,*) 'ERR=',ERR
      WRITE(*,*) 'C=',C
      END ! SUBROUTINE MACCORMAK
      
      
C====================================================================
C             SUBROUTINE PREDICTOR -- luzuo
C====================================================================
      SUBROUTINE PREDICTOR
      USE MAIN
      INTEGER I
      
      DO I = 2 , NX - 1
          DO J = 1 ,3
              DU_DT(J,I) = - ( F(J,I+1) - F(J,I) ) / X_STEP - H(J,I)
              S(J,I) = NCX * ABS( P(I+1) - 2*P(I) + P(I-1) ) / 
     1        ( P(I+1) + 2*P(I) + P(I-1) ) * 
     1        ( U(J,I+1) - 2*U(J,I) + U(J,I-1) )
              U_PRE(J,I) = U(J,I) + DU_DT(J,I) * T_STEP(I) + S(J,I)
              !WRITE(*,*) DU_DT(J,I)
          END DO ! DO J = 1 ,3
          
          RHO_PRE(I) = U_PRE(1,I)
          V_PRE(I) = U_PRE(2,I) / U_PRE(1,I)
          P_PRE(I) = ( GAMMA - 1 ) * ( U_PRE(3,I) - 
     1    0.5 * RHO_PRE(I) * V_PRE(I)**2 )
          
      END DO ! DO I = 2 , NX - 1
      
      V_PRE(1) = 2 * V_PRE(2) - V_PRE(3)
      RHO_PRE(1) = ( ( GAMMA / ( GAMMA - 1 ) * P_0 / RHO_0 
     1    - 0.5 * V_PRE(1)**2 ) * ( ( GAMMA - 1) / GAMMA 
     1    * RHO_0**GAMMA / P_0 ) )**( 1 / ( GAMMA - 1 ) )
      P_PRE(1) = P_0 * ( RHO_PRE(1) / RHO_0 )**GAMMA
      
      E = P_PRE(1)/( RHO_PRE(1) * ( GAMMA-1 ) ) + 0.5 * V_PRE(1)**2
      U_PRE(1,1) = RHO_PRE(1)
      U_PRE(2,1) = RHO_PRE(1) * V_PRE(1)
      U_PRE(3,1) = RHO_PRE(1) * E
      
      P_PRE(NX) = P_E
      !P_PRE(NX) = 2*P_PRE(NX-1)-P_PRE(NX-2)
      RHO_PRE(NX) = ( P_PRE(NX)/P_PRE(NX-1) )**( 1/GAMMA )*RHO_PRE(NX-1)
      V_PRE(NX) = 2 * V_PRE(NX-1) - V_PRE(NX-2)
      
      E = P_PRE(NX)/( RHO_PRE(NX) * ( GAMMA-1 ) ) + 0.5 * V_PRE(NX)**2
      U_PRE(1,NX) = RHO_PRE(NX)
      U_PRE(2,NX) = RHO_PRE(NX) * V_PRE(NX)
      U_PRE(3,NX) = RHO_PRE(NX) * E
      
      DO I = 1 , NX
          AX_A = 2.0 * X(I) / ( 0.5 + X(I)**2 )
          E = P_PRE(I)/( RHO_PRE(I) * ( GAMMA -1 ) ) 
     1    + 0.5 * V_PRE(I)**2
          F_PRE(1,I) = RHO_PRE(I) * V_PRE(I)
          F_PRE(2,I) = RHO_PRE(I) * V_PRE(I)**2 + P_PRE(I)
          F_PRE(3,I) = ( RHO_PRE(I) * E + P_PRE(I) ) * V_PRE(I)
          H_PRE(1,I) = AX_A * RHO_PRE(I) * V_PRE(I)
          H_PRE(2,I) = AX_A * RHO_PRE(I) * V_PRE(I)**2
          H_PRE(3,I) = AX_A * ( RHO_PRE(I) * E + P_PRE(I) ) * V_PRE(I)
          !WRITE(*,*) 'RHO_pre(I) =',RHO_PRE(I)
          !WRITE(*,*) 'P_PRE(I)=',P_PRE(I)
          !WRITE(*,*) 'V_PRE(I)=',V_PRE(I)
      END DO ! DO I = 1 , NX
      
      END ! SUBROUTINE PREDICTOR

C====================================================================
C             SUBROUTINE Corrector -- luzuo
C====================================================================
      SUBROUTINE CORRECTOR
      USE MAIN
      
      DO I = 2 , NX - 1
          DO J = 1 ,3
              DU_DT_PRE(J,I) = - ( F_PRE(J,I) - F_PRE(J,I-1) ) 
     1        / X_STEP - H_PRE(J,I)
              
              S_PRE(J,I) = NCX * 
     1        ABS( P_PRE(I+1) - 2*P_PRE(I) + P_PRE(I-1) ) / 
     1        ( P_PRE(I+1) + 2*P_PRE(I) + P_PRE(I-1) ) * 
     1        ( U_PRE(J,I+1) - 2*U_PRE(J,I) + U_PRE(J,I-1) )
              
              DU_DT_AV(J,I) = 0.5 * ( DU_DT_PRE(J,I) + DU_DT(J,I) ) 
              
              U(J,I) = U(J,I) + DU_DT_AV(J,I) * T_STEP(I) + S_PRE(J,I)
              
              
              !WRITE(*,*) 'DU_DT_AV(J,I) =',DU_DT_AV(J,I)
          END DO ! DO J = 1 ,3
          ERR = ERR + ABS(DU_DT_AV(1,I)) / ( NX - 1 )
          RHO(I) = U(1,I)
          V(I) = U(2,I) / U(1,I)
          P(I) = ( GAMMA - 1 ) * ( U(3,I) - 
     1    0.5 * RHO(I) * V(I)**2 )
          !WRITE(*,*) 'RHO(I) =',RHO(I)
          !WRITE(*,*) 'P(I)=',P(I)
      END DO ! DO I = 2 , NX - 1
      
      V(1) = 2 * V(2) - V(3)
      RHO(1) = ( ( GAMMA / ( GAMMA - 1 ) * P_0 / RHO_0 
     1    - 0.5 * V(1)**2 ) * ( ( GAMMA - 1) / GAMMA 
     1    * RHO_0**GAMMA / P_0 ) )**( 1 / ( GAMMA - 1 ) )
      P(1) = P_0 * ( RHO(1) / RHO_0 )**GAMMA
      
      E = P(1)/( RHO(1) * ( GAMMA-1 ) ) + 0.5 * V(1)**2
      U(1,1) = RHO(1)
      U(2,1) = RHO(1) * V(1)
      U(3,1) = RHO(1) * E
      
      P(NX) = P_E
      !P(NX) = 2*P(NX-1)-P(NX-2)
      RHO(NX) = ( P(NX)/P(NX-1) )**( 1/GAMMA )*RHO(NX-1)
      V(NX) = 2 * V(NX-1) - V(NX-2)
      
      E = P(NX)/( RHO(NX) * ( GAMMA-1 ) ) + 0.5 * V(NX)**2
      U(1,NX) = RHO(NX)
      U(2,NX) = RHO(NX) * V(NX)
      U(3,NX) = RHO(NX) * E
      
      DO I = 1 , NX
          AX_A = 2.0 * X(I) / ( 0.5 + X(I)**2 )
          E = P(I)/( RHO(I) * ( GAMMA -1 ) ) 
     1    + 0.5 * V(I)**2
          F(1,I) = RHO(I) * V(I)
          F(2,I) = RHO(I) * V(I)**2 + P(I)
          F(3,I) = ( RHO(I) * E + P(I) ) * V(I)
          H(1,I) = AX_A * RHO(I) * V(I)
          H(2,I) = AX_A * RHO(I) * V(I)**2
          H(3,I) = AX_A * ( RHO(I) * E + P(I) ) * V(I)
          
          A = SQRT( GAMMA * P(I) / RHO(I) )
          T_STEP(I) = NCFL * X_STEP / ( V(I) + A )
          !WRITE(*,*) 'RHO(I) =',RHO(I)
          !WRITE(*,*) 'P(I)=',P(I)
          !WRITE(*,*) 'V(I)=',V(I)
          M(I) = V(I) / A
          FLUX(i) = U(2,I) * (0.5+V(I)**2)
      END DO ! DO I = 1 , NX
      
      END ! SUBROUTINE CORRECTOR