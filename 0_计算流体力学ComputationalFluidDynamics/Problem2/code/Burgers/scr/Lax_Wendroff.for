C===================================================================
C     Lax-Wendroff --LUZUO
C===================================================================
      SUBROUTINE LAX_WENDROFF
      USE MAIN
      INTEGER :: I , J , K
      
      DO I = 1 , NT
          
          DO J = 2 , NX - 1
              DU_DT(J) = MU * ( U(J+1)-2*U(J)+U(J-1) )/(DELTA_X**2)
     1        - ( F(J+1)-F(J-1) )/( 2*DELTA_X )
              
              IF ( J == NX - 1 ) THEN
                  DU3_DX2DT = MU * ( ( 2*U(J+1)-U(J) ) - 4 * U(J+1) + 
     1            6 * U(J) -4 * U(J-1) + U(J-2) )/(DELTA_X**4) - 
     1            ( ( 2*F(J+1)-F(J) ) - 2*F(J+1) + 2*F(J-1) - F(J-2) )
     1            / (2 * DELTA_X**3 )
                  DU2_DXDT = MU * ( ( 2*U(J+1)-U(J) ) - 2 * U(J+1) + 
     1            2 * U(J-1) - U(J-2) )/( 2 * DELTA_X**3 ) - 
     1            ( ( F(J+1)) - 2*F(J) +F(J-1) )/(DELTA_X**2)
              ELSE IF ( J == 2 ) THEN
                  DU3_DX2DT = MU * ( ( 2*U(J-1)-U(J) ) - 4 * U(J-1) + 
     1            6 * U(J) -4 * U(J+1) + U(J+2) )/(DELTA_X**4) - 
     1            ( -( 2*F(J-1)-F(J) ) + 2*F(J-1) - 2*F(J+1) + F(J+2) )
     1            / (2 * DELTA_X**3 )
                  DU2_DXDT = MU * ( -( 2*U(J-1)-U(J) ) + 2 * U(J-1) - 
     1            2 * U(J+1) + U(J+2) )/( 2 * DELTA_X**3 ) - 
     1            ( ( F(J+1)) - 2*F(J) +F(J-1) )/(DELTA_X**2)
              ELSE
                  DU3_DX2DT = MU * ( U(J-2) - 4 * U(J-1) + 
     1            6 * U(J) -4 * U(J+1) + U(J+2) )/(DELTA_X**4) - 
     1            ( F(J+2) - 2*F(J+1) + 2*F(J-1) - F(J-2) )
     1            / (2 * DELTA_X**3 )
                  DU2_DXDT = MU * ( -U(J-2) + 2 * U(J-1) - 
     1            2 * U(J+1) + U(J+2) )/( 2 * DELTA_X**3 ) - 
     1            ( ( F(J+1)) - 2*F(J) +F(J-1) )/(DELTA_X**2)
              END IF ! IF ( J = NX - 1 )
              
              DF2_DTDX = ( (U(J+1) - U(J-1))/( 2*DELTA_X ) )*(DU_DT(J))
     1        + U(J) * DU2_DXDT
              DU2_DT2 = MU * DU3_DX2DT - DF2_DTDX
              !DU2_DT2 = 2.0*U(J)*(DU_DX(J)**2) + 
     1        !(U(J)**2)*( U(J+1)-2*U(J)+U(J-1) )/(DELTA_X**2)
              U_NEXT(J) = U(J) + DU_DT(J) * DELTA_T + 
     1        DU2_DT2 * ( 0.5 * DELTA_T**2 )
          END DO ! DO J = 2 , NX - 1
          
          DO J = 2 , NX - 1
              U(J) = U_NEXT(J)
              F(J) = 0.5 * U(J)**2
          END DO ! DO J = 2 , NX - 1
          
          T = T + DELTA_T
          WRITE(*,*) 'T=',T
          CALL OUTPUT
      END DO ! DO I = 1 , NT
      
      
      END ! SUBROUTINE LAX_WENDROFF 