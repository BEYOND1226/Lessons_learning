C====================================================================
C             GRIDGENERATION -- luzuo
C====================================================================
      SUBROUTINE GRIDGENERATION
      USE MAIN
      
      INTEGER :: I , J
      REAL(KIND = 8) :: X_TEMP , Y_TEMP , ERR_TEMP
      ! ³õÊ¼»¯
      DO I = 1,NX
          X(I,1) =  0.5*(1+cos(2 * PI / (NX-1) * ( I-1 )))
          X(I,NY) = 0.5 + 0.5 * NFAR * cos(2 * PI / (NX-1) * ( I-1 ))
          
          Y(I,1) = 0.6*(-0.1015*X(I,1)**4 + 0.2843*X(I,1)**3 - 
     1    0.3576*X(I,1)**2 - 0.1221*X(I,1) + 0.2969*X(I,1)**0.5)
          Y(I,NY)  = ((NFAR/2)**2 - (X(I,NY)-0.5)**2)**0.5

      END DO ! DO I = 1,NX
      
      DO I = 1 , NX / 2
          Y(I,1) = -Y(I,1)
          Y(I,NY) = -Y(I,NY)
      END DO ! DO I = 1 , NX / 2
      
      DO J = 2 , NY-1
          DO I = 1 , NX
              X(I,J) = (X(I,NY) - X(I,1))/(NY-1) * (J - 1) + X(I,1)
              Y(I,J) = (Y(I,NY) - Y(I,1))/(NY-1) * (J - 1) + Y(I,1)
              !WRITE(*,*) X(I,J)
          END DO ! DO I = 1 , NX
      END DO ! DO J = 2 , NY-1
      
      !CALL OUTPUT
      
      !µü´ú
      ERR = 1000.D0
      C = 0
      DO WHILE( ERR > ERRMAX .AND. C < CMAX)
          ERR = 0.D0
          DO I = 1 , NX - 1
              
              DO J = 2 , NY - 1
                  IF( I == 1 )THEN
                      X_ETA = ( X(I,J+1)-X(I,J-1) )/ ( 2 * ETA_STEP )
                      Y_ETA = ( Y(I,J+1)-Y(I,J-1) )/ ( 2 * ETA_STEP )
                      X_XI = ( X(I+1,J)-X(NX-1,J) )/ ( 2 * XI_STEP )
                      Y_XI = ( Y(I+1,J)-Y(NX-1,J) )/ ( 2 * XI_STEP )
                  
                      ALPHA = X_ETA**2 + Y_ETA**2
                      BETA = X_ETA * X_XI + Y_ETA * Y_XI
                      GAMMA = X_XI**2 + Y_XI**2
                      
                      X_TEMP = X(I,J)
                      Y_TEMP = Y(I,J)
                  
                      X(I,J) =  0.5 * ( ALPHA * ( X(I+1,J) + X(NX-1,J) )
     1                + GAMMA * ( X(I,J+1) + X(I,J-1) ) -
     1                0.5 * BETA * ( X(I+1,J+1) + X(NX-1,J-1) -
     1                X(NX-1,J+1) - X(I+1,J-1) ) ) / ( ALPHA + GAMMA )
                      Y(I,J) =  0.5 * ( ALPHA * ( Y(I+1,J) + Y(NX-1,J) )
     1                + GAMMA * ( Y(I,J+1) + Y(I,J-1) ) -
     1                0.5 * BETA * ( Y(I+1,J+1) + Y(NX-1,J-1) -
     1                Y(NX-1,J+1) - Y(I+1,J-1) ) ) / ( ALPHA + GAMMA )
                      
                      ERR = MAX( ABS( X_TEMP - X(I,J) )
     1                 , ABS( Y_TEMP - Y(I,J) ) , ERR )
                      
                      X(NX,J) = X(I,J)
                      Y(NX,J) = Y(I,J)
                  ELSE
                      X_ETA = ( X(I,J+1)-X(I,J-1) )/ ( 2 * ETA_STEP )
                      Y_ETA = ( Y(I,J+1)-Y(I,J-1) )/ ( 2 * ETA_STEP )
                      X_XI = ( X(I+1,J)-X(I-1,J) )/ ( 2 * XI_STEP )
                      Y_XI = ( Y(I+1,J)-Y(I-1,J) )/ ( 2 * XI_STEP )
                  
                      ALPHA = X_ETA**2 + Y_ETA**2
                      BETA = X_ETA * X_XI + Y_ETA * Y_XI
                      GAMMA = X_XI**2 + Y_XI**2
                      !WRITE(*,*) X_ETA,Y_ETA,X_XI,Y_XI
                      X_TEMP = X(I,J)
                      Y_TEMP = Y(I,J)
                  
                      X(I,J) =  0.5 * ( ALPHA * ( X(I+1,J) + X(I-1,J) )
     1                + GAMMA * ( X(I,J+1) + X(I,J-1) ) -
     1                0.5 * BETA * ( X(I+1,J+1) + X(I-1,J-1) -
     1                X(I-1,J+1) - X(I+1,J-1) ) ) / ( ALPHA + GAMMA )
                      Y(I,J) =  0.5 * ( ALPHA * ( Y(I+1,J) + Y(I-1,J) )
     1                + GAMMA * ( Y(I,J+1) + Y(I,J-1) ) -
     1                0.5 * BETA * ( Y(I+1,J+1) + Y(I-1,J-1) -
     1                Y(I-1,J+1) - Y(I+1,J-1) ) ) / ( ALPHA + GAMMA )
                      
                      ERR = MAX( ABS( X_TEMP - X(I,J) )
     1                 , ABS( Y_TEMP - Y(I,J) ) , ERR )
                  END IF ! IF( I == 1 )THEN
              END DO  ! DO J = 2 , NY - 1
              
          END DO ! DO I = 2 , NX - 1
          C = C + 1    
      END DO ! DO WHILE( ERR > ERRMAX .AND. C < CMAX)
      WRITE(*,*) 'END GRIDGENERATION'
      WRITE(*,*) 'ERR=',ERR,'C=',C
      END !SUBROUTINE GRIDGENERATION
      
C====================================================================
C             SOLVER -- luzuo
C====================================================================
      SUBROUTINE SOLVER
      USE MAIN
      INTEGER :: I , J
      REAL(KIND = 8) :: PHI_TEMP
      
      WRITE(*,*) 'BEGIN SOLVER'
      
      DO I = 1 , NX
          DO J = 1 , NY
              PHI(I,J) = V_FAR * X(I,J)
          END DO
      END DO ! DO I = 1 , NX
      
      !DO J = 2 , NY-1
      !    DO I = 1 , NX
      !        PHI(I,J) = (PHI(I,NY)-PHI(I,1))/(NY-1)*(J - 1) + PHI(I,1)
      !        !WRITE(*,*) X(I,J)
      !    END DO ! DO I = 1 , NX
      !END DO ! DO J = 2 , NY-1 
      
      ERR = 1.D0
      C = 0
      DO WHILE( ERR > ERRMAX .AND. C < CMAX )
          ERR = 0.D0
          DO I = 1 , NX - 1
              DO J = 2 , NY - 1
                  IF( I == 1 )THEN
                      !PHI_ETA = ( PHI(I,J+1)-PHI(I,J-1) )/( 2*ETA_STEP )
                      !PHI_XI = ( PHI(I+1,J)-PHI(NX-1,J) )/( 2*XI_STEP )
                  
                      !ALPHA = PHI_ETA**2
                      !BETA = PHI_ETA * PHI_XI
                      !GAMMA = PHI_XI**2
                      
                      X_ETA = ( X(I,J+1)-X(I,J-1) )/ ( 2 * ETA_STEP )
                      Y_ETA = ( Y(I,J+1)-Y(I,J-1) )/ ( 2 * ETA_STEP )
                      X_XI = ( X(I+1,J)-X(NX-1,J) )/ ( 2 * XI_STEP )
                      Y_XI = ( Y(I+1,J)-Y(NX-1,J) )/ ( 2 * XI_STEP )
                  
                      ALPHA = X_ETA**2 + Y_ETA**2
                      BETA = X_ETA * X_XI + Y_ETA * Y_XI
                      GAMMA = X_XI**2 + Y_XI**2
                      
                      PHI_TEMP = PHI(I,J)
                  
                      PHI(I,J) =  0.5 * ( ALPHA * 
     1                ( PHI(I+1,J) + PHI(NX-1,J) )
     1                + GAMMA * ( PHI(I,J+1) + PHI(I,J-1) ) -
     1                0.5 * BETA * ( PHI(I+1,J+1) + PHI(NX-1,J-1) -
     1                PHI(NX-1,J+1) - PHI(I+1,J-1) ) )/( ALPHA + GAMMA )
                      
                      ERR = MAX( ABS( PHI_TEMP - PHI(I,J) ) , ERR )
                      
                      PHI(NX,J) = PHI(I,J)
                  ELSE
                      !PHI_ETA = ( PHI(I,J+1)-PHI(I,J-1) )/( 2*ETA_STEP )
                      !PHI_XI = ( PHI(I+1,J)-PHI(I-1,J) )/( 2*XI_STEP )
                  
                      !ALPHA = PHI_ETA**2
                      !BETA = PHI_ETA * PHI_XI
                      !GAMMA = PHI_XI**2
                      X_ETA = ( X(I,J+1)-X(I,J-1) )/ ( 2 * ETA_STEP )
                      Y_ETA = ( Y(I,J+1)-Y(I,J-1) )/ ( 2 * ETA_STEP )
                      X_XI = ( X(I+1,J)-X(I-1,J) )/ ( 2 * XI_STEP )
                      Y_XI = ( Y(I+1,J)-Y(I-1,J) )/ ( 2 * XI_STEP )
                  
                      ALPHA = X_ETA**2 + Y_ETA**2
                      BETA = X_ETA * X_XI + Y_ETA * Y_XI
                      GAMMA = X_XI**2 + Y_XI**2 
                  
                      PHI_TEMP = PHI(I,J)
                  
                      PHI(I,J) =  0.5 * ( ALPHA * 
     1                ( PHI(I+1,J) + PHI(I-1,J) )
     1                + GAMMA * ( PHI(I,J+1) + PHI(I,J-1) ) -
     1                0.5 * BETA * ( PHI(I+1,J+1) + PHI(I-1,J-1) -
     1                PHI(I-1,J+1) - PHI(I+1,J-1) ) )/( ALPHA + GAMMA )
                      
                      ERR = MAX( ABS( PHI_TEMP - PHI(I,J) ) , ERR )
                  END IF ! IF( I == 1 )THEN
              END DO  ! DO J = 2 , NY - 1
              
              PHI_TEMP = PHI(I,1)
              IF( I == 1 )THEN
                  X_ETA = ( -3*X(I,1)+4*X(I,2)-X(I,3) )/( 2 * ETA_STEP )
                  Y_ETA = ( -3*y(I,1)+4*Y(I,2)-Y(I,3) )/( 2 * ETA_STEP )
                  X_XI = ( X(I+1,1)-X(NX-1,1) )/ ( 2 * XI_STEP )
                  Y_XI = ( Y(I+1,1)-Y(NX-1,1) )/ ( 2 * XI_STEP )
              
                  BETA = X_ETA * X_XI + Y_ETA * Y_XI
                  GAMMA = X_XI**2 + Y_XI**2
                  
                  PHI(I,1) = -1.D0/3.D0 * (
     1            ( beta/gamma * ETA_STEP/xi_STEP) *
     1            (phi(i+1,1)-phi(NX-1,1)) - 4*phi(i,2) + phi(i,3) )
                  PHI(NX,1) = PHI(I,1) 
              ELSE
                  X_ETA = ( -3*X(I,1)+4*X(I,2)-X(I,3))/ ( 2 * ETA_STEP )
                  Y_ETA = ( -3*y(I,1)+4*Y(I,2)-Y(I,3))/ ( 2 * ETA_STEP )
                  X_XI = ( X(I+1,1)-X(I-1,1) )/ ( 2 * XI_STEP )
                  Y_XI = ( Y(I+1,1)-Y(I-1,1) )/ ( 2 * XI_STEP )
              
                  BETA = X_ETA * X_XI + Y_ETA * Y_XI
                  GAMMA = X_XI**2 + Y_XI**2
                  
                  PHI(I,1) = -1.D0/3.D0 * (
     1            ( beta/gamma * ETA_STEP/xi_STEP) *
     1            (phi(i+1,1)-phi(I-1,1)) - 4*phi(i,2) + phi(i,3) )
              END IF !  IF( I == 1 )THEN
              ERR = MAX( ABS( PHI_TEMP - PHI(I,1) ) , ERR )
              
          END DO ! DO I = 2 , NX - 1
          C = C + 1
          
      END DO ! DO WHILE( ERR > ERRMAX .AND. C < CMAX )
      
      
      WRITE(*,*) 'ERR=',ERR,'C=',C
      
      DO I = 1 ,NX - 1
          DO J = 2 , NY-1
              IF ( I==1 ) THEN
                  X_ETA = ( X(I,J+1)-X(I,J-1) )/ ( 2 * ETA_STEP )
                  Y_ETA = ( Y(I,J+1)-Y(I,J-1) )/ ( 2 * ETA_STEP )
                  X_XI = ( X(I+1,J)-X(NX-1,J) )/ ( 2 * XI_STEP )
                  Y_XI = ( Y(I+1,J)-Y(NX-1,J) )/ ( 2 * XI_STEP )
                  PHI_ETA = ( PHI(I,J+1)-PHI(I,J-1) )/( 2*ETA_STEP )
                  PHI_XI = ( PHI(I+1,J)-PHI(NX-1,J) )/( 2*XI_STEP )
              ELSE
                  X_ETA = ( X(I,J+1)-X(I,J-1) )/ ( 2 * ETA_STEP )
                  Y_ETA = ( Y(I,J+1)-Y(I,J-1) )/ ( 2 * ETA_STEP )
                  X_XI = ( X(I+1,J)-X(I-1,J) )/ ( 2 * XI_STEP )
                  Y_XI = ( Y(I+1,J)-Y(I-1,J) )/ ( 2 * XI_STEP )
                  PHI_ETA = ( PHI(I,J+1)-PHI(I,J-1) )/( 2*ETA_STEP )
                  PHI_XI = ( PHI(I+1,J)-PHI(I-1,J) )/( 2*XI_STEP )
              END IF ! IF ( I==1 ) THEN
              
              JAC = x_xi*y_eta - x_eta*y_xi
              U(I,J) = 1.D0/JAC * (phi_xi*y_eta - phi_eta*y_xi)
              v(i,j) = 1.D0/JAC * (phi_eta*x_xi - phi_xi*x_eta)
              cp(i,j) = 1 - (u(i,j)**2 + v(i,j)**2)/(v_far**2)
              
          END DO ! DO J = 2 , NY - 1
          
          IF( I == 1 ) THEN
              X_ETA = ( -3*X(I,1)+4*X(I,2)-X(I,3) )/ ( 2 * ETA_STEP )
              Y_ETA = ( -3*y(I,1)+4*Y(I,2)-Y(I,3) )/ ( 2 * ETA_STEP )
              X_XI = ( X(I+1,1)-X(NX-1,1) )/ ( 2 * XI_STEP )
              Y_XI = ( Y(I+1,1)-Y(NX-1,1) )/ ( 2 * XI_STEP )
              PHI_ETA = ( -3*PHI(I,3)+4*PHI(I,2)-PHI(I,1) )/(2*ETA_STEP)
              PHI_XI = ( PHI(I+1,1)-PHI(NX-1,1) )/( 2*XI_STEP )
          ELSE
              X_ETA = ( -3*X(I,1)+4*X(I,2)-X(I,3) )/ ( 2 * ETA_STEP )
              Y_ETA = ( -3*y(I,1)+4*Y(I,2)-Y(I,3) )/ ( 2 * ETA_STEP )
              X_XI = ( X(I+1,1)-X(I-1,1) )/ ( 2 * XI_STEP )
              Y_XI = ( Y(I+1,1)-Y(I-1,1) )/ ( 2 * XI_STEP )
              PHI_ETA = ( -3*PHI(I,3)+4*PHI(I,2)-PHI(I,1) )/(2*ETA_STEP)
              PHI_XI = ( PHI(I+1,1)-PHI(I-1,1) )/( 2*XI_STEP )
          END IF !  IF( I == 1 ) THEN
          
          JAC = x_xi*y_eta - x_eta*y_xi
          U(I , 1) = 1.D0/JAC * (phi_xi*y_eta - phi_eta*y_xi)
          v(i,1) = 1.D0/JAC * (phi_eta*x_xi - phi_xi*x_eta)
          cp(i,1) = 1 - (u(i,1)**2 + v(i,1)**2)/(v_far**2)
          
          IF( I == 1 ) THEN
              X_ETA = ( 3*X(I,NY)-4*X(I,NY-1)+X(I,NY-2) )/( 2*ETA_STEP )
              Y_ETA = ( 3*Y(I,NY)-4*Y(I,NY-1)+Y(I,NY-2) )/( 2*ETA_STEP )
              X_XI = ( X(I+1,NY)-X(NX-1,NY) )/ ( 2 * XI_STEP )
              Y_XI = ( Y(I+1,NY)-Y(NX-1,NY) )/ ( 2 * XI_STEP )
              PHI_ETA = ( 3*PHI(I,NX)-4*PHI(I,NX-1)+PHI(I,NX-2) )
     1         /(2*ETA_STEP)
              PHI_XI = ( PHI(I+1,NY)-PHI(NX-1,NY) )/( 2*XI_STEP )
          ELSE
              X_ETA = ( 3*X(I,NY)-4*X(I,NY-1)+X(I,NY-2) )/( 2*ETA_STEP )
              Y_ETA = ( 3*Y(I,NY)-4*Y(I,NY-1)+Y(I,NY-2) )/( 2*ETA_STEP )
              X_XI = ( X(I+1,NY)-X(I-1,NY) )/ ( 2 * XI_STEP )
              Y_XI = ( Y(I+1,NY)-Y(I-1,NY) )/ ( 2 * XI_STEP )
              PHI_ETA = ( 3*PHI(I,NX)-4*PHI(I,NX-1)+PHI(I,NX-2) )
     1         /(2*ETA_STEP)
              PHI_XI = ( PHI(I+1,NY)-PHI(I-1,NY) )/( 2*XI_STEP )
          END IF !  IF( I == 1 ) THEN
          
          JAC = x_xi*y_eta - x_eta*y_xi
          U(I,NY) = 1.D0/JAC * (phi_xi*y_eta - phi_eta*y_xi)
          v(i,NY) = 1.D0/JAC * (phi_eta*x_xi - phi_xi*x_eta)
          cp(i,NY) = 1 - (u(i,NX)**2 + v(i,NX)**2)/(v_far**2)
                  
      END DO !  DO I = 1 ,NX-1
      
      DO J = 1 , NY
          U(NX,J) = U(1,J)
          v(NX,J) = V(1,J)
          cp(NX,J) = CP(1,J)
      END DO ! DO J = 1 , NY
      
      WRITE(*,*) 'END SOLVER'
      END  ! SUBROUTINE SOLVER
      