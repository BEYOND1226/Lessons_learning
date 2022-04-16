C====================================================================
C             SUBROUTINE INPUT -- luzuo
C====================================================================
      SUBROUTINE INPUT
      USE MAIN
      LOGICAL ALIVE
      
      INQUIRE( FILE = '.\indata\Lavar.in' , EXIST = ALIVE )
      IF ( ALIVE ) THEN
          OPEN(10 , FILE = '.\indata\Lavar.in' )
          READ(10,*)
          READ(10,*) P_0 , RHO_0 , P_E ,M1 
          READ(10,*)
          READ(10,*) NX , NCFL , NCX , ERRMAX , CMAX
          CLOSE(10)
          WRITE(*,*) 'P0=',P_0 ,'RHO_0=',RHO_0,'P_E=',P_E
          WRITE(*,*) 'NX=',NX ,'NCFL=',NCFL,'NCX=',NCX,'ERRMAX=',ERRMAX
          WRITE(*,*) 'CMAX=',CMAX
      ELSE
          WRITE(*,*) 'CAN NOT FIND LAVAR.IN !!!'
      END IF ! IF ( ALIVE ) THEN
      WRITE(*,*) 'END INPUT'
      
      END ! SUBROUTINE INPUT
      
C====================================================================
C             SUBROUTINE INITIALIZATION -- luzuo
C====================================================================
      SUBROUTINE INITIALIZATION
      USE MAIN
      INTEGER :: I
      REAL(KIND=8) :: E
      
      ALLOCATE( RHO(NX) , P(NX) , V(NX) , X(NX) , T_STEP(NX) )
      ALLOCATE( U(3,NX) , F(3,NX) , H(3,NX) , S(3,NX) )
      ALLOCATE( RHO_PRE(NX) , P_PRE(NX) , V_PRE(NX) )
      ALLOCATE( U_PRE(3,NX) , F_PRE(3,NX) , H_PRE(3,NX) , S_PRE(3,NX) )
      ALLOCATE( DU_DT(3,NX) , DU_DT_PRE(3,NX) , DU_DT_AV(3,NX) )
      ALLOCATE( M(NX),FLUX(NX) )
      
      GAMMA = 1.4D0
      X_STEP = ( 1 - (0.1) ) / (NX-1)
      WRITE(*,*) 'X_STEP=',X_STEP
      DO I = 1 , NX
          
          X(I) = 0.1 + (I-1) * X_STEP
          AX_A = 2.0 * X(I) / ( 0.5 + X(I)**2 )
          
          IF ( I == 1 ) THEN
              P(1) = P_0 / ( 1+ (GAMMA-1)/2 * M1**2 )**(GAMMA/(GAMMA-1))
              RHO(1) = RHO_0 / (P_0/P(1))**(1/GAMMA)
              V(1) = M1 * SQRT(GAMMA * P(1)/RHO(1))
          ELSE IF ( I == NX )THEN
              P(NX) = P_E
              !P(NX) = 2*P(NX-1)-P(NX-2)
              RHO(NX) = ( P(NX)/P(NX-1) )**( 1/GAMMA )*RHO(NX-1)
              V(NX) = 2 * V(NX-1) - V(NX-2)
          ELSE
              P(I) = P(1) - (I-1) * ( P(1) - P_E )/(NX-1)
              RHO(I) = ( P(I)/P_0 )**( 1.0/GAMMA ) * RHO_0
              V(I) = SQRT( 2*GAMMA/( GAMMA-1 )*(P_0/RHO_0-P(I)/RHO(I)))
          END IF
          
          
          E = P(I)/( RHO(I) * ( GAMMA-1 ) ) + 0.5 * V(I)**2 
          U(1,I) = RHO(I)
          U(2,I) = RHO(I) * V(I)
          U(3,I) = RHO(I) * E
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
      END DO ! DO I = 1 , NX
      WRITE(*,*) 'END INITIALIZATION'
      
      END ! SUBROUTINE INITIALIZATION
      
C====================================================================
C             SUBROUTINE OUTPUT -- luzuo
C====================================================================
      SUBROUTINE OUTPUT
      USE MAIN
      
      OPEN(100,FILE = '.\outdata\Laval.dat',STATUS = 'replace' )
      WRITE(100,*) 'Variables=' , 'X, ','Ma,','u,','P,','Rho,','Mflux' 
      WRITE(100,*) 'ZONE'
      DO I = 1 , NX
          WRITE(100,"(6F16.6)") X(I),M(I),V(I),P(I),RHO(I),FLUX(I)
      END DO
      CLOSE(100)
      
      
      
      END ! SUBROUTINE OUTPUT