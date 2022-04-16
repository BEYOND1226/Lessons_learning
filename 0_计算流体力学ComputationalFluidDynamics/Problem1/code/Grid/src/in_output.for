C====================================================================
C             SUBROUTINE INPUT -- luzuo
C====================================================================
      SUBROUTINE INPUT
      USE MAIN
      LOGICAL ALIVE
      
      INQUIRE( FILE = '.\indata\Grid.in' , EXIST = ALIVE )
      IF ( ALIVE ) THEN
          OPEN(10 , FILE = '.\indata\Grid.in' )
          READ(10,*)
          READ(10,*) NX , NY , NFAR  , V_FAR
          READ(10,*)
          READ(10,*)  ERRMAX ,CMAX
          CLOSE(10)
          WRITE(*,*) 'NX=',NX ,'NY=',NY,'NFAR',NFAR
          WRITE(*,*) 'V_FAR=',V_FAR
          WRITE(*,*) 'CMAX=',CMAX
      ELSE
          WRITE(*,*) 'CAN NOT FIND LAVAR.IN !!!'
      END IF ! IF ( ALIVE ) THEN
      
      INQUIRE( FILE = '.\output\Grid.dat' , EXIST = ALIVE )
      IF( ALIVE ) THEN
          OPEN(100,FILE = '.\output\Grid.dat',STATUS = 'replace' )
          WRITE(100,*) 'Variables=X,Y,U,V,CP'
          WRITE(100,*) 'ZONE  I=',NX,'J=',NY
          CLOSE(100)
      ELSE
          OPEN(100,FILE = '.\output\Grid.dat',STATUS = 'NEW' )
          WRITE(100,*) 'Variables=X,Y,U,V,CP' 
          WRITE(100,*) 'ZONE','I=',NX,'J=',NY
          CLOSE(100)
      END IF ! IF( ALIVE ) THEN
      
      WRITE(*,*) 'END INPUT'
      
      END ! SUBROUTINE INPUT
      
C====================================================================
C             SUBROUTINE INITIALIZATION -- luzuo
C====================================================================
      SUBROUTINE INITIALIZATION
      USE MAIN
      INTEGER :: I
      
      ALLOCATE( X(NX,NY) , Y(NX,NY) )
      ALLOCATE( PHI(NX,NY) , U(NX,NY) ,V(NX,NY),CP(NX,NY) )
      
      XI_STEP = 1.D0 / ( NX -1 )
      ETA_STEP = 1.D0 / ( NY - 1 )
      
      WRITE(*,*) 'END INITIALIZATION'
      
      END ! SUBROUTINE INITIALIZATION
      
C====================================================================
C             SUBROUTINE OUTPUT -- luzuo
C====================================================================
      SUBROUTINE OUTPUT
      USE MAIN
      
      OPEN(100,FILE='.\output\Grid.dat',STATUS='OLD',POSITION='APPEND')
      DO J = 1 , NY
          DO I = 1 , NX 
              WRITE(100,"(5F16.6)") X(I,J),Y(I,J),U(I,J),V(I,J),CP(I,J)
          END DO 
      END DO
      CLOSE(100)
      
      
      
      END ! SUBROUTINE OUTPUT