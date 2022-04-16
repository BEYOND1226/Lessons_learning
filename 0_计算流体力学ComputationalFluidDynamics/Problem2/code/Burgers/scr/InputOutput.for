C===================================================================
C     INPUT --LUZUO
C===================================================================
      SUBROUTINE INPUT
      USE MAIN
      LOGICAL ALIVE
      
      INQUIRE( FILE = '.\indata\Burgers.in',EXIST = ALIVE )
      IF ( ALIVE ) THEN
          OPEN(10 , FILE = '.\indata\Burgers.in')
          READ(10,*)
          READ(10,*) MU , NTECHNIQUE ,DELTA_X , DELTA_T
          CLOSE(10)
      ELSE
          WRITE(*,*) 'CAN NOT FIND Burgers.in !!!'
      END IF ! IF ( ALIVE ) THEN
      
      IF( NTECHNIQUE == 0 )THEN
          INQUIRE( FILE = '.\outdata\Lax_Wendroff.dat',EXIST = ALIVE )
          IF ( ALIVE ) THEN
              OPEN( 100 ,FILE = '.\outdata\Lax_Wendroff.dat',
     1        STATUS='replace')
          ELSE
              OPEN( 100 ,FILE = '.\outdata\Lax_Wendroff.dat',
     1        STATUS='NEW ')
          END IF ! IF ( ALIVE ) THEN
          WRITE(100,*) 'variables = X , U'
          CLOSE(100)
      ELSE
          INQUIRE( FILE = '.\outdata\MarCormak.dat',EXIST = ALIVE )
          IF ( ALIVE ) THEN
              OPEN( 101 ,FILE = '.\outdata\MarCormak.dat',
     1        STATUS='replace')
          ELSE
              OPEN( 101 ,FILE = '.\outdata\MarCormak.dat',
     1        STATUS='NEW ')
          END IF ! IF ( ALIVE ) THEN
          WRITE(101,*) 'variables = X , U'
          CLOSE(101)
          
c          INQUIRE( FILE = '.\outdata\MarCormack_t0.8.dat',EXIST = ALIVE)
c          IF ( ALIVE ) THEN
c              OPEN( 108 ,FILE = '.\outdata\MarCormack_t0.8.dat',
c     1        STATUS='replace')
c          ELSE
c              OPEN( 108 ,FILE = '.\outdata\MarCormack_t0.8.dat',
c     1        STATUS='NEW ')
c          END IF ! IF ( ALIVE ) THEN
c          WRITE(108,*) 'variables = X , U'
c          CLOSE(108)
c          
c          INQUIRE( FILE = '.\outdata\MarCormack_t1.6.dat',EXIST = ALIVE)
c          IF ( ALIVE ) THEN
c              OPEN( 116 ,FILE = '.\outdata\MarCormack_t1.6.dat',
c     1        STATUS='replace')
c          ELSE
c              OPEN( 116 ,FILE = '.\outdata\MarCormack_t1.6.dat',
c     1        STATUS='NEW ')
c          END IF ! IF ( ALIVE ) THEN
c          WRITE(116,*) 'variables = X , U'
c          CLOSE(116)
c          
c          INQUIRE( FILE = '.\outdata\MarCormack_t2.0.dat',EXIST = ALIVE)
c          IF ( ALIVE ) THEN
c              OPEN( 120 ,FILE = '.\outdata\MarCormack_t2.0.dat',
c     1        STATUS='replace')
c          ELSE
c              OPEN( 120 ,FILE = '.\outdata\MarCormack_t2.0.dat',
c     1        STATUS='NEW ')
c          END IF ! IF ( ALIVE ) THEN
c          WRITE(120,*) 'variables = X , U'
c          CLOSE(120)
          
      END IF ! IF( NTECHNIQUE == 0 )THEN

      END ! SUBROUTINE INPUT
      
C===================================================================
C     OUTPUT --LUZUO
C===================================================================
      SUBROUTINE OUTPUT
      USE MAIN
      
      IF( NTECHNIQUE == 0 )THEN
          OPEN( 100 ,FILE = '.\outdata\Lax_Wendroff.dat',STATUS='OLD'
     1                  ,POSITION = 'APPEND'   )
          WRITE(100,*) 'ZONE  '
          DO I = 1 , NX
              WRITE(100,*) X(I) , U(I)
          END DO ! DO I = 1 , NX
          CLOSE(100)
      ELSE
          OPEN( 101 ,FILE = '.\outdata\MarCormak.dat',STATUS='OLD'
     1                  ,POSITION = 'APPEND'   )
          WRITE(101,*) 'ZONE  '
          DO I = 1 , NX
              WRITE(101,*) X(I) , U(I)
          END DO ! DO I = 1 , NX
          CLOSE(101)
          
c          IF ( ABS(T - 0.8) <= 0.0005  ) THEN
c              OPEN( 108 ,FILE = '.\outdata\MarCormack_t0.8.dat',
c     1                  STATUS='OLD',POSITION = 'APPEND'   )
c              WRITE(108,*) 'ZONE  '
c              DO I = 1 , NX
c                  WRITE(108,*) X(I) , U(I)
c              END DO ! DO I = 1 , NX
c          END IF ! IF ( T == 0.8 ) THEN
c          CLOSE(108)
c          
c          IF ( ABS(T - 1.6) <= 0.0005 ) THEN
c              OPEN( 116 ,FILE = '.\outdata\MarCormack_t1.6.dat',
c     1                  STATUS='OLD',POSITION = 'APPEND'   )
c              WRITE(116,*) 'ZONE  '
c              DO I = 1 , NX
c                  WRITE(108,*) X(I) , U(I)
c              END DO ! DO I = 1 , NX
c          END IF ! IF ( T == 1.6 ) THEN
c          CLOSE(116)
c          
c          IF ( ABS(T - 2.0) <= 0.0005  ) THEN
c              OPEN( 120 ,FILE = '.\outdata\MarCormack_t2.0.dat',
c     1                  STATUS='OLD',POSITION = 'APPEND'   )
c              WRITE(120,*) 'ZONE  '
c              DO I = 1 , NX
c                  WRITE(120,*) X(I) , U(I)
c              END DO ! DO I = 1 , NX
c          END IF ! IF ( T == 2.0 ) THEN
c          CLOSE(120)
          
      END IF ! IF( NTECHNIQUE == 0 )THEN
      END ! SUBROUTINE OUTOUT