










CCCCC        DEFINITIONEN FUER C-PREPROZESSOR
C            MASCHINE

C              MESSAGE PASSING INTERFACE

C            RAEUMLICHE DISKRETISIERUNG
C
***********************************************************************
C                       KOMPAKT-UPWIND in X-RICHTUNG (3-TER ORDNUNG) fuer
C                       die U-Komponente (V und W bleiben Kompakt 4ter ordnung)
C                       falls im Stroemungsfeld eine Koerper vorhanden ist
C
***********************************************************************
C                       KOMPAKTVERF. in XYZ-RICHTUNG (4-TER ORDNUNG)


**********************************************************************
C                      Preprocessing with ADM for 2nd Order Central
***********************************************************************
CCC                     Bei periodischen Randbedingungen in X-Richtung
CCC                     ist eine hoehere O
C
CCC                     Bei periodischen Randbedingungen in Y-Richtung
CCC                     ist eine hoehere Ordnung moeglich
C
CCC                     Bei periodischen Randbedingungen in Y-Richtung:
CCC                     Zentraldiff. 4-ter Ordnung (Parallelisieung moeglich)
C
***********************************************************************
C
C
C            INTERPOLATION AN GITTER-GRENZEN

C            BEHANDLUNG DER TOPAR-RANDBEDINGUNG

C            ZEITLICHE DISKRETISIERUNG

C            FEINSTRUKTURMODELL

C            NICHT-NEWTONSCHE SPANNUNGEN


C            TRANSPORT UND ORIENTIERUNG VON PARTIKELN


C            SCALAR HEAT/TEMPERATURE TRANSPORT
C            APROXIMATE DECONVOLUTION FOR SCALAR
C            TURBULENT PRANDTL NUMBER 
C            PLOT LOCAL SCALAR CONVECTION DIFFUSION EXTREMES 
C            ANALYZE AND PLOT NEAR WALL GRID RESLOLUTION 
C            WRITE SPECIAL 1D LINE FOR TMIX FOR SPECTRA AND PDF
C            SCALAR TIME ADVANCEMENT/DISCRETISATION METHOD

C            STROEMUNG NACH OBEN?


C            BEHANDLUNG DER FLUKTUATIONS-RANDBEDINGUNG

C            BEHANDLUNG VON PPHYS ALS QUELLTERM IN TSTLE2

C            POSITIONIERUNG DER EINSTROEMPROFILE

C           STATISTIK


C                 GEMISCHTES

C                GRDFMI WIRD NICHT VERWENDET, DAHER EINSPARUNG DER FELDER
C                IGRHF UND GRHF

C                VERGROESSERN DER KJI-RICHTUNG BEI FORTSETZUNGSLAUF ERLAUBT!

C                RAUSSCHREIBEN VON ZEITRECORDS

C                FUER KORRELATIONEN UND HAEFIGKEITSVERTEILUNGEN
C$SET SEPARATE                                                           
C$RESET FREE                                                             
      SUBROUTINE FIXRL(DATA,N,NREM,ISIGN,IFORM,NDAT)                    
      DIMENSION DATA(NDAT)                                              
      TWOPI=6.2831853071795865D0*FLOAT(ISIGN)                           
      IP0=2                                                             
      IP1=IP0*(N/2)                                                     
      IP2=IP1*NREM                                                      
      IF (IFORM) 1,7,7                                                  
C PACK THE REAL INPUT VALUES (TWO PER COLUMN)                           
    1     J1=IP1+1                                                          
      DATA(2)=DATA(J1)                                                  
      IF (NREM-1) 7,7,2                                                 
    2     J1=J1+IP0                                                         
      I2MIN=IP1+1                                                       
      DO 6 I2=I2MIN,IP2,IP1                                             
         DATA(I2)=DATA(J1)                                              
         J1=J1+IP0                                                      
         IF (N-2) 5,5,3                                                 
    3        I1MIN=I2+IP0                                                   
         I1MAX=I2+IP1-IP0                                               
         DO 4 I1=I1MIN,I1MAX,IP0                                        
            DATA(I1)=DATA(J1)                                           
            DATA(I1+1)=DATA(J1+1)                                       
    4        J1=J1+IP0                                                      
    5        DATA(I2+1)=DATA(J1)                                            
    6     J1=J1+IP0                                                         
    7     DO 9 I2=1,IP2,IP1                                                 
         TEMPR=DATA(I2)                                                 
         DATA(I2)=DATA(I2)+DATA(I2+1)                                   
         DATA(I2+1)=TEMPR-DATA(I2+1)                                    
         IF (IFORM) 8,9,9                                               
    8        DATA(I2)=DATA(I2)/2.                                           
         DATA(I2+1)=DATA(I2+1)/2.                                       
    9     CONTINUE                                                          
      IF (N-2) 18,18,10                                                 
   10    THETA=TWOPI/FLOAT(N)                                              
      SINTH=SIN(THETA/2.)                                               
      ZSTPR=-2.*SINTH*SINTH                                             
      ZSTPI=SIN(THETA)                                                  
      ZR=(1.-ZSTPI)/2.                                                  
      ZI=(1.+ZSTPR)/2.                                                  
      IF (IFORM) 11,12,12                                               
   11    ZR=1.-ZR                                                          
      ZI=-ZI                                                            
   12    I1MIN=IP0+1                                                       
      I1MAX=IP0*(N/4)+1                                                 
      DO 17 I1=I1MIN,I1MAX,IP0                                          
         DO 16 I2=I1,IP2,IP1                                            
            I2CNJ=N+IP0-2*I1+I2                                         
            IF (I2-I2CNJ) 15,13,13                                      
   13          IF (ISIGN*(2*IFORM+1)) 14,16,16                             
   14          DATA(I2+1)=-DATA(I2+1)                                      
         GO TO 16                                                       
   15          DIFR=DATA(I2)-DATA(I2CNJ)                                   
            DIFI=DATA(I2+1)+DATA(I2CNJ+1)                               
            TEMPR=DIFR*ZR-DIFI*ZI                                       
            TEMPI=DIFR*ZI+DIFI*ZR                                       
            DATA(I2)=DATA(I2)-TEMPR                                     
            DATA(I2+1)=DATA(I2+1)-TEMPI                                 
            DATA(I2CNJ)=DATA(I2CNJ)+TEMPR                               
            DATA(I2CNJ+1)=DATA(I2CNJ+1)-TEMPI                           
   16    CONTINUE                                                          
         TEMPR=ZR-.5                                                    
         ZR=ZSTPR*TEMPR-ZSTPI*ZI+ZR                                     
   17    ZI=ZSTPR*ZI+ZSTPI*TEMPR+ZI                                        
   18    IF (IFORM) 25,19,19                                               
   19    I2=IP2+1                                                          
      I1=I2                                                             
      J1=IP0*(N/2+1)*NREM+1                                             
      GO TO 23                                                          
   20    DATA(J1)=DATA(I1)                                                 
      DATA(J1+1)=DATA(I1+1)                                             
      I1=I1-IP0                                                         
      J1=J1-IP0                                                         
   21    IF (I2-I1) 20,22,22                                               
   22    DATA(J1)=DATA(I1)                                                 
      DATA(J1+1)=0.                                                     
   23    I2=I2-IP1                                                         
      J1=J1-IP0                                                         
      DATA(J1)=DATA(I2+1)                                               
      DATA(J1+1)=0.                                                     
      I1=I1-IP0                                                         
      J1=J1-IP0                                                         
      IF (I2-1) 24,24,21                                                
   24    DATA(2)=0.                                                        
   25    RETURN                                                            
      END                                                               
C$SET SEPARATE                                                           
C$RESET FREE                                                             
      SUBROUTINE PREPO(BR,X,K1,K2,ISIG,IFORM,N1MAX)                     
       DIMENSION BR(1),X(1)                                             
      K12=2*N1MAX                                                       
      K21=K2+1                                                          
      IA=K1*K21                                                         
      IF (IFORM.EQ.0) KX=K1                                             
      IF (IFORM.EQ.-1) KX=K1+2                                          
      IS=1-K12                                                          
      DO 1 J=1,K21                                                      
      IS=IS+K12                                                         
      CALL SELFFT(X,KX,IS,BR,IA)                                        
      CALL FIXRL(X,K1,1,ISIG,IFORM,KX)                                  
      CALL BACSET(X,K1,IS,BR,IA)                                        
    1     CONTINUE                                                          
      RETURN                                                            
      END                                                               
C$SET SEPARATE                                                           
C$RESET FREE                                                             
      SUBROUTINE BACSET(AS,M2,IS,A,IA)                                  
      DIMENSION AS(M2),A(IA)                                            
C VECTORS AS ARE COLUMNS OF ARRAY A                                     
      IS1=IS-1                                                          
      DO 1 I=1,M2                                                       
    1     A(IS1+I)=AS(I)                                                    
      RETURN                                                            
      END                                                               
C     C$SET SEPARATE                                                         
C     C$RESET FREE                                                           
      SUBROUTINE FFTGUT(AR,AI,NCLBT,NSTRN,EX,WORK,LWORK)                
      DIMENSION AR(NCLBT*NSTRN*NPTS*NSEP),
     1          AI(NCLBT*NSTRN*NPTS*NSEP),
     2          EX(2*NPTS),
     3          WORK(LWORK)
      COMMON /FFTCOM/ NPTS,NSEP,SIGNEX,NA1,NA2,NA3,ISCRAM,LWRK          
      DATA IZR/0/                                                       
      DATA E1D3R/-.5/,                                                  
     1     E1D3I/ .866025403784439/                                     
      DATA E1D5R/ .309016994374947/,                                    
     1     E1D5I/ .951056516295154/,                                    
     2     E2D5R/-.809016994374947/,                                    
     3     E2D5I/ .587785252292473/                                     
      MA1 = NSEP*NA1                                                    
      MA1NA3 = MA1*NA3                                                  
      MA1NA2 = MA1*NA2                                                  
      IA3MAX = (NA3-1)*MA1NA2                                           
      IA2MAX = (NA2-1)*MA1                                              
      IF (ISCRAM.EQ.+1) GO TO 150                                       
      IF (NA1.LT.2) GO TO 150                                           
      INCINC = 2*NA3                                                    
      IXINC = 0                                                         
      IA1MIN = NSEP+1                                                   
      DO 100 IA1=IA1MIN,MA1,NSEP                                        
      IXINC = IXINC+INCINC                                              
      MIN = IA1+MA1                                                     
      MAX = IA1+IA2MAX                                                  
      IX = 1                                                            
      DO 100 IAMIN=MIN,MAX,MA1                                          
      IX = IX+IXINC                                                     
      EXR = EX(IX)                                                      
      EXI = SIGNEX*EX(IX+1)                                             
      IAMAX = IAMIN+IA3MAX                                              
      DO 100 IA=IAMIN,IAMAX,MA1NA2                                      
      DO 100 IT=1,NSTRN                                                 
      IB=IA+(IT-1)*NCLBT                                                
      ERAR=EXR*AR(IB)                                                   
      EIAR=EXI*AR(IB)                                                   
      ERAI = EXR*AI(IB)                                                 
      EIAI = EXI*AI(IB)                                                 
      AR(IB)= ERAR-EIAI                                                 
      AI(IB)= EIAR+ERAI                                                 
  100   CONTINUE                                                          
  150   IGO = MIN0(NA2,6)                                                 
      IGO = IGO-1                                                       
      GO TO (200,300,400,500,700),IGO                                   
  200   CONTINUE                                                          
C  TRANSFORM OF LENGTH 2                                                
      DO 250 IA3=IZR,IA3MAX,MA1NA2                                      
      DO 250 J1=1,MA1,NSEP                                              
!OCL NOVREC(AR,AI)
      DO 250 IT=1,NSTRN                                                 
      IIA0 = IA3+J1                                                     
      IIA0=IIA0+(IT-1)*NCLBT                                            
      IIA1 = IIA0+MA1                                                   
      S0R = AR(IIA0)+AR(IIA1)                                           
      D0R = AR(IIA0)-AR(IIA1)                                           
      AR(IIA0) = S0R                                                    
      AR(IIA1) = D0R                                                    
      S0I = AI(IIA0)+AI(IIA1)                                           
      D0I = AI(IIA0)-AI(IIA1)                                           
      AI(IIA0) = S0I                                                    
      AI(IIA1) = D0I                                                    
  250   CONTINUE                                                          
      GO TO 900                                                         
  300   CONTINUE                                                          
C  TRANSFORM OF LENGTH 3                                                
      E1D3I = SIGN(E1D3I,SIGNEX)                                        
      DO 350 IA3=IZR,IA3MAX,MA1NA2                                      
      DO 350 J1=1,MA1,NSEP                                              
!OCL NOVREC(AR,AI)
      DO 350 IT=1,NSTRN                                                 
      IIA0 = IA3+J1                                                     
      IIA0=IIA0+(IT-1)*NCLBT                                            
      IIA1 = IIA0+MA1                                                   
      IIA2 = IIA1+MA1                                                   
      S0R = AR(IIA0)                                                    
      S0I = AI(IIA0)                                                    
      S1R = AR(IIA1)+AR(IIA2)                                           
      D1R = AR(IIA1)-AR(IIA2)                                           
      S1I = AI(IIA1)+AI(IIA2)                                           
      D1I = AI(IIA1)-AI(IIA2)                                           
      AR(IIA0) = S0R+S1R                                                
      AI(IIA0) = S0I+S1I                                                
      BS1R = S0R+E1D3R*S1R                                              
      BS1I = S0I+E1D3R*S1I                                              
      BD1I = E1D3I*D1R                                                  
      BD1R = -E1D3I*D1I                                                 
      AR(IIA1) = BS1R+BD1R                                              
      AR(IIA2) = BS1R-BD1R                                              
      AI(IIA1) = BS1I+BD1I                                              
      AI(IIA2) = BS1I-BD1I                                              
  350   CONTINUE                                                          
      GO TO 900                                                         
  400   CONTINUE                                                          
C  TRANSFORM OF LENGTH 4                                                
      DO 450 IA3=IZR,IA3MAX,MA1NA2                                      
      DO 450 J1=1,MA1,NSEP                                              
!OCL NOVREC(AR,AI)
      DO 450 IT=1,NSTRN                                                 
      IIA0 = IA3+J1                                                     
      IIA0=IIA0+(IT-1)*NCLBT                                            
      IIA1 = IIA0+MA1                                                   
      IIA2 = IIA1+MA1                                                   
      IIA3 = IIA2+MA1                                                   
      S0R = AR(IIA0)+AR(IIA2)                                           
      D0R = AR(IIA0)-AR(IIA2)                                           
      S0I = AI(IIA0)+AI(IIA2)                                           
      D0I = AI(IIA0)-AI(IIA2)                                           
      S1R = AR(IIA1)+AR(IIA3)                                           
      D1R = SIGNEX*(AR(IIA1)-AR(IIA3))                                  
      S1I = AI(IIA1)+AI(IIA3)                                           
      D1I = SIGNEX*(AI(IIA1)-AI(IIA3))                                  
      AR(IIA0) = S0R+S1R                                                
      AR(IIA2) = S0R-S1R                                                
      AR(IIA1) = D0R-D1I                                                
      AR(IIA3) = D0R+D1I                                                
      AI(IIA0) = S0I+S1I                                                
      AI(IIA2) = S0I-S1I                                                
      AI(IIA1) = D0I+D1R                                                
      AI(IIA3) = D0I-D1R                                                
  450   CONTINUE                                                          
      GO TO 900                                                         
  500   CONTINUE                                                          
C  TRANSFORM OF LENGTH 5                                                
      E1D5I = SIGN(E1D5I,SIGNEX)                                        
      E2D5I = SIGN(E2D5I,SIGNEX)                                        
      DO 550 IA3=IZR,IA3MAX,MA1NA2                                      
      DO 550 J1=1,MA1,NSEP                                              
!OCL NOVREC(AR,AI)
      DO 550 IT=1,NSTRN                                                 
      IIA0 = IA3+J1                                                     
      IIA0=IIA0+(IT-1)*NCLBT                                            
      IIA1 = IIA0+MA1                                                   
      IIA2 = IIA1+MA1                                                   
      IIA3 = IIA2+MA1                                                   
      IIA4 = IIA3+MA1                                                   
      S0R = AR(IIA0)                                                    
      S0I = AI(IIA0)                                                    
      S1R = AR(IIA1)+AR(IIA4)                                           
      D1R = AR(IIA1)-AR(IIA4)                                           
      S1I = AI(IIA1)+AI(IIA4)                                           
      D1I = AI(IIA1)-AI(IIA4)                                           
      S2R = AR(IIA2)+AR(IIA3)                                           
      D2R = AR(IIA2)-AR(IIA3)                                           
      S2I = AI(IIA2)+AI(IIA3)                                           
      D2I = AI(IIA2)-AI(IIA3)                                           
      AR(IIA0) = S0R+S1R+S2R                                            
      AI(IIA0) = S0I+S1I+S2I                                            
      BS1R = S0R+E1D5R*S1R+E2D5R*S2R                                    
      BS1I = S0I+E1D5R*S1I+E2D5R*S2I                                    
      BD1I =     E1D5I*D1R+E2D5I*D2R                                    
      BD1R =    -E1D5I*D1I-E2D5I*D2I                                    
      AR(IIA1) = BS1R+BD1R                                              
      AR(IIA4) = BS1R-BD1R                                              
      AI(IIA1) = BS1I+BD1I                                              
      AI(IIA4) = BS1I-BD1I                                              
      BS2R = S0R+E2D5R*S1R+E1D5R*S2R                                    
      BS2I = S0I+E2D5R*S1I+E1D5R*S2I                                    
      BD2I =     E2D5I*D1R-E1D5I*D2R                                    
      BD2R =    -E2D5I*D1I+E1D5I*D2I                                    
      AR(IIA2) = BS2R+BD2R                                              
      AR(IIA3) = BS2R-BD2R                                              
      AI(IIA2) = BS2I+BD2I                                              
      AI(IIA3) = BS2I-BD2I                                              
  550   CONTINUE                                                          
      GO TO 900                                                         
  700   CONTINUE                                                          
C  TRANSFORM OF ANY ODD LENGTH                                          
      IF (2*NA2.GT.LWRK) GO TO 9000                                     
      IXINC2 = 2*NA1*NA3                                                
      N2 = 2*NPTS                                                       
      IA2HMX = MA1*NA2/2                                                
      IWHMX = NA2                                                       
      IWCM = 2*NA2+2                                                    
      DO 850 IA3=IZR,IA3MAX,MA1NA2                                      
      DO 850 IA1=1,MA1,NSEP                                             
      DO 850 IT=1,NSTRN                                                 
      IIA0 = IA3+IA1                                                    
      IIA0=IIA0+(IT-1)*NCLBT                                            
      BSR = AR(IIA0)                                                    
      WORK(1) = BSR                                                     
      BSI = AI(IIA0)                                                    
      WORK(2) = BSI                                                     
      IACM = MA1NA2+2*IIA0                                              
      IAMIN = IIA0+MA1                                                  
      IAMAX = IIA0+IA2MAX                                               
      IW = 1                                                            
      DO 750 IA=IAMIN,IAMAX,MA1                                         
      IW = IW+2                                                         
      WORK(IW) = AR(IA)                                                 
      BSR = BSR+AR(IA)                                                  
      WORK(IW+1) = AI(IA)                                               
      BSI = BSI+AI(IA)                                                  
  750   CONTINUE                                                          
      AR(IIA0) = BSR                                                    
      AI(IIA0) = BSI                                                    
      IAMAX = IIA0+IA2HMX                                               
      IXINC = 0                                                         
      DO 850 IA=IAMIN,IAMAX,MA1                                         
      IXINC = IXINC+IXINC2                                              
      BSR = WORK(1)                                                     
      BSI = WORK(2)                                                     
      BDR = 0.                                                          
      BDI = 0.                                                          
      IX = 1                                                            
      DO 800 IW=3,IWHMX,2                                               
      IX = IX+IXINC                                                     
      IF (IX.GT.N2) IX=IX-N2                                            
      EXR = EX(IX)                                                      
      EXI = SIGNEX*EX(IX+1)                                             
      IWC = IWCM-IW                                                     
      BSR = BSR+EXR*(WORK(IW  )+WORK(IWC  ))                            
      BDI = BDI+EXI*(WORK(IW  )-WORK(IWC  ))                            
      BDR = BDR-EXI*(WORK(IW+1)-WORK(IWC+1))                            
      BSI = BSI+EXR*(WORK(IW+1)+WORK(IWC+1))                            
  800   CONTINUE                                                          
      IAC = IACM-IA                                                     
      AR(IA ) = BSR+BDR                                                 
      AR(IAC) = BSR-BDR                                                 
      AI(IAC) = BSI-BDI                                                 
      AI(IA ) = BSI+BDI                                                 
  850   CONTINUE                                                          
  900   IF (ISCRAM.NE.+1) RETURN                                          
      IF (NA1.LT.2) RETURN                                              
      INCINC = 2*NA3                                                    
      IXINC = 0                                                         
      IA1MIN = NSEP+1                                                   
      DO 950 IA1=IA1MIN,MA1,NSEP                                        
      IXINC = IXINC+INCINC                                              
      MIN = IA1+MA1                                                     
      MAX = IA1+IA2MAX                                                  
      IX = 1                                                            
      DO 950 IAMIN=MIN,MAX,MA1                                          
      IX = IX+IXINC                                                     
      EXR = EX(IX)                                                      
      EXI = SIGNEX*EX(IX+1)                                             
      IAMAX = IAMIN+IA3MAX                                              
      DO 950 IA=IAMIN,IAMAX,MA1NA2                                      
!OCL NOVREC(AR,AI)
      DO 950 IT=1,NSTRN                                                 
      IB=IA+(IT-1)*NCLBT                                                
      ERAR = EXR*AR(IB)                                                 
      EIAR = EXI*AR(IB)                                                 
      ERAI = EXR*AI(IB)                                                 
      EIAI = EXI*AI(IB)                                                 
      AR(IB)=ERAR-EIAI                                                  
      AI(IB)=EIAR+ERAI                                                  
  950   CONTINUE                                                          
      RETURN                                                            
 9000  PRINT 6900,NA2,LWRK                                               
 6900  FORMAT (' INSUFFICIENT WORKSPACE PRIME FACTOR AND LWORK =',2I10)   
      STOP                                                              
      END                                                               
C      C$SET SEPARATE                                                        
C      C$RESET FREE                                                          
C                                                                       
      SUBROUTINE FFTSET(N,NFACT,IPERMU,EXPON)                           
C  FACTORS N INTO PRIMES AND SETS UP ARRAYS NFACT, EXPON AND IPERMU FOR 
C         USE BY FFT ROUTIMES.                                          
C                                                                       
C  ON INPUT                                                             
C                                                                       
C  N          LENGTH OF TRANSFORM TO BE PERFORMED.                      
C                                                                       
C  ON OUTPUT                                                            
C                                                                       
C  NFACT      NFACT(1) IS M, THE NUMBER OF PRIME FACTORS OF N.          
C             NFACT(2),...,NFACT(M+1) ARE THOSE PRIME FACTORS.  FOR     
C         PURPOSES OF THE FFT IT IS USEFUL TO CONSIDER 4 A PRIME, AND   
C         ALL FACTORS OF 4 APPEAR 1ST IN THE LIST.                      
C    CAUTION - NFACT MUST BE DIMENSIONED SUFFICIENTLY IN THE CALLING    
C         PROGRAM.  A DIMENSION OF 16 WILL SUFFICE FOR ANY N UP TO ABOUT
C         9 MILLION.                                                    
C  IPERMU     AN ARRAY CONTAINING THE PERMUTATION INDICES FOR USE       
C         IN REORDERING TRANSFORMED ARRAYS.                             
C    CAUTION - IPERMU MUST BE DIMENSIONED SUFFICIENTLY IN THE CALLING   
C         PROGRAM.  A DIMENSION OF N WILL SUFFICE.                      
C  EXPON      AN ARRAY CONTAINING CEXP(2*PI*I*(J-1)/N) FOR J=1,N STORED 
C         IN THE USUAL COMPLEX CONVENTION.                              
C    CAUTION - EXPON MUST BE DIMENSIONED SUFFICIENTLY IN THE CALLING    
C         PROGRAM.  A DIMENSION OF 2*N WILL SUFFICE.                    
C                                                                       
C                                                                       
C  NOTE - THE DIMENSIONS OF THE FOLLOWING LOCAL VARIABLES LIMIT THE     
C         NUMBER OF FACTORS N MAY HAVE.  16 IS SUFFICIENT FOR ANY N     
C         UP TO ABOUT 9 MILLION.                                        
      DIMENSION IDIG(16),NPROD(16)   ,NFACT(16),EXPON(2*N),IPERMU(2*N)  
      DATA MMAX/16/                                                     
      DATA PI/3.141592653589793/                                        
      IF (N-1) 900,10,20                                                
   10    NFACT(1) = 1                                                      
      NFACT(2) = 1                                                      
      EXPON(1) = 1.                                                     
      EXPON(2) = 0.                                                     
      GO TO 400                                                         
   20    N2 = 2*N                                                          
      M = 0                                                             
      NL = N                                                            
      NTRY = 4                                                          
   50    NQ = NL/NTRY                                                      
      NR = NL-NTRY*NQ                                                   
      IF (NR.EQ.0) GO TO 80                                             
      IF (NTRY.GT.4) GO TO 70                                           
      NGO = NTRY-1                                                      
      GO TO (62,63,64),NGO                                              
   62    NTRY = 3                                                          
      GO TO 50                                                          
   63    NTRY = 5                                                          
      GO TO 50                                                          
   64    NTRY = 2                                                          
      GO TO 50                                                          
   70    NTRY = NTRY+2                                                     
      IF (NTRY.LE.NQ) GO TO 50                                          
      NTRY = NL                                                         
      NQ = 1                                                            
   80    M = M+1                                                           
      NFACT(M+1) = NTRY                                                 
      NL = NQ                                                           
      IF (NL.NE.1) GO TO 50                                             
      NFACT(1) = M                                                      
      ARG = 2.*PI/N                                                     
      DO 500 J=1,N                                                      
      ARGJ = (J-1)*ARG                                                  
      EXPON(2*J-1) = COS(ARGJ)                                          
      EXPON(2*J  ) = SIN(ARGJ)                                          
  500   CONTINUE                                                          
      IF (M.EQ.1) GO TO 400                                             
      IF (M.GT.MMAX) GO TO 910                                          
      MP1 = M+1                                                         
      MP2 = M+2                                                         
      IP = 1                                                            
      DO 150 I=1,M                                                      
      IP = IP*NFACT(MP2-I)                                              
      NPROD(MP1-I) = IP                                                 
      IDIG(I) = 0                                                       
  150   CONTINUE                                                          
      NM1 = N-1                                                         
      MM1 = M-1                                                         
      JREV = 0                                                          
      IPERMU(1) = 1                                                     
      DO 300 J=2,NM1                                                    
C  INCREMENT DIGIT REVERSED COUNTER                                     
      DO 200 I=1,MM1                                                    
      JREV = JREV+NPROD(I+1)                                            
      IDIG(I) = IDIG(I)+1                                               
      IF (IDIG(I).LT.NFACT(I+1)) GO TO 250                              
      IDIG(I) = 0                                                       
      JREV = JREV-NPROD(I)                                              
  200   CONTINUE                                                          
      JREV = JREV+1                                                     
  250   IPERMU(J) = JREV+1                                                
  300   CONTINUE                                                          
      IPERMU(N) = N                                                     
      RETURN                                                            
  400   DO 450 I=1,N                                                      
      IPERMU(I) = I                                                     
  450   CONTINUE                                                          
      RETURN                                                            
  900   PRINT 6900,N                                                      
 6900  FORMAT ("0N = ",I20," IS 0 OR NEGATIVE")                            
      STOP                                                              
  910   PRINT 6910                                                        
 6910  FORMAT (" INSUFFICIENT LOCAL STORAGE")                            
      RETURN                                                            
      END                                                               
C      C$SET SEPARATE                                                        
C      C$RESET FREE                                                          
      SUBROUTINE FTCCSC(N,AR,AI,NCLBP,NCLBT,NSTRN,SGNEX,NFACT,EXPON,    
     $                  WORK,LWORK)                                     
C. FOURIER TRANSFORM - COMPLEX TO COMPLEX SCRAMBLED                     
C                                                                       
C  ON INPUT                                                             
C                                                                       
C  N      LENGTH OF TRANSFORM                                           
C  AR     REAL PART OF DATA TO BE TRANSFORMED                           
C  AI     IMAGINARY PART OF DATA TO BE TRANSFORMED                      
C  NCLBP  NUMBER OF CORE LOCATIONS BETWEEN POINTS TO BE TRANSFORMED.    
C         I.E. AR AND AI MAY BE CONSIDERED TO HAVE DIMENSION (NCLBP,N)  
C         AND THE TRANSFORM IS PERFORMED ALONG THE 2ND DIMENSION.       
C         THIS CONFIGURATION IS USED BOTH FOR INPUT AND OUTPUT.         
C  NCLBT  NUMBER OF CORE LOCATIONS BETWEEN TRANSFORMATIONS.             
C  NSTRN  NUMBER OF SIMULTANEOUS TRANSFORMATIONS.                       
C  SGNEX  SIGN OF THE ARGUMENT OF THE COMPLEX EXPONENTIAL USED IN THE   
C         TRANSFORM                                                     
C  NFACT  AN ARRAY CONTAINING A PRIME FACTORIZATION OF N.  IT MUST BE   
C         SET UP BY CALLING FFTSET.                                     
C  EXPON  AN ARRAY OF SINES AND COSINES.  IT MUST BE SET UP BY CALLING  
C         FFTSET.                                                       
C  WORK   A WORKSPACE FOR USE BY THE ROUTINE.  IT IS NOT USED IF N HAS  
C         NO PRIME FACTORS LARGER THAN 5.                               
C  LWORK  THE LENGTH OF THE WORKSPACE.  IT MUST BE .GE. THE LARGEST     
C         PRIME FACTOR OF N UNLESS N HAS NO PRIME FACTORS LARGER THAN 5.
C         IN THAT CASE IT MAY BE 0                                      
C                                                                       
C  ON OUTPUT                                                            
C                                                                       
C  AR     REAL PART OF TRANSFORM RESULT IN SCRAMBLED ORDER              
C  AI     IMAGINARY PART OF TRANSFORM RESULT IN SCRAMBLED ORDER         
      DIMENSION AR(NCLBP,N),AI(NCLBP,N),
     1          NFACT(16),EXPON(2,N),WORK(LWORK)
      COMMON /FFTCOM/ NPTS,NSEP,SIGNEX,NA1,NA2,NA3,ISCRAM,LWRK          
      ISCRAM = +1                                                       
      NPTS = N                                                          
      IF (NPTS.EQ.1) RETURN                                             
      NSEP = NCLBP                                                      
      SIGNEX = SIGN(1.,SGNEX)                                           
      LWRK = LWORK                                                      
      M = NFACT(1)                                                      
      NA1 = N                                                           
      NA3 = 1                                                           
      DO 200 I=1,M                                                      
      NA2 = NFACT(M+2-I)                                                
      NA1 = NA1/NA2                                                     
      CALL FFTGUT(AR,AI,NCLBT,NSTRN,EXPON,WORK,LWORK)                   
      NA3 = NA2*NA3                                                     
  200   CONTINUE                                                          
      RETURN                                                            
      END                                                               
C      C$SET SEPARATE                                                        
C      C$RESET FREE                                                          
      SUBROUTINE FTSCCC(N,AR,AI,NCLBP,NCLBT,NSTRN,SGNEX,NFACT,EXPON,    
     $                  WORK,LWORK)                                     
C. FOURIER TRANSFORM - SCRAMBLED COMPLEX TO COMPLEX                     
C  THE ARGUMENTS FOR THIS ROUTINE ARE IDENTICAL TO THOSE FOR FTCCSC     
C         EXCEPT THAT THE INPUT IS ASSUMED SCRAMBLED AND THE OUTPUT IS  
C         IN THE CORRECT ORDER.                                         
      DIMENSION AR(NCLBP,N),AI(NCLBP,N),
     1          NFACT(16),EXPON(2,N),WORK(LWORK) 
      COMMON /FFTCOM/ NPTS,NSEP,SIGNEX,NA1,NA2,NA3,ISCRAM,LWRK          
      ISCRAM = -1                                                       
      NPTS = N                                                          
      IF (NPTS.EQ.1) RETURN                                             
      NSEP = NCLBP                                                      
      SIGNEX = SIGN(1.,SGNEX)                                           
      LWRK = LWORK                                                      
      M = NFACT(1)                                                      
      NA1 = 1                                                           
      NA3 = N                                                           
      DO 200 I=1,M                                                      
      NA2 = NFACT(I+1)                                                  
      NA3 = NA3/NA2                                                     
      CALL FFTGUT(AR,AI,NCLBT,NSTRN,EXPON,WORK,LWORK)                   
      NA1 = NA2*NA1                                                     
  200   CONTINUE                                                          
      RETURN                                                            
      END                                                               
C$SET SEPARATE                                                           
C$RESET FREE                                                             
      SUBROUTINE SELFFT(AS,M2,IS,A,IA)                                  
      DIMENSION AS(M2),A(IA)                                            
C ONE ROW AFTER THE OTHER OF ARRAY A IS ASSIGNED TO THE VECTOR S        
      IS1=IS-1                                                          
      DO 1 I=1,M2                                                       
    1     AS(I)=A(IS1+I)                                                    
      RETURN                                                            
      END                                                               
C      C$SET SEPARATE                                                        
C      C$RESET FREE                                                          
      SUBROUTINE SCRAMB(NTYP,N,AR,AI,NCLBP,NCLBT,NSTRN,IPERMU,WORK)     
C. UNSCRAMBLES OR SCRAMBLES (DEPENDING ON WHETHER NTYP IS -1 OR +1) THE 
C.        ARRAYS AR AND AI.                                             
C  THE SCRAMBLING IS THE DIGIT REVERSAL REORDERING WHICH OCCURS WITH    
C         THE USE OF THE FAST FOURIER TRANSFORM ROUTINE FTCCSC.         
C  SIMILAR TO THE OTHER ROUTINES IN THIS PACKAGE, A NUMBER OF           
C         SIMULTANEOUS REORDERINGS MAY BE ACCOMPLISHED IN A SINGLE CALL.
C                                                                       
C  ON INPUT                                                             
C                                                                       
C  NTYP       A FLAG TO INDICATE THE TYPE OF REORDERING.                
C         NTYP = -1 MEANS THE ARRAYS ARE TO BE UNSCRAMBLED.             
C         NTYP = +1 MEANS THE ARRAYS ARE TO BE SCRAMBLED.               
C  N          THE NUMBER OF POINTS TO BE REORDERED.                     
C  AR         REAL PART OF TRANSFORM RESULT TO BE REORDERED.            
C  AI         IMAGINARY PART OF TRANSFORM RESULT TO BE REORDERED.       
C  NCLBP      THE NUMBER OF CORE LOCATIONS BETWEEN POINTS TO BE         
C         REORDERED.  I.E. A MAY BE CONSIDERED TO HAVE DIMENSION        
C         (NCLBP.N) AND THE REORDERING IS TO PERFORMED ALONG THE 2ND    
C         DIMENSION.                                                    
C         THIS CONFIGURATION IS USED BOTH FOR INPUT AND OUTPUT.         
C  NCLBT      THE NUMBER OF CORE LOCATIONS BETWEEN TRANSFORMATIONS.     
C  NSTRN      THE NUMBER OF SIMULTANEOUS TRANSFORMATIONS.               
C  IPERMU     AN ARRAY CONTAINING THE PERMUTATION INDEX FOR REORDERING. 
C         IT MUST BE SET UP BY CALLING FFTSET.                          
C  WORK       A WORKSPACE OF LENGTH 2*N FOR USE BY THE ROUTINE          
C                                                                       
C  ON OUTPUT                                                            
C                                                                       
C  AI         IMAGINARY PART OF THE REORDERED ARRAY.                    
C  AR         REAL PART OF REORDERED ARRAY.                             
      DIMENSION AR(N*NSTRN*2),AI(N*NSTRN*2),IPERMU(N),WORK(2,N)
C   CHECKING ARRAY BOUNDS (M.M. 170998)
      IF (NCLBT .GT. N*2) CALL ERRR(501,'SCRAMB')
C  NOTE - THE FIRST AND LAST ELEMENTS ARE NEVER SCRAMBLED               
      IF (N.LE.2) RETURN                                                
      NM2 = N-2                                                         
      NSEP = NCLBP                                                      
!OCL NOVREC(AR,AI)
      DO 500 IT=1,NSTRN                                                 
      IBASE=(IT-1)*NCLBT+1                                              
      JMIN = NSEP+IBASE                                                 
      JMAX = NM2*NSEP+IBASE                                             
      IF (NTYP .EQ. -1) GO TO 250                                       
      JW = 1                                                            
      DO 100 J=JMIN,JMAX,NSEP                                           
      JW = JW+1                                                         
      WORK(1,JW) = AR(J)                                                    
      WORK(2,JW) = AI(J)                                                    
  100 CONTINUE                                                          
      JW = 1                                                            
      DO 200 J=JMIN,JMAX,NSEP                                           
      JW = JW+1                                                         
      JREV = IPERMU(JW)                                                 
      AR(J) = WORK(1,JREV)                                                  
      AI(J) = WORK(2,JREV)                                                  
  200 CONTINUE                                                          
      GO TO 500                                                         
  250 CONTINUE                                                          
      JW = 1                                                            
      DO 300 J=JMIN,JMAX,NSEP                                           
      JW = JW+1                                                         
      JREV = IPERMU(JW)                                                 
      WORK(1,JREV) = AR(J)                                                  
      WORK(2,JREV) = AI(J)                                                  
  300 CONTINUE                                                          
      JW = 1                                                            
      DO 400 J=JMIN,JMAX,NSEP                                           
      JW = JW+1                                                         
      AR(J) = WORK(1,JW)                                                    
      AI(J) = WORK(2,JW)                                                    
  400 CONTINUE                                                          
  500 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE SCRAMB_S(NTYP,N,AR,AI,NCLBP,NCLBT,NSTRN,IPERMU,WORK,
     $     H2D)     
C. UNSCRAMBLES OR SCRAMBLES (DEPENDING ON WHETHER NTYP IS -1 OR +1) THE 
C.        ARRAYS AR AND AI.                                             
C  THE SCRAMBLING IS THE DIGIT REVERSAL REORDERING WHICH OCCURS WITH    
C         THE USE OF THE FAST FOURIER TRANSFORM ROUTINE FTCCSC.         
C  SIMILAR TO THE OTHER ROUTINES IN THIS PACKAGE, A NUMBER OF           
C         SIMULTANEOUS REORDERINGS MAY BE ACCOMPLISHED IN A SINGLE CALL.
C                                                                       
C  ON INPUT                                                             
C                                                                       
C  NTYP       A FLAG TO INDICATE THE TYPE OF REORDERING.                
C         NTYP = -1 MEANS THE ARRAYS ARE TO BE UNSCRAMBLED.             
C         NTYP = +1 MEANS THE ARRAYS ARE TO BE SCRAMBLED.               
C  N          THE NUMBER OF POINTS TO BE REORDERED.                     
C  AR         REAL PART OF TRANSFORM RESULT TO BE REORDERED.            
C  AI         IMAGINARY PART OF TRANSFORM RESULT TO BE REORDERED.       
C  NCLBP      THE NUMBER OF CORE LOCATIONS BETWEEN POINTS TO BE         
C         REORDERED.  I.E. A MAY BE CONSIDERED TO HAVE DIMENSION        
C         (NCLBP.N) AND THE REORDERING IS TO PERFORMED ALONG THE 2ND    
C         DIMENSION.                                                    
C         THIS CONFIGURATION IS USED BOTH FOR INPUT AND OUTPUT.         
C  NCLBT      THE NUMBER OF CORE LOCATIONS BETWEEN TRANSFORMATIONS.     
C  NSTRN      THE NUMBER OF SIMULTANEOUS TRANSFORMATIONS.               
C  IPERMU     AN ARRAY CONTAINING THE PERMUTATION INDEX FOR REORDERING. 
C         IT MUST BE SET UP BY CALLING FFTSET.                          
C  WORK       A WORKSPACE OF SAME SIZE AS AR FOR USE BY THE ROUTINE          
C                                                                       
C  ON OUTPUT                                                            
C                                                                       
C  AI         IMAGINARY PART OF THE REORDERED ARRAY.                    
C  AR         REAL PART OF REORDERED ARRAY.                             
C----------------------------------------------------------------
C    VERSION WITH INVERTED LOOP-ORDERING OF SCRAMB
C    CAN IMPROVE PERFORMANCE, DEPENDING ON NCLBP
C      M.M. 23.10.00
C----------------------------------------------------------------
      DIMENSION AR(N*NSTRN*2),AI(N*NSTRN*2),IPERMU(N),WORK(2,N)
     1     ,H2D(NSTRN,N)
C   CHECKING ARRAY BOUNDS (M.M. 170998)
      IF (NCLBT .GT. N*2) CALL ERRR(501,'SCRAMB')
C  NOTE - THE FIRST AND LAST ELEMENTS ARE NEVER SCRAMBLED               
      IF (N.LE.2) RETURN                                                
      NM2 = N-2                                                         
      NSEP = NCLBP
C----------------------------------------
      IF (NTYP .NE. -1) THEN
C---------------------------------------- Case I
         DO JT=1,NM2
            DO IT=1,NSTRN                                                 
               IBASE=(IT-1)*NCLBT+1
               J  = NSEP+IBASE + NSEP*(JT-1)
               JW = JT+1
               H2D(IT,JW) = AR(J)
            ENDDO
         ENDDO
         DO JT=1,NM2
            DO IT=1,NSTRN                                                 
               IBASE=(IT-1)*NCLBT+1
               J  = NSEP+IBASE + NSEP*(JT-1)
               JW = JT+1
               JREV = IPERMU(JW)
               AR(J) = H2D(IT,JREV)
            ENDDO
         ENDDO
         DO JT=1,NM2
            DO IT=1,NSTRN                                                 
               IBASE=(IT-1)*NCLBT+1
               J  = NSEP+IBASE + NSEP*(JT-1)
               JW = JT+1
               H2D(IT,JW) = AI(J)
            ENDDO
         ENDDO
         DO JT=1,NM2
            DO IT=1,NSTRN                                                 
               IBASE=(IT-1)*NCLBT+1
               J  = NSEP+IBASE + NSEP*(JT-1)
               JW = JT+1
               JREV = IPERMU(JW)
               AI(J) = H2D(IT,JREV)
            ENDDO
         ENDDO
C---------------------------------------- 
      ELSE
C---------------------------------------- Case II
         DO JT=1,NM2
            DO IT=1,NSTRN                                                 
               IBASE=(IT-1)*NCLBT+1
               J  = NSEP+IBASE + NSEP*(JT-1)
               JW = JT+1
               JREV = IPERMU(JW)
C               WORK(1,JREV) = AR(J)
               H2D(IT,JREV) = AR(J)
            ENDDO
         ENDDO
         DO JT=1,NM2
            DO IT=1,NSTRN                                                 
               IBASE=(IT-1)*NCLBT+1
               J  = NSEP+IBASE + NSEP*(JT-1)
               JW = JT+1
C               AR(J) = WORK(1,JW)
               AR(J) = H2D(IT,JW)
            ENDDO
         ENDDO
         DO JT=1,NM2
            DO IT=1,NSTRN                                                 
               IBASE=(IT-1)*NCLBT+1
               J  = NSEP+IBASE + NSEP*(JT-1)
               JW = JT+1
               JREV = IPERMU(JW)
               H2D(IT,JREV) = AI(J)
            ENDDO
         ENDDO
         DO JT=1,NM2
            DO IT=1,NSTRN                                                 
               IBASE=(IT-1)*NCLBT+1
               J  = NSEP+IBASE + NSEP*(JT-1)
               JW = JT+1
               AI(J) = H2D(IT,JW)
            ENDDO
         ENDDO
C---------------------------------------- 
      ENDIF
C---------------------------------------- 
      RETURN                                                            
      END                                                               
