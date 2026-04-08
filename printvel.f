










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
      SUBROUTINE PRINTVEL (KK,JJ,II,U,V,W,P,X,Y,Z,TIMEPH)

      REAL         U(KK,JJ,II),   V(KK,JJ,II),   W(KK,JJ,II),
     $             P(KK,JJ,II)
      REAL         X(II), Y(JJ),Z(KK)

     
C      REAL       Uex (KK,JJ,II),Vex (KK,JJ,II),Wex (KK,JJ,II),
C     $           Pex (KK,JJ,II)

       REAL 	tm1,umax,vmax,wmax,pmax,time1,time2,Uex,Vex,Wex,Pex
     $          ue1,ve1,we1,pe1,ue2,ve2,we2,pe2,rue,rve,rpe,maxdPe,dp
       INTEGER 	I,J,K,Imu,Jmu,Imv,Jmv,Imp,Jmp,Imdp,Jmdp


      K  = 3

      hx=X(3)-X(2)
      hy=hx
C      DO I=3,II
C         DO J=3,JJ
C	 write(98,*) 'Uex(',J-2,',',I-2,')=',Uex(K,J,I),';'
C	 write(101,*) (x(i)+0.5*hx),y(j),Uex(K,J,I)
C	 write(103,*) (X(I)+0.5*hx),Y(J),U(K,J,I)
C	 write(104,*) (X(I)),Y(J)+0.5*hy,V(K,J,I)
C	 write(105,*) (X(I)),Y(J),P(K,J,I)
	 
C	 ENDDO
C      ENDDO


C-----------------print max
      umax=0.0
      vmax=0.0
      wmax=0.0
      pmax=0.0
      tm1=0.0

      K=3
     
      DO J=3,JJ-2
         DO I=3,II-2
	 
	 IF (  abs(U(K,J,I)) .GT. umax) THEN
	 umax=abs(U(K,J,I))
	 ENDIF
	 
	 IF (  abs(V(K,J,I)) .GT. vmax) THEN
	 vmax= abs(V(K,J,I))
         ENDIF
 
	 IF (  abs(W(K,J,I)) .GT. wmax) THEN
	 wmax= abs(W(K,J,I))
	 ENDIF
	 
	 IF (  abs(P(K,J,I)) .GT. pmax) THEN
	 pmax= abs(P(K,J,I))
	 ENDIF
	 
	 ENDDO
c	 write(6,*) P(K,J,I)
      ENDDO
	write(6,*) 'oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO'
	write(6,*) 'MAXU',umax,'MAXV',vmax
	write(6,*) 'MAXW',wmax,'MAXP',pmax




C  compute maximum error


 
	  ue1=0.0
	  ve1=0.0
	  We1=0.0
	  pe1=0.0
	  ue2=0.0
	  ve2=0.0
	  we2=0.0
	  pe2=0.0
	  rue=0.0
	  rve=0.0
	  rpe=0.0
	  dpmax=0.0
      
      time1=exp(-2.0*TIMEPH/100.0)
      time2=exp(-4.0*TIMEPH/100.0)
     
      hx=X(2)-X(1)	
      hy=hx
	DO I=3,II-2
	   DO J=3,JJ-2

      Uex= time1*(-4.0*cos(X(I)+hx/2.)*sin(0.5*hx)*sin(0.5*hy)
     &             *sin(Y(J))/(hx*hy))
      Vex= time1*4.0*sin(X(I))*sin(0.5*hx)*sin(0.5*hy)
     &             *cos(Y(J)+hy/2.)/(hx*hy)
      Pex= time2*(0.125*(sin(2.*X(I)-hx)*hy+sin(2.*Y(J)-hy)*hx
     &          -sin(2.*X(I)+hx)*hy-sin(2.*Y(J)+hy)*hx)/(hx*hy))
      Wex= 0.0
      
      dpex= (0.125*(sin(2.*X(I+1)-hx)*hy+sin(2.*Y(J)-hy)*hx
     &          -sin(2.*X(I+1)+hx)*hy-sin(2.*Y(J)+hy)*hx)/(hx*hy))-
     &     (0.125*(sin(2.*X(I-1)-hx)*hy+sin(2.*Y(J)-hy)*hx
     &          -sin(2.*X(I-1)+hx)*hy-sin(2.*Y(J)+hy)*hx)/(hx*hy))
     
      dp= P(K,J,I+1)-P(K,J,I-1)

           tm1=abs( U(K,J,I)- Uex )
	   IF (  tm1 .GT. ue1 ) THEN
	   ue1=tm1
	   rue=tm1/Uex
	   Imu=I
	   Jmu=J
	   ENDIF
	   ue2=ue2+tm1**2
	   
	   tm1=abs( V(K,J,I)- Vex )
	   IF (  tm1 .GT. ve1 ) THEN
	   ve1=tm1
	   rve=tm1/Vex
	   Imv=I
	   Jmv=J
	   ENDIF
	   ve2=ve2+tm1**2
	   
	   tm1=abs( W(K,J,I)- Wex )
	   IF (  tm1 .GT. we1 ) THEN
	   we1=tm1
C no use
C	   rwe=tm1/Wex
	   ENDIF
	   we2=we2+tm1**2

	   tm1=abs( P(K,J,I)- Pex )
	   IF (  tm1 .GT. pe1 ) THEN
	   pe1=tm1
	   rpe=tm1/Pex
   	   Imp=I
	   Jmp=J
	   ENDIF
	   pe2=pe2+tm1**2
	   
	   tm1=abs( dp- dpex )	   
   	   IF (  tm1 .GT. dpmax ) THEN
	   dpmax=tm1
	   rpe=tm1/Pex
   	   Imdp=I
	   Jmdp=J
	   ENDIF      
	
	   ENDDO
	ENDDO
    
         Ue = 0.0 
         Ve = 0.0
         We = 0.0

    	energy=0.0
	  K=KK/2
      DO I=3,II-2
        DO J=3,JJ-2
         Ue = U(K,J,I) 
         Ve = V(K,J,I)
         We = W(K,J,I)
	    energy=energy+0.5*(Ue**2+Ve**2+We**2)

        ENDDO
      ENDDO

      write (6,*) 'KINETIC_ENERGEY=',energy
      write (6,*) 'UEi',ue1,Imu,Jmu
      write (6,*) 'VEi',ve1,Imv,Jmv
      write (6,*) 'WEi',we1
      write (6,*) 'PEi',pe1,Imp,Jmp
      write (6,*) 'dPEi',dpmax,Imdp,Jmdp
      write (6,*) 'UEel',ue2,'VEel',ve2
      write (6,*) 'WEel',we2,'PEel',pe2
      write (6,*) 'Rue=',rue,'Rve=',rve
      write (6,*) 'Rpe=',rpe
      write (6,*) 'Time=',TIMEPH




      RETURN
      END


