










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
      subroutine init_filter_periodic(N,GQ_N2)
C------------------------------------------------------
C     Constructs desired transfere funktion
C     Copy from Holger Foysi
C     by Florian Schwertfirm
C     17.06.04
C     number of nterms not yet checked
C     4 or 6 ?? 
C-----------------------------------------------------
 
      integer N, i, nu, nterms
      real    alpha, beta, pi, dw
      real    a, b, c, d
      real    w, r
      real    GQ_N2(N+1),GW(N+1),Q_N(N+1),GQ_N(N+1)
      real    w1, w2, T1, T2, a11, a22, a21, a12, DET, R1, R2

      pi = ACOS(-1.0)
      w1 = 0.72
      w2 = (((w1-0.4)/(1.0-0.4)-0.5)*0.2+1.0)*(1.0-w1)/2.0 + w1
      w2 = w2*pi
      w1 = w1*pi
      T1 = 0.95*16.0
      T2 = 0.5*16.0

      a11 = 12.0 + (16.0 -2.0*T1)*cos(w1)+4.0*cos(2.0*w1)
      a12 = 10.0*cos(w1) + (16.0 - 2.0*T1)*cos(2.0*w1)+6.0*cos(3.0*w1)

      a21 = 12.0 + (16.0 -2.0*T2)*cos(w2)+4.0*cos(2.0*w2)
      a22 = 10.0*cos(w2)+(16.0-2.0*T2)*cos(2.0*w2)+6.0*cos(3.0*w2)

      DET = a11*a22 - a21*a12
      R1 = T1 - 8.0 - 9.0*cos(w1) + cos(3.0*w1)
      R2 = T2 - 8.0 - 9.0*cos(w2) + cos(3.0*w2)

      alpha = (R1*a22 -R2*a12)/DET
      beta = (a11*R2 -a21*R1)/DET

      nterms = 0

C---- Transfere Funktion  C.2.10.b acc. Lelel

      a = 1.0/4.0*(2.0 + 3.0*alpha)
      b = 1.0/16.0*(9.0 + 16.0*alpha + 10.0*beta)
      c = 1.0/4.0*(alpha + 4.0*beta)
      d = 1.0/16.0*(6.0*beta - 1.0)

C---- Construct filter response funktion G(w) formula C.2.2

      dw = pi/float(N)
 
      do i=1,N+1
       w = (i-1)*dw
       GW(i) = (a + b*cos(w) + c*cos(2*w) + d*cos(3*w))/
     $          (1.0 + 2.0*alpha*cos(w)+2.0*beta*cos(2*w))
      enddo

C---- Construct approx. dec. op. Q_N

      do i=1,N+1
       r = 1.0-GW(i)
       Q_N(i) = 1.0
       do nu = 1,nterms
        Q_N(i) = Q_N(i)*r + 1.0
       enddo
      enddo

C---- Construct G*Q_N and (G*Q_N)**2
      do i=1,n+1
       if(GW(i) .GT. 1.0) GW(i) = 1.0
       if(GW(i) .LT. 0.0) GW(i) = 0.0
       GQ_N2(i) = GW(i)
      enddo

C---- Output for plot

C      do i=1,n+1
C       write(30,*)float(i-1)*dw/pi,GW(i),GQ_N(i),GQ_N2(i)
C      enddo

      RETURN
      END
