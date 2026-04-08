










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
      SUBROUTINE PREPROC_FIELD(kk,jj,ii,PHI,H3D1,H3D2,H3D3,H3D4,H3D5,
     $                         DZ,KON)

C..... Prerpocessing the field to achieve better resolution
C      26.07.04
C      by Florian Schwertfirm
C..........................................................

      IMPLICIT NONE
      integer i,j,k
      integer kk,jj,ii
      

      logical KON
      real PHI(kk,jj,ii),H3D1(kk,jj,ii),H3D2(kk,jj,ii),H3D3(kk,jj,ii),
     $     H3D4(kk,jj,ii),H3D5(kk,jj,ii)
      real DZ(kk)
c..... filter coefficients
      real a1,b1,a2,b2,a11,b11

      if(KON) THEN
      a1 = 9.0/140.0+sqrt(265.0)/140.0
      b1 = 23.0/35.0
      a2 = -21.0/92.0+sqrt(265.0)*7.0/276.0
      b2 = 4.0/6.0
C      a1 = 1.0/6.0
C      b1 = 4.0/6.0
C      a2 = a1
C      b2 = b1

      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D1)
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D2)
c......................................................... x-direction
c..... first filter step
        
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
       H3D1(k,j,i) = a1*PHI(k,j,i-1) + b1*PHI(k,j,i) + 
     $               a1*PHI(k,j,i+1)
        enddo
       enddo
      enddo

C---- Preiodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D1(k,j,ii-1) = H3D1(k,j,3)
        H3D1(k,j,2) = H3D1(k,j,ii-2)
       enddo
      enddo
C---- Periodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D1(k,2,i) = H3D1(k,jj-2,i)
        H3D1(k,jj-1,i) = H3D1(k,3,i)
       enddo
      enddo
C

c..... second filter step
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
       H3D2(k,j,i) = a2*H3D1(k,j,i-1) + b2*H3D1(k,j,i) +
     $                  a2*H3D1(k,j,i+1)
        enddo
       enddo
      enddo

C---- Preiodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D2(k,j,ii-1) = H3D2(k,j,3)
        H3D2(k,j,2) = H3D2(k,j,ii-2)
       enddo
      enddo
C---- Periodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D2(k,2,i) = H3D2(k,jj-2,i)
        H3D2(k,jj-1,i) = H3D2(k,3,i)
       enddo
      enddo

c      write(6,*)'PR: ',H3D2(65,70,60)

c..... deconvolution
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
        H3D3(k,j,i)= 3.*PHI(k,j,i) 
     $            - 3.*H3D1(k,j,i) + H3D2(k,j,i)
        enddo
       enddo
      enddo
C.... Preiodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D3(k,j,ii-1) = H3D3(k,j,3)
        H3D3(k,j,2) = H3D3(k,j,ii-2)
       enddo
      enddo
C---- Periodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D3(k,2,i) = H3D3(k,jj-2,i)
        H3D3(k,jj-1,i) = H3D3(k,3,i)
       enddo
      enddo


c 
c....................................................... y-direction

c..... first filter step
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D1)
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D2)
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
       H3D1(k,j,i) = a1*PHI(k,j-1,i) + b1*PHI(k,j,i) +
     $                  a1*PHI(k,j+1,i)
        enddo
       enddo
      enddo

C---- Periodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D1(k,j,ii-1) = H3D1(k,j,3)
        H3D1(k,j,2) = H3D1(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D1(k,jj-1,i) = H3D1(k,3,i)
        H3D1(k,2,i) = H3D1(k,jj-2,i)
       enddo
      enddo

c..... second filter step
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
       H3D2(k,j,i) = a2*H3D1(k,j-1,i) + b2*H3D1(k,j,i) +
     $                  a2*H3D1(k,j+1,i)
        enddo
       enddo
      enddo

C---- Preiodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D2(k,j,ii-1) = H3D2(k,j,3)
        H3D2(k,j,2) = H3D2(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D2(k,jj-1,i) = H3D2(k,3,i)
        H3D2(k,2,i) = H3D2(k,jj-2,i)
       enddo
      enddo

c..... deconvolution
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
        H3D4(k,j,i)= 3.*PHI(k,j,i) 
     $               - 3.*H3D1(k,j,i) + H3D2(k,j,i)
        enddo
       enddo
      enddo

C.... Preiodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D4(k,j,ii-1) = H3D4(k,j,3)
        H3D4(k,j,2) = H3D4(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D4(k,jj-1,i) = H3D4(k,3,i)
        H3D4(k,2,i) = H3D4(k,jj-2,i)
       enddo
      enddo
c....................................................... z-direction
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D1)
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D2)

CC Problem der Randbedingung
C    Moegl. 1: nicht aequidistanter Filter
C    Moegl. 2: Randpunkt unveraendert
C    Moegl. 3: Geschwindigkeit an der Wand gespiegelt
C..... first filter step
C       a1 = 1.0/6.0
C       b1 = 4.0/6.0
C       a2 = 1.0/6.0
C       b2 = 4.0/6.0

      do i=3,ii-2
       do j=3,jj-2
Cc...... unterer Rand
C.......................................................................     Moegl. 1
C      H3D1(3,j,i) = 1.0/3.0*DZ(3)/(0.5*DZ(2)+DZ(3))*PHI(2,j,i) +
C     $              2.0/3.0*PHI(3,j,i) +
C     $              1.0/3.0*1.0/(1.0+DZ(3)*2.0/DZ(2))*PHI(4,j,i)
C.......................................................................     Moegl. 2
C       H3D1(3,j,i) = PHI(3,j,i)
C.......................................................................     Moegl. 3
       H3D1(3,j,i) = -a1*PHI(3,j,i) + b1*PHI(3,j,i) +
     $               a1*PHI(4,j,i)
Cc..... im Gebiet
        do k=4,kk-3
C      H3D1(k,j,i) = DZ(k)/(3.0*(DZ(k) + DZ(k-1)))*PHI(k-1,j,i) +
C     $               2.0/3.0                     *PHI(k,j,i) +
C     $               1.0/3.0/(DZ(k)/DZ(k-1)+1.0) *PHI(k+1,j,i)
       H3D1(k,j,i) = a1*PHI(k-1,j,i) + b1*PHI(k,j,i) +
     $               a1*PHI(k+1,j,i)
c       H3D1(k,j,i) = a1*H3D2(k-1,j,i) + b1*H3D2(k,j,i) +
c     $               a1*H3D2(k+1,j,i)
        enddo
Cc...... oberer Rand
C      H3D1(kk-2,j,i)=
C     $    1.0/3.0*(0.5*DZ(kk-2))/(0.5*DZ(kk-2)+DZ(kk-3))*PHI(kk-3,j,i) +
C     $    2.0/3.0*PHI(kk-2,j,i) +
C     $    1.0/3.0*1.0/(1.0+DZ(kk-2)/(2*DZ(kk-3)))*PHI(kk-1,j,i)
C       H3D1(kk-2,j,i) = PHI(kk-2,j,i)
       H3D1(kk-2,j,i) = a1*PHI(kk-3,j,i) + b1*PHI(kk-2,j,i) 
     $               -a1*PHI(kk-2,j,i)
C       H3D1(kk-2,j,i) = -1./6.*PHI(kk-3,j,i) + 5.*PHI(kk-2,j,i) +
C     $               11./3.*PHI(kk-2,j,i)
       enddo
      enddo

C---- Periodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D1(k,j,ii-1) = H3D1(k,j,3)
        H3D1(k,j,2) = H3D1(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D1(k,jj-1,i) = H3D1(k,3,i)
        H3D1(k,2,i) = H3D1(k,jj-2,i)
       enddo
      enddo

CCc..... second filter step
C
      do i=3,ii-2
       do j=3,jj-2
CCc...... unterer Rand
C      H3D2(3,j,i) = 1.0/3.0*DZ(3)/(0.5*DZ(2)+DZ(3))*H3D1(2,j,i) +
C     $              2.0/3.0*H3D1(3,j,i) +
C     $              1.0/3.0*1.0/(1.0+DZ(3)*2.0/DZ(2))*H3D1(4,j,i)

C       H3D2(3,j,i) = PHI(3,j,i)
       H3D2(3,j,i) = -a2*H3D1(3,j,i) +b2*H3D1(3,j,i) +
C       H3D2(3,j,i) = -1./6.*PHI(3,j,i) +4./6.*H3D1(3,j,i) +
     $               a2*H3D1(4,j,i)
CCCCc..... im Gebiet
        do k=4,kk-3
C      H3D2(k,j,i) = DZ(k)/(3.0*(DZ(k) + DZ(k-1)))*H3D1(k-1,j,i) +
C     $               2.0/3.0                     *H3D1(k,j,i) +
C     $               1.0/3.0/(DZ(k)/DZ(k-1)+1.0)*H3D1(k+1,j,i)
       H3D2(k,j,i) = a2*H3D1(k-1,j,i) +b2*H3D1(k,j,i) +
     $               a2*H3D1(k+1,j,i)
        enddo
CCCc...... oberer Rand
C      H3D2(kk-2,j,i)=
C     $   1.0/3.0*(0.5*DZ(kk-2))/(0.5*DZ(kk-2)+DZ(kk-3))*H3D1(kk-3,j,i) +
C     $   2.0/3.0*H3D1(kk-2,j,i) +
C     $   1.0/3.0*1.0/(1.0+DZ(kk-2)/(2*DZ(kk-3)))*H3D1(kk-1,j,i)



C       H3D2(kk-2,j,i) = PHI(kk-2,j,i)
       H3D2(kk-2,j,i) = a2*H3D1(k-3,j,i) +b2*H3D1(kk-2,j,i) 
C     $               -1./6.*PHI(kk-2,j,i)
     $               -a2*H3D1(kk-2,j,i)
C
       enddo
      enddo
CC
C---- Preiodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D2(k,j,ii-1) = H3D2(k,j,3)
        H3D2(k,j,2) = H3D2(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D2(k,jj-1,i) = H3D2(k,3,i)
        H3D2(k,2,i) = H3D2(k,jj-2,i)
       enddo
      enddo

CCc..... deconvolution
      do i=2,ii-1
       do j=2,jj-1
        do k=3,kk-2
C        H3D5(k,j,i)= 3.*PHI(k,j,i)
C     $               - 3.*H3D1(k,j,i) + H3D2(k,j,i)
         H3D5(k,j,i) = PHI(k,j,i)
        enddo
       enddo
      enddo

C---- Periodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D5(k,j,ii-1) = H3D5(k,j,3)
        H3D5(k,j,2) = H3D5(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D5(k,jj-1,i) = H3D5(k,3,i)
        H3D5(k,2,i) = H3D5(k,jj-2,i)
       enddo
      enddo

CCC.......................................................................
C............ Here the same for diffusive term .........................
C.......................................................................
      else
      a1 = 23.0/780+1.0/260*sqrt(345.0)
      b1 = 161.0/195.0
      a2 = -13.0/168.0+13.0/1288.0*sqrt(345.0)
      b2 = 10.0/12.0
C      a1 = 1.0/12.0
C      b1 = 10.0/12.0
C      a2 = 1.0/12.0
C      b2 = 10.0/12.0

      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D1)
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D2)
Cc..... x-direction
Cc..... first filter step
C  
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
       H3D1(k,j,i) = a1*PHI(k,j,i-1) + b1*PHI(k,j,i) + 
     $               a1*PHI(k,j,i+1)
        enddo
       enddo
      enddo

C---- Preiodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D1(k,j,ii-1) = H3D1(k,j,3)
        H3D1(k,j,2) = H3D1(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D1(k,jj-1,i) = H3D1(k,3,i)
        H3D1(k,2,i) = H3D1(k,jj-2,i)
       enddo
      enddo


c..... second filter step
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
       H3D2(k,j,i) = a2*H3D1(k,j,i-1) + b2*H3D1(k,j,i) +
     $                  a2*H3D1(k,j,i+1)
        enddo
       enddo
      enddo

CC---- Preiodic boundary in x
      do j=1,jj
       do k=1,kk
        H3D2(k,j,ii-1) = H3D2(k,j,3)
        H3D2(k,j,2) = H3D2(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D2(k,jj-1,i) = H3D2(k,3,i)
        H3D2(k,2,i) = H3D2(k,jj-2,i)
       enddo
      enddo


c..... deconvolution
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
        H3D3(k,j,i)= 3.*PHI(k,j,i) 
     $            - 3.*H3D1(k,j,i) + H3D2(k,j,i)
        enddo
       enddo
      enddo
c...... Periodic in x
      do j=1,jj
       do k=1,kk
        H3D3(k,j,ii-1) = H3D3(k,j,3)
        H3D3(k,j,2) = H3D3(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D3(k,jj-1,i) = H3D3(k,3,i)
        H3D3(k,2,i) = H3D3(k,jj-2,i)
       enddo
      enddo

Cc 
Cc..... y-direction
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D1)
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D2)
Cc..... first filter step
C
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
       H3D1(k,j,i) = a1*PHI(k,j-1,i) + b1*PHI(k,j,i) +
     $                  a1*PHI(k,j+1,i)
        enddo
       enddo
      enddo

c...... Periodic in x
      do j=1,jj
       do k=1,kk
        H3D1(k,j,ii-1) = H3D1(k,j,3)
        H3D1(k,j,2) = H3D1(k,j,ii-2)
       enddo
      enddo
c.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D1(k,jj-1,i) = H3D1(k,3,i)
        H3D1(k,2,i) = H3D1(k,jj-2,i)
       enddo
      enddo

c..... second filter step
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
       H3D2(k,j,i) = a2*H3D1(k,j-1,i) + b2*H3D1(k,j,i) +
     $                  a2*H3D1(k,j+1,i)
        enddo
       enddo
      enddo
c...... Periodic in x
      do j=1,jj
       do k=1,kk
        H3D2(k,j,ii-1) = H3D2(k,j,3)
        H3D2(k,j,2) = H3D2(k,j,ii-2)
       enddo
      enddo
Cc.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D2(k,jj-1,i) = H3D2(k,3,i)
        H3D2(k,2,i) = H3D2(k,jj-2,i)
       enddo
      enddo

Cc..... deconvolution
      do i=3,ii-2
       do j=3,jj-2
        do k=3,kk-2
        H3D4(k,j,i)= 3.*PHI(k,j,i) 
     $               - 3.*H3D1(k,j,i) + H3D2(k,j,i)
        enddo
       enddo
      enddo
c...... Periodic in x
      do j=1,jj
       do k=1,kk
        H3D4(k,j,ii-1) = H3D4(k,j,3)
        H3D4(k,j,2) = H3D4(k,j,ii-2)
       enddo
      enddo
Cc.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D4(k,jj-1,i) = H3D4(k,3,i)
        H3D4(k,2,i) = H3D4(k,jj-2,i)
       enddo
      enddo

c............................................................................. z-direction

c..... first filter step
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D1)
      CALL DPHI0(kk,jj,ii,kk,jj,ii,H3D2)

      do i=3,ii-2
       do j=3,jj-2
c--------------------------------------------------------------------
C      H3D1(3,j,i) = 1.0/6.0*DZ(3)/(0.5*DZ(2)+DZ(3))*PHI(2,j,i) +
C     $              10.0/12.0*PHI(3,j,i) +
C     $              1.0/6.0*1.0/(1.0+DZ(3)*2.0/DZ(2))*PHI(4,j,i)
C--------------------------------------------------------------------
C       H3D1(3,j,i) = PHI(3,j,i)
c--------------------------------------------------------------------
       H3D1(3,j,i) = -a1*PHI(3,j,i) + b1*PHI(3,j,i) +
     $               a1*PHI(4,j,i)
c--------------------------------------------------------------------
C       H3D1(3,j,i) = -1./6.*PHI(3,j,i) + 5./6.*PHI(3,j,i) +
C     $               1./6.*PHI(4,j,i)
c--------------------------------------------------------------------
CCc..... im Gebiet
        do k=4,kk-3
C        H3D1(k,j,i) = DZ(k)/(6.0*(DZ(k) + DZ(k-1)))*PHI(k-1,j,i) +
C     $               10.0/12.0                    *PHI(k,j,i) +
C     $               1.0/(6.0*(DZ(k)/DZ(k-1)+1.0))*PHI(k+1,j,i)
       H3D1(k,j,i) = a1*PHI(k-1,j,i) + b1*PHI(k,j,i) +
     $               a1*PHI(k+1,j,i)
        enddo
CCc...... oberer Rand
c-------------------------------------------------------------------
C      H3D1(kk-2,j,i)=
C     $    1.0/6.0*(0.5*DZ(kk-2))/(0.5*DZ(kk-2)+DZ(kk-3))*PHI(kk-3,j,i) +
C     $    10.0/12.0*PHI(kk-2,j,i) +
C     $    1.0/6.0*1.0/(1.0+DZ(kk-2)/(2*DZ(kk-3)))*PHI(kk-1,j,i)
c-------------------------------------------------------------------
C       H3D1(kk-2,j,i) = PHI(kk-2,j,i)
c-------------------------------------------------------------------
       H3D1(kk-2,j,i) = a1*PHI(kk-3,j,i) + b1*PHI(kk-2,j,i) 
     $               -a1*PHI(kk-2,j,i)
c-------------------------------------------------------------------
C       H3D1(kk-2,j,i) = 1./6.*PHI(kk-3,j,i) + 5./6.*PHI(kk-2,j,i) 
C     $               -1./6.*PHI(kk-2,j,i)

       enddo
      enddo

c...... Periodic in x
      do j=1,jj
       do k=1,kk
        H3D1(k,j,ii-1) = H3D1(k,j,3)
        H3D1(k,j,2) = H3D1(k,j,ii-2)
       enddo
      enddo
Cc.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D1(k,jj-1,i) = H3D1(k,3,i)
        H3D1(k,2,i) = H3D1(k,jj-2,i)
       enddo
      enddo
CC
CCc..... second filter step
      do i=3,ii-2
       do j=3,jj-2
Cc...... unterer Rand
C      H3D2(3,j,i) = 1.0/3.0*DZ(3)/(0.5*DZ(2)+DZ(3))*H3D1(2,j,i) +
C     $              2.0/3.0*H3D1(3,j,i) +
C     $              1.0/3.0*1.0/(1.0+DZ(3)*2.0/DZ(2))*H3D1(4,j,i)
c------------------------------------------------------------------
C       H3D2(3,j,i) = PHI(3,j,i)
c------------------------------------------------------------------
       H3D2(3,j,i) = -a2*H3D1(3,j,i) + b2*H3D1(3,j,i) +
C       H3D2(3,j,i) = -1./12.*PHI(3,j,i) + 10./12.*H3D1(3,j,i) +
     $               a2*H3D1(4,j,i)
c------------------------------------------------------------------
C       H3D2(3,j,i) = -1./6.*H3D1(3,j,i) + 5./6.*H3D1(3,j,i) +
C     $               1./6.*H3D1(4,j,i)

CCc..... im Gebiet
        do k=4,kk-3
C      H3D2(k,j,i) = DZ(k)/(6.0*(DZ(k) + DZ(k-1)))*H3D1(k-1,j,i) +
C     $               10.0/12.0*H3D1(k,j,i) +
C     $               1.0/(6.0*(DZ(k)/DZ(k-1)+1.0))*H3D1(k+1,j,i)
       H3D2(k,j,i) = a2*H3D1(k-1,j,i) + b2*H3D1(k,j,i) + 
     $               a2*H3D1(k+1,j,i)
        enddo
Cc...... oberer Rand

C      H3D2(kk-2,j,i)=
C     $   1.0/6.0*(0.5*DZ(kk-2))/(0.5*DZ(kk-2)+DZ(kk-3))*H3D1(kk-3,j,i) +
C     $   10.0/12.0*H3D1(kk-2,j,i) +
C     $   1.0/6.0*1.0/(1.0+DZ(kk-2)/(2*DZ(kk-3)))*H3D1(kk-1,j,i)
C-------------------------------------------------------------------
C       H3D2(kk-2,j,i) = PHI(kk-2,j,i)
c-------------------------------------------------------------------
       H3D2(kk-2,j,i) = a2*H3D1(kk-3,j,i) + b2*H3D1(kk-2,j,i)
C     $               -1./12.*PHI(kk-2,j,i)
     $               -a2*H3D1(kk-2,j,i)
c-------------------------------------------------------------------
C       H3D2(kk-2,j,i) = 1./6.*H3D1(kk-3,j,i) + 5./6.*H3D1(kk-2,j,i) 
C     $               -1./6.*H3D1(kk-2,j,i)
       enddo
      enddo

c...... Periodic in x
      do j=1,jj
       do k=1,kk
        H3D2(k,j,ii-1) = H3D2(k,j,3)
        H3D2(k,j,2) = H3D2(k,j,ii-2)
       enddo
      enddo
Cc.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D2(k,jj-1,i) = H3D2(k,3,i)
        H3D2(k,2,i) = H3D2(k,jj-2,i)
       enddo
      enddo

CCc..... deconvolution
      do i=2,ii-1
       do j=2,jj-1
        do k=3,kk-2
C        H3D5(k,j,i)= 3.*PHI(k,j,i)
C     $               - 3.*H3D1(k,j,i) + H3D2(k,j,i)
         H3D5(k,j,i)=PHI(k,j,i)
        enddo
       enddo
      enddo

c...... Periodic in x
      do j=1,jj
       do k=1,kk
        H3D5(k,j,ii-1) = H3D5(k,j,3)
        H3D5(k,j,2) = H3D5(k,j,ii-2)
       enddo
      enddo
Cc.... Preiodic boundary in y
      do i=1,ii
       do k=1,kk
        H3D5(k,jj-1,i) = H3D5(k,3,i)
        H3D5(k,2,i) = H3D5(k,jj-2,i)
       enddo
      enddo
      endif
C
      return
      end

