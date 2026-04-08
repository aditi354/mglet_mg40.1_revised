










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
      SUBROUTINE PREPROC_PRESSURE(kk,jj,ii,PHI,H3D1,H3D2,H3D3)

C..... Prerpocessing the field to achieve better resolution
C      26.07.04
C      by Florian Schwertfirm
C..........................................................

      IMPLICIT NONE
      integer i,j,k
      integer kk,jj,ii
      

      logical KON
      real PHI(kk,jj,ii),H3D1(kk,jj,ii),H3D2(kk,jj,ii),H3D3(kk,jj,ii)
      real DZ(kk)
c..... filter coefficients
      real a1,b1,a2,b2,a11,b11

      a1 = 29.0/2240+1.0/560.0*sqrt(565.0)
      b1 = 911.0/1120.0
      a2 = -203.0/5466.0+14.0/2733.0*sqrt(565.0)
      b2 = 4.0/6.0

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
        H3D2(k,j,i)= 3.*PHI(k,j,i) 
     $            - 3.*H3D1(k,j,i) + H3D2(k,j,i)
        enddo
       enddo
      enddo
C.... Preiodic boundary in x
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

      return
      end

