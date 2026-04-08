










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
      SUBROUTINE XTRACTXY (KMXC,JMXC,IMXC,KMX,Z,DZ,DDZ,ZC,DZC,DDZC,
     $                   IN,OUT,
     $                   ISTA,ISTO,JSTA,JSTO,
     $                   IPOS,JPOS,KPOS,KSTAG,IORDER)
C*MGLET*****************************************************************
C     X T R A C T X Y    EXTRACTS AN X-Y-PLANE OUT OF A 3D-FIELD
C*MGLET*****************************************************************
C
C PARAM: 
C        KMXC, JMXC, IMXC  - DIMENSIONEN DES GITTERS 
C       IN(KMXC,JMXC,IMXC) + INPUT-3D-FIELD
C      OUT(    JMXC,IMXC) + OUTPUT-2D-FIELD
C        ISTA,ISTO,JSTA,JSTO
C        KSTART,KSTOP   - BEREICH IN DEM RANDBEDINGUNG GESETZT WIRD
C        IPOS,JPOS,KPOS - POSITION IM NACHBARGITTER, AUF DER DER PUNKT
C                         MIT DEN INDIZES (3,3,3) ZU LIEGEN KOMMT
C
C VERS:  16.09.96 (MM)  : ORIGINAL AUS CONTOPAR ABGELEITET
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************

C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $    IN(KMXC,JMXC,IMXC),  OUT(JMXC,IMXC), 
     $       ZC(KMXC),            Z(KMX),
     $      DZC(KMXC),           DZ(KMX),
     $     DDZC(KMXC),          DDZ(KMX)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                         FINE-GRID: K-INDEX=KMX -1  (TOP)
      KF = KMX - 1
C                         K-INDEX OF COARSE-GRID
C                                   (=PARENT, CALLED NEIGHBOUR)
      KC = KPOS - 1 + (KF-1)/2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                         RUNNING INIDICES
C
      ICSTA = IPOS - 1 + (ISTA-1)/2
      ICSTO = IPOS - 1 + (ISTO-1)/2
      JCSTA = JPOS - 1 + (JSTA-1)/2
      JCSTO = JPOS - 1 + (JSTO-1)/2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                                VARIABLE NOT STAGGERED IN FIRST DIR.
      IF (KSTAG .EQ. 0) THEN

C                                Interpolation of first order
         IF (IORDER .EQ. 1) THEN
            DO IC = ICSTA,ICSTO
               DO JC = JCSTA,JCSTO
                  
                  OUT(JC,IC) = IN(KC,JC,IC)
                  
               ENDDO
            ENDDO
         ENDIF

C                                Interpolation of second order
         IF (IORDER .EQ. 2) THEN
            DO IC = ICSTA,ICSTO
               DO JC = JCSTA,JCSTO
                  
                  OUT(JC,IC)  = IN(KC  ,JC,IC)*0.75 +
     +                          IN(KC-1,JC,IC)*0.25
                  
               ENDDO
            ENDDO
         ENDIF
C                                VARIABLE STAGGERED IN FIRST DIR.
      ELSEIF  (KSTAG .EQ. 1) THEN
C                                Interpolation of first order
         IF (IORDER .EQ. 1) THEN
            
            DO IC = ICSTA,ICSTO
               DO JC = JCSTA,JCSTO
                  
                  OUT(JC,IC) = 0.5*(IN(KC-1,JC,IC) + IN(KC  ,JC,IC))
                  
               ENDDO
            ENDDO
         ENDIF
C                                Interpolation of second order
         IF (IORDER .EQ. 2) THEN
            DO IC = ICSTA,ICSTO
               DO JC = JCSTA,JCSTO
                  
                   OUT(JC,IC)  = 1./DDZC(KC) *
     +                          (IN(KC  ,JC,IC)*DDZ(KF) +
     +                           IN(KC-1,JC,IC)*DDZ(KF+1))
                  
               ENDDO
            ENDDO
         ENDIF

C                                VARIABLE NEGATIV STAGGERED IN FIRST DIR.
      ELSEIF  (KSTAG .EQ. -1) THEN

C                                No  Interpolation necessary
         DO IC = ICSTA,ICSTO
            DO JC = JCSTA,JCSTO
               
               OUT(JC,IC) = IN(KC-1,JC,IC)
               
            ENDDO
         ENDDO

      ENDIF


      RETURN
      END
      
