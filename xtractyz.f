










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
      SUBROUTINE XTRACTYZ (KMXC,JMXC,IMXC,IMX,X,DX,DDX,XC,DXC,DDXC,
     $                   IN,OUT,
     $                   JSTA,JSTO,KSTA,KSTO,
     $                   IPOS,JPOS,KPOS,ISTAG,IORDER)
C*MGLET*****************************************************************
C     X T R A C T Y Z    EXTRACTS AN Y-Z-PLANE OUT OF A 3D-FIELD
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
C VERS:  27.10.97 (MM)  : ORIGINAL AUS XTRACTXY ABGELEITET
C
C DEFINE-DIREKTIVEN     : QUD 
C
C*STARLET***************************************************************

C
      CHARACTER (LEN=1)  ITYP
C
      REAL
     $    IN(KMXC,JMXC,IMXC),  OUT(KMXC,JMXC), 
     $       XC(IMXC),            X(IMX),
     $      DXC(IMXC),           DX(IMX),
     $     DDXC(IMXC),          DDX(IMX)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                         FINE-GRID: I-INDEX=IMX -1  (BAC)
      IIF = IMX - 1
C                         I-INDEX OF COARSE-GRID
C                                   (=PARENT, CALLED NEIGHBOUR)
      IC = IPOS - 1 + (IIF-1)/2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                         RUNNING INIDICES
C
      KCSTA = KPOS - 1 + (KSTA-1)/2
      KCSTO = KPOS - 1 + (KSTO-1)/2
      JCSTA = JPOS - 1 + (JSTA-1)/2
      JCSTO = JPOS - 1 + (JSTO-1)/2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                                VARIABLE NOT STAGGERED IN LAST DIR.
      IF (ISTAG .EQ. 0) THEN

C                                Interpolation of first order
         IF (IORDER .EQ. 1) THEN
            DO JC = JCSTA,JCSTO
               DO KC = KCSTA,KCSTO
                  
                  OUT(KC,JC) = IN(KC,JC,IC)
                  
               ENDDO
            ENDDO
         ENDIF

C                                Interpolation of second order
         IF (IORDER .EQ. 2) THEN
            DO JC = JCSTA,JCSTO
               DO KC = KCSTA,KCSTO
                  
                  OUT(KC,JC)  = IN(KC,JC,IC  )*0.75 +
     +                          IN(KC,JC,IC-1)*0.25
                  
               ENDDO
            ENDDO
         ENDIF
C                                VARIABLE STAGGERED IN LAST DIR.
      ELSEIF  (ISTAG .EQ. 1) THEN
C                                Interpolation of first order
         IF (IORDER .EQ. 1) THEN
            
            DO JC = JCSTA,JCSTO
               DO KC = KCSTA,KCSTO
                  
                  OUT(KC,JC) = 0.5*(IN(KC,JC,IC-1) + IN(KC,JC,IC-1))
                  
               ENDDO
            ENDDO
         ENDIF
C                                Interpolation of second order
         IF (IORDER .EQ. 2) THEN
            DO JC = JCSTA,JCSTO
               DO KC = KCSTA,KCSTO
                  
                   OUT(KC,JC)  = 1./DDXC(IC) *
     +                          (IN(KC,JC,IC  )*DDX(IIF) +
     +                           IN(KC,JC,IC-1)*DDX(IIF+1))
                  
               ENDDO
            ENDDO
         ENDIF

C                                VARIABLE NEGATIV STAGGERED IN FIRST DIR.
      ELSEIF  (ISTAG .EQ. -1) THEN

C                                No  Interpolation necessary
         DO JC = JCSTA,JCSTO
            DO KC = KCSTA,KCSTO
               
               OUT(KC,JC) = IN(KC,JC,IC-1)
               
            ENDDO
         ENDDO

      ENDIF


      RETURN
      END
      
