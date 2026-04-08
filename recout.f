










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
       SUBROUTINE recout (KMX,JMX,IMX,K1,K2,J1,J2,I1,I2,
     $                   X,Y,Z,A,KANGEO,KANREC,CIDREC,TIMEPH,IINIT)
C*MGLET*****************************************************************
C     R E C O U T     Schreiben von Felder im RECORD-Format, simultan
C
C*MGLET*****************************************************************
C
C PARAM: 
C        KMX, JMX, IMX  - DIMENSIONEN DES GITTERS 
C        K1,K2,J1,J2,I1,I2 - INDEX-GRENZEN DES OUTPUTS
C        A(KMX,JMX,IMX) - FELD, DAS RAUSGESCHRIEBEN WIRD
C        KANGEO         - KANAL D. GEOMETRIE-FILES
C        KANREC         - KANAL DER RECORDS
C        CIDREC         - CHARACTER (LEN=10)-VARIABLE, KENNSTRING
C        TIMEPH         - ZEITPUNKT, DER RAUSGESCHRIEBEN WIRD
C        IINIT          - STEUERT OUTPUT:
C                         IINIT=0: NUR FELD WIRD GESCHRIEBEN
C                         IINIT=1: RECORD WIRD GESCHREIBEN MIT TIMEPH
C                         IINIT=2: GEOMETRIE-FILE WIRD GESCHRIEBEN, RECORD
C                                  KANAL WIRD INITIALISIERT, MIT TIMEPH
C                         IINIT=3: GEOMETRIE-FILE UND RECORD-FILE WERDEN
C                                  INITIALISIERT UND GESCHRIEBEN
C                         IINIT=4: TIMEPH WIRD GESCHRIEBEN
C
C VERS:   1.10.96 (MM)  : ORIGINAL
C
C
C*MGLET************************************************************

C
C
      REAL A(KMX,JMX,IMX),
     $       X(IMX),            Y(JMX),            Z(KMX)

       CHARACTER (LEN=16) CIDREC
C
	write (44,*) 'recout called',kangeo,kanrec,timeph,iinit
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
 6050  FORMAT (6(E12.5E3,1X))
 6060  FORMAT (4(I9,1X))
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                        INITIALISING GEOMETRY-FILE

       IF (IINIT .EQ. 2 .OR. IINIT .EQ. 3) THEN
          open  (KANGEO,form='FORMATTED')
          write (KANGEO,'(A)') CIDREC
          write (KANGEO,'(A)') '   X   '
          write (KANGEO,6050 ) (x(i),i=i1,i2)
          write (KANGEO,'(A)') '   Y   '
          write (KANGEO,6050 ) (y(j),j=j1,j2)
          write (KANGEO,'(A)') '   Z   '
          write (KANGEO,6050 ) (z(k),k=k1,k2)
          write (KANGEO,'(A)') ' IARR '
          write (KANGEO,6060 ) (i,i=i1,i2)
          write (KANGEO,'(A)') ' JARR '
          write (KANGEO,6060 ) (j,j=j1,j2)
          write (KANGEO,'(A)') ' KARR '
          write (KANGEO,6060 ) (k,k=k1,k2)
          write (KANGEO,'(A)') ' NTREC'
          write (KANGEO,6060 ) 1
       ENDIF
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                        INITIALISING RECORD-FILE
C
       IF (IINIT .EQ. 2 .OR. IINIT .EQ. 3) THEN
          open  (KANREC,form='UNFORMATTED')
          write (KANREC) CIDREC
       ENDIF

       IF  (IINIT .EQ. 1 .OR. IINIT .EQ. 2 
     $ .OR. IINIT .EQ. 3 .OR. IINIT .EQ. 4) THEN
          write (KANREC) timeph
       ENDIF

       IF (IINIT .EQ. 0 .OR. IINIT .EQ. 1 .OR. IINIT .EQ. 3) THEN
          write (KANREC) (((A(K,J,I),K=K1,K2),J=J1,J2),I=I1,I2)
       ENDIF

      RETURN
      END

