










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
      SUBROUTINE ADJDPDX(KK,JJ,II,NBND,ITSTEP,ITMIT,
     $     X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,U,V,W,UO,VO,WO,B,
     $     WALLSSX,WALLSSY,WALLSSZ,GRADPX,
     $     UBULK,UBULKX,GMOL,RHO,UGRID,DT,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP,BP,GRADPXOLD)
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     BERECHNUNG DER GEMITTELTEN WANDSCHUBSPANNUNG
C
C     VERSION VOM 25.1.1998 (MM) 
C                           (NICHT ALLE SCHUBSPANNUNGEN IMPLEMENTIERT!)
C     AENDERUNG : N.Peller 08.09.03 BP Feld hinzugefuegt
C                                   ifdef 
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      REAL U(KK,JJ,II), V(KK,JJ,II), W(KK,JJ,II)
      REAL UO(KK,JJ,II), VO(KK,JJ,II), WO(KK,JJ,II)
      REAL B(KK,JJ,II),BP(KK,JJ,II),GRADPXTMP

      REAL DDX(II),     DDY(JJ),     DDZ(KK)
      REAL  DX(II),      DY(JJ),      DZ(KK)
      REAL   X(II),       Y(JJ),       Z(KK)

      REAL WALLSSX(6), WALLSSY(6), WALLSSZ(6)

      REAL LX,LY,LZ
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      CALL WIFAK  (ITSTEP,ITMIT,WPHI,WKON,WDIF,WSOR,1)

      CALL CHANNELD(KK,JJ,II,NBND,
     $     DX,DY,DZ,DDX,DDY,DDZ,UO,VO,WO,
     $     WALLSSX,WALLSSY,WALLSSZ,GMOL,RHO,UGRID,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)

      CALL CAL_UBULK(KK,JJ,II,NBND,DDX,
     $                    DDY,DDZ,UO,WO,B,
     $                    UGRID,UBULK,BP)

C      WRITE (6,*) 'adj, wphi,wdif,wsor',wphi,wdif,wsor
      SHEARSTRESS = WALLSSX(3)+WALLSSX(4)+WALLSSX(5)+WALLSSX(6)


      GRADPXTMP = GRADPX

      GRADPX = GRADPX + 0.1*(UBULK - UBULKX)
     $                - 0.05*abs(UBULK - UBULKX)*(GRADPXOLD)
      GRADPXOLD = GRADPXTMP 

      RETURN

CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72


      CALL WIFAK  (ITSTEP,ITMIT,WPHI,WKON,WDIF,WSOR,2)


      CALL CHANNELD(KK,JJ,II,NBND,
     $     DX,DY,DZ,DDX,DDY,DDZ,U,V,W,
     $     WALLSSX,WALLSSY,WALLSSZ,GMOL,RHO,UGRID,
     $     NFRO,NBAC,NRGT,NLFT,NBOT,NTOP)

      CALL CAL_UBULK(KK,JJ,II,NBND,DDX,
     $                    DDY,DDZ,U,W,B,
     $                    UGRID,UBULK,BP)

C      WRITE (6,*) 'adj, wphi,wdif',wphi,wdif


      GRADPX = GRADPX + 0.5*(UBULK - UBULKX)


C      write (6,*)'adjusting pressure gradient to: ',gradpx
CCCCC CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC72

      RETURN
      END

