










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
      SUBROUTINE EXPAN1  (FVT, IDIMF, DELTO, DELT, N, NLARGE)
C*STARLET***************************************************************
C        E X P A N 1      DIE AUF EINER LAENGE VON  N * DELTO  AEQUI-
C                         DISTANT VERTEILTEN FUNKTIONSWERTE  FVT  MIT
C                         DEN STUETZSTELLEN  T_I = FLOAT(I) * DELTO
C                         , WOBEI I = 0,1,2, ... , N   WERDEN
C                         AUF EINER LAENGE VON  (N+1) * DELT = N * DELTO
C                         NEU VERTEILT. (LINEARE INTERPOLATION)
C                         DIES IST BEISPIELSWEISE BEI EINER FOURIER-
C                         ANALYSE NOTWENDIG, WENN   N   U N G E R A D -
C                         Z A H L I G  IST, DIE FFT JEDOCH EINE GERADE
C                         ANZAHL VON STUETZPUNKTEN VERLANGT.
C*STARLET***************************************************************
C
C PARAM: FVT (0:IDIMF)  + ENTHAELT DIE FUNKTIONSWERTE AN DEN ALTEN
C                         STUETZSTELLEN (BEIM EINTRITT IN DIE SUBR.)
C                         UND DIE INTERPOLIERTEN FUNKTIONSWERTE AN DEN
C                         NEUEN STUETZSTELLEN BEIM VERLASSEN DER SUBR.
C        IDIMF          - DIMENSION DES FVT-FELDES
C        DELTO          - ORIGINALABSTAND DER STUETZSTELLEN, AN DENEN
C                         FVT GEGEBEN IST
C        DELT           + ABSTAND DER NEUEN STUETZSTELLEN :
C                         DELT = N/(N+1)*DELTO
C        N              - ANZAHL DER STUETZSTELLEN (UNGERADZAHLIG !)
C        NLARGE         + NLARGE = N + 1 (GERADZAHLIG)
C
C UPROG                 : ERRR
C
C DEFINE-DIREKTIVEN     : KEINE
C
C        03.11.88 (HW)  : ORIGINAL
C
C*STARLET***************************************************************
C
      REAL           FVT (0:IDIMF)
C
C                                 TEST, OB  N  WIRKLICH UNGERADZAHLIG
C                                 IST
C
      NHALF  = N / 2
      IF((N - 2*NHALF) .NE. 1) THEN
         WRITE (6,6010)
         WRITE (6,6020) N
         CALL ERRR (101,' EXPAN1   ')
      END IF
C
      IF(N .LT. 1) THEN
         WRITE (6,6010)
         WRITE (6,6025) N
         CALL ERRR (501,' EXPAN1   ')
      END IF
C
      NLARGE = N + 1
C
      IF(NLARGE .GT. IDIMF) THEN
         WRITE (6,6010)
         WRITE (6,6030) IDIMF, NLARGE
         CALL ERRR (502,' EXPAN1   ')
      END IF
C
C                                 STRIKTE ERFUELLUNG DER PERIODIZITAET
C                                 FVT(0) =! FVT(N)
C
      FVT(0) = 0.5 * (FVT(0) + FVT(N))
      FVT(N) = FVT(0)
C
      FVT0   = FVT(0)
      DELT   = DELTO * FLOAT (N) / FLOAT (NLARGE)
C
C                                 ABARBEITUNG DER 0. STUETZSTELLE
C
      TINT   = (FLOAT(0) + 0.5) * DELT
      TOP    = (FLOAT(0) + 0.5) * DELTO
      TOM    = TOP - DELTO
      FVTI   = (FVT0 - FVT(N-1)) / (TOP - TOM) * (TINT - TOM) + FVT(N-1)
      FVTIO  = FVTI
C
C                                 ABARBEITUNG DER STUETZSTELLEN
C                                 1,2,3, ... N-1
C
      IF(N .GE. 2) THEN
         DO 100 I = 1,N-1
            TINT     = (FLOAT(I) + 0.5) * DELT
            TOP      = (FLOAT(I) + 0.5) * DELTO
            TOM      = TOP - DELTO
            FVTI     = (FVT(I) - FVT(I-1)) / (TOP - TOM) * (TINT - TOM)
     $               +  FVT(I-1)
            FVT(I-1) = FVTIO
            FVTIO    = FVTI
  100    CONTINUE
      END IF
C
C                                 ABARBEITUNG DER (N). STUETZSTELLE
C
      TINT     = (FLOAT(N) + 0.5) * DELT
      TOP      = (FLOAT(N) + 0.5) * DELTO
      TOM      = TOP - DELTO
      FVTI     = (FVT0 - FVT(N-1)) / (TOP - TOM) * (TINT - TOM)
     $         +  FVT(N-1)
      FVT(N-1) =  FVTIO
      FVT(N)   =  FVTI
C
      FVT(N+1) =  FVT(0)
C
      RETURN
 6010 FORMAT (/,' ********** FEHLERMELDUNG AUS SUBR. EXPAN1',
     $        ' **********')
 6020 FORMAT (' N = ',I6,' IST GERADZAHLIG, SOLLTE ABER EIGENTLICH',
     $        ' UNGERADZAHLIG SEIN !')
 6025 FORMAT (' N = ',I6,' IST KLEINER ALS 1 ??? !!!')
 6030 FORMAT (' FELDUEBERSCHREITUNG !! ISTWERT : FVT (',I6,
     $        ')    SOLLWERT : FVT (',I6,')')
      END
