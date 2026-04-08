










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
      SUBROUTINE ITINF(ITTOT,TIME,IPCORR,DIVGMX,IGRID,CPSEC
     $                ,MTURB,EPSU,EPSV,EPSW,ESUMG,ESUMS
     $                ,IPRINT_WSS,WSSX,WSSY,WSSZ
     $                ,IPRINT_WNS,WNSX,WNSY,WNSZ
     $                ,UBLK,GRDPX
     $                ,WALLSSX,WALLSSY,WALLSSZ 
     $                ,IDIVMAX,JDIVMAX,KDIVMAX
     $                ,XDIVMAX,YDIVMAX,ZDIVMAX,GRDIVMAX    
     $                )


      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS         =128 )
      PARAMETER ( MAXBOCONDS       = 10 )

      
      COMMON /COGRDPRO/ LEVEL,LCHILD,XMIN,YMIN,ZMIN,XTOT,YTOT,ZTOT
      COMMON /COGRDPRO/ XHOMOG,YHOMOG,ZHOMOG,NXGRAE,NYGRAE,NZGRAE
      COMMON /COGRDPRO/ GRADPX,UBULKX,LTST,LVP,LSCAI,LPLEVEL,LPOISSONDIR
      COMMON /COGRDPRO/ LSLICE,NXSLICE,NYSLICE,NZSLICE,NVPGRIDS
      COMMON /COGRDPRO/ CONV1SANF,CONV1SEND,TRANSLES1,TRANSLES2
      COMMON /COGRDPRO/ GRADPXOLD

      INTEGER LEVEL(MAXGRIDS),NVPGRIDS(MAXGRIDS)
      INTEGER NXGRAE(MAXGRIDS),NYGRAE(MAXGRIDS),NZGRAE(MAXGRIDS)
      INTEGER NXSLICE(MAXGRIDS),NYSLICE(MAXGRIDS),NZSLICE(MAXGRIDS)

      REAL XTOT(MAXGRIDS),YTOT(MAXGRIDS),ZTOT(MAXGRIDS)
      REAL XMIN(MAXGRIDS),YMIN(MAXGRIDS),ZMIN(MAXGRIDS)
      REAL GRADPX(MAXGRIDS),UBULKX(MAXGRIDS)
      REAL CONV1SANF(MAXGRIDS),CONV1SEND(MAXGRIDS)
      REAL TRANSLES1(MAXGRIDS),TRANSLES2(MAXGRIDS),GRADPXOLD(MAXGRIDS)

      LOGICAL LCHILD(MAXGRIDS),LSLICE(MAXGRIDS),LPOISSONDIR(MAXGRIDS)
      LOGICAL XHOMOG(MAXGRIDS),YHOMOG(MAXGRIDS),ZHOMOG(MAXGRIDS)
      LOGICAL LTST(MAXGRIDS),LVP(MAXGRIDS),LSCAI(MAXGRIDS)
      LOGICAL LPLEVEL(MAXGRIDS)
     

C***********************************************************************
C  I T I N F     ITERATIONSINFORMATION
C*************************************************** F.BAETKE 16.10.81 *
C                                         GEAENDERT: J.FERSTL 17.11.82
C                                                  : F.BAETKE 04.10.84 *
C                                                  : M.MANHART22. 4.92 *
C                                                  : T.BRUNNER 07.05.03
C
C  PARAMETER  ITTOT     - DURCHGEFUEHRTER ZEITSCHRITT
C             TIME      - ERREICHTE DIMENSIONSLOSE ZEIT
C             IPCORR    - ZAHL DER DURCHGEFUEHRTEN PCORR-ZYKLEN
C             DIVGMX    - MAXIMALE DIVERGENZ NACH KORREKTUR
C             CPSEC     - CP-ZEIT BEI LETZTER ABFRAGE
C                       + AKTUELLE CP-ZEIT
C             MTURB     - =0 ES WIRD INFORMATION UEBER EPSU.. AUSGEG.
C                          1 KEINE INFORMATION UEBER EPSU.. AUSGEG.
C             EPSU      - NORMIERTE ABWEICHUNG U-U0
C             EPSV      - NORMIERTE ABWEICHUNG V-V0
C             EPSW      - NORMIERTE ABWEICHUNG W-W0
C             WSSX,Y,Z- WALL-SHEAR-STRESS IN X-, Y-, Z-DIRECTION
C             WNSX,Y,Z- WALL-NORMAL-STRESS IN X-, Y-, Z-DIRECTION
C
C  VERS: 04.02.86 (HW)  : FORMATE Z.T. GEAENDERT
C        07.02.86 (FB)  : JETZT AUF GETSEC UMGESTELLT
C        07.05.03 (TB)  : SCALAR OUTPUT IMPLEMENTED
C
C***********************************************************************
C
      REAL WALLSSX(6), WALLSSY(6), WALLSSZ(6)
C
      INTEGER GRDIVMAX
C      

      CPSECA = GETSEC(0)
      CPDIFF = CPSECA - CPSEC

      WRITE(6,6000) IGRID,ITTOT,TIME,IPCORR,DIVGMX,CPDIFF
C
CTBC  140703: CHANGED FOR MORE DETAILED DIV OUTPUT
      IF (LSLICE(IGRID)) THEN
        WRITE(6,6014) GRDIVMAX,XDIVMAX,YDIVMAX,ZDIVMAX,
     $                         IDIVMAX,JDIVMAX,KDIVMAX 	
      ELSE
        WRITE(6,6015) IGRID,XDIVMAX,YDIVMAX,ZDIVMAX,
     $                      IDIVMAX,JDIVMAX,KDIVMAX 		
      ENDIF            
C
      IF(MTURB .EQ. 0) THEN
           WRITE(6,6001) EPSU,EPSV,EPSW
      ELSE
           WRITE(6,6002) IGRID,ITTOT,TIME,ESUMG,0.0
      ENDIF
      IF (IPRINT_WSS .EQ. 1) THEN
         WRITE (6,6003) IGRID,ITTOT,TIME,WSSX,WSSY,WSSZ
      ENDIF
      IF (IPRINT_WNS .EQ. 1) THEN
         WRITE (6,6004) IGRID,ITTOT,TIME,WNSX,WNSY,WNSZ
      ENDIF

      
 6000 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' IPC= ',I3,' DIV=',E14.6,
     $       ' CPS= ',F8.2)

 6001 FORMAT(1X,'* EPSU = ',G13.6,' EPSV = ',G13.6,
     $       '  EPSW = ',G13.6)

 6002 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' ESUMG=',E14.6,' ESUMS=',E14.6)
 6102 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' ESUMG=',E20.12)

 6003 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' WSSX=',E12.5,' WSSY=',E12.5,' WSSZ=',E12.5)

 6004 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' WNSX=',E14.6,' WNSY=',E14.6,' WNSZ=',E14.6)
 6005 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' TOPWALLX=',E12.5,' BOTWALLX=',E12.5)
 6035 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' QBWALLX= ',E12.5,' QTWALLX= ',E12.5)
 6006 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' TOPWALLY=',E12.5,' BOTWALLY=',E12.5)
 6007 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' RGTWALLX=',E12.5,' LFTWALLX=',E12.5)
 6027 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' RGTWALLZ=',E12.5,' LFTWALLZ=',E12.5)
 6028 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' FROWALLZ=',E12.5,' BACWALLZ=',E12.5)
 6008 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       ' UBULK= ',E12.5,' GRADPX= ',E12.5)

 6009 FORMAT(1X,'GRID:',I3,' IT= ',I8,' TI= ',F9.3,
     $       '     T=',F14.6,'    DT=',E14.6)
 6010 FORMAT(1X,'SGRD:',I3,'(X,Y,Z): (',F11.5,',',F11.5,',',F11.5,')'
     $       ,' (IJK):(',I4,',',I4,',',I4,')') 
 6011 FORMAT(1X,'GRID:',I3,'(X,Y,Z): (',F11.5,',',F11.5,',',F11.5,')'
     $       ,' (IJK):(',I4,',',I4,',',I4,')') 
 6012 FORMAT(1X,'GRID:',I3,'(X,Y,Z): (',F12.6,',',F12.6,',',F12.6,')',
     $       'PRTMAX=',E13.6) 
 6013 FORMAT(1X,'GRID:',I3,'(X,Y,Z): (',F12.6,',',F12.6,',',F12.6,')',
     $       'PRTMIN=',E13.6)                  
C
 6014 FORMAT(1X,'SGRD:',I3,'(X,Y,Z): (',F11.5,',',F11.5,',',F11.5,')'
     $       ,' (IJK):(',I4,',',I4,',',I4,')')
 6015 FORMAT(1X,'GRID:',I3,'(X,Y,Z): (',F11.5,',',F11.5,',',F11.5,')'
     $       ,' (IJK):(',I4,',',I4,',',I4,')')
 6016 FORMAT(1X,'QC= ',E14.6,'                     QD= ',E14.6,
     $          '  QDC=',E14.6)
 6017 FORMAT(1X,'QCE=',E14.6,' QCW=',E14.6,'  QDE=',E14.6,'  QDW=',
     $                  E14.6)                     
 6018 FORMAT(1X,'QCN=',E14.6,' QCS=',E14.6,'  QDN=',E14.6,'  QDS=',
     $                  E14.6)
 6019 FORMAT(1X,'QCT=',E14.6,' QCB=',E14.6,'  QDT=',E14.6,'  QDB=',
     $                  E14.6)   
 6020 FORMAT(1X,'GRID:',I3,' DDS: ',E14.6,'  <=  SM: ',E14.6,
     $       '     SPLUS=',F14.6)
 6021 FORMAT(1X,'GRID:',I3,' DDS: ',E14.6,'  >  SM: ',E14.6,
     $       '     SPLUS=',F14.6)            
 6024 FORMAT(1X,'----------------------------------------------------')
 6025 FORMAT(1X,'====================================================')
      RETURN
      END
      FUNCTION GETSEC(NULL)
C***********************************************************************
C   G E T S E C    GETS THE ACCUMULATED CPU-TIME IN SECONDS
C*************************************************** F.BAETKE 21.03.85 *
C
C     THIS FUNCTION SHOULD CONTAIN THE ONLY NON-ANSI CALL
C     IN THE PROGRAM  A B O X 7 7  AND HAS TO BE ADAPTED
C     ACCORDING TO THE HOST-COMPUTER.
C
C  DEFINE-DIREKTIVEN    : B7800,   CRAY,    CYBER
C
C  VERS: 20.02.89 (HW)  : DEFINE-DIREKTIVEN B7800, CRAY, CYBER
C                         EINGEFUEHRT.
C***********************************************************************
C
      DATA  SECINI /10000./
C
C                                 CPU-TIME FOR LINUX (G77)
       CALL CPU_TIME(GETSEC)
C      GETSEC = 0.0
C                                 CPU-TIME FOR  SIEMENS 7800
*     CALL CLOCKM(ICPII)
*     GETSEC = FLOAT(ICPII)/1000.
C                                  CPU-TIME FOR MVS (IBM)
*     GETSEC = SECINI - ITIME(0)/100.
C
      RETURN
      END
      FUNCTION GETDAT(NULL)
C***********************************************************************
C   G E T D A T    GETS THE CURRENT DATE (CHARACTER (LEN=8))
C*************************************************** H.WERNER 17.12.86 *
C
C***********************************************************************
C
C
C
C
C
C
C                                  WHEN USING  LINUX (IFC)
C
      CHARACTER (LEN=8)  GETDAT
      CALL DATE_AND_TIME(GETDAT)

C
      RETURN
      END
      FUNCTION GETAET(NULL)
C***********************************************************************
C   G E T A E T    GETS THE A_CTUAL E_LAPSED T_IME (CHARACTER * 8)
C***********************************************************************
C
C     THIS FUNCTION CONTAINS A NON-ANSI CALL
C     IN THE PROGRAM  S T A R L E T  AND HAS TO BE ADAPTED
C     ACCORDING TO THE HOST-COMPUTER.
C
C  DEFINE-DIREKTIVEN    : CRAY,    CYBER
C
C  VERS: 21.02.89 (HW)  : ORIGINAL
C
C***********************************************************************
C
C
C
C                      AUF LINUX (G77)
       CHARACTER (LEN=8) GETAET
       GETAET = '        '
C
C
      RETURN
      END
