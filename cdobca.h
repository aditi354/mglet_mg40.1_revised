C
C                                 START- UND STOP-KOORDINATEN DER
C                                 LINIEN- BZW. FLAECHENMITTELUNG
C
#ifdef _XHOMOG_
      XSLM   = X(ISTALM)
      XELM   = X(ISTPLM)
C                                 FUER STAGGERED VARIABLE:
      XSLMST = XSLM + 0.5*DX(ISTALM)
      XELMST = XELM + 0.5*DX(ISTPLM)
#else
      XSLM   = 0.0
      XELM   = 0.0
C                                 FUER STAGGERED VARIABLE:
      XSLMST = 0.0
      XELMST = 0.0
#endif
#ifdef _YHOMOG_
      YSLM   = Y(JSTALM)
      YELM   = Y(JSTPLM)
      YSLMST = YSLM + 0.5*DY(JSTALM)
      YELMST = YELM + 0.5*DY(JSTPLM)
#else
      YSLM   = 0.0
      YELM   = 0.0
      YSLMST = 0.0
      YELMST = 0.0
#endif
C
C                                 IN Z-RICHTUNG IST KEINE LINIENMITTE-
C                                 LUNG VORGESEHEN, DAHER WIRD FORMAL EIN
C                                 GEFUEHRT:
C
      ZSLM   = 0.0
      ZELM   = 0.0
      ZSLMST = 0.0
      ZELMST = 0.0
C
C                                 PHYSIKALISCHE ZEIT ZU DER DIE LETZTE
C                                 STICHPROBE FUER DIE BILDUNG VON
C                                 ENSEMBLE-MITTELWERTEN GENOMMEN WURDE
C
      TEEM   = FLOAT(ITTOT - MOD(ITTOT,ITFLUC)) * DT
C
C                                 STICHPROBEN F. D. BILD. V. ENSEMBLE-M.
C                                 WERDEN JEWEILS NACH ABLAUF FOLGENDER
C                                 PHYSIKAL. ZEIT ENTNOMMEN:
C
      DTEM   = FLOAT(ITFLUC) * DT
