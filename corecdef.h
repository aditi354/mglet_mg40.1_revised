
         PARAMETER   ( NPREC= 16500000 , NMREC= 1 )

      COMMON /CORECDEF/
     &                 NVREC,    ITREC,    MAXREC,
     &                 NGREC,
     &                 CIDREC,   IVEREC,   IVOREC,
     &                 KANREC,   KANGEO,
     &                 INXREC,   INYREC,   INZREC,
     &                 XUGREC,      XOGREC,
     &                 YUGREC,      YOGREC,
     &                 ZUGREC,      ZOGREC

       INTEGER
     +         NGREC(NMREC)
     +        ,KANREC(NMREC),KANGEO(NMREC)
     +        ,IVEREC(NMREC),IVOREC(NMREC)
     +        ,INXREC(NMREC),INYREC(NMREC),INZREC(NMREC)
     +        ,NTREC(NMREC)


       REAL
     +      XUGREC(NMREC),XOGREC(NMREC)
     +     ,YUGREC(NMREC),YOGREC(NMREC)
     +     ,ZUGREC(NMREC),ZOGREC(NMREC)
C
       CHARACTER (LEN=16) CIDREC(NMREC),CIDRE2
C
