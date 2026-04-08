
      PARAMETER (NPPHYS_MAX=1)
      COMMON /COPHYSPAR/
     $                  MTURB,  TU_LEVEL,   RHO,    GMOL,   UGRID,
     $                  VREF,   EXPON,  UTAUX,
     $                  CIDUFR, UFRCON, UFRFREQ, DELTA,  XREF,
     $                  IPP,    JPP,    KPP,
     $                  XPER,   YPER,   ZPER, 
     $                  NXPER,  NYPER,  NZPER,
     $                  NPPHYS, PPHYS, XPPHYS
#ifdef _TSCAL_
     $                  ,LAMDA, PRMOL,  PRTURB, TREF, EXPONT, 
     $                  CIDTFR, TFRCON, DELTAT
#endif

      INTEGER           IPP,  JPP,  KPP,   NPPHYS
     
      REAL              TU_LEVEL,  RHO,   GMOL,   UGRID,  VREF,  
     $                  EXPON, UFRCON, 
     $                  PPHYS(NPPHYS_MAX), XPPHYS(NPPHYS_MAX)
#ifdef _TSCAL_
     $               ,LAMDA, PRMOL, PRTURB, TREF, EXPONT, TFRCON, DELTAT
#endif     

#ifdef _TSCAL_
      CHARACTER (LEN=16)      CIDUFR, CIDTFR
#else
      CHARACTER (LEN=16)      CIDUFR
#endif

