
      PARAMETER ( NBSP  =  8 )
      
      COMMON /COGRDDEF/
     &                 NBX,      NBY,      NBZ,
     &                 NTXB,     NTYB,     NTZB,
     &                 DLX,      DLY,      DLZ,
     &                 NLXB,     NLYB,     NLZB,
     &                 SX,       SY,       SZ,
     &                 NRXB,     NRYB,     NRZB,
     &                 DRX,      DRY,      DRZ

      INTEGER
     &       NBX(MAXGRIDS),       NBY(MAXGRIDS),       NBZ(MAXGRIDS),
     & NTXB(NBSP,MAXGRIDS), NTYB(NBSP,MAXGRIDS), NTZB(NBSP,MAXGRIDS),
     & NLXB(NBSP,MAXGRIDS), NLYB(NBSP,MAXGRIDS), NLZB(NBSP,MAXGRIDS),
     & NRXB(NBSP,MAXGRIDS), NRYB(NBSP,MAXGRIDS), NRZB(NBSP,MAXGRIDS)

      REAL
     &  DLX(NBSP,MAXGRIDS),  DLY(NBSP,MAXGRIDS),  DLZ(NBSP,MAXGRIDS),
     &   SX(NBSP,MAXGRIDS),   SY(NBSP,MAXGRIDS),   SZ(NBSP,MAXGRIDS),
     &  DRX(NBSP,MAXGRIDS),  DRY(NBSP,MAXGRIDS),  DRZ(NBSP,MAXGRIDS)

