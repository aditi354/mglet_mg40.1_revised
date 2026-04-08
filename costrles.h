
      COMMON /COSTRLES/
     $                  VERS,   NRRUN,
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,    LGRIDREMOVE,
     $                  NPRNEU, FPRNEU, MTSTEP,
     $                  DT,     ITPRIN, IPINF,  ITINT,
     $                  IPRINT_WSS, IPRINT_WNS, ITFLUC, ITMIT,  
     $                  MPCORR, EPCORR, EPFAK,  MPCVOR, MPCNACH,
     $                  IVPINF, 
     $                  OMG,    LDIMLO, ISETRE,
     $                  LINPRN, IWRB,   LREC
#ifdef _TSCAL_
     $                 ,LTXSTART,TXSTART
#endif
#ifdef _SCALGRIDDIAGNOSE_
     $                 ,LLIMITSPDIAG,SPXSTART,SPXSTOP,SPYSTART,SPYSTOP,
     $                  SPZSTART,SPZSTOP
#endif

      INTEGER 
     $                  NRRUN,
     $                  MTURB,  NPRNEU, MTSTEP,
     $                          ITPRIN, IPINF,  ITINT,
     $                  ITFLUC, ITMIT,  MPCORR,        
     $                                          ISETRE,
     $                          MSLIN,  IWRB          

      LOGICAL
     $                  DREAD,  DWRITE, DCONT,  DCUB,   
     $                  LBODYINSERT,
     $                  LGRIDNEW,       LGRIDREMOVE,
     $                  LDIMLO, LINPRN, LREC
#ifdef _TSCAL_
     $                 ,LTXSTART
#endif 
#ifdef _SCALGRIDDIAGNOSE_
     $                 ,LLIMITSPDIAG
#endif    

      REAL
     $                  DT,     EPCORR,  OMG, FPRNEU
#ifdef _TSCAL_
     $                 ,TXSTART
#endif
#ifdef _SCALGRIDDIAGNOSE_
     $               ,SPXSTART,SPXSTOP,SPYSTART,SPYSTOP,SPZSTART,SPZSTOP
#endif     

      CHARACTER (LEN=8)       VERS

