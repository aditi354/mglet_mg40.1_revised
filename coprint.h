
      COMMON /COPRINT/
     &                IPRGRID,    PRINTFORMAT,
     &                ISELPREC,
#ifdef _TSCAL_
     &                LLU,  LLV,  LLW,  LLT, LLP,  LLG,  LLB,
#else
     &                LLU,  LLV,  LLW,  LLP,  LLG,  LLB,
#endif
     &                NPRINTPOS,
     &                ISTPR,   JSTPR,   KSTPR,
     &                PRINTEBE,   IPRINTEBE,
     &                IP1,  IP2,  IPS,
     &                JP1,  JP2,  JPS,
     &                KP1,  KP2,  KPS,
     &                CIP1

      INTEGER
     &        IPRGRID,    ISELPRINT,
     &        LLU,  LLV,  LLW,  LLP,  LLG,  LLB,
     &        NPRINTPOS,
     &        ISTPR(8),   JSTPR(8),   KSTPR(8),
     &        IPRINTEBE,
     &        IP1,  IP2,  IPS,
     &        JP1,  JP2,  JPS,
     &        KP1,  KP2,  KPS,
     &        CIP1

      CHARACTER (LEN=16) PRINTFORMAT
      CHARACTER (LEN=8)  PRINTEBE

