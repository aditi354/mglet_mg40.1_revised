C
C                                 REZIPROKE BEZUGSGROESSEN
C
      RUREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(75))))
     $       *       SIGN(1.0,RIDENT(75))
      RLREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(76))))
     $       *       SIGN(1.0,RIDENT(76))
      ROMREF =       (AMAX1((10.0*SMALL),ABS(RIDENT(77))))
     $       *       SIGN(1.0,RIDENT(77))
      REREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(78))))
     $       *       SIGN(1.0,RIDENT(78))
      RGREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(79))))
     $       *       SIGN(1.0,RIDENT(79))
      RPREF  = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(80))))
     $       *       SIGN(1.0,RIDENT(80))
      RTAURE = 1.0 / (AMAX1((10.0*SMALL),ABS(RIDENT(81))))
     $       *       SIGN(1.0,RIDENT(81))
      RO2REF = ROMREF**2
      RHEREF = RUREF**2 / RLREF
      RWAVEN = 1.0      / RLREF
      RASDUI = 2.0 * RUREF**2 * RLREF
      RASDOM = 2.0 * RUREF**2 / RLREF
