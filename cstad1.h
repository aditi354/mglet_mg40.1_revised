
      REAL     FVT    (0:IDIMF),           HVK    (KLASS,JJ2L)
      INTEGER  ISELEA (2,752),             ISELEP (2,752),
     $         ISLINP (ISLIDI),            ISLINA (ISLIDI),
     $         KSTW   (JJ  ,II  )
      REAL     RIDENT (100)
      REAL     HILF   (KK  ,JJ  ,II     )
C
      REAL     UFG    (KK  ,JJ  ,II     ), AURG   (KKA ,JJA ,IIA    ),
     $         VFG    (KK  ,JJ  ,II     ), AVRG   (KKA ,JJA ,IIA    ),
     $         WFG    (KK  ,JJ  ,II     ), AWRG   (KKA ,JJA ,IIA    ),
     $         OX     (KK  ,JJ  ,II     ), AOX    (KKA ,JJA ,IIA    ),
     $         OXFG   (KK  ,JJ  ,II     ), AOXRG  (KKA ,JJA ,IIA    ),
     $         OY     (KK  ,JJ  ,II     ), AOY    (KKA ,JJA ,IIA    ),
     $         OYFG   (KK  ,JJ  ,II     ), AOYRG  (KKA ,JJA ,IIA    ),
     $         OZ     (KK  ,JJ  ,II     ), AOZ    (KKA ,JJA ,IIA    ),
     $         OZFG   (KK  ,JJ  ,II     ), AOZRG  (KKA ,JJA ,IIA    )
      REAL     ARXUU  (KKXL,JJ2L,ILIMXP ), SRXUU  (KKXL,JJ1L,ILIMXP ),
     $         ARYUU  (KKYL,JJ2L,ILIMXP ), SRYUU  (KKYL,JJ1L,ILIMXP ),
     $         ARZUU  (KKZL,JJ2L,ILIMXP ), SRZUU  (KKZL,JJ1L,ILIMXP ),
     $         ASXUU  (KSXL,JJ1L,ILIMXP ), SSXUU  (KSXL,JJ1L,ILIMXP ),
     $         ASYUU  (KSYL,JJ1L,ILIMXP ), SSYUU  (KSYL,JJ1L,ILIMXP ),
     $         ASZUU  (KSZL,JJ1L,ILIMXP ), SSZUU  (KSZL,JJ1L,ILIMXP )
      REAL     ARXVV  (KKXL,JJ2L,ILIMXP ), SRXVV  (KKXL,JJ1L,ILIMXP ),
     $         ARYVV  (KKYL,JJ2L,ILIMXP ), SRYVV  (KKYL,JJ1L,ILIMXP ),
     $         ARZVV  (KKZL,JJ2L,ILIMXP ), SRZVV  (KKZL,JJ1L,ILIMXP ),
     $         ASXVV  (KSXL,JJ1L,ILIMXP ), SSXVV  (KSXL,JJ1L,ILIMXP ),
     $         ASYVV  (KSYL,JJ1L,ILIMXP ), SSYVV  (KSYL,JJ1L,ILIMXP ),
     $         ASZVV  (KSZL,JJ1L,ILIMXP ), SSZVV  (KSZL,JJ1L,ILIMXP )
      REAL     ARXWW  (KKXL,JJ2L,ILIMXP ), SRXWW  (KKXL,JJ1L,ILIMXP ),
     $         ARYWW  (KKYL,JJ2L,ILIMXP ), SRYWW  (KKYL,JJ1L,ILIMXP ),
     $         ARZWW  (KKZL,JJ2L,ILIMXP ), SRZWW  (KKZL,JJ1L,ILIMXP ),
     $         ASXWW  (KSXL,JJ1L,ILIMXP ), SSXWW  (KSXL,JJ1L,ILIMXP ),
     $         ASYWW  (KSYL,JJ1L,ILIMXP ), SSYWW  (KSYL,JJ1L,ILIMXP ),
     $         ASZWW  (KSZL,JJ1L,ILIMXP ), SSZWW  (KSZL,JJ1L,ILIMXP )
C
C                                  KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------
C
      REAL     ARXUV  (KKXL,JJ2L,ILIMXP ), SRXUV  (KKXL,JJ1L,ILIMXP ),
     $         ARYUV  (KKYL,JJ2L,ILIMXP ), SRYUV  (KKYL,JJ1L,ILIMXP ),
     $         ARXUW  (KKXL,JJ2L,ILIMXP ), SRXUW  (KKXL,JJ1L,ILIMXP ),
     $         ARYUW  (KKYL,JJ2L,ILIMXP ), SRYUW  (KKYL,JJ1L,ILIMXP ),
     $         ARYVW  (KKYL,JJ2L,ILIMXP ), SRYVW  (KKYL,JJ1L,ILIMXP )
C
C                                  KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------
C
      REAL     ACXUW  (KKXL,JJ2L,ILIMXP ), SCXUW  (KKXL,JJ1L,ILIMXP ),
     $         ACZUW  (KKZL,JJ2L,ILIMXP ), SCZUW  (KKZL,JJ1L,ILIMXP )
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
      REAL     ARXOXX (KKXL,JJ2L,ILIMXP ), SRXOXX (KKXL,JJ1L,ILIMXP ),
     $         ARXOXY (KKXL,JJ2L,ILIMXP ), SRXOXY (KKXL,JJ1L,ILIMXP ),
     $         ARXOXZ (KKXL,JJ2L,ILIMXP ), SRXOXZ (KKXL,JJ1L,ILIMXP ),
     $         ARXOYY (KKXL,JJ2L,ILIMXP ), SRXOYY (KKXL,JJ1L,ILIMXP ),
     $         ARXOYZ (KKXL,JJ2L,ILIMXP ), SRXOYZ (KKXL,JJ1L,ILIMXP ),
     $         ARXOZZ (KKXL,JJ2L,ILIMXP ), SRXOZZ (KKXL,JJ1L,ILIMXP )
      REAL     ARYOXX (KKYL,JJ2L,ILIMXP ), SRYOXX (KKYL,JJ1L,ILIMXP ),
     $         ARYOXY (KKYL,JJ2L,ILIMXP ), SRYOXY (KKYL,JJ1L,ILIMXP ),
     $         ARYOXZ (KKYL,JJ2L,ILIMXP ), SRYOXZ (KKYL,JJ1L,ILIMXP ),
     $         ARYOYY (KKYL,JJ2L,ILIMXP ), SRYOYY (KKYL,JJ1L,ILIMXP ),
     $         ARYOYZ (KKYL,JJ2L,ILIMXP ), SRYOYZ (KKYL,JJ1L,ILIMXP ),
     $         ARYOZZ (KKYL,JJ2L,ILIMXP ), SRYOZZ (KKYL,JJ1L,ILIMXP )
      REAL     ARZOXX (KKZL,JJ2L,ILIMXP ), SRZOXX (KKZL,JJ1L,ILIMXP ),
     $         ARZOXY (KKZL,JJ2L,ILIMXP ), SRZOXY (KKZL,JJ1L,ILIMXP ),
     $         ARZOXZ (KKZL,JJ2L,ILIMXP ), SRZOXZ (KKZL,JJ1L,ILIMXP ),
     $         ARZOYY (KKZL,JJ2L,ILIMXP ), SRZOYY (KKZL,JJ1L,ILIMXP ),
     $         ARZOYZ (KKZL,JJ2L,ILIMXP ), SRZOYZ (KKZL,JJ1L,ILIMXP ),
     $         ARZOZZ (KKZL,JJ2L,ILIMXP ), SRZOZZ (KKZL,JJ1L,ILIMXP )
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
      REAL     ASXOXX (KSXL,JJ1L,ILIMXP ), SSXOXX (KSXL,JJ1L,ILIMXP ),
     $         ASXOYY (KSXL,JJ1L,ILIMXP ), SSXOYY (KSXL,JJ1L,ILIMXP ),
     $         ASXOZZ (KSXL,JJ1L,ILIMXP ), SSXOZZ (KSXL,JJ1L,ILIMXP ),
     $         ASYOXX (KSYL,JJ1L,ILIMXP ), SSYOXX (KSYL,JJ1L,ILIMXP ),
     $         ASYOYY (KSYL,JJ1L,ILIMXP ), SSYOYY (KSYL,JJ1L,ILIMXP ),
     $         ASYOZZ (KSYL,JJ1L,ILIMXP ), SSYOZZ (KSYL,JJ1L,ILIMXP )
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
      REAL     AHOZOY (KH0L,JJ2L,ILIMXP ), SHOZOY (KH0L,JJ2L,ILIMXP ),
     $         AHOZOX (KH0L,JJ2L,ILIMXP ), SHOZOX (KH0L,JJ2L,ILIMXP ),
     $         AHOYOX (KH0L,JJ2L,ILIMXP ), SHOYOX (KH0L,JJ2L,ILIMXP )
