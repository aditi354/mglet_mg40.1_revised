
      INTEGER  ISELEA (2,752)      , ISELEP (2,752),
     $         ISAMPA (752)        , ISAMPP (752),
     $         ISLINP (ISLIDI)     , ISLINA (ISLIDI)
      REAL     RIDENT (100)
      REAL     HILF   (KK ,JJ ,II ), B      (KK ,JJ ,II )
C
      REAL     U      (KK ,JJ ,II ), UO     (KK ,JJ ,II ),
     $         AU     (KKA,JJA,IIA), SU     (KKA,JJA,IIA),
     $         UFG    (KK ,JJ ,II ),
     $         AURG   (KKA,JJA,IIA), SURG   (KKA,JJA,IIA),
     $         AUSKG  (KKA,JJA,IIA), SUSKG  (KKA,JJA,IIA),
     $         AUFLG  (KKA,JJA,IIA), SUFLG  (KKA,JJA,IIA)
      REAL     V      (KK ,JJ ,II ), VO     (KK ,JJ ,II ),
     $         AV     (KKA,JJA,IIA), SV     (KKA,JJA,IIA),
     $         VFG    (KK ,JJ ,II ),
     $         AVRG   (KKA,JJA,IIA), SVRG   (KKA,JJA,IIA),
     $         AVSKG  (KKA,JJA,IIA), SVSKG  (KKA,JJA,IIA),
     $         AVFLG  (KKA,JJA,IIA), SVFLG  (KKA,JJA,IIA)
      REAL     W      (KK ,JJ ,II ), WO     (KK ,JJ ,II ),
     $         AW     (KKA,JJA,IIA), SW     (KKA,JJA,IIA),
     $         WFG    (KK ,JJ ,II ),
     $         AWRG   (KKA,JJA,IIA), SWRG   (KKA,JJA,IIA),
     $         AWSKG  (KKA,JJA,IIA), SWSKG  (KKA,JJA,IIA),
     $         AWFLG  (KKA,JJA,IIA), SWFLG  (KKA,JJA,IIA)
      REAL     P      (KK ,JJ ,II ),
     $         AP     (KKA,JJA,IIA), SP     (KKA,JJA,IIA),
     $         PFG    (KK ,JJ ,II ),
     $         APRG   (KKA,JJA,IIA), SPRG   (KKA,JJA,IIA),
     $         APSKG  (KKA,JJA,IIA), SPSKG  (KKA,JJA,IIA),
     $         APFLG  (KKA,JJA,IIA), SPFLG  (KKA,JJA,IIA)
      REAL     G      (KK ,JJ ,II ),
     $         AEFG   (KKA,JJA,IIA), SEFG   (KKA,JJA,IIA),
     $         AEFS   (KKA,JJA,IIA), SEFS   (KKA,JJA,IIA),
     $         ADFG   (KKA,JJA,IIA), SDFG   (KKA,JJA,IIA),
     $         ADUDX2 (KKA,JJA,IIA), SDUDX2 (KKA,JJA,IIA),
     $         ADUDY2 (KKA,JJA,IIA), SDUDY2 (KKA,JJA,IIA),
     $         ADUDZ2 (KKA,JJA,IIA), SDUDZ2 (KKA,JJA,IIA),
     $         ADVDX2 (KKA,JJA,IIA), SDVDX2 (KKA,JJA,IIA),
     $         ADVDY2 (KKA,JJA,IIA), SDVDY2 (KKA,JJA,IIA),
     $         ADVDZ2 (KKA,JJA,IIA), SDVDZ2 (KKA,JJA,IIA),
     $         ADWDX2 (KKA,JJA,IIA), SDWDX2 (KKA,JJA,IIA),
     $         ADWDY2 (KKA,JJA,IIA), SDWDY2 (KKA,JJA,IIA),
     $         ADWDZ2 (KKA,JJA,IIA), SDWDZ2 (KKA,JJA,IIA),
     $         AUFWFG (KKA,JJA,IIA), SUFWFG (KKA,JJA,IIA),
     $         AUFWFS (KKA,JJA,IIA), SUFWFS (KKA,JJA,IIA),
     $         AUFWFM (KKA,JJA,IIA), SUFWFM (KKA,JJA,IIA)
      REAL     AVFWFG (KKA,JJA,IIA), SVFWFG (KKA,JJA,IIA),
     $         AVFWFS (KKA,JJA,IIA), SVFWFS (KKA,JJA,IIA),
     $         AVFWFM (KKA,JJA,IIA), SVFWFM (KKA,JJA,IIA),
     $         AUFVFG (KKA,JJA,IIA), SUFVFG (KKA,JJA,IIA),
     $         AUFVFS (KKA,JJA,IIA), SUFVFS (KKA,JJA,IIA),
     $         AUFVFM (KKA,JJA,IIA), SUFVFM (KKA,JJA,IIA)
      REAL     OX     (KK ,JJ ,II ), OXFG   (KK ,JJ ,II ),
     $         AOX    (KKA,JJA,IIA), SOX    (KKA,JJA,IIA),
     $         AOXRG  (KKA,JJA,IIA), SOXRG  (KKA,JJA,IIA)
      REAL     OY     (KK ,JJ ,II ), OYFG   (KK ,JJ ,II ),
     $         AOY    (KKA,JJA,IIA), SOY    (KKA,JJA,IIA),
     $         AOYRG  (KKA,JJA,IIA), SOYRG  (KKA,JJA,IIA)
      REAL     OZ     (KK ,JJ ,II ), OZFG   (KK ,JJ ,II ),
     $         AOZ    (KKA,JJA,IIA), SOZ    (KKA,JJA,IIA),
     $         AOZRG  (KKA,JJA,IIA), SOZRG  (KKA,JJA,IIA)
      REAL                           O2FG   (KK ,JJ ,II ),
     $         AO2    (KKA,JJA,IIA), SO2    (KKA,JJA,IIA),
     $         AO2RG  (KKA,JJA,IIA), SO2RG  (KKA,JJA,IIA)
      REAL                           HEFG   (KK ,JJ ,II ),
     $         AHE    (KKA,JJA,IIA), SHE    (KKA,JJA,IIA),
     $         AHERG  (KKA,JJA,IIA), SHERG  (KKA,JJA,IIA)
#ifdef _TAU_NN_
C                             NICHT-NEWTONSCHER SPANNUNGSTENSOR
      REAL       TAU11(KK ,JJ ,II ),   TAU12(KK ,JJ ,II ), 
     $           TAU13(KK ,JJ ,II ),   TAU22(KK ,JJ ,II ),
     $           TAU23(KK ,JJ ,II ),   TAU33(KK ,JJ ,II ),
     $         ATAU11(KKA,JJA,IIA), STAU11(KKA,JJA,IIA),
     $         ATAU12(KKA,JJA,IIA), STAU12(KKA,JJA,IIA),
     $         ATAU13(KKA,JJA,IIA), STAU13(KKA,JJA,IIA),
     $         ATAU22(KKA,JJA,IIA), STAU22(KKA,JJA,IIA),
     $         ATAU23(KKA,JJA,IIA), STAU23(KKA,JJA,IIA),
     $         ATAU33(KKA,JJA,IIA), STAU33(KKA,JJA,IIA)
#endif
#if (defined _FRED_BODY_)
      REAL BP(KK,JJ,II),BU(KK,JJ,II),BV(KK,JJ,II),BW(KK,JJ,II)
#endif
#ifdef _TSCAL_     
      REAL     T      (KK ,JJ ,II ), TO     (KK ,JJ ,II ),
     $         AT     (KKA,JJA,IIA), ST     (KKA,JJA,IIA),
     $         TFG    (KK ,JJ ,II ),
     $         ATRG   (KKA,JJA,IIA), STRG   (KKA,JJA,IIA), 
     $         ATSKG  (KKA,JJA,IIA), STSKG  (KKA,JJA,IIA),
     $         ATFLG  (KKA,JJA,IIA), STFLG  (KKA,JJA,IIA),         
     $         AUFTFG (KKA,JJA,IIA), SUFTFG (KKA,JJA,IIA),
     $         AUFTFS (KKA,JJA,IIA), SUFTFS (KKA,JJA,IIA),
     $         AUFTFM (KKA,JJA,IIA), SUFTFM (KKA,JJA,IIA),
     $         AVFTFG (KKA,JJA,IIA), SVFTFG (KKA,JJA,IIA),
     $         AVFTFS (KKA,JJA,IIA), SVFTFS (KKA,JJA,IIA),
     $         AVFTFM (KKA,JJA,IIA), SVFTFM (KKA,JJA,IIA),  
     $         AWFTFG (KKA,JJA,IIA), SWFTFG (KKA,JJA,IIA),
     $         AWFTFS (KKA,JJA,IIA), SWFTFS (KKA,JJA,IIA),
     $         AWFTFM (KKA,JJA,IIA), SWFTFM (KKA,JJA,IIA)        
#endif
#if defined _NEWSTAT_
      REAL     HILF3D1(KK ,JJ ,II ),HILF3D2(KK ,JJ ,II ),
     $         HILF3D3(KK ,JJ ,II )
#endif
#if (defined _STAT11_)
      REAL       AUUM(KKA,JJA,IIA),    SUUM(KKA,JJA,IIA),
     $           AVVM(KKA,JJA,IIA),    SVVM(KKA,JJA,IIA),
     $           AWWM(KKA,JJA,IIA),    SWWM(KKA,JJA,IIA),
     $           APPM(KKA,JJA,IIA),    SPPM(KKA,JJA,IIA)
#endif
#if (defined _STAT12_)
      REAL       AUVM(KKA,JJA,IIA),    SUVM(KKA,JJA,IIA),
     $           AUWM(KKA,JJA,IIA),    SUWM(KKA,JJA,IIA),
     $           AVWM(KKA,JJA,IIA),    SVWM(KKA,JJA,IIA)
#endif
#if (defined _STAT13_)
      REAL     AUXUXM(KKA,JJA,IIA),  SUXUXM(KKA,JJA,IIA),
     $         AUYUYM(KKA,JJA,IIA),  SUYUYM(KKA,JJA,IIA),
     $         AUZUZM(KKA,JJA,IIA),  SUZUZM(KKA,JJA,IIA),
     $         AVXVXM(KKA,JJA,IIA),  SVXVXM(KKA,JJA,IIA),
     $         AVYVYM(KKA,JJA,IIA),  SVYVYM(KKA,JJA,IIA),
     $         AVZVZM(KKA,JJA,IIA),  SVZVZM(KKA,JJA,IIA),
     $         AWXWXM(KKA,JJA,IIA),  SWXWXM(KKA,JJA,IIA),
     $         AWYWYM(KKA,JJA,IIA),  SWYWYM(KKA,JJA,IIA),
     $         AWZWZM(KKA,JJA,IIA),  SWZWZM(KKA,JJA,IIA)
#endif
#if (defined _STAT14_)
      REAL      AUUUM(KKA,JJA,IIA),  SUUUM(KKA,JJA,IIA),
     $          AVVVM(KKA,JJA,IIA),  SVVVM(KKA,JJA,IIA),
     $          AWWWM(KKA,JJA,IIA),  SWWWM(KKA,JJA,IIA),
     $          AUUVM(KKA,JJA,IIA),  SUUVM(KKA,JJA,IIA),
     $          AUUWM(KKA,JJA,IIA),  SUUWM(KKA,JJA,IIA),
     $          AVVUM(KKA,JJA,IIA),  SVVUM(KKA,JJA,IIA),
     $          AVVWM(KKA,JJA,IIA),  SVVWM(KKA,JJA,IIA),
     $          AWWUM(KKA,JJA,IIA),  SWWUM(KKA,JJA,IIA),
     $          AWWVM(KKA,JJA,IIA),  SWWVM(KKA,JJA,IIA),
     $          AUVWM(KKA,JJA,IIA),  SUVWM(KKA,JJA,IIA),
     $           AUPM(KKA,JJA,IIA),   SUPM(KKA,JJA,IIA),
     $           AVPM(KKA,JJA,IIA),   SVPM(KKA,JJA,IIA),
     $           AWPM(KKA,JJA,IIA),   SWPM(KKA,JJA,IIA),
     $          AUXPM(KKA,JJA,IIA),  SUXPM(KKA,JJA,IIA),
     $          AVYPM(KKA,JJA,IIA),  SVYPM(KKA,JJA,IIA),
     $          AWZPM(KKA,JJA,IIA),  SWZPM(KKA,JJA,IIA),
     $         AUYVXP(KKA,JJA,IIA), SUYVXP(KKA,JJA,IIA),
     $         AUZWXP(KKA,JJA,IIA), SUZWXP(KKA,JJA,IIA),
     $         AVZWYP(KKA,JJA,IIA), SVZWYP(KKA,JJA,IIA),
     $         AUXVXM(KKA,JJA,IIA), SUXVXM(KKA,JJA,IIA),
     $         AUYVYM(KKA,JJA,IIA), SUYVYM(KKA,JJA,IIA),
     $         AUZVZM(KKA,JJA,IIA), SUZVZM(KKA,JJA,IIA),
     $         AUXWXM(KKA,JJA,IIA), SUXWXM(KKA,JJA,IIA),
     $         AUYWYM(KKA,JJA,IIA), SUYWYM(KKA,JJA,IIA),
     $         AUZWZM(KKA,JJA,IIA), SUZWZM(KKA,JJA,IIA),
     $         AVXWXM(KKA,JJA,IIA), SVXWXM(KKA,JJA,IIA),
     $         AVYWYM(KKA,JJA,IIA), SVYWYM(KKA,JJA,IIA),
     $         AVZWZM(KKA,JJA,IIA), SVZWZM(KKA,JJA,IIA)
#endif
#if (defined _STAT15_)
      REAL  AUTM(KKA,JJA,IIA),SUTM(KKA,JJA,IIA)
     $    ,AVTM(KKA,JJA,IIA),SVTM(KKA,JJA,IIA)
     $    ,AWTM(KKA,JJA,IIA),SWTM(KKA,JJA,IIA)
     $    ,ATTM(KKA,JJA,IIA),STTM(KKA,JJA,IIA)
     $    ,ASSM(KKA,JJA,IIA),SSSM(KKA,JJA,IIA)
#endif
#if (defined _STAT17_)
      REAL AUTTM(KKA,JJA,IIA),SUTTM(KKA,JJA,IIA),
     $     AVTTM(KKA,JJA,IIA),SVTTM(KKA,JJA,IIA),
     $     AWTTM(KKA,JJA,IIA),SWTTM(KKA,JJA,IIA),
     $     ATXTXM(KKA,JJA,IIA),STXTXM(KKA,JJA,IIA),
     $     ATYTYM(KKA,JJA,IIA),STYTYM(KKA,JJA,IIA),
     $     ATZTZM(KKA,JJA,IIA),STZTZM(KKA,JJA,IIA)
#endif
#if (defined _STAT20_)
      REAL     AUUUUM (KKA,JJA,IIA), SUUUUM (KKA,JJA,IIA),
     &         AVVVVM (KKA,JJA,IIA), SVVVVM (KKA,JJA,IIA),
     &         AWWWWM (KKA,JJA,IIA), SWWWWM (KKA,JJA,IIA),
     &         APPPPM (KKA,JJA,IIA), SPPPPM (KKA,JJA,IIA)
#ifdef _TSCAL_
      REAL     ATTTTM (KKA,JJA,IIA), STTTTM (KKA,JJA,IIA),
     $         ATTTM (KKA,JJA,IIA), STTTM (KKA,JJA,IIA)
#endif
#endif
