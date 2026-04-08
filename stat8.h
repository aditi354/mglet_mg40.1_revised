C
C           DIMENSIONIERUNGEN FUER DIE STATISTISCHE AUSWERTUNG
C
C                  LEVEL: <U>, <V>, <W> UND <P>
C                         UFG, VFG, WFG, PFG
C                         <U-RMS> G, <V.RMS> G, <W-RMS> G, <P-RMS> G
C                    <U-RMS> G+S, <V-RMS> G+S, <W-RMS> G+S,
C                    <UV-F> G,      <UW-F> G,     <VW-F> G,    <E-F> G
C                    <UV-F> S,      <UW-F> S,     <VW-F> S,    <E-F> S
C                    <UV-F> G+S,    <UW-F> G+S,   <VW-F> G+S,  <E-F> G+S
C                    <UV-F> M,      <UW-F> M,     <VW-F> M,
C                    <UV-F> G+S+M,  <UW-F> G+S+M, <VW-F> G+S+M,
C                    <U-SKE> G, <V-SKE> G, <W-SKE> G, <P-SKE> G
C                    <U-FLA> G, <V-FLA> G, <W-FLA> G, <P-FLA> G
C
#ifdef _TAU_NN_
C                   <TAU11>, <TAU12>, <TAU13>, <TAU22>, <TAU23>, <TAU33>
#endif
C
C
      REAL     AU     (   IDIMA         ), SU     (   IDIMA         ),
***   REAL     AU     (        1        ), SU     (        1        ),
     $         UFG    (   IDIM3D     ),
***  $         UFG    (       1      ),
     $         AURG   (   IDIMA         ), SURG   (   IDIMA         ),
***  $         AURG   (        1        ), SURG   (        1        ),
     $         AUSKG  (   IDIMA         ), SUSKG  (   IDIMA         ),
***  $         AUSKG  (        1        ), SUSKG  (        1        ),
     $         AUFLG  (   IDIMA         ), SUFLG  (   IDIMA         ),
***  $         AUFLG  (        1        ), SUFLG  (        1        ),
*    $         ARXUU  (   IDIM2L        ), SRXUU  (     IDIM1L      ),
     $         ARXUU  (        1        ), SRXUU  (        1        ),
*    $         ARYUU  (   IDIM2L        ), SRYUU  (     IDIM1L      ),
     $         ARYUU  (        1        ), SRYUU  (        1        ),
*    $         ARZUU  (   IDIM2L        ), SRZUU  (     IDIM1L      ),
     $         ARZUU  (        1        ), SRZUU  (        1        ),
*    $         ASXUU  (     IDIM1L      ), SSXUU  (     IDIM1L      ),
     $         ASXUU  (        1        ), SSXUU  (        1        ),
*    $         ASYUU  (     IDIM1L      ), SSYUU  (     IDIM1L      ),
     $         ASYUU  (        1        ), SSYUU  (        1        ),
*    $         ASZUU  (     IDIM1L      ), SSZUU  (     IDIM1L      )
     $         ASZUU  (        1        ), SSZUU  (        1        )
      REAL     AV     (   IDIMA         ), SV     (   IDIMA         ),
***   REAL     AV     (        1        ), SV     (        1        ),
     $         VFG    (  IDIM3D      ),
***  $         VFG    (       1      ),
     $         AVRG   (   IDIMA         ), SVRG   (   IDIMA         ),
***  $         AVRG   (        1        ), SVRG   (        1        ),
     $         AVSKG  (   IDIMA         ), SVSKG  (   IDIMA         ),
***  $         AVSKG  (        1        ), SVSKG  (        1        ),
     $         AVFLG  (   IDIMA         ), SVFLG  (   IDIMA         ),
***  $         AVFLG  (        1        ), SVFLG  (        1        ),
*    $         ARXVV  (   IDIM2L        ), SRXVV  (     IDIM1L      ),
     $         ARXVV  (        1        ), SRXVV  (        1        ),
*    $         ARYVV  (   IDIM2L        ), SRYVV  (     IDIM1L      ),
     $         ARYVV  (        1        ), SRYVV  (        1        ),
*    $         ARZVV  (   IDIM2L        ), SRZVV  (     IDIM1L      ),
     $         ARZVV  (        1        ), SRZVV  (        1        ),
*    $         ASXVV  (     IDIM1L      ), SSXVV  (     IDIM1L      ),
     $         ASXVV  (        1        ), SSXVV  (        1        ),
*    $         ASYVV  (     IDIM1L      ), SSYVV  (     IDIM1L      ),
     $         ASYVV  (        1        ), SSYVV  (        1        ),
*    $         ASZVV  (     IDIM1L      ), SSZVV  (     IDIM1L      )
     $         ASZVV  (        1        ), SSZVV  (        1        )
      REAL     AW     (   IDIMA         ), SW     (   IDIMA         ),
***   REAL     AW     (        1        ), SW     (        1        ),
     $         WFG    (  IDIM3D      ),
***  $         WFG    (       1      ),
     $         AWRG   (   IDIMA         ), SWRG   (   IDIMA         ),
***  $         AWRG   (        1        ), SWRG   (        1        ),
     $         AWSKG  (   IDIMA         ), SWSKG  (   IDIMA         ),
***  $         AWSKG  (        1        ), SWSKG  (        1        ),
     $         AWFLG  (   IDIMA         ), SWFLG  (   IDIMA         ),
***  $         AWFLG  (        1        ), SWFLG  (        1        ),
*    $         ARXWW  (   IDIM2L        ), SRXWW  (     IDIM1L      ),
     $         ARXWW  (        1        ), SRXWW  (        1        ),
*    $         ARYWW  (   IDIM2L        ), SRYWW  (     IDIM1L      ),
     $         ARYWW  (        1        ), SRYWW  (        1        ),
*    $         ARZWW  (   IDIM2L        ), SRZWW  (     IDIM1L      ),
     $         ARZWW  (        1        ), SRZWW  (        1        ),
*    $         ASXWW  (     IDIM1L      ), SSXWW  (     IDIM1L      ),
     $         ASXWW  (        1        ), SSXWW  (        1        ),
*    $         ASYWW  (     IDIM1L      ), SSYWW  (     IDIM1L      ),
     $         ASYWW  (        1        ), SSYWW  (        1        ),
*    $         ASZWW  (     IDIM1L      ), SSZWW  (     IDIM1L      )
     $         ASZWW  (        1        ), SSZWW  (        1        )
      REAL     AP     (   IDIMA         ), SP     (   IDIMA         ),
***   REAL     AP     (        1        ), SP     (        1        ),
     $         PFG    (  IDIM3D      ),
***  $         PFG    (       1      ),
     $         APRG   (   IDIMA         ), SPRG   (   IDIMA         ),
***  $         APRG   (        1        ), SPRG   (        1        ),
     $         APSKG  (   IDIMA         ), SPSKG  (   IDIMA         ),
***  $         APSKG  (        1        ), SPSKG  (        1        ),
     $         APFLG  (   IDIMA         ), SPFLG  (   IDIMA         )
***  $         APFLG  (        1        ), SPFLG  (        1        )
      REAL     AEFG   (   IDIMA         ), SEFG   (   IDIMA         ),
***   REAL     AEFG   (        1        ), SEFG   (        1        ),
     $         AEFS   (   IDIMA         ), SEFS   (   IDIMA         )
***  $         AEFS   (        1        ), SEFS   (        1        )
      REAL     ADFG   (   IDIMA         ), SDFG   (   IDIMA         ),
*     REAL     ADFG   (        1        ), SDFG   (        1        ),
     $         ADUDX2 (   IDIMA         ), SDUDX2 (   IDIMA         ),
*    $         ADUDX2 (        1        ), SDUDX2 (        1        ),
     $         ADUDY2 (   IDIMA         ), SDUDY2 (   IDIMA         ),
*    $         ADUDY2 (        1        ), SDUDY2 (        1        ),
     $         ADUDZ2 (   IDIMA         ), SDUDZ2 (   IDIMA         ),
*    $         ADUDZ2 (        1        ), SDUDZ2 (        1        ),
     $         ADVDX2 (   IDIMA         ), SDVDX2 (   IDIMA         ),
*    $         ADVDX2 (        1        ), SDVDX2 (        1        ),
     $         ADVDY2 (   IDIMA         ), SDVDY2 (   IDIMA         ),
*    $         ADVDY2 (        1        ), SDVDY2 (        1        ),
     $         ADVDZ2 (   IDIMA         ), SDVDZ2 (   IDIMA         ),
*    $         ADVDZ2 (        1        ), SDVDZ2 (        1        ),
     $         ADWDX2 (   IDIMA         ), SDWDX2 (   IDIMA         ),
*    $         ADWDX2 (        1        ), SDWDX2 (        1        ),
     $         ADWDY2 (   IDIMA         ), SDWDY2 (   IDIMA         ),
*    $         ADWDY2 (        1        ), SDWDY2 (        1        ),
     $         ADWDZ2 (   IDIMA         ), SDWDZ2 (   IDIMA         )
*    $         ADWDZ2 (        1        ), SDWDZ2 (        1        )
      REAL     AUFWFG (   IDIMA         ), SUFWFG (   IDIMA         ),
***   REAL     AUFWFG (        1        ), SUFWFG (        1        ),
     $         AUFWFS (   IDIMA         ), SUFWFS (   IDIMA         ),
***  $         AUFWFS (        1        ), SUFWFS (        1        ),
     $         AUFWFM (   IDIMA         ), SUFWFM (   IDIMA         )
***  $         AUFWFM (        1        ), SUFWFM (        1        )
      REAL     AVFWFG (   IDIMA         ), SVFWFG (   IDIMA         ),
***   REAL     AVFWFG (        1        ), SVFWFG (        1        ),
     $         AVFWFS (   IDIMA         ), SVFWFS (   IDIMA         ),
***  $         AVFWFS (        1        ), SVFWFS (        1        ),
     $         AVFWFM (   IDIMA         ), SVFWFM (   IDIMA         ),
***  $         AVFWFM (        1        ), SVFWFM (        1        ),
     $         AUFVFG (   IDIMA         ), SUFVFG (   IDIMA         ),
***  $         AUFVFG (        1        ), SUFVFG (        1        ),
     $         AUFVFS (   IDIMA         ), SUFVFS (   IDIMA         ),
***  $         AUFVFS (        1        ), SUFVFS (        1        ),
     $         AUFVFM (   IDIMA         ), SUFVFM (   IDIMA         )
***  $         AUFVFM (        1        ), SUFVFM (        1        )
*     REAL     OX     (  IDIM3D      ),
      REAL     OX     (       1      ),
*    $         OXFG   (  IDIM3D      ),
     $         OXFG   (       1      ),
*    $         AOX    (   IDIMA         ), SOX    (   IDIMA         ),
     $         AOX    (        1        ), SOX    (        1        ),
*    $         AOXRG  (   IDIMA         ), SOXRG  (   IDIMA         )
     $         AOXRG  (        1        ), SOXRG  (        1        )
*     REAL     OY     (  IDIM3D      ),
      REAL     OY     (       1      ),
*    $         OYFG   (  IDIM3D      ),
     $         OYFG   (       1      ),
*    $         AOY    (   IDIMA         ), SOY    (   IDIMA         ),
     $         AOY    (        1        ), SOY    (        1        ),
*    $         AOYRG  (   IDIMA         ), SOYRG  (   IDIMA         )
     $         AOYRG  (        1        ), SOYRG  (        1        )
*     REAL     OZ     (  IDIM3D      ),
      REAL     OZ     (       1      ),
*    $         OZFG   (  IDIM3D      ),
     $         OZFG   (       1      ),
*    $         AOZ    (   IDIMA         ), SOZ    (   IDIMA         ),
     $         AOZ    (        1        ), SOZ    (        1        ),
*    $         AOZRG  (   IDIMA         ), SOZRG  (   IDIMA         )
     $         AOZRG  (        1        ), SOZRG  (        1        )
*     REAL                                 O2FG   (  IDIM3D      ),
      REAL                                 O2FG   (       1      ),
*    $         AO2    (   IDIMA         ), SO2    (   IDIMA         ),
     $         AO2    (        1        ), SO2    (        1        ),
*    $         AO2RG  (   IDIMA         ), SO2RG  (   IDIMA         )
     $         AO2RG  (        1        ), SO2RG  (        1        )
*     REAL                                 HEFG   (  IDIM3D      ),
      REAL                                 HEFG   (       1      ),
*    $         AHE    (   IDIMA         ), SHE    (   IDIMA         ),
     $         AHE    (        1        ), SHE    (        1        ),
*    $         AHERG  (   IDIMA         ), SHERG  (   IDIMA         )
     $         AHERG  (        1        ), SHERG  (        1        )
C
#ifdef _TAU_NN_
C                                  NICHT-NEWTONSCHE SPANNUNGEN
C                                  ---------------------------C
      REAL
     $         ATAU11 (   IDIMA         ), STAU11 (   IDIMA         ),
*    $         ATAU11 (        1        ), STAU11 (        1        ),
     $         ATAU12 (   IDIMA         ), STAU12 (   IDIMA         ),
*    $         ATAU12 (        1        ), STAU12 (        1        ),
     $         ATAU13 (   IDIMA         ), STAU13 (   IDIMA         ),
*    $         ATAU13 (        1        ), STAU13 (        1        ),
     $         ATAU22 (   IDIMA         ), STAU22 (   IDIMA         ),
*    $         ATAU22 (        1        ), STAU22 (        1        ),
     $         ATAU23 (   IDIMA         ), STAU23 (   IDIMA         ),
*    $         ATAU23 (        1        ), STAU23 (        1        ),
     $         ATAU33 (   IDIMA         ), STAU33 (   IDIMA         )
*    $         ATAU33 (        1        ), STAU33 (        1        )
#else
      REAL
     $         ATAU11 (        1        ), STAU11 (        1        ),
     $         ATAU12 (        1        ), STAU12 (        1        ),
     $         ATAU13 (        1        ), STAU13 (        1        ),
     $         ATAU22 (        1        ), STAU22 (        1        ),
     $         ATAU23 (        1        ), STAU23 (        1        ),
     $         ATAU33 (        1        ), STAU33 (        1        )
#endif
C
C                                  KREUZKORRELATIONSFUNKTIONEN
C                                  ---------------------------
C
*     REAL     ARXUV  (   IDIM2L        ), SRXUV  (     IDIM1L      ),
      REAL     ARXUV  (        1        ), SRXUV  (        1        ),
*    $         ARYUV  (   IDIM2L        ), SRYUV  (     IDIM1L      ),
     $         ARYUV  (        1        ), SRYUV  (        1        ),
*    $         ARXUW  (   IDIM2L        ), SRXUW  (     IDIM1L      ),
     $         ARXUW  (        1        ), SRXUW  (        1        ),
*    $         ARYUW  (   IDIM2L        ), SRYUW  (     IDIM1L      ),
     $         ARYUW  (        1        ), SRYUW  (        1        ),
*    $         ARYVW  (   IDIM2L        ), SRYVW  (     IDIM1L      )
     $         ARYVW  (        1        ), SRYVW  (        1        )
C
C                                  KORRELATIONSKOEFFIZIENTEN
C                                  -------------------------
C
*     REAL     ACXUW  (   IDIM2L        ), SCXUW  (     IDIM1L      ),
      REAL     ACXUW  (        1        ), SCXUW  (        1        ),
*    $         ACZUW  (   IDIM2L        ), SCZUW  (     IDIM1L      )
     $         ACZUW  (        1        ), SCZUW  (        1        )
C
C                                  OBERE DREIECKSMATRIX DES TENSORS
C                                  DER KORRELATIONSFUNKTIONEN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  --------------------------------
C
*     REAL     ARXOXX (   IDIM2L        ), SRXOXX (     IDIM1L      ),
      REAL     ARXOXX (        1        ), SRXOXX (        1        ),
*    $         ARXOXY (   IDIM2L        ), SRXOXY (     IDIM1L      ),
     $         ARXOXY (        1        ), SRXOXY (        1        ),
*    $         ARXOXZ (   IDIM2L        ), SRXOXZ (     IDIM1L      ),
     $         ARXOXZ (        1        ), SRXOXZ (        1        ),
*    $         ARXOYY (   IDIM2L        ), SRXOYY (     IDIM1L      ),
     $         ARXOYY (        1        ), SRXOYY (        1        ),
*    $         ARXOYZ (   IDIM2L        ), SRXOYZ (     IDIM1L      ),
     $         ARXOYZ (        1        ), SRXOYZ (        1        ),
*    $         ARXOZZ (   IDIM2L        ), SRXOZZ (     IDIM1L      )
     $         ARXOZZ (        1        ), SRXOZZ (        1        )
*     REAL     ARYOXX (   IDIM2L        ), SRYOXX (     IDIM1L      ),
      REAL     ARYOXX (        1        ), SRYOXX (        1        ),
*    $         ARYOXY (   IDIM2L        ), SRYOXY (     IDIM1L      ),
     $         ARYOXY (        1        ), SRYOXY (        1        ),
*    $         ARYOXZ (   IDIM2L        ), SRYOXZ (     IDIM1L      ),
     $         ARYOXZ (        1        ), SRYOXZ (        1        ),
*    $         ARYOYY (   IDIM2L        ), SRYOYY (     IDIM1L      ),
     $         ARYOYY (        1        ), SRYOYY (        1        ),
*    $         ARYOYZ (   IDIM2L        ), SRYOYZ (     IDIM1L      ),
     $         ARYOYZ (        1        ), SRYOYZ (        1        ),
*    $         ARYOZZ (   IDIM2L        ), SRYOZZ (     IDIM1L      )
     $         ARYOZZ (        1        ), SRYOZZ (        1        )
*     REAL     ARZOXX (   IDIM2L        ), SRZOXX (     IDIM1L      ),
      REAL     ARZOXX (        1        ), SRZOXX (        1        ),
*    $         ARZOXY (   IDIM2L        ), SRZOXY (     IDIM1L      ),
     $         ARZOXY (        1        ), SRZOXY (        1        ),
*    $         ARZOXZ (   IDIM2L        ), SRZOXZ (     IDIM1L      ),
     $         ARZOXZ (        1        ), SRZOXZ (        1        ),
*    $         ARZOYY (   IDIM2L        ), SRZOYY (     IDIM1L      ),
     $         ARZOYY (        1        ), SRZOYY (        1        ),
*    $         ARZOYZ (   IDIM2L        ), SRZOYZ (     IDIM1L      ),
     $         ARZOYZ (        1        ), SRZOYZ (        1        ),
*    $         ARZOZZ (   IDIM2L        ), SRZOZZ (     IDIM1L      )
     $         ARZOZZ (        1        ), SRZOZZ (        1        )
C
C                                  LEISTUNGSDICHTESPEKTREN DER
C                                  FLUKTUATIONEN DER VORTICITY
C                                  ----------------------------
C
*     REAL     ASXOXX (     IDIM1L      ), SSXOXX (     IDIM1L      ),
      REAL     ASXOXX (        1        ), SSXOXX (        1        ),
*    $         ASXOYY (     IDIM1L      ), SSXOYY (     IDIM1L      ),
     $         ASXOYY (        1        ), SSXOYY (        1        ),
*    $         ASXOZZ (     IDIM1L      ), SSXOZZ (     IDIM1L      ),
     $         ASXOZZ (        1        ), SSXOZZ (        1        ),
*    $         ASYOXX (     IDIM1L      ), SSYOXX (     IDIM1L      ),
     $         ASYOXX (        1        ), SSYOXX (        1        ),
*    $         ASYOYY (     IDIM1L      ), SSYOYY (     IDIM1L      ),
     $         ASYOYY (        1        ), SSYOYY (        1        ),
*    $         ASYOZZ (     IDIM1L      ), SSYOZZ (     IDIM1L      )
     $         ASYOZZ (        1        ), SSYOZZ (        1        )
C
C                                  HAEUFIGKEITSVERTEILUNG DER
C                                  INKLINATIONSWINKEL DER VORTICITY-
C                                  VEKTOREN
C                                  ---------------------------------
C
*     REAL     AHOZOY (     IDIM2L      ), SHOZOY (     IDIM2L      ),
      REAL     AHOZOY (        1        ), SHOZOY (        1        ),
*    $         AHOZOX (     IDIM2L      ), SHOZOX (     IDIM2L      ),
     $         AHOZOX (        1        ), SHOZOX (        1        ),
*    $         AHOYOX (     IDIM2L      ), SHOYOX (     IDIM2L      )
     $         AHOYOX (        1        ), SHOYOX (        1        )
