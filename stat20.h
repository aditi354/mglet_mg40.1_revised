
C
C           DIMENSIONIERUNGEN FUER DIE STATISTISCHE AUSWERTUNG
C
C                LEVEL:
C
c               <UUUU> G   <VVVU> G   <WWWU> G
c               <PPPP> G   <TTTT> G
C
c
C               OUTPUT LEVEL(DEFINE OPTION "STAT20FLA"):
C
C               <U-FLG> G  <V-FLG> G  <W-FLG> G <P-FLG> G
C               <T-FLG> G
C
c     <X-FLG> G =
C
C <XXXX>G-<X>^4-6<X>^2(<XX>G-<X>^2)-3<X>(<XXX>G-<X>^3-3<X>(<XX>G-<X>^2))
C ----------------------------------------------------------------------
C                         (<XX> G - <X>^2)^2 

      REAL     AUUUUM (   IDIMA         ), SUUUUM (   IDIMA         ),
     &         AVVVVM (   IDIMA         ), SVVVVM (   IDIMA         ),
     &         AWWWWM (   IDIMA         ), SWWWWM (   IDIMA         ),
     &         APPPPM (   IDIMA         ), SPPPPM (   IDIMA         )
#ifdef _TSCAL_
      REAL     ATTTTM (   IDIMA         ), STTTTM (   IDIMA         ),
     $          ATTTM  (   IDIMA         ), STTTM  (   IDIMA         )
#endif

