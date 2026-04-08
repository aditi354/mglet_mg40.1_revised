C
C           DIMENSIONIERUNGEN FUER DIE STATISTISCHE AUSWERTUNG
C
C                LEVEL:
C
c               <UUU> G   <VVV> G   <WWW> G
c               <UUV> G   <UUW> G   <VVU> G
c               <VVW> G   <WWU> G   <WWV> G   <UVW> G
c               <UP> G   <VP> G   <WP> G
c               <(dU/dX)P> G   <(dV/dY)P> G   <(dW/dZ)P> G
c               <(dU/dY+dV/dX)P> G   <(dU/dZ+dW/dX)P> G
c               <(dV/dZ+dW/dY)P> G
c               <dU/dX dV/dX > G   <dU/dY dV/dY > G   <dU/dZ dV/dZ > G
c               <dU/dX dW/dX > G   <dU/dY dW/dY > G   <dU/dZ dW/dZ > G
c               <dV/dX dW/dX > G   <dV/dY dW/dY > G   <dV/dZ dW/dZ > G

      REAL     AUUUM  (   IDIMA         ), SUUUM  (   IDIMA         ),
     &         AVVVM  (   IDIMA         ), SVVVM  (   IDIMA         ),
     &         AWWWM  (   IDIMA         ), SWWWM  (   IDIMA         ),
     &         AUUVM  (   IDIMA         ), SUUVM  (   IDIMA         ),
     &         AUUWM  (   IDIMA         ), SUUWM  (   IDIMA         ),
     &         AVVUM  (   IDIMA         ), SVVUM  (   IDIMA         ),
     &         AVVWM  (   IDIMA         ), SVVWM  (   IDIMA         ),
     &         AWWUM  (   IDIMA         ), SWWUM  (   IDIMA         ),
     &         AWWVM  (   IDIMA         ), SWWVM  (   IDIMA         ),
     &         AUVWM  (   IDIMA         ), SUVWM  (   IDIMA         )
c
      REAL     AUPM   (   IDIMA         ), SUPM   (   IDIMA         ),
     &         AVPM   (   IDIMA         ), SVPM   (   IDIMA         ),
     &         AWPM   (   IDIMA         ), SWPM   (   IDIMA         )
c
      REAL     AUXPM  (   IDIMA         ), SUXPM  (   IDIMA         ),
     &         AVYPM  (   IDIMA         ), SVYPM  (   IDIMA         ),
     &         AWZPM  (   IDIMA         ), SWZPM  (   IDIMA         ),
     &         AUYVXP (   IDIMA         ), SUYVXP (   IDIMA         ),
     &         AUZWXP (   IDIMA         ), SUZWXP (   IDIMA         ),
     &         AVZWYP (   IDIMA         ), SVZWYP (   IDIMA         )
c
      REAL     AUXVXM (   IDIMA         ), SUXVXM (   IDIMA         ),
     &         AUYVYM (   IDIMA         ), SUYVYM (   IDIMA         ),
     &         AUZVZM (   IDIMA         ), SUZVZM (   IDIMA         ),
     &         AUXWXM (   IDIMA         ), SUXWXM (   IDIMA         ),
     &         AUYWYM (   IDIMA         ), SUYWYM (   IDIMA         ),
     &         AUZWZM (   IDIMA         ), SUZWZM (   IDIMA         ),
     &         AVXWXM (   IDIMA         ), SVXWXM (   IDIMA         ),
     &         AVYWYM (   IDIMA         ), SVYWYM (   IDIMA         ),
     &         AVZWZM (   IDIMA         ), SVZWZM (   IDIMA         )
