
	PARAMETER ( NPART_MAX =  1, NFAM_MAX = 1 )
C      PARAMETER ( NPART_MAX =   1, NFAM_MAX = 1 )

      COMMON /COPART/  NPART, NFAM, SHAPE_PART,
     $                 ASPECT, VOL_FRAC, IRANDOM, FRANDOM, 
     $                 ITSKIP, 
     $                 IPART_OUT,
     $                 CONF, UPART, DUPART, XPART, RRANDOM,
     $                 RR, RRRREE, TAU,
     $                 FXPART,FYPART,FZPART

      INTEGER IPART_OUT, IRANDOM, ITSKIP
      INTEGER NPART, NPART_MAX, NFAM, NFAM_MAX,
     $     IINDEX( NFAM_MAX ),
     $     JINDEX( NFAM_MAX ),
     $     KINDEX( NFAM_MAX )

      REAL   SHAPE_PART, ASPECT, VOL_FRAC, FRANDOM
      REAL    CONF  ( NFAM_MAX * NPART_MAX * 3 )
      REAL    UPART ( NFAM_MAX * 3 ),
     $       DUPART ( NFAM_MAX * 9 ),
     $        XPART ( NFAM_MAX * 3 ),
     $      RRANDOM ( NPART_MAX* 3 ),
     $           RR ( NFAM_MAX*9 ),
     $       RRRREE ( NFAM_MAX*9 ),
     $          TAU ( NFAM_MAX*9 )

      REAL FXPART(NFAM_MAX,4),FYPART(NFAM_MAX,4),FZPART(NFAM_MAX,4)

C--------------------------------------------------------
C       X,Y,Z DER FAMILIES BEFINDEN SICH DAMIT AUF:
C     XPART(       1), XPART(NFAM + 1) UND XPART(2*NFAM + 1)
