
      PARAMETER ( NPART_MAX = 20000, NFAM_MAX = 20000 )

      COMMON /COPART2/  NPART, NFAM, RUN,WLOP,
     $                 IPART_OUT,NXNPART,PARREAD,
     $                 IPART, UPART, DUPART,
     $                 XPART,FXPART,FYPART,FZPART,
     $                 START_VEC,DS,PFLOWCASE
      LOGICAL PARREAD

      INTEGER IPART_OUT,WLOP,PFLOWCASE
      INTEGER NPART, NPART_MAX, NFAM, NFAM_MAX,
     $     NXNPART,
     $        IPART (    NPART_MAX ),
     $     IINDEX( NFAM_MAX ),IINDEX_N(NFAM_MAX),
     $     JINDEX( NFAM_MAX ),JINDEX_N(NFAM_MAX),
     $     KINDEX( NFAM_MAX ),KINDEX_N(NFAM_MAX),
     $     RUN   ( NFAM_MAX )

      REAL   SHAPE_PART
      REAL   DISSI(NFAM_MAX),
     $     DISTURB(NFAM_MAX),
     $     PARTSCA(NFAM_MAX),
     $     START_VEC ( NFAM_MAX*3 ),
     $     XPART (NFAM_MAX*3),
     $     UPART (NFAM_MAX*3),
     $     DUPART (NFAM_MAX*9),
     $     DS(NFAM_MAX)
      REAL    FXPART(NFAM_MAX,4),FYPART(NFAM_MAX,4),FZPART(NFAM_MAX,4)

C--------------------------------------------------------
C       X,Y,Z DER FAMILIES BEFINDEN SICH DAMIT AUF:
C     XPART(       1), XPART(NFAM + 1) UND XPART(2*NFAM + 1)
