
      PARAMETER (NMBODOLD = 100)
      
      COMMON /COBODYOLD/
     &               NBODOLD,  CTYPO,
     &               IB1O,   IB2O,   JB1O,   JB2O,   KB1O,   KB2O,
     &               XB1O,   XB2O,   YB1O,   YB2O,   ZB1O,   ZB2O,
     &               NCOUNO, XMITO,  HEIGHTO, ALPHAO, CDIRO

      INTEGER
     &       NBODOLD,
     &       IB1O(NMBODOLD),   IB2O(NMBODOLD),
     &       JB1O(NMBODOLD),   JB2O(NMBODOLD),
     &       KB1O(NMBODOLD),   KB2O(NMBODOLD),
     &       NCOUNO(NMBODOLD) 

      REAL
     &       XB1O(NMBODOLD),   XB2O(NMBODOLD),
     &       YB1O(NMBODOLD),   YB2O(NMBODOLD),
     &       ZB1O(NMBODOLD),   ZB2O(NMBODOLD),
     &      XMITO(NMBODOLD),HEIGHTO(NMBODOLD),ALPHAO(NMBODOLD)

      CHARACTER (LEN=16) CTYPO(NMBODOLD),CDIRO(NMBODOLD)
 
