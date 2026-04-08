SUBROUTINE bubvbw(kk,jj,ii,bu,bv,bw,bp)
  ! ======================================================================
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !
  !   copyright (c) 1999-2002
  !   all rights reserved TU-Muenchen and F.Tremblay. No part of this publication
  !   may be reproduced,
  !   stored in a retrieval system ( e.g., in memory, disk, or core)
  !   or be transmitted by any means, electronic, mechanical, photocopy,
  !   recording, or otherwise, without written permission from the
  !   author.
  !
  !
  ! **********************************************************************
  IMPLICIT NONE
  INTEGER kk,jj,ii
  REAL bu(kk,jj,ii),bv(kk,jj,ii),bw(kk,jj,ii),bp(kk,jj,ii)

  INTEGER i,j,k

  DO i = 1,ii-1
     DO j = 1,jj
        DO k = 1,kk
           bu(k,j,i) = bp(k,j,i)*bp(k,j,i+1)
        END DO
     END DO
  END DO
  DO j=1,jj
     DO k = 1,kk
        bu(k,j,ii) = bp(k,j,ii)
     END DO
  END DO

  DO i = 1,ii
     DO j = 1,jj-1
        DO k = 1,kk
           bv(k,j,i) = bp(k,j,i)*bp(k,j+1,i)
        END DO
     END DO
  END DO
  DO i=1,ii
     DO k = 1,kk
        bv(k,jj,i) = bp(k,jj,i)
     END DO
  END DO


  DO i = 1,ii
     DO j = 1,jj
        DO k = 1,kk-1
           bw(k,j,i) = bp(k,j,i)*bp(k+1,j,i)
        END DO
     END DO
  END DO
  DO i=1,ii
     DO j = 1,jj
        bw(kk,j,i) = bp(kk,j,i)
     END DO
  END DO

END SUBROUTINE bubvbw


