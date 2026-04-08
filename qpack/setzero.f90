SUBROUTINE setzero (kk,jj,ii,u,v,w,p,bp,bu,bv,bw)
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

  INTEGER kk,jj,ii,i,j,k
  REAL u(kk,jj,ii),v(kk,jj,ii),w(kk,jj,ii),p(kk,jj,ii),bp(kk,jj,ii)
  REAL bu(kk,jj,ii)
  REAL  bv(kk,jj,ii),bw(kk,jj,ii)

  DO i=2,ii-1
     DO j=2,jj-1
        DO k=2,kk-1
           p(k,j,i) = p(k,j,i)*bp(k,j,i)
        END DO
     END DO
  END DO

  DO i=2,ii-1
     DO j=2,jj-1
        DO k=2,kk-1
           u(k,j,i) =  u(k,j,i)*bp(k,j,i)*bp(k,j,i+1)
        END DO
     END DO
  END DO

  DO i=2,ii-1
     DO j=2,jj-1
        DO k=2,kk-1
           v(k,j,i) = v(k,j,i)*bp(k,j,i)*bp(k,j+1,i)
        END DO
     END DO
  END DO

  DO i=2,ii-1
     DO j=2,jj-1
        DO k=2,kk-1
           w(k,j,i) = w(k,j,i)*bp(k,j,i)*bp(k+1,j,i)
        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE setzero

