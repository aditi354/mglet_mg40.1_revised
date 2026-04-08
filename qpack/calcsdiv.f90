SUBROUTINE calcsdiv(kk,jj,ii,u,v,w,sdiv, &
     rddx,rddy,rddz,bp)
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
  REAL u(kk,jj,ii),v(kk,jj,ii),w(kk,jj,ii),sdiv(kk,jj,ii),rddx(ii),rddy(jj),rddz(kk),bp(kk,jj,ii)


  INTEGER i,j,k

  DO i = 3,ii-2
     DO j = 3,jj-2
        DO k = 3,kk-2
           sdiv(k,j,i) = sdiv(k,j,i) + bp(k,j,i)*(1-bp(k,j,i+1))*  &
                (rddx(i))*u(k,j,i)
           sdiv(k,j,i) = sdiv(k,j,i) + bp(k,j,i)*(1-bp(k,j+1,i))*  &
                (rddy(j))*v(k,j,i)
           sdiv(k,j,i) = sdiv(k,j,i) + bp(k,j,i)*(1-bp(k+1,j,i))*  &
                (rddz(k))*w(k,j,i)
           sdiv(k,j,i) = sdiv(k,j,i) - bp(k,j,i)*(1-bp(k,j,i-1))*  &
                (rddx(i))*u(k,j,i-1)
           sdiv(k,j,i) = sdiv(k,j,i) - bp(k,j,i)*(1-bp(k,j-1,i))*  &
                (rddy(j))*v(k,j-1,i)
           sdiv(k,j,i) = sdiv(k,j,i) - bp(k,j,i)*(1-bp(k-1,j,i))*  &
                (rddz(k))*w(k-1,j,i)
        END DO
     END DO
  END DO
END SUBROUTINE calcsdiv
