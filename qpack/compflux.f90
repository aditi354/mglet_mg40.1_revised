SUBROUTINE compflux(kk,jj,ii,u,v,w,ncountu,ncountv,ncountw,listu,listv,listw,facu,facv,facw,flux,option)
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

  INTEGER kk,jj,ii,ncountu,ncountv,ncountw,listu(ncountu),listv(ncountv),listw(ncountw),option
  REAL flux, facu(ncountu),facv(ncountv),facw(ncountw)
  REAL u(kk*jj*ii),v(kk*jj*ii),w(kk*jj*ii)


  INTEGER i, ip

  IF (option.EQ.0) THEN

     DO i = 1,ncountu
        ip = listu(i)
        flux = flux+u(ip)*facu(i)
     END DO
     DO i = 1,ncountv
        ip = listv(i)
        flux = flux+v(ip)*facv(i)
     END DO
     DO i = 1,ncountw
        ip = listw(i)
        flux = flux+w(ip)*facw(i)
     END DO

  ELSE
     DO i = 1,ncountu
        ip = listu(i)
        u(ip) = u(ip)+flux*facu(i)
     END DO
     DO i = 1,ncountv
        ip = listv(i)
        v(ip) = v(ip)+flux*facv(i)
     END DO
     DO i = 1,ncountw
        ip = listw(i)
        w(ip) = w(ip)+flux*facw(i)
     END DO

  END IF
END SUBROUTINE compflux
