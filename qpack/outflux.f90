SUBROUTINE outflux(kk,jj,ii,u,v,w,bp,bpc,stot,flux,ddx,ddy,ddz, &
     rdx,rdy,rdz,option,corr,gfluxt)
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
  REAL u(kk,jj,ii),v(kk,jj,ii),w(kk,jj,ii),bp(kk,jj,ii),flux
  REAL ddx(ii),ddy(jj),ddz(kk), rdx(ii), rdy(jj),rdz(kk),bpc(kk,jj,ii), &
       stot,corr,gfluxt
  LOGICAL, SAVE :: firstcall 
  INTEGER i,j,k,option, ip1,ip2,ip3,ip4,ip5,ip6,ip7
  INTEGER, SAVE :: listu,listv,listw,ncountu,ncountv,ncountw,maxpts
  REAL , SAVE :: facu,facv,facw,stot2,gflux 
  REAL areaxy,areaxz,areayz
  ALLOCATABLE :: facu(:,:),facv(:,:),facw(:,:),listu(:),listv(:),listw(:)
  DATA firstcall /.TRUE./


  IF (option.EQ.0) THEN
     IF (firstcall) THEN
        firstcall = .FALSE.
        OPEN (unit=13,file="QPack.dat",form='formatted')
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*) 
        READ(13,*)
        READ(13,*)
        READ(13,*) 
        READ(13,*) 
        READ(13,*) maxpts
        READ(13,*)       
        READ(13,*) gflux
        CLOSE(13)
        ALLOCATE (facu(maxpts,2),facv(maxpts,2),facw(maxpts,2),listu(maxpts),listv(maxpts),listw(maxpts))
        ncountu = 0
        ncountv = 0
        ncountw = 0
        stot2 = 0.0
        flux = 0.0
        DO i=3,ii-2
           DO j=3,jj-2
              IF (ncountu.GT.maxpts.OR.ncountv.GT.maxpts.OR.ncountw.GT.maxpts) THEN
                 WRITE(*,*) 'outflux', ncountu,ncountv,ncountw
                 STOP
              END IF
              DO k=3,kk-2
                 ip1 = k+(j-1)*kk+(i-1)*kk*jj
                 ip2 = k+(j-1)*kk+(i-1+1)*kk*jj
                 ip3 = k+(j-1)*kk+(i-1-1)*kk*jj
                 ip4 = k+(j-1+1)*kk+(i-1)*kk*jj
                 ip5 = k+(j-1-1)*kk+(i-1)*kk*jj
                 ip6 = k+1+(j-1)*kk+(i-1)*kk*jj
                 ip7 = k-1+(j-1)*kk+(i-1)*kk*jj
                 areayz = ddy(j)*ddz(k)
                 areaxz = ddx(i)*ddz(k)
                 areaxy = ddx(i)*ddy(j)
                 IF ((1.0-bp(k,j,i+1))*bp(k,j,i).EQ.1) THEN
                    ncountu = ncountu+1
                    listu(ncountu) = ip1
                    facu(ncountu,1) = areayz  
                    stot2 = stot2 + rdx(i)*areayz
                    facu(ncountu,2) = rdx(i)
                 END IF
                 IF ((1.0-bp(k,j,i-1))*bp(k,j,i).EQ.1) THEN
                    ncountu = ncountu+1
                    listu(ncountu) = ip3
                    facu(ncountu,1) = -areayz  
                    stot2 = stot2 + rdx(i-1)*areayz
                    facu(ncountu,2) = -rdx(i-1)
                 END IF
                 IF ((1.0-bp(k,j+1,i))*bp(k,j,i).EQ.1) THEN
                    ncountv = ncountv+1
                    listv(ncountv) = ip1
                    facv(ncountv,1) = areaxz  
                    stot2 = stot2 + rdy(j)*areaxz
                    facv(ncountv,2) = rdy(j)
                 END IF
                 IF ((1.0-bp(k,j-1,i))*bp(k,j,i).EQ.1) THEN
                    ncountv = ncountv+1
                    listv(ncountv) = ip5
                    facv(ncountv,1) = -areaxz
                    stot2 = stot2 + rdy(j-1)*areaxz
                    facv(ncountv,2) = -rdy(j-1)
                 END IF
                 IF ((1.0-bp(k+1,j,i))*bp(k,j,i).EQ.1) THEN
                    ncountw = ncountw+1
                    listw(ncountw) = ip1
                    facw(ncountw,1) = areaxy
                    stot2 = stot2 + rdz(k)*areaxy
                    facw(ncountw,2) = rdz(k)
                 END IF
                 IF ((1.0-bp(k-1,j,i))*bp(k,j,i).EQ.1) THEN
                    ncountw = ncountw+1
                    listw(ncountw) = ip7
                    facw(ncountw,1) = -areaxy
                    stot2 = stot2 + rdz(k-1)*areaxy
                    facw(ncountw,2) = - rdz(k-1)
                 END IF

              END DO
           END DO
        END DO
     END IF

     stot = stot2
     flux = 0.0
     gfluxt = gflux


     CALL compflux(kk,jj,ii,u,v,w,ncountu,ncountv,ncountw,listu,listv,listw,facu(1,1),facv(1,1),facw(1,1),flux,option)

  ELSE

     CALL compflux(kk,jj,ii,u,v,w,ncountu,ncountv,ncountw,listu,listv,listw,facu(1,2),facv(1,2),facw(1,2),corr,option)
  END IF

END SUBROUTINE outflux
