SUBROUTINE blocking(ni,nj,nk,x,y,z,ntopol,topol,&
     blocked,tmptabi,ntrimax,&
     xc,yc,zc,maccur,itermax,ilim,option)
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
  !   input: 
  !   x(ni+1),y(nj+1),z(nk+1) : Koordinaten der Zellgrenzen
  !   xc(ni),yc(nj),zc(nk) : Koordinaten der Zellmitelpunkte
  !   topol(4,3,ntopol)    : \vec x_1, \vec x_2, \vec x_3, \vec u
  !                          der Dreiecke 
  !
  !   output:
  !   blocked(ni,nj,nk)    : -Anzahl der Dreiecke in einer Zelle
  !   tmptabi(ni,nj,nk,ntrimax) : Nummer des Dreiecks
  !
  !   modifications:
  !   13.02.2001  jk  changed loops: first loop over all triangels, then
  !                   loop only over cells near this triangels
  !   13. 7.2005  jk  ntrilim neu bestimmt, wesentlich kleinerer Suchbereich
  !                   Fortschrittszaehler
  ! **********************************************************************

  IMPLICIT NONE
  INTEGER ni,nj,nk,ntopol,&
       tmptabi(ni,nj,nk,ntrimax),ntrimax,&
       blocked(ni,nj,nk),option
  REAL topol(4,3,ntopol),&
       x(ni+1),y(nj+1),z(nk+1),xc(ni),yc(nj),zc(nk)

  INTEGER i,j,k,itri, nlist, nbloc,ibloc,&
       ibloc2,i2,j2,k2, ntri,iitri, iter,itermax,ncount,&
       icount,ic,jc, sum1, iii, ilim(3,2), &
       blockedint(ni,nj,nk),ntrilim(3,2),&
       imod, nmoda, nnmoda, nmodb, nnmodb,&
       imin,imax,jmin,jmax,kmin,kmax,ncount2,zaehler
  REAL xp,yp,zp,xx,yy,zz,a,b,c,d,x1,y1,z1,&
       x2,y2,z2,x3,y3,z3,jac,invjac,&
       px,py,pz,s1,s2,s3,abx,aby,abz,tiny, temp,eps,&
       maxsid, a1,b1,c1,d1,e1,f1,t1,found,nturn,&
       small,&
       eps1,eps2,eps3, small2, maccur, &
       xtrilim(3,2)
  LOGICAL check

  WRITE(*,*) '   ALLOCATING INTEGERS FOR PAINTING ALGO :',ni*nj*nk 

  ! Initialisierung einiger Variablen

  tmptabi = 0
  blockedint = 0
  tiny = maccur
  eps = maccur


  ! Schleife .bŽüber alle Dreiecke
  ! um deren Schnittpunkte mit dem
  ! kartesischen Gitter zu bestimmen
  zaehler = 0
  DO itri = 1, ntopol
     ! Die Koordinaten der drei
     ! Punkte eines Dreiecks
     x1 = topol(1,1,itri)
     x2 = topol(2,1,itri)
     x3 = topol(3,1,itri)
     y1 = topol(1,2,itri)
     y2 = topol(2,2,itri)
     y3 = topol(3,2,itri)
     z1 = topol(1,3,itri)
     z2 = topol(2,3,itri)
     z3 = topol(3,3,itri)

     ! Normalenvektor
     a = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
     b = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
     c = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
     ! d = x1*a + y1*b + z1*c

     ! Ausdehnung des Dreiecks im
     ! kartesischen Gitter --> 
     ! Koordinaten
     DO iii=1,3
        xtrilim(iii,1) = MIN (topol(1,iii,itri),&
             topol(2,iii,itri),topol(3,iii,itri))
        xtrilim(iii,2) = MAX (topol(1,iii,itri),&
             topol(2,iii,itri),topol(3,iii,itri))
     END DO

     ! Indices der Ausdehnung des Dreiecks
     ! bestimmen mit Hilfe der Koordinaten
     DO i = 1,3
        ntrilim(i,1) = -999
        ntrilim(i,2) = -999
     END DO

     ! Fuer x-Richtung
     DO i = ilim(1,1),ilim(1,2)
           IF(xtrilim(1,1).GT.x(i)) ntrilim(1,1) = i
           IF(xtrilim(1,2).GT.x(i)) ntrilim(1,2) = i
     END DO

     ! Fuer y-Richtung
     DO j = ilim(2,1),ilim(2,2)
           IF(xtrilim(2,1).GT.y(j)) ntrilim(2,1) = j
           IF(xtrilim(2,2).GT.y(j)) ntrilim(2,2) = j
     END DO

     ! Fuer z-Richtung
     DO k = ilim(3,1),ilim(3,2)
           IF(xtrilim(3,1).GT.z(k)) ntrilim(3,1) = k
           IF(xtrilim(3,2).GT.z(k)) ntrilim(3,2) = k
     END DO


     ! Wenn keine Grenzen bestimmt werden konnten,
     ! dann ilim untere Grenzen verwenden.
     IF (ntrilim(1,1).EQ.-999) ntrilim(1,1) = ilim(1,1)
     IF (ntrilim(2,1).EQ.-999) ntrilim(2,1) = ilim(2,1)
     IF (ntrilim(3,1).EQ.-999) ntrilim(3,1) = ilim(3,1)
     IF (ntrilim(1,2).EQ.-999) ntrilim(1,2) = ilim(1,1)
     IF (ntrilim(2,2).EQ.-999) ntrilim(2,2) = ilim(2,1)
     IF (ntrilim(3,2).EQ.-999) ntrilim(3,2) = ilim(3,1)


     ! Zellen suchen, die das Dreieck schneiden
     ! Innerhalb der zuvor bestimmten Index-Grenzen
     ! wird nach Schnittpunkten gesucht
     DO k = ntrilim(3,1),ntrilim(3,2)
        DO j = ntrilim(2,1),ntrilim(2,2)
           DO i = ntrilim(1,1),ntrilim(1,2)
              IF ((blocked(i,j,k).EQ.0).OR.&
                   (option.EQ.2)) THEN
                 !              (option.EQ.2)   THEN
                 !            sum1 = 0
                 eps1 = eps*ABS(x(i+1)-x(i))
                 eps2 = eps*ABS(y(j+1)-y(j))
                 eps3 = eps*ABS(z(k+1)-z(k))

                 found = 0
                 !               xp = 0.5*(x(i)+x(i+1))
                 !               yp = 0.5*(y(j)+y(j+1))
                 !               zp = 0.5*(z(k)+z(k+1))
                 xp = xc(i)
                 yp = yc(j)
                 zp = zc(k)
                 ! checking for surface
                 ! in the x-direction
                 IF (a.NE.0) THEN
                    xx = (-c*(zp-z1)-b*(yp-y1))/a + x1
                    yy = yp
                    zz = zp
                    IF ((xx.GE.x(i)-eps1).AND.(xx.LE.x(i+1)+eps1)) THEN
                       ! d1 Vektor berechnen = Gerade von Dreiecksseitenmittelpkt. 
                       ! zu Schnittpunkt xx,yy,zz
                       px = xx-0.5*(x2+x1)
                       py = yy-0.5*(y2+y1)
                       pz = zz-0.5*(z2+z1)

                       ! Hier nur ein paar AbkŽürzungen
                       abx= (x2-x1)
                       aby= (y2-y1)
                       abz= (z2-z1)

                       ! Jetzt Skalarprodukt d1 * ( P1P2 x P1P3 ) x P1P2
                       ! Vorzeichen vertauscht, um lt.tiny zu schreiben
                       s1=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       px = xx-0.5*(x3+x2)
                       py = yy-0.5*(y3+y2)
                       pz = zz-0.5*(z3+z2)

                       abx= (x3-x2)
                       aby= (y3-y2)
                       abz= (z3-z2)

                       s2=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       px = xx-0.5*(x1+x3)
                       py = yy-0.5*(y1+y3)
                       pz = zz-0.5*(z1+z3)

                       abx= (x1-x3)
                       aby= (y1-y3)
                       abz= (z1-z3)

                       ! Zum Schluss noch schauen, ob der Schnittpunkt
                       ! wirklich im Dreieck liegt.
                       s3=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       IF ((s1.LT.tiny).AND.(s2.LT.tiny).AND.(s3.LT.tiny)) THEN
                          found = 1
                          ! write(199,*) "pressure x-direction",xx,yy,zz
                       END IF
                    END IF
                 END IF
                 IF (found.NE.1.AND.b.NE.0) THEN
                    xx = xp
                    yy = (-c*(zp-z1)-a*(xp-x1))/b + y1
                    zz = zp
                    IF ((yy.GE.y(j)-eps2).AND.(yy.LE.y(j+1)+eps2)) THEN
                       px = xx-0.5*(x2+x1)
                       py = yy-0.5*(y2+y1)
                       pz = zz-0.5*(z2+z1)

                       abx= (x2-x1)
                       aby= (y2-y1)
                       abz= (z2-z1)

                       s1=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       px = xx-0.5*(x3+x2)
                       py = yy-0.5*(y3+y2)
                       pz = zz-0.5*(z3+z2)

                       abx= (x3-x2)
                       aby= (y3-y2)
                       abz= (z3-z2)

                       s2=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       px = xx-0.5*(x1+x3)
                       py = yy-0.5*(y1+y3)
                       pz = zz-0.5*(z1+z3)

                       abx= (x1-x3)
                       aby= (y1-y3)
                       abz= (z1-z3)

                       s3=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       IF ((s1.LT.tiny).AND.(s2.LT.tiny).AND.(s3.LT.tiny)) THEN
                          found = 1
                          ! write(199,*) "pressure y-direction",xx,yy,zz
                       END IF
                    END IF
                 END IF
                 IF (found.NE.1.AND.c.NE.0) THEN
                    xx = xp
                    yy = yp
                    zz = (-b*(yp-y1)-a*(xp-x1))/c + z1
                    IF ((zz.GE.z(k)-eps3).AND.(zz.LE.z(k+1)+eps3)) THEN
                       px = xx-0.5*(x2+x1)
                       py = yy-0.5*(y2+y1)
                       pz = zz-0.5*(z2+z1)

                       abx= (x2-x1)
                       aby= (y2-y1)
                       abz= (z2-z1)

                       s1=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       px = xx-0.5*(x3+x2)
                       py = yy-0.5*(y3+y2)
                       pz = zz-0.5*(z3+z2)

                       abx= (x3-x2)
                       aby= (y3-y2)
                       abz= (z3-z2)

                       s2=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       px = xx-0.5*(x1+x3)
                       py = yy-0.5*(y1+y3)
                       pz = zz-0.5*(z1+z3)

                       abx= (x1-x3)
                       aby= (y1-y3)
                       abz= (z1-z3)

                       s3=px*(aby*c-abz*b)+py*(abz*a-abx*c)+pz*(abx*b-aby*a)
                       IF ((s1.LT.tiny).AND.(s2.LT.tiny).AND.(s3.LT.tiny)) THEN
                          found = 1
                          ! write(199,*) "pressure z-direction",xx,yy,zz
                       END IF
                    END IF
                 END IF

                 IF (found.NE.0) THEN
                    blockedint(i,j,k) = blockedint(i,j,k) - 1 
                    IF (blockedint(i,j,k).LT.-ntrimax) THEN
                       WRITE(*,*) '   Too many triangles in one cell, enlarge ntrimax'
                       WRITE(*,*) '   ',ntrimax, blockedint(i,j,k), i,j,k
                       STOP
                    END IF
                    tmptabi(i,j,k,-blockedint(i,j,k)) = itri  
                 END IF
              END IF
           END DO
        END DO
     END DO
     IF (real(itri)/real(ntopol) .ge. zaehler*0.1) THEN
        WRITE(*,'(a17,f6.2,a1)') '       blocking: ',real(itri)/real(ntopol)*100.,'%'
        zaehler = zaehler+1
     END IF
  END DO


  ! Fuer das P-Blockfeld wird blockedint
  ! direkt in blocked geschrieben. Fuer die
  ! Geschwindigkeiten nicht.
  IF (option.EQ.2) THEN
     blocked = blockedint
  ELSE 
     DO k = ilim(3,1),ilim(3,2)
        DO j = ilim(2,1),ilim(2,2)
           DO i = ilim(1,1),ilim(1,2)
              IF (blocked(i,j,k).EQ.0) THEN
                 blocked(i,j,k)=blockedint(i,j,k)
              END IF
           END DO
        END DO
     END DO
  END IF



END SUBROUTINE blocking
