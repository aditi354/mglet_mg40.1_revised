SUBROUTINE blockquad(ni,nj,nk,nvirtual,x,y,z,ntopol,topol,maxblock,norder,&
     ndepth,gridnum,dx,dy,dz,ddx,ddy,ddz, ntrimax,xb,yb,zb,maccur, &
     itermax,sump1,sump2,xe,ye,ze,xmaxg,whatipol,bconds,narea,redord,bp,istenc)

  ! **********************************************************************
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !
  !   modifications:
  !   04.09.2003  modifications N.Peller 
  !   11. 8.2005  jk Moeglichkeit mehrerer paint-Gebiete (flow) wiedereingefuehrt
  !               offx, offy, offz abgeschafft
  ! **********************************************************************

  IMPLICIT NONE
  INTEGER ni,nj,nk,ntopol,maxblock,norder,&
       ndepth,bconds(6),gridnum,nvirtual,&
       ntrimax,whatipol,narea,redord,istenc
  REAL x(ni+2*nvirtual),y(nj+2*nvirtual),&
       z(nk+2*nvirtual),topol(4,3,ntopol),&
       dx(ni+2*nvirtual),dy(nj+2*nvirtual),&
       dz(nk+2*nvirtual),ddx(ni+2*nvirtual),&
       ddy(nj+2*nvirtual),ddz(nk+2*nvirtual),&
       xb(narea),yb(narea),zb(narea),maccur,&
       xe(narea),ye(narea),ze(narea),xmaxg(3,2)
  ! lokale Variablen
  INTEGER idnum,nip,njp,nkp,i,j,k,option,&
       blocked(ni,nj,nk),& 
       ib,jb,kb,ie,je,ke,itermax,sump1,sump2,&
       ilim(3,2),blockedp(ni,nj,nk),flow(3,2,narea)
  REAL xp(ni+1+2*(nvirtual)),yp(nj+1+2*(nvirtual)),&
       zp(nk+1+2*(nvirtual)),xc(ni+2*nvirtual),&
       yc(nj+2*nvirtual),zc(nk+2*nvirtual), &
       bp(ni,nj,nk)
  LOGICAL foundb,founde



  WRITE(*,*) 'ALLOCATING ',ni*nj*nk*2,' INTEGERS FOR TEMP BLOCKING'


  ! Herausfinden der Start und Endpunkte fuer
  ! den Paint-Algorithmus aus den Stroemungs-
  ! punkten in QPack.dat
  CALL blockquad_searcharea(flow,x,y,z,&
       ni,nj,nk,nvirtual,xb,yb,zb,xe,ye,ze,narea,gridnum)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !first the P cells
  if (istenc.ne.-1) then
     blockedp = nint(bp)
  else
     blockedp = 1
  end if

  option = 2
  idnum = 10*gridnum+4

  nip = ni+ 2*(nvirtual)
  njp = nj+ 2*(nvirtual)
  nkp = nk+ 2*(nvirtual)

  xp(1) = x(1)+0.5*dx(1)-ddx(1)
  DO i = 1,nip
     xp(i+1) = x(i)+0.5*dx(i)
  END DO
  yp(1) = y(1)+0.5*dy(1)-ddy(1)
  DO j = 1,njp
     yp(j+1) = y(j)+0.5*dy(j)
  END DO
  zp(1) = z(1)+0.5*dz(1)-ddz(1)
  DO k = 1,nkp
     zp(k+1) = z(k)+0.5*dz(k)
  END DO
  DO i = 1,nip
     xc(i) = x(i)
  END DO
  DO j = 1,njp
     yc(j) = y(j)
  END DO
  DO k = 1,nkp
     zc(k) = z(k)
  END DO

  WRITE(*,*) ' '
  WRITE(*,*) '   IDNUM',idnum
  WRITE(*,*) '   (blocking+stencils for pressure)'

  ! Herausfinden des Anfangs und des Endes des Suchgebietes
  CALL blockquad_geomlim(xc,yc,zc,ni,nj,nk,nvirtual,&
       ilim,xmaxg,nip,njp,nkp)
  WRITE(*,*) '   START SEARCH AREA AT :'
  WRITE(*,*) '   (Kb,Jb,Ib)',ilim(3,1),ilim(2,1),ilim(1,1)
  WRITE(*,*) '   END   SEARCH AREA AT :'
  WRITE(*,*) '   (Ke,Je,Ie)',ilim(3,2),ilim(2,2),ilim(1,2)

  ! Aufruf der Blocking Routine zum Erstellen 
  ! der stencils in Abhaengigkeit 
  ! von ib, ie, ... xc, ...
  ! Blocking for the pressure
  CALL blocking_xtrpol(nip,njp,nkp,xp,yp,zp,&
       ntopol,topol,maxblock,norder,&
       ndepth,option,idnum,xc,yc,zc,&
       ntrimax, blockedp,&
       4,maccur,itermax,sump1,sump2,&
       ilim,whatipol,gridnum,bconds,redord,narea,flow,istenc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! istenc =< 0: nur BP-Feld erzeugen
  if (istenc.le.0) then
     bp = real(blockedp)
     return
  end if
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  jk 29.7.2005 auskommentiert, um sicher zu gehen
!!$  !now the Scalar test 
  blocked = 1.0
  DO k = 1,nk
   DO j = 1,nj
    DO i = 2,ni
     blocked(i,j,k) = blockedp(i,j,k)
    END DO
   END DO
  END DO
!!$
  option = 3
!!$  offx = 0
!!$  offy = 0
!!$  offz = 0
  idnum = 10*gridnum+5
  nip = ni+ 2*(nvirtual)
  njp = nj+ 2*(nvirtual)
  nkp = nk+ 2*(nvirtual)
!!$
  xp(1) = x(1)+0.5*dx(1)-ddx(1)
  DO i = 1,nip
     xp(i+1) = x(i)+0.5*dx(i)
  END DO
  yp(1) = y(1)+0.5*dy(1)-ddy(1)
  DO j = 1,njp
     yp(j+1) = y(j)+0.5*dy(j)
  END DO
  zp(1) = z(1)+0.5*dz(1)-ddz(1)
  DO k = 1,nkp
     zp(k+1) = z(k)+0.5*dz(k)
  END DO
  DO i = 1,nip
     xc(i) = x(i)
  END DO
  DO j = 1,njp
     yc(j) = y(j)
  END DO
  DO k = 1,nkp
     zc(k) = z(k)
  END DO

  WRITE(*,*) ' '
  WRITE(*,*) '   IDNUM',idnum
  WRITE(*,*) '   (blocking+stencils for scalar)'
!!$
!!$  ! Herausfinden des Anfangs und des Endes des Suchgebietes
  CALL blockquad_geomlim(xc,yc,zc,ni,nj,nk,nvirtual,&
       ilim,xmaxg,nip,njp,nkp)
  WRITE(*,*) '   START SEARCH AREA AT :'
  WRITE(*,*) '   (Kb,Jb,Ib)',ilim(3,1),ilim(2,1),ilim(1,1)
  WRITE(*,*) '   END   SEARCH AREA AT :'
  WRITE(*,*) '   (Ke,Je,Ie)',ilim(3,2),ilim(2,2),ilim(1,2)
!!$
!!$  ! Aufruf der Blocking Routine zum Erstellen
!!$  ! der stencils in Abhaengigkeit
!!$  ! von ib, ie, ... xc, ...
!!$  ! Blocking for the pressure
!  CALL blocking_xtrpol(nip,njp,nkp,xp,yp,zp,&
!      ntopol,topol,maxblock,norder,&
!      ndepth,option,idnum,xc,yc,zc,&
!      ntrimax, blockedp, offx, offy, offz,&
!      1,ib,jb,kb,maccur,itermax,sump1,sump2,&
!      ie,je,ke,ilim,whatipol,gridnum,bconds,redord)
  CALL blocking_xtrpol(nip,njp,nkp,xp,yp,zp,&
       ntopol,topol,maxblock,norder,&
       ndepth,option,idnum,xc,yc,zc,&
       ntrimax, blockedp,&
       1,maccur,itermax,sump1,sump2,&
       ilim,whatipol,gridnum,bconds,redord,narea,flow,istenc)
!!$
!!$
!!$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! first the U1 grid
  ! Hier wird das Blocking-Feld fuer die
  ! U Geschwindigkeitszellen geschrieben
  ! und das entsprechende Gitter erzeugt.


  blocked = 1
  DO k = 1,nk
     DO j = 1,nj
        DO i = 2,ni
           ! blocking of the right u velocity of p-cell
           blocked(i,j,k) = blockedp(i,j,k)*blocked(i,j,k)
           ! blocking of the left u velocity of p-cell
           blocked(i-1,j,k) = blockedp(i,j,k)*blocked(i-1,j,k)             
        END DO
     END DO
  END DO

  option = 1
  idnum = 10*gridnum+1
  nip = ni + 2*(nvirtual)! - 1
  njp = nj + 2*(nvirtual)
  nkp = nk + 2*(nvirtual)

  DO i = 1,nip
     xp(i) = x(i)
  END DO
  xp(nip+1) = x(nip-1)+2*ddx(nip)
  yp(1) = y(1)-0.5*ddy(1)
  DO j = 1,njp
     yp(j+1) = y(j)+0.5*dy(j)
  END DO
  zp(1) = z(1)-0.5*ddz(1)
  DO k = 1,nkp
     zp(k+1) = z(k)+0.5*dz(k)
  END DO

  DO i = 1,nip
     xc(i) = x(i)+0.5*dx(i)
  END DO
  DO j = 1,njp
     yc(j) = y(j)
  END DO
  DO k = 1,nkp
     zc(k) = z(k)
  END DO

  WRITE(*,*) ' '
  WRITE(*,*) '   IDNUM',idnum
  WRITE(*,*) '   (blocking+stencils for u-velocity)'

  ! Herausfinden des Anfangs und des Endes des Suchgebietes
  CALL blockquad_geomlim(xc,yc,zc,ni,nj,nk,nvirtual,&
       ilim,xmaxg,nip,njp,nkp)
  WRITE(*,*) '   START SEARCH AREA AT :'
  WRITE(*,*) '   (Kb,Jb,Ib)',ilim(3,1),ilim(2,1),ilim(1,1)
  WRITE(*,*) '   END   SEARCH AREA AT :'
  WRITE(*,*) '   (Ke,Je,Ie)',ilim(3,2),ilim(2,2),ilim(1,2)

  CALL blocking_xtrpol(ni,nj,nk,xp,yp,zp,&
       ntopol,topol,maxblock,norder,&
       ndepth,option,idnum,xc,yc,zc,&
       ntrimax, blocked, &
       1,maccur,itermax,sump1,sump2,&
       ilim,whatipol,gridnum,bconds,redord,narea,flow,istenc)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! now the U2 cells
  ! Hier wird das Blocking-Feld fuer die
  ! V Geschwindigkeitszellen geschrieben
  ! und das entsprechende Gitter erzeugt.

  blocked = 1
  DO k = 1,nk
     DO j = 2,nj
        DO i = 1,ni
           blocked(i,j,k) = blockedp(i,j,k)*blocked(i,j,k)
           blocked(i,j-1,k) = blockedp(i,j,k)*blocked(i,j-1,k)             
        END DO
     END DO
  END DO
  option = 1
  idnum = 10*gridnum+2
  nip = ni + 2*(nvirtual)
  njp = nj + 2*(nvirtual)
  nkp = nk+ 2*(nvirtual)

  xp(1) = x(1)+0.5*dx(1)-ddx(1)
  DO i = 1,nip
     xp(i+1) = x(i)+0.5*dx(i)
  END DO
  DO j = 1,njp
     yp(j) = y(j)
  END DO
  yp(njp+1) = y(njp-1) + 2*ddy(njp)
  zp(1) = z(1)+0.5*dz(1)-ddz(1)
  DO k = 1,nkp
     zp(k+1) = z(k)+0.5*dz(k)
  END DO

  DO i = 1,nip
     xc(i) = x(i)
  END DO
  DO j = 1,njp
     yc(j) = y(j)+0.5*dy(j)
  END DO
  DO k = 1,nkp
     zc(k) = z(k)
  END DO

  WRITE(*,*) ' '
  WRITE(*,*) '   IDNUM',idnum
  WRITE(*,*) '   (blocking+stencils for v-velocity)'

  ! Herausfinden des Anfangs und des Endes des Suchgebietes
  CALL blockquad_geomlim(xc,yc,zc,ni,nj,nk,nvirtual,&
       ilim,xmaxg,nip,njp,nkp)
  WRITE(*,*) '   START SEARCH AREA AT :'
  WRITE(*,*) '   (Kb,Jb,Ib)',ilim(3,1),ilim(2,1),ilim(1,1)
  WRITE(*,*) '   END   SEARCH AREA AT :'
  WRITE(*,*) '   (Ke,Je,Ie)',ilim(3,2),ilim(2,2),ilim(1,2)

  CALL blocking_xtrpol(ni,nj,nk,xp,yp,zp,&
       ntopol,topol,maxblock,norder,&
       ndepth,option,idnum,xc,yc,zc,& 
       ntrimax, blocked,&
       2,maccur,itermax,sump1,sump2,&
       ilim,whatipol,gridnum,bconds,redord,narea,flow,istenc)








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! now the U3 cells
  ! Hier wird das Blocking-Feld fuer die 
  ! W Geschwindigkeitszellen geschrieben
  ! und das entsprechende Gitter erzeugt. 

  blocked = 1
  DO k = 2,nk
     DO j = 1,nj
        DO i = 1,ni
           blocked(i,j,k) = blockedp(i,j,k)*blocked(i,j,k)
           blocked(i,j,k-1) = blockedp(i,j,k)*blocked(i,j,k-1)             
        END DO
     END DO
  END DO
  option = 1
  idnum = 10*gridnum+3
  nip = ni + 2*(nvirtual)
  njp = nj + 2*(nvirtual)
  nkp = nk + 2*(nvirtual)

  xp(1) = x(1)+0.5*dx(1)-ddx(1)
  DO i = 1,nip
     xp(i+1) = x(i)+0.5*dx(i)
  END DO
  yp(1) = y(1)+0.5*dy(1)-ddy(1)
  DO j = 1,njp
     yp(j+1) = y(j)+0.5*dy(j)
  END DO
  DO k = 1,nkp
     zp(k) = z(k)
  END DO
  zp(nkp+1) = z(nkp-1) + 2*ddz(nkp)

  DO i = 1,nip
     xc(i) = x(i)
  END DO
  DO j = 1,njp
     yc(j) = y(j)
  END DO
  DO k = 1,nkp
     zc(k) = z(k)+0.5*dz(k)
  END DO

  WRITE(*,*) ' '
  WRITE(*,*) '   IDNUM',idnum
  WRITE(*,*) '   (blocking+stencils for w-velocity)'

  ! Herausfinden des Anfangs und des Endes des Suchgebietes
  CALL blockquad_geomlim(xc,yc,zc,ni,nj,nk,nvirtual,&
       ilim,xmaxg,nip,njp,nkp)
  WRITE(*,*) '   START SEARCH AREA AT :'
  WRITE(*,*) '   (Kb,Jb,Ib)',ilim(3,1),ilim(2,1),ilim(1,1)
  WRITE(*,*) '   END   SEARCH AREA AT :'
  WRITE(*,*) '   (Ke,Je,Ie)',ilim(3,2),ilim(2,2),ilim(1,2)

  CALL blocking_xtrpol(ni,nj,nk,xp,yp,zp,&
       ntopol,topol,maxblock,norder,&
       ndepth,option,idnum,xc,yc,zc,&
       ntrimax, blocked,&
       3,maccur,itermax,sump1,sump2,&
       ilim,whatipol,gridnum,bconds,redord,narea,flow,istenc)

  WRITE(*,*) ' '
  WRITE(*,*) 'DEALLOCATING ',ni*nj*nk*2,' REALS FOR TEMP BLOCKING'




END SUBROUTINE blockquad




