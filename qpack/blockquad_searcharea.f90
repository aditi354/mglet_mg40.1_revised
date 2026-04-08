SUBROUTINE blockquad_searcharea(flow,x,y,z,&
     ni,nj,nk,nvirtual,xb,yb,zb,xe,ye,ze,narea,gridnum)

  !   modifications:
  !   11. 8.2005  jk Moeglichkeit mehrerer paint-Gebiete (flow) wiedereingefuehrt

  IMPLICIT NONE
  INTEGER i,j,k,ni,nj,nk,nvirtual,&
       nip,njp,nkp,narea,gridnum, &
       flow(3,2,narea)
  REAL x(ni+2*nvirtual),y(nj+2*nvirtual),&
       z(nk+2*nvirtual),xb(narea),yb(narea),&
       zb(narea),xe(narea),ye(narea),ze(narea)
  ! lokale Variablen
  INTEGER iarea,foundarea, &
       ib,jb,kb,ie,je,ke
  LOGICAL foundb,founde

  nip = ni+ 2*(nvirtual)
  njp = nj+ 2*(nvirtual)
  nkp = nk+ 2*(nvirtual)


  ! Herausfinden der Start und Endpunkte 
  ! für den Paint-Algorithmus; das sind die
  ! Zellen, die als erste eine Geometrie
  ! enthalten

  iarea = 1
  foundarea = 0
  DO WHILE ((iarea.LE.narea))
!!$  DO WHILE ((foundarea.NE.1).AND.&
!!$       (iarea.LE.narea))
     ib = -998
     jb = -998
     kb = -998
     ie = -999
     je = -999
     ke = -999

     foundb = .FALSE.
     founde = .FALSE.
     IF (xe(iarea).GE.x(1)) THEN
        DO i = 1,nip
           !xc(i) = x(i)
           IF (.NOT.foundb) THEN
              IF(xb(iarea).LE.x(i)) THEN
                 foundb = .TRUE.
                 ib = i
              END IF
           END IF
           IF (.NOT.founde) THEN
              IF(xe(iarea).LE.x(i)) THEN
                 founde = .TRUE.
                 ie = i
              END IF
           END IF
        END DO
     END IF

     foundb = .FALSE.
     founde = .FALSE.
     IF (ye(iarea).GE.y(1)) THEN
        DO j = 1,njp
           !yc(j) = y(j)
           IF (.NOT.foundb) THEN
              IF(yb(iarea).LE.y(j)) THEN
                 foundb = .TRUE.
                 jb = j
              END IF
           END IF
           IF (.NOT.founde) THEN
              IF(ye(iarea).LE.y(j)) THEN
                 founde = .TRUE.
                 je = j
              END IF
           END IF
        END DO
     END IF

     foundb = .FALSE.
     founde = .FALSE.
     IF (ze(iarea).GE.z(1)) THEN
        DO k = 1,nkp
           !zc(k) = z(k)
           IF (.NOT.foundb) THEN
              IF(zb(iarea).LE.z(k)) THEN
                 foundb = .TRUE.
                 kb = k
              END IF
           END IF
           IF (.NOT.founde) THEN
              IF(ze(iarea).LE.z(k)) THEN
                 founde = .TRUE.
                 ke = k
              END IF
           END IF
        END DO
     END IF

     IF ((ib.NE.-998).AND.&
          (jb.NE.-998).AND.&
          (kb.NE.-998)) THEN        
        IF (ie.EQ.-999) ie = nip
        IF (je.EQ.-999) je = njp
        IF (ke.EQ.-999) ke = nkp
        foundarea = 1
     END IF


     WRITE(*,*) 'START PAINT AREA AT (Kb,Jb,Ib) : ',kb,jb,ib
     WRITE(*,*) 'END   PAINT AREA AT (Ke,Je,Ie) : ',ke,je,ie

     flow(1,1,iarea) = ib
     flow(1,2,iarea) = ie
     flow(2,1,iarea) = jb
     flow(2,2,iarea) = je
     flow(3,1,iarea) = kb
     flow(3,2,iarea) = ke

     iarea = iarea + 1
  END DO

  IF (foundarea.EQ.0) THEN
     WRITE(*,*)
     WRITE(*,*) 'ERROR: cannot find painting area.',gridnum
     WRITE(*,*)
     STOP
  END IF




END SUBROUTINE blockquad_searcharea
