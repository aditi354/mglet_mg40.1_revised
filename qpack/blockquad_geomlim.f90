SUBROUTINE blockquad_geomlim(xc,yc,zc,ni,nj,nk,&
              nvirtual,ilim,xmaxg,nip,njp,nkp)
              
  !,xc,yc,zc)

  IMPLICIT NONE
  INTEGER i,j,k,ni,nj,nk,nvirtual,&
          ilim(3,2),&
          nip,njp,nkp
  REAL  xc(ni+2*nvirtual),yc(nj+2*nvirtual),&
        zc(nk+2*nvirtual),xmaxg(3,2)
  LOGICAL foundb,founde

  
 
  DO i = 1,3
     ilim(i,1) = -999
     ilim(i,2) = -999
  END DO


  foundb = .FALSE.
  founde = .FALSE.
  DO i = 1,nip
     IF (.NOT.foundb) THEN
        IF(xmaxg(1,1).LE.xc(i)) THEN
           foundb = .TRUE.
           ilim(1,1) = i
        END IF
     END IF
     IF (.NOT.founde) THEN
        IF(xmaxg(1,2).LE.xc(i)) THEN
           founde = .TRUE.
           ilim(1,2) = i
        END IF
     END IF
  END DO
  foundb = .FALSE.
  founde = .FALSE.
  DO j = 1,njp
     IF (.NOT.foundb) THEN
        IF(xmaxg(2,1).LE.yc(j)) THEN
           foundb = .TRUE.
           ilim(2,1) = j
        END IF
     END IF
     IF (.NOT.founde) THEN
        IF(xmaxg(2,2).LE.yc(j)) THEN
           founde = .TRUE.
           ilim(2,2) = j
        END IF
     END IF
  END DO
  foundb = .FALSE.
  founde = .FALSE.
  DO k = 1,nkp
     IF (.NOT.foundb) THEN
        IF(xmaxg(3,1).LE.zc(k)) THEN
           foundb = .TRUE.
           ilim(3,1) = k
        END IF
     END IF
     IF (.NOT.founde) THEN
        IF(xmaxg(3,2).LE.zc(k)) THEN
           founde = .TRUE.
           ilim(3,2) = k
        END IF
     END IF
  END DO
  IF (ilim(1,1).EQ.-999) ilim(1,1) = nip
  IF (ilim(2,1).EQ.-999) ilim(2,1) = njp
  IF (ilim(3,1).EQ.-999) ilim(3,1) = nkp
  IF (ilim(1,2).EQ.-999) ilim(1,2) = nip
  IF (ilim(2,2).EQ.-999) ilim(2,2) = njp
  IF (ilim(3,2).EQ.-999) ilim(3,2) = nkp
  ilim(1,1) = MAX(1, ilim(1,1) - 5)
  ilim(1,2) = MIN(nip, ilim(1,2) + 5)
  ilim(2,1) = MAX(1, ilim(2,1) - 5)
  ilim(2,2) = MIN(njp, ilim(2,2) + 5)
  ilim(3,1) = MAX(1, ilim(3,1) - 5)
  ilim(3,2) = MIN(nkp, ilim(3,2) + 5)


END SUBROUTINE
