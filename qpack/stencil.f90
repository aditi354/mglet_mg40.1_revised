SUBROUTINE stencil(ni,nj,nk,x,y,z,tmpabi,blocked,ntopol, topol,i,j,k,&
     nstenc,nused,const,stenci,stencr,ndepth,xc,yc,zc, ntrimax,im,jm,km,&
     used,iproxi2,compon,maccur,whatipol,gridnum,bconds,redord,option,&
     nredord,nresc1,nresc2)

  ! **********************************************************************
  !
  !        programmer     Nikolaus Peller
  !        version        1.0   
  !        date           04-09-2003
  !
  !   modifications:
  !   04.09.2003  derived from stencil4 subroutine 
  !   11. 8.2005  jk neuer rescue-stencil, s.u.
  ! **********************************************************************
  ! TODO : jk 11. 8.2005 maccur unsauber eingesetzt, absolut, nicht relativ
  ! TODO : jk 15. 8.2005 die Sache mit der Nachbarzelle

  IMPLICIT NONE
  INTEGER ntrimax,ntopol,im,jm,km,compon,i,j,k,&
       nused, ndepth, used(6),ni,nj,nk,&
       nstenc,stenci(3,6,nstenc+1),&
       tmpabi(ni,nj,nk,ntrimax),iproxi2(6),&
       whatipol,gridnum,bconds(6),&
       blocked(ni,nj,nk),redord,option,last, &
       nredord, nresc1, nresc2
  REAL x(ni+1),y(nj+1),z(nk+1),const(6),&
       stencr(6,nstenc+1),tmpcst(6),&
       xc(ni),yc(nj),zc(nk),&
       topol(4,3,ntopol),maccur
  ! lokal neu eingefuehrte Variablen
  INTEGER direction, position,jmpdir,&
       tmpsteni(nstenc+1,3,6),&
       counter,counter2,error(6,6),&
       nstencil, nstenc_reduced(6),&
       errout(6),errout2,casep,&
       it,jt,kt, nvalue, itri
  REAL tmpstenr(nstenc+1,6),distance(6),&
       sum,beta(6),sum_beta,product,gamma(6),value
!!$  REAL cdist(6),cgamma(6)

  LOGICAL firstcall, lopen, lcoarse
  DATA firstcall /.TRUE./
  SAVE :: firstcall

  ! check if it is really a blocked cell
  IF (blocked(i,j,k).EQ.1) THEN
     WRITE(*,*) 'ERROR, THE CELL IS NOT BLOCKED'
     STOP
  END IF

  ! For pressure interpolation - change indices
  IF(option.eq.2) THEN
     casep = 1
  ELSE
     casep = 0
  END IF


  const = 0.0
  iproxi2 = 0
  tmpcst = 0
  tmpsteni = -999
  tmpstenr = 0.0
  stencr = 0.0
  stenci = 0
  error = 0
  nstenc_reduced = 0
  distance = 0.0
  do direction = 1,6
     errout(direction) = 0
  enddo
  errout2 = 0

! Create stencil information
!---------------------------
 

     position = (ni*nj)*(k-1)+(ni)*(j-1)+i
     direction = 1
     jmpdir = 1
     CALL stencil_determine(ni,nj,nk,x,y,z,tmpabi,blocked,ntopol, topol,i,j,k,&
          nstenc,const,stenci,stencr,ndepth,xc,yc,zc, ntrimax,im,jm,km,&
          iproxi2,compon,maccur,tmpcst,tmpsteni,tmpstenr,distance,&
          direction,i,ni,im,x,xc,xc(i),jmpdir,position,&
          nstenc_reduced,error,whatipol,gridnum,bconds,redord,option,nredord)
     direction = 2
     jmpdir = ni
     CALL stencil_determine(ni,nj,nk,x,y,z,tmpabi,blocked,ntopol, topol,i,j,k,&
          nstenc,const,stenci,stencr,ndepth,xc,yc,zc, ntrimax,im,jm,km,&
          iproxi2,compon,maccur,tmpcst,tmpsteni,tmpstenr,distance,&
          direction,j,nj,jm,y,yc,yc(j),jmpdir,position,&
          nstenc_reduced,error,whatipol,gridnum,bconds,redord,option,nredord)
     direction = 3
     jmpdir = ni*nj
     CALL stencil_determine(ni,nj,nk,x,y,z,tmpabi,blocked,ntopol, topol,i,j,k,&
          nstenc,const,stenci,stencr,ndepth,xc,yc,zc, ntrimax,im,jm,km,&
          iproxi2,compon,maccur,tmpcst,tmpsteni,tmpstenr,distance,&
          direction,k,nk,km,z,zc,zc(k),jmpdir,position,&
          nstenc_reduced,error,whatipol,gridnum,bconds,redord,option,nredord)


  sum = 0.0
  used = 0
  nused = 0


! Wenn eine Nachbarzelle Stroemung enthaelt, 
! darf der Stencil nicht in die andere Richtung gehen.
! jk 15. 8.2005
  if (blocked(i-1,j,k).eq.1) nstenc_reduced(1) = 0
  if (blocked(i+1,j,k).eq.1) nstenc_reduced(2) = 0
  if (blocked(i,j-1,k).eq.1) nstenc_reduced(3) = 0
  if (blocked(i,j+1,k).eq.1) nstenc_reduced(4) = 0
  if (blocked(i,j,k-1).eq.1) nstenc_reduced(5) = 0
  if (blocked(i,j,k+1).eq.1) nstenc_reduced(6) = 0
 
! Die Sache mit der Nachbarzelle
! Wenn eine Nachbarzelle Stroemung enthaelt, 
! darf der Stencil nicht durch eine Geblockte gehen.
! Wichtig, damit nicht von falsche Koerperseite Daten genommen werden,
! verschlechtert allerdings das Verfahren bei hohen Seitenverhaeltnissen
! des kartesischen Gitters
! Empfehlung: bei "dicker" Geometrie auskommentieren
! jk 15. 8.2005
  lopen = .false.
  IF (blocked(min(i+1,ni),j,k).EQ.1) lopen = .true.
  IF (blocked(max(i-1, 1),j,k).EQ.1) lopen = .true.
  IF (blocked(i,min(j+1,nj),k).EQ.1) lopen = .true.
  IF (blocked(i,max(j-1, 1),k).EQ.1) lopen = .true.
  IF (blocked(i,j,min(k+1,nk)).EQ.1) lopen = .true.
  IF (blocked(i,j,max(k-1, 1)).EQ.1) lopen = .true.
  IF (lopen) THEN
     if (blocked(i+1,j,k).le.0) nstenc_reduced(1) = 0
     if (blocked(i-1,j,k).le.0) nstenc_reduced(2) = 0
     if (blocked(i,j+1,k).le.0) nstenc_reduced(3) = 0
     if (blocked(i,j-1,k).le.0) nstenc_reduced(4) = 0
     if (blocked(i,j,k+1).le.0) nstenc_reduced(5) = 0
     if (blocked(i,j,k-1).le.0) nstenc_reduced(6) = 0
  END IF

! Wenn zwei stencils in +/- Richtung dann
! eleminiere den stencil, dessen Randbedingung
! weiter vom zu interpolierenden Punkt weg ist
! For sharp edges to ensure one-sided interpolation
!---------------------------------------------
!  x-direction
  IF ( nstenc_reduced(1).NE.0 .AND. nstenc_reduced(2).NE.0 )  THEN
     IF((tmpsteni(1,1,1)-i).gt.(i-tmpsteni(1,1,2))) THEN
        nstenc_reduced(1) = 0
     ELSE
        nstenc_reduced(2) = 0
     END IF
  END IF

!  y-direction
  IF ( nstenc_reduced(3).NE.0 .AND. nstenc_reduced(4).NE.0 )  THEN
     IF((tmpsteni(1,2,3)-j).gt.(j-tmpsteni(1,2,4))) THEN
       nstenc_reduced(3) = 0
     ELSE
       nstenc_reduced(4) = 0
     END IF
  END IF

!  z-direction
  IF ( nstenc_reduced(5).NE.0 .AND. nstenc_reduced(6).NE.0 )  THEN
     IF((tmpsteni(1,3,5)-k).gt.(k-tmpsteni(1,3,6))) THEN
       nstenc_reduced(5) = 0
     ELSE
       nstenc_reduced(6) = 0
     END IF
  END IF



!!$! Weighting
!!$!--------------------------------------------------
!!$
!!$  ! CDIST itroduced for better wheighting ... let's see
!!$  cdist = 0.0
!!$  DO direction = 1,2
!!$     IF (nstenc_reduced(direction).NE.0) THEN
!!$        it = tmpsteni(1,1,direction)
!!$        cdist(direction) = abs(xc(i)-xc(it))
!!$     END IF
!!$  END DO
!!$  DO direction = 3,4
!!$     IF (nstenc_reduced(direction).NE.0) THEN
!!$        jt = tmpsteni(1,2,direction)
!!$        cdist(direction) = abs(yc(j)-yc(jt))
!!$     END IF
!!$  END DO
!!$  DO direction = 5,6
!!$     IF (nstenc_reduced(direction).NE.0) THEN
!!$        kt = tmpsteni(1,3,direction)
!!$        cdist(direction) = abs(zc(k)-zc(kt))
!!$     END IF
!!$  END DO
!!$
!!$  DO direction = 1,6
!!$     IF (nstenc_reduced(direction).NE.0) THEN
!!$        nused = nused + 1
!!$        used(nused) = direction
!!$        sum = sum + cdist(direction)
!!$     END IF
!!$  END DO
!!$
!!$  IF (sum.NE.0) THEN
!!$     DO direction = 1,6
!!$        cdist(direction) = cdist(direction) / sum
!!$     END DO
!!$  END IF
!!$
!!$  DO direction = 1,6
!!$     IF (nstenc_reduced(direction).NE.0) THEN
!!$        cdist(direction) = MAX(maccur,cdist(direction))
!!$     END IF
!!$  END DO
!!$
!!$  beta = 0.0
!!$  DO counter2 = 1,nused
!!$     product = 1.0
!!$     DO counter = 1,nused
!!$        IF (counter2.NE.counter) THEN
!!$           product = product * cdist(used(counter))
!!$        END IF
!!$     END DO
!!$     beta(counter2) = product / (cdist(used(counter2)))
!!$  END DO
!!$
!!$  sum_beta = 0.0
!!$  DO counter2 = 1,nused
!!$     sum_beta = sum_beta + beta(counter2)
!!$  END DO
!!$
!!$  cgamma = 0.0
!!$  DO counter2 = 1,nused
!!$     cgamma(counter2) = beta(counter2) / (sum_beta)
!!$  END DO




! Here wheigting for stencils
!------------------------------------------


  nused = 0
  sum = 0.0
  DO direction = 1,6
     IF (nstenc_reduced(direction).NE.0) THEN
        nused = nused + 1
        used(nused) = direction
        sum = sum + distance(direction)
     END IF
  END DO

  IF (sum.GT.0.0) THEN
     DO direction = 1,6
        distance(direction) = distance(direction) / sum
     END DO
  END IF

  DO direction = 1,6
     IF (nstenc_reduced(direction).NE.0) THEN
        distance(direction) = MAX(maccur,distance(direction))
     END IF
  END DO

  beta = 0.0
  DO counter2 = 1,nused
     product = 1.0
     DO counter = 1,nused
        IF (counter2.NE.counter) THEN
           product = product * distance(used(counter))
        END IF
     END DO
     beta(counter2) = product / (distance(used(counter2)))
  END DO

  sum_beta = 0.0
  DO counter2 = 1,nused
     sum_beta = sum_beta + beta(counter2)
  END DO

  gamma = 0.0
  DO counter2 = 1,nused
     gamma(counter2) = beta(counter2) / (sum_beta)
  END DO


! Combination of gamma and cgamma wheighting
! Only for velocities. For the pressure interpolation
! the wheighting factors gamma and cgamma are equivalent
!-------------------------------------------------------

!IF (option.NE.2) THEN
!  sum = 0.0
!  DO direction = 1,nused
!     sum = sum + gamma(direction) * cgamma(direction)
!  END DO
!  DO direction = 1,nused
!     gamma(direction) = gamma(direction) * cgamma(direction) / sum
!  END DO
!ENDIF



  ! Jetzt nur noch die berechneten Stencils mit gamma
  ! multiplizieren und in die endgueltigen Varaiablen schreiben.
  const = 0.0
  nstencil = 0
  DO counter2 = 1,nused
     nstencil = used(counter2)
     const(counter2) = gamma(counter2) * tmpcst(nstencil)
     DO counter = 1, (nstenc+casep)
        stenci(1,counter2,counter) = tmpsteni(counter,1,nstencil)
        stenci(2,counter2,counter) = tmpsteni(counter,2,nstencil)
        stenci(3,counter2,counter) = tmpsteni(counter,3,nstencil)
        stencr(counter2,counter) = gamma(counter2) & 
             * tmpstenr(counter,nstencil)
     END DO
  END DO



! Writing of rescue stencil when there has been the
! problem that no boundary condition could be found
! Die Geschwindigkeiten werden gesetzt:
! a. Geschwindigkeit der Dreiecke, die in der Zelle liegne
! b. Geschwindigkeit der umgebenden Dreiecke
! c. Geschwindigkeit Null
!-------------------------------------------------- 


  lopen = .false.
!jk test     IF ((nused.EQ.0).AND.(option.NE.2)) THEN
  lcoarse = .false.
  IF (x(i+1)-x(i).gt.10.0) lcoarse=.true.
  IF (y(j+1)-y(j).gt.10.0) lcoarse=.true.
  IF (z(k+1)-z(k).gt.10.0) lcoarse=.true.
  IF (((nused.EQ.0).AND.(option.NE.2)).OR.lcoarse) THEN
        ! kein Stencil gefunden und Geschwindigkeit
        IF (blocked(min(i+1,ni),j,k).EQ.1) lopen = .true.
        IF (blocked(max(i-1, 1),j,k).EQ.1) lopen = .true.
        IF (blocked(i,min(j+1,nj),k).EQ.1) lopen = .true.
        IF (blocked(i,max(j-1, 1),k).EQ.1) lopen = .true.
        IF (blocked(i,j,min(k+1,nk)).EQ.1) lopen = .true.
        IF (blocked(i,j,max(k-1, 1)).EQ.1) lopen = .true.
        IF (lopen) THEN
           ! Es gibt eine offene Nachbarzelle
           value = 0.0
           nvalue = 0
           IF (blocked(i,j,k).lt.0) then
              ! a. In der Zelle gibt es RB
              DO counter = 1,-blocked(i,j,k)
                 itri = tmpabi(i,j,k,counter)
                 value = value + topol(4,compon,itri)/ real(-blocked(i,j,k))
              END DO
              nvalue = -1
           ELSE
              ! b. Nachbarzellen abfragen
              DO it = max(i-1, 1),min(i+1,ni)
                 DO jt = max(j-1, 1),min(j+1,nj)
                    DO kt = max(k-1, 1),min(k+1,nk)
                       DO counter = 1,-blocked(it,jt,kt)
                          itri = tmpabi(it,jt,kt,counter)
                          value = value +  &
                               topol(4,compon,itri)/real(-blocked(it,jt,kt))
                       END DO
                       IF (blocked(it,jt,kt).le.-1) nvalue = nvalue + 1
                    END DO
                 END DO
              END DO
              IF (nvalue.ne.0) THEN
                 value = value / real(nvalue)
                 nresc1 = nresc1 + 1
              ELSE
                 ! c.
                 value = 0.0
!!$                 write(*,*) 'stencil: Geschw. auf Null, obwohl keine Nachbar RB gefunden'
!!$                 write(*,*) '         in Pos.:', i,j,k
                 nresc2 = nresc2 + 1
              END IF
           END IF

           ! Stencil bilden
           nused = 1
           IF (i.le.ni-nstenc-casep) THEN
              used(1) = 1
           ELSE IF (i.gt.nstenc+casep) THEN
              used(1) = 2
           ELSE IF (j.le.nj-nstenc-casep) THEN
              used(1) = 3
           ELSE IF (j.gt.nstenc+casep) THEN
              used(1) = 4
           ELSE IF (k.le.nk-nstenc-casep) THEN
              used(1) = 5
           ELSE IF (k.gt.nstenc+casep) THEN
              used(1) = 6
           END IF
           DO counter = 1, (nstenc+casep)
              stenci(1,1,counter) = i
              stenci(2,1,counter) = j
              stenci(3,1,counter) = k
           END DO

           stencr(1,:) = 0.0
           const(1) = value
  
           ! fuer folgende Ausgabe
           errout2 = 1
        END IF

        IF (errout2.GE.1) THEN

!!$           if (compon.eq.1) open(unit=31,file='nostenc_u.dat',form='FORMATTED',position='APPEND')
!!$           if (compon.eq.2) open(unit=31,file='nostenc_v.dat',form='FORMATTED',position='APPEND')
!!$           if (compon.eq.3) open(unit=31,file='nostenc_w.dat',form='FORMATTED',position='APPEND')
!!$           write(31,fmt="(3e16.7,I4,2X,e16.7,3(I7,2X))") x(i),y(j),z(k),nvalue,value,i,j,k
!!$           close(31)
           CALL stencil_output(blocked,ni,nj,nk,&
                i,j,k,error,direction,nused,used,&
                iproxi2,nstenc,nstenc_reduced,&
                stenci,stencr,const,distance,gridnum,&
                xc,yc,zc,compon)
        END IF
     END IF
!  END IF





END SUBROUTINE stencil

