SUBROUTINE blockingadd(blocked,nflow,flow,&
     ni,nj,nk,itermax,bconds,redord)

  ! **********************************************************************
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !
  !   modifications:
  !   04.09.2003  modifications N.Peller - routine
  !               derived from blocking.f90
  !jk 11. 8.2005  Fluid : blocked=1
  !               RB    : blocked < 0
  !               additional blocking: blocked = -999
  !jk 11. 8.2005  Alle zusaetzlichen blocked-Felder gespart
  ! **********************************************************************
  ! TODO: jk 11. 8.2005 zweimal Luecken schliessen ist ueberfluessig.

  IMPLICIT NONE
  INTEGER ni,nj,nk,blocked(ni,nj,nk),&
       nflow,flow(3,2,nflow),bconds(6),redord
  ! lokale Variablen
  INTEGER iter,&
       ncount,icount,nnmoda,nnmodb,nnmodc,&
       i,j,k,kmin,kmax,imin,imax,&
       jmin,jmax,nmoda,nmodb,nmodc,imod,&
       itermax
!!$!!! Rausschreiben des Fuellens zur Fehlersuche
!!$  character ctext*80, filen*17


  ! now painting the exterior of the body
  ! Den ersten Teil des Gebietes ausserhalb
  ! des Koerpers schreiben
  DO icount = 1,nflow
     DO k = flow(3,1,icount),flow(3,2,icount)
        DO j = flow(2,1,icount),flow(2,2,icount)
           DO i = flow(1,1,icount),flow(1,2,icount)
              IF (blocked(i,j,k).LT.0) &
                   STOP 'blockingadd: painting area schneidet RB.'
              blocked(i,j,k) = 1
           END DO
        END DO
     END DO
  END DO



  ! Ab hier werden ausgehend vom zuvor bestimmten Blocked
  ! Feld mit Geometrieschnittpunkten die Regionen 
  ! innerhalb und ausserhalb des Koerpers 
  ! bestimmt. --> Fuellalgorithmus 

  ! zuerst: kleine Luecken schliessen, damit nichts rausfliesst

  nnmoda = 0
  DO imod = 1,100
     nmoda = 0
     DO k = 2,nk-1
        DO j = 1,nj
           DO i = 1,ni
              IF (blocked(i,j,k).EQ.0) THEN
                 IF (blocked(i,j,k+1).LT.0.AND.&
                      blocked(i,j,k-1).LT.0) THEN
                    blocked(i,j,k) = -999
                    nmoda = nmoda+1
                 END IF
              END IF
           END DO
        END DO
     END DO
     DO k = 1,nk
        DO j = 2,nj-1
           DO i = 1,ni
              IF (blocked(i,j,k).EQ.0) THEN
                 IF (blocked(i,j+1,k).LT.0.AND.&
                      blocked(i,j-1,k).LT.0) THEN
                    blocked(i,j,k) = -999
                    nmoda = nmoda+1
                 END IF
              END IF
           END DO
        END DO
     END DO
     DO k = 1,nk
        DO j = 1,nj
           DO i = 2,ni-1
              IF (blocked(i,j,k).EQ.0) THEN
                 IF (blocked(i-1,j,k).LT.0.AND.&
                      blocked(i+1,j,k).LT.0) THEN
                    blocked(i,j,k) = -999
                    nmoda = nmoda+1
                 END IF
              END IF
           END DO
        END DO
     END DO

     IF (nmoda.EQ.0 .AND. nmodb.EQ.0) EXIT
     nnmoda = nnmoda + nmoda
     IF (imod.EQ.100) THEN
        WRITE(*,*) 'ERROR: In 100 loops it was not possible '
        WRITE(*,*) '       to close all cells lying between'
        WRITE(*,*) '       closed cells (A)'
     END IF
  END DO


  ! Schau in Sternform um die nicht-geblockte Zelle
  ! Wenn die betrachtete Zelle nicht geblockt
  ! ist, dann setzte auch die Nachbarzelle
  ! auf nicht-geblockt, wenn diese keine
  ! Geometrie enthaelt 

  iter = 0
100 CONTINUE
  iter = iter + 1
  ncount = 0

  DO k = 1,nk-1
  DO j = 1,nj
  DO i = 1,ni
     IF (blocked(i,j,k).EQ.1) THEN
        IF (blocked(i,j,k+1).EQ.0) THEN
           blocked(i,j,k+1) = 1
           ncount = ncount+1
        END IF
     END IF
  END DO
  END DO
  END DO

  DO k = 2,nk
  DO j = 1,nj
  DO i = 1,ni
     IF (blocked(i,j,k).EQ.1) THEN
        IF (blocked(i,j,k-1).EQ.0) THEN
           blocked(i,j,k-1) = 1
           ncount = ncount+1
        END IF
     END IF
  END DO
  END DO
  END DO

  DO k = 1,nk
  DO j = 1,nj-1
  DO i = 1,ni
     IF (blocked(i,j,k).EQ.1) THEN
        IF (blocked(i,j+1,k).EQ.0) THEN
           blocked(i,j+1,k) = 1
           ncount = ncount+1
        END IF
     END IF
  END DO
  END DO
  END DO


  DO k = 1,nk
  DO j = 2,nj
  DO i = 1,ni
     IF (blocked(i,j,k).EQ.1) THEN
        IF (blocked(i,j-1,k).EQ.0) THEN
           blocked(i,j-1,k) = 1
           ncount = ncount+1
        END IF
     END IF
  END DO
  END DO
  END DO

  DO k = 1,nk
  DO j = 1,nj
  DO i = 1,ni-1
     IF (blocked(i,j,k).EQ.1) THEN
        IF (blocked(i+1,j,k).EQ.0) THEN
           blocked(i+1,j,k) = 1
           ncount = ncount+1
        END IF
     END IF
  END DO
  END DO
  END DO

  DO k = 1,nk
  DO j = 1,nj
  DO i = 2,ni
     IF (blocked(i,j,k).EQ.1) THEN
        IF (blocked(i-1,j,k).EQ.0) THEN
           blocked(i-1,j,k) = 1
           ncount = ncount+1
        END IF
     END IF
  END DO
  END DO
  END DO

!!$!!! Rausschreiben des Fuellens zur Fehlersuche
!!$!  if (modulo(iter,1).eq.1) then
!!$     !   write(filen,'(a12,i5)') 'blocking.BP_00000'
!!$     write(filen,'(a12,i5.5)') 'blocking.BP_',3*iter
!!$     open(19,file=filen)
!!$     ctext =  'BP'
!!$     write(19,'(a80)') ctext
!!$     ctext =  'part'
!!$     write(19,'(a80)') ctext
!!$     write(19,'(i10)') 1
!!$     ctext = 'block'
!!$     write(19,'(a80)') ctext
!!$     
!!$     DO k = 1,nk
!!$        DO j = 1,nj
!!$           DO i = 1,ni
!!$              write (19,'(e12.5)') real(blocked (i,j,k))
!!$           END DO
!!$        END DO
!!$     END DO
!!$     close(19)
!!$!  end if
!!$!!! Ende Rausschreiben

  IF (ncount.NE.0) THEN
     IF (iter.GT.itermax) THEN
        WRITE(*,*) '   ERROR: could not fill the body'
        STOP
     ELSE
        GOTO 100
     END IF
  END IF


  
  ! performing additional blocking
  ! Nachdem der Stern Algorithmus das 
  ! Stroemungsfeld bestimmt hat muessen jetzt
  ! noch Luecken aufgefuellt werden
  WRITE(*,*) '   START ADDITIONAL BLOCKING'
  WRITE(*,*) '      CASE A: in one direction both'
  WRITE(*,*) '              neighbours are blocked'
  WRITE(*,*) '      CASE B: in one direction only'
  WRITE(*,*) '              one neighbourcell is open,'
  WRITE(*,*) '              the next is blocked'

  nnmodb = 0
  nnmodc = 0


  ! Kleine Luecken fuellen, um Zwickel zu vermeiden
  DO imod = 1,100
     nmoda = 0
     nmodb = 0
     nmodc = 0

     ! einser Luecken
     DO k = 2,nk-1
        DO j = 1,nj
           DO i = 1,ni
              IF (blocked(i,j,k).EQ.1) THEN
                 IF (blocked(i,j,k+1).LE.0.AND.&
                      blocked(i,j,k-1).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmoda = nmoda+1
                 END IF
              END IF
           END DO
        END DO
     END DO
     DO k = 1,nk
        DO j = 2,nj-1
           DO i = 1,ni
              IF (blocked(i,j,k).EQ.1) THEN
                 IF (blocked(i,j+1,k).LE.0.AND.&
                      blocked(i,j-1,k).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmoda = nmoda+1
                 END IF
              END IF
           END DO
        END DO
     END DO
     DO k = 1,nk
        DO j = 1,nj
           DO i = 2,ni-1
              IF (blocked(i,j,k).EQ.1) THEN
                 IF (blocked(i-1,j,k).LE.0.AND.&
                      blocked(i+1,j,k).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmoda = nmoda+1
                 END IF
              END IF
           END DO
        END DO
     END DO

     ! zweier Luecken
     DO k = 2,nk-2
        DO j = 1,nj
           DO i = 1,ni
              IF (blocked(i,j,k).EQ.1) THEN
                 IF (blocked(i,j,k+2).LE.0.AND.&
                      blocked(i,j,k-1).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmodb = nmodb+1
                 END IF
                 IF (blocked(i,j,k+1).LE.0.AND.&
                      blocked(i,j,k-2).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmodb = nmodb+1
                 END IF
              END IF
           END DO
        END DO
     END DO
     DO k = 1,nk
        DO j = 2,nj-2
           DO i = 1,ni
              IF (blocked(i,j,k).EQ.1) THEN
                 IF (blocked(i,j+2,k).LE.0.AND.&
                      blocked(i,j-1,k).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmodb = nmodb+1
                 END IF
                 IF (blocked(i,j+1,k).LE.0.AND.&
                      blocked(i,j-2,k).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmodb = nmodb+1
                 END IF
              END IF
           END DO
        END DO
     END DO
     DO k = 1,nk
        DO j = 1,nj
           DO i = 2,ni-2
              IF (blocked(i,j,k).EQ.1) THEN
                 IF (blocked(i+2,j,k).LE.0.AND.&
                      blocked(i-1,j,k).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmodb = nmodb+1
                 END IF
                 IF (blocked(i+1,j,k).LE.0.AND.&
                      blocked(i-2,j,k).LE.0) THEN
                    blocked(i,j,k) = -999
                    nmodb = nmodb+1
                 END IF
              END IF
           END DO
        END DO
     END DO

     ! Umgebung der Randbedingungen:
     
     ! OP1
     IF (bconds(5).eq.3) THEN
        DO j = 1,nj
           DO i = 1,ni
              if (blocked(i,j,5).LE.0) then
                 if (blocked(i,j,4).gt.0) then
                    blocked(i,j,4) = -999
                    nmodc = nmodc+1
                 end if
              end if
              if (blocked(i,j,4).LE.0) then
                 if (blocked(i,j,3).gt.0) then
                    blocked(i,j,3) = -999
                    nmodc = nmodc+1
                 end if
              end if
              if (blocked(i,j,3).LE.0) then
                 ! fuer Steuerung in BLOCKBP
                 if (blocked(i,j,2).gt.0) then
                    blocked(i,j,2) = -999
                 end if
              end if
           END DO
        END DO
     END IF
     IF (bconds(6).eq.3) THEN
        DO j = 1,nj
           DO i = 1,ni
              if (blocked(i,j,nk-4).LE.0) then
                 if (blocked(i,j,nk-3).gt.0) then
                    blocked(i,j,nk-3) = -999
                    nmodc = nmodc+1
                 end if
              end if
              if (blocked(i,j,nk-3).LE.0) then
                 if (blocked(i,j,nk-2).gt.0) then
                    blocked(i,j,nk-2) = -999
                    nmodc = nmodc+1
                 end if
              end if
              if (blocked(i,j,nk-2).LE.0) then
                 ! fuer Steuerung in BLOCKBP
                 if (blocked(i,j,nk-1).gt.0) then
                    blocked(i,j,nk-1) = -999
                 end if
              end if
           END DO
        END DO
     END IF
     IF (bconds(2).eq.3) THEN
        DO k = 1,nk
           DO j = 1,nj
              if (blocked(ni-4,j,k).LE.0) then
                 if (blocked(ni-3,j,k).gt.0) then
                    blocked(ni-3,j,k) = -999
                    nmodc = nmodc+1
                 end if
              end if
              if (blocked(ni-3,j,k).LE.0) then
                 if (blocked(ni-2,j,k).gt.0) then
                    blocked(ni-2,j,k) = -999
                    nmodc = nmodc+1
                 end if
              end if
              if (blocked(ni-2,j,k).LE.0) then
                 ! fuer Steuerung in BLOCKBP
                 if (blocked(ni-1,j,k).gt.0) then
                    blocked(ni-1,j,k) = -999
                 end if
              end if
           END DO
        END DO
     END IF
     IF (bconds(1).eq.3) STOP 'blockingadd, OP1 FR nicht implementiert'
     IF (bconds(3).eq.3) STOP 'blockingadd, OP1 RI nicht implementiert'
     IF (bconds(4).eq.3) STOP 'blockingadd, OP1 LE nicht implementiert'

     IF (nmoda.EQ.0 .AND. nmodb.EQ.0) EXIT
     nnmoda = nnmoda + nmoda
     nnmodb = nnmodb + nmodb
     nnmodc = nnmodc + nmodc
     IF (imod.EQ.100) THEN
        WRITE(*,*) 'ERROR: In 100 loops it was not possible '
        WRITE(*,*) '       to close all cells lying between'
        WRITE(*,*) '       closed cells'
     END IF
  END DO
  WRITE(*,*) '      LOOPS FOR ADDITIONAL BLOCKING: ',imod
  WRITE(*,*) '      CELLS BLOCKED: CASE A:', nnmoda, 'CASE B', nnmodb, 'CASE C', nnmodc




  ! COUNTING HOW MANY BLOCKED CELLS
  imin = 999e6
  imax = -999e6
  jmin = imin
  jmax = imax
  kmin = imin
  kmax = imax
  ncount = 0
  DO k = 1,nk
     DO j = 1,nj
        DO i = 1,ni
           IF (blocked(i,j,k).LE.0) THEN
              ncount = ncount + 1
              IF (i.LT.imin) imin = i
              IF (i.GT.imax) imax = i
              IF (j.LT.jmin) jmin = j
              IF (j.GT.jmax) jmax = j
              IF (k.LT.kmin) kmin = k
              IF (k.GT.kmax) kmax = k
           END IF
        END DO
     END DO
  END DO
  WRITE(*,*) '   THE PRESSURE GRID HAS ',ncount,' BLOCKED CELLS'
  WRITE(*,*) '      FROM I = ',imin, ' TO ', imax
  WRITE(*,*) '      FROM J = ',jmin, ' TO ', jmax
  WRITE(*,*) '      FROM K = ',kmin, ' TO ', kmax
  WRITE(*,*) '   RATIO BLOCKED/ALL   :   ', float(ncount)/float(nk*nj*ni)




END SUBROUTINE blockingadd
