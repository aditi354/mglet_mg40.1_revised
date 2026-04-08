SUBROUTINE blocking_xtrpol(ni,nj,nk,x,y,z,ntopol,topol,maxblock,norder,&
     ndepth,option,idnum,xc,yc,zc,ntrimax, blocked,&
     compon,maccur,itermax,sump1,sump2,ilim,&
     whatipol,gridnum,bconds,redord,nflow,flow,istenc)

  ! **********************************************************************
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !
  !   modifications:
  !   13.02.2001  JK  loop over cells given to subroutine blocking2
  !   26.02.2001  JK  additional blocking modified:
  !   04.09.2003  modifications N.Peller 
  !   15. 7.2005  JK Schleife fuer Geschw. Stencils von nsetnc+1 auf nsetnc reduziert
  !   11. 8.2005  jk Moeglichkeit mehrerer paint-Gebiete (flow) wiedereingefuehrt
  ! **********************************************************************
  ! TODO : jk 11. 8.2005 Skalar muss noch implementiert werden

  IMPLICIT NONE
  INTEGER ni,nj,nk,ntopol,maxblock,face,norder,ndepth,&
          gridnum,nxpol,infxpo(2,maxblock),xpoli(12*maxblock),&
          maxsurf, option,idnum, ntrimax, used(6),&
          compon,ib,jb,kb,itermax,whatipol,&
          bconds(6),blocked(ni,nj,nk),redord,nflow,flow(3,2,nflow),istenc
  REAL x(ni+1),y(nj+1),z(nk+1),topol(4,3,ntopol),&
       blockedr(ni,nj,nk),xpolr(6*(norder+2)*maxblock),&
       xc(ni),yc(nj),zc(nk),maccur
  !lokale Variablen
  INTEGER nblocked,i,j,k,nused,stenci(3,6,norder+1), ncount,&
          nutot,tmpabi(ni,nj,nk,ntrimax),&
          isten, pointer1,dummy,&
          nlist(6), nlistot,sump1,sump2,icount2,icount3,&
          nchar, tabli, tablr, iproxi2(6),&
          ko,jo,io, pointer2,&
          ncount2, ie,je,ke, ilim(3,2),&
          itri, npex ,nredord, nresc1, nresc2
  REAL const(6),stencr(6,norder+1),cratio 
  CHARACTER(len=10) :: prefix
  CHARACTER(len=4) :: suffix
  CHARACTER(len=4) :: number
  CHARACTER (len=18) :: filout
  prefix = 'BLOCK_XTR_'
  suffix = '.DAT'


  WRITE(*,*) '   ALLOCATING INTEGERS FOR STENCILS :',&
             ni*nj*nk*ntrimax+14*maxblock 
  WRITE(*,*) '   ALLOCATING REALS FOR STENCILS    :',&
             6*(norder+2)*maxblock


  nutot = 0
  nxpol = 0
  dummy = 0
  
  nredord = 0
  nresc1  = 0
  nresc2  = 0

  IF (option.EQ.2) THEN

     !WRITE(*,*) '---BLOCKING FOR PRESSURE---'
     ! Here blocking for the pressure
     ! Hier wird zuerst nur das Blocking Feld
     ! mit den Geometriezellen geschrieben.
     if (istenc.eq.-1) then
        ! istenc = 0,-1: nur BP-Feld erzeugen
        CALL blocking(ni,nj,nk,x,y,z,ntopol,topol,&
             blocked,tmpabi,ntrimax,&
             xc,yc,zc,maccur,itermax,ilim,option)
     end if
     if (istenc.le.0) then
        ! Jetzt ausgehend von den Zellen mit
        ! Geometrie das Innere und Aeussere
        ! des Koerpers auffuellen.
        CALL blockingadd(blocked,nflow,flow,&
             ni,nj,nk,itermax,bconds,redord)   
       return
     end if
     pointer1 = 0
     pointer2 = 0
     DO k = max(3,ilim(3,1)),min(nk-2,ilim(3,2))
        DO j = max(3,ilim(2,1)),min(nj-2,ilim(2,2))
           DO i = max(3,ilim(1,1)),min(ni-2,ilim(1,2))
              blocked(i,j,k) = MAX(blocked(i,j,k), 0)
              IF (blocked(i,j,k).EQ.0) THEN
                 CALL stencil(ni,nj,nk,x,y,z,tmpabi,blocked,&
                      ntopol,topol,i,j,k,norder,nused,const,&
                      stenci,stencr,ndepth,xc,yc,zc, ntrimax,&
                      ni,nj,nk, used, iproxi2,&
                      compon,maccur,whatipol,gridnum,bconds,&
                      redord,option,nredord, nresc1, nresc2)
          
                 IF (nused.NE.0) THEN
                    nxpol = nxpol + 1
                    IF (nxpol.GT.maxblock) THEN
                       WRITE(*,*) 'MAXBLOCK EXCEEDED'
                       STOP
                    END IF
                    infxpo(1,nxpol) = 1+(k-1)+(j-1)*nk+(i-1)*nj*nk
                    infxpo(2,nxpol) = nused    
                    DO ncount = 1,nused
                       ko = stenci(3,ncount,1)
                       jo = stenci(2,ncount,1)
                       io = stenci(1,ncount,1)
                       pointer1 = pointer1 + 1
                       IF (used(ncount).LE.2) THEN
                          xpoli(pointer1) = nj*nk
                       ELSE IF (used(ncount).LE.4) THEN
                          xpoli(pointer1) = nk
                       ELSE IF (used(ncount).LE.6) THEN
                          xpoli(pointer1) = 1
                       END IF
                       xpoli(pointer1) = xpoli(pointer1)*iproxi2(used(ncount))
                       pointer1 = pointer1 + 1
                       xpoli(pointer1) = 1+(ko-1)+(jo-1)*nk+(io-1)*nj*nk
                       nutot = nutot + 1                      
                       DO isten = 1,norder+1  
                          pointer2 = pointer2 + 1
                          xpolr(pointer2) = stencr(ncount,isten)
                       END DO
                    END DO
                 END IF

              END IF
           END DO
        END DO
     END DO
  END IF
 


  ! Dieser Abschnitt erstellt das Blocked-Feld
  ! fuer die Geschwindigkeitszellen. Das Feld wird
  ! aus dem Blocked-Feld der Druckzellen
  ! abgeleitet.
  IF (option.NE.2) THEN

     tmpabi = 0
     !WRITE(*,*) '---BLOCKING FOR VELOCITIES---'
     CALL blocking(ni,nj,nk,x,y,z,ntopol,topol,&
          blocked,tmpabi,ntrimax,&
          xc,yc,zc,maccur,itermax,ilim,option)

     pointer1 = 0
     pointer2 = 0
     DO k = max(3,ilim(3,1)),min(nk-2,ilim(3,2))
        DO j = max(3,ilim(2,1)),min(nj-2,ilim(2,2))
           DO i = max(3,ilim(1,1)),min(ni-2,ilim(1,2))
               IF (blocked(i,j,k).NE.1) THEN
                 CALL stencil(ni,nj,nk,x,y,z,tmpabi,blocked,&
                      ntopol,topol,i,j,k,norder,nused,const,&
                      stenci,stencr,ndepth,xc,yc,zc, ntrimax,&
                      ni,nj,nk, used, iproxi2,&
                      compon,maccur,whatipol,gridnum,bconds,&
                      redord,option,nredord, nresc1, nresc2) 

                 IF (nused.NE.0) THEN
                    nxpol = nxpol + 1
                    IF (nxpol.GT.maxblock) THEN
                       WRITE(*,*) 'MAXBLOCK EXCEEDED'
                       STOP
                    END IF
                    infxpo(1,nxpol) = 1+(k-1)+(j-1)*nk+(i-1)*nj*nk
                    infxpo(2,nxpol) = nused    
                    DO ncount = 1,nused
                       ko = stenci(3,ncount,1)
                       jo = stenci(2,ncount,1)
                       io = stenci(1,ncount,1)
                       pointer1 = pointer1 + 1
                       ! Used(ncount) gibt an welcher der 6 moeglichen
                       ! verwendet wird = used(nused)
                       IF (used(ncount).LE.2) THEN
                          xpoli(pointer1) = nj*nk
                       ELSE IF (used(ncount).LE.4) THEN
                          xpoli(pointer1) = nk
                       ELSE IF (used(ncount).LE.6) THEN
                          xpoli(pointer1) = 1
                       END IF
                       ! multipliziert mit 1 oder -1 je nach iproxi2
                       ! iproxi2 ist also nur die Stencil-Richtung
                       xpoli(pointer1) = xpoli(pointer1)*iproxi2(used(ncount))
                       pointer1 = pointer1 + 1
                       ! Umschreiben in andere Laufvariable
                       xpoli(pointer1) = 1+(ko-1)+(jo-1)*nk+(io-1)*nj*nk
                       ! Anzahl der used stencils um 1 erhoehen
                       nutot = nutot + 1                      
                       DO isten = 1,norder
                          pointer2 = pointer2 + 1
                          xpolr(pointer2) = stencr(ncount,isten)
                       END DO
                       pointer2 = pointer2 + 1
                       xpolr(pointer2) = const(ncount)                        
                       !write(100,*) ncount,const(ncount) 
                    END DO
                 END IF
              END IF
           END DO
        END DO
     END DO
  END IF

  WRITE(*,*) '   USING REDUCED ORDER                : ',nredord
  WRITE(*,*) '   USING VELOCITY OF NEAREST TRIANGLE : ',nresc1
  WRITE(*,*) '   USING VELOCITY = 0                 : ',nresc2

  sump1 = sump1 + pointer1 + 2*nxpol
  sump2 = sump2 + pointer2
  WRITE(*,*) '   POINTS TO EXTRAPOLATE :',nxpol
  WRITE(*,*) '   REQUIERS INTEGER NUMBERS : ', pointer1 + 2*nxpol
  WRITE(*,*) '   REQUIERS REAL NUMBERS    : ', pointer2
  WRITE(*,*) '   FOR A TOTAL OF : ',sump1,' AND ',sump2



  ! Hier kann das momentane Blocked Feld 
  ! in ein file ausgegeben werden, das
  ! von gnuplot gelesen werden kann
  ! call blocked_gnuplot(blocked,ni,nj,nk,x,&
  !                       y,z,xc,yc,zc,maccur) 



  ! In dieser Schleife werden alle 
  ! negativen blocking Werte auf Null 
  ! gessetzt. Fuer MGlet ist das wichtig
  DO k = 1,nk
     DO j = 1,nj
        DO i = 1,ni
           blocked(i,j,k) = MAX(0,blocked(i,j,k))
        END DO
     END DO
  END DO


  ! Hier wird der Filename in Abhängigkeit der idnum
  ! bestimmt. Die idnum ist bestimmt durch
  ! idnum = 10 x gridnum + 2
  WRITE(number(1:4),'(i4)') idnum
  filout = prefix//number//suffix
  DO icount2 = 1,8
     IF (filout(10+1:10+1).EQ.' ') THEN
        DO icount3 = 1,8-icount2
           filout(10+icount3:10+icount3)=filout(10+icount3+1:10+icount3+1)
        END DO
     END IF
  END DO
  ! Hier wird die Anzahl der Stellen fuer den
  ! Namen der Datei bestimmt. 10 --> 2
  ! 100 --> 3 Stellen + 4 Stellen fuer '.dat'
  nchar = INT(LOG10(float(idnum))+1) + 4


  ! Jetzt wird das File mit dem zuvor
  ! bestimmten Filenamen geoeffnet
  ! und die stencil-Koeffizienten 
  ! geschrieben. Das blocked-Feld muss 
  ! als real Feld gespeichert werden.
  blockedr = REAL(blocked)
  OPEN (unit=13,file=filout(1:10+nchar),form='UNFORMATTED')
  WRITE(13) option
  WRITE(13) nk,nj,ni
  WRITE(13) nxpol
  WRITE(13) pointer1
  WRITE(13) pointer2
  WRITE(13) norder
  WRITE(13) ((infxpo(j,i),j=1,2),i=1,nxpol),1.0
  WRITE(13) (xpoli(i),i=1,pointer1),1.0
  WRITE(13) (xpolr(i),i=1,pointer2) ,1.0
  IF (option.EQ.2) THEN
     WRITE(13) (((blockedr(i,j,k),k=1,nk),j=1,nj),i=1,ni)
  END IF
  CLOSE(13)




  write(*,*) 'DEALLOCATING ',ni*nj*nk*ntrimax+14*maxblock ,' INTEGERS FOR STENCILS'
  write(*,*) 'DEALLOCATING ',6*(norder+1)*maxblock,' REALS FOR STENCILS'



END SUBROUTINE blocking_xtrpol













