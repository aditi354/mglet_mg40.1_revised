SUBROUTINE stencil_output(blocked,ni,nj,nk,i,j,k,error,direction,&
     nused,used,iproxi2,nstenc,nstenc_reduced,&
     stenci,stencr,const,distance,gridnum,&
     xc,yc,zc,compon)

  IMPLICIT NONE 
  INTEGER direction,ni,nj,nk,i,j,k,error(6,6),zaehler,nused
  INTEGER ausgabe,fehler_ja,counter,used(6),iproxi2(6),gridnum,ll
  INTEGER nstenc,nstenc_reduced,stenci(3,6,nstenc+1),&
          blocked(ni,nj,nk),compon
  REAL stencr(6,nstenc+1),const(6) 
  REAL distance(6),xc(ni),yc(nj),zc(nk)
  LOGICAL firstcall
  DATA firstcall /.TRUE./
  SAVE :: firstcall


  ! Fehlerbehandlung.
!  IF (firstcall) THEN
     OPEN(unit=10,file="QPack.error.in.stencil",&
          STATUS="REPLACE",POSITION="REWIND")
     WRITE(10,*) " "
     WRITE(10,*) "---Error in stencil creation---"
     WRITE(10,*) "   Error description:"
     WRITE(10,*) "      Nr.1 : limit_neg   / limit_pos     == 0"
     WRITE(10,*) "      Nr.2 : gdepth_neg  / gdepth_pos    == 0"
     WRITE(10,*) "      Nr.3 : nstenc_reduced_neg / ..pos  << nstenc"
     WRITE(10,*) "      Nr.4 : sbegin_neg  / sbegin_pos    == 0"
     WRITE(10,*) "      Nr.5 : xb_found                    == 0"
     WRITE(10,*) "      Nr.6 : space between xb too big        "
     WRITE(10,*) "---Error in stencil creation---"
     WRITE(10,*) " "
     CLOSE(10)
!  ENDIF
  firstcall = .FALSE.

  OPEN(unit=10,file="QPack.error.in.stencil",position="append")
   
  DO LL = 1,nused

     WRITE (unit=10,fmt='(A,I3,I2,I2,1X,A,1X,3(G10.4),6(I4))') &
          'GRID,DIR,COMPON',gridnum,used(LL),compon, &
          'COORDINATES', xc(i),yc(j),zc(k),&
          ERROR(LL,1),ERROR(LL,2),ERROR(LL,3),&
          ERROR(LL,4),ERROR(LL,5),ERROR(LL,6)
     
  END DO

  CLOSE(10)





!!$  ! Ausgabe eines stencils zum Vergleich,
!!$  ! wenn die Variable ausgabe = 1 ist
!!$  ausgabe=1
!!$  ! if ((i.eq.32).and.(j.eq.1).and.(k.eq.41)) then ausgabe=1
!!$  IF (ausgabe.EQ.1) THEN
!!$     WRITE(199,*) 'Stencil-Routine'
!!$     WRITE(199,*) '   i,j,k,blocked(i,j,k)',i,j,k,blocked(i,j,k)
!!$     DO direction = 1,nused
!!$        !        nstencil = used(direction)
!!$        DO counter = 1, nstenc
!!$           WRITE(199,*) '   i', stenci(1,direction,counter)
!!$           WRITE(199,*) '   j',  stenci(2,direction,counter)
!!$           WRITE(199,*) '   k',  stenci(3,direction,counter)
!!$           WRITE(199,*) '   koeff',  stencr(direction,counter)
!!$        END DO
!!$        WRITE(199,*) '   add konst',  const(direction) 
!!$        WRITE(199,*) '   distance',distance(used(direction))
!!$        WRITE(199,*) '   iproxi2, nused',iproxi2(used(direction)),nused
!!$        WRITE(199,*) '   ---'
!!$     END DO
!!$     WRITE(199,*) '   ' 
!!$     ! Fehlerbehandlung.
!!$     fehler_ja=0
!!$     IF(fehler_ja.EQ.1) THEN
!!$        WRITE(199,*) '-------------------'
!!$        WRITE(199,*) 'Konnte keine Stencils erstellen fuer:'
!!$        WRITE(199,*) 'Blocked(i,j,k)',  blocked(i,j,k)
!!$        WRITE(199,*) '   nstenc_reduced',nstenc_reduced
!!$        WRITE(199,*) '   ---' 
!!$        direction = 0
!!$        DO WHILE (direction.LT.6)
!!$           WRITE(199,*) '   Directions:', direction+1, direction+2
!!$           DO zaehler=1, 6
!!$              WRITE(199,*) '   Error :', zaehler
!!$              WRITE(199,*) '     +Richtung :',error(direction+1,zaehler),&
!!$                   '     = error(1,zaehler2)'
!!$              WRITE(199,*) '     -Richtung :',error(direction+2,zaehler),&
!!$                   '     = error(2,zaehler2)'
!!$           END DO
!!$           direction = direction + 2
!!$        END DO
!!$        WRITE(199,*) '-------------------'
!!$     END IF
!!$  END IF
!!$



   END SUBROUTINE stencil_output

