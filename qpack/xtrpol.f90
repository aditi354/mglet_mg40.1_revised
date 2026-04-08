SUBROUTINE xtrpol(ni,nj,nk,nxpol,nutot,&
           norder,nlist,POINTER, pointer2,&
           infxpo,xpoli,xpolr,field,value)

  ! ********************************************************************** 
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !        modifications  07.11.03 (NP)
  !                       15.7.05  Schleife ueber iorder reduziert jk
  ! **********************************************************************

  IMPLICIT NONE
  INTEGER POINTER, pointer2
  INTEGER ni,nj,nk,nxpol,nutot,norder,&
       infxpo(2,nxpol),xpoli(POINTER),&
       nlist(6)
  REAL xpolr(pointer2),field(ni*nj*nk), value

  INTEGER i,j,k,ilist,ipont,iorder,&
          ip,jp,kp,npoint,npoint2,&
          idir,pntr,is,point,point2,ilist2
  REAL coef


  point = 0
  point2 = 0
  DO ilist = 1,nxpol
     i = infxpo(1,ilist)
     field(i) = 0.0
     ilist2 = infxpo(2,ilist)
     DO ipont = 1, ilist2
        point = point + 1
        idir = xpoli(point)
        point = point + 1
        is = xpoli(point)
        ip = is
        DO iorder = 1, norder
           point2 = point2 + 1
           coef = xpolr(point2)
           field(i) = field(i) + coef*field(ip)
           ip = ip + idir
        END DO
        point2 = point2 + 1
        coef = xpolr(point2)
        field(i) = field(i) + coef
     END DO
  END DO


END SUBROUTINE xtrpol








