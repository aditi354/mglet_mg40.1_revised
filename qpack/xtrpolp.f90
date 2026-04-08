SUBROUTINE xtrpolp(ni,nj,nk,nxpol,nutot,&
           norder,nlist,POINTER,pointer2,&
           infxpo,xpoli,xpolr,field)

  ! ********************************************************************** 
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !        modifications  07.11.03 (NP)
  !
  ! **********************************************************************

  IMPLICIT NONE
  INTEGER POINTER, pointer2
  INTEGER ni,nj,nk,nxpol,nutot,&
       norder,infxpo(2,nxpol),&
       xpoli(POINTER),&
       nlist(6)
  REAL xpolr(pointer2),field(ni*nj*nk)


  INTEGER i,j,k,ilist, ipont,&
       iorder,ip, jp, kp, npoint, npoint2,&
       idir, pntr, point, point2, is, ilist2
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
        DO iorder = 1, norder+1
           point2 = point2 + 1
           coef = xpolr(point2)
           field(i) = field(i) + coef*field(ip)
           ip = ip + idir
        END DO
     END DO
  END DO


END SUBROUTINE xtrpolp








