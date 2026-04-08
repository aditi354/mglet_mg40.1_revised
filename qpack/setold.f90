SUBROUTINE setold(ni,nj,nk,nxpol,nutot,&
     norder,nlist,POINTER, pointer2,infxpo,&
     xpoli,xpolr,field,oldsol,&
     maxpts)

  ! ********************************************************************** 
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !        modifications  07.11.03 (NP)
  !
  ! **********************************************************************

  IMPLICIT NONE
  INTEGER POINTER, pointer2, maxpts
  INTEGER ni,nj,nk,nxpol,nutot,norder,&
       infxpo(2,nxpol),xpoli(POINTER),&
       nlist(6)
  REAL xpolr(pointer2),field(ni*nj*nk),&
        oldsol(maxpts)

  INTEGER i,j,k,ilist,ipont,iorder,&
       ip,jp,kp,npoint, npoint2,&
       idir, pntr
  REAL coef

  IF (nxpol.GT.maxpts) THEN
     WRITE(*,*) 'more than maxpts extrapolated points', nxpol, maxpts
     STOP
  END IF

  DO ilist = 1,nxpol
     i = infxpo(1,ilist)
     oldsol(ilist) = field(i) 
  END DO

END SUBROUTINE setold








