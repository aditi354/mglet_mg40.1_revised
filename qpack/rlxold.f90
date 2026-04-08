SUBROUTINE rlxold(ni,nj,nk,nxpol,nutot,&
     norder,nlist,POINTER, pointer2,&
     infxpo,xpoli,xpolr,field,&
     oldsol,rlx,maxpts)

  ! **********************************************************************
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !        modifications  07.11.03 (NP)
  !
  ! **********************************************************************

  IMPLICIT NONE
  INTEGER POINTER, pointer2,maxpts
  INTEGER ni,nj,nk,nxpol,nutot,norder,&
       infxpo(2,nxpol),xpoli(POINTER),&
       nlist(6)
  REAL xpolr(pointer2),field(ni*nj*nk),&
       oldsol(maxpts), rlx


  INTEGER i,j,k,ilist,ipont,iorder,&
       ip,jp,kp,npoint,npoint2,&
       idir, pntr
  REAL coef

  IF (nxpol.GT.maxpts) THEN
     WRITE(*,*) 'more than maxpts extrapolated points', nxpol,maxpts
     STOP
  END IF

  IF (rlx.NE.0) THEN
     ! Calculating the solution with relaxation
     DO ilist = 1,nxpol
        field(infxpo(1,ilist)) = &
        (1-rlx)*oldsol(ilist) + rlx*field(infxpo(1,ilist)) 
     END DO
  ELSE
     DO ilist = 1,nxpol
        i = infxpo(1,ilist)
        field(i) = oldsol(ilist) 
     END DO
  END IF

END SUBROUTINE rlxold








