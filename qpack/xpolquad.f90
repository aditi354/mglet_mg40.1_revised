SUBROUTINE xpolquad(ni,nj,nk,nvirtual,&
              gridnum,comp,field, relaxkind)

  ! ********************************************************************* 
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !        modifications  07.11.03 (NP)
  !
  ! **********************************************************************

  IMPLICIT NONE
  INTEGER gridnum,comp,ni,nj,nk,&
          nvirtual, memax
  REAL field(ni+2*nvirtual,nj+2*nvirtual,nk+2*nvirtual),&
       value, relax

  INTEGER, PARAMETER :: maxmax = 100
  INTEGER nlistmax,maxpts ,maxgrd,narea,&
       i,j,k,nip,njp,nkp,nslist,slist,&
       size,POINTER,infxpo,xpoli,ipoint1,&
       ipoint2,ipoint3,ipoint4,memax1,memax2,memax3,memax4,&
       option,nutot,nxpol,norder,nlist(6),nchar,icount,&
       icount2,memu1,memu2,memu3,memu4,list,ipoint5,icount3,&
       off1,off2,off3,offset,xpoli2,xpoli3,xpoli4,&
       infxpo2,infxpo3,infxpo4,pointer3,pointer2,iold
  REAL xpolr,dummy,oldsol,xpolr2,&
       xpolr3,xpolr4,relaxkind,rlxfac
  LOGICAL firstcall
  CHARACTER(len=10) :: prefix
  CHARACTER(len=4) :: suffix
  CHARACTER(len=4) :: number
  CHARACTER (len=18) :: filout
  ALLOCATABLE xpoli(:), xpolr(:), oldsol(:,:),&
              slist(:,:),SIZE(:,:),POINTER(:,:)
  DATA firstcall /.TRUE./
  SAVE :: nslist, slist, firstcall, memu1,&
          POINTER, size, xpoli, xpolr, memu2, memu3,&
          memax3,memu4,memax4,list,oldsol,nlistmax,&
          rlxfac,maxpts,memax1,memax2


  offset = 0
  IF (firstcall) THEN
     OPEN (unit=13,file="QPack.dat",form='formatted')
     READ(13,*)
     READ(13,*) narea
     DO i=1,narea
        READ(13,*)
        READ(13,*)
     END DO
     READ(13,*)
     READ(13,*)
     READ(13,*)
     READ(13,*)
     READ(13,*)
     READ(13,*)
     READ(13,*)
     READ(13,*) nlistmax
     READ(13,*) memax1
     READ(13,*) memax2
     READ(13,*) maxpts
     READ(13,*) rlxfac
     CLOSE(13)

     ALLOCATE (xpoli(memax1))
     ALLOCATE (xpolr(memax2))
     WRITE(*,*) 'ALLOCATING ',memax1,memax2,' integers and reals in XpolQuad'
     ALLOCATE (oldsol(maxpts,nlistmax))
     ALLOCATE (slist(5,maxmax),SIZE(10,maxmax),POINTER(5,maxmax))
     firstcall = .FALSE.
     size = 0
     nslist = 0
     slist = 0
     POINTER = 0
     memu1 = 0
     memu2 = 0
     memu3 = 0
     memu4 = 0
     oldsol = 0
  END IF

  ! check if in memory
  IF (gridnum.GT.maxmax) THEN
     WRITE(*,*) 'igrid is greater than maxmax', gridnum, maxmax
     STOP
  END IF

  ipoint1 = slist(comp,gridnum)
  IF (ipoint1.NE.0) THEN
     nip = SIZE(1,ipoint1)
     njp = SIZE(2,ipoint1)
     nkp = SIZE(3,ipoint1)
     option = SIZE(4,ipoint1)
     nxpol = SIZE(5,ipoint1)
     pointer3 = SIZE(6,ipoint1)
     norder = SIZE(7,ipoint1)
     pointer2 = SIZE(10, ipoint1)
     ipoint3 = SIZE(8,ipoint1)
     ipoint4 = SIZE(9,ipoint1)
  ELSE
     ! not in memory
     prefix = 'BLOCK_XTR_'
     suffix = '.DAT'
     WRITE(number(1:4),'(i4)') 10*gridnum+comp
     filout = prefix//number//suffix
     DO icount2 = 1,8
        IF (filout(10+1:10+1).EQ.' ') THEN
           DO icount3 = 1,8-icount2
              filout(10+icount3:10+icount3) = &
              filout(10+icount3+1:10+icount3+1)
           END DO
        END IF
     END DO
     nchar = INT(LOG10(float(10*gridnum+comp))+1) + 4
     OPEN (unit=13,file=filout(1:10+nchar),form='UNFORMATTED')

     READ(13) option
     READ(13) nip,njp,nkp
     READ(13) nxpol
     READ(13) pointer3
     READ(13) pointer2
     READ(13) norder
     nslist = nslist + 1
     IF (nslist.GT.nlistmax) THEN
        WRITE(*,*) 'maximum number of grid exceeded'
        STOP
     END IF
     slist(comp,gridnum) = nslist
     ipoint1 = nslist

     IF (memu1+2*nxpol+pointer3.LE.memax1) THEN  
        ipoint3 = memu1+1
        READ(13) (xpoli(i),i=memu1+1,memu1+2*nxpol)
        READ(13) (xpoli(i),i=memu1+2*nxpol+1,memu1+2*nxpol+pointer3)
        memu1 = memu1 + 2*nxpol+pointer3
     ELSE
        WRITE(*,*) 'memax1 too small', memu1+2*nxpol+pointer3
        WRITE(*,*) 'idnum : ',10*gridnum+comp
        STOP
     END IF
     IF (memu2+pointer2.LE.memax2) THEN  
        ipoint4 = memu2+1
        READ(13) (xpolr(i),i=memu2+1,memu2+pointer2) 
        memu2 = memu2 +pointer2
     ELSE
        WRITE(*,*) 'memax2 too small', memu2+pointer2
        STOP
     END IF
     CLOSE(13)

     SIZE(1,nslist) = nip
     SIZE(2,nslist) = njp
     SIZE(3,nslist) = nkp
     SIZE(4,nslist) = option
     SIZE(5,nslist) = nxpol
     SIZE(6,nslist) = pointer3
     SIZE(10,nslist) = pointer2
     SIZE(7,nslist) = norder 
     SIZE(8,nslist) = ipoint3 
     SIZE(9,nslist) = ipoint4
     WRITE(*,*) 'BLOCKING READ , ',nslist 
        
  END IF


  relax = relaxkind
  IF (relax.EQ.1.0) relax = rlxfac


  ! Storing the old solution
  IF (relax.LT.0) THEN
     CALL setold(nip,njp,nkp,nxpol,nutot,&
          norder,nlist,pointer3, pointer2,&
          xpoli(ipoint3),xpoli(ipoint3+2*nxpol),&
          xpolr(ipoint4),field, &
          oldsol(1,ipoint1),maxpts)
     RETURN
  END IF

  ! Inter-/ Extrapolation the new solution
  IF (relax.NE.0) THEN
     IF ((option.EQ.1).OR.(option.EQ.3)) THEN
        CALL xtrpol(nip,njp,nkp,nxpol,nutot,&
             norder,nlist,pointer3, pointer2,&
             xpoli(ipoint3),xpoli(ipoint3+2*nxpol),&
             xpolr(ipoint4),field, value)
     ELSE
        CALL xtrpolp(nip,njp,nkp,nxpol,nutot,&
             norder,nlist,pointer3, pointer2,&
             xpoli(ipoint3),xpoli(ipoint3+2*nxpol),&
             xpolr(ipoint4),field)
     END IF
  END IF

  ! Calculating the relaxed solution
  IF (relax.NE.1.0) THEN 
     CALL rlxold(nip,njp,nkp,nxpol,nutot,&
          norder,nlist,pointer3, pointer2,&
          xpoli(ipoint3),xpoli(ipoint3+2*nxpol),&
          xpolr(ipoint4),field, &
          oldsol(1,ipoint1),relax,maxpts)
  END IF

  ! deallocate(xpoli,xpolr,oldsol,slist,size,pointer)

END SUBROUTINE xpolquad





