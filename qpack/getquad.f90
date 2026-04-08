SUBROUTINE getquad(ni,nj,nk,nvirtual,blocked,gridnum,comp)
  ! ======================================================================
  !
  !        programmer     Frederic Tremblay
  !        version 1.0    date   02-05-1999
  !
  !   copyright (c) 1999-2002
  !   all rights reserved TU-Muenchen and F.Tremblay. No part of this publication
  !   may be reproduced,
  !   stored in a retrieval system ( e.g., in memory, disk, or core)
  !   or be transmitted by any means, electronic, mechanical, photocopy,
  !   recording, or otherwise, without written permission from the
  !   author.
  !
  !
  ! **********************************************************************

  IMPLICIT NONE
  INTEGER ni,nj,nk,nvirtual,gridnum,comp
  REAL blocked(ni*nj*nk)


  INTEGER, PARAMETER :: maxmax = 200
  INTEGER nlistmax ,narea
  INTEGER slist, i,j,k,nip,njp,nkp, memu, POINTER,dum1,&
       icount,icount2,icount3,nchar, size, ipoint, ipoint2, nslist, dummy, offset
  ALLOCATABLE :: slist(:,:),POINTER(:),SIZE(:,:) 
  LOGICAL firstcall, trashold
  CHARACTER(len=10) :: prefix
  CHARACTER(len=4) :: suffix
  CHARACTER(len=4) :: number
  CHARACTER (len=18) :: filout
  DATA firstcall /.TRUE./

  SAVE :: nslist, slist, firstcall, memu, POINTER, size, nlistmax
  !   trashold = .TRUE.
  IF (firstcall) THEN
     OPEN (unit=13,file="QPack.dat")
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
     READ(13,*)
     READ(13,*) nlistmax 
     CLOSE(13)

     ALLOCATE (POINTER(maxmax),slist(4,maxmax),SIZE(3,maxmax))

     !      firstcall = .FALSE.
     nslist = 0
     slist = 0

     POINTER = 0
     memu = 0
  END IF

  ! check if in memory
  IF (gridnum.GT.maxmax) THEN

     WRITE(*,*) 'error, gridnum is greater than maxmax', gridnum, maxmax
     STOP
  END IF
  ipoint = slist(comp,gridnum)
  IF (ipoint.NE.0) THEN
     blocked = 0.0
     ipoint2 = POINTER(ipoint)
     nip = SIZE(1,ipoint)
     njp = SIZE(2,ipoint)
     nkp = SIZE(3,ipoint)
  ELSE
     prefix = 'BLOCK_XTR_'
     suffix = '.DAT'
     WRITE(number(1:4),'(i4)') 10*gridnum+comp
     filout = prefix//number//suffix
     DO icount2 = 1,8
        IF (filout(10+1:10+1).EQ.' ') THEN
           DO icount3 = 1,8-icount2
              filout(10+icount3:10+icount3)=filout(10+icount3+1:10+icount3+1)
           END DO
        END IF
     END DO
     nchar = INT(LOG10(float(10*gridnum+comp))+1) + 4
     OPEN (unit=12,file=filout(1:10+nchar),form='UNFORMATTED')
     READ(12) dummy
     READ(12) nip,njp,nkp
     READ(12) dummy
     !      if (dummy.eq.0) dum1 = 999
     READ(12) dummy
     READ(12) dummy
     READ(12) dummy
     !      if (dum1.eq.999) goto 111
     READ(12) dummy
     READ(12) dummy
     READ(12) dummy
111  CONTINUE       
     READ(12) (blocked(i),i=1,nip*njp*nkp)
     CLOSE(12)
     nslist = nslist + 1
     IF (nslist.GT.nlistmax) THEN
        WRITE(*,*) 'nslist exceeded'
        STOP
     END IF
     slist(comp,gridnum) = nslist
     POINTER(nslist) = memu + 1
     memu = memu + nip*njp*nkp
     SIZE(1,nslist) = nip
     SIZE(2,nslist) = njp
     SIZE(3,nslist) = nkp
     ipoint2 = POINTER(nslist)
  END IF

  !   if (trashold) then
  !      firstcall = .TRUE.
  !   end if
  DEALLOCATE (POINTER, slist, size)
END SUBROUTINE getquad





