subroutine stencil_lagrange(koeffizient,nstenc,xin,xb,xp,&
     nstenc_reduced,vel,tmpcst,alldir,option)

  implicit none
  integer nstenc,alldir,option,nstenc_reduced
  real koeffizient(4),xin(nstenc+1),xb,xp,vel,tmpcst(6)
  ! lokale Variablen
  integer counter,counter2,casep 
  real sum

  koeffizient = 0.0
  ! Unterschied fuer Druck- bzw. UVW Interpol.
  if (option.eq.2) then
     casep = 1
     ! Fuer Druck Center Stencil Null
     tmpcst(alldir) = 0.0

  else
     casep = 0
     ! Der additive Center Stencil
     sum = 1
     do counter = 2,nstenc_reduced+1
        sum = sum*((xp-xin(counter))/&
             (xin(1)-xin(counter)))
     end do
     tmpcst(alldir) = sum * vel

  end if

  ! Die multiplikativen Stencils
  do counter=(2-casep),nstenc_reduced+1-casep
     sum = 1   
     do counter2=1,nstenc_reduced+1-casep
        if (counter2.ne.counter) then
           sum = sum*((xp-xin(counter2))/&
                (xin(counter)-xin(counter2)))
        end if
     end do
     koeffizient(counter-1+casep) = sum
  end do


! ! Die multiplikativen Stencils
!  do counter=2,nstenc_reduced+1
!     sum = 1
!     do counter2=(1+casep),nstenc_reduced+1
!        if (counter2.ne.counter) then
!           sum = sum*((xp-xin(counter2))/&
!                (xin(counter)-xin(counter2)))
!        end if
!     end do
!     koeffizient(counter-1) = sum
!  end do


end subroutine stencil_lagrange
