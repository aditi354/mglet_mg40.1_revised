subroutine stencil_spline3(klsm,nstenc,xin,xx,xp,nstenc_reduced,&
              vel,cst,alldir,option,part)

! natural spline with polynom third grade
! funktioniert nur mit nstenc=3 !!!
! Es werden drei Polynome = parts erzeugt
! xin(1) = xx

  implicit none
  integer nstenc,part,alldir,option,nstenc_reduced
  real klsm(4),xin(nstenc+1),&
       xx,xp,vel,cst(6)
  !local variables
  !integer 
  real h(3),f(4),bb(12),k(4),&
       kk(8),dd(12),g(12),dx

  !h(1) = x2 - x1
  !h(2) = x3 - x2
  !h(3) = x4 - x3
  write(33,*) "jep"
 
  h(1) = xin(2) - xin(1)
  h(2) = xin(3) - xin(2)
  h(3) = xin(4) - xin(3)

  f(1) = 2 * (h(1)+h(2))
  f(2) = h(2)
  f(3) = h(2)
  f(4) = 2 * (h(2)+h(3))

  k(1) = +f(4) / (f(1)*f(4)-f(2)*f(3))
  k(2) = -f(2) / (f(1)*f(4)-f(2)*f(3)) 
  k(3) = +f(3) / (f(2)*f(3)-f(1)*f(4))
  k(4) = -f(1) / (f(2)*f(3)-f(1)*f(4))

  kk(1) = 3*k(1)/h(1)
  kk(2) = 3*k(2)/h(2) - 3*k(1)/h(1) - 3*k(1)/h(2)
  kk(3) = 3*k(1)/h(2) - 3*k(2)/h(2) - 3*k(2)/h(3)
  kk(4) = 3*k(2)/h(3)

  kk(5) = 3*k(3)/h(1)
  kk(6) = 3*k(4)/h(2) - 3*k(3)/h(1) - 3*k(3)/h(2)
  kk(7) = 3*k(3)/h(2) - 3*k(4)/h(2) - 3*k(4)/h(3)
  kk(8) = 3*k(4)/h(3) 

  bb(1) = -1/h(1) - kk(1)*h(1)/3
  bb(2) =  1/h(1) - kk(2)*h(1)/3
  bb(3) =          -kk(3)*h(1)/3
  bb(4) =          -kk(4)*h(1)/3

  bb(5) =          -2*h(2)*kk(1)/3 - h(2)*kk(5)/3
  bb(6) = -1/h(2) - 2*h(2)*kk(2)/3 - h(2)*kk(6)/3
  bb(7) =  1/h(2) - 2*h(2)*kk(3)/3 - h(2)*kk(7)/3
  bb(8) =          -2*h(2)*kk(4)/3 - h(2)*kk(8)/3

  bb(9)  =         -2*h(3)*kk(5)/3
  bb(10) =         -2*h(3)*kk(6)/3
  bb(11) = -1/h(3) -2*h(3)*kk(7)/3
  bb(12) =  1/h(3) -2*h(3)*kk(8)/3
  !write(*,*) bb(11),bb(12),kk(7),kk(8),h(3)

  dd(1) = kk(1)/(3*h(1))
  dd(2) = kk(2)/(3*h(1))
  dd(3) = kk(3)/(3*h(1))
  dd(4) = kk(4)/(3*h(1))

  dd(5) = (kk(5)-kk(1)) / (3*h(2))
  dd(6) = (kk(6)-kk(2)) / (3*h(2))
  dd(7) = (kk(7)-kk(3)) / (3*h(2))
  dd(8) = (kk(8)-kk(4)) / (3*h(2))

  dd(9)  = -kk(5) / (3*h(3))
  dd(10) = -kk(6) / (3*h(3))
  dd(11) = -kk(7) / (3*h(3))
  dd(12) = -kk(8) / (3*h(3))

  dx = xp - xin(1)
  
  g(1) = 1 + bb(1)*dx + dd(1)*dx**3
  g(2) =     bb(2)*dx + dd(2)*dx**3
  g(3) =     bb(3)*dx + dd(3)*dx**3
  g(4) =     bb(4)*dx + dd(4)*dx**3

  dx = xp - xin(2)
 
  g(5) =     bb(5)*dx + kk(1)*dx**2 + dd(5)*dx**3 
  g(6) = 1 + bb(6)*dx + kk(2)*dx**2 + dd(6)*dx**3 
  g(7) =     bb(7)*dx + kk(3)*dx**2 + dd(7)*dx**3 
  g(8) =     bb(8)*dx + kk(4)*dx**2 + dd(8)*dx**3
 
  dx = xp - xin(3)
 
  g(9)  =     bb(9) *dx + kk(5)*dx**2 + dd(9) *dx**3 
  g(10) =     bb(10)*dx + kk(6)*dx**2 + dd(10)*dx**3
  g(11) = 1 + bb(11)*dx + kk(7)*dx**2 + dd(11)*dx**3
  g(12) =     bb(12)*dx + kk(8)*dx**2 + dd(12)*dx**3
  !write(*,*) bb(11),bb(12),kk(7),kk(8),dd(11),dd(12)


  if (part.eq.1) then 
     cst(alldir) = g(1) * vel 
     klsm(1) = g(2)
     klsm(2) = g(3)
     klsm(3) = g(4)
  elseif (part.eq.2) then
     cst(alldir) = g(5) * vel 
     klsm(1) = g(6)
     klsm(2) = g(7)
     klsm(3) = g(8)
  elseif (part.eq.3) then
     cst(alldir) = g(9) * vel 
     klsm(1) = g(10)
     klsm(2) = g(11)
     klsm(3) = g(12)
     !write(*,*) g(10),g(11),g(12)
  end if



end subroutine stencil_spline3


