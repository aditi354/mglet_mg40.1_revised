subroutine stencil_least_square(KLSM,nstenc,&
              xin,xx,xp,nstenc_reduced,VEL,tmpcst,&
              alldir,option,maccur)
 
  implicit none
  integer ncount,nstenc,nstenc_reduced,&
          alldir,option
  real xin(nstenc+1), xx, xp, VEL, maccur,tmpcst(6)

! *********************************************************************
! N.Peller : Variablendefinitionen fuer least square stencils
! VEL   : Geschwindigkeit des Starrkoerpers
! SA    : Summe A (1 - 4)
! SV    : Summe V (1,2)
! DXMR  : Differenz von Xm - Xp  (m = Laufindex für Pkt 1 - 4)
! DQXMR : Differenz von Xm^2 - Xp^2
! DXPR  : Differenz von Xp - Xr  (r = Randpunkt)
! DQXPR : Differenz von Xp^2 - Xr^2
! KLSM  : Koeffizient Least Square Method,der berechnet wird; K (1 - 4)
! CONST : Konstante auf Grund Drehgeschwindigkeit der Geometrie
  real SA(4), DXMR(4), DQXMR(4), DXPR, DQXPR, KLSM(4),SV(2), nenner
  real shift,zaehler
! *********************************************************************

  ! Initialisierung der Variablen und arrays
  SV = 0.0
  SA = 0.0
  DXMR = 0.0
  DQXMR = 0.0
  DXPR = 0.0
  DQXPR = 0.0
  KLSM = 0.0
  shift = 1e6

  !Berechnung der Summenterme
  do ncount=2,nstenc_reduced+1
     SA(1) = SA(1) + (xin(ncount) - xx)**2 
     SA(2) = SA(2) + (xin(ncount)**2 - xx**2) * (xin(ncount) - xx)
     SA(3) = SA(2)
     SA(4) = SA(4) + (xin(ncount)**2 - xx**2)**2 
     
     SV(1) = SV(1) + (xin(ncount)) - xx
     SV(2) = SV(2) + (xin(ncount))**2 - xx**2
 
     DXMR(ncount-1)  = xin(ncount) - xx
     DQXMR(ncount-1) = xin(ncount)**2 - xx**2
  end do


  !Berechnung der Differenzen zw. Xp und dem Randwert xx=xr
  DXPR  = (xp - xx)
  DQXPR = (xp**2 - xx**2)
  

  !Schliesslich Einsetzen der Summen und Differenzen in die Formel
  !für die least square Koeffizienten
  !Nenner ist immer gleich
  SA(1) = SA(1) * shift 
  SA(2) = SA(2) * shift 
  SA(3) = SA(3) * shift 
  SA(4) = SA(4) * shift  
  nenner = SA(2)*SA(3)-SA(1)*SA(4)
  nenner = nenner / shift
  if (abs(nenner).gt.maccur) then  
     !Least Square stencils ohne Geschwindigkeit
     do ncount=2,nstenc_reduced+1
        zaehler =  (SA(2) * DQXMR(ncount-1) - SA(4) * DXMR(ncount-1) ) * DXPR +&
                   (SA(3) * DXMR(ncount-1) - SA(1) * DQXMR(ncount-1) ) * DQXPR
        KLSM(ncount-1) = zaehler /  nenner
     end do
     !Center stencil aus Geschw. Ransbedingung
     zaehler =  ( -SA(2) * SV(2) + SA(4) * SV(1) ) * DXPR +&
                ( -SA(3) * SV(1) + SA(1) * SV(2) ) * DQXPR
     tmpcst(alldir) = VEL * ( 1 + ( zaehler / nenner ))
  else
     write(76,*) '      Nenner Null!!'  
  end if
 

end subroutine 
