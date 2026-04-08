SUBROUTINE  stencil_determine(ni,nj,nk,x,y,z,tmpabi,blocked,ntopol, topol,i,j,k,&
     nstenc,const,stenci,stencr,ndepth,xc,yc,zc, ntrimax,im,jm,km,&
     iproxi2,compon,maccur,tmpcst,tmpsteni,tmpstenr,distance,&
     direction,i1d,ni1d,im1d,x1d,xc1d,xp1d,jmpdir,position,&
     nstenc_reduced,error,whatipol,gridnum,bconds,&
     redord,option,nredord)

  ! Fred Tremblay
  ! Nick Peller
  !
  ! modifications:
  ! jk 11. 8.2005 : Einfuehren von reduced order
  !                 Verhindern, das stencils aus dem Gebiet ragen

  IMPLICIT NONE 
  INTEGER ntrimax, ntopol, im,jm,km, compon,&
       i,j,k, ndepth, &
       ni,nj,nk,nstenc,stenci(3,6,nstenc+1),&
       tmpabi(ni,nj,nk,ntrimax),iproxi2(6),&
       tmpsteni(nstenc+1,3,6),error(6,6),&
       blocked(ni*nj*nk),whatipol,nstenc_reduced(6),&
       gridnum,bconds(6),im1d,position,jmpdir,&
       i1d,ni1d,direction,redord,option,nredord
  REAL x(ni+1),y(nj+1),z(nk+1),const(6),&
       stencr(6,nstenc+1), xc(ni),yc(nj),&
       zc(nk),topol(4,3,ntopol),maccur,&
       distance(6),xb_pos,xb_neg,xp,yp,zp,&
       tmpcst(6),tmpstenr(nstenc+1,6),&
       x1d(ni1d+1),xc1d(ni1d),xp1d
  ! lokale Variablen
  INTEGER limit_pos,limit_neg,gdepth_pos,&
       gdepth_neg,sbegin_pos,sbegin_neg,&
       ipol_dir,ibound_pos,ibound_neg,&
       xb_found,ntriangles,itri,icount2,&
       it,jt,kt, alldir,ndepth_pos, ndepth_neg,&
       zaehler,zaehler2,it1d,sbegin,inc,dir1,dir2,&
       nstenc_reduced_pos, nstenc_reduced_neg,&
       stopit,casep,nocopy(6)
  REAL minu,value,x1,y1,z1,x2,y2,z2,&
       x3,y3,z3,a,b,c,d, velocity,xb_tmp,&
       eps,small,koeffizient(4),&
       xin1d(nstenc+1),vel_pos, vel_neg



  dir1 = 1+(direction-1)*2
  dir2 = 2+(direction-1)*2
  xb_found = 0
  limit_pos = 0
  limit_neg = 0
  ipol_dir = 0
  gdepth_pos = 0
  gdepth_neg = 0
  ndepth_pos = ndepth
  ndepth_neg = ndepth
  nstenc_reduced_pos = 0
  nstenc_reduced_neg = 0
  sbegin_pos = 0
  sbegin_neg = 0
  ibound_pos = 0
  ibound_neg = 0 
  vel_pos = 0.0
  vel_neg = 0.0
  xp = xc(i)
  yp = yc(j)
  zp = zc(k)
  small = maccur
  it = i
  jt = j
  kt = k
  stopit = 0
  nocopy = 0
  im1d = ni1d
  ! For pressure interpolation - change indices
  IF(option.eq.2) THEN
     casep = 1
  ELSE 
     casep = 0
  END IF


  ! *** check_domain_limits ***
  ! Hier werden die maximalen und minimalen Inkremente bestimmt, so dass
  ! der Algorithmus spaeter nicht ueber die domain Raender hinweglaeuft
  ! Ausserdem wird gleich hier geprueft, ob genuegend Zellen fuer einen
  ! stencil vorhanden sind. 
  ! limit_pos/neg
  limit_pos = im1d - i1d
  IF ((limit_pos).GT.(ndepth+nstenc+casep+1)) THEN
     limit_pos = (ndepth+nstenc+casep+1)
  END IF
  IF (limit_pos.EQ.0) error(dir1,1) = 1

  ! Noch einmal für die negative Richtung
  limit_neg = i1d - 1
  IF ((limit_neg).GT.(ndepth+nstenc+casep+1)) THEN
     limit_neg = (ndepth+nstenc+casep+1)
  END IF
  IF (limit_neg.EQ.0) error(dir2,1) = 1 

  ! *** interpolation_direction ***
  ! Hier wird gesucht in welcher Richtung der Körper sich befindet 
  ! ipol_dir: 1 = positive Richtung, 2 = negative Richtung, 3 = beide
  DO zaehler = -(limit_neg), (limit_pos)
     IF ((blocked(position+(jmpdir*zaehler)).EQ.1).AND.(zaehler.LT.0)) THEN
        ipol_dir = 2
     END IF
     IF ((ipol_dir.EQ.0).AND.(blocked(position+(jmpdir*zaehler)).EQ.1).AND.&
          (zaehler.GT.0)) THEN 
        ipol_dir = 1
     END IF
     IF ((ipol_dir.EQ.2).AND.(blocked(position+(jmpdir*zaehler)).EQ.1).AND.&
          (zaehler.GT.0)) THEN
        ipol_dir = 3
     END IF
  END DO

  ! *** check_geometry_depth ***
  ! Hier wird bestimmt wie viele Nachbarzellen ebenfalls 
  ! geblockt sind  ! WICHTIG: Es wird nur
  ! aufsummiert, wenn der unmittelbare Nachbar geblockt ist.
  ! Eine Zelle mit blocked > 0 stoppt die Summation.
  DO zaehler = 1, (limit_pos)
     IF ((blocked(position+(jmpdir*zaehler)).LT.1).AND.&
          (zaehler.EQ.(gdepth_pos+1))) THEN
        gdepth_pos = gdepth_pos + 1
     END IF
  END DO
  IF (gdepth_pos.LT.ndepth) ndepth_pos = gdepth_pos
  IF (gdepth_pos.EQ.0) error(dir1,2) = 1
  ! Noch einmal fuer die negative Richtung
  DO zaehler = 1, (limit_neg) 
     IF ((blocked(position+(jmpdir*(-zaehler))).LT.1).AND.&
          (zaehler.EQ.(gdepth_neg+1))) THEN
        gdepth_neg = gdepth_neg + 1
     END IF
  END DO
  IF (gdepth_neg.LT.ndepth) ndepth_neg = gdepth_neg
  IF (gdepth_neg.EQ.0) error(dir2,2) = 1
  ! nachher wird von ndepth_neg bis ndepth_pos nach RBen gesucht

  ! *** check_stencil_cells ***
  ! Hier wird berechnet wie viele Zellen fuer den stencil
  ! vorhanden sind und bei welcher Zelle er beginnt
  nstenc_reduced_pos = 0

  DO zaehler = gdepth_pos + 1, min(limit_pos,gdepth_pos+nstenc+casep)
     sbegin_pos = gdepth_pos + 1
     IF ((blocked(position+jmpdir*zaehler).NE.1)) THEN
        EXIT
     ELSE  
        nstenc_reduced_pos = nstenc_reduced_pos + 1
     END IF
  END DO

  IF (nstenc_reduced_pos.GT.nstenc+casep) nstenc_reduced_pos = nstenc+casep
  IF (nstenc_reduced_pos.LT.nstenc+casep-redord) error(dir1,3) = 1
  IF (sbegin_pos.EQ.0) error(dir1,4) = 1
  ! Noch einmal fuer die negative Richtung
  nstenc_reduced_neg = 0

  DO zaehler = gdepth_neg + 1, min(limit_neg,gdepth_neg+nstenc+casep)
     sbegin_neg = gdepth_neg + 1
     IF ((blocked(position+jmpdir*(-zaehler)).NE.1)) THEN
        EXIT
     ELSE 
        nstenc_reduced_neg = nstenc_reduced_neg + 1
     END IF
  END DO

  IF (nstenc_reduced_neg.GT.nstenc+casep) nstenc_reduced_neg = nstenc+casep 
  IF (nstenc_reduced_neg.LT.nstenc+casep-redord) error(dir2,3) = 1
  IF (sbegin_neg.EQ.0) error(dir2,4) = 1


  ! If these errors have already occured
  ! most of the rest can be skipped
  IF ((error(dir1,4).EQ.1).AND.&
      (error(dir2,4).EQ.1)) THEN
!      (nstenc.NE.0)) THEN
     stopit = 1
  END IF



  ! Only neccessary for velocity interpolation
  IF (option.NE.2) THEN
     ! *** find_geom_boundary ***
     ! Durchsuchen der Zellen von -gdepth_neg bis gdepth_pos nach
     ! einer passenden Geometrie Randbedingung. Es wird die am weitesten
     ! links und rechts liegende Geometrierandbedingung gespeichert. Gibt
     ! es nur eine RB sind xb_neg und xb_pos identisch.
     xb_pos = -9e6
     xb_neg = 9e6  
     zaehler = -ndepth_neg  
     minu = 1e6
     ! Nur Berechnen, wenn zuvor kein Fehler
     ! aufgetreten ist. 
     IF (stopit.NE.1) THEN
        DO WHILE (zaehler.LT.(ndepth_pos+1))
           IF (direction.EQ.1)  it = i + zaehler
           IF (direction.EQ.2)  jt = j + zaehler  
           IF (direction.EQ.3)  kt = k + zaehler
           ntriangles = -blocked(position+jmpdir*zaehler)
           IF (zaehler.EQ.0) minu = 1e6 
           DO icount2 = 1, ntriangles
              xb_tmp = 0    
              itri = tmpabi(it,jt,kt,icount2)
              value = topol(4,compon,itri)
              x1 = topol(1,1,itri)
              x2 = topol(2,1,itri)
              x3 = topol(3,1,itri)
              y1 = topol(1,2,itri)
              y2 = topol(2,2,itri)
              y3 = topol(3,2,itri)
              z1 = topol(1,3,itri)
              z2 = topol(2,3,itri)
              z3 = topol(3,3,itri)
              a = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
              b = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
              c = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
              d = x1*a + y1*b + z1*c
              IF ((direction.EQ.1).AND.(a.NE.0)) THEN 
                 xb_tmp = (-c*(zp-z1)-b*(yp-y1))/a + x1
              END IF
              IF ((direction.EQ.2).AND.(b.NE.0)) THEN
                 xb_tmp = (-c*(zp-z1)-a*(xp-x1))/b + y1
              END IF
              IF ((direction.EQ.3).AND.(c.NE.0)) THEN 
                 xb_tmp = (-b*(yp-y1)-a*(xp-x1))/c + z1
              END IF
              it1d = i1d + zaehler
              eps = small*ABS(x1d(it1d+1)-x1d(it1d))            
              IF ((xb_tmp.LT.(x1d(it1d)-eps)).OR.&
                   (xb_tmp.GT.(x1d(it1d+1)+eps))) THEN 
                 xb_tmp = 0  
              ELSE  
                 IF (ABS(xb_tmp-xp1d).LT.minu) THEN
                    minu = ABS(xb_tmp-xp1d)     
                    IF ((zaehler).LE.0) THEN
                       ibound_neg = ABS(it1d - i1d)
                       xb_neg = xb_tmp
                       vel_neg = value
                    ELSE IF(zaehler.NE.0) THEN
                       ibound_pos = ABS(it1d - i1d)
                       xb_pos = xb_tmp  
                       vel_pos = value
                    END IF
                    xb_found = 1                
                 END IF
              END IF
           END DO
           zaehler = zaehler + 1   
        END DO
        IF (xb_neg.GT.1e6)  THEN
           xb_neg = xb_pos
           vel_neg = vel_pos
        ELSE IF (xb_pos.LT.-1e6) THEN
           xb_pos = xb_neg
           vel_pos = vel_neg
        END IF

     END IF

     ! Abstaende und xb-Fehlermoeglichkeit
     ! Wenn keine RB gefunden wurde muss der
     ! Abstand Null gesetzt werden
     IF (xb_found.EQ.0) THEN
        distance(dir1) = 0.0
        error(dir1,5) = 1
     ELSE
        distance(dir1) =  ABS(xb_pos-xp1d)
     END IF
     IF (xb_found.EQ.0) THEN
        distance(dir2) = 0.0
        error(dir2,5) = 1
     ELSE 
        distance(dir2) =  ABS(xb_neg-xp1d)
     END IF

  ELSE
     ! For pressure int.; distances for wheighting
     xb_pos = xc1d(i1d+sbegin_pos)
     xb_neg = xc1d(i1d-sbegin_neg)
     distance(dir1) =  ABS(xb_pos-xp1d)
     distance(dir2) =  ABS(xb_neg-xp1d)
  END IF



  ! Eine Zelle Zwischenraum maximal zw. Randbed.Zelle und 
  ! Anfang stencil erste Zelle
  error(dir1,6) = 1
  IF ((sbegin_pos - ibound_pos).LE.2) error(dir1,6) = 0
  IF ((ibound_neg.LE.2).AND.(sbegin_pos.LE.2)) error(dir1,6) = 0
  error(dir2,6) = 1
  IF ((sbegin_neg - ibound_neg).LE.2) error(dir2,6) = 0
  IF ((ibound_pos.LE.2).AND.(sbegin_neg.LE.2)) error(dir2,6) = 0 



  ! Ist zuvor ein Fehler aufgetreten?
  ! Dann kann kein stencil erstellt werden
  IF ((error(dir1,1).EQ.1).OR.(error(dir1,3).EQ.1).OR.&
       (error(dir1,4).EQ.1).OR.(error(dir1,5).EQ.1).OR.&
       (error(dir1,6).EQ.1)) THEN
     nstenc_reduced_pos = 0
     nocopy(dir1) = 1
  END IF
  IF ((error(dir2,1).EQ.1).OR.(error(dir2,3).EQ.1).OR.&
       (error(dir2,4).EQ.1).OR.(error(dir2,5).EQ.1).OR.&
       (error(dir2,6).EQ.1)) THEN
     nstenc_reduced_neg = 0 
     nocopy(dir2) = 1
  END IF




  ! *** calculate_stencils ***
  ! Jetzt endlich die eigentliche Berechnung 
  ! der stencils in beide Richtungen...
  ! write(76,*) '   subroutine calculate_stencils'
  !IF (stopit.NE.1) THEN
  it1d = 0
  xp1d = xc1d(i1d)
  DO zaehler = 1 ,2
     ! Die 2*3 Richtungen zusammenfassen
     alldir=zaehler+(direction-1)*2
     IF (zaehler.EQ.1) THEN
        sbegin = sbegin_pos
        inc = 1
        xin1d(1) = xb_pos
        velocity = vel_pos
        iproxi2(alldir) = 1
        nstenc_reduced(alldir) = nstenc_reduced_pos 
     ELSE
        sbegin = -sbegin_neg
        inc = -1
        xin1d(1) = xb_neg
        velocity = vel_neg
        iproxi2(alldir) = -1 
        nstenc_reduced(alldir) = nstenc_reduced_neg
     END IF

     ! Schreibt Positionen der stencils in anderes array um.
     ! Dadurch wird die Richtungs-Abhaengigkeit beseitigt!!
     ! -->xin1d(nstenc+1): Positionen der stencils
     it1d = i1d + sbegin
     DO zaehler2 = (2-casep), (nstenc+1)
        IF ((zaehler2-1+casep).LE.nstenc_reduced(alldir)) THEN
           xin1d(zaehler2) = xc1d(it1d) 
           it1d = it1d + inc         
        ELSE
           xin1d(zaehler2) = 0.0
        END IF
     END DO

     ! Wenn ein stencil vorhanden, dann berechnen
     IF (nstenc_reduced(alldir).NE.0) THEN
        IF (option.EQ.2) THEN ! For pressure interpolation
           CALL stencil_lagrange(koeffizient,nstenc,&
                xin1d,xin1d(1),xp1d,nstenc_reduced(alldir),&
                velocity,tmpcst,alldir,option)
        ELSEIF (option.EQ.3) THEN ! Stencil for adiabatic wall - Scalar
           CALL stencil_least_square_adiabat(koeffizient,nstenc,&
                xin1d,xin1d(1),xp1d,nstenc_reduced(alldir),&
                velocity,tmpcst,alldir,option,maccur)
        ELSEIF (whatipol.EQ.0) THEN
           IF (redord.eq.0) THEN
              CALL stencil_lagrange(koeffizient,nstenc,&
                   xin1d,xin1d(1),xp1d,nstenc_reduced(alldir),&
                   velocity,tmpcst,alldir,option)
           ELSE
              IF (nstenc_reduced(alldir).eq.nstenc+casep) THEN
                 CALL stencil_lagrange(koeffizient,nstenc,&
                      xin1d,xin1d(1),xp1d,nstenc_reduced(alldir),&
                      velocity,tmpcst,alldir,option)
              ELSE IF (nstenc_reduced(alldir).lt.nstenc+casep) THEN
                 ! Lineare Interpolation
                 CALL stencil_lagrange(koeffizient,nstenc,&
                      xin1d,xin1d(1),xp1d,nstenc_reduced(alldir)-redord-1,&
                      velocity,tmpcst,alldir,option)
              END IF
           END IF
        ELSEIF (whatipol.EQ.1) THEN
           IF (redord.eq.0) THEN
              CALL stencil_least_square(koeffizient,nstenc,&
                   xin1d,xin1d(1),xp1d,nstenc_reduced(alldir),&
                   velocity,tmpcst,alldir,option,maccur)
           ELSE
              IF (nstenc_reduced(alldir).eq.nstenc+casep) THEN
                 CALL stencil_least_square(koeffizient,nstenc,&
                      xin1d,xin1d(1),xp1d,nstenc_reduced(alldir),&
                      velocity,tmpcst,alldir,option,maccur)
              ELSE IF (nstenc_reduced(alldir).lt.nstenc+casep) THEN
                 ! Lineare Interpolation
                 CALL stencil_lagrange(koeffizient,nstenc,&
                      xin1d,xin1d(1),xp1d,nstenc_reduced(alldir)-redord,&
                      velocity,tmpcst,alldir,option)
!                 write(*,*) 'stencil_determine, redord:',i,j,k
                 nredord = nredord + 1
!!$                 write(*,*) koeffizient
!!$                 write(*,*) tmpcst(alldir),compon
!!$                 write(*,*) nstenc_reduced(alldir)
!!$                 write(*,*) xin1d, xp1d
!!$                 stop
              END IF
           END IF
        ELSEIF (whatipol.EQ.2) THEN
           CALL stencil_spline3(koeffizient,nstenc,&
                xin1d,xin1d(1),xp1d,nstenc_reduced(alldir),&
                velocity,tmpcst,alldir,option,1)
        ELSE
           WRITE(*,*) '***************************************'
           WRITE(*,*) 'ERROR IN QPACK: COULD NOT RECOGNIZE THE'
           WRITE(*,*) '                CHOICE OF INTERPOLATION' 
           WRITE(*,*) '                IN THE QPack.dat file' 
           WRITE(*,*) '***************************************'
           stop
        END IF

        ! erster Punkt des Stencils
        it1d = i1d + sbegin    

        ! Verschieben, falls stencil (nstenc+casep ) ueber Gebiet hinausragt
        DO
           IF ( ((zaehler.eq.1).and.(it1d+(nstenc+casep-1).gt.im1d))  &
                .OR. ((zaehler.eq.2).and.(it1d-(nstenc+casep-1).lt.1)) ) THEN
              
              IF ((abs(koeffizient(nstenc+casep)).gt.2*epsilon(1.0)) &
                   .and. (option.ne.2)) THEN
                 WRITE(*,*) 'stencil_determine, Fehler: Stencil ragt hinaus', &
                      koeffizient,i,j,k,nstenc_reduced(alldir)
                 STOP
              END IF
              it1d = it1d - inc
              DO zaehler2=(nstenc+casep),2,-1
                 koeffizient(zaehler2) = koeffizient(zaehler2-1)
              END DO
              koeffizient(1) = 0.0
           ELSE
              exit
           END IF
        END DO

        DO zaehler2=1, (nstenc+casep) 
           tmpsteni(zaehler2,1,alldir) = i
           tmpsteni(zaehler2,2,alldir) = j
           tmpsteni(zaehler2,3,alldir) = k 
           IF (zaehler2.LE.nstenc_reduced(alldir)) THEN
              ! Hier wird einfach die "direction" überschrieben
              tmpsteni(zaehler2,direction,alldir) = it1d
              tmpstenr(zaehler2,alldir) = koeffizient(zaehler2) 
           ELSE
              tmpstenr(zaehler2,alldir) = 0.0 
           END IF
           it1d = it1d + inc
        END DO



     ELSE
        IF((nstenc.eq.0).AND.(option.ne.2)) then 
           ! Wenn Treppchengitter
           tmpcst(alldir) = velocity  
           IF(nocopy(alldir).NE.1) THEN 
              nstenc_reduced(alldir) = 1
           ENDIF
        ELSE
           tmpcst(alldir) = 0.0
        END IF
        DO zaehler2=1, (nstenc+casep)
           tmpsteni(zaehler2,1,alldir) = 1
           tmpsteni(zaehler2,2,alldir) = 1
           tmpsteni(zaehler2,3,alldir) = 1
           tmpstenr(zaehler2,alldir) = 0.0     
        END DO
     END IF

!     IF (nstenc_reduced(alldir).EQ.2) THEN
!        write(99,*) "nstenc_red.:",i,j,k
!     END IF



  END DO





END SUBROUTINE stencil_determine
