!######## DEFINE STATEMENTS ###########


!######################################

SUBROUTINE solint(IDIM3D,IDIM2D,IDIM1D,IDIMA,&
     X,Y,Z,DX,DY,DZ,DDX,DDY,DDZ,&

     U,V,W,P,TIME,GMOL)




  IMPLICIT NONE

  INTEGER idim3d,idim2d,idim1d,idima,kk,jj,ii,&
          i,ilevel,igrid,ibb,ibu,ib3,ip1,ip2,ip3



  REAL u (idim3d),v (idim3d),w (idim3d),p (idim3d),&
       uav(idim2d),vav(idim2d),wav(idim2d),pav(idim2d)

  REAL x (idim1d),y (idim1d),z (idim1d),&
       dx (idim1d),dy (idim1d),dz (idim1d),&
       ddx(idim1d),ddy(idim1d),ddz(idim1d)
  REAL time,gmol

  CHARACTER (LEN=10) action,forma,homog
  CHARACTER (LEN=3) doit

  LOGICAL,SAVE:: first = .TRUE.



      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS =128 )
      PARAMETER ( MAXBOCONDS = 10 )



      INTEGER MAXLVLMAX,MINLVLMIN,MAXGRDSOFLVL

      PARAMETER ( MAXLVLMAX = 4, MINLVLMIN = -4 )
      PARAMETER ( MAXGRDSOFLVL =260)

      COMMON /COLEVEL/ MAXLEVEL,MINLEVEL
      COMMON /COLEVEL/ NOFLEVEL,IGRDOFLEVEL,NOFPSDIR,IGRDOFPSDIR
      COMMON /COLEVEL/ NOFVPIT ,IGRDOFVPIT ,NOFTST, IGRDOFTST
      COMMON /COLEVEL/ NOFSLCED,IGRDOFSLCED,NOFSLICE,IGRDOFSLICE
      COMMON /COLEVEL/ NOFSCAI, IGRDOFSCAI,NGRIDPERLEVEL,NOFTHISGRID

      INTEGER MAXLEVEL,MINLEVEL

      INTEGER NOFLEVEL ( MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFLEVEL(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFVPIT ( MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFVPIT (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFTST ( MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFTST (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSCAI ( MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSCAI (MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSLCED ( MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSLCED(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFSLICE ( MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFSLICE(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NOFPSDIR ( MINLVLMIN:MAXLVLMAX)
      INTEGER IGRDOFPSDIR(MAXGRDSOFLVL,MINLVLMIN:MAXLVLMAX)
      INTEGER NGRIDPERLEVEL
      INTEGER NOFTHISGRID(MAXGRDSOFLVL)

  INTEGER,PARAMETER :: IDIMP = 5000

  COMMON /COSOLINT/ NLOCP,LOCPNTS,DXLOC,&
       SOLOUT,DIST,STAG,SHIFT

  INTEGER NLOCP(MAXGRIDS),SHIFT

  LOGICAL STAG

  REAL LOCPNTS(3,IDIMP,MAXGRIDS),SOLOUT(4,4,IDIMP,MAXGRIDS),&
       DXLOC(3,IDIMP,MAXGRIDS),DIST


      COMMON /COGRDPRO/ LEVEL,LCHILD,XMIN,YMIN,ZMIN,XTOT,YTOT,ZTOT
      COMMON /COGRDPRO/ XHOMOG,YHOMOG,ZHOMOG,NXGRAE,NYGRAE,NZGRAE
      COMMON /COGRDPRO/ GRADPX,UBULKX,LTST,LVP,LSCAI,LPLEVEL,LPOISSONDIR
      COMMON /COGRDPRO/ LSLICE,NXSLICE,NYSLICE,NZSLICE,NVPGRIDS
      COMMON /COGRDPRO/ CONV1SANF,CONV1SEND,TRANSLES1,TRANSLES2
      COMMON /COGRDPRO/ GRADPXOLD

      INTEGER LEVEL(MAXGRIDS),NVPGRIDS(MAXGRIDS)
      INTEGER NXGRAE(MAXGRIDS),NYGRAE(MAXGRIDS),NZGRAE(MAXGRIDS)
      INTEGER NXSLICE(MAXGRIDS),NYSLICE(MAXGRIDS),NZSLICE(MAXGRIDS)

      REAL XTOT(MAXGRIDS),YTOT(MAXGRIDS),ZTOT(MAXGRIDS)
      REAL XMIN(MAXGRIDS),YMIN(MAXGRIDS),ZMIN(MAXGRIDS)
      REAL GRADPX(MAXGRIDS),UBULKX(MAXGRIDS)
      REAL CONV1SANF(MAXGRIDS),CONV1SEND(MAXGRIDS)
      REAL TRANSLES1(MAXGRIDS),TRANSLES2(MAXGRIDS),GRADPXOLD(MAXGRIDS)

      LOGICAL LCHILD(MAXGRIDS),LSLICE(MAXGRIDS),LPOISSONDIR(MAXGRIDS)
      LOGICAL XHOMOG(MAXGRIDS),YHOMOG(MAXGRIDS),ZHOMOG(MAXGRIDS)
      LOGICAL LTST(MAXGRIDS),LVP(MAXGRIDS),LSCAI(MAXGRIDS)
      LOGICAL LPLEVEL(MAXGRIDS)
! Fuer Postprocessing offline ohne mglet







  WRITE(*,*) 'Aufruf solint.f90!!'

  OPEN (unit = 100,file='solint.dat')

  READ(100,*) doit
  READ(100,*) action, forma
  READ(100,*) homog
  READ(100,*) dist
  CLOSE(100)

  ! Exit the routine when no values
  ! should be calculated

  IF (doit(1:2).EQ.'no') THEN
     RETURN
  END IF



  DO ILEVEL = MINLEVEL,MAXLEVEL
     DO I = 1,NOFTST(ILEVEL)
        IGRID = IGRDOFTST(I,ILEVEL)




! NUR, FALLS GITTER NICHT GEBIETSZERLEGT

           IF ( .NOT. LSLICE(IGRID) ) THEN

              CALL MGDIMS (KK,JJ,II,IGRID)
              CALL MGPOINT (IP3,IP2,IP1,IBB,IB3,IBU,IGRID)

              CALL homoy(kk,jj,ii,x(ip1),z(ip1),&
                   u(ip3),v(ip3),w(ip3),p(ip3),&
                   uav(ip2),vav(ip2),wav(ip2),pav(ip2),forma,igrid)

              CALL interpxz(kk,ii,x(ip1),z(ip1),dx(ip1),dz(ip1),&
                   uav(ip2),vav(ip2),wav(ip2),pav(ip2),time,forma,igrid)

              CALL cfcdxz(gmol,time,igrid)

           ENDIF





     ENDDO
  ENDDO
END SUBROUTINE solint


  SUBROUTINE interpxz(kk,ii,x,z,dx,dz,uav,vav,wav,pav,time,forma,igrid)

  IMPLICIT NONE



      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS =128 )
      PARAMETER ( MAXBOCONDS = 10 )


  INTEGER,PARAMETER :: IDIMP = 5000

  COMMON /COSOLINT/ NLOCP,LOCPNTS,DXLOC,&
       SOLOUT,DIST,STAG,SHIFT

  INTEGER NLOCP(MAXGRIDS),SHIFT

  LOGICAL STAG

  REAL LOCPNTS(3,IDIMP,MAXGRIDS),SOLOUT(4,4,IDIMP,MAXGRIDS),&
       DXLOC(3,IDIMP,MAXGRIDS),DIST

  INTEGER kk,ii,npoint,i,k,ipoint,nbody,ibody,igrid,ilocp

  REAL uav(kk,ii),vav(kk,ii),wav(kk,ii),pav(kk,ii),&
       phi(idimp),dphidx(idimp),dphidz(idimp),&
       x(ii),xp(ii),dx(ii),z(kk),zp(kk),dz(kk)

! REAL locx,locz,phi,dphidx,dphidz,time,minx,minz,maxx,maxz
  REAL locx,locz,time,minx,minz,maxx,maxz


  REAL points(3,idimp),norm

  CHARACTER (LEN=10) forma,homog

  LOGICAL,SAVE :: firstcall

  DATA firstcall /.TRUE./

  IF (FIRSTCALL) THEN
     write(*,*) ii,kk
     minx = x(2) + 0.5 * dx(2)
     minz = z(2) + 0.5 * dz(2)
     maxx = x(ii-2) + 0.5 * dx(ii-2)
     maxz = z(kk-2) + 0.5 * dz(kk-2)

     ilocp = 0
     solout = 0.0

     OPEN (unit=85,file='solint.pnts')
     READ(85,*) npoint
     IF (npoint .GT. idimp) THEN
        WRITE(6,*)'idimp smaller than npoint in interpxz'
! CALL ERRR ( 505,' SOLINT ' )
     ENDIF

     DO i = 1, npoint
        READ(85,*) points(1,i),points(2,i),points(3,i)
        IF (points(1,i) .LT. minx .OR. points(1,i) .GT. maxx .OR. &
             points(3,i) .LT. minz .OR. points(3,i) .GT. maxz) CYCLE
        ilocp = ilocp + 1
        locpnts(:,ilocp,igrid) = points(:,i)
     END DO

     nlocp(igrid) = 2 * ilocp ! Aenderung
     shift = ilocp
     write(*,*) nlocp(igrid),igrid,minx,minz,maxx,maxz

     DO i = 1, npoint-1
        dxloc(1,i,igrid) = 0.5*(locpnts(1,i+1,igrid) - locpnts(1,i,igrid))
        dxloc(2,i,igrid) = 0.5*(locpnts(2,i+1,igrid) - locpnts(2,i,igrid))
        dxloc(3,i,igrid) = 0.5*(locpnts(3,i+1,igrid) - locpnts(3,i,igrid))
        !write(21,*) dxloc(1,i,igrid),dxloc(2,i,igrid),dxloc(3,i,igrid)
     END DO

     OPEN(UNIT=10,FILE='solint.oldpnts',FORM='FORMATTED')
     DO i = 1, npoint-1
        ! Writing old location of points
        write(10,*) locpnts(1,i,igrid),&
                    locpnts(2,i,igrid),&
                    locpnts(3,i,igrid)
     END DO
     CLOSE(10)


     DO i = 1, npoint-1
        locpnts(1,i,igrid) = locpnts(1,i,igrid) + dxloc(1,i,igrid)
        !locpnts(2,i,igrid) = ...bleibt gleich
        locpnts(3,i,igrid) = locpnts(3,i,igrid) + dxloc(3,i,igrid)
     END DO


     DO i = 1, npoint-1
        norm = sqrt((dxloc(1,i,igrid))**2+(dxloc(3,i,igrid))**2)
        locpnts(1,i+shift,igrid) = (locpnts(1,i,igrid)) - &
                                   (dxloc(3,i,igrid)/norm)*dist
        locpnts(2,i+shift,igrid) = locpnts(2,i,igrid)
        locpnts(3,i+shift,igrid) = (locpnts(3,i,igrid)) + &
                                   (dxloc(1,i,igrid)/norm)*dist
    END DO


    OPEN(UNIT=10,FILE='solint.newpnts',FORM='FORMATTED')
    OPEN(UNIT=11,FILE='solint.perpnts',FORM='FORMATTED')
    OPEN(UNIT=12,FILE='solint.normals',FORM='FORMATTED')
    DO i = 1, npoint-1
        write(10,*) locpnts(1,i,igrid),&
                    locpnts(2,i,igrid),&
                    locpnts(3,i,igrid)
        write(11,*) locpnts(1,i+shift,igrid),&
                    locpnts(2,i+shift,igrid),&
                    locpnts(3,i+shift,igrid)
        write(12,*) locpnts(1,i,igrid),&
                    locpnts(2,i,igrid),&
                    locpnts(3,i,igrid)
        write(12,*) locpnts(1,i+shift,igrid),&
                    locpnts(2,i+shift,igrid),&
                    locpnts(3,i+shift,igrid)
        write(12,*)
        write(12,*)
     END DO
     CLOSE(10)
     CLOSE(11)
     CLOSE(12)


  END IF

! first the U1 grid

  DO i = 1, ii
     xp(i) = x(i) + 0.5 * dx(i)
  END DO

  CALL findlocxz(kk,ii,xp,z,uav,phi,dphidx,dphidz,igrid)

  do i = 1,nlocp(igrid)
     solout(1,1,i,igrid) = phi(i)
     solout(1,2,i,igrid) = dphidx(i)
     solout(1,4,i,igrid) = dphidz(i)
  enddo
  !now the U2 Grid

  CALL findlocxz(kk,ii,x,z,vav,phi,dphidx,dphidz,igrid)

  do i = 1,nlocp(igrid)
     solout(2,1,i,igrid) = phi(i)
     solout(2,2,i,igrid) = dphidx(i)
     solout(2,4,i,igrid) = dphidz(i)
  enddo

  !now the U3 Grid

  DO k = 1, kk
     zp(k) = z(k)+ 0.5 * dz(k)
  END DO

  CALL findlocxz(kk,ii,x,zp,wav,phi,dphidx,dphidz,igrid)

  do i = 1,nlocp(igrid)
     solout(3,1,i,igrid) = phi(i)
     solout(3,2,i,igrid) = dphidx(i)
     solout(3,4,i,igrid) = dphidz(i)
  enddo

  !now the P Grid

  CALL findlocxz(kk,ii,x,z,pav,phi,dphidx,dphidz,igrid)

  do i = 1,nlocp(igrid)
     solout(4,1,i,igrid) = phi(i)
     solout(4,2,i,igrid) = dphidx(i)
     solout(4,4,i,igrid) = dphidz(i)
  enddo

!!$ IF(forma(1:5).EQ.'ascii') THEN

  if (igrid .eq. 2) then
     OPEN (unit=86, file='ft.86.2',form = 'formatted')
  else
     OPEN (unit=86, file='ft.86.3',form = 'formatted')
  endif


! WRITE(86,*) time
     DO i = 1,nlocp(igrid)
        WRITE(86,*) i,locpnts(1,i,igrid),&
             solout(4,1,i,igrid),solout(4,2,i,igrid),solout(4,4,i,igrid)
     ENDDO

!!$ ELSE
!!$ OPEN (unit=86,file='ft.86', form = 'unformatted' )
!!$ WRITE(86) time
!!$ WRITE(86) (solout(:,:,i),i=1,nlocp)
!!$ END IF
!!$
  CLOSE (85)
  CLOSE (86)

END SUBROUTINE interpxz

  SUBROUTINE findlocxz(kk,ii,x,z,field,phi,dphidx,dphidz,igrid)

  IMPLICIT NONE



      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS =128 )
      PARAMETER ( MAXBOCONDS = 10 )


  INTEGER,PARAMETER :: IDIMP = 5000

  COMMON /COSOLINT/ NLOCP,LOCPNTS,DXLOC,&
       SOLOUT,DIST,STAG,SHIFT

  INTEGER NLOCP(MAXGRIDS),SHIFT

  LOGICAL STAG

  REAL LOCPNTS(3,IDIMP,MAXGRIDS),SOLOUT(4,4,IDIMP,MAXGRIDS),&
       DXLOC(3,IDIMP,MAXGRIDS),DIST

  INTEGER ip,kp,i,k,ii,kk,ipoint,igrid

  REAL x(ii),z(kk),field(kk,ii)

  REAL facx,facz,ufield(4),psi(4),dfacx,dfacz,&
       phi(idimp),dphidx(idimp),dphidz(idimp)

  phi = 0.0 ; dphidx = 0.0 ; dphidz = 0.0

! if (igrid.eq.2) then
! write(6,*)'homoy f',igrid,field(172,163),field(172,164)
! endif

  DO ipoint = 1,nlocp(igrid)

     DO i = 3, ii-1
        IF (x(i) .ge. locpnts(1,ipoint,igrid)) THEN
           ip = i-1
           exit
        END IF
     END DO

     DO k = 3,kk-1
        IF (z(k) .ge. locpnts(3,ipoint,igrid)) THEN
           kp = k-1
           exit
        END IF
     END DO



        if (igrid.eq.2) then
           if (ipoint.eq.44) then
              write(6,*)'loco',igrid,ipoint,kp,ip,field(kp,ip),field(kp,ip+1)
           endif
        else
           if (ipoint.eq.1) then
              write(6,*)'loco',igrid,ipoint,kp,ip,field(kp,ip),field(kp,ip+1)
           endif
        endif

     ufield(1) = field(kp, ip )
     ufield(2) = field(kp ,ip+1)
     ufield(3) = field(kp+1,ip )
! ufield(2) = field(kp+1,ip)
! ufield(3) = field(kp,ip+1 )
     ufield(4) = field(kp+1,ip+1)

     facx = (locpnts(1,ipoint,igrid)-x(ip))/(x(ip+1)-x(ip))
     facz = (locpnts(3,ipoint,igrid)-z(kp))/(z(kp+1)-z(kp))

     IF (facx .LT. 0.0 .OR. facx .GT. 1.0 .OR. &
         facz .LT. 0.0 .OR. facz .GT. 1.0 ) THEN
        WRITE(6,*)'interpxz failed for ipoint',igrid,ipoint
! CALL ERRR (506,' SOLINT ')
     ENDIF

     dfacx = 1.0/(x(ip+1)-x(ip))
     dfacz = 1.0/(z(kp+1)-z(kp))

     psi(1) = (1.0-facx)*(1.0-facz)
     psi(2) = facx *(1.0-facz)
     psi(3) = (1.0-facx)* facz
     psi(4) = facx * facz

     phi(ipoint) = psi(1) * ufield(1) &
                   + psi(2) * ufield(2) &
                   + psi(3) * ufield(3) &
                   + psi(4) * ufield(4)

     psi(1) = -dfacx*(1.0-facz)
     psi(2) = dfacx*(1.0-facz)
     psi(3) = -dfacx* facz
     psi(4) = dfacx* facz

     dphidx(ipoint) = psi(1) * ufield(1) &
                     + psi(2) * ufield(2) &
                     + psi(3) * ufield(3) &
                     + psi(4) * ufield(4)

     psi(1) = (1.0-facx)* (-dfacz)
     psi(2) = facx * (-dfacz)
     psi(3) = (1.0-facx)* dfacz
     psi(4) = facx * dfacz

     dphidz(ipoint) = psi(1) * ufield(1) &
                     + psi(2) * ufield(2) &
                     + psi(3) * ufield(3) &
                     + psi(4) * ufield(4)

  END DO ! IPOINT

END SUBROUTINE findlocxz

SUBROUTINE homoy(kk,jj,ii,x,z,u,v,w,p,uav,vav,wav,pav,forma,igrid)

  IMPLICIT NONE
  INTEGER kk,jj,ii,k,j,i,ierr,ifile,ios,igrid,ja,je

  REAL u(kk,jj,ii),v(kk,jj,ii),w(kk,jj,ii),p(kk,jj,ii),&
       uav(kk,ii),vav(kk,ii),wav(kk,ii),pav(kk,ii),rnhmg,x(ii),z(kk)

  CHARACTER (LEN=10) forma
  CHARACTER (LEN=20) outfile

  LOGICAL firstcall,ex

  DATA firstcall /.TRUE./
  SAVE :: firstcall,ifile

  IF (jj.EQ.1) THEN
     rnhmg = 1
     ja = 1
     je = 1
  ELSE
     ja = 3
     je = jj-2
     rnhmg = 1.0/real(jj-4)
  END IF


  uav = 0.0
  vav = 0.0
  wav = 0.0
  pav = 0.0

  DO i = 1,ii
     DO k=1,kk
        DO j = ja,je

           uav(k,i) = uav(k,i) + u(k,j,i)
           vav(k,i) = vav(k,i) + v(k,j,i)
           wav(k,i) = wav(k,i) + w(k,j,i)
           pav(k,i) = pav(k,i) + p(k,j,i)
           !write(11,*) i,k,uav(k,i),pav(k,i)

        END DO

        uav(k,i) = uav(k,i) * rnhmg
        vav(k,i) = vav(k,i) * rnhmg
        wav(k,i) = wav(k,i) * rnhmg
        pav(k,i) = pav(k,i) * rnhmg

     END DO
  END DO


END SUBROUTINE homoy

SUBROUTINE cfcdxz(gmol,time,igrid)

  IMPLICIT NONE



      INTEGER MAXGRIDS,MAXBOCONDS
      PARAMETER ( MAXGRIDS =128 )
      PARAMETER ( MAXBOCONDS = 10 )


  INTEGER,PARAMETER :: IDIMP = 5000

  COMMON /COSOLINT/ NLOCP,LOCPNTS,DXLOC,&
       SOLOUT,DIST,STAG,SHIFT

  INTEGER NLOCP(MAXGRIDS),SHIFT

  LOGICAL STAG

  REAL LOCPNTS(3,IDIMP,MAXGRIDS),SOLOUT(4,4,IDIMP,MAXGRIDS),&
       DXLOC(3,IDIMP,MAXGRIDS),DIST

  INTEGER i,j,igrid

  REAL tx,tz,nx,nz,dudx,dudz,dwdx,dwdz,area,areasum,cx,cz,&
       tauw,tauw2,p,rho,gmol,time,cf(nlocp(igrid)),&
       cp(nlocp(igrid)),rstau,cx2,cz2

  LOGICAL found
  character (len=11) outfile

  rho = 1.0

  cx = 0.0 ; cz = 0.0 ; areasum = 0.0

! All Points exept the last one
! write(*,*) nlocp(igrid)
  DO i = 1, (nlocp(igrid)/2)-1
     !tx = locpnts(1,i,igrid)-locpnts(1,i+1,igrid)
     !tz = locpnts(3,i,igrid)-locpnts(3,i+1,igrid)
     tx = -2*dxloc(1,i,igrid)
     tz = -2*dxloc(3,i,igrid)
     area = SQRT(tx**2.0+tz**2.0)

     tx = tx/area
     tz = tz/area

     nx = -tz
     nz = tx

     ! dudx = 0.5*(solout(1,2,i,igrid)+solout(1,2,i+1,igrid))
     ! dudz = 0.5*(solout(1,4,i,igrid)+solout(1,4,i+1,igrid))
     ! dwdx = 0.5*(solout(3,2,i,igrid)+solout(3,2,i+1,igrid))
     ! dwdz = 0.5*(solout(3,4,i,igrid)+solout(3,4,i+1,igrid))
     ! p = 0.5*(solout(4,1,i,igrid)+solout(4,1,i+1,igrid))

     dudx = solout(1,2,i,igrid)
     dudz = solout(1,4,i,igrid)
     dwdx = solout(3,2,i,igrid)
     dwdz = solout(3,4,i,igrid)
     p = solout(4,1,i,igrid)


     tauw = rho * gmol * ( (dudx*tx+dwdx*tz) * nx &
                           +(dudz*tx+dwdz*tz) * nz )

     tauw2 = -rho*gmol*((solout(1,1,i+shift,igrid))*tx + &
                        (solout(3,1,i+shift,igrid))*tz) / dist

     !write(10,*) locpnts(1,i,igrid),solout(1,1,i+shift,igrid)/20.0,tauw,tauw2
     !write(15,*) locpnts(1,i,igrid),solout(4,1,i,igrid),solout(1,1,i,igrid)
     ! tauw = tauw2

     cf(i) = tauw
     cp(i) = p
     cx = cx + (tauw*nz - p*nx) * area
     cz = cz - (tauw*nx + p*nz) * area
     cx2 = cx2 + (tauw2*nz - p*nx) * area
     cz2 = cz2 - (tauw2*nx + p*nz) * area
     areasum = areasum + area
     !write(16,*) locpnts(1,i,igrid),cx,cz,areasum

  END DO

  rstau = 2.0 / (rho * (0.2**2.0) * 0.10)

  write(*,*) 'cx,cz',cx,cz,areasum,rstau
  write(*,*) 'cd,cl',cx * rstau ,cz * rstau
  write(*,*) 'cd2,cl2',cx2 * rstau ,cz2 * rstau

  outfile = 'cfcd.grd.00'

  WRITE(outfile(10:11),'(i2.2)') igrid

  open (unit = 11, file=outfile, form = 'formatted')

  WRITE(11,*)

  DO i=1,nlocp(igrid)-1
     WRITE(11,*) i,locpnts(1,i,igrid),locpnts(3,i,igrid),cf(i),cp(i)
  ENDDO
  WRITE(11,*)

  CLOSE(11)
  RETURN

END SUBROUTINE cfcdxz
