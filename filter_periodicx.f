










CCCCC        DEFINITIONEN FUER C-PREPROZESSOR
C            MASCHINE

C              MESSAGE PASSING INTERFACE

C            RAEUMLICHE DISKRETISIERUNG
C
***********************************************************************
C                       KOMPAKT-UPWIND in X-RICHTUNG (3-TER ORDNUNG) fuer
C                       die U-Komponente (V und W bleiben Kompakt 4ter ordnung)
C                       falls im Stroemungsfeld eine Koerper vorhanden ist
C
***********************************************************************
C                       KOMPAKTVERF. in XYZ-RICHTUNG (4-TER ORDNUNG)


**********************************************************************
C                      Preprocessing with ADM for 2nd Order Central
***********************************************************************
CCC                     Bei periodischen Randbedingungen in X-Richtung
CCC                     ist eine hoehere O
C
CCC                     Bei periodischen Randbedingungen in Y-Richtung
CCC                     ist eine hoehere Ordnung moeglich
C
CCC                     Bei periodischen Randbedingungen in Y-Richtung:
CCC                     Zentraldiff. 4-ter Ordnung (Parallelisieung moeglich)
C
***********************************************************************
C
C
C            INTERPOLATION AN GITTER-GRENZEN

C            BEHANDLUNG DER TOPAR-RANDBEDINGUNG

C            ZEITLICHE DISKRETISIERUNG

C            FEINSTRUKTURMODELL

C            NICHT-NEWTONSCHE SPANNUNGEN


C            TRANSPORT UND ORIENTIERUNG VON PARTIKELN


C            SCALAR HEAT/TEMPERATURE TRANSPORT
C            APROXIMATE DECONVOLUTION FOR SCALAR
C            TURBULENT PRANDTL NUMBER 
C            PLOT LOCAL SCALAR CONVECTION DIFFUSION EXTREMES 
C            ANALYZE AND PLOT NEAR WALL GRID RESLOLUTION 
C            WRITE SPECIAL 1D LINE FOR TMIX FOR SPECTRA AND PDF
C            SCALAR TIME ADVANCEMENT/DISCRETISATION METHOD

C            STROEMUNG NACH OBEN?


C            BEHANDLUNG DER FLUKTUATIONS-RANDBEDINGUNG

C            BEHANDLUNG VON PPHYS ALS QUELLTERM IN TSTLE2

C            POSITIONIERUNG DER EINSTROEMPROFILE

C           STATISTIK


C                 GEMISCHTES

C                GRDFMI WIRD NICHT VERWENDET, DAHER EINSPARUNG DER FELDER
C                IGRHF UND GRHF

C                VERGROESSERN DER KJI-RICHTUNG BEI FORTSETZUNGSLAUF ERLAUBT!

C                RAUSSCHREIBEN VON ZEITRECORDS

C                FUER KORRELATIONEN UND HAEFIGKEITSVERTEILUNGEN
      subroutine filter_periodicx(kk,jj,ii,T,HELP1,HELP2,wsavex,feldx,
     $           filterx)
C-----------------------------------------------------------------
C      Periodischer filter der die FFT-Coefficients mit
C      der in filterx gespeicherten transfer-funktion multipliziert
C
C      FFT aus FFTPACK mit dependecies kopiert
C      17.06.04   Florian Schwertfirm
C-----------------------------------------------------------------
       implicit none
       integer nx,ny,nz,ii,jj,kk
       integer l,m,NCOUNT,k,j,i
       real T(kk,jj,ii),pi

       real HELP1((kk)*(jj)*(ii)),HELP2((kk)*(jj)*(ii))
       real filterx((ii-4)/2+1)
       real wsavex(2*(ii-4)+15)
       real feldx((ii-4))
C       real*8 feldxdp((ii-4))
       integer,parameter::singleprec = SELECTED_REAL_KIND(4)
       integer,parameter::doubleprec = SELECTED_REAL_KIND(8)

       nx=ii-4
       ny=jj-4
       nz=kk-4
       NCOUNT = 0
       pi = acos(-1.0)

C--- DEBUG
C       do i=1,nx/2+1
C        filterx(i) = 1.0
C       enddo
C       do k=1,nz
C        do j=1,ny
C         do i=1,nx
C          T(k+2,j+2,i+2) = cos(2.0*pi*float((i-1))/float(nx)) + 
C     $     2.0*cos(4.0*pi*float((i-1))/float(nx))
C     $   + 3.0*cos(6.0*pi*float((i-1))/float(nx))
C     $   + sin(2.0*pi*float((i-1))/float(nx))
C     $   + 2.0*sin(4.0*pi*float((i-1))/float(nx))
C     $   + 3.0*sin(6.0*pi*float((i-1))/float(nx))
C        T(k+2,j+2,i+2) = 0.0
C         enddo
C        enddo
C       enddo
C       do k=1,nz
C        do j=1,ny
C         T(k+2,j+2,3)=1.0
C        enddo
C       enddo     
C      do i=1,nx
C         write(20,*)i,T(10,10,i+2)
C       enddo


C------ i-direction must be innermost
       do k=1,nz
        do j=1,ny
         do i=1,nx
         NCOUNT = NCOUNT + 1 
         HELP1(NCOUNT) = T(k+2,j+2,i+2)
         enddo
        enddo
       enddo
C------ fft initialize

        call rffti(nx,wsavex)

C----- forward fft in nx-stripes
        do l=1,(nz*ny)
         do i = 1,nx
C          feldxdp(i) = dble(HELP1((l-1)*nx+i))
          feldx(i) = HELP1((l-1)*nx+i)
         enddo
C         call dfftf(nx,feldxdp,wsavex)
         call rfftf(nx,feldx,wsavex)
         do i = 1,nx
C          HELP1((l-1)*nx+i) = real(feldxdp(i),singleprec)
          HELP1((l-1)*nx+i) = feldx(i)
         enddo
        enddo

C----- Multiplikation with transfere funktion         
 
       do l=1,(nz*ny)
        i = 1
C----- Behandlung des nullten mode
        HELP2((l-1)*nx+i) = HELP1((l-1)*nx+i)*filterx(i)
C----- restliche mode
        do i = 2,nx/2
            HELP2((l-1)*nx+(2*i-2)) = 
     $                           HELP1((l-1)*nx+(2*i-2))*filterx(i)
            HELP2((l-1)*nx+(2*i-1)) = 
     $                           HELP1((l-1)*nx+(2*i-1))*filterx(i)
        enddo
C----- letzter mode
        i = nx
        HELP2((l-1)*nx+(i)) = 
     $                           HELP1((l-1)*nx+(i))*filterx(i/2+1)
       enddo

C----- backward fft
        do l=1,(ny*nz)
         do i=1,nx
C          feldxdp(i) = dble(HELP2(l+(nz*ny)*(i-1)))
          feldx(i) = HELP2((l-1)*nx+i)
         enddo
C         call dfftb(nx,feldxdp,wsavex)
         call rfftb(nx,feldx,wsavex)
         do i=1,nx
C          HELP2(l+(nz*ny)*(i-1)) = real(feldxdp(i),singleprec)/float(nx)
          HELP2((l-1)*nx+i) = feldx(i)/float(nx)
         enddo
        enddo
C----- Transpose bak into T-Field
       NCOUNT = 0
       do k=1,nz
        do j=1,ny
         do i=1,nx
         NCOUNT = NCOUNT + 1
         T(k+2,j+2,i+2) = HELP2(NCOUNT)
         enddo
        enddo
       enddo
 
       RETURN
       END 
