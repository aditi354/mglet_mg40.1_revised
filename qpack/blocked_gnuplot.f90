subroutine blocked_gnuplot(blocked,ni,nj,nk,x,y,z,xc,yc,zc,maccur)

  implicit none
  ! aus subroutine
  integer ni,nj,nk
  real  x(ni+1),y(nj+1),z(nk+1) ,blocked(ni,nj,nk)
  real  xc(ni),yc(nj),zc(nk), maccur
  ! lokale variablen
  integer direction,a,b,c,ip,tmpip,itriangle
integer i,j,k, npoints, ntriangles, idummy, ncells
integer number, nedges, ed_number, tri_point, ndoubles
integer :: ncalls=1
real rdummy,xmax_body(3,2)
character*20  cdummy
type point
  integer :: number
  integer :: flag,position
  real    :: x, y, z
end type point
type triangle
  integer :: number
  integer :: p1, p2, p3
end type triangle
type edge
  integer :: number
  integer :: p1, p2
end type edge
type surface
  integer :: number
  integer :: ed1, ed2, ed3
end type surface
type(point) points, sorted_points 
type(triangle) triangles
type(edge) edges
type(surface) surfaces
allocatable :: points(:), triangles(:), edges(:), surfaces(:)
allocatable :: sorted_points(:)


if (ncalls.eq.3) then
  
direction = 1
!call blocked_plot(blocked,ni,nj,nk,x,y,z,xc,yc,zc,direction)
direction = 2
!call blocked_plot(blocked,ni,nj,nk,x,y,z,xc,yc,zc,direction)
direction = 3
!call blocked_plot(blocked,ni,nj,nk,x,y,z,xc,yc,zc,direction)
  
end if

if (ncalls.eq.4) then

ncells = 0
do i = 1, ni
   do j = 1, nj
      do k=1, nk
      ! Wenn eine Geometriezelle
      if ((blocked(i,j,k)).lt.0) then
      ncells = ncells + 1
      end if
      end do
   end do
end do
write(130,*) ncells


ntriangles = 12 * ncells
npoints = 8 * ncells
!nsurfaces = 12 * ncells
!nedges =  3 * ntriangles * ncells

! Speicher allocieren
allocate(points(npoints),triangles(ntriangles))
!         edges(nedges), surfaces(nedges))

ip = 1
itriangle = 1
tmpip = 0
ncells = 1
do i = 1, ni
   do j = 1, nj
      do k = 1, nk
         write(135,*) x(i), y(j), z(k)
         ! Wenn eine Geometriezelle
         if ((blocked(i,j,k)).lt.0) then
            write(130,*) ip  
            ! 8 Quadereckpunkte erstellen
            
            do a = 1,2
               
               do b = 1,2
                   
                  do c = 1,2
                       points(ip) % number = ip
                       points(ip) % position = (ni*nj)*((k+c-1)-1)+(ni)*((j+b-1)-1)+(i+a-1)
                       points(ip) % x = x(i+a-1) 
                       points(ip) % y = y(j+b-1) 
                       points(ip) % z = z(k+c-1) 
                       if (ncells.lt.9) write(130,*) x(i+a),(i+a), y(j+b),(j+b), z(k+c),(k+c)
                       if (c.lt.2) ip = ip + 1
                       
                  if (ncells.lt.9) write(130,*) 'c',ip, c
                  end do 
                  if (b.lt.2) ip = ip + 1
                  
                  if (ncells.lt.9) write(130,*) 'b', ip, b
               end do  
               if (a.lt.2) ip = ip + 1
               
               if (ncells.lt.9) write(130,*) 'a', ip, a
            end do
            ! Die Dreiecke erstellen 
            ! 1. Dreieck
            ip = (ip - 7)
            if (ncells.lt.9) write(130,*) ip  
            ncells = ncells + 1 
   
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip) % number  
            triangles(itriangle) % p2 = points(ip+1) % number  
            triangles(itriangle) % p3 = points(ip+2) % number  
            ! 2. Dreieck
            itriangle = itriangle + 1
   
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip+1) % number  
            triangles(itriangle) % p2 = points(ip+2) % number  
            triangles(itriangle) % p3 = points(ip+3) % number  
            ! 3. Dreieck      
            itriangle = itriangle + 1
    
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip+4) % number  
            triangles(itriangle) % p2 = points(ip+5) % number  
            triangles(itriangle) % p3 = points(ip+6) % number  
            ! 4. Dreieck
            itriangle = itriangle + 1
   
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip+5) % number  
            triangles(itriangle) % p2 = points(ip+6) % number  
            triangles(itriangle) % p3 = points(ip+7) % number  
            ! 5. Dreieck
            itriangle = itriangle + 1
   
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip) % number  
            triangles(itriangle) % p2 = points(ip+2) % number  
            triangles(itriangle) % p3 = points(ip+4) % number  
            ! 6. Dreieck
            itriangle = itriangle + 1
   
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip+2) % number  
            triangles(itriangle) % p2 = points(ip+4) % number  
            triangles(itriangle) % p3 = points(ip+6) % number  
            ! 7. Dreieck
            itriangle = itriangle + 1
   
            triangles(itriangle) % number = itriangle
            triangles(itriangle) % p1 = points(ip+1) % number  
            triangles(itriangle) % p2 = points(ip+3) % number  
            triangles(itriangle) % p3 = points(ip+5) % number  
            ! 8. Dreieck
            itriangle = itriangle + 1
    
            triangles(itriangle) % number = itriangle
            triangles(itriangle) % p1 = points(ip+3) % number  
            triangles(itriangle) % p2 = points(ip+5) % number  
            triangles(itriangle) % p3 = points(ip+7) % number  
            ! 9. Dreieck
            itriangle = itriangle + 1
   
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip+2) % number  
            triangles(itriangle) % p2 = points(ip+3) % number  
            triangles(itriangle) % p3 = points(ip+6) % number  
            ! 10. Dreieck
            itriangle = itriangle + 1
   
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip+3) % number  
            triangles(itriangle) % p2 = points(ip+6) % number  
            triangles(itriangle) % p3 = points(ip+7) % number  
            ! 11. Dreieck
            itriangle = itriangle + 1
   
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip) % number  
            triangles(itriangle) % p2 = points(ip+1) % number  
            triangles(itriangle) % p3 = points(ip+4) % number  
            ! 12. Dreieck
            itriangle = itriangle + 1
    
            triangles(itriangle) % number = itriangle 
            triangles(itriangle) % p1 = points(ip+1) % number  
            triangles(itriangle) % p2 = points(ip+4) % number  
            triangles(itriangle) % p3 = points(ip+5) % number  
          
            itriangle = itriangle + 1
            ip = ip + 8        
         
         end if
      end do
   end do
end do 
 
! Doppelte Punkte Verweis aendern
! Schoen aber langsam wegen vergleichen
!write(*,*) 'Verweise in Punkten aendern'
!do i=1,npoints
!   write(*,*) i
!   points(i) % flag = 1
!   do j=1,i  
!      if (j.ne.i) then 
!         if ((points(j) % flag).ne.0) then
!            if ((((points(j) % x)-maccur).le.(points(i) % x)).and.&
!                (((points(j) % x)+maccur).ge.(points(i) % x))) then
!               if ((((points(j) % y)-maccur).le.(points(i) % y)).and.&
!                  (((points(j) % y)+maccur).ge.(points(i) % y))) then
!                  if ((((points(j) % z)-maccur).le.(points(i) % z)).and.&
!                     (((points(j) % z)+maccur).ge.(points(i) % z))) then
!                     points(i) % number = points(j) % number
!                     points(i) % flag = 0
!                  end if
!               end if
!            end if
!         end if
!      end if
!   end do
!end do        


! Doppelte Punkte Verweis aendern
write(*,*) 'Verweise in Punkten aendern'
write(*,*) 'npoints', npoints
do i=1,npoints
   ! write(*,*) i
   points(i) % flag = 1
   do j=1,i  
      if (j.ne.i) then 
         if ((points(j) % flag).ne.0) then
            if ((points(j) % position).eq.(points(i) % position)) then
               points(i) % number = points(j) % number
               points(i) % flag = 0
            end if
         end if
      end if
   end do
end do        

! Zaehlen der doppelten Punkte        
ndoubles = 0
do i=1,npoints
   if (points(i) % flag.eq.0) then
      ndoubles = ndoubles + 1
   end if
end do

      write(*,*)    ndoubles, npoints 

! Punktnummern in den Dreiecken entsprechend umschreiben
write(*,*) 'Punkte in Dreiecken umschreiben'
do i=1, ntriangles
   
   if ((points(triangles(i) % p1) % number).ne.&
       (triangles(i) % p1)) then
      triangles(i) % p1 = points(triangles(i) % p1) % number
   end if
   if ((points(triangles(i) % p2) % number).ne.&
       (triangles(i) % p2)) then
      triangles(i) % p2 = points(triangles(i) % p2) % number
   end if
   if ((points(triangles(i) % p3) % number).ne.&
       (triangles(i) % p3)) then
      triangles(i) % p3 = points(triangles(i) % p3) % number
   end if
end do
 

! Speicher allocieren
allocate(sorted_points(npoints-ndoubles)) 


! Jetzt koennen doppelte Elemente 
! in eine neue Liste geschrieben werden
write(*,*) 'Neue Liste schreiben'
j = 0
do i=1, npoints
   if (points(i) % flag.ne.0) then
      j = j + 1
      sorted_points(j) % number = points(i) %  number

      sorted_points(j) % x = points(i) %  x
      sorted_points(j) % y = points(i) %  y
      sorted_points(j) % z = points(i) %  z
      sorted_points(j) % flag = 1
   end if
end do
   


! Dreiecks-Verweise auf die neue Liste umschreiben
write(*,*) 'Dreiecks Verweise auf neue Liste umschreiben'
do i=1, ntriangles
   do j=1, npoints-ndoubles
      if ((triangles(i) % p1).eq.&
         (sorted_points(j) % number)) then
         triangles(i) % p1 = j
      end if
      if ((triangles(i) % p2).eq.&
         (sorted_points(j) % number)) then
         triangles(i) % p2 = j
      end if
      if ((triangles(i) % p3).eq.&
         (sorted_points(j) % number)) then
         triangles(i) % p3 = j
      end if    
   end do
end do

! Punktnummern in der neuen Liste noch anpassen
write(*,*) 'Punktnummern in neuer Liste anpassen'
do i=1, (npoints-ndoubles)
   sorted_points(i) % number = i
end do


! Zusaetzlich: Funktionalitaet um
! gts file fuer meshviewer zu erstellen

! Anzahl der Punkte, Kanten und Flaechen
! Flaechen sind auch hier Dreiecke...
write(125,*) (npoints-ndoubles), (ntriangles*3), ntriangles

! Liste mit den Punkten (angefangen bei "1")
!do i=1, ntriangles
!   number = triangles(i) % p1
!   write(125,*) points(number) % x, &
!                points(number) % y, &
!                points(number) % z
!   number = triangles(i) % p2
!   write(125,*) points(number) % x, &
!                points(number) % y, &
!                points(number) % z
!   number = triangles(i) % p3
!   write(125,*) points(number) % x, &
!                points(number) % y, &
!                points(number) % z
!end do
! Neu: weil sonst so grosses file entsteht
do i=1, (npoints-ndoubles)
   write(125,*) sorted_points(i) % x, &
                sorted_points(i) % y, &
                sorted_points(i) % z
end do


! Liste mit den Kanten
!i = 1
!do while (i.le.(ntriangles*3))
!   write(125,*) (i),  (i+1)
!   write(125,*) (i+1),(i+2)
!   write(125,*) (i+2),(i)
!   i = i + 3
! end do
! Neu wegen filegroesse ...
do i=1, ntriangles
   write(125,*) triangles(i) % p1, &
                triangles(i) % p2
   write(125,*) triangles(i) % p2, &
                triangles(i) % p3 
   write(125,*) triangles(i) % p3, &
                triangles(i) % p1
end do

i = 1
! Liste mit Verknüfung der Kanten...
do while (i.le.(ntriangles*3))
   write(125,*) (i),(i+1),(i+2)
   i = i + 3
end do


end if
ncalls = ncalls + 1

   
contains

subroutine blocked_plot(blocked,ni,nj,nk,x,y,z,xc,yc,zc,direction)


  implicit none
  ! aus subroutine
  integer ni,nj,nk
  real  x(ni+1),y(nj+1),z(nk+1) ,blocked(ni,nj,nk)
  real  xc(ni),yc(nj),zc(nk)
  ! lokale variablen
  integer i,j,k,a,b, number, counter, nchar
  integer i1, i2, i3, ni1, ni2, ni3, direction
  !integer :: ncalls=1
  character(len=10) name
  character(len=4) ext,exe, plane
  character(len=4) num, format
  character(len=4) info
  character(len=6) script
  character(len=24) fileptr
  name = 'bd_gnuplot'
  ext = '.dat'
  exe = '.exe'
  info = 'info'
  script = 'script'
  plane = '_xz_'

  select case (direction)  
     case(1)  
        ni1 = ni
        ni2 = nj
        ni3 = nk   
        plane = '_yz_'
     case(2)  
        ni1 = nj
        ni2 = ni
        ni3 = nk
        plane = '_xz_'
     case(3)  
        ni1 = nk
        ni2 = ni
        ni3 = nj   
        plane = '_xy_' 
  end select
    

  !if (ncalls.eq.1) then
    
     do i1 = 1, ni1

        number = i1
        write(num(1:4),'(i4)') number
        fileptr = name//plane//num//ext
        do a = 1,8
           if (fileptr(14+1:14+1).eq.' ') then
              do b = 1, (9-a)
                 fileptr(14+b:14+b)=fileptr(14+b+1:14+b+1)
              end do
           end if
        end do
        ! Berechnet die benötigten Stellen fuer den Namen
        nchar = int(log10(float(number))+1) + 4
        open (unit=103,file=fileptr(1:14+nchar),form='FORMATTED')
        ! open (unit=103,file=filout(1:10+nchar))
        ! write(103,*) 'test', x(i)

        do i2 = 1, ni2
           do i3 = 1, ni3
              ! Wenn eine Geometriezelle
              select case (direction)
                 case (1)
                    i = i1
                    j = i2
                    k = i3
                 case (2)
                    i = i2
                    j = i1
                    k = i3
                 case (3)
                    i = i2
                    j = i3
                    k = i1
              end select
              if ((blocked(i,j,k)).lt.0) then
                 ! Quader zeichnen
                 ! erstes Viereck
                 write(103,*) x(i)  , y(j)  , z(k)
                 write(103,*) x(i+1), y(j)  , z(k)
                 write(103,*) x(i+1), y(j+1), z(k)
                 write(103,*) x(i)  , y(j+1), z(k)
                 write(103,*) x(i)  , y(j)  , z(k)
                 ! Zwei Leerzeilen
                 write(103,*) ' '
                 write(103,*) ' '
                 ! Zweites Viereck
                 write(103,*) x(i)  , y(j)  , z(k+1)
                 write(103,*) x(i+1), y(j)  , z(k+1)
                 write(103,*) x(i+1), y(j+1), z(k+1)
                 write(103,*) x(i)  , y(j+1), z(k+1)
                 write(103,*) x(i)  , y(j)  , z(k+1)
                 ! Zwei Leerzeilen
                 write(103,*) ' '
                 write(103,*) ' '
                 ! 4 Verbindungslinien
                 write(103,*) x(i)  , y(j)  , z(k)
                 write(103,*) x(i)  , y(j)  , z(k+1)
                 write(103,*) ' '
                 write(103,*) ' '
                 write(103,*) x(i+1), y(j)  , z(k)
                 write(103,*) x(i+1), y(j)  , z(k+1)
                 write(103,*) ' '
                 write(103,*) ' '
                 write(103,*) x(i+1), y(j+1), z(k)
                 write(103,*) x(i+1), y(j+1), z(k+1)
                 write(103,*) ' '
                 write(103,*) ' '
                 write(103,*) x(i)  , y(j+1), z(k)
                 write(103,*) x(i)  , y(j+1), z(k+1)
                 write(103,*) ' '
                 write(103,*) ' '
                 ! Noch einmal fuer 3d-Isometrie Plot
                 ! Quader zeichnen
                 ! erstes Viereck
                 open(104,file='bd_3d_gnuplot')
                 write(104,*) x(i)  , y(j)  , z(k)
                 write(104,*) x(i+1), y(j)  , z(k)
                 write(104,*) x(i+1), y(j+1), z(k)
                 write(104,*) x(i)  , y(j+1), z(k)
                 write(104,*) x(i)  , y(j)  , z(k)
                 ! Zwei Leerzeilen
                 write(104,*) ' '
                 write(104,*) ' '
                 ! Zweites Viereck
                 write(104,*) x(i)  , y(j)  , z(k+1)
                 write(104,*) x(i+1), y(j)  , z(k+1)
                 write(104,*) x(i+1), y(j+1), z(k+1)
                 write(104,*) x(i)  , y(j+1), z(k+1)
                 write(104,*) x(i)  , y(j)  , z(k+1)
                 ! Zwei Leerzeilen
                 write(104,*) ' '
                 write(104,*) ' '
                 ! 4 Verbindungslinien
                 write(104,*) x(i)  , y(j)  , z(k)
                 write(104,*) x(i)  , y(j)  , z(k+1)
                 write(104,*) ' '
                 write(104,*) ' '
                 write(104,*) x(i+1), y(j)  , z(k)
                 write(104,*) x(i+1), y(j)  , z(k+1)
                 write(104,*) ' '
                 write(104,*) ' '
                 write(104,*) x(i+1), y(j+1), z(k)
                 write(104,*) x(i+1), y(j+1), z(k+1)
                 write(104,*) ' '
                 write(104,*) ' '
                 write(104,*) x(i)  , y(j+1), z(k)
                 write(104,*) x(i)  , y(j+1), z(k+1)
                 write(104,*) ' '
                 write(104,*) ' '
              end if
           end do
        end do
        close(103)
     end do 
  !end if

  !ncalls = ncalls + 1


  ! Ausgabe in file "bd_gnuplot_xz_info.dat"
  fileptr = name//plane//info//ext
  open (unit=103,file=fileptr(1:22),form='FORMATTED')
  write(103,FMT="(A10)") '#!/bin/csh'
  write(103,*) '#Script fuer das automatische Erstellen'
  write(103,*) '#der Graphiken mit gnuplot'
  write(103,*) 'setenv ni' ,ni
  write(103,*) 'setenv nj' ,nj
  write(103,*) 'setenv nk' ,nk
  close(103)

  ! Ausgabe in file "bd_gnuplot_xz_script.exe"
  fileptr = name//plane//script//exe
  open (unit=103,file=fileptr(1:24),form='FORMATTED')
  write(103,FMT="(A10)") '#!/bin/csh'
  write(103,*) ' '
  write(103,*) '#Schleife fuer die einzelnen Bilder '
  write(103,*) ' '
  write(103,fmt="(A17,A4,A8)") 'source bd_gnuplot',plane,'info.dat'
  write(103,*) ' '
  select case (direction)
     case (1)
        write(103,*) '@ i = 1'
        write(103,*) 'while ( $i <= $ni ) '
        write(103,*) ' '
!        write(103,*) '   setenv dateiname  plot_yz_$i.ps'
!        write(103,*) '   setenv geometry   bd_gnuplot_yz_$i.dat'
        write(103,*) '   echo $i '
!        write(103,*) '   bd_gnuplot_yz_script.dat'
        write(103,*) '   @ i++'
     case (2)
        write(103,*) '@ j = 1'
        write(103,*) 'while ( $j <= $nj ) '
        write(103,*) ' '
!        write(103,*) '   setenv dateiname  plot_xz_$j.ps'
!        write(103,*) '   setenv geometry   bd_gnuplot_xz_$j.dat'
        write(103,*) '   echo $j '
!        write(103,*) '   bd_gnuplot_xz_script.dat'
        write(103,*) '   @ j++'
     case (3)
        write(103,*) '@ k = 1'
        write(103,*) 'while ( $k <= $nk ) '
        write(103,*) ' '
!        write(103,*) '   setenv dateiname  plot_xz_$k.ps'
!        write(103,*) '   setenv geometry   bd_gnuplot_xz_$k.dat'
        write(103,*) '   echo $k '
!        write(103,*) '   bd_gnuplot_xz_script.dat'
        write(103,*) '   @ k++'
  end select
  write(103,*) ' '
  write(103,*) 'end'
  write(103,*) ' '
  write(103,*) 'gnuplot <<EOF '
  write(103,*) '   set terminal postscript' 
  write(103,fmt="(A21,A4,A4)") '   set output "slices',plane,'.ps"'
!  write(103,fmt="(A19,A4,A23)") '   plot "bd_gnuplot',plane,'1.dat" u 1:3 with lines'
  do j=1,ni1
  if ((j.ge.0).and.(j.lt.10)) write(103,FMT="(A19,A4,I1,A22)") &
   '   plot "bd_gnuplot',plane,j,'.dat" u 1:3 with lines'
  if ((j.ge.10).and.(j.lt.100)) write(103,FMT="(A19,A4,I2,A22)") &
   '   plot "bd_gnuplot',plane,j,'.dat" u 1:3 with lines'
  if ((j.ge.100).and.(j.lt.1000)) write(103,FMT="(A19,A4,I3,A22)") &
   '   plot "bd_gnuplot',plane,j,'.dat" u 1:3 with lines'
  if ((j.ge.1000).and.(j.lt.10000)) write(103,FMT="(A19,A4,I4,A22)") &
   '   plot "bd_gnuplot',plane,j,'.dat" u 1:3 with lines'
  if ((j.ge.10000).and.(j.lt.100000)) write(103,FMT="(A19,A4,I5,A22)") &
   '   plot "bd_gnuplot',plane,j,'.dat" u 1:3 with lines'
  end do
  write(103,*) '   quit'
  write(103,fmt="(A3)") 'EOF'
  write(103,*) ' '
  write(103,fmt="(A16,A4,A6)") 'mkdir bd_gnuplot',plane,'slices'
  write(103,fmt="(A13,A4,A1,A13,A4,A6)") 'mv bd_gnuplot',plane,'*',' ./bd_gnuplot',plane,'slices'
  write(103,*) ' '
  write(103,*) ' '
  close(103)

  ! Ausgabe in file "bd_gnuplot_xz_script.dat"
  fileptr = name//plane//script//ext
  open (unit=103,file=fileptr(1:24),form='FORMATTED')
  write(103,FMT="(A10)") '#!/bin/csh'
  write(103,*) ' '
  write(103,*) 'echo $dateiname'
  write(103,*) 'gnuplot << EOF'
  write(103,*) 'set terminal postscript'
  write(103,*) 'set output "$dateiname"'
  write(103,*) 'plot "$geometry" u 1:3 with lines'
  write(103,*) 'quit'
  write(103,*) 'EOF'
  close(103) 

  ! Abschlussscript, das alle anderen startet
  if (direction.eq.3) then
     fileptr = name//script//exe
     open (unit=103,file=fileptr(1:24),form='FORMATTED')
     write(103,FMT="(A10)") '#!/bin/csh'
     write(103,*) ' '
     write(103,*) 'bd_gnuplot_xy_script.exe'
     write(103,*) 'bd_gnuplot_xz_script.exe'
     write(103,*) 'bd_gnuplot_yz_script.exe'
     write(103,*) ' '
     
  end if 


end subroutine

end subroutine
