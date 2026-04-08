










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
        subroutine brenner(rmu0,rmu1,rmu2,rmu3,rmu4,rnu,phi,asprat)
c
c       this subroutine serves to compute the particle stress 
c     coefficients and shape factor of the rigid, ellipsoidal particles 
c       from their aspect ratio.
c
c       coefficients in constitutive law:
c               
c               asprat                  aspect ratio of particles
c               rmu0,rmu1,rmu2,rmu3     particle stress coefficients
c               G                       shape factor of particles
c               rnu                     Brownian rotory diffusivity
c               phi                     particle volume fraction
c                                  (of order <= 1/asprat**2 for dilute)
c
c
c       CAUTION: THE ENTROPIC STRESS COEFFICIENT (rmu1) IS
c       CORRECT ONLY FOR ASPECT RATIO = 5 !!!!
c
C        implicit double precision (a-h,o-z)
c
c
c       first, based on the aspect ratio of the particles, 
c       select the correct qofr(r) = capital theta (aspect ratio)
c
        r = asprat
c------------------------------------------(r_p, S. 210 Brenner, 1974)
c------------------------------------------(qofr :: beta, S. 210 Br.)
        r1 = log(r + sqrt(r**2 - 1.0))

        if (r.le.1.0) then

                beta = r1/(r*sqrt(r**2 + 1.0))

        else

                beta = r1/(r*sqrt(r**2 - 1.0))

        end if

c        write (6,*) 'beta::::::::::::::::::' , beta
c
c
c       now compute the particle stress coefficients
c
        an = r**2    / (r**2 - 1.0) * (1.0 - beta)
        ap = 2.0 / (r**2 - 1.0) * (r**2 * beta - 1.0)

        an1 = r / (r**2 - 1.0)**2
     $       * (r**2 + 2.0 - 3.0*r**2*beta)

        ap1 = r**2  / (4.0*(r**2 - 1.0)**2)
     $       * (3.0*beta + 2.0*r**2 - 5.0)

        an2 = r**2  / (r**2 - 1.0)**2
     $       * ((2.0 * r**2 + 1.0)*beta - 3.0)

        ap2 = r**2  / (4.0*(r**2 - 1.0)**2)
     $       * (2.0*r**2 + 1.0 - 
     $             (4.0 * r**2 - 1.0) * beta)


        rk = 2.0*(r**2 + 1.0) / 
     $       (3.0*(r**2 * ap + an ))

        rn = 2.0*(r**2 - 1.0) / 
     $       (5.0*(r**2 * ap + an ))

        q1 = 0.2 / ap1

        q2 = 2.0/(15.0*ap1)
     $       * (1.0 - (ap2/an2))

        q3 = 1.0/(5.0*ap1)
     $       * (( (r*(ap + an)) / (r**2 * an + ap) * (ap1 / an1))
     $            - 1.0)

        q30 = 1.0/(5.0*ap1)
     $       * (( (2.0*r) / (r**2  + 1.0) * (ap1 / an1))
     $            - 1.0)

        rmu0 = 5.0*q1
        rmu1 = 5.0*q2
        rmu2 = 5.0*(-3.0*q2 - 4.0 * q30)
        rmu3 = 5.0*q30
        rmu4 = 5.0*rn*2.0*rnu

c        write (6,*) 'Q1, Q2, Q3, Q30: ',Q1, Q2, Q3, Q30

c        write (6,*) ' intrinsic viscosity for Pe = 0.0 : ', 
c     $       5.0*q1 - q2 + 2.0 * q3

        q3test = 0.5*( 6.0 - 5.0 * q1 + q2)
c        write (6,*) 'q3test:',q3test

c        write (6,*) ' intrinsic viscosity for Pe ->  + infty: ', 
c     $       5.0*(q1 - q2)

c        write (6,*) ' intrinsic viscosity for Pe ->  - infty: ', 
c     $       5.0*q1 - 1.2*q2
C
C---------------------------- Introduction of Volume fraction
C

        rmu0 = phi*rmu0
        rmu1 = phi*rmu1
        rmu2 = phi*rmu2
        rmu3 = phi*rmu3
        rmu4 = phi*rmu4


        return
        end
