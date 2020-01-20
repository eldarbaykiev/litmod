! Subroutine SUB_trigo.f 
!Geographical coordinates in º (internally tranformed into rads)


c***********************************************
      SUBROUTINE trigo(a,b,c,d,delta,dist,acim)
       real(4) a,b,c,d,pi,R_T,long,delta,dist,indic,acim,indic_ac
       parameter (pi=3.14159,R_T=6370)
       integer j
!Profile coordinates (long,lat): (c,a)-->(d,b)
!Output angular distance(rads), distance (km), acimuth(rads)
!delta,dist,acim
       a=a*pi/180
       b=b*pi/180
       c=c*pi/180
       d=d*pi/180

C j es un parámetro de control de flujo del programa
       j=0
       
      if (c.GE.0) then
        if (d.GE.0) then
         long=abs(c-d)  
         delta=acos(cos(pi/2-a)*cos(pi/2-b)+sin(pi/2-a)*sin(pi/2-b)*co
     *s(long))
         dist=R_T*delta
C         acim=asin(cos(b)*sin(long)/sin(delta))
         indic=0
c         write(*,*) 'Los dos puntos quedan al E del meridiano 0'
        else
c El punto 2 queda al Oeste del punto 1
        
c          write(*,*) 'El punto 2 queda al Oeste del punto 1' 
          long=c+abs(d)
          delta=acos(cos(pi/2-a)*cos(pi/2-b)+sin(pi/2-a)*sin(pi/2-b)*co
     *s(long))
          dist=R_T*delta

          if ((cos(b)*sin(long)/sin(delta))>1) then
          acim=3*pi/2
          j=1
          else 
          acim=asin(cos(b)*sin(long)/sin(delta))
          endif


           indic_ac=1
          indic=1
        endif


       else
        if (d.LE.0) then
c         write(*,*) 'Los dos puntos quedan al Oeste del meridiano 0'       
         long=abs(c-d)
         delta=acos(cos(pi/2-a)*cos(pi/2-b)+sin(pi/2-a)*sin(pi/2-b)*co
     *s(long))
         dist=R_T*delta
C         acim=asin(cos(b)*sin(long)/sin(delta))
         indic=0
        else
c  El punto 2 queda al Este del punto 1
c          write(*,*) 'El punto 2 queda al Este del punto 1'
         long=d+abs(c)
         delta=acos(cos(pi/2-a)*cos(pi/2-b)+sin(pi/2-a)*sin(pi/2-b)*co
     *s(long))
         dist=R_T*delta

          if ((cos(b)*sin(long)/sin(delta))>1) then
          acim=pi/2
          j=1 
          else
          acim=asin(cos(b)*sin(long)/sin(delta))
          endif


         indic=1
         indic_ac=0
        endif

       endif
     
      if (indic.EQ.0) then
        if (c.GE.d) then
c         write(*,*) 'El punto 2 queda al Oeste del punto 1'
          if ((cos(b)*sin(long)/sin(delta))>1) then
          acim=3*pi/2
          j=1
          else
          acim=asin(cos(b)*sin(long)/sin(delta))
          endif

         indic_ac=1

        else
c        write(*,*) 'El punto 2 queda al Este del punto 1' 
          if ((cos(b)*sin(long)/sin(delta))>1) then
          acim=pi/2
          j=1
          else
          acim=asin(cos(b)*sin(long)/sin(delta))
          endif


         indic_ac=0
        endif
      endif
C*****************************************
C El programa mira en que cuadrante se encuentra el punto 2
C tomando el origen en el punto 1, para calcular el acimut
C correctamente      

      if (j==1) then
C La situación corresponde a los casos en 
C que es (cos(b)*sin(long)/sin(delta))>1, es decir, cuando 
C se impone acim=pi/2 ó 3*pi/2. Estos valores no necesitan 
C ser modificados según el cuadrante.
      else 

      if (indic_ac.EQ.1) then
      
       if (cos(pi/2-b)-cos(pi/2-a)*cos(delta).LE.0) then
c Punto en el 3º cuadrante
         acim=pi+acim
       else
c  Punto en el 4º cuadrante
         acim=2*pi-acim
       endif

      else
       
       if (cos(pi/2-b)-cos(pi/2-a)*cos(delta).LE.0) then
c  Punto en el 2º cuadrante
         acim=pi-acim
       endif
c  Punto en el 1º cuadrante
      endif  

      endif
      
       a=a*180/pi
       b=b*180/pi
       c=c*180/pi
       d=d*180/pi       


      return
      END  SUBROUTINE
