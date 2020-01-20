! Subroutine SUB_in_out_pol.f to determine if the points
! of a regular grid are inside/outside a x-ordered or
! y-ordered polygon. One polygon is x(y)-ordered if the
! x(y) coordinate increases monotonally going through connected 
! vertex in the 'upper branch' from x(y)_min to x(y)_max and 
! decreases monotonally from x(y)_max to x(y)_min going
! through the 'lower branch'. x(y)-ordered is equivalent to say
! x(y)-convexe

      SUBROUTINE in_out_pol (opt)
      USE M_imag
      USE M_profile
      USE M_pol
      USE M_material
      USE M_sublitho_an
      REAL(4) xp_min,xp_max,yp_min,yp_max,dx_l,dx_r,dy_l,dy_r
      REAL(4) csx,csy,x_0,y_0,p,b,x,y,xp,yp
      REAL(4) vec_px,vec_py,deph_nl
      REAL(4) rho_low,rho_up,rho_c
      REAL(4), allocatable:: mark(:,:),v_dum(:,:)
      CHARACTER(16) UPDOWN,rot 
      INTEGER idummy,kx,ky  
      CHARACTER(1):: opt
      CHARACTER(2):: chdm
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       IF (ins>9) THEN
        WRITE(chdm,'(I2)')ins
       ELSE
        WRITE(chdm,'(I1)')ins
       ENDIF
       fil=TRIM(string)//TRIM(chdm)//'.xyz'
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      open(30,file=fil,status='UNKNOWN')
      x_0=lon_min
      y_0=lat_min
      ALLOCATE (mark(IDIM,JDIM))
      mark=0
      DO i=1,n_pol
       write(*,*)'POLYGON NUM-->',i 
       xp_min=minval(A_pol(2*i-1,1:NP_POL(i)))
       xp_max=maxval(A_pol(2*i-1,1:NP_POL(i)))
       yp_min=minval(A_pol(2*i,1:NP_POL(i)))
       yp_max=maxval(A_pol(2*i,1:NP_POL(i)))
       Ixp_min=1+nint((xp_min-x_0)/p_m_x)
       Ixp_max=1+nint((xp_max-x_0)/p_m_x) 
       Iyp_min=1+nint((yp_min-y_0)/p_m_y)
       Iyp_max=1+nint((yp_max-y_0)/p_m_y)
       IF(Ixp_min<1) Ixp_min=1
       IF(Ixp_max>IDIM) Ixp_max=IDIM 
       IF(Iyp_min<1) Iyp_min=1
       IF(Iyp_max>JDIM) Iyp_max=JDIM
       csx=0
       csy=0
!Check if the polygon x(y)-ordered
       DO j=3,NP_POL(i)
        dx_l=A_pol(2*i-1,j-1)-A_pol(2*i-1,j-2)
        dx_r=A_pol(2*i-1,j)-A_pol(2*i-1,j-1)
        dy_l=A_pol(2*i,j-1)-A_pol(2*i,j-2)
        dy_r=A_pol(2*i,j)-A_pol(2*i,j-1)
        IF ((dx_l>0.AND.dx_r<0).OR.(dx_l<0.AND.dx_r>0)) THEN
        csx=csx+1
        vec_px=-dx_l*dy_r+dy_l*dx_r
!         write(*,*)'vec_px',vec_px 
        ENDIF

        IF ((dy_l>0.AND.dy_r<0).OR.(dy_l<0.AND.dy_r>0)) THEN
        csy=csy+1
        vec_py=-dy_l*dx_r+dx_l*dy_r
!        write(*,*)'vec_p y',vec_py
        ENDIF 

        
       ENDDO
       IF (csx>2.AND.csy<=2) write(*,*)'THE POLYGON IS NOT X-ORDERED'
       IF (csy>2.AND.csx<=2) write(*,*)'THE POLYGON IS NOT Y-ORDERED'
       IF (csx>2.AND.csy>2) THEN
       write(*,*)'THE POLYGON IS NEITHER X-ORDERED NOR Y-ORDERED'
       write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@'
       mark(Ixp_min:Ixp_max,Iyp_min:Iyp_max)=3
       CYCLE
       ENDIF
!       write(*,*)'csx,csy',csx,csy
       IF (csx<=2) THEN
        if (vec_px>0) then
         rot='CW'
         else
         rot='CCW'
         endif
       ELSE
         if (vec_py>0) then
         rot='CCW'
         else
         rot='CW'
         endif
       ENDIF
       write(*,*)'ROTATION---->  ',rot
!If rotation is counter clock-wise re-arrange its coordinates 
       IF (rot=='CCW') THEN
       write(*,*)'REORDERING POLYGON VERTEX'
       ALLOCATE (v_dum(2,NP_POL(i))) 
       v_dum(1,:)=A_pol(2*i-1,1:NP_POL(i))
       v_dum(2,:)=A_pol(2*i,1:NP_POL(i))
        do il=1,NP_POL(i)
        A_pol(2*i-1,il)=v_dum(1,NP_POL(i)-il+1)
        A_pol(2*i,il)=v_dum(2,NP_POL(i)-il+1)
        enddo
       DEALLOCATE (v_dum)
       ENDIF
       write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@'
!********************************************************************
!Determine grid point in/outside of polygon 
         DO j=2,NP_POL(i)

         IF (csx<=2) THEN
!         write(*,*)'calculo por x'
!x-ordered polygon
!/////////////////////////////////////////////////////////
          Ix_l=1+nint((A_pol(2*i-1,j-1)-x_0)/p_m_x)
          Ix_r=1+nint((A_pol(2*i-1,j)-x_0)/p_m_x)
           IF (Ix_l==Ix_r) CYCLE
           IF (Ix_r>Ix_l) THEN
!From x_max to x_min
           p=(A_pol(2*i,j)-A_pol(2*i,j-1))/(A_pol(2*i-1,j)-
     *        A_pol(2*i-1,j-1))
           b=A_pol(2*i,j-1)-A_pol(2*i-1,j-1)*p 
           UPDOWN='up'
           ELSE
!From x_min to x_max
           p=(A_pol(2*i,j-1)-A_pol(2*i,j))/(A_pol(2*i-1,j-1)-
     *        A_pol(2*i-1,j))
           b=A_pol(2*i,j)-A_pol(2*i-1,j)*p
           idummy=Ix_l
           Ix_l=Ix_r
           Ix_r=idummy
           UPDOWN='down'
           ENDIF

          IF(Ix_l<1) Ix_l=1
          IF(Ix_r>IDIM) Ix_r=IDIM

          do kx=Ix_l,Ix_r
           x=(kx-1)*p_m_x+x_0
           yp=p*x+b
           do ky=Iyp_min,Iyp_max
            y=(ky-1)*p_m_y+y_0
           IF (y>=yp.AND.UPDOWN=='up') mark(kx,ky)=mark(kx,ky)-3
           IF (y<=yp.AND.UPDOWN=='down') mark(kx,ky)=mark(kx,ky)-3
            mark(kx,ky)=mark(kx,ky)+1 
           enddo
          enddo

        ELSEIF (csy<=2) THEN
!        write(*,*)'calculo por y'
!y-ordered polygon
!/////////////////////////////////////////////////////////
          Iy_l=1+nint((A_pol(2*i,j-1)-y_0)/p_m_y)
          Iy_r=1+nint((A_pol(2*i,j)-y_0)/p_m_y)
           IF (Iy_l==Iy_r) CYCLE
           IF (Iy_r>Iy_l) THEN
!From y_max to y_min
           p=(A_pol(2*i-1,j)-A_pol(2*i-1,j-1))/(A_pol(2*i,j)-
     *        A_pol(2*i,j-1))
           b=A_pol(2*i-1,j-1)-A_pol(2*i,j-1)*p
           UPDOWN='up'
           ELSE
!From y_min to y_max
           p=(A_pol(2*i-1,j-1)-A_pol(2*i-1,j))/(A_pol(2*i,j-1)-
     *        A_pol(2*i,j))
           b=A_pol(2*i-1,j)-A_pol(2*i,j)*p
           idummy=Iy_l
           Iy_l=Iy_r
           Iy_r=idummy
           UPDOWN='down'
           ENDIF
 
          IF(Iy_l<1) Iy_l=1
          IF(Iy_r>JDIM) Iy_r=JDIM 

          do ky=Iy_l,Iy_r
           y=(ky-1)*p_m_y+y_0
           xp=p*y+b
           do kx=Ixp_min,Ixp_max
           x=(kx-1)*p_m_x+x_0
           IF (x<=xp.AND.UPDOWN=='up') mark(kx,ky)=mark(kx,ky)-3
           IF (x>=xp.AND.UPDOWN=='down') mark(kx,ky)=mark(kx,ky)-3
           mark(kx,ky)=mark(kx,ky)+1
           enddo
          enddo
        ENDIF
!/////////////////////////////////////////////////////////
        
        ENDDO
!*****************************************************************      
      ENDDO 
       IF (opt=='A') THEN
!CALL FROM new_lay subroutine
         do ky=JDIM,1,-1
          y=(ky-1)*p_m_y+y_0
          do kx=1,IDIM
           x=(kx-1)*p_m_x+x_0
          IF (mark(kx,ky)>=2) THEN
           IF (ins==1)  THEN
            deph_nl=layers(kx,ky,ins)/2
           ELSE
            deph_nl=(layers(kx,ky,ins)+layers(kx,ky,ins-1))/2
           ENDIF
          layers(kx,ky,ins)=deph_nl
          write (30,*)x,y,deph_nl
          ELSE
            IF (i_top==0)  THEN
            deph_nl=-E(kx,ky)*1e-3
            ELSE
            deph_nl=layers(kx,ky,i_top)
            ENDIF
          layers(kx,ky,ins)=deph_nl
          write (30,*)x,y,deph_nl
          ENDIF
          enddo
         enddo
!CALL FROM two2three_lay subroutine
       ELSEIF (opt=='T') THEN
         open(1,file='./layers/layer1.xyz',status='UNKNOWN')
         open(2,file='./layers/layer2.xyz',status='UNKNOWN')
         open(3,file='./layers/layer3.xyz',status='UNKNOWN')
         rho_up=prop(2)%rho
         rho_low=prop(3)%rho
         rho_c=prop(1)%rho
          do ky=JDIM,1,-1
           y=(ky-1)*p_m_y+y_0
           do kx=1,IDIM
            x=(kx-1)*p_m_x+x_0
            IF ((mark(kx,ky)>=2.AND.cr_type=='O').OR.
     *          (mark(kx,ky)<2.AND.cr_type=='C')) THEN
             layers(kx,ky,2)=layers(kx,ky,1)
             layers(kx,ky,3)=layers(kx,ky,1) 
            ELSEIF ((mark(kx,ky)<2.AND.cr_type=='O').OR.
     *              (mark(kx,ky)>=2.AND.cr_type=='C')) THEN 
             layers(kx,ky,2)=((E(kx,ky)*1e-3+layers(kx,ky,1))*
     *                       rho_c-E(kx,ky)*1e-3*rho_up-
     *                       layers(kx,ky,1)*rho_low)/(rho_up-rho_low)
             layers(kx,ky,3)=layers(kx,ky,1)
             layers(kx,ky,1)=-E(kx,ky)*1e-3
            ENDIF      
            write (1,*)x,y,layers(kx,ky,1)
            write (2,*)x,y,layers(kx,ky,2)
            write (3,*)x,y,layers(kx,ky,3)
           enddo
          enddo
         CLOSE(1)
         CLOSE(2)
!CALL FROM join_lay subroutine
       ELSEIF (opt=='J') THEN
         do ky=JDIM,1,-1
          y=(ky-1)*p_m_y+y_0
          do kx=1,IDIM
           x=(kx-1)*p_m_x+x_0
           IF ((mark(kx,ky)>=2.AND.yin_yang).OR.
     *         (mark(kx,ky)<2.AND..NOT.yin_yang)) THEN
            IF (i_top==0) THEN
             dm_lay(kx,ky,ins)=-E(kx,ky)*1e-3
            ELSE
             dm_lay(kx,ky,ins)=layers(kx,ky,i_top) 
            ENDIF 
           ENDIF 
          CALL lay_order(kx,ky,ins)
          layers(kx,ky,:)=dm_lay(kx,ky,:)
          write (30,*)x,y,layers(kx,ky,ins)
          enddo
         enddo
!CALL FROM sublitho_an subroutine
       ELSEIF (opt=='U') THEN
        do ky=JDIM,1,-1
          y=(ky-1)*p_m_y+y_0
          do kx=1,IDIM
           x=(kx-1)*p_m_x+x_0
           IF (mark(kx,ky)>=2) THEN
            WRITE(400,'(A1,1X,3I3,1X,3F9.2)')'/',I_z_an,kx,JDIM-ky,
     *      z_sublit_an*1D-3,x,y 
           ENDIF
          enddo
         enddo

       ENDIF
 
      close(30)
      DEALLOCATE (mark)

      END SUBROUTINE
