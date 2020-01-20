! two2three_lay.f 



      SUBROUTINE two2three_lay
      use M_imag
      use M_profile
      use M_material 
      real(4) r_c,x,y,r1,r2,r3,r4,E_d,E_u,rho_low,rho_up,z_c,rho_m
      real(4) dE,numt,den,E_fil

      IF (N_lay>2) RETURN
!      rho_up=2700
!      rho_low=2920
!      rho_m=prop(2)%rho
!      rho_m=3245
      rho_c=prop(1)%rho
      call pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               -3) 
      CALL SYSTEM('mv ./layers/layer2.xyz ./layers/layer4.xyz')
      prop(4)=prop(2)

      prop(1)%name_mat='oceanic_crust' 
      prop(1)%K=2.1
      prop(1)%A=2E-7
      prop(1)%alpha=0E0
      prop(1)%beta=0E0

      prop(2)%name_mat='upper_crust'
      prop(2)%rho=2700
      prop(2)%K=2.5
      prop(2)%A=1e-6
      prop(2)%alpha=0E0
      prop(2)%beta=0E0

      prop(3)%name_mat='lower_crust'
      prop(3)%rho=2920
      prop(3)%K=2.1
      prop(3)%A=2E-7
      prop(3)%alpha=0E0
      prop(3)%beta=0E0       

      E_d=-2000
      E_u=-200  
      dE=E_u-E_d

      dm_lay=layers
      IF(ALLOCATED(layers)) DEALLOCATE(layers)
      ALLOCATE (layers(IDIM,JDIM,N_lay+2))
      layers(:,:,1:N_lay)=dm_lay
      DEALLOCATE(dm_lay,GLUE_m)
      ALLOCATE (dm_lay(IDIM,JDIM,N_lay+2))
      ALLOCATE (GLUE_m(IDIM,JDIM,N_lay+2,N_lay+2)) 

      layers(:,:,4)=layers(:,:,2)
      N_lay=N_lay+2
!Rewrite layers.info with the new layer
      open(15,file='layers.info',status='OLD')
      DO i=0,N_lay
       IF (i==0) THEN
       write(15,2)"'BODY NUMBER | DENSITY | ALPHA | THERMAL CONDUCTIVITY
     * | HEAT PRODUCTION RATE | PRESSURE COEFFICIENT'"
       ELSE
       write(15,3)i,prop(i)%name_mat,prop(i)%rho,prop(i)%alpha,prop(i)%K
     *            ,prop(i)%A,prop(i)%beta
       ENDIF
      ENDDO
  2   FORMAT (A100)
  3   FORMAT (I2,'  ',A20,F6.0,E12.3,F5.1,E10.2,E10.2)
      close (15)

      CALL pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               -8)
      CALL PGSLCT(ISTAT)
      CALL clear_scr  (0)
      CALL get_pol
      CALL in_out_pol ('T') 

!      open(1,file='./layers/layer1.xyz',status='UNKNOWN')
!      open(2,file='./layers/layer2.xyz',status='UNKNOWN')
!      open(3,file='./layers/rho_m_TER.xyz',status='OLD') 
! Filtered ELEVATION
!      open(20,file='elev_fil.xyz',status='OLD')
!      DO j=JDIM,1,-1
!       y=(j-1)*p_m_y+lat_min
!        DO i=1,IDIM
!         x=(i-1)*p_m_x+lon_min
!         z_c=layers(i,j,1)
!         IF (E(i,j)<E_d) THEN
!         layers(i,j,1)=-E(i,j)*1e-3
!         layers(i,j,2)=(z_c*(rho_m-rho_c)+E(i,j)*1e-3*(rho_low-rho_c))/
!     *                 (rho_m-rho_low)
!!         ELSEIF (E(i,j)<=E_u) THEN
!!         numt=(z_c+E(i,j)*1e-3)*rho_c-z_c*rho_m-E(i,j)*1e-3*rho_low
!!         den=((E(i,j)-E_d)*rho_up+(2*E_u-E_d-E(i,j))*rho_low
!!     *       -2*dE*rho_m)/(2*dE)
!!         layers(i,j,2)=numt/den
!!         layers(i,j,1)=layers(i,j,2)*(E(i,j)-E_d)/(2*dE)-E(i,j)*1e-3
         
!         r_c=prop(2)%rho-((prop(2)%rho-2820)/(abs(E_d-E_u)))*
!     *       (E(i,j)-E_d)
!         layers(i,j,1)=(layers(i,j,2)*(r_c-prop(2)%rho)+E(i,j)*1e-3*
!     *                 (r_c-prop(1)%rho))/(prop(1)%rho-prop(2)%rho)
!         ELSE
!!         numt=(z_c+E(i,j)*1e-3)*rho_c-z_c*rho_m-E(i,j)*1e-3*rho_up 
!!         den=(rho_up+rho_low-2*rho_m)/2
!!         layers(i,j,2)=numt/den
!!         layers(i,j,1)=layers(i,j,2)/2 
!          layers(i,j,2)=z_c
!          layers(i,j,1)=((E(i,j)*1e-3+z_c)*rho_c-E(i,j)*1e-3*rho_up-
!     *                  z_c*rho_low)/(rho_up-rho_low)
!         layers(i,j,1)=(layers(i,j,2)*(2820-prop(2)%rho)+E(i,j)*1e-3*
!     *                 (2820-prop(1)%rho))/(prop(1)%rho-prop(2)%rho)
!         ENDIF
!         write (1,*)x,y,layers(i,j,1)
!         write (2,*)x,y,layers(i,j,2) 
!        ENDDO
!      ENDDO
!      close(1)
!      close(2)
!      close(3)
!      close(20)
      
      END SUBROUTINE
