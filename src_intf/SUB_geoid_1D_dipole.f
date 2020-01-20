! SUB_geoid_1D_dipole.f

      SUBROUTINE geoid_1D_dipole(z_lims,n,d_N,E_2D,T_mh)
      USE M_profile
      USE M_material
      USE M_ctes
      REAL(4) z_lims(0:N_lay),d_N,rho_av,rxz,E_2D,T_mh
      REAL(8) z_c,z_l,E_pt
      INTEGER n
!Z increases and is positive downwards 
      z_c=z_lims(N_lay-1)
      z_l=z_lims(N_lay)
      d_N=0.
      rxz=0.
      DO i=0,n-1
       IF (z_lims(i+1)<=z_lims(i)) CYCLE 
       IF (i==0.AND.z_lims(0)>0) THEN
        rho_av=prop(i+1)%rho+prop(i+1)%beta*(-z_lims(0)
     *         +(z_lims(i)+z_lims(i+1))/2) 
        d_N=rho_av*(z_lims(i+1)**2-z_lims(i)**2)+rho_w*z_lims(0)**2
        rxz=rho_av*(z_lims(i+1)-z_lims(i))
       ELSEIF (i==N_lay-1) THEN
        E_pt=-z_lims(0) 
        IF (prop(i+1)%rho>0) THEN
         d_N=d_N+prop(i+1)%rho*(z_l**2-z_c**2)
         rxz=rxz+prop(i+1)%rho*(z_l-z_c)
        ELSE
!!         CALL geoterm_1D(z_c,z_l,E_pt,T_mh)  
         d_N=d_N+rho_a*(z_bot**2-z_c**2+prop(N_lay)%alpha*(T_a-T_mh)*
     *       ((z_l-z_c)*(z_l+2*z_c))/3)
         rxz=rxz+rho_a*(1+prop(N_lay)%alpha*(T_a-T_mh)/2)*(z_l-z_c)
        ENDIF
       ELSE
        rho_av=prop(i+1)%rho+prop(i+1)%beta*(-z_lims(0)
     *         +(z_lims(i)+z_lims(i+1))/2)
        d_N=d_N+rho_av*(z_lims(i+1)**2-z_lims(i)**2)
        rxz=rxz+rho_av*(z_lims(i+1)-z_lims(i))
!        write(*,*)'crt',i,z_lims(i+1)*1e-3,z_lims(i)*1e-3,rho_av,
!     * rho_av*(z_lims(i+1)-z_lims(i))
       ENDIF
      ENDDO
      d_N=-d_N*pi*G_u/g
      E_2D=z_l+E_pt-2320-rxz/rho_a 
!      write(*,*)'z_l,E_pt,rxz',z_l,E_pt,rxz
      IF (E_pt<0) E_2D=E_2D*rho_a/(rho_a-rho_w)

      END SUBROUTINE  
