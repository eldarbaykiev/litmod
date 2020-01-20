! Subroutine  SUB_geoterm_1D.f


      SUBROUTINE geoterm_1D(z_c,z_l,E,T_mh)
      USE M_ctes
      REAL(4) theta,k_c,T_s,H_s,h_r,delta,k_m,ter1
      REAL(4) z_c,z_l,T_mh,E
      k_m=3.2
      h_r=15D3
      T_s=0
      H_s=2.50D-6
      k_c=2.5       
!      write(*,*)'k_c,T_s,H_s,z_c,'
!      z_c=-z_c
!      z_l=-z_l
      theta=k_c*T_s+H_s*(1-exp(-(z_c+E)/h_r))*h_r**2
      delta=k_m*T_a*(z_c+E)
      ter1=z_c*(k_m-k_c)+k_c*z_l+E*k_m
      T_mh=(theta*(z_l-z_c)+delta)/ter1 
            
 
      END SUBROUTINE    
