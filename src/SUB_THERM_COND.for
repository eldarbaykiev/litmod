!Subroutine THERM_COND

      SUBROUTINE THERM_COND 
      USE M_temperatures
      USE M_grid_parameters
      USE M_columns
      USE M_material
      REAL(8) k_rad,k_pres,k_exp,T_kel
      REAL(8) a,b,c,d

      a=2.85D-5
      b=1.008D-8
      c=-.384
      d=0D0  
      DO I_x=1,N_x
       DO I_y=1,N_y
        DO I_z=1,N_z
        z=-(E_max-(I_z-1)*d_z)
        T_kel=T(I_x,I_y,I_z)+273D0
        IF (mat(I_x,I_y,I_z)<=0) CYCLE
        IF (prop(mat(I_x,I_y,I_z))%gmma_T==0) THEN
         k_T_P(I_x,I_y,I_z)=prop(mat(I_x,I_y,I_z))%K
         CYCLE          
        ENDIF

!!   The T and P dependence follows the model of Hofmeister (1999)
!!   Refer to Afonso et al. G3 (2008)

!!   Radiative contribution (polynomial)

        k_rad=.0175-(1.037D-4*T_kel)+(2.245D0*T_kel**2/1D7)-
     *       (3.407*T_kel**3/1D11)

!!  Pressure effect
        k_pres=1+(prop(mat(I_x,I_y,I_z))%K_0p*g_av*3350D0*z*1D-9)/
     *            prop(mat(I_x,I_y,I_z))%K_T 
!!  Exponential factor in main formula

        k_exp=-(a*(T_kel-298D0)+b*.5*(T_kel**2-88804D0)+
     *        c*(3.3557D-3-1/T_kel)+d*(T_kel**5-2.35D12)/5)*
     *        (prop(mat(I_x,I_y,I_z))%gmma_T*4+1/3) 

!!   Final form of k(T,P)...

        k_T_P(I_x,I_y,I_z)=(prop(mat(I_x,I_y,I_z))%K*(298D0/T_kel)**.45)
     *  *EXP(k_exp)*k_pres+k_rad
        ENDDO
       ENDDO
      ENDDO 
      



      END SUBROUTINE
