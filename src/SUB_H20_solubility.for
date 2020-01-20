       SUBROUTINE H20_solubility_ol(T_cel,P_Pa)
       USE M_fugacity_ol_solu
       USE M_conductivity, ONLY: R_uni,X_Fe_ol
       REAL(8) T_kel,P_GPa,T_cel,P_Pa
       T_kel=T_cel+273D0
       P_GPa=P_Pa*1D-9 
!Initial guess for the molar volume (cm3/mol)...?
      rho_mol=1D0/5
      T_vec_C(1)=1/T_kel**4D0
      T_vec_C(2)=1/T_kel**2D0
      T_vec_C(3)=1/T_kel
      T_vec_C(4)=1
      T_vec_C(5)=T_kel
      T_vec_C(6)=T_kel**2D0 
!Initial molar density (in mol/cm3)
      cte=P_GPa/(R_uni*T_kel) 
!Compute the coefficients     
      c_EoS(:)=0
      DO i=1,10
       DO j=1,6
        c_EoS(i)=c_EoS(i)+c_T(i,j)*T_vec_C(j)
       ENDDO
      ENDDO
      
      del_rho=rho_mol/10     
      DO
       rho_0=rho_mol    
       CALL molar_dens
       F0=F

       rho_mol=rho_0+del_rho
       CALL molar_dens
       Ffwd=F
 
       rho_mol=rho_0-del_rho
       CALL molar_dens
       Fbwd=F

!Check if the function crosses zero in the interval
       IF ((Ffwd>0.AND.Fbwd>0).OR.(Ffwd<0.AND.Fbwd<0)) THEN       
        IF(F0>0) THEN
         IF(Ffwd>Fbwd) THEN
          rho_mol=rho_0-del_rho
         ELSE
          rho_mol=rho_0+del_rho
         ENDIF         
  
        ELSE
         IF(Ffwd>Fbwd) THEN
          rho_mol=rho_0+del_rho
         ELSE
          rho_mol=rho_0-del_rho
         ENDIF
        ENDIF
       ELSE
!Crossing  zero
        del_rho=del_rho/2
        IF(F0>0) THEN
         IF(Ffwd>0) THEN
          rho_mol=rho_0-del_rho
         ELSE
          rho_mol=rho_0+del_rho
         ENDIF 
        ELSE
         IF(Ffwd>0) THEN
          rho_mol=rho_0+del_rho
         ELSE
          rho_mol=rho_0-del_rho
         ENDIF
        ENDIF
       ENDIF

!       WRITE(*,*)'RHO_0,RHO,del_rho',rho_0-rho,del_rho 
!       WRITE(*,*)'molar volumen (cm3/mol)',1/rho
!       pause
      IF(rho_mol<1D-8.OR.rho_mol>1D8) THEN
       WRITE(*,*) 'WATER FUGACITY EoS'
       WRITE(*,*)'CONVERGENCE NOT ACHIEVED, TRY OTHER INITIAL MOLAR 
     *VOLUME!'
       write(*,*)'T(K), P(GPa)',T_kel,P_GPa 
       STOP
      ENDIF   
      IF(ABS(rho_mol-rho_0)<.00001) EXIT

      ENDDO
!      WRITE(*,*)'molar volumen (cm3/mol)',1/rho_mol

!Calculate water fugacity
      Ares_adim=c_EoS(1)*rho_mol+(1/(c_EoS(2)+c_EoS(3)*rho_mol+
     *c_EoS(4)*rho_mol**2D0+c_EoS(5)*rho_mol**3D0+
     *c_EoS(6)*rho_mol**4D0)-1/c_EoS(2))-(c_EoS(7)/c_EoS(8))*(exp
     *(-c_EoS(8)*rho_mol)-1)-(c_EoS(9)/c_EoS(10))*(exp
     *(-c_EoS(10)*rho_mol)-1)
       
      fH2O=rho_mol*R_uni*T_kel*exp(Ares_adim+cte/rho_mol-1)
!      WRITE(*,*)'water fugacity (GPa)',fH2O 

!Calculate water solubility in olivine according to Zhao et al, 2004
!      WRITE(*,*)'Olivine iron content X_Fe_ol...?'
!      READ(*,*)X_Fe_ol 
      C_OH=Azhao*fH2O*exp((-(E_zhao+P_GPa*dV_zhao)+alpha_zhao*X_Fe_ol)/
     *    (R_uni*T_kel))
!      WRITE(*,*)'H+ solubility in olivine (wt ppm): ',C_OH 
      C_OH=C_OH*1D-4
!Calculate water solubility in opx according to Mierdel et al., 2007
      C_OH_opx_0=Amie*fH2O*exp(-(E_mie+P_GPa*dV_mie)/(R_uni*T_kel))
      C_OH_opx_Al=Amie_al*SQRT(fH2O)*exp(-(E_mie_al+P_GPa*dV_mie_al)
     +         /(R_uni*T_kel))
      C_OH_opx=C_OH_opx_0+C_OH_opx_Al
!      WRITE(*,*)'H+ solubility in opx (wt ppm): ',C_OH_opx_0,
!     *           C_OH_opx_Al,C_OH_opx
       C_OH_opx=C_OH_opx*1D-4
!Calculate water solubility in cpx according to Gavrilenko 2008 PhD
      C_OH_cpx_0=Agav*fH2O*exp(-(E_gav+P_GPa*dV_gav)/(R_uni*T_kel))
      C_OH_cpx=C_OH_cpx_0+C_OH_opx_Al
!      WRITE(*,*)'H+ solubility in cpx (wt ppm): ',C_OH_cpx_0,
!     *           C_OH_opx_Al,C_OH_cpx      
      C_OH_cpx=C_OH_cpx*1D-4     
!Calculate water solubility in garnet according to Lu and Keppler,1997
      C_OH_gt=Alu*SQRT(fH2O)*exp(-(E_lu+P_GPa*dV_lu)/
     *    (R_uni*T_kel))
!      WRITE(*,*)'H+ solubility in olivine (wt ppm): ',C_OH 
!      WRITE(*,*)C_OH_gt,Alu*SQRT(fH2O)*exp(-(P_GPa*dV_lu)/(R_uni*T_kel))
!      pause
      C_OH_gt=C_OH_gt*1D-4
      END

      SUBROUTINE molar_dens
      USE M_fugacity_ol_solu
      F=rho_0+(rho_mol**2D0*(c_EoS(1)-((c_EoS(3)+2*c_EoS(4)*rho_mol+
     *  3*c_EoS(5)*rho_mol**2D0+4*c_EoS(6)*rho_mol**3D0)/(c_EoS(2)+
     *  c_EoS(3)*rho_mol+c_EoS(4)*rho_mol**2D0+
     *  c_EoS(5)*rho_mol**3D0+c_EoS(6)*rho_mol**4D0)**2))+
     *  c_EoS(7)*rho_mol**2D0*exp(-c_EoS(8)*rho_mol)+
     *  c_EoS(9)*rho_mol**2D0*exp(-c_EoS(10)*rho_mol))-cte  
      END SUBROUTINE 
