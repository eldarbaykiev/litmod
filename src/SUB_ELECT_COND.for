!*********** SUB_ELECT_COND.for******************!


      SUBROUTINE ELECT_COND(T_cel,P_Pa)
!      USE M_columns
!      USE M_grid_parameters_1D
      USE M_phases
      USE M_conductivity
      USE M_mant
      USE M_fugacity_ol_solu, ONLY:C_OH,C_OH_cpx,C_OH_opx,C_OH_gt 
      USE M_columns, ONLY:mat_z,I_z
      USE M_SIGMELTS
      REAL(8),ALLOCATABLE:: sg_ph(:,:)
      REAL(8) T_kel,P_GPa,T_cel,P_Pa
       INTEGER i_sg
      
      T_kel=T_cel+273D0 
      P_GPa=P_Pa*1D-9
!      z=-E_max+(I_z-1)*d_z 
      i_sg=0
      C_w=0D0 
      sigma=0D0
      res=0D0
!Olivine
       IF(I_phases(i_LK,i_near,j_near,1)/=0) THEN
!Compute water solubility in olivine
       CALL H20_solubility_ol(T_cel,P_Pa)
!Iron number of olivine
       X_Fe_ol=wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,1))
     */(wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,1))+
     * wt_MgO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,1)))
!Water content
       IF(C_w_onoff(mat_z)==1 ) THEN
        C_w=C_w_ol(I_z)
!      write(*,*)'C_w,I_z, C_OH',C_w,I_z, C_OH
!      pause
!Limit the water content to the solubility of H+ in olivine
        IF(C_w>C_OH) C_w=C_OH
       ENDIF
!!       IF(I_phases(i_LK,i_near,j_near,5)/=0) C_w=C_w/4
!Different conductivity models
       IF(ol_mod==1) THEN
!Formula by Hirsch et al, 1993
       sigma(1)=(10.**6.54/T_kel)*X_Fe_ol**1.81*
     *             EXP(-1.35/(k_bolt*T_kel))
       ELSEIF(ol_mod==2) THEN
!Formula after Yoshino et al 2009
       sigma(1)=10.**4.73*EXP(-2.31/(k_bolt*T_kel))+
     * 10.**2.98*EXP(-1.71/(k_bolt*T_kel))+10.**1.9*C_w*
     * EXP(-(.92-.16*C_w**(1D0/3D0))/(k_bolt*T_kel))
       ELSEIF(ol_mod==3) THEN
!Formula by Vacher and Verhoeven 2007
       sigma(1)=10.**sg_ref_ol_VV*(X_Fe_ol/.1)**al_ol_VV*
     * EXP(-(U_ref_ol_VV+bt_ol_VV*(X_Fe_ol-.1))/(k_bolt*T_kel))   
      ELSEIF(ol_mod==4) THEN
! SO2 (Constable et al. 1992)
       sigma(1)=10.**2.402*EXP(-1.6/(k_bolt*T_kel))+
     *                    10.**9.17*EXP(-4.25/(k_bolt*T_kel))
      ELSEIF(ol_mod==5) THEN
!SEO3 (Constable 2006)
!Oxygen fugacity (Myers and Eugster 1983 )
!QMF (Pa)
       fo2=10**(-2.44419D4/T_kel+13.296)
!IW (Pa)
!       fo2=10**(-2.68347D4/T_kel+11.477)
       sigma(1)=10.**.9952*EXP(-1.407/(k_bolt*T_kel))+
     *                    10.**2.601*EXP(-1.842/(k_bolt*T_kel))+
     *                   (10.**.8135*EXP(-1.07/(k_bolt*T_kel))+
     *                    10.**6.733*EXP(-2.92/(k_bolt*T_kel)))*
     *                    fo2**(1D0/6D0)             
      ELSEIF(ol_mod==6) THEN
!Preferred formula (Based on Wang et al., 2006 Nature))
!Activation enthalpy
       d_H=1.642+.246*X_Fe_ol-4.85*X_Fe_ol**2+3.259*X_Fe_ol**3
       d_H=d_H+(P_GPa*0.68)*k_bolt/R_uni

       sigma(1)=10.**3*EXP(-d_H/(k_bolt*T_kel))+10.**3*
     * C_w**.62*EXP(-.9/(k_bolt*T_kel)) 
     *  +10.**4.73*EXP(-2.31/(k_bolt*T_kel))
!       write(*,*)'ol',d_H
      ELSEIF(ol_mod==7) THEN
!Formula by Du Frane et al., 05
!Oxygen fugacity (Myers and Eugster 1983 )
!QMF (Pa)
       fo2=10**(-2.44419D4/T_kel+8.29)
!IW (Pa)
!       fo2=10**(-2.68347D4/T_kel+6.471)         
       sigma(1)=(2.51D0*fo2**(2D0/11D0)+6.53D-2)
     *                    *EXP(-.531/(k_bolt*T_kel))
       ELSEIF(ol_mod==8) THEN
!Formula after Yoshino et al., 2009 with the water content dependence 
!according to Poe et al., 2010
!      s100=10.**2.57*C_w*EXP(-(1.26-1.18*C_w**(1D0/3D0))/(k_bolt*T_kel))
!      s010=10.**3.46*C_w*EXP(-(1.5-1.43*C_w**(1D0/3D0))/(k_bolt*T_kel))
!      s001=10.**1.02*C_w*EXP(-(.81-.7*C_w**(1D0/3D0))/(k_bolt*T_kel))
!sg_H_Poe is the geometric average of the conductivity along a,b and c
!axes, sg_H_Poe=(s100*s010*s001)**(1/3)
      sg_H_Poe=10.**2.35*C_w*EXP(-(1.19-1.104*C_w**(1D0/3D0))
     *          /(k_bolt*T_kel))
       sigma(1)=10.**4.73*EXP(-2.31/(k_bolt*T_kel))+
     * 10.**2.98*EXP(-1.71/(k_bolt*T_kel))+sg_H_Poe
      ENDIF
       i_sg=i_sg+1
       ELSE
        X_Fe_ol=0D0 
       ENDIF

!Garnet
      IF(I_phases(i_LK,i_near,j_near,4)/=0) THEN       
       X_Fe_gt=wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,4))
     */(wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,4))+
     * wt_MgO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,4)))
!Water content
        IF(C_w_onoff(mat_z)==1) THEN
         C_w=C_w_gt(I_z)
!Limit the water content in garnet according the solubility 
!(Lu and Keppler 1997)
         IF(C_w>C_OH_gt) C_w=C_OH_gt
        ENDIF 
!Jones et al. 2008 Lithos formula
       IF(gt_mod==1) THEN
       sigma(4)=10.**(4.26-12.26*X_Fe_gt)*
     *                    EXP(-(2.4-6*X_Fe_gt)/(k_bolt*T_kel))
!Formula with dependence on T, P but no X_Fe (Dai-Karato_09)
       ELSEIF(gt_mod==2) THEN
       sigma(4)=A_0_Dai_d*(1-B_Dai_d*P_GPa)*
     *                    EXP(-(E_dai_d+P_GPa*V_Dai_d)/(R_uni*T_kel))
     * +A_0_Dai_w*C_w**.63*EXP(-(E_dai_w+P_GPa*V_Dai_w)/(R_uni*T_kel))
       ELSEIF(gt_mod==3) THEN
!Formula by Xu and Shankland JGR 99
       sigma(4)=2239D0*EXP(-1.66/(k_bolt*T_kel))
       ELSEIF(gt_mod==4) THEN 
!Preferred formula 1 (based on Dai and Karato 2009 (PEPI) for small
!polaron and proton conduction, + Xfe from Romano et al., 2006 (Am. Miner))
!Activation energy (eV) at 10 Gpa (Romano et al., 06)
       d_H=2.6-15.33*X_Fe_gt+80.4*X_Fe_gt**2-194.6*X_Fe_gt**3+
     *202.6*X_Fe_gt**4-75*X_Fe_gt**5
!Convert to kJ/mol
       d_H=(d_H-1.55D0)*R_uni/k_bolt
       sigma(4)=A_0_Dai_d*(1-B_Dai_d*P_GPa)*
!!     *  EXP(-(d_H+(P_GPa-10D0)*V_Dai_d)/(R_uni*T_kel))
     *  EXP(-(d_H+E_dai_d+P_GPa*V_Dai_d)/(R_uni*T_kel))
     * +A_0_Dai_w*C_w**.63*EXP(-(E_dai_w+P_GPa*V_Dai_w)/(R_uni*T_kel))
!       sigma(4,I_z-N_z_c)=10.**3*EXP(-d_H/(k_bolt*T_kel))+10.**3.29*
!     * C_w**.63*EXP(-.725/(k_bolt*T_kel))
!       write(*,*)'X_Fe_gt',sigma(4),X_Fe_gt
!       pause
      ELSEIF(gt_mod==5) THEN 
!Preferred formula 2 (based on Yoshino et al. 2008 (PEPI) for small
!polaron and Xfe + proton conduction from Poe et al. 2010
!Activation energy (eV) from Yoshino et al. 2008 
       d_H=1.47-.87*(X_Fe_gt-0.08)
        sg_H_Poe=10.**2.35*C_w*EXP(-(1.19-1.104*C_w**(1D0/3D0))
     *          /(k_bolt*T_kel))
        sigma(4)=522*EXP(-d_H/(k_bolt*T_kel))+sg_H_Poe
      ENDIF
      i_sg=i_sg+1 
!D      WRITE(308,'(F10.3,1X,F9.1)')(E_max-(I_z-1)*d_z)*1D-3,
!D    * C_w*1D4
      ELSE
       X_Fe_gt=0D0
      ENDIF
!Orthopyroxene and C2/c HP clinopyroxene
      IF(I_phases(i_LK,i_near,j_near,2)/=0.OR.
     *   I_phases(i_LK,i_near,j_near,7)/=0) THEN
       IF (I_phases(i_LK,i_near,j_near,2)/=0) THEN
        iopx=2
       ELSE
        iopx=7
       ENDIF
      X_Fe_opx=wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,
     *iopx))/(wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,
     *iopx))+ wt_MgO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,
     *iopx)))
!Water content
       IF(C_w_onoff(mat_z)==1) THEN
        C_w=C_w_opx(I_z)
!Limit the water content to the solubility of H+ in opx
        IF(C_w>C_OH_opx) C_w=C_OH_opx
       ENDIF
       IF(opx_mod==1) THEN 
!Formula by Xu and Shankland 99 JGR
       sigma(2)=10.**3.72*EXP(-1.8/(k_bolt*T_kel))       
       ELSEIF(opx_mod==2) THEN
!Formula by Vacher and Verhoeven 2007
       sigma(2)=10.**sg_ref_opx_VV*(X_Fe_opx/.1)**al_opx_VV*
     * EXP(-(U_ref_opx_VV+bt_opx_VV*(X_Fe_opx-.1))/(k_bolt*T_kel))
       ELSEIF(opx_mod==3) THEN
!Preferred formula 1 based on Dai and Karato 09 P Jpn A SB for small
! polaron and proton conduction, + from XFe Seifert et al.,82
       d_H=1.9-2.77*X_Fe_opx+2.61*X_Fe_opx**2-1.09*X_Fe_opx**3
       sigma(2)=10.**3* EXP(-d_H/(k_bolt*T_kel))+
     *10.**2.6*C_w**.62*EXP(-.85/(k_bolt*T_kel))
       ELSEIF(opx_mod==4) THEN
!Preferred formula 2 based on Dai and Karato 09 P Jpn A SB for small
! polaron + proton conduction from Poe et al. 2010 + XFe from Seifert et al.,82
        d_H=1.9-2.77*X_Fe_opx+2.61*X_Fe_opx**2-1.09*X_Fe_opx**3
        sg_H_Poe=10.**2.35*C_w*EXP(-(1.19-1.104*C_w**(1D0/3D0))
     *          /(k_bolt*T_kel))
        sigma(2)=10.**3*EXP(-d_H/(k_bolt*T_kel))+sg_H_Poe
       ENDIF 
       i_sg=i_sg+1 
!       write(*,*)'opx',d_H 
!D       WRITE(307,'(F10.3,1X,F9.1)')(E_max-(I_z-1)*d_z)*1D-3,
!D    * C_w*1D4
      ENDIF 
!Clinopyroxene
      IF(I_phases(i_LK,i_near,j_near,3)/=0) THEN
      X_Fe_cpx=wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,3))
     */(wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,3))+
     * wt_MgO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,3)))
!Water content
       IF(C_w_onoff(mat_z)==1) THEN
        C_w=C_w_cpx(I_z)
!Limit the water content to the solubility of H+ in opx
        IF(C_w>C_OH_cpx) C_w=C_OH_cpx
       ENDIF
        IF(opx_mod==1) THEN 
! Xu and Shankland 99 JGR
         sigma(3)=10.**3.25*EXP(-1.87/(k_bolt*T_kel))
        ELSEIF(opx_mod==2) THEN 
!Formula by Vacher and Verhoeven 2007
        sigma(3)=10.**sg_ref_cpx_VV*(X_Fe_cpx/.1)**al_cpx_VV*
     *  EXP(-(U_ref_cpx_VV+bt_cpx_VV*(X_Fe_cpx-.1))/(k_bolt*T_kel))
        ELSEIF(opx_mod==3) THEN
!Preferred formula 1 based on Xu and Shankland 99 JGR for small
! polaron and Dai and Karato 09 P Jpn A SB for proton conduction, 
!+ XFe from Seifert et al.,82
       d_H=2.075-2.77*X_Fe_cpx+2.61*X_Fe_cpx**2-1.09*X_Fe_cpx**3
       sigma(3)=10.**3.25*EXP(-d_H/(k_bolt*T_kel))+
     *10.**2.6*C_w**.62*EXP(-.85/(k_bolt*T_kel))
       ELSEIF(opx_mod==4) THEN
!Preferred formula 2 based on Xu and Shankland 99 JGR for small
! polaron + proton conduction from Poe et al. 2010 + XFe from Seifert et al.,82
       d_H=2.075-2.77*X_Fe_cpx+2.61*X_Fe_cpx**2-1.09*X_Fe_cpx**3
       sg_H_Poe=10.**2.35*C_w*EXP(-(1.19-1.104*C_w**(1D0/3D0))
     *          /(k_bolt*T_kel))
       sigma(3)=10.**3.25*EXP(-d_H/(k_bolt*T_kel))+sg_H_Poe
       ENDIF

       i_sg=i_sg+1
!D       WRITE(309,'(F10.3,1X,F9.1)')(E_max-(I_z-1)*d_z)*1D-3,
!D    * C_w*1D4
      ENDIF

!Spinel
      IF(I_phases(i_LK,i_near,j_near,5)/=0) THEN
!      X_Fe_sp=wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,5))
!     */(wt_FeO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,5))+
!     * wt_MgO(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,5)))
!Preferred formula from (Nell and Wood 89 and 91, Am. Min.)
       d_H=.385
       sigma(5)=EXP(11.635D0)*EXP(-d_H/(k_bolt*T_kel))/T_kel
!       write(*,*)'sp',sigma(5,I_z-N_z_c) 
       i_sg=i_sg+1
      ENDIF
!      write(*,*)sigma(:),X_Fe_ol,X_Fe_gt,X_Fe_opx,X_Fe_cpx
!      pause
!Partial melt
      IF(melt(I_z)>0D0.AND.melt_onoff(mat_z)==1) THEN
       IF(mlt_mod==1) THEN
!Carbonatitic melt
!Yoshino et al., 2010 EPSL
       sigma_mlt=955.3*EXP(-.38/(k_bolt*T_kel))
      ELSEIF(mlt_mod==2) THEN
!Gaillard et al., 2008:
       sigma_mlt=3442*EXP(-.33/(k_bolt*T_kel))   
       ELSEIF(mlt_mod==3) THEN
!Rhyolitic melt (Gaillard 2004 EPSL) 
       sigma_mlt=850*EXP(-(70.5D0+P_GPa*5.4D0)/(R_uni*T_kel))
       ELSEIF(mlt_mod==4) THEN
!Silicate (basaltic and rhyolitic) melt (Pommier et al, 2008, Pommier and LeTrong 2011)
!        CALL SIGMELTS (T_kel,P_GPa)
!        sigma_mlt=s_0_pom
          WRITE(*,*)'SIGMELTS subroutine not avialable in this version'
          WRITE(*,*)'Aborting...'
          STOP
       ENDIF
!      write(*,*)'mlt_mod,sigma_mlt',mlt_mod,sigma_mlt
!!      pause
      ENDIF
!D      WRITE(300,'(F10.3,4F8.3)')(E_max-(I_z-1)*d_z)*1D-3,
!D    *                           X_Fe_ol,X_Fe_gt,X_Fe_cpx,X_Fe_opx

      IF (ALLOCATED(sg_ph)) DEALLOCATE(sg_ph)
      ALLOCATE(sg_ph(i_sg,2))
      
      i_sg=0
       IF(I_phases(i_LK,i_near,j_near,1)/=0) THEN
        i_sg=i_sg+1
        sg_ph(i_sg,1)=sigma(1)
        sg_ph(i_sg,2)=
     *  vol_frac(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,1))
!D       WRITE(301,'(F10.3,1X,F7.2,E12.3)')(E_max-(I_z-1)*d_z)*1D-3,
!D    *  sg_ph(i_sg,2),1/sg_ph(i_sg,1)
       ENDIF 
       IF(I_phases(i_LK,i_near,j_near,4)/=0) THEN
        i_sg=i_sg+1
        sg_ph(i_sg,1)=sigma(4)
        sg_ph(i_sg,2)=
     *  vol_frac(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,4))
!D       WRITE(302,'(F10.3,1X,F7.2,E12.3)')(E_max-(I_z-1)*d_z)*1D-3,
!D    *  sg_ph(i_sg,2),1/sg_ph(i_sg,1)
       ENDIF 
       IF(I_phases(i_LK,i_near,j_near,2)/=0.OR.
     *    I_phases(i_LK,i_near,j_near,7)/=0) THEN
        i_sg=i_sg+1
        sg_ph(i_sg,1)=sigma(2)
        IF (I_phases(i_LK,i_near,j_near,2)/=0) THEN
         iopx=2
        ELSE
         iopx=7
        ENDIF
        sg_ph(i_sg,2)=
     *  vol_frac(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,iopx))
!D       WRITE(303,'(F10.3,1X,F7.2,E12.3)')(E_max-(I_z-1)*d_z)*1D-3,
!D    *  sg_ph(i_sg,2),1/sg_ph(i_sg,1)
       ENDIF
       IF(I_phases(i_LK,i_near,j_near,3)/=0) THEN
        i_sg=i_sg+1 
        sg_ph(i_sg,1)=sigma(3)
        sg_ph(i_sg,2)=
     *  vol_frac(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,3))
!D       WRITE(304,'(F10.3,1X,F7.2,E12.3)')(E_max-(I_z-1)*d_z)*1D-3,
!D    *  sg_ph(i_sg,2),1/sg_ph(i_sg,1)
       ENDIF
       IF(I_phases(i_LK,i_near,j_near,5)/=0) THEN
        i_sg=i_sg+1
        sg_ph(i_sg,1)=sigma(5)
        sg_ph(i_sg,2)=
     *  vol_frac(i_LK,i_near,j_near,I_phases(i_LK,i_near,j_near,5))
!D       WRITE(305,'(F10.3,1X,F7.2,E12.3)')(E_max-(I_z-1)*d_z)*1D-3,
!D    *sg_ph(i_sg,2),1/sg_ph(i_sg,1) 
       ENDIF
!Averging schemes
!Hashin-Shtrikman extremal bounds
!      write(*,*) i_sg,sg_ph(:,1)
!      pause
      sig_min=MINVAL(sg_ph(:,1))
      sig_max=MAXVAL(sg_ph(:,1))
      sum_min=0D0
      sum_max=0D0
      sum_s=0D0
      sum_p=0D0 
      DO i=1,i_sg
!LK   
       sum_min=sum_min+1D-2*sg_ph(i,2)/(sg_ph(i,1)+2*sig_min)
       sum_max=sum_max+1D-2*sg_ph(i,2)/(sg_ph(i,1)+2*sig_max)
       sum_s=sum_s+1D-2*sg_ph(i,2)*sg_ph(i,1)
       sum_p=sum_p+1D-2*sg_ph(i,2)/sg_ph(i,1)
      ENDDO
      sum_p=1D0/sum_p
!sigma(1,:)-->HS min conductivity (max resistivity)
!sigma(2,:)-->HS max conductivity (min resistivity)
!sigma(3,:)-->parellel resistivity (min resistivity)
!sigma(4,:)-->series resistivity (max resistivity)
      sigma(1)=1/sum_min-2*sig_min
      sigma(2)=1/sum_max-2*sig_max
      sigma(3)=sum_s
      sigma(4)=sum_p
!      z_lay_sg(I_z-N_z_c+1)=-E_max+(I_z-1)*d_z
!Add partial melt to the electrical conductivity
!Archie's law
       IF(mlt_mod==1) THEN
        C_arch=.97D0
        n_arch=1.14D0
       ELSE
        C_arch=.67D0
        n_arch=.89D0
       ENDIF
       IF(melt(I_z)>0D0.AND.melt_onoff(mat_z)==1) THEN
!D       WRITE(311,'(F10.3,1X,E12.3)')(E_max-(I_z-1)*d_z)*1D-3,
!D    *  1/sigma_mlt(I_z-N_z_c)
        IF(melt_mix==1) THEN
         sigma(2)=C_arch*melt(I_z)**n_arch*sigma_mlt
        ELSEIF(melt_mix==2) THEN           
!Cubic grains (Waff, 1974)
       sigma(2)=(1-(1-melt(I_z))**(2D0/3D0))*sigma_mlt 
!HS max
        ELSEIF(melt_mix==3) THEN 
       sigma(2)=sigma_mlt+(1-melt(I_z))/
     *  (1/(sigma(1)-sigma_mlt)+melt(I_z)/
     *  (3*sigma_mlt))
        ENDIF
!      write(*,*)I_z,sigma_mlt,melt(I_z),P_GPa
!      pause
       ENDIF
       res(:)=1/sigma(1:4)
!D       WRITE(310,'(F10.3,1X,F11.4)')(E_max-(I_z-1)*d_z)*1D-3,
!D    * melt(I_z)   
!Electrical asthenosphere
!      dz_elas=50D3
!      rho_elas=10D0
!      IF ((I_z-1)*d_z-z_lay(N_lay)<dz_elas*1D3.AND.
!     *    (I_z-1)*d_z>z_lay(N_lay)) THEN
!       sigma(:,I_z-N_z_c)=1D0/rho_elas
!      ENDIF
!      write(*,'(F7.2,1X,F7.2,1X,2F7.3)')(E_max-(I_z-1)*d_z)
!     *      *1d-3,T(I_z),LOG10(1/sigma(5:6,I_z-N_z_c))
      END SUBROUTINE   
