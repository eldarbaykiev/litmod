!Subroutine  SUB_cross_sect.f


      SUBROUTINE cross_sect
      USE M_imag
      USE M_profile
      USE M_material
      USE M_ctes
      INTEGER::PGOPEN,ICS
      REAL(4) L_prf,d_x,FA_av,GEOID_av,a_obs,b_obs,a,b,N_0
      REAL(4) d_N,z_lims(0:N_lay)
      REAL(8)DELTA,XB(2),YB(2),ZB(2),ZBG(2),RHO_B,rho_red
      REAL(8) rho_up,rho_down,T_geo,T_mh,XP,YP,GR,GEO,BG,z_up,z_down
      REAL(8) E_ref,z_cref,z_lref,rho_cref,rho_mref,E_pt,rho_ref 
      CHARACTER(1) chdm 
      
      IF(ALLOCATED(FA_2D)) DEALLOCATE(FA_2D) 
      IF(ALLOCATED(GEOID_2D)) DEALLOCATE(GEOID_2D)
      IF(ALLOCATED(BGA_2D)) DEALLOCATE(BGA_2D)
      IF(ALLOCATED(GEOID_1D)) DEALLOCATE(GEOID_1D)
      IF(ALLOCATED(E_2D_v)) DEALLOCATE(E_2D_v)
      ALLOCATE (FA_2D(I_pr_end),GEOID_2D(I_pr_end),BGA_2D(I_pr_end),
     *          GEOID_1D(I_pr_end),E_2D_v(I_pr_end))
      rho_red=2670.

      FA_2D=0.
      GEOID_2D=0.
      BGA_2D=0.
 
!!      IF(LINUX) THEN
!!         ICS=PGOPEN('/XWINDOW')
!!      ELSE
!!         ICS=PGOPEN('/WX')
!!      ENDIF
!!      IF (ICS.LE.0) STOP
!!      CALL pgscr(0,1.,1.,1.)
!!      CALL pgscr(1,0.,0.,0.)
!!      CALL pgask(.false.)      
 
       H_dist=H_dist*1E3
       bod_pr=-bod_pr*1E3      
       L_prf=MAXVAL(H_dist(:))
       d_x=L_prf/(I_pr_end-1)
       DELTA=1D9
!Transversal limits (m)
       YB(1)=-1D6
       YB(2)=1D6

       write(*,*)'START 2D CALCULATION...'
 
      DO
!      write(*,*)'DELTA,YB(1),YB(2) (m)'
!      read(*,*)DELTA,YB(1),YB(2) 
      FA_2D=0.
      GEOID_2D=0.
      BGA_2D=0.
!!            bod_pr(:,0)=0. !QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
!!            bod_pr(:,1)=-31.4E3
!!            bod_pr(:,2)= -151.2E3
!Loop for columns
       DO i=1,I_pr_end
        XB(1)=H_dist(i)-d_x*.5D0
        XB(2)=H_dist(i)+d_x*.5D0
!Extend the borders
         IF (i==1) THEN
          XB(1)=XB(1)-DELTA
         ELSEIF (i==I_pr_end) THEN
          XB(2)=XB(2)+DELTA
         ENDIF
!Loop for bodies
        DO k=0,N_lay
         IF (k==0.AND.bod_pr(i,k)<0D0) THEN
         rho_B=(rho_red-rho_w)
         ZBG(2)=bod_pr(i,k) 
         ZBG(1)=0D0
         ELSEIF (k==0.AND.bod_pr(i,k)>=0D0) THEN
         rho_B=-rho_red
         ZBG(2)=0D0
         ZBG(1)=bod_pr(i,k)
         ENDIF
!Densities and vertical limits of the prisms
         IF (k==N_lay-1) THEN
          ZB(1)=bod_pr(i,k)
          ZB(2)=bod_pr(i,k+1)
          IF (ZB(2)>=ZB(1)) CYCLE
          E_pt=bod_pr(i,0)
!!          CALL geoterm_1D(-ZB(1),-ZB(2),E_pt,T_mh)
          
          IF (prop(k+1)%rho>0) THEN
           RHO_UP=prop(k+1)%rho*(1+prop(k+1)%beta*(bod_pr(i,0)-ZB(1)))
           RHO_DOWN=prop(k+1)%rho*(1+prop(k+1)%beta*(bod_pr(i,0)-ZB(2)))
          ELSE
           RHO_UP=-prop(k+1)%rho*(1+prop(k+1)%alpha*(T_a-T_mh_prf(i)))+
     *           prop(k+1)%beta*(bod_pr(i,0)-ZB(1))
           RHO_DOWN=rho_a
          ENDIF

         ELSEIF (k==N_lay) THEN
          ZB(1)=bod_pr(i,N_lay)
          ZB(2)=z_bot
          RHO_UP=rho_a
          RHO_DOWN=rho_a
         ELSE
          ZB(1)=bod_pr(i,k)
          ZB(2)=bod_pr(i,k+1)
          T_geo=0.
          RHO_UP=prop(k+1)%rho*(1-prop(k+1)%alpha*T_geo)+
     *          prop(k+1)%beta*(bod_pr(i,0)-ZB(1))
          RHO_DOWN=prop(k+1)%rho*(1-prop(k+1)%alpha*T_geo)+
     *            prop(k+1)%beta*(bod_pr(i,0)-ZB(2))

         ENDIF
!Loop for calculation points
!---------------------------------------
        YP=0D0
        DO j=1,I_pr_end
         XP=H_dist(j)
         E_pt=bod_pr(j,0)
!Onshore calc point
         IF (E_pt>0) THEN
          CALL GRAV_GRAD3D(XB,YB,ZB,RHO_up,RHO_down,XP,YP,
     *                     E_pt,GR)
          FA_2D(j)=FA_2D(j)+GR
          CALL GEO_GRAD3D(XB,YB,ZB-E_pt,RHO_up,RHO_down,XP,YP,0D0,
     *                    GEO)
          GEOID_2D(j)=GEOID_2D(j)+GEO 

          IF (k==0) THEN
!Water layer
           IF (ZB(1)<0) THEN
            CALL GRAV_GRAD3D(XB,YB,ZBG,rho_w,rho_w,XP,YP,
     *                       E_pt,GR)
            FA_2D(j)=FA_2D(j)+GR
            CALL GEO_GRAD3D(XB,YB,ZBG-E_pt,rho_w,rho_w,XP,YP,0D0,
     *                      GEO)
            GEOID_2D(j)=GEOID_2D(j)+GEO
           ENDIF
           CALL GRAV_GRAD3D(XB,YB,ZBG,RHO_B,RHO_B,XP,
     *                      YP,E_pt,BG)
           BGA_2D(j)=BGA_2D(j)+BG
          ENDIF
!Offshore calc point 
         ELSE
!          write(*,*)'j,ZB,RHO_up,RHO_down',
!     +    j,ZB(1),ZB(2),RHO_up,RHO_down   
!          pause
          CALL GEOGRAV_GRAD3D(XB,YB,ZB,RHO_up,RHO_down,XP,YP,0D0,
     *                        GEO,GR)
          FA_2D(j)=FA_2D(j)+GR
          GEOID_2D(j)=GEOID_2D(j)+GEO
           IF (k==0) THEN
!Water layer&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            IF (ZB(1)<0) THEN
             CALL GEOGRAV_GRAD3D(XB,YB,ZBG,rho_w,rho_w,XP,YP,0D0,
     *                           GEO,GR)
             FA_2D(j)=FA_2D(j)+GR
             GEOID_2D(j)=GEOID_2D(j)+GEO
            ENDIF
!Water layer&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            CALL GRAV_GRAD3D(XB,YB,ZBG,RHO_B,RHO_B,XP,
     *                       YP,0D0,BG)
            BGA_2D(j)=BGA_2D(j)+BG
           ENDIF

         ENDIF
        
!!         write(10,'(9F12.2)')XP*1e-3,XB(1)*1e-3,XB(2)*1e-3,
!!     *zB(1)*1D-3,zB(2)*1D-3,GEO,GR,RHO_up,RHO_down
         ENDDO
!---------------------------------------
        ENDDO
       ENDDO


      if(1==1) then   
!Reference model
      E_ref=0D0
      z_cref=-31.4D3
      z_lref=-151.2D3
      rho_cref=2780D0
      rho_mref=3245D0
      XB(1)=-d_x*.5-DELTA
      XB(2)=L_prf+d_x*.5+DELTA
!Initial value for 1D geoid reference level
      N_0=0
      write(*,*)'REFERENCE MODEL SUBSTRACTED'
      DO j=1,I_pr_end
       XP=H_dist(j)
       E_pt=bod_pr(j,0) 
         DO il=1,3
         IF (IL==1) THEN
         ZB(1)=0D0
         ZB(2)=z_cref
         rho_ref=rho_cref
         ELSEIF (IL==2) THEN
         ZB(1)=z_cref
         ZB(2)=z_lref
         rho_ref=rho_mref
         ELSE
         ZB(1)=z_lref
         ZB(2)=z_bot
         rho_ref=rho_a
         ENDIF
!Calculate 1D geoid reference level
         IF (j==1) THEN
          N_0=N_0+(pi*G_u/g)*((ZB(2)**2-ZB(1)**2)*rho_ref)
         ENDIF
!Calculate 1D geoid reference level
        IF (E_pt>0) THEN
          call GRAV_GRAD3D(XB,YB,ZB,rho_ref,rho_ref,XP,YP,
     *                     E_pt,GR)
          FA_2D(j)=FA_2D(j)-GR
          call GEO_GRAD3D(XB,YB,ZB-E_pt,rho_ref,rho_ref,XP,YP,0D0,
     *                    GEO)
          GEOID_2D(j)=GEOID_2D(j)-GEO
         ELSE
          call GEOGRAV_GRAD3D(XB,YB,ZB,rho_ref,rho_ref,XP,YP,0D0,
     *                        GEO,GR)
          FA_2D(j)=FA_2D(j)-GR
          GEOID_2D(j)=GEOID_2D(j)-GEO
         ENDIF
        ENDDO
      ENDDO
      endif
!1D Geoid anomaly--dipole moment++++++++++++++++++
!Loop for columns
       DO i=1,I_pr_end
        z_lims=bod_pr(i,:)
        CALL geoid_1D_dipole(-z_lims,N_lay,d_N,E_2D_v(i),T_mh_prf(i)) 
        GEOID_1D(i)=d_N+N_0
       ENDDO
       write(*,*)'E_2D_v',minval(E_2D_v),maxval(E_2D_v)
!1D Geoid anomaly--dipole moment+++++++++++++++++

!      FA_av=(SUM(FA_2D)-SUM(OBS_prf(:,6)))/(I_pr_end)
!      FA_2D=FA_2D-FA_av
!      BGA_2D=FA_2D+BGA_2D 
!      GEOID_av=(SUM(GEOID_2D)-SUM(OBS_prf(:,5)))/(I_pr_end)
!      GEOID_2D=GEOID_2D-GEOID_av 
!      BGA_av=(SUM(BGA_2D)-SUM(OBS_prf(:,7)))/(I_pr_end)
!      write(*,*)'FA_av,GEOID_av',FA_av,GEOID_av 
!Eliminate tilting
       CALL LNR_REG(H_dist,FA_2D,I_pr_end,a,b)
       CALL LNR_REG(H_dist,OBS_prf(:,6),I_pr_end,a_obs,b_obs)
       IF (a_obs==0.AND.b_obs==0) THEN 
        CALL LNR_REG(H_dist,OBS_prf(:,2),I_pr_end,a_obs,b_obs)
       ENDIF
       write(*,*)'FA_2D',a,b,a_obs,b_obs
       DO ii=1,I_pr_end
        FA_2D(ii)=FA_2D(ii)-(a-a_obs)*(H_dist(ii))-(b-b_obs)
       ENDDO

       CALL LNR_REG(H_dist,GEOID_2D,I_pr_end,a,b)
       CALL LNR_REG(H_dist,OBS_prf(:,5),I_pr_end,a_obs,b_obs)
       IF (a_obs==0.AND.b_obs==0) THEN
!For syn models|||||||||||||||||||||||||||||||||||
         IF(OBS_prf(INT(I_pr_end/2),1)>0) THEN
          GEOID_1D=GEOID_1D-minval(GEOID_1D)
          GEOID_2D=GEOID_2D-minval(GEOID_2D)
          OBS_prf(:,1)=OBS_prf(:,1)-minval(OBS_prf(:,1))
         ELSE
          GEOID_1D=GEOID_1D-maxval(GEOID_1D)
          GEOID_2D=GEOID_2D-maxval(GEOID_2D)
          OBS_prf(:,1)=OBS_prf(:,1)-maxval(OBS_prf(:,1))
         ENDIF
!For syn models|||||||||||||||||||||||||||||||||||
       ELSE
        GEOID_1D=GEOID_1D-(SUM(GEOID_1D)-SUM(OBS_prf(:,5)))/I_pr_end
        GEOID_2D=GEOID_2D-(a-a_obs)*H_dist-(b-b_obs) 
       ENDIF
       write(*,*)'GE,1D,2D,3D',maxval(GEOID_1D),maxval(GEOID_2D),
     *            maxval(OBS_prf(:,1))
!!        write(*,*)'mean geoid',SUM(OBS_prf(:,1))/I_pr_end

       BGA_2D=FA_2D+BGA_2D
       CALL LNR_REG(H_dist,BGA_2D,I_pr_end,a,b)
       CALL LNR_REG(H_dist,OBS_prf(:,7),I_pr_end,a_obs,b_obs)
       IF (a_obs==0.AND.b_obs==0) THEN
!For syn models|||||||||||||||||||||||||||||||||||
         IF(OBS_prf(INT(I_pr_end/2),3)>0) THEN
          BGA_2D=BGA_2D-minval(BGA_2D)
          OBS_prf(:,3)=OBS_prf(:,3)-minval(OBS_prf(:,3))
         ELSE
          write(*,*)'2D,3D',minval(BGA_2D),maxval(BGA_2D),
     *    minval(OBS_prf(:,3)),maxval(OBS_prf(:,3))
          BGA_2D=BGA_2D-maxval(BGA_2D)
          OBS_prf(:,3)=OBS_prf(:,3)-maxval(OBS_prf(:,3))
         ENDIF
!For syn models|||||||||||||||||||||||||||||||||||
       ELSE
        BGA_2D=BGA_2D-(a-a_obs)*H_dist-(b-b_obs)  
       ENDIF 
       write(*,*)'BG',a,b,a_obs,b_obs

       bod_pr=-bod_pr*1E-3
       RETURN

!AHORA NO SE EJECUTA
!Plot 2D calc
!      CALL PGSWIN (0.,MAXVAL(H_dist(:)*1E-3),MINVAL(GEOID_2D),
!     *             MAXVAL(GEOID_2D))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL PGSVP(.07,.91,.03,.25)
      CALL PGSWIN (0.,MAXVAL(H_dist(:)*1E-3),-150.,
     *             150.)
      CALL PGSCH(.7)
      CALL PGBOX('bntsi',0.0,0,'bctsivmn',0.0,0)
      CALL PGLAB('HORIZONTAL DISTANCE (km) ', 'DEPTH','')
      CALL PGSCI(3)
      CALL PGLINE(I_pr_end,H_dist*1E-3,BGA_2D)
      CALL PGSCI(1) 
      CALL PGLINE(I_pr_end,H_dist*1E-3,
     *             OBS_prf(:,3))
      CALL PGSLS(4)
      CALL PGLINE(I_pr_end,H_dist*1E-3,
     *             OBS_prf(:,7)) 
       CALL PGSLS(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GEOID_2D
      CALL PGSVP(.07,.91,.25,.5)
      CALL PGSWIN (0.,MAXVAL(H_dist(:)*1E-3),-10.,
     *             15.)
      CALL PGSCH(.7)
      CALL PGBOX('bntsi',0.0,0,'bctsivmn',0.0,0)
      CALL PGLAB('HORIZONTAL DISTANCE (km) ', 'DEPTH','')
      CALL PGSCI(3)
      CALL PGLINE(I_pr_end,H_dist*1E-3,GEOID_2D)
      CALL PGSCI(2)
      CALL PGLINE(I_pr_end,H_dist*1E-3,GEOID_1D) 
      CALL PGSCI(1)
      CALL PGLINE(I_pr_end,H_dist*1E-3,
     *             OBS_prf(:,1))
      CALL PGSLS(4)
      CALL PGLINE(I_pr_end,H_dist*1E-3,
     *             OBS_prf(:,5))
       CALL PGSLS(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL PGSVP(.07,.91,.5,.75)
      CALL PGSWIN (0.,MAXVAL(H_dist(:)*1E-3),-150.,
     *             150.)
      CALL PGSCH(.7)
      CALL PGBOX('bntsi',0.0,0,'bctsivmn',0.0,0)
      CALL PGLAB('HORIZONTAL DISTANCE (km) ', 'DEPTH','')
      CALL PGSCI(3)
      CALL PGLINE(I_pr_end,H_dist*1E-3,FA_2D)
      CALL PGSCI(1)
      CALL PGLINE(I_pr_end,H_dist*1E-3,
     *             OBS_prf(:,2))
      CALL PGSLS(4)
      CALL PGLINE(I_pr_end,H_dist*1E-3,
     *             OBS_prf(:,6))
       CALL PGSLS(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL PGSVP(.07,.91,.75,1.)
      CALL PGSWIN (0.,MAXVAL(H_dist(:)*1E-3),-3500.,
     *             2000.)
      CALL PGSCH(.7)
      CALL PGBOX('bntsi',0.0,0,'bctsivmn',0.0,0)
      CALL PGLAB('HORIZONTAL DISTANCE (km) ', 'DEPTH','')
      CALL PGSCI(3)
      CALL PGLINE(I_pr_end,H_dist*1E-3,E_2D_v)
      CALL PGSCI(1)
      CALL PGLINE(I_pr_end,H_dist*1E-3,
     *             OBS_prf(:,4))
      CALL PGSLS(4)
      CALL PGLINE(I_pr_end,H_dist*1E-3,
     *             OBS_prf(:,8))
       CALL PGSLS(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1



      CALL pgband(0,1,xdm,ydm,xdm,ydm,chdm)
      CALL PGERAS 
      IF (ICHAR(chdm)==27) EXIT
      ENDDO
      END SUBROUTINE
