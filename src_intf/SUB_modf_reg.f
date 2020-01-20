! Subroutine SUB_modf_reg.f

      SUBROUTINE modf_reg
      USE M_imag
      USE M_profile
      USE M_material
      USE M_label 
      
      CHARACTER(1) le
      CHARACTER(35) lbl
      CHARACTER(80) TEXT
      CHARACTER(80) TEXT_ud
      REAL(4) xdm,ydm,x1,y1,x2,y2,z1,z2
      REAL(4) x,y,z_obs_max,z_obs_min,ytick,ysub
      REAL(4) slp,dm,dm_p
      REAL(4), ALLOCATABLE:: bod_PG(:,:)
      REAL(4) R,G,B
      INTEGER kx,ky,ii_dm,i_pr_old,hh
      INTEGER,ALLOCATABLE:: pr_pos_dum(:)
      REAL(4), ALLOCATABLE:: xy_pr_p_dum(:)
      REAL(4), ALLOCATABLE:: lat_vec1(:,:),lat_vec2(:,:)
      INTEGER N_lat,i_lat,low_fin,up_fin,left_fin,right_fin,i_del
      INTEGER PGOPEN
      CHARACTER(2) chdm


      
!Bodies info
      open(15,file='layers.info',status='OLD')
      open(16,file='color.rgb',status='UNKNOWN')
      DO i=1,N_lay+1
       IF (i==1) THEN
       read(15,*)
       ELSE
       read(15,*)j,prop(j)%name_mat,prop(j)%rho,
     *                      prop(j)%alpha,prop(j)%K,
     *                      prop(j)%A,prop(j)%beta,prop(j)%gmma_T,
     *                      prop(j)%K_0p,prop(j)%K_T,prop(j)%litho 
!Create color file
       R=1.
!       G=(rho_max-prop(j)%rho)/(rho_max-rho_min)
       G=1-((ABS(prop(j)%rho)-rho_min)/(rho_max-rho_min))**2
       B=0.
       write(16,*)R,G,B
       ENDIF
      ENDDO
      close (15)
      close(16)
!'SELECT LAYER TO MODIFY...'
      CALL lay_sel
      IF (LBL_SGN==1) RETURN
      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *               5)
!!       IF (updown_lay=='U') THEN
!!       TEXT_ud='PUSH UP LAYER'
!!       ELSE
!!        IF (ICHAR(updown_lay)==27) RETURN
!!       TEXT_ud='PUSH DOWN LAYER'
!!       ENDIF
!      write(*,*)'NUMBER OF PROFILES : '
!      write(*,*)'(If = 0 -->PROFILES ARE SELECTED BY MOUSE)' 
!      read(*,*) N_pr
       write(*,*)'N_pr: ',N_pr
!      write(*,*)'East-West----> E'
!      write(*,*)'North-South--> N'
!      read(*,*) ort
       IF (ort=='E') THEN
        lbl='TRANSECT W-E'
!Band width: project reference marks 
        bnd_lbl=1*p_m_y       
!Interpolation(longitudinal) parameter
        N_lat=NINT(IDIM*.04)
        IF (N_lat<1) N_lat=1 
       ELSE
!Band width: project reference marks 
        bnd_lbl=1*p_m_x
        lbl='TRANSECT S-N'
!Interpolation(longitudinal) parameter
        N_lat=NINT(JDIM*.04)
        IF (N_lat<1) N_lat=1
       ENDIF
      write(*,*)'PROFILE ORIENTATION:  ',lbl
      write(*,*)'SELECT REGION...'
      write(*,*)'BANDWIDTH OF LONGITUDINAL INTERPOLATION IN Nº OF 
     *PROFILES: ',N_lat
      IF (glue_l) THEN
       CALL glue('C',0,0)
       WRITE(*,*)'GLUE ON'
      ELSE
       WRITE(*,*)'GLUE OFF'
      ENDIF
      write(*,*)'Z: CHANGE VERTICAL SCALE'
!****************************************
       IF(LINUX) THEN
          IPG_LAY=PGOPEN('/XWINDOW')
       ELSE
          IPG_LAY=PGOPEN('/WX')
       ENDIF
       IF (IPG_LAY<=0) STOP
       CALL pgask(.false.)
       CALL plot
       CALL pgqvp(0,r1,r2,r3,r4) 
       CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *               lat_max,7)
!**************************************** 
      CALL pgask(.false.) 
      CALL PGSFS(3)
      xdm=(lon_max+lon_min)/2
      ydm=(lat_max+lat_min)/2
      CALL pgband(0,1,0.,0.,xdm,ydm,le)
      IF (le=='x'.OR.ICHAR(le)==27) RETURN
      x1=xdm
      y1=ydm
      CALL pgband(2,1,xdm,ydm,xdm,ydm,le)
      IF (le=='x'.OR.ICHAR(le)==27) RETURN
      x2=xdm
      y2=ydm
      if (x2>x1) then
      lnz_min=x1
      lnz_max=x2
      else
      lnz_min=x2
      lnz_max=x1
      endif
      if (y2>y1) then
      ltz_min=y1
      ltz_max=y2
      else
      ltz_min=y2
      ltz_max=y1
      endif

      IF (lnz_min<lon_min) lnz_min=lon_min
      IF (lnz_max>lon_max) lnz_max=lon_max 
      IF (ltz_min<lat_min) ltz_min=lat_min
      IF (ltz_max>lat_max) ltz_max=lat_max

      i_mreg_min=nint((lnz_min-lon_min)*IDIM/(lon_max-lon_min))
      i_mreg_max=nint((lnz_max-lon_min)*IDIM/(lon_max-lon_min))
      IF (i_mreg_min<1) i_mreg_min=1
      IF (i_mreg_max>IDIM) i_mreg_max=IDIM
      i_mreg=i_mreg_max-i_mreg_min+1

      j_mreg_min=nint((ltz_min-lat_min)*JDIM/(lat_max-lat_min))
      j_mreg_max=nint((ltz_max-lat_min)*JDIM/(lat_max-lat_min))
      IF (j_mreg_min<1) j_mreg_min=1
      IF (j_mreg_max>JDIM) j_mreg_max=JDIM
      j_mreg=j_mreg_max-j_mreg_min+1
      IF (N_pr<=0) THEN
       MOUSE='on' 
       IF (ort=='E') THEN
        N_pr=j_mreg
       ELSE
        N_pr=i_mreg
       ENDIF 
       WRITE(*,*)'NUMBER OF POINTS IN PGPOLY',N_pr
!      write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
!      write(*,*)'KEYS & ACTIONS'
!      write(*,*)'s---------->SAVE,GO TO NEXT '
!      write(*,*)'x--->INTERPOLATE,GO TO NEXT '
!      write(*,*)'m--------->SAVE,SELECT NEXT '
!      write(*,*)'n-->INTERPOLATE,SELECT NEXT '
!      write(*,*)'d-------> DELETE LAST POINT '
!      write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      ELSE
       MOUSE='off'
!      write(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@'
!      write(*,*)'KEYS & ACTIONS'
!      write(*,*)'s--------> NEXT PROFILE'
!      write(*,*)'d---> DELETE LAST POINT'
!      write(*,*)'@@@"@@@@@@@@@@@@@@@@@@@@'
      ENDIF
      IF(ALLOCATED(pr_pos)) DEALLOCATE(pr_pos)
      IF(ALLOCATED(xy_pr_p)) DEALLOCATE(xy_pr_p)
      IF(ALLOCATED(pr_pos_dum)) DEALLOCATE(pr_pos_dum)
      IF(ALLOCATED(xy_pr_p_dum)) DEALLOCATE(xy_pr_p_dum) 
      ALLOCATE (pr_pos(0:N_pr+1),xy_pr_p(0:N_pr+1),pr_pos_dum(0:N_pr+1),
     *          xy_pr_p_dum(0:N_pr+1))
      pr_pos=0
! Open a window with a legend of the layers
       CALL legend
      
!Loop for transects      
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      dm_lay=layers
      lems='X'
      fw_pr=1
      i_pr=0
                     
      DO WHILE (i_pr<=N_pr)
      i_pr_old=i_pr 
      i_pr=i_pr+fw_pr

!E-W oriented profiles
       IF (ort=='E') THEN
       ifl=i_mreg
       ifl_p=j_mreg
       IF (i_pr==1) THEN
         IF(ALLOCATED(E_pr)) DEALLOCATE(E_pr)
         IF(ALLOCATED(bod_pr)) DEALLOCATE(bod_pr)
         IF(ALLOCATED(bod_sf)) DEALLOCATE(bod_sf)
         IF(ALLOCATED(bod_pr_bak)) DEALLOCATE(bod_pr_bak)
!cambio
         ALLOCATE (E_pr(ifl),bod_pr(ifl,0:N_lay),bod_sf(ifl,0:ifl_p+1,
     *             0:N_lay),bod_pr_bak(ifl))
         bod_sf=0
       ENDIF
        IF (MOUSE=='on') THEN
         d_pr=p_m_y
        ELSE
         d_pr=(ltz_max-ltz_min)/(N_pr-1)
        ENDIF
       y1=ltz_min+d_pr*(i_pr-1)
       y_pl(:)=y1
       xy_pr_min=lnz_min
       xy_pr_max=lnz_max
       xy_pr_p_min=ltz_min
       xy_pr_p_max=ltz_max
       xy_pr_p(i_pr)=y1
       ijdm=nint((y1-lat_min)*JDIM/(lat_max-lat_min))
        IF (ijdm<1) ijdm=1
        IF (ijdm>JDIM) ijdm=JDIM
       E_pr=E(i_mreg_min:i_mreg_max,ijdm)
!cambio
        bod_pr(:,0)=-E_pr*1D-3
        DO i=1,N_lay
        bod_pr(:,i)=layers(i_mreg_min:i_mreg_max,ijdm,i)
        ENDDO
       L=xy_pr_max-xy_pr_min
       L_p=xy_pr_p_max-xy_pr_p_min
       I_p=nint ((y1-xy_pr_p_min)*ifl_p/L_p)
       if (I_p>ifl_p) I_p=ifl_p

!N-S oriented profiles
       ELSE
       ifl=j_mreg
       ifl_p=i_mreg
       IF (i_pr==1) THEN 
         IF(ALLOCATED(E_pr)) DEALLOCATE(E_pr)
         IF(ALLOCATED(bod_pr)) DEALLOCATE(bod_pr)
         IF(ALLOCATED(bod_sf)) DEALLOCATE(bod_sf)
         IF(ALLOCATED(bod_pr_bak)) DEALLOCATE(bod_pr_bak)
!cambio
         ALLOCATE (E_pr(ifl),bod_pr(ifl,0:N_lay),bod_sf(ifl,0:ifl_p+1,
     *             0:N_lay),bod_pr_bak(ifl))
         bod_sf=0
       ENDIF
        IF (MOUSE=='on') THEN
         d_pr=p_m_x
        ELSE
         d_pr=(lnz_max-lnz_min)/N_pr
        ENDIF
       x1=lnz_min+d_pr*(i_pr-1)
       x_pl(:)=x1
       xy_pr_min=ltz_min
       xy_pr_max=ltz_max
       xy_pr_p_min=lnz_min
       xy_pr_p_max=lnz_max
       xy_pr_p(i_pr)=x1
       ijdm=nint((x1-lon_min)*IDIM/(lon_max-lon_min))
        IF (ijdm<1) ijdm=1
        IF (ijdm>IDIM) ijdm=IDIM
        E_pr=E(ijdm,j_mreg_min:j_mreg_max) 
        bod_pr(:,0)=-E_pr*1D-3 
        DO i=1,N_lay
        bod_pr(:,i)=layers(ijdm,j_mreg_min:j_mreg_max,i)
        ENDDO
       L=xy_pr_max-xy_pr_min
       L_p=xy_pr_p_max-xy_pr_p_min
       I_p=nint ((x1-xy_pr_p_min)*ifl_p/L_p)
       if (I_p>ifl_p) I_p=ifl_p
       ENDIF

       IF (MOUSE=='on') I_p=i_pr

       IF (i_pr==1) THEN 
        IF(ALLOCATED(xy_pr)) DEALLOCATE(xy_pr) 
        IF(ALLOCATED(bod_PG)) DEALLOCATE(bod_PG) 
        IF(ALLOCATED(lat_vec1)) DEALLOCATE (lat_vec1)
        IF(ALLOCATED(lat_vec2)) DEALLOCATE (lat_vec2)
        ALLOCATE (xy_pr(ifl))
        ALLOCATE (bod_PG(2,2*ifl))
        ALLOCATE (lat_vec1(0:N_lat,0:ifl_p+1),
     *            lat_vec2(0:N_lat,0:ifl_p+1))
        lat_vec1=0D0
        lat_vec2=0D0
!!        DO j=1,ifl
!!        xy_pr(j)=xy_pr_min+(j-1)*p_m
!!       ENDDO
!Save extreme profiles
        IF (ort=='E') THEN
         DO j=1,ifl
          xy_pr(j)=xy_pr_min+(j-1)*p_m_x
         ENDDO
         in_min=j_mreg_min-1
         in_max=j_mreg_max+1
         IF (in_min<1) THEN
          low_fin=1
          in_min=1
          xy_pr_p(0)=ltz_min
         ELSE
          low_fin=0
          xy_pr_p(0)=ltz_min-d_pr
         ENDIF 
         IF (in_max>JDIM) THEN
          up_fin=1
          in_max=JDIM
          xy_pr_p(N_pr+1)=ltz_max
         ELSE
          up_fin=0
          xy_pr_p(N_pr+1)=ltz_max+d_pr
         ENDIF
         i_del=in_max-in_min+1+(low_fin-1)
!····························
         IF (i_mreg_min==1) THEN
         left_fin=1 
         in_lat=1
         ELSE
         left_fin=0
         in_lat=i_mreg_min-1
         ENDIF 
         IF (i_lay==0) THEN
          lat_vec1(0,low_fin:i_del)=-E(in_lat,in_min:in_max)*1D-3
         ELSE
          lat_vec1(0,low_fin:i_del)=layers(in_lat,in_min:in_max,i_lay)
         ENDIF
         IF (i_mreg_max==IDIM) THEN
         right_fin=1
         in_lat=IDIM
         ELSE
         right_fin=0
         in_lat=i_mreg_max+1
         ENDIF
         IF (i_lay==0) THEN
          lat_vec2(0,low_fin:i_del)=-E(in_lat,in_min:in_max)*1D-3
         ELSE
          lat_vec2(0,low_fin:i_del)=layers(in_lat,in_min:in_max,i_lay) 
         ENDIF
!····························
         bod_sf(:,0,0)=-E(i_mreg_min:i_mreg_max,in_min)*1D-3
         bod_sf(:,ifl_p+1,0)=-E(i_mreg_min:i_mreg_max,in_max)*1D-3   
         DO i=1,N_lay
         bod_sf(:,0,i)=layers(i_mreg_min:i_mreg_max,in_min,i)   
         bod_sf(:,ifl_p+1,i)=layers(i_mreg_min:i_mreg_max,in_max,i)
         ENDDO

        ELSE
         DO j=1,ifl
          xy_pr(j)=xy_pr_min+(j-1)*p_m_y
         ENDDO
         in_min=i_mreg_min-1
         in_max=i_mreg_max+1
         IF (in_min<1) THEN
          low_fin=1
          in_min=1
          xy_pr_p(0)=lnz_min
         ELSE
          low_fin=0
          xy_pr_p(0)=lnz_min-d_pr 
         ENDIF
         IF (in_max>IDIM) THEN
          up_fin=1
          in_max=IDIM
          xy_pr_p(N_pr+1)=lnz_max
          ELSE
          up_fin=0
          xy_pr_p(N_pr+1)=lnz_max+d_pr
         ENDIF
         i_del=in_max-in_min+1+(low_fin-1)
!··················
         IF (j_mreg_min==1) THEN
         left_fin=1
         in_lat=1
         ELSE
         left_fin=0
         in_lat=j_mreg_min-1
         ENDIF
         IF (i_lay==0) THEN
          lat_vec1(0,low_fin:i_del)=-E(in_min:in_max,in_lat)*1D-3
         ELSE
          lat_vec1(0,low_fin:i_del)=layers(in_min:in_max,in_lat,i_lay)
         ENDIF
         IF (j_mreg_max==JDIM) THEN
         right_fin=1
         in_lat=JDIM
         ELSE
         right_fin=0
         in_lat=j_mreg_max+1
         ENDIF
         IF (i_lay==0) THEN
          lat_vec2(0,low_fin:i_del)=-E(in_min:in_max,in_lat)*1D-3
         ELSE
          lat_vec2(0,low_fin:i_del)=layers(in_min:in_max,in_lat,i_lay)
         ENDIF
!····························
         bod_sf(:,0,0)=-E(in_min,j_mreg_min:j_mreg_max)*1D-3
         bod_sf(:,ifl_p+1,0)=-E(in_max,j_mreg_min:j_mreg_max)*1D-3 
         DO i=1,N_lay
         bod_sf(:,0,i)=layers(in_min,j_mreg_min:j_mreg_max,i)
         bod_sf(:,ifl_p+1,i)=layers(in_max,j_mreg_min:j_mreg_max,i)
         ENDDO
        ENDIF
       ENDIF

!Save intermediate profiles
       IF ((fw_pr>1.AND.lems=='X')) THEN
        DO i=i_pr_old+1,i_pr-1
         IF (ort=='E') THEN
         xy_pr_p(i)=ltz_min+d_pr*(i-1)
         ijdm=nint((xy_pr_p(i)-lat_min)*JDIM/(lat_max-lat_min))
         IF (i_lay==0) THEN
          bod_pr(:,i_lay)=-E(i_mreg_min:i_mreg_max,ijdm)*1D-3  
         ELSE
          bod_pr(:,i_lay)=layers(i_mreg_min:i_mreg_max,ijdm,i_lay)
         ENDIF
         pr_pos(i)=i 
         ELSE
         xy_pr_p(i)=lnz_min+d_pr*(i-1)
         ijdm=nint((xy_pr_p(i)-lon_min)*IDIM/(lon_max-lon_min))
         IF (i_lay==0) THEN
          bod_pr(:,i_lay)=-E(ijdm,j_mreg_min:j_mreg_max)*1D-3
         ELSE
          bod_pr(:,i_lay)=layers(ijdm,j_mreg_min:j_mreg_max,i_lay)
         ENDIF
         pr_pos(i)=i
         ENDIF
         bod_sf(:,pr_pos(i),i_lay)=bod_pr(:,i_lay)
        ENDDO
       ENDIF 
!Erase profiles       
       IF (fw_pr<0) pr_pos(i_pr:N_pr)=0
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!GRAPHICAL MODIFICATION
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ch_v_scl=.FALSE.
      z_max_mdfr=MAXVAL(bod_pr(:,i_lay))
      Do
!Display profile location 
!#########################################################
      DO iw_pl=1,2
      IF(iw_pl==1) THEN
       CALL PGSLCT(ISTAT)
      ELSE
       CALL PGSLCT(IPG_LAY) 
      ENDIF 
      CALL PGSFS(3) 
      CALL PGRECT(lnz_min,lnz_max,ltz_min,ltz_max) 
      CALL PGSLW(10)     
      IF (ort=='E') THEN
      x_pl(1)=xy_pr_min
      x_pl(2)=xy_pr_max
      CALL PGLINE(2,x_pl,y_pl)     
      ELSE
      y_pl(1)=xy_pr_min
      y_pl(2)=xy_pr_max
      CALL PGLINE(2,x_pl,y_pl)
      ENDIF  
      CALL pgscr(17,.5,.5,.5)
      CALL PGSCI(17)
      CALL PGSLW(3)
      CALL PGRECT(lnz_min,lnz_max,ltz_min,ltz_max)
      ENDDO
!#########################################################

      IF(LINUX) THEN
         IPR=PGOPEN('/XWINDOW') 
      ELSE
         IPR=PGOPEN('/WX') 
      ENDIF
      IF (IPR.LE.0) STOP
!Vertical scale
      IF (i_lay==0) THEN
       IF (z_max_mdfr>400.) z_max_mdfr=MAXVAL(bod_pr(:,i_lay))+1.
       z_max_mdfr=z_max_mdfr+1.
      ELSE
       IF (z_max_mdfr>400.) z_max_mdfr=MAXVAL(bod_pr(:,i_lay))+10.
       z_max_mdfr=z_max_mdfr+10.
      ENDIF 

      CALL pgscr(0,1.,1.,1.)
      CALL pgscr(1,0.,0.,0.)
      CALL pgask(.false.)
      CALL PGPAP(8.0,.8)
      CALL PGSVP(.075,.92,.1,.9)
!      write(TEXT,11)TRIM(prop(i_lay)%name_mat),i_lay
! 11   FORMAT (A,'--> ',I2)
      CALL PGQVP(0,r1,r2,r3,r4)
!/////////////Plot Observable profile////////////////////////
      IF(ALLOCATED(OBS_pr)) DEALLOCATE(OBS_pr)
      ALLOCATE(OBS_pr(ifl))
      IF (ort=='E') THEN
       OBS_pr=OBS(i_obs,i_mreg_min:i_mreg_max,ijdm)
      ELSE
       OBS_pr=OBS(i_obs,ijdm,j_mreg_min:j_mreg_max)
      ENDIF
      z_obs_max=MAXVAL(OBS_pr)
      z_obs_min=MINVAL(OBS_pr)
      IF (z_obs_max>0.AND.z_obs_min>0) z_obs_min=0
      z_obs_max=z_obs_max+(MAXVAL(OBS_pr)-MINVAL(OBS_pr))/10
      z_obs_min=z_obs_min-(MAXVAL(OBS_pr)-MINVAL(OBS_pr))/10
      IF (z_obs_max>0.AND.z_obs_min>0) z_obs_min=0
      IF (z_obs_max<0.AND.z_obs_min<0) z_obs_max=0 
      IF(z_obs_max==z_obs_min) z_obs_min=-z_obs_max
      CALL PGSVP(r1,r2,.75,.9)
      CALL PGSWIN (xy_pr_min,xy_pr_max,z_obs_min,z_obs_max)
      CALL PGBOX('ab',0.0,0,'bctivmn',0.0,0) 
      CALL PGLINE(ifl,xy_pr,OBS_pr)
!/////////////Plot Observable profile////////////////////////
      WRITE(title,10) i_pr,TRIM(lbl),i_lay
 10   FORMAT ('PROFILE: ',I3,', ',A,', ','MODIFYING LAYER--> ',
     *        ' ',I2)
      CALL text_plot(IPR,r1,r2,r3,r4,xy_pr_min,xy_pr_max,z_obs_min,
     *              z_obs_max, 6)
      CALL PGSVP(r1,r2,r3,.70)
      z_min_mdfr=1.-MINVAL(bod_pr(:,0))
      CALL PGSWIN (xy_pr_min,xy_pr_max,-z_max_mdfr,z_min_mdfr)
!      ytick=NINT((z_max_mdfr-z_min_mdfr)/10.)
      ysub=5 
      CALL PGBOX('bntsi',0.,0,'bctsivmn',0.,ysub)
      CALL PGLAB('HORIZONTAL DISTANCE ', 'DEPTH','')
!Horizontal coordinates
      bod_PG(1,1:ifl)=xy_pr(ifl:1:-1)
      bod_PG(1,ifl+1:2*ifl)=xy_pr
!Water layer
      CALL PGSCR (50,.0,.97,.97)
      bod_PG(2,1:ifl)=0.
      bod_PG(2,ifl+1:2*ifl)=-8.
      CALL PGSCI(50)
      CALL PGPOLY(2*ifl,bod_PG(1,1:2*ifl),bod_PG(2,1:2*ifl))
!Body layers
!Smooth layers
!!      DO i=1,N_lay
!!       CALL smooth1D(bod_pr(:,i),ifl)
!!      ENDDO
      open(16,file='color.rgb',status='UNKNOWN')
      DO i=0,N_lay-1
      read(16,*) R,G,B
      CALL PGSCR (15+i,R,G,B)
      CALL PGSCI (15+i)
      IF (i==0) THEN
      bod_PG(2,1:ifl)=E_pr(ifl:1:-1)*1e-3
      bod_PG(2,ifl+1:2*ifl)=-bod_pr(:,1)
      ELSE  
      bod_PG(2,1:ifl)=-bod_pr(ifl:1:-1,i)
      bod_PG(2,ifl+1:2*ifl)=-bod_pr(:,i+1)
      ENDIF
      CALL PGPOLY(2*ifl,bod_PG(1,1:2*ifl),bod_PG(2,1:2*ifl))
      ENDDO 
      close (16)
!Asthenospheric layer
      CALL PGSCR (52,.97,.02,.95)
      bod_PG(2,1:ifl)=-bod_pr(ifl:1:-1,N_lay)
      bod_PG(2,ifl+1:2*ifl)=-z_max_mdfr
      CALL PGSCI(52)
      CALL PGPOLY(2*ifl,bod_PG(1,1:2*ifl),bod_PG(2,1:2*ifl)) 
!Draw line between bodies
      IF (i_lay==0) THEN
       CALL PGSLW(7)
       CALL PGSLS(4)
      ELSE
       CALL PGSLW(3)
       CALL PGSCI(1)
      ENDIF
      CALL PGLINE(ifl,xy_pr,E_pr*1e-3)
      DO i=1,N_lay
       CALL PGSCI(1)
       IF (i==i_lay) THEN
       CALL PGSLW(7)
       CALL PGSLS(4)
       ELSE
       CALL PGSLW(3)
       CALL PGSLS(1)
       ENDIF 
       CALL PGLINE(ifl,xy_pr,-bod_pr(:,i))
      ENDDO
      CALL PGSLW(3)
      CALL PGSLS(1)
      CALL PGIDEN
!Plot "backup line" (last profile modified)
       IF (i_pr>1.OR.fw_pr<0) THEN
        CALL PGSLW(5)
        CALL PGSLS(3)
        CALL PGSCI(3)
        CALL PGLINE(ifl,xy_pr,-bod_pr_bak(:))
        CALL PGSLW(1)
        CALL PGSLS(1)
       ENDIF
       CALL PGSCI(1)
!Plot the available resolution in longitudinal nodes
       bod_PG(2,ifl+1:2*ifl)=-z_max_mdfr+z_max_mdfr/15
       CALL PGSLW(3)
       CALL PGPT(ifl,bod_PG(1,1:ifl),bod_PG(2,ifl+1:2*ifl),5)
       CALL PGSLW(1)
!Plot projected labels
      CALL PGSVP(r1,r2,r3,.70)
      CALL label_proj
!Plot projected reference points
      CALL PGSVP(r1,r2,r3,r4)
      CALL ref_points_proj
      CALL PGSVP(r1,r2,r3,.70)
!      CALL PGLAB('','',TRIM(obs_file(i_obs)))
      CALL PGPTEXT((xy_pr_min+xy_pr_max)/2,5.,0.,0.5,
     *             TRIM(obs_file(i_obs)))
!      write(*,*)'$$$$$$$$$$$$$$$$$$$$'
!      write(*,*)'PROFILE Nº--> ',i_pr
!      write(*,*)'$$$$$$$$$$$$$$$$$$$$'
!Modify layers
       CALL new_coord
       IF (exit_mdreg=='y') THEN
        CALL PGSLCT(I_LEG)
        CALL PGCLOS
        CALL PGSLCT(IPR)
        CALL PGCLOS
        CALL PGSLCT(ISTAT)
        CALL clear_scr(0)
        EXIT
       ENDIF
       IF (abort_mdreg=='y') THEN
        IF(ALLOCATED(xy_pr)) DEALLOCATE(xy_pr)
        IF(ALLOCATED(bod_sf)) DEALLOCATE(bod_sf)
        IF(ALLOCATED(E_pr)) DEALLOCATE(E_pr)
        IF(ALLOCATED(bod_pr)) DEALLOCATE(bod_pr)
        IF(ALLOCATED(bod_PG)) DEALLOCATE(bod_PG)
        IF(ALLOCATED(bod_pr_bak)) DEALLOCATE(bod_pr_bak) 
        CALL PGSLCT(I_LEG)
        CALL PGCLOS
        CALL PGSLCT(IPG_LAY)
        CALL PGCLOS 
        CALL PGSLCT(IPR) 
        CALL PGCLOS
        CALL PGSLCT(ISTAT)
        CALL clear_scr(0)
        RETURN
       ENDIF
      CALL PGSLCT(IPR)
      CALL PGCLOS 
      CALL PGSLCT(IPG_LAY)
      CALL plot
      CALL PGSLCT(ISTAT)
      CALL clear_scr(0)
      IF (.NOT.ch_v_scl) EXIT
      
      ENDDo
       IF (exit_mdreg=='y') EXIT
      ENDDO
!Close the legend window in the case of exiting by saving the last profile
      IF (exit_mdreg=='n') THEN
       CALL PGSLCT(I_LEG)
       CALL PGCLOS 
      ENDIF 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!Redimension pr_pos and redefine N_pr for the interpolation
       IF (MOUSE=='on') THEN
        i_nonzero=0
        DO kh=1,N_pr
         IF(pr_pos(kh)>0) i_nonzero=i_nonzero+1
        ENDDO
        pr_pos_dum=pr_pos
        xy_pr_p_dum=xy_pr_p
        DEALLOCATE(pr_pos,xy_pr_p)
        ALLOCATE(pr_pos(0:i_nonzero+1),xy_pr_p(0:i_nonzero+1))
        ii_dm=1
        DO i=1,N_pr
         IF (pr_pos_dum(i)>0) THEN
          pr_pos(ii_dm)=pr_pos_dum(i)
          xy_pr_p(ii_dm)=xy_pr_p_dum(i)
          ii_dm=ii_dm+1
         ENDIF
        ENDDO
        pr_pos(0)=0
        pr_pos(i_nonzero+1)=N_pr+1   
        xy_pr_p(0)=xy_pr_p_dum(0)
        xy_pr_p(i_nonzero+1)=xy_pr_p_dum(N_pr+1)
        
        N_pr=i_nonzero
       ENDIF 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DO i=0,N_pr
       IF (i==0.AND.low_fin==1) CYCLE
       IF (i==N_pr.AND.up_fin==1) CYCLE
       IF (ort=='E') THEN
        p_m=p_m_y
       ELSE
        p_m=p_m_x
       ENDIF
       DO j=1,ifl
        z1=bod_sf(j,pr_pos(i),i_lay)
        z2=bod_sf(j,pr_pos(i+1),i_lay) 
         d_pr=(pr_pos(i+1)-pr_pos(i))*p_m
         slp=(z2-z1)/d_pr
         dm=(j-1)*L/ifl+xy_pr_min
        DO k=pr_pos(i),pr_pos(i+1)
         bod_sf(j,k,i_lay)=slp*((k-pr_pos(i))*p_m)+z1
        ENDDO

       ENDDO
      ENDDO
      IF (left_fin==1) lat_vec1(0,:)=bod_sf(1,:,i_lay)
      IF (right_fin==1)lat_vec2(0,:)=bod_sf(ifl,:,i_lay)
      DO i_lat=1,N_lat
      lat_vec1(i_lat,:)=lat_vec1(i_lat-1,:)*(N_lat-i_lat)/N_lat+
     *                  bod_sf(i_lat,:,i_lay)*i_lat/N_lat 
      lat_vec2(i_lat,:)=lat_vec2(i_lat-1,:)*(N_lat-i_lat)/N_lat+
     *                  bod_sf(ifl-i_lat+1,:,i_lay)*i_lat/N_lat
      ENDDO 
      bod_sf(1:N_lat,:,i_lay)=lat_vec1(1:N_lat,:)
      bod_sf(ifl:ifl-N_lat+1:-1,:,i_lay)=lat_vec2(1:N_lat,:)
      DO j=1,ifl
       DO k=0,ifl_p+1
        IF (k==0.AND.low_fin==1) CYCLE
        IF (k==ifl_p+1.AND.up_fin==1) CYCLE

        IF (ort=='E') THEN
        i_gl=i_mreg_min+(j-1)
        j_gl=j_mreg_min+(k-1)
        IF (j_gl<1) j_gl=1
        IF (j_gl>JDIM) j_gl=JDIM
        ELSE
        i_gl=i_mreg_min+(k-1)
        j_gl=j_mreg_min+(j-1)
        IF (i_gl<1) i_gl=1
        IF (i_gl>IDIM) i_gl=IDIM
        ENDIF
        IF (i_lay==0) THEN
         E(i_gl,j_gl)=bod_sf(j,k,i_lay)
         DO hh=1,N_lay
          IF (dm_lay(i_gl,j_gl,hh)<E(i_gl,j_gl))
     *        dm_lay(i_gl,j_gl,hh)=E(i_gl,j_gl)
         ENDDO
         E(i_gl,j_gl)=-E(i_gl,j_gl)*1D3
        ELSE
         dm_lay(i_gl,j_gl,i_lay)=bod_sf(j,k,i_lay)
         CALL lay_order(i_gl,j_gl,i_lay)
        ENDIF
        layers(i_gl,j_gl,1:N_lay)=dm_lay(i_gl,j_gl,:)
        IF (glue_l) CALL glue('U',i_gl,j_gl)
       ENDDO
      ENDDO 
!Write the layer files
      DO i=1,N_lay
       IF (i>9) THEN
       WRITE(chdm,'(I2)')i
       ELSE
       WRITE(chdm,'(I1)')i
       ENDIF
       fil=TRIM(string)//TRIM(chdm)//'.xyz'
       open(20,file=fil,status='UNKNOWN')
       do ky=JDIM,1,-1
        y=(ky-1)*p_m_y+lat_min
         do kx=1,IDIM
          x=(kx-1)*p_m_x+lon_min
          write (20,*)x,y,layers(kx,ky,i)
         enddo
       enddo
       close(20)
      ENDDO
!cambio
      OBS(1,:,:)=E(:,:)
      open(20,file='./GEO_DATA/elev_dat.xyz',status='UNKNOWN')
       do ky=JDIM,1,-1
        y=(ky-1)*p_m_y+lat_min
         do kx=1,IDIM
          x=(kx-1)*p_m_x+lon_min
          write (20,*)x,y,E(kx,ky)
         enddo
       enddo
       close(20)

      IF(ALLOCATED(xy_pr)) DEALLOCATE(xy_pr)
      IF(ALLOCATED(bod_sf)) DEALLOCATE(bod_sf)
      IF(ALLOCATED(E_pr)) DEALLOCATE(E_pr)
      IF(ALLOCATED(bod_pr)) DEALLOCATE(bod_pr)
      IF(ALLOCATED(bod_PG)) DEALLOCATE(bod_PG)
      IF(ALLOCATED(lat_vec1)) DEALLOCATE (lat_vec1)
      IF(ALLOCATED(lat_vec2)) DEALLOCATE (lat_vec2)
      IF(ALLOCATED(bod_pr_bak)) DEALLOCATE(bod_pr_bak)

!      CALL PGSLCT(I_LEG)
!      CALL PGCLOS
       
      CALL PGSLCT(IPG_LAY)
      CALL plot
      CALL pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *               lat_max,12)
      CALL pgband(0,1,0.,0.,xdm,ydm,le)
      CALL PGCLOS
      CALL PGSLCT(ISTAT)
 
 
      END SUBROUTINE 
