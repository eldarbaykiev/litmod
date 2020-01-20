! Subroutine  SUB_profile_screen.f


      SUBROUTINE profile_screen
      USE M_imag
      USE M_profile 
      USE M_material
      USE M_ctes,ONLY: pi,T_a,rho_a,Rd_T
      INTEGER:: i_prf,j_prf,PGOPEN,rho_m_av,n_prf_max
      INTEGER,ALLOCATABLE:: INDEX_prf(:,:)
      REAL(4):: a,lon,lat,dist,dmy,z_prf_max,z_obs_max(5),z_obs_min(5),
     *          margn
      REAL(4), ALLOCATABLE:: bod_PG(:,:),E_fil(:)
      CHARACTER(2) dmch,z_prf_lb
      CHARACTER(7) ort_lbl
      CHARACTER(50) txt_lbl,UCASE,name_bod,name_obs,T_txt 
      REAL(4) xbox(4),ybox(4)
!Isoterm parameters
      T0_iso=600
      dT_iso=200
      NC=4

      IF (ABS(y_pl(2)-y_pl(1))>ABS(x_pl(2)-x_pl(1))) THEN
!Reorder the extrems of the profile to satisfy: x_pl(2)>x_pl(1)
      IF (y_pl(2)<y_pl(1)) THEN
       dmy=y_pl(2)
       y_pl(2)=y_pl(1)
       y_pl(1)=dmy
       dmy=x_pl(2)
       x_pl(2)=x_pl(1)
       x_pl(1)=dmy
      ENDIF
!Profile slope
       a=(x_pl(2)-x_pl(1))/(y_pl(2)-y_pl(1))
!Number of longitudinal nodes in the profile
       I_pr_end=NINT(ABS((y_pl(2)-y_pl(1))/p_m_y))+1
      ELSE      
!Reorder the extrems of the profile to satisfy: x_pl(2)>x_pl(1)
      IF (x_pl(2)<x_pl(1)) THEN
       dmy=x_pl(2)
       x_pl(2)=x_pl(1)
       x_pl(1)=dmy
       dmy=y_pl(2)
       y_pl(2)=y_pl(1)
       y_pl(1)=dmy
      ENDIF      
!Profile slope
       a=(y_pl(2)-y_pl(1))/(x_pl(2)-x_pl(1))
!Number of longitudinal nodes in the profile
       I_pr_end=NINT(ABS((x_pl(2)-x_pl(1))/p_m_x))+1
      ENDIF

      IF(ALLOCATED(bod_PG)) DEALLOCATE(bod_PG) 
      IF(ALLOCATED(INDEX_prf)) DEALLOCATE(INDEX_prf)
      IF(ALLOCATED(bod_pr)) DEALLOCATE(bod_pr)
      IF(ALLOCATED(OBS_prf)) DEALLOCATE(OBS_prf)
      IF(ALLOCATED(H_dist)) DEALLOCATE(H_dist)
      IF(ALLOCATED(T_mh_prf)) DEALLOCATE(T_mh_prf)
      IF(ALLOCATED(z_T_1D)) DEALLOCATE(z_T_1D)
      IF(ALLOCATED(E_fil)) DEALLOCATE(E_fil) 
      ALLOCATE (bod_PG(2,2*I_pr_end),INDEX_prf(2,I_pr_end),
     *          bod_pr(I_pr_end,0:N_lay),OBS_prf(I_pr_end,10),
     *          H_dist(I_pr_end),T_mh_prf(I_pr_end),z_T_1D(I_pr_end,NC),
     *          E_fil(I_pr_end))

!Loop for profile points
       DO k=1,I_pr_end
       IF (ABS(y_pl(2)-y_pl(1))>ABS(x_pl(2)-x_pl(1))) THEN
         lat=y_pl(1)+(k-1)*p_m_y
         lon=x_pl(1)+a*(lat-y_pl(1))
!Index of matrices in the profile
         INDEX_prf(1,k)=NINT(ABS((lon-lon_min)/p_m_x))+1
         INDEX_prf(2,k)=NINT(ABS((lat-lat_min)/p_m_y))+1
       ELSE
         lon=x_pl(1)+(k-1)*p_m_x
         lat=y_pl(1)+a*(lon-x_pl(1))
!Index of matrices in the profile
         INDEX_prf(1,k)=NINT(ABS((lon-lon_min)/p_m_x))+1
         INDEX_prf(2,k)=NINT(ABS((lat-lat_min)/p_m_y))+1
       ENDIF

!Check borders
        IF (INDEX_prf(1,k)>IDIM) INDEX_prf(1,k)=IDIM
        IF (INDEX_prf(1,k)<1) INDEX_prf(1,k)=1
        IF (INDEX_prf(2,k)>JDIM) INDEX_prf(2,k)=JDIM
        IF (INDEX_prf(2,k)<1) INDEX_prf(2,k)=1
!Horizontal distance
        IF (k==1) THEN
         dist=0
         H_dist(k)=0
        ELSE
!Dist es la distancia calculada con trigonometría esférica
         CALL trigo(y_pl(1),lat,x_pl(1),lon,delta,dist,acim)
!H_dist es input para cross_sect (los elementos de los extremos
!tienen que tener la misma longitud, si no, hay efectos de borde)
         H_dist(k)=SQRT((lat-y_pl(1))**2+(lon-x_pl(1))**2)*pi*Rd_T/180
!!         dist=dist*cos((lat+y_pl(1))*pi/360)
        ENDIF
        bod_PG(1,I_pr_end+k)=dist
!Layer profiles
        bod_pr(k,0)=-E(INDEX_prf(1,k),INDEX_prf(2,k))*1e-3
        bod_pr(k,1:N_lay)=layers(INDEX_prf(1,k),INDEX_prf(2,k),:) 
        OBS_prf(k,1:4)=OBS_calc(:,INDEX_prf(1,k),INDEX_prf(2,k))
        OBS_prf(k,5:8)=OBS_dat(:,INDEX_prf(1,k),INDEX_prf(2,k))
        OBS_prf(k,9)=OBS(12,INDEX_prf(1,k),INDEX_prf(2,k))
        OBS_prf(k,10)=OBS_calc(5,INDEX_prf(1,k),INDEX_prf(2,k))
        T_mh_prf(k)=T_moho(INDEX_prf(1,k),INDEX_prf(2,k))  
!Generate 1D isoterms
!!      CALL geoterm_1D(bod_pr(k,N_lay-1)*1E3,bod_pr(k,N_lay)*1E3,
!!     * -bod_pr(k,0)*1E3,T_mh)
!!       T_mh_prf(k)=T_mh
         DO i_iso=1,NC
          T_iso=T0_iso+(i_iso-1)*dT_iso
          z_T_1D(k,i_iso)=(T_iso-T_mh_prf(k))*(bod_pr(k,N_lay)-
     *                     bod_pr(k,N_lay-1))/
     *                    (T_a-T_mh_prf(k))+bod_pr(k,N_lay-1)
         ENDDO 
!Generate 1D isoterms
       ENDDO
       bod_PG(1,1:I_pr_end)=bod_PG(1,2*I_pr_end:I_pr_end+1:-1)
!Smooth layers
!!       bod_pr(:,0)=-OBS_prf(:,4)*1e-3
!!       bod_pr(:,0)=-OBS_prf(:,8)*1e-3
       DO i=11,N_lay 
        CALL smooth1D(bod_pr(:,i),I_pr_end)
       ENDDO
!Filtered elevation
       E_fil=OBS_prf(:,8)
!!!       CALL smooth1D(E_fil(:),I_pr_end)


      IF(LINUX) THEN
         IPRF=PGOPEN('/XWINDOW')
      ELSE
         IPRF=PGOPEN('/WX')
      ENDIF
      IF (IPRF.LE.0) STOP
      CALL pgscr(0,1.,1.,1.)
      CALL pgscr(1,0.,0.,0.)
      CALL pgask(.false.)
      CALL PGPAP(8.0,.8)
!Title label
       IF (y_pl(2)>y_pl(1)) THEN
        IF (x_pl(2)>x_pl(1)) THEN
         ort_lbl='SW-->NE'
        ELSE
         ort_lbl='SE-->NW'
        ENDIF
       ELSE
        IF (x_pl(2)>x_pl(1)) THEN
         ort_lbl='NW-->SE'
        ELSE
         ort_lbl='NE-->SW'
        ENDIF
       ENDIF
       IF (ABS(x_pl(2)-x_pl(1))<.001) ort_lbl='S-->N'
       IF (ABS(y_pl(2)-y_pl(1))<.001) ort_lbl='W-->E'


       r1=.07
       r2=.91 
       call pgsclp(0)
       CALL PGSVP(r1,r2,0.,1.)
       CALL PGSWIN(0.,1.,0.,1.)
       CALL PGPTEXT(.5,.97,0.,0.5,TRIM(ort_lbl))
!Call to cross_sect**********************
       IF (act=='K') THEN
        CALL cross_sect
       ENDIF
!Call to cross_sect**********************

!/////////////Plot Observable profiles////////////////////////
      DO i=1,9   
       IF (i==1.OR.i==5) THEN
        txt_lbl='Geoid(m)'
        r4=.88 
        r3=.78
       ELSEIF(i==2.OR.i==6) THEN  
        txt_lbl='FA(mGal)'
        r4=.76
        r3=.66
       ELSEIF(i==3.OR.i==7) THEN
        txt_lbl='Bouguer(mGal)'
        r4=.64
        r3=.54
       ELSEIF(i==4.OR.i==8) THEN
        txt_lbl='Topo(m)'
        r4=.52
        r3=.42
       ELSEIF(i==9) THEN
        txt_lbl='HF(mW/m2)'
        r4=.98
        r3=.9
       ENDIF    
       CALL PGSVP(r1,r2,r3,r4)
       IF(i<5) THEN
        IF (act=='K'.AND.i==1) THEN
         z_obs_max(i)=MAX(MAXVAL(OBS_prf(:,i)),MAXVAL(OBS_prf(:,i+4)),
     *                    MAXVAL(GEOID_2D(:)),MAXVAL(GEOID_1D(:)))
         z_obs_min(i)=MIN(MINVAL(OBS_prf(:,i)),MINVAL(OBS_prf(:,i+4)),
     *                    MINVAL(GEOID_2D(:)),MINVAL(GEOID_1D(:)))
        ELSEIF (act=='K'.AND.i==2) THEN        
         z_obs_max(i)=MAX(MAXVAL(OBS_prf(:,2)),MAXVAL(OBS_prf(:,6)),
     *                    MAXVAL(FA_2D(:)))
         z_obs_min(i)=MIN(MINVAL(OBS_prf(:,2)),MINVAL(OBS_prf(:,6)),
     *                    MINVAL(FA_2D(:)))
        ELSEIF (act=='K'.AND.i==3) THEN
          z_obs_max(i)=MAX(MAXVAL(OBS_prf(:,i)),MAXVAL(OBS_prf(:,i+4)),
     *                    MAXVAL(BGA_2D(:)))
          z_obs_min(i)=MIN(MINVAL(OBS_prf(:,i)),MINVAL(OBS_prf(:,i+4)),
     *                    MINVAL(BGA_2D(:)))
        ELSEIF (act=='K'.AND.i==4) THEN
         z_obs_max(i)=MAX(MAXVAL(OBS_prf(:,i)),MAXVAL(OBS_prf(:,i+4)),
     *                    MAXVAL(E_2D_v(:)))
         z_obs_min(i)=MIN(MINVAL(OBS_prf(:,i)),MINVAL(OBS_prf(:,i+4)),
     *                    MINVAL(E_2D_v(:)))
        ELSE 
         z_obs_max(i)=MAX(MAXVAL(OBS_prf(:,i)),MAXVAL(OBS_prf(:,i+4)))
         z_obs_min(i)=MIN(MINVAL(OBS_prf(:,i)),MINVAL(OBS_prf(:,i+4)))
        ENDIF
        margn=(z_obs_max(i)-z_obs_min(i))/10
        IF (z_obs_max(i)==z_obs_min(i)) THEN
         IF (z_obs_max(i)==0D0) THEN
          z_obs_min(i)=1
          z_obs_max(i)=-1
         ELSE
          margn=(z_obs_max(5))/10
         ENDIF
        ENDIF
        z_obs_max(i)=z_obs_max(i)+margn
        z_obs_min(i)=z_obs_min(i)-margn
        CALL PGSCH(.8)
        CALL PGSWIN(0.,1.,0.,1.)
        CALL PGPTEXT(-.07,.5,90.,0.5,TRIM(txt_lbl))
        CALL PGSWIN (0.,MAXVAL(bod_PG(1,:)),z_obs_min(i),z_obs_max(i))
        CALL PGSCH(.7)
        CALL PGBOX('ab',0.0,0,'bctivmn',0.0,0)
        CALL PGLAB('','' ,'')
        CALL PGSLS(1)
       ELSEIF (i<9) THEN
        CALL PGSWIN (0.,MAXVAL(bod_PG(1,:)),z_obs_min(i-4),
     *               z_obs_max(i-4))
        CALL PGSLS(4) 
       ELSE
!For heat flow
        z_obs_max(5)=MAXVAL(OBS_prf(:,i))
        z_obs_min(5)=MINVAL(OBS_prf(:,i))
        margn=(z_obs_max(5)-z_obs_min(5))/10
        IF (z_obs_max(5)==z_obs_min(5))margn=(z_obs_max(5))/10
        z_obs_max(5)=z_obs_max(5)+margn
        z_obs_min(5)=z_obs_min(5)-margn
        CALL PGSLS(1)
        CALL PGSCH(.8)
        CALL PGSWIN(0.,1.,0.,1.)
        CALL PGPTEXT(-.07,.5,90.,0.5,TRIM(txt_lbl))
        CALL PGSWIN (0.,MAXVAL(bod_PG(1,:)),z_obs_min(5),z_obs_max(5))
        CALL PGSCH(.7)
        CALL PGBOX('ab',0.0,0,'bctivmn',0.0,0)
        CALL PGLAB('','' ,'')
        CALL PGSLS(1)
       ENDIF  
!Write observable files for post script%%%%%%%%%%%%%%%%%
       IF (prf_act=='P') THEN
        IF (i==1) name_obs='geoid_calc.xy'
        IF (i==2) name_obs='FA_calc.xy'
        IF (i==3) name_obs='Boug_calc.xy'
        IF (i==4) name_obs='elevation_calc.xy'
        IF (i==5) name_obs='geoid_dat.xy'
        IF (i==6) name_obs='FA_dat.xy'
        IF (i==7) name_obs='Boug_dat.xy'
        IF (i==8) name_obs='elevation_dat.xy' 
        IF (i==9) name_obs='HF_calc.xy'
        open(50,file=TRIM(name_obs),status='UNKNOWN')
        DO ll=1,I_pr_end
         IF (act=='K'.AND.i<5) THEN
          IF (i==1) THEN
           write(50,*)bod_PG(1,I_pr_end+ll),OBS_prf(ll,i), 
     *                GEOID_2D(ll),GEOID_1D(ll)
          ELSEIF (i==2) THEN
           write(50,*)bod_PG(1,I_pr_end+ll),OBS_prf(ll,i),FA_2D(ll)
          ELSEIF (i==3) THEN
           write(50,*)bod_PG(1,I_pr_end+ll),OBS_prf(ll,i),BGA_2D(ll) 
          ELSEIF (i==4) THEN
           write(50,*)bod_PG(1,I_pr_end+ll),OBS_prf(ll,i),E_2D_v(ll)
          ENDIF
         ELSE
          IF (i==4) THEN
           write(50,*)bod_PG(1,I_pr_end+ll),OBS_prf(ll,i),OBS_prf(ll,10)
          ELSE
           write(50,*)bod_PG(1,I_pr_end+ll),OBS_prf(ll,i) 
          ENDIF

         ENDIF
        ENDDO
        close(50)
       ENDIF
!Write observable files for post script%%%%%%%%%%%%%%%%%
      
       CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             OBS_prf(I_pr_end:1:-1,i))
       IF (i==4) THEN
!Plot observed elevation
        CALL PGSCI(4)
        CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             E_fil(I_pr_end:1:-1))
!Plot flexural isostasy elevation
        CALL PGSCR (70,.0,.9,.0)
        CALL PGSCI (70) 
        CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             OBS_prf(I_pr_end:1:-1,10))
        CALL PGSCI(1)
       ENDIF
!Write 2D/1D values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IF (act=='K'.AND.i==1) THEN
         CALL PGSCI(3)
         CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             GEOID_2D(I_pr_end:1:-1))
         CALL PGSCI(4)
         CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             GEOID_1D(I_pr_end:1:-1)) 
         CALL PGSCI(1)
        ELSEIF (act=='K'.AND.i==2) THEN
         CALL PGSCI(3)
         CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             FA_2D(I_pr_end:1:-1))
         CALL PGSCI(1)
        ELSEIF (act=='K'.AND.i==3) THEN
         CALL PGSCI(3)
         CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             BGA_2D(I_pr_end:1:-1))
         CALL PGSCI(1)
        ELSEIF (act=='K'.AND.i==4) THEN
         CALL PGSCI(3)
         CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             E_2D_v(I_pr_end:1:-1))
         CALL PGSCI(1)
        ENDIF      
!Write 2D/1D values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ENDDO   
!Write a file containing the min/max values for Post Script plot
      IF (prf_act=='P') THEN
       open(51,file='obs_scale.dat',status='UNKNOWN')
       WRITE(51,'(10F9.2)')z_obs_min(:),z_obs_max(:) 
       close(51)
      ENDIF
      CALL PGSLS(1)
!/////////////Plot Observable profile////////////////////////

!@@@@@@@@@@@@@@Plot transect@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      z_prf_lb='C'
      IN_prf_ly=0
      n_prf_max=N_lay+1
      DO
!Erase transect
      CALL PGSVP(0.,1.,.07,.4)
      CALL PGSCI(0)
      CALL PGSWIN(0.,1.,0.,1.)
      CALL PGRECT(-.1,1.1,-.1,1.)
      CALL PGSCI(1)
!!      CALL PGSCH(1.)
      CALL PGSCH(.7) 
      CALL PGPTEXT(.2,-.18,0.,0.5,'CHANGE SCALE: c (DEF=NO) ')
      CALL PGSVP(r1,r2,.07,.4)

      IF (z_prf_lb=='C') THEN
        n_prf_max=n_prf_max-1
        IF (n_prf_max<0) n_prf_max=N_lay
        z_prf_max=-(MAXVAL(bod_pr(:,n_prf_max))+10) 
      ELSE
       EXIT 
      ENDIF
      CALL PGSLW(1)
      CALL PGSWIN (0.,MAXVAL(bod_PG(1,:)),z_prf_max,5.)
      CALL PGSCH(.7)
      CALL PGBOX('bntsi',0.0,0,'bctsivmn',0.0,0)
      CALL PGLAB('HORIZONTAL DISTANCE (km) ', 'DEPTH','') 
!Water layer
      CALL PGSCR (50,.0,.97,.97)
      bod_PG(2,1:I_pr_end)=0.
      bod_PG(2,I_pr_end+1:2*I_pr_end)=-8.
      CALL PGSCI(50)
      CALL PGPOLY(2*I_pr_end,bod_PG(1,:),bod_PG(2,:))
!Body layers
      open(16,file='color.rgb',status='UNKNOWN')
       DO i=0,N_lay-1
       read(16,*) R,G,B
       R=R/255
       G=G/255
       B=B/255 
       CALL PGSCR (15+i,R,G,B)
       CALL PGSCI (15+i)
       bod_PG(2,I_pr_end:1:-1)=-bod_pr(:,i)
       bod_PG(2,I_pr_end+1:2*I_pr_end)=-bod_pr(:,i+1)
!Write body files for post script%%%%%%%%%%%%%%%%%
      IF (prf_act=='P') THEN
       IF (i>8) THEN
        WRITE(dmch,'(I2)')i+1
       ELSE
        WRITE(dmch,'(I1)')i+1
       ENDIF
       name_bod='bod'//TRIM(dmch)//'.xy'
       open(50,file=name_bod,status='UNKNOWN')
       DO ll=1,2*I_pr_end
        write(50,*)bod_PG(1,ll),bod_PG(2,ll)
       ENDDO
       close(50)
      ENDIF
!Write body files for post script%%%%%%%%%%%%%%%%%
       CALL PGPOLY(2*I_pr_end,bod_PG(1,:),bod_PG(2,:))
       ENDDO
      close (16)
!Asthenospheric layer
      CALL PGSCR (52,.97,.02,.95)
      bod_PG(2,1:I_pr_end)=-bod_pr(I_pr_end:1:-1,N_lay) 
      bod_PG(2,I_pr_end+1:2*I_pr_end)=-400.
      CALL PGSCI(52)
      CALL PGPOLY(2*I_pr_end,bod_PG(1,:),bod_PG(2,:))
!Draw line between bodies
      CALL PGSLW(3)
      CALL PGSCI(1)
      DO i=0,N_lay
       CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             -bod_pr(I_pr_end:1:-1,i))
      ENDDO
!IF rho_m(T) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!CONDICIÓN TEMPORAL, CAMBIAR PARA HACER CORRER EL INDEX DE PROFILES
!EN Y CUANDO PROCEDA (Y2-Y1)>(X2-X1)
       IF (prop(N_lay)%rho<0.and.I_pr_end>35) THEN
!Plot lithospheric mantle density labels
       i_cont=0
       DO i_lb=25,I_pr_end-10,NINT(I_pr_end/9.)
       i_cont=i_cont+1
       rho_m_av=NINT(rho_a*(1+prop(N_lay)%alpha*(T_a-T_mh_prf(i_lb))/2))
       write(T_txt,*)rho_m_av,'kg/m3'
       z_lb=-(bod_pr(i_lb,N_lay)+bod_pr(i_lb,N_lay-1))/2
       CALL PGQTXT(bod_PG(1,I_pr_end-i_lb),z_lb,
     *              0.,0.5,TRIM(T_txt),xbox,ybox)
       CALL PGSCI(14)
       CALL PGRECT(xbox(1),xbox(3),ybox(1),ybox(2))
       CALL PGSCI(1)
       CALL PGPTEXT(bod_PG(1,I_pr_end-i_lb),z_lb,
     *             0.,0.5,TRIM(T_txt))
!Write lithospheric mantle density file for post script%%%%%%%%%%%%%%%%%
       IF (prf_act=='P') THEN
        IF (i_cont>9) THEN
         WRITE(dmch,'(I2)')i_cont
        ELSE
         WRITE(dmch,'(I1)')i_cont
        ENDIF
        name_bod='rho_m_av'//TRIM(dmch)//'.box'
        open(50,file='rho_m_av.lbl',status='UNKNOWN')  
        open(52,file=name_bod,status='UNKNOWN')
        write(50,'(2F12.2,1X,3I2,2A)') (xbox(1)+xbox(3))/2,
     *           (ybox(1)+ybox(3))/2,8,0,0,' CM',TRIM(T_txt)
        write(52,*) xbox(1),ybox(1)
        write(52,*) xbox(2),ybox(2)
        write(52,*) xbox(3),ybox(3)
        write(52,*) xbox(4),ybox(4)
        close(52)
       ENDIF
!Write lithospheric mantle density file for post script%%%%%%%%%%%%%%%%%
       ENDDO
!Plot lithospheric mantle density labels

!Plot isoterms +++++++++++++++++++++++++++++++++++++++++++++       
       DO i_iso=1,NC
       CALL PGSLS(4)
       CALL PGSLW(5) 
       CALL PGLINE(I_pr_end,bod_PG(1,1:I_pr_end),
     *             -z_T_1D(I_pr_end:1:-1,i_iso))  
       write(T_txt,*) INT(T0_iso+(i_iso-1)*dT_iso),'ºC'
       CALL PGSLS(1)
       CALL PGSLW(3)
       CALL PGQTXT(bod_PG(1,I_pr_end-i_iso*5),-z_T_1D(i_iso*5,i_iso)-5,
     *              0.,0.5,TRIM(T_txt),xbox,ybox)
        CALL PGSCI(0)
        CALL PGRECT(xbox(1),xbox(3),ybox(1),ybox(2))
        CALL PGSCI(1)
       CALL PGPTEXT(bod_PG(1,I_pr_end-i_iso*5),-z_T_1D(i_iso*5,i_iso)-5,
     *             0.,0.5,TRIM(T_txt))
!Write isoterms files for post script%%%%%%%%%%%%%%%%%
       IF (prf_act=='P') THEN
        close(50)
        IF (INT(T0_iso+(i_iso-1)*dT_iso)<1000) THEN
         write(T_txt,'(I3)') INT(T0_iso+(i_iso-1)*dT_iso)
        ELSE
         write(T_txt,'(I4)') INT(T0_iso+(i_iso-1)*dT_iso)
        ENDIF
        name_bod='T_'//TRIM(T_txt)//'.xy'
        open(50,file=name_bod,status='UNKNOWN')
        DO ll=1,I_pr_end
         write(50,*) bod_PG(1,I_pr_end-ll+1),-z_T_1D(ll,i_iso)
        ENDDO
        close(50)
        name_bod='T_'//TRIM(T_txt)//'.lbl'
        write(T_txt,*) INT(T0_iso+(i_iso-1)*dT_iso),'ºC'
        open(50,file=name_bod,status='UNKNOWN')
        write(50,'(2F12.2,1X,3I2,2A)') (xbox(1)+xbox(3))/2,
     *           (ybox(1)+ybox(3))/2,8,0,0,' CM',TRIM(T_txt)
        write(50,*) xbox(1),ybox(1)
        write(50,*) xbox(2),ybox(2)
        write(50,*) xbox(3),ybox(3)
        write(50,*) xbox(4),ybox(4)
        close(50)
       ENDIF
!Write isoterms files for post script%%%%%%%%%%%%%%%%%

       ENDDO 
!Plot isoterms +++++++++++++++++++++++++++++++++++++++++++++
      ENDIF
!IF rho_m(T) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      CALL pgband(0,1,xdm,ydm,xdm,ydm,z_prf_lb)
      z_prf_lb=UCASE(z_prf_lb)
      ENDDO
!@@@@@@@@@@@@@@Plot transect@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      IF (prf_act=='P') THEN
       open(7,file='qq',status='UNKNOWN')
       write(7,5) x_pl(1),y_pl(1),x_pl(2),y_pl(2),N_prf,N_lay,
     *            z_prf_max,MAXVAL(bod_PG(1,:))
       close (7) 
       CALL SYSTEM ('profile.job')
      ENDIF
   5     FORMAT (4F7.2,I3,1X,I3,2F12.2)
  
      CALL PGSLW(3)
      CALL PGSLS(1)
      CALL PGSCH(1.)
!!      CALL pgband(0,1,xdm,ydm,xdm,ydm,dmch)
      CALL PGCLOS 
      

      END SUBROUTINE 
