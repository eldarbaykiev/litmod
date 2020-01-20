!  SUB_visualize.f


      SUBROUTINE visualize
      USE M_imag
      USE M_costa
      USE M_pol
      USE M_label
      USE M_profile, ONLY:x_pl,y_pl
      USE M_ref_points
      USE M_visualize
      LOGICAL STAT
      REAL(4)dx,dy,dxy1,dxy2,Topo_max,sl_dep,xy1,xy2,
     *       xy1_min,xy1_max,xy2_min,xy2_max
      REAL(4) xdm,ydm,L_x,L_y,TR_SLC(6)
      CHARACTER(1) le,le2,ps
      CHARACTER(50) UCASE
      INTEGER h,l
      INTEGER PGOPEN
      INTEGER,ALLOCATABLE:: INDEX_prf(:,:)
      REAL(4) z_lst,FLD_MIN,FLD_MAX
      REAL(4), ALLOCATABLE:: FLD(:,:,:),dist(:,:)
      CHARACTER(50) FLD_txt(6),WDG_txt(6),xy1_txt,xy2_txt,title_slc
      REAL(4):: a,lon,lat,dmy

      INQUIRE(FILE='LITMOD3D.info',EXIST=STAT)
      IF (STAT) THEN
       OPEN(15,file='LITMOD3D.info',STATUS='OLD')
      ELSE
       WRITE(*,*) 'FILE LITMOD3D.info NOT FOUND!'
       RETURN
      ENDIF

      READ(15,*)L_x,L_y,Topo_max,Nx,Ny,Nz,dz
      CLOSE(15) 
      FLD_txt(1)='Temperature'
      FLD_txt(2)='Density'
      FLD_txt(5)='Pressure'
      FLD_txt(6)='Thermal cond'
      WDG_txt(1)='C'
      WDG_txt(2)='kg/m3'
      WDG_txt(5)='Mpa'
      WDG_txt(6)='(W/m K)'

      L_x=L_x*1E3
      L_y=L_y*1E3
      dx=L_x/(Nx-1)
      dy=L_y/(Ny-1)
      dz=dz*1E3
      Topo_max=Topo_max*1E3
      z_lst=Topo_max-(Nz-1)*dz

      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -12) 
       CALL pgband(0,1,xdm,ydm,xdm,ydm,ps)
       ps=UCASE(ps)


       IF(a_v=='V') THEN
!Seismic velocities
        IF (att_y/='N') THEN
         FLD_txt(3)='Vp velocity (atten)'
         FLD_txt(4)='Vs velocity (atten)'
        ELSE
         FLD_txt(3)='Vp velocity (anharmonic)'
         FLD_txt(4)='Vs velocity (anharmonic)'
        ENDIF
        WDG_txt(3)='km/s'
        WDG_txt(4)='km/s'
        
        IF (fst_call_v.AND.att_y/='N') THEN
!Attenuated velocities
         INQUIRE(FILE='./OUT_xyz/V_attenuation.xyz',EXIST=STAT) 
         IF (.NOT.STAT) THEN
          WRITE(*,*)'FILE OUT_xyz/V_attenuation.xyz NOT FOUND!'
          WRITE(*,*)'RUN atten PROGRAM'
          RETURN
         ENDIF  
         
         OPEN(31,file='OUT_xyz/V_attenuation.xyz',STATUS='OLD')
         WRITE(*,*)'LOADING SEISMIC VELOCITIES WITH ATTENUATION...'
         ALLOCATE (vp_at(Nx,Ny,Nz),vs_at(Nx,Ny,Nz))
          DO j=1,Ny
           DO i=1,Nx
            DO k=1,Nz
             READ(31,*)dm,dm,dm,vp_at(i,j,k),vs_at(i,j,k)
            ENDDO
           ENDDO
          ENDDO
          CLOSE(31)
          fst_call_v=.FALSE.
        ENDIF

       ELSE
!Velocity anomalies
        INQUIRE(FILE='./OUT_xyz/V_anomaly.xyz',EXIST=STAT)
        IF (.NOT.STAT) THEN
         WRITE(*,*)'FILE OUT_xyz/V_anomaly.xyz NOT FOUND!'
         WRITE(*,*)'RUN vel2anomaly PROGRAM'
         RETURN
        ENDIF
         FLD_txt(3)='Vp anomaly'
         FLD_txt(4)='Vs anomaly'
         WDG_txt(3)='%'
         WDG_txt(4)='%'
         IF (fst_call_a) THEN
          OPEN(31,file='OUT_xyz/V_anomaly.xyz',STATUS='OLD')
          WRITE(*,*)'LOADING VELOCITY ANOMALIES...'
          ALLOCATE (a_vp(Nx,Ny,Nz),a_vs(Nx,Ny,Nz))
          DO j=1,Ny
           DO i=1,Nx
            DO k=1,Nz
             READ(31,*)dm,dm,dm,a_vp(i,j,k),a_vs(i,j,k)
            ENDDO
           ENDDO
          ENDDO
          CLOSE(31)
          fst_call_a=.FALSE.
         ENDIF
       ENDIF


      IF (fst_call) THEN
      NN=MAX(Nx,Ny,Nz)
      INQUIRE(FILE='./OUT_xyz/Trhopvels.xyz',EXIST=STAT)
       IF (.NOT.STAT) THEN
        WRITE(*,*)'FILE OUT_xyz/Trhopvels.xyz NOT FOUND!'
        RETURN
       ENDIF
       OPEN(30,file='./OUT_xyz/Trhopvels.xyz',STATUS='OLD')
       OPEN(31,file='./OUT_xyz/Thermal_cond.xyz',STATUS='OLD')
       
      ALLOCATE (T(Nx,Ny,Nz),rho(Nx,Ny,Nz),press(Nx,Ny,Nz),
     *          vp(Nx,Ny,Nz),vs(Nx,Ny,Nz),K_therm(Nx,Ny,Nz)) 
      WRITE(*,*)'LOADING T, RHO, P, K AND ELASTIC SEISMIC VELOCITIES...'
       DO j=1,Ny
        DO i=1,Nx
         DO k=1,Nz
          READ(30,*)dm,dm,dm,T(i,j,k),rho(i,j,k),press(i,j,k),
     *              vp(i,j,k),vs(i,j,k)
           READ(31,*)dm,dm,dm,K_therm(i,j,k)
         ENDDO
        ENDDO
       ENDDO
       CLOSE(30)
       CLOSE(31)
       fst_call=.FALSE.
      ENDIF 
 

      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -11)

      xdm=(lon_max+lon_min)/2
      ydm=(lat_max+lat_min)/2
      CALL pgband(0,1,xdm,ydm,xdm,ydm,le)
      le=UCASE(le) 
      IF (ICHAR(le)==27) RETURN
      IF (le=='A') THEN
       x1=xdm
       y1=ydm
       CALL pgband(1,1,xdm,ydm,xdm,ydm,le2)
       IF (ICHAR(le2)==27) RETURN
       x2=xdm
       y2=ydm
       IF (x1<lon_min) x1=lon_min
       IF (x2<lon_min) x2=lon_min
       IF (x1>lon_max) x1=lon_max
       IF (x2>lon_max) x2=lon_max

       IF (y1<lat_min) y1=lat_min
       IF (y2<lat_min) y2=lat_min
       IF (y1>lat_max) y1=lat_max
       IF (y2>lat_max) y2=lat_max
       x_pl(1)=x1
       x_pl(2)=x2
       y_pl(1)=y1
       y_pl(2)=y2 
      ENDIF 

      IF (le=='A') THEN 
!Vertical profile
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
       Npr_v=NINT(ABS((y_pl(2)-y_pl(1))/p_m_y))+1
       IF (Npr_v>Ny) Npr_v=Ny
       IF (Npr_v<1) Npr_v=1
       IF(ALLOCATED(INDEX_prf)) DEALLOCATE(INDEX_prf)
       IF(ALLOCATED(dist)) DEALLOCATE(dist)
       ALLOCATE(INDEX_prf(2,Npr_v),dist(3,Npr_v))      
!Loop for profile points
       DO k=1,Npr_v
        lat=y_pl(1)+(k-1)*p_m_y
        lon=x_pl(1)+a*(lat-y_pl(1))
!Index of matrices in the profile
        INDEX_prf(1,k)=NINT(ABS((lon-lon_min)/p_m_x))+1
        INDEX_prf(2,k)=NINT(ABS((lat-lat_min)/p_m_y))+1
        dist(1,k)=lon
        dist(2,k)=lat
!Check borders
        IF (INDEX_prf(1,k)>IDIM) INDEX_prf(1,k)=IDIM
        IF (INDEX_prf(1,k)<1) INDEX_prf(1,k)=1
        IF (INDEX_prf(2,k)>JDIM) INDEX_prf(2,k)=JDIM
        IF (INDEX_prf(2,k)<1) INDEX_prf(2,k)=1
       ENDDO
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
       Npr_v=NINT(ABS((x_pl(2)-x_pl(1))/p_m_x))+1
       IF (Npr_v>Nx) Npr_v=Nx
       IF (Npr_v<1) Npr_v=1
       IF(ALLOCATED(INDEX_prf)) DEALLOCATE(INDEX_prf)
       IF(ALLOCATED(dist)) DEALLOCATE(dist)
       ALLOCATE(INDEX_prf(2,Npr_v),dist(3,Npr_v))
!Loop for profile points
       DO k=1,Npr_v 
        lon=x_pl(1)+(k-1)*p_m_x
        lat=y_pl(1)+a*(lon-x_pl(1))
!Index of matrices in the profile
        INDEX_prf(1,k)=NINT(ABS((lon-lon_min)/p_m_x))+1
        INDEX_prf(2,k)=NINT(ABS((lat-lat_min)/p_m_y))+1
        dist(1,k)=lon
        dist(2,k)=lat
!Check borders
       IF (INDEX_prf(1,k)>IDIM) INDEX_prf(1,k)=IDIM
       IF (INDEX_prf(1,k)<1) INDEX_prf(1,k)=1
       IF (INDEX_prf(2,k)>JDIM) INDEX_prf(2,k)=JDIM
       IF (INDEX_prf(2,k)<1) INDEX_prf(2,k)=1
       ENDDO 
      ENDIF
!Horizontal distance
       DO k=1,Npr_v
        IF (k==1) THEN
         dist(3,k)=0
        ELSE
         CALL trigo(y_pl(1),dist(2,k),x_pl(1),dist(1,k),
     *              delta,dist(3,k),acim)        
        ENDIF
       ENDDO     
       i_or=1
       NN1=Npr_v
       NN2=Nz
       dxy2=dz
!Define the transformation matrix
       TR_SLC(1)=0
       TR_SLC(4)=Topo_max*1E-3
       TR_SLC(2)=dist(3,Npr_v)/Npr_v
       TR_SLC(6)=(Topo_max-z_lst)*1E-3/Nz
       TR_SLC(3)=0
       TR_SLC(5)=0       
       xy1_min=0
       xy1_max=dist(3,Npr_v)
       xy2_min=-z_lst*1E-3
       xy2_max=Topo_max*1E-3
       xy1_txt='Horizontal distance (km)'
       xy2_txt='Deph (km)'
      ELSE
       i_or=3
       CALL lay_sel
       Npr_v=NINT((Topo_max+z_slice)*Nz/((Nz-1)*dz))
       write(*,*)'Npr_v,dz',Npr_v,dz
       IF (Npr_v>Nz) Npr_v=Nz
       IF (Npr_v<1) Npr_v=1
       NN1=Nx
       NN2=Ny
       dxy1=dx
       dxy2=dy
       k=Npr_v 
       sl_dep=(((Npr_v-1)*dz)-Topo_max) *1E-3
       IF (sl_dep>400) sl_dep=400
       WRITE(title_slc,15) sl_dep
 15    FORMAT('Slice at depth', F6.1,'km')
       write(*,*)'Slice at depth: ',sl_dep, ' km'
!Define the transformation matrix
       TR_SLC(1)=lon_min
       TR_SLC(4)=lat_min
       TR_SLC(2)=(lon_max-lon_min)/IDIM
       TR_SLC(6)=(lat_max-lat_min)/JDIM
       TR_SLC(3)=0
       TR_SLC(5)=0
       xy1_min=lon_min
       xy1_max=lon_max
       xy2_min=lat_min
       xy2_max=lat_max
       xy1_txt='Latitude '
       xy2_txt='Longitude'       
      ENDIF 
      OPEN(30,file='./OUT_xyz/T.sl',STATUS='UNKNOWN')
      OPEN(31,file='./OUT_xyz/dens.sl',STATUS='UNKNOWN') 
      OPEN(32,file='./OUT_xyz/Vp.sl',STATUS='UNKNOWN')
      OPEN(33,file='./OUT_xyz/Vs.sl',STATUS='UNKNOWN')
 
      ALLOCATE(FLD(NN1,NN2,5)) 

      DO h=1,NN1
        IF (i_or==1) THEN
         xy1=dist(3,h)*1D3
         i=INDEX_prf(1,h)
         j=JDIM-INDEX_prf(2,h)+1
        ELSE
         xy1=(h-1)*dxy1
         i=h
        ENDIF

       DO l=1,NN2
        IF (i_or==1) THEN
         k=l
         xy2=Topo_max-(l-1)*dz
        ELSE
         j=Ny-l+1
         xy2=(l-1)*dxy2
        ENDIF
        
        FLD(h,l,1)=T(i,j,k)
        FLD(h,l,2)=rho(i,j,k)
        FLD(h,l,5)=press(i,j,k)
        WRITE(30,*)xy1*1E-3,xy2*1E-3,T(i,j,k)
        WRITE(31,*)xy1*1E-3,xy2*1E-3,rho(i,j,k)
        IF(a_v/='V') THEN
         FLD(h,l,3)=a_vp(i,j,k)
         FLD(h,l,4)=a_vs(i,j,k)
         WRITE(32,*)xy1*1E-3,xy2*1E-3,a_vp(i,j,k)
         WRITE(33,*)xy1*1E-3,xy2*1E-3,a_vs(i,j,k)
        ELSE
         IF(att_y=='N') THEN
          FLD(h,l,3)=vp(i,j,k)
          FLD(h,l,4)=vs(i,j,k)
          WRITE(32,*)xy1*1E-3,xy2*1E-3,vp(i,j,k)
          WRITE(33,*)xy1*1E-3,xy2*1E-3,vs(i,j,k)
         ELSE
          FLD(h,l,3)=vp_at(i,j,k)
          FLD(h,l,4)=vs_at(i,j,k)
          WRITE(32,*)xy1*1E-3,xy2*1E-3,vp_at(i,j,k)
          WRITE(33,*)xy1*1E-3,xy2*1E-3,vs_at(i,j,k)
         ENDIF
        ENDIF

       ENDDO
      ENDDO 
      CLOSE(30)
      CLOSE(31)
      CLOSE(32)
      CLOSE(33) 

!Display profile location
      CALL PGSLCT(ISTAT)
       IF(i_or<3) THEN
        CALL PGSLW(10)
        CALL PGLINE(2,x_pl,y_pl)
      ENDIF
      IF (ps=='Y') THEN
!Call the script to plot the slices
       OPEN(15,file='qq',STATUS='UNKNOWN')
       WRITE(15,10)NINT(MINVAL(T)),NINT(MAXVAL(T)),NINT(MINVAL(rho)),
     *             NINT(MAXVAL(rho)),MINVAL(vp),MAXVAL(vp),MINVAL(vs),
     *             MAXVAL(vs),xy1*1E-3,xy2*1E-3,dxy1*1E-3,dxy2*1E-3,i_or
 10    FORMAT(4I6,4F7.3,2F9.0,1X,2F7.3,I2)
       CLOSE(15)
       CALL SYSTEM ('visualize.job')  
      ELSE
       i=1
       DO
        IF(LINUX) THEN
         IPG_SLC=PGOPEN('/XWINDOW')
        ELSE
         IPG_SLC=PGOPEN('/WX')
        ENDIF
        IF (IPG_SLC<=0) STOP
        CALL PALETT(2, CONTRA, BRIGHT)
        CALL pgscr(0,1.,1.,1.)
        CALL pgscr(1,0.,0.,0.)
        CALL pgask(.false.)
        CALL PGPAP(8.0,.8)
        CALL PGSCI(1)
        CALL PGSVP(.3,.9,.1,.9)
        CALL text_plot(IPG_SLC,.3,.9,.1,.9,lon_min,lon_max,lat_min,
     *                lat_max,-13)
        IF (i>5)i=1
         IF (i_or<3.AND.i==3.AND.a_v=='V') THEN
          FLD_MIN=6
          FLD_MAX=MAXVAL(FLD(:,:,i)) 
         ELSEIF (i_or<3.AND.i==4.AND.a_v=='V') THEN 
          FLD_MIN=4
          FLD_MAX=MAXVAL(FLD(:,:,i))
         ELSEIF (i_or<3.AND.i==3.AND.a_v/='V') THEN
          FLD_MIN=-4
          FLD_MAX=4
         ELSEIF (i_or<3.AND.i==4.AND.a_v/='V') THEN
          FLD_MIN=-6
          FLD_MAX=6
         ELSEIF (i_or<3.AND.i==2) THEN
          FLD_MIN=2000
          FLD_MAX=MAXVAL(FLD(:,:,i))
         ELSE
          FLD_MIN=MINVAL(FLD(:,:,i))
          FLD_MAX=MAXVAL(FLD(:,:,i))
         ENDIF
        CALL PGSWIN(xy1_min,xy1_max,xy2_min,xy2_max)
        CALL PGBOX('bctsin',0.0,0,'bctsivn',0.0,0)
        CALL PGIMAG(FLD(:,:,i),NN1,NN2,1,NN1,1,NN2,
     *              FLD_MIN,FLD_MAX,TR_SLC)
        CALL PGLAB(TRIM(xy1_txt),TRIM(xy2_txt),TRIM(FLD_txt(i)))
        CALL PGWEDG('RI',2.5,3.,FLD_MIN,FLD_MAX,TRIM(WDG_txt(i)))
      IF (i_or==3) THEN
!Dibujar los puntos de la linea de costa
       DO ii=1,coast_lin
        CALL PGLINE(npt_coast(ii),L_coast(2*ii-1,1:npt_coast(ii)),
     *              L_coast(2*ii,1:npt_coast(ii)))
       ENDDO
!Plot labels
       DO jj=1,NMARK
        CALL PGPT1(labels(jj)%long,labels(jj)%lat,18)
       ENDDO
!Plot reference points
       CALL PGSFS(3)
       DO j=1,N_RPTS
       CALL PGPOLY(N_pol_RPTS(j),pol_RPTS(2*j-1,:),pol_RPTS(2*j,:))
       CALL PGLINE(N_pol_RPTS(j)+1,pol_RPTS(2*j-1,:),pol_RPTS(2*j,:))
       ENDDO
       CALL PGSFS(1)
       CALL PGSVP(.01,.3,.1,.9)
       CALL PGSWIN(0.,1.,0.,1.)
       CALL PGPTEXT(.01,.1,0.,0.0,TRIM(title_slc))

      ENDIF

         CALL pgband(0,1,xdm,ydm,xdm,ydm,ps)
         IF (ICHAR(ps)==27) EXIT
         i=i+1 
         CALL PGSLCT(IPG_SLC)
         CALL PGCLOS
       ENDDO
        CALL PGSLCT(IPG_SLC)
        CALL PGCLOS
      ENDIF
      CALL PGSLCT(ISTAT)
      CALL clear_scr(0)

      END SUBROUTINE
