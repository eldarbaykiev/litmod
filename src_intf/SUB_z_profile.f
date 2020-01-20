! Subroutine  SUB_z_profile.f


      SUBROUTINE z_profile
      USE M_imag
      USE M_visualize
      USE M_profile, ONLY: layers,N_lay
      REAL(4) xdm,ydm
      LOGICAL STAT
      CHARACTER(1)le,ps
      CHARACTER(50) UCASE
      INTEGER i_zp,j_zp
      INTEGER PGOPEN
      REAL(4) z_lst,FLD_MIN,FLD_MAX,L_x,L_y,Topo_max
      REAL(4), ALLOCATABLE:: Z_KR(:),FLD(:,:)
      CHARACTER(50) FLD_txt(6)
       


      INQUIRE(FILE='LITMOD3D.info',EXIST=STAT)
      IF (STAT) THEN
       OPEN(15,file='LITMOD3D.info',STATUS='OLD')
      ELSE
       WRITE(*,*) 'FILE LITMOD3D.info NOT FOUND!'
       RETURN
      ENDIF

      READ(15,*)L_x,L_y,Topo_max,Nx,Ny,Nz,dz
      CLOSE(15)

      FLD_txt(1)='Temperature (C)'
      FLD_txt(2)='Density (kg/m3)'
      FLD_txt(5)='Pressure (Mpa)'
      FLD_txt(6)='Thermal cond (W/m K)'

      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -12)
      CALL pgband(0,1,xdm,ydm,xdm,ydm,ps)
      ps=UCASE(ps)

      IF(a_v=='V') THEN
!Seismic velocities
       IF (att_y/='N') THEN
        FLD_txt(3)='Vp velocity (atten) (km/s)'
        FLD_txt(4)='Vs velocity (atten) (km/s)'
       ELSE
        FLD_txt(3)='Vp velocity (elastic) (km/s)'
        FLD_txt(4)='Vs velocity (elastic) (km/s)'
       ENDIF
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
       FLD_txt(3)='Vp anomaly (%)'
       FLD_txt(4)='Vs anomaly (%)'
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

      INQUIRE(FILE='./OUT_xyz/Trhopvels.xyz',EXIST=STAT)
      IF (.NOT.STAT) THEN
       WRITE(*,*)'FILE OUT_xyz/Trhopvels.xyz NOT FOUND!'
       RETURN
      ENDIF
      INQUIRE(FILE='./OUT_xyz/Thermal_cond.xyz',EXIST=STAT)
      IF (.NOT.STAT) THEN
       WRITE(*,*)'FILE OUT_xyz/Thermal_cond.xyz NOT FOUND!'
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
     *                -3)

       xdm=(lon_max+lon_min)/2
       ydm=(lat_max+lat_min)/2
       CALL pgband(0,1,xdm,ydm,xdm,ydm,le)
       IF (ICHAR(le)==27) RETURN

       i_zp=NINT((xdm-lon_min)*Nx/(lon_max-lon_min))
       j_zp=Ny-NINT((ydm-lat_min)*Ny/(lat_max-lat_min))
       IF (i_zp<1) i_zp=1
       IF (i_zp>Nx) i_zp=Nx 
       IF (j_zp<1) j_zp=1
       IF (j_zp>Ny) j_zp=Ny
       ALLOCATE(Z_KR(Nz),FLD(Nz,6))

       OPEN(30,file='./OUT_xyz/T.z',STATUS='UNKNOWN')
       OPEN(31,file='./OUT_xyz/dens.z',STATUS='UNKNOWN')
       OPEN(32,file='./OUT_xyz/Vp.z',STATUS='UNKNOWN')
       OPEN(33,file='./OUT_xyz/Vs.z',STATUS='UNKNOWN')
       OPEN(34,file='./OUT_xyz/press.z',STATUS='UNKNOWN')
       OPEN(35,file='./OUT_xyz/Thermal_cond.z',STATUS='UNKNOWN')
       DO k=1,Nz
        z=Topo_max-(k-1)*dz
        Z_KR(k)=z
        FLD(k,1)=T(i_zp,j_zp,k)
        FLD(k,2)=rho(i_zp,j_zp,k)
        FLD(k,5)=press(i_zp,j_zp,k)
        FLD(k,6)=K_therm(i_zp,j_zp,k)
        WRITE(30,*) z,T(i_zp,j_zp,k)
        WRITE(31,*) z,rho(i_zp,j_zp,k)
        WRITE(34,*) z,press(i_zp,j_zp,k)        
        WRITE(35,*) z,K_therm(i_zp,j_zp,k)
        
        IF(a_v/='V') THEN  
         FLD(k,3)=a_vp(i_zp,j_zp,k)
         FLD(k,4)=a_vs(i_zp,j_zp,k) 
         WRITE(32,*) z,a_vp(i_zp,j_zp,k)
         WRITE(33,*) z,a_vs(i_zp,j_zp,k)
        ELSE
         IF(att_y=='N') THEN
          FLD(k,3)=vp(i_zp,j_zp,k)
          FLD(k,4)=vs(i_zp,j_zp,k)
          WRITE(32,*) z,vp(i_zp,j_zp,k)
          WRITE(33,*) z,vs(i_zp,j_zp,k)
         ELSE
          FLD(k,3)=vp_at(i_zp,j_zp,k)
          FLD(k,4)=vs_at(i_zp,j_zp,k)
          WRITE(32,*) z,vp_at(i_zp,j_zp,k)
          WRITE(33,*) z,vs_at(i_zp,j_zp,k)
         ENDIF
        ENDIF
       
       ENDDO
       CLOSE(30)
       CLOSE(31)
       CLOSE(32)
       CLOSE(33)
       CLOSE(34)
       CLOSE(35)

       IF (ps=='Y') THEN
        CALL SYSTEM ('z_profile.job')
       ELSE
        i=1
        z_lst=Topo_max-(Nz-1)*dz
        DO
        IF(LINUX) THEN
         IPG_Z=PGOPEN('/XWINDOW')
        ELSE
         IPG_Z=PGOPEN('/WX')
        ENDIF
        IF (IPG_Z<=0) STOP
        CALL pgscr(0,1.,1.,1.)
        CALL pgscr(1,0.,0.,0.)
        CALL pgask(.false.)
        CALL PGPAP(8.0,.8)
        CALL pgsvp(.3,.9,.1,.9)
       CALL text_plot(IPG_Z,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -13)
        
         IF (i>6)i=1 
         IF (i==3.AND.a_v=='V') THEN
          FLD_MIN=6
          FLD_MAX=MAXVAL(FLD(:,i))
         ELSEIF (i==4.AND.a_v=='V') THEN
          FLD_MIN=4
          FLD_MAX=MAXVAL(FLD(:,i))
         ELSEIF (i==3.AND.a_v/='V') THEN
          FLD_MIN=-4
          FLD_MAX=4
         ELSEIF (i==4.AND.a_v/='V') THEN
          FLD_MIN=-6
          FLD_MAX=6
         ELSEIF (i==2) THEN
          FLD_MIN=2000
          FLD_MAX=MAXVAL(FLD(:,i))
         ELSEIF (i==6) THEN
          z_lst=-layers(i_zp,j_zp,N_lay) 
          FLD_MIN=1.5
          FLD_MAX=MAXVAL(FLD(:,i))+.5
!        write(*,*)'z_lst',z_lst,N_lay,layers(i_zp,j_zp,N_lay),Topo_max 
         ELSE
          z_lst=Topo_max-(Nz-1)*dz 
          FLD_MIN=MINVAL(FLD(:,i))
          FLD_MAX=MAXVAL(FLD(:,i))
         ENDIF
         CALL PGSWIN (FLD_MIN,FLD_MAX,z_lst,5.)
         CALL PGBOX('bctsinm',0.0,0,'bctsivm',0.0,0)
         CALL PGLAB('', 'DEPTH (km)',TRIM(FLD_txt(i)))
         CALL PGLINE(Nz,FLD(:,i),Z_KR(:))
         CALL pgband(0,1,xdm,ydm,xdm,ydm,ps)
         IF (ICHAR(ps)==27) EXIT
         i=i+1
         CALL PGSLCT(IPG_Z)
         CALL PGCLOS
        ENDDO
        CALL PGSLCT(IPG_Z)
        CALL PGCLOS
       ENDIF



        END SUBROUTINE
