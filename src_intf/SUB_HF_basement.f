!! SUB_HF_basement.f


      SUBROUTINE HF_basement
      USE M_imag
      USE M_material
      USE M_visualize 
      USE M_profile, ONLY: i_lay,layers
      USE M_costa
      USE M_ref_points
      USE M_label
      LOGICAL STAT
      INTEGER PGOPEN,IPG_HF,K_BASp
      REAL(4) xdm,ydm,L_x,L_y,TR_SLC(6)
      REAL(4) xy1_min,xy1_max,xy2_min,xy2_max
      REAL(4) k_thT,k_thB,dz_BAS
      CHARACTER(50) FLD_txt(5),WDG_txt(5),xy1_txt,xy2_txt,title_slc
      CHARACTER(1) le,le2,ps

      INQUIRE(FILE='LITMOD3D.info',EXIST=STAT)
      IF (STAT) THEN
       OPEN(15,file='LITMOD3D.info',STATUS='OLD')
      ELSE
       WRITE(*,*) 'FILE LITMOD3D.info NOT FOUND!'
       RETURN
      ENDIF
      READ(15,*)L_x,L_y,Topo_max,Nx,Ny,Nz,dz
      CLOSE(15)
      FLD_txt(1)='Basal Heat Flow'
      WDG_txt(1)='mW/m2'

!Load Temperature file if needed
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
       
! Select basement layer to compute the HF at
       CALL lay_sel
!Compute HF at the selected horizon
        k_thT=prop(i_lay)%K
        k_thB=prop(i_lay+1)%K
        IF(k_thT==0) k_thT=3.2
        IF(k_thB==0) k_thB=3.2
!        a_h1=prop(i_top)%A
!        a_h2=prop(i_bot)%A
!        FLD(:,:,i)
!        layers(ii,jj,i_lay)
        write(*,*)'i_top,i_bot,i_lay',i_top,i_bot,i_lay 
        write(*,*)'k_thT,k_thB',k_thT,k_thB
! Open file to write basal HF output 
        OPEN(30,file='./OUT_xyz/basal_HF.xyz',STATUS='UNKNOWN')
        WRITE(30,*)'Basal heat flow between layers: ',i_lay,i_lay+1
        WRITE(30,*) TRIM(prop(i_lay)%name_mat),' - ',
     *              TRIM(prop(i_lay+1)%name_mat)
        DO jj=JDIM,1,-1
         DO ii=1,IDIM
         K_BAS=NINT((Topo_max+layers(ii,jj,i_lay))/dz+1)
         IF (T(ii,jj,K_BAS)==T(ii,jj,K_BAS+1)) THEN
          K_BASp=K_BAS+2
          dz_BAS=2*dz
         ELSE
          K_BASp=K_BAS+1
          dz_BAS=dz 
         ENDIF
         HF(ii,jj)=-2E3*(T(ii,jj,K_BAS)-T(ii,jj,K_BASp))*k_thT*k_thB/
     *            ((k_thT+k_thB)*dz_BAS*1E3)
!        write(*,*)'HF',HF(ii,jj)
!        write(*,*) 'T',T(ii,jj,K_BAS),T(ii,jj,K_BAS+1),dz,k_thT,k_thB
         WRITE(30,*) coord(1:2,ii,jj),HF(ii,jj)
        ENDDO
       ENDDO
       CLOSE(30)
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
       FLD_MIN=10
       FLD_MAX=100
       write(*,*)minval(HF), maxval(HF)
       FLD_MIN=minval(HF)
       FLD_MAX=maxval(HF)
!Plot HF
       IF(LINUX) THEN
        IPG_HF=PGOPEN('/XWINDOW')
       ELSE
        IPG_HF=PGOPEN('/WX')
       ENDIF
       IF (IPG_HF<=0) STOP
       CALL PALETT(2, CONTRA, BRIGHT)
       CALL pgscr(0,1.,1.,1.)
       CALL pgscr(1,0.,0.,0.)
       CALL pgask(.false.)
       CALL PGPAP(8.0,.8)
       CALL PGSCI(1)
       CALL PGSVP(.3,.9,.1,.9)
       CALL text_plot(IPG_HF,.3,.9,.1,.9,lon_min,lon_max,lat_min,
     *                lat_max,-18)
       CALL PGSWIN(xy1_min,xy1_max,xy2_min,xy2_max)    
       CALL PGBOX('bctsin',0.0,0,'bctsivn',0.0,0)
       CALL PGIMAG(HF(:,:),Nx,Ny,1,Nx,1,Ny,
     *              FLD_MIN,FLD_MAX,TR_SLC)
       CALL PGLAB(TRIM(xy1_txt),TRIM(xy2_txt),TRIM(FLD_txt(1)))
       CALL PGWEDG('RI',2.5,3.,FLD_MIN,FLD_MAX,TRIM(WDG_txt(1)))
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
!       CALL PGPTEXT(.01,.1,0.,0.0,TRIM(title_slc))
       CALL pgband(0,1,xdm,ydm,xdm,ydm,ps)
       CALL PGCLOS  
       CALL PGSLCT(ISTAT)
       CALL clear_scr(0)


      END SUBROUTINE

 
