!SUB_sublitho_an.f


      SUBROUTINE sublitho_an
      USE M_imag
      USE M_pol
      USE M_visualize
      USE M_costa
      USE M_sublitho_an
      USE M_label, ONLY:LBL_SGN
      CHARACTER(1) dmch
      CHARACTER(50) xy1_txt,xy2_txt
      REAL(4) TR_SLC(6),xy1_min,xy1_max,xy2_min,xy2_max
      REAL(4), ALLOCATABLE:: FLD(:,:)
      LOGICAL STAT
      INTEGER PGOPEN
      IF (fst_call) THEN
       INQUIRE(FILE='LITMOD3D.info',EXIST=STAT)
       IF (STAT) THEN
        OPEN(15,file='LITMOD3D.info',STATUS='OLD')
       ELSE
        WRITE(*,*) 'FILE LITMOD3D.info NOT FOUND!'
        RETURN
       ENDIF

     
      READ(15,*)L_x,L_y,Topo_max,Nx,Ny,Nz,dz
      CLOSE(15)
      dz=dz*1E3
      Topo_max=Topo_max*1E3
       NN=MAX(Nx,Ny,Nz)
       INQUIRE(FILE='./OUT_xyz/Trhopvels.xyz',EXIST=STAT)
       IF (.NOT.STAT) THEN
        WRITE(*,*)'FILE OUT_xyz/Trhopvels.xyz NOT FOUND!'
        RETURN
       ENDIF
       OPEN(30,file='./OUT_xyz/Trhopvels.xyz',STATUS='OLD')

      ALLOCATE (T(Nx,Ny,Nz),rho(Nx,Ny,Nz),press(Nx,Ny,Nz),
     *          vp(Nx,Ny,Nz),vs(Nx,Ny,Nz))
      WRITE(*,*)'LOADING T, RHO, P AND ANHARMONIC SEISMIC VELOCITIES...'
       DO j=1,Ny
        DO i=1,Nx
         DO k=1,Nz
          READ(30,*)dm,dm,dm,T(i,j,k),rho(i,j,k),press(i,j,k),
     *              vp(i,j,k),vs(i,j,k)
         ENDDO
        ENDDO
       ENDDO
       CLOSE(30)
       CLOSE(31)
       fst_call=.FALSE.
      ENDIF
    
      
      ALLOCATE(FLD(Nx,Ny))  
      n_pol_lst=0 
      OPEN(400,file='sublith_an.xyz',status='UNKNOWN')
!Start writing at the end of the file
       DO
        READ(400,*,IOSTAT=i_status)
        IF (i_status/=0) EXIT
       ENDDO      
!Input D_T (K) and litho
!       WRITE(*,*)'D_T (ºC), litho number of anomalous sublith body...?'
!       READ(*,*)D_T_an,litho_an
      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -14)
       WRITE(400,'(A1,1X,F7.2,1X,I3)') '*',D_T_an,litho_an
!       CALL pgband(0,1,xdm,ydm,xdm,ydm,dmch)
!       IF (ICHAR(dmch)==27) RETURN
       CALL lay_sel
       IF (LBL_SGN==1) RETURN
       write(*,*)z_an_up,z_an_dw
       IF(z_an_up>z_an_dw) THEN
        dm=z_an_up
        z_an_up=z_an_dw
        z_an_dw=dm
        ENDIF
        Npr_v_up=NINT((Topo_max+z_an_up)/dz+1)
        Npr_v_dw=NINT((Topo_max+z_an_dw)/dz+1)
        IF(Npr_v_up<1)Npr_v_up=1
        IF(Npr_v_up>Nz)Npr_v_up=Nz
        IF(Npr_v_dw<1)Npr_v_dw=1
        IF(Npr_v_dw>Nz)Npr_v_dw=Nz

!        write(*,*)'z_an_up,z_an_dw',z_an_up,z_an_dw
!        write(*,*)'Npr_v_up,Npr_v_dw',Npr_v_up,Npr_v_dw
       i_z_pas=INT((Npr_v_dw-Npr_v_up)/10)
      DO I_z_an=Npr_v_up,Npr_v_dw,1 
       z_sublit_an=Topo_max-(I_z_an-1)*dz
       WRITE(title_slc_an,15) -z_sublit_an*1E-3
 15    FORMAT('Slice at depth', F6.1,'km')
       DO ih=Ny,1,-1
        DO il=1,Nx
         FLD(il,Ny-ih+1)=T(il,ih,I_z_an)
        ENDDO
       ENDDO 
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

!Open the window
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
        
!Plot the temperatures
        T_an_min=MINVAL(T(:,:,I_z_an))
        T_an_max=MAXVAL(T(:,:,I_z_an))
       CALL PGSWIN(xy1_min,xy1_max,xy2_min,xy2_max)
        CALL PGBOX('bctsin',0.0,0,'bctsivn',0.0,0)
        CALL PGIMAG(FLD(:,:),Nx,Ny,1,Nx,1,Ny,T_an_min,T_an_max,
     *              TR_SLC)
        CALL PGLAB(TRIM(xy1_txt),TRIM(xy2_txt),'Temperature (C)')
        CALL PGWEDG('RI',2.5,3.,T_an_min,T_an_max,'ºC')
!Dibujar los puntos de la linea de costa
       DO ii=1,coast_lin
        CALL PGLINE(npt_coast(ii),L_coast(2*ii-1,1:npt_coast(ii)),
     *              L_coast(2*ii,1:npt_coast(ii)))
       ENDDO
      

      CALL get_pol
      CALL in_out_pol ('U')

      CALL PGSLCT(IPG_SLC)
      CALL PGCLOS

      ENDDO

      CLOSE(400)

      END SUBROUTINE



