!Para compilar programas que usan las librerias pgplot:


!g77 [].f -o [] -L/home/jfullea/pgplot -lpgplot -L/usr/X11R6/lib -lX11
!ifc GEO3Dmod_INTF.f SUB* -o GEO3Dmod_INTF -L/home/jfullea/pgplot -lpgplot -L/usr/X11R6/lib -lX11 -lg2c -Vaxlib -C -w
! For compilation:
!ifort LITMOD3D_INTF.f SUB* -watch source -o LITMOD3D_INTF  -L/home/jfullea/pgplot -lpgplot -L/usr/lib -lX11 -w
! gfortran LITMOD3D_INTF.f SUB* -o LITMOD3D_INTF -lpgplot -L/usr/lib -lX11 -fd-lines-as-comments -lpng
!ifort LITMOD3D_INTF.f SUB* -watch source -o LITMOD3D_INTF -L/home/jfullea/pgplot -lpgplot -L/usr/lib -lX11 -w -lgfortran
!PGOPEN must always be declared as INTEGER
      MODULE M_imag
      INTEGER:: IDIM,JDIM,ISTAT,MODE,IPS,IPR,IPG_LAY,IS_LAY,i_geo,i_obs
      INTEGER:: I_LEG,N_cr
      INTEGER:: i_top,i_bot,ins,I_BODY,i_up_h,i_dw_h
      REAL(4):: lon_min,lon_max,lat_min,lat_max,E_max,E_min,p_m_x,
     *          p_m_y,p_m,TR(6)
      REAL(4):: z_min,z_max,CONTRA,BRIGHT,z_min_pl,z_max_pl
      REAL(4):: rho_max,rho_min
      REAL(4), ALLOCATABLE:: E(:,:),OBS(:,:,:),OBS_calc(:,:,:),
     *                       OBS_dat(:,:,:),bod_tk(:,:),IMG(:,:),
     *                       T_moho(:,:),HF(:,:),coord(:,:,:)
      CHARACTER(LEN=*), PARAMETER::num='0123456789'
      CHARACTER(LEN=14):: string
      CHARACTER(50) :: fil,obs_file(12),IMG_txt
      CHARACTER(2):: act
      CHARACTER(1):: cr_type,hi_ch
      LOGICAL LINUX,out_win,fst_call,fst_call_a,fst_call_v
      END MODULE M_imag

      MODULE M_costa
      REAL(4), ALLOCATABLE:: L_coast(:,:)
      INTEGER num_geo,coast_max,coast_lin
      INTEGER,ALLOCATABLE:: npt_coast(:)
      REAL(4) x_coast,y_coast
      CHARACTER(50) ch_coast_x,ch_coast_y
      END MODULE M_costa

      MODULE M_profile
      REAL(4) d_pr,xy_pr_min,xy_pr_max,x_pl(2),y_pl(2),
     *        xy_pr_p_min,xy_pr_p_max,z_max_mdfr,z_min_mdfr
      REAL(4) L,L_p,T_iso,T0_iso,dT_iso
      REAL(4), ALLOCATABLE:: xy_pr(:),xy_pr_p(:),E_pr(:),bod_pr(:,:),
     *                      OBS_pr(:),H_dist(:),OBS_prf(:,:),z_T_1D(:,:)
     *                      ,T_mh_prf(:),bod_pr_bak(:)
      REAL(4), ALLOCATABLE:: bod_sf(:,:,:),layers(:,:,:),dm_lay(:,:,:)
      REAL(4), ALLOCATABLE:: FA_2D(:),GEOID_2D(:),BGA_2D(:),GEOID_1D(:),
     *                       E_2D_v(:)
      INTEGER N_lay,ifl,ifl_p,N_pr,i_pr,i_lay,I_p,IPRF,N_prf,I_pr_end,NC
      INTEGER, ALLOCATABLE:: pr_pos(:),GLUE_m(:,:,:,:)
      CHARACTER(1) ort,updown_lay,bod_ok
      CHARACTER(3) MOUSE
      CHARACTER(1) exit_mdreg,abort_mdreg,lems,ms_flow
      INTEGER fw_pr,i_gl,j_gl,ms_pr,i_nonzero
      REAL(4) lnz_min,lnz_max,ltz_min,ltz_max
      CHARACTER(80) TEXT_BDSEL,title
      CHARACTER(1), ALLOCATABLE::i_bod_plot(:)
      CHARACTER(1)prf_act
      LOGICAL glue_l,ch_v_scl
      END MODULE M_profile

      MODULE M_material
      TYPE MATERIAL
!K -->Thermal conductivity
!A -->Heat production rate
!rho-->Density at T=Ts,P=0 (except for materials with
!      rho<0, in which case rho is defined at
!      T=Ta)
!alpha-->Volumetric expansion coefficient
!beta-->Pressure coefficient
!name_mat-->Short material description
!gmma_T-->Thermal Gr�neisen parameter
!K_0p-->Pressure derivative of isothermal bulk modulus
!K_T-->isothermal bulk modulus
!litho-->Index describing the type of mantle material
! litho=0-->crustal material
! litho=11-19-->oceanic mantle
! litho=21-29-->Phenerozoic mantle
! litho=31-39-->Proterozoic mantle
! litho=41-49-->Archean mantle
! litho=91-99-->Sublithospheric mantle
      REAL(4) K,A,rho,alpha,beta,gmma_T,K_0p,K_T
      INTEGER litho
      CHARACTER(LEN=50) name_mat
      END TYPE MATERIAL
      TYPE (MATERIAL) prop (-1:100)
      INTEGER ins_mat
      CHARACTER(1) abort_body
      END MODULE M_material

      MODULE M_pol
      INTEGER n_pol,n_pol_lst
      INTEGER,ALLOCATABLE:: NP_POL(:),NP_POL_lst(:)
      REAL(4),ALLOCATABLE:: A_pol(:,:),A_pol_lst(:,:)
      LOGICAL yin_yang
      END MODULE M_pol

      MODULE M_label
      INTEGER, PARAMETER:: NMARK_MAX=1500
      TYPE ANNOT
      REAL(4):: long,lat,z_lbl
      CHARACTER(80):: text_mark
      END TYPE ANNOT
      TYPE (ANNOT) labels(NMARK_MAX)
      INTEGER:: LBL_SGN,NMARK,NTICK
      CHARACTER(1):: lbl_act
      REAL(4):: bnd_lbl
      END MODULE M_label

      MODULE M_ref_points
      INTEGER, PARAMETER:: N_pol_RPTS_MAX=100,NP_pol_RPTS_MAX=100
      INTEGER:: N_pol_RPTS(N_pol_RPTS_MAX),N_RPTS
      REAL(4):: pol_RPTS(2*N_pol_RPTS_MAX,NP_pol_RPTS_MAX)
      CHARACTER(1):: RPTS_act
      END MODULE M_ref_points

      MODULE M_ctes
      REAL(4) pi,G_u,g,T_a,rho_a,z_bot,Rd_T
      REAL(8) rho_w
      END MODULE M_ctes

      MODULE M_visualize
      REAL(4), ALLOCATABLE:: T(:,:,:),rho(:,:,:),press(:,:,:),
     *                       vp(:,:,:),vs(:,:,:),a_vp(:,:,:),
     *                       a_vs(:,:,:),vp_at(:,:,:),vs_at(:,:,:),
     *                       K_therm(:,:,:)
      REAL(4) z_slice,dz
      INTEGER Nx,Ny,Nz,Npr_v
      CHARACTER(1) a_v,att_y
      END MODULE M_visualize

      MODULE M_sublitho_an
      REAL(4) D_T_an
      INTEGER litho_an,IPG_SLC,I_z_an,Npr_v_up
      REAL(4)z_sublit_an,z_an_up, z_an_dw
      CHARACTER(50) title_slc_an
      END MODULE M_sublitho_an

      USE M_imag
      USE M_costa
      USE M_profile
      USE M_material
      USE M_label
      USE M_ref_points
      USE M_ctes

D     USE DFLIB
D     TYPE (QWINFO) qinfo
C      CHARACTER GETCHARQQ
      INTEGER PGOPEN
      REAL(4) x,y,dm_dx,lon_min_dm,lon_max_dm,lat_min_dm,lat_max_dm,zdd
      CHARACTER le,KEY,GETCHARQQ
      CHARACTER(50) UCASE
      LOGICAL LINUX_CHECK,STAT,STAT1
      integer MAC
      CHARACTER(2) chdm
      CHARACTER(100) CWD
      INTEGER N_xdmm,N_ydmm,syn
D     qinfo%TYPE = QWIN$MAX
D     STAT=SETWSIZEQQ(QWIN$FRAMEWINDOW, qinfo)
D     STAT=SETWSIZEQQ(0, qinfo)


!Palette parameters
      CONTRA=1.3
      BRIGHT=0.5

      rho_max=3600
      rho_min=2100
!Constants
      rho_a=3200.
      z_bot=-400D3
      T_a=1350.
      rho_w=1030.
      pi=3.14159
      G_u=6.67e-11
      g=9.8
      Rd_T=6370.



      LINUX=LINUX_CHECK()
      IF(LINUX) THEN
        MAC=MAC_CHECK()
        if(MAC.eq.1) then
          write(*,*)'OS: MACOS'
          string='./layers/layer'
        else
          write(*,*)'OS: LINUX'
          string='./layers/layer'
        end if

      ELSE
       write(*,*)'OS: WINDOWS'
       string='.\layers\layer'
      ENDIF

!Count the number of bodies,N_lay
!/////////////////////////////////////////////////
      open(15,file='layers.info',status='UNKNOWN')
      jfl=0
      N_lay=0
      read(15,*)
      DO
       read(15,*,IOSTAT=jfl)j,prop(j)%name_mat,prop(j)%rho,
     *                      prop(j)%alpha,prop(j)%K,
     *                      prop(j)%A,prop(j)%beta,prop(j)%gmma_T,
     *                      prop(j)%K_0p,prop(j)%K_T,prop(j)%litho
       IF (jfl/=0) EXIT
       IF(j>1) THEN
        IF(prop(j)%litho/=0.AND.prop(j-1)%litho==0) then
          N_cr=j-1
        ENDIF
       ENDIF
       N_lay=N_lay+1
      ENDDO
      close (15)
!Delete empty files
       CALL SYSTEM ('find -type f -size 0 -print0 | xargs -0 rm ')
!Call to the script which generates the coast line and resamples
!the  elevation and layer files
      INQUIRE(FILE='./GEO_DATA/elev_dat.xyz',EXIST=STAT)
      IF (.NOT.STAT) THEN
       WRITE(*,*)'./GEO_DATA/elev_dat.xyz DOES NOT EXIST'
       WRITE(*,*)'CALL LITMOD_3D.job SETTING pre_pro=1'
       STOP
      ENDIF
!Call to the script LITMOD_3D.job to obtain the coordinates of the region
!and the number of nodes
      OPEN(1,file='dmmy',status='UNKNOWN')
      WRITE(1,*) 1
      CLOSE(1)
      if(MAC.eq.1) then
        CALL SYSTEM('./LITMOD_3D.job')
      else
        CALL SYSTEM('./LITMOD_3D.job')
      end if
!      CALL SYSTEM('./LITMOD_3D_mac.job')
      OPEN(1,file='info.dat',status='OLD')
      READ(1,*) lon_min,lon_max,lat_min,lat_max,IDIM,JDIM,syn
      CLOSE (1)
!Check if the number of nodes or region coordinates have changed since last
!run
      OPEN(200,file='./layers/layer1.xyz',status='UNKNOWN')
      y_a=0
      i_con=0
      N_xdmm=0
      DO
       read(200,*,IOSTAT=jfl) x,y,zdd
       IF (jfl/=0) EXIT
       IF (i_con>0.AND.N_xdmm==0.AND.(y-y_a)/=0) N_xdmm=i_con
       y_a=y
       i_con=i_con+1
      ENDDO
      CLOSE(200)
      N_ydmm=i_con/N_xdmm

      !current way to call minmax (gmt independent)
      CALL SYSTEM('python3 minmax.py ./layers/layer1.xyz>qq')

      OPEN(1,file='qq',status='OLD')
      READ(1,*)lon_min_dm,lon_max_dm,lat_min_dm,lat_max_dm
      WRITE(*,*)lon_min_dm,lon_max_dm,lat_min_dm,lat_max_dm
      CLOSE(1)
      CALL SYSTEM('rm qq')
      IF(IDIM/=N_xdmm.OR.JDIM/=N_ydmm.OR.ABS(lon_min_dm-lon_min)>1D-3
     *  .OR.ABS(lon_max_dm-lon_max)>1D-3.OR.ABS(lat_min_dm-lat_min)>
     *   1D-3.OR.ABS(lat_max_dm-lat_max)>1D-3) THEN
        WRITE(*,*)"NUMBER OF NODES OR REGION COORDINATES CHANGED
     * SINCE LAST RUN"
        WRITE(*,*)"RESAMPLING LAYER AND ELEVATION FILES "
!Rewrite the layer files
      IF(((lon_min<lon_min_dm.AND.lon_max<lon_min_dm).OR.
     *   lon_min>lon_max_dm.AND.lon_max>lon_max_dm).OR.
     *   ((lat_min<lat_min_dm.AND.lat_max<lat_min_dm).OR.
     *   lat_min>lat_max_dm.AND.lat_max>lat_max_dm)) THEN
       DO i=1,N_lay
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        IF (i>9) THEN
        WRITE(chdm,'(I2)')i
        ELSE
        WRITE(chdm,'(I1)')i
        ENDIF
        fil=TRIM(string)//TRIM(chdm)//'.xyz'
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       open(20,file=fil,status='OLD')
        DO jj=1,JDIM
         c_lat=lat_min+(lat_max-lat_min)*jj/(JDIM-1)
         DO ii=1,IDIM
          c_lon=lon_min+(lon_max-lon_min)*ii/(IDIM-1)
         WRITE(20,*) c_lon,c_lat,150*i/(N_lay+1)
         ENDDO
        ENDDO
        close(20)
       ENDDO
      ENDIF
        IF(LINUX) THEN
!Call to the script which generates the coast line and resamples
!the  elevation and layer files
            CALL SYSTEM ('./elev.job')
         ELSE
            CALL SYSTEM ('elev.bat')
         ENDIF
      ENDIF

      INQUIRE(FILE='./GEO_DATA/coast_prev.xy',EXIST=STAT1)

      IF(.NOT. STAT1.AND.syn==0) THEN
         PRINT '('' Coast file does not exist.''/
     *           '' Do you want to create it? (N=def, Y): '',$)'
!         KEY=GETCHARQQ()
         READ(*,*)KEY
         IF(ICHAR(KEY) .EQ. 13) THEN
            KEY='N'
         ELSE
            KEY=UCASE(KEY)
         ENDIF
      ELSE
         KEY='N'
      ENDIF

      IF(KEY .EQ. 'Y') THEN
        STAT1=.TRUE.
         IF((LINUX).OR.(MAC.eq.1)) THEN
!Call to the script which generates the coast line and resamples
!the  elevation and layer files
            CALL SYSTEM ('./elev.job')
         ELSE
            CALL SYSTEM ('elev.bat')
         ENDIF
      ENDIF

      !write(*,*)"Line 317"
!Abre el device for graphic
!Abrir dos devices: fichero postscript y una xwin. Para pasar de uno
!a otro a de usarse el comando PGSLCT(device)
!      IPS=PGOPEN('interfaz.ps/vcps')
!      ISTAT=PGOPEN('?')
      IF((LINUX).OR.(MAC.eq.1)) THEN
         ISTAT=PGOPEN('/XWINDOW')
      ELSE
         ISTAT=PGOPEN('/WX')
      ENDIF
      IF (ISTAT<=0) THEN
         PRINT '('' Cannot open graphics window, program stops'')'
         STOP
      ENDIF

!Matrix dimensions
      p_m_x=(lon_max-lon_min)/(IDIM-1)
      p_m_y=(lat_max-lat_min)/(JDIM-1)
!!        p_m=(lon_max-lon_min)/IDIM
      ALLOCATE (E(IDIM,JDIM),OBS(12,IDIM,JDIM),OBS_calc(5,IDIM,JDIM),
     *          OBS_dat(4,IDIM,JDIM),IMG(IDIM,JDIM),T_moho(IDIM,JDIM),
     *          HF(IDIM,JDIM),coord(2,IDIM,JDIM) )
!The transformation matrix TR is used to calculate the world
!coordinates of the center of the "cell" that represents each
!array element. The world coordinates of the center of the cell
!corresponding to array element A(I,J) are given by:
!                X = TR(1) + TR(2)*I + TR(3)*J
!                Y = TR(4) + TR(5)*I + TR(6)*J

!Usually TR(3) and TR(5) are zero -- unless the coordinate
!transformation involves a rotation or shear.
      TR(1)=lon_min
      TR(4)=lat_min
      TR(2)=(lon_max-lon_min)/IDIM
      TR(6)=(lat_max-lat_min)/JDIM
      TR(3)=0
      TR(5)=0
!Establece el color del fondo (0) y de los ejes, etiquetas etc. (1)
!Utiliza escala Red Green Blue normalizada [0.,1.]
      call pgscr(0,1.,1.,1.)
      call pgscr(1,0.,0.,0.)
      call pgscr(17,.5,.5,.5)


!LOAD COAST LINE***************************************
      IF (STAT1.AND.syn==0) THEN
      open(200,file='./GEO_DATA/coast_prev.xy',status='UNKNOWN')
!Almacena vectorialmente las coordenadas de los puntos de la
!l�nea de costa en los vectores L_coast(:,:) y npt_coast(:)
      open(1,file='./GEO_DATA/coast.xy',status='UNKNOWN')
      DO i=1,2
       read(200,*)
      ENDDO
      jfl=0
      coast_lin=1
      j_coast=0
      num_geo=1
      coast_max=0

      DO
      read(200,*,IOSTAT=jfl) ch_coast_x,ch_coast_y
      IF (jfl/=0) EXIT
      IF ((ICHAR(ch_coast_x(1:1))<48.OR.ICHAR(ch_coast_x(1:1))>=58)
     *     .AND.ICHAR(ch_coast_x(1:1))/=45) THEN
       IF (j_coast>coast_max) coast_max=j_coast
       coast_lin=coast_lin+1
       j_coast=0
       write(1,*)'500 500'
      ELSE
       write(1,*)TRIM(ch_coast_x),' ',TRIM(ch_coast_y)
       j_coast=j_coast+1
      ENDIF
      num_geo=num_geo+1
      ENDDO

      write(1,*)'500 500'
      REWIND (1)
      IF(coast_lin==1) coast_max=j_coast

      write(*,*)2*coast_lin,coast_max
      ALLOCATE (L_coast(2*coast_lin,coast_max),npt_coast(coast_lin))
      i_coast=1
      j_coast=1
!      write(*,*)num_geo
      DO i=1,num_geo
        read(1,*) x_coast,y_coast
        if (x_coast/=500) then
!La siguiente l�nea asegura que la longitud se continua
!en el meridiano de Greenwich
          if (x_coast>180) x_coast=x_coast-360

!          write(*,*)2*i_coast-1,j_coast
!          write(*,*)2*coast_lin,coast_max
!ERROR HERE
!ELDAR: quick fix
          if ((2*i_coast-1).le.(2*coast_lin)) then
            if ((2*i_coast).le.(2*coast_lin)) then
              if (j_coast.le.coast_max) then
                L_coast(2*i_coast-1,j_coast)=x_coast
                L_coast(2*i_coast,j_coast)=y_coast
              end if
            end if
          end if
          j_coast=j_coast+1
        else
          npt_coast(i_coast)=j_coast-1
          i_coast=i_coast+1
          j_coast=1
        endif
      ENDDO
      close(1)
      close(200)

      ENDIF
!LOAD COAST LINE***************************************
!LOAD LABEL FILE+++++++++++++++++++++++++++++++++++++++
       INQUIRE(FILE='labels.dat',EXIST=STAT)
       IF (STAT) THEN
        OPEN (911,file='labels.dat',status='OLD')
        jfl=0
        NMARK=0
        DO
         read(911,*,IOSTAT=jfl) labels(NMARK+1)%long,
     *   labels(NMARK+1)%lat,labels(NMARK+1)%z_lbl,
     *   labels(NMARK+1)%text_mark
         IF (jfl/=0) EXIT
         NMARK=NMARK+1
        ENDDO
       ENDIF
!LOAD LABEL FILE+++++++++++++++++++++++++++++++++++++++

!LOAD REFERENCE POINTS FILE&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!Reference points are suposed to be grouped in one or more
!polygons
!ref_points.info contains one register for each polygon
!with its number of vertex. This info is loaded in
!N_pol_RPTS (INTEGER)
!ref_points.dat contains the x, y coordinates of all
!polygon vertexs ordered according to ref_points.info.
!This info is loaded in pol_RPTS (REAL)

      INQUIRE(FILE='ref_points.dat',EXIST=STAT)
      INQUIRE(FILE='ref_points.info',EXIST=STAT1)
      STAT=STAT.AND.STAT1
      IF (STAT) THEN
       OPEN (40,file='ref_points.info',status='OLD')
        jfl=0
        N_RPTS=0
        DO
         READ(40,*,IOSTAT=jfl) N_pol_RPTS(N_RPTS+1)
         IF (jfl/=0) EXIT
         N_RPTS=N_RPTS+1
        ENDDO
       OPEN (45,file='ref_points.dat',status='OLD')
        DO i=1,N_RPTS
         DO j=1,N_pol_RPTS(i)
         READ(45,*)pol_RPTS(2*i-1,j),pol_RPTS(2*i,j)
         ENDDO
         pol_RPTS(2*i-1,N_pol_RPTS(i)+1)=pol_RPTS(2*i-1,1)
         pol_RPTS(2*i,N_pol_RPTS(i)+1)=pol_RPTS(2*i,1)
        ENDDO

      ENDIF
!LOAD REFERENCE POINTS FILE&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!i_geo=1 indicates the first run (i.e. is necessary to load
!the layer and observable files), otherwise i_geo>1
      i_geo=1
!When fst_call=.true. geophysical observables are load in
!memmory in SUB_visualize.f,SUB_sublitho_an.f, and SUB_z_profile.for
!(only needed once in running time)
      fst_call=.true.
      fst_call_v=.true.
      fst_call_a=.true.
      i_obs=1
!      write(*,*)"Line 480"
      IF(LINUX) THEN
         CALL SYSTEM('rm ./bak/layer*')
         CALL SYSTEM('cp ./layers/layer*.xyz layers.info ./bak')
      ELSE
         CALL SYSTEM('copy .\layers\layer*.* .\bak')
      ENDIF
!      write(*,*)"Line 487"
!Starts the endless loop, waiting for options
!################################################################
      DO
!Count the number of bodies,N_lay
!/////////////////////////////////////////////////
      open(15,file='layers.info',status='UNKNOWN')
      jfl=0
      N_lay=0
      read(15,*)
      DO
       read(15,*,IOSTAT=jfl)j,prop(j)%name_mat,prop(j)%rho,
     *                      prop(j)%alpha,prop(j)%K,
     *                      prop(j)%A,prop(j)%beta,prop(j)%gmma_T,
     *                      prop(j)%K_0p,prop(j)%K_T,prop(j)%litho
       IF (jfl/=0) EXIT
       N_lay=N_lay+1
      ENDDO
      close (15)

      write(*,'(A30,I2)')'NUMBER OF BODIES DEFINED--> ',N_lay

!ONLY FOR THE FIRST RUN
      IF (i_geo==1) THEN
      i_geo=i_geo+1
      ALLOCATE (layers(IDIM,JDIM,N_lay),dm_lay(IDIM,JDIM,N_lay))
      ALLOCATE (GLUE_m(IDIM,JDIM,N_lay,N_lay))
!ELEVATION
       open(2,file='./GEO_DATA/elev_dat.xyz',status='OLD')
       DO j=JDIM,1,-1
        DO i=1,IDIM
         read(2,*) coord(1,i,j),coord(2,i,j),E(i,j)
        ENDDO
       ENDDO
       close(2)
       E_min=minval(E)-200
       E_max=maxval(E)+200
!CALCULATED OBSERVABLES
      INQUIRE(FILE='./GEO_DATA/geoid_calc.xyz',EXIST=STAT)
      INQUIRE(FILE='./GEO_DATA/FA_calc.xyz',EXIST=STAT1)
      STAT=STAT.AND.STAT1
      INQUIRE(FILE='./GEO_DATA/Boug_calc.xyz',EXIST=STAT1)
      STAT=STAT.AND.STAT1
      INQUIRE(FILE='./GEO_DATA/elev_calc.xyz',EXIST=STAT1)
      STAT=STAT.AND.STAT1
      IF (STAT) THEN
      DO i_obs=1,4
       IF (i_obs==1) open(30,file='./GEO_DATA/geoid_calc.xyz',
     *                    status='OLD')
       IF (i_obs==2) open(30,file='./GEO_DATA/FA_calc.xyz',
     *                    status='OLD')
       IF (i_obs==3) open(30,file='./GEO_DATA/Boug_calc.xyz',
     *                    status='OLD')
       IF (i_obs==4) open(30,file='./GEO_DATA/elev_calc.xyz',
     *                    status='OLD')
       DO j=JDIM,1,-1
        DO i=1,IDIM
         read(30,*) dum,dum,OBS_calc(i_obs,i,j)
        ENDDO
       ENDDO
       close(30)
      ENDDO
      ELSE
       OBS_calc=0.
      ENDIF
!Regional (flexural) isostasy
      INQUIRE(FILE='./GEO_DATA/elev_calc_flex.xyz',EXIST=STAT)
      IF (STAT) THEN
       open(31,file='./GEO_DATA/elev_calc_flex.xyz',status='OLD')
        DO j=JDIM,1,-1
         DO i=1,IDIM
          read(31,*) dum,dum,OBS_calc(5,i,j)
         ENDDO
        ENDDO
        close(31)
      ENDIF
!Moho temperature
!      INQUIRE(FILE='./GEO_DATA/T_moho.xyz',EXIST=STAT)
!      IF (STAT) THEN
!       open(31,file='./GEO_DATA/T_moho.xyz',status='OLD')
!       DO j=JDIM,1,-1
!         DO i=1,IDIM
!         read(31,*) dum,dum,T_moho(i,j)
!         ENDDO
!        ENDDO
!        close(31)
!       ENDIF
!Heat flow
      INQUIRE(FILE='./GEO_DATA/HF.xyz',EXIST=STAT)
      IF (STAT) THEN
       open(20,file='./GEO_DATA/HF.xyz',status='OLD')
         DO j=JDIM,1,-1
          DO i=1,IDIM
          read(20,*) dum,dum,OBS(12,i,j)
          ENDDO
         ENDDO
         close(20)
       ENDIF
!DATA
      INQUIRE(FILE='./GEO_DATA/geoid_dat.xyz',EXIST=STAT)
      INQUIRE(FILE='./GEO_DATA/FA_dat.xyz',EXIST=STAT1)
      STAT=STAT.AND.STAT1
      INQUIRE(FILE='./GEO_DATA/Boug_dat.xyz',EXIST=STAT1)
      STAT=STAT.AND.STAT1
      INQUIRE(FILE='./GEO_DATA/elev_dat.xyz',EXIST=STAT1)
      STAT=STAT.AND.STAT1
      IF (STAT) THEN
      DO i_obs=1,4
       IF (i_obs==1) open(30,file='./GEO_DATA/geoid_dat.xyz',
     *                    status='OLD')
       IF (i_obs==2) open(30,file='./GEO_DATA/FA_dat.xyz',
     *                    status='OLD')
       IF (i_obs==3) open(30,file='./GEO_DATA/Boug_dat.xyz',
     *                    status='OLD')
       IF (i_obs==4) open(30,file='./GEO_DATA/elev_dat.xyz',
     *                    status='OLD')
       DO j=JDIM,1,-1
        DO i=1,IDIM
         read(30,*) dum,dum,OBS_dat(i_obs,i,j)
        ENDDO
       ENDDO
       close(30)
      ENDDO
      ELSE
      OBS_dat=0.
      OBS_dat(4,:,:)=E(:,:)
      ENDIF
! Observable file names
      obs_file(2)=' Res Geoid'
      obs_file(3)=' Res Free air'
      obs_file(4)=' Res Bouguer'
      obs_file(5)=' Res Elevation'
      obs_file(6)=' Res Elevation (flex)'
      obs_file(7)=' Calc Geoid'
      obs_file(8)=' Calc Free air'
      obs_file(9)=' Calc Bouguer'
      obs_file(10)=' Calc Elevation'
      obs_file(11)=' Calc Elevation (flex)'
      obs_file(12)=' Heat flow'
      OBS(2:5,:,:)= OBS_calc(:,:,:)-OBS_dat(:,:,:)
      OBS(6,:,:)= OBS_calc(5,:,:)-OBS_dat(4,:,:)
      OBS(7:10,:,:)= OBS_calc(:,:,:)
      OBS(11,:,:)= OBS_calc(5,:,:)
!Starts the graphic with elevation
      call pgask(.false.)
      OBS(1,:,:)=E(:,:)
      obs_file(1)=' Elevation'
      z_min=minval(OBS(1,:,:))-200
      z_max=maxval(OBS(1,:,:))+200
      i_obs=1
      IMG=OBS(1,:,:)
      IMG_txt=obs_file(1)
      CALL PALETT(2, CONTRA, BRIGHT)
      CALL PGPAP(8.0,.8)
      CALL clear_scr  (0)

!Layers in memmory
      DO i=1,N_lay
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       IF (i>9) THEN
       WRITE(chdm,'(I2)')i
       ELSE
       WRITE(chdm,'(I1)')i
       ENDIF
       fil=TRIM(string)//TRIM(chdm)//'.xyz'
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      open(20,file=fil,status='OLD')
       DO jj=JDIM,1,-1
        DO ii=1,IDIM
        read(20,*) dumm,dumm,layers(ii,jj,i)
        dm_lay(ii,jj,i)=layers(ii,jj,i)
        ENDDO
       ENDDO
      close(20)
      ENDDO

      ENDIF

!/////////////////////////////////////////////////
!Write actions and obtain act value
      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *               1)
      act=UCASE(act)
      SELECT CASE (act)
      CASE ('P')
      CALL profile
      CASE ('K')
      CALL profile
      CASE ('M')
      CALL modf_reg
      CASE ('A')
      CALL new_lay
      CASE ('J')
      CALL joint_lay
      CASE ('B')
      CALL body_plot
      CASE ('S')
      CALL body_plot_ps
      CASE ('T')
      CALL two2three_lay
      CASE ('I')
      CALL smooth_lay
      CASE ('L')
      CALL label
      CASE ('C')
      CALL ch_obs
      CASE ('R')
      CALL ref_points
      CASE ('H')
      CALL hierarchy
      CASE ('V')
      CALL visualize
      CASE ('Z')
      CALL z_profile
      CASE ('D')
      CALL delete_lay
      CASE ('U')
      CALL sublitho_an
      CASE ('O')
      CALL histogram
      CASE ('F')
      CALL HF_basement
      CASE ('E')
      STOP
      CASE DEFAULT
      write(*,*)'OPTION NOT RECOGNISED, PLEASE, RE-ENTER'
      END SELECT
      ENDDO
!################################################################

!Finaliza el device for graphic
       CALL PGCLOS
       CLOSE(911)
       CLOSE(40)
       CLOSE(45)
      STOP
      END




      SUBROUTINE PALETT(PTYPE, CONTRA, BRIGHT)
C-----------------------------------------------------------------------
C Set a "palette" of colors in the range of color indices used by
C PGIMAG.
C-----------------------------------------------------------------------
      INTEGER PTYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
      REAL TL(10), TR(10), TG(10), TB(10)
C
      DATA GL /0.0, 1.0/
       DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
CC
      DATA TL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.0, 1.7/
      DATA TR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0, 1.0/
      DATA TG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0, 1.0/
      DATA TB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/

      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (PTYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (PTYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (PTYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (PTYPE.EQ.4) THEN
C        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (PTYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      ELSE IF (PTYPE.EQ.6) THEN
C        -- Thickness
         CALL PGCTAB(TL, TR, TG, TB, 10, CONTRA, BRIGHT)
      END IF
      END SUBROUTINE

      SUBROUTINE zoom (ik)
      USE M_imag
      USE M_costa
      USE M_pol
      REAL(4) lnz_min,lnz_max,ltz_min,ltz_max
      REAL(4) x1,y1,x2,y2,xdm,ydm
      INTEGER IXMIN,IXMAX,IYMIN,IYMAX
      CHARACTER le
      BRIGHT = 0.5
      CONTRA  = 1.3
      CALL PALETT(2, CONTRA, BRIGHT)
!Establece el color del fondo (0) y de los ejes, etiquetas etc. (1)
!Utiliza escala Red Green Blue normalizada [0.,1.]
      call pgscr(0,1.,1.,1.)
      call pgscr(1,0.,0.,0.)
      call pgscr(17,.5,.5,.5)
!pgask controls new page prompting
      call pgask(.false.)
      CALL PGSCI(1)
      xdm=(lon_max+lon_min)/2
      ydm=(lat_max+lat_min)/2
      call pgband(0,1,0.,0.,xdm,ydm,le)
      x1=xdm
      y1=ydm
      call pgband(2,1,xdm,ydm,xdm,ydm,le)
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
!      CALL PGENV(lnz_min,lnz_max,ltz_min,ltz_max,1,0)
      call pgsvp(.3,.9,.1,.9)
      call pgswin(lnz_min,lnz_max,ltz_min,ltz_max)
      CALL PGBOX('bctsim',0.0,0,'bctsivm',0.0,0)
       IXMIN=1+nint((lnz_min-lon_min)/p_m_x)
       IXMAX=IXMIN+nint((lnz_max-lnz_min)/p_m_x)
       IYMIN=1+nint((ltz_min-lat_min)/p_m_y)
       IYMAX= IYMIN+nint((ltz_max-ltz_min)/p_m_y)
      CALL PGIMAG(E,IDIM,JDIM,IXMIN,IXMAX,IYMIN,IYMAX,E_MIN,E_MAX,TR)
      CALL PGLAB('Long ', 'Lat', 'ZOOM')
      CALL PGIDEN
!Dibujar los puntos de la linea de costa
       DO ii=1,coast_lin
        CALL PGLINE(npt_coast(ii),L_coast(2*ii-1,1:npt_coast(ii)),
     *              L_coast(2*ii,1:npt_coast(ii)))
      ENDDO
      if (ik==0) return
      DO j=1,ik
      CALL PGPOLY(NP_POL(j),A_pol(2*j-1,:),A_pol(2*j,:))
      CALL PGLINE(NP_POL(j)+1,A_pol(2*j-1,:),A_pol(2*j,:))
      ENDDO


      END SUBROUTINE


      LOGICAL FUNCTION LINUX_CHECK()
C     ============================
C Function checks whether the operating system is Linux (.TRUE.) or
C Windows (.FALSE.)
!     USE DFPORT
      CHARACTER(100) CWD
      CWD=''
!!      I=GETDRIVEDIRQQ(CWD)
      CALL GETCWD(CWD)
      L=LEN_TRIM(CWD)
      I=INDEX(CWD,'/')
      IF(I .EQ. 0) THEN
         LINUX_CHECK=.FALSE.
      ELSE
         LINUX_CHECK=.TRUE.
      ENDIF
      END

      integer function MAC_CHECK()
        character(100) uname_output
        call system("uname > uname.out")
        open(800,file='uname.out',status='UNKNOWN')
        read(800,*)uname_output
        write(*,*)uname_output
        if(trim(uname_output).eq.'Darwin') then
          MAC_CHECK=1
        else
          MAC_CHECK=0
        end if
      end


      CHARACTER(*) FUNCTION UCASE(TEXT)
C     ======================
C
C***********************************************************************
C
C Subroutine changes all letters in TEXT to upper case
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      CHARACTER(*) TEXT
      UCASE=' '
      DO I=1,LEN_TRIM(TEXT)
       KEY=ICHAR(TEXT(I:I))
       IF(KEY.GT.96 .AND. KEY.LT. 123) THEN
          UCASE(I:I)=CHAR(KEY-32)
       ELSE
          UCASE(I:I)=CHAR(KEY)
       ENDIF
      ENDDO
      END
