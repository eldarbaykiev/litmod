!Subroutine MATER


      SUBROUTINE MATER
      USE M_temperatures
      USE M_grid_parameters
      USE M_layers
      USE M_material
!LK
      USE M_mant
      USE M_phases
      USE M_sublit
      USE M_conductivity


      REAL(8) z_top,z_boto
      REAL(8) C_w_up,C_w_dw,slp_cw,z_moho,melt_up,melt_dw,slp_melt
      REAL(8) Cw_melt_up,Cw_melt_dw,slp_Cw_melt,Na_up,Na_dw,slp_Na
      REAL(8) Si_up,Si_dw,slp_Si
      LOGICAL C_w_onoff_l,melt_onoff_l
      CHARACTER(2) chdm,idmm_c
      INTEGER N_mant
!LK
      CHARACTER(len=255) :: path_mnt
      REAL(8) P_dm_lst,P_dm,T_dm_lst,T_dm


      REAL(8) E_ice

      !ICE PARAM FROM WINTERC



      ALLOCATE (layers(N_x,N_y,0:N_lay+1),c_zs(N_x,N_y),
     *          d_z_topo(N_x,N_y))
      N_mant_max=0
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
      open(20,file=fil,status='UNKNOWN')
       DO I_y=1,N_y
        DO I_x=1,N_x
         read(20,*) dumm,dumm,layers(I_x,I_y,i)
         layers(I_x,I_y,i)=-layers(I_x,I_y,i)*1D3
        ENDDO
       ENDDO
      close(20)
      ENDDO
!Bottom layer
      layers(:,:,N_lay+1)=E_max-(N_z-1)*d_z

      IF (topo==1) THEN
       OPEN(1,file='./layers_xy/elev_c.xyz',status='UNKNOWN')
      ELSE
       OPEN(1,file='./layers_xy/elev_fil_c.xyz',status='UNKNOWN')
      ENDIF
      OPEN(7,file='OUT_xyz/NODOS',status='UNKNOWN')
      OPEN(777,file='E_real-discr.xyz',status='UNKNOWN')
      ALLOCATE (mat(N_x,N_y,N_z),DUM(N_x,N_y,N_z),
     *          DUM_sublit(N_x,N_y,N_z))
      ALLOCATE (T(N_x,N_y,N_z),E(N_x,N_y))
      ALLOCATE (T_old(N_x,N_y,N_z),k_T_P(N_x,N_y,N_z))
      mat=0
      DUM=0
      k_T_P=0D0
      T_sub=0D0
      N_mant=0
      DO I_y=1,N_y
       y=d_y*(I_y-1)
       DO I_x=1,N_x
        x=d_x*(I_x-1)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
       read(1,*) dumm,dumm,E(I_x,I_y)

       !assign ice



       IF (E(I_x,I_y)>-1D-2.AND.E(I_x,I_y)<0D0) E(I_x,I_y)=0D0
       layers(I_x,I_y,0)=E(I_x,I_y)
!The Moho deph, z_c, and the Litho-asthenosphere boundary deph, z_l,
!must be ALWAYS NEGATIVE
!The deph of the base of sediments,z_s, can be either >0 or <0
        iiq=0
        DO i_lay=0,N_lay-1
         z_top=layers(I_x,I_y,i_lay)
         z_boto=layers(I_x,I_y,i_lay+1)
        IF(z_boto>z_top) THEN
         layers(I_x,I_y,i_lay+1)=layers(I_x,I_y,i_lay)
         IF(z_boto-z_top>d_z) THEN
          WRITE(*,*)'WARNING, HIERARCHY VIOLATION (>d_z) BETWEEN'
          WRITE(*,*)'LAYER ',i_lay,'AND LAYER ',i_lay+1
          WRITE(*,*) 'AT X= ',x,'AND Y= ',y
          WRITE(*,*)'AUTOMATICALLY CORRECTED'
         ENDIF
         z_boto=z_top
        ENDIF
        difff=ABS(z_boto-z_top)
        IF(prop(i_lay+1)%litho>0.AND.difff<3*d_z) THEN
          !for mantle-mantle only
          if(prop(i_lay+1)%crst==0)then
            layers(I_x,I_y,i_lay+1)=layers(I_x,I_y,i_lay)
          endif
        ENDIF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        DO I_z=1,N_z
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        z=E_max-(I_z-1)*d_z
        !E_ice = E(I_x,I_y)

         IF (i_lay==0) THEN
          IF (z>0.AND.z>z_top) THEN
            mat(I_x,I_y,I_z)=-100
            T(I_x,I_y,I_z)=Ts
          ELSEIF (z<=0.AND.z>z_top) THEN
            mat(I_x,I_y,I_z)=-200
            T(I_x,I_y,I_z)=Ts
            !REMOVEDBUG
          !ELSEIF (z<=0.AND.-z>E(I_x,I_y)) THEN
          !  mat(I_x,I_y,I_z)=-300
          !  T(I_x,I_y,I_z)=Ts
          !  IF(E_ice==E(I_x,I_y)) mat(I_x,I_y,I_z)=-100
!          write(*,*)'ice',z
          ELSEIF (z<layers(I_x,I_y,N_lay)) THEN
            mat(I_x,I_y,I_z)=0
            T(I_x,I_y,I_z)=Ta
          ENDIF
         ENDIF
         IF (z<=z_top.AND.z>=z_boto) THEN
          mat(I_x,I_y,I_z)=i_lay+1
          T(I_x,I_y,I_z)=-1000D0
          IF (i_lay==0.AND.iiq==0) THEN
           write(777,*)x,y,z,z_top
           iiq=1
          ENDIF
         ENDIF
!Maximum number of nodes in the mantle(s)
         IF (mat(I_x,I_y,I_z)<0) THEN
          DUM(I_x,I_y,I_z)=-1
         ELSE
          DUM(I_x,I_y,I_z)=mat(I_x,I_y,I_z)
         ENDIF
!Check for the presence of mantle-crustal-mantle stratification
!(not allowed)
          IF(I_z>1) THEN
           IF(prop(DUM(I_x,I_y,I_z))%litho==0.AND.
     *        prop(DUM(I_x,I_y,I_z-1))%litho/=0) THEN
            WRITE(*,*)'MANTLE-CRUST-MANTLE STRATIFICATION NOT ALLOWED!!'
            WRITE(*,*)'CHECK THE ORDER OF THE LAYERS'
            WRITE(*,*)'LitMod3D ABORTS...'
            STOP
           ENDIF
          ENDIF

!

          IF(prop(DUM(I_x,I_y,I_z))%litho/=0) THEN
           IF (prop(DUM(I_x,I_y,I_z-1))%crst==1) THEN
            N_mant=1
           ELSE
            N_mant=N_mant+1
           ENDIF
           IF(N_mant_max<N_mant) N_mant_max=N_mant
          ENDIF
!Initial value for the T-P-dependent thermal conductivity
          IF (mat(I_x,I_y,I_z)>0) THEN
           k_T_P(I_x,I_y,I_z)=prop(mat(I_x,I_y,I_z))%K
          ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       ENDDO
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
       ENDDO
      ENDDO
      write(*,*)'MAX NUMBER OF MANTLE NODES-->  ',
     *           N_mant_max

!Boundary conditions for the thermal calculation
      DUM=mat
      DUM_sublit=mat
      DO I_x=1,N_x
       DO I_y=1,N_y
        DO I_z=2,N_z-1
         IF(mat(I_x,I_y,I_z)<=0) CYCLE
         ki=0
         ik=0
         DO ii=-1,1
         DO jj=-1,1
         DO kk=-1,1
         IF(I_x==1.AND.ii==-1.OR.I_x==N_x.AND.ii==1)CYCLE
         IF(I_y==1.AND.jj==-1.OR.I_y==N_y.AND.jj==1)CYCLE
         IF(((ii==0).AND.(jj==0).AND.(kk==0))) CYCLE
         IF(ABS(ii)+ABS(jj)+ABS(kk)>1) CYCLE
         IF(ki>0 .OR. ik>0) EXIT
         IF (mat(I_x+ii,I_y+jj,I_z+kk)<0) THEN
          DUM(I_x,I_y,I_z)=-DUM(I_x,I_y,I_z)-100
          T(I_x,I_y,I_z)=Ts
          ki=1
          EXIT
           ELSEIF (mat(I_x+ii,I_y+jj,I_z+kk)==0) THEN
          DUM(I_x,I_y,I_z)=-DUM(I_x,I_y,I_z)
          T(I_x,I_y,I_z)=Ta
          ik=1
          EXIT
          ENDIF
          ENDDO
          IF(ki>0 .OR. ik>0) EXIT
          ENDDO
          IF(ki>0 .OR. ik>0) EXIT
          ENDDO
        ENDDO
       ENDDO
      ENDDO

!Load sublithospheric mantle bodies
      IF (N_sublit_t>0) THEN
!Change material number
      DO i=1,N_sublit_t
       z=E_max-(IND_SUB(i,3)-1)*d_z
       IF(IND_SUB(i,1)==-1) CYCLE
       IF(z>layers(IND_SUB(i,1),IND_SUB(i,2),N_lay)) CYCLE
       mat(IND_SUB(i,1),IND_SUB(i,2),IND_SUB(i,3))=I_sublit(i)+N_lay
      ENDDO
!Change litho number
      DO i=N_lay+1,N_lay+N_sublit
       prop(i)%litho=litho_sub(i-N_lay)
      ENDDO
!Check if the sublithospheric materials are already present in
!layers.info
      RON=.TRUE.
      I_subRON=N_sublit
      DO i=1,N_sublit
       DO j=0,N_lay
        IF(litho_sub(i)==prop(j)%litho)THEN
         IF(RON(i)) I_subRON=I_subRON-1
         RON(i)=.FALSE.
        ENDIF
       ENDDO
      ENDDO

      ENDIF
      write(*,*)'LOADING MANTLE FILES...'
!LK Load look-up tables information
       open(15,file='mnt.info',status='OLD',IOSTAT=II)
       IF(II==0) THEN
        jfl=0
        N_pm=1
        DO
         READ(15,*,IOSTAT=jfl) I_mant_file(N_pm),mant_file(N_pm)
         IF (jfl/=0) EXIT
         mant_file_sys(N_pm)=TRIM(mant_file(N_pm))//'_sys.dat'
         mant_file_ph(N_pm)=TRIM(mant_file(N_pm))//'_phases.dat'
         N_pm=N_pm+1
        ENDDO
       ELSE
        WRITE(*,*)'FILE mnt.info DOES NOT EXIST'
       STOP
       ENDIF
       close (15)

!Load look-up tables
       OPEN(111,file='Perplex_failure.dat',status='UNKNOWN')
       CALL GETCWD(path_mnt)
       CALL CHDIR (TRIM(path_mnt)//'/mant_data')
       IF (N_sublit_t>0) I_mnt=I_mnt+I_subRON
       ALLOCATE(LK_index(I_mnt),d_P_LK(I_mnt),d_T_LK(I_mnt),
     *          T_LK_0(I_mnt),P_LK_0(I_mnt),N_LK(I_mnt),
     *          N_LK_T(I_mnt),N_LK_P(I_mnt),litho2LK(I_mnt))
       I_mnt=0
       DO i_lay=0,N_lay+N_sublit
        IF (prop(i_lay)%litho==0) CYCLE
        IF (i_lay>N_lay) THEN
         IF(.NOT.RON(i_lay-N_lay)) CYCLE
        ENDIF
        I_mnt=I_mnt+1
        DO i=1,N_pm
         IF(I_mant_file(i)==prop(i_lay)%litho) THEN
          litho2LK(I_mnt)=prop(i_lay)%litho
          OPEN(100,file=mant_file_sys(i),status='OLD')
          READ(100,*)N_LK(I_mnt)
          READ(100,*)T_LK_0(I_mnt),P_LK_0(I_mnt)
          T_dm_lst=T_LK_0(I_mnt)
          P_dm_lst=P_LK_0(I_mnt)
!Find the number of nodes in the T-P space and the grid steps
          DO j=2,N_LK(I_mnt)
           READ(100,*)T_dm,P_dm
           IF(j==2) d_T_LK(I_mnt)=ABS(T_dm-T_dm_lst)
           IF(P_dm_lst/=P_dm) THEN
            N_LK_T(I_mnt)=j-1
            N_LK_P(I_mnt)=N_LK(I_mnt)/N_LK_T(I_mnt)
            d_P_LK(I_mnt)=ABS(P_dm-P_dm_lst)
            EXIT
           ENDIF
           P_dm_lst=P_dm
          ENDDO
          CLOSE(100)
          LK_index(I_mnt)=i
         ENDIF
        ENDDO
       ENDDO


      !!! Convert units.
      !!! The Temperature in the Perple_X file is in Kelvin.
      !!! The Pressure in the Perple_X file is in Bar.
      !!! write(*,*)'T_LK_0',T_LK_0
      !!! stop !!! MATTIA
      T_LK_0=T_LK_0-273.15D0
      P_LK_0=P_LK_0*1D5
      d_P_LK=d_P_LK*1D5
!Number of crustal bodies
       I_crst=N_lay-I_mnt+1
       N_LK_T_max=MAXVAL(N_LK_T(:))
       N_LK_P_max=MAXVAL(N_LK_P(:))
       ALLOCATE(LK_vec(I_mnt,N_LK_T_max,N_LK_P_max,5),dm_vec(5))
       WRITE(*,*)'READING MANTLE FILE: '
       DO k=1,I_mnt
        OPEN(100,file=mant_file_sys(LK_index(k)),status='OLD')
        WRITE(*,*)TRIM(mant_file(LK_index(k)))
        READ(100,*)
        LK_vec(k,:,:,:)=0
        DO
         READ(100,*,IOSTAT=jfl)T_dm,P_dm,dm_vec(1:5)
         IF (jfl/=0) EXIT
         IF((T_dm>1823D0.AND.T_LK_0(k)<300D0).OR.P_dm>145100D0
     *      .OR.T_dm>2200D0) CYCLE
            i_tb=NINT((T_dm-T_LK_0(k)-273.15D0)/d_T_LK(k)+1)
            j_tb=NINT((P_dm*1D5-P_LK_0(k))/d_P_LK(k)+1)
            LK_vec(k,i_tb,j_tb,1:5)=dm_vec(1:5)
        ENDDO
!In points where the Gibbs energy minimization failed, copy last point
        DO j=1,N_LK_P(k)
         DO i=1,N_LK_T(k)
          IF(LK_vec(k,i,j,1)==0) THEN
          write(111,*)'PERPLEX MINIMIZATION FAILURE IN THE SYSTEM FILE'
          WRITE(111,*)TRIM(mant_file(LK_index(k)))
          write(111,*)'T(K),P(kbar)= ',(i-1)*d_T_LK(k)+273.5D0,
     *                       ((j-1)*d_P_LK(k)+P_LK_0(k))*1D-5
           IF(i==1) THEN
            IF(j==1) THEN
             WRITE(*,*)'PERPLEX MINIMIZATION FAILURE IN THE FIRST NODE
     * OF THE TABLE',TRIM(mant_file_sys(LK_index(k)))
             WRITE(*,*)'TRY ANOTHER BULK COMPOSITION'
             write(*,*)'LitMod3D aborts...'
             STOP
            ENDIF
            LK_vec(k,i,j,1:5)=LK_vec(k,i,j-1,1:5)
            write(111,*)'COPYING THE PROPERTIES OF THE CLOSEST
     *SUCCESFUL NODE, T(K),P(kbar)=',(i-1)*d_T_LK(k)+273.5D0,
     *                       ((j-2)*d_P_LK(k)+P_LK_0(k))*1D-5
           ELSE
            LK_vec(k,i,j,1:5)=LK_vec(k,i-1,j,1:5)
            write(111,*)'COPYING THE PROPERTIES OF THE CLOSEST
     *SUCCESFUL NODE, T(K),P(kbar)=',(i-2)*d_T_LK(k)+273.5D0,
     *                       ((j-1)*d_P_LK(k)+P_LK_0(k))*1D-5
           ENDIF
          ENDIF
         ENDDO
        ENDDO

       CLOSE(100)
       ENDDO
!Read phases info
      ALLOCATE(phases(I_mnt,N_LK_T_max,N_LK_P_max,N_tot_ph),
     *         vol_frac(I_mnt,N_LK_T_max,N_LK_P_max,N_tot_ph),
     *         wt_MgO(I_mnt,N_LK_T_max,N_LK_P_max,N_tot_ph),
     *         wt_FeO(I_mnt,N_LK_T_max,N_LK_P_max,N_tot_ph),
     *         I_phases(I_mnt,N_LK_T_max,N_LK_P_max,7))
      I_phases=0
      phases=''
      vol_frac=0D0
      wt_MgO=0D0
      wt_FeO=0D0
      DO k=1,I_mnt
        OPEN(100,file=mant_file_ph(LK_index(k)),status='OLD')
        jfl=0
        N_phases=0
        N_pt_LK=0
        DO
          READ(100,'(I2)',ADVANCE='NO',IOSTAT=jfl) isys
          IF (jfl/=0) EXIT
          IF(isys==0) THEN
            N_phases=N_phases+1
            IF(T_dm>1823D0.AND.T_LK_0(k)<300D0) THEN
              !!! write(*,*)'T_dm',T_dm
              !!! write(*,*)'T_LK_0(k)',T_LK_0(k)
              write(*,*)'WARNING: T in the phases file is out of range.'
              T_dm=1823D0
            endif
         IF(P_dm>145100D0) THEN
	  !RETURN IT BACK ELDAR BAYKIEV
          !WRITE(*,*)'p IN THE PHASES FILE OUT OF RANGE!!'
          P_dm=145000D0
         ENDIF
         i=NINT((T_dm-T_LK_0(k)-273.15D0)/d_T_LK(k)+1)
         j=NINT((P_dm*1D5-P_LK_0(k))/d_P_LK(k)+1)
         READ(100,*)phases(k,i,j,N_phases),dm,vol_frac(k,i,j,N_phases),
     *   dm,dm,dm,
     *   wt_FeO(k,i,j,N_phases), wt_MgO(k,i,j,N_phases)
         DO iph=1,N_tot_ph
          IF (phases(k,i,j,iph)(1:2)=='O(')
     *       I_phases(k,i,j,1)=iph
          IF (phases(k,i,j,iph)(1:3)=='Opx'.OR.
     *        phases(k,i,j,iph)(1:5)=='CrOpx')
     *       I_phases(k,i,j,2)=iph
          IF (phases(k,i,j,iph)(1:3)=='Cpx')
     *       I_phases(k,i,j,3)=iph
          IF (phases(k,i,j,iph)(1:2)=='Gt'.OR.
     *        phases(k,i,j,iph)(1:4)=='CrGt')
     *       I_phases(k,i,j,4)=iph
          IF (phases(k,i,j,iph)(1:2)=='Sp'.OR.
     *        phases(k,i,j,iph)(1:4)=='CrSp')
     *       I_phases(k,i,j,5)=iph
          IF (phases(k,i,j,iph)(1:2) =='an'
     *        .OR.phases(k,i,j,iph)(1:2)=='Pl')
     *       I_phases(k,i,j,6)=iph
          IF (phases(k,i,j,iph)(1:2)=='C2'.OR.
     *        phases(k,i,j,iph)(1:2)=='c2')
     *       I_phases(k,i,j,7)=iph
         ENDDO
        ELSEIF (isys==-2) THEN
         N_pt_LK=N_pt_LK+1
         N_phases=0
         READ(100,*)
        ELSE
         N_pt_LK=N_pt_LK+1
         N_phases=0
         READ(100,*)T_dm,P_dm
        ENDIF
       ENDDO
       CLOSE(100)
       DO j=1,N_LK_P(k)
        DO i=1,N_LK_T(k)
        IF(phases(k,i,j,1)=='') THEN
        write(111,*)'GIBBS ENERGY MINIMIZATION FAILURE IN PHASES FILE: '
        WRITE(111,*)TRIM(mant_file(LK_index(k)))
         IF(i==1) THEN
          IF(j==1) THEN
           WRITE(*,*)'PERPLEX MINIMIZATION FAILURE IN THE FIRST NODE
     * OF THE TABLE',TRIM(mant_file_sys(LK_index(k)))
           WRITE(*,*)'TRY ANOTHER BULK COMPOSITION'
           write(*,*)'LitMod3D aborts...'
            STOP
           ENDIF
           phases(k,i,j,1:N_phases)=phases(k,i,j-1,1:N_phases)
           vol_frac(k,i,j,1:N_phases)=vol_frac(k,i,j-1,1:N_phases)
           wt_FeO(k,i,j,1:N_phases)=wt_FeO(k,i,j-1,1:N_phases)
           wt_MgO(k,i,j,1:N_phases)=wt_MgO(k,i,j-1,1:N_phases)
           I_phases(k,i,j,1:7)=I_phases(k,i,j-1,1:7)
           write(111,*)'COPYING THE PROPERTIES OF THE CLOSEST
     *SUCCESFUL NODE, T(K),P(kbar)=',(i-1)*d_T_LK(k)+273.5D0,
     *                      ((j-2)*d_P_LK(k)+P_LK_0(k))*1D-5
          ELSE
           phases(k,i,j,1:N_phases)=phases(k,i-1,j,1:N_phases)
           vol_frac(k,i,j,1:N_phases)=vol_frac(k,i-1,j,1:N_phases)
           wt_FeO(k,i,j,1:N_phases)=wt_FeO(k,i-1,j,1:N_phases)
           wt_MgO(k,i,j,1:N_phases)=wt_MgO(k,i-1,j,1:N_phases)
           I_phases(k,i,j,1:7)=I_phases(k,i-1,j,1:7)
           write(111,*)'COPYING THE PROPERTIES OF THE CLOSEST
     *SUCCESFUL NODE, T(K),P(kbar)=',(i-2)*d_T_LK(k)+273.5D0,
     *                       ((j-1)*d_P_LK(k)+P_LK_0(k))*1D-5
          ENDIF
        ENDIF
        ENDDO
       ENDDO

      ENDDO
      CLOSE(111)
      CALL CHDIR(path_mnt)
!!Read mantle water contents
      ALLOCATE(C_w_ol(N_z),C_w_opx(N_z),C_w_cpx(N_z),C_w_gt(N_z),
     *         C_w_onoff(0:N_lay),melt_onoff(0:N_lay),
     *         melt(N_z),Na_melt(N_z),Si_melt(N_z),Cw_melt_melt(N_z))
      C_w_onoff_l=.FALSE.
      C_w_onoff(0)=0
      melt_onoff(0)=0
      C_w_ol=0D0
      C_w_opx=0D0
      C_w_cpx=0D0
      C_w_gt=0D0
      melt=0D0
      Na_melt=0D0
      Si_melt=0D0
      Cw_melt_melt=0D0
      INQUIRE(FILE='C_w_mant_3D.info',EXIST=C_w_onoff_l)
      INQUIRE(FILE='melt_mant_3D.info',EXIST=melt_onoff_l)

      IF (C_w_onoff_l) THEN
      OPEN(700,file='C_w_mant_3D.info',status='UNKNOWN')
      DO kk=1,3
       READ(700,*)
      ENDDO
       jfl=0
       N_cw=0
       DO
        READ(700,*,IOSTAT=jfl) dm
        IF (jfl/=0) EXIT
        N_cw=N_cw+1
       ENDDO
       REWIND(700)
       N_cw=N_cw-1
       ALLOCATE(C_w_mod(N_cw,4,2),z_cw(N_cw))
       READ(700,*)
!Dry/wet layers
       READ(700,*)C_w_onoff(1:N_lay)
       WRITE(*,*)'ELECTRICAL CONDUCTIVITY: WET/DRY MANTLE LAYERS'
       DO l=1,N_lay
        IF(prop(l)%litho/=0) THEN
         IF(C_w_onoff(l)==0) THEN
          WRITE(*,*) 'LAYER:',l,TRIM(prop(l)%name_mat),'--> DRY'
         ELSEIF(C_w_onoff(l)==1) THEN
          WRITE(*,*) 'LAYER:',l,TRIM(prop(l)%name_mat),'--> WET'
         ENDIF
        ENDIF
       ENDDO
       READ(700,*)
       READ(700,*) z_moho
       z_moho=z_moho*1D3
       DO i=1,N_cw
        READ(700,*)z_cw(i),C_w_mod(i,1,:),C_w_mod(i,2,:),C_w_mod(i,3,:),
     *             C_w_mod(i,4,:)
        z_cw(i)=z_cw(i)*1D3
       ENDDO

       DO i=1,N_cw
        DO j=1,4
         z_bot=z_cw(i)
         C_w_up=C_w_mod(i,j,1)
         C_w_dw=C_w_mod(i,j,2)
         IF(i==1) THEN
          z_top=z_moho
         ELSE
          z_top=z_cw(i-1)
         ENDIF
         I_z1=NINT((z_top+E_max)/d_z)+1
         I_z2=NINT((z_bot+E_max)/d_z)+1
         IF(I_z2>N_z)I_z2=N_z
         slp_cw=(C_w_dw-C_w_up)/(I_z2-I_z1)
         DO I_zz=I_z1+1,I_z2
          IF(j==1) THEN
           C_w_ol(I_zz)=C_w_up+slp_cw*(I_zz-I_z1)
          ELSEIF(j==2) THEN
           C_w_opx(I_zz)=C_w_up+slp_cw*(I_zz-I_z1)
          ELSEIF(j==3) THEN
           C_w_cpx(I_zz)=C_w_up+slp_cw*(I_zz-I_z1)
          ELSEIF(j==4) THEN
           C_w_gt(I_zz)=C_w_up+slp_cw*(I_zz-I_z1)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
       CLOSE(700)
       DEALLOCATE(z_cw)
      ELSE
       C_w_onoff=0
       C_w_ol=0D0
       C_w_opx=0D0
       C_w_cpx=0D0
       C_w_gt=0D0
       WRITE(*,*)'ELECTRICAL CONDUCTIVITY: FILE C_w_mant_3D.info DOES
     *NOT EXIST--> ALL MANTLE LAYERS ARE ASSUMED TO BE DRY '
      ENDIF


!Read melt distribution
      IF (melt_onoff_l) THEN
      OPEN(700,file='melt_mant_3D.info',status='UNKNOWN')
      DO kk=1,3
       READ(700,*)
      ENDDO
       jfl=0
       N_cw=0
       DO
        READ(700,*,IOSTAT=jfl) dm
        IF (jfl/=0) EXIT
        N_cw=N_cw+1
       ENDDO
       REWIND(700)
       N_cw=N_cw-1
       ALLOCATE(melt_mod(N_cw,8),z_cw(N_cw))
       READ(700,*)
!Melt in layers
       READ(700,*)melt_onoff(1:N_lay)
       WRITE(*,*)'ELECTRICAL CONDUCTIVITY: MELT IN LAYERS'
       DO l=1,N_lay
        IF(melt_onoff(l)==0) THEN
        ELSEIF(melt_onoff(l)==1) THEN
         WRITE(*,*) 'LAYER:',l,TRIM(prop(l)%name_mat),'--> MELT'
        ENDIF
       ENDDO
       READ(700,*)
       READ(700,*) z_moho
       z_moho=z_moho*1D3
       DO i=1,N_cw
        READ(700,*)z_cw(i),melt_mod(i,1:8)
        z_cw(i)=z_cw(i)*1D3
       ENDDO
       DO i=1,N_cw
         z_bot=z_cw(i)
         melt_up=melt_mod(i,1)
         melt_dw=melt_mod(i,2)
         Na_up=melt_mod(i,3)
         Na_dw=melt_mod(i,4)
         Si_up=melt_mod(i,5)
         Si_dw=melt_mod(i,6)
         Cw_melt_up=melt_mod(i,7)
         Cw_melt_dw=melt_mod(i,8)
         IF(i==1) THEN
          z_top=z_moho
         ELSE
          z_top=z_cw(i-1)
         ENDIF
         I_z1=NINT((z_top+E_max)/d_z)+1
         I_z2=NINT((z_bot+E_max)/d_z)+1
         IF(I_z2>N_z)I_z2=N_z
         slp_melt=(melt_dw-melt_up)/(I_z2-I_z1)
         slp_Na=(Na_dw-Na_up)/(I_z2-I_z1)
         slp_Si=(Si_dw-Si_up)/(I_z2-I_z1)
         slp_Cw_melt=(Cw_melt_dw-Cw_melt_up)/(I_z2-I_z1)
         DO I_zz=I_z1+1,I_z2
          melt(I_zz)=melt_up+slp_melt*(I_zz-I_z1)
          Na_melt(I_zz)=Na_up+slp_Na*(I_zz-I_z1)
          Si_melt(I_zz)=Si_up+slp_Si*(I_zz-I_z1)
          Cw_melt_melt(I_zz)=Cw_melt_up+slp_Cw_melt*(I_zz-I_z1)
         ENDDO
       ENDDO
       CLOSE(700)
       DEALLOCATE(z_cw)
      ELSE
       melt_onoff=0
       melt=0D0
       WRITE(*,*)'ELECTRICAL CONDUCTIVITY: FILE melt_mant_3D.info DOES
     *NOT EXIST--> NO MELT WILL BE CONSIDERED IN THE MODEL '
      ENDIF


!Write NODOS file (info about node material)

      DO I_x=1,N_x
       x=d_x*(I_x-1)
       DO I_y=1,N_y
        y=d_y*(I_y-1)
         icon=1
        DO I_z=1,N_z
         z=E_max-(I_z-1)*d_z
       write(7,'(3F13.3,I5)')x,y,z,mat(I_x,I_y,I_z)
        ENDDO
       ENDDO
      ENDDO
      close (7)
      close(777)

      END SUBROUTINE
