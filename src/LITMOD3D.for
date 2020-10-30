



!!!=============================================================================
!!!*****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   LitMod   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     Version: 4.0    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!=============================================================================


!!! To compile (on macOS Catalina 10.15.5, Xcode 11.5):
!!! gfortran -o litmod  modules.for LITMOD3D.for SUB* -Ofast -fno-automatic -fd-lines-as-code -ffixed-line-length-none -std=legacy


      program litmod


      use M_grid_parameters
      use M_material
      use M_temperatures
      use M_layers
      use M_L_EQ_SYS
      use M_columns
      use M_lsq_plane
      use M_mant
      use M_sublit
      use M_conductivity
      use M_fugacity_ol_solu, ONLY:c_T
      use m_itherc

!MINEOS
      use M_surf_waves
      use M_column_1D
      use M_d

      REAL(8),ALLOCATABLE::dmmy_vec(:)

      !!! Coefficients for pure water EoS from Pitzer and Sterner, 1994.
      c_T(:,:)=0
      c_T(1,3)=.24657688D6
      c_T(1,4)=.51359951D2
      c_T(2,3)=.58638965
      c_T(2,4)=-.28646939D-2
      c_T(2,5)=.31375577D-4
      c_T(3,3)=-.6278384D1
      c_T(3,4)=0.14791599D-1
      c_T(3,5)=.35779579D-3
      c_T(3,6)=.15432925D-7
      c_T(4,4)=-.42719875
      c_T(4,5)=-.16325155D-4
      c_T(5,3)=.56654978D4
      c_T(5,4)=-.16580167D2
      c_T(5,5)=.76560762D-1
      c_T(6,4)=.10917883
      c_T(7,1)=.38878656D13
      c_T(7,2)=-.13494878D9
      c_T(7,3)=.30916564D6
      c_T(7,4)=.75591105D1
      c_T(8,3)=-.65537898D5
      c_T(8,4)=.18810675D3
      c_T(9,1)=-.14182435D14
      c_T(9,2)=.18165390D9
      c_T(9,3)=-.19769068D6
      c_T(9,4)=-.23530318D2
      c_T(10,3)=.92093375D5
      c_T(10,4)=.12246777D3


!Definition of materials:
!////////////////////////////////////////////////
!Material 0: Sublithospheric mantle
      prop(0)%name_mat='Sub_litho_mantle'
      prop(0)%rho=3390D0
      prop(0)%beta=7.D-12
      prop(0)%alpha=3D-5
      prop(0)%litho=99
!Material -1: sea water
      prop(-1) % name_mat='Water'
      prop(-1) % rho=1030D0
      prop(-1)%litho=0
!////////////////////////////////////////////////
!Gravity atracction,g, gradient
      d_gz =(g_400-g_s )/400D3
!Bodies info
      open(15,file='layers.info',status='UNKNOWN')
      jfl=0
      N_lay=0
!LK
      I_mnt=1
      read(15,*)
      DO
       read(15,*,IOSTAT=jfl)j,prop(j)%name_mat,prop(j)%rho,
     *                      prop(j)%alpha,prop(j)%K,
     *                      prop(j)%A,prop(j)%beta,prop(j)%gmma_T,
     *                      prop(j)%K_0p,prop(j)%K_T,prop(j)%litho,
     *                      prop(j)%Vp,prop(j)%Vs,prop(j)%crst
      ! *                      prop(j)%cla , prop(j)%inpfil
       IF (jfl/=0) EXIT
       N_lay=N_lay+1
       IF (prop(j)%litho/=0) THEN
        I_mnt=I_mnt+1
!
       ENDIF

       IF (prop(j)%crst==0) N_lay_cr=j-1





      ENDDO
      close (15)
      IF (prop(N_lay)%rho<0.OR.MAXVAL(prop(:)%litho)==0) THEN
       WRITE(*,*)'layers.info INPUT FILE SEEMS TO BE SUITABLE FOR
     * GEO3Dmod. CHECK THE INPUT PARAMETER classic IN LITMOD3D.job  '
       WRITE(*,*)'LitMod3D aborts...'
       STOP
      ENDIF
!Read Sublithospheric mantle bodies (sublith_an.xyz)
      INQUIRE(FILE='sublith_an.xyz',EXIST=sub_onoff)
      i_status=0
      N_sublit=0
      N_sublit_t=0
      IF(sub_onoff) THEN

       OPEN(400,file='sublith_an.xyz',status='UNKNOWN')
       DO
        READ(400,'(A1)',IOSTAT=i_status) chdm_s
        IF (i_status/=0) EXIT
        N_sublit_t=N_sublit_t+1
        IF (ICHAR(chdm_s)==42) THEN
         N_sublit=N_sublit+1
        ENDIF
       ENDDO
       REWIND(400)
!       write(*,*)'N_sublit,N_sublit_t',N_sublit,N_sublit_t
       IF(N_sublit_t>0) THEN
        ALLOCATE(D_T_sub(N_sublit),litho_sub(N_sublit),
     *           I_sublit(N_sublit_t),IND_SUB(N_sublit_t,3),
     *           RON(N_sublit))
        N_sublit=0
        DO i=1,N_sublit_t
         READ(400,'(A1)',ADVANCE='NO',IOSTAT=i_status) chdm_s
         IF (ICHAR(chdm_s)==42) THEN
          N_sublit=N_sublit+1
          READ(400,*)D_T_sub(N_sublit),litho_sub(N_sublit)
          IND_SUB(i,:)=-1
!         write(*,'(A,1X,F7.2,1X,A,1X,I3)')'D_T(�C)',
!     *         D_T_sub(N_sublit),'litho n�: ',litho_sub(N_sublit)
         ELSE
          READ(400,*)IND_SUB(i,3),IND_SUB(i,1),IND_SUB(i,2)
         ENDIF
         I_sublit(i)=N_sublit
        ENDDO
       ENDIF
      ENDIF
!Solving the equation system
      ITMAX=150000
      ERROR=1.D-7

!Input parameters
      OPEN(15,FILE='LITMOD3D.info',STATUS='OLD',IOSTAT=II)
      IF(II==0) THEN
      READ(15,*) L_x,L_y,E_max,N_x,N_y,N_z,d_z,topo,temp_calc,Ts,Ta,
     *           rho_red,dm,E_cali,ol_mod,gt_mod,opx_mod,mlt_mod,
     *           melt_mix,Z_Usec,volexp,d_size


       volexp=volexp/1000000.0
!Elevation calibration
       IF (E_cali==1) THEN
!Stixrude et al., 05 (5 oxides) (d_z=2 km)
        rho_m_MOR=3421D0
       ELSEIF (E_cali==2) THEN
!H&P98 modified by JC Afonso (2010) (d_z=2 km)
        rho_m_MOR=3408D0
       ELSEIF (E_cali==3) THEN
!Stixrude et al., 08 (6 oxides) (d_z=2 km)
        rho_m_MOR=3424.5D0
       ENDIF
       CLOSE(15)


!	stop
      ELSE
      WRITE(*,*)'FILE LITMOD3D.info DOES NOT EXIST'
      STOP
      ENDIF

!Unit conversion km-->m
      E_max=E_max*1D3
      L_x=L_x*1D3
      L_y=L_y*1D3
      d_x=L_x/(N_x-1)
      d_y=L_y/(N_y-1)
      d_z=d_z*1D3
      Z_Usec=Z_Usec*1D3
!Check if the model enters in the transition zone
      IF (-E_max+(N_z-1)*d_z-400D3>d_z) THEN
       WRITE(*,*)'THE VERTICAL EXTENSION OF THE MODEL EXCEDES THE
     *DEPTH OF THE TRANSITION ZONE: REDEFINING N_z'
       N_z=INT((400D3+E_max)/d_z)+1
       WRITE(*,*)'NEW N_z= ',N_z
       OPEN(15,FILE='LITMOD3D.info',STATUS='OLD')
       WRITE(15,'(2F12.2,F4.1,3I4,F6.2,2I2,3F8.1,5I2)')
     *L_x*1D-3,L_y*1D-3,
     *E_max*1D-3,N_x,N_y,N_z,d_z*1D-3,topo,temp_calc,Ts,Ta,
     *rho_red,1,E_cali,ol_mod,gt_mod,opx_mod

      ENDIF

!Total number of effective nodes
      NOD_TOT=(N_x-2)*(N_y-2)*(N_z-2)

      write(*,*)''
      write(*,*)'##################### LitMod3D ####################'
      write(*,*)'#                    Version 4.0                  #'
      write(*,*)'#                     (itherc)                    #'
      write(*,*)'#                    written by                   #'
      write(*,*)'#              J. Fullea and J.C Afonso           #'
      write(*,*)'#                                                 #'
      write(*,*)'#               in collaboration with             #'
      write(*,*)'#                   Hermann Zeyen                 #'
      write(*,*)'#                                                 #'
      write(*,*)'# Reference:                                      #'
      write(*,*)'# -Fullea et al., 2009, G-cubed,                  #'
      write(*,*)'#  doi:10.1029/2009GC002391                       #'
      write(*,*)'#                                                 #'
      write(*,*)'# Copyright J. Fullea & J.Afonso                  #'
      write(*,*)'# modified by Ben Mather (thermal solver)         #'
      write(*,*)'#             Mattia Guerri (crustal composition) #'
      write(*,*)'#             Eldar Baykiev (parallel gravity,    #'
      write(*,*)'#                            dispersion curves)   #'
      write(*,*)'###################################################'
      write(*,*)''
      write(*,*)''
      write(*,*)''
      write(*,*)''
      write(*,*)'Grid spacing:'
      write(*,*)'-------------------------------------'
      write(*,*)'d_x (km) -->',d_x*1e-3
      write(*,*)'d_y (km) -->',d_y*1e-3
      write(*,*)'-------------------------------------'
      write(*,*)'Vertical nodes -->',N_z
      write(*,*)'-------------------------------------'
      write(*,*)'Total number of nodes -->',NOD_TOT
      write(*,*)'-------------------------------------'
      write(*,*)''
      write(*,*)'Number of layers: ',N_lay
      write(*,*)'BODY NUMBER | DENSITY | ALPHA | THERMAL CONDUCTIVITY |
     *HEAT PRODUCTION RATE | PRESSURE COEFFICIENT | GR�NEISEN PARAMETER
     *|PRESSURE DERIVATIVE OF ISOTHERMAL BULK MODULUS | ISOTHERMAL BULK
     *MODULUS| MANTLE INDEX | Vp | Vs'
      DO i=-1,N_lay
       write(*,'(A)')prop(i)%name_mat
       write(*,3)prop(i)%rho,prop(i)%alpha,prop(i)%K,
     *  prop(i)%A,prop(i)%beta,prop(i)%gmma_T,prop(i)%K_0p,prop(i)%K_T,
     *  prop(i)%litho,prop(i)%Vp,prop(i)%Vs
      ENDDO
   3   FORMAT (F6.0,E12.3,F5.1,2E12.3,2F6.2,1X,F5.1,1X,I2,F5.1,F5.1)
      write(*,*)'*************************************'

      IF(N_sublit_t>0) THEN
!       WRITE(*,*)'NUMBER OF SUBLITHOSPHERIC MANTLE BODIES: ',N_sublit
       WRITE(*,*)'SUB-LITHOSPHERIC BODIES'
       WRITE(*,*)'----------------------'
         DO i=1,N_sublit
          WRITE(*,'(A,1X,F7.2,1X,A,1X,I3)')'D_T(K)',
     *         D_T_sub(N_sublit),'litho number: ',litho_sub(N_sublit)
         ENDDO
        WRITE(*,*)'----------------------'
!        pause
      ENDIF

!      call cpu_time(time1)
      CALL MATER

      write(*,*)'MATERIAL ASSIGNATION COMPLETED'
      write(*,*)'*************************************'
      write(*,*)''
!Average gravity atracction for the thermal calculation
      g_av=g_s+d_gz*ABS(MINVAL(layers(:,:,N_lay)))*.5

      IF (temp_calc==1) THEN
       WRITE(*,*)'THERMAL CALCULATION ON'
       WRITE(*,*)'*************************'
!       CALL TEMPERATURE_3D
!       !Thermal iteration
!       DO
!        IF (MAXVAL(prop(:)%gmma_T)==0) THEN
!         WRITE(*,*)'ALL THERMAL CONDUCTIVITIES ARE CONSTANT'
!         WRITE(*,*)'SKIPPING THE K-T LOOP'
!         EXIT
!        ENDIF
!        CALL THERM_COND
!        T_old=T
!        call cpu_time(time1)
!        CALL TEMPERATURE_3D
!        call cpu_time(time2)
!         WRITE(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
!        write(*,*)'No iteration k(T), DTemperature(K),time--> ',I_TEM,
!     *MAXVAL(ABS(T-T_old)),time2-time1
!        WRITE(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
!        IF (MAXVAL(ABS(T-T_old))<0.1) THEN
!         DEALLOCATE(T_old)
!         EXIT
!        ENDIF
!        I_TEM=I_TEM+1
!       ENDDO
!
!       OPEN(55,FILE='OUT_xyz/Temperature_old.xyz', STATUS='UNKNOWN')
       WRITE(*,*)'New solver...'
!Thermal calculation requires thermal properties to be written
       OPEN(69,FILE='OUT_xyz/Material_index.xyz',STATUS='UNKNOWN')
       DO I_y=1,N_y
        y=d_y*(I_y-1)*1D-3
        DO I_x=1,N_x
         x=d_x*(I_x-1)*1D-3
         DO I_z=1,N_z
          z=(E_max-(I_z-1)*d_z)*1D-3
!          WRITE(55,'(3F13.3,F10.2)')x,y,z,T(I_x,I_y,I_z)
          WRITE(69,'(3F13.3,I6)')x,y,z,mat(I_x,I_y,I_z)
         ENDDO
        ENDDO
       ENDDO
!       CLOSE(55)
       CLOSE(69)
       CALL SYSTEM('python3 temperature_solver.py')

!Finished temperature solve
!Read temperature and conductivity files
       OPEN(80,FILE='OUT_xyz/Thermal_cond.xyz',STATUS='UNKNOWN')
       OPEN(82,FILE='OUT_xyz/Temperature.xyz',STATUS='UNKNOWN')
       DO I_y=1,N_y
        DO I_x=1,N_x
         DO I_z=1,N_z
          READ(80,*)x,y,z,k_T_P(I_x,I_y,I_z)
          READ(82,*)x,y,z,T(I_x,I_y,I_z)
         ENDDO
        ENDDO
       ENDDO
       CLOSE(80)
       CLOSE(82)
       WRITE(*,*)'min/max T', MINVAL(T), MAXVAL(T)
       !T = T+12D0 Why it is here??

      ELSE
       WRITE(*,*)'THERMAL CALCULATION OFF'
      ENDIF



      ! !!! Read an external file with density.
      ! allocate(inpdenarr(N_x,N_y,N_z))
      ! open(222,file='xyden_car_fli.xyz',status='UNKNOWN')
      ! do I_y = 1,N_y
        ! do I_x = 1,N_x
          ! do I_z = 1,N_z
          ! read(222,*) dum1, dum1, inpdenarr(I_x,I_y,I_z)
          ! enddo
        ! enddo
      ! enddo
      ! close(222)






! Calculate elevation without any sublithospheric contribution
      ! CALL COLUMNS(0)
      call columns_itherc(0)
!Calculate density, pressure,gravity and geoid anomaly considering
!sublithospheric contribution
      ! CALL COLUMNS(1)
      call columns_itherc(1)
!      call cpu_time(time2)
!      write(*,*)'Total time:', time2-time1
!Observable residuals
      ALLOCATE(D_OBS(N_x*N_y),D_calc(N_x*N_y))
       I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=E(I_x,I_y)
        D_calc(I_c)=E_calc(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, ELEVATION, Av, Std dev',AV,std_dev

      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=FA_OBS(I_x,I_y)
        D_calc(I_c)=FA(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, FA, Av, Std dev',AV,std_dev

      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=BGA_OBS(I_x,I_y)
        D_calc(I_c)=BGA(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, BOUGUER, Av, Std dev',AV,std_dev


      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=GEOID_OBS(I_x,I_y)
        D_calc(I_c)=GEOID_vec(I_c)
!!         D_calc(I_c)=GEOID(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, GEOID, Av, Std dev',AV,std_dev

      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=Uzz_OBS(I_x,I_y)
        D_calc(I_c)=Uzz(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, Uzz, Av, Std dev',AV,std_dev

      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=Uxx_OBS(I_x,I_y)
        D_calc(I_c)=Uxx(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, Uxx, Av, Std dev',AV,std_dev

      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=Uyy_OBS(I_x,I_y)
        D_calc(I_c)=Uyy(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, Uyy, Av, Std dev',AV,std_dev

      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=Uzx_OBS(I_x,I_y)
        D_calc(I_c)=Uzx(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, Uzx, Av, Std dev',AV,std_dev

      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=Uzy_OBS(I_x,I_y)
        D_calc(I_c)=Uzy(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, Uzy, Av, Std dev',AV,std_dev

      I_c=1
      DO I_x=1,N_x
       DO I_y=1,N_y
        D_OBS(I_c)=Uxy_OBS(I_x,I_y)
        D_calc(I_c)=Uxy(I_x,I_y)
        I_c=I_c+1
       ENDDO
      ENDDO
      CALL STATT(D_OBS(1:N_x*N_y),D_calc(1:N_x*N_y),N_x*N_y,AV,std_dev)
      write(*,*) 'RESIDUAL, Uxy, Av, Std dev',AV,std_dev



! EXTRACT COLUMNS LITMOD1D
      ALLOCATE (anis_mnt(N_z))

      OPEN(890,file='mineos_lonlat.xy',status='UNKNOWN')

      read(890, *) d_coordnum
      do d_coordind = 1,d_coordnum
        read(890, *) d_coordx,d_coordy,
     *     name_ray,name_love,I_anis,name_anis

        write(*,*) 'POINT: ', d_coordx, d_coordy

	      d_distmax = sqrt((L_x)**2 + (L-y)**2)
        write(*,*) 'd_distmax: ', d_distmax
	      do I_y = 1,N_y
          y = d_y*(I_y-1)
          do I_x = 1,N_x
            x = d_x*(I_x-1)
!           write(*,*) 'x, y: ', x, y
	          d_dist = sqrt((x - d_coordx)**2 + (y - d_coordy)**2)
            if (d_dist < d_distmax) then
              d_Ix = I_x
              d_Iy = I_y
              d_distmax = d_dist
            endif
          enddo
        enddo
        write(*,*) 'NEAREST POINT ON THE GRID FOUND:',
     *          d_Ix, d_Iy, d_x*(d_Ix-1), d_y*(d_Iy-1)

        write (indexstring, "(I6.6)") d_coordind
        write(*,*) indexstring

	      call read_DAT

        call LITMOD2MINEOS

        write(*,*)"mineos finished"

        IF(SW_RAY_on) THEN
          v_ray_calc=0D0
          write(*,*)"rayleighC_fund"
          CALL rayleighC_fund
          write(*,*)"rayleighC_fund finished"
          OPEN(100,file=indexstring//'_out_phase_ray.dat',STATUS='unknown')
          do i=1,N_ray
            write(100,*) v_ray_dat(i,1),v_ray_calc(i)
          enddo
          CLOSE(100)
!Misfit
          ALLOCATE(dmmy_vec(N_ray))
          dm_vec=v_ray_dat(1:N_ray,2)-v_ray_calc(1:N_ray)
          CALL MISFIT(N_ray,dm_vec(1:N_ray),v_ray_dat(1:N_ray,2),
     *          v_ray_dat(1:N_ray,3),rms,xi_sq)
          write(*,*)'Misfit Rayleigh rms (m/s), Xi squared: ', rms,xi_sq
        ENDIF

        IF(SW_LOVE_on) THEN
          v_love_calc=0D0
          write(*,*)"isoloveC_fund"
          CALL isoloveC_fund
          write(*,*)"isoloveC_fund finished"
          OPEN(101,file=indexstring//'_out_phase_love.dat',STATUS='unknown')
          do i=1,N_love
            write(101,*) v_love_dat(i,1),v_love_calc(i)
          enddo
          CLOSE(101)
!Misfit
          CALL MISFIT(N_love,v_love_dat(:,2)-v_love_calc(:),v_love_dat(:,2),
     *     v_love_dat(:,3),rms,xi_sq)
          write(*,*) 'Misfit Love rms (m/s), Xi squared: ', rms,xi_sq

        ENDIF

!deallocate variables in SUB_read_DAT
        DEALLOCATE (MDL,MDL_dm)
        DEALLOCATE(v_ray_dat,v_ray_calc)
        DEALLOCATE(v_love_dat,v_love_calc)


        DEALLOCATE(dmmy_vec)




      enddo

      close(890)



      STOP
      END
