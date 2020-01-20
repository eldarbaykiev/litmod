

!!!=============================================================================
!!!!!!!!!! columns_itherc
!!!=============================================================================


      subroutine columns_itherc(mode_col)


      use M_temperatures
      use M_grid_parameters
      use M_layers
      use M_material
      use M_columns
      use M_lsq_plane
      use M_phases
      !LK
      use M_mant
      use M_sublit
      use M_conductivity
      use M_SIGMELTS, ONLY: Cw_eps
      use m_itherc


      use M_d


      use ieee_arithmetic


      REAL(8) slp,T_int,P_int,L_lit,P_1,P_2,dens_av,a_h1,a_h2
!     REAL(8) vel_p,vel_s,ZBCV(3),rhoCV(3),wg_up,wg_dwn
      REAL(8) ZBCV(3),rhoCV(3),wg_up,wg_dwn
      REAL(8) U_sec(6),rr, GR_w,GEO_w,GR_mnt,GEO_mnt,ZB_mnt(2),
     *        ZB_mnt_g(2),rho_up_mnt,rho_down_mnt,FA_av_cr,
     *        U_sec_mnt(6),U_sec_w(6)
      INTEGER Iz1,Iz2,mat_zz,mode_col
      LOGICAL STAT,STAT1,onoff
      CHARACTER(len=255) :: path
      REAL(8) top_bot_ph(10,4),av_P_CL
      INTEGER N_mant

      integer numcom
      real(8) pre, tem
      integer check
      integer num1
      real(8) dum1

      !added by ELDAR++++++++++++++++
      integer i_counter
      integer n_elem
      integer compute_internally
      logical checknan
      real(8) xnan

      REAL(8), PARAMETER :: AA1=750, alfa3=0.26, energi=424.0E3
      REAL(8), PARAMETER :: T_zero=50
      REAL(8) melt_fac,melt_scalar,T_solidus,T_liquidus,T_sol_minus,d_T_sol
       REAL(8) melt_crit
      REAL(8) parexp, auxtop, auxbot, auxexp
      !++++++++++++++++++++++++++++++

      ! ! Datum for the second derivatives (elevation above the sea level in m).
      ! Z_usec=250D3
      compute_internally = 0

      ! Small number to avoid zero in SIGMELTS.
      Cw_eps = EPSILON(Cw_eps)


		call GETCWD(path)
		! open(14,file='OUT_xyz/vol_frac.xyz',status='UNKNOWN')
		! open(15,file='OUT_xyz/phases.xyz',status='UNKNOWN')


      open(11,file='OUT_xyz/Trhopvels.xyz',status='UNKNOWN')
      open(333,file='xoutput0.xyz',status='unknown')
      open(555,file='xoutput1.xyz',status='unknown')
      open(777,file='xoutput2.xyz',status='unknown')
      open(12,file='OUT_xyz/T_adiab.xyz',status='UNKNOWN')
      open(13,file='OUT_xyz/HF.xyz',status='UNKNOWN')
      open(14,file='OUT_xyz/T_400.xyz',status='UNKNOWN')


		if (temp_calc==1.and.mode_col==0) then


      ! write unperturbed temperatures in a file.
			do I_y = 1,N_y
				y = d_y*(I_y-1)
				do I_x = 1,N_x
					x = d_x*(I_x-1)
					do I_z = 1,N_z
						z = E_max-(I_z-1)*d_z
						write(12,'(3F13.3,F10.3)')x,y,z,T(I_x,I_y,I_z)
						if (I_z>1) then
							if(T(I_x,I_y,I_z-1)==Ts.and.T(I_x,I_y,I_z)>Ts) then
								if(mat(I_x,I_y,I_z-1)<1) then
                  if (I_x==1) mat(I_x,I_y,I_z-1) = mat(2,I_y,I_z-1)
                  if (I_x==N_x) mat(I_x,I_y,I_z-1)=mat(N_x-1,I_y,I_z-1)
                  if (I_y==1) mat(I_x,I_y,I_z-1) = mat(I_x,2,I_z-1)
						if (I_y==N_y) mat(I_x,I_y,I_z-1) = mat(I_x,N_y-1,I_z-1)
								endif
								if(mat(I_x,I_y,I_z)<-1.OR.mat(I_x,I_y,I_z-1)<-1) then
									a_h1 = prop(mat(I_x,I_y,I_z))%A
									a_h2 = prop(mat(I_x,I_y,I_z-1))%A
								else
									a_h1 = prop(mat(I_x,I_y,I_z))%A
									a_h2 = prop(mat(I_x,I_y,I_z-1))%A
								endif
								HF = (k_T_P(I_x,I_y,I_z)+k_T_P(I_x,I_y,I_z-1))*.5D0*
     * (T(I_x,I_y,I_z)-T(I_x,I_y,I_z-1))/
     * d_z+d_z*(a_h2*.5D0)
								write(13,'(2F13.3,1X,F9.5)')x,y,HF
							endif
            		endif
						if (I_z==N_z) write(14,'(2F13.3,F10.2)')x,y,T(I_x,I_y,I_z)
         		enddo
				enddo
			enddo


		elseif (temp_calc==0.and.mode_col==0) then


      ! Read unperturbed temperatures from file T_adiab.xyz.
			write(*,*)''
			write(*,*)'Reading unperturbed T from T_adiab.xyz'
			write(*,*)''


			do I_y = 1,N_y
				y = d_y*(I_y-1)
				do I_x = 1,N_x
					x = d_x*(I_x-1)
					do I_z = 1,N_z
						z = E_max-(I_z-1)*d_z
						read(12,*) dm, dm, dm, T(I_x,I_y,I_z)
						if(I_z>1) then
							if(T(I_x,I_y,I_z-1)==Ts.and.T(I_x,I_y,I_z)>Ts) then
								if(mat(I_x,I_y,I_z-1)<1) then
									if (I_x==1) mat(I_x,I_y,I_z-1)=mat(2,I_y,I_z-1)
									if (I_x==N_x) mat(I_x,I_y,I_z-1)=mat(N_x-1,I_y,I_z-1)
									if (I_y==1) mat(I_x,I_y,I_z-1)=mat(I_x,2,I_z-1)
									if (I_y==N_y) mat(I_x,I_y,I_z-1)=mat(I_x,N_y-1,I_z-1)
								endif
								if(mat(I_x,I_y,I_z)<-1.OR.mat(I_x,I_y,I_z-1)<-1) then
									a_h1 = prop(mat(I_x,I_y,I_z))%A
									a_h2 = prop(mat(I_x,I_y,I_z-1))%A
								else
									a_h1 = prop(mat(I_x,I_y,I_z))%A
									a_h2 = prop(mat(I_x,I_y,I_z-1))%A
								endif
								HF = (k_T_P(I_x,I_y,I_z)+k_T_P(I_x,I_y,I_z-1))*.5D0*
     * (T(I_x,I_y,I_z)-T(I_x,I_y,I_z-1))/
     * d_z+d_z*(a_h2*.5D0)
								write(13,'(2F13.3,1X,F9.5)')x,y,HF
           				endif
						endif
						if (I_z==N_z) write(14,'(2F13.3,F10.2)')x,y,T(I_x,I_y,I_z)
					enddo
        enddo
      enddo


      endif


		close(13)
		close(14)


		if(mode_col==0) then


      ALLOCATE (DENS(N_x,N_y,N_z),P(N_x,N_y,N_z),P_CL(N_x,N_y))
      ALLOCATE (BGA(N_x,N_y),FA(N_x,N_y),GEOID(N_x,N_y),
     * PH_TR(N_x,N_y,N_mant_max,3),E_calc(N_x,N_y))
      ALLOCATE (GEOID_OBS(N_x,N_y),FA_OBS(N_x,N_y),BGA_OBS(N_x,N_y))
      ALLOCATE(Uxx(N_x,N_y),Uyy(N_x,N_y),Uzz(N_x,N_y),Uzx(N_x,N_y),
     * Uzy(N_x,N_y),Uxy(N_x,N_y))
      ALLOCATE (GEOID_cr(N_x,N_y),FA_cr(N_x,N_y),BGA_cr(N_x,N_y))
      ALLOCATE(Uxx_cr(N_x,N_y),Uyy_cr(N_x,N_y),Uzz_cr(N_x,N_y),
     *         Uzx_cr(N_x,N_y),Uzy_cr(N_x,N_y),Uxy_cr(N_x,N_y))
      ALLOCATE(Uxx_OBS(N_x,N_y),Uyy_OBS(N_x,N_y),Uzz_OBS(N_x,N_y),
     *         Uzx_OBS(N_x,N_y),Uzy_OBS(N_x,N_y),Uxy_OBS(N_x,N_y))
      ALLOCATE(I_sub_cont(N_x,N_y),I_sub_dis(N_x,N_y,50))
      ALLOCATE (Gt(N_mant_max),Sp(N_mant_max),an(N_mant_max))

      ALLOCATE(vel_p(N_x,N_y,N_z), vel_s(N_x,N_y,N_z))

      ALLOCATE(Q_s(N_x,N_y,N_z))


      OPEN(16,file='OUT_xyz/elev_calc_c.xyz',status='UNKNOWN')
      OPEN(17,file='OUT_xyz/FA_calc_c.xyz',status='UNKNOWN')
      OPEN(18,file='OUT_xyz/Boug_calc_c.xyz',status='UNKNOWN')
      OPEN(19,file='OUT_xyz/geoid_calc_c.xyz',status='UNKNOWN')
      OPEN(24,file='OUT_xyz/U_sec_deriv_calc_c.xyz',status='UNKNOWN')
      OPEN(20,file='OUT_xyz/wadsleyite.xyz',status='UNKNOWN')
      OPEN(44,file='OUT_xyz/P_400_c.xyz',status='UNKNOWN')
      open(33,file='OUT_xyz/moho.xyz',status='UNKNOWN')
      OPEN(191,file='OUT_xyz/geoid_calc_c_1D.xyz',status='UNKNOWN')
      OPEN(50,file='OUT_xyz/FA_calc_cr_c.xyz',status='UNKNOWN')
      OPEN(51,file='OUT_xyz/Boug_calc_cr_c.xyz',status='UNKNOWN')
      OPEN(52,file='OUT_xyz/geoid_calc_cr_c.xyz',status='UNKNOWN')
      OPEN(53,file='OUT_xyz/U_sec_deriv_calc_cr_c.xyz',status='UNKNOWN')


      else


		! Perturbed thermal field due to sublithospheric mantle bodies.
      do i=1,N_sublit_t
        if(IND_SUB(i,1)==-1) cycle
          T(IND_SUB(i,1),IND_SUB(i,2),IND_SUB(i,3)) =
     * T(IND_SUB(i,1),IND_SUB(i,2),IND_SUB(i,3))+D_T_sub(I_sublit(i))
        enddo
      endif

      n_elem = 0

      DENS=0D0
      P=0D0
      FA=0D0
      Uzz=0D0
      Uxx=0D0
      Uyy=0D0
      Uzx=0D0
      Uzy=0D0
      Uxy=0D0
      BGA=0D0
      P_CL=0D0
      GEOID=0D0
      FA_OBS=0D0
      BGA_OBS=0D0
      GEOID_OBS=0D0
      DELTA=1D9
      I_sub_cont=0
      I_sub_dis=0
      PH_TR(:,:,:,1)=1000
      PH_TR(:,:,:,2)=0
      PH_TR(:,:,:,3)=0
      rho_a=prop(0)%rho
      rho_w=prop(-1)%rho
      w1=rho_a/(rho_a-rho_w)
      w2=rho_bot/(rho_bot-rho_w)
      z_bot=-E_max+(N_z-1)*d_z


		! Assumed error for pressure (Pa).
		eps = 1D5


		! key_m=1 reads thermodynamic data.
		N_TOT = N_x*N_y*N_z
		N_cont = 0
		N_mant = 0


		! Set the compensation depth 2 nodes up the 400km-discontinuty
		! (to avoid entering into the transition zone).


      if (ABS(z_bot-400D3)<d_z) z_bot = z_bot-2*d_z
      ! Elevation calibration (Mid Oceanic Ridge).
      rho_MOR = ((z_bot-E_MOR-dzc_MOR)*rho_m_MOR+dzc_MOR*rho_c_MOR)/
     * (z_bot-E_MOR)

      pi_par = z_bot-(rho_MOR*(z_bot-E_MOR)+E_MOR*rho_w)/rho_bot
      write(*,*)rho_MOR, rho_m_MOR, pi_par

! 3401.0363036303629        3408.0000000000000        24792.803145681624


		z_bot = -z_bot


      !!! Average gravity attraction.
      g_av = g_s-d_gz*z_bot * .5


      if(mode_col==0) then
      	write(*,*)''
      	write(*,*)'Calculating elevation.'
      	write(*,*)''
      else
      	write(*,*)''
      	write(*,*)'Calculating density and seismic velocities.'
      	write(*,*)''
      endif


      call cpu_time(time1)


      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ! Initialize Perple_X project.
      ! numcom = 12
      ! call initia(numcom)
      ! numcom = 45
      ! call initia_II(numcom)
      ! numcom = 99
      ! call initia_III(numcom)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ! This is necessary because somehow calling Preple_X
      ! ! interacts with the opening of the file.
      ! open(11,file='OUT_xyz/Trhopvels.xyz',status='UNKNOWN')
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do I_y = 1,N_y
        y = d_y*(I_y-1)
        do I_x = 1,N_x
          x = d_x*(I_x-1)
          rho = 0D0
          i_st = 0
          D_G = 0D0
          do I_z = 1,N_z


            ! Initialize seismic velocities.
            vel_p(I_x,I_y,I_z) = 0D0
            vel_s(I_x,I_y,I_z) = 0D0

					  z = E_max-(I_z-1)*d_z
					  g = g_s-z*d_gz


            if (mat(I_x,I_y,I_z)==-100) then ! Air node.
              if(mode_col==1) then
							  write(11,'(3F13.3,3F10.3,2F7.4,I4,A)')x,y,z,T(I_x,I_y,I_z),
     * DENS(I_x,I_y,I_z),0D0,0D0,0D0,-1,''
              endif
						  cycle
					  endif


            mat_z = mat(I_x,I_y,I_z)




            if (I_z<N_z) then
              mat_z_plus = mat(I_x,I_y,I_z+1)
					  else
						  mat_z_plus = mat(I_x,I_y,I_z)
            endif


            if (I_z>1) mat_z_minus = mat(I_x,I_y,I_z-1)
            if (mat_z==-200) then	! Water node.
              DENS(I_x,I_y,I_z) = rho_w
              P(I_x,I_y,I_z) = g*rho_w*abs(z)

              !correction copied from SUB_COLUMNS_1D
              vel_p(I_x,I_y,I_z)=1450D-3
              vel_s(I_x,I_y,I_z)=0D0

              if(mode_col==1) then
      write(11,'(3F13.3,3F10.3,2F7.4,I4,A)')x,y,z,T(I_x,I_y,I_z),
     * DENS(I_x,I_y,I_z),P(I_x,I_y,I_z)*1D-6,0D0,0D0,-1,''
              endif
              cycle
            endif

            !correction copied from SUB_COLUMNS_1D
            if (mat_z==-300) then
              !Ice node
              DENS(I_x,I_y,I_z)=910D0
              P(I_x,I_y,I_z)=g*910D0*abs(z)
              vel_p(I_x,I_y,I_z)=3.6D0
              vel_s(I_x,I_y,I_z)=1.8D0
              !!        write(*,*)'sub mater, icy'

             cycle

            endif


            if (i_st==0) then ! First node different from air or water.


              if (mat(I_x,I_y,I_z-1)==-200.OR.E(I_x,I_y)<=0) then
      P(I_x,I_y,I_z)=(rho_w*ABS(E(I_x,I_y))+prop(mat_z)%rho*
     * (E(I_x,I_y)-z))*g
                D_G=rho_w*E(I_x,I_y)**2*.5
              else
                P(I_x,I_y,I_z)=prop(mat_z)%rho*g*(E(I_x,I_y)-z)
							if (P(I_x,I_y,I_z)<0) P(I_x,I_y,I_z) = 0
							  D_G=-prop(mat_z)%rho*E(I_x,I_y)**2*.5
              endif


					     DENS(I_x,I_y,I_z) = prop(mat_z)%rho*
     * (1-prop(mat_z)%alpha*T(I_x,I_y,I_z))
					     ! rho=DENS(I_x,I_y,I_z)*(d_z*.5D0+(E(I_x,I_y)-z))
					     if (mode_col==1) then
      write(11,'(3F13.3,3F10.3,2F7.4,I4,A)')x,y,z,T(I_x,I_y,I_z),
     * DENS(I_x,I_y,I_z),P(I_x,I_y,I_z)*1D-6,prop(mat_z)%Vp,prop(mat_z)%Vs,-1,
     * ''
               endif


            else ! Any other node inside the model.


              ! Weigths for layer interfaces.
              if (mat_z/=mat_z_minus.and.mat_z<=N_lay
     * .and.mat_z_minus<=N_lay) then
                wg_up=(d_z-layers(I_x,I_y,mat_z_minus)+z)/d_z
                wg_dwn=(layers(I_x,I_y,mat_z_minus)-z)/d_z
              else
                wg_up = .5
                wg_dwn = .5
              endif


              ! Depth increment for the geoid anomaly.
              if (mat_z/=mat_z_plus) then
                d_z_geoid=d_z*.5D0+(z-layers(I_x,I_y,mat_z))
              elseif (mat_z/=mat_z_minus) then
                d_z_geoid=d_z*.5D0+(-z+layers(I_x,I_y,mat_z_minus))
              else
                d_z_geoid=d_z
              endif


              if ((prop(mat_z)%litho==0).and.(prop(mat_z)%crst==1)) then
                ! The following expression is the explicit formulation of the first order
                ! iteration for a typical eq. of state rho=rho_0*(1+alpha*dT+beta*dP)
                ! if the pressure coefficient, prop(mat_z)%beta, is <0, then it is
                ! assumed that it is referred to deph rather than pressure (i.e.
                ! rho = rho_0*(1+alpha*dT)+beta*z)
                if (mat_z/=mat_z_minus) then
      dens_av = (wg_up*DENS(I_x,I_y,I_z-1)+wg_dwn*prop(mat_z)%rho)
                else
                  dens_av = DENS(I_x,I_y,I_z-1)
                endif

                if (prop(mat_z)%beta>0) then
       P(I_x,I_y,I_z)=P(I_x,I_y,I_z-1)+g*dens_av*d_z
     * *(1+g*prop(mat_z)%rho*prop(mat_z)%beta*d_z*.5D0)-g*d_z*
     * prop(mat_z)%rho*prop(mat_z)%alpha*
     * (T(I_x,I_y,I_z)-T(I_x,I_y,I_z-1))*.5D0
       DENS(I_x,I_y,I_z)=prop(mat_z)%rho*
     * (1-prop(mat_z)%alpha*T(I_x,I_y,I_z)+prop(mat_z)%beta*
     * P(I_x,I_y,I_z))
                else
                  P(I_x,I_y,I_z)=P(I_x,I_y,I_z-1)+g*dens_av*d_z
     * -g*d_z*prop(mat_z)%rho*prop(mat_z)%alpha*
     * (T(I_x,I_y,I_z)-T(I_x,I_y,I_z-1))*.5D0
                  DENS(I_x,I_y,I_z)=prop(mat_z)%rho*
     * (1-prop(mat_z)%alpha*T(I_x,I_y,I_z))-prop(mat_z)%beta*
     * (E(I_x,I_y)-z)
                endif

                !GET VP VS FOR CRUST!!!!!!!!!! INSTER HERE LIKE IN LITMOD1D.info

                vel_p(I_x,I_y,I_z) = prop(mat_z)%Vp
                vel_s(I_x,I_y,I_z) = prop(mat_z)%Vs
              endif


              if ((prop(mat_z)%litho)>1) then
              ! Avoid compositional heterogeneities when calculating elevation.
              if (mat_z>N_lay.and.mode_col==0) mat_z=0
              P_1 = P(I_x,I_y,I_z-1)+DENS(I_x,I_y,I_z-1)*g*d_z
              call LK_interpol(T(I_x,I_y,I_z),P_1,prop(mat_z)%litho,
     * DENS(I_x,I_y,I_z),vel_p(I_x,I_y,I_z),vel_s(I_x,I_y,I_z))
              endif


      P_1=P(I_x,I_y,I_z-1)+g*(wg_up*DENS(I_x,I_y,I_z-1)+
     * wg_dwn*DENS(I_x,I_y,I_z))*d_z
      i_af=1
      do
        if (prop(mat_z)%litho==0) exit
        call LK_interpol(T(I_x,I_y,I_z),P_1,prop(mat_z)%litho,
     * DENS(I_x,I_y,I_z),vel_p(I_x,I_y,I_z),vel_s(I_x,I_y,I_z))
      ! !------------------------------
      ! ! Launch Perple_X.
      ! if(mode_col==1) then
        ! numcom = prop(mat_z)%litho
        ! pre = P_1 * 1D-5
        ! tem = T(I_x,I_y,I_z) + 273.15
        ! if (numcom==12) then
          ! call meemum4_sub2(pre,tem)
        ! elseif (numcom==45) then
          ! call meemum4_sub2_II(pre,tem)
        ! elseif (numcom==99) then
          ! call meemum4_sub2_III(pre,tem)
        ! endif
        ! DENS(I_x,I_y,I_z) = perden
        ! vel_p = pervp
        ! vel_s = pervs
        ! endif
      !------------------------------
        P_2=P(I_x,I_y,I_z-1)+g*(wg_up*DENS(I_x,I_y,I_z-1)+
     * wg_dwn*DENS(I_x,I_y,I_z))*d_z
!LK
        d_prs=ABS(P_2-P_1)
        if (d_prs<eps)then
          P(I_x,I_y,I_z)=P_2
          !!!---------------------------------------------
          !!! Read the density from a file.
            !!! if ( prop(mat_z)%litho==28 ) then
              !!! DENS(I_x,I_y,I_z) = inpdenarr(I_x,I_y,I_z)
            !!!endif
          !!!---------------------------------------------
          exit
        endif
          P_1=P_2
          i_af=i_af+1
        enddo

!+++++++phase changes+++++++++++
!       write(*,*)"energi,volexp", energi, volexp
!       write(*,*)"P(I_x, I_y, I_z)", P(I_x, I_y, I_z)
!       write(*,*)"R_uni", R_uni
!       WRITE(*,*)"T(I_x, I_y, I_z)", T(I_x, I_y, I_z)

	auxtop=energi+volexp*P(I_x, I_y, I_z)
!       write(*,*)"auxtop", auxtop
        auxbot=R_uni*1D3*(T(I_x, I_y, I_z)+273.0D0)
!       write(*,*)"auxbot", auxbot
        auxexp=-(auxtop)/(auxbot)
!       write(*,*)"auxexp", auxexp
	parexp=DEXP(auxexp)
!       write(*,*)"parexp",parexp

         Q_s(I_x, I_y, I_z)=AA1*(((T_zero/(d_size*1000.0D0))*parexp))**alfa3  !done
!       write(*,*)"Q_s", Q_s(I_x, I_y, I_z)
         Q_s(I_x, I_y, I_z)=1/Q_s(I_x, I_y, I_z)   !done
!       write(*,*)"Q_s", Q_s(I_x, I_y, I_z)



!!      Q_s(I_z)=Q_best(I_z)*.7

         IF(Q_s(I_x, I_y, I_z)>300D0)Q_s(I_x, I_y, I_z)=300D0    !done
!!         IF((Q_s>100.OR.T(I_z)<T_a).AND.-z<z_lay(N_lay)) THEN
!!          Q_s=v_lay(4,mat_z(N_z_c))+Q_slp*(T(I_z)-T(N_z_c))
!!         ENDIF

!Check if there is melt and correct Vp and Vs !////////////////////
          T_solidus=1080d0+134.2D0*(P(I_x, I_y, I_z)/1d9)
     *              -6.581D0*(P(I_x, I_y, I_z)/1d9)
     *              **2D0+0.1054D0*(P(I_x, I_y, I_z)/1d9)**3D0   !done
          T_liquidus=1762d0+57.46D0*P(I_x, I_y, I_z)/1d9
     *              -3.487D0*(P(I_x, I_y, I_z)/1d9)
     *              **2D0+0.077D0*(P(I_x, I_y, I_z)/1d9)**3D0  !done
          d_T_sol=T_solidus*(1-.98D0)  !done
!          d_T_sol=20D0
          melt_crit=7D-1  !done
!          T_sol_minus=T_solidus-d_T_sol
          T_sol_minus=T_solidus*.98D0  !done
          IF((T(I_x, I_y, I_z)>T_sol_minus).AND.
     *        (prop(mat_z)%litho/=0)) THEN
!!          IF(T(I_z)>T_solidus) THEN
           IF(T(I_x, I_y, I_z)>T_solidus) THEN
            melt_scalar=100*(T(I_x, I_y, I_z)-T_solidus)/
     *         (T_liquidus-T_solidus)+melt_crit
            IF(melt_scalar>1) THEN
             melt_fac=1
            ELSE
             melt_fac=melt_scalar
            ENDIF
            WRITE(*,*)'Melt fraction (%), depth (km): ',
     *         melt_scalar,z*1D-3,melt_fac
!After Hammond	and Humphreys 2000 JGR  valid for melt<1%
! Tubules dLVs/dF=-2.7; org cuspate relaxed: -7.9 (melt<1%)
!If melt>1% dLVs/dF=-14.5 org cuspate relaxed:
!            vel_s(I_z)=(-melt_fac*7.9*vel_s(I_z))/100D0+vel_s(I_z)
!            vel_p(I_z)=(-melt_fac*3.6*vel_p(I_z))/100D0+vel_p(I_z)
!After Chantel et al., 2016 Science
!             vel_s(I_z)=vel_s(I_z)+0.065*melt_fac**2D0-0.5565*melt_fac
!             vel_p(I_z)=vel_p(I_z)+0.07*melt_fac**2D0-0.5566*melt_fac
!            Q_s_melt=100D0/(2.4063*log(melt_fac)+5.9284)
!            IF(Q_s_melt<0) Q_s_melt=1D3
!            Q_s(I_z)=MIN(Q_s(I_z),Q_s_melt)
           ELSE
!Buffer melt zone T_solidus > T >T_solidus-d_T_sol
           melt_fac=melt_crit*(T(I_x, I_y, I_z)-T_solidus+d_T_sol)/d_T_sol

           ENDIF
          vel_s(I_x, I_y, I_z)=vel_s(I_x, I_y, I_z)+
     *                  0.065*melt_fac**2D0-0.5565*melt_fac
          vel_p(I_x, I_y, I_z)=vel_p(I_x, I_y, I_z)+
     *                  0.07*melt_fac**2D0-0.5566*melt_fac
!           vel_s(I_z)=(-melt_fac*2.7*vel_s(I_z))/100D0+vel_s(I_z)
!           vel_p(I_z)=(-melt_fac*1.2*vel_p(I_z))/100D0+vel_p(I_z)
          Q_s_melt=100D0/(2.4063*log(melt_fac)+5.9284)
          IF(Q_s_melt<0) Q_s_melt=1D3
          Q_s(I_x, I_y, I_z)=MIN(Q_s(I_x, I_y, I_z),Q_s_melt)
          ELSE
           melt_scalar=0
           melt_fac=0
          ENDIF


      endif
! Include attenuation

!          write(*,*) "vel_s", vel_s(I_x, I_y, I_z)
!          write(*,*) "Q_s", Q_s(I_x, I_y, I_z)
!          write(*,*) "alfa3", alfa3
!          write(*,*) "_", (1-2D0/(9D0*DTAN(pi*alfa3/
!     *          2.0D0)*Q_s(I_x, I_y, I_z)))
          vel_p(I_x, I_y, I_z)=vel_p(I_x, I_y, I_z)*
     *   (1-2D0/(9D0*DTAN(pi*alfa3/2.0D0)*Q_s(I_x, I_y, I_z)))
          vel_s(I_x, I_y, I_z)=vel_s(I_x, I_y, I_z)*
     *   (1-1D0/(2D0*DTAN(pi*alfa3/2.0D0)*Q_s(I_x, I_y, I_z)))
!          write(*,*) "vel_s after", vel_s(I_x, I_y, I_z)




		! Iteration P-rho
		! Compute mantle electrical resistivity.
!      if (prop(mat_z)%litho/=0)
!     * CALL ELECT_COND(T(I_x,I_y,I_z),P(I_x,I_y,I_z))
      if(mode_col==1) then
        if ((prop(mat_z)%litho==0).and.(prop(mat_z)%crst==1)) then
          !to remove stupid NaN bug
          xnan = vel_p(I_x,I_y,I_z)
          checknan = ((xnan*2d0==xnan).AND.
     * ((xnan<-epsilon(1.d0)).OR.(xnan>+epsilon(1.d0))))
          if(checknan.eqv..true.) then
            continue
          else
            write(11,'(3F13.3,3F10.3,2F7.4,I4,A)')x,y,z,T(I_x,I_y,I_z),
     * DENS(I_x,I_y,I_z),P(I_x,I_y,I_z)*1D-6,
     * vel_p(I_x,I_y,I_z),vel_s(I_x,I_y,I_z),-1,''
          endif



        else
          write(11,'(3F13.3,3F10.3,2F7.4,4E12.3,A)')x,y,z,T(I_x,I_y,I_z),
     * DENS(I_x,I_y,I_z),P(I_x,I_y,I_z)*1D-6,
     * vel_p(I_x,I_y,I_z),vel_s(I_x,I_y,I_z),res(1:4),''
        endif
		! Write Moho temparature and mantle velocities
!        if (prop(mat(I_x,I_y,I_z-1))%litho==0.and.
!     * prop(mat(I_x,I_y,I_z))%litho/=0)then
!          write(33,'(3F13.3,F10.3,2F7.4)')x,y,z,T(I_x,I_y,I_z),
!     * vel_p(I_x,I_y,I_z),vel_s(I_x,I_y,I_z)
!        endif
      endif





               i_st=i_st+1
               if (i_st==1) cycle
               ! Search for abrupt density changes in the mantle(s).
               if(prop(mat_z)%litho/=0.and.mode_col==1) then
                  if (mat(I_x,I_y,I_z-2)==-200 .OR.
     * mat(I_x,I_y,I_z-2)==-100) then
                     mat_zz=-1
                  else
                     mat_zz=mat(I_x,I_y,I_z-1)
                  endif
                  ! Search for bodies in the sublithospheric mantle
                  if((ABS(mat_zz-mat_z)>N_lay.and.mat_zz/=-1).OR.
     * (mat_zz>N_lay.and.mat_z>N_lay.and.mat_zz/=mat_z))then
                     I_sub_cont(I_x,I_y)=I_sub_cont(I_x,I_y)+1
                     I_sub_dis(I_x,I_y,I_sub_cont(I_x,I_y))=I_z
                  endif


                  if (prop(mat_zz)%litho==0) then
                     N_mant=1
                  elseif ((prop(mat_zz)%litho/=0).and.(prop(mat_zz)%crst==1)) then
                     N_mant=1
                  else
                     N_mant=N_mant+1
                  endif
                  Gt(N_mant)=0
                  Sp(N_mant)=0
                  an(N_mant)=0
                  PH_TR(I_x,I_y,N_mant,1)=mat_z
                  PH_TR(I_x,I_y,N_mant,2)=I_z
                  PH_TR(I_x,I_y,N_mant,3)=0
                  do i=1,N_tot_ph
                     if (phases(i_LK,i_near,j_near,i)(1:3)
     * =='Wad') wads_tr=.TRUE.
                     if (phases(i_LK,i_near,j_near,i)(1:2)
     * =='Gt'.OR.phases(i_LK,i_near,j_near,i)(1:4)=='CrGt')
     * Gt(N_mant)=1
                     if (phases(i_LK,i_near,j_near,i)(1:2)
     * =='Sp'.OR. phases(i_LK,i_near,j_near,i)(1:4)=='CrSp')
     * Sp(N_mant)=1
                     if (phases(i_LK,i_near,j_near,i)(1:2)
     * =='an'.OR.phases(i_LK,i_near,j_near,i)(1:2)
     * =='Pl') an(N_mant)=1
                  enddo
                  if (N_mant>1) then
                     if (Gt(N_mant)+Sp(N_mant)+an(N_mant)>1) then
                        PH_TR(I_x,I_y,N_mant,3)=1
                     elseif(((Gt(N_mant)/=Gt(N_mant-1)).OR.
     * (Sp(N_mant)/=Sp(N_mant-1)).OR.
     * (an(N_mant)/=an(N_mant-1))).and.
     * PH_TR(I_x,I_y,N_mant-1,3)==0) then
                        PH_TR(I_x,I_y,N_mant,3)=-1
                     endif
                  endif


               endif


               ! ! Percentage of thermodyn execution
               ! if (N_cont==NINT(N_tot*.1)) then
                  ! write(*,*)'**10%**'
               ! elseif (N_cont==NINT(N_tot*.3)) then
                  ! write(*,*)'****30%****'
               ! elseif (N_cont==NINT(N_tot*.5)) then
                  ! write(*,*)'********50%********'
               ! elseif (N_cont==NINT(N_tot*.9)) then
                  ! write(*,*)'***********90%************'
               ! endif
               ! N_cont=N_cont+1


               ! Check if the code enters into the transition zone.
               if (wads_tr.and.mode_col==1) then
                  write(20,'(3F13.3)')x,y,z
                  wads_tr=.FALSE.
               endif
               D_G=D_G+DENS(I_x,I_y,I_z)*ABS(z)*d_z_geoid


            enddo ! End the loop over the 'z'.


            ! Elevation calculation integrating along the whole litho-sublithospheric column.
            if(mode_col==0) then
               L_lit=layers(I_x,I_y,0)-z_bot
               if (ABS(layers(1,1,N_lay+1)+400D3)<d_z) then
                  P_CL(I_x,I_y)=P(I_x,I_y,N_z-2)
               else
                  P_CL(I_x,I_y)=P(I_x,I_y,N_z)
               endif

               if(mat_z>N_lay)then
                  write(*,*)z_bot,P_CL(I_x,I_y),L_lit
               endif
               E_calc(I_x,I_y)=L_lit-P_CL(I_x,I_y)/(rho_bot*g_av)-pi_par
               write(16,'(2F13.3,F10.2)')x,y,E_calc(I_x,I_y)
               D_G=G_ref_abs-D_G*2*pi*G_gr/g_av
               write(191,'(2F13.3,F10.2)')x,y,D_G
            endif


         !! exit !solo hace una columna. End loop on the 'x'.
         enddo


      !! exit !solo hace una columna. End loop on the 'y'.
      enddo



      call cpu_time(time2)


      ! write(*,*) 'Thermodynamic time (s)... ',time2-time1


      if(mode_col==0) then
         av_P_CL=SUM(P_CL)/(N_x*N_y)
         write(*,*)'Av. press. at the bottom of the model assuming
     * calculated elevation (GPa)',
     * av_P_CL*1D-9
         av_P_CL=-(z_bot+pi_par)*rho_bot*g_av
         write(*,*)'Av. press. at the bottom of the model assuming
     * observed elevation (GPa)',
     * av_P_CL*1D-9


         do I_y=1,N_y
            y=d_y*(I_y-1)
            do I_x=1,N_x
               x=d_x*(I_x-1)
      write(44,'(2F13.3,1X,E14.6)')x*1D-3,((N_y-1)*d_y-y)*1D-3,
     * P_CL(I_x,I_y)-av_P_CL
            enddo
         enddo


         close(16)
         close(44)
         close(12)
         close(191)


         RETURN


      else


         close(11)
         close(333)
         close(555)
         close(777)
         close(33)


         I_sub_max=MAXVAL(I_sub_cont)
         layers(:,:,N_lay+1)=z_bot
         allocate(rho_lay(N_x,N_y,0:(N_lay+I_sub_max)*2+3),
     *          z_lay(N_x,N_y,0:(N_lay+I_sub_max)*2+3))
         allocate(N_ph(N_x,N_y,(N_lay+I_sub_max)+1),
     *          z_ph(N_x,N_y,(N_lay+I_sub_max)+1,20),
     *          rho_ph(N_x,N_y,(N_lay+I_sub_max)+1,20))
      endif


      ! Set layer limits and densities for the potential field calculation.
      do I_y=1,N_y
!       y=d_y*(I_y-1)
       do I_x=1,N_x
!        x=d_x*(I_x-1)
        rho=0D0
        I_ss=I_sub_cont(I_x,I_y)
        I_sub_dis(I_x,I_y,I_sub_cont(I_x,I_y)+1)=
     *  NINT((-z_bot+E_max)/d_z+1)
        do i_lay=0,N_lay+I_ss
         if (i_lay<=N_lay) then
          z_up=layers(I_x,I_y,i_lay)
          if ((i_lay==N_lay.and.I_ss==0).OR.i_lay<N_lay) then
           z_down=layers(I_x,I_y,i_lay+1)
          elseif (i_lay==N_lay.and.I_ss>0) then
           z_down=E_max-(I_sub_dis(I_x,I_y,1)-1)*d_z
          endif
         else
           z_up=E_max-(I_sub_dis(I_x,I_y,i_lay-N_lay)-1)*d_z
           z_down=E_max-(I_sub_dis(I_x,I_y,i_lay-N_lay+1)-1)*d_z
         endif
         z_lay(I_x,I_y,2*i_lay)=z_up
         z_lay(I_x,I_y,2*i_lay+1)=z_down
      if(z_down>z_up)write(*,*)'zd>z_p',I_x,I_y,i_lay,(z_down-z_up)*1d-3
!z_up
         Iz1=INT((E_max-z_up)/d_z)+1
         Iz2=NINT((E_max-z_up)/d_z)+1
         if (Iz2==Iz1) Iz2=Iz2+1
         slp=-(T(I_x,I_y,Iz2)-T(I_x,I_y,Iz1))/d_z
         T_int=T(I_x,I_y,Iz1)+slp*(z_up+(Iz1-1)*d_z-E_max)
         slp=-(P(I_x,I_y,Iz2)-P(I_x,I_y,Iz1))/d_z
         P_int=P(I_x,I_y,Iz1)+slp*(z_up+(Iz1-1)*d_z-E_max)

         if (i_lay<=N_lay) then
!Lithosphere
          if (prop(i_lay+1)%litho/=0.OR.(i_lay==N_lay.and.I_ss==0))then
!Lithospheric mantle
           if (mat(I_x,I_y,Iz1)==i_lay+1.OR.i_lay==N_lay
     *         .and.mat(I_x,I_y,Iz1)==0) then
            rho_up=DENS(I_x,I_y,Iz1)
           else
            rho_up=DENS(I_x,I_y,Iz2)
           endif
          else
!Crust
           if (prop(i_lay+1)%beta>0) then
            rho_up=prop(i_lay+1)%rho*(1-prop(i_lay+1)%alpha*T_int+
     *          prop(i_lay+1)%beta*P_int)
           else
            rho_up=prop(i_lay+1)%rho*(1-prop(i_lay+1)%alpha*T_int)-
     *            prop(i_lay+1)%beta*(E(I_x,I_y)-z_up)
           endif
          endif
         else
!Sublithospheric bodies
           rho_up=DENS(I_x,I_y,I_sub_dis(I_x,I_y,i_lay-N_lay))
          endif

         Iz1=INT((E_max-z_down)/d_z)+1
         Iz2=NINT((E_max-z_down)/d_z)+1
         if (Iz2==Iz1) Iz1=Iz1-1
         slp=-(T(I_x,I_y,Iz2)-T(I_x,I_y,Iz1))/d_z
         T_int=T(I_x,I_y,Iz1)+slp*(z_down+(Iz1-1)*d_z-E_max)
         slp=-(P(I_x,I_y,Iz2)-P(I_x,I_y,Iz1))/d_z
         P_int=P(I_x,I_y,Iz1)+slp*(z_down+(Iz1-1)*d_z-E_max)
          if (i_lay<=N_lay) then
!Lithosphere
           if (prop(i_lay+1)%litho/=0.OR.(i_lay==N_lay.and.I_ss==0))then
!Lithospheric mantle
            if (mat(I_x,I_y,Iz1)==i_lay+1.OR.(i_lay==N_lay.and.
     *         mat(I_x,I_y,Iz1)==0)) then
             rho_down=DENS(I_x,I_y,Iz1)
            else
             rho_down=DENS(I_x,I_y,Iz2)
            endif
           else
!Crust
            if (prop(i_lay+1)%beta>0) then
            rho_down=prop(i_lay+1)%rho*(1-prop(i_lay+1)%alpha*T_int+
     *            prop(i_lay+1)%beta*P_int)
            else
            rho_down=prop(i_lay+1)%rho*(1-prop(i_lay+1)%alpha*T_int)-
     *             prop(i_lay+1)%beta*(E(I_x,I_y)-z_down)
            endif
           endif
!Sublithospheric mantle
          else
           rho_down=DENS(I_x,I_y,I_sub_dis(I_x,I_y,i_lay-N_lay+1)-1)
          endif
!Check for the density of prisms with a thickness < d_z
         if(ABS(z_up-z_down)<d_z.and.z_up/=z_down) then
          if(prop(i_lay+1)%litho/=0) then
           rho_up=DENS(I_x,I_y,NINT((E_max-z_up)/d_z)+1)
           rho_down=DENS(I_x,I_y,NINT((E_max-z_up)/d_z)+1)
          else
           rho_down=rho_up
          endif
         endif
         rho_lay(I_x,I_y,2*i_lay)=rho_up
         rho_lay(I_x,I_y,2*i_lay+1)=rho_down
!+++++++phase changes+++++++++++
         i_pha=0
         top_bot_ph=0
         onoff=.FALSE.
         do i=1,N_mant_max-1
!Skip very thin mantle layers
!because new crust is thin
!          if(ABS(z_up-z_down)<3*d_z) exit
!Select the appropriate mantle body
          if (i_lay+1==PH_TR(I_x,I_y,i,1).OR.(i_lay==N_lay.and.
     *        PH_TR(I_x,I_y,i,1)==0 )) then
!
           if (PH_TR(I_x,I_y,i,3)==0.and.PH_TR(I_x,I_y,i+1,3)==1) then
            i_pha=i_pha+1
            top_bot_ph(i_pha,1)=E_max-(PH_TR(I_x,I_y,i,2)-1)*d_z
            top_bot_ph(i_pha,3)=DENS(I_x,I_y,PH_TR(I_x,I_y,i,2))
            onoff=.TRUE.
           elseif (PH_TR(I_x,I_y,i,3)==0.and.PH_TR(I_x,I_y,i+1,3)==-1)
     *     then
            i_pha=i_pha+1
            top_bot_ph(i_pha,1)=E_max-(PH_TR(I_x,I_y,i,2)-1)*d_z
            top_bot_ph(i_pha,3)=DENS(I_x,I_y,PH_TR(I_x,I_y,i,2))
            top_bot_ph(i_pha,2)=E_max-(PH_TR(I_x,I_y,i+1,2)-1)*d_z
            top_bot_ph(i_pha,4)=DENS(I_x,I_y,PH_TR(I_x,I_y,i+1,2))
            if(top_bot_ph(i_pha,2)<z_down) then
             top_bot_ph(i_pha,2)=z_down
            endif
           endif
           if (PH_TR(I_x,I_y,i,3)==1.and.PH_TR(I_x,I_y,i+1,3)==0) then
            if(.not.onoff) then
             i_pha=i_pha+1
             top_bot_ph(i_pha,1:2)=E_max-(PH_TR(I_x,I_y,i+1,2)-1)*d_z
             top_bot_ph(i_pha,3:4)=DENS(I_x,I_y,PH_TR(I_x,I_y,i+1,2))
            else
             top_bot_ph(i_pha,2)=E_max-(PH_TR(I_x,I_y,i+1,2)-1)*d_z
             top_bot_ph(i_pha,4)=DENS(I_x,I_y,PH_TR(I_x,I_y,i+1,2))
             onoff=.FALSE.
            endif
           endif
          endif
         enddo
!if i_pha=0 no phase changes are present and z_ph reproduces z_lay
!(same thing for rho_ph and rho_lay)
!Check if multiple phases crosses the boundaries of the mantle layers

          if(i_pha>0) then
           if(top_bot_ph(i_pha,2)==0) then
            top_bot_ph(i_pha,2)=z_down
            top_bot_ph(i_pha,4)=rho_down
           endif
          endif
!          if (i_lay==3) THEN
!            write(*,*)i_pha, z_up, z_down
!          endif
          N_ph(I_x,I_y,i_lay+1)=i_pha
          z_ph(I_x,I_y,i_lay+1,1)=z_up
          z_ph(I_x,I_y,i_lay+1,2*i_pha+2)=z_down
          rho_ph(I_x,I_y,i_lay+1,1)=rho_up
          rho_ph(I_x,I_y,i_lay+1,2*i_pha+2)=rho_down
          do i=1,i_pha
           z_ph(I_x,I_y,i_lay+1,2*i)=top_bot_ph(i,1)
           z_ph(I_x,I_y,i_lay+1,2*i+1)=top_bot_ph(i,2)
           rho_ph(I_x,I_y,i_lay+1,2*i)=top_bot_ph(i,3)
           rho_ph(I_x,I_y,i_lay+1,2*i+1)=top_bot_ph(i,4)
           if (top_bot_ph(i,2)>top_bot_ph(i,1)) then
!This is for very thin mantle layers
            z_ph(I_x,I_y,i_lay+1,2*i)=top_bot_ph(i,2)
            rho_ph(I_x,I_y,i_lay+1,2*i)=top_bot_ph(i,4)
           endif
          enddo







        enddo
       enddo
      enddo


!ELDAR: POTENTIAL FIELD CALCULATION LOOP
      close(21)
      write(*,*)'*******************************************'
      write(*,*)'CALCULATING POTENTIAL FIELDS...'
      write(*,*)'*******************************************'
      call cpu_time(time1)

      open(678,file='grav_params.dat',status='UNKNOWN')
      write(678, *) N_x, N_y, d_x, d_y, N_lay, N_lay_cr, rho_w, Z_usec
      close(678)

      open(998,file='grav_E_array.dat',status='UNKNOWN')
      do i_counter=1,N_x
        write(998,*) E(i_counter,:)
      end do
      close(998)

      open(999,file='grav_input.dat',status='UNKNOWN')



      do I_x=1,N_x
       x=d_x*(I_x-1)
       do I_y=1,N_y
        y=d_y*(I_y-1)
        XB(1)=x-d_x*.5D0
        XB(2)=x+d_x*.5D0
        YB(1)=y-d_y*.5D0
        YB(2)=y+d_y*.5D0
        if (I_x==1) then
         XB(1)=XB(1)-DELTA
        elseif (I_x==N_x) then
         XB(2)=XB(2)+DELTA
        endif
        if (I_y==1) then
         YB(1)=YB(1)-DELTA
        elseif(I_y==N_y) then
         YB(2)=YB(2)+DELTA
        endif
!Loop over layers
        I_ss=I_sub_cont(I_x,I_y)
        do i_lay=0,N_lay+I_ss
        if (i_lay==0.and.E(I_x,I_y)<0D0) then
         rho_B=(rho_red-rho_w)
         ZBG(2)=E(I_x,I_y)
         ZBG(1)=0D0
         elseif (i_lay==0.and.E(I_x,I_y)>=0D0) then
         rho_B=-rho_red
         ZBG(2)=0D0
         ZBG(1)=E(I_x,I_y)
         endif
!Loop over phase changes
         do i=1,2*N_ph(I_x,I_y,i_lay+1)+1
         ZB(1)=z_ph(I_x,I_y,i_lay+1,i)
         ZB(2)=z_ph(I_x,I_y,i_lay+1,i+1)
         rho_up=rho_ph(I_x,I_y,i_lay+1,i)
         rho_down=rho_ph(I_x,I_y,i_lay+1,i+1)
!Check for curvature in the density distribution within the stability
!field
          r1=rho_up
          r2=rho_down
          Iz1=NINT((E_max-ZB(1))/d_z+1)
          Iz2=NINT((E_max-ZB(2))/d_z+1)
          Iz3=NINT((Iz2+Iz1)/2D0)
          j_CV=1
          if ((Iz2-Iz1)>6.and.ABS(DENS(I_x,I_y,Iz3)-(r1+r2)*.5)>5D0
     *        .and.prop(i_lay+1)%litho/=0) then
           ZBCV(1)=ZB(1)
           ZBCV(3)=ZB(2)
           rhoCV(1)=rho_up
           rhoCV(3)=rho_down
           ZBCV(2)=E_max-(Iz3-1)*d_z
           rhoCV(2)=DENS(I_x,I_y,Iz3)
           j_CV=2
          endif
          do j=1,j_CV
           if (j_CV==2) then
            ZB(1)=ZBCV(j)
            ZB(2)=ZBCV(j+1)
            rho_up=rhoCV(j)
            rho_down=rhoCV(j+1)
           endif
              rr=(rho_up+rho_down)/2
!---------------------------------------

!LOOP FOR COMPUTATIONAL GRID
      if (compute_internally.eq.1) then
      do I_xx=1,N_x
       XP=d_x*(I_xx-1)
       if(ZB(2)>ZB(1)) exit
        do I_yy=1,N_y
         YP=d_y*(I_yy-1)
         ZB_g(:)=ZB(:)-E(I_xx,I_yy)
         ZBG_g=ZBG-E(I_xx,I_yy)
         GR=0
         GEO=0
         GR_w=0
         GEO_w=0
         if (E(I_xx,I_yy)>0) then
      !Calculation point E>0
          call GRAV_GRAD3D(XB,YB,ZB,RHO_up,RHO_down,XP,YP,
     *                             E(I_xx,I_yy),GR)
          FA(I_xx,I_yy)=FA(I_xx,I_yy)+GR
      call GEO_GRAD3D(XB,YB,ZB_g,RHO_up,RHO_down,XP,YP,0D0,GEO)
          GEOID(I_xx,I_yy)=GEOID(I_xx,I_yy)+GEO
      call U_SECOND_DER(XB,YB,ZB,RHO_up,RHO_down,XP,YP,Z_usec,U_sec)
          Uzz(I_xx,I_yy)=Uzz(I_xx,I_yy)+U_sec(1)
          Uxx(I_xx,I_yy)=Uxx(I_xx,I_yy)+U_sec(2)
          Uyy(I_xx,I_yy)=Uyy(I_xx,I_yy)+U_sec(3)
          Uzx(I_xx,I_yy)=Uzx(I_xx,I_yy)+U_sec(4)
          Uzy(I_xx,I_yy)=Uzy(I_xx,I_yy)+U_sec(5)
          Uxy(I_xx,I_yy)=Uxy(I_xx,I_yy)+U_sec(6)

          if (i_lay==0) then
           if (ZBG(2)<0) then
      !Effect of water
            call GRAV_GRAD3D(XB,YB,ZBG,rho_w,rho_w,XP,YP,
     *                              E(I_xx,I_yy),GR_w)
            FA(I_xx,I_yy)=FA(I_xx,I_yy)+GR_w
      call GEO_GRAD3D(XB,YB,ZBG_g,rho_w,rho_w,XP,YP,0D0,
     *                               GEO_w)
             GEOID(I_xx,I_yy)=GEOID(I_xx,I_yy)+GEO_w
      call U_SECOND_DER(XB,YB,ZBG,rho_w,rho_w,XP,YP,Z_usec,U_sec_w)
          Uzz(I_xx,I_yy)=Uzz(I_xx,I_yy)+U_sec_w(1)
          Uxx(I_xx,I_yy)=Uxx(I_xx,I_yy)+U_sec_w(2)
          Uyy(I_xx,I_yy)=Uyy(I_xx,I_yy)+U_sec_w(3)
          Uzx(I_xx,I_yy)=Uzx(I_xx,I_yy)+U_sec_w(4)
          Uzy(I_xx,I_yy)=Uzy(I_xx,I_yy)+U_sec_w(5)
          Uxy(I_xx,I_yy)=Uxy(I_xx,I_yy)+U_sec_w(6)
           endif
          call GRAV_GRAD3D(XB,YB,ZBG,RHO_B,RHO_B,XP,
     *                             YP,E(I_xx,I_yy),BG)
          BGA(I_xx,I_yy)=BGA(I_xx,I_yy)+BG
          endif
      if (isnan(U_sec(2))) stop

         else
      !Calculation point E<0
         call GEOGRAV_GRAD3D(XB,YB,ZB,RHO_up,RHO_down,XP,YP,0D0,
     *                               GEO,GR)
         FA(I_xx,I_yy)=FA(I_xx,I_yy)+GR
         GEOID(I_xx,I_yy)=GEOID(I_xx,I_yy)+GEO
         call U_SECOND_DER(XB,YB,ZB,RHO_up,RHO_down,XP,YP,
     *                             Z_usec,U_sec)
         Uzz(I_xx,I_yy)=Uzz(I_xx,I_yy)+U_sec(1)
         Uxx(I_xx,I_yy)=Uxx(I_xx,I_yy)+U_sec(2)
         Uyy(I_xx,I_yy)=Uyy(I_xx,I_yy)+U_sec(3)
         Uzx(I_xx,I_yy)=Uzx(I_xx,I_yy)+U_sec(4)
         Uzy(I_xx,I_yy)=Uzy(I_xx,I_yy)+U_sec(5)
         Uxy(I_xx,I_yy)=Uxy(I_xx,I_yy)+U_sec(6)
          if (i_lay==0) then
           if (ZBG(2)<0) then
            call GEOGRAV_GRAD3D(XB,YB,ZBG,rho_w,rho_w,XP,YP,0D0,
     *                                 GEO_w,GR_w)
            FA(I_xx,I_yy)=FA(I_xx,I_yy)+GR_w
            GEOID(I_xx,I_yy)=GEOID(I_xx,I_yy)+GEO_w
            call U_SECOND_DER(XB,YB,ZBG,rho_w,rho_w,XP,YP,
     *                                Z_usec,U_sec_w)
          Uzz(I_xx,I_yy)=Uzz(I_xx,I_yy)+U_sec_w(1)
          Uxx(I_xx,I_yy)=Uxx(I_xx,I_yy)+U_sec_w(2)
          Uyy(I_xx,I_yy)=Uyy(I_xx,I_yy)+U_sec_w(3)
          Uzx(I_xx,I_yy)=Uzx(I_xx,I_yy)+U_sec_w(4)
          Uzy(I_xx,I_yy)=Uzy(I_xx,I_yy)+U_sec_w(5)
          Uxy(I_xx,I_yy)=Uxy(I_xx,I_yy)+U_sec_w(6)
           endif
          call GRAV_GRAD3D(XB,YB,ZBG,RHO_B,RHO_B,XP,
     *                             YP,0D0,BG)
          BGA(I_xx,I_yy)=BGA(I_xx,I_yy)+BG
          endif

         endif

      !write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$4
         if(prop(i_lay+1)%crst==1.and.i_lay<N_lay) then
      !Add homogenous mantle to compute the crustal gravity effect
          if (i_lay==N_lay_cr-1) then
           ZB_mnt(1)=layers(I_x,I_y,N_lay_cr)
           ZB_mnt(2)=-176.D3
           rho_up_mnt=3257.8D0
           rho_down_mnt=3257.8D0
           ZB_mnt_g(:)=ZB_mnt(:)-E(I_xx,I_yy)
           if (E(I_xx,I_yy)>0) then
            call GRAV_GRAD3D(XB,YB,ZB_mnt,RHO_up_mnt,
     *                   RHO_down_mnt,XP,YP,E(I_xx,I_yy),GR_mnt)
            call GEO_GRAD3D(XB,YB,ZB_mnt_g,RHO_up_mnt,
     *                             RHO_down_mnt,XP,YP,0D0,GEO_mnt)
           else
            call GEOGRAV_GRAD3D(XB,YB,ZB_mnt,RHO_up_mnt,
     *                                  RHO_down_mnt,XP,YP,0D0,GEO_mnt,
     *                                  GR_mnt)
           endif
           call U_SECOND_DER(XB,YB,ZB_mnt,RHO_up_mnt,
     *                              RHO_down_mnt,XP,YP,Z_usec,U_sec_mnt)

           FA_cr(I_xx,I_yy)=FA_cr(I_xx,I_yy)+GR_mnt
           GEOID_cr(I_xx,I_yy)=GEOID_cr(I_xx,I_yy)+GEO_mnt
           Uzz_cr(I_xx,I_yy)=Uzz_cr(I_xx,I_yy)+U_sec_mnt(1)
           Uxx_cr(I_xx,I_yy)=Uxx_cr(I_xx,I_yy)+U_sec_mnt(2)
           Uyy_cr(I_xx,I_yy)=Uyy_cr(I_xx,I_yy)+U_sec_mnt(3)
           Uzx_cr(I_xx,I_yy)=Uzx_cr(I_xx,I_yy)+U_sec_mnt(4)
           Uzy_cr(I_xx,I_yy)=Uzy_cr(I_xx,I_yy)+U_sec_mnt(5)
           Uxy_cr(I_xx,I_yy)=Uxy_cr(I_xx,I_yy)+U_sec_mnt(6)
          endif

         FA_cr(I_xx,I_yy)=FA_cr(I_xx,I_yy)+GR+GR_w
         GEOID_cr(I_xx,I_yy)=GEOID_cr(I_xx,I_yy)+GEO+GEO_w
         Uzz_cr(I_xx,I_yy)=Uzz_cr(I_xx,I_yy)+U_sec(1)+U_sec_w(1)
         Uxx_cr(I_xx,I_yy)=Uxx_cr(I_xx,I_yy)+U_sec(2)+U_sec_w(2)
         Uyy_cr(I_xx,I_yy)=Uyy_cr(I_xx,I_yy)+U_sec(3)+U_sec_w(3)
         Uzx_cr(I_xx,I_yy)=Uzx_cr(I_xx,I_yy)+U_sec(4)+U_sec_w(4)
         Uzy_cr(I_xx,I_yy)=Uzy_cr(I_xx,I_yy)+U_sec(5)+U_sec_w(5)
         Uxy_cr(I_xx,I_yy)=Uxy_cr(I_xx,I_yy)+U_sec(6)+U_sec_w(6)


          endif
      !write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$
       enddo
      enddo
      end if
              write(999,*) i_lay,prop(i_lay+1)%litho,prop(i_lay+1)%crst,
     *                     rho_B,XB,YB,ZB,ZBG,
     *                     RHO_up,RHO_down, layers(I_x,I_y,N_lay_cr)
              n_elem = n_elem+1

!              pause

!output
!GR=0
!GEO=0
!GR_w=0
!GEO_w=0
!CUTOUT

!---------------------------------------
          enddo
         enddo
        enddo
       enddo
      enddo
      close(999)

      if (compute_internally.eq.0) then
!        call system("pwd")
        call system("python3 gravity_calculator.py")


        open(200,file='grav_input_BGA.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(200,*) BGA(i_counter,:)
        end do
        close(200)

        open(201,file='grav_input_FA.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(201,*) FA(i_counter,:)
        end do
        close(201)

        open(202,file='grav_input_GEOID.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(202,*) GEOID(i_counter,:)
        end do
        close(202)

        open(203,file='grav_input_Uxx.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(203,*) Uxx(i_counter,:)
        end do
        close(203)

        open(204,file='grav_input_Uyy.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(204,*) Uyy(i_counter,:)
        end do
        close(204)

        open(205,file='grav_input_Uzz.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(205,*) Uzz(i_counter,:)
        end do
        close(205)

        open(206,file='grav_input_Uzx.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(206,*) Uzx(i_counter,:)
        end do
        close(206)

        open(207,file='grav_input_Uzy.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(207,*) Uzy(i_counter,:)
        end do
        close(207)

        open(208,file='grav_input_Uxy.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(208,*) Uxy(i_counter,:)
        end do
        close(208)


        open(300,file='grav_input_BGA_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(300,*) BGA_cr(i_counter,:)
        end do
        close(300)

        open(301,file='grav_input_FA_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(301,*) FA_cr(i_counter,:)
        end do
        close(301)

        open(302,file='grav_input_GEOID_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(302,*) GEOID_cr(i_counter,:)
        end do
        close(302)

        open(303,file='grav_input_Uxx_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(303,*) Uxx_cr(i_counter,:)
        end do
        close(303)

        open(304,file='grav_input_Uyy.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(304,*) Uyy_cr(i_counter,:)
        end do
        close(304)

        open(305,file='grav_input_Uzz_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(305,*) Uzz_cr(i_counter,:)
        end do
        close(305)

        open(306,file='grav_input_Uzx_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(306,*) Uzx_cr(i_counter,:)
        end do
        close(306)

        open(307,file='grav_input_Uzy_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(307,*) Uzy_cr(i_counter,:)
        end do
        close(307)

        open(308,file='grav_input_Uxy_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          read(308,*) Uxy_cr(i_counter,:)
        end do
        close(308)

      else
        open(200,file='grav_BGA.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(200,*) BGA(i_counter,:)
        end do
        close(200)

        open(201,file='grav_FA.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(201,*) FA(i_counter,:)
        end do
        close(201)

        open(202,file='grav_GEOID.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(202,*) GEOID(i_counter,:)
        end do
        close(202)

        open(203,file='grav_Uxx.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(203,*) Uxx(i_counter,:)
        end do
        close(203)

        open(204,file='grav_Uyy.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(204,*) Uyy(i_counter,:)
        end do
        close(204)

        open(205,file='grav_Uzz.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(205,*) Uzz(i_counter,:)
        end do
        close(205)

        open(206,file='grav_Uzx.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(206,*) Uzx(i_counter,:)
        end do
        close(206)

        open(207,file='grav_Uzy.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(207,*) Uzy(i_counter,:)
        end do
        close(207)

        open(208,file='grav_Uxy.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(208,*) Uxy(i_counter,:)
        end do
        close(208)


        open(300,file='grav_BGA_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(300,*) BGA_cr(i_counter,:)
        end do
        close(300)

        open(301,file='grav_FA_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(301,*) FA_cr(i_counter,:)
        end do
        close(301)

        open(302,file='grav_GEOID_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(302,*) GEOID_cr(i_counter,:)
        end do
        close(302)

        open(303,file='grav_Uxx_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(303,*) Uxx_cr(i_counter,:)
        end do
        close(303)

        open(304,file='grav_Uyy_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(304,*) Uyy_cr(i_counter,:)
        end do
        close(304)

        open(305,file='grav_Uzz_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(305,*) Uzz_cr(i_counter,:)
        end do
        close(305)

        open(306,file='grav_Uzx_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(306,*) Uzx_cr(i_counter,:)
        end do
        close(306)

        open(307,file='grav_Uzy_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(307,*) Uzy_cr(i_counter,:)
        end do
        close(307)

        open(308,file='grav_Uxy_cr.dat',status='UNKNOWN')
        do i_counter=1,N_x
          write(308,*) Uxy_cr(i_counter,:)
        end do
        close(308)
      endif
!      resulta=MAXLOC(FA)
!Reference model
      E_ref=0D0
      z_cref=-30.0D3
      z_lref=-176.D3
      rho_cref=2780D0
      rho_mref=3257.8D0
      rho_sublref=3573.93D0
      XB(1)=-d_x*.5-DELTA
      XB(2)=L_x+d_x*.5+DELTA
      YB(1)=-d_y*.5-DELTA
      YB(2)=L_y+d_y*.5+DELTA
      if (1==1) then
      write(*,*)'REFERENCE MODEL SUBSTRACTED'
      do I_x=1,N_x
       xp=d_x*(I_x-1)
       do I_y=1,N_y
        yp=d_y*(I_y-1)
         do il=1,3
         if (IL==1) then
         ZB(1)=0D0
         ZB(2)=z_cref
         rho_ref=rho_cref
         elseif (IL==2) then
         ZB(2)=z_lref
         ZB(1)=z_cref
         rho_ref=rho_mref
         else
         ZB(1)=z_lref
         ZB(2)=layers(1,1,N_lay+1)
         rho_ref=rho_sublref
         endif
         ZB_g=ZB-E(I_x,I_y)
         if (E(I_x,I_y)>0) then
           call GRAV_GRAD3D(XB,YB,ZB,rho_ref,rho_ref,XP,YP,
     *                      E(I_x,I_y),GR)
           FA(I_x,I_y)=FA(I_x,I_y)-GR

           call GEO_GRAD3D(XB,YB,ZB_g,rho_ref,rho_ref,XP,YP,0D0,
     *                     GEO)
           GEOID(I_x,I_y)=GEOID(I_x,I_y)-GEO
          else
           call GEOGRAV_GRAD3D(XB,YB,ZB,rho_ref,rho_ref,XP,YP,0D0,
     *                         GEO,GR)
           FA(I_x,I_y)=FA(I_x,I_y)-GR
           GEOID(I_x,I_y)=GEOID(I_x,I_y)-GEO
          endif

          call U_SECOND_DER(XB,YB,ZB,rho_ref,rho_ref,XP,YP,
     *                      Z_usec,U_sec)
          Uzz(I_x,I_y)=Uzz(I_x,I_y)-U_sec(1)
          Uxx(I_x,I_y)=Uxx(I_x,I_y)-U_sec(2)
          Uyy(I_x,I_y)=Uyy(I_x,I_y)-U_sec(3)
          Uzx(I_x,I_y)=Uzx(I_x,I_y)-U_sec(4)
          Uzy(I_x,I_y)=Uzy(I_x,I_y)-U_sec(5)
          Uxy(I_x,I_y)=Uxy(I_x,I_y)-U_sec(6)

!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$4
          if (IL<3) then
           FA_cr(I_x,I_y)=FA_cr(I_x,I_y)-GR
           GEOID_cr(I_x,I_y)=GEOID_cr(I_x,I_y)-GEO
           Uzz_cr(I_x,I_y)=Uzz_cr(I_x,I_y)-U_sec(1)
           Uxx_cr(I_x,I_y)=Uxx_cr(I_x,I_y)-U_sec(2)
           Uyy_cr(I_x,I_y)=Uyy_cr(I_x,I_y)-U_sec(3)
           Uzx_cr(I_x,I_y)=Uzx_cr(I_x,I_y)-U_sec(4)
           Uzy_cr(I_x,I_y)=Uzy_cr(I_x,I_y)-U_sec(5)
           Uxy_cr(I_x,I_y)=Uxy_cr(I_x,I_y)-U_sec(6)
          endif
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$4
        enddo
       enddo
      enddo
      endif

      call cpu_time(time2)

      write(*,*) 'Potential field time (s)... ',time2-time1
      write(*,*) 'DATUM FOR THE GRAVITY GRADIENTS (km) ',Z_usec*1D-3
! read measured observables, if needed
      call CHDIR (TRIM(path))
      INQUIRE(FILE='GEO_DATA/geoid_dat_c.xyz',EXIST=STAT)
      INQUIRE(FILE='GEO_DATA/FA_dat_c.xyz',EXIST=STAT1)
      STAT=STAT.and.STAT1
      INQUIRE(FILE='GEO_DATA/Boug_dat_c.xyz',EXIST=STAT1)
      STAT=STAT.and.STAT1

      if (STAT) then
       write(*,*)'readING OBSERVED GEOPHYSICAL DATA'
       write(*,*)'(GEOID, FA and BOUGUER ANOMALIES)'
       open(21,file='GEO_DATA/geoid_dat_c.xyz',status='OLD')
       open(22,file='GEO_DATA/FA_dat_c.xyz',status='OLD')
       open(23,file='GEO_DATA/Boug_dat_c.xyz',status='OLD')
       do I_y=1,N_y
        do I_x=1,N_x
         read(23,*) q,q,BGA_OBS(I_x,I_y)
         read(22,*) q,q,FA_OBS(I_x,I_y)
         read(21,*) q,q,GEOID_OBS(I_x,I_y)
        enddo
       enddo
       close(21)
       close(22)
       close(23)
      endif
!Gravity gradients
      INQUIRE(FILE='GEO_DATA/U_sec_deriv_dat_c.xyz',EXIST=STAT1)
      if (STAT1) then
       write(*,*)'readING OBSERVED GEOPHYSICAL DATA'
       write(*,*)'(GRAVITY GRADIENTS Uzz Uxx Uyy Uzx Uzy Uxy (MRF))'
       open(21,file='GEO_DATA/U_sec_deriv_dat_c.xyz',status='OLD')
       do I_y=1,N_y
        do I_x=1,N_x
         read(21,*) q,q,Uzz_OBS(I_x,I_y),Uxx_OBS(I_x,I_y),
     *              Uyy_OBS(I_x,I_y),Uzx_OBS(I_x,I_y),Uzy_OBS(I_x,I_y),
     *              Uxy_OBS(I_x,I_y)
        enddo
       enddo
       close(21)
      endif

      write(*,*)'MAX MIN BOUGUER ANOMALY DATA',
     *maxval(BGA_OBS),minval(BGA_OBS)
      write(*,*)'MAX MIN FA ANOMALY DATA',maxval(FA_OBS),minval(FA_OBS)
      write(*,*)'MAX MIN GEOID ANOMALY DATA',
     *maxval(GEOID_OBS),minval(GEOID_OBS)
      write(*,*)'MAX MIN Uzz DATA',maxval(Uzz_OBS),minval(Uzz_OBS)
      write(*,*)'MAX MIN Uxx DATA',maxval(Uxx_OBS),minval(Uxx_OBS)
      write(*,*)'MAX MIN Uyy DATA',maxval(Uyy_OBS),minval(Uyy_OBS)
      write(*,*)'MAX MIN Uzx DATA',maxval(Uzx_OBS),minval(Uzx_OBS)
      write(*,*)'MAX MIN Uzy DATA',maxval(Uzy_OBS),minval(Uzy_OBS)
      write(*,*)'MAX MIN Uxy DATA',maxval(Uxy_OBS),minval(Uxy_OBS)

      if (STAT) then
       FA_av=(sum(FA)-sum(FA_OBS))/(N_x*N_y)
       FA_av_cr=(sum(FA_cr)-sum(FA_OBS))/(N_x*N_y)
      else
       FA_av=sum(FA)/(N_x*N_y)
      endif
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$4
      FA_cr=FA_cr-FA_av_cr
      BGA_cr=FA_cr+BGA
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$4
      FA=FA-FA_av
      BGA=FA+BGA
      if (STAT) then
       BGA_av=(sum(BGA)-sum(BGA_OBS))/(N_x*N_y)
      else
       BGA_av=sum(BGA)/(N_x*N_y)
      endif
       call LSQ_PLANE
! Uzy and Uxy have a explicit sign change to
! account for the reference frame implicitly assumed by LitMod in which
! x--> E and y-->S (with respect to the LNOF,x-->N and y-->W,
! in which GOCE gravity gradients are usually provided)
!       Uzx=-Uzx
       Uzy=-Uzy
       Uxy=-Uxy
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$
       Uzy_cr=-Uzy_cr
       Uxy_cr=-Uxy_cr
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$

       if (STAT1) then
        Uzz=Uzz-(sum(Uzz)-sum(Uzz_OBS))/(N_x*N_y)
        Uxx=Uxx-(sum(Uxx)-sum(Uxx_OBS))/(N_x*N_y)
        Uyy=Uyy-(sum(Uyy)-sum(Uyy_OBS))/(N_x*N_y)
        Uzx=Uzx-(sum(Uzx)-sum(Uzx_OBS))/(N_x*N_y)
        Uzy=Uzy-(sum(Uzy)-sum(Uzy_OBS))/(N_x*N_y)
        Uxy=Uxy-(sum(Uxy)-sum(Uxy_OBS))/(N_x*N_y)
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$
        Uzz_cr=Uzz_cr-(sum(Uzz_cr)-sum(Uzz_OBS))/(N_x*N_y)
        Uxx_cr=Uxx_cr-(sum(Uxx_cr)-sum(Uxx_OBS))/(N_x*N_y)
        Uyy_cr=Uyy_cr-(sum(Uyy_cr)-sum(Uyy_OBS))/(N_x*N_y)
        Uzx_cr=Uzx_cr-(sum(Uzx_cr)-sum(Uzx_OBS))/(N_x*N_y)
        Uzy_cr=Uzy_cr-(sum(Uzy_cr)-sum(Uzy_OBS))/(N_x*N_y)
        Uxy_cr=Uxy_cr-(sum(Uxy_cr)-sum(Uxy_OBS))/(N_x*N_y)
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$
       else
        Uzz=Uzz-sum(Uzz)/(N_x*N_y)
        Uxx=Uxx-sum(Uxx)/(N_x*N_y)
        Uyy=Uyy-sum(Uyy)/(N_x*N_y)
        Uzx=Uzx-sum(Uzx)/(N_x*N_y)
        Uzy=Uzy-sum(Uzy)/(N_x*N_y)
        Uxy=Uxy-sum(Uxy)/(N_x*N_y)
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$
        Uzz_cr=Uzz_cr-(sum(Uzz_cr))/(N_x*N_y)
        Uxx_cr=Uxx_cr-(sum(Uxx_cr))/(N_x*N_y)
        Uyy_cr=Uyy_cr-(sum(Uyy_cr))/(N_x*N_y)
        Uzx_cr=Uzx_cr-(sum(Uzx_cr))/(N_x*N_y)
        Uzy_cr=Uzy_cr-(sum(Uzy_cr))/(N_x*N_y)
        Uxy_cr=Uxy_cr-(sum(Uxy_cr))/(N_x*N_y)
!write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$

       endif
      write(*,*)''
      write(*,*)'**********************************************'
      write(*,*)'FREE AIR, BOUGUER and GEOID ANOMALY CALCULATED'
      write(*,*)'**********************************************'
      I_c=1
      do I_x=1,N_x
       x=d_x*(I_x-1)
       do I_y=1,N_y
        y=d_y*(I_y-1)
         write(17,'(2F13.3,F14.2)')x,y,FA(I_x,I_y)
         write(18,'(2F13.3,F14.2)')x,y,BGA(I_x,I_y)
         write(19,'(2F13.3,F14.3)')x,y,GEOID_vec(I_c)
         write(24,'(2F13.3,6E12.4)')x,y,Uzz(I_x,I_y),Uxx(I_x,I_y),
     *                              Uyy(I_x,I_y),Uzx(I_x,I_y),
     *                              Uzy(I_x,I_y),Uxy(I_x,I_y)
        write(50,'(2F13.3,F14.2)')x,y,FA_cr(I_x,I_y)
        write(51,'(2F13.3,F14.2)')x,y,BGA_cr(I_x,I_y)
        write(52,'(2F13.3,F14.3)')x,y,GEOID_vec_cr(I_c)
        write(53,'(2F13.3,6E12.4)')x,y,Uzz_cr(I_x,I_y),Uxx_cr(I_x,I_y),
     *                              Uyy_cr(I_x,I_y),Uzx_cr(I_x,I_y),
     *                              Uzy_cr(I_x,I_y),Uxy_cr(I_x,I_y)
         I_c=I_c+1
       enddo
      enddo
      close(17)
      close(18)
      close(19)
      close(24)
      close(50)
      close(51)
      close(52)
      close(53)
		close(0001)
      DEallocate (rho_lay,z_lay)

      end SUBROUTINE
