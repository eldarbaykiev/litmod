!Version for LITMOD_1D_MINEOS

      SUBROUTINE LITMOD2MINEOS
      USE M_surf_waves
!     MODULE M_surf_waves
!      INTEGER(4) N_ray,N_love,I_anis
!      CHARACTER(300)name_ray,name_love
!      REAL(8),ALLOCATABLE::v_ray_dat(:,:),v_ray_calc(:),
!    *                      v_love_dat(:,:),v_love_calc(:),srf_wv(:,:)
!      REAL(8),ALLOCATABLE::v_ray_calc_best(:),v_love_calc_best(:)
!      LOGICAL MT_on
!      LOGICAL SW_RAY_on,SW_LOVE_on
!      REAL(8) r_mineos(100),rho_mineos(100),vpv_mineos(100),
!    *         vsv_mineos(100),qkappa_mineos(100),qshear_mineos(100),
!    *         vph_mineos(100),vsh_mineos(100)
!      REAL(8) R_L2M(300),rho_L2M(300),vpv_L2M(300),vsv_L2M(300),
!    *         Q_kappa_L2M(300),Q_s_L2M(300),vph_L2M(300),vsh_L2M(300)
!      REAL(8) z_anis(11),anis(11)
!     INTEGER N_L2M
!     REAL(8)vs_660,vs_410,an_660,an_410,anis_tran
!     END MODULE


      USE M_grid_parameters
      USE M_column_1D
!     MODULE M_column_1D
!     REAL(8),ALLOCATABLE:: anis_mnt(:),MDL(:,:),MDL_dm(:,:)
!     END MODULE

      USE M_columns
!     module M_columns
!     real(8), allocatable:: DENS(:,:,:),BGA(:,:),FA(:,:),GEOID(:,:),
!    * Uxx(:,:),Uyy(:,:),Uzz(:,:),Uzx(:,:),Uzy(:,:),Uxy(:,:)
!     real(8), allocatable:: BGA_cr(:,:),FA_cr(:,:),GEOID_cr(:,:),
!    * Uxx_cr(:,:),Uyy_cr(:,:),Uzz_cr(:,:),Uzx_cr(:,:),Uzy_cr(:,:),
!    * Uxy_cr(:,:)
!     real(8), allocatable:: GEOID_OBS(:,:),FA_OBS(:,:),E_calc(:,:)
!     real(8), allocatable:: BGA_OBS(:,:),P(:,:,:),Uxx_OBS(:,:)
!    *,Uyy_OBS(:,:),Uzz_OBS(:,:),Uzx_OBS(:,:),Uzy_OBS(:,:),Uxy_OBS(:,:)
!     real(8), allocatable:: rho_lay(:,:,:),z_lay(:,:,:)
!     real(8), allocatable:: vel_p(:,:,:),vel_s(:,:,:)
!     INTEGER(4),ALLOCATABLE::PH_TR(:,:,:,:),N_ph(:,:,:),
!    *                        Gt(:),Sp(:),an(:)
!     REAL(8),ALLOCATABLE::z_ph(:,:,:,:),rho_ph(:,:,:,:),P_CL(:,:)
!     REAL(8),ALLOCATABLE:: D_OBS(:),D_calc(:)
!     integer(2) mat_z,mat_z_plus,mat_z_minus
!     integer(4) i_st,I_z
!     real(8) rho,rho_L,AV,std_dev,Z_usec
!     real(8) z_bot,rho_up,rho_down,z_up,z_down,rho_a,rho_red,rho_B
!     real(8) w1,w2,rho_w,L,z_1,rho_1,xp,yp,GR,GEO,BG,BGA_av,FA_av,
!    *        U_sec_av
!     real(8) rho_ref,ra,RHO_AV,DELTA
!     real(8) E_ref,z_cref,z_lref,rho_cref,rho_mref,rho_sublref
!     real(8),dimension(2):: XB,YB,ZB,ZBG,ZB_g,ZBG_g
!     real(8),parameter:: g_s=9.81D0,g_400=9.97D0
!     REAL(8) g,d_gz,g_av
!     REAL(8),PARAMETER:: E_MOR=2.6D3,dzc_MOR=6.5D3
!!      REAL(8),PARAMETER:: rho_m_MOR=3420.D0,rho_c_MOR=2930D0,
!     REAL(8) rho_m_MOR
!     INTEGER E_cali
!     REAL(8),PARAMETER:: rho_c_MOR=2986D0,rho_bot=3611.3D0
!     REAL(8) rho_MOR,pi_par
!     REAL(8),PARAMETER::G_gr=6.67D-11,pi=3.14159,
!    *                   G_ref_abs=11912.4
!     end module

      use M_layers
      use M_material

      use M_d

      INTEGER rtio_up,rtio_dw
      REAL(8) an1,an2,slp_a,z1,z2
      INTEGER I_z1,I_z2,I_zz,index_ani


! Deep part of Earth model z> lower mantle
!que pasa si ianis=0!!!!!!!!!!!!!!!!!!!!!!!!
        R_L2M(1:100)=r_mineos(1:100)		!M_surf_waves
        rho_L2M(1:100)=rho_mineos(1:100)	!M_surf_waves
        vpv_L2M(1:100)=vpv_mineos(1:100)	!M_surf_waves
        vsv_L2M(1:87)=vsv_mineos(1:87)		!M_surf_waves
        vph_L2M(1:100)=vph_mineos(1:100)	!M_surf_waves
        vsh_L2M(1:87)=vsh_mineos(1:87)		!M_surf_waves
        Q_kappa_L2M(1:100)=qkappa_mineos(1:100)	!M_surf_waves
        Q_s_L2M(1:100)=qshear_mineos(1:100)	!M_surf_waves
!Transition zone in the mantle
        vs_660=5569.85				!M_surf_waves
        an_660=0				!M_surf_waves
        vs_410=5041.52				!M_surf_waves
        an_410=0				!M_surf_waves
        vsv_L2M(88)=vs_660*(1-an_660/3D0)	!M_surf_waves
        vsh_L2M(88)=vs_660*(1+2D0*an_660/3D0)	!M_surf_waves

        vsv_L2M(100)=vs_410*(1-an_410/3D0)	!M_surf_waves
        vsh_L2M(100)=vs_410*(1+2D0*an_410/3D0)	!M_surf_waves
!Linearly interpolate in the mantle transition zone
        DO i=89,99
          vsv_L2M(i)=vsv_L2M(88)-(vsv_L2M(88)-vsv_L2M(100))*
     *               (R_L2M(i)-R_L2M(88))/(250D3)			!M_surf_waves
          vsh_L2M(i)=vsh_L2M(88)-(vsh_L2M(88)-vsh_L2M(100))*
     *               (R_L2M(i)-R_L2M(88))/(250D3)			!M_surf_waves
!      write(*,*)'i,vsv,vsh',i,vsv_L2M(i), vsh_L2M(i)
        ENDDO
!Upper mantle discontinuity with LITMOD properties
        R_L2M(101)=r_mineos(100)					!M_surf_waves
        rho_L2M(101)=DENS(d_Ix, d_Iy, N_z)				!DONE
        vpv_L2M(101)=vel_p(d_Ix, d_Iy, N_z)*1D3					!M_surf_waves
        vsv_L2M(101)=vel_s(d_Ix, d_Iy, N_z)*1D3*(1-anis_mnt(N_z)/3D0)		!M_surf_waves M_column_1D
        vph_L2M(101)=vel_p(d_Ix, d_Iy, N_z)*1D3					!M_surf_waves
        vsh_L2M(101)=vel_s(d_Ix, d_Iy, N_z)*1D3*
     *          (1+2D0*anis_mnt(N_z)/3D0)		!M_surf_waves M_column_1D
        Q_kappa_L2M(101)=qkappa_mineos(100)				!M_surf_waves
        Q_s_L2M(101)=Q_s(d_Ix, d_Iy, N_z)						!M_surf_waves  M_grid_parameters


!       WRITE(*,*)"101", R_L2M(101),rho_L2M(101),
!     *             vpv_L2M(101),vsv_L2M(101),
!     *             Q_kappa_L2M(101),Q_s_L2M(101),
!     *             (-R_L2M(101)+6371D3)*1d-3
!LitMod section, mantle part
       rtio_up=NINT(2D3/d_z)						!M_grid_parameters
       rtio_dw=NINT(10D3/d_z)						!M_grid_parameters
       ii=0
!Mantle loop for z>120 km   KEEEP
       write(*,*) 'E_max',E_max
       write(*,*) 'N_moho',N_moho

!DO LOOP until 120 km
       DO i=N_z,N_moho,-rtio_dw
        z=E_max-(i-1)*d_z
        IF (z>-120D3) EXIT
        ii=ii+1
        i_L2M=101+ii
!      write(*,*)'i,ii,i_L2M 1 loop dentro ',i,i_L2M,z
        R_L2M(i_L2M)=6371D3+z
        rho_L2M(i_L2M)=DENS(d_Ix, d_Iy, i)
        vpv_L2M(i_L2M)=vel_p(d_Ix, d_Iy, i)*1D3
        vsv_L2M(i_L2M)=vel_s(d_Ix, d_Iy, i)*1D3*(1-anis_mnt(i)/3D0)
        vph_L2M(i_L2M)=vel_p(d_Ix, d_Iy, i)*1D3
        vsh_L2M(i_L2M)=vel_s(d_Ix, d_Iy, i)*1D3*(1+2D0*anis_mnt(i)/3D0)
        Q_kappa_L2M(i_L2M)=qkappa_mineos(1)
        Q_s_L2M(i_L2M)=Q_s(d_Ix, d_Iy, i)
!        WRITE(*,*)"l141",i_L2M, R_L2M(i_L2M),rho_L2M(i_L2M),
!     *             vpv_L2M(i_L2M),vsv_L2M(i_L2M),
!     *             Q_kappa_L2M(i_L2M),Q_s_L2M(i_L2M),
!     *             (-R_L2M(i_L2M)+6371D3)*1d-3
       ENDDO     !END KEEP
!      write(*,*)'(i_L2 1 loop',i_L2M
!      pause
!Mantle loop for z_moho< z <120 km  should go to the top
!       write(*,*)'E_max', E_max
!       write(*,*)'E(d_Ix, d_Iy)', E(d_Ix, d_Iy)
       DO i=N_z,1,-rtio_up
         z=E_max-(i-1)*d_z
!        write(*,*)"z", z
         IF (z<-120D3) CYCLE
         IF (z>E(d_Ix, d_Iy)) CYCLE

         IF (prop(mat(d_Ix,d_Iy,i))%litho.NE.0) THEN
           i_L2M=i_L2M+1
!      write(*,*)'i,ii,i_L2M 2 loop dentro ',i,i_L2M,z
           R_L2M(i_L2M)=6371D3+z
           rho_L2M(i_L2M)=DENS(d_Ix, d_Iy, i)
           vpv_L2M(i_L2M)=vel_p(d_Ix, d_Iy, i)*1D3 !m/s
           vsv_L2M(i_L2M)=vel_s(d_Ix, d_Iy, i)*1D3*(1-anis_mnt(i)/3D0)
           vph_L2M(i_L2M)=vel_p(d_Ix, d_Iy, i)*1D3
           vsh_L2M(i_L2M)=vel_s(d_Ix, d_Iy, i)*1D3*(1+2D0*anis_mnt(i)/3D0)
           Q_kappa_L2M(i_L2M)=qkappa_mineos(1)
           Q_s_L2M(i_L2M)=Q_s(d_Ix, d_Iy, i)
!           WRITE(*,*)"m1",i, z, i_L2M, R_L2M(i_L2M),rho_L2M(i_L2M),
!     *             vpv_L2M(i_L2M),vsv_L2M(i_L2M),
!     *             Q_kappa_L2M(i_L2M),Q_s_L2M(i_L2M),
!     *             (-R_L2M(i_L2M)+6371D3)*1d-3
!Check for very sttep gradients
           IF((ABS(vsv_L2M(i_L2M)-vsv_L2M(i_L2M-1))>100D0.OR.
     *       ABS(vsh_L2M(i_L2M)-vsh_L2M(i_L2M-1))>100D0).AND.
     *       i_L2M<213 ) THEN
!            write(*,*)'steep!',vel_s(i),vsv_L2M(i_L2M),vsv_L2M(i_L2M-1)
!            R_L2M(i_L2M)=R_L2M(i_L2M-1)

             R_L2M(i_L2M+1)=6371D3+z+d_z
             rho_L2M(i_L2M+1)=rho_L2M(i_L2M)
             vpv_L2M(i_L2M+1)=vpv_L2M(i_L2M)
             vsv_L2M(i_L2M+1)=vsv_L2M(i_L2M)
             vph_L2M(i_L2M+1)=vph_L2M(i_L2M)
             vsh_L2M(i_L2M+1)=vsh_L2M(i_L2M)
             Q_kappa_L2M(i_L2M+1)=Q_kappa_L2M(i_L2M)
             Q_s_L2M(i_L2M+1)=Q_s_L2M(i_L2M)
!             WRITE(*,*)"mg",i, z, i_L2M+1, R_L2M(i_L2M+1),rho_L2M(i_L2M+1),
!     *             vpv_L2M(i_L2M+1),vsv_L2M(i_L2M+1),
!     *             Q_kappa_L2M(i_L2M+1),Q_s_L2M(i_L2M+1),
!     *             (-R_L2M(i_L2M+1)+6371D3)*1d-3
             i_L2M=i_L2M+1
           ENDIF

           IF(prop(mat(d_Ix,d_Iy,i))
     *         %litho.NE.prop(mat(d_Ix,d_Iy,i-1))%litho) THEN
!            write(*,*)'steep!',vel_s(i),vsv_L2M(i_L2M),vsv_L2M(i_L2M-1)
!            R_L2M(i_L2M)=R_L2M(i_L2M-1)

             R_L2M(i_L2M+1)=6371D3+z+d_z
             rho_L2M(i_L2M+1)=rho_L2M(i_L2M)
             vpv_L2M(i_L2M+1)=vpv_L2M(i_L2M)
             vsv_L2M(i_L2M+1)=vsv_L2M(i_L2M)
             vph_L2M(i_L2M+1)=vph_L2M(i_L2M)
             vsh_L2M(i_L2M+1)=vsh_L2M(i_L2M)
             Q_kappa_L2M(i_L2M+1)=Q_kappa_L2M(i_L2M)
             Q_s_L2M(i_L2M+1)=Q_s_L2M(i_L2M)
!             WRITE(*,*)"m2",i, z, i_L2M+1, R_L2M(i_L2M+1),rho_L2M(i_L2M+1),
!     *             vpv_L2M(i_L2M+1),vsv_L2M(i_L2M+1),
!     *             Q_kappa_L2M(i_L2M+1),Q_s_L2M(i_L2M+1),
!     *             (-R_L2M(i_L2M+1)+6371D3)*1d-3
             i_L2M=i_L2M+1
           ENDIF

         ELSEIF (prop(mat(d_Ix,d_Iy,i))%litho.EQ.0) THEN
           i_L2M=i_L2M+1
!      write(*,*)'i,ii,i_L2M 2 loop dentro ',i,i_L2M,z
           R_L2M(i_L2M)=6371D3+z
           rho_L2M(i_L2M)=DENS(d_Ix, d_Iy, i)
           if(vel_p(d_Ix, d_Iy, i).eq.0) then
             vpv_L2M(i_L2M) = vpv_L2M(i_L2M-1)
             vph_L2M(i_L2M) = vph_L2M(i_L2M-1)
           else
             vpv_L2M(i_L2M)=vel_p(d_Ix, d_Iy, i)*1D3 !m/s
             vph_L2M(i_L2M)=vel_p(d_Ix, d_Iy, i)*1D3
           endif

	   if(vel_s(d_Ix, d_Iy, i).eq.0) then
             vsv_L2M(i_L2M) = vsv_L2M(i_L2M-1)
             vsh_L2M(i_L2M) = vsh_L2M(i_L2M-1)
           else
             vsv_L2M(i_L2M)=vel_s(d_Ix, d_Iy, i)*1D3*(1-anis_mnt(i)/3D0)
             vsh_L2M(i_L2M)=vel_s(d_Ix, d_Iy, i)*1D3*(1+2D0*anis_mnt(i)/3D0)
           endif

           Q_kappa_L2M(i_L2M)=qkappa_mineos(1)

 	   if(Q_s(d_Ix, d_Iy, i).eq.0) then
             Q_s_L2M(i_L2M) = Q_s_L2M(i_L2M-1)
           else
             Q_s_L2M(i_L2M)=Q_s(d_Ix, d_Iy, i)
           endif

!           WRITE(*,*)"c1",i, z, i_L2M, R_L2M(i_L2M),rho_L2M(i_L2M),
!     *             vpv_L2M(i_L2M),vsv_L2M(i_L2M),
!     *             Q_kappa_L2M(i_L2M),Q_s_L2M(i_L2M),
!     *             (-R_L2M(i_L2M)+6371D3)*1d-3

           IF((prop(mat(d_Ix,d_Iy,i))%rho.NE.prop(mat(d_Ix,d_Iy,i-1))%rho)
     *         ) THEN
!            write(*,*)'steep!',vel_s(i),vsv_L2M(i_L2M),vsv_L2M(i_L2M-1)
!            R_L2M(i_L2M)=R_L2M(i_L2M-1)

             if((6371D3+z+d_z).gt.(6371D3+E(d_Ix, d_Iy))) then
               R_L2M(i_L2M+1)=6371D3+E(d_Ix, d_Iy)
             else
               R_L2M(i_L2M+1)=6371D3+z+d_z
             endif

             rho_L2M(i_L2M+1)=rho_L2M(i_L2M)
             vpv_L2M(i_L2M+1)=vpv_L2M(i_L2M)
             vsv_L2M(i_L2M+1)=vsv_L2M(i_L2M)
             vph_L2M(i_L2M+1)=vph_L2M(i_L2M)
             vsh_L2M(i_L2M+1)=vsh_L2M(i_L2M)
             Q_kappa_L2M(i_L2M+1)=Q_kappa_L2M(i_L2M)
             Q_s_L2M(i_L2M+1)=Q_s_L2M(i_L2M)
!             WRITE(*,*)"c2",i,z, i_L2M+1, R_L2M(i_L2M+1),rho_L2M(i_L2M+1),
!     *             vpv_L2M(i_L2M+1),vsv_L2M(i_L2M+1),
!     *             Q_kappa_L2M(i_L2M+1),Q_s_L2M(i_L2M+1),
!     *             (-R_L2M(i_L2M+1)+6371D3)*1d-3
             i_L2M=i_L2M+1
           ENDIF
         ENDIF
       ENDDO

       N_L2M = i_L2M

!      write(*,*)"z_lay(0)", z_lay(d_Ix, d_Iy, 0)

       IF(E(d_Ix, d_Iy)<0) THEN
         write(*,*) "water layer"
         N_L2M = i_L2M+1
         R_L2M(N_L2M)=R_L2M(N_L2M-1)
         R_L2M(N_L2M+1)=6371D3
         rho_L2M(N_L2M:N_L2M+1)=1030D0
         vpv_L2M(N_L2M:N_L2M+1)=1450D0
         vsv_L2M(N_L2M:N_L2M+1)=0D0
         vph_L2M(N_L2M:N_L2M+1)=1450D0
         vsh_L2M(N_L2M:N_L2M+1)=0D3
         Q_kappa_L2M(N_L2M:N_L2M+1)=qkappa_mineos(1)
         Q_s_L2M(N_L2M:N_L2M+1)=300D0
         N_L2M = N_L2M+1
       endif
!        N_L2M = i_L2M+1
!        R_L2M(N_L2M)=6371D3+z-d_z*2

!        R_L2M(N_L2M+1)=6371D3+E(d_Ix, d_Iy)
!        rho_L2M(N_L2M:N_L2M+1)=1030D0
!        vpv_L2M(N_L2M:N_L2M+1)=1450D0
!        vsv_L2M(N_L2M:N_L2M+1)=0D0
!        vph_L2M(N_L2M:N_L2M+1)=1450D0
!        vsh_L2M(N_L2M:N_L2M+1)=0D3
!        Q_kappa_L2M(N_L2M:N_L2M+1)=qkappa_mineos(1)
!        Q_s_L2M(N_L2M:N_L2M+1)=300D0


!        WRITE(*,*)"w",N_L2M, R_L2M(N_L2M),rho_L2M(N_L2M),
!     *             vpv_L2M(N_L2M),vsv_L2M(N_L2M),
!     *             Q_kappa_L2M(N_L2M),Q_s_L2M(N_L2M),
!     *             (-R_L2M(N_L2M)+6371D3)*1d-3

!        WRITE(*,*)"w",N_L2M+1, R_L2M(N_L2M+1),rho_L2M(N_L2M+1),
!     *             vpv_L2M(N_L2M+1),vsv_L2M(N_L2M+1),
!     *             Q_kappa_L2M(N_L2M+1),Q_s_L2M(N_L2M+1),
!     *             (-R_L2M(N_L2M+1)+6371D3)*1d-3
!        N_L2M = N_L2M+1
!       ENDIF

!      pause
!Write input model for MINEOS
       OPEN(1000,file='ME01_litmod_ray_newver',STATUS='unknown')
       DO i=1,N_L2M
       WRITE(1000,97)R_L2M(i),rho_L2M(i),vpv_L2M(i),vsv_L2M(i),
     *             Q_kappa_L2M(i),Q_s_L2M(i),(-R_L2M(i)+6371D3)*1d-3
       ENDDO
       close(1000)

       OPEN(1000,file='ME01_litmod_love_newver',STATUS='unknown')
       DO i=1,N_L2M
        WRITE(1000,97)R_L2M(i),rho_L2M(i),vph_L2M(i),vsh_L2M(i),
     *             Q_kappa_L2M(i),Q_s_L2M(i),(-R_L2M(i)+6371D3)*1d-3
       ENDDO
       close(1000)
!       pause
 97    FORMAT(F9.0,1X,F9.2,1X,F9.2,1X,F9.2,1X,F9.1,1X,F9.1,F8.1)

       END SUBROUTINE
