!!!!!!!!!!!!!!!!SUB_read_DAT.for!!!!!!!!!!!!!!!!!!!!!!!!
!version for LITMOD_1D_MINEOS

      SUBROUTINE read_DAT

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

!     MODULE M_grid_parameters_1D
!     INTEGER  N_z,N_z_lit,N_z_c,N_top,I_z,isg_crt,N_moho
!     REAL(8) d_z
!     INTEGER N_lay,temp_calc,N_lay_cr
!     LOGICAL mltfrst,frt_crt_mlt
!     END MODULE


!     module M_grid_parameters
!     IMPLICIT INTEGER (I-N)
!     IMPLICIT real(8) (A-H,O-Z)
!     integer(4) N_x,N_y,N_z
!     integer(4) N_B
!     real(8) d_x,d_y,d_z,x,y,z
!     real(8) c_za,c_x,c_y,c_z
!     real(8), allocatable:: c_zs(:,:),d_z_topo(:,:)
!     real (8) L_x,L_y,L_z
!     real (8) E_max,E_min,z_lmax,z_lmin
!     integer(4) N_NOD,N_NODic,NOD_TOT,N_diag
!     integer(2) topo,temp_calc
!     end module M_grid_parameters

      USE M_column_1D,ONLY: anis_mnt,MDL,MDL_dm

      REAL(8) an1,an2,slp_a,z1,z2
      INTEGER I_z1,I_z2,I_zz,N_ref
      LOGICAL STAT
!Read radial anisotropy
!       write(*,*) 'I_anis',I_anis
       IF(I_anis==1) THEN 					!M_surf_waves
        WRITE(*,*)'ANISOTROPIC CASE'
        INQUIRE(FILE=name_anis,EXIST=STAT)
        IF (STAT) THEN
         OPEN(19,file=name_anis,STATUS='OLD')
         OPEN(18,file=indexstring//'_anis.z',STATUS='UNKNOWN')
        ELSE
         WRITE(*,*) 'FILE '//trim(name_anis)//' NOT FOUND!'
         WRITE(*,*) 'seting ianis=0 for isotropic case!'
         anis_mnt(:)=0D0 					!M_column_1D
         ianis=0 						!??????
!         STOP
        ENDIF
!Read and redefine the anisotropic column to have N_z nodes
       N_ref=0
       DO
        READ (19,*,IOSTAT=ios)
        IF (IOS<0) EXIT
        N_ref=N_ref+1
       ENDDO
       REWIND(19)
!      write(*,*)'N_ref',N_ref
      ALLOCATE (MDL(N_ref,2),MDL_dm(N_ref,2))  			!M_column_1D
      MDL=0 							!M_column_1D
      MDL_dm=0 							!M_column_1D
      DO i=1,N_ref
       READ (19,*)MDL(i,:)					!M_column_1D
      ENDDO
!Reorder, if necessary
      IF (MDL(1,1)>MDL(N_ref,1))THEN				!M_column_1D
       MDL_dm(:,:)=MDL(:,:)					!M_column_1D
       DO ik=1,N_ref
        MDL(ik,:)=MDL_dm(N_ref-ik+1,:)				!M_column_1D
       ENDDO
      ENDIF
!Redefine the reference column to have N_z nodes
      DO i=1,N_ref-1
       an1=MDL(i,2)						!M_column_1D
       an2=MDL(i+1,2)						!M_column_1D
       I_z1=NINT((MDL(i,1)*1D3+E_max)/d_z)			!M_column_1D
       I_z2=NINT((MDL(i+1,1)*1D3+E_max)/d_z)			!M_column_1D
       IF(I_z1<0)I_z1=1
       IF(I_z2>N_z)I_z2=N_z
       z1=((I_z1-1)*d_z)
       z2=((I_z2-1)*d_z)
       IF(z2>z1) THEN
        slp_a=(an2-an1)/(z2-z1)
        DO I_zz=I_z1+1,I_z2
         z=(I_zz-1)*d_z-E_max					!M_column_1D
         anis_mnt(I_zz)=(an1+slp_a*(z-z1))/1D2			!M_column_1D
         write(18,*)anis_mnt(I_zz)*1D2 ,-z*1D-3			!M_column_1D
        ENDDO
       ENDIF
      ENDDO

       ELSE
        WRITE(*,*)'ISOTROPIC CASE'
        anis_mnt(:)=0D0						!M_column_1D
       ENDIF
       !      z_lay=z_lay*1D3





!Surface wave data
      OPEN(3,file=name_ray,STATUS='OLD',IOSTAT=II)		!M_surf_waves
      IF(II==0) THEN
!       READ(3,*)name_ray
       READ(3,*)N_ray						!M_surf_waves
       ALLOCATE(v_ray_dat(N_ray,3),v_ray_calc(N_ray))		!M_surf_waves
       DO i=1,N_ray						!M_surf_waves
        READ(3,*)v_ray_dat(i,1),v_ray_dat(i,2),v_ray_dat(i,3)	!M_surf_waves
        v_ray_dat(i,2)=v_ray_dat(i,2)*1D3			!M_surf_waves
        v_ray_dat(i,3)=v_ray_dat(i,3)*1D3			!M_surf_waves
       ENDDO
       CLOSE(3)
       SW_RAY_on=.TRUE.						!M_surf_waves
       WRITE(*,*)'RAYLEIGH SURFACE WAVES INVERSION SWITCHED ON'
       write(*,*)'INPUT FILE FOR RAYLEIGH PHASE VEL: ',trim(name_ray)	!M_surf_waves
      ELSE
       WRITE(*,*)'RAYLEIGH SURFACE WAVES INVERSION SWITCHED OFF'
       SW_RAY_on=.FALSE.					!M_surf_waves
      ENDIF
      OPEN(3,file=name_love,STATUS='OLD',IOSTAT=II)		!M_surf_waves
      IF(II==0) THEN
       READ(3,*)N_love						!M_surf_waves
       ALLOCATE(v_love_dat(N_love,3),v_love_calc(N_love))	!M_surf_waves

       DO i=1,N_love						!M_surf_waves
        READ(3,*)v_love_dat(i,1),v_love_dat(i,2),v_love_dat(i,3)	!M_surf_waves
         v_love_dat(i,2)=v_love_dat(i,2)*1D3			!M_surf_waves
         v_love_dat(i,3)=v_love_dat(i,3)*1D3			!M_surf_waves
       ENDDO
       CLOSE(3)
       SW_LOVE_on=.TRUE.					!M_surf_waves
       WRITE(*,*)'LOVE SURFACE WAVES INVERSION SWITCHED ON'
       write(*,*)'INPUT FILE FOR love PHASE VEL: ',trim(name_LOVE)	!M_surf_waves
      ELSE
       WRITE(*,*)'LOVE SURFACE WAVES INVERSION SWITCHED OFF'
       SW_LOVE_on=.FALSE. 					!M_surf_waves
      ENDIF
!Read input reference model for surface waves
      IF(I_anis==0) THEN					!M_surf_waves
!Isotropic case
       OPEN(3,file='ME01',STATUS='OLD',IOSTAT=II)
!Read header
       DO i=1,103
        j=i-3
        IF (i<4) THEN
         READ(3,*)
        ELSE
        READ(3,*) r_mineos(j),rho_mineos(j),vpv_mineos(j),vsv_mineos(j),	!M_surf_waves
     *            qkappa_mineos(j),qshear_mineos(j)            			!M_surf_waves
          vph_mineos(j)=vpv_mineos(j)						!M_surf_waves
          vsh_mineos(j)=vsv_mineos(j)						!M_surf_waves
        ENDIF
       ENDDO
      ELSE
!Anisotropic case
       OPEN(3,file='ME01b',STATUS='OLD',IOSTAT=II)
!Read header
       DO i=1,103
         j=i-3
        IF (i<4) THEN
         READ(3,*)
        ELSE
        READ(3,*) r_mineos(j),rho_mineos(j),vpv_mineos(j),vsv_mineos(j),	!M_surf_waves
     *             qkappa_mineos(j),qshear_mineos(j),				!M_surf_waves
     *             vph_mineos(j),vsh_mineos(j)					!M_surf_waves
        ENDIF
       ENDDO

      ENDIF

      CLOSE(3)

      WRITE(*,*)'read_DAT FINISHED!'
      END SUBROUTINE
