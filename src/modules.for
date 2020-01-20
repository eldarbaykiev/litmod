

!##################################################
      ! Module with variables peculiar of the itherc LitMod version.
      module m_itherc
      real(8) perden, pervp, pervs
      real(8), allocatable :: inpdenarr(:,:,:)
      end module
!##################################################


!##################################################
      module m_pervar_III
      include 'perplex_parameters.h'
      integer i,ier,idead
      logical nodata, bulk
      character amount*6, yes*1
      integer itri(4),jtri(4),ijpt
      double precision wt(3), num
      integer iwt
      common/ cst209 /iwt
      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
      double precision atwt
      common/ cst45 /atwt(k0)
      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)
      character*5 cname
      common/ csta4 /cname(k5)
      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
      integer io3,io4,io9
      common / cst41 /io3,io4,io9
      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn
      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
      integer iam
      common/ cst4 /iam
		end module
!##################################################


!##################################################
      module m_pervar_II
      include 'perplex_parameters.h'
      integer i,ier,idead
      logical nodata, bulk
      character amount*6, yes*1
      integer itri(4),jtri(4),ijpt
      double precision wt(3), num
      integer iwt
      common/ cst209 /iwt
      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
      double precision atwt
      common/ cst45 /atwt(k0)
      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)
      character*5 cname
      common/ csta4 /cname(k5)
      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
      integer io3,io4,io9
      common / cst41 /io3,io4,io9
      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn
      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
      integer iam
      common/ cst4 /iam
      end module
!##################################################


!##################################################
      module m_pervar
      include 'perplex_parameters.h'
      integer i,ier,idead
      logical nodata, bulk
      character amount*6, yes*1
      integer itri(4),jtri(4),ijpt
      double precision wt(3), num
      integer iwt
      common/ cst209 /iwt
      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
      double precision atwt
      common/ cst45 /atwt(k0)
      double precision v,tr,pr,r,ps
      common/ cst5  /v(l2),tr,pr,r,ps
      integer ipot,jv,iv
      common / cst24 /ipot,jv(l2),iv(l2)
      character*8 vname,xname
      common/ csta2  /xname(k5),vname(l2)
      character*5 cname
      common/ csta4 /cname(k5)
      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp
      integer io3,io4,io9
      common / cst41 /io3,io4,io9
      logical gflu,aflu,fluid,shear,lflu,volume,rxn
      common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn
      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
      integer iam
      common/ cst4 /iam
      end module
!##################################################





!##################################################
      module M_grid_parameters

      IMPLICIT INTEGER (I-N)
      IMPLICIT real(8) (A-H,O-Z)
      ! Grid steps :d_x,d_y,d_z1,d_z2
      ! Number of nodes in the x, y and z axis: N_X,N_Y,N_Z
      integer(4) N_x,N_y,N_z
      integer(4) N_B
      real(8) d_x,d_y,d_z,x,y,z
      real(8) c_za,c_x,c_y,c_z
      real(8), allocatable:: c_zs(:,:),d_z_topo(:,:)
      real (8) L_x,L_y,L_z
      real (8) E_max,E_min,z_lmax,z_lmin
      integer(4) N_NOD,N_NODic,NOD_TOT,N_diag
      integer(2) topo,temp_calc

      end module M_grid_parameters
!##################################################




!##################################################
      module M_material


      TYPE MATERIAL
      ! K -->Thermal conductivity
      ! A -->Heat production rate
      ! rho-->Density at T=Ts,P=0
      ! alpha-->Volumetric expansion coefficient
      ! beta-->Pressure coefficient
      ! name_mat-->Short material description
      ! gmma_T-->Thermal Grüneisen parameter
      ! K_0p-->Pressure derivative of isothermal bulk modulus
      ! K_T-->isothermal bulk modulus
      ! litho-->Index describing the type of mantle material
      ! litho=0-->crustal material
      ! litho>1-->Lithospheric mantle
      ! litho=99-->Sublithospheric mantle
      real(8) K,A,rho,alpha,beta,gmma_T,K_0p,K_T
      INTEGER litho
      INTEGER crst
      character(LEN=50)name_mat

      integer inp             ! Define whether the properties are input or not.
      character(len=50)inpfil ! Define the name of the file to read the properties from.
      real(8), allocatable::inpden(:,:,:)

      real(8) Vp,Vs

      end TYPE MATERIAL


      TYPE (MATERIAL) prop (-1:100)
      integer(2), allocatable:: mat(:,:,:),DUM(:,:,:),
     *                          DUM_sublit(:,:,:)
      INTEGER N_mant_max


      integer checkinp


      end module M_material
!##################################################




!##################################################
      module M_temperatures
      real(8), allocatable:: T(:,:,:),T_old(:,:,:),k_T_P(:,:,:)
      real(8), allocatable:: T_vec(:),T_vec_dm(:)
      real(8) Ts,Ta,HF,z_lit,z_lit_buf
      INTEGER sw_clm
      end module M_temperatures
!##################################################



!##################################################
      module M_layers
      INTEGER N_lay,N_lay_cr
      CHARACTER(LEN=*), PARAMETER::string='./layers_xy/layer',
     *                              num="0123456789"
      CHARACTER(50) :: fil
      real(8), allocatable:: layers(:,:,:)
      real(8), allocatable:: E(:,:)
      end module M_layers
!##################################################



!##################################################
      module M_L_EQ_SYS
      real(8), allocatable:: B(:,:)
      real(8), allocatable:: A_RHS(:)
      real(8) UD_x,UD_y,UD_z,dgn_x,dgn_y,dgn_z
      REAL(8) LD_x,LD_y,LD_z
      integer(4) ITMAX,IT
      real(8) R1,RLAST,error
      integer(4),allocatable:: IDCOL(:,:),IDUM(:,:),IROW_node(:)
      integer(4) IROW_comp,NEQ
      end module M_L_EQ_SYS
!##################################################



!##################################################
      module M_columns
      real(8), allocatable:: DENS(:,:,:),BGA(:,:),FA(:,:),GEOID(:,:),
     * Uxx(:,:),Uyy(:,:),Uzz(:,:),Uzx(:,:),Uzy(:,:),Uxy(:,:)
      real(8), allocatable:: BGA_cr(:,:),FA_cr(:,:),GEOID_cr(:,:),
     * Uxx_cr(:,:),Uyy_cr(:,:),Uzz_cr(:,:),Uzx_cr(:,:),Uzy_cr(:,:),
     * Uxy_cr(:,:), Q_s(:,:,:)
      real(8), allocatable:: GEOID_OBS(:,:),FA_OBS(:,:),E_calc(:,:)
      real(8), allocatable:: BGA_OBS(:,:),P(:,:,:),Uxx_OBS(:,:)
     *,Uyy_OBS(:,:),Uzz_OBS(:,:),Uzx_OBS(:,:),Uzy_OBS(:,:),Uxy_OBS(:,:)
      real(8), allocatable:: rho_lay(:,:,:),z_lay(:,:,:)
      real(8), allocatable:: vel_p(:,:,:),vel_s(:,:,:)
      INTEGER(4),ALLOCATABLE::PH_TR(:,:,:,:),N_ph(:,:,:),
     *                        Gt(:),Sp(:),an(:)
      REAL(8),ALLOCATABLE::z_ph(:,:,:,:),rho_ph(:,:,:,:),P_CL(:,:)
      REAL(8),ALLOCATABLE:: D_OBS(:),D_calc(:)
      integer(2) mat_z,mat_z_plus,mat_z_minus
      integer(4) i_st,I_z
      real(8) rho,rho_L,AV,std_dev,Z_usec
      real(8) z_bot,rho_up,rho_down,z_up,z_down,rho_a,rho_red,rho_B
      real(8) w1,w2,rho_w,L,z_1,rho_1,xp,yp,GR,GEO,BG,BGA_av,FA_av,
     *        U_sec_av
      real(8) rho_ref,ra,RHO_AV,DELTA
      real(8) E_ref,z_cref,z_lref,rho_cref,rho_mref,rho_sublref
      real(8),dimension(2):: XB,YB,ZB,ZBG,ZB_g,ZBG_g
      real(8),parameter:: g_s=9.81D0,g_400=9.97D0
      REAL(8) g,d_gz,g_av
      REAL(8),PARAMETER:: E_MOR=2.6D3,dzc_MOR=6.5D3
!!      REAL(8),PARAMETER:: rho_m_MOR=3420.D0,rho_c_MOR=2930D0,
      REAL(8) rho_m_MOR
      INTEGER E_cali
      REAL(8),PARAMETER:: rho_c_MOR=2986D0,rho_bot=3611.3D0
      REAL(8) rho_MOR,pi_par
      REAL(8),PARAMETER::G_gr=6.67D-11,pi=3.14159,
     *                   G_ref_abs=11912.4
      end module
!##################################################




!##################################################
      module M_lsq_plane
      real(8) det_D
      real(8), allocatable:: Dt(:,:),GEOID_vec(:),GEOID_OBS_vec(:),
     *                       Z_pl(:),Z_pl_cr(:),GEOID_vec_cr(:)
      real(8),dimension(3,3):: DtD_1
      real(8) S_x,S_y,S_x2,S_y2,S_xy
      real(8) ddmy,ddmy_o,DMY(3),DMY_O(3),CF_pl(3),CF_pl_O(3)
      real(8) ddmy_cr,DMY_cr(3),CF_pl_cr(3)

      end module
!##################################################




!##################################################
      module M_phases
      CHARACTER(14),ALLOCATABLE:: phases(:,:,:,:)
      REAL(8),ALLOCATABLE::vol_frac(:,:,:,:),wt_MgO(:,:,:,:),
     *                     wt_FeO(:,:,:,:)
      INTEGER,ALLOCATABLE:: I_phases(:,:,:,:)
      LOGICAL wads_tr
      INTEGER,PARAMETER:: N_tot_ph=15
      REAL(8),PARAMETER:: mm_FeO=71.844D0,mm_MgO=40.3044
      CHARACTER(14) str_ph
      INTEGER N_phases,N_pt_LK
      end module
!##################################################




!##################################################
      module M_mant
!LK
      CHARACTER(120) mant_file(99),mant_file_sys(99),mant_file_ph(99)
      INTEGER I_mnt,I_crst,I_mant_file(99),i_near,j_near,i_LK
      INTEGER,ALLOCATABLE:: LK_index(:),N_LK(:),N_LK_T(:),N_LK_P(:),
     *                      litho2LK(:)
      REAL(8),ALLOCATABLE::LK_vec(:,:,:,:),dm_vec(:)
      REAL(8),ALLOCATABLE:: d_P_LK(:),d_T_LK(:),T_LK_0(:),P_LK_0(:)
      end module
!##################################################




!##################################################
      module M_sublit
      INTEGER N_sublit,N_sublit_t
      INTEGER,ALLOCATABLE::I_sublit(:),IND_SUB(:,:),litho_sub(:)
      INTEGER,ALLOCATABLE::I_sub_cont(:,:),I_sub_dis(:,:,:)
      REAL(8),ALLOCATABLE::D_T_sub(:)
      LOGICAL,ALLOCATABLE::RON(:)
      LOGICAL sub_onoff
      CHARACTER(1)chdm_s
      end module
!##################################################




!##################################################
       MODULE M_conductivity
       INTEGER N_sg_cr,pt_dec,N_w,N_lay_sg
       REAL(8),PARAMETER::k_bolt=8.61734315D-5,R_uni=8.314472D-3
!Parameters for garnet Dai-Karato_09 PEPI
!Dry
       REAL(8),PARAMETER::A_0_Dai_d=1036,B_Dai_d=0.044,E_dai_d=128,
     *                    V_Dai_d=2.5
!Wet
       REAL(8),PARAMETER::A_0_Dai_w=195,E_dai_w=70,V_Dai_w=-.57
!Parameters for Vacher and Verhoeven 2007
!Olivine
       REAL(8),PARAMETER::sg_ref_ol_VV=2.69,al_ol_VV=2.42,
     * U_ref_ol_VV=1.62,bt_ol_VV=-.48
!Opx
       REAL(8),PARAMETER::sg_ref_opx_VV=3.66,al_opx_VV=-1.28,
     * U_ref_opx_VV=1.79,bt_opx_VV=-1.41
!Cpx
      REAL(8),PARAMETER::sg_ref_cpx_VV=3.19,al_cpx_VV=-1.28,
     * U_ref_cpx_VV=1.86,bt_cpx_VV=-1.41
       REAL(8) X_Fe_ol,X_Fe_gt,X_Fe_cpx,X_Fe_opx
       REAL(8),ALLOCATABLE::C_w_ol(:),C_w_opx(:),C_w_cpx(:),
     *                      C_w_gt(:),C_w_mod(:,:,:),z_cw(:)
       REAL(8),ALLOCATABLE:: melt(:),melt_mod(:,:),Na_melt(:),
     *                       Si_melt(:),Cw_melt_melt(:)
       REAL(8) sig_min,sig_max,sum_min,sum_max,sum_s,sum_p
       REAL(8) sigma(5),res(4),sigma_mlt


       REAL(8) rho_HS,dz_elas,rho_elas,C_arch,n_arch
       INTEGER mlt_mod,melt_mix,ol_mod,gt_mod,opx_mod,mix_rule
       INTEGER,ALLOCATABLE:: C_w_onoff(:),melt_onoff(:)

       end MODULE
!##################################################




!##################################################
      MODULE M_fugacity_ol_solu
      REAL(8) rho_0,rho_mol,c_EoS(10),c_T(10,6),T_vec_C(6),cte
      REAL(8) F,F0,Ffwd,Fbwd,del_rho,Ares_adim,fH2O,C_OH
      REAL(8) C_OH_cpx,C_OH_cpx_0,C_OH_opx_Al,C_OH_opx,C_OH_opx_0,
     *        C_OH_gt
      REAL(8),PARAMETER::Azhao=5.62D3,E_zhao=50,dV_zhao=10,alpha_zhao=97
      REAL(8),PARAMETER::Amie=13.54D1,E_mie=-4.563,dV_mie=12.1
      REAL(8),PARAMETER::Agav=18.5D1,E_gav=-11.117,dV_gav=14.62
      REAL(8),PARAMETER::Amie_al=4.2D0,E_mie_al=-79.685,dV_mie_al=11.3
!      REAL(8),PARAMETER::Alu=67.9D0,E_lu=-14D0,dV_lu=5.71
      REAL(8),PARAMETER::Alu=18.09D0,E_lu=-14D0,dV_lu=5.71
      end MODULE
!##################################################




!##################################################
      MODULE M_SIGMELTS
      REAL(8),PARAMETER::rNa=.95D-10,rO=1.4D-10,Navog=6.02214179D23
      REAL(8) gma,E_b,G_sh,E_s,E_pom,V_melt,H_pom,s_0_pom
      REAL(8) intc,C_Si,C_T,C_p,C_LH2O,C_H2O,C_Na,C_Ha,Cw_eps
      end MODULE
!##################################################




!##################################################
      MODULE M_TARANTOLA
      INTEGER N_GTOT,N_DAT,N_PAR
      REAL(8),ALLOCATABLE:: G_vec(:)
      INTEGER,ALLOCATABLE:: IG_COL(:),IG_ROW(:)
      REAL(8),ALLOCATABLE:: D_OBS(:),m_PAR(:),C_OBS(:),C_PAR(:)
      REAL(8),ALLOCATABLE:: Qdmm(:),V_dummy(:)
      end MODULE
!##################################################





!MINEOS MODULES
!##################################################
      MODULE M_surf_waves
       INTEGER(4) N_ray,N_love,I_anis
       CHARACTER(300)name_ray,name_love,name_anis
       CHARACTER(3)indexstring
       REAL(8),ALLOCATABLE::v_ray_dat(:,:),v_ray_calc(:),
     *                      v_love_dat(:,:),v_love_calc(:),srf_wv(:,:)
       REAL(8),ALLOCATABLE::v_ray_calc_best(:),v_love_calc_best(:)
       LOGICAL MT_on
       LOGICAL SW_RAY_on,SW_LOVE_on
       REAL(8) r_mineos(100),rho_mineos(100),vpv_mineos(100),
     *         vsv_mineos(100),qkappa_mineos(100),qshear_mineos(100),
     *         vph_mineos(100),vsh_mineos(100)
       REAL(8) R_L2M(300),rho_L2M(300),vpv_L2M(300),vsv_L2M(300),
     *         Q_kappa_L2M(300),Q_s_L2M(300),vph_L2M(300),vsh_L2M(300)
       REAL(8) z_anis(11),anis(11)
      INTEGER N_L2M
      REAL(8)vs_660,vs_410,an_660,an_410,anis_tran
      END MODULE
!##################################################

!##################################################
!     MODULE M_grid_parameters_1D
!     INTEGER  N_z,N_z_lit,N_z_c,N_top,I_z,isg_crt,N_moho
!     REAL(8) d_z
!     INTEGER N_lay,temp_calc,N_lay_cr
!     LOGICAL mltfrst,frt_crt_mlt
!     END MODULE
!##################################################


!##################################################
      MODULE M_column_1D
      REAL(8),ALLOCATABLE:: anis_mnt(:),MDL(:,:),MDL_dm(:,:)
      END MODULE
!##################################################


      module M_d
      integer d_coordnum, d_coordind
      real(8) d_coordx, d_coordy
      real(8) d_dist, d_distmax
      integer d_Ix, d_Iy
      REAL(8) d_size, volexp
      end module
