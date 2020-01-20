!Subroutine:  TEMPERATURE_3D
!Sistema de ecs:   B*T_vec=A_RHS


      subroutine TEMPERATURE_3D 
      use M_temperatures
      use M_grid_parameters
      use M_L_EQ_SYS
      use M_material
      
      allocate (B(NOD_TOT,4)) 
      allocate (A_RHS(NOD_TOT))
      allocate (T_vec(NOD_TOT),T_vec_dm(NOD_TOT))
!!      open(11,file='OUT_xyz/T.xyz',status='UNKNOWN')
      open(33,file='B_solve.inf',status='UNKNOWN')
!      open(13,file='OUT_xyz/HF.xyz',status='UNKNOWN')
!      open(50,file='OUT_xyz/T_400.xyz',status='UNKNOWN')  
      write(*,*)'THERMAL PROBLEM -->  B*T_vec=A_RHS'
      c_zs=Ts/(2*d_z*d_z)
      c_za=Ta/(2*d_z*d_z)
      c_x=1/(2*d_x*d_x)
      c_y=1/(2*d_y*d_y) 
      c_z=1/(2*d_z*d_z)
      N_diag=1
      DO I_z=2,N_z-1
       DO I_y=2,N_y-1
        DO I_x=2,N_x-1
         IF (DUM(I_x,I_y,I_z)<1) THEN
          B(N_diag,1:4)=0D0
          A_RHS(N_diag)=0D0
         ELSE
         call UpperD_x(k_T_P(I_x,I_y,I_z),
     *                 k_T_P(I_x+1,I_y,I_z),UD_x)
         call UpperD_y(k_T_P(I_x,I_y,I_z),
     *                 k_T_P(I_x,I_y+1,I_z),UD_y)
         call UpperD_z(k_T_P(I_x,I_y,I_z),
     *                 k_T_P(I_x,I_y,I_z+1),UD_z)
          B(N_diag,2)=UD_x
          B(N_diag,3)=UD_y
          B(N_diag,4)=UD_z
          A_RHS(N_diag)=-prop(mat(I_x,I_y,I_z))%A          

!Boundary conditions: Fixed temperature
!Eje X*********************************************
          if (DUM(I_x-1,I_y,I_z)<-99) then
          A_RHS(N_diag)=A_RHS(N_diag)-c_x*
     *                  Ts*(k_T_P(I_x-1,I_y,I_z)+
     *                  k_T_P(I_x,I_y,I_z))
          endif
          if (DUM(I_x+1,I_y,I_z)<-99) then
          A_RHS(N_diag)=A_RHS(N_diag)-Ts*B(N_diag,2)
          B(N_diag,2)=0D0
          endif 
          if ((DUM(I_x-1,I_y,I_z)<1).AND.(DUM(I_x-1,I_y,I_z)>-100)) then
          A_RHS(N_diag)=A_RHS(N_diag)-c_x*
     *                  Ta*(k_T_P(I_x-1,I_y,I_z)+
     *                  k_T_P(I_x,I_y,I_z))
          endif
          if ((DUM(I_x+1,I_y,I_z)<1).AND.(DUM(I_x+1,I_y,I_z)>-100)) then
          A_RHS(N_diag)=A_RHS(N_diag)-Ta*B(N_diag,2)
          B(N_diag,2)=0D0
          endif
!Eje Y*********************************************         
          if (DUM(I_x,I_y-1,I_z)<-99) then
          A_RHS(N_diag)=A_RHS(N_diag)-c_y*
     *                  Ts*(k_T_P(I_x,I_y-1,I_z)+
     *                  k_T_P(I_x,I_y,I_z))
          endif
          if (DUM(I_x,I_y+1,I_z)<-99) then
          A_RHS(N_diag)=A_RHS(N_diag)-Ts*B(N_diag,3)
          B(N_diag,3)=0D0
          endif
          if ((DUM(I_x,I_y-1,I_z)<1).AND.(DUM(I_x,I_y-1,I_z)>-100)) then
          A_RHS(N_diag)=A_RHS(N_diag)-c_y*
     *                  Ta*(k_T_P(I_x,I_y-1,I_z)+
     *                  k_T_P(I_x,I_y,I_z))
          endif
          if ((DUM(I_x,I_y+1,I_z)<1).AND.(DUM(I_x,I_y+1,I_z)>-100)) then
          A_RHS(N_diag)=A_RHS(N_diag)-Ta*B(N_diag,3)
          B(N_diag,3)=0D0
          endif 
!Eje Z*********************************************
          if (DUM(I_x,I_y,I_z-1)<-99) then                              
          A_RHS(N_diag)=A_RHS(N_diag)-c_zs(I_x,I_y)*
     *                  (k_T_P(I_x,I_y,I_z-1)+
     *                  k_T_P(I_x,I_y,I_z))                      
          elseif ((DUM(I_x,I_y,I_z+1)<1).AND.(DUM(I_x,I_y,I_z+1)>-100))
     *    then
          A_RHS(N_diag)=A_RHS(N_diag)-Ta*B(N_diag,4)
          B(N_diag,4)=0D0
          endif

          call Diag_z(k_T_P(I_x,I_y,I_z-1),
     *               k_T_P(I_x,I_y,I_z),
     *               k_T_P(I_x,I_y,I_z+1),dgn_z)
!Boundary conditions:lateral heat flow = 0
          if (I_x==2) then
          call UpperD_x(k_T_P(I_x,I_y,I_z),
     *                 k_T_P(I_x+1,I_y,I_z),UD_x)
          dgn_x=-UD_x
          elseif (I_x==N_x-1) then
          call n_Diag_x(k_T_P(I_x-1,I_y,I_z),
     *                 k_T_P(I_x,I_y,I_z),dgn_x )
          B(N_diag,2)=0D0
          else
          call Diag_x(k_T_P(I_x-1,I_y,I_z),
     *                k_T_P(I_x,I_y,I_z),
     *                k_T_P(I_x+1,I_y,I_z),dgn_x)
          endif

          if (I_y==2) then
          call UpperD_y(k_T_P(I_x,I_y,I_z),
     *                 k_T_P(I_x,I_y+1,I_z),UD_y)
          dgn_y=-UD_y
          elseif (I_y==N_y-1) then
          call n_Diag_y(k_T_P(I_x,I_y-1,I_z),
     *                 k_T_P(I_x,I_y,I_z),dgn_y )
          B(N_diag,3)=0D0
          else
          call Diag_y(k_T_P(I_x,I_y-1,I_z),
     *                k_T_P(I_x,I_y,I_z),
     *                k_T_P(I_x,I_y+1,I_z),dgn_y)
          endif
          B(N_diag,1)=dgn_x+dgn_y+dgn_z
         ENDIF
        N_diag=N_diag+1
        ENDDO
       ENDDO
      ENDDO 

!Rewrite matrix B compresed whithout 0's (constant temperature nodes)
      allocate (IDCOL(NOD_TOT,2))
      allocate (IDUM(NOD_TOT,2))
      allocate (IROW_node(NOD_TOT))
      if(1==2) then
      IDCOL(:,1)=N_x-2
      IDCOL(:,2)=(N_x-2)*(N_y-2)
      NEQ= NOD_TOT
      else
      IROW_comp=0
      IROW_node=0
!SOLO ES NECESARIO DEFINIR ESTAS COLUMNAS DE IDCOL, LA
!1º ES 0 Y LA 2º 1.IDCOL REPRESENTA DISTANCIA A LA DIAGONAL
! DE LOS ELEMENTOS NO NULOS DE B
      IDCOL(:,1)=N_x-2
      IDCOL(:,2)=(N_x-2)*(N_y-2)  
      DO N_diag=1,NOD_TOT
       IF (B(N_diag,1).ne.0) ICONT=N_diag
      ENDDO
      DO N_diag=1,ICONT
       IF (B(N_diag,1).NE.0) THEN
          IROW_comp=IROW_comp+1
          B(IROW_comp,:)=B(N_diag,:)
          A_RHS(IROW_comp)=A_RHS(N_diag)
          IDCOL(IROW_comp,:)=IDCOL(N_diag,:)  
          IDUM(IROW_comp,:)=IDCOL(N_diag,:)
          IROW_node(IROW_comp)=N_diag
       ELSEIF (IROW_comp>0) THEN
          do i=IROW_comp,1,-1
           if((IDCOL(i,1)+IROW_node(i)).LT.N_diag ) EXIT
           IDUM(i,1)=IDUM(i,1)-1
          enddo
          do i=IROW_comp,1,-1
           if((IDCOL(i,2)+IROW_node(i)).LT.N_diag ) EXIT
           IDUM(i,2)=IDUM(i,2)-1
          enddo
       ENDIF
      ENDDO
      IDCOL=IDUM
      NEQ=IROW_comp      
      deallocate (IDUM,IROW_node)
!!!      write(*,*)'Nodes on the compressed B -->',IROW_comp
      endif
!Resolución del sistema de ecuaciones para la Temperatura
!¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬¬ 
!      IDCOL(1)=0
!      IDCOL(2)=1
!      IDCOL(3)=N_x-2
!      IDCOL(4)=(N_x-2)*(N_y-2)
      call LSQR_FD(33)
      write(*,*)'THERMAL PROBLEM SOLVED AFTER ',IT,'ITERATIONS'     
!From vector to matrix
      N_B=1
      DO I_z=2,N_z-1
       DO I_y=2,N_y-1
        DO I_x=2,N_x-1
         IF (DUM(I_x,I_y,I_z)>0) THEN
          T(I_x,I_y,I_z)=T_vec(N_B)
          N_B=N_B+1
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      T(1,:,:)=T(2,:,:)
      T(N_x,:,:)=T(N_x-1,:,:)
      T(:,1,:)=T(:,2,:)
      T(:,N_y,:)=T(:,N_y-1,:)
      T(:,:,1)=Ts
      T(:,:,N_z)=1520D0

      DO I_y=1,N_y
      y=d_y*(I_y-1)
       DO I_x=1,N_x
       x=d_x*(I_x-1) 
        sw_clm=0
        DO I_z=1,N_z
         z=E_max-(I_z-1)*d_z
         I_zplus1=I_z+1 
         IF (I_z==N_z) I_zplus1=I_z
!"""""""""""""ADIABATIC GRADIENT""""""""""""""""""""""""
       IF (DUM_sublit(I_x,I_y,I_zplus1)==0) THEN
!Get the depth of the lithosphere and the buffer relaxation zone
        IF (sw_clm==0) THEN
         z_lit=z
         z_lit_buf=z-40D3
!For thin lithospheres (<60 km thick) set the depth of the buffer 
!to 100 km
         IF (z_lit>=-60D3) z_lit_buf=-100D3
         sw_clm=1          
        ENDIF

         IF(z>z_lit_buf) THEN
!Inside the buffer
          T(I_x,I_y,I_z)=Ta+(z-z_lit)*(1400D0-Ta)/(z_lit_buf-z_lit) 
         ELSE  
!Adiabatic mantle
           grd_ad=(1520D0-1400D0)/(400D3+z_lit_buf)
!Limit the adiabatic gradient (note this may chacnge the temperature
!at the bottom of the model)
           IF (grd_ad<.35D-3) grd_ad=.35D-3
           IF (grd_ad>.6D-3) grd_ad=.6D-3
          T(I_x,I_y,I_z)=1400D0-(z-z_lit_buf)*grd_ad
         ENDIF

        ENDIF 
!"""""""""""""ADIABATIC GRADIENT""""""""""""""""""""""""
!      write(11,'(3F13.3,F10.2)')x,y,z,T(I_x,I_y,I_z)
!!       IF (I_z==N_z) write(50,'(2F13.3,F10.2)')x,y,T(I_x,I_y,I_z) 
!!         IF(I_z>1) THEN
!!          if (T(I_x,I_y,I_z-1)==Ts.AND.T(I_x,I_y,I_z)>Ts) then
!        HF=k_T_P(I_x,I_y,I_z+1)*(T(I_x,I_y,I_z+1)-T(I_x,I_y,I_z))/
!     *      d_z+d_z*prop(mat(I_x,I_y,I_z+1))%A+
!     *      prop(mat(I_x,I_y,I_z))%A*d_z_topo(I_x,I_y)
!!        HF=k_T_P(I_x,I_y,I_z)*(T(I_x,I_y,I_z)-T(I_x,I_y,I_z-1))/
!!     *      d_z+d_z*(prop(mat(I_x,I_y,I_z))%A+
!!     *      prop(mat(I_x,I_y,I_z-1))%A*.5D0)
!        HF=HF+d_z*(prop(mat(I_x,I_y,I_z))%A+prop(mat(I_x,I_y,I_z+1))%A+
!     *     prop(mat(I_x,I_y,I_z+2))%A)  
!         if ((I_z<N_z).AND.(T(I_x,I_y,I_z+1)==Ta).AND.
!     *      (T(I_x,I_y,I_z)<Ta)) then
!       HF=k_T_P(I_x,I_y,I_z)*(T(I_x,I_y,I_z)-T(I_x,I_y,I_z-1))/d_z
!     *      +d_z*prop(mat(I_x,I_y,I_z))%A/2
!!         write(13,1)x,y,HF
!!1     FORMAT (2F13.3,' ',F7.5)  
!!         endif
!!        ENDIF 
        ENDDO
       ENDDO
      ENDDO 
!!      close(11) 
      CLOSE(33)
      close(13)
      CLOSE(50)
      deallocate (T_vec,T_vec_dm,A_RHS,B,IDCOL)
      end subroutine
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
!      SUBROUTINE LSQR_FD(NEQ,ANZ,BIND,TSOL,IDCOL,ITMAX,IT
!     *      ,ERROR,IO,R1, RLAST,NOD_TOT,IROW_node)
       SUBROUTINE LSQR_FD(IO)
C     ==================
C
C**********************************************************************
C
C Subroutine solves iteratively a system of linear equations for 3D
C finite difference calculations. It is supposed that the matrix is
C symmetric and has potentially non-zero valus only on 4 diagonals in
C the upper triangle (plus 3 in the lower one). In order to reduce
C storage needs, these values are stored in a matrix of order (NEQ,4).
C
C**********************************************************************
C
C NEQ:   Number of equations (=number of parameters to be calculated)
C B  :   Matrix containing all potentially non-zero values in the matrix
C        the data are stored by diagonals, only the upper triangle is
C        considered;
C A_RHS:  Vector of independent terms
C T_vec:  Final result
C IDCOL: Vector with four entries describing the distance of the
C        diagonals from the main diagonal. Usually these values are:
C        IDCOL(1)=0, i.e. the main diagonal
C        IDCOL(2)=1, (3)=number of nodes in the innermost loop
C             (4)=product of number of nodes in the two innermost loops
C ITMAX: Maximum number of iterations allowed
C IT:    Number of iterations effectively used
C ERROR: Maximum error allowed for the solution
C IO:    Logical number of output device
C R1, RLAST: Initial and final error in solution
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use M_L_EQ_SYS
      use M_temperatures
      IMPLICIT REAL*8 (A-H,O-Z)
C      integer row(N*4),column(N*4)
      ALLOCATABLE :: DUM1(:),DUM2(:)
      IF(.NOT. ALLOCATED(DUM1)) ALLOCATE (DUM1(NEQ),DUM2(NEQ))
  
!       do N_diag=1,NEQ
!      write(13,*)IDCOL(N_diag,1),IDCOL(N_diag,2),IROW_node(N_diag)
!      write(8,*)B(N_diag,1),B(N_diag,2),B(N_diag,3),B(N_diag,4)
!      write(10,*) A_RHS(N_diag)  
!      enddo
!      close(13)
!      close(10)
!      close(8)
!      stop 
      T_vec=0.D0
      DUM1=0.D0
      DERROR=ERROR*0.01
c
      CALL NORMLZ(NEQ,A_RHS,BETA)
c
      B1=BETA
c
C      call atupv(m,n,u,v,anz,row,column)
      CALL ATUPV(A_RHS,DUM1)
      CALL NORMLZ(NEQ,DUM1,ALFA)
c
      RHOBAR=ALFA
      PHIBAR=BETA
c
      DUM2=DUM1
c
!!      DO ITER=1,ITMAX
      ITER=1
      DO
       T_vec_dm=T_vec
       A=-ALFA
       A_RHS=A*A_RHS
C       do i=1,m
C        u(i)=a*u(i)
C       end do
c
C     common /blk1/ anz,row,column,nnz
C       call avpu(m,n,u,v,anz,row,column)
       CALL AVPU(A_RHS,DUM1)
       CALL NORMLZ(NEQ,A_RHS,BETA)
c
       BB=-BETA
c
       DUM1=BB*DUM1
C       do i=1,n
C        v(i)=b*v(i)
C       end do
c
C       call atupv(m,n,u,v,anz,row,column)
       CALL ATUPV(A_RHS,DUM1)
       CALL NORMLZ(NEQ,DUM1,ALFA)
c
       RHO=DSQRT(RHOBAR*RHOBAR+BETA*BETA)
       C=RHOBAR/RHO
       S=BETA/RHO
       TETA=S*ALFA
       RHOBAR=-C*ALFA
       PHI=C*PHIBAR
       PHIBAR=S*PHIBAR
       T1=PHI/RHO
       T2=-TETA/RHO
c
       DO I=1,NEQ
        T_vec(I)=T1*DUM2(I)+T_vec(I)
        DUM2(I)=T2*DUM2(I)+DUM1(I)
       ENDDO
c
       R=PHIBAR/B1
c    
       IF(ITER .EQ. 1) THEN
          WRITE(IO,'(I7,D15.6,2D15.6,F10.3)') ITER,PHIBAR,R,0.,0D0
          R1=R
       ELSE
          WRITE(IO,'(I7,D15.6,3D15.6)') ITER,PHIBAR,R,RPREV-R,
     * MAXVAL(ABS(T_vec-T_vec_dm) )
          RLAST=R
       ENDIF
      
      IF (ITER>1) THEN 
!       IF(ABS(R-RPREV) .LT. DERROR) THEN
        IF(ABS(R)<1D-6.AND.ABS(R-RPREV)<DERROR) THEN
          WRITE(IO,*) 'Number of iterations in LSQR: ',iter
          IT=ITER
          RETURN
       ENDIF
      ENDIF
      IF (ITER==150000) THEN
       WRITE(*,*)'CONVERGENCE IS REALLY SLOW... TRY TO REDUCE THE
     *NUMBER OF NODES'
       WRITE(*,*)'LitMod3D aborts...'
       STOP
      ENDIF 
c
       RPREV=R
      ITER=ITER+1
      ENDDO
c
      IT=ITMAX
c
      END
c
c     ----------------------------------------------------------------
c
      SUBROUTINE NORMLZ(N,X,S)
c
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N)
c
      S=0.D0
c
      DO I=1,N
       S=S+X(I)*X(I)
      ENDDO
c
      S=DSQRT(S)
      SS=1.D0/S
c
      DO I=1,N
       X(I)=X(I)*SS
      ENDDO
c
      END
c
c     ----------------------------------------------------------------
c
C      subroutine avpu(m,n,u,v,anz,row,column)
      SUBROUTINE AVPU(U,V)
      use M_L_EQ_SYS 
C      include 'ray.par'
c
C      parameter(nnzmax=2500000, ndmax=25000)
C      parameter(nsconst=3*nximax*nzimax,
C     +          nszero=7*nximax*nzimax,
C     +          nmmax=nximax*nzimax)
c
      IMPLICIT REAL*8 (A-H,O-Z)
C      integer row(nnzmax+nszero),column(nnzmax+nszero)
C      real anz(nnzmax+nszero),u(m),v(n),sum(ndmax+nsconst)
C      real*8 anz(nnzmax+nszero),u(ndmax+nsconst),v(n)
      DIMENSION U(NEQ),V(NEQ),IA(4)
C      COMMON /blk1/ NNZERO,NNZMAX,NDMAX,NSZERO,NSCONST,NMMAX
      ALLOCATABLE :: SUM(:)
      data IA /0,1,0,0/
!      write(*,*)'dimA_RHS,B,T,IDCOL',size(u),size(ANZ(:,1)),
!     * size(IDCOL(:,2))
C      common /blk1/ anz,row,column,nnzero
C      IF(.NOT. ALLOCATED(SUM)) ALLOCATE(SUM(NNZMAX+NSZERO))
      IF(.NOT. ALLOCATED(SUM)) ALLOCATE(SUM(NEQ))
c
      SUM=0.D0
c
C      do i=1,nnzero
      DO I=1,NEQ
C       IF(COLUMN(I) .EQ. 0) CYCLE
C       sum(row(i))=sum(row(i))+anz(i)*v(column(i))
       IA(3)=IDCOL(I,1)
       IA(4)=IDCOL(I,2)
       SUM(I)=SUM(I)+B(I,1)*V(I)
       IF (IA(4)==0) CYCLE
       IF(IA(3)==0) THEN
          J1=4
!       ELSEIF (IA(3)==1 .AND. B(I,3).NE.0.D0) THEN
!          J1=3
       ELSE
          J1=2
       ENDIF
       DO J=J1,4
        K=I+IA(J)
        IF(K .LE. NEQ) THEN
           SUM(I)=SUM(I)+B(I,J)*V(K)
           SUM(K)=SUM(K)+B(I,J)*V(I)
        ENDIF
       ENDDO
      ENDDO
c
      DO I=1,NEQ
       U(I)=U(I)+SUM(I)
      ENDDO
c
      END
c
c     ----------------------------------------------------------------
c
      SUBROUTINE ATUPV(U,V)
      use M_L_EQ_SYS 
C      include 'ray.par'
c
C      parameter(nnzmax=2500000, ndmax=25000)
C      parameter(nsconst=3*nximax*nzimax,
C     +          nszero=7*nximax*nzimax,
C     +          nmmax=nximax*nzimax)
c
      IMPLICIT REAL*8 (A-H,O-Z)
C      integer row(nnzmax+nszero),column(nnzmax+nszero)
C      real anz(nnzmax+nszero),u(m),v(n),sum(nmmax)
C      DIMENSION ANZ(nnzmax+nszero),u(ndmax+nsconst),v(n)
      DIMENSION U(NEQ),V(NEQ),IA(4)
C      COMMON /blk1/ NNZERO,NNZMAX,NDMAX,NSZERO,NSCONST,NMMAX
      ALLOCATABLE :: SUM(:)
      data IA /0,1,0,0/
c
C      IF(.NOT. ALLOCATED(SUM)) ALLOCATE(SUM(nnzmax+nszero))
      IF(.NOT. ALLOCATED(SUM)) ALLOCATE(SUM(NEQ))
c
C      common /blk1/ anz,row,column,nnzero
c
      SUM=0.D0
c
      DO I=1,NEQ
C       IF(COLUMN(I) .LT. 1) CYCLE
       IA(3)=IDCOL(I,1)
       IA(4)=IDCOL(I,2)
       SUM(I)=SUM(I)+B(I,1)*U(I)
       IF (IA(4)==0) CYCLE
       IF(IA(3)==0) THEN
          J1=4
!       ELSEIF (IA(3)==1 .AND. B(I,3).NE.0.D0) THEN
!          J1=3
       ELSE
          J1=2
       ENDIF
       DO J=J1,4
        K=I+IA(J)
        IF(K .LE. NEQ) THEN
           SUM(K)=SUM(K)+B(I,J)*U(I)
           SUM(I)=SUM(I)+B(I,J)*U(K)
        ENDIF
C       sum(column(i))=sum(column(i))+anz(i)*u(row(i))
       ENDDO
      ENDDO
c
      DO J=1,NEQ
       V(J)=V(J)+SUM(J)
      ENDDO
c
      END

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




