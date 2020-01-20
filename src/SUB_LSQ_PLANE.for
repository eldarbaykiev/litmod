!Subroutine:  LSQ_PLANE

    
      subroutine LSQ_PLANE ()
      use M_grid_parameters
      use M_columns
      use M_lsq_plane
      allocate (Dt(3,N_x*N_y),GEOID_vec(N_x*N_y),GEOID_OBS_vec(N_x*N_y),
     *          Z_pl(N_x*N_y),GEOID_vec_cr(N_x*N_y),Z_pl_cr(N_x*N_y))
      S_x=d_x*(N_x-1)*N_x*.5D0
      S_y=d_y*(N_y-1)*N_y*.5D0
      S_x2=d_x*d_x*(N_x*N_x*(N_x-2)+N_x*(N_x+1)*.5D0)/3D0
      S_y2=d_y*d_y*(N_y*N_y*(N_y-2)+N_y*(N_y+1)*.5D0)/3D0 
      S_xy=S_x*S_y
      CF_pl=0D0
      Ind=1
      DO I_x=1,N_x
       x=d_x*(I_x-1)
       DO I_y=1,N_y
        y=d_y*(I_y-1)
        Dt(1,Ind)=x
        Dt(2,Ind)=y
        GEOID_vec(Ind)=GEOID(I_x,I_y)
        GEOID_OBS_vec(Ind)=GEOID_OBS(I_x,I_y)
!Write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$
       GEOID_vec_cr(Ind)=GEOID_cr(I_x,I_y) 
!Write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$
        Ind=Ind+1
       ENDDO
      ENDDO
      Dt(3,:)=1D0
      DtD_1(1,1)=N_x*N_x*(N_y*S_y2-S_y*S_y)
      DtD_1(1,2)=N_x*N_y*(S_x*S_y-S_xy)
      DtD_1(1,3)=N_x*S_xy*S_y-N_x*N_y*S_y2*S_x
      DtD_1(2,2)=N_y*N_y*(N_x*S_x2-S_x*S_x)
      DtD_1(2,3)=N_y*(S_xy*S_x-N_x*S_x2*S_y) 
      DtD_1(3,3)=N_x*N_y*S_x2*S_y2-S_xy*S_xy
      DtD_1(2,1)=DtD_1(1,2)
      DtD_1(3,1)=DtD_1(1,3)
      DtD_1(3,2)=DtD_1(2,3)
    
      det_D=N_x*N_y*(N_x*N_y*S_x2*S_y2+S_x*S_y*S_xy-N_y*S_y2*
     *      S_x*S_x-N_x*S_x2*S_y*S_y) 
      DtD_1=DtD_1/det_D
      DO i=1,3
       ddmy=0D0
       ddmy_o=0D0
       ddmy_cr=0D0
       DO j=1,N_x*N_y
        ddmy=ddmy+Dt(i,j)*GEOID_vec(j)
        ddmy_o=ddmy_o+Dt(i,j)*GEOID_OBS_vec(j)
        ddmy_cr=ddmy_cr+Dt(i,j)*GEOID_vec_cr(j)
       ENDDO
       DMY(i)=ddmy      
       DMY_O(i)=ddmy_o
       DMY_cr(i)=ddmy_cr
      ENDDO
      DO i=1,3
       DO j=1,3
        CF_pl(i)=CF_pl(i)+DtD_1(i,j)*DMY(j)
        CF_pl_O(i)=CF_pl_O(i)+DtD_1(i,j)*DMY_O(j)
        CF_pl_cr(i)=CF_pl_cr(i)+DtD_1(i,j)*DMY_cr(j)
       ENDDO   
      ENDDO
!      write(*,*)'CF_pl,CF_pl_O',CF_pl,CF_pl_O
      DO j=1,N_x*N_y
       ddmy=0D0
       ddmy_cr=0D0
       DO i=1,3
       ddmy=ddmy+Dt(i,j)*(CF_pl(i)-CF_pl_O(i))
!!        ddmy=ddmy+Dt(i,j)*CF_pl_O(i)
       ddmy_cr=ddmy_cr+Dt(i,j)*(CF_pl_cr(i)-CF_pl_O(i)) 
       ENDDO
       Z_pl(j)=ddmy
       Z_pl_cr(j)=ddmy_cr
       GEOID_vec(j)=GEOID_vec(j)-Z_pl(j)
!Write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$
       GEOID_vec_cr(j)=GEOID_vec_cr(j)-Z_pl_cr(j)
!Write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$
      ENDDO

      end subroutine  
    
