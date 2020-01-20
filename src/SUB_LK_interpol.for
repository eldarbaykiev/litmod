! Subroutine SUB_LK_interpol.for
 

      SUBROUTINE LK_interpol(T_pt,P_pt,litho_LK,rho_LK,v_p_LK,v_s_LK)
      USE M_mant
      REAL(8)T_pt,P_pt,rho_LK,v_p_LK,v_s_LK
      REAL(8)T_1,T_2,P_1,P_2,pmt_LK
      INTEGER litho_LK
!Index of mantle material in the look up table vector
      DO i_LK=1,I_mnt
       IF(litho2LK(i_LK)==litho_LK) EXIT
      ENDDO  
!      write(*,*)'i_LK',i_LK,litho2LK(i_LK),T_pt
!      i_LK=mat_z_LK-I_crst+1
!      IF(mat_z_LK==0) i_LK=1
!Find the 4 nearest neighbours
      i_1=INT((T_pt-T_LK_0(i_LK))/d_T_LK(i_LK)+1)
      j_1=INT((P_pt-P_LK_0(i_LK))/d_P_LK(i_LK)+1)
      i_2=i_1+1
      j_2=j_1+1
      T_1=(i_1-1)*d_T_LK(i_LK)+T_LK_0(i_LK)
      T_2=(i_2-1)*d_T_LK(i_LK)+T_LK_0(i_LK) 
      P_1=(j_1-1)*d_P_LK(i_LK)+P_LK_0(i_LK)
      P_2=(j_2-1)*d_P_LK(i_LK)+P_LK_0(i_LK)
! Determine the nearest node of the four neighbours
      IF (ABS(T_2-T_pt)<ABS(T_1-T_pt)) THEN
       i_near=i_2
      ELSE
       i_near=i_1
      ENDIF
      IF (ABS(P_2-P_pt)<ABS(P_1-P_pt)) THEN
       j_near=j_2
      ELSE
       j_near=j_1
      ENDIF

!Check for T and P values out of the range of the compositional files

      IF(i_1<0) THEN
       write(*,*)''
       write(*,*)'TRYING TO USE AN ANOMALOUS SUB-LITHOSPHERIC MANTLE COM
     *POSITIONAL FILE IN THE LITHOSPHERE: '
       write(*,*) TRIM(mant_file_sys(LK_index(i_LK)))
      ELSEIF(i_2>N_LK_T(i_LK)) THEN
       write(*,*)''
       write(*,*)'TRYING TO RETRIEVE TOO HIGH T FROM THE COMPOSITIONAL 
     *FILE: ' 
       write(*,*) TRIM(mant_file_sys(LK_index(i_LK)))
       write(*,*)'T(C) = ',T_pt
       write(*,*)'T max COMPOSITIONAL FILE(C) =', 
     *           (N_LK_T(i_LK)-1)*d_T_LK(i_LK)+T_LK_0(i_LK)
      ELSEIF(ABS(j_2)>N_LK_P(i_LK)) THEN 
       write(*,*)''
       write(*,*)'TRYING TO RETRIEVE TOO HIGH P FROM THE COMPOSITIONAL 
     *FILE: '
       write(*,*) TRIM(mant_file_sys(LK_index(i_LK)))
       write(*,*)'P(GPa) = ',P_pt*1D-9
       write(*,*)'P max COMPOSITIONAL FILE(GPa) =',
     *           ((N_LK_P(i_LK)-1)*d_P_LK(i_LK)+P_LK_0(i_LK))*1D-9
       write(*,*)'CHECK THE PRESSURE COEFFICIENT OF THE CRUSTAL BODIES'
      ENDIF

      IF(i_1<0.OR.i_2>N_LK_T(i_LK).OR.ABS(j_2)>N_LK_P(i_LK)) THEN
       write(*,*)'LitMod3D aborts...'
       write(*,*)''
       STOP
      ENDIF 

      DO k=3,5 
       pmt_LK=(LK_vec(i_LK,i_1,j_1,k)*(T_2-T_pt)*(P_2-P_pt)+
     *        LK_vec(i_LK,i_2,j_1,k)*(T_pt-T_1)*(P_2-P_pt)+
     *        LK_vec(i_LK,i_1,j_2,k)*(T_2-T_pt)*(P_pt-P_1)+
     *        LK_vec(i_LK,i_2,j_2,k)*(T_pt-T_1)*(P_pt-P_1))/
     *        ((T_2-T_1)*(P_2-P_1))
       IF(k==3) rho_LK=pmt_LK
       IF(k==4) v_p_LK=pmt_LK
       IF(k==5) v_s_LK=pmt_LK
      ENDDO
      END SUBROUTINE
