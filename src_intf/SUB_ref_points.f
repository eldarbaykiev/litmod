! Subroutine  SUB_ref_points.f

      SUBROUTINE ref_points
      USE M_imag
      USE M_pol
      USE M_ref_points

      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -7)
      IF (RPTS_act=='R') THEN 
      OPEN (40,file='ref_points.info',status='UNKNOWN')
      OPEN (45,file='ref_points.dat',status='UNKNOWN')
      N_RPTS=0
      ENDIF

      CALL get_pol
      pol_RPTS(2*N_RPTS+1:2*N_RPTS+2*n_pol,:)=A_pol(:,:)
      N_pol_RPTS(N_RPTS+1:N_RPTS+n_pol)=NP_POL(:)

!Write new reference points in ref_points.info and ref_points.dat
      DO i=N_RPTS+1,(N_RPTS+1)+(n_pol-1)
       WRITE(40,*)N_pol_RPTS(i)
       DO j=1,N_pol_RPTS(i)
         WRITE(45,*)pol_RPTS(2*i-1,j),pol_RPTS(2*i,j)
       ENDDO
       pol_RPTS(2*i-1,N_pol_RPTS(i)+1)=pol_RPTS(2*i-1,1)
       pol_RPTS(2*i,N_pol_RPTS(i)+1)=pol_RPTS(2*i,1) 
       ENDDO
       N_RPTS=N_RPTS+n_pol
      END SUBROUTINE
