! SUB_smooth_lay.f
      


      SUBROUTINE smooth_lay
      use M_imag
      use M_profile
      INTEGER::I_win
      real(4) CR(4),DG(4),w_in,w_cr,w_dg
      real(4) x,y
      character(2) chdm
!      CALL lay_sel 
      
      I_win=3
      DO k=1,N_lay
       DO j=JDIM-I_win,1+I_win,-1
        DO i=I_win+1,IDIM-I_win
         layers(i,j,k)=SUM(layers(i-I_win:i+I_win,j-I_win:j+I_win,k))
     *                    /(2*I_win+1)**2
        ENDDO
       ENDDO
      ENDDO
   
      DO k=0,N_lay-1
       DO j=JDIM,1,-1
        DO i=1,IDIM
         IF (k==0) THEN
          IF (ABS(-E(i,j)*1e-3-layers(i,j,k+1))<.5) THEN
           layers(i,j,k+1)=-E(i,j)*1e-3 
          ENDIF
         ELSE  
          IF (ABS(layers(i,j,k)-layers(i,j,k+1))<.5) THEN
           layers(i,j,k+1)=layers(i,j,k)      
          ENDIF
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      DO i=1,N_lay
      IF (i>9) THEN
       WRITE(chdm,'(I2)')i
       ELSE
       WRITE(chdm,'(I1)')i
       ENDIF
       fil=TRIM(string)//TRIM(chdm)//'.xyz'
       open(20,file=fil,status='UNKNOWN') 
      do ky=JDIM,1,-1
        y=(ky-1)*p_m_y+lat_min
         do kx=1,IDIM
          x=(kx-1)*p_m_x+lon_min
          write (20,*)x,y,layers(kx,ky,i)
         enddo
       enddo
      close(20)
      ENDDO
      END SUBROUTINE
