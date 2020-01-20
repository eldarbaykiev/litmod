! SUB_body_plot_ps.f

      SUBROUTINE body_plot_ps
      use M_imag
      use M_profile
      use M_material
      real(4) x_0,y_0,x,y
      
!Bodies info

      open(15,file='layers.info',status='UNKNOWN')
      DO i=1,N_lay+1
       IF (i==1) THEN
       read(15,*)
       ELSE
       read(15,*)j,prop(j)%name_mat,prop(j)%rho,prop(j)%alpha,prop(j)%K,
     *           prop(j)%A,prop(j)%beta
       ENDIF
      ENDDO
      close (15)

      IF( ALLOCATED(i_bod_plot)) DEALLOCATE(i_bod_plot)
      ALLOCATE (i_bod_plot(N_lay))
!      DO i=1,N_lay
!      write(*,*)'PLOT BODY (Y/N)... ?'
!      write(*,*) prop(i)%name_mat
!      read(*,*) i_bod_plot(i)
      call pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               -4)  
!      ENDDO

      x_0=lon_min
      y_0=lat_min
      IF(ALLOCATED(bod_tk)) DEALLOCATE (bod_tk)
      ALLOCATE (bod_tk(IDIM,JDIM))
      DO i_lay=1,N_lay
       open(27,file='bod_thick.xyz',status='UNKNOWN')
       IF (i_bod_plot(i_lay)/='Y') CYCLE
       write(*,*)'BODY---> ',prop(i_lay)%name_mat
       write(*,*)'*****************'
       IF (i_lay==1) THEN
          bod_tk=layers(:,:,i_lay)+E*1e-3
       ELSE 
          bod_tk=layers(:,:,i_lay)-layers(:,:,i_lay-1)
       ENDIF
       DO jj=JDIM,1,-1
        y=(jj-1)*p_m_y+y_0
        DO ii=1,IDIM
         x=(ii-1)*p_m_x+x_0
         write(27,*) x,y,bod_tk(ii,jj)
        ENDDO
       ENDDO
       close (27)
       open(7,file='qq',status='UNKNOWN')
       HMIN=FLOOR(MINVAL(BOD_TK))
       HMAX=CEILING(MAXVAL(BOD_TK))
       write(7,6) HMIN,HMAX,prop(i_lay)%name_mat
  6    FORMAT (2F7.2,' ',A15)
       close (7)

!       write(*,*)'bodplot',minval(bod_tk),maxval(bod_tk),
!     *           prop(i_lay)%name_mat
       IF(LINUX) THEN
          CALL SYSTEM ('body_plot.job')
       ELSE
          CALL BODY_PREP(I_LAY)
       ENDIF
      ENDDO
      DEALLOCATE (bod_tk)

      END SUBROUTINE
      

