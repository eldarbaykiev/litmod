!Subroutine  SUB_delete_lay.f


      SUBROUTINE delete_lay
      USE M_imag
      USE M_material
      USE M_label
      USE M_material  
      USE M_profile
      TYPE (MATERIAL) dum (-1:100)
      CHARACTER(100) :: command
      CHARACTER(2) chdm1,chdm2
      open(15,file='layers.info',status='UNKNOWN')
!Bodies info
      DO i=1,N_lay+1
       IF (i==1) THEN
       read(15,*)
       ELSE
       read(15,*)j,prop(j)%name_mat,prop(j)%rho,
     *                      prop(j)%alpha,prop(j)%K,
     *                      prop(j)%A,prop(j)%beta,prop(j)%gmma_T,
     *                      prop(j)%K_0p,prop(j)%K_T,prop(j)%litho
       ENDIF
      ENDDO
      rewind (15)
      call pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               -3)
      CALL PGSLCT(ISTAT)
      CALL lay_sel
      IF (LBL_SGN==1) RETURN     
      dum=prop
      IF (i_top>i_bot) THEN
       prop(i_lay:N_lay-1)=dum(i_lay+1:N_lay)
      ELSE
       prop(i_lay+1:N_lay-1)=dum(i_lay+2:N_lay)
      ENDIF
      dm_lay=layers
      IF( ALLOCATED(layers)) DEALLOCATE(layers)
      ALLOCATE (layers(IDIM,JDIM,N_lay-1))
      IF (i_lay>1) layers(:,:,1:i_lay-1)=dm_lay(:,:,1:i_lay-1)
      layers(:,:,i_lay:N_lay-1)=dm_lay(:,:,i_lay+1:N_lay)
      DEALLOCATE(dm_lay,GLUE_m)
      ALLOCATE (dm_lay(IDIM,JDIM,N_lay-1),
     *            GLUE_m(IDIM,JDIM,N_lay-1,N_lay-1))
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      DO i=i_lay+1,N_lay
       IF (i>10) THEN
       WRITE(chdm1,'(I2)')i
       WRITE(chdm2,'(I2)')i-1
       ELSEIF (i==10) THEN
       WRITE(chdm1,'(I2)')i
       WRITE(chdm2,'(I1)')i-1
       ELSE
       WRITE(chdm1,'(I1)')i
       WRITE(chdm2,'(I1)')i-1
       ENDIF
       command='mv '//TRIM(string)//TRIM(chdm1)//'.xyz '//
     *         TRIM(string)//TRIM(chdm2)//'.xyz'
!       write(*,*)'command ',command
       call SYSTEM(TRIM(command))
      ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
      N_lay=N_lay-1
!Rewrite layers.info without the erased layer 
      DO i=0,N_lay
       IF (i==0) THEN
       write(15,2)"'BODY NUMBER | DENSITY | ALPHA | THERMAL CONDUCTIVITY
     * | HEAT PRODUCTION RATE | PRESSURE COEFFICIENT | GRÜNEISEN
     * PARAMETER | PRESSURE DERIVATIVE OF ISOTHERMAL BULK MODULUS |
     * ISOTHERMAL BULK MODULUS | MANTLE INDEX'"
       ELSE
       write(15,3)i,(prop(i)%name_mat),prop(i)%rho,prop(i)%alpha,
     *            prop(i)%K,prop(i)%A,prop(i)%beta,prop(i)%gmma_T,
     *            prop(i)%K_0p,prop(i)%K_T,prop(i)%litho
       ENDIF
      ENDDO
  2   FORMAT (A215)
  3   FORMAT (I2,' ',A20,1X,F5.0,E12.3,F5.1,2E12.3,2F6.2,1X,F5.1,1X,I2)
      close (15)      
      CALL PGSLCT(ISTAT)
      CALL clear_scr (0) 

      END SUBROUTINE 


 
