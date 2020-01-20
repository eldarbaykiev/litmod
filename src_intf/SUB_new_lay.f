! Subroutine  SUB_new_lay.f

      SUBROUTINE new_lay
      use M_imag
      use M_profile
      use M_material
      use M_pol
      use M_label
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
!      write(*,*)'INSERT LAYER...'
      CALL PGSLCT(ISTAT)
      CALL lay_sel
      IF (LBL_SGN==1) RETURN
!      write(*,*)'From layer --->'
!      read(*,*)i_top
!      write(*,*)'to layer ----->'
!      read(*,*)i_bot 
      IF (i_bot<i_top) THEN 
      ins=i_top
      ins_mat=ins+1
      ELSE
      ins=i_bot
      ins_mat=ins
      ENDIF
      dum=prop
      dum(ins_mat+1:N_lay+1)=prop(ins_mat:N_lay)
      
      CALL body_prop
      IF (abort_body=='y') RETURN

      prop(ins_mat+1:N_lay+1)=dum(ins_mat+1:N_lay+1)

!      write(*,*)'Body properties'
!      write(*,*)'***************'
!      write(*,*)'BODY NAME'
!      read(*,*) prop(ins_mat)%name_mat
!      write(*,*)'DENSITY (kg/m³)'
!      read(*,*) prop(ins_mat)%rho
!      write(*,*)'ALPHA (1/k)'
!      read(*,*) prop(ins_mat)%alpha
!      write(*,*)'THERMAL CONDUCTIVITY'
!      read(*,*) prop(ins_mat)%K
!      write(*,*)'HEAT PRODUCTION RATE'
!      read(*,*) prop(ins_mat)%A
!      write(*,*)' PRESSURE COEFFICIENT'
!      read(*,*) prop(ins_mat)%beta
      dm_lay=layers
      IF( ALLOCATED(layers)) DEALLOCATE(layers)
      ALLOCATE (layers(IDIM,JDIM,N_lay+1))
      layers(:,:,1:N_lay)=dm_lay
      DEALLOCATE(dm_lay,GLUE_m)
      ALLOCATE (dm_lay(IDIM,JDIM,N_lay+1),
     *            GLUE_m(IDIM,JDIM,N_lay+1,N_lay+1))
      

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      DO i=N_lay,ins,-1
       IF (i>9) THEN
       WRITE(chdm1,'(I2)')i
       WRITE(chdm2,'(I2)')i+1
       ELSEIF (i==9) THEN
       WRITE(chdm1,'(I1)')i
       WRITE(chdm2,'(I2)')i+1
       ELSE
       WRITE(chdm1,'(I1)')i
       WRITE(chdm2,'(I1)')i+1
       ENDIF
       command='mv '//TRIM(string)//TRIM(chdm1)//'.xyz '//
     *         TRIM(string)//TRIM(chdm2)//'.xyz'
!       write(*,*)'command ',command
       call SYSTEM(TRIM(command))
      layers(:,:,i+1)=layers(:,:,i)
      ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      

      N_lay=N_lay+1

!Rewrite layers.info with the new layer 
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
!  3   FORMAT (I2,'  ',A20,F6.0,E12.3,F5.1,E10.2,E10.2,2F4.2,F5.1,I2)
   3   FORMAT (I2,' ',A20,1X,F6.0,E12.3,F5.1,2E12.3,2F6.2,1X,F5.1,1X,I2)
      close (15)
      
      CALL PGSLCT(ISTAT)
      CALL clear_scr  (0)
      CALL get_pol
      CALL in_out_pol ('A')    
      dm_lay=layers
 
      END SUBROUTINE
