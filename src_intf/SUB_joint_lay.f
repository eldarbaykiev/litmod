! Subroutine SUB_joint_lay.f


      SUBROUTINE joint_lay
      use M_imag
      use M_profile
      use M_material
      use M_pol 
      CHARACTER(1):: chdm,GETCHARQQ
      CHARACTER(50) UCASE
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
      dm_lay=layers 
      

      call pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               3)
!!      write(*,*)'PUT LAYER A AT DEPTH OF LAYER B'
!!      write(*,*)'LAYER A --->'
!!      read(*,*) ins
!!      write(*,*)'LAYER B---->'
!!      read(*,*)i_top
!!      write(*,*)'JOIN LAYERS INSIDE POLYGON,CHANGE (N=def, Y)?'
!!      READ(*,*)chdm
!!!      chdm=GETCHARQQ()
!!      IF(ICHAR(chdm)== 13) THEN
!!       chdm='N'
!!      ELSE
!!       chdm=UCASE(chdm)
!!      ENDIF

!!      IF (TRIM(chdm)=='Y') THEN
!!       yin_yang=.FALSE.
!!      ELSE
!!       yin_yang=.TRUE.
!!      ENDIF  
      IF (ins==-1000.OR.i_top==-1000) RETURN

      IF (i_top==0) THEN
       write(*,11)'Elevation  ',prop(ins)%name_mat
      ELSE 
       write(*,11)prop(ins)%name_mat,prop(i_top)%name_mat
      ENDIF
      IF (yin_yang) THEN
       WRITE(*,*)'JOIN LAYERS INSIDE POLYGON(S)' 
      ELSE
       WRITE(*,*)'JOIN LAYERS OUTSIDE POLYGON(S)'
      ENDIF
  11  FORMAT('PUT LAYER: ',A20,'AT DEPTH OF LAYER: ',A20) 
      CALL get_pol
      CALL in_out_pol ('J')
       

      END SUBROUTINE
