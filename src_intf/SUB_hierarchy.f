! SUB_hierarchy.f


      SUBROUTINE hierarchy
      USE M_imag
      USE M_profile
      USE M_material
!      INTEGER:: i_up,i_dw
      CHARACTER(50) UCASE
      CHARACTER(2) chdm !,GETCHARQQ,hi_ch
      REAL(4)::L_top,L_bot,x,y 
 

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
      
      call pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               0)


!       write(*,*)'Upper layer number: '
!       read(*,*) i_up_h
!       write(*,*)'Lower layer number: '
!       read(*,*) i_dw_h
!       write(*,*)'For hierarchy violations, set both layers'
!       write(*,*)'to the depth of upper (u) or lower (l) layer...?'
!!       hi_ch=GETCHARQQ()
!       READ(*,*) hi_ch
!       hi_ch=UCASE(hi_ch) 
!       WRITE(*,*)'GLUE OPTION OFF, CHANGE ? (N==def, Y)?'
!        READ(*,*)chdm
!        IF(ICHAR(chdm)== 13) THEN
!         chdm='N'
!        ELSE
!         chdm=UCASE(chdm)
!        ENDIF
!       chdm=GETCHARQQ()
!       chdm=UCASE(chdm)
 
       IF (i_up_h==-1000.OR.i_dw_h==-1000) RETURN       

       IF (i_up_h==0.AND.hi_ch=='L') THEN
        WRITE(*,*)'THAT WOULD MEAN CHANGING ELEVATION LAYER!!'
        WRITE(*,*)'IF YOU REALLY WANT TO DO THAT, USE MODIFY_REG OPTION'
        RETURN  
       ENDIF
       
!       IF (chdm=='Y') THEN
!        glue_l=.TRUE.
!       ELSE
!        glue_l=.FALSE.
!       ENDIF

       IF (i_up_h>i_dw_h.OR.i_dw_h>N_lay) THEN
        WRITE(*,*)'ILLEGAL ACTION, NOTHING WILL BE DONE'
        RETURN
       ENDIF

       IF(hi_ch/='U'.AND.hi_ch/='L') THEN
        WRITE(*,*)'ILLEGAL ACTION, NOTHING WILL BE DONE'
        RETURN
       ELSEIF (hi_ch=='U') THEN
        WRITE(*,*)'SET BOTH LAYERS TO THE DEPTH OF LAYER:',i_up_h,
     *            TRIM(prop(i_up_h)%name_mat)
       ELSEIF (hi_ch=='L') THEN
        WRITE(*,*)'SET BOTH LAYERS TO THE DEPTH OF LAYER:',i_dw_h,
     *            TRIM(prop(i_dw_h)%name_mat) 
       ENDIF


       IF (glue_l)  CALL glue('C',0,0)
       DO j=JDIM,1,-1
        DO i=1,IDIM
         IF (i_up_h==0) THEN
          L_top=-E(i,j)*1e-3
         ELSE
          L_top=layers(i,j,i_up_h)
         ENDIF
         L_bot=layers(i,j,i_dw_h)
         IF (L_bot<L_top) THEN
          IF (hi_ch=='U') THEN
           IF (i_up_h==0) THEN
            layers(i,j,i_dw_h)=-E(i,j)*1e-3
           ELSE
            layers(i,j,i_dw_h)=layers(i,j,i_up_h)
           ENDIF
           i_lay=i_dw_h
          ELSEIF(hi_ch=='L') THEN
           layers(i,j,i_up_h)=layers(i,j,i_dw_h)
           i_lay=i_up_h 
          ENDIF
         IF (glue_l) CALL glue('U',i,j)
         ENDIF
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
