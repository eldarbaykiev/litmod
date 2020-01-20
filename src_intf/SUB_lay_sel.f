! SUB_lay_sel.f


      SUBROUTINE lay_sel
      USE M_imag
      USE M_profile
      USE M_material
      USE M_pol
      USE M_label
      USE M_visualize, ONLY:z_slice,Nz,dz  
      USE M_sublitho_an, ONLY: z_an_up, z_an_dw
      real(4), allocatable:: z_lay (:)
      character(1) le
      character(80) TEXT  
      real(4) xdm,ydm,zdm,z_text,x_pos_txt,z_lst
      real(4),dimension(2) :: x_2p,y_2p
      integer ipt,jpt
      INTEGER PGOPEN,NMARK_old 

      out_win=.FALSE.
      LBL_SGN=0
      x_2p(1)=0
      x_2p(2)=1
      NMARK_old=NMARK 
      IF(ALLOCATED(z_lay)) DEALLOCATE(z_lay) 
      ALLOCATE (z_lay (0:N_lay)) 
      DO 
      CALL PGSLCT(ISTAT)
       CALL clear_scr(0)
       CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -3)
       xdm=(lon_max+lon_min)/2
       ydm=(lat_max+lat_min)/2 
       CALL pgband(0,1,xdm,ydm,xdm,ydm,le) 
       IF (ICHAR(le)==27) THEN
         LBL_SGN=1
         IF (act=='L') THEN
          NMARK=NMARK_old-1
         ENDIF
        RETURN
       ELSEIF (ICHAR(le)==13.AND.act/='M'.AND.act/='B') THEN
         LBL_SGN=1
        RETURN
       ENDIF
       ipt=nint((xdm-lon_min)*IDIM/(lon_max-lon_min))
       jpt=nint((ydm-lat_min)*JDIM/(lat_max-lat_min))
       IF (ipt<1.OR.ipt>IDIM.OR.jpt<1.OR.jpt>JDIM) THEN
        LBL_SGN=1 
        RETURN
       ELSE
        out_win=.TRUE.   
       ENDIF

       IF (act=='L') THEN
        NTICK=NMARK+1
        labels(NTICK)%long=xdm
        labels(NTICK)%lat=ydm
       ENDIF
 
       DO i=0,N_lay
        IF (i==0) THEN
           z_lay(i)=E(ipt,jpt)*1e-3
        ELSE
           z_lay(i)=-layers(ipt,jpt,i)
        ENDIF
       ENDDO
       IF(LINUX) THEN
          IS_LAY=PGOPEN('/XWINDOW')
       ELSE
          IS_LAY=PGOPEN('/WX')
       ENDIF
       IF (IS_LAY<=0) STOP
       call pgscr(0,1.,1.,1.)
       call pgscr(1,0.,0.,0.)
       call pgask(.false.)
!       CALL PGENV(0.,1.,-180.,5.,0,-2)
       CALL pgsvp(.3,.9,.1,.9)
       IF(act=='V'.OR.act=='U') THEN
        z_lst=-(Nz-1)*dz*1E-3
       ELSEIF (act=='F') THEN
        z_lst=-MAXVAL(layers(:,:,N_cr))
       ELSE
        z_lst=-MAXVAL(layers(:,:,N_lay)) 
       ENDIF
       CALL pgswin(0.,1.,z_lst,5.)
       CALL PGBOX('bi',0.0,0,'bctsivm',0.0,0)
       IF(act=='L') THEN
       CALL PGLAB('', 'DEPTH','INTRODUCE REFERENCE LABELS')
       ELSE 
       CALL PGLAB('', 'DEPTH','GET LAYER')
       CALL text_plot(IS_LAY,.3,.9,.1,.9,0.,1.,z_lst,5.,-10)   
       ENDIF
 
!Plot layers
       DO i=0,N_lay
        y_2p= z_lay(i)
        IF (i==0) THEN
           CALL PGPTEXT(1.,y_2p(2)+5,0.,1.,'ELEVATION')
           CALL PGSLS(5)
        ELSE
           IF (z_lay(i-1)-z_lay(i)<5e-3) CYCLE
           write(TEXT,11)TRIM(prop(i)%name_mat),i
 11        FORMAT (A,'--> ',I2)
           z_text=(z_lay(i)+z_lay(i-1))/2 
           IF(act=='L') THEN
            x_pos_txt=.8
           ELSE
            x_pos_txt=.5
           ENDIF 
           CALL PGPTEXT(x_pos_txt,z_text,0.,0.5,TEXT)
           CALL PGSLS(1)
        ENDIF
        CALL PGLINE(2,x_2p,y_2p)
       ENDDO 
!CALL FROM SUBROUTINE label ******************** 
       IF(act=='L') THEN
       CALL text_plot(IS_LAY,.3,.9,.1,.9,0.,1.,z_lst,5.,-15)
       DEALLOCATE (z_lay)
       CALL PGCLOS
       RETURN
       ENDIF
!*********************************************
       xdm=.5
       zdm=-70
       call pgband(5,0,xdm,zdm,xdm,zdm,le)

!CALL FROM  SUBROUTINE visualize ******************** 
       IF(act=='V') THEN
       z_slice=-zdm*1E3
       DEALLOCATE (z_lay)
       CALL PGCLOS
       RETURN
       ENDIF
!*********************************************

!CALL FROM  SUBROUTINE sublitho_an ******************** 
       IF(act=='U') THEN
       IF (ICHAR(le)==27) THEN
        LBL_SGN=1
        RETURN
       ENDIF
       z_an_up=-zdm*1E3
       CALL pgband(5,0,xdm,zdm,xdm,zdm,le)
       IF (ICHAR(le)==27) THEN
        LBL_SGN=1
        RETURN
       ENDIF
       z_an_dw=-zdm*1E3
!       CALL text_plot(IS_LAY,.3,.9,.1,.9,0.,1.,z_lst,5.,-16)
       DEALLOCATE (z_lay)
       CALL PGCLOS
       RETURN
       ENDIF
!*********************************************
      
       DO i=0,N_lay-1
        IF (z_lay(i)>zdm.AND.z_lay(i+1)<zdm) THEN
           IF ((z_lay(i)-zdm)<zdm-z_lay(i+1)) THEN
              i_lay=i
              i_top=i
              i_bot=i+1
              IF(act=='B') THEN   
                 write(TEXT_BDSEL,14) TRIM(prop(i+1)%name_mat)
              ELSE
                 write(TEXT_BDSEL,12) TRIM(prop(i+1)%name_mat)
              ENDIF
 12           FORMAT ('BODY: ',A,' --> Upper limit')
 14           FORMAT ('BODY: ',A) 
              EXIT
           ELSE
              i_lay=i+1
              i_top=i+1
              i_bot=i 
              IF(act=='B') THEN
                 write(TEXT_BDSEL,14) TRIM(prop(i+1)%name_mat)
              ELSE
                 write(TEXT_BDSEL,13) TRIM(prop(i+1)%name_mat)
              ENDIF
 13           FORMAT ('BODY: ',A,' --> Lower limit')
              EXIT
           ENDIF
        ELSEIF (zdm<z_lay(N_lay)) THEN
           i_lay=N_lay 
           i_top=N_lay
           i_bot=N_lay-1
           IF(act=='B') THEN
              TEXT_BDSEL='NO BODY SELECTED !!'
              i_lay=-1
              EXIT
           ELSE
              write(TEXT_BDSEL,13) TRIM(prop(N_lay)%name_mat)
           ENDIF
            EXIT
        ELSEIF (zdm>z_lay(0)) THEN
!        write(*,*)'!!!!!!!!!!!!!!!!!'
!        write(*,*)'NO LAYER SELECTED'
!        write(*,*)'!!!!!!!!!!!!!!!!!'
           TEXT_BDSEL='NO LAYER SELECTED !!' 
           i_lay=-1 
           EXIT 
        ENDIF
       ENDDO
       CALL text_plot(IS_LAY,.3,.9,.1,.9,0.,1.,z_lst,5.,4)
       IF (i_lay==-1) THEN
         CALL PGCLOS
         CYCLE
       ENDIF
!      write(*,*)'OK ...? (y/n)' 
!      read(*,*) le
       IF (bod_ok=='Y') THEN
          EXIT
       ELSE
          CALL PGCLOS
       ENDIF

      ENDDO 
      DEALLOCATE (z_lay)
      CALL PGCLOS
      END SUBROUTINE  
