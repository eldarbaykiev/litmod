! SUB_body_plot.f


      SUBROUTINE body_plot
      USE M_imag
      USE M_profile
      USE M_material
      USE M_label
      INTEGER PGOPEN
      CHARACTER LE
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

!SELECT BODY
      CALL pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               -3)
      CALL PGSLCT(ISTAT)
       
      out_win=.FALSE.
      DO
       IF (out_win) EXIT
       CALL lay_sel
       IF (LBL_SGN==1) RETURN
      ENDDO

!!      IF(LINUX) THEN
!!         IPG_LAY=PGOPEN('/XWINDOW')
!!      ELSE
!!         IPG_LAY=PGOPEN('/WX')
!!      ENDIF
!!      IF (IPG_LAY<=0) STOP
!!      CALL pgask(.false.)
      CALL plot

!!      CALL pgqvp(0,r1,r2,r3,r4) 
!!      CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max
!!     *              ,12)
!!      CALL pgband(0,1,0.,0.,xdm,ydm,le)
!!      CALL PGCLOS
      CALL PGSLCT(ISTAT)
!Erase the title of the plot
       CALL PGSCI(0)
       CALL PGSVP(.27,.9,.93,1.)
       CALL PGSWIN(0.,1.,0.,1.)
       CALL PGRECT(0.,1.,0.,1.)
!Erase the palette
       CALL PGSVP(.9,1.,.1,1.)
       CALL PGSWIN(0.,1.,0.,1.)
       CALL PGRECT(0.,1.,0.,1.)
      IMG=bod_tk
      IMG_txt=prop(I_BODY)%name_mat 
!      z_min=minval(bod_tk)
      z_min=0.
      z_max=maxval(bod_tk)
      CALL PALETT(6,-1.,.5)
      CALL clear_scr(0) 
!Plot the palette
       CALL PGWEDG('RI',2.5,3.,0.,z_max,'(km)')

      END SUBROUTINE
