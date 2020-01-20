! Subroutine SUB_profile.f

      SUBROUTINE profile
      USE M_imag
      USE M_profile
      USE M_material
      INTEGER il,sense
      INTEGER PGOPEN
      REAL(4) xdm,ydm,x1,y1,x2,y2
      REAL(4) R,G,B
      CHARACTER(2) chdm
      CHARACTER(1) le
D     CHARACTER*50 START,FINISH
!Bodies info
      open(15,file='layers.info',status='UNKNOWN')
      open(16,file='color.rgb',status='UNKNOWN')
      DO i=1,N_lay+1
       IF (i==1) THEN
       read(15,*)
       ELSE
       read(15,*)j,prop(j)%name_mat,prop(j)%rho,prop(j)%alpha,prop(j)%K,
     *           prop(j)%A,prop(j)%beta
!Create color file
       R=1.
!       G=(rho_max-prop(j)%rho)/(rho_max-rho_min)
       G=1-((ABS(prop(j)%rho)-rho_min)/(rho_max-rho_min))**2
       B=0.
       write(16,*)nint(R*255),nint(G*255),nint(B*255)
       ENDIF
      ENDDO
      close (15)
      close(16)

      MODE=0
      le='a'
      i=0
      DO WHILE (ICHAR(le)/=27)
D      OPEN(15,FILE='prof_plot.bat')
!      write(*,*)'***********************'
!      write(*,*)'Enter the extreme coordinates of the profile...'
       IF(i==0) CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,
     *                   lat_min,lat_max,2)
       i=i+1
       IF (i>=1) THEN
        CALL PGSLCT(ISTAT)
        CALL clear_scr(0)
       ENDIF
       xdm=(lon_max+lon_min)/2
       ydm=(lat_max+lat_min)/2
       call pgband(0,1,0.,0.,xdm,ydm,le)
       IF (ICHAR(le)==27) EXIT
       x1=xdm
       y1=ydm
       call pgband(1,1,xdm,ydm,xdm,ydm,le)
       IF (ICHAR(le)==27) EXIT
       x2=xdm
       y2=ydm
!Display profile location
!//////////////////////////////////////
      IF (x1<lon_min) x1=lon_min
      IF (x2<lon_min) x2=lon_min
      IF (x1>lon_max) x1=lon_max
      IF (x2>lon_max) x2=lon_max

      IF (y1<lat_min) y1=lat_min
      IF (y2<lat_min) y2=lat_min
      IF (y1>lat_max) y1=lat_max
      IF (y2>lat_max) y2=lat_max

       CALL PGSLW(10)
       CALL PGSCI(1)
       x_pl(1)=x1
       x_pl(2)=x2
       y_pl(1)=y1
       y_pl(2)=y2
       CALL PGLINE(2,x_pl,y_pl)
       CALL PGSLW(3)
!//////////////////////////////////////

!CHECK SCREEN/PS OUTPUT**********
!       IF (prf_act=='S') THEN
        N_prf=i
        CALL profile_screen
!        CYCLE
!       ENDIF
!CHECK SCREEN/PS OUTPUT**********

!       sense=0
D      WRITE(START,'(''-C'',F0.4,''/'',F0.4)') X1,Y1
D      WRITE(FINISH,'(''-E'',F0.4,''/'',F0.4,'' -G5'')') X2,Y2
D      CALL SYSTEM('GMT project '//TRIM(START)//' '//TRIM(FINISH)//
D    *            ' -Q > profile_project.xyp')
D      OPEN(55,FILE='profile_project.xyp')
D      DO
D       READ(55,*,IOSTAT=II) DUM1,DUM2,PER_L
D       IF(II .NE. 0) EXIT
D      ENDDO
D      CLOSE(55)
!       DO il=1,N_lay
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!       IF (il>9) THEN
!       WRITE(chdm,'(I2)')il
!       ELSE
!       WRITE(chdm,'(I1)')il
!       ENDIF
!       fil=TRIM(string)//TRIM(chdm)//'.xyz'
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        sense=1-sense
!        open(7,file='qq',status='UNKNOWN')
C        IF(LINUX) THEN
!         write(7,5) x1,y1,x2,y2,i,il,N_lay,TRIM(fil),sense
C        ELSE
C
C        ENDIF
!        close (7)
!  5     FORMAT (4F7.2,3I2,' ',A,' ',I1)
!        IF(LINUX) THEN
!           CALL SYSTEMQQ ('profile.job')
!        ELSE
!           CALL PROF_PREP(5.,170.,PER_L)
!        ENDIF
!      ENDDO
D      CLOSE(15)
D      CALL SYSTEMQQ('prof_plot.bat')
!       call pgband(0,1,0.,0.,dm,dm,le)
      ENDDO

      END SUBROUTINE
