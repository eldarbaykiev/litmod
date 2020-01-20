! SUB_histogram.f


      SUBROUTINE histogram
      USE M_imag
      REAL(4), ALLOCATABLE:: VALUES(:)
      REAL(4) hst_min,hst_max,AV,std_dev,dm
      INTEGER::PGOPEN,N_hist
      CHARACTER(1) ps
      CHARACTER(50) UCASE
      CHARACTER(80) TEXT


      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -17)
      CALL pgband(0,1,xdm,ydm,xdm,ydm,ps)
      ps=UCASE(ps)
      N_hist=IDIM*JDIM
      IF(ALLOCATED(VALUES)) DEALLOCATE (VALUES)
      ALLOCATE (VALUES(N_hist))
!      write(*,*) N_hist
      VALUES=0D0
      i_con=2
      IF (ps=='Y') THEN
       DO i_con=2,5
        OPEN(30,file='./OUT_xyz/hist_obs.dat',STATUS='UNKNOWN')
! Compute average and standad deviation
        AV=SUM(OBS(i_con,:,:))/N_hist
        dm=0
         DO i=1,IDIM
          DO j=1,JDIM
           dm=dm+(OBS(i_con,i,j)-AV)**2
!           WRITE(30,*)OBS(2:5,i,j)
          ENDDO
         ENDDO
         std_dev=SQRT(dm/N_hist) 
        WRITE(30,*)AV,std_dev

       ENDDO
         DO i=1,IDIM
          DO j=1,JDIM
           WRITE(30,*)OBS(2:5,i,j)
          ENDDO
         ENDDO
         CLOSE(30)
         CALL SYSTEM ('histogram.job')
         CALL PGSLCT(ISTAT)
         CALL clear_scr(0) 
         RETURN
      ENDIF

      DO
      IF (i_con>5)i_con=2
      DO i=1,IDIM
       DO j=1,JDIM
          VALUES((i-1)*JDIM+j)=OBS(i_con,i,j)
       ENDDO
      ENDDO
! Compute average and standad deviation
      AV=SUM(VALUES(:))/N_hist
      dm=0
      DO i=1,N_hist
       dm=dm+(VALUES(i)-AV)**2
      ENDDO
      std_dev=SQRT(dm/N_hist)
      hst_min=MINVAL(VALUES)
      hst_max=MAXVAL(VALUES)


      IF(LINUX) THEN
         IHIST=PGOPEN('/XWINDOW')
!       IHIST=PGOPEN('histo.ps/vcps')
      ELSE
         IHIST=PGOPEN('/WX')
      ENDIF
      IF (IHIST.LE.0) STOP
      CALL pgscr(0,1.,1.,1.)
      CALL pgscr(1,0.,0.,0.)
      CALL pgask(.false.)
      CALL PGPAP(8.0,.8)
      CALL PGSCI(1)
      CALL PGSVP(.3,.9,.1,.9)
      rh_min=0D0
      rh_max=N_hist/8
      CALL text_plot(IHIST,.3,.9,.1,.9,hst_min, hst_max,rh_min,
     *                rh_max,-13)
!      write(*,*)hst_min,hst_max,MAXVAL(OBS(i_con,:,:))
!PGFLAG (input)  : if PGFLAG = 1, the histogram is plotted in the
!                   current window and viewport; if PGFLAG = 0,
!                   PGENV is called automatically by PGHIST to start
!                   a new plot (the x-limits of the window will be
!                   DATMIN and DATMAX; the y-limits will be chosen
!                   automatically.
!                   IF PGFLAG = 2,3 the histogram will be in the same
!                   window and viewport but with a filled area style.
!                   If pgflag=4,5 as for pgflag = 0,1, but simple
!                   line drawn as for PGBIN
      CALL PGHIST(N_hist, VALUES, hst_min, hst_max, 50,1)
!      call pgqvp(0,r1,r2,r3,r4)
!      write(*,*)r1,r2,r3,r4
!       call pgqwin(r1,r2,r3,r4)
!      write(*,*)r1,r2,r3,r4
      CALL PGBOX('bctsinm',0.0,0,'bctsivm',0.0,0)
      IF(i_con==2.OR.i_con==5) THEN
        TEXT='Calculated-observed (m)'
      ELSE 
         TEXT='Calculated-observed (mGal)'
      ENDIF 
      CALL PGLAB(TRIM(TEXT),
     1           'Number of data points',
     2           TRIM(obs_file(i_con)) )
       WRITE(TEXT,12) AV
 12    FORMAT('Average: ',F10.3)  
       CALL PGPTEXT(.0,rh_max-15,.9,0.5,TRIM(TEXT))
       WRITE(TEXT,13) std_dev
 13    FORMAT('Std dev: ',F10.3)
       CALL PGPTEXT(.0,rh_max-30,.7,0.5,TRIM(TEXT)) 
       CALL pgband(0,1,xdm,ydm,xdm,ydm,ps)  
       IF (ICHAR(ps)==27) EXIT
       i_con=i_con+1
       CALL PGSLCT(IHIST)
       CALL PGCLOS
       ENDDO
       CALL PGSLCT(IHIST)
       CALL PGCLOS


       CALL PGSLCT(ISTAT)
       CALL clear_scr(0)
      END SUBROUTINE
