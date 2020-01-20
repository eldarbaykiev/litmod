! Subroutine SUB_ch_obs.f

      SUBROUTINE ch_obs
      USE M_imag
      USE M_costa
      USE M_pol
      USE M_label 
      USE M_ref_points
      LOGICAL STAT
      REAL(4) xdm,ydm
      CHARACTER(1) :: le
      CHARACTER(50) UCASE
       
!      IF (fst_call) THEN
!       INQUIRE(FILE='obs_name.dat',EXIST=STAT)
!       IF(.NOT.STAT) THEN
!        WRITE(*,*)'obs_name.dat DOES NOT EXIST' 
!        RETURN
!       ENDIF
!       open(1,file='obs_name.dat',status='OLD')
!Load observable files
!       write(*,*)'LOADING OBSERVABLE FILES'  
!       DO i_obs=2,6 
!        READ(1,*)obs_file(i_obs)
!        write(*,*)obs_file(i_obs)
!        open(20,file=obs_file(i_obs),status='OLD')
!        DO j=JDIM,1,-1
!         DO i=1,IDIM
!         read(20,*) dum,dum,OBS(i_obs,i,j)
!         ENDDO
!        ENDDO
!        close(20)
!       ENDDO
!       close (1) 
!       fst_call=.false.
!      ENDIF 

      CALL PGSLCT(ISTAT)
      i_obs=1
      DO
      z_min=minval(OBS(i_obs,:,:))
      z_max=maxval(OBS(i_obs,:,:))
      CALL pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               -3) 
       xdm=(lon_max+lon_min)/2
       ydm=(lat_max+lat_min)/2
!Erase the title of the plot
       CALL PGSCI(0)
       CALL PGSVP(.27,.9,.93,1.)
       CALL PGSWIN(0.,1.,0.,1.)
       CALL PGRECT(0.,1.,0.,1.)
!Erase the palette
       CALL PGSVP(.9,1.,.1,1.)
       CALL PGSWIN(0.,1.,0.,1.)
       CALL PGRECT(0.,1.,0.,1.)
!Plot Observable
       IMG=OBS(i_obs,:,:)
       IMG_txt=obs_file(i_obs) 
       CALL PALETT(2, CONTRA, BRIGHT)
       CALL clear_scr(0)
!Plot the palette
       CALL PGWEDG('RI',2.5,3.,z_min,z_max,'')

       CALL pgband(0,1,xdm,ydm,xdm,ydm,le)
       IF (ICHAR(le)==13)EXIT
!CHANGE PALETTE parameters*******************************
       IF (ICHAR(le)==9) THEN
!Read min value for the palette
       CALL pgqvp(0,r1,r2,r3,r4)
       CALL text_plot(ISTAT,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               -9)
       DO
!Erase the palette
       CALL PGSCI(0)
       CALL PGSVP(.9,1.,.1,1.)
       CALL PGSWIN(0.,1.,0.,1.)
       CALL PGRECT(0.,1.,0.,1.) 
!Plot Observable
       CALL PALETT(2, CONTRA, BRIGHT)
       CALL clear_scr(0)
!Plot the palette
       CALL PGWEDG('RI',2.5,3.,z_min,z_max,'')
       CALL pgband(0,1,xdm,ydm,xdm,ydm,le)
       IF (ICHAR(le)==13) EXIT
       le=UCASE(le)
        IF(le=='B')THEN
        BRIGHT=BRIGHT+.1
         IF(BRIGHT>1.1) BRIGHT=.0
        ELSEIF(le=='C')THEN
        CONTRA=CONTRA+.1
         IF(CONTRA>3) CONTRA=.0  
        ENDIF
       ENDDO
       ENDIF
!CHANGE PALETTE parameters*******************************
       IF (ICHAR(le)==13) EXIT
       le=UCASE(le)
       IF (ICHAR(le)==82) THEN
        i_obs=i_obs+1
        IF(i_obs>6) i_obs=1
       ELSEIF (ICHAR(le)==67) THEN
        i_obs=i_obs+1
        IF(i_obs<7.OR.i_obs>12) i_obs=7
       ENDIF
      ENDDO
      


      END SUBROUTINE
