!Subroutine SUB_new_coord.f

 

      SUBROUTINE new_coord 
      USE M_imag
      USE M_profile
      character(1) le,ok,UCASE
      character(12) TEXT 
      real(4) xy,z,xy_2p(2),z_2p(2),xy_vec(0:100),z_vec(0:100)
      REAL(4) xy_vec_dm(0:100),z_vec_dm(0:100)
      real(4) slp,L_fr,z_vec_ini,d_click,d_m
      real(4) xms,yms
      real(4) ms_min,ms_max,ms 
      INTEGER NPT,I_fr_l,I_fr_r,I_L_fr

      exit_mdreg='n'
      abort_mdreg='n'
      ch_v_scl=.FALSE.
      NPT=0
      xy=(xy_pr_max+xy_pr_min)/2
      z=-70
      le='a'
      ok='a' 
      IF (ort=='E') THEN
       d_m=p_m_x
      ELSE
       d_m=p_m_y
      ENDIF

      DO
!//////////First node/////////////////////////
       IF (NPT==0) THEN
        MODE=7
        call pgband(MODE,1,xy,z,xy,z,le)
        LE=UCASE(LE)
        IF (le=='S'.OR.le=='X'.OR.le=='M'.OR.le=='N') EXIT
        IF (le=='C') THEN
         WRITE(*,*)'COPY LAST PROFILE'
         call pgband(0,0,dm,dm,dm,dm,ok)
         EXIT
        ENDIF
        
        IF (le=='Z') THEN
         WRITE(*,*)'CHANGING VERTICAL SCALE'
!!         call pgband(0,0,dm,dm,dm,dm,ok)
         ch_v_scl=.TRUE.      
         RETURN
        ENDIF
        xy=xy_pr_min+nint((xy-xy_pr_min)/d_m)*d_m  
        CALL PGPT1(xy,z,18)
        xy_vec(NPT)=xy
        z_vec(NPT)=z
        NPT=NPT+1  
        cycle
!//////////First node/////////////////////////
       ELSE
        MODE=7
       ENDIF 

       call pgband(MODE,1,xy,z,xy,z,le)
       LE=UCASE(LE)
       IF (le=='S'.OR.le=='X'.OR.le=='M'.OR.le=='N') EXIT
       IF (le=='C') THEN
        call pgband(0,0,dm,dm,dm,dm,ok) 
        EXIT
       ENDIF
!Check for confirmation
       call pgband(0,0,dm,dm,dm,dm,ok)
       OK=UCASE(ok)
       IF (ok=='S'.OR.ok=='X'.OR.ok=='M'.OR.ok=='N') EXIT

        IF (ok=='D') THEN
        xy=xy_vec(NPT-1)
        z=z_vec(NPT-1)
        CYCLE
        ENDIF
!Check if clicks are of shorter wavelength than the grid step
        d_click=ABS((xy-xy_vec(NPT-1))/d_m)
!      write(*,*)'d_click',d_click,z
        IF (d_click<1) CYCLE
        IF (xy>xy_vec(NPT-1)) THEN
         xy=xy_vec(NPT-1)+NINT(d_click)*d_m 
        ELSE
         xy=xy_vec(NPT-1)-NINT(d_click)*d_m        
        ENDIF
                 
 
       CALL PGSCI(1)
       CALL PGSLW(5)

       CALL PGPT1(xy,z,18)
       xy_vec(NPT)=xy
       z_vec(NPT)=z
       xy_2p=xy_vec(NPT-1:NPT)
       z_2p=z_vec(NPT-1:NPT)
       CALL PGLINE(2,xy_2p,z_2p)
       NPT=NPT+1
      ENDDO  


      NPT=NPT-1
      OK=UCASE(OK)
      IF (le=='S'.OR.le=='X') ok=''
 
      IF (ok=='S'.OR.ok=='X'.OR.le=='S'.OR.le=='X') fw_pr=1 
       
      IF (ok=='S'.OR.ok=='M'.OR.le=='S'.OR.le=='M') THEN
!Reorder, if necessary
        IF (xy_vec(1)-xy_vec(0)<0) THEN
         WRITE(*,*)'RE-ORDER PROFILE' 
         xy_vec_dm(:)=xy_vec(:)
         z_vec_dm(:)=z_vec(:)
         DO ik=0,NPT
          xy_vec(ik)=xy_vec_dm(NPT-ik)
          z_vec(ik)=z_vec_dm(NPT-ik)
         ENDDO
        ENDIF    
  
        DO k=0,NPT-1
         
         L_fr=xy_vec(k+1)-xy_vec(k)
         IF (ABS(z_vec(k+1)-z_vec(k))<.3.AND.i_lay>0)
     *       z_vec(k+1)=z_vec(k)
         I_fr_l=nint((xy_vec(k)-xy_pr_min)/d_m)+1
         I_fr_r=nint((xy_vec(k+1)-xy_pr_min)/d_m)+1
         IF (I_fr_r<I_fr_l) THEN
          I_fr_dm=I_fr_l
          I_fr_l=I_fr_r
          I_fr_r=I_fr_dm
         ENDIF
         IF (I_fr_r>ifl) I_fr_r=ifl
         IF (I_fr_l<1) I_fr_l=1
         I_L_fr=I_fr_r-I_fr_l
         slp=(z_vec(k+1)-z_vec(k))/I_L_fr 
         z_vec_ini=z_vec(k)

         IF (I_fr_l==I_fr_r-1) THEN
!Click of length equal to 1 node
          bod_pr(I_fr_l,i_lay)=-z_vec(k)
          bod_pr(I_fr_r,i_lay)=-z_vec(k+1)
         ELSE
!Click of length greater than 1 node
          IF (I_fr_l==1) bod_pr(1,i_lay)=-z_vec(0)
          DO j=I_fr_l,I_fr_r-1
           bod_pr(j+1,i_lay)=-1*(slp*(j+1-I_fr_l)+z_vec_ini)
          ENDDO
         ENDIF

        ENDDO
       IF (le=='C') THEN
        bod_pr(:,i_lay)=bod_pr_bak
       ENDIF 
       bod_pr_bak=bod_pr(:,i_lay)
       bod_sf(:,I_p,i_lay)=bod_pr(:,i_lay) 
       pr_pos(i_pr)=I_p
      ENDIF

      IF (ok=='M'.OR.ok=='N'.OR.le=='M'.OR.le=='N') THEN
      CALL PGSLCT(IPG_LAY)
      call pgqvp(0,r1,r2,r3,r4)
      CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *               lat_max,8)
!      write(*,*)'MOUSE INPUT' 
!      write(*,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' 
       IF (ort=='E') THEN
!       write(*,*)'MOVING EAST:' 
!       write(*,*)'  LEFT CLICK--->INTERPOLATE'
!       write(*,*)'  RIGHT CLICK-------->SAVE '
!       write(*,*)'MOVING WEST:'
!       write(*,*)'  ERASE PROFILES EAST OF CURSOR'
!       write(*,*)'ENTER: EXIT' 
!       write(*,*)'Esc:  ABORT'
       ms_min=ltz_min
       ms_max=ltz_max
       ELSE
!       write(*,*)'MOVING NORTH:'
!       write(*,*)'  LEFT CLICK: INTERPOLATE'
!       write(*,*)'  RIGHT CLICK: SAVE '
!       write(*,*)'MOVING SOUTH:'
!       write(*,*)'  ERASE PROFILES NORTH OF CURSOR'
!       write(*,*)'ENTER: EXIT'
!       write(*,*)'Esc:  ABORT'
       ms_min=lnz_min
       ms_max=lnz_max
       ENDIF
      write(*,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      DO 
       call pgband(0,1,0.,0.,xms,yms,lems)
       LEMS=UCASE(LEMS)
       IF(LEMS .EQ. 'S') LEMS='X'
       IF (ort=='E') THEN
          ms=yms 
       ELSE
          ms=xms
       ENDIF 

       IF (ICHAR(lems)==27) THEN
          CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *               lat_max,9) 
!       write(*,*)'REALLY ABORT ? (y/n)' 
!       read(*,*) abort_mdreg
          IF (abort_mdreg=='y') EXIT
          CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *               lat_max,8)
          CYCLE
       ENDIF
       IF (ICHAR(lems)==13) THEN
          CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *                 lat_max,10)
!        write(*,*)'REALLY EXIT ? (y/n)'
!        read(*,*) exit_mdreg
          IF (exit_mdreg=='y') EXIT
          CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *               lat_max,8)
          CYCLE
       ENDIF
!Clicking outside the map nothing happens
!       IF (xms<lon_min.OR.xms>lon_max.OR.yms<lat_min.OR.yms>lat_max)
!     *     CYCLE

       IF (lems=='A') THEN
          TEXT='INTERPOLATE'
       ELSE
          TEXT='SAVE'
       ENDIF

       ms_pr=nint ((ms-xy_pr_p_min)*ifl_p/L_p)+1
       fw_pr=ms_pr-i_pr 
       write(*,*)TEXT,'  INTERMEDIATE PROFILES : '
       write(*,9) i_pr+1,ms_pr-1
       IF (ms>xy_pr_p(i_pr)) THEN
          IF (ms>=ms_max) fw_pr=N_pr-i_pr
          EXIT 
       ELSE
          IF (ms<=ms_min) fw_pr=-i_pr+1
          CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *               lat_max,11)        
!          write(*,*)'REALLY ERASE MODIFICATIONS ? (y/n)'
!          read(*,*) ms_flow
          IF (ms_flow=='y') EXIT
          CALL text_plot(IPG_LAY,r1,r2,r3,r4,lon_min,lon_max,lat_min,
     *               lat_max,8)
       ENDIF

      ENDDO

      ENDIF
  9   FORMAT ('(',I3,',',I3,')')
      END SUBROUTINE
