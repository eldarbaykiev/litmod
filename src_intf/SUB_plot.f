!SUB_layer_plot.f


      SUBROUTINE plot
      USE M_imag
      USE M_costa
      USE M_profile
      USE M_material
!      USE M_pol
      USE M_label
      USE M_ref_points
      INTEGER i_up,i_dwn
      REAL(4) CONTRA_pl,BRIGHT_pl
      IF(ALLOCATED(bod_tk)) DEALLOCATE (bod_tk)
      ALLOCATE (bod_tk(IDIM,JDIM))
      IF (i_top>i_bot) THEN
       i_up=i_top
       i_dwn=i_bot
      ELSE
       i_up=i_bot
       i_dwn=i_top
      ENDIF
      I_BODY=i_up
      IF (i_dwn==0) THEN
        bod_tk=layers(:,:,1)+E*1e-3
       ELSE
        bod_tk=layers(:,:,i_up)-layers(:,:,i_dwn)
       ENDIF
      z_min_pl=minval(bod_tk)
      z_max_pl=maxval(bod_tk)
      IF (act=='B') THEN
       RETURN
      ENDIF
       
C Set up the color map. A negative sign in CONTRA reverses
C the color scale
      CALL PGSLCT(IPG_LAY)
      BRIGHT_pl = 0.5
      CONTRA_pl  = -1.
      CALL PALETT(6,CONTRA_pl,BRIGHT_pl) 
      call pgscr(0,1.,1.,1.)
      call pgscr(1,0.,0.,0.)
!      call pgscr(17,.5,.5,.5)
      call pgslw(1)
      CALL PGSCI(1)
!      CALL PGENV(lon_min,lon_max,lat_min,lat_max,1,0)
      call pgsvp(.17,.7,.1,.9)
      call pgswin(lon_min,lon_max,lat_min,lat_max) 
      CALL PGBOX('bctsinm',0.0,0,'bctsivnm',0.0,0)
      CALL PGIMAG(bod_tk,IDIM,JDIM,1,IDIM,1,JDIM,
     *            0,z_max_pl,TR)
      CALL PGLAB('Long ','Lat', TRIM(prop(i_up)%name_mat))

!!      IF (i_pr>=0.OR.act=='B')  THEN
      IF (i_pr>=0)  THEN
      CALL  PGWEDG('LI',4.,5.,0,z_max_pl,'Thickness (km)')
      ENDIF 

      CALL PGIDEN
!Dibujar los puntos de la linea de costa
      DO ii=1,coast_lin
        CALL PGLINE(npt_coast(ii),L_coast(2*ii-1,1:npt_coast(ii)),
     *              L_coast(2*ii,1:npt_coast(ii)))
      ENDDO
!Plot labels
      DO i=1,NMARK
       CALL PGPT1(labels(i)%long,labels(i)%lat,18)
      ENDDO
!Plot reference points
      CALL PGSFS(3)
      DO j=1,N_RPTS
      CALL PGPOLY(N_pol_RPTS(j),pol_RPTS(2*j-1,:),pol_RPTS(2*j,:))
      CALL PGLINE(N_pol_RPTS(j)+1,pol_RPTS(2*j-1,:),pol_RPTS(2*j,:))
      ENDDO
      CALL PGSFS(1)      

      END SUBROUTINE 
