! SUB_clear_scr.f

      SUBROUTINE clear_scr  (ik)
      USE M_imag
      USE M_costa
      USE M_pol
      USE M_label 
      USE M_ref_points
      
!!       CALL PALETT(2, CONTRA, BRIGHT)
!Establece el color del fondo (0) y de los ejes, etiquetas etc. (1)
!Utiliza escala Red Green Blue normalizada [0.,1.]
      CALL pgscr(0,1.,1.,1.)
      CALL pgscr(1,0.,0.,0.)
!      CALL pgscr(17,.5,.5,.5)
!pgask controls new page prompting
      CALL pgask(.false.)
!      CALL PGPAP(8.0,.8)
      CALL pgslw(1)
      CALL PGSCI(1)
      CALL pgsvp(.3,.9,.1,.9)
!      CALL PGENV(lon_min,lon_max,lat_min,lat_max,1,0)
      CALL pgswin(lon_min,lon_max,lat_min,lat_max)
      CALL PGBOX('bctsinm',0.0,0,'bctsivnm',0.0,0)
!!      IF (act=='B'.AND.out_win) THEN
!!      CALL PGIMAG(bod_tk(:,:),IDIM,JDIM,1,IDIM,1,JDIM,z_min_pl,
!!     *            z_max_pl,TR)
!!      CALL PGLAB('Long ', 'Lat',TRIM(prop(I_BODY)%name_mat))
!!      ELSE
!!      CALL PGIMAG(OBS(i_obs,:,:),IDIM,JDIM,1,IDIM,1,JDIM,z_min,z_max,TR)
      CALL PGIMAG(IMG(:,:),IDIM,JDIM,1,IDIM,1,JDIM,z_min,z_max,TR)
      CALL PGLAB('Long ', 'Lat',TRIM(IMG_txt))
           CALL pgqvp(0,r1,r2,r3,r4)
!!      ENDIF
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
      IF (ik==0) RETURN

      DO j=1,ik
      CALL PGPOLY(NP_POL(j),A_pol(2*j-1,:),A_pol(2*j,:))
      CALL PGLINE(NP_POL(j)+1,A_pol(2*j-1,:),A_pol(2*j,:))
      ENDDO

       END SUBROUTINE
