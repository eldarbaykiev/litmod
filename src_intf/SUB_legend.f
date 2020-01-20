! SUB_legend.f


      SUBROUTINE legend
      USE M_imag
      USE M_profile
      USE M_material 
      INTEGER PGOPEN
      REAL(4) xbox(4),ybox(4),BOX_E(4),vxdm(2),vydm(2)
      CHARACTER(2) chdm
      CHARACTER(LEN=*), PARAMETER::string_l='LAYER'
      CHARACTER(50) :: fil_l 
      IF(LINUX) THEN
         I_LEG=PGOPEN('/XWINDOW')
      ELSE
         I_LEG=PGOPEN('/WX')
      ENDIF
      IF (I_LEG<=0) STOP
      call pgask(.false.)
      CALL PGPAP(5.0,1.318)
      CALL PGSVP(.0,1.,.0,1.)
      CALL PGSWIN(.0,1.,.0,2.)
      CALL PGSCR(0,1.,1.,1.)
      CALL PGSCR(1,0.,0.,0.)
      CALL PGSCH(1.6)
     
       CALL PGPTEXT(.5,1.95,0.,0.5,'LAYERS INFO')
!Water layer
!       CALL PGSCI(50)
!       CALL PGQTXT(.5,1.85,0.,0.5,'',xbox,ybox)
       open(16,file='color.rgb',status='UNKNOWN')
       CALL PGQTXT(.5,1.95,0.,0.5,'LAYER',xbox,ybox) 
       write(*,*) xbox,ybox
       xbox(1)=xbox(1)-.04
       xbox(2)=xbox(2)-.04
       xbox(3)=xbox(3)+.04  
       xbox(4)=xbox(4)+.04
       xbox=xbox-.3
       BOX_E(1)=xbox(1)
       BOX_E(2)=xbox(3)

       ydmm1=ybox(1)
       ydmm2=ybox(2)
       ydmm2=ydmm2+.05 
       DO i=1,N_lay
        IF (i>9) THEN
         WRITE(chdm,'(I2)')i
        ELSE
         WRITE(chdm,'(I1)')i
        ENDIF
        fil_l=TRIM(string_l)//TRIM(chdm)
        ybox(1)=ydmm1-.06-i*(1.8/N_lay)
        ybox(2)=ydmm1+.06-i*(1.8/N_lay)
        ybox(3)=ydmm1+.06-i*(1.8/N_lay)   
        ybox(4)=ydmm1-.06-i*(1.8/N_lay)
        BOX_E(3)=ybox(1)
        BOX_E(4)=ybox(2)
        read(16,*) R,G,B
        CALL PGSCR (15+i,R,G,B)
        CALL PGSCI (15+i)
        CALL PGRECT(BOX_E(1),BOX_E(2),BOX_E(3),BOX_E(4))
        CALL PGSCI(1)
        CALL PGLINE(4,xbox,ybox)
        CALL PGPTEXT((BOX_E(1)+BOX_E(2))*.5,(BOX_E(3)+BOX_E(4))*.5,
     *               0.,0.5,fil_l)
        CALL PGPTEXT(.5,(BOX_E(3)+BOX_E(4))*.5,0.,0.,
     *               TRIM(prop(i)%name_mat))
        vxdm(1)=xbox(1)
        vxdm(2)=xbox(4)
        vydm(1)=ybox(1)
        vydm(2)=ybox(4)
        IF(i==i_lay) THEN
         CALL PGSLW(7)
         CALL PGSLS(4)
        ENDIF
        CALL PGLINE(2,vxdm,vydm)
        CALL PGSLW(3)
        CALL PGSLS(1)
       ENDDO

       
       CLOSE(16)
      END SUBROUTINE 
