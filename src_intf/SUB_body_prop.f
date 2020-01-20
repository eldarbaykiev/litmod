! SUB_body_prop.f


      SUBROUTINE body_prop
      use M_imag
      use M_material
      INTEGER PGOPEN,N_CH,I_PROP
      REAL(4) xbox(4),ybox(4),BOX_E(4),BOX_R(4),BOX_A(4),xdm,ydm,ydf
      REAL(4) BOX_P(4,10) 
      CHARACTER(1) dmch
      REAL(4) dfts(10)
      CHARACTER(20) STRG(10),dfts_name
      CHARACTER(20) format_txt(10)
      LOGICAL EXIT_MN
       EXIT_MN=.FALSE.
       abort_body='n'
       dfts_name='new_mat'
       dfts(2)=2750
       dfts(3)=0.350D-04
       dfts(4)=2.3
       dfts(5)=1D-7
       dfts(6)=7.69D-12
       dfts(7)=1.25
       dfts(8)=4.3
       dfts(9)=130
       dfts(10)=21
       format_txt(1)='(A)'
       format_txt(2)='(F5.0)'
       format_txt(3)='(E12.3)'
       format_txt(4)='(F3.1)'
       format_txt(5)='(E12.3)'
       format_txt(6)='(E12.3)'
       format_txt(7)='(F4.2)'
       format_txt(8)='(F4.2)'
       format_txt(9)='(F5.1)'
       format_txt(10)='(F3.0)'

       prop(ins_mat)%name_mat=dfts_name
       prop(ins_mat)%rho=dfts(2)
       prop(ins_mat)%alpha=dfts(3)
       prop(ins_mat)%K=dfts(4)
       prop(ins_mat)%A=dfts(5)
       prop(ins_mat)%beta=dfts(6)
       prop(ins_mat)%gmma_T=dfts(7)
       prop(ins_mat)%K_0p=dfts(8)
       prop(ins_mat)%K_T=dfts(9)
       prop(ins_mat)%litho=dfts(10)
       
      IF(LINUX) THEN
         I_PROP=PGOPEN('/XWINDOW')
      ELSE
         I_PROP=PGOPEN('/WX')
      ENDIF
      IF (I_PROP<=0) STOP
      call pgask(.false.)
      CALL PGPAP(5.0,1.318)
      CALL PGSVP(.0,1.,.0,1.)
      CALL PGSWIN(.0,1.,.0,2.)
      CALL PGSCR(0,1.,1.,1.)
      CALL PGSCR(1,0.,0.,0.)
      CALL PGSCH(1.6)
!Exit button
      CALL PGPTEXT(.5,1.95,0.,0.5,'BODY PROPERTIES')
      CALL PGSCH(1.3)
      CALL PGQTXT(.25,.05,0.,0.5,'EXIT',xbox,ybox)
      BOX_E(1)=xbox(1)-.01
      BOX_E(2)=xbox(3)+.01
      BOX_E(3)=ybox(1)-.01
      BOX_E(4)=ybox(2)+.01
      CALL PGSCI(3)
      CALL PGRECT(BOX_E(1),BOX_E(2),BOX_E(3),BOX_E(4))
      CALL PGSCI(1) 
      CALL PGPTEXT(.25,.05,0.,0.5,'EXIT')
!Reset button
      CALL PGQTXT(.5,.05,0.,0.5,'RESET',xbox,ybox)
      BOX_R(1)=xbox(1)-.01
      BOX_R(2)=xbox(3)+.01
      BOX_R(3)=ybox(1)-.01
      BOX_R(4)=ybox(2)+.01
      CALL PGSCI(5)
      CALL PGRECT(BOX_R(1),BOX_R(2),BOX_R(3),BOX_R(4))
      CALL PGSCI(1)
      CALL PGPTEXT(.5,.05,0.,0.5,'RESET')
!Abort button
      CALL PGQTXT(.75,.05,0.,0.5,'ABORT',xbox,ybox)
      BOX_A(1)=xbox(1)-.01
      BOX_A(2)=xbox(3)+.01
      BOX_A(3)=ybox(1)-.01
      BOX_A(4)=ybox(2)+.01
      CALL PGSCI(2)
      CALL PGRECT(BOX_A(1),BOX_A(2),BOX_A(3),BOX_A(4))
      CALL PGSCI(1)
      CALL PGPTEXT(.75,.05,0.,0.5,'ABORT')
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Body properties buttons
      CALL PGQTXT(.5,1.85,0.,0.5,'BODY NAME',xbox,ybox)
      BOX_P(1,1)=.25
      BOX_P(2,1)=.75
      BOX_P(3,1)=ybox(1)-.08
      BOX_P(4,1)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,1),BOX_P(2,1),BOX_P(3,1),BOX_P(4,1))
      CALL PGPTEXT(.5,1.85,0.,0.5,'BODY NAME')

      CALL PGQTXT(.5,1.7,0.,0.5,'DENSITY (kg/m³)',xbox,ybox)
      BOX_P(1,2)=.25
      BOX_P(2,2)=.75
      BOX_P(3,2)=ybox(1)-.08
      BOX_P(4,2)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,2),BOX_P(2,2),BOX_P(3,2),BOX_P(4,2))
      CALL PGPTEXT(.5,1.7,0.,0.5,'DENSITY (kg/m³)')

      CALL PGQTXT(.5,1.55,0.,0.5,'ALPHA (1/k)',xbox,ybox)
      BOX_P(1,3)=.25
      BOX_P(2,3)=.75
      BOX_P(3,3)=ybox(1)-.08
      BOX_P(4,3)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,3),BOX_P(2,3),BOX_P(3,3),BOX_P(4,3))
      CALL PGPTEXT(.5,1.55,0.,0.5,'ALPHA (1/K)')

      CALL PGQTXT(.5,1.4,0.,0.5,'THERMAL CONDUCTIVITY (W/K)',xbox,ybox)
      BOX_P(1,4)=.25
      BOX_P(2,4)=.75
      BOX_P(3,4)=ybox(1)-.08
      BOX_P(4,4)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,4),BOX_P(2,4),BOX_P(3,4),BOX_P(4,4))
      CALL PGPTEXT(.5,1.4,0.,0.5,'THERMAL CONDUCTIVITY (W/K)')

      CALL PGQTXT(.5,1.25,0.,0.5,'HEAT PRODUCTION RATE (W/m³)',
     *            xbox,ybox)
      BOX_P(1,5)=.25
      BOX_P(2,5)=.75
      BOX_P(3,5)=ybox(1)-.08
      BOX_P(4,5)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,5),BOX_P(2,5),BOX_P(3,5),BOX_P(4,5))
      CALL PGPTEXT(.5,1.25,0.,0.5,'HEAT PRODUCTION RATE (W/m³)')

      CALL PGQTXT(.5,1.1,0.,0.5,'PRESSURE COEFFICIENT',xbox,ybox)
      BOX_P(1,6)=.25
      BOX_P(2,6)=.75
      BOX_P(3,6)=ybox(1)-.08
      BOX_P(4,6)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,6),BOX_P(2,6),BOX_P(3,6),BOX_P(4,6))
      CALL PGPTEXT(.5,1.1,0.,0.5,'PRESSURE COEFFICIENT') 

      CALL PGQTXT(.5,.95,0.,0.5,'GRÜNEISEN PARAM',xbox,ybox)
      BOX_P(1,7)=.25
      BOX_P(2,7)=.75
      BOX_P(3,7)=ybox(1)-.08
      BOX_P(4,7)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,7),BOX_P(2,7),BOX_P(3,7),BOX_P(4,7))
      CALL PGPTEXT(.5,.95,0.,0.5,'GRÜNEISEN PARAM')
      
      CALL PGQTXT(.5,.8,0.,0.5,'PRES DER ISOTH BULK MODULII',xbox,ybox)
      BOX_P(1,8)=.25
      BOX_P(2,8)=.75
      BOX_P(3,8)=ybox(1)-.08
      BOX_P(4,8)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,8),BOX_P(2,8),BOX_P(3,8),BOX_P(4,8))
      CALL PGPTEXT(.5,.8,0.,0.5,'PRES DER ISOTH BULK MODULUS')
       
      CALL PGQTXT(.5,.65,0.,0.5,'ISOTH BULK MODULUS',xbox,ybox)
      BOX_P(1,9)=.25
      BOX_P(2,9)=.75
      BOX_P(3,9)=ybox(1)-.08
      BOX_P(4,9)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,9),BOX_P(2,9),BOX_P(3,9),BOX_P(4,9))
      CALL PGPTEXT(.5,.65,0.,0.5,'ISOTH BULK MODULUS')

      CALL PGQTXT(.5,.5,0.,0.5,'MANTLE INDEX',xbox,ybox)
      BOX_P(1,10)=.25
      BOX_P(2,10)=.75
      BOX_P(3,10)=ybox(1)-.08
      BOX_P(4,10)=ybox(2)-.055
      CALL PGSCI(15)
      CALL PGRECT(BOX_P(1,10),BOX_P(2,10),BOX_P(3,10),BOX_P(4,10))
      CALL PGPTEXT(.5,.5,0.,0.5,'MANTLE INDEX')
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Write default values
      DO i =1,10
      ydf=(BOX_P(3,i)+BOX_P(4,i))/2-.01
       IF (i==1) THEN
       write(STRG(i),format_txt(i))dfts_name 
       ELSE
       write(STRG(i),format_txt(i)) dfts(i)
       ENDIF
!       write(*,format_txt(i)) dfts(i)
      CALL PGSCI(1) 
      CALL PGPTEXT(.5,ydf,0.,0.5,TRIM(STRG(i)))
      ENDDO
!Require action
      xdm=.5
      ydm=1.
      DO
      DO i =1,10
!      CALL PGBAND(0,1,xdm,ydm,xdm,ydm,dmch)
!       IF (xdm>=BOX_P(1,i).AND.xdm<=BOX_P(2,i).AND.
!     *    ydm>=BOX_P(3,i).AND.ydm<=BOX_P(4,i)) THEN  

       CALL PGSCI(2)
       CALL PGSFS(2)
       CALL PGRECT(BOX_P(1,i),BOX_P(2,i),BOX_P(3,i),
     *             BOX_P(4,i))
       CALL PGSFS(1)
       CALL PGSCI(1)
       DO

       CALL PGBAND(0,0,xdm,ydm,xdm,ydm,dmch)
!Abort
       IF (xdm>=BOX_A(1).AND.xdm<=BOX_A(2).AND.
     *     ydm>=BOX_A(3).AND.ydm<=BOX_A(4)) THEN

       abort_body='y'
       CALL PGCLOS
       RETURN
       ENDIF
!Exit
       IF (xdm>=BOX_E(1).AND.xdm<=BOX_E(2).AND.
     *     ydm>=BOX_E(3).AND.ydm<=BOX_E(4)) THEN
       EXIT_MN=.TRUE.
       EXIT
       ENDIF
!Reset
      IF (xdm>=BOX_R(1).AND.xdm<=BOX_R(2).AND.
     *    ydm>=BOX_R(3).AND.ydm<=BOX_R(4)) THEN
        DO k =1,10
         IF (k==1) THEN
         STRG(k)=dfts_name
         ELSE
         write(STRG(k),format_txt(k)) dfts(k)
         ENDIF
         CALL PGSCI(15)
         CALL PGSFS(1)
         CALL PGRECT(BOX_P(1,k),BOX_P(2,k),BOX_P(3,k),BOX_P(4,k))
         ydf=(BOX_P(3,k)+BOX_P(4,k))/2-.01
         CALL PGSCI(1)
         CALL PGPTEXT(.5,ydf,0.,0.5,TRIM(STRG(k)))
        ENDDO
       xdm=.5
       ydm=1
       EXIT
      ENDIF
         N_CH=LEN_TRIM(STRG(i))
!         DO
         IF (ICHAR(dmch)==13.OR.ICHAR(dmch)==9) THEN
         CALL PGSCI(1)
         CALL PGSFS(2) 
         CALL PGRECT(BOX_P(1,i),BOX_P(2,i),BOX_P(3,i),
     *               BOX_P(4,i))
         CALL PGSFS(1) 
         EXIT
         ENDIF
!Delete
         IF (ICHAR(dmch)==127.OR.ICHAR(dmch)==8) THEN
          IF (N_CH==0) CYCLE
          STRG(i)(N_CH:N_CH)=''
          N_CH=N_CH-1
         ELSE 
          IF ((ICHAR(dmch)<45.OR.ICHAR(dmch)>57.OR.ICHAR(dmch)==47)
     *        .AND.i>1.AND.ICHAR(dmch)/=100.AND.ICHAR(dmch)/=101
     *        .AND.ICHAR(dmch)/=68.AND.ICHAR(dmch)/=69) CYCLE     
         N_CH=N_CH+1
         STRG(i)(N_CH:N_CH)=dmch
         ENDIF 
         CALL PGSCI(15)
         CALL PGRECT(BOX_P(1,i),BOX_P(2,i),BOX_P(3,i),BOX_P(4,i))
         CALL PGSCI(1)
         ydf=(BOX_P(3,i)+BOX_P(4,i))/2-.01
         CALL PGPTEXT(.5,ydf,0.,0.5,STRG(i)(1:N_CH))
         CALL PGSCI(2)
         CALL PGSFS(2)
         CALL PGRECT(BOX_P(1,i),BOX_P(2,i),BOX_P(3,i),
     *               BOX_P(4,i))
         CALL PGSFS(1)
         ENDDO
        IF (EXIT_MN) EXIT
!       ENDIF
       ENDDO 
       
      IF (EXIT_MN) EXIT 
      ENDDO
       
      DO i =1,10
       IF (i==1) THEN
       prop(ins_mat)%name_mat=STRG(i)
       ELSE
       read(STRG(i),*) dfts(i)
       ENDIF
      ENDDO
      prop(ins_mat)%rho=dfts(2)
      prop(ins_mat)%alpha=dfts(3)
      prop(ins_mat)%K=dfts(4)
      prop(ins_mat)%A=dfts(5)
      prop(ins_mat)%beta=dfts(6)
      prop(ins_mat)%gmma_T=dfts(7)
      prop(ins_mat)%K_0p=dfts(8)
      prop(ins_mat)%K_T=dfts(9)
      prop(ins_mat)%litho=NINT(dfts(10))

      CALL PGSCH(1.)
      CALL PGCLOS 
      END SUBROUTINE  
