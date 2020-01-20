! SUB_text_plot.f

      SUBROUTINE text_plot(IPLOT,xvp_min,xvp_max,yvp_min,yvp_max,
     *                     xw_min,xw_max,yw_min,yw_max,ICALL )

      USE M_imag
      USE M_profile
      USE M_material
      USE M_pol
      USE M_label
      USE M_ref_points
      USE M_visualize, ONLY:a_v,att_y
      USE M_sublitho_an, ONLY:title_slc_an,D_T_an,litho_an 
      INTEGER IPLOT,ICALL
      REAL(4) xvp_min,xvp_max,yvp_min,yvp_max
      REAL(4) xw_min,xw_max,yw_min,yw_max,p1,p2 
      CHARACTER(80) TEXT 
      CHARACTER(1) dmch
      CHARACTER(20) str_j
      REAL(4) xbox(2),zbox(2)
      REAL(4) xdm,ydm,zdm
      LOGICAL STAT,STAT1
      CHARACTER(50) UCASE 

      CALL PGSLCT(IPLOT)  
      CALL PGSCI(0)
      CALL PGSFS(1)
      CALL PGSLW(1)
!ICALL<6-->print text at the left of the screen
      IF (ICALL<6) THEN
         CALL PGSCH(.7)
         CALL PGSVP(.05,.29,.1,.9)
         CALL PGSWIN(0.,1.,0.,1.)
         CALL PGRECT(0.,1.,0.,.94)
         CALL PGSCI(1)
      ELSEIF (ICALL>=7) THEN
         CALL PGSVP(.7,1.,.1,.9)
         CALL PGSWIN(0.,1.,0.,1.)
         CALL PGRECT(0.15,1.,0.,1.)
         CALL PGSCI(1)
         CALL PGSCH(.7)
      ENDIF

      IF (ICALL==1) THEN
         write(TEXT,12) N_lay
 12      FORMAT ('Nº OF BODIES DEFINED--> ',I2)
         CALL PGPTEXT(.01,1.,0.,0.0,TRIM(TEXT))
         CALL PGPTEXT(.01,.9,0.,0.0,'EXTRACT PROFILE: p ')
!         CALL PGPTEXT(.01,.85,0.,0.0,'CROSS SECTION: k ')
         CALL PGPTEXT(.01,.85,0.,0.0,'MODIFY REGION: m ')
         CALL PGPTEXT(.01,.8,0.,0.0,'ADD NEW LAYER: a ')
         CALL PGPTEXT(.01,.75,0.,0.0,'JOIN TWO LAYERS: j ')
         CALL PGPTEXT(.01,.7,0.,0.0,'BODY PLOT: b ')
         CALL PGPTEXT(.01,.65,0.,0.0,'BODY PLOT (PS): s')
         CALL PGPTEXT(.01,.6,0.,0.0,'CHANGE OBS: c')
         CALL PGPTEXT(.01,.55,0.,0.0,'LABELS: l')
         CALL PGPTEXT(.01,.5,0.,0.0,'REFERENCE POINTS: r')
!         CALL PGPTEXT(.01,.4,0.,0.0,'OCEANIC/CONT CRUST: t')
         CALL PGPTEXT(.01,.45,0.,0.0,'DELETE LAYER: d')
         CALL PGPTEXT(.01,.4,0.,0.0,'SMOOTH LAY: i')
         CALL PGPTEXT(.01,.35,0.,0.0,'LAYER HIERARCHY: h')
         CALL PGPTEXT(.01,.3,0.,0.0,'VISUALIZE SLICE: v')
         CALL PGPTEXT(.01,.25,0.,0.0,'Z-PROFILE: z')
         CALL PGPTEXT(.01,.2,0.,0.0,'ADD SUBLITH BODY: u')
         CALL PGPTEXT(.01,.15,0.,0.0,'HISTOGRAM RESIDUAL: o')
         CALL PGPTEXT(.01,.10,0.,0.0,'HEAT FLOW BASEMENT: f') 
         CALL PGPTEXT(.01,.05,0.,0.0,'EXIT: e')
         CALL PGPTEXT(.01,.03,0.,0.0,'************************')
         CALL PGPTEXT(.01,.01,0.,0.0,'SELECT OPTION ...'   )
         call pgband(0,0,xdm,ydm,xdm,ydm,act)
      ELSEIF (ICALL==2) THEN
         CALL PGPTEXT(.01,.9,0.,0.0,'PS OUTPUT: p (DEF=NO)')
         CALL pgband(0,0,xdm,ydm,xdm,ydm,prf_act)
         prf_act=UCASE(prf_act)
         CALL PGPTEXT(.01,.85,0.,0.0,'ENTER PROFILE ') 
         CALL PGPTEXT(.01,.83,0.,0.0,'COORDINATES...')
         IF(LINUX)THEN
            CALL PGPTEXT(.01,.8,0.,0.0,'OR')
            CALL PGPTEXT(.01,.7,0.,0.0,'PRESS Esc TO EXIT')
         ELSE
            CALL PGPTEXT(0.01,.85,0.,0.0,'Place cursor on point')
            CALL PGPTEXT(0.01,.8,0.,0.0,'Press "a" to accept point')
            CALL PGPTEXT(.01,.75,0.,0.0,'OR')
            CALL PGPTEXT(.01,.7,0.,0.0,'PRESS Esc TO EXIT')
         ENDIF

      ELSEIF (ICALL==3) THEN
         CALL PGPTEXT(.01,.9,0.,0.0,'PUT LAYER i AT DEPTH OF LAYER j')
         CALL PGPTEXT(.01,.85,0.,0.0,'IF i<j --> CREATE  BODY i ')
         CALL PGPTEXT(.01,.8,0.,0.0,'IF i>j --> REMOVE  BODY i ')
         IF (N_lay>9) THEN
          CALL PGPTEXT(.01,.77,0.,0.0,'N_lay>9, USE 2 DIGIT FORMAT ')
          CALL PGPTEXT(.01,.745,0.,0.0,'(i.e. LAYER 1-->01)')
         ENDIF
         CALL PGPTEXT(.01,.7,0.,0.0,'LAYER i ...? ')
         str_j=''
         DO i=1,2
          CALL pgband(0,1,0.,0.,xdm,ydm,dmch)
          dmch=UCASE(dmch)
          IF (ICHAR(dmch)>47.AND.ICHAR(dmch)<58) THEN
           str_j(i:i)=dmch
           IF (N_lay<10) EXIT
          ELSE
           WRITE(*,*)'INVALID NUMBER, ABORTING'
           ins=-1000
           RETURN 
          ENDIF
         ENDDO
         READ(str_j,*) ins
         CALL PGPTEXT(.5,.7,0.,0.0,TRIM(str_j))
         IF (ins>N_lay) THEN
          WRITE(*,*)'LAYER DOES NOT EXIST (i>N_lay), ABORTING'
          ins=-1000
          RETURN
         ENDIF

         str_j=''
         CALL PGPTEXT(.01,.65,0.,0.0,'LAYER j ...? ')
         DO i=1,2 
          CALL pgband(0,1,0.,0.,xdm,ydm,dmch)
          dmch=UCASE(dmch)
          IF (ICHAR(dmch)>47.AND.ICHAR(dmch)<58) THEN
           str_j(i:i)=dmch
           IF (N_lay<10) EXIT
          ELSE
            WRITE(*,*)'INVALID NUMBER, ABORTING'
            i_top=-1000
            RETURN
          ENDIF
         ENDDO
         READ(str_j,*)i_top
         CALL PGPTEXT(.5,.65,0.,0.0,TRIM(str_j))
         IF (i_top>N_lay) THEN
          WRITE(*,*)'LAYER DOES NOT EXIST (j>N_lay), ABORTING'
          i_top=-1000
          RETURN
         ENDIF

         CALL PGPTEXT(.01,.6,0.,0.0,'JOIN LAYERS INSIDE POLYGON(S)')
         CALL PGPTEXT(.01,.55,0.,0.0,'CHANGE (N=def, Y)? ')
         CALL pgband(0,1,0.,0.,xdm,ydm,dmch)
         IF(ICHAR(dmch)== 13) THEN
          dmch='N'
         ELSE
          dmch=UCASE(dmch)
         ENDIF  
         
         IF (TRIM(dmch)=='Y') THEN
          yin_yang=.FALSE.
         ELSE
          yin_yang=.TRUE.
         ENDIF
      ELSEIF (ICALL==0) THEN

       CALL PGPTEXT(.01,.9,0.,0.0,'CHECK HIERARCHY VIOLATIONS')
       CALL PGPTEXT(.01,.85,0.,0.0,'BETWEEN LAYER i AND j (i<j) ')
       IF (N_lay>9) THEN
        CALL PGPTEXT(.01,.82,0.,0.0,'N_lay>9, USE 2 DIGIT FORMAT ')
        CALL PGPTEXT(.01,.795,0.,0.0,'(i.e. LAYER 1-->01)')
       ENDIF
       str_j=''
       CALL PGPTEXT(.01,.75,0.,0.0,'LAYER i ...? ')
       DO i=1,2 
        CALL pgband(0,1,0.,0.,xdm,ydm,dmch)
        dmch=UCASE(dmch)
        IF (ICHAR(dmch)>47.AND.ICHAR(dmch)<58) THEN
         str_j(i:i)=dmch
         IF (N_lay<10) EXIT
!         READ(dmch,'(I3)') i_up_h
!        C ALL PGPTEXT(.5,.75,0.,0.0,TRIM(dmch))
        ELSE
         WRITE(*,*)'INVALID NUMBER, ABORTING'
         i_up_h=-1000
         RETURN
        ENDIF
       ENDDO
       READ(str_j,*)i_up_h
       CALL PGPTEXT(.5,.75,0.,0.0,TRIM(str_j))
       IF (i_up_h>N_lay) THEN
        WRITE(*,*)'LAYER DOES NOT EXIST (i>N_lay), ABORTING'
        i_up_h=-1000
        RETURN
       ENDIF

       str_j=''
       CALL PGPTEXT(.01,.7,0.,0.0,'LAYER j ...? ')
       DO i=1,2 
        CALL pgband(0,1,0.,0.,xdm,ydm,dmch)
        dmch=UCASE(dmch)
        IF (ICHAR(dmch)>47.AND.ICHAR(dmch)<58) THEN
!         READ(dmch,'(I3)') i_dw_h
!         CALL PGPTEXT(.5,.7,0.,0.0,TRIM(dmch))
         str_j(i:i)=dmch
         IF (N_lay<10) EXIT
        ELSE
         WRITE(*,*)'INVALID NUMBER, ABORTING'
         i_dw_h=-1000
         RETURN
        ENDIF
       ENDDO
       READ(str_j,*)i_dw_h
       CALL PGPTEXT(.5,.7,0.,0.0,TRIM(str_j))
       IF (i_dw_h>N_lay) THEN
        WRITE(*,*)'LAYER DOES NOT EXIST (j>N_lay), ABORTING'
        i_dw_h=-1000
        RETURN
       ENDIF

       CALL PGPTEXT(.01,.65,0.,0.0,'GLUE OPTION OFF,')
       CALL PGPTEXT(.01,.6,0.,0.0,'CHANGE (N=def, Y)? ')
       CALL pgband(0,1,0.,0.,xdm,ydm,dmch)
       IF(ICHAR(dmch)== 13) THEN
        dmch='N'
       ELSE
        dmch=UCASE(dmch)
       ENDIF
       IF (TRIM(dmch)=='Y') THEN
        glue_l=.TRUE.
        CALL PGPTEXT(.01,.55,0.,0.0,'GLUE ON,')
       ELSE
        glue_l=.FALSE.
        CALL PGPTEXT(.01,.55,0.,0.0,'GLUE OFF,') 
       ENDIF

       CALL PGPTEXT(.01,.5,0.,0.0,'FOR HIERARCHY VIOLATIONS ')
       CALL PGPTEXT(.01,.45,0.,0.0,'SET BOTH LAYERS TO THE DEPTH')
       CALL PGPTEXT(.01,.4,0.,0.0,'OF LAYER i (U) OR j (L) ...? ') 
       CALL pgband(0,1,0.,0.,xdm,ydm,hi_ch) 
       hi_ch=UCASE(hi_ch)

      ELSEIF (ICALL==-2) THEN

       INQUIRE(FILE='labels.dat',EXIST=STAT)
       IF(STAT) THEN
        CALL PGPTEXT(.0,.77,0.,0.0,'labels.dat EXIST')
        CALL PGPTEXT(.0,.67,0.,0.0,'ADD NEW MARKS (a)')
        CALL PGPTEXT(.0,.57,0.,0.0,'OR REWRITE ALL (r)? (a=def) ')
        CALL pgband(0,0,xdm,ydm,xdm,ydm,lbl_act) 
        IF(ICHAR(lbl_act)==13) THEN
            lbl_act='A'
         ELSE
            lbl_act=UCASE(lbl_act)
        ENDIF
        IF (lbl_act=='R') CLOSE(911)
      ELSE
       CALL PGPTEXT(.0,.77,0.,0.0,'labels.dat DOES NOT EXIST')
       lbl_act='R'
      ENDIF

      ELSEIF (ICALL==-7) THEN
!Call from SUB_ref_points.f 
       INQUIRE(FILE='ref_points.dat',EXIST=STAT)
       INQUIRE(FILE='ref_points.info',EXIST=STAT1)
       STAT=STAT.AND.STAT1
      IF(STAT) THEN
        CALL PGPTEXT(.0,.77,0.,0.0,'Ref points files EXISTS')
        CALL PGPTEXT(.0,.67,0.,0.0,'ADD NEW REF POINTS (a)')
        CALL PGPTEXT(.0,.57,0.,0.0,'OR REWRITE ALL (r)? (a=def) ')
        CALL pgband(0,0,xdm,ydm,xdm,ydm,RPTS_act)
        IF(ICHAR(RPTS_act)==13) THEN
            RPTS_act='A'
        ELSE
            RPTS_act=UCASE(RPTS_act)
        ENDIF
        IF (RPTS_act=='R') THEN
         CLOSE(40)
         CLOSE(45)
        ENDIF
      ELSE
        CALL PGPTEXT(.0,.77,0.,0.0,'Ref points files DO NOT EXIST')
        RPTS_act='R'
      ENDIF
      ELSEIF (ICALL==-3) THEN
         IF(LINUX.AND.act=='A'.OR.act=='M'.OR.act=='B'.OR.act=='I'.OR.
     *      act=='D') THEN
            CALL PGPTEXT(.01,.9,0.,0.0,'CLICK ON THE MAP') 
            CALL PGPTEXT(.01,.85,0.,0.0,'TO GET A COLUMN')
         ELSEIF(act=='A') THEN
            CALL PGPTEXT(.01,.9,0.,0.0,'Get column for choosing layer:')
            CALL PGPTEXT(.01,.85,0.,0.0,'    Place cursor on map')
            CALL PGPTEXT(.01,.8,0.,0.0,'    Press "a" to validate')
         ELSEIF(LINUX.AND.act=='L') THEN      
            CALL PGPTEXT(.01,.9,0.,0.0,'SELECT POINT FOR LABELING')
            CALL PGPTEXT(.01,.8,0.,0.0,'ENTER: EXIT')
            CALL PGPTEXT(.01,.7,0.,0.0,'ESC: ABORT')
         ELSEIF(LINUX.AND.act=='C') THEN
            CALL PGPTEXT(.01,.9,0.,0.0,'R--> RES. OBSERVABLES')
            CALL PGPTEXT(.01,.85,0.,0.0,'C--> CALC. OBSERVABLES') 
            CALL PGPTEXT(.01,.8,0.,0.0,'ENTER: EXIT')
            CALL PGPTEXT(.01,.75,0.,0.0,'TAB: CHANGE PALETTE')
         ELSEIF(LINUX.AND.act=='Z') THEN
            CALL PGPTEXT(.01,.9,0.,0.0,'CLICK ON THE MAP')
            CALL PGPTEXT(.01,.85,0.,0.0,'TO GET A Z-PROFILE')
         ELSEIF(LINUX.AND.act=='V') THEN
            CALL PGPTEXT(.01,.9,0.,0.0,'CLICK ON THE MAP')
            CALL PGPTEXT(.01,.85,0.,0.0,'TO GET THE DEPTH OF THE SLICE')
         ELSEIF(LINUX.AND.act=='U') THEN
            CALL PGPTEXT(.01,.9,0.,0.0,'CLICK ON THE MAP')
            CALL PGPTEXT(.01,.85,0.,0.0,'TO GET THE DEPTH RANGE OF')
            CALL PGPTEXT(.01,.8,0.,0.0,'THE ANOMALOUS SUBLITH. BODIES')
         ELSEIF(LINUX.AND.act=='F') THEN 
           CALL PGPTEXT(.01,.9,0.,0.0,'CLICK ON THE MAP')
           CALL PGPTEXT(.01,.87,0.,0.0,'TO SELECT LAYER HORIZON FOR')
           CALL PGPTEXT(.01,.85,0.,0.0,'BASAL HEAT FLOW')   
         ENDIF
         IF (act=='A'.OR.act=='T'.OR.act=='D') THEN
            CALL PGSVP(.0,.26,.9,1.)
            CALL PGSCI(0)
            CALL PGRECT(0.,1.,0.,1.)
            CALL PGSCI(1)
            CALL PGSVP(.05,.23,.1,.9) 
         ELSEIF (act=='C') THEN
!            CALL PGSVP(.0,1.,.9,1.)
!            CALL PGSCI(0)
!            CALL PGRECT(0.,1.,0.,1.)
!            CALL PGSCI(1)
!            CALL PGSVP(.05,.23,.1,.9)
         ENDIF
      ELSEIF (ICALL==-11) THEN
        CALL PGPTEXT(.01,.9,0.,0.0,'CLICK ON THE MAP')
        CALL PGPTEXT(.01,.85,0.,0.0,'TO GET A SLICE')
        CALL PGPTEXT(.01,.8,0.,0.0,'LEFT CLICK: VERTICAL PROFILE')
!        CALL PGPTEXT(.01,.8,0.,0.0,'LEFT CLICK: E-W PROFILE')
!        CALL PGPTEXT(.01,.75,0.,0.0,'RIGHT CLICK: N-S PROFILE')
        CALL PGPTEXT(.01,.75,0.,0.0,'OTHERWISE: Z-SLICE')     
      ELSEIF (ICALL==-12) THEN
        CALL PGSCH(.67) 
        CALL PGPTEXT(.01,.9,0.,0.0,'PLOT SEISMIC VELOCITIES: v')    
        CALL PGPTEXT(.01,.87,0.,0.0,'PLOT VELOCITY ANOMALIES: a (=DEF)')
        CALL pgband(0,1,xdm,ydm,xdm,ydm,a_v)
        a_v=UCASE(a_v)
        CALL PGPTEXT(.01,.83,0.,0.0,'ATTENUATION: y/n (DEF=y) ?')
        CALL pgband(0,1,xdm,ydm,xdm,ydm,att_y)
        att_y=UCASE(att_y)
        CALL PGPTEXT(.01,.78,0.,0.0,'PS OUTPUT: y/n (DEF=n) ?')
      ELSEIF (ICALL==-13) THEN 
        CALL PGPTEXT(.01,.9,0.,0.0,'CYCLE BETWEEN PLOTS: ANY KEY')
        CALL PGPTEXT(.01,.85,0.,0.0,'EXIT: ESC')  
      ELSEIF (ICALL==-14) THEN
         CALL PGSCI(0)
         CALL PGSVP(.05,.26,.88,.93)
         CALL PGSWIN(0.,1.,0.,1.)
         CALL PGRECT(0.,1.,.0,1.)
         CALL PGSVP(.05,.23,.1,.9)
         CALL PGSWIN(0.,1.,0.,1.)
         CALL PGSCI(1)
         CALL PGPTEXT(.01,.9,0.,0.0,'ANOMALOUS TEMP, D_T (K)')
         CALL PGPTEXT(.01,.85,0.,0.0,'FOR SUBLITHOSPHERIC BODY...?')
         CALL lec_char(.83,.79,TEXT) 
         READ(TEXT,*)D_T_an
         WRITE(*,*)'ANOMALOUS TEMP, D_T (K) FOR SUBLITHOSPHERIC BODY',
     *D_T_an
         CALL PGPTEXT(.01,.7,0.,0.0,'LITHO NUMBER (MANT COMPOSITION)')
         CALL PGPTEXT(.01,.65,0.,0.0,'FOR SUBLITHOSPHERIC BODY...?')
         CALL lec_char(.63,.59,TEXT) 
         READ(TEXT,'(I2)')litho_an
         WRITE(*,*)'LITHO NUMBER (MANT COMPOSITION) FOR SUBLITHOSPHERIC
     * BODY',litho_an 
      ELSEIF (ICALL==-17) THEN 
        CALL PGPTEXT(.01,.78,0.,0.0,'PS OUTPUT: y/n (DEF=n) ?') 
!Read D_T and litho
      ELSEIF (ICALL==-18) THEN
       CALL PGPTEXT(.01,.9,0.,0.0,'HEAT FLOW AT THE INTERFACE')
       CALL PGPTEXT(.01,.85,0.,0.0,'BETWEEN LAYERS: ')
       CALL PGPTEXT(.01,.8,0.,0.0,TRIM(prop(i_lay)%name_mat)) 
       CALL PGPTEXT(.01,.77,0.,0.0,TRIM(prop(i_lay+1)%name_mat))
      ELSEIF (ICALL==-19) THEN
!       CALL PGSCH(.67)
!       CALL PGPTEXT(.01,.9,0.,0.0,'BASAL HEAT FLOW')
!       CALL PGPTEXT(.01,.87,0.,0.0,'SELECT LAYER HORIZON')
      ELSEIF (ICALL==4) THEN
         CALL PGPTEXT(.01,.9,0.,0.0,TRIM(TEXT_BDSEL))
         IF (i_lay/=-1) THEN
            CALL PGPTEXT(.01,.8,0.,0.0,'OK ...? (y/n)')
            call pgband(0,0,0.,0.,xdm,ydm,bod_ok)
            bod_ok=UCASE(bod_ok)
         ENDIF
      ELSEIF (ICALL==-16) THEN
       CALL PGPTEXT(.0,.77,0.,0.0,'ENTER: SELECT DEPH')
       CALL PGPTEXT(.0,.67,0.,0.0,'ESC: ABORT')
       CALL pgband(5,0,0.,0.,xdm,zdm,dmch)
       IF (ICHAR(dmch)==27) RETURN
       z_an_up=-zdm*1E3  
       CALL pgband(5,0,0.,0.,xdm,zdm,dmch)
       IF (ICHAR(dmch)==27) RETURN
       z_an_dw=-zdm*1E3
      ELSEIF (ICALL==-15) THEN 
       CALL PGPTEXT(.0,.77,0.,0.0,'ENTER: FINISH LABEL WRITING / EXIT')
       CALL PGPTEXT(.0,.67,0.,0.0,'ESC: ABORT')
       CALL PGSVP(xvp_min,xvp_max,yvp_min,yvp_max)
       CALL PGSWIN (xw_min,xw_max,yw_min,yw_max)
       DO
        CALL pgband(5,0,0.,0.,xdm,zdm,dmch)
        write(*,*)'zdm',zdm
        CALL PGPT1(xdm,zdm,18)
        xbox(1)=xdm+.01
        xbox(2)=xdm+.3
        zbox(1)=zdm-2
        zbox(2)=zdm-10
        IF (ICHAR(dmch)==27) THEN
         NMARK=NTICK-1 
         RETURN
        ELSEIF (ICHAR(dmch)==13) THEN
         EXIT
        ENDIF
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4
         ILET=0
         TEXT=''
         CALL PGSCI(15)
         CALL PGRECT(xbox(1),xbox(2),zbox(1),zbox(2))
         CALL PGSCH(1.)
         DO
          CALL PGBAND(0,0,0.,0.,ydm,ydm,dmch)
          I_LAB=ICHAR(dmch)
          IF (I_LAB==13) THEN
             IF(ILET>0) EXIT
             CYCLE
          ELSEIF (I_LAB==27) THEN
           NMARK=NTICK
           RETURN
          ELSEIF (I_LAB==127.OR.I_LAB==8) THEN
             IF(I_LAB<1) CYCLE
             TEXT(ILET:ILET)=''
             ILET=ILET-1
             IF (ILET==0) THEN
                ILET=1
                TEXT=''
             ENDIF
          ELSE
            ILET=ILET+1
            TEXT(ILET:ILET)=dmch
          ENDIF
          CALL PGSCI(15)
          CALL PGRECT(xbox(1),xbox(2),zbox(1),zbox(2))
          CALL PGSCI(1)
          xdm=xbox(1)
          ydm=((zbox(2)+zbox(1))/2)-.7
          CALL PGPTEXT(xdm,ydm,0.,0.0,TRIM(TEXT))
         ENDDO
!         write(*,*)'NMARK,NTICK',NMARK,NTICK
         labels(NMARK+1)%long=labels(NTICK)%long
         labels(NMARK+1)%lat=labels(NTICK)%lat
         labels(NMARK+1)%z_lbl=zdm
         labels(NMARK+1)%text_mark=TRIM(TEXT) 
         NMARK=NMARK+1
         CALL PGSCH(.7)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
       ENDDO
      ELSEIF (ICALL==-10) THEN
         IF(act=='M') THEN
            IF(LINUX) THEN
               CALL PGPTEXT(.01,.9,0.,0.0,'SELECT LAYER TO MODIFY')
            ELSE
               CALL PGPTEXT(.01,.9,0.,0.0,'Place cursor near')
               CALL PGPTEXT(.01,.85,0.,0.0,'  layer to modify')
               CALL PGPTEXT(.01,.8,0.,0.0,'  Press "a" to accept')
            ENDIF
         ELSEIF (act=='B') THEN
            IF(LINUX) THEN
               CALL PGPTEXT(.01,.9,0.,0.0,'SELECT BODY TO PLOT')
            ELSE
               CALL PGPTEXT(.01,.9,0.,0.0,'Place cursor into body')
               CALL PGPTEXT(.01,.85,0.,0.0,'   Press "a" to accept')
            ENDIF
         ELSEIF (act=='A') THEN 
            IF(LINUX) THEN
               CALL PGPTEXT(.01,.9,0.,0.0,'SELECT THE BASE')
               CALL PGPTEXT(.01,.85,0.,0.0,'OF NEW LAYER')
            ELSE
               CALL PGPTEXT(.01,.9,0.,0.0,'Place cursor near base')
               CALL PGPTEXT(.01,.85,0.,0.0,'          of new layer')
               CALL PGPTEXT(.01,.8,0.,0.0,'   Press "a" to accept')
            ENDIF
         ELSEIF (act=='I') THEN
          IF(LINUX) THEN
               CALL PGPTEXT(.01,.9,0.,0.0,'SELECT LAYER TO SMOOTH')
          ELSE
          ENDIF
        ELSEIF (act=='D') THEN
          IF(LINUX) THEN
               CALL PGPTEXT(.01,.9,0.,0.0,'SELECT LAYER TO DELETE')
          ELSE
          ENDIF 
        ELSEIF (act=='U') THEN
          IF(LINUX) THEN
               CALL PGPTEXT(.01,.9,0.,0.0,'SELECT DEPH RANGE')
               CALL PGPTEXT(.01,.85,0.,0.0,'FOR THE ANOMALOUS')
               CALL PGPTEXT(.01,.8,0.,0.0,'SUBLITHOSPHERIC BODIES')
          ELSE
          ENDIF 
         
         ENDIF
      ELSEIF (ICALL==-4) THEN
         DO i=1,N_lay
          CALL PGSCI(0)
          CALL PGSWIN(0.,1.,0.,1.)
          CALL PGRECT(0.,1.,0.,.94)
          CALL PGSCI(1)
!          p1=.9-(i-1)*.12
!          p2=.85-(i-1)*.12
          CALL PGPTEXT(.01,.9,0.,0.0,TRIM(prop(i)%name_mat))
          CALL PGPTEXT(.01,.85,0.,0.0,'PLOT BODY (y/n)... ?')
          CALL pgband(0,1,0.,0.,xdm,ydm,i_bod_plot(i))
          i_bod_plot(i)=UCASE(i_bod_plot(i))
         ENDDO
      ELSEIF (ICALL==-8) THEN 
         CALL PGPTEXT(.01,.9,0.,0.0,'SELECT OCEANIC/TRANS CRUST: o ')
         CALL PGPTEXT(.01,.85,0.,0.0,'SELECT CONTINENT CRUST: c ')
         CALL pgband(0,1,0.,0.,xdm,ydm,cr_type)
         cr_type=UCASE(cr_type)
      ELSEIF (ICALL==-9) THEN
!!!!!!!
!Read min value for the palette
         CALL PGPTEXT(.01,.9,0.,0.0,'MIN VALUE FOR THE PALETTE')  
         ILET=0
         TEXT=''
         CALL PGSCI(15)
         CALL PGRECT(0.,.8,.53,.56)
         DO
          CALL PGBAND(0,1,0.,0.,xdm,ydm,dmch)
          z_min=ICHAR(dmch)
          IF (z_min==13) THEN
             IF(ILET>0) EXIT
             CYCLE
          ELSEIF (z_min>44.AND.z_min<58.AND.z_min/=47.OR.z_min==68.OR.
     *            z_min==69.OR.z_min==100.OR.z_min==101) THEN
             ILET=ILET+1
             TEXT(ILET:ILET)=dmch
          ELSEIF (z_min==127.OR.z_min==8) THEN
             IF(ILET<1) CYCLE
             TEXT(ILET:ILET)=''
             ILET=ILET-1
             IF (ILET==0) THEN
                ILET=1
                TEXT=''
             ENDIF
          ELSE
             write(*,*)'INVALID NUMBER'
          ENDIF
          CALL PGSCI(15)
          CALL PGRECT(0.,0.8,.53,.56)
          CALL PGSCI(1)
          CALL PGPTEXT(.0,.54,0.,0.0,TRIM(TEXT))
         ENDDO
         CALL PGSCI(1)
         READ(TEXT,*) z_min
!Read max value for the palette
         CALL PGPTEXT(.01,.8,0.,0.0,'MAX VALUE FOR THE PALETTE')
         ILET=0
         TEXT=''
         CALL PGSCI(15)
         CALL PGRECT(0.,.7,.53,.56)
         DO
          CALL PGBAND(0,1,0.,0.,xdm,ydm,dmch)
          z_max=ICHAR(dmch)
          IF (z_max==13) THEN
             IF(ILET>0) EXIT
             CYCLE
          ELSEIF (z_max>44.AND.z_max<58.AND.z_max/=47.OR.z_max==68.OR.
     *            z_max==69.OR.z_max==100.OR.z_max==101) THEN
             ILET=ILET+1
             TEXT(ILET:ILET)=dmch
          ELSEIF (z_max==127.OR.z_max==8) THEN
             IF(ILET<1) CYCLE
             TEXT(ILET:ILET)=''
             ILET=ILET-1
             IF (ILET==0) THEN
                ILET=1
                TEXT=''
             ENDIF
          ELSE
             write(*,*)'INVALID NUMBER'
          ENDIF
          CALL PGSCI(15)
          CALL PGRECT(0.,0.7,.53,.56)
          CALL PGSCI(1)
          CALL PGPTEXT(.0,.54,0.,0.0,TRIM(TEXT))
         ENDDO
         CALL PGSCI(1)
         READ(TEXT,*) z_max
         CALL PGPTEXT(.01,.4,0.,0.0,'b: CHANGE BRIGHT')
         CALL PGPTEXT(.01,.35,0.,0.0,'c: CHANGE CONTRAST')
         CALL PGPTEXT(.01,.3,0.,0.0,'ENTER: EXIT')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      ELSEIF (ICALL==5) THEN
!!         CALL PGPTEXT(.01,.9,0.,0.0,'PUSH UP LAYER: u ') 
!!         CALL PGPTEXT(.01,.85,0.,0.0,'PUSH DOWN LAYER: d ')
!!         CALL pgband(0,1,0.,0.,xdm,ydm,updown_lay)
!!         updown_lay=UCASE(updown_lay)
!!          IF (updown_lay=='U') THEN
!!           CALL PGPTEXT(.01,.80,0.,0.0,'ACT SELECTED:')
!!           CALL PGPTEXT(.01,.77,0.,0.0,'PUSH UP LAYER')
!!          ELSE
!!           CALL PGPTEXT(.01,.80,0.,0.0,'ACT SELECTED:')
!!           CALL PGPTEXT(.01,.77,0.,0.0,'PUSH DOWN LAYER') 
!!          ENDIF
         CALL PGPTEXT(.01,.65,0.,0.0,'NUMBER OF PROFILES : ')
         CALL PGPTEXT(.01,.6,0.,0.0,'(If = 0 ->PROFILES ')
         CALL PGPTEXT(.01,.57,0.,0.0,'ARE SELECTED BY MOUSE)')
         ILET=0
         TEXT=''
         CALL PGSCI(15)
         CALL PGRECT(0.,0.3,.53,.56)        
         DO
          CALL PGBAND(0,1,0.,0.,xdm,ydm,dmch)
          N_pr=ICHAR(dmch)
          IF (N_pr==13) THEN
             IF(ILET>0) EXIT
             CYCLE   
          ELSEIF (N_pr>47.AND.N_pr<58) THEN
             ILET=ILET+1
             TEXT(ILET:ILET)=dmch
          ELSEIF (ICHAR(dmch)==127.OR.ICHAR(dmch)==8) THEN
             IF(ILET<1) CYCLE
             TEXT(ILET:ILET)=''
             ILET=ILET-1
             IF (ILET==0) THEN
                ILET=1
                TEXT=''
             ENDIF
          ELSE
             write(*,*)'INVALID NUMBER OF PROFILES'
          ENDIF
          CALL PGSCI(15)
          CALL PGRECT(0.,0.3,.53,.56)
          CALL PGSCI(1)
          CALL PGPTEXT(.0,.54,0.,0.0,TRIM(TEXT))
         ENDDO
         CALL PGSCI(1)
         READ(TEXT,'(I3)') N_pr
         CALL PGPTEXT(.01,.45,0.,0.0,'East-West----> E')
         CALL PGPTEXT(.01,.4,0.,0.0,'North-South--> N')
         CALL pgband(0,1,0.,0.,xdm,ydm,ort)
         ort=UCASE(ort)
         CALL PGPTEXT(.01,.35,0.,0.0,'GLUE: ON, CHANGE? (y/n, def=n)')
         CALL pgband(0,1,0.,0.,xdm,ydm,dmch)
         dmch=UCASE(dmch)
         IF (dmch=='Y') THEN
          glue_l=.false.
         ELSE
          glue_l=.true. 
         ENDIF
      ELSEIF (ICALL==-5) THEN
         CALL PGSCI(0)
         CALL PGSVP(.05,.26,.88,.93)
         CALL PGSWIN(0.,1.,0.,1.)
         CALL PGRECT(0.,1.,.0,1.)
         CALL PGSVP(.05,.23,.1,.9)
         CALL PGSWIN(0.,1.,0.,1.)
         CALL PGSCI(1)
         IF (act=='U') THEN
          CALL PGPTEXT(.01,.25,0.,0.0,TRIM(title_slc_an))
          CALL PGPTEXT(.01,.2,0.,0.0,'DEFINE THE BORDERS OF')
          CALL PGPTEXT(.01,.16,0.,0.0,'THE ANOMALOUS BODIES AT') 
          CALL PGPTEXT(.01,.12,0.,0.0,'EACH DEPTH USING POLYGONS')
          CALL PGPTEXT(.01,.08,0.,0.0,'C: COPY LAST POLYGON(S)')
         ENDIF
         CALL PGPTEXT(.01,.9,0.,0.0,'NUMBER OF POLYGONS...?')
         ILET=0
         TEXT=''
         CALL PGSCI(15)
         CALL PGRECT(0.,0.3,.83,.86)
         DO
          call pgband(0,1,0.,0.,xdm,ydm,dmch)
          n_pol=ichar(dmch)
          IF (n_pol==13) THEN
             IF(ILET>0) EXIT
             CYCLE
          ELSEIF (n_pol>47.AND.n_pol<58) THEN
             ILET=ILET+1
             TEXT(ILET:ILET)=dmch
          ELSEIF (ICHAR(dmch)==127.OR.ICHAR(dmch)==8) THEN
             IF(ILET<1) CYCLE
             TEXT(ILET:ILET)=''
             ILET=ILET-1
             IF (ILET==0) THEN
                ILET=1
                TEXT=''
             ENDIF
          ELSE
             IF(act=='U'.AND.(ICHAR(dmch)==67.OR.ICHAR(dmch)==99))THEN
              TEXT='-1'
              EXIT
             ELSE
              write(*,*)'INVALID NUMBER OF POLYGONS !'
             ENDIF
          ENDIF
          CALL PGSCI(15)
          CALL PGRECT(0.,0.3,.83,.86)
          CALL PGSCI(1)
          CALL PGPTEXT(.0,.84,0.,0.0,TRIM(TEXT))
         ENDDO
         CALL PGSCI(1) 
         write(*,*)TEXT
         READ(TEXT,'(I3)') n_pol
         CALL PGSCH(.67)
         CALL PGPTEXT(.01,.8,0.,0.0,'**************************')
         CALL PGPTEXT(.0,.77,0.,0.0,'ENTER: LAST POLYGON NODE')
         CALL PGPTEXT(.01,.74,0.,0.0,'**************************')
      ELSEIF (ICALL==-6) THEN
         CALL PGSVP(.0,.26,.9,1.) 
         CALL PGSCI(0)
         CALL PGRECT(0.,1.,0.,1.)
         CALL PGSCI(1)
         CALL PGSVP(.05,.23,.1,.9)
         CALL PGPTEXT(.01,.9,0.,0.0,'INSERT LAYER...')
      ELSEIF (ICALL==6) THEN
         CALL PGSVP(.0,1.,.8,1.)
         CALL PGSWIN(0.,1.,0.,1.)
         CALL PGSCI(1)
         CALL PGSCH(.74)
         IF (MOUSE=='on') THEN
            CALL PGPTEXT(.1,.92,0.,0.0,'s: SAVE, GO TO NEXT ')
            CALL PGPTEXT(.4,.92,0.,0.0,'x: DO NOT SAVE, GO TO NEXT ')
            CALL PGPTEXT(.7,.92,0.,0.0,'m: SAVE, SELECT NEXT ')
            CALL PGPTEXT(.1,.82,0.,0.0,'n: DO NOT SAVE, SELECT NEXT ')
            CALL PGPTEXT(.4,.82,0.,0.0,'d: DELETE LAST POINT ')
            CALL PGSCF(2)
            CALL PGSCH(1.)     
            CALL PGPTEXT(.5,.70,0.,0.5,TRIM(title))
            CALL PGSCF(1)
         ELSE
            CALL PGPTEXT(.1,.9,0.,0.0,'s: NEXT PROFILE ')
            CALL PGPTEXT(.7,.9,0.,0.0,'d: DELETE LAST POINT ') 
         ENDIF 
      ELSEIF (ICALL==7) THEN
         IF(LINUX) THEN
            CALL PGPTEXT(.2,.9,0.,0.0,'SELECT REGION...')
            CALL PGPTEXT(.2,.8,0.,0.0,'CLICK Esc TO EXIT')
         ELSE
            CALL PGPTEXT(.2,.9,0.,0.0,'SELECT REGION...')
            CALL PGPTEXT(.2,.87,0.,0.0,'Place cursor on first corner')
            CALL PGPTEXT(.2,.84,0.,0.0,'   Press "a"')
            CALL PGPTEXT(.2,.81,0.,0.0,'Place cursor on second corner')
            CALL PGPTEXT(.2,.78,0.,0.0,'   Press "a"')
            CALL PGPTEXT(.2,.73,0.,0.0,'CLICK Esc TO EXIT')
         ENDIF
      ELSEIF (ICALL==8) THEN
        CALL PGQLW(illl)
         write(*,*)illl 
         IF(LINUX) THEN
            CALL PGPTEXT(.2,.95,0.,0.0,'MOUSE INPUT') 
            IF (ort=='n'.OR.ort=='N') THEN 
               CALL PGPTEXT(.2,.9,0.,0.0,'MOVING EAST:')
               CALL PGPTEXT(.2,.7,0.,0.0,'MOVING WEST:')
               CALL PGPTEXT(.2,.65,0.,0.0,
     *              'ERASE PROFILES EAST OF CURSOR')
            ELSE
               CALL PGPTEXT(.2,.9,0.,0.0,'MOVING NORTH:')
               CALL PGPTEXT(.2,.7,0.,0.0,'MOVING SOUTH:')
               CALL PGPTEXT(.2,.65,0.,0.0,
     *              'ERASE PROFILES NORTH OF CURSOR')  
            ENDIF 
            CALL PGPTEXT(.2,.85,0.,0.0,'LEFT CLICK: INTERPOLATE')
            CALL PGPTEXT(.2,.8,0.,0.0,'RIGHT CLICK: SAVE ')
            CALL PGPTEXT(.2,.55,0.,0.0,'ENTER: EXIT')
            CALL PGPTEXT(.2,.45,0.,0.0,'Esc: ABORT')
         ELSE
            CALL PGPTEXT(.2,.95,0.,0.0,'Keyboard input') 
            IF (ort=='n'.OR.ort=='N') THEN 
               CALL PGPTEXT(.2,.9,0.,0.0,'MOVING EAST:')
               CALL PGPTEXT(.2,.7,0.,0.0,'MOVING WEST:')
               CALL PGPTEXT(.2,.65,0.,0.0,
     *              'ERASE PROFILES EAST OF CURSOR')
            ELSE
               CALL PGPTEXT(.2,.9,0.,0.0,'MOVING NORTH:')
               CALL PGPTEXT(.2,.7,0.,0.0,'MOVING SOUTH:')
               CALL PGPTEXT(.2,.65,0.,0.0,
     *              'ERASE PROFILES NORTH OF CURSOR') 
            ENDIF 
            CALL PGPTEXT(.2,.85,0.,0.0,'"a": INTERPOLATE')
            CALL PGPTEXT(.2,.8,0.,0.0,'"s": SAVE ')
            CALL PGPTEXT(.2,.55,0.,0.0,'ENTER: EXIT')
            CALL PGPTEXT(.2,.45,0.,0.0,'Esc: ABORT')
         ENDIF
      ELSEIF (ICALL==9) THEN
         CALL PGSCI(2)
         CALL PGSLW(1)
         CALL PGPTEXT(.2,.9,0.,0.0,'REALLY ABORT ? (y/n)')
         call pgband(0,1,0.,0.,xdm,ydm,abort_mdreg)
         CALL PGSCI(1)
      ELSEIF (ICALL==10) THEN
         CALL PGSCI(3)
         CALL PGPTEXT(.2,.9,0.,0.0,'REALLY EXIT ? (y/n)')
         CALL PGSCI(1)
         call pgband(0,1,0.,0.,xdm,ydm,exit_mdreg)
      ELSEIF (ICALL==11) THEN
         CALL PGSCI(2)
         CALL PGPTEXT(.17,.9,0.,0.0,'REALLY ERASE MODIFICATIONS? (y/n)')
         CALL PGSCI(1)
         call pgband(0,1,0.,0.,xdm,ydm,ms_flow) 
      ELSEIF (ICALL==12) THEN
         CALL PGPTEXT(.17,.9,0.,0.0,'PRESS ANY KEY TO CONTINUE...')
      ENDIF

      CALL PGSVP(xvp_min,xvp_max,yvp_min,yvp_max)
      CALL PGSWIN (xw_min,xw_max,yw_min,yw_max) 
      CALL PGSCH(1.)

      END SUBROUTINE
