! Subroutine SUB_label_proj.f

      SUBROUTINE label_proj
      USE M_profile
      USE M_label
      REAL(4) xy_lbl,xy_lbl_p
      REAL(4) xbox(4),ybox(4) 

!Find marks within the band +-bnd_lbl
      lbl_proy=0 
      DO i=1,NMARK
       IF (ort=='E') THEN
       xy_lbl=labels(i)%lat
       xy_lbl_p=labels(i)%long
       ELSE
       xy_lbl=labels(i)%long
       xy_lbl_p=labels(i)%lat
       ENDIF 
        IF (ABS(xy_pr_p(i_pr)-xy_lbl)<=bnd_lbl.AND.xy_lbl_p>xy_pr_min
     *      .AND.xy_lbl_p<xy_pr_max ) THEN
        IF (ort=='E') THEN
         xy_lbl=labels(i)%long
        ELSE
         xy_lbl=labels(i)%lat
        ENDIF
        CALL PGPT1(xy_lbl,labels(i)%z_lbl,18)
        CALL PGQTXT(xy_lbl+.03,labels(i)%z_lbl,0.,0.,
     *              TRIM(labels(i)%text_mark),xbox,ybox)
        CALL PGSCI(0)
        CALL PGRECT(xbox(1),xbox(3),ybox(1),ybox(2))
        CALL PGSCI(1)
        CALL PGPTEXT(xy_lbl+.03,labels(i)%z_lbl,0.,0.,
     *               TRIM(labels(i)%text_mark))
        ENDIF
      ENDDO

      END SUBROUTINE 
