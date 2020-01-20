! Subroutine  SUB_label.f

      SUBROUTINE label
      USE M_imag
      USE M_label
      INTEGER NMARK_old
      CHARACTER KEY,GETCHARQQ,dmch


      CALL text_plot(ISTAT,.3,.9,.1,.9,lon_min,lon_max,lat_min,lat_max,
     *                -2)

      IF (lbl_act=='R') THEN
       NMARK=0
       open(911,file='labels.dat',status='UNKNOWN')
      ENDIF
      NMARK_old=NMARK+1 

      LBL_SGN=0
      DO 
        IF (LBL_SGN==1) EXIT
       CALL lay_sel
      ENDDO 
!Write new labels in labels.dat
      IF(NMARK_old<=NMARK) THEN
       DO i=NMARK_old,NMARK
       WRITE(911,*)labels(i)%long,labels(i)%lat,
     *                         labels(i)%z_lbl,
     *                         TRIM(labels(i)%text_mark)
       ENDDO
      ENDIF

      END SUBROUTINE
