CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE BODY_PREP(IBODY)
C     ====================
C
C***********************************************************************
C
C Subroutine prepares a batch file to plot the bodies chosen in
C SUB_body_plot_ps
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
D      USE DFLIB
      REAL*4 LON_MIN,LON_MAX,LAT_MIN,LAT_MAX,LON_INI,LAT_INI,LON_FIN,
     *       LAT_FIN
      CHARACTER FILE*100,SCALE*10,TEXT*200,PSFILE*100
      OPEN(55,FILE='info.dat',STATUS='OLD')
      READ(55,*) IDUM,LON_MIN,LON_MAX,LAT_MIN,LAT_MAX
      CLOSE(55)
C #!/bin/ksh 
C # body_plot.job
C set lon_min=-12
C set lon_max=0
C set lat_min=28
C set lat_max=39

C awk '{print $0}' qq |
C read min max name
C echo $min $max $name
      OPEN(55,FILE='qq',STATUS='OLD')
      READ(55,*) HMIN,HMAX,FILE
      CLOSE(55)
      OPEN(55,FILE='body_plot.bat')
      WRITE(PSFILE,'(''body'',I0,''_'',A,''.ps'')') IBODY,TRIM(FILE)
C if [ $min -lt 0 ] ; then
C ((int=($max+$min)/10))
C else
C ((int=($max-$min)/10))
C fi
      IF(HMIN .LT. 0) THEN
         INT=CEILING((HMAX+HMIN)/10)
      ELSE
         INT=CEILING((HMAX-HMIN)/10)
      ENDIF
      NINT=(HMAX-HMIN)/INT
      HMAX=HMIN+INT*NINT
C #min=0
C #int=5
C echo $name $min $max $int
C #Proyección
C #p=m1.25
C echo $lon_min $lon_max $lat_min $lat_max |
C awk '{print ($2-$1),($4-$3)}' |
C read delta_lon delta_lat
      DELTA_LON=LON_MAX-LON_MIN
      DELTA_LAT=LAT_MAX-LAT_MIN
C if [ $delta_lon -le $delta_lat ]
C  then
C echo $delta_lat |
C awk '{print (10.24/$1)}' |
C read p
C p=m$p
C  else
C echo $delta_lon |
C awk '{print (15.36/$1)}' |
C read p
C p=m$p
C fi
      IF(DELTA_LON .LE. DELTA_LAT) THEN
         P=10.24/DELTA_LAT
      ELSE
         P=15.36/DELTA_LON
      ENDIF
      WRITE(SCALE,'(''-Jm'',F0.3)') P
C echo "ESCALA ..."
C echo $p

C #surface bod_thick.xyz -Gbod.grd  -R$lon_min/$lon_max/$lat_min/$lat_max -I2m

C xyz2grd bod_thick.xyz -Gbod.grd  -R$lon_min/$lon_max/$lat_min/$lat_max  -I3m
      WRITE(TEXT,'(''-R'',3(F0.2,''/''),F0.2)') LON_MIN,LON_MAX,
     *      LAT_MIN,LAT_MAX
      WRITE(55,'(''xyz2grd bod_thick.xyz -Gbod.grd '',A,'' -I3m'')')
     *      TRIM(TEXT)
      
C grd2cpt bod.grd -Chaxby -S$min/$max/$int -I >bod.cpt
C      WRITE(55,'(''grd2cpt bod.grd -Chaxby -S'',2(F0.0,''/''),I0,
      WRITE(55,'(''makecpt -Chaxby -T'',2(F0.0,''/''),I0,
     *      '' -I > bod.cpt'')') HMIN,HMAX,INT

C grdimage bod.grd -R$lon_min/$lon_max/$lat_min/$lat_max -J$p -Cbod.cpt  -K -B2f1:." $name ": -E150 -Y1  >bod.ps
      WRITE(55,'(''grdimage bod.grd '',A,'' '',A,'' -Cbod.cpt -K '',
     *      ''-B2f1:."'',A,''": -E150 -Y1 > '',A)') TRIM(TEXT),
     *      TRIM(SCALE),TRIM(FILE),TRIM(PSFILE)
C psscale -Cbod.cpt -D18/7.9/17/0.9 -V -K -O -B::/:"(km)": >>bod.ps
      WRITE(55,'(''psscale -Cbod.cpt -D18/7.9/17/0.9 -V -K -O -B::/'',
     *      '':"(km)": >> '',A)') TRIM(PSFILE)
C grdcontour bod.grd -R$lon_min/$lon_max/$lat_min/$lat_max -J$p -K -O -W2/0/0/0  -Cbod.cpt -Af12a0 >>bod.ps
      WRITE(55,'(''grdcontour bod.grd '',A,'' '',A,'' -K -O '',
     *      ''-W2/0/0/0 -Cbod.cpt -Af12a0 >> '',A)') TRIM(TEXT),
     *      TRIM(SCALE),TRIM(PSFILE)
C pscoast -R$lon_min/$lon_max/$lat_min/$lat_max -J$p  -Dh -W2p -V -O -A500 >>bod.ps
      WRITE(55,'(''pscoast '',A,'' '',A,'' -Dh -W2p -V -O -A500 '',
     *      ''>> '',A)') TRIM(TEXT),TRIM(SCALE),TRIM(PSFILE)
C gv bod.ps &
      WRITE(55,'(''gsview32 '',A)') TRIM(PSFILE)
C rm qq bod.cpt bod.grd 
      WRITE(55,'(''DEL qq bod.cpt bod.grd'')')
      CLOSE(55)
      CALL SYSTEM('body_plot.bat')
      END
