CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PROF_PREP(D,Z_LMAX,PER_L)
C     ===================
C
C***********************************************************************
C
C Subroutine performs the GMT commands necessary to plot a profile out
C of the 3D data file under Windows (same as profile.job for Linux)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      REAL*4 LON_MIN,LON_MAX,LAT_MIN,LAT_MAX,LON_INI,LAT_INI,LON_FIN,
     *       LAT_FIN
      ALLOCATABLE :: XP(:),YP(:)
      CHARACTER PASO_MALLA*3,NAME*100,COL*11,REGION*50,GRID*20,BD*10,
     *          START*20,FINISH*30,TEXT*100,FILE*20,PROJECT*10,T*200,
     *          F1*5,F2*5
C      CALL SYSTEMQQ('GMT gmtset D_FORMAT %lg')
      OPEN(55,FILE='info.dat',STATUS='OLD')
      READ(55,*) IDUM,LON_MIN,LON_MAX,LAT_MIN,LAT_MAX
      CLOSE(55)
      paso_malla='3m'

C Coordenadas de los extremos del perfil
C gawk '{print $0}' qq |
C read lon_ini lat_ini lon_fin lat_fin i i_lay N_lay name sense
      OPEN(55,FILE='qq',STATUS='OLD')
      READ(55,*) LON_INI,LAT_INI,LON_FIN,LAT_FIN,I1,I_LAY,N_LAY,NAME,
     *           I_SENSE
!!      CLOSE(55,DISP='DELETE')
      BD='bod'
      WRITE(BD(4:10),'(I0)') I_LAY
      F1='l'
      WRITE(F1(2:5),'(I0)') I_LAY

C Color
C awk '{if (NR=='$i_lay') print $0}' color.rgb |
C read R G B
      OPEN(55,FILE='color.rgb',STATUS='OLD')
      DO I=1,I_LAY
       READ(55,*) IR,IG,IB
      ENDDO
      CLOSE(55)
      WRITE(COL,'(I0,''/'',I0,''/'',I0)') IR,IG,IB
      WRITE(REGION,'(''-R'',3(F0.3,''/''),F0.3)') LON_MIN,LON_MAX,
     *      LAT_MIN,LAT_MAX
      WRITE(GRID,'(''-I'',A)') TRIM(PASO_MALLA)

C col=$R/$G/$B
C surface $name -Glayer.grd  -R$lon_min/$lon_max/$lat_min/$lat_max -I3m
      T='surface '//TRIM(NAME)//' -Glayer.grd '//TRIM(REGION)//' '//
     *  TRIM(GRID)
      WRITE(15,'(A)') TRIM(T)

C project -C$lon_ini/$lat_ini -E$lon_fin/$lat_fin -G$d -Q  >profile_project.xyp
C grdtrack profile_project.xyp -Glayer.grd >$name'p'
      T='grdtrack profile_project.xyp -Glayer.grd > '//TRIM(F1)//'p'
      WRITE(15,'(A)') TRIM(T)

C minmax -C $name'p' |awk '{print $6}' |
C read  L_per

C if [ $i_lay -eq 1 ] ; then
      WRITE(REGION,'(''-R0/'',F0.1,''/-'',F0.1,''/10'')') PER_L,Z_LMAX
      WRITE(PROJECT,'(''-JX25/15'')')
      IF(I_LAY .EQ. 1) THEN

C #FIRST CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C grdtrack profile_project.xyp -Gelev.grd |
C    awk '{print $3,$4*1E-3}' >profile_elev.xypz
         T='grdtrack profile_project.xyp -Gelev.grd | '//
     *     'gawk "{print $3,$4*1E-3}" > profile_elev.xypz'
         WRITE(15,'(A)') TRIM(T)
C awk '{print $3,-$4}' $name'p'  >old
C sort -r -n  old -o old
         T='gawk "{print $3,-$4}" '//TRIM(F1)//'p > '//TRIM(F1)
         WRITE(15,'(A)') TRIM(T)
         WRITE(15,'(''inverse '',A)') TRIM(F1)
C #Water layer

C psxy  -JX25/15 -R0/$L_per/-$z_lmax/10 -L -K -G0/250/250 -W1p/0/0/0
C -Bf50g$L_per:" (km)":/f10g300:" Depth (km)":WS:."Profile $lon_ini ,
C $lat_ini ; $lon_fin , $lat_fin ": << END> perfil_$i.ps
C 0 0
C $L_per 0
C $L_per -10
C 0 -10
C END
         OPEN(55,FILE='water')
         WRITE(55,'(''0 0'')')
         WRITE(55,'(F0.1,'' 0'')') PER_L
         WRITE(55,'(F0.1,'' -10'')') PER_L
         WRITE(55,'(''0 -10'')')
         CLOSE(55)
         WRITE(FILE,'(''perfil_'',I0,''.ps'')') I1
         T='psxy water '//TRIM(REGION)//' '//TRIM(PROJECT)//
     *     ' -L -K -G0/250/250 -W1p > '//TRIM(FILE)
         WRITE(15,'(A)') TRIM(T)

C cat profile_elev.xypz old > bod
         T='COPY profile_elev.xypz+'//TRIM(F1)//' bod'
         WRITE(15,'(A)') TRIM(T)

C psxy bod  -JX25/15 -R0/$L_per/-$z_lmax/10 -K -O -L -G$col -W1p/0/0/0 >> perfil_$i.ps
         T='psxy bod '//TRIM(REGION)//' '//TRIM(PROJECT)//
     *     ' -K -O -L -G'//TRIM(COL)//' -W1p >> '//TRIM(FILE)
         WRITE(15,'(A)') TRIM(T)
         F2=F1

C #FIRST CALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C else

C if [ $i_lay -eq $N_lay ] ; then
      ELSEIF(I_LAY .EQ. N_LAY) THEN

C #LAST CALL*********************************************************
C #+++++++++++++++++++++++++++++++++++++++
C if [ $sense -eq 0 ] ; then
         IF(I_SENSE .EQ. 0) THEN
C awk '{print $3,-$4}' $name'p' >new.c
            T='gawk "{print $3,-$4}" '//TRIM(F1)//'p > '//TRIM(F1)
            WRITE(15,'(A)') TRIM(T)
C cat  old  new.c > bod
            WRITE(15,'(''Copy '',A,''+'',A,'' bod'')') TRIM(F2),
     *            TRIM(F1)
            F2=F1
C awk '{print $1,$2}'<<END>asth
C $L_per -$z_lmax
C 0 -$z_lmax
C END
            OPEN(55,FILE='asth')
            WRITE(55,'(F0.3,''  '',F0.3)') PER_L,-Z_LMAX
            WRITE(55,'(''0  '',F0.3)') -Z_LMAX
            CLOSE(55)
C cat new.c asth > asth.c
            WRITE(15,'(''COPY '',A,''+asth asth_c'')') TRIM(F1)
C else
         ELSE
C awk '{print $3,-$4}' $name'p' >new.c
C sort -r -n new.c -o new.c
            T='gawk "{print $3,-$4}" '//TRIM(F1)//'p > '//TRIM(F1)
            WRITE(15,'(A)') TRIM(T)
C cat  old  new.c > bod
            WRITE(15,'(''Copy '',A,''+'',A,'' bod'')') TRIM(F2),
     *            TRIM(F1)
            F2=F1
C awk '{print $1,$2}'<<END>asth
C 0 -$z_lmax
C $L_per -$z_lmax
C END
            OPEN(55,FILE='asth')
            WRITE(55,'(''0  '',F0.3)') -Z_LMAX
            WRITE(55,'(F0.3,''  '',F0.3)') PER_L,-Z_LMAX
            CLOSE(55)
C cat new.c asth > asth.c
            WRITE(15,'(''COPY '',A,''+asth bod'')') TRIM(F1)
C fi
         ENDIF
C #+++++++++++++++++++++++++++++++++++++++
C #cccccccccccccccccccccccccccccc
C #cccccccccccccccccccccccccccccc
C psxy bod  -JX25/15 -R0/$L_per/-$z_lmax/10 -W1p/0/0/0 -L -G$col -K -O >>perfil_$i.ps
         T='psxy bod '//TRIM(REGION)//' '//TRIM(PROJECT)//
     *     ' -W1p -L -G'//TRIM(COL)//' -O -K >> '//TRIM(FILE)
         WRITE(15,'(A)') TRIM(T)
C psxy asth.c  -JX25/15 -R0/$L_per/-$z_lmax/10 -W1p/0/0/0 -L -G250/5/250 -O >>perfil_$i.ps
         WRITE(TEXT,'(''-B50g100:"Distance (km)":/10g20:"Depth (km)":'',
     *         '':".Profile '',3(F0.2,'', ''),F0.2,''":WSe > '')')
     *         LON_INI,LAT_INI,LON_FIN,LAT_FIN
         T='psxy asth_c '//TRIM(REGION)//' '//TRIM(PROJECT)//' '//
     *     TRIM(TEXT)//' -W1p -L -G250/5/250 -O >> '//TRIM(FILE)
         WRITE(15,'(A)') TRIM(T)
C rm *.xypz *.xyp *xyzp layer.grd old bod asth asth.c new.c
         T='DEL *.xypz *.xyp *.xyzp layer.grd old bod asth asth_c '//
     *     'new_c water scratch'
         WRITE(15,'(A)') TRIM(T)
C gv perfil_$i.ps &
         WRITE(15,'(''gsview32 '',A)') TRIM(FILE)
C #LAST CALL*********************************************************

C else
      ELSE
C #INTERMEDIATE CALLS/////////////////////////////////////////////////
C #+++++++++++++++++++++++++++++++++++++++
C if [ $sense -eq 0 ] ; then
         IF(I_SENSE .EQ. 0) THEN
C awk '{print $3,-$4}' $name'p' >new.c
            T='gawk "{print $3,-$4}" '//TRIM(F1)//'p > '//TRIM(F1)
            WRITE(15,'(A)') TRIM(T)
            F2=F1
C else
         ELSE
C awk '{print $3,-$4}' $name'p' >new.c
C sort -r -n new.c -o new.c
            T='gawk "{print $3,-$4}"'//TRIM(F1)//'p > '//TRIM(F1)
            WRITE(15,'(A)') TRIM(T)
            WRITE(15,'(''inverse '',A)') TRIM(F1)
C cat  old  new.c > bod
             WRITE(15,'(''COPY '',A,''+'',A,'' bod'')') TRIM(F2),
     *             TRIM(F1)
             F2=F1
C fi
         ENDIF
C #+++++++++++++++++++++++++++++++++++++++
C #cccccccccccccccccccccccccccccc
C psxy bod  -JX25/15 -R0/$L_per/-$z_lmax/10 -W1p/0/0/0 -L -G$col  -K -O >>perfil_$i.ps
         T='psxy bod '//TRIM(REGION)//' '//TRIM(PROJECT)//
     *     ' -W1p -L -G'//TRIM(COL)//' -K -O >> '//TRIM(FILE)
         WRITE(15,'(A)') TRIM(T)
C mv new.c old
C #INTERMEDIATE CALLS/////////////////////////////////////////////////
C fi
      ENDIF

C fi

C rm qq
C      CALL SYSTEMQQ('DEL qq')
C #awk '{print $1,-$2}' /home/jfullea/Tesis/datos/terremotos/terre.dh |
C #psxy -V -Jx$esc_x/.05 -R0/$L_per/-$z_lmax/10 -O -W3/0/0/0 -Sc0.25 -G0/0/255>>perfil_cget.ps

C STILL TO BE IMPLEMENTED

      END
