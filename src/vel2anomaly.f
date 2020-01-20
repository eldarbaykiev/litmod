! vel2anomaly.f
      
      LOGICAL STAT
      CHARACTER(100)name_ref,name_out
      CHARACTER(1)at
      REAL(8),ALLOCATABLE::MDL(:,:),MDL_dm(:,:),VREF(:,:),
     *                     vp_vec(:,:,:),vs_vec(:,:,:) 
      INTEGER N_x,N_y,N_z,N_ref,I_z1,I_z2,I_zz
      REAL(8) d_z,vp1,vp2,vs1,vs2,z1,z2,z,slp_p,slp_s
      REAL(8) L_x,L_y,Topo_max
      
      WRITE(*,*)'*****************************************************'
      WRITE(*,*)'    PROGRAM TO SUBSTRACT A REFERENCE SEISMIC MODEL'
      WRITE(*,*)'      FROM THE VELOCITIES OBTAINED WITH LITMOD3D    ' 
      WRITE(*,*)''
      WRITE(*,*)'  A reference 1D model is necessary (e.g. PREM)      '
       WRITE(*,*)'Two reference models (¨smooth PREM" and ak135) are'
       WRITE(*,*)        '  provided with this program '
       write(*,*)' '
       WRITE(*,*)'If you want to use another 1D model, use this format:'
      WRITE(*,*)'FORMAT: name.xy; 1 COLUMN: DEPTH (KM) 2-3 COLUMNS: Vp,
     *Vs'
      WRITE(*,*)'*****************************************************'
      WRITE(*,*)'name of the 1D reference model(file name without .xy)?'
      READ(*,*)name_ref
      WRITE(*,*)'CONSIDER SEISMIC ATTENUATION...(y/n DEF=y)?'
      WRITE(*,*)'(Note: remember that comparing elastic velocities'
      WRITE(*,*)'with a reference model such as PREM is physically'
      WRITE(*,*)'incorrect. Recommendation: Run atten.f before)'

      READ(*,*)at
      name_out=TRIM(name_ref)//'_ref.xy'
      name_ref=TRIM(name_ref)//'.xy'

      INQUIRE(FILE='LITMOD3D.info',EXIST=STAT)
      IF (STAT) THEN
       OPEN(15,file='LITMOD3D.info',STATUS='OLD')
      ELSE
       WRITE(*,*) 'FILE LITMOD3D.info NOT FOUND!'
       STOP
      ENDIF
      READ(15,*)L_x,L_y,Topo_max,N_x,N_y,N_z,d_z
      CLOSE(15)

!      WRITE(*,*)'NUMBER OF NODES AND VERTICAL GRID STEP'
!      WRITE(*,*)'N_x,N_y,N_z,d_z (km)...?'
!      READ(*,*) N_x,N_y,N_z,d_z

      OPEN(1,file=TRIM(name_ref),status='OLD')
      OPEN(2,file=TRIM(name_out),status='unknown') 
      N_ref=0
      DO 
       READ (1,*,IOSTAT=ios) 
       IF (IOS<0) EXIT
       N_ref=N_ref+1
      ENDDO
      REWIND(1)
      ALLOCATE (MDL(N_ref,3),MDL_dm(N_ref,3),VREF(N_z,2))
      VREF=0
      MDL=0
      MDL_dm=0
      DO i=1,N_ref
       READ (1,*)MDL(i,1),MDL(i,2),MDL(i,3)
      ENDDO 
      CLOSE(1)
!Reorder, if necessary
      IF (MDL(1,1)>MDL(N_ref,1))THEN
       MDL_dm(:,:)=MDL(:,:)
       DO ik=1,N_ref
        MDL(ik,:)=MDL_dm(N_ref-ik+1,:)
       ENDDO
      ENDIF

!Avoid sharp discontinuities
      DO i=1,N_ref-1
       IF (MDL(i,1)==MDL(i+1,1)) THEN
        MDL(i,1)=MDL(i,1)-d_z
        MDL(i+1,1)=MDL(i+1,1)+d_z
       ENDIF
      ENDDO

!Redefine the reference column to have N_z nodes     
      DO i=1,N_ref-1
       vp1=MDL(i,2)
       vp2=MDL(i+1,2)
       vs1=MDL(i,3)
       vs2=MDL(i+1,3)
       I_z1=NINT((MDL(i,1)+Topo_max)/d_z) 
       I_z2=NINT((MDL(i+1,1)+Topo_max)/d_z)
       IF(I_z1<0)I_z1=1
       IF(I_z2>N_z)I_z2=N_z
       z1=((I_z1-1)*d_z)
       z2=((I_z2-1)*d_z)
       slp_p=(vp2-vp1)/(z2-z1)
       slp_s=(vs2-vs1)/(z2-z1)
       DO I_zz=I_z1+1,I_z2
        z=(I_zz-1)*d_z
        VREF(I_zz,1)=vp1+slp_p*(z-z1)
        VREF(I_zz,2)=vs1+slp_s*(z-z1) 
       ENDDO
      ENDDO 
!Smooth the profile
!      I_win=5
!      DO i=I_win+1,N_z-I_win 
!       VREF(i)=SUM(VREF(i-I_win:i+I_win))/(2*I_win+1)
!      ENDDO
!Write the reference column
      DO i=1,N_z
       z=-Topo_max+(i-1)*d_z
       WRITE(2,*)z,VREF(i,1),VREF(i,2)
      ENDDO
      CLOSE(2)

      IF(at=='n'.OR.at=='N') THEN
       INQUIRE(FILE='./OUT_xyz/Trhopvels.xyz',EXIST=STAT)
       IF (.NOT.STAT) THEN
        WRITE(*,*)'FILE OUT_xyz/Trhopvels.xyz NOT FOUND!'
        STOP
       ENDIF
       OPEN(30,file='./OUT_xyz/Trhopvels.xyz',STATUS='OLD')
      ELSE

       INQUIRE(FILE='./OUT_xyz/V_attenuation.xyz',EXIST=STAT)
       IF (.NOT.STAT) THEN
        WRITE(*,*)'FILE OUT_xyz/V_attenuation.xyz NOT FOUND!'
        STOP
       ENDIF
       OPEN(30,file='OUT_xyz/V_attenuation.xyz',STATUS='OLD')
      ENDIF 

      OPEN(31,file='OUT_xyz/V_anomaly.xyz',STATUS='UNKNOWN')
      ALLOCATE (vp_vec(N_x,N_y,N_z),vs_vec(N_x,N_y,N_z))      
       DO j=1,N_y
        DO i=1,N_x
         DO k=1,N_z
          IF(at=='n'.OR.at=='N') THEN
           READ(30,*)x,y,z,dm,dm,dm,vp_vec(i,j,k),vs_vec(i,j,k)
          ELSE
           READ(30,*)x,y,z,vp_vec(i,j,k),vs_vec(i,j,k)
          ENDIF
          IF(vp_vec(i,j,k)/=0) THEN
           vp_vec(i,j,k)=100*(vp_vec(i,j,k)-VREF(k,1))/VREF(k,1)
           vs_vec(i,j,k)=100*(vs_vec(i,j,k)-VREF(k,2))/VREF(k,2)
          ENDIF
          WRITE(31,'(3F13.3,F8.4,1X,F8.4)')x,y,z,vp_vec(i,j,k),
     *                                     vs_vec(i,j,k) 
         ENDDO
        ENDDO
       ENDDO
       CLOSE(30)
       CLOSE(31)

      END

