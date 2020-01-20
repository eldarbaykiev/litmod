! Subroutine SUB_get_pol.f

      SUBROUTINE get_pol
      USE M_imag
      USE M_pol
      USE M_sublitho_an
!      USE M_imag, ONLY:act
      INTEGER NPT,COL
      INTEGER, PARAMETER::MAXPT=100
      real(4),dimension(MAXPT) :: x_vec,y_vec
      real(4),dimension(2) :: x_2p,y_2p
      real(4) x,y
      character le
      !Introduce the new body limits
      COL=0
      write(*,*)'MÁX NUMBER OF POINTS-->',MAXPT
      call pgqvp(0,r1,r2,r3,r4)
      IF (act=='U') THEN
       i_pan=IPG_SLC
      ELSE
       i_pan=ISTAT 
       CALL clear_scr(0)
      ENDIF
!Plot last polygon
      IF (act=='U'.AND.I_z_an>Npr_v_up) THEN
       write(*,*)'n_pol_lst',n_pol_lst,MAXVAL(NP_POL_lst(:)),
     *  MAXVAL(A_pol_lst(:,:))
       DO i=1,n_pol_lst
        CALL PGSLW(5)
        CALL PGSLS(3)
        CALL PGSCI(3)
        CALL PGPOLY(NP_POL_lst(i),A_pol_lst((2*i-1),:),A_pol_lst(2*i,:))
        CALL PGLINE(NP_POL_lst(i),A_pol_lst((2*i-1),:),A_pol_lst(2*i,:))
        CALL PGSLW(1)
        CALL PGSLS(1)
       ENDDO
      ENDIF
      CALL text_plot(i_pan,r1,r2,r3,r4,lon_min,lon_max,lat_min,lat_max,
     *               -5)

!      write(*,*)'NUMBER OF POLYGONS...?'
!      read(*,*) n_pol
!      write(*,*)'******************************'
!      write(*,*)'/ or * --> last polygon node'
!      write(*,*)'******************************'
      write(*,*)'n_pol',n_pol

      IF(n_pol==-1.AND.I_z_an>Npr_v_up) THEN
       A_pol=A_pol_lst
       NP_POL=NP_POL_lst 
       n_pol=n_pol_lst
       WRITE(*,*)'COPYING POLYGON(S) FROM LAST SLICE'
      ELSE

      IF(ALLOCATED(A_pol)) DEALLOCATE(A_pol)
      IF(ALLOCATED(NP_POL)) DEALLOCATE(NP_POL)
      ALLOCATE (A_pol(n_pol*2,MAXPT),NP_POL(n_pol))
      IF(ALLOCATED(A_pol_lst)) DEALLOCATE(A_pol_lst)
      IF(ALLOCATED(NP_POL_lst)) DEALLOCATE(NP_POL_lst)
      ALLOCATE (A_pol_lst(n_pol*2,MAXPT),NP_POL_lst(n_pol))
      DO i=1,n_pol
      write(*,*)'POLYGON NUM-->',i
      MODE=0
      COL=COL+1
      NPT = 0
      x=(lon_max+lon_min)/2
      y=(lat_max+lat_min)/2
      le='a'
      CALL PGSCI(0)
      CALL PGSLW(4)
      DO 
       IF (NPT==0) MODE=0
       call pgband(MODE,1,x,y,x,y,le)
       IF (ICHAR(le)==13) THEN
        IF (NPT>0) EXIT
        CYCLE
       ENDIF
!Check for zoom
       ik=i-1
       IF (le.EQ.'z') call zoom(ik)
       IF (le.EQ.'u') call clear_scr (ik)
       MODE=1
      IF (le.NE.'z'.AND.le.NE.'u') THEN
!PGPT1 dibuja los puntos
       CALL PGSCI(1)
       call PGPT1(x,y,18)
       CALL PGSCI(0)
       write(*,*)'COORDINATES',x,y
       NPT=NPT+1
       x_vec(NPT)=x
       y_vec(NPT)=y

        IF (NPT>=2) THEN
        x_2p=x_vec(NPT-1:NPT)
        y_2p=y_vec(NPT-1:NPT)
        CALL PGLINE(2,x_2p,y_2p)
        ENDIF
      
      ENDIF

      ENDDO

      IF (NPT.GE.3) THEN
           x_vec(NPT+1)=x_vec(1)
           y_vec(NPT+1)=y_vec(1)

           CALL PGSCI(COL)
           CALL PGSFS(3)
           CALL PGPOLY(NPT,x_vec,y_vec)
           CALL PGLINE(NPT+1,x_vec,y_vec)
!           CALL PGSLCT(IPS)
!           CALL PGSCI(COL)
!           CALL PGPOLY(NPT,x_vec,y_vec)
!           CALL PGLINE(NPT+1,x_vec,y_vec)
!           CALL PGSLCT(ISTAT)
      ENDIF
       A_pol((2*i-1),:)=x_vec
       A_pol(2*i,:)= y_vec
       NP_POL(i)=NPT
       write(*,*)'****************************'

      ENDDO
      CALL PGSLW(1)

      NP_POL=NP_POL+1

      ENDIF
!Copy last polygon defined
      A_pol_lst=A_pol
      NP_POL_lst=NP_POL
      n_pol_lst=n_pol

      END SUBROUTINE
