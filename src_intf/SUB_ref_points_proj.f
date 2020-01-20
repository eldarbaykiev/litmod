! SUB_ref_points_proj.f

      SUBROUTINE ref_points_proj
      USE M_profile
      USE M_ref_points
      REAL(4) xp_min,xp_max,yp_min,yp_max,p
      REAL(4) xv_min,xv_max,yv_min,yv_max,x_dm,x_0,y_0
      REAL(4) xy_crss(2),z_crss(2),xy_pr_crss
      CALL PGSLW(10)
      DO i=1,N_RPTS
       xp_min=minval(pol_RPTS(2*i-1,1:N_pol_RPTS(i)))
       xp_max=maxval(pol_RPTS(2*i-1,1:N_pol_RPTS(i)))
       yp_min=minval(pol_RPTS(2*i,1:N_pol_RPTS(i)))
       yp_max=maxval(pol_RPTS(2*i,1:N_pol_RPTS(i)))      
       IF (xp_min>lnz_min.OR.xp_max<lnz_max.OR.
     *     yp_min>ltz_min.OR.yp_max<ltz_max) THEN
        DO j=1,N_pol_RPTS(i)
        xv_min=pol_RPTS(2*i-1,j)
        yv_min=pol_RPTS(2*i,j)
        xv_max=pol_RPTS(2*i-1,j+1)
        yv_max=pol_RPTS(2*i,j+1) 
!Rearrange
        IF (xv_min>xv_max) THEN
         x_dm=xv_min
         xv_min=xv_max
         xv_max=x_dm
        ENDIF
        IF (yv_min>yv_max) THEN
         x_dm=yv_min
         yv_min=yv_max
         yv_max=x_dm
        ENDIF
!Slope
        p=(yv_max-yv_min)/(xv_max-xv_min)

        IF (ort=='E') THEN
         y_0=xy_pr_p(i_pr) 
         IF(y_0>yv_min.AND.y_0<yv_max)THEN
          xy_pr_crss=xv_min+(y_0-yv_min)/p
          xy_crss(1)=xy_pr_crss
          xy_crss(2)=xy_pr_crss
          z_crss(1)=-z_max_mdfr
          z_crss(2)=5.
          CALL PGSLW(5)
          CALL PGSLS(2)
          CALL PGLINE(2,xy_crss,z_crss)
          CALL PGSLW(3)
          CALL PGSLS(1)
         ENDIF    
        ELSE
         x_0=xy_pr_p(i_pr)
         IF(x_0>xv_min.AND.x_0<xv_max)THEN
          xy_pr_crss=p*(x_0-xv_min)+yv_min 
          xy_crss(1)=xy_pr_crss
          xy_crss(2)=xy_pr_crss
          z_crss(1)=-z_max_mdfr
          z_crss(2)=5.
          write(*,*)'xy_crss,i',xy_crss,i
          CALL PGSLW(5)
          CALL PGSLS(2)
          CALL PGLINE(2,xy_crss,z_crss)
          CALL PGSLW(3)
       CALL PGSLS(1)
         ENDIF
        ENDIF

        ENDDO
       ENDIF

      ENDDO 
      CALL PGSLW(3)
      END SUBROUTINE
