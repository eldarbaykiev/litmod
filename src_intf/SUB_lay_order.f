! SUB_lay_order.f
!Checks the hierarchy for layer number i_modf with respect to rest of the 
!layers in the point i,j.
!It actuates over the matrix dm_lay, e.g: 
!CALL lay_order(i_gl,j_gl,i_lay)
!layers(i_gl,j_gl,:)=dm_lay(i_gl,j_gl,:)
      SUBROUTINE lay_order(i,j,i_modf)
      use M_imag
      use M_profile
      INTEGER :: h,i_modf 
!Check for proper layer ordering

!!      IF (updown_lay=='U') THEN

!       IF (layers(i,j,i_modf)<dm_lay(i,j,i_modf).OR.
!     *     ABS(layers(i,j,i_modf)-dm_lay(i,j,i_modf))<1) THEN
!        dm_lay(i,j,i_modf)=layers(i,j,i_modf)
!       ENDIF
!Control hierarchy upwards
       DO h=0,i_modf-1
        IF (h==0) THEN
         IF (-E(i,j)*1e-3>dm_lay(i,j,i_modf)) THEN
          dm_lay(i,j,i_modf)=-E(i,j)*1e-3
         ENDIF
        ELSE 
         IF (abs(layers(i,j,h))-abs(dm_lay(i,j,i_modf))>0)
     *   THEN
          dm_lay(i,j,h)=dm_lay(i,j,i_modf)
         ENDIF
        ENDIF
       ENDDO
!!      ELSE
!       IF (layers(i,j,i_modf)>dm_lay(i,j,i_modf).OR.
!     *    ABS(layers(i,j,i_modf)-dm_lay(i,j,i_modf))<1 ) THEN
!        dm_lay(i,j,i_modf)=layers(i,j,i_modf)
!       ENDIF
!Control hierarchy downwards
       DO h=i_modf+1,N_lay
        IF (abs(layers(i,j,h))-abs(dm_lay(i,j,i_modf))<0) 
     *  THEN
         dm_lay(i,j,h)=dm_lay(i,j,i_modf)
        ENDIF
       ENDDO
!!      ENDIF
      END SUBROUTINE
