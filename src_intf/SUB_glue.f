! SUB_glue.f

      SUBROUTINE glue(glue_act,i_glue,j_glue)
      USE M_imag
      USE M_profile
      INTEGER lay,i_glue,j_glue
      REAL(4)::dz
      CHARACTER(1)::glue_act
!Minimum distance between layers(km)     
      dz=.5 
!CREATE THE GLUE_m MATRIX真真真真真真真真真真真真真真真真真真真  
      IF (glue_act=='C') THEN
       DO j=JDIM,1,-1
        DO i=1,IDIM
         DO k=1,N_lay
          DO m=k,N_lay
          IF (ABS(layers(i,j,k)-layers(i,j,m))<dz) THEN
           GLUE_m(i,j,k,m)=1     
          ELSE
           GLUE_m(i,j,k,m)=0
          ENDIF  
          GLUE_m(i,j,m,k)=GLUE_m(i,j,k,m)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF 
!CREATE THE GLUE_m MATRIX真真真真真真真真真真真真真真真真真真真
 
!USE THE GLUE_m MATRIX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF (glue_act=='U') THEN
       DO k=0,N_lay
!        write(*,*)'i_lay,k',i_lay,k,GLUE_m(i_glue,j_glue,i_lay,k)
!        pause
        IF (k==0) THEN
         IF (ABS(-E(i_glue,j_glue)*1e-3-layers(i_glue,j_glue,i_lay))<dz)
     *   THEN
          layers(i_glue,j_glue,i_lay)=-E(i_glue,j_glue)*1e-3 
         ENDIF
        ELSE 
         IF (GLUE_m(i_glue,j_glue,i_lay,k)==1) THEN
          layers(i_glue,j_glue,k)=layers(i_glue,j_glue,i_lay)
         ENDIF
        ENDIF
       ENDDO 
      ENDIF
!USE THE GLUE_m MATRIX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END SUBROUTINE
