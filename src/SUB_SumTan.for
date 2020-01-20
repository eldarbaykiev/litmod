CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SUMTAN(X,Y,SUM)
C     =================
C
C***********************************************************************
C
C Subroutine prepares the coefficients of an two ATAN to be summed in
C order to calculate only once the ATAN. For subtraction:
C ATAN(X)-ATAN(Y) = ATAN(X)+ATAN(-Y) !!
C
C***********************************************************************
C
C ATAN(X)+ATAN(Y)=ATAN((X+Y)/(1-XY))+/-PI
C
C The result (X+Y)/(1-XY) is stored in X
C
C SUM is the multiple of PI which has to be added
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      IMPLICIT REAL*8 (A-H,O-Z)
      DATA PI2/1.5707963267949D0/,PI/3.1415926535898D0/
      SAVE
      F1=X+Y
      F2=X*Y
      IF(F2 .EQ. 1.) THEN
         IF(X .GT. 0.) THEN
            SUM=SUM+PI2
         ELSE
            SUM=SUM-PI2
         ENDIF
         X=0.
      ELSE
         IF(F2 .GT. 1.) THEN
            IF(X .GT. 0.) THEN
               SUM=SUM+PI
            ELSE
               SUM=SUM-PI
            ENDIF
         ENDIF
         X=F1/(1.-F2)
      ENDIF
      RETURN
      END
