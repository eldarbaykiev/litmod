CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE GEOGRAV_GRAD3D(XB,YB,ZB,RHO_UP,RHO_LOW,XP,YP,ZP,GEO,GR)
C     =========================
C
C***********************************************************************
C
C Subroutine calculates the geoid and gravity effects of a rectangular
C prism with its sides parallel to the coordinate system and a density
C that may present a vertical linear gradient. Geoid formula derived by
C Hermann Zeyen, gravity by Gallardo-Delgado et al., Geophysics, 68
C (2002), p.949.
C
C Condition for using this subroutine: The topography has to at or
C below sea level, i.e. ZP=0.
C
C***********************************************************************
C
C XB,YB,ZB: Coordinates of the prisms
C           Each of them is a vector of two elements ordered by
C           increasing coordinate.
C RHO_UP, RHO_LOW : Densities at ZB(1) and ZB(2) respectively (kg/m3)
C XP,YP,ZP: Position of the observation point (X=E-W, Y=N-S)
C GEO:      Geoid effect (m)
C GR:       Gravity effect (mGal)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XB(2),YB(2),ZB(2),X(2),Y(2),Z(2)
      PARAMETER (G=6.67D-11)
      PARAMETER (GEO_FAC=G/9.81D-9)
      PARAMETER (GR_FAC=G*1.D11)
C      PARAMETER (GEO_FAC=G/9.81)
C      PARAMETER (GR_FAC=G*1.D5)
!      write(*,*)'GeoGrav_Grad3D'
!      write(*,*)'RHO_UP,RHO_LOW',RHO_UP,RHO_LOW
!      pause 
      IF(ZP .NE. 0.D0) THEN
         PRINT '('' Subroutine GeoGrav_Grad3D called with a non-zero '',
     *           ''topography!''/'' Program stopped.'')'
         STOP
      ENDIF

      IF (ZB(1)==ZB(2)) THEN
       GEO=0D0
       GR=0D0 
       RETURN
      ENDIF

      GEO=0.D0
      GR=0.D0
      GAMMA=(RHO_LOW-RHO_UP)/(ZB(2)-ZB(1))
      G6=GAMMA/6.D0
      RHO0=(RHO_UP-GAMMA*ZB(1))*1.D-3
      RHO_Z0=RHO0+GAMMA*ZP*1.D-3
C      RHO0=RHO_UP-GAMMA*ZB(1)
C      RHO_Z0=RHO0+GAMMA*ZP
      X(1)=(XB(1)-XP)*1.D-3
      X(2)=(XB(2)-XP)*1.D-3
      Y(1)=(YB(1)-YP)*1.D-3
      Y(2)=(YB(2)-YP)*1.D-3
      Z(1)=(ZB(1)-ZP)*1.D-3
      Z(2)=(ZB(2)-ZP)*1.D-3
C      X(1)=XB(1)-XP
C      X(2)=XB(2)-XP
C      Y(1)=YB(1)-YP
C      Y(2)=YB(2)-YP
C      Z(1)=ZB(1)-ZP
C      Z(2)=ZB(2)-ZP
      X12=X(1)*X(1)
      X22=X(2)*X(2)
      Y12=Y(1)*Y(1)
      Y22=Y(2)*Y(2)
      Z12=Z(1)*Z(1)
      Z123=3.D0*Z12
      Z22=Z(2)*Z(2)
      Z223=3.D0*Z22
      R111=DSQRT(X12+Y12+Z12)
      R112=DSQRT(X12+Y12+Z22)
      R121=DSQRT(X12+Y22+Z12)
      R122=DSQRT(X12+Y22+Z22)
      R211=DSQRT(X22+Y12+Z12)
      R212=DSQRT(X22+Y12+Z22)
      R221=DSQRT(X22+Y22+Z12)
      R222=DSQRT(X22+Y22+Z22)
C
      IF(DABS(X(1)).LT.1.D-10 .OR. DABS(Y(1)).LT.1.D-10) THEN
         DX1Y1=0.D0
      ELSE
         DX1Y1=X(1)*Y(1)*DLOG((Z(2)+R112)/(Z(1)+R111))
      ENDIF
      IF(DABS(X(2)).LT.1.D-10 .OR. DABS(Y(1)).LT.1.D-10) THEN
         DX2Y1=0.D0
      ELSE
         DX2Y1=X(2)*Y(1)*DLOG((Z(1)+R211)/(Z(2)+R212))
      ENDIF
      IF(DABS(X(1)).LT.1.D-10 .OR. DABS(Y(2)).LT.1.D-10) THEN
         DX1Y2=0.D0
      ELSE
         DX1Y2=X(1)*Y(2)*DLOG((Z(1)+R121)/(Z(2)+R122))
      ENDIF
      IF(DABS(X(2)).LT.1.D-10 .OR. DABS(Y(2)).LT.1.D-10) THEN
         DX2Y2=0.D0
      ELSE
         DX2Y2=X(2)*Y(2)*DLOG((Z(2)+R222)/(Z(1)+R221))
      ENDIF
      IF(DABS(X(1)) .LT. 1.D-10) THEN
         DX1Z1=0.D0
         DX1Z2=0.D0
      ELSE
         DX1Z1=X(1)*DLOG((Y(2)+R121)/(Y(1)+R111))
         DX1Z2=X(1)*DLOG((Y(1)+R112)/(Y(2)+R122))
      ENDIF
      IF(DABS(X(2)) .LT. 1.D-10) THEN
         DX2Z1=0.D0
         DX2Z2=0.D0
      ELSE
         DX2Z1=X(2)*DLOG((Y(1)+R211)/(Y(2)+R221))
         DX2Z2=X(2)*DLOG((Y(2)+R222)/(Y(1)+R212))
      ENDIF
      IF(DABS(Y(1)) .LT. 1.D-10) THEN
         DY1Z1=0.D0
         DY1Z2=0.D0
      ELSE
         DY1Z1=Y(1)*DLOG((X(2)+R211)/(X(1)+R111))
         DY1Z2=Y(1)*DLOG((X(1)+R112)/(X(2)+R212))
      ENDIF
      IF(DABS(Y(2)) .LT. 1.D-10) THEN
         DY2Z1=0.D0
         DY2Z2=0.D0
      ELSE
         DY2Z1=Y(2)*DLOG((X(1)+R121)/(X(2)+R221))
         DY2Z2=Y(2)*DLOG((X(2)+R222)/(X(1)+R122))
      ENDIF
      SLOG_GEO=RHO_Z0*(DX1Y1+DX2Y2+DX1Y2+DX2Y1+
     *                 Z(1)*(DX1Z1+DX2Z1+DY1Z1+DY2Z1)+
     *                 Z(2)*(DX2Z2+DX1Z2+DY2Z2+DY1Z2))+
     *         G6*((X12+Z123)*DX1Z1+(X12+Z223)*DX1Z2+
     *             (X22+Z223)*DX2Z2+(X22+Z123)*DX2Z1+
     *             (Y12+Z123)*DY1Z1+(Y12+Z223)*DY1Z2+
     *             (Y22+Z223)*DY2Z2+(Y22+Z123)*DY2Z1)
      SLOG_GRA=-RHO_Z0*(DX1Z2+DX2Z1+DY1Z2+DY2Z1+
     *                 DX1Z1+DX2Z2+DY1Z1+DY2Z2)+
     *         GAMMA*(DX1Y1+DX2Y2+DX1Y2+DX2Y1)
C
C Calculate ATAN terms which are multiplied with X1.
C
      SUMPI=0.D0
      IF(DABS(X(1)) .LT. 1.D-10) THEN
         SUMX1=0.D0
      ELSE
         SUMPI1=0.D0
         T1=Y(1)*Z(2)/(X(1)*R112)
         T2=Y(2)*Z(1)/(X(1)*R121)
         CALL SUMTAN(T1,T2,SUMPI1)
         SUMPI2=0.D0
         T3=Y(1)*Z(1)/(X(1)*R111)
         T4=Y(2)*Z(2)/(X(1)*R122)
         CALL SUMTAN(T3,T4,SUMPI2)
         CALL SUMTAN(T1,-T3,SUMPI)
         SUMX1=DATAN(T1)+SUMPI+SUMPI1-SUMPI2
      ENDIF
C
C Calculate ATAN terms which are multiplied with X2.
C
      SUMPI=0.D0
      IF(DABS(X(2)) .LT. 1.D-10) THEN
         SUMX2=0.D0
      ELSE
         SUMPI2=0.D0
         T3=Y(1)*Z(2)/(X(2)*R212)
         T4=Y(2)*Z(1)/(X(2)*R221)
         CALL SUMTAN(T3,T4,SUMPI2)
         SUMPI1=0.D0
         T1=Y(1)*Z(1)/(X(2)*R211)
         T2=Y(2)*Z(2)/(X(2)*R222)
         CALL SUMTAN(T1,T2,SUMPI1)
         CALL SUMTAN(T1,-T3,SUMPI)
         SUMX2=DATAN(T1)+SUMPI+SUMPI1-SUMPI2
      ENDIF
C
C Calculate ATAN terms which are multiplied with Y1.
C
      SUMPI=0.D0
      IF(DABS(Y(1)) .LT. 1.D-10) THEN
         SUMY1=0.D0
      ELSE
         SUMPI1=0.D0
         T1=X(1)*Z(2)/(Y(1)*R112)
         T2=X(2)*Z(1)/(Y(1)*R211)
         CALL SUMTAN(T1,T2,SUMPI1)
         SUMPI2=0.D0
         T3=X(1)*Z(1)/(Y(1)*R111)
         T4=X(2)*Z(2)/(Y(1)*R212)
         CALL SUMTAN(T3,T4,SUMPI2)
         CALL SUMTAN(T1,-T3,SUMPI)
         SUMY1=DATAN(T1)+SUMPI+SUMPI1-SUMPI2
      ENDIF
C
C Calculate ATAN terms which are multiplied with Y2.
C
      SUMPI=0.D0
      IF(DABS(Y(2)) .LT. 1.D-10) THEN
         SUMY2=0.D0
      ELSE
         SUMPI2=0.D0
         T3=X(1)*Z(2)/(Y(2)*R122)
         T4=X(2)*Z(1)/(Y(2)*R221)
         CALL SUMTAN(T3,T4,SUMPI2)
         SUMPI1=0.D0
         T1=X(1)*Z(1)/(Y(2)*R121)
         T2=X(2)*Z(2)/(Y(2)*R222)
         CALL SUMTAN(T1,T2,SUMPI1)
         CALL SUMTAN(T1,-T3,SUMPI)
         SUMY2=DATAN(T1)+SUMPI+SUMPI1-SUMPI2
      ENDIF
C
C Calculate ATAN terms which are multiplied with Z1.
C
      SUMPI=0.D0
      IF(DABS(Z(1)) .LT. 1.D-10) THEN
         SUMZ1=0.D0
      ELSE
         SUMPI1=0.D0
         T1=X(1)*Y(2)/(Z(1)*R121)
         T2=X(2)*Y(1)/(Z(1)*R211)
         CALL SUMTAN(T1,T2,SUMPI1)
         SUMPI2=0.D0
         T3=X(1)*Y(1)/(Z(1)*R111)
         T4=X(2)*Y(2)/(Z(1)*R221)
         CALL SUMTAN(T3,T4,SUMPI2)
         CALL SUMTAN(T1,-T3,SUMPI)
         SUMZ1=DATAN(T1)+SUMPI+SUMPI1-SUMPI2
      ENDIF
C
C Calculate ATAN terms which are multiplied with Z2.
C
      SUMPI=0.D0
      IF(DABS(Z(2)) .LT. 1.D-10) THEN
         SUMZ2=0.D0
      ELSE
         SUMPI2=0.D0
         T3=X(1)*Y(2)/(Z(2)*R122)
         T4=X(2)*Y(1)/(Z(2)*R212)
         CALL SUMTAN(T3,T4,SUMPI2)
         SUMPI1=0.D0
         T1=X(1)*Y(1)/(Z(2)*R112)
         T2=X(2)*Y(2)/(Z(2)*R222)
         CALL SUMTAN(T1,T2,SUMPI1)
         CALL SUMTAN(T1,-T3,SUMPI)
         SUMZ2=DATAN(T1)+SUMPI+SUMPI1-SUMPI2
      ENDIF
      STAN_GEO=RHO_Z0*(X12*SUMX1+X22*SUMX2+Y12*SUMY1+Y22*SUMY2)+
     *        (RHO_Z0+4.D0*G6*Z(1))*Z12*SUMZ1+
     *        (RHO_Z0+4.D0*G6*Z(2))*Z22*SUMZ2
      SUM_3=(X(1)*(Y(1)*(R112-R111)+Y(2)*(R121-R122))+
     *       X(2)*(Y(1)*(R211-R212)+Y(2)*(R222-R221)))*G6*2.D0
      GEO=GEO_FAC*(STAN_GEO*0.5-(SLOG_GEO+SUM_3))
      STAN_GRA=(RHO_Z0+GAMMA*Z(1)*0.5D0)*Z(1)*SUMZ1+
     *         (RHO_Z0+GAMMA*Z(2)*0.5D0)*Z(2)*SUMZ2-
     *         GAMMA*0.5D0*(X12*SUMX1+X22*SUMX2+Y12*SUMY1+Y22*SUMY2)
      GR=GR_FAC*(SLOG_GRA+STAN_GRA)
      END

