CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE U_SECOND_DER(XB,YB,ZB,RHO_UP,RHO_LOW,XP,YP,ZP,U_sec)
C     ======================
C
C***********************************************************************
C
C Subroutine calculates the second derivatives of the gravity potential
C of a rectangular
C prism with its sides parallel to the coordinate system and a density
C that may present a vertical linear gradient. 
C See appendix in Fullea et al. 2015 JAG
C***********************************************************************
C
C XB,YB,ZB: Coordinates of the prisms
C           Each of them is a vector of two elements ordered by
C           increasing coordinate.
C RHO_UP, RHO_LOW : Densities at ZB(1) and ZB(2) respectively (kg/m3)
C XP,YP,ZP: Position of the observation point (X=E-W, Y=N-S)
C U_sec(6):     Second derivatives matrix (Eotvos) Uzz Uxx Uyy Uxz Uyz Uxy
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
!      IMPLICIT REAL*8 (A-H,O-Z)
      REAL(8) XB(2),YB(2),ZB(2),X(2),Y(2),Z(2),U_sec(6),v(6)
      REAL(8) r,GAMMA,G_uni,RHO_UP,RHO_LOW,XP,YP,ZP,rho_0
      INTEGER i,j,k,sgn
      LOGICAL NO_GRAD_FLAG
      PARAMETER (G_uni=6.67D-11)
!      PARAMETER (GRAV_FAC=6.67D0)
!      write(*,*)'GRAV_GRAD3D' 
!      write(*,*)'RHO_UP,RHO_LOW',RHO_UP,RHO_LOW
!      pause
      IF (ZB(1)==ZB(2)) THEN
       U_sec=0D0
       RETURN
      ENDIF
      X(1)=(XB(1)-XP)
      X(2)=(XB(2)-XP)
      Y(1)=(YB(1)-YP)
      Y(2)=(YB(2)-YP)
      Z(1)=(ZB(1)-ZP)
      Z(2)=(ZB(2)-ZP)
      IF(DABS(RHO_UP-RHO_LOW)<1D-3) THEN
         NO_GRAD_FLAG=.TRUE.
         GAMMA=0.D0
      ELSE
         NO_GRAD_FLAG=.FALSE.
         GAMMA=(RHO_LOW-RHO_UP)/(ZB(2)-ZB(1))
      ENDIF
!      write(*,*)'Z, ZB, ZP,gamma',Z, ZB, ZP,GAMMA
!      write(*,*)'X, XB, XP',X, XB, XP
!      write(*,*)'Y, YB, YP',Y, YB, YP
      U_sec=0D0
      v=0D0
      rho_0=RHO_UP-GAMMA*Z(1)
!      rho_0=(RHO_UP+RHO_LOW)/2d0
!      write(*,*)'rho_0',rho_0
      do i=1,2
       do j=1,2
        do k=1,2

        sgn=(-1)**(i+j+k)
        r=SQRT(x(i)**2+y(j)**2+z(k)**2)

        IF(r==0D0) THEN
         U_sec=0D0
         CYCLE
        ENDIF

        IF(x(i)==0D0) THEN
         U_sec(2)=0D0       
         U_sec(4)=-RHO_0*DLOG(y(j)+r)+GAMMA*(y(j)*DLOG(z(k)+r)) 
        ELSE
         U_sec(2)=RHO_0*DATAN(z(k)*y(j)/(x(i)*r))-
     *            GAMMA*x(i)*DLOG(y(j)+r)
         U_sec(4)=-RHO_0*LOG(y(j)+r)+
     *         GAMMA*(y(j)*DLOG(z(k)+r)-x(i)*DATAN(z(k)*y(j)/(x(i)*r)))
        ENDIF

        IF(y(j)==0D0) THEN
         U_sec(3)=0D0
         U_sec(5)=-RHO_0*DLOG(x(i)+r)+GAMMA*(x(i)*DLOG(z(k)+r))
        ELSE
         U_sec(3)=RHO_0*DATAN(z(k)*x(i)/(y(j)*r))-
     *            GAMMA*y(j)*DLOG(x(i)+r)
         U_sec(5)=-RHO_0*DLOG(x(i)+r)+
     *          GAMMA*(x(i)*DLOG(z(k)+r)-y(j)*DATAN(z(k)*x(i)/(y(j)*r)))
        ENDIF 

        IF(z(k)==0D0) THEN
         U_sec(1)=0D0
        ELSE
         U_sec(1)=DATAN(x(i)*y(j)/(z(k)*r))*RHO_0+
     *   GAMMA*(x(i)*DLOG(y(j)+r)+y(j)*DLOG(x(i)+r))
        ENDIF 

        U_sec(6)=-RHO_0*DLOG(z(k)+r)-GAMMA*r

        U_sec=U_sec*sgn
         
        v=v+U_sec
        enddo
       enddo
      enddo

      U_sec=v*G_uni*1D9
!      write(*,*)'Uzz,Uxx,Uyy,Uzx,Uzy',U_sec(1:5)
!      write(*,*)'Uxx+Uxx+Uzz',sum(U_sec(1:3))




      END SUBROUTINE

