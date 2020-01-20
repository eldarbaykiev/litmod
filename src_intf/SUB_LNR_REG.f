! SUB_LNR_REG.f


       SUBROUTINE LNR_REG(X,Y,n,a,b)
       REAL(4) X(n),Y(n),a,b,Sxy,Sx2
       INTEGER n

       X_m=SUM(X)/n
       Y_m=SUM(Y)/n
       DO i=1,n
        Sx2=SUM(X**2)/n-X_m**2
        Sxy=SUM(X*Y)/n-X_m*Y_m
       ENDDO
       a=Sxy/Sx2
       b=Y_m-a*X_m 


      END SUBROUTINE
