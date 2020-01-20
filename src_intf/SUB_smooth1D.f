! Subroutine  SUB_smooth1D.f


      SUBROUTINE smooth1D(A,I_n)
      INTEGER::I_n,I_win
      REAL(4):: A(I_n)
      I_win=3
      DO i=I_win+1,I_n-I_win
       A(i)=SUM(A(i-I_win:i+I_win))/(2*I_win+1)
      ENDDO
 



      END SUBROUTINE 
