!Subroutine: SUB_STAT.for 


      SUBROUTINE STATT(V1,V2,n,AV,std_dev)
      REAL(8) V1(n),V2(n),Vdiff(n),AV,std_dev,dm
      INTEGER n  
      Vdiff(:)=V1(:)-V2(:)
!      write(*,*)n,minval(Vdiff),maxval(Vdiff)
      AV=SUM(Vdiff(:))/n
      dm=0
      DO i=1,n
       dm=dm+(Vdiff(i)-AV)**2
      ENDDO
      std_dev=SQRT(dm/n)

      END SUBROUTINE 
