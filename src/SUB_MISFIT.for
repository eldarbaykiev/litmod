!SUB_MISFIT.for
!Subroutine to calculate the chi squared and standard deviation 
!of data + error bars   

      SUBROUTINE MISFIT(N_dat,mft,dat,err_bar,rms,xi_sq) 
      INTEGER N_dat
      REAL(8) Err(N_dat,2),mft(N_dat),dat(N_dat),err_bar(N_dat)
      REAL(8) E_av,MIS,rms,xi_sq
      LOGICAL STAT
! Read the name of the measured and calculated data files and the number
! of points     
      Err(:,1)=mft(:)
      Err(:,2)=err_bar(:)
!Compute the mean
      E_av=SUM(Err(:,1))/N_dat
!Compute the standard deviation
      xi_sq=0D0
      rms=0D0
      DO i=1,N_dat
!       STD=STD+(Err(i,1)-E_av)**2
       xi_sq=xi_sq+(Err(i,1)/Err(i,2))**2
       rms=rms+(Err(i,1)/Err(i,2))**2/dat(i)**2
!       STD_w=STD_w+EXP(-Err(i,2))*(Err(i,1)-E_av)**2
      ENDDO
!      STD=SQRT(STD/N_dat) 
!      STD_w=SQRT(STD_w/SUM(EXP(-Err(:,2))))
      rms=SQRT(rms/N_dat)
      xi_sq=xi_sq/N_dat
      MIS=rms
!      MIS=xi_sq
!      write(*,*)'Chi squared: ',xi_sq
!      write(*,*)'rms',rms 

      END SUBROUTINE
