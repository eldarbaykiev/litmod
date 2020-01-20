


      subroutine meemum4_sub1(numcom,pre,tem)


		use m_pervar


		integer numcom
		real(8) pre, tem


		! Assign name to Perple_X project.
		call num2nam(numcom)


c                                 iam is a flag indicating the Perple_X program
      iam = 2
! c                                 version info
      ! call vrsion (6)
c                                 initialization, read files etc. 
      rxn = .false.


      call iniprp


      ! write (*,1000) 
      ! read (*,'(a)') yes

      ! if (yes.eq.'y'.or.yes.eq.'Y') then 
! c                                 bulk is true, user enters composition and p-t conditions
         ! bulk = .true.

      ! else 
! c                                 else user enters only p-t and composition read from input file.
         bulk = .false.

      ! end if 

c                                 iwt is set by input, it is only used below to determine
c                                 whether to convert weight compositions to molar. the 
c                                 computations are done solely in molar units. 
      amount = 'molar '

      if (iwt.eq.1) amount = 'weight'

c                                 computational loop


		! ! Start the iterations.
		! do

			! ! Write the number of iterations in a file.
			! con = con+1
			! open(0001,file='con.xyz',status='UNKNOWN')
			! write(0001,*)con
			! close(0001)


			! ! read potential variable values    
			! ! v(1) is P(bar), v(2) is T(K) the pointer jv used 
         ! ! for general problems but can be eliminated for calculations 
			! ! simply as a f(P,T)       
         ! write (*,1070) (vname(jv(i)), i = 1, ipot)
         ! read (*,*,iostat=ier) (v(jv(i)), i = 1, ipot)
			v(jv(2)) = pre ! Pressure.
			v(jv(1)) = tem  ! Temperature.
			ier = 0


         ! if (ier.ne.0) cycle
         ! if (v(jv(1)).eq.0d0) exit 


         if (bulk) then 
c                                 load the composition into b, the component names are  
c                                 in cname, if iwt = 1 the composition is in mass fractions
c                                 otherwise in molar units.

				! do 
					! write (*,1060) amount
               ! write (*,'(12(a,1x))') (cname(i),i=1,jbulk)
               ! read (*,*,iostat=ier) (cblk(i),i=1,jbulk)
               ! if (ier.eq.0) exit
            ! end do  
         
            if (iwt.eq.1) then 
c                                 convert mass to molar 
               do i = 1, jbulk
                  cblk(i) = cblk(i)/atwt(i)
               end do 

            end if
c                                 normalize the composition vector, this 
c                                 is necessary for reasons of stupidity (lpopt0). 
            ctotal = 0d0

            do i = 1, icp
               ctotal = ctotal + cblk(i)
            end do 

            do i = 1, icp
               b(i) = cblk(i)/ctotal
            end do

         end if 
c                                 set dependent variables
         call incdp0
c                                 lpopt does the minimization and outputs
c                                 the results to the print file.
         call lpopt0 (idead)

         if (idead.eq.0) then
c                                 compute derivative properties
            call getloc (itri,jtri,ijpt,wt,nodata)
c                                 print summary to LUN 6
            call calpr0 (6)

            if (io3.eq.0) call calpr0 (n3)

         end if 

         if (goodc(1)+badc(1).gt.0d0) then
            num = badc(1)/(badc(1)+goodc(1))*1d2
            if (num.gt.1d-1) call warn (53,num,i,'MEEMUM')
            goodc(1) = 0d0
            badc(1) = 0d0 
         end if


		! ! End of the iterations.
      ! end do


      end
    


















