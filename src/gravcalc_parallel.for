!!!=============================================================================
!!!!!!!!!! gravcalc_parallel
!!!=============================================================================

      program gravcalc_parallel

      use M_grid_parameters, only: N_x,N_y
      use M_columns, only: XB, YB, ZB,
     *    ZBG, ZB_g, ZBG_g,
     *    xp, yp,
     *    rho_w, rho_B, GR, GEO,
     *    rho_up, rho_down,
     *    FA, GEOID, BGA, Uxx, Uyy, Uzz, Uzx, Uzy, Uxy,
     *    Z_usec, BG, BGA,
     *    FA,GEOID,BGA,
     *    Uxx,Uyy,Uzz,Uzx,Uzy,Uxy,
     *    FA_cr, GEOID_cr, BGA_cr,
     *    Uxx_cr, Uyy_cr, Uzz_cr, Uzx_cr, Uzy_cr, Uxy_cr
      use M_layers, only: E, N_lay, N_lay_cr



      integer output_yn
      integer i_counter, i_a
      character(len=32) :: arg
      character(len=32) :: filename_grav_input
      character(len=32) :: filename_grav_input_ne
      character(len=32) :: filename_nlines
      character(len=64) :: filename_grid
      integer filename_len_grav_input
      character(len=255):: command_nlines
      real(4) pr_b
      integer n_elem

      real(8) lrs

!from SUB_COLUMNS
      integer i_lay
      integer litho
      integer crst
      real(8) GR_w, GEO_w
      real(8) U_sec(6), U_sec_w(6), U_sec_mnt(6)
      real(8) ZB_mnt(2), ZB_mnt_g(2), GR_mnt, GEO_mnt

      call cpu_time(time1)

!###############################################################################

      if (iargc().eq.0) then
        output_yn = 1
      else
        output_yn = 0
      end if

      if (output_yn.gt.0) then
        write(*, *)'EXTERNAL GRAVITY COMPUTATION'
        write(*, *)'   Eldar  Baykiev'
        write(*, *)'   Javier Fullea'
        write(*, *)'                    2018'
        write(*, *)''
      end if

      if (iargc().eq.0) then
        write(*,*)'Input filename ([filename].dat):'
        read(*,*)filename_grav_input
      else
        if (output_yn.gt.0) then
          write(*, *)'Arguments (', iargc(), ')'
          do i = 1, iargc()
            call getarg(i, arg)
            write (*,*) i, arg
          end do
        end if

        call getarg(1, filename_grav_input)
      end if


!load parameters
      open(99, file='grav_params.dat', status='UNKNOWN')
      read(99, *) N_x, N_y, d_x, d_y, N_lay, N_lay_cr, rho_w, Z_Usec
      close(99)

      if (output_yn.gt.0) then
        write(*, *) 'Parameters (grav_params.dat):'
        write(*, *) 'N_x N_y d_x d_y N_lay N_lay_cr rho_w Z_Usec'
        write(*, *) N_x, N_y, d_x, d_y, N_lay, N_lay_cr, rho_w, Z_Usec
      end if

!load the elevation array
      open(100, file='grav_E_array.dat', status='UNKNOWN')
      allocate (E(N_x, N_y))
      do i_counter = 1, N_x
        read(100, *) E(i_counter,:)
      end do

      close(100)

!count number of lines in the input file
      command_nlines = "sed -n '$=' " //  filename_grav_input
      command_nlines = trim(command_nlines) // '>'

      filename_len_grav_input=len(trim(filename_grav_input))
!     write(*,*)filename_len_grav_input
!     write(*,*)filename_grav_input(1:filename_len_grav_input-4)
      filename_grav_input_ne =
     *      trim(filename_grav_input(1:filename_len_grav_input-4))
      filename_nlines =
     *     trim(filename_grav_input_ne) //
     *     "_nlines.dat"

      command_nlines = trim(command_nlines) // filename_nlines
      call system(command_nlines)
!     call system('awk 'END { print NR }' grav_input.dump > grav_calc_nlines.dat')

      open(101, file=filename_nlines,status='UNKNOWN')
      read(101, *) n_elem
      close(101)

      if (output_yn.gt.0) then
        write(*, *)'Number of the elements: ', n_elem
      end if


!allocate grids
      allocate (BGA(N_x,N_y),FA(N_x,N_y),GEOID(N_x,N_y))
      allocate (Uxx(N_x,N_y),Uyy(N_x,N_y),Uzz(N_x,N_y),Uzx(N_x,N_y),
     *          Uzy(N_x,N_y),Uxy(N_x,N_y))
      allocate (GEOID_cr(N_x,N_y),FA_cr(N_x,N_y),BGA_cr(N_x,N_y))
      allocate (Uxx_cr(N_x,N_y),Uyy_cr(N_x,N_y),Uzz_cr(N_x,N_y),
     *          Uzx_cr(N_x,N_y),Uzy_cr(N_x,N_y),Uxy_cr(N_x,N_y))

!loop over all stations
      if (output_yn.gt.0) then
        write(*,*)'Starting loop over all elements...'
      end if

1000  format(I3,A)
!read elements
      open(102, file=filename_grav_input,status='UNKNOWN')
!     write(*, *) 'i_lay litho rho_B XB YB ZB ZBG RHO_up RHO_down l'
      do i_counter=1,n_elem
        read(102, *)i_lay,litho,crst,rho_B,XB,YB,
     *              ZB,ZBG,RHO_up,RHO_down,lrs
!       write(*, *)i_lay,litho,rho_B,XB,YB,ZB,ZBG,RHO_up,RHO_down,lrs

        pr_b = i_counter/(n_elem/100.0)
        if ((pr_b-floor(pr_b)).le.0.001) then
          if (mod(int(pr_b), 5).eq.0) then
            if (output_yn.gt.0) then
              write(*,1000) int(pr_b), '%'
            end if
          end if
        end if

!from SUB_COLUMNS
        do I_xx=1,N_x
          XP=d_x*(I_xx-1)
          if(ZB(2)>ZB(1)) exit

          do I_yy=1,N_y
            YP=d_y*(I_yy-1)
            ZB_g(:)=ZB(:)-E(I_xx,I_yy)
            ZBG_g=ZBG-E(I_xx,I_yy)
            GR=0
            GEO=0
            GR_w=0
            GEO_w=0
            if (E(I_xx,I_yy)>0) then
                !Calculation point E>0
              call GRAV_GRAD3D(XB,YB,ZB,RHO_up,RHO_down,XP,YP,
     *                           E(I_xx,I_yy),GR)
              FA(I_xx,I_yy)=FA(I_xx,I_yy)+GR
              call GEO_GRAD3D(XB,YB,ZB_g,RHO_up,RHO_down,XP,YP,
     *                        0D0,GEO)

              GEOID(I_xx,I_yy)=GEOID(I_xx,I_yy)+GEO
              call U_SECOND_DER(XB,YB,ZB,RHO_up,RHO_down,XP,YP,
     *                        Z_usec,U_sec)
              Uzz(I_xx,I_yy)=Uzz(I_xx,I_yy)+U_sec(1)
              Uxx(I_xx,I_yy)=Uxx(I_xx,I_yy)+U_sec(2)
              Uyy(I_xx,I_yy)=Uyy(I_xx,I_yy)+U_sec(3)
              Uzx(I_xx,I_yy)=Uzx(I_xx,I_yy)+U_sec(4)
              Uzy(I_xx,I_yy)=Uzy(I_xx,I_yy)+U_sec(5)
              Uxy(I_xx,I_yy)=Uxy(I_xx,I_yy)+U_sec(6)

              if (i_lay==0) then
                if (ZBG(2)<0) then
                  !Effect of water
                  call GRAV_GRAD3D(XB,YB,ZBG,rho_w,rho_w,XP,YP,
     *                               E(I_xx,I_yy),GR_w)
                  FA(I_xx,I_yy)=FA(I_xx,I_yy)+GR_w
                  call GEO_GRAD3D(XB,YB,ZBG_g,rho_w,rho_w,XP,YP,0D0,
     *                              GEO_w)
                  GEOID(I_xx,I_yy)=GEOID(I_xx,I_yy)+GEO_w
                  call U_SECOND_DER(XB,YB,ZBG,rho_w,rho_w,
     *                                XP,YP,Z_usec,U_sec_w)
                  Uzz(I_xx,I_yy)=Uzz(I_xx,I_yy)+U_sec_w(1)
                  Uxx(I_xx,I_yy)=Uxx(I_xx,I_yy)+U_sec_w(2)
                  Uyy(I_xx,I_yy)=Uyy(I_xx,I_yy)+U_sec_w(3)
                  Uzx(I_xx,I_yy)=Uzx(I_xx,I_yy)+U_sec_w(4)
                  Uzy(I_xx,I_yy)=Uzy(I_xx,I_yy)+U_sec_w(5)
                  Uxy(I_xx,I_yy)=Uxy(I_xx,I_yy)+U_sec_w(6)
                end if
                call GRAV_GRAD3D(XB,YB,ZBG,RHO_B,RHO_B,XP,
     *                             YP,E(I_xx,I_yy),BG)
                BGA(I_xx,I_yy)=BGA(I_xx,I_yy)+BG
              end if
              if (isnan(U_sec(2))) stop

            else
              !Calculation point E<0
!             write(*,*)'before GEOGRAV_GRAD3D'
!             write(*,*)XB,YB,ZB,RHO_up,RHO_down,XP,YP,0D0
              call GEOGRAV_GRAD3D(XB,YB,ZB,RHO_up,RHO_down,XP,YP,0D0,
     *                                GEO,GR)
              FA(I_xx,I_yy)=FA(I_xx,I_yy)+GR
              GEOID(I_xx,I_yy)=GEOID(I_xx,I_yy)+GEO
              call U_SECOND_DER(XB,YB,ZB,RHO_up,RHO_down,XP,YP,
     *                              Z_usec,U_sec)
              Uzz(I_xx,I_yy)=Uzz(I_xx,I_yy)+U_sec(1)
              Uxx(I_xx,I_yy)=Uxx(I_xx,I_yy)+U_sec(2)
              Uyy(I_xx,I_yy)=Uyy(I_xx,I_yy)+U_sec(3)
              Uzx(I_xx,I_yy)=Uzx(I_xx,I_yy)+U_sec(4)
              Uzy(I_xx,I_yy)=Uzy(I_xx,I_yy)+U_sec(5)
              Uxy(I_xx,I_yy)=Uxy(I_xx,I_yy)+U_sec(6)
              if (i_lay==0) then
                if (ZBG(2)<0) then
                  call GEOGRAV_GRAD3D(XB,YB,ZBG,rho_w,rho_w,XP,YP,0D0,
     *                                    GEO_w,GR_w)
                  FA(I_xx,I_yy)=FA(I_xx,I_yy)+GR_w
                  GEOID(I_xx,I_yy)=GEOID(I_xx,I_yy)+GEO_w
                  call U_SECOND_DER(XB,YB,ZBG,rho_w,rho_w,XP,YP,
     *                                  Z_usec,U_sec_w)
                  Uzz(I_xx,I_yy)=Uzz(I_xx,I_yy)+U_sec_w(1)
                  Uxx(I_xx,I_yy)=Uxx(I_xx,I_yy)+U_sec_w(2)
                  Uyy(I_xx,I_yy)=Uyy(I_xx,I_yy)+U_sec_w(3)
                  Uzx(I_xx,I_yy)=Uzx(I_xx,I_yy)+U_sec_w(4)
                  Uzy(I_xx,I_yy)=Uzy(I_xx,I_yy)+U_sec_w(5)
                  Uxy(I_xx,I_yy)=Uxy(I_xx,I_yy)+U_sec_w(6)
                end if
                call GRAV_GRAD3D(XB,YB,ZBG,RHO_B,RHO_B,XP,
     *                               YP,0D0,BG)
                BGA(I_xx,I_yy)=BGA(I_xx,I_yy)+BG
              end if
            endif

            !write crust/mantle output !STATOIL$$$$$$$$$$$$$$$$$$$$$$$$4
            if(crst==1.and.i_lay<N_lay) then
              !Add homogenous mantle to compute the crustal gravity effect
              if (i_lay==N_lay_cr-1) then
                ZB_mnt(1)=lrs
                ZB_mnt(2)=-176.D3
                rho_up_mnt=3257.8D0
                rho_down_mnt=3257.8D0
                ZB_mnt_g(:)=ZB_mnt(:)-E(I_xx,I_yy)
                if (E(I_xx,I_yy)>0) then
                  call GRAV_GRAD3D(XB,YB,ZB_mnt,RHO_up_mnt,
     *                     RHO_down_mnt,XP,YP,E(I_xx,I_yy),GR_mnt)
                  call GEO_GRAD3D(XB,YB,ZB_mnt_g,RHO_up_mnt,
     *                     RHO_down_mnt,XP,YP,0D0,GEO_mnt)
                else
                  call GEOGRAV_GRAD3D(XB,YB,ZB_mnt,RHO_up_mnt,
     *                     RHO_down_mnt,XP,YP,0D0,GEO_mnt,
     *                     GR_mnt)
                end if
                call U_SECOND_DER(XB,YB,ZB_mnt,RHO_up_mnt,
     *                   RHO_down_mnt,XP,YP,Z_usec,U_sec_mnt)

                FA_cr(I_xx,I_yy)=FA_cr(I_xx,I_yy)+GR_mnt
                GEOID_cr(I_xx,I_yy)=GEOID_cr(I_xx,I_yy)+GEO_mnt
                Uzz_cr(I_xx,I_yy)=Uzz_cr(I_xx,I_yy)+U_sec_mnt(1)
                Uxx_cr(I_xx,I_yy)=Uxx_cr(I_xx,I_yy)+U_sec_mnt(2)
                Uyy_cr(I_xx,I_yy)=Uyy_cr(I_xx,I_yy)+U_sec_mnt(3)
                Uzx_cr(I_xx,I_yy)=Uzx_cr(I_xx,I_yy)+U_sec_mnt(4)
                Uzy_cr(I_xx,I_yy)=Uzy_cr(I_xx,I_yy)+U_sec_mnt(5)
                Uxy_cr(I_xx,I_yy)=Uxy_cr(I_xx,I_yy)+U_sec_mnt(6)
              end if

              FA_cr(I_xx,I_yy)=FA_cr(I_xx,I_yy)+GR+GR_w
              GEOID_cr(I_xx,I_yy)=GEOID_cr(I_xx,I_yy)+GEO+GEO_w
              Uzz_cr(I_xx,I_yy)=Uzz_cr(I_xx,I_yy)+U_sec(1)+U_sec_w(1)
              Uxx_cr(I_xx,I_yy)=Uxx_cr(I_xx,I_yy)+U_sec(2)+U_sec_w(2)
              Uyy_cr(I_xx,I_yy)=Uyy_cr(I_xx,I_yy)+U_sec(3)+U_sec_w(3)
              Uzx_cr(I_xx,I_yy)=Uzx_cr(I_xx,I_yy)+U_sec(4)+U_sec_w(4)
              Uzy_cr(I_xx,I_yy)=Uzy_cr(I_xx,I_yy)+U_sec(5)+U_sec_w(5)
              Uxy_cr(I_xx,I_yy)=Uxy_cr(I_xx,I_yy)+U_sec(6)+U_sec_w(6)

            end if
          end do
        end do
      end do

      close(102)

      if (output_yn.gt.0) then
        write(*,*) '    done'
      end if
!loop finished

!deallocate elevation
      if (allocated (E)) deallocate (E)

!save grids to files
      if (output_yn.gt.0) then
        write(*,*)'Saving grids to files...'
      end if


      filename_grid = trim(filename_grav_input_ne) // '_BGA.dat'
      open(200,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(200,*) BGA(i_counter,:)
      end do
      close(200)
      if (allocated (BGA)) deallocate (BGA)

      filename_grid = trim(filename_grav_input_ne) // '_FA.dat'
      open(201,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(201,*) FA(i_counter,:)
      end do
      close(201)
      if (allocated (FA)) deallocate (FA)

      filename_grid = trim(filename_grav_input_ne) // '_GEOID.dat'
      open(202,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(202,*) GEOID(i_counter,:)
      end do
      close(202)
      if (allocated (GEOID)) deallocate (GEOID)

      filename_grid = trim(filename_grav_input_ne) // '_Uxx.dat'
      open(203,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(203,*) Uxx(i_counter,:)
      end do
      close(203)
      if (allocated (Uxx)) deallocate (Uxx)

      filename_grid = trim(filename_grav_input_ne) // '_Uyy.dat'
      open(204,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(204,*) Uyy(i_counter,:)
      end do
      close(204)
      if (allocated (Uyy)) deallocate (Uyy)

      filename_grid = trim(filename_grav_input_ne) // '_Uzz.dat'
      open(205,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(205,*) Uzz(i_counter,:)
      end do
      close(205)
      if (allocated (Uzz)) deallocate (Uzz)

      filename_grid = trim(filename_grav_input_ne) // '_Uzx.dat'
      open(206,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(206,*) Uzx(i_counter,:)
      end do
      close(206)
      if (allocated (Uzx)) deallocate (Uzx)

      filename_grid = trim(filename_grav_input_ne) // '_Uzy.dat'
      open(207,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(207,*) Uzy(i_counter,:)
      end do
      close(207)
      if (allocated (Uzy)) deallocate (Uzy)

      filename_grid = trim(filename_grav_input_ne) // '_Uxy.dat'
      open(208,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(208,*) Uxy(i_counter,:)
      end do
      close(208)
      if (allocated (Uxy)) deallocate (Uxy)

      filename_grid = trim(filename_grav_input_ne) // '_BGA_cr.dat'
      open(300,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(300,*) BGA_cr(i_counter,:)
      end do
      close(300)
      if (allocated (BGA_cr)) deallocate (BGA_cr)

      filename_grid = trim(filename_grav_input_ne) // '_FA_cr.dat'
      open(301,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(301,*) FA_cr(i_counter,:)
      end do
      close(301)
      if (allocated (FA_cr)) deallocate (FA_cr)

      filename_grid = trim(filename_grav_input_ne) // '_GEOID_cr.dat'
      open(302,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(302,*) GEOID_cr(i_counter,:)
      end do
      close(302)
      if (allocated (GEOID_cr)) deallocate (GEOID_cr)

      filename_grid = trim(filename_grav_input_ne) // '_Uxx_cr.dat'
      open(303,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(303,*) Uxx_cr(i_counter,:)
      end do
      close(303)
      if (allocated (Uxx_cr)) deallocate (Uxx_cr)

      filename_grid = trim(filename_grav_input_ne) // '_Uyy_cr.dat'
      open(304,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(304,*) Uyy_cr(i_counter,:)
      end do
      close(304)
      if (allocated (Uyy_cr)) deallocate (Uyy_cr)

      filename_grid = trim(filename_grav_input_ne) // '_Uzz_cr.dat'
      open(305,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(305,*) Uzz_cr(i_counter,:)
      end do
      close(305)
      if (allocated (Uzz_cr)) deallocate (Uzz_cr)

      filename_grid = trim(filename_grav_input_ne) // '_Uzx_cr.dat'
      open(306,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(306,*) Uzx_cr(i_counter,:)
      end do
      close(306)
      if (allocated (Uzx_cr)) deallocate (Uzx_cr)

      filename_grid = trim(filename_grav_input_ne) // '_Uzy_cr.dat'
      open(307,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(307,*) Uzy_cr(i_counter,:)
      end do
      close(307)
      if (allocated (Uzy_cr)) deallocate (Uzy_cr)

      filename_grid = trim(filename_grav_input_ne) // '_Uxy_cr.dat'
      open(308,file=filename_grid,status='UNKNOWN')
      do i_counter=1,N_x
        write(308,*) Uxy_cr(i_counter,:)
      end do
      close(308)
      if (allocated (Uxy_cr)) deallocate (Uxy_cr)

      call cpu_time(time2)

      if (output_yn.gt.0) then
        write(*,*) '    done'
        write(*,*) 'program finished in ',time2-time1, ' sec'
      end if

      stop
      end
