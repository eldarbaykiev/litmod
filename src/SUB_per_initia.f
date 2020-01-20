

      subroutine initia(numcom)

      use m_pervar

      ! Assign name to Perple_X project.
      call num2nam(numcom)

c                                 iam is a flag indicating the Perple_X program
      iam = 2
! c                                 version info
      ! call vrsion (6)
c                                 initialization, read files etc.
      rxn = .false.

      call iniprp

      end subroutine








      subroutine initia_II(numcom)

      use m_pervar_II

      ! Assign name to Perple_X project.
      call num2nam(numcom)

c                                 iam is a flag indicating the Perple_X program
      iam = 2
! c                                 version info
      ! call vrsion (6)
c                                 initialization, read files etc.
      rxn = .false.

      call iniprp

      end subroutine








      subroutine initia_III(numcom)

      use m_pervar_III

      ! Assign name to Perple_X project.
      call num2nam(numcom)

c                                 iam is a flag indicating the Perple_X program
      iam = 2
! c                                 version info
      ! call vrsion (6)
c                                 initialization, read files etc.
      rxn = .false.

      call iniprp

      end subroutine
