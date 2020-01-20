

      subroutine Diag_x(n_kx,kx,kx_n,diagonal_x)
      use M_grid_parameters
      real(8) n_kx,kx,kx_n
      real(8) diagonal_x
      diagonal_x=-c_x*(n_kx+2*kx+kx_n)
      end subroutine

      subroutine Diag_y(n_ky,ky,ky_n,diagonal_y)
      use M_grid_parameters
      real(8) n_ky,ky,ky_n
      real(8) diagonal_y
      diagonal_y=-c_y*(n_ky+2*ky+ky_n)
      end subroutine

      subroutine Diag_z(n_kz,kz,kz_n,diagonal_z)
      use M_grid_parameters
      real(8) n_kz,kz,kz_n
      real(8) diagonal_z
      diagonal_z=-c_z*(n_kz+2*kz+kz_n)
      end subroutine


      subroutine n_Diag_x(n_kx,kx,diagonal_x)
      use M_grid_parameters
      real(8) n_kx,kx
      real(8) diagonal_x
      diagonal_x=-c_x*(n_kx+kx)
      end subroutine

      subroutine n_Diag_y(n_ky,ky,diagonal_y)
      use M_grid_parameters
      real(8) n_ky,ky
      real(8) diagonal_y
      diagonal_y=-c_y*(n_ky+ky)
      end subroutine


      subroutine UpperD_x(kx,kx_n,UD_x)
      use M_grid_parameters
      real(8) kx,kx_n
      real(8) UD_x
      UD_x=c_x*(kx+kx_n)
      end subroutine

      subroutine UpperD_y(ky,ky_n,UD_y)
      use M_grid_parameters
      real(8) ky,ky_n
      real(8) UD_y
      UD_y=c_y*(ky+ky_n)
      end subroutine

      subroutine UpperD_z(kz,kz_n,UD_z)
      use M_grid_parameters
      real(8) kz,kz_n
      real(8) UD_z
      UD_z=c_z*(kz+kz_n)
      end subroutine

      SUBROUTINE LowerD_x(kx,kx_n,LD_x)
      USE M_grid_parameters
      real(8) kx,kx_n
      real(8) LD_x
      LD_x=c_x*(kx+kx_n)
      END SUBROUTINE

      SUBROUTINE LowerD_y(ky,ky_n,LD_y)
      USE M_grid_parameters
      real(8) ky,ky_n
      real(8) LD_y
      LD_y=c_y*(ky+ky_n)
      END SUBROUTINE

      SUBROUTINE LowerD_z(kz,kz_n,LD_z)
      USE M_grid_parameters
      real(8) kz,kz_n
      real(8) LD_z
      LD_z=c_z*(kz+kz_n)
      END SUBROUTINE
