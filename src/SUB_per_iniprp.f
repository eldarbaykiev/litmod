      subroutine iniprp
c----------------------------------------------------------------------
c iniprp - read data files and initialization for meemum
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical first, output, err

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer iemod,kmod
      logical smod,pmod
      double precision emod
      common/ cst319 /emod(k15,k10),smod(h9),pmod(k10),iemod(k10),kmod
c-----------------------------------------------------------------------
      first = .true.
      output = .false.
      err = .false.
c                                 elastic modulii flag
      kmod = 0
c                                 -------------------------------------------
c                                 open statements for units n1-n5 and n9
c                                 are in subroutine input1
      call input1 (first,output,err)
c                                 for meemum turn auto_refine OFF
      iopt(6) = 0
c                                 read thermodynamic data on unit n2:
      call input2 (first)
c                                 allow reading of auto-refine data
      call setau1 (output)
c                                 read data for solution phases on n9:
      call input9 (first,output)
c                                 call initlp to initialize arrays
c                                 for optimization.
      call initlp

      end
