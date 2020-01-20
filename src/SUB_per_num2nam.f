
      subroutine num2nam(inpnumcom)

      character*100 prject,tfname
		common/ cst228 /prject,tfname

		integer inpnumcom,numcom
		character(len=50)namcom


      ! Link the numcom to the namcom.
      ! namcom will be the name for the Perple_X project.
      open(24,file='mnt.info')
      do
        read(24,*,iostat=ios)numcom,namcom
        if(ios.eq.0)then
          if(numcom.eq.inpnumcom) prject = namcom
        else
          exit
        endif
      enddo
      close(24)


      end subroutine
