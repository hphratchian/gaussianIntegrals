      module memory_utils
!
!     This is a simple utility module that uses c intrinsics to figure out how
!     much memory is currently being used by a fortran program.
!
!     H. P. Hratchian, 2025.
!
      use iso_c_binding
      implicit none
!
!
      interface
        function c_getpid() bind(C, name="getpid")
          use iso_c_binding, only: c_int
          integer(c_int) :: c_getpid
        end function c_getpid
      end interface
!
!
      CONTAINS
!
!PROCEDURE print_memory_usage
      subroutine print_memory_usage(iOut,label)
      use iso_c_binding, only: c_int
      implicit none
      integer::iOut,ios,unit,pid
      real::memVal
      character(len=*),intent(in)::label
      character(len=512)::line,filename,memType,memUnits
!
      pid = c_getpid()
      write(filename,'(A,I0,A)') '/proc/',pid,'/status'
!
      open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
      if(ios.ne.0) then
        print *, "Error opening ", trim(filename)
        return
      endIf
      do
        read(unit,'(A)',iostat=ios) line
        if(ios.ne.0) exit
        if(index(line,"VmRSS:") == 1) then
          read(line,*) memType,memVal,memUnits
          if(memVal.gt.1023.0) then
            memVal = memVal/1024
            memUnits = 'mb'
          endIf
          if(memVal.gt.1023.0) then
            memVal = memVal/1024
            memUnits = 'gb'
          endIf
          write(iOut,'(A,1X,A,": ",f6.1,1x,A)') '*** Memory usage at ',  &
            Trim(label),memVal,Trim(memUnits)
          exit
        endIf
      endDo
      close(unit)
!
      return
      end subroutine print_memory_usage


      end module memory_utils
