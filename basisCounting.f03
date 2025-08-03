      program basisCounting
!
!     This program takes an orbital angular momentum quantum number, l, at the
!     command line and evaluates the number of functions in a shell with orbital
!     angular momentum quantum number l. The program also tests MQCPack code
!     that outputs the list of (lx,ly,lz) values for a full set of AO basis
!     functions in the type of shell requested by the user.
!
!     Hrant P. Hratchian, 2023.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      USE iso_fortran_env
      USE MQC_Integrals
      implicit none
      integer(kind=int64),parameter::iOut=6
      integer(kind=int64)::nArguments,lInput,nCartesians
      character(len=256)::commandLineArg
      logical::fail=.false.
!
!     Format statements.
!
 2000 format(1x,'l: ',i2,3x,'N Cartesians: ',i5)
 9000 format(1x,'Expected 1 (and only 1) command line argument.')
 9010 format(1x,'Invalid l. l must be >= 0.')
!
!
!     Begin by reading the command line to get lInput.
!
      nArguments = COMMAND_ARGUMENT_COUNT()
      if(nArguments.ne.1) then
        write(iOut,9000)
        fail = .true.
      endIf
      if(fail) goto 999
      call GET_COMMAND_ARGUMENT(1,commandLineArg)
      read(commandLineArg,'(I1)') lInput
      if(lInput.lt.0) then
        write(iOut,9010)
        fail = .true.
      endIf
      if(fail) goto 999
!
!     Using lInput, determine nCartesians and report the result.
!
      nCartesians = (lInput+1)*(lInput+2)
      nCartesians = nCartesians/2
      write(iOut,2000) lInput,nCartesians
!
!     The end of the program.
!
  999 Continue
      end program basisCounting

