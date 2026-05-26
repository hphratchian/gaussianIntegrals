      program unitTest2
!
!     This program tests the lab-frame grid builders used for simple rotational
!     averaging models in PAD.
!
!     Command line arguments:
!           1. lab-frame model flag (optional; default=0)
!              0 = 3 Cartesian lab-frame orientations
!              1 = sphere-grid lab-frame orientations
!           2. number of lab-frame theta points (optional; default=5)
!           3. number of lab-frame phi points (optional; default=8)
!
!
!     H. P. Hratchian, 2026.
!
      use pad_mod
      implicit none
      integer(kind=int64),parameter::nOMP=1
      integer(kind=int64)::i,nLabFrames,nExpected
      real(kind=real64),parameter::tolGeom=1.0d-10,tolWeight=1.0d-10
      real(kind=real64)::weightSum,expectedWeightSum
      real(kind=real64),dimension(3)::weightedMoment
      real(kind=real64),dimension(:),allocatable::labFrameWeights
      real(kind=real64),dimension(:,:),allocatable::epsilonVectors,  &
        kPlaneVectors
      character(len=8),dimension(:),allocatable::planeLabels
      type(pad_options)::options
!
!     Format statements.
!
 1000 format(1x,'Starting unitTest2')
 1100 format(1x,'Lab-frame model flag   = ',i3)
 1110 format(1x,'nLabFrameTheta        = ',i6)
 1120 format(1x,'nLabFramePhi          = ',i6)
 2000 format(1x,'No. of lab frames     = ',i8)
 2010 format(1x,'Sum of weights        = ',f20.10,3x,'expected = ',f20.10)
 2020 format(1x,'Weighted first moment = ',3f20.10)
 2999 format(/,1x,'unitTest2 complete.')
!
!     Start the unit test program.
!
      call parseUnitTest2CommandLine(options)
      write(iOut,1000)
      call omp_set_num_threads(nOMP)
      call mqc_version_print(iOut)
      write(iOut,1100) options%labFrameType
      write(iOut,1110) options%nLabFrameTheta
      write(iOut,1120) options%nLabFramePhi
!
!     Build the requested lab-frame grid and verify its geometry.
!
      call buildPADLabFrames(options,epsilonVectors,kPlaneVectors,  &
        planeLabels,labFrameWeights)
      nLabFrames = Size(labFrameWeights)
      weightSum = SUM(labFrameWeights)
      write(iOut,2000) nLabFrames
!
      do i = 1,nLabFrames
        call assertNear(sqrt(dot_product(epsilonVectors(:,i),epsilonVectors(:,i))),  &
          mqc_float(1),tolGeom,'epsilon vector norm')
        call assertNear(sqrt(dot_product(kPlaneVectors(:,i),kPlaneVectors(:,i))),  &
          mqc_float(1),tolGeom,'k-plane vector norm')
        call assertNear(dot_product(epsilonVectors(:,i),kPlaneVectors(:,i)),  &
          mqc_float(0),tolGeom,'epsilon/k-plane orthogonality')
      endDo
!
!     Run model-specific tests.
!
      if(options%labFrameType.eq.PAD_LAB_FRAMES_CARTESIAN) then
        nExpected = 3
        expectedWeightSum = mqc_float(3)
        call assertNear(weightSum,expectedWeightSum,tolWeight,  &
          'Cartesian weight sum')
        if(nLabFrames.ne.nExpected) then
          write(iOut,*)' nLabFrames = ',nLabFrames
          call mqc_error('unitTest2: unexpected number of Cartesian lab frames.')
        endIf
        call assertNear(labFrameWeights(1),mqc_float(1),tolWeight,  &
          'Cartesian weight 1')
        call assertNear(labFrameWeights(2),mqc_float(1),tolWeight,  &
          'Cartesian weight 2')
        call assertNear(labFrameWeights(3),mqc_float(1),tolWeight,  &
          'Cartesian weight 3')
        if(TRIM(planeLabels(1)).ne.'xy'.or.TRIM(planeLabels(2)).ne.'yz'.or.  &
          TRIM(planeLabels(3)).ne.'zx') then
          write(iOut,*)' planeLabels = ',planeLabels
          call mqc_error('unitTest2: unexpected Cartesian plane labels.')
        endIf
!
      elseIf(options%labFrameType.eq.PAD_LAB_FRAMES_SPHERE) then
        nExpected = (options%nLabFrameTheta-2)*options%nLabFramePhi+2
        expectedWeightSum = mqc_float(4)*Pi
        call assertNear(weightSum,expectedWeightSum,tolWeight,  &
          'sphere weight sum')
        if(nLabFrames.ne.nExpected) then
          write(iOut,*)' nLabFrames = ',nLabFrames
          write(iOut,*)' nExpected  = ',nExpected
          call mqc_error('unitTest2: unexpected number of sphere lab frames.')
        endIf
        if(MINVAL(labFrameWeights).le.mqc_float(0)) then
          write(iOut,*)' min weight = ',MINVAL(labFrameWeights)
          call mqc_error('unitTest2: sphere weights must all be positive.')
        endIf
        weightedMoment = mqc_float(0)
        do i = 1,nLabFrames
          weightedMoment = weightedMoment+  &
            labFrameWeights(i)*epsilonVectors(:,i)
        endDo
        write(iOut,2020) weightedMoment
        call assertNear(weightedMoment(1),mqc_float(0),tolWeight,  &
          'sphere weighted x moment')
        call assertNear(weightedMoment(2),mqc_float(0),tolWeight,  &
          'sphere weighted y moment')
        call assertNear(weightedMoment(3),mqc_float(0),tolWeight,  &
          'sphere weighted z moment')
!
      else
        call mqc_error('unitTest2 only supports Cartesian and sphere lab-frame models.')
      endIf
!
      write(iOut,2010) weightSum,expectedWeightSum
      write(iOut,2999)
!
      contains

!PROCEDURE assertNear
      subroutine assertNear(value,target,tolerance,label)
!
!     This routine checks whether two real values are equal within a specified
!     tolerance.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),intent(in)::value,target,tolerance
      character(len=*),intent(in)::label
!
      if(abs(value-target).gt.tolerance) then
        write(iOut,*)' label     = ',TRIM(label)
        write(iOut,*)' value     = ',value
        write(iOut,*)' target    = ',target
        write(iOut,*)' tolerance = ',tolerance
        call mqc_error('unitTest2: scalar comparison failed.')
      endIf
!
      return
      end subroutine assertNear

!
!PROCEDURE parseUnitTest2CommandLine
      subroutine parseUnitTest2CommandLine(options)
!
!     This routine processes the unitTest2 command line and fills the PAD
!     options object used to build lab-frame grids.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      type(pad_options),intent(out)::options
      integer(kind=int64)::nCommandLineArgs
!
      options = pad_options()
      options%printResults = .false.
      options%printThetaTable = .false.
      options%nOMP = nOMP
      nCommandLineArgs = command_argument_count()
      if(nCommandLineArgs.gt.3)  &
        call mqc_error('unitTest2 expects 0-3 command line arguments.')
      if(nCommandLineArgs.ge.1) then
        call mqc_get_command_argument_integer(1,options%labFrameType)
      endIf
      if(nCommandLineArgs.ge.2) then
        call mqc_get_command_argument_integer(2,options%nLabFrameTheta)
      endIf
      if(nCommandLineArgs.ge.3) then
        call mqc_get_command_argument_integer(3,options%nLabFramePhi)
      endIf
      if(options%labFrameType.eq.PAD_LAB_FRAMES_CUSTOM) then
        call mqc_error('unitTest2 does not support custom lab-frame input.')
      endIf
!
      return
      end subroutine parseUnitTest2CommandLine

      end program unitTest2
