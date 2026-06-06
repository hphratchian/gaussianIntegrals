      program unitTest2
!
!     This program tests the lab-frame grid builders used for simple rotational
!     averaging models in PAD.
!
!     Command line arguments:
!           1. lab-frame model flag (optional; default=0)
!             -1 = custom lab-frame orientations filled by this test
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
      real(kind=real64),dimension(3,2)::expectedEpsilon,expectedKPlane
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
      if(options%labFrameType.eq.PAD_LAB_FRAMES_CUSTOM) then
        nExpected = 2
        expectedWeightSum = mqc_float(5)
        call assertNear(weightSum,expectedWeightSum,tolWeight,  &
          'custom weight sum')
        if(nLabFrames.ne.nExpected) then
          write(iOut,*)' nLabFrames = ',nLabFrames
          call mqc_error('unitTest2: unexpected number of custom lab frames.')
        endIf
        expectedEpsilon(:,1) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
        expectedEpsilon(:,2) = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
        expectedKPlane(:,1) = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
        expectedKPlane(:,2) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
        call assertVectorNear(epsilonVectors(:,1),expectedEpsilon(:,1),  &
          tolGeom,'custom epsilon 1')
        call assertVectorNear(epsilonVectors(:,2),expectedEpsilon(:,2),  &
          tolGeom,'custom epsilon 2')
        call assertVectorNear(kPlaneVectors(:,1),expectedKPlane(:,1),  &
          tolGeom,'custom k-plane 1')
        call assertVectorNear(kPlaneVectors(:,2),expectedKPlane(:,2),  &
          tolGeom,'custom k-plane 2')
        call assertNear(labFrameWeights(1),mqc_float(2),tolWeight,  &
          'custom weight 1')
        call assertNear(labFrameWeights(2),mqc_float(3),tolWeight,  &
          'custom weight 2')
        if(TRIM(planeLabels(1)).ne.'custx'.or.  &
          TRIM(planeLabels(2)).ne.'custz') then
          write(iOut,*)' planeLabels = ',planeLabels
          call mqc_error('unitTest2: unexpected custom plane labels.')
        endIf
!
      elseIf(options%labFrameType.eq.PAD_LAB_FRAMES_CARTESIAN) then
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
        call mqc_error('unitTest2: unsupported lab-frame model.')
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
!PROCEDURE assertVectorNear
      subroutine assertVectorNear(value,target,tolerance,label)
!
!     This routine checks whether two real vectors are equal within a specified
!     component-wise tolerance.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::value,target
      real(kind=real64),intent(in)::tolerance
      character(len=*),intent(in)::label
!
      if(Size(value).ne.Size(target)) then
        write(iOut,*)' label       = ',TRIM(label)
        write(iOut,*)' Size(value) = ',Size(value)
        write(iOut,*)' Size(target)= ',Size(target)
        call mqc_error('unitTest2: vector shape comparison failed.')
      endIf
      if(maxval(abs(value-target)).gt.tolerance) then
        write(iOut,*)' label     = ',TRIM(label)
        call mqc_print(value,iOut,header='value')
        call mqc_print(target,iOut,header='target')
        write(iOut,*)' tolerance = ',tolerance
        call mqc_error('unitTest2: vector comparison failed.')
      endIf
!
      return
      end subroutine assertVectorNear

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
      if(options%labFrameType.eq.PAD_LAB_FRAMES_CUSTOM)  &
        call setupUnitTest2CustomLabFrames(options)
!
      return
      end subroutine parseUnitTest2CommandLine

!
!PROCEDURE setupUnitTest2CustomLabFrames
      subroutine setupUnitTest2CustomLabFrames(options)
!
!     This routine fills a small custom lab-frame set directly in the options
!     object to exercise the programmatic PAD_LAB_FRAMES_CUSTOM path.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      type(pad_options),intent(inout)::options
!
      Allocate(options%labEpsilonVector(3,2),  &
        options%labKPlaneVector(3,2),options%labFrameWeights(2),  &
        options%labFrameLabels(2))
      options%labEpsilonVector(:,1) =  &
        [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      options%labEpsilonVector(:,2) =  &
        [ mqc_float(0),mqc_float(0),mqc_float(1) ]
      options%labKPlaneVector(:,1) =  &
        [ mqc_float(0),mqc_float(1),mqc_float(0) ]
      options%labKPlaneVector(:,2) =  &
        [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      options%labFrameWeights = [ mqc_float(2),mqc_float(3) ]
      options%labFrameLabels(1) = 'custx'
      options%labFrameLabels(2) = 'custz'
!
      return
      end subroutine setupUnitTest2CustomLabFrames

      end program unitTest2
