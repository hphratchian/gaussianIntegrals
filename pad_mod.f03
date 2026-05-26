      module pad_mod
!
!     This module supports the PAD program.
!
      use iso_fortran_env
      use omp_lib
      use mqc_general
      use mqc_integrals1
      use mqc_gaussian
      use memory_utils
      use gbs_mod
!
      implicit none
      integer(kind=int64),parameter::PAD_LAB_FRAMES_CARTESIAN=0_int64
      integer(kind=int64),parameter::PAD_LAB_FRAMES_SPHERE=1_int64
      integer(kind=int64),parameter::PAD_LAB_FRAMES_CUSTOM=2_int64
!
!
!     The pad_options object collects user-facing and model-control options for
!     a PAD calculation. The photon and binding energies are supplied in eV.
      type pad_options
        integer(kind=int64)::dysonMOIndex=1_int64
        integer(kind=int64)::nGridPointsTheta=5_int64
        integer(kind=int64)::nGridPointsM=101_int64
        integer(kind=int64)::nChi=36_int64
        integer(kind=int64)::iPEType=0_int64
        integer(kind=int64)::lMax=6_int64
        integer(kind=int64)::nOMP=1_int64
        integer(kind=int64)::labFrameType=PAD_LAB_FRAMES_CARTESIAN
        integer(kind=int64)::nLabFrameTheta=5_int64
        integer(kind=int64)::nLabFramePhi=8_int64
        logical::printResults=.true.
        logical::printThetaTable=.true.
        real(kind=real64)::photonEnergyEV=0.0_real64
        real(kind=real64)::bindingEnergyEV=0.0_real64
        real(kind=real64),dimension(:,:),allocatable::labEpsilonVector
        real(kind=real64),dimension(:,:),allocatable::labKPlaneVector
        real(kind=real64),dimension(:),allocatable::labFrameWeights
        character(len=8),dimension(:),allocatable::labFrameLabels
      end type pad_options
!
!     The pad_results object stores the central numerical results from a PAD
!     calculation so that other codes can call the driver as a routine.
      type pad_results
        integer(kind=int64)::nLabFrames=0_int64
        integer(kind=int64)::nTheta=0_int64
        integer(kind=int64)::nChi=0_int64
        real(kind=real64)::photoelectronEnergyEV=0.0_real64
        real(kind=real64)::photoelectronEnergyHartree=0.0_real64
        real(kind=real64)::kMag=0.0_real64
        real(kind=real64)::dysonSelfOverlap=0.0_real64
        real(kind=real64)::labFrameWeightSum=0.0_real64
        real(kind=real64)::chiWeightSum=0.0_real64
        real(kind=real64)::averageIntensity0=0.0_real64
        real(kind=real64)::averageIntensity90=0.0_real64
        real(kind=real64)::averageIntegratedIntensity=0.0_real64
        real(kind=real64)::averageSolidAngleIntegratedIntensity=0.0_real64
        real(kind=real64)::averageBetaParaPerp=0.0_real64
        real(kind=real64)::averageBetaFit=0.0_real64
        real(kind=real64)::averageRSquared=0.0_real64
        real(kind=real64)::meanBetaParaPerp=0.0_real64
        real(kind=real64)::meanBetaFit=0.0_real64
        real(kind=real64),dimension(:),allocatable::theta
        real(kind=real64),dimension(:),allocatable::thetaWeights
        real(kind=real64),dimension(:),allocatable::thetaSolidAngleWeights
        real(kind=real64),dimension(:),allocatable::chi
        real(kind=real64),dimension(:),allocatable::chiWeights
        real(kind=real64),dimension(:),allocatable::labFrameWeights
        real(kind=real64),dimension(:),allocatable::integratedIntensity
        real(kind=real64),dimension(:),allocatable::solidAngleIntegratedIntensity
        real(kind=real64),dimension(:),allocatable::intensity0
        real(kind=real64),dimension(:),allocatable::intensity90
        real(kind=real64),dimension(:),allocatable::betaValsParaPerp
        real(kind=real64),dimension(:),allocatable::betaValsFit
        real(kind=real64),dimension(:),allocatable::rSquared
        real(kind=real64),dimension(:),allocatable::averageIntensityTheta
        real(kind=real64),dimension(:,:),allocatable::intensityTheta
        real(kind=real64),dimension(:,:),allocatable::epsilonVector
        real(kind=real64),dimension(:,:),allocatable::kPlaneVector
        character(len=8),dimension(:),allocatable::intPlaneLabels
      end type pad_results
!
!
      CONTAINS
!
!
!PROCEDURE padCommandLine
      subroutine padCommandLine(options,fafName)
!
!     This routine processes the command line arguments and fills the PAD
!     options object and FAF filename.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      type(pad_options),intent(out)::options
      character(len=256),intent(out)::fafName
!
      integer(kind=int64)::nCommandLineArgs
!
!     Walk through the command line arguments to fill required input parameters
!     and other option flags.
!
      options = pad_options()
      nCommandLineArgs = command_argument_count()
      if(nCommandLineArgs.lt.4.or.nCommandLineArgs.gt.11)  &
        call mqc_error('PAD expects 4-11 command line arguments.')
      call get_command_argument(1,fafName)
      call mqc_get_command_argument_integer(2,options%dysonMOIndex)
      call mqc_get_command_argument_real(3,options%photonEnergyEV)
      call mqc_get_command_argument_real(4,options%bindingEnergyEV)
      if(nCommandLineArgs.ge.5) then
        call mqc_get_command_argument_integer(5,options%nGridPointsTheta)
      endIf
      if(nCommandLineArgs.ge.6) then
        call mqc_get_command_argument_integer(6,options%nGridPointsM)
      endIf
      if(nCommandLineArgs.ge.7) then
        call mqc_get_command_argument_integer(7,options%iPEType)
      endIf
      if(nCommandLineArgs.ge.8) then
        call mqc_get_command_argument_integer(8,options%labFrameType)
      endIf
      if(nCommandLineArgs.ge.9) then
        call mqc_get_command_argument_integer(9,options%nLabFrameTheta)
      endIf
      if(nCommandLineArgs.ge.10) then
        call mqc_get_command_argument_integer(10,options%nLabFramePhi)
      endIf
      if(nCommandLineArgs.ge.11) then
        call mqc_get_command_argument_integer(11,options%nChi)
      endIf
!
      return
      end subroutine padCommandLine

!
!PROCEDURE padComputePhotoelectronEnergyAndK
      subroutine padComputePhotoelectronEnergyAndK(options,  &
        photoelectronEnergyEV,  &
        photoelectronEnergyHartree,kMag)
!
!     This routine computes the photoelectron kinetic energy and plane-wave
!     magnitude from the user-facing photon and binding energies.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      type(pad_options),intent(in)::options
      real(kind=real64),intent(out)::photoelectronEnergyEV,  &
        photoelectronEnergyHartree,kMag
!
!     Convert the user-facing energy inputs to the plane-wave magnitude used by
!     the transition-moment routines.
!
      if(options%photonEnergyEV.le.mqc_float(0))  &
        call mqc_error('PAD: photon energy must be positive.')
      if(options%bindingEnergyEV.le.mqc_float(0))  &
        call mqc_error('PAD: binding energy must be positive.')
      photoelectronEnergyEV = options%photonEnergyEV-options%bindingEnergyEV
      if(photoelectronEnergyEV.le.mqc_float(0))  &
        call mqc_error('PAD: photon energy must exceed binding energy.')
      photoelectronEnergyHartree = photoelectronEnergyEV/evPHartree
      kMag = sqrt(mqc_float(2)*photoelectronEnergyHartree)
!
      return
      end subroutine padComputePhotoelectronEnergyAndK


!PROCEDURE runPADCalculation
      subroutine runPADCalculation(faf,options,results)
!
!     This routine carries out the PAD calculation for a selected Dyson orbital
!     from an already-loaded Gaussian FAF object. The molecule is kept fixed and
!     the lab-frame polarization vectors are varied. For each polarization
!     direction, the scan plane can also be rotated uniformly about epsilon
!     through a periodic chi quadrature.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      type(mqc_gaussian_unformatted_matrix_file),intent(inout)::faf
      type(pad_options),intent(in)::options
      type(pad_results),intent(out)::results
!
      integer(kind=int64)::i,j,iChi,nIntPlanes
      real(kind=real64)::tStart1,tEnd1,stepSizeIntM,stepSizeTheta,  &
        thetaStart,MSquared0,MSquared90
      real(kind=real64),dimension(3)::cartStart,cartEnd,epsilonVec,  &
        epsilonUnit,uBasis,vBasis,uChi
      real(kind=real64),dimension(:),allocatable::lWeights0,lWeights90,  &
        lWeights0Tmp,lWeights90Tmp
      real(kind=real64),dimension(:),allocatable::quadWeightsM,  &
        MSquaredList
      real(kind=real64),dimension(:,:),allocatable::quadGridM,moCoeffs
      logical::found
      type(mqc_basisset)::basisSet
      type(MQC_Variable)::tmp,quadTmp
!
 1100 format(/,1x,'Lab frame ',i4,': epsilon = (',f7.3,',',f7.3,',',f7.3,')',  &
        '  reference k-plane vector = (',f7.3,',',f7.3,',',f7.3,')')
 1200 format(1x,'Photoelectron model flag = ',i3)
 1210 format(1x,'Photon energy          = ',f12.6,' eV')
 1220 format(1x,'Electron binding energy= ',f12.6,' eV')
 1230 format(1x,'Electron kinetic energy= ',f12.6,' eV = ',es14.6,' Eh')
 1240 format(1x,'Photoelectron k        = ',f12.6,' a.u.')
 1250 format(1x,'Lab-frame model flag   = ',i3)
 1260 format(1x,'Lab-frame orientations = ',i8,3x,'weight sum = ',es14.6)
 1270 format(1x,'Chi quadrature points  = ',i8,3x,'weight sum = ',es14.6)
  2000 format(1x,'Data for the ',A,' grid:',/,  &
        3x,'nPoints = ',i15,' step size = ',f20.5)
 2010 format(1x,'Data for the ',A,' grid:',/,3x,'nPoints = ',i15)
 2500 format(/,1x,'Dyson orbital self-overlap = ',f12.8)
 3000 format(/,1x,55('-'),/,  &
        12x,'Intensity as a Function of Theta',/,  &
        17x,'(k = ',f10.6,' a.u.)',/,  &
        9x,'theta (rad)',15x,'I(theta)',/,  &
        1x,55('='))
 3010 format(9x,f7.3,3x,es25.12)
 3020 format(1x,55('='))
 3100 format(1x,'I(0) = ',es15.8,3x,'I(90) = ',es15.8,3x,'beta = ',f12.6)
 3200 format(1x,'l weights at theta = ',A,':',*(1x,f9.5))
 3500 format(/,1x,87('-'),/,  &
        31x,'Summary of PAD Calculation',/,  &
        34x,'(k = ',f10.6,' a.u.)',/,  &
        1x,87('='))
  3510 format(1x,'Plane: ',A,'  |  epsilon: (',f5.2,',',f5.2,',',f5.2,  &
        ')',/,18x,'theta-integrated intensity       = ',es14.6,/,  &
        18x,'solid-angle integrated intensity = ',es14.6,/,  &
        18x,'beta(ratio) = ',f10.6,'  |  beta(fit) = ',f10.6,  &
        ' (R**2 = ',f8.5,')')
  3520 format(1x,87('='))
  3530 format(1x,'Weighted mean beta(ratio) = ',f10.6,/,  &
        1x,'Weighted mean beta(fit)   = ',f10.6)
 3540 format(1x,'Averaged theta-integrated intensity = ',es14.6,/,  &
        1x,'Averaged solid-angle integrated intensity = ',es14.6,/,  &
        1x,'Beta from averaged PAD(ratio) = ',f10.6,/,  &
        1x,'Beta from averaged PAD(fit)   = ',f10.6,  &
        ' (R**2 = ',f8.5,')',/)
 8998 format(1x,'Time for ',A,': ',f15.1,' s')
!
!     Validate user/model options and compute the photoelectron kinematics.
!
      if(options%dysonMOIndex.lt.1)  &
        call mqc_error('PAD: requested Dyson orbital index is out of range.')
      if(options%nGridPointsTheta.lt.2)  &
        call mqc_error('PAD: nGridPointsTheta must be at least 2.')
      if(options%nGridPointsM.lt.2)  &
        call mqc_error('PAD: nGridPointsM must be at least 2.')
      if(options%nChi.lt.1)  &
        call mqc_error('PAD: nChi must be positive.')
      if(options%iPEType.lt.0.or.options%iPEType.gt.2)  &
        call mqc_error('PAD: invalid photoelectron model flag.')
      if(options%lMax.lt.0)  &
        call mqc_error('PAD: lMax must be nonnegative.')
      if(options%nLabFrameTheta.lt.2)  &
        call mqc_error('PAD: nLabFrameTheta must be at least 2.')
      if(options%nLabFramePhi.lt.1)  &
        call mqc_error('PAD: nLabFramePhi must be positive.')
      call padComputePhotoelectronEnergyAndK(options,  &
        results%photoelectronEnergyEV,results%photoelectronEnergyHartree,  &
        results%kMag)
      if(options%printResults) then
        write(iOut,1200) options%iPEType
        write(iOut,1210) options%photonEnergyEV
        write(iOut,1220) options%bindingEnergyEV
        write(iOut,1230) results%photoelectronEnergyEV,  &
          results%photoelectronEnergyHartree
        write(iOut,1240) results%kMag
      endIf
!
!     Set up the lab-frame polarizations. These are rotations of the lab frame
!     relative to a fixed molecular frame.
!
      call buildPADLabFrames(options,results%epsilonVector,  &
        results%kPlaneVector,results%intPlaneLabels,  &
        results%labFrameWeights)
      nIntPlanes = Size(results%epsilonVector,2)
      results%labFrameWeightSum = SUM(results%labFrameWeights)
      if(results%labFrameWeightSum.le.mqc_small)  &
        call mqc_error('PAD: lab-frame weights sum to zero.')
      if(options%printResults) then
        write(iOut,1250) options%labFrameType
        write(iOut,1260) nIntPlanes,results%labFrameWeightSum
      endIf
      results%nLabFrames = nIntPlanes
      results%nTheta = options%nGridPointsTheta
      results%nChi = options%nChi
      call buildChiQuadrature(results%chi,results%chiWeights,results%nChi)
      results%chiWeightSum = SUM(results%chiWeights)
      if(results%chiWeightSum.le.mqc_small)  &
        call mqc_error('PAD: chi quadrature weights sum to zero.')
      if(options%printResults) write(iOut,1270) results%nChi,  &
        results%chiWeightSum
      Allocate(results%integratedIntensity(nIntPlanes),  &
        results%solidAngleIntegratedIntensity(nIntPlanes),  &
        results%intensity0(nIntPlanes),results%intensity90(nIntPlanes),  &
        results%betaValsParaPerp(nIntPlanes),  &
        results%betaValsFit(nIntPlanes),results%rSquared(nIntPlanes),  &
        results%intensityTheta(options%nGridPointsTheta,nIntPlanes),  &
        results%averageIntensityTheta(options%nGridPointsTheta),  &
        results%theta(options%nGridPointsTheta),  &
        results%thetaWeights(options%nGridPointsTheta),  &
        results%thetaSolidAngleWeights(options%nGridPointsTheta))
      results%integratedIntensity = mqc_float(0)
      results%solidAngleIntegratedIntensity = mqc_float(0)
      results%intensity0 = mqc_float(0)
      results%intensity90 = mqc_float(0)
      results%betaValsParaPerp = mqc_float(0)
      results%betaValsFit = mqc_float(0)
      results%rSquared = mqc_float(0)
      results%intensityTheta = mqc_float(0)
      results%averageIntensityTheta = mqc_float(0)
!
!     Read the basis set and MO coefficients from the FAF.
!
      call loadGaussianBasisSet(faf,basisSet)
      call faf%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmp)
      moCoeffs = tmp
      if(options%dysonMOIndex.gt.Size(moCoeffs,2))  &
        call mqc_error('PAD: requested Dyson orbital index is out of range.')
!
!     Use the Gaussian XC quadrature grid from the FAF when available. Otherwise
!     fall back to a simple Cartesian trapezoid grid.
!
      call CPU_TIME(tStart1)
      call faf%getArray('3D Quadrature Grid',mqcVarOut=quadTmp,foundOut=found)
      if(found) then
        Allocate(quadWeightsM(SIZE(quadTmp,2)),quadGridM(3,SIZE(quadTmp,2)))
        do i = 1,SIZE(quadTmp,2)
          quadWeightsM(i) = quadTmp%getVal([ 1,i ])
          quadGridM(1,i) = quadTmp%getVal([ 2,i ])
          quadGridM(2,i) = quadTmp%getVal([ 3,i ])
          quadGridM(3,i) = quadTmp%getVal([ 4,i ])
        endDo
        if(options%printResults) write(iOut,2010)  &
          'FAF M quadrature',SIZE(quadWeightsM)
      else
        cartStart = [ -mqc_float(6),-mqc_float(6),-mqc_float(6) ]
        cartEnd = [ mqc_float(6),mqc_float(6),mqc_float(6) ]
        stepSizeIntM = (cartEnd(1)-cartStart(1))/  &
          mqc_float(options%nGridPointsM-1)
        Allocate(quadGridM(3,options%nGridPointsM**3),  &
          quadWeightsM(options%nGridPointsM**3))
        call setup_quadrature_trapezoid3d(options%nGridPointsM,  &
          stepSizeIntM,cartStart,quadGridM,quadWeightsM)
        if(options%printResults) write(iOut,2000) 'fallback M quadrature',  &
          options%nGridPointsM**3,stepSizeIntM
      endIf
      call CPU_TIME(tEnd1)
      if(options%printResults) write(iOut,8998)  &
        'M quadrature setup',tEnd1-tStart1
      if(options%printResults) flush(iOut)
!
!     Fill the theta grid for the PAD scan.
!
      thetaStart = mqc_float(0)
      stepSizeTheta = Pi/mqc_float(options%nGridPointsTheta-1)
      call setup_quadrature_trapezoid1d(options%nGridPointsTheta,  &
        stepSizeTheta,thetaStart,results%theta,results%thetaWeights)
      results%thetaSolidAngleWeights = sin(results%theta)*results%thetaWeights
      if(options%printResults) write(iOut,2000)  &
        'theta',options%nGridPointsTheta,stepSizeTheta
      if(options%printResults) flush(iOut)
!
!     Test that the chosen MO is normalized on the quadrature grid.
!
      call CPU_TIME(tStart1)
      results%dysonSelfOverlap = moInnerProductNumericalIntegration(  &
        moCoeffs(:,options%dysonMOIndex),quadGridM,quadWeightsM,basisSet)
      if(options%printResults) write(iOut,2500) results%dysonSelfOverlap
      call CPU_TIME(tEnd1)
      if(options%printResults) write(iOut,8998) 'norm test',tEnd1-tStart1
      if(options%printResults) flush(iOut)
!
!     Evaluate the PAD for each lab-frame orientation.
!
      Allocate(lWeights0(0:options%lMax),lWeights90(0:options%lMax),  &
        lWeights0Tmp(0:options%lMax),lWeights90Tmp(0:options%lMax),  &
        MSquaredList(options%nGridPointsTheta))
      do i = 1,nIntPlanes
        call buildTransverseBasis(results%epsilonVector(:,i),  &
          results%kPlaneVector(:,i),epsilonUnit,uBasis,vBasis)
        results%epsilonVector(:,i) = epsilonUnit
        results%kPlaneVector(:,i) = uBasis
        if(options%printResults) write(iOut,1100) i,  &
          results%epsilonVector(:,i),results%kPlaneVector(:,i)
        results%intensityTheta(:,i) = mqc_float(0)
        results%intensity0(i) = mqc_float(0)
        results%intensity90(i) = mqc_float(0)
        lWeights0 = mqc_float(0)
        lWeights90 = mqc_float(0)
        call CPU_TIME(tStart1)
!
        do iChi = 1,results%nChi
          epsilonVec = results%epsilonVector(:,i)
          uChi = cos(results%chi(iChi))*uBasis+sin(results%chi(iChi))*vBasis
          call mqc_normalizeVector(uChi)
          if(options%iPEType.eq.0) then
            MSquared0 = dysonPlaneWaveMatrixElementSquared(mqc_float(0),  &
              results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM)
            MSquared90 = dysonPlaneWaveMatrixElementSquared(Pi/mqc_float(2),  &
              results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM)
            MSquaredList = dysonPlaneWaveMatrixElementSquaredThetaList(  &
              results%theta,results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM)
          else
            call dysonMatrixElement1Angle(options%iPEType,options%lMax,  &
              mqc_float(0),results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM,MSquared0,lWeights0Tmp)
            call dysonMatrixElement1Angle(options%iPEType,options%lMax,  &
              Pi/mqc_float(2),results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM,MSquared90,lWeights90Tmp)
            call dysonMatrixElementThetaList(options%iPEType,options%lMax,  &
              results%theta,results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM,MSquaredList)
            if(options%iPEType.eq.2) then
              lWeights0 = lWeights0+results%chiWeights(iChi)*lWeights0Tmp
              lWeights90 = lWeights90+results%chiWeights(iChi)*lWeights90Tmp
            endIf
          endIf
          results%intensityTheta(:,i) = results%intensityTheta(:,i)+  &
            results%chiWeights(iChi)*MSquaredList
          results%intensity0(i) = results%intensity0(i)+  &
            results%chiWeights(iChi)*MSquared0
          results%intensity90(i) = results%intensity90(i)+  &
            results%chiWeights(iChi)*MSquared90
        endDo
        results%intensityTheta(:,i) = results%intensityTheta(:,i)/  &
          results%chiWeightSum
        results%intensity0(i) = results%intensity0(i)/results%chiWeightSum
        results%intensity90(i) = results%intensity90(i)/results%chiWeightSum
        if(options%iPEType.eq.2) then
          lWeights0 = lWeights0/results%chiWeightSum
          lWeights90 = lWeights90/results%chiWeightSum
        endIf
!
        call CPU_TIME(tEnd1)
        if(options%printResults) then
          if(results%nChi.gt.1) then
            write(iOut,8998) 'chi-averaged MSquared(theta) list',tEnd1-tStart1
          else
            write(iOut,8998) 'MSquared(theta) list',tEnd1-tStart1
          endIf
        endIf
        if(options%iPEType.eq.2.and.options%printResults) then
          write(iOut,3200) '0 ',lWeights0
          write(iOut,3200) '90',lWeights90
        endIf
!
!       Compute beta using both the parallel/perpendicular ratio and the
!       least-squares fit to the PAD shape.
!
        call betaLeastSquares(results%intensityTheta(:,i),results%theta,  &
          results%betaValsFit(i),results%rSquared(i))
        results%integratedIntensity(i) =  &
          dot_product(results%intensityTheta(:,i),results%thetaWeights)
        results%solidAngleIntegratedIntensity(i) = results%chiWeightSum*  &
          dot_product(results%intensityTheta(:,i),  &
          results%thetaSolidAngleWeights)
        results%betaValsParaPerp(i) = betaParaPerp(results%intensity0(i),  &
          results%intensity90(i))
!
        if(options%printResults.and.options%printThetaTable) then
          write(iOut,3000) results%kMag
          do j = 1,options%nGridPointsTheta
            write(iOut,3010) results%theta(j),results%intensityTheta(j,i)
          endDo
          write(iOut,3020)
          write(iOut,3100) results%intensity0(i),results%intensity90(i),  &
            results%betaValsParaPerp(i)
          flush(iOut)
        endIf
      endDo
!
!     Build the weighted orientation-averaged PAD and associated beta values.
!
      do i = 1,nIntPlanes
        results%averageIntensityTheta = results%averageIntensityTheta +  &
          results%labFrameWeights(i)*results%intensityTheta(:,i)
      endDo
      results%averageIntensityTheta = results%averageIntensityTheta/  &
        results%labFrameWeightSum
      results%averageIntensity0 = dot_product(results%labFrameWeights,  &
        results%intensity0)/results%labFrameWeightSum
      results%averageIntensity90 = dot_product(results%labFrameWeights,  &
        results%intensity90)/results%labFrameWeightSum
      results%averageIntegratedIntensity =  &
        dot_product(results%averageIntensityTheta,results%thetaWeights)
      results%averageSolidAngleIntegratedIntensity = results%chiWeightSum*  &
        dot_product(results%averageIntensityTheta,  &
        results%thetaSolidAngleWeights)
      results%averageBetaParaPerp = betaParaPerp(results%averageIntensity0,  &
        results%averageIntensity90)
      call betaLeastSquares(results%averageIntensityTheta,results%theta,  &
        results%averageBetaFit,results%averageRSquared)
      results%meanBetaParaPerp = dot_product(results%labFrameWeights,  &
        results%betaValsParaPerp)/results%labFrameWeightSum
      results%meanBetaFit = dot_product(results%labFrameWeights,  &
        results%betaValsFit)/results%labFrameWeightSum
!
!     Print the summary.
!
      if(options%printResults) then
        write(iOut,3500) results%kMag
        do i = 1,nIntPlanes
          write(iOut,3510) TRIM(results%intPlaneLabels(i)),  &
            results%epsilonVector(:,i),results%integratedIntensity(i),  &
            results%solidAngleIntegratedIntensity(i),  &
            results%betaValsParaPerp(i),results%betaValsFit(i),  &
            results%rSquared(i)
        endDo
        write(iOut,3520)
        write(iOut,3530) results%meanBetaParaPerp,results%meanBetaFit
        write(iOut,3540) results%averageIntegratedIntensity,  &
          results%averageSolidAngleIntegratedIntensity,  &
          results%averageBetaParaPerp,results%averageBetaFit,  &
          results%averageRSquared
      endIf
!
      return
      end subroutine runPADCalculation


!PROCEDURE buildPADLabFrames
      subroutine buildPADLabFrames(options,epsilonVectors,kPlaneVectors,  &
        planeLabels,labFrameWeights)
!
!     This routine builds the lab-frame polarization and reference k-plane
!     vector arrays requested in the PAD options object. The calculation driver
!     only loops over these arrays and does not need to know how they were
!     generated.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      type(pad_options),intent(in)::options
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        epsilonVectors,kPlaneVectors
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        labFrameWeights
      character(len=8),dimension(:),allocatable,intent(out)::planeLabels
!
      integer(kind=int64)::i,nLabFrames
!
!     Build the requested lab-frame vector set.
!
      select case(options%labFrameType)
      case(PAD_LAB_FRAMES_CARTESIAN)
        call buildCartesianLabFrames(epsilonVectors,kPlaneVectors,  &
          planeLabels,labFrameWeights)
!
      case(PAD_LAB_FRAMES_SPHERE)
        call buildSphereLabFrames(epsilonVectors,kPlaneVectors,  &
          planeLabels,labFrameWeights,options%nLabFrameTheta,  &
          options%nLabFramePhi)
!
      case(PAD_LAB_FRAMES_CUSTOM)
        if(.not.Allocated(options%labEpsilonVector))  &
          call mqc_error('buildPADLabFrames: custom epsilon vectors missing.')
        if(.not.Allocated(options%labKPlaneVector))  &
          call mqc_error('buildPADLabFrames: custom k-plane vectors missing.')
        if(Size(options%labEpsilonVector,1).ne.3)  &
          call mqc_error('buildPADLabFrames: epsilon vectors must be 3 x N.')
        if(Size(options%labKPlaneVector,1).ne.3)  &
          call mqc_error('buildPADLabFrames: k-plane vectors must be 3 x N.')
        if(Size(options%labEpsilonVector,2).ne.  &
          Size(options%labKPlaneVector,2))  &
          call mqc_error('buildPADLabFrames: custom vector counts differ.')
        nLabFrames = Size(options%labEpsilonVector,2)
        if(nLabFrames.lt.1)  &
          call mqc_error('buildPADLabFrames: no custom lab frames supplied.')
        epsilonVectors = options%labEpsilonVector
        kPlaneVectors = options%labKPlaneVector
        if(Allocated(options%labFrameLabels)) then
          if(Size(options%labFrameLabels).ne.nLabFrames)  &
            call mqc_error('buildPADLabFrames: custom label count differs.')
          planeLabels = options%labFrameLabels
        else
          Allocate(planeLabels(nLabFrames))
          do i = 1,nLabFrames
            write(planeLabels(i),'("c",i7.7)') i
          endDo
        endIf
        if(Allocated(options%labFrameWeights)) then
          if(Size(options%labFrameWeights).ne.nLabFrames)  &
            call mqc_error('buildPADLabFrames: custom weight count differs.')
          if(MINVAL(options%labFrameWeights).lt.mqc_float(0))  &
            call mqc_error('buildPADLabFrames: negative custom weight.')
          labFrameWeights = options%labFrameWeights
        else
          Allocate(labFrameWeights(nLabFrames))
          labFrameWeights = mqc_float(1)
        endIf
!
      case default
        call mqc_error('buildPADLabFrames: invalid lab-frame type.')
      end select
!
      return
      end subroutine buildPADLabFrames


!PROCEDURE buildCartesianLabFrames
      subroutine buildCartesianLabFrames(epsilonVectors,kPlaneVectors,  &
        planeLabels,labFrameWeights)
!
!     Build the default fixed-orientation lab-frame set. The molecule stays in
!     the input molecular frame; each column gives a lab electric-field
!     polarization vector and one perpendicular reference vector used to define
!     the chi-rotated k-vector scan plane.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        epsilonVectors,kPlaneVectors
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        labFrameWeights
      character(len=8),dimension(:),allocatable,intent(out)::planeLabels
!
      Allocate(epsilonVectors(3,3),kPlaneVectors(3,3),planeLabels(3),  &
        labFrameWeights(3))
      labFrameWeights = mqc_float(1)
!
      planeLabels(1) = 'xy'
      epsilonVectors(:,1) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      kPlaneVectors(:,1)  = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
!
      planeLabels(2) = 'yz'
      epsilonVectors(:,2) = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
      kPlaneVectors(:,2)  = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
!
      planeLabels(3) = 'zx'
      epsilonVectors(:,3) = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
      kPlaneVectors(:,3)  = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
!
      return
      end subroutine buildCartesianLabFrames


!PROCEDURE buildSphereLabFrames
      subroutine buildSphereLabFrames(epsilonVectors,kPlaneVectors,  &
        planeLabels,labFrameWeights,nLabFrameTheta,nLabFramePhi)
!
!     This routine builds a lab-frame vector set from points on the unit sphere.
!     The sphere points are used as electric-field polarization vectors, and the
!     associated orthogonal vectors define reference directions for the
!     chi-rotated k-vector scan planes.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        epsilonVectors,kPlaneVectors
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        labFrameWeights
      character(len=8),dimension(:),allocatable,intent(out)::planeLabels
      integer(kind=int64),intent(in)::nLabFrameTheta,nLabFramePhi
!
      integer(kind=int64)::i,nLabFrames
!
!     Build the vectors and labels.
!
      if(nLabFrameTheta.lt.3)  &
        call mqc_error('buildSphereLabFrames: theta grid must be at least 3.')
      call buildSphereGrid(epsilonVectors,kPlaneVectors,nLabFrameTheta,  &
        nLabFramePhi,weights=labFrameWeights)
      nLabFrames = Size(epsilonVectors,2)
      Allocate(planeLabels(nLabFrames))
      do i = 1,nLabFrames
        write(planeLabels(i),'("s",i7.7)') i
      endDo
!
      return
      end subroutine buildSphereLabFrames


!PROCEDURE buildChiQuadrature
      subroutine buildChiQuadrature(chiVals,chiWeights,nChi)
!
!     This routine builds a periodic uniform quadrature over chi from 0 to
!     2*pi. The point at 2*pi is omitted because it is equivalent to the point
!     at 0 for periodic sampling.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::nChi
      real(kind=real64),dimension(:),allocatable,intent(out)::chiVals,  &
        chiWeights
!
      integer(kind=int64)::i
      real(kind=real64)::chiStep
!
      if(nChi.lt.1) call mqc_error('buildChiQuadrature: nChi must be positive.')
      chiStep = mqc_float(2)*Pi/mqc_float(nChi)
      Allocate(chiVals(nChi),chiWeights(nChi))
      do i = 1,nChi
        chiVals(i) = mqc_float(i-1)*chiStep
        chiWeights(i) = chiStep
      endDo
!
      return
      end subroutine buildChiQuadrature


!PROCEDURE buildTransverseBasis
      subroutine buildTransverseBasis(epsilonVector,kPlaneReference,  &
        epsilonUnit,uBasis,vBasis)
!
!     This routine normalizes the input polarization vector and constructs an
!     orthonormal transverse basis around it. The input reference k-plane vector
!     is projected into the transverse space when possible and is otherwise
!     replaced by a simple Cartesian fallback direction.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(3),intent(in)::epsilonVector,  &
        kPlaneReference
      real(kind=real64),dimension(3),intent(out)::epsilonUnit,uBasis,vBasis
!
      real(kind=real64),dimension(3)::fallbackVector
!
      epsilonUnit = epsilonVector
      if(dot_product(epsilonUnit,epsilonUnit).le.mqc_small)  &
        call mqc_error('buildTransverseBasis: epsilon vector has zero norm.')
      call mqc_normalizeVector(epsilonUnit)
!
!     Start from the supplied reference vector and project out any component
!     parallel to epsilon.
!
      uBasis = kPlaneReference-dot_product(kPlaneReference,epsilonUnit)*  &
        epsilonUnit
      if(dot_product(uBasis,uBasis).le.mqc_small) then
        if(abs(epsilonUnit(1)).le.abs(epsilonUnit(2)).and.  &
          abs(epsilonUnit(1)).le.abs(epsilonUnit(3))) then
          fallbackVector = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
        elseIf(abs(epsilonUnit(2)).le.abs(epsilonUnit(3))) then
          fallbackVector = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
        else
          fallbackVector = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
        endIf
        uBasis = fallbackVector-dot_product(fallbackVector,epsilonUnit)*  &
          epsilonUnit
      endIf
      if(dot_product(uBasis,uBasis).le.mqc_small)  &
        call mqc_error('buildTransverseBasis: could not form transverse basis.')
      call mqc_normalizeVector(uBasis)
!
!     Complete the right-handed transverse basis and re-orthogonalize uBasis.
!
      vBasis = mqc_crossProduct3D_real(epsilonUnit,uBasis)
      if(dot_product(vBasis,vBasis).le.mqc_small)  &
        call mqc_error('buildTransverseBasis: degenerate transverse basis.')
      call mqc_normalizeVector(vBasis)
      uBasis = mqc_crossProduct3D_real(vBasis,epsilonUnit)
      call mqc_normalizeVector(uBasis)
!
      return
      end subroutine buildTransverseBasis



!PROCEDURE betaParaPerp
      function betaParaPerp(IPara,IPerp) result(beta)
!
!     This function evaluates the photoelectron anisotropy parameter, beta,
!     using the analytic form based on the angular intensities when the photon
!     electric field is parallel and perpendicular to the detected detachment
!     direction.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),intent(in)::IPara,IPerp
      real(kind=real64)::beta
!
!     Do the work...
!
      beta = IPara-IPerp
      if(abs(beta).lt.mqc_small) then
        beta = mqc_float(0)
      else
        beta = beta/((IPara/mqc_float(2))+IPerp)
      endIf
!
      return
      end function betaParaPerp

!
!PROCEDURE betaLeastSquares
      subroutine betaLeastSquares(MSquaredList,thetaList,beta,rSquared)
!
!     This subroutine is given a list of intensities as a function of theta and
!     a list of theta values. Then, using a least squares fitting scheme, the
!     routine returns the fitted value of beta.
!
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::MSquaredList,thetaList
      real(kind=real64)::beta,rSquared
!
      integer(kind=int64)::n
      real(kind=real64)::slope,intercept
      real(kind=real64),dimension(:),allocatable::x,y
!
!     Before doing any work, check to see if all intensities are zero. If they
!     are, then beta is zero.
!
      beta = mqc_float(0)
      rSquared = mqc_float(0)
      if(maxval(MSquaredList).lt.mqc_small) return
      if((maxval(MSquaredList)-minval(MSquaredList)).lt.  &
        mqc_small*max(mqc_float(1),maxval(abs(MSquaredList)))) then
        rSquared = mqc_float(1)
        return
      endIf
!
!     Set n and then allocate the temporary arrays.
!
      n = SIZE(MSquaredList)
      Allocate(x(n),y(n))
!
!     Fill the x array by linearizing the anisotropy intensity equation's
!     independent variable.
!
      x = mqc_float(3)*(Cos(thetaList))**2-mqc_float(1)
      x = x/mqc_float(2)
      y = MSquaredList/maxval(MSquaredList)
!
!     Call the least squares fitting routine to get slope and intercept values.
!     Then, compuate the value of beta.
!
      call MQC_leastSquaresFit(x,y,slope,intercept,rSquared)
      beta = slope/intercept
!
      DeAllocate(x)
      return
      end subroutine betaLeastSquares
!
!PROCEDURE buildSphereGrid
      subroutine buildSphereGrid(xyz,orthogVector,nThetaInput,nPhiInput,  &
        weights)
!
!     This routine builds a set of Cartesian coordinates on the surface of a
!     unit sphere. The points are evenly spaced in theta and phi. The output
!     from this routine is array xyz, which is allocatable. The output array
!     orthogVector gives the projection onto the xy plance of each point. At the
!     poles, orthogVector is set to the positive x-direction.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),dimension(:,:),allocatable::xyz,orthogVector
      integer(kind=int64),intent(in),optional::nThetaInput,nPhiInput
      real(kind=real64),dimension(:),allocatable,intent(out),optional::weights
      integer(kind=int64)::i,j,nThetaPoints,nPhiPoints,ixyz
      real(kind=real64)::thetaStep,phiStep,theta,phi,x,y,z,  &
        thetaLower,thetaUpper,pointWeight
!
 1000 format(1x,i4,': theta,phi (pi)=',f5.2,',',f5.2,' | x,y,z=',f8.4,',',f8.4,',',f8.4)
!
!     Start by setting nThetaPoints and  nPhiPoints, and figuring out thetaStep
!     and phiStep. Then, allocate xyz.
!
      nThetaPoints = 5
      nPhiPoints = 8
      if(PRESENT(nThetaInput)) nThetaPoints = nThetaInput
      if(PRESENT(nPhiInput)) nPhiPoints = nPhiInput
      if(nThetaPoints.lt.2)  &
        call mqc_error('buildSphereGrid: nThetaPoints must be at least 2.')
      if(nPhiPoints.lt.1)  &
        call mqc_error('buildSphereGrid: nPhiPoints must be positive.')
      thetaStep = Pi/mqc_float(nThetaPoints-1)
      phiStep = mqc_float(2)*Pi/mqc_float(nPhiPoints)
!
!     Loop over theta and phi to build the array xyz.
!
      Allocate(xyz(3,(nThetaPoints-2)*nPhiPoints+2),  &
        orthogVector(3,(nThetaPoints-2)*nPhiPoints+2))
      if(PRESENT(weights))  &
        Allocate(weights((nThetaPoints-2)*nPhiPoints+2))
      theta = 0
      ixyz = 1
      do i = 1,nThetaPoints
        theta = mqc_float(i-1)*thetaStep
        thetaLower = max(mqc_float(0),theta-(thetaStep/mqc_float(2)))
        thetaUpper = min(Pi,theta+(thetaStep/mqc_float(2)))
        if(i.eq.1) then
          phi = mqc_float(0)
          x = sin(theta)*cos(phi)
          y = sin(theta)*sin(phi)
          z = cos(theta)
          xyz(1,ixyz) = x
          xyz(2,ixyz) = y
          xyz(3,ixyz) = z
          orthogVector(1,ixyz) = mqc_float(1)
          orthogVector(2,ixyz) = mqc_float(0)
          orthogVector(3,ixyz) = mqc_float(0)
          if(PRESENT(weights)) then
            pointWeight = mqc_float(2)*Pi*(cos(thetaLower)-cos(thetaUpper))
            weights(ixyz) = pointWeight
          endIf
          ixyz = ixyz+1
        elseIf(i.eq.nThetaPoints) then
          phi = mqc_float(0)
          x = sin(theta)*cos(phi)
          y = sin(theta)*sin(phi)
          z = cos(theta)
          xyz(1,ixyz) = x
          xyz(2,ixyz) = y
          xyz(3,ixyz) = z
          orthogVector(1,ixyz) = mqc_float(-1)
          orthogVector(2,ixyz) = mqc_float(0)
          orthogVector(3,ixyz) = mqc_float(0)
          if(PRESENT(weights)) then
            pointWeight = mqc_float(2)*Pi*(cos(thetaLower)-cos(thetaUpper))
            weights(ixyz) = pointWeight
          endIf
          ixyz = ixyz+1
        else
          do j = 1,nPhiPoints
            phi = mqc_float(j-1)*phiStep
            x = sin(theta)*cos(phi)
            y = sin(theta)*sin(phi)
            z = cos(theta)
            xyz(1,ixyz) = x
            xyz(2,ixyz) = y
            xyz(3,ixyz) = z
            x = sin(theta+0.5*pi)*cos(phi)
            y = sin(theta+0.5*pi)*sin(phi)
            z = cos(theta+0.5*pi)
            orthogVector(1,ixyz) = x
            orthogVector(2,ixyz) = y
            orthogVector(3,ixyz) = z
            if(PRESENT(weights)) then
              pointWeight = phiStep*(cos(thetaLower)-cos(thetaUpper))
              weights(ixyz) = pointWeight
            endIf
            ixyz = ixyz+1
          endDo
        endIf
      endDo
!
      return
      end subroutine buildSphereGrid

!
!PROCEDURE compute_legendre
      subroutine compute_legendre(l,m,x,Plm)
!
!     Evaluates the associated Legendre polynomial P_l^m(x) using upward
!     recursion from P_m^m(x).
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      integer(kind=int64),intent(in)::l,m
      real(kind=real64),intent(in)::x
      real(kind=real64),intent(out)::Plm
!
      integer(kind=int64)::i,abs_m
      real(kind=real64)::pmm,pmmp1,pll,somx2
!
!     Check for valid input range.
!
      if(l.lt.0.or.abs(m).gt.l) call mqc_error('compute_legendre: Invalid l or m.')
!
!     Compute P_m^m(x)
!
      abs_m = abs(m)
      somx2 = sqrt(max(mqc_float(0),mqc_float(1)-x*x))
      pmm = mqc_float(1)
      if(abs_m.gt.0) then
        pmm = (-mqc_float(1))**abs_m
        do i=1,abs_m
          pmm = pmm*somx2*mqc_float(2*i-1)
        endDo
      endIf
      if(l.eq.abs_m) then
        Plm = pmm
        return
      endIf
!
!     Compute P_{m+1}^m(x)
!
      pmmp1 = x*pmm*mqc_float(2*abs_m+1)
      if(l.eq.abs_m + 1) then
        Plm = pmmp1
        return
      endIf
!
!     Upward recurrence for l > m + 1
!
      do i=abs_m+2,l
        pll = ((2*i-1)*x*pmmp1-(i+abs_m-1)*pmm)/(i-abs_m)
        pmm = pmmp1
        pmmp1 = pll
      endDo
      Plm = pll
!
!     Apply Condon-Shortley phase if m < 0
!
      if(m.lt.0) then
        Plm = (-mqc_float(1))**abs_m * Plm
      endIf
!
      return
      end subroutine compute_legendre

!
!PROCEDURE Ylm_complex
      function Ylm_complex(l,m,theta,phi) result(Y)
!
!     Evaluates the complex spherical harmonic Y_l^m(theta,phi) using associated
!     Legendre polynomials and exp(i*m*phi).
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      integer(kind=int64),intent(in)::l,m
      real(kind=real64),intent(in)::theta,phi
      complex(kind=real64)::Y
!
      real(kind=real64)::P_lm,x,norm
      complex(kind=real64)::eimphi
      integer(kind=int64)::abs_m
!
!     Compute the key input values.
!
      x = cos(theta)
      abs_m = abs(m)
!
!     Compute associated Legendre polynomial.
!
      call compute_legendre(l,m,x,P_lm)
!
!     Figure out the normalization factor.
!
      norm = sqrt((2*l + 1)/(mqc_float(4)*acos(-mqc_float(1))) *  &
        mqc_float(factorial(l-abs_m))/mqc_float(factorial(l+abs_m)))
!
!     Compute e^(i m phi)
!
      eimphi = cmplx(mqc_float(0),mqc_float(m)*phi,kind=real64)
      eimphi = exp(eimphi)
!
!     Final result
!
      Y = norm*P_lm*eimphi
!
      return
      end function Ylm_complex

!
!PROCEDURE sph_bessel_j
      function sph_bessel_j(l,x) result(jl)
!
!     Computes the spherical Bessel function j_l(x) using direct expressions for
!     l=0,1 and recurrence for l>1.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      integer(kind=int64),intent(in)::l
      real(kind=real64),intent(in)::x
      real(kind=real64)::jl
!
      integer(kind=int64)::n
      real(kind=real64)::j0,j1,jn
!
!     Handle small-x regime explicitly...
!
      if(abs(x).lt.mqc_small) then
        if(l.eq.0) then
          jl = mqc_float(1)
        else
          jl = mqc_float(0)
        endIf
        return
      endIf
!
!     Handle the base cases...
!
      if(l.eq.0) then
        jl = sin(x)/x
        return
      endIf
      if(l.eq.1) then
        jl = (sin(x)/x**2)-(cos(x)/x)
        return
      endIf
!
!     Recursion for l > 1...
!
      j0 = sin(x)/x
      j1 = (sin(x)/x**2)-(cos(x)/x)
!
      do n=1,l-1
        jn = ((2*n+1)/x)*j1-j0
        j0 = j1
        j1 = jn
      endDo
      jl = jn
!
      return
      end function sph_bessel_j

!PROCEDURE dysonMatrixElement1Angle
subroutine dysonMatrixElement1Angle(iPEType,lMax,theta,kMag,  &
  photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,  &
  quadraturePoints,quadratureWeights,MSquared,lWeights)
!
  implicit none
  integer(kind=int64),intent(in) :: iPEType,lMax
  real(kind=real64),intent(in) :: theta,kMag
  real(kind=real64),dimension(3),intent(inOut) :: photonVector, orthogPlaneVector
  real(kind=real64),dimension(:),intent(in) :: dysonCoeffs, quadratureWeights
  real(kind=real64),dimension(:,:),intent(in) :: quadraturePoints
  real(kind=real64),intent(out) :: MSquared
  real(kind=real64),dimension(0:),intent(out) :: lWeights
  class(mqc_basisSet),intent(in) :: aoBasisSet
!
  integer(kind=int64) :: i, l, m
  real(kind=real64) :: rVec(3), r, thetaVal, phiVal, w
  real(kind=real64) :: j_l, dysonVal, dysonNorm, epsilonDotR
  complex(kind=real64) :: Ylm, muVal
  complex(kind=real64) :: dipoleAmp
  complex(kind=real64),dimension(0:lMax,-lMax:lMax) :: cLM
  real(kind=real64),dimension(0:lMax) :: W_l
  real(kind=real64),dimension(:),allocatable :: aoBasisValues,        &
       MValuesReal, MValuesImaginary, dysonNormTest
  real(kind=real64),dimension(3) :: kVector
  complex(kind=real64) :: iunit
!
  iunit = (0.0_real64, 1.0_real64)
  call mqc_normalizeVector(photonVector)
  call mqc_normalizeVector(orthogPlaneVector)
  kVector = cos(theta)*photonVector + sin(theta)*orthogPlaneVector
  kVector = kMag * kVector
!
  allocate(MValuesReal(size(quadratureWeights)),  &
           MValuesImaginary(size(quadratureWeights)),  &
           dysonNormTest(size(quadratureWeights)))
!
  cLM = (0.0, 0.0)
!
!$omp parallel do private(i, aoBasisValues, rVec, r, dysonVal, thetaVal, phiVal, &
!$omp& l, m, j_l, Ylm, muVal, epsilonDotR, w)                                     &
!$omp& shared(MValuesReal, MValuesImaginary, dysonNormTest, cLM) schedule(dynamic)
  do i = 1, size(quadratureWeights)
    rVec = quadraturePoints(:,i)
    r = sqrt(dot_product(rVec,rVec))
    call basisSetValuesList1(aoBasisSet, rVec, aoBasisValues)
    dysonVal = dot_product(dysonCoeffs, aoBasisValues)
    dysonNormTest(i) = dysonVal*dysonVal
!
    select case(iPEType)
    case(1)  ! Plane wave
      w = dot_product(kVector, rVec)
      epsilonDotR = dot_product(photonVector, rVec)
      MValuesReal(i) = cos(w) * epsilonDotR * dysonVal
      MValuesImaginary(i) = -sin(w) * epsilonDotR * dysonVal
!
    case(2)  ! Spherical Bessel
      if(r > MQC_Small) then
        thetaVal = acos(rVec(3)/r)
        phiVal = atan2(rVec(2),rVec(1))
      else
        thetaVal = 0.0_real64
        phiVal = 0.0_real64
      end if
!
      epsilonDotR = dot_product(photonVector, rVec)
!
      do l = 0, lMax
        j_l = sph_bessel_j(l, kMag * r)
        do m = -l, l
          Ylm = Ylm_complex(l, m, thetaVal, phiVal)
          muVal = j_l * Ylm * epsilonDotR * dysonVal
!$omp atomic
          cLM(l,m) = cLM(l,m) + quadratureWeights(i) * muVal
        end do
      end do
!
    case default
      call mqc_error('ERROR: Invalid iPEType in dysonMatrixElement1Angle.')
    end select
  end do
!$omp end parallel do
!
  dysonNorm = dot_product(quadratureWeights, dysonNormTest)
!
  if(iPEType == 2) then
!
!   Compute l-weights and MSquared
!
    do l = 0, lMax
      W_l(l) = 0.0_real64
      do m = -l, l
        W_l(l) = W_l(l) + abs(cLM(l,m))**2
      end do
    end do
!
    if(sum(W_l) > MQC_Small) then
      lWeights(0:lMax) = W_l / sum(W_l)
    else
      lWeights(0:lMax) = 0.0_real64
    end if
!
    dipoleAmp = (0.0, 0.0)
    do l = 0, lMax
      do m = -l, l
        dipoleAmp = dipoleAmp + cLM(l,m)
      end do
    end do
    MSquared = abs(dipoleAmp)**2
!
  else
    MSquared = dot_product(quadratureWeights, MValuesReal)**2 +  &
               dot_product(quadratureWeights, MValuesImaginary)**2
    lWeights(0) = 1.0_real64
    lWeights(1:lMax) = 0.0_real64
  end if
!
  return
end subroutine dysonMatrixElement1Angle

!
!PROCEDURE dysonMatrixElementThetaList
      subroutine dysonMatrixElementThetaList(iPEType,lMax,thetaVals,  &
        kMag,photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,  &
        quadraturePoints,quadratureWeights,Itheta,lWeights,  &
        lWeightsTheta)
!
!     Computes I(theta) at each angle in thetaVals using the specified
!     outgoing wave representation. Also returns optional partial-wave
!     angular momentum weights: intensity-averaged (lWeights) and
!     per-theta resolved (lWeightsTheta).
!
!     H. P. Hratchian, 2025.
!
      implicit none
!
      integer(kind=int64),intent(in)::iPEType,lMax
      real(kind=real64),dimension(:),intent(in)::thetaVals
      real(kind=real64),intent(in)::kMag
      real(kind=real64),dimension(3),intent(in)::photonVector,           &
        orthogPlaneVector
      real(kind=real64),dimension(:),intent(in)::dysonCoeffs,            &
        quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      real(kind=real64),dimension(:),intent(out)::Itheta
      real(kind=real64),dimension(0:),optional,intent(out)::lWeights
      real(kind=real64),dimension(0:,:),optional,intent(out)::lWeightsTheta
      class(mqc_basisSet),intent(in)::aoBasisSet
!
      integer(kind=int64)::iTh,nTh,l
      real(kind=real64)::theta,MSq,Itot
      real(kind=real64),dimension(0:lMax)::lWTemp
      real(kind=real64),dimension(0:lMax)::Wsum
      real(kind=real64),dimension(3)::pVec,oVec
!
!     Initialization
!
      nTh = size(thetaVals)
      pVec = photonVector
      oVec = orthogPlaneVector
      Wsum = mqc_float(0)
      Itot = mqc_float(0)
!
!     Loop through theta values
!
      do iTh=1,nTh
        theta = thetaVals(iTh)
        call dysonMatrixElement1Angle(iPEType,lMax,theta,kMag, pVec,  &
          oVec,dysonCoeffs,aoBasisSet,quadraturePoints,  &
          quadratureWeights,MSq,lWTemp)
        Itheta(iTh) = MSq
        Itot = Itot + MSq
        if(present(lWeightsTheta)) then
          do l=0,lMax
            lWeightsTheta(l,iTh) = lWTemp(l)
          end do
        end if
        if(present(lWeights)) then
          do l=0,lMax
            Wsum(l) = Wsum(l) + MSq * lWTemp(l)
          end do
        end if
      end do
!
      if(present(lWeights)) then
        if(Itot.gt.mqc_small) then
          do l=0,lMax
            lWeights(l) = Wsum(l) / Itot
          end do
        else
          do l=0,lMax
            lWeights(l) = mqc_float(0)
          end do
        end if
      end if
!
      return
      end subroutine dysonMatrixElementThetaList

!
!PROCEDURE generate_sph_grid
      subroutine generate_sph_grid(nTheta,nPhi,thetaVals,phiVals,  &
        weights)
      implicit none
      integer(kind=int64),intent(in)::nTheta,nPhi
      real(real64),intent(out)::thetaVals(nTheta),phiVals(nPhi)
      real(real64),intent(out)::weights(nTheta,nPhi)
      integer::i,j
      real(real64)::dtheta,dphi
!
!     Set dtheta and dphi.
!
      dtheta = pi/mqc_float(nTheta-1)
      dphi   = mqc_float(2)*pi/mqc_float(nPhi)
!
!     Create theta grid...
!
      do i = 1,nTheta
        thetaVals(i) = dtheta*mqc_float(i-1)
      end do
!
!     Create phi grid...
!
      do j = 1,nPhi
        phiVals(j) = dphi*mqc_float(j-1)
      end do
!
!     Compute weights: sin(theta)*dtheta*dphi...
!
      do i = 1,nTheta
        do j = 1,nPhi
          weights(i,j) = sin(thetaVals(i))*dtheta*dphi
        endDo
      endDo
!
      return
      end subroutine generate_sph_grid



      end module pad_mod
