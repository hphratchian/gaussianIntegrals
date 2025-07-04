Include "gbs_mod.f03"
      program pad
!
!     This test program evaluates the photoelectron angular distribution,
!     intensities as a function of theta, for a given MO extracted from a FAF.
!     There is 1 required command line argument and 3 optional command line
!     arguments:
!           1. FAF filename (required)
!           2. (alpha) MO number to serve as the Dyson orbital (optional;
!              default=1)
!           3. magnitude of the k-vector in eV0
!           4. number of spatial grid points (optional; default=101)
!           5. number of angles to evaluate from 0 --> pi (optional; default=15)
!
!     CODE STATUS: At this time, it seems that the program works numerically.
!     However, some tests result in very large PAD intensities. I'm not sure if
!     these are correct or if there's an error or bug in the code somewhere. I
!     still need to run that to ground. Additionally, the code's output is not
!     printed in a convenient way yet. Furthermore, there's a lot of room for
!     exploring quadrature choices and also parallel coding.
!
!
!     Hrant P. Hratchian, 2025.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      USE OMP_LIB
      USE gbs_mod
      implicit none
      integer(kind=int64),parameter::nOMP=12,nIntPlanes=3
      logical,parameter::extraPrint=.false.
      integer(kind=int64)::i,j,iMODyson,nGridPointsM,nGridPointsTheta
      real(kind=real64)::tStart,tEnd,tstart1,tEnd1,stepSizeIntM,  &
        stepSizeTheta,thetaStart,moVal1,moVal2,kMag,MSquared0,  &
        MSquared90,betaTmp
      real(kind=real64),dimension(3)::cartStart,cartEnd
      real(kind=real64),dimension(nIntPlanes)::integratedIntensity,  &
        betaValsParaPerp,betaValsFit
      real(kind=real64),dimension(3,nIntPlanes)::laserVector,orthogPlaneVector
      real(kind=real64),dimension(:),allocatable::quadGridTheta,  &
        quadWeightsTheta,quadWeightsM,basisValues,MSquaredList
      real(kind=real64),dimension(:,:),allocatable::quadGridM,moCoeffs
      logical::fail=.false.
      character(len=256)::fafName
      character(len=2),dimension(nIntPlanes)::intPlaneLabels
      type(mqc_gaussian_unformatted_matrix_file)::faf
      type(mqc_basisset)::basisSet
      type(MQC_Variable)::tmp
!
      real(kind=real64),dimension(:,:),allocatable::basisIntegrals,  &
        overlapMatrix
      real(kind=real64),dimension(:),allocatable::quadValues
!
!     Format statements.
!
 1000 format(1x,'Program PAD.')
 1100 format(1x,'Current laser field: (',f6.2,',',f6.2,',',f6.2,')')
 2000 format(1x,'Data for the ',A,' grid:',/,  &
        3x,'nPoints = ',i15,' step size = ',f20.5)
 2500 format(/,1x,'Dyson orbital self-overlap = ',f12.8)
 3000 format(/,1x,55('-'),/,  &
        12x,'Intensity as a Function of Theta',/,  &
        18x,'(kMag = ',f8.3,' eV)',/,  &
        9x,'theta (rad)',15x,'I(theta)',/,  &
        1x,55('='),/)
 3010 format(9x,f7.3,3x,f25.8)
 3020 format(1x,55('='))
 3100 format(1x,'I(0) = ',f12.6,3x,'I(90) = ',f12.6,3x,'beta = ',f12.6)
 3500 format(/,1x,87('-'),/,  &
        31x,'Summary of PAD Calculation',/,  &
        34x,'(kMag = ',f8.3,' eV)',/,  &
        1x,87('='))
 3510 format(1x,'Integration Plane: ',A,'  |  Laser field: (',  &
        f5.2,',',f5.2,',',f5.2,')  |  Intensity = ',f10.6,/,  &
        20x,'beta(analytic) = ',f10.6,'  |  beta(fit) = ',f10.6)
 3520 format(1x,87('='))
 8998 format(1x,'Time for ',A,': ',f15.1,' s')
 8999 format(/,1x,'Job Time: ',f15.1,' s',/,1x,'PAD Complete.')
!
!
!     Begin the program.
!
      call CPU_TIME(tStart)
      call omp_set_num_threads(nOMP)
      fail = .false.
      write(iOut,1000)
      call mqc_version_print(iOut)
!
!     Memory check...
!
      if(MEMChecks) call print_memory_usage(iOut,'At top of PAD.')
!
!     Read the FAF name and MO number (the one we use for the Dyson orbtial)
!     from the command line. If the MO number isn't included, we set it to zero
!     here and later default it to the HOMO.
!
      if(command_argument_count().lt.1.or.  &
        command_argument_count().gt.5)  &
        call mqc_error('PAD expects 1-5 command line arguments.')
      call get_command_argument(1,fafName)
      if(command_argument_count().ge.2) then
        call mqc_get_command_argument_integer(2,iMODyson)
      else
        iMODyson = 1
      endIf
      if(command_argument_count().ge.3) then
        call mqc_get_command_argument_real(3,kMag)
      else
        kMag = mqc_float(1)/mqc_float(500)
      endIf
      if(command_argument_count().ge.4) then
        call mqc_get_command_argument_integer(4,nGridPointsM)
      else
        nGridPointsM = 101
      endIf
      if(command_argument_count().ge.5) then
        call mqc_get_command_argument_integer(5,nGridPointsTheta)
      else
        nGridPointsTheta = 15
      endIf
!
!     Load the FAF and set the MO number if it wasn't provided on the command
!     line.
      call faf%load(fafName)
!
!     Set the laser electric field vector and the orthogonal vector defining the
!     integration plane to be used with each electric field. At present we
!     hardwire 3 experiments with the electric field vector and planar slice of
!     photoelectron transition dipole matrices:
!           1. eVector: (1,0,0)  |  plane: yx
!           2. eVector: (0,1,0)  |  plane: yz
!           3. eVector: (0,0,1)  |  plane: xz
!
      intPlaneLabels(1) = 'xy'
      laserVector(:,1) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      orthogPlaneVector(:,1) = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
!
      intPlaneLabels(2) = 'yz'
      laserVector(:,2) = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
      orthogPlaneVector(:,2) = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
!
      intPlaneLabels(3)      = 'zx'
      laserVector(:,3)       = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
      orthogPlaneVector(:,3) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
!
!     Read the basis set and MO coefficients from faf.
!
      call loadGaussianBasisSet(faf,basisSet)
      call faf%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmp)
      moCoeffs = tmp
!
!     Fill in the grid of theta points.
!
      if(MEMChecks) call print_memory_usage(iOut,'Before building grids.')
      thetaStart = mqc_float(0)
      stepSizeTheta = Pi/mqc_float(nGridPointsTheta-1)
      Allocate(quadGridTheta(nGridPointsTheta),  &
        quadWeightsTheta(nGridPointsTheta))
      call setup_quadrature_trapezoid1d(nGridPointsTheta,stepSizeTheta,  &
        thetaStart,quadGridTheta,quadWeightsTheta)
      write(iOut,2000) 'theta',nGridPointsTheta,stepSizeTheta
      flush(iOut)
!
!     Prepare the integration grid and quadrature weights for the M evaluations.
!     There is one M per theta.
!
      cartStart = [ -6.0,-6.0,-6.0 ]
      cartEnd = [ 6.0,6.0,6.0 ]
      stepSizeIntM = (cartEnd(1)-cartStart(1))/mqc_float(nGridPointsM-1)
      Allocate(quadGridM(3,nGridPointsM**3),quadWeightsM(nGridPointsM**3),  &
        quadValues(nGridPointsM**3))
      call CPU_TIME(tStart1)
      call setup_quadrature_trapezoid3d(nGridPointsM,stepSizeIntM,  &
        cartStart,quadGridM,quadWeightsM)
      write(iOut,2000) 'M quadrature',nGridPointsM**3,stepSizeIntM
      call CPU_TIME(tEnd1)
      write(iOut,8998) '3D grid setup',tEnd1-tStart1
      flush(iOut)
!
!     Memory check...
!
      if(MEMChecks) call print_memory_usage(iOut,'After building grids.')
!
!     Test that the chosen MO is normalized using quadrature.
!
      if(extraPrint)  &
        call mqc_print(moCoeffs(:,iMODyson),iOut,header='Dyson coefficients')
      call CPU_TIME(tStart1)
      write(iOut,2500)  &
        moInnerProductNumericalIntegration(moCoeffs(:,iMODyson),  &
        quadGridM,quadWeightsM,basisSet)
      call CPU_TIME(tEnd1)
      write(iOut,8998) 'norm test',tEnd1-tStart1
      flush(iOut)
!
!     Integrate the Dyson/plane wave matrix elements as a function of theta.
!
      do i = 1,nIntPlanes
        write(iOut,1100) laserVector(:,i)
        call CPU_TIME(tStart1)
        MSquared0 = dysonPlaneWaveMatrixElementSquared(mqc_float(0),kMag,  &
          laserVector(:,i),orthogPlaneVector(:,i),moCoeffs(:,iMODyson),  &
          basisSet,quadGridM,quadWeightsM)
        call CPU_TIME(tEnd1)
        write(iOut,8998) 'MSquared0 ',tEnd1-tStart1
        call CPU_TIME(tStart1)
        MSquared90 = dysonPlaneWaveMatrixElementSquared(Pi/mqc_float(2),  &
          kMag,laserVector(:,i),orthogPlaneVector(:,i),  &
          moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
        call CPU_TIME(tEnd1)
        write(iOut,8998) 'MSquared90',tEnd1-tStart1
        call CPU_TIME(tStart1)
        flush(iOut)
        MSquaredList = dysonPlaneWaveMatrixElementSquaredThetaList(  &
          quadGridTheta,kMag,  &
          laserVector(:,i),orthogPlaneVector(:,i),moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
        call CPU_TIME(tEnd1)
        write(iOut,8998) 'MSquared(theta) list',tEnd1-tStart1
        flush(iOut)
!
!       Compute beta using a least squares fitting scheme.
!
        call betaLeastSquares(MSquaredList,quadGridTheta,betaTmp)
        betaValsFit(i) = betaTmp
!
!       Print the intensity data table.
!
        write(iOut,3000) kMag/evPHartree
        do j = 1,nGridPointsTheta
          write(iOut,3010) quadGridTheta(j),MSquaredList(j)
        endDo
        write(iOut,3020)
        integratedIntensity(i) = dot_product(MSquaredList,quadWeightsTheta)
!
!       Evaluate the anisotropy parameter, beta, using the computed intensities
!       at 0 and 90 degrees.
!
        betaValsParaPerp(i) = betaParaPerp(MSquared0,MSquared90)
        write(iOut,3100) MSquared0,MSquared90,  &
          betaParaPerp(MSquared0,MSquared90)
        flush(iOut)
      endDo
!
!     Print out the integrated planar intensity for each laser vector. Then
!     compute an analytic beta.
!
      write(iOut,3500) kMag
      do i = 1,nIntPlanes
        write(iOut,3510) TRIM(intPlaneLabels(i)),laserVector(:,i),  &
          integratedIntensity(i),betaValsParaPerp(i),betaValsFit(i)
      endDo
      write(iOut,3520)
!
!     The end of the program.
!
  999 Continue
      call CPU_TIME(tEnd)
      write(iOut,8999) tEnd-tStart
      if(MEMChecks) call print_memory_usage(iOut,'End of PAD.')
      end program pad
