Include "gbs_mod.f03"
      program pad
!
!     This test program evaluates the photoelectron angular distribution,
!     intensities as a function of theta, for a given MO extracted from a FAF.
!     There are up to 2 command line arguments. The first command line argument
!     is a FAF from a Gaussian job, which is required. The second command line
!     arguement is the (alpha) MO number to use as the Dyson orbtial, which is
!     optional and defaults to MO 1.
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
      integer(kind=int64),parameter::nOMP=48
      integer(kind=int64)::i,j,iMODyson,nGridPointsM,nGridPointsTheta
      real(kind=real64)::tStart,tEnd,tstart1,tEnd1,stepSizeIntM,  &
        stepSizeTheta,thetaStart,moVal1,moVal2,kMag,MSquared0,MSquared90
      real(kind=real64),dimension(3)::cartStart,cartEnd,laserVector
      real(kind=real64),dimension(:),allocatable::quadGridTheta,  &
        quadWeightsTheta,quadWeightsM,basisValues,MSquaredList
      real(kind=real64),dimension(:,:),allocatable::quadGridM,moCoeffs
      logical::fail=.false.
      character(len=256)::fafName
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
 2000 format(1x,'Data for the ',A,' grid:',/,  &
        3x,'nPoints = ',i15,' step size = ',f20.5)
 3000 format(/,1x,55('-'),/,  &
        12x,'Intensity as a Function of Theta',/,  &
        18x,'(kMag = ',f8.3,' eV)',/,  &
        9x,'theta (rad)',15x,'I(theta)',/,  &
        1x,55('='))
 3010 format(9x,f7.3,3x,f25.8)
 3020 format(1x,55('='),/)
 3100 format(/,1x,'I(0) = ',f25.8,3x,'I(90) = ',f25.8,/,  &
        1x,'beta = ',f25.8,/)
 8999 format(1x,'Job Time: ',f15.1,' s')
!
!
!     Begin the program.
!
      call CPU_TIME(tStart)
      call omp_set_num_threads(nOMP)
      fail = .false.
      write(iOut,1000)

!hph+
!!$OMP PARALLEL
!PRINT *, "Thread number:", OMP_GET_THREAD_NUM()
!!$OMP END PARALLEL
!hph-

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
        command_argument_count().gt.2)  &
        call mqc_error('PAD expects 1 or 2 command line arguments.')
      call get_command_argument(1,fafName)
      iMODyson = 1
      if(command_argument_count().ge.2)  &
        call mqc_get_command_argument_integer(2,iMODyson)
!
!     Load the FAF and set the MO number if it wasn't provided on the command
!     line.
      call faf%load(fafName)
!
!     Set the laser polarization vector.
!
      laserVector = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
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
      nGridPointsTheta = 5
      stepSizeTheta = Pi/mqc_float(nGridPointsTheta-1)
      Allocate(quadGridTheta(nGridPointsTheta),  &
        quadWeightsTheta(nGridPointsTheta))
      call setup_quadrature_trapezoid1d(nGridPointsTheta,stepSizeTheta,  &
        thetaStart,quadGridTheta,quadWeightsTheta)
      write(iOut,2000) 'theta',nGridPointsTheta,stepSizeTheta
!
!     Prepare the integration grid and quadrature weights for the M evaluations.
!     There is one M per theta.
!
      cartStart = [ -6.0,-6.0,-6.0 ]
      cartEnd = [ 6.0,6.0,6.0 ]
      nGridPointsM = 501
      stepSizeIntM = (cartEnd(1)-cartStart(1))/mqc_float(nGridPointsM-1)
      Allocate(quadGridM(3,nGridPointsM**3),quadWeightsM(nGridPointsM**3),  &
        quadValues(nGridPointsM**3))
      call CPU_TIME(tStart1)
      call setup_quadrature_trapezoid3d(nGridPointsM,stepSizeIntM,  &
        cartStart,quadGridM,quadWeightsM)
      write(iOut,2000) 'M quadrature',nGridPointsM**3,stepSizeIntM
      call CPU_TIME(tEnd1)
      write(iOut,*)' Time for 3D grid setup = ',tEnd1-tStart1
!
!     Memory check...
!
      if(MEMChecks) call print_memory_usage(iOut,'After building grids.')
!
!     Test that the chosen MO is normalized using quadrature.
!
      call mqc_print(moCoeffs(:,iMODyson),iOut,header='Dyson coefficients')
      call CPU_TIME(tStart1)
      write(iOut,*)
      write(iOut,*)' Test of <dyson|dyson> = ',  &
        moInnerProductNumericalIntegration(moCoeffs(:,iMODyson),  &
        quadGridM,quadWeightsM,basisSet)
      write(iOut,*)
      call CPU_TIME(tEnd1)
      write(iOut,*)' Time for norm test = ',tEnd1-tStart1
!
!     Try out the Dyson Transition Dipole magnitude function. For now, we
!     calculate the value at 0, 90, and 180 degrees.
!
      kMag = mqc_float(1)/mqc_float(100)
      laserVector = [ 0.0,0.0,1.0 ]
      call CPU_TIME(tStart1)
      MSquared0 = dysonPlaneWaveMatrixElementSquared(mqc_float(0),kMag,  &
        laserVector,moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
      call CPU_TIME(tEnd1)
      write(iOut,*)' Time for MSquared0  = ',tEnd1-tStart1
      call CPU_TIME(tStart1)
      MSquared90 = dysonPlaneWaveMatrixElementSquared(Pi/mqc_float(2),  &
        kMag,laserVector,moCoeffs(:,iMODyson),basisSet,quadGridM,  &
        quadWeightsM)
      call CPU_TIME(tEnd1)
      write(iOut,*)' Time for MSquared90 = ',tEnd1-tStart1
!hph+
!      MSquaredList = dysonPlaneWaveMatrixElementSquaredThetaList(  &
!        [ mqc_float(0),Pi/mqc_float(4),Pi/mqc_float(2),mqc_float(3)*Pi/mqc_float(4),Pi  ],kMag,  &
!        laserVector,moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
      call CPU_TIME(tStart1)
      MSquaredList = dysonPlaneWaveMatrixElementSquaredThetaList(  &
        quadGridTheta,kMag,  &
        laserVector,moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
!hph-
      call CPU_TIME(tEnd1)
      write(iOut,*)' Time for MSquared(theta) list = ',tEnd1-tStart1
      write(iOUt,*)
!
!     Print the intensity data table.
!
      write(iOut,3000) kMag*evPHartree
      do i = 1,nGridPointsTheta
        write(iOut,3010) quadGridTheta(i),MSquaredList(i)
      endDo
      write(iOut,3020)
!
!     Evaluate the anisotropy parameter, beta, using the computed intensities at
!     0 and 90 degrees.
!
      write(iOut,3100) MSquared0,MSquared90,  &
        (MSquared0-MSquared90)/((MSquared0/mqc_float(2))+MSquared90)
!
!     The end of the program.
!
  999 Continue
      call CPU_TIME(tEnd)
      write(iOut,8999) tEnd-tStart
      if(MEMChecks) call print_memory_usage(iOut,'End of PAD.')
      end program pad
