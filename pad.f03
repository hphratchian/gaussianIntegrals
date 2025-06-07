Include "gbs_mod.f03"
      program pad
!
!     This test program evaluates the photoelectron angular distribution,
!     intensities as a function of theta, for a given MO extracted from a FAF.
!     There is one command line argument, a FAF from a Gaussian job. The program
!     has a parameter iMODyson declared below. That number corresponds to the
!     (alpha) MO that is used in the program as the Dyson orbital.
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
      integer(kind=int64),parameter::nOMP=1,iMODyson=4
      integer(kind=int64)::i,j,nGridPointsM,nGridPointsTheta
      real(kind=real64)::tStart,tEnd,stepSizeIntM,stepSizeTheta,  &
        thetaStart,moVal1,moVal2,kMag,MSquared
      real(kind=real64),dimension(3)::cartStart,cartEnd,laserVector
      real(kind=real64),dimension(:),allocatable::quadGridTheta,  &
        quadWeightsTheta,quadWeightsM,basisValues
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
 1000 format(1x,'Program gbs.')
!
!
!     Begin the program.
!
      call CPU_TIME(tStart)
      call omp_set_num_threads(nOMP)
      fail = .false.
      write(iOut,1000)
!$OMP PARALLEL
PRINT *, "Thread number:", OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
!
!     Read the FAF name from the command line and load the object faf.
!
      call get_command_argument(1,fafName)
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
!     Prepare the integration grid and quadrature weights for the intensity
!     vector, I(theta).
!
      thetaStart = mqc_float(0)
      nGridPointsTheta = 50
      stepSizeTheta = Pi/mqc_float(nGridPointsTheta-1)
      write(iOut,*)' nGridPointsTheta = ',nGridPointsTheta
      write(iOut,*)' stepSizeTheta    = ',stepSizeTheta
      write(iOut,*)' thetaStart       = ',thetaStart
      Allocate(quadGridTheta(nGridPointsTheta),  &
        quadWeightsTheta(nGridPointsTheta))
      call setup_quadrature_trapezoid1d(nGridPointsTheta,stepSizeTheta,  &
        thetaStart,quadGridTheta,quadWeightsTheta)
      write(iOut,*)' max theta grid point: ',maxval(quadGridTheta)
!
!     Prepare the integration grid and quadrature weights for the M evaluations.
!     There is one M per theta.
!
      cartStart = [ -5.0,-5.0,-5.0 ]
      cartEnd = [ 5.0,5.0,5.0 ]
      nGridPointsM = 201
      stepSizeIntM = (cartEnd(1)-cartStart(1))/mqc_float(nGridPointsM-1)
      write(iOut,*)' nGridPointsM = ',nGridPointsM
      write(iOut,*)' stepSizeIntM = ',stepSizeIntM
      Allocate(quadGridM(3,nGridPointsM**3),quadWeightsM(nGridPointsM**3),  &
        quadValues(nGridPointsM**3))
      call setup_quadrature_trapezoid3d(nGridPointsM,stepSizeIntM,  &
        cartStart,quadGridM,quadWeightsM)
      write(iOut,*)' max grid point: ',maxval(quadGridM)
!
!     Test that the chosen MO is normalized using quadrature.
!
      write(iOut,*)
      write(iOut,*)' Test of <dyson|dyson> = ',  &
        moInnerProductNumericalIntegration(moCoeffs(:,iMODyson),  &
        quadGridM,quadWeightsM,basisSet)
      write(iOut,*)
!
!     Try out the Dyson Transition Dipole magnitude function. For now, we
!     calculate the value at 0, 90, and 180 degrees.
!
      kMag = mqc_float(1)/mqc_float(100)
      laserVector = [ 0.0,0.0,1.0 ]
      MSquared = dysonTransitionDipole(mqc_float(0),kMag,  &
        laserVector,moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
      MSquared = dysonTransitionDipole(Pi/mqc_float(2),kMag,  &
        laserVector,moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
      MSquared = dysonTransitionDipole(Pi,kMag,  &
        laserVector,moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
!
!     The end of the program.
!
  999 Continue
      call CPU_TIME(tEnd)
      write(iOut,*)' TIME = ',tEnd-tStart
      write(iOut,*)
      end program pad
