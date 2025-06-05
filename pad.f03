Include "gbs_mod.f03"
      program pad
!
!     This program evaluates the photoelectron angular distribution.
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
      integer(kind=int64),parameter::nOMP=1
      integer(kind=int64)::i,j,nGridPointsM,nGridPointsTheta
      real(kind=real64)::tStart,tEnd,stepSizeIntM,stepSizeTheta,  &
        thetaStart,moVal1,moVal2,totalIntegral
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
!     The end of the program.
!
  999 Continue
      call CPU_TIME(tEnd)
      write(iOut,*)' TIME = ',tEnd-tStart
      write(iOut,*)
      end program pad
