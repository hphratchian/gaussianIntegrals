include "memory_utils.f03"
      module gbs_mod
!
!     This module supports the GBS test program.
!
      use iso_fortran_env
      use mqc_general
      use mqc_integrals1
      use mqc_gaussian
      use memory_utils
!
      implicit none
      integer(kind=int64),parameter::iOut=6
      logical::MEMChecks=.false.
!
!
      CONTAINS
!
!PROCEDURE
      subroutine loadGaussianBasisSet(faf,basisSet)
!
!     This routine reads a basis set from the Gaussian Fortran Array File (FAF)
!     sent as <faf> and loads it into argument <basisSet>. The basisSet argument
!     should be an MQC_BasisSet object.
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      type(mqc_gaussian_unformatted_matrix_file)::faf
      type(MQC_basisSet)::basisSet
      integer(kind=int64)::i,iCurrentPrim,nBasis,nShells,nPrims,  &
        nShellsFull
      integer(kind=int64),dimension(:),allocatable::shellToAtomMap,  &
        shellTypes,nPrimsPerShell
      real(kind=real64),dimension(:),allocatable::primitiveExponents,  &
        contractionCoefficients,contractionCoefficientsP,coordinates
      type(MQC_Variable)::tmp
!
 1000 format(/,1x,'Data in loadGaussianBasisSet',/,  &
        3x,'nBasis=',i6,3x,'nShells=',i6,3x,'nPrims=',i6,/)
!
!     Get the key parameters of the basis set on the FAF.
!
      nBasis = faf%getVal('nBasis')
      nShells = faf%getVal('nShlAO')
      nPrims = faf%getVal('nPrmAO')
      write(iOut,1000) nBasis,nShells,nPrims
!
!     Read in the shell info from faf so we can load basisSet.
!
      call faf%getArray('shell to atom map',mqcVarOut=tmp)
      shellToAtomMap = tmp
      call faf%getArray('shell types',mqcVarOut=tmp)
      shellTypes = tmp
      call faf%getArray('number of primitives per shell',mqcVarOut=tmp)
      nPrimsPerShell = tmp
      call faf%getArray('primitive exponents',mqcVarOut=tmp)
      primitiveExponents = tmp
      call faf%getArray('contraction coefficients',mqcVarOut=tmp)
      contractionCoefficients = tmp
      call faf%getArray('P(S=P) CONTRACTION COEFFICIENTS',mqcVarOut=tmp)
      contractionCoefficientsP = tmp
      call faf%getArray('coordinates of each shell',mqcVarOut=tmp)
      coordinates = tmp
!
!     Figure out the value of nShellsFull...
!
      nShellsFull = 0
      do i = 1,nShells
        if(shellTypes(i).ge.0) then
          nShellsFull = nShellsFull+1
        elseIf(shellTypes(i).eq.-1) then
          nShellsFull = nShellsFull+4
        else
          write(iOut,*)' shellTypes(i) = ',shellTypes(i)
          call mqc_error('Invalid shell type in loadGaussianBasisSet.')
        endIf
      endDo
      call basisSet%init(nShellsFull)
!
!     Initiate the basis set object and then fill it.
!
      iCurrentPrim = 1
      do i = 1,nShells
        if(shellTypes(i).ge.0) then
          nShellsFull = nShellsFull+1
          call MQC_basisSet_addShell(basisSet,shellTypes(i),  &
            coordinates((3*i-2):3*i),  &
            contractionCoefficients(iCurrentPrim:iCurrentPrim+nPrimsPerShell(i)-1),  &
            primitiveExponents(iCurrentPrim:iCurrentPrim+nPrimsPerShell(i)-1))
        elseIf(shellTypes(i).eq.-1) then
          nShellsFull = nShellsFull+1
          call MQC_basisSet_addShell(basisSet,0,coordinates((3*i-2):3*i),  &
            contractionCoefficients(iCurrentPrim:iCurrentPrim+nPrimsPerShell(i)-1),  &
            primitiveExponents(iCurrentPrim:iCurrentPrim+nPrimsPerShell(i)-1))
          nShellsFull = nShellsFull+3
          call MQC_basisSet_addShell(basisSet,1,coordinates((3*i-2):3*i),  &
            contractionCoefficientsP(iCurrentPrim:iCurrentPrim+nPrimsPerShell(i)-1),  &
            primitiveExponents(iCurrentPrim:iCurrentPrim+nPrimsPerShell(i)-1))
        else
          call mqc_error('Invalid shell type in loadGaussianBasisSet.')
        endIf
        iCurrentPrim = iCurrentPrim+nPrimsPerShell(i)
      endDo
!
      return
      end subroutine loadGaussianBasisSet

!
!PROCEDURE setup_quadrature_trapezoid1d
      subroutine setup_quadrature_trapezoid1d(nPoints,stepSize,origin,  &
        quadGrid,quadWeights)
!
!     This subroutine generates a 1D trapezoidal quadrature grid and weights
!     over a grid with uniform spacing.
!
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      integer(kind=int64),intent(in)::nPoints
      real(kind=real64),intent(in)::stepSize,origin
      real(kind=real64),dimension(nPoints),intent(out)::quadGrid,  &
        quadWeights
      integer(kind=int64)::i,idx
      real(kind=real64)::wx,wy,wz
      real(kind=real64)::x
!
      quadGrid(1) = origin
      quadWeights(1) = stepSize/mqc_float(2)
      do i = 2,nPoints-1
        quadGrid(i) = quadGrid(i-1)+stepSize
        quadWeights(i) = stepSize
      endDo
      quadGrid(nPoints) = quadGrid(nPoints-1)+stepSize
      quadWeights(nPoints) = stepSize/mqc_float(2)
!
      return
      end subroutine setup_quadrature_trapezoid1d

!
!PROCEDURE setup_quadrature_trapezoid3d
      subroutine setup_quadrature_trapezoid3d(nPoints,stepSize,origin,  &
        quadGrid,quadWeights)
!
!     This subroutine generates a 3D trapezoidal quadrature grid and weights
!     over a cube with uniform spacing.
!
!
!     H. P. Hratchian, 2025.
!
!
      implicit none
      integer(kind=int64),intent(in)::nPoints
      real(kind=real64),intent(in)::stepSize
      real(kind=real64),dimension(3),intent(in)::origin
      real(kind=real64),dimension(3,nPoints**3),intent(out)::quadGrid
      real(kind=real64),dimension(nPoints**3),intent(out)::quadWeights
      integer(kind=int64)::i,j,k,idx
      real(kind=real64)::wx,wy,wz
      real(kind=real64)::x,y,z
!
      idx = 0
      do k = 1, nPoints
        z = origin(3) + (k - 1) * stepSize
        wz = stepSize
        if (k == 1 .or. k == nPoints) wz = 0.5d0 * stepSize
        do j = 1, nPoints
          y = origin(2) + (j - 1) * stepSize
          wy = stepSize
          if (j == 1 .or. j == nPoints) wy = 0.5d0 * stepSize
          do i = 1, nPoints
            x = origin(1) + (i - 1) * stepSize
            wx = stepSize
            if (i == 1 .or. i == nPoints) wx = 0.5d0 * stepSize
            idx = idx + 1
            quadGrid(1,idx) = x
            quadGrid(2,idx) = y
            quadGrid(3,idx) = z
            quadWeights(idx) = wx * wy * wz
          end do
        end do
      end do
!
      return
      end subroutine setup_quadrature_trapezoid3d

!
!PROCEDURE dysonPlaneWaveMatrixElementSquared
      function dysonPlaneWaveMatrixElementSquared(theta,kMag,  &
        photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,  &
        quadraturePoints,quadratureWeights) result(MSquared)
!
!     This function computes the Dyson transition dipole squared for a given
!     photon electric field vector, <photonVector>, the angle between that
!     vector and the outgoing plane wave, <theta>, and the magnitude of the
!     outgoing plane wave, <kMag>. This function assumes the outgoing plane wave
!     travels in the xz plane.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),intent(in)::theta,kMag
      real(kind=real64),dimension(3),intent(inOut)::photonVector,  &
        orthogPlaneVector
      real(kind=real64),dimension(:),intent(in)::dysonCoeffs,quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      real(kind=real64)::MSquared
      class(mqc_basisSet),intent(in)::aoBasisSet
!
      integer(kind=int64)::i
      real(kind=real64)::dysonVal,dysonNorm,epsilonDotMu,kVectorA,  &
        kVectorB,MReal,MImaginary,w
      real(kind=real64),dimension(3)::kVector
      real(kind=real64),dimension(:),allocatable::aoBasisValues,  &
        MValuesReal,MValuesImaginary,dysonNormTest
!
!     Ensure the photon electric field vector and the orthogonal vector defining
!     the integration plane are both normalized. Then, set-up the k-vector.
!
      call mqc_normalizeVector(photonVector)
      call mqc_normalizeVector(orthogPlaneVector)
      kVector = cos(theta)*photonVector+  &
        sin(theta)*orthogPlaneVector
      kVector = kMag*kVector
!
!     Loop through the quadrature points to evaluate integrand values.
!
      Allocate(MValuesReal(SIZE(quadratureWeights)),  &
        MValuesImaginary(SIZE(quadratureWeights)),  &
        dysonNormTest(SIZE(quadratureWeights)))
!$omp parallel do private(i, aoBasisValues, w, epsilonDotMu, dysonVal)  &
!$omp&  shared(MValuesReal, MValuesImaginary, dysonNormTest)  &
!$omp&  schedule(dynamic)
      do i = 1,SIZE(quadratureWeights)
        w = kVector(1)*quadraturePoints(1,i)  &
          + kVector(2)*quadraturePoints(2,i)  &
          + kVector(3)*quadraturePoints(3,i)
        epsilonDotMu = photonVector(1)*quadraturePoints(1,i)  &
          + photonVector(2)*quadraturePoints(2,i)  &
          + photonVector(3)*quadraturePoints(3,i)
        call basisSetValuesList1(aoBasisSet,  &
          quadraturePoints(:,i),aoBasisValues)
        dysonVal = dot_product(dysonCoeffs,aoBasisValues)
        dysonNormTest(i) = dysonVal*dysonVal
        MValuesReal(i) = cos(w)*epsilonDotMu*dysonVal
        MValuesImaginary(i) = -sin(w)*epsilonDotMu*dysonVal
      endDo
!$omp end parallel do
      dysonNorm = dot_product(quadratureWeights,dysonNormTest)
      MReal = dot_product(quadratureWeights,MValuesReal)
      MImaginary = dot_product(quadratureWeights,MValuesImaginary)
      MSquared = MReal**2 + MImaginary**2
      if(abs(MSquared).lt.mqc_small) MSquared = mqc_float(0)
      return
      end function dysonPlaneWaveMatrixElementSquared

!
!PROCEDURE dysonPlaneWaveMatrixElementSquaredThetaList
      function dysonPlaneWaveMatrixElementSquaredThetaList(thetaList,  &
        kMag,photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,quadraturePoints,&
        quadratureWeights) result(MSquared)
!
!     This function computes the Dyson transition dipole squared for a given
!     photon electric field vector, <photonVector>, the angle between that
!     vector and the outgoing plane wave, <theta>, and the magnitude of the
!     outgoing plane wave, <kMag>. This function assumes the outgoing plane wave
!     travels in the xz plane.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),intent(in)::kMag
      real(kind=real64),dimension(3),intent(inOut)::photonVector,orthogPlaneVector
      real(kind=real64),dimension(:),intent(in)::thetaList,dysonCoeffs,quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      class(mqc_basisSet),intent(in)::aoBasisSet
      real(kind=real64),dimension(:),allocatable::MSquared
!
      integer(kind=int64)::i,j,nTheta,nGrid
      real(kind=real64)::kVectorA,kVectorB,w,dysonVal,epsilonDotMu,thetaTest
      real(kind=real64),dimension(:),allocatable::aoBasisValues,dysonNormTest
      real(kind=real64),dimension(:),allocatable::MReal,MImag
      real(kind=real64),dimension(3)::gridPoint,kVector
!
!     Allocate memory and initialize variables.
!
      nTheta = SIZE(thetaList)
      nGrid = SIZE(quadratureWeights)
      Allocate(MSquared(nTheta))
      Allocate(MReal(nTheta),MImag(nTheta))
      Allocate(dysonNormTest(nGrid))
      MSquared = mqc_float(0)
      MReal = mqc_float(0)
      MImag = mqc_float(0)
      call mqc_normalizeVector(photonVector)
      if(MEMChecks) call print_memory_usage(iOut,'dysonPlaneWaveMatrixElementSquaredThetaList before OMP loop.')
!
!     Ensure the photon electric field vector is normalize and also set up
!     orthogPlaneVector.
!
      if(dot_product(photonVector,photonVector).gt.MQC_Small)  &
        call mqc_normalizeVector(photonVector)
      if(dot_product(orthogPlaneVector,orthogPlaneVector).gt.MQC_Small)  &
        call mqc_normalizeVector(orthogPlaneVector)
!
!     Loop through the quadrature points to evaluate integrand values.
!
!hph !$omp parallel default(shared) private(i,j,w,dysonVal,epsilonDotMu,aoBasisValues,gridPoint,kVector) &
!hph !$omp& reduction(+:MReal, MImag, dysonNormTest)
      allocate(aoBasisValues(SIZE(dysonCoeffs)))
!hph !$omp do schedule(dynamic)
      do i = 1, nGrid
        gridPoint = quadraturePoints(:,i)
        call basisSetValuesList1(aoBasisSet, gridPoint, aoBasisValues)
        dysonVal = dot_product(dysonCoeffs, aoBasisValues)
        dysonNormTest(i) = dysonVal * dysonVal
        epsilonDotMu = dot_product(photonVector, gridPoint)
        do j = 1, nTheta
          kVector = cos(thetaList(j))*photonVector+  &
            sin(thetaList(j))*orthogPlaneVector
          thetaTest = vectorAngle(kVector,photonVector)
!hph          write(iOut,'(1x,"theta=",f8.4," | thetaTest=",f8.4)') thetaList(j),thetaTest
          if(abs(thetaList(j)-thetaTest).gt.0.001)  &
            write(iOut,'(4x,"PROBLEM  ",f6.3,5(",",f6.3))') kVector,photonVector
!hph write(iOut,'(A,5(3x,f10.3))')' angle = ',vectorAngle(kVector,photonVector),thetaList(j),kVectorA,kVectorB,dot_product(kVector,photonVector)
          w = kMag * dot_product(kVector, gridPoint)
          MReal(j) = MReal(j) + cos(w) * epsilonDotMu * dysonVal * quadratureWeights(i)
          MImag(j) = MImag(j) - sin(w) * epsilonDotMu * dysonVal * quadratureWeights(i)
        end do
      end do
!hph !$omp end do
      deallocate(aoBasisValues)
!hph !$omp end parallel
      MSquared = MReal**2 + MImag**2
      call mqc_vectorTrimZero(MSquared)
      if(MEMChecks) call print_memory_usage(iOut,'dysonPlaneWaveMatrixElementSquaredThetaList after OMP loop.')
!
      return
      end function dysonPlaneWaveMatrixElementSquaredThetaList

!
!PROCEDURE moInnerProductNumericalIntegration
      function moInnerProductNumericalIntegration(moCoeffsBra,  &
        quadraturePoints,quadratureWeights,aoBasisSet,moCoeffsKet)  &
        result(integralValue)
!
!     This function computes the inner-product of two MOs numerically using the
!     quadrature grid and weights sent as input dummy arguments. Argument
!     moCoeffsKet is optional; if it is NOT sent, the bra and ket are both taken
!     to be moCoeffsBra.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::moCoeffsBra,quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      class(mqc_basisSet),intent(in)::aoBasisSet
      real(kind=real64),dimension(:),intent(in),optional::moCoeffsKet
      real(kind=real64)::integralValue,localValue,minValue
!
      integer(kind=int64)::i
      real(kind=real64),dimension(:),allocatable::aoBasisValues,valuesGrid
!
!     Loop through the quadrature points to evaluate integrand values.
!
      Allocate(valuesGrid(SIZE(quadratureWeights)))
      if(MEMChecks) call print_memory_usage(iOut,'moInnerProductNumericalIntegration before OMP loop.')
!$omp parallel do private(i,aoBasisValues,localValue) shared(valuesGrid) schedule(dynamic)
      do i = 1,SIZE(quadratureWeights)
!        aoBasisValues = basisSetValuesList(aoBasisSet,  &
!          quadraturePoints(:,i))
        call basisSetValuesList1(aoBasisSet,  &
          quadraturePoints(:,i),aoBasisValues)
        localValue = dot_product(moCoeffsBra,aoBasisValues)
        if(PRESENT(moCoeffsKet)) then
          localValue = localValue*dot_product(moCoeffsKet,aoBasisValues)
        else
          localValue = localValue*localValue
        endIf
        valuesGrid(i) = localValue
      endDo
!$omp end parallel do
      minValue = MinVal(valuesGrid)
      integralValue = dot_product(quadratureWeights,valuesGrid)
      if(abs(integralValue).lt.mqc_small) integralValue = mqc_float(0)
      return
      end function moInnerProductNumericalIntegration

!
!PROCEDURE vectorAngle
      function vectorAngle(v1,v2,degrees) result(angle)
!
!     This function computes the angle between two vectors, v1 and v2. The dummy
!     argument degrees is optional, defaulting to .FALSE., and indicates if the
!     angle should be reported in degrees (.TRUE.) or radians (.FALSE.).
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::v1,v2
      logical,optional::degrees
      real(kind=real64)::angle
      real(kind=real64)::denominator
!
!     Do the work...
!
      denominator = vectorMagnitude(v1)*vectorMagnitude(v2)
      if(denominator.gt.MQC_Small) then
        angle = acos(dot_product(v1,v2)/denominator)
      else
        angle = mqc_float(0)
      endIf
      if(PRESENT(degrees)) then
        if(degrees) angle = mqc_float(180)*angle/Pi
      endIf
!
      return
      end function vectorAngle

!
!PROCEDURE vectorMagnitude
      function vectorMagnitude(v) result(vMagnitude)
!
!     This function computes the magnitude of vector v.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::v
      real(kind=real64)::vMagnitude
!
!     Do the work...
!
      vMagnitude = sqrt(dot_product(v,v))
      return
      end function vectorMagnitude

!
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
      subroutine betaLeastSquares(MSquaredList,thetaList,beta)
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
      real(kind=real64)::beta
!
      integer(kind=int64)::n
      real(kind=real64)::slope,intercept,rSquared
      real(kind=real64),dimension(:),allocatable::x,y
!
!     Before doing any work, check to see if all intensities are zero. If they
!     are, then beta is zero.
!
      beta = mqc_float(0)
      if(maxval(MSquaredList).lt.mqc_small) return
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


      end module gbs_mod
