      module gbs_mod
!
!     This module supports the GBS test program.
!
      use iso_fortran_env
      use mqc_general
      use mqc_integrals1
      use mqc_gaussian
!
      implicit none
      integer(kind=int64),parameter::iOut=6
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
!PROCEDURE dysonTransitionDipole
      function dysonTransitionDipole(theta,kMag,  &
        photonVector,dysonCoeffs,aoBasisSet,quadraturePoints,  &
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
      real(kind=real64),intent(in)::theta,kMag
      real(kind=real64),dimension(3),intent(in)::photonVector
      real(kind=real64),dimension(:),intent(in)::dysonCoeffs,quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      real(kind=real64)::MSquared
      class(mqc_basisSet),intent(in)::aoBasisSet
!
      integer(kind=int64)::i
      real(kind=real64)::dysonVal,dysonNorm,epsilonDotMu,MReal,  &
        MImaginary,w
      real(kind=real64),dimension(3)::kVector
      real(kind=real64),dimension(:),allocatable::aoBasisValues,  &
        MValuesReal,MValuesImaginary,dysonNormTest
!
!     Set up the elements of kVector. For now, we only consider the k-vector in
!     the xz plane.
!
      kVector(1) = sin(theta)
      kVector(2) = mqc_float(0)
      kVector(3) = cos(theta)
      kVector = kMag*kVector
!
!     Loop through the quadrature points to evaluate integrand values.
!
      Allocate(MValuesReal(SIZE(quadratureWeights)),  &
        MValuesImaginary(SIZE(quadratureWeights)),  &
        dysonNormTest(SIZE(quadratureWeights)))
      do i = 1,SIZE(quadratureWeights)
        w = kVector(1)*quadraturePoints(1,i)  &
          + kVector(2)*quadraturePoints(2,i)  &
          + kVector(3)*quadraturePoints(3,i)
        epsilonDotMu = photonVector(1)*quadraturePoints(1,i)  &
          + photonVector(2)*quadraturePoints(2,i)  &
          + photonVector(3)*quadraturePoints(3,i)
        aoBasisValues = basisSetValuesList(aoBasisSet,  &
          quadraturePoints(:,i))
        dysonVal = dot_product(dysonCoeffs,aoBasisValues)
        dysonNormTest(i) = dysonVal*dysonVal
        MValuesReal(i) = cos(w)*epsilonDotMu*dysonVal
        MValuesImaginary(i) = -sin(w)*epsilonDotMu*dysonVal
      endDo
      dysonNorm = dot_product(quadratureWeights,dysonNormTest)
      MReal = dot_product(quadratureWeights,MValuesReal)
      MImaginary = dot_product(quadratureWeights,MValuesImaginary)
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Hrant - MReal      = ',MReal
      write(iOut,*)'         MImaginary = ',MImaginary
      MSquared = MReal**2 + MImaginary**2
      write(iOut,'(A,f6.4,3x,f20.12,3x,f20.12)')' Hrant - theta, MSquared, dysonNorm = ',theta,MSquared,dysonNorm
      return
      end function dysonTransitionDipole

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
      real(kind=real64)::integralValue
!
      integer(kind=int64)::i
      real(kind=real64),dimension(:),allocatable::aoBasisValues,valuesGrid
!
!     Loop through the quadrature points to evaluate integrand values.
!
      Allocate(valuesGrid(SIZE(quadratureWeights)))
      do i = 1,SIZE(quadratureWeights)
        aoBasisValues = basisSetValuesList(aoBasisSet,  &
          quadraturePoints(:,i))
        valuesGrid(i) = dot_product(moCoeffsBra,aoBasisValues)
        if(PRESENT(moCoeffsKet)) then
          valuesGrid(i) = valuesGrid(i)*dot_product(moCoeffsKet,aoBasisValues)
        else
          valuesGrid(i) = valuesGrid(i)*valuesGrid(i)
        endIf
      endDo
      integralValue = dot_product(quadratureWeights,valuesGrid)
      return
      end function moInnerProductNumericalIntegration

      end module gbs_mod
