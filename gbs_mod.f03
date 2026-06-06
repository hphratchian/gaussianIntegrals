      module gbs_mod
!
!     This module supports Gaussian basis-set loading, simple quadrature grids,
!     and numerical MO integration used by the GBS and PAD programs.
!
      use iso_fortran_env
      use mqc_general
      use mqc_integrals1
      use mqc_gaussian
      use memory_utils
!
      implicit none
      integer(kind=int64),parameter::iOut=6
      integer(kind=int64),parameter::QUAD_RADIAL_TRAPEZOID=0_int64
      integer(kind=int64),parameter::QUAD_RADIAL_GAUSS_LEGENDRE=1_int64
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
!PROCEDURE setup_quadrature_gauss_legendre1d
      subroutine setup_quadrature_gauss_legendre1d(nPoints,lowerBound,  &
        upperBound,quadGrid,quadWeights)
!
!     This subroutine generates a Gauss-Legendre quadrature grid and weights
!     on the interval [lowerBound,upperBound].
!
!
!     H. P. Hratchian, 2026.
!
!
      implicit none
      integer(kind=int64),intent(in)::nPoints
      real(kind=real64),intent(in)::lowerBound,upperBound
      real(kind=real64),dimension(nPoints),intent(out)::quadGrid,  &
        quadWeights
      integer(kind=int64)::i,j,m
      real(kind=real64)::p1,p2,p3,pp,tol,xMid,xHalf,z,zOld
!
      if(nPoints.lt.1)  &
        call mqc_error('setup_quadrature_gauss_legendre1d: nPoints must be positive.')
      if(upperBound.le.lowerBound)  &
        call mqc_error('setup_quadrature_gauss_legendre1d: invalid interval.')
!
      tol = 1.0e-14_real64
      m = (nPoints+1)/2
      xMid = (upperBound+lowerBound)/mqc_float(2)
      xHalf = (upperBound-lowerBound)/mqc_float(2)
      do i = 1,m
        z = cos(Pi*(mqc_float(i)-mqc_float(1)/mqc_float(4))/  &
          (mqc_float(nPoints)+mqc_float(1)/mqc_float(2)))
        do
          p1 = mqc_float(1)
          p2 = mqc_float(0)
          do j = 1,nPoints
            p3 = p2
            p2 = p1
            p1 = ((mqc_float(2*j-1)*z*p2)-mqc_float(j-1)*p3)/mqc_float(j)
          endDo
          pp = mqc_float(nPoints)*(z*p1-p2)/(z*z-mqc_float(1))
          zOld = z
          z = zOld-p1/pp
          if(abs(z-zOld).le.tol) exit
        endDo
        quadGrid(i) = xMid-xHalf*z
        quadGrid(nPoints+1-i) = xMid+xHalf*z
        quadWeights(i) = mqc_float(2)*xHalf/  &
          ((mqc_float(1)-z*z)*pp*pp)
        quadWeights(nPoints+1-i) = quadWeights(i)
      endDo
!
      return
      end subroutine setup_quadrature_gauss_legendre1d


!PROCEDURE setup_lebedev_angular_grid
      subroutine setup_lebedev_angular_grid(lebedevOrder,angularGrid,  &
        angularWeights)
!
!     This routine builds selected low-order Lebedev angular grids on the unit
!     sphere. The returned angular weights include the 4*pi solid-angle factor.
!
!
!     H. P. Hratchian, 2026.
!
!
      implicit none
      integer(kind=int64),intent(in)::lebedevOrder
      real(kind=real64),dimension(:,:),allocatable,intent(out)::angularGrid
      real(kind=real64),dimension(:),allocatable,intent(out)::angularWeights
      integer(kind=int64)::idx
      real(kind=real64)::a,wAxis,wEdge,wCube
!
      select case(lebedevOrder)
      case(6)
        Allocate(angularGrid(3,6),angularWeights(6))
        wAxis = mqc_float(4)*Pi/mqc_float(6)
        idx = 0
        call addLebedevPoint(angularGrid,angularWeights,idx,  &
          [ mqc_float(1),mqc_float(0),mqc_float(0) ],wAxis)
        call addLebedevPoint(angularGrid,angularWeights,idx,  &
          [ -mqc_float(1),mqc_float(0),mqc_float(0) ],wAxis)
        call addLebedevPoint(angularGrid,angularWeights,idx,  &
          [ mqc_float(0),mqc_float(1),mqc_float(0) ],wAxis)
        call addLebedevPoint(angularGrid,angularWeights,idx,  &
          [ mqc_float(0),-mqc_float(1),mqc_float(0) ],wAxis)
        call addLebedevPoint(angularGrid,angularWeights,idx,  &
          [ mqc_float(0),mqc_float(0),mqc_float(1) ],wAxis)
        call addLebedevPoint(angularGrid,angularWeights,idx,  &
          [ mqc_float(0),mqc_float(0),-mqc_float(1) ],wAxis)
!
      case(14)
        Allocate(angularGrid(3,14),angularWeights(14))
        wAxis = mqc_float(4)*Pi/mqc_float(15)
        wCube = mqc_float(4)*Pi*mqc_float(3)/mqc_float(40)
        idx = 0
        call addLebedevAxes(angularGrid,angularWeights,idx,wAxis)
        a = mqc_float(1)/sqrt(mqc_float(3))
        call addLebedevCube(angularGrid,angularWeights,idx,a,wCube)
!
      case(26)
        Allocate(angularGrid(3,26),angularWeights(26))
        wAxis = mqc_float(4)*Pi/mqc_float(21)
        wEdge = mqc_float(4)*Pi*mqc_float(4)/mqc_float(105)
        wCube = mqc_float(4)*Pi*mqc_float(9)/mqc_float(280)
        idx = 0
        call addLebedevAxes(angularGrid,angularWeights,idx,wAxis)
        a = mqc_float(1)/sqrt(mqc_float(2))
        call addLebedevEdges(angularGrid,angularWeights,idx,a,wEdge)
        a = mqc_float(1)/sqrt(mqc_float(3))
        call addLebedevCube(angularGrid,angularWeights,idx,a,wCube)
!
      case default
        call mqc_error('setup_lebedev_angular_grid: unsupported Lebedev order.')
      end select
!
      return
      end subroutine setup_lebedev_angular_grid


!PROCEDURE setup_quadrature_lebedev_spherical
      subroutine setup_quadrature_lebedev_spherical(originBohr,rMaxBohr,  &
        nRadial,lebedevOrder,radialRule,quadGrid,quadWeights)
!
!     This routine builds a one-center spherical product quadrature grid using
!     selected Lebedev angular grids and either trapezoid or Gauss-Legendre
!     radial quadrature on 0 <= r <= rMaxBohr. Coordinates are in bohr and
!     weights are in bohr**3.
!
!
!     H. P. Hratchian, 2026.
!
!
      implicit none
      real(kind=real64),dimension(3),intent(in)::originBohr
      real(kind=real64),intent(in)::rMaxBohr
      integer(kind=int64),intent(in)::nRadial,lebedevOrder,radialRule
      real(kind=real64),dimension(:,:),allocatable,intent(out)::quadGrid
      real(kind=real64),dimension(:),allocatable,intent(out)::quadWeights
      integer(kind=int64)::iAngular,iRadial,idx,nAngular
      real(kind=real64)::stepSize
      real(kind=real64),dimension(:),allocatable::angularWeights,  &
        radialGrid,radialWeights
      real(kind=real64),dimension(:,:),allocatable::angularGrid
!
      if(rMaxBohr.le.mqc_float(0))  &
        call mqc_error('setup_quadrature_lebedev_spherical: rMaxBohr must be positive.')
      if(nRadial.lt.1)  &
        call mqc_error('setup_quadrature_lebedev_spherical: nRadial must be positive.')
      if(radialRule.eq.QUAD_RADIAL_TRAPEZOID.and.nRadial.lt.2)  &
        call mqc_error('setup_quadrature_lebedev_spherical: trapezoid radial grid needs nRadial >= 2.')
!
      Allocate(radialGrid(nRadial),radialWeights(nRadial))
      select case(radialRule)
      case(QUAD_RADIAL_TRAPEZOID)
        stepSize = rMaxBohr/mqc_float(nRadial-1)
        call setup_quadrature_trapezoid1d(nRadial,stepSize,mqc_float(0),  &
          radialGrid,radialWeights)
      case(QUAD_RADIAL_GAUSS_LEGENDRE)
        call setup_quadrature_gauss_legendre1d(nRadial,mqc_float(0),  &
          rMaxBohr,radialGrid,radialWeights)
      case default
        call mqc_error('setup_quadrature_lebedev_spherical: unsupported radial rule.')
      end select
!
      call setup_lebedev_angular_grid(lebedevOrder,angularGrid,  &
        angularWeights)
      nAngular = SIZE(angularWeights)
      Allocate(quadGrid(3,nRadial*nAngular),quadWeights(nRadial*nAngular))
!
      idx = 0
      do iRadial = 1,nRadial
        do iAngular = 1,nAngular
          idx = idx+1
          quadGrid(:,idx) = originBohr+radialGrid(iRadial)*  &
            angularGrid(:,iAngular)
          quadWeights(idx) = radialWeights(iRadial)*  &
            radialGrid(iRadial)*radialGrid(iRadial)*angularWeights(iAngular)
        endDo
      endDo
!
      return
      end subroutine setup_quadrature_lebedev_spherical


!PROCEDURE addLebedevPoint
      subroutine addLebedevPoint(angularGrid,angularWeights,idx,point,weight)
!
!     This small helper appends one angular point and weight to a Lebedev grid.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),dimension(3),intent(in)::point
      real(kind=real64),intent(in)::weight
!
      idx = idx+1
      angularGrid(:,idx) = point
      angularWeights(idx) = weight
!
      return
      end subroutine addLebedevPoint


!PROCEDURE addLebedevAxes
      subroutine addLebedevAxes(angularGrid,angularWeights,idx,weight)
!
!     This helper appends the six Cartesian axis points used in low-order
!     Lebedev grids.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
      integer(kind=int64),intent(inout)::idx
      real(kind=real64),intent(in)::weight
!
      call addLebedevPoint(angularGrid,angularWeights,idx,  &
        [ mqc_float(1),mqc_float(0),mqc_float(0) ],weight)
      call addLebedevPoint(angularGrid,angularWeights,idx,  &
        [ -mqc_float(1),mqc_float(0),mqc_float(0) ],weight)
      call addLebedevPoint(angularGrid,angularWeights,idx,  &
        [ mqc_float(0),mqc_float(1),mqc_float(0) ],weight)
      call addLebedevPoint(angularGrid,angularWeights,idx,  &
        [ mqc_float(0),-mqc_float(1),mqc_float(0) ],weight)
      call addLebedevPoint(angularGrid,angularWeights,idx,  &
        [ mqc_float(0),mqc_float(0),mqc_float(1) ],weight)
      call addLebedevPoint(angularGrid,angularWeights,idx,  &
        [ mqc_float(0),mqc_float(0),-mqc_float(1) ],weight)
!
      return
      end subroutine addLebedevAxes


!PROCEDURE addLebedevCube
      subroutine addLebedevCube(angularGrid,angularWeights,idx,a,weight)
!
!     This helper appends the eight cubic diagonal points (a,a,a) with all
!     sign combinations.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
      integer(kind=int64),intent(inout)::idx
      integer(kind=int64)::iSign,jSign,kSign
      real(kind=real64),intent(in)::a,weight
!
      do kSign = -1,1,2
        do jSign = -1,1,2
          do iSign = -1,1,2
            call addLebedevPoint(angularGrid,angularWeights,idx,  &
              [ mqc_float(iSign)*a,mqc_float(jSign)*a,  &
                mqc_float(kSign)*a ],weight)
          endDo
        endDo
      endDo
!
      return
      end subroutine addLebedevCube


!PROCEDURE addLebedevEdges
      subroutine addLebedevEdges(angularGrid,angularWeights,idx,a,weight)
!
!     This helper appends the twelve points with two nonzero equal-magnitude
!     coordinates used in the 26-point Lebedev grid.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),intent(inout)::angularGrid
      real(kind=real64),dimension(:),intent(inout)::angularWeights
      integer(kind=int64),intent(inout)::idx
      integer(kind=int64)::iSign,jSign
      real(kind=real64),intent(in)::a,weight
!
      do jSign = -1,1,2
        do iSign = -1,1,2
          call addLebedevPoint(angularGrid,angularWeights,idx,  &
            [ mqc_float(iSign)*a,mqc_float(jSign)*a,mqc_float(0) ],  &
            weight)
          call addLebedevPoint(angularGrid,angularWeights,idx,  &
            [ mqc_float(iSign)*a,mqc_float(0),mqc_float(jSign)*a ],  &
            weight)
          call addLebedevPoint(angularGrid,angularWeights,idx,  &
            [ mqc_float(0),mqc_float(iSign)*a,mqc_float(jSign)*a ],  &
            weight)
        endDo
      endDo
!
      return
      end subroutine addLebedevEdges


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
      integralValue = dot_product(quadratureWeights,valuesGrid)
      if(abs(integralValue).lt.mqc_small) integralValue = mqc_float(0)
      return
      end function moInnerProductNumericalIntegration

      end module gbs_mod
