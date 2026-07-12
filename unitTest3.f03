      program unitTest3
!
!     This program tests spherical Lebedev product quadrature utilities used
!     for future PAD matrix-element integration grids.
!
!     Command line arguments:
!           No command line arguments are required.
!
!     Usage:
!           ./unitTest3.exe
!
!
!     H. P. Hratchian, 2026.
!
      use iso_fortran_env
      use mqc_general
      use gbs_mod
      implicit none
      real(kind=real64),parameter::tolAngular=1.0e-10_real64
!
!     Format statements.
!
 1000 format(1x,'Starting unitTest3')
 2000 format(1x,'Lebedev order ',i4,3x,'angular weight sum = ',f20.10)
 2100 format(1x,a,3x,'nPoints = ',i8,3x,'volume = ',f20.10)
 2999 format(/,1x,'unitTest3 complete.')
!
!     Start the unit test program.
!
      write(iOut,1000)
      call mqc_version_print(iOut)
!
      call runLebedevAngularCase(6_int64)
      call runLebedevAngularCase(14_int64)
      call runLebedevAngularCase(26_int64)
      call runSphericalProductCase('trapezoid',QUAD_RADIAL_TRAPEZOID,  &
        801_int64,26_int64,1.0e-5_real64)
      call runSphericalProductCase('gauss-legendre',  &
        QUAD_RADIAL_GAUSS_LEGENDRE,6_int64,26_int64,1.0e-10_real64)
      call runRadialSizingCase()
!
      write(iOut,2999)
!
      contains

!PROCEDURE runLebedevAngularCase
      subroutine runLebedevAngularCase(lebedevOrder)
!
!     This routine checks angular Lebedev grid normalization and low-order
!     spherical moments.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::lebedevOrder
      integer(kind=int64)::i
      real(kind=real64)::weightSum
      real(kind=real64),dimension(3)::firstMoment,secondMoment
      real(kind=real64),dimension(:),allocatable::angularWeights
      real(kind=real64),dimension(:,:),allocatable::angularGrid
!
      call setup_lebedev_angular_grid(lebedevOrder,angularGrid,  &
        angularWeights)
      if(Size(angularWeights).ne.lebedevOrder) then
        write(iOut,*)' Size(angularWeights) = ',Size(angularWeights)
        write(iOut,*)' lebedevOrder         = ',lebedevOrder
        call mqc_error('unitTest3: unexpected Lebedev grid size.')
      endIf
      weightSum = SUM(angularWeights)
      firstMoment = mqc_float(0)
      secondMoment = mqc_float(0)
      do i = 1,Size(angularWeights)
        call assertNear(sqrt(dot_product(angularGrid(:,i),angularGrid(:,i))),  &
          mqc_float(1),tolAngular,'Lebedev point norm')
        firstMoment = firstMoment+angularWeights(i)*angularGrid(:,i)
        secondMoment = secondMoment+angularWeights(i)*  &
          angularGrid(:,i)*angularGrid(:,i)
      endDo
      write(iOut,'(1x,''Lebedev order '',i4,3x,''angular weight sum = '',f20.10)')  &
        lebedevOrder,weightSum
      call assertNear(weightSum,mqc_float(4)*Pi,tolAngular,  &
        'Lebedev angular weight sum')
      call assertVectorNear(firstMoment,[ mqc_float(0),mqc_float(0),  &
        mqc_float(0) ],tolAngular,'Lebedev angular first moment')
      call assertVectorNear(secondMoment,[ mqc_float(4)*Pi/mqc_float(3),  &
        mqc_float(4)*Pi/mqc_float(3),mqc_float(4)*Pi/mqc_float(3) ],  &
        tolAngular,'Lebedev angular second moment')
!
      return
      end subroutine runLebedevAngularCase


!PROCEDURE runSphericalProductCase
      subroutine runSphericalProductCase(label,radialRule,nRadial,  &
        lebedevOrder,tolerance)
!
!     This routine checks one-center spherical product quadrature volume and
!     central Cartesian moments.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::label
      integer(kind=int64),intent(in)::radialRule,nRadial,lebedevOrder
      real(kind=real64),intent(in)::tolerance
      integer(kind=int64)::i
      real(kind=real64)::rMaxBohr,volume,expectedVolume,  &
        expectedSecondMoment
      real(kind=real64),dimension(3)::centralPoint,firstMoment,originBohr,  &
        secondMoment
      real(kind=real64),dimension(:),allocatable::quadWeights
      real(kind=real64),dimension(:,:),allocatable::quadGrid
!
      rMaxBohr = mqc_float(8)
      originBohr = [ mqc_float(2)/mqc_float(10),  &
        -mqc_float(3)/mqc_float(10),mqc_float(4)/mqc_float(10) ]
      call setup_quadrature_lebedev_spherical(originBohr,rMaxBohr,nRadial,  &
        lebedevOrder,radialRule,quadGrid,quadWeights)
      if(Size(quadWeights).ne.nRadial*lebedevOrder) then
        write(iOut,*)' Size(quadWeights) = ',Size(quadWeights)
        write(iOut,*)' expected          = ',nRadial*lebedevOrder
        call mqc_error('unitTest3: unexpected spherical grid size.')
      endIf
      volume = SUM(quadWeights)
      expectedVolume = mqc_float(4)*Pi*rMaxBohr**3/mqc_float(3)
      expectedSecondMoment = mqc_float(4)*Pi*rMaxBohr**5/mqc_float(15)
      firstMoment = mqc_float(0)
      secondMoment = mqc_float(0)
      do i = 1,Size(quadWeights)
        centralPoint = quadGrid(:,i)-originBohr
        firstMoment = firstMoment+quadWeights(i)*centralPoint
        secondMoment = secondMoment+quadWeights(i)*centralPoint*centralPoint
      endDo
      write(iOut,'(1x,a,3x,''nPoints = '',i8,3x,''volume = '',f20.10)')  &
        label,Size(quadWeights),volume
      call assertNear(volume,expectedVolume,tolerance*expectedVolume,  &
        TRIM(label)//' spherical volume')
      call assertVectorNear(firstMoment,[ mqc_float(0),mqc_float(0),  &
        mqc_float(0) ],tolerance*expectedVolume,  &
        TRIM(label)//' central first moment')
      call assertVectorNear(secondMoment,  &
        [ expectedSecondMoment,expectedSecondMoment,  &
          expectedSecondMoment ],tolerance*expectedSecondMoment,  &
        TRIM(label)//' central second moment')
!
      return
      end subroutine runSphericalProductCase


!PROCEDURE runRadialSizingCase
      subroutine runRadialSizingCase()
!
!     This routine checks automatic radial cutoff and radial point-count helper
!     formulas for future Lebedev PAD quadrature selection.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64)::nRadial
      real(kind=real64)::alphaMin,expectedRMaxBohr,kMag,  &
        pointsPerWavelength,rMaxBohr,tailTol
      real(kind=real64),dimension(3)::originBohr
      real(kind=real64),dimension(3,2)::centerCoordsBohr
      real(kind=real64),dimension(3)::primitiveExponents
!
      originBohr = [ mqc_float(0),mqc_float(0),mqc_float(0) ]
      centerCoordsBohr(:,1) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      centerCoordsBohr(:,2) = [ mqc_float(0),mqc_float(3),mqc_float(4) ]
      primitiveExponents = [ mqc_float(2),mqc_float(1),  &
        mqc_float(1)/mqc_float(4) ]
      alphaMin = mqc_float(1)/mqc_float(4)
      tailTol = 1.0e-8_real64
      expectedRMaxBohr = mqc_float(5)+sqrt(-log(tailTol)/alphaMin)
      rMaxBohr = lebedev_spherical_rmax_from_basis_tail(originBohr,  &
        centerCoordsBohr,primitiveExponents,tailTol)
      call assertNear(rMaxBohr,expectedRMaxBohr,1.0e-10_real64,  &
        'Lebedev radial tail rMax')
!
      kMag = Pi
      pointsPerWavelength = mqc_float(8)
      nRadial = lebedev_spherical_nradial_from_k(mqc_float(3),kMag,  &
        pointsPerWavelength,5_int64)
      if(nRadial.ne.13_int64) then
        write(iOut,*)' nRadial = ',nRadial
        call mqc_error('unitTest3: unexpected k-derived nRadial.')
      endIf
      nRadial = lebedev_spherical_nradial_from_k(mqc_float(3),  &
        mqc_float(0),pointsPerWavelength,7_int64)
      if(nRadial.ne.7_int64) then
        write(iOut,*)' nRadial = ',nRadial
        call mqc_error('unitTest3: zero-k nRadial should use minimum.')
      endIf
!
      return
      end subroutine runRadialSizingCase


!PROCEDURE assertNear
      subroutine assertNear(value,target,tolerance,label)
!
!     This routine checks whether two real values are equal within a specified
!     tolerance.
!
!
!     H. P. Hratchian, 2026.
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
        call mqc_error('unitTest3: scalar comparison failed.')
      endIf
!
      return
      end subroutine assertNear


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
        call mqc_error('unitTest3: vector shape comparison failed.')
      endIf
      if(maxval(abs(value-target)).gt.tolerance) then
        write(iOut,*)' label     = ',TRIM(label)
        call mqc_print(value,iOut,header='value')
        call mqc_print(target,iOut,header='target')
        write(iOut,*)' tolerance = ',tolerance
        call mqc_error('unitTest3: vector comparison failed.')
      endIf
!
      return
      end subroutine assertVectorNear

      end program unitTest3
