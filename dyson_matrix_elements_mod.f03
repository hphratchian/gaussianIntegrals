      module dyson_matrix_elements_mod
!
!     This module contains reusable Dyson matrix-element kernels and spherical
!     helper routines used by PAD calculations.
!
      use iso_fortran_env
      use omp_lib
      use mqc_general
      use mqc_integrals1
      use mqc_gaussian
      use memory_utils
!
      implicit none
      private
      public::dysonPlaneWaveMatrixElementSquared
      public::dysonPlaneWaveMatrixElementSquaredThetaList
      public::dysonPlaneWaveMatrixElementSquaredPhiThetaList
      public::compute_legendre
      public::Ylm_complex
      public::sph_bessel_j
      public::dysonMatrixElement1Angle
      public::dysonMatrixElementThetaList
      public::generate_sph_grid
      integer(kind=int64),parameter::iOut=6_int64
      logical::MEMChecks=.false.
!
!
      CONTAINS
!
!PROCEDURE dysonPlaneWaveMatrixElementSquared
      function dysonPlaneWaveMatrixElementSquared(theta,kMag,  &
        photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,  &
        quadraturePoints,quadratureWeights) result(MSquared)
!
!     Computes the plane-wave Dyson transition dipole squared for a single
!     theta value. The photoelectron momentum magnitude, kMag, is in atomic
!     units.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      real(kind=real64),intent(in)::theta,kMag
      real(kind=real64),dimension(3),intent(inOut)::photonVector,  &
        orthogPlaneVector
      real(kind=real64),dimension(:),intent(in)::dysonCoeffs,  &
        quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      class(mqc_basisSet),intent(in)::aoBasisSet
      real(kind=real64)::MSquared
!
      integer(kind=int64)::i
      real(kind=real64)::dysonVal,epsilonDotMu,MReal,MImaginary,w
      real(kind=real64),dimension(3)::kVector
      real(kind=real64),dimension(:),allocatable::aoBasisValues,  &
        MValuesReal,MValuesImaginary
!
      call mqc_normalizeVector(photonVector)
      call mqc_normalizeVector(orthogPlaneVector)
      kVector = cos(theta)*photonVector+sin(theta)*orthogPlaneVector
      kVector = kMag*kVector
      Allocate(MValuesReal(SIZE(quadratureWeights)),  &
        MValuesImaginary(SIZE(quadratureWeights)))
!$omp parallel do private(i,aoBasisValues,w,epsilonDotMu,dysonVal)  &
!$omp& shared(MValuesReal,MValuesImaginary) schedule(dynamic)
      do i = 1,SIZE(quadratureWeights)
        w = dot_product(kVector,quadraturePoints(:,i))
        epsilonDotMu = dot_product(photonVector,quadraturePoints(:,i))
        call basisSetValuesList1(aoBasisSet,  &
          quadraturePoints(:,i),aoBasisValues)
        dysonVal = dot_product(dysonCoeffs,aoBasisValues)
        MValuesReal(i) = cos(w)*epsilonDotMu*dysonVal
        MValuesImaginary(i) = -sin(w)*epsilonDotMu*dysonVal
      endDo
!$omp end parallel do
      MReal = dot_product(quadratureWeights,MValuesReal)
      MImaginary = dot_product(quadratureWeights,MValuesImaginary)
      MSquared = MReal**2+MImaginary**2
!
      return
      end function dysonPlaneWaveMatrixElementSquared

!
!PROCEDURE dysonPlaneWaveMatrixElementSquaredThetaList
      function dysonPlaneWaveMatrixElementSquaredThetaList(thetaList,  &
        kMag,photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,  &
        quadraturePoints,quadratureWeights) result(MSquared)
!
!     Computes the plane-wave Dyson transition dipole squared for a list of
!     theta values. The photoelectron momentum magnitude, kMag, is in atomic
!     units.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      real(kind=real64),intent(in)::kMag
      real(kind=real64),dimension(3),intent(inOut)::photonVector,  &
        orthogPlaneVector
      real(kind=real64),dimension(:),intent(in)::thetaList,dysonCoeffs,  &
        quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      class(mqc_basisSet),intent(in)::aoBasisSet
      real(kind=real64),dimension(:),allocatable::MSquared
!
      integer(kind=int64)::i,j,nTheta,nGrid
      real(kind=real64)::w,dysonVal,epsilonDotMu,thetaTest
      real(kind=real64),dimension(:),allocatable::aoBasisValues
      real(kind=real64),dimension(:),allocatable::MReal,MImag,  &
        MRealThread,MImagThread
      real(kind=real64),dimension(3)::gridPoint,kVector
!
      nTheta = SIZE(thetaList)
      nGrid = SIZE(quadratureWeights)
      Allocate(MSquared(nTheta),MReal(nTheta),MImag(nTheta))
      MSquared = mqc_float(0)
      MReal = mqc_float(0)
      MImag = mqc_float(0)
      if(dot_product(photonVector,photonVector).gt.mqc_small)  &
        call mqc_normalizeVector(photonVector)
      if(dot_product(orthogPlaneVector,orthogPlaneVector).gt.mqc_small)  &
        call mqc_normalizeVector(orthogPlaneVector)
      if(MEMChecks) call print_memory_usage(iOut,  &
        'dysonPlaneWaveMatrixElementSquaredThetaList before OMP loop.')
!
!$omp parallel default(shared) private(i,j,w,dysonVal,epsilonDotMu,  &
!$omp& aoBasisValues,gridPoint,kVector,MRealThread,MImagThread,thetaTest)
      allocate(aoBasisValues(SIZE(dysonCoeffs)))
      allocate(MRealThread(nTheta),MImagThread(nTheta))
      MRealThread = mqc_float(0)
      MImagThread = mqc_float(0)
!$omp do schedule(static)
      do i = 1,nGrid
        gridPoint = quadraturePoints(:,i)
        call basisSetValuesList1(aoBasisSet,gridPoint,aoBasisValues)
        dysonVal = dot_product(dysonCoeffs,aoBasisValues)
        epsilonDotMu = dot_product(photonVector,gridPoint)
        do j = 1,nTheta
          kVector = cos(thetaList(j))*photonVector+  &
            sin(thetaList(j))*orthogPlaneVector
          thetaTest = vectorAngle(kVector,photonVector)
          if(abs(thetaList(j)-thetaTest).gt.0.001_real64) then
            write(iOut,'(4x,"PROBLEM  ",f6.3,5(",",f6.3))')  &
              kVector,photonVector
            call mqc_print(photonVector,iOut,header='photonVector')
            call mqc_print(orthogPlaneVector,iOut,header='orthogPlaneVector')
            call mqc_print(kVector,iOut,header='kVector')
            write(iOut,'(1x,"theta=",f8.4," | thetaTest=",f8.4)')  &
              thetaList(j),thetaTest
            write(iOut,'(A,3(3x,f10.3))')' angle = ',  &
              vectorAngle(kVector,photonVector),thetaList(j),  &
              dot_product(kVector,photonVector)
            call mqc_error('STOP: ERROR!')
          endIf
          w = kMag*dot_product(kVector,gridPoint)
          MRealThread(j) = MRealThread(j)+cos(w)*epsilonDotMu*  &
            dysonVal*quadratureWeights(i)
          MImagThread(j) = MImagThread(j)-sin(w)*epsilonDotMu*  &
            dysonVal*quadratureWeights(i)
        endDo
      endDo
!$omp end do
!$omp critical
      MReal = MReal+MRealThread
      MImag = MImag+MImagThread
!$omp end critical
      deallocate(aoBasisValues)
      deallocate(MRealThread,MImagThread)
!$omp end parallel
      MSquared = MReal**2+MImag**2
      if(MEMChecks) call print_memory_usage(iOut,  &
        'dysonPlaneWaveMatrixElementSquaredThetaList after OMP loop.')
!
      return
      end function dysonPlaneWaveMatrixElementSquaredThetaList

!
!PROCEDURE dysonPlaneWaveMatrixElementSquaredPhiThetaList
      function dysonPlaneWaveMatrixElementSquaredPhiThetaList(phiList,  &
        thetaList,kMag,photonVector,orthogPlaneVector,dysonCoeffs,  &
        aoBasisSet,quadraturePoints,quadratureWeights) result(MSquared)
!
!     Computes the plane-wave Dyson transition dipole squared for lists of phi
!     and theta values. The photoelectron momentum magnitude, kMag, is in
!     atomic units.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      real(kind=real64),intent(in)::kMag
      real(kind=real64),dimension(3),intent(inOut)::photonVector,  &
        orthogPlaneVector
      real(kind=real64),dimension(:),intent(in)::phiList,thetaList,  &
        dysonCoeffs,quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      class(mqc_basisSet),intent(in)::aoBasisSet
      real(kind=real64),dimension(:,:),allocatable::MSquared
!
      integer(kind=int64)::i,j,k,nPhi,nTheta,nGrid
      real(kind=real64)::w,dysonVal,epsilonDotMu,thetaTest
      real(kind=real64),dimension(:),allocatable::aoBasisValues
      real(kind=real64),dimension(:,:),allocatable::MReal,MImag
      real(kind=real64),dimension(3)::orthog2PlaneVector,planeVector,  &
        gridPoint,kVector
!
      nPhi = SIZE(phiList)
      nTheta = SIZE(thetaList)
      nGrid = SIZE(quadratureWeights)
      Allocate(MSquared(nPhi,nTheta),MReal(nPhi,nTheta),MImag(nPhi,nTheta))
      MSquared = mqc_float(0)
      MReal = mqc_float(0)
      MImag = mqc_float(0)
      if(dot_product(photonVector,photonVector).gt.mqc_small)  &
        call mqc_normalizeVector(photonVector)
      if(dot_product(orthogPlaneVector,orthogPlaneVector).gt.mqc_small)  &
        call mqc_normalizeVector(orthogPlaneVector)
      orthog2PlaneVector = mqc_crossProduct3D_real(photonVector,  &
        orthogPlaneVector)
      if(MEMChecks) call print_memory_usage(iOut,  &
        'dysonPlaneWaveMatrixElementSquaredPhiThetaList before OMP loop.')
!
      allocate(aoBasisValues(SIZE(dysonCoeffs)))
      do i = 1,nGrid
        gridPoint = quadraturePoints(:,i)
        call basisSetValuesList1(aoBasisSet,gridPoint,aoBasisValues)
        dysonVal = dot_product(dysonCoeffs,aoBasisValues)
        epsilonDotMu = dot_product(photonVector,gridPoint)
        do j = 1,nTheta
          do k = 1,nPhi
            planeVector = cos(phiList(k))*orthogPlaneVector+  &
              sin(phiList(k))*orthog2PlaneVector
            kVector = cos(thetaList(j))*photonVector+  &
              sin(thetaList(j))*planeVector
            thetaTest = vectorAngle(kVector,photonVector)
            if(abs(thetaList(j)-thetaTest).gt.0.001_real64)  &
              write(iOut,'(4x,"PROBLEM  ",f6.3,5(",",f6.3))')  &
              kVector,photonVector
            w = kMag*dot_product(kVector,gridPoint)
            MReal(k,j) = MReal(k,j)+cos(w)*epsilonDotMu*dysonVal*  &
              quadratureWeights(i)
            MImag(k,j) = MImag(k,j)-sin(w)*epsilonDotMu*dysonVal*  &
              quadratureWeights(i)
          endDo
        endDo
      endDo
      deallocate(aoBasisValues)
      MSquared = MReal**2+MImag**2
      call mqc_matrixTrimZero(MSquared)
      if(MEMChecks) call print_memory_usage(iOut,  &
        'dysonPlaneWaveMatrixElementSquaredPhiThetaList after OMP loop.')
!
      return
      end function dysonPlaneWaveMatrixElementSquaredPhiThetaList

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
      if(l.lt.0.or.abs(m).gt.l)  &
        call mqc_error('compute_legendre: Invalid l or m.')
!
      abs_m = abs(m)
      somx2 = sqrt(max(mqc_float(0),mqc_float(1)-x*x))
      pmm = mqc_float(1)
      if(abs_m.gt.0) then
        pmm = (-mqc_float(1))**abs_m
        do i = 1,abs_m
          pmm = pmm*somx2*mqc_float(2*i-1)
        endDo
      endIf
      if(l.eq.abs_m) then
        Plm = pmm
        return
      endIf
!
      pmmp1 = x*pmm*mqc_float(2*abs_m+1)
      if(l.eq.abs_m+1) then
        Plm = pmmp1
        return
      endIf
!
      do i = abs_m+2,l
        pll = ((2*i-1)*x*pmmp1-(i+abs_m-1)*pmm)/(i-abs_m)
        pmm = pmmp1
        pmmp1 = pll
      endDo
      Plm = pll
!
      if(m.lt.0) Plm = (-mqc_float(1))**abs_m*Plm
!
      return
      end subroutine compute_legendre

!
!PROCEDURE Ylm_complex
      function Ylm_complex(l,m,theta,phi) result(Y)
!
!     Evaluates the complex spherical harmonic Y_l^m(theta,phi) using
!     associated Legendre polynomials and exp(i*m*phi).
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      integer(kind=int64),intent(in)::l,m
      real(kind=real64),intent(in)::theta,phi
      complex(kind=real64)::Y
!
      integer(kind=int64)::abs_m
      real(kind=real64)::P_lm,x,norm
      complex(kind=real64)::eimphi
!
      x = cos(theta)
      abs_m = abs(m)
      call compute_legendre(l,m,x,P_lm)
      norm = sqrt((2*l+1)/(mqc_float(4)*Pi)*  &
        mqc_float(factorial(l-abs_m))/mqc_float(factorial(l+abs_m)))
      eimphi = cmplx(mqc_float(0),mqc_float(m)*phi,kind=real64)
      Y = norm*P_lm*exp(eimphi)
!
      return
      end function Ylm_complex

!
!PROCEDURE sph_bessel_j
      function sph_bessel_j(l,x) result(jl)
!
!     Computes the spherical Bessel function j_l(x) using direct expressions
!     for l=0,1 and recurrence for l>1.
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
      if(abs(x).lt.mqc_small) then
        if(l.eq.0) then
          jl = mqc_float(1)
        else
          jl = mqc_float(0)
        endIf
        return
      endIf
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
      j0 = sin(x)/x
      j1 = (sin(x)/x**2)-(cos(x)/x)
      do n = 1,l-1
        jn = ((2*n+1)/x)*j1-j0
        j0 = j1
        j1 = jn
      endDo
      jl = jn
!
      return
      end function sph_bessel_j

!
!PROCEDURE dysonMatrixElement1Angle
      subroutine dysonMatrixElement1Angle(iPEType,lMax,theta,kMag,  &
        photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,  &
        quadraturePoints,quadratureWeights,MSquared,lWeights)
!
!     Computes the Dyson transition matrix element squared for one theta value
!     using the requested outgoing-electron representation. The photoelectron
!     momentum magnitude, kMag, is in atomic units.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::iPEType,lMax
      real(kind=real64),intent(in)::theta,kMag
      real(kind=real64),dimension(3),intent(inOut)::photonVector,  &
        orthogPlaneVector
      real(kind=real64),dimension(:),intent(in)::dysonCoeffs,  &
        quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      real(kind=real64),intent(out)::MSquared
      real(kind=real64),dimension(0:),intent(out)::lWeights
      class(mqc_basisSet),intent(in)::aoBasisSet
!
      integer(kind=int64)::i,l,m
      real(kind=real64)::r,thetaVal,phiVal,w,j_l,dysonVal,  &
        epsilonDotR
      real(kind=real64),dimension(3)::rVec,kVector
      real(kind=real64),dimension(0:lMax)::W_l
      real(kind=real64),dimension(:),allocatable::aoBasisValues,  &
        MValuesReal,MValuesImaginary
      complex(kind=real64)::Ylm,muVal,dipoleAmp
      complex(kind=real64),dimension(0:lMax,-lMax:lMax)::cLM
!
      call mqc_normalizeVector(photonVector)
      call mqc_normalizeVector(orthogPlaneVector)
      kVector = cos(theta)*photonVector+sin(theta)*orthogPlaneVector
      kVector = kMag*kVector
      Allocate(MValuesReal(size(quadratureWeights)),  &
        MValuesImaginary(size(quadratureWeights)))
      MValuesReal = mqc_float(0)
      MValuesImaginary = mqc_float(0)
      cLM = cmplx(mqc_float(0),mqc_float(0),kind=real64)
!
!$omp parallel do private(i,aoBasisValues,rVec,r,dysonVal,thetaVal,phiVal,  &
!$omp& l,m,j_l,Ylm,muVal,epsilonDotR,w)  &
!$omp& shared(MValuesReal,MValuesImaginary,cLM) schedule(dynamic)
      do i = 1,size(quadratureWeights)
        rVec = quadraturePoints(:,i)
        r = sqrt(dot_product(rVec,rVec))
        call basisSetValuesList1(aoBasisSet,rVec,aoBasisValues)
        dysonVal = dot_product(dysonCoeffs,aoBasisValues)
!
        select case(iPEType)
        case(1)
          w = dot_product(kVector,rVec)
          epsilonDotR = dot_product(photonVector,rVec)
          MValuesReal(i) = cos(w)*epsilonDotR*dysonVal
          MValuesImaginary(i) = -sin(w)*epsilonDotR*dysonVal
!
        case(2)
          if(r.gt.mqc_small) then
            thetaVal = acos(rVec(3)/r)
            phiVal = atan2(rVec(2),rVec(1))
          else
            thetaVal = mqc_float(0)
            phiVal = mqc_float(0)
          endIf
          epsilonDotR = dot_product(photonVector,rVec)
          do l = 0,lMax
            j_l = sph_bessel_j(l,kMag*r)
            do m = -l,l
              Ylm = Ylm_complex(l,m,thetaVal,phiVal)
              muVal = j_l*Ylm*epsilonDotR*dysonVal
!$omp atomic
              cLM(l,m) = cLM(l,m)+quadratureWeights(i)*muVal
            endDo
          endDo
!
        case default
          call mqc_error('ERROR: Invalid iPEType in dysonMatrixElement1Angle.')
        end select
      endDo
!$omp end parallel do
!
      if(iPEType.eq.2) then
        do l = 0,lMax
          W_l(l) = mqc_float(0)
          do m = -l,l
            W_l(l) = W_l(l)+abs(cLM(l,m))**2
          endDo
        endDo
        if(sum(W_l).gt.mqc_small) then
          lWeights(0:lMax) = W_l/sum(W_l)
        else
          lWeights(0:lMax) = mqc_float(0)
        endIf
        dipoleAmp = cmplx(mqc_float(0),mqc_float(0),kind=real64)
        do l = 0,lMax
          do m = -l,l
            dipoleAmp = dipoleAmp+cLM(l,m)
          endDo
        endDo
        MSquared = abs(dipoleAmp)**2
      else
        MSquared = dot_product(quadratureWeights,MValuesReal)**2+  &
          dot_product(quadratureWeights,MValuesImaginary)**2
        lWeights(0) = mqc_float(1)
        lWeights(1:lMax) = mqc_float(0)
      endIf
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
!     outgoing-electron representation. Optional partial-wave angular momentum
!     weights are returned either intensity-averaged or per-theta.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::iPEType,lMax
      real(kind=real64),dimension(:),intent(in)::thetaVals
      real(kind=real64),intent(in)::kMag
      real(kind=real64),dimension(3),intent(in)::photonVector,  &
        orthogPlaneVector
      real(kind=real64),dimension(:),intent(in)::dysonCoeffs,  &
        quadratureWeights
      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
      real(kind=real64),dimension(:),intent(out)::Itheta
      real(kind=real64),dimension(0:),optional,intent(out)::lWeights
      real(kind=real64),dimension(0:,:),optional,intent(out)::lWeightsTheta
      class(mqc_basisSet),intent(in)::aoBasisSet
!
      integer(kind=int64)::iTh,l,nTh
      real(kind=real64)::theta,MSq,Itot
      real(kind=real64),dimension(0:lMax)::lWTemp,Wsum
      real(kind=real64),dimension(3)::pVec,oVec
!
      nTh = size(thetaVals)
      pVec = photonVector
      oVec = orthogPlaneVector
      Wsum = mqc_float(0)
      Itot = mqc_float(0)
      do iTh = 1,nTh
        theta = thetaVals(iTh)
        call dysonMatrixElement1Angle(iPEType,lMax,theta,kMag,pVec,  &
          oVec,dysonCoeffs,aoBasisSet,quadraturePoints,  &
          quadratureWeights,MSq,lWTemp)
        Itheta(iTh) = MSq
        Itot = Itot+MSq
        if(present(lWeightsTheta)) then
          do l = 0,lMax
            lWeightsTheta(l,iTh) = lWTemp(l)
          endDo
        endIf
        if(present(lWeights)) then
          do l = 0,lMax
            Wsum(l) = Wsum(l)+MSq*lWTemp(l)
          endDo
        endIf
      endDo
!
      if(present(lWeights)) then
        if(Itot.gt.mqc_small) then
          do l = 0,lMax
            lWeights(l) = Wsum(l)/Itot
          endDo
        else
          do l = 0,lMax
            lWeights(l) = mqc_float(0)
          endDo
        endIf
      endIf
!
      return
      end subroutine dysonMatrixElementThetaList

!
!PROCEDURE generate_sph_grid
      subroutine generate_sph_grid(nTheta,nPhi,thetaVals,phiVals,  &
        weights)
!
!     Generates a simple theta/phi spherical quadrature grid with
!     sin(theta) dtheta dphi weights.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::nTheta,nPhi
      real(kind=real64),dimension(nTheta),intent(out)::thetaVals
      real(kind=real64),dimension(nPhi),intent(out)::phiVals
      real(kind=real64),dimension(nTheta,nPhi),intent(out)::weights
!
      integer(kind=int64)::i,j
      real(kind=real64)::dtheta,dphi
!
      if(nTheta.lt.2)  &
        call mqc_error('generate_sph_grid: nTheta must be at least 2.')
      if(nPhi.lt.1)  &
        call mqc_error('generate_sph_grid: nPhi must be positive.')
      dtheta = Pi/mqc_float(nTheta-1)
      dphi = mqc_float(2)*Pi/mqc_float(nPhi)
      do i = 1,nTheta
        thetaVals(i) = dtheta*mqc_float(i-1)
      endDo
      do j = 1,nPhi
        phiVals(j) = dphi*mqc_float(j-1)
      endDo
      do i = 1,nTheta
        do j = 1,nPhi
          weights(i,j) = sin(thetaVals(i))*dtheta*dphi
        endDo
      endDo
!
      return
      end subroutine generate_sph_grid

!
!PROCEDURE vectorAngle
      function vectorAngle(v1,v2,degrees) result(angle)
!
!     Computes the angle between two vectors. The optional degrees argument
!     requests output in degrees rather than radians.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::v1,v2
      logical,optional::degrees
      real(kind=real64)::angle,denominator
!
      denominator = vectorMagnitude(v1)*vectorMagnitude(v2)
      if(denominator.gt.mqc_small) then
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
!     Computes the magnitude of vector v.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::v
      real(kind=real64)::vMagnitude
!
      vMagnitude = sqrt(dot_product(v,v))
!
      return
      end function vectorMagnitude

      end module dyson_matrix_elements_mod
