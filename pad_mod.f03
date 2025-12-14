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
!
!
      CONTAINS
!
!
!PROCEDURE padCommandLine
      subroutine padCommandLine(iPEType,iMODyson,nGridPointsTheta,  &
        nGridPointsM,kMag,fafName)
!
!     This routine processes the command line arguments and fills key variables
!     and option flags.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      integer(kind=int64)::iPEType,iMODyson,nGridPointsTheta,  &
        nGridPointsM
      real(kind=real64)::kMag
      character(len=256)::fafName
!
      integer(kind=int64)::i,nCommandLineArgs
!
!     Walk through the command line arguments to fill required input parameters
!     and other option flags.
!
      nCommandLineArgs = command_argument_count()
      if(nCommandLineArgs.lt.1.or.nCommandLineArgs.gt.6)  &
        call mqc_error('PAD expects 1-6 command line arguments.')
      if(.true.) then
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
          call mqc_get_command_argument_integer(4,nGridPointsTheta)
        else
          nGridPointsTheta = 15
          nGridPointsTheta = 5
        endIf
        if(command_argument_count().ge.5) then
          call mqc_get_command_argument_integer(5,nGridPointsM)
        else
          nGridPointsM = 101
        endIf
        if(command_argument_count().ge.6) then
          call mqc_get_command_argument_integer(6,iPEType)
        else
          iPEType = 0
        endIf
      else
        do i = 1,nCommandLineArgs

        endDo
      endIf
!
      return
      end subroutine padCommandLine


!hph+
!
!PROCEDURE pad
!
!     Allocate arrays used for the number of integration planes. Then fill the
!     arrays laserVector and orthogPlaneVector. Currently, there are two methods
!     in the program to do this. The first method simply sets up 3 integration
!     plane along the primary axes. The second methods sets up a set of equally
!     spaced vectors going around a unit sphere.
!
      if(.false.) then
        nIntPlanes = 3
        Allocate(integratedIntensity(nIntPlanes),  &
          betaValsParaPerp(nIntPlanes),betaValsFit(nIntPlanes),  &
          rSquared(nIntPlanes))
        Allocate(laserVector(3,nIntPlanes),  &
          orthogPlaneVector(3,nIntPlanes))
        Allocate(intPlaneLabels(nIntPlanes))
!
!       Set the laser electric field vector and the orthogonal vector defining
!       the integration plane to be used with each electric field. At present we
!       hardwire 3 experiments with the electric field vector and planar slice
!       of photoelectron transition dipole matrices:
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
      else
        call buildSphereGrid(laserVector,orthogPlaneVector)
        nIntPlanes = Size(laserVector,2)
        write(iOut,*)' nIntPlanes = ',nIntPlanes
        Allocate(integratedIntensity(nIntPlanes),  &
          betaValsParaPerp(nIntPlanes),betaValsFit(nIntPlanes),  &
          rSquared(nIntPlanes))
        Allocate(intPlaneLabels(nIntPlanes))
        intPlaneLabels = 'a'
      endIf
!hph-





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
!
!PROCEDURE buildSphereGrid
      subroutine buildSphereGrid(xyz,orthogVector)
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
      integer(kind=int64)::i,j,nThetaPoints,nPhiPoints,ixyz
      real(kind=real64)::thetaStep,phiStep,theta,phi,x,y,z
!
 1000 format(1x,i4,': theta,phi (pi)=',f5.2,',',f5.2,' | x,y,z=',f8.4,',',f8.4,',',f8.4)
!
!     Start by setting nThetaPoints and  nPhiPoints, and figuring out thetaStep
!     and phiStep. Then, allocate xyz.
!
      nThetaPoints = 5
      nPhiPoints = 8
      nThetaPoints = 2
      nPhiPoints = 3
      thetaStep = Pi/mqc_float(nThetaPoints-1)
      phiStep = mqc_float(2)*Pi/mqc_float(nPhiPoints)
      write(*,*)' thetaStep = ',thetaStep
      write(*,*)'           = ',thetaStep/Pi
      write(*,*)
      write(*,*)' phiStep   = ',phiStep
      write(*,*)'           = ',phiStep/Pi
      write(*,*)
!
!     Loop over theta and phi to build the array xyz.
!
      Allocate(xyz(3,(nThetaPoints-2)*nPhiPoints+2),  &
        orthogVector(3,(nThetaPoints-2)*nPhiPoints+2))
      theta = 0
      ixyz = 1
      do i = 1,nThetaPoints
        theta = mqc_float(i-1)*thetaStep
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
          write(iOut,1000) ixyz,theta/pi,phi/pi,xyz(:,ixyz)
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
          write(iOut,1000) ixyz,theta/pi,phi/pi,xyz(:,ixyz)
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
            write(iOut,1000) ixyz,theta/pi,phi/pi,xyz(:,ixyz)
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

!hph+
!!
!!PROCEDURE dysonMatrixElement1Angle
!      subroutine dysonMatrixElement1Angle(iPEType,lMax,theta,kMag,  &
!        photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,  &
!        quadraturePoints,quadratureWeights,MSquared,lWeights)
!!
!!     Computes the Dyson transition dipole squared at angle theta using
!!     the specified photoelectron representation. If iPEType==2, this
!!     also computes the normalized partial wave ℓ contributions.
!!
!!
!!     H. P. Hratchian, 2025.
!!
!      implicit none
!      integer(kind=int64),intent(in)::iPEType,lMax
!      real(kind=real64),intent(in)::theta,kMag
!      real(kind=real64),dimension(3),intent(inOut)::photonVector,       &
!        orthogPlaneVector
!      real(kind=real64),dimension(:),intent(in)::dysonCoeffs,           &
!        quadratureWeights
!      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
!      real(kind=real64),intent(out)::MSquared
!      real(kind=real64),dimension(0:),intent(out)::lWeights
!      class(mqc_basisSet),intent(in)::aoBasisSet
!!
!      integer(kind=int64)::i,l,m,mu
!      real(kind=real64)::rVec(3),r,thetaVal,phiVal,w
!      real(kind=real64)::j_l,dysonVal,dysonNorm
!      real(kind=real64)::epsilonDotMu
!      complex(kind=real64)::Ylm,psiF,muVal
!      complex(kind=real64)::dipX,dipY,dipZ
!      complex(kind=real64),dimension(0:lMax,-lMax:lMax,3)::cLM
!      real(kind=real64),dimension(0:lMax)::W_l
!      real(kind=real64),dimension(:),allocatable::aoBasisValues,        &
!        MValuesReal,MValuesImaginary,dysonNormTest
!      real(kind=real64),dimension(3)::kVector
!      complex(kind=real64)::iunit
!!
!!     Start by setting up a set of initial variables/arrays.
!!
!      iunit = (0.0,1.0)
!      call mqc_normalizeVector(photonVector)
!      call mqc_normalizeVector(orthogPlaneVector)
!      kVector = cos(theta)*photonVector + sin(theta)*orthogPlaneVector
!      kVector = kMag*kVector
!!
!!     Allocate some memory.
!!
!      Allocate(MValuesReal(SIZE(quadratureWeights)),  &
!        MValuesImaginary(SIZE(quadratureWeights)),    &
!        dysonNormTest(SIZE(quadratureWeights)))
!!
!!     Do the work...
!!
!      cLM = (0.0,0.0)
!!
!!$omp parallel do private(i, aoBasisValues, rVec, r, dysonVal, thetaVal, phiVal, &
!!$omp& l, m, mu, j_l, Ylm, muVal, w, epsilonDotMu, psiF)                          &
!!$omp& shared(MValuesReal, MValuesImaginary, dysonNormTest, cLM) schedule(dynamic)
!      do i = 1,SIZE(quadratureWeights)
!        rVec = quadraturePoints(:,i)
!        r = sqrt(dot_product(rVec,rVec))
!        call basisSetValuesList1(aoBasisSet, rVec, aoBasisValues)
!        dysonVal = dot_product(dysonCoeffs, aoBasisValues)
!        dysonNormTest(i) = dysonVal*dysonVal
!!
!        select case(iPEType)
!        case(1)
!          w = dot_product(kVector, rVec)
!          epsilonDotMu = dot_product(photonVector, rVec)
!          MValuesReal(i) = cos(w)*epsilonDotMu*dysonVal
!          MValuesImaginary(i) = -sin(w)*epsilonDotMu*dysonVal
!!
!        case(2)
!          if(r.gt.1.0d-10) then
!            thetaVal = acos(rVec(3)/r)
!            phiVal = atan2(rVec(2),rVec(1))
!          else
!            thetaVal = mqc_float(0)
!            phiVal = mqc_float(0)
!          endIf
!!
!!         Loop over l and m to build c_{lm}^{(mu)} terms
!          do l=0,lMax
!            j_l = sph_bessel_j(l,kMag*r)
!            do m=-l,l
!              Ylm = Ylm_complex(l,m,thetaVal,phiVal)
!              do mu=1,3
!                epsilonDotMu = rVec(mu)
!                muVal = j_l*conjg(Ylm)*epsilonDotMu*dysonVal
!!$omp atomic
!                cLM(l,m,mu) = cLM(l,m,mu) + quadratureWeights(i)*muVal
!              endDo
!            endDo
!          endDo
!!
!        case default
!          call mqc_error('ERROR: Invalid iPEType in dysonMatrixElement1Angle.')
!        end select
!      end do
!!$omp end parallel do
!!
!      dysonNorm = dot_product(quadratureWeights,dysonNormTest)
!!
!!     If we're working with partial waves, compute the normalized weights by
!!     angular momentum type.
!!
!      if(iPEType.eq.2) then
!        do l=0,lMax
!          W_l(l) = mqc_float(0)
!          do m=-l,l
!            do mu=1,3
!              W_l(l) = W_l(l) + abs(cLM(l,m,mu))**2
!            endDo
!          endDo
!        endDo
!        if(sum(W_l).gt.mqc_small) then
!          lWeights(0:lMax) = W_l / sum(W_l)
!        else
!          lWeights(0:lMax) = mqc_float(0)
!        endIf
!        MValuesReal = mqc_float(0)
!        MValuesImaginary = mqc_float(0)
!!
!!       Reconstruct dipole projection from cLM
!!
!        dipX = (0.0,0.0)
!        dipY = (0.0,0.0)
!        dipZ = (0.0,0.0)
!        do l=0,lMax
!          do m=-l,l
!            dipX = dipX + cLM(l,m,1)
!            dipY = dipY + cLM(l,m,2)
!            dipZ = dipZ + cLM(l,m,3)
!          endDo
!        endDo
!!
!        MSquared = abs( dipX*photonVector(1) +  &
!                        dipY*photonVector(2) +  &
!                        dipZ*photonVector(3) )**2
!!
!      else
!        MSquared = dot_product(quadratureWeights,MValuesReal)**2 +  &
!                   dot_product(quadratureWeights,MValuesImaginary)**2
!        lWeights = mqc_float(-1)
!      endIf
!!
!      return
!      end subroutine dysonMatrixElement1Angle
!hph-

!hph+
!!PROCEDURE dysonMatrixElement1Angle
!      subroutine dysonMatrixElement1Angle(iPEType,lMax,theta,kMag,  &
!        photonVector,orthogPlaneVector,dysonCoeffs,aoBasisSet,  &
!        quadraturePoints,quadratureWeights,MSquared,lWeights)
!!
!      implicit none
!      integer(kind=int64),intent(in)::iPEType,lMax
!      real(kind=real64),intent(in)::theta,kMag
!      real(kind=real64),dimension(3),intent(inOut)::photonVector,       &
!        orthogPlaneVector
!      real(kind=real64),dimension(:),intent(in)::dysonCoeffs,           &
!        quadratureWeights
!      real(kind=real64),dimension(:,:),intent(in)::quadraturePoints
!      real(kind=real64),intent(out)::MSquared
!      real(kind=real64),dimension(0:),intent(out)::lWeights
!      class(mqc_basisSet),intent(in)::aoBasisSet
!!
!      integer(kind=int64)::i,l,m,mu
!      real(kind=real64)::rVec(3),r,thetaVal,phiVal,w
!      real(kind=real64)::j_l,dysonVal,dysonNorm
!      real(kind=real64)::epsilonDotMu
!      complex(kind=real64)::Ylm,psiF,muVal
!      complex(kind=real64)::dipX,dipY,dipZ
!      complex(kind=real64),dimension(0:lMax,-lMax:lMax,3)::cLM
!      real(kind=real64),dimension(0:lMax)::W_l
!      real(kind=real64),dimension(:),allocatable::aoBasisValues,        &
!        MValuesReal,MValuesImaginary,dysonNormTest
!      real(kind=real64),dimension(3)::kVector
!      complex(kind=real64)::iunit
!!
!      iunit = (0.0,1.0)
!      call mqc_normalizeVector(photonVector)
!      call mqc_normalizeVector(orthogPlaneVector)
!      kVector = cos(theta)*photonVector + sin(theta)*orthogPlaneVector
!      kVector = kMag*kVector
!!
!      Allocate(MValuesReal(SIZE(quadratureWeights)),  &
!        MValuesImaginary(SIZE(quadratureWeights)),    &
!        dysonNormTest(SIZE(quadratureWeights)))
!!
!      cLM = (0.0,0.0)
!!
!!$omp parallel do private(i, aoBasisValues, rVec, r, dysonVal, thetaVal, phiVal, &
!!$omp& l, m, mu, j_l, Ylm, muVal, w, epsilonDotMu, psiF)                          &
!!$omp& shared(MValuesReal, MValuesImaginary, dysonNormTest, cLM) schedule(dynamic)
!      do i = 1,SIZE(quadratureWeights)
!        rVec = quadraturePoints(:,i)
!        r = sqrt(dot_product(rVec,rVec))
!        call basisSetValuesList1(aoBasisSet, rVec, aoBasisValues)
!        dysonVal = dot_product(dysonCoeffs, aoBasisValues)
!        dysonNormTest(i) = dysonVal*dysonVal
!!
!        select case(iPEType)
!        case(1)
!          w = dot_product(kVector, rVec)
!          epsilonDotMu = dot_product(photonVector, rVec)
!          MValuesReal(i) = cos(w)*epsilonDotMu*dysonVal
!          MValuesImaginary(i) = -sin(w)*epsilonDotMu*dysonVal
!!
!        case(2)
!          if(r.gt.1.0d-10) then
!            thetaVal = acos(rVec(3)/r)
!            phiVal = atan2(rVec(2),rVec(1))
!          else
!            thetaVal = mqc_float(0)
!            phiVal = mqc_float(0)
!          endIf
!!
!          do l=0,lMax
!            j_l = sph_bessel_j(l,kMag*r)
!            do m=-l,l
!              Ylm = Ylm_complex(l,m,thetaVal,phiVal)
!              do mu=1,3
!                epsilonDotMu = rVec(mu)
!                muVal = j_l * Ylm * epsilonDotMu * dysonVal
!!$omp atomic
!                cLM(l,m,mu) = cLM(l,m,mu) + quadratureWeights(i)*muVal
!              endDo
!            endDo
!          endDo
!!
!        case default
!          call mqc_error('ERROR: Invalid iPEType in dysonMatrixElement1Angle.')
!        end select
!      end do
!!$omp end parallel do
!!
!      dysonNorm = dot_product(quadratureWeights,dysonNormTest)
!!
!      if(iPEType.eq.2) then
!        do l=0,lMax
!          W_l(l) = mqc_float(0)
!          do m=-l,l
!            do mu=1,3
!              W_l(l) = W_l(l) + abs(cLM(l,m,mu))**2
!            endDo
!          endDo
!          write(iOut,*)'l,W_l(l) = ',l,W_l(l)
!        endDo
!        if(sum(W_l).gt.MQC_Small) then
!          lWeights(0:lMax) = W_l / sum(W_l)
!        else
!          lWeights(0:lMax) = mqc_float(0)
!        endIf
!!
!        dipX = (0.0,0.0)
!        dipY = (0.0,0.0)
!        dipZ = (0.0,0.0)
!        do l=0,lMax
!          do m=-l,l
!            dipX = dipX + cLM(l,m,1)
!            dipY = dipY + cLM(l,m,2)
!            dipZ = dipZ + cLM(l,m,3)
!          endDo
!        endDo
!!
!!       Debug print: inspect dipole components and projected MSquared
!        write(*,*) 'DEBUG: dipX = ', dipX
!        write(*,*) 'DEBUG: dipY = ', dipY
!        write(*,*) 'DEBUG: dipZ = ', dipZ
!        write(*,*) 'DEBUG: photonVector = ', photonVector
!        write(*,*) 'DEBUG: dot(dip, epsilon) = ',  &
!                    dipX*photonVector(1) +  &
!                    dipY*photonVector(2) +  &
!                    dipZ*photonVector(3)
!
!!
!        MSquared = abs( dipX*photonVector(1) +  &
!                        dipY*photonVector(2) +  &
!                        dipZ*photonVector(3) )**2
!!
!      else
!        MSquared = dot_product(quadratureWeights,MValuesReal)**2 +  &
!                   dot_product(quadratureWeights,MValuesImaginary)**2
!        lWeights(0) = mqc_float(1)
!        do l=1,lMax
!          lWeights(l) = mqc_float(0)
!        endDo
!      endIf
!!
!      return
!      end subroutine dysonMatrixElement1Angle
!hph-

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
      write(*,*)' l,W_l(l) = ',l,W_l(l)
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
      write(*,*)
      write(*,*)' Hrant - pi     = ',pi
      write(*,*)' Hrant - nTheta = ',nTheta
      write(*,*)' Hrant - nPhi   = ',nPhi
      write(*,*)
      write(*,*)' Hrant - dTheta = ',dTheta
      write(*,*)' Hrant - dPhi   = ',dPhi
      write(*,*)
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
