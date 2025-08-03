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
      subroutine padCommandLine(iMODyson,nGridPointsTheta,nGridPointsM,  &
        kMag,fafName)
!
!     This routine processes the command line arguments and fills key variables
!     and option flags.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      integer(kind=int64)::iMODyson,nGridPointsTheta,nGridPointsM
      real(kind=real64)::kMag
      character(len=256)::fafName
!
      integer(kind=int64)::i,nCommandLineArgs
!
!     Walk through the command line arguments to fill required input parameters
!     and other option flags.
!
      nCommandLineArgs = command_argument_count()
      if(nCommandLineArgs.lt.1.or.nCommandLineArgs.gt.5)  &
        call mqc_error('PAD expects 1-5 command line arguments.')
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
      else
        do i = 1,nCommandLineArgs

        endDo
      endIf

!
      return
      end subroutine padCommandLine

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
      somx2 = sqrt(max(0.0_real64,1.0_real64 - x*x))
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


      end module pad_mod
