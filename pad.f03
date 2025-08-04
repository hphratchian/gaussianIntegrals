      program pad
!
!     This test program evaluates the photoelectron angular distribution,
!     intensities as a function of theta, for a given MO extracted from a FAF.
!     There is 1 required command line argument and 3 optional command line
!     arguments:
!           1. FAF filename (required)
!           2. (alpha) MO number to serve as the Dyson orbital (optional;
!              default=1)
!           3. magnitude of the k-vector in eV
!           4. number of angles to evaluate from 0 --> pi (optional; default=15)
!           5. number of spatial grid points (optional; default=101)
!           6. iPEType flag (optional; default=0)
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
      use pad_mod
      implicit none
      integer(kind=int64),parameter::nOMP=1,lMax=2
      logical,parameter::extraPrint=.false.
      integer(kind=int64)::i,j,k,iMODyson,nGridPointsM,  &
        nGridPointsTheta,nGridPointsPhi,nIntPlanes,iPEType
      real(kind=real64)::tStart,tEnd,tstart1,tEnd1,stepSizeIntM,  &
        stepSizeTheta,thetaStart,stepSizePhi,moVal1,moVal2,kMag,  &
        MSquared0,MSquared90
      real(kind=real64),dimension(3)::cartStart,cartEnd
      real(kind=real64),dimension(0:lMax)::lWeights0,lWeights90,lWeights
      real(kind=real64),dimension(:),allocatable::integratedIntensity,  &
        betaValsParaPerp,betaValsFit,rSquared
      real(kind=real64),dimension(:,:),allocatable::laserVector,  &
        orthogPlaneVector
      real(kind=real64),dimension(:),allocatable::quadGridTheta,  &
        quadWeightsTheta,quadWeightsM,basisValues,MSquaredList,  &
        quadGridPhi,MSquaredList1
      real(kind=real64),dimension(:,:),allocatable::quadGridM,moCoeffs,  &
        MSquaredList2
      logical::fail=.false.,found
      character(len=256)::fafName
      character(len=2),dimension(:),allocatable::intPlaneLabels
      type(mqc_gaussian_unformatted_matrix_file)::faf
      type(mqc_basisset)::basisSet
      type(MQC_Variable)::tmp
      real(kind=real64),dimension(3)::v1,v2,v3
!
      real(kind=real64)::tmpBeta,tmpR2
      real(kind=real64),dimension(:,:),allocatable::basisIntegrals,  &
        overlapMatrix
      real(kind=real64),dimension(:),allocatable::quadValues
!
!     Format statements.
!
 1000 format(1x,'Program PAD.')
 1100 format(1x,'Current laser field: (',f6.2,',',f6.2,',',f6.2,')')
 2000 format(1x,'Data for the ',A,' grid:',/,  &
        3x,'nPoints = ',i15,' step size = ',f20.5)
 2500 format(/,1x,'Dyson orbital self-overlap = ',f12.8)
 3000 format(/,1x,55('-'),/,  &
        12x,'Intensity as a Function of Theta',/,  &
        18x,'(kMag = ',f8.3,' eV)',/,  &
        9x,'theta (rad)',15x,'I(theta)',/,  &
        1x,55('='))
 3010 format(9x,f7.3,3x,f25.8)
 3020 format(1x,55('='))
 3100 format(1x,'I(0) = ',f12.6,3x,'I(90) = ',f12.6,3x,'beta = ',f12.6)
 3500 format(/,1x,87('-'),/,  &
        31x,'Summary of PAD Calculation',/,  &
        34x,'(kMag = ',f8.3,' eV)',/,  &
        1x,87('='))
 3510 format(1x,'Integration Plane: ',A,'  |  Laser field: (',  &
        f5.2,',',f5.2,',',f5.2,')  |  Intensity = ',f10.6,/,  &
        20x,'beta(ratio) = ',f10.6,'  |  beta(fit) = ',f10.6,  &
        ' (R**2 = ',f8.5,')')
 3520 format(1x,87('='))
 3530 format(1x,'Average Beta (ratio) = ',f10.6,/,  &
        1x,'Average Beta (fit)   = ',f10.6,/)
 8998 format(1x,'Time for ',A,': ',f15.1,' s')
 8999 format(/,1x,'Job Time: ',f15.1,' s',/,1x,'PAD Complete.')
!
!
!     Begin the program.
!
      call CPU_TIME(tStart)
      call omp_set_num_threads(nOMP)
      fail = .false.
      write(iOut,1000)
      call mqc_version_print(iOut)
!
!     Memory check...
!
      if(MEMChecks) call print_memory_usage(iOut,'At top of PAD.')
!
!     Call padCommandLine to parse the command line and fill input and other job
!     option flags.
!
      call padCommandLine(iPEType,iMODyson,nGridPointsTheta,  &
        nGridPointsM,kMag,fafName)
!
!     Load the FAF.
!
      call faf%load(fafName)
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
!
!     Read the basis set and MO coefficients from faf.
!
      call loadGaussianBasisSet(faf,basisSet)
      call faf%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmp)
      moCoeffs = tmp
!
!     Prepare the integration grid and quadrature weights for the M evaluations.
!     There are two quadratures coded here: (1) 3D trapezoid quadrature; and (2)
!     XC quadrature grid saved on the Gaussian FAF.
!
      if(.false.) then
        cartStart = [ -6.0,-6.0,-6.0 ]
        cartEnd = [ 6.0,6.0,6.0 ]
        stepSizeIntM = (cartEnd(1)-cartStart(1))/mqc_float(nGridPointsM-1)
        Allocate(quadGridM(3,nGridPointsM**3),quadWeightsM(nGridPointsM**3),  &
          quadValues(nGridPointsM**3))
        call CPU_TIME(tStart1)
        call setup_quadrature_trapezoid3d(nGridPointsM,stepSizeIntM,  &
          cartStart,quadGridM,quadWeightsM)
        write(iOut,2000) 'M quadrature',nGridPointsM**3,stepSizeIntM
        call CPU_TIME(tEnd1)
        write(iOut,8998) '3D grid setup',tEnd1-tStart1
        flush(iOut)
      else
        write(iOut,*)
        write(iOut,*)' Hrant - calling MQC get_quad routine...'
        call MQC_Gaussian_FAF_Get_3DQuadratureGrid(faf,  &
          quadWeightsM,quadGridM,found)
        write(iOut,*)' Hrant - back from MQC get_quad routine...'
        write(iOut,*)'         found               = ',found
        write(iOut,*)'         Alloc(quadWeightsM) = ',Allocated(quadWeightsM)
        write(iOut,*)'         Alloc(quadGridM)    = ',Allocated(quadGridM)
        write(iOut,*)'         No. of grid points  = ',SIZE(quadWeightsM)
        write(iOut,*)
      endIf
!
!     Fill in the grid of I(theta) points.
!
      if(MEMChecks) call print_memory_usage(iOut,'Before building grids.')
      thetaStart = mqc_float(0)
      stepSizeTheta = Pi/mqc_float(nGridPointsTheta-1)
      Allocate(quadGridTheta(nGridPointsTheta),  &
        quadWeightsTheta(nGridPointsTheta))
      call setup_quadrature_trapezoid1d(nGridPointsTheta,stepSizeTheta,  &
        thetaStart,quadGridTheta,quadWeightsTheta)
      write(iOut,2000) 'theta',nGridPointsTheta,stepSizeTheta
      call mqc_print(quadGridTheta,iOut,header='Theta Grid')
      call mqc_print(quadGridTheta/Pi,iOut,header='Theta Grid/Pi')
      flush(iOut)
      nGridPointsPhi = 2*nGridPointsTheta-1
      stepSizePhi = mqc_float(2)*Pi/mqc_float(nGridPointsPhi-1)
      Allocate(quadGridPhi(nGridPointsPhi))
      quadGridPhi(1) = mqc_float(0)
      do i = 2,nGridPointsPhi
        quadGridPhi(i) = quadGridPhi(i-1) + stepSizeTheta
      endDo
      call mqc_print(quadGridPhi,iOut,header='phi grid')
      call mqc_print(quadGridPhi/Pi,iOut,header='phi grid/Pi')
      write(iOut,*)
      write(iOut,*)' Hrant - stepSizeTheta = ',stepSizeTheta
      write(iOut,*)'         stepSizePhi   = ',stepSizePhi
      write(iOut,*)
!
!     Memory check...
!
      if(MEMChecks) call print_memory_usage(iOut,'After building grids.')
!
!     Test that the chosen MO is normalized using quadrature.
!
      if(extraPrint)  &
        call mqc_print(moCoeffs(:,iMODyson),iOut,header='Dyson coefficients')
      call CPU_TIME(tStart1)
      write(iOut,2500)  &
        moInnerProductNumericalIntegration(moCoeffs(:,iMODyson),  &
        quadGridM,quadWeightsM,basisSet)
      call CPU_TIME(tEnd1)
      write(iOut,8998) 'norm test',tEnd1-tStart1
      flush(iOut)
!
!     Integrate the Dyson/plane wave matrix elements as a function of theta.
!

!hph+
      write(iOut,*)
      write(iOut,*)' Hrant - iPEType = ',iPEType
      lWeights0 = -1.0
      lWeights90 = -1.0
!hph-

      do i = 1,nIntPlanes
        write(iOut,*)
        write(iOut,*)' Hrant - i = ',i
        write(iOut,1100) laserVector(:,i)
        call CPU_TIME(tStart1)

!hph+
        if(iPEType.eq.0) then
          MSquared0 = dysonPlaneWaveMatrixElementSquared(mqc_float(0),kMag,  &
            laserVector(:,i),orthogPlaneVector(:,i),moCoeffs(:,iMODyson),  &
            basisSet,quadGridM,quadWeightsM)
        else
          call dysonMatrixElement1Angle(iPEType,lMax,mqc_float(0),kMag,  &
            laserVector(:,i),orthogPlaneVector(:,i),moCoeffs(:,iMODyson),  &
            basisSet,quadGridM,quadWeightsM,MSquared0,lWeights0)
        endIf
        write(iOut,*)' Hrant - MSquared0  = ',mSquared0
        call mqc_print(lWeights0,iOut,header='l weights (0)')
!hph-

        call CPU_TIME(tEnd1)
        write(iOut,8998) 'MSquared0 ',tEnd1-tStart1
        call CPU_TIME(tStart1)

!hph+
        if(iPEType.eq.0) then
          MSquared90 = dysonPlaneWaveMatrixElementSquared(Pi/mqc_float(2),  &
            kMag,laserVector(:,i),orthogPlaneVector(:,i),  &
            moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
        else
          call dysonMatrixElement1Angle(iPEType,lMax,Pi/mqc_float(2),kMag,  &
            laserVector(:,i),orthogPlaneVector(:,i),moCoeffs(:,iMODyson),  &
            basisSet,quadGridM,quadWeightsM,MSquared90,lWeights90)
        endIf
        write(iOut,*)' Hrant - MSquared90 = ',mSquared90
        call mqc_print(lWeights90,iOut,header='l weights (90)')
!hph-

        call CPU_TIME(tEnd1)
        write(iOut,8998) 'MSquared90',tEnd1-tStart1
        call CPU_TIME(tStart1)
        flush(iOut)

!hph+
        goto 999
!hph-

        if(iPEType.eq.0) then
          MSquaredList = dysonPlaneWaveMatrixElementSquaredThetaList(  &
            quadGridTheta,kMag,laserVector(:,i),orthogPlaneVector(:,i),  &
            moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
        else

          if(Allocated(mSquaredList)) DeAllocate(mSquaredList)
          Allocate(mSquaredList(nGridPointsTheta))

!          call dysonMatrixElementThetaList(iPEType,lMax,quadGridTheta,kMag,  &
!            laserVector(:,i),orthogPlaneVector(:,i),moCoeffs(:,iMODyson),  &
!            basisSet,quadGridM,quadWeightsM,MSquaredList,lWeights)
          call dysonMatrixElementThetaList(iPEType,lMax,quadGridTheta,kMag,  &
            laserVector(:,i),orthogPlaneVector(:,i),moCoeffs(:,iMODyson),  &
            basisSet,quadGridM,quadWeightsM,MSquaredList)
        endIf
        call mqc_print(MSquaredList,iOut,header='MSquaredList',blank_at_top=.true.)


!hph+
!        MSquaredList2 = dysonPlaneWaveMatrixElementSquaredPhiThetaList(  &
!          quadGridPhi,quadGridTheta,kMag,laserVector(:,i),  &
!          orthogPlaneVector(:,i),moCoeffs(:,iMODyson),basisSet,  &
!          quadGridM,quadWeightsM)
!        call mqc_print(MSquaredList2,iOut,header='MSquaredList2',blank_at_top=.true.)
!        do j = 1,nGridPointsPhi
!          call betaLeastSquares(MSquaredList2(j,:),quadGridTheta,  &
!            tmpBeta,tmpR2)
!          write(iOut,*) quadGridPhi(j)/Pi,tmpBeta,tmpR2
!        endDo
!hph-

        call CPU_TIME(tEnd1)
        write(iOut,8998) 'MSquared(theta) list',tEnd1-tStart1
        flush(iOut)
!
!       Compute beta using a least squares fitting scheme.
!
        call betaLeastSquares(MSquaredList,quadGridTheta,betaValsFit(i),  &
          rSquared(i))
        write(iOut,*)' Hrant - Saving betaValsFit(i) and rSquared(i): ',  &
          betaValsFit(i),rSquared(i)
        write(iOut,*)
!
!       Print the intensity data table.
!
        write(iOut,3000) kMag/evPHartree
        do j = 1,nGridPointsTheta
          write(iOut,3010) quadGridTheta(j),MSquaredList(j)
        endDo
        write(iOut,3020)
        integratedIntensity(i) = dot_product(MSquaredList,quadWeightsTheta)
!
!       Evaluate the anisotropy parameter, beta, using the computed intensities
!       at 0 and 90 degrees.
!
        betaValsParaPerp(i) = betaParaPerp(MSquared0,MSquared90)
        write(iOut,3100) MSquared0,MSquared90,  &
          betaParaPerp(MSquared0,MSquared90)
        flush(iOut)

!hph+
!        if(.not.Allocated(MSquaredList1))  &
!          Allocate(MSquaredList1(Size(MSquaredList)))
!        MSquaredList1 = mqc_float(0)
!        do j = 1,nGridPointsTheta
!          do k = 1,nGridPointsPhi
!            MSquaredList1(j) = MSquaredList1(j)+MSquaredList2(k,j)*stepSizeTheta
!          endDo
!        endDo
!        call mqc_print(MSquaredList1,iOut,header='MSquaredList1')
!        call betaLeastSquares(MSquaredList1,quadGridTheta,tmpBeta,tmpR2)
!        write(iOut,*)' beta,R**2:',tmpBeta,tmpR2
!        write(iOut,*)
!hph-

      endDo
!
!     Print out the integrated planar intensity for each laser vector. Then
!     compute an analytic beta.
!
      write(iOut,3500) kMag
      do i = 1,nIntPlanes
        write(iOut,3510) TRIM(intPlaneLabels(i)),laserVector(:,i),  &
          integratedIntensity(i),betaValsParaPerp(i),betaValsFit(i),  &
          rSquared(i)
      endDo
      write(iOut,3520)
      write(iOut,3530) SUM(betaValsParaPerp)/mqc_float(nIntPlanes),  &
        SUM(betaValsFit)/mqc_float(nIntPlanes)

!hph+
      goto 999
!hph-

!
!     The end of the program.
!
  999 Continue
      call CPU_TIME(tEnd)
      write(iOut,8999) tEnd-tStart
      if(MEMChecks) call print_memory_usage(iOut,'End of PAD.')
      end program pad
