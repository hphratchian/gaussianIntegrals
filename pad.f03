Include "gbs_mod.f03"
      program pad
!
!     This test program evaluates the photoelectron angular distribution,
!     intensities as a function of theta, for a given MO extracted from a FAF.
!     There is 1 required command line argument and 3 optional command line
!     arguments:
!           1. FAF filename (required)
!           2. (alpha) MO number to serve as the Dyson orbital (optional;
!              default=1)
!           3. number of spatial grid points (optional; default=501)
!           4. number of angles to evaluate from 0 --> pi (optional; default=31)
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
      integer(kind=int64),parameter::nOMP=48
      integer(kind=int64)::i,j,iMODyson,nGridPointsM,nGridPointsTheta
      real(kind=real64)::tStart,tEnd,tstart1,tEnd1,stepSizeIntM,  &
        stepSizeTheta,thetaStart,moVal1,moVal2,kMag,MSquared0,MSquared90
      real(kind=real64),dimension(3)::cartStart,cartEnd,integratedIntensity
      real(kind=real64),dimension(3,3)::laserVector
      real(kind=real64),dimension(:),allocatable::quadGridTheta,  &
        quadWeightsTheta,quadWeightsM,basisValues,MSquaredList
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
 1000 format(1x,'Program PAD.')
 1100 format(/,1x,'Current laser field: (',f6.2,',',f6.2,',',f6.2,')   Intensity = ',f25.8)
 2000 format(1x,'Data for the ',A,' grid:',/,  &
        3x,'nPoints = ',i15,' step size = ',f20.5)
 3000 format(/,1x,55('-'),/,  &
        12x,'Intensity as a Function of Theta',/,  &
        18x,'(kMag = ',f8.3,' eV)',/,  &
        9x,'theta (rad)',15x,'I(theta)',/,  &
        1x,55('='))
 3010 format(9x,f7.3,3x,f25.8)
 3020 format(1x,55('='),/)
 3100 format(/,1x,'I(0) = ',f25.8,3x,'I(90) = ',f25.8,/,  &
        1x,'beta = ',f25.8,/)
 3500 format(1x,'Laser field: (',f5.2,',',f5.2,',',f5.2,')   Intensity = ',f10.6)
 3510 format(1x,'Beta(z|',A,') = ',f10.6)
 8999 format(/,1x,'Job Time: ',f15.1,' s')
!
!
!     Begin the program.
!
      call CPU_TIME(tStart)
      call omp_set_num_threads(nOMP)
      fail = .false.
      write(iOut,1000)
!
!     Memory check...
!
      if(MEMChecks) call print_memory_usage(iOut,'At top of PAD.')
!
!     Read the FAF name and MO number (the one we use for the Dyson orbtial)
!     from the command line. If the MO number isn't included, we set it to zero
!     here and later default it to the HOMO.
!
      iMODyson = 0
      nGridPointsM = 0
      nGridPointsTheta = 0
      if(command_argument_count().lt.1.or.  &
        command_argument_count().gt.4)  &
        call mqc_error('PAD expects 1-4 command line arguments.')
      call get_command_argument(1,fafName)
      if(command_argument_count().ge.2)   &
        call mqc_get_command_argument_integer(2,iMODyson)
      if(iMODyson.eq.0) iMODyson = 1
      if(command_argument_count().ge.3)  &
        call mqc_get_command_argument_integer(3,nGridPointsM)
      if(nGridPointsM.eq.0) nGridPointsM = 501
      if(command_argument_count().ge.4)  &
        call mqc_get_command_argument_integer(4,nGridPointsTheta)
      if(nGridPointsTheta.eq.0) nGridPointsTheta = 31
!
!     Load the FAF and set the MO number if it wasn't provided on the command
!     line.
      call faf%load(fafName)
!
!     Set the laser polarization vector.
!
      laserVector(:,1) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      laserVector(:,2) = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
      laserVector(:,3) = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
!
!     Read the basis set and MO coefficients from faf.
!
      call loadGaussianBasisSet(faf,basisSet)
      call faf%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmp)
      moCoeffs = tmp
!
!     Fill in the grid of theta points.
!
      if(MEMChecks) call print_memory_usage(iOut,'Before building grids.')
      thetaStart = mqc_float(0)
      stepSizeTheta = Pi/mqc_float(nGridPointsTheta-1)
      Allocate(quadGridTheta(nGridPointsTheta),  &
        quadWeightsTheta(nGridPointsTheta))
      call setup_quadrature_trapezoid1d(nGridPointsTheta,stepSizeTheta,  &
        thetaStart,quadGridTheta,quadWeightsTheta)
      write(iOut,2000) 'theta',nGridPointsTheta,stepSizeTheta
      flush(iOut)
!
!     Prepare the integration grid and quadrature weights for the M evaluations.
!     There is one M per theta.
!
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
      write(iOut,*)' Time for 3D grid setup = ',tEnd1-tStart1
      flush(iOut)
!
!     Memory check...
!
      if(MEMChecks) call print_memory_usage(iOut,'After building grids.')
!
!     Test that the chosen MO is normalized using quadrature.
!
      call mqc_print(moCoeffs(:,iMODyson),iOut,header='Dyson coefficients')
      call CPU_TIME(tStart1)
      write(iOut,*)
      write(iOut,*)' Test of <dyson|dyson> = ',  &
        moInnerProductNumericalIntegration(moCoeffs(:,iMODyson),  &
        quadGridM,quadWeightsM,basisSet)
      write(iOut,*)
      call CPU_TIME(tEnd1)
      write(iOut,*)' Time for norm test = ',tEnd1-tStart1
      flush(iOut)
!
!     Integrate the Dyson/plane wave matrix elements as a function of theta.
!
      kMag = mqc_float(1)/mqc_float(500)
      do i = 1,3
        write(iOut,1100) laserVector(:,i)
        call CPU_TIME(tStart1)
        MSquared0 = dysonPlaneWaveMatrixElementSquared(mqc_float(0),kMag,  &
          laserVector(:,i),moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
        call CPU_TIME(tEnd1)
        write(iOut,*)' Time for MSquared0  = ',tEnd1-tStart1
        call CPU_TIME(tStart1)
        MSquared90 = dysonPlaneWaveMatrixElementSquared(Pi/mqc_float(2),  &
          kMag,laserVector(:,i),moCoeffs(:,iMODyson),basisSet,quadGridM,  &
          quadWeightsM)
        call CPU_TIME(tEnd1)
        write(iOut,*)' Time for MSquared90 = ',tEnd1-tStart1
        call CPU_TIME(tStart1)
        call mqc_print(laserVector(:,i),iOut,  &
          header='electric field vector...',blank_at_top=.true.)
        flush(iOut)
        MSquaredList = dysonPlaneWaveMatrixElementSquaredThetaList(  &
          quadGridTheta,kMag,  &
          laserVector(:,i),moCoeffs(:,iMODyson),basisSet,quadGridM,quadWeightsM)
        call CPU_TIME(tEnd1)
        write(iOut,*)' Time for MSquared(theta) list = ',tEnd1-tStart1
        write(iOUt,*)
        flush(iOut)
!
!       Print the intensity data table.
!
        write(iOut,3000) kMag*evPHartree
        do j = 1,nGridPointsTheta
          write(iOut,3010) quadGridTheta(j),MSquaredList(j)
        endDo
        write(iOut,3020)
        write(iOut,*)' Integrated I = ',dot_product(quadGridTheta,quadWeightsTheta)
        integratedIntensity(i) = dot_product(quadGridTheta,quadWeightsTheta)
        write(iOut,*)
!
!       Evaluate the anisotropy parameter, beta, using the computed intensities
!       at 0 and 90 degrees.
!
        write(iOut,3100) MSquared0,MSquared90,  &
          (MSquared0-MSquared90)/((MSquared0/mqc_float(2))+MSquared90)
        flush(iOut)
      endDo
!
!     Print out the integrated planar intensity for each laser vector. Then
!     compute an analytic beta.
!
      do i = 1,3
        write(iOut,3500) laserVector(:,i),integratedIntensity(i)
      endDo
      write(iOut,3510) 'x',(integratedIntensity(3)-integratedIntensity(1))/  &
        (integratedIntensity(3)+integratedIntensity(1)/mqc_float(2))
      write(iOut,3510) 'y',(integratedIntensity(3)-integratedIntensity(2))/  &
        (integratedIntensity(3)+integratedIntensity(2)/mqc_float(2))
!
!     The end of the program.
!
  999 Continue
      call CPU_TIME(tEnd)
      write(iOut,8999) tEnd-tStart
      if(MEMChecks) call print_memory_usage(iOut,'End of PAD.')
      end program pad
