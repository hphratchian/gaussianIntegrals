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
      use dyson_matrix_elements_mod
!
      implicit none
      integer(kind=int64),parameter::PAD_LAB_FRAMES_CARTESIAN=0_int64
      integer(kind=int64),parameter::PAD_LAB_FRAMES_SPHERE=1_int64
      integer(kind=int64),parameter::PAD_LAB_FRAMES_AXISYMMETRIC=2_int64
      integer(kind=int64),parameter::PAD_LAB_FRAMES_CUSTOM=-1_int64
!
!
!     The pad_options object collects user-facing and model-control options for
!     a PAD calculation. The photon and binding energies are supplied in eV.
      type pad_options
        integer(kind=int64)::dysonMOIndex=1_int64
        integer(kind=int64)::nGridPointsTheta=5_int64
        integer(kind=int64)::nGridPointsM=101_int64
        integer(kind=int64)::nChi=36_int64
        integer(kind=int64)::iPEType=0_int64
        integer(kind=int64)::lMax=6_int64
        integer(kind=int64)::nOMP=1_int64
        integer(kind=int64)::labFrameType=PAD_LAB_FRAMES_CARTESIAN
        integer(kind=int64)::nLabFrameTheta=5_int64
        integer(kind=int64)::nLabFramePhi=8_int64
        logical::printResults=.true.
        logical::printThetaTable=.true.
        real(kind=real64)::photonEnergyEV=0.0_real64
        real(kind=real64)::bindingEnergyEV=0.0_real64
        real(kind=real64)::labFrameAlignment=0.0_real64
        real(kind=real64),dimension(:,:),allocatable::labEpsilonVector
        real(kind=real64),dimension(:,:),allocatable::labKPlaneVector
        real(kind=real64),dimension(:),allocatable::labFrameWeights
        character(len=8),dimension(:),allocatable::labFrameLabels
      end type pad_options
!
!     The pad_results object stores the central numerical results from a PAD
!     calculation so that other codes can call the driver as a routine.
      type pad_results
        integer(kind=int64)::nLabFrames=0_int64
        integer(kind=int64)::nTheta=0_int64
        integer(kind=int64)::nChi=0_int64
        real(kind=real64)::photoelectronEnergyEV=0.0_real64
        real(kind=real64)::photoelectronEnergyHartree=0.0_real64
        real(kind=real64)::kMag=0.0_real64
        real(kind=real64)::dysonSelfOverlap=0.0_real64
        real(kind=real64)::labFrameWeightSum=0.0_real64
        real(kind=real64)::chiWeightSum=0.0_real64
        real(kind=real64)::averageIntensity0=0.0_real64
        real(kind=real64)::averageIntensity90=0.0_real64
        real(kind=real64)::averageThetaIntegratedIntensity=0.0_real64
        real(kind=real64)::averageSolidAngleIntegratedIntensity=0.0_real64
        real(kind=real64)::averageBetaParaPerp=0.0_real64
        real(kind=real64)::averageBetaFit=0.0_real64
        real(kind=real64)::averageRSquared=0.0_real64
        real(kind=real64)::meanBetaParaPerp=0.0_real64
        real(kind=real64)::meanBetaFit=0.0_real64
        real(kind=real64),dimension(:),allocatable::theta
        real(kind=real64),dimension(:),allocatable::thetaWeights
        real(kind=real64),dimension(:),allocatable::thetaSolidAngleWeights
        real(kind=real64),dimension(:),allocatable::chi
        real(kind=real64),dimension(:),allocatable::chiWeights
        real(kind=real64),dimension(:),allocatable::labFrameWeights
        real(kind=real64),dimension(:),allocatable::thetaIntegratedIntensity
        real(kind=real64),dimension(:),allocatable::solidAngleIntegratedIntensity
        real(kind=real64),dimension(:),allocatable::intensity0
        real(kind=real64),dimension(:),allocatable::intensity90
        real(kind=real64),dimension(:),allocatable::betaValsParaPerp
        real(kind=real64),dimension(:),allocatable::betaValsFit
        real(kind=real64),dimension(:),allocatable::rSquared
        real(kind=real64),dimension(:),allocatable::averageIntensityTheta
        real(kind=real64),dimension(:,:),allocatable::intensityTheta
        real(kind=real64),dimension(:,:),allocatable::epsilonVector
        real(kind=real64),dimension(:,:),allocatable::kPlaneVector
        character(len=8),dimension(:),allocatable::intPlaneLabels
      end type pad_results
!
!
      CONTAINS
!
!
!PROCEDURE padCommandLine
      subroutine padCommandLine(options,fafName)
!
!     This routine processes the command line arguments and fills the PAD
!     options object and FAF filename.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      type(pad_options),intent(out)::options
      character(len=256),intent(out)::fafName
!
      integer(kind=int64)::nCommandLineArgs
!
!     Walk through the command line arguments to fill required input parameters
!     and other option flags.
!
      options = pad_options()
      nCommandLineArgs = command_argument_count()
      if(nCommandLineArgs.gt.0) then
        call padDispatchCommandLine(options,fafName,nCommandLineArgs)
      else
        call padPrintUsage()
        call mqc_error('PAD expects command line arguments.')
      endIf
!
      return
      end subroutine padCommandLine


!PROCEDURE padDispatchCommandLine
      subroutine padDispatchCommandLine(options,fafName,nCommandLineArgs)
!
!     This routine dispatches command-line parsing to the named-option parser
!     or to the legacy positional parser.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      type(pad_options),intent(inout)::options
      character(len=256),intent(out)::fafName
      integer(kind=int64),intent(in)::nCommandLineArgs
!
      character(len=256)::arg
!
      call get_command_argument(1,arg)
      if(TRIM(arg).eq.'-help'.or.TRIM(arg).eq.'--help'.or.  &
        TRIM(arg).eq.'-h') then
        call padPrintUsage()
        stop
      endIf
      if(arg(1:1).eq.'-') then
        call padNamedCommandLine(options,fafName,nCommandLineArgs)
      else
        call padLegacyCommandLine(options,fafName,nCommandLineArgs)
      endIf
!
      return
      end subroutine padDispatchCommandLine


!PROCEDURE padLegacyCommandLine
      subroutine padLegacyCommandLine(options,fafName,nCommandLineArgs)
!
!     This routine processes the legacy positional command-line arguments.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      type(pad_options),intent(inout)::options
      character(len=256),intent(out)::fafName
      integer(kind=int64),intent(in)::nCommandLineArgs
!
      if(nCommandLineArgs.lt.4.or.nCommandLineArgs.gt.11)  &
        call mqc_error('PAD expects 4-11 command line arguments.')
      call get_command_argument(1,fafName)
      call mqc_get_command_argument_integer(2,options%dysonMOIndex)
      call mqc_get_command_argument_real(3,options%photonEnergyEV)
      call mqc_get_command_argument_real(4,options%bindingEnergyEV)
      if(nCommandLineArgs.ge.5) then
        call mqc_get_command_argument_integer(5,options%nGridPointsTheta)
      endIf
      if(nCommandLineArgs.ge.6) then
        call mqc_get_command_argument_integer(6,options%nGridPointsM)
      endIf
      if(nCommandLineArgs.ge.7) then
        call mqc_get_command_argument_integer(7,options%iPEType)
      endIf
      if(nCommandLineArgs.ge.8) then
        call mqc_get_command_argument_integer(8,options%labFrameType)
      endIf
      if(nCommandLineArgs.ge.9) then
        call mqc_get_command_argument_integer(9,options%nLabFrameTheta)
      endIf
      if(nCommandLineArgs.ge.10) then
        call mqc_get_command_argument_integer(10,options%nLabFramePhi)
      endIf
      if(nCommandLineArgs.ge.11) then
        call mqc_get_command_argument_integer(11,options%nChi)
      endIf
!
      return
      end subroutine padLegacyCommandLine


!PROCEDURE padNamedCommandLine
      subroutine padNamedCommandLine(options,fafName,nCommandLineArgs)
!
!     This routine processes named command-line options supplied as
!     -option value pairs.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      type(pad_options),intent(inout)::options
      character(len=256),intent(out)::fafName
      integer(kind=int64),intent(in)::nCommandLineArgs
!
      integer(kind=int64)::i,equalPos
      logical::gotFaf,gotMO,gotPhotonEV,gotBindingEV
      character(len=256)::arg,value
!
      gotFaf = .false.
      gotMO = .false.
      gotPhotonEV = .false.
      gotBindingEV = .false.
      fafName = ''
      i = 1
      do while(i.le.nCommandLineArgs)
        call get_command_argument(i,arg)
        if(TRIM(arg).eq.'-help'.or.TRIM(arg).eq.'--help'.or.  &
          TRIM(arg).eq.'-h') then
          call padPrintUsage()
          stop
        endIf
        if(arg(1:1).ne.'-')  &
          call mqc_error('PAD: expected named option beginning with "-".')
        equalPos = INDEX(arg,'=')
        if(equalPos.gt.0) then
          value = arg(equalPos+1:)
          arg = arg(1:equalPos-1)
        else
          if(i.eq.nCommandLineArgs)  &
            call mqc_error('PAD: missing value after option '//TRIM(arg)//'.')
          i = i+1
          call get_command_argument(i,value)
        endIf
        call padSetNamedOption(options,fafName,arg,value,gotFaf,gotMO,  &
          gotPhotonEV,gotBindingEV)
        i = i+1
      endDo
!
      if(.not.gotFaf) call mqc_error('PAD: missing required -faf option.')
      if(.not.gotMO) call mqc_error('PAD: missing required -dyson-mo option.')
      if(.not.gotPhotonEV)  &
        call mqc_error('PAD: missing required -photon-ev option.')
      if(.not.gotBindingEV)  &
        call mqc_error('PAD: missing required -binding-ev option.')
!
      return
      end subroutine padNamedCommandLine


!PROCEDURE padSetNamedOption
      subroutine padSetNamedOption(options,fafName,arg,value,gotFaf,gotMO,  &
        gotPhotonEV,gotBindingEV)
!
!     This routine assigns one parsed named command-line option.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      type(pad_options),intent(inout)::options
      character(len=256),intent(inout)::fafName
      character(len=*),intent(in)::arg,value
      logical,intent(inout)::gotFaf,gotMO,gotPhotonEV,gotBindingEV
!
      character(len=64)::key
!
      call padOptionKey(arg,key)
      select case(TRIM(key))
      case('faf','faffile','file')
        fafName = TRIM(value)
        gotFaf = .true.
      case('dysonmo','mo','moindex','dysonmoindex')
        call padReadIntegerOption(arg,value,options%dysonMOIndex)
        gotMO = .true.
      case('photonev','photonenergy','photonenergyev')
        call padReadRealOption(arg,value,options%photonEnergyEV)
        gotPhotonEV = .true.
      case('bindingev','bindingenergy','bindingenergyev')
        call padReadRealOption(arg,value,options%bindingEnergyEV)
        gotBindingEV = .true.
      case('ntheta','ngridpointstheta','theta')
        call padReadIntegerOption(arg,value,options%nGridPointsTheta)
      case('ngrid','ngridpointsm','mgrid')
        call padReadIntegerOption(arg,value,options%nGridPointsM)
      case('petype','ipetype','photoelectronmodel')
        call padReadIntegerOption(arg,value,options%iPEType)
      case('labframe','labframetype')
        call padReadLabFrameOption(arg,value,options%labFrameType)
      case('labtheta','nlabtheta','nlabframetheta')
        call padReadIntegerOption(arg,value,options%nLabFrameTheta)
      case('labphi','nlabphi','nlabframephi')
        call padReadIntegerOption(arg,value,options%nLabFramePhi)
      case('labalignment','alignment','labalign','labalignp2',  &
        'labalignmentp2','alignmentp2')
        call padReadRealOption(arg,value,options%labFrameAlignment)
      case('nchi','chi')
        call padReadIntegerOption(arg,value,options%nChi)
      case('lmax')
        call padReadIntegerOption(arg,value,options%lMax)
      case('nomp','omp','threads')
        call padReadIntegerOption(arg,value,options%nOMP)
      case default
        call mqc_error('PAD: unknown option '//TRIM(arg)//'.')
      end select
!
      return
      end subroutine padSetNamedOption


!PROCEDURE padReadIntegerOption
      subroutine padReadIntegerOption(arg,value,iValue)
!
!     This routine reads an integer command-line option value.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::arg,value
      integer(kind=int64),intent(out)::iValue
!
      integer::iStat
!
      read(value,*,iostat=iStat) iValue
      if(iStat.ne.0)  &
        call mqc_error('PAD: invalid integer value for '//TRIM(arg)//'.')
!
      return
      end subroutine padReadIntegerOption


!PROCEDURE padReadRealOption
      subroutine padReadRealOption(arg,value,rValue)
!
!     This routine reads a real command-line option value.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::arg,value
      real(kind=real64),intent(out)::rValue
!
      integer::iStat
!
      read(value,*,iostat=iStat) rValue
      if(iStat.ne.0)  &
        call mqc_error('PAD: invalid real value for '//TRIM(arg)//'.')
!
      return
      end subroutine padReadRealOption


!PROCEDURE padReadLabFrameOption
      subroutine padReadLabFrameOption(arg,value,labFrameType)
!
!     This routine reads a lab-frame model option as text or as an integer.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::arg,value
      integer(kind=int64),intent(out)::labFrameType
!
      character(len=64)::key
!
      call padOptionKey(value,key)
      select case(TRIM(key))
      case('cartesian','cart','0')
        labFrameType = PAD_LAB_FRAMES_CARTESIAN
      case('sphere','spherical','spheregrid','1')
        labFrameType = PAD_LAB_FRAMES_SPHERE
      case('axisymmetric','axial','aligned','alignment','2')
        labFrameType = PAD_LAB_FRAMES_AXISYMMETRIC
      case('custom','user','userprovided','-1')
        labFrameType = PAD_LAB_FRAMES_CUSTOM
      case default
        call padReadIntegerOption(arg,value,labFrameType)
      end select
!
      return
      end subroutine padReadLabFrameOption


!PROCEDURE padOptionKey
      subroutine padOptionKey(text,key)
!
!     This routine normalizes command-line option names for matching.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::text
      character(len=*),intent(out)::key
!
      integer::i,j,code
      character(len=1)::c
!
      key = ''
      j = 0
      do i = 1,LEN_TRIM(text)
        c = text(i:i)
        if(c.eq.'-'.or.c.eq.'_'.or.c.eq.' ') cycle
        code = IACHAR(c)
        if(code.ge.IACHAR('A').and.code.le.IACHAR('Z'))  &
          c = ACHAR(code+IACHAR('a')-IACHAR('A'))
        if(j.lt.LEN(key)) then
          j = j+1
          key(j:j) = c
        endIf
      endDo
!
      return
      end subroutine padOptionKey


!PROCEDURE padPrintUsage
      subroutine padPrintUsage()
!
!     This routine prints a short command-line usage summary.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
!
      write(iOut,'(1x,A)')  &
        'Usage: ./pad.exe -faf FILE -dyson-mo N -photon-ev EV -binding-ev EV'
      write(iOut,'(8x,A)')  &
        '[-n-theta N] [-n-grid N] [-pe-type N]'
      write(iOut,'(8x,A)')  &
        '[-lab-frame cartesian|sphere|axisymmetric] [-lab-theta N] [-lab-phi N]'
      write(iOut,'(8x,A)')  &
        '[-lab-alignment A] [-n-chi N] [-lmax N] [-threads N]'
      write(iOut,'(1x,A)')  &
        'Legacy positional arguments are still accepted.'
!
      return
      end subroutine padPrintUsage


!PROCEDURE padPrintOpenMPSettings
      subroutine padPrintOpenMPSettings(options)
!
!     This routine prints the OpenMP settings visible to the PAD driver.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      type(pad_options),intent(in)::options
!
 1000 format(1x,'OpenMP requested threads = ',i8)
 1010 format(1x,'OpenMP max threads       = ',i8)
!
      write(iOut,1000) options%nOMP
      write(iOut,1010) omp_get_max_threads()
!
      return
      end subroutine padPrintOpenMPSettings


!PROCEDURE padPrintReproducibleCommand
      subroutine padPrintReproducibleCommand(fafName,options)
!
!     This routine prints a command line that reproduces the current PAD CLI
!     calculation options.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      character(len=*),intent(in)::fafName
      type(pad_options),intent(in)::options
!
 1000 format(/,1x,'Reproducible command line:',/)
 1010 format(3x,'./pad.exe -faf ',A,' \')
 1020 format(5x,'-dyson-mo ',i0,' \')
 1030 format(5x,'-photon-ev ',es24.16,' \')
 1040 format(5x,'-binding-ev ',es24.16,' \')
 1050 format(5x,'-n-theta ',i0,' -n-grid ',i0,' -pe-type ',i0,' \')
 1060 format(5x,'-lab-frame ',i0,' -lab-theta ',i0,' -lab-phi ',i0,' \')
 1070 format(5x,'-lab-alignment ',es24.16,' -n-chi ',i0,' \')
 1080 format(5x,'-lmax ',i0,' -threads ',i0,/)
!
      write(iOut,1000)
      write(iOut,1010) TRIM(fafName)
      write(iOut,1020) options%dysonMOIndex
      write(iOut,1030) options%photonEnergyEV
      write(iOut,1040) options%bindingEnergyEV
      write(iOut,1050) options%nGridPointsTheta,options%nGridPointsM,  &
        options%iPEType
      write(iOut,1060) options%labFrameType,options%nLabFrameTheta,  &
        options%nLabFramePhi
      write(iOut,1070) options%labFrameAlignment,options%nChi
      write(iOut,1080) options%lMax,options%nOMP
!
      return
      end subroutine padPrintReproducibleCommand

!
!PROCEDURE padComputePhotoelectronEnergyAndK
      subroutine padComputePhotoelectronEnergyAndK(options,  &
        photoelectronEnergyEV,  &
        photoelectronEnergyHartree,kMag)
!
!     This routine computes the photoelectron kinetic energy and plane-wave
!     magnitude from the user-facing photon and binding energies.
!
!
!     H. P. Hratchian, 2025.
!
      implicit none
      type(pad_options),intent(in)::options
      real(kind=real64),intent(out)::photoelectronEnergyEV,  &
        photoelectronEnergyHartree,kMag
!
!     Convert the user-facing energy inputs to the plane-wave magnitude used by
!     the transition-moment routines.
!
      if(options%photonEnergyEV.le.mqc_float(0))  &
        call mqc_error('PAD: photon energy must be positive.')
      if(options%bindingEnergyEV.le.mqc_float(0))  &
        call mqc_error('PAD: binding energy must be positive.')
      photoelectronEnergyEV = options%photonEnergyEV-options%bindingEnergyEV
      if(photoelectronEnergyEV.le.mqc_float(0))  &
        call mqc_error('PAD: photon energy must exceed binding energy.')
      photoelectronEnergyHartree = photoelectronEnergyEV/evPHartree
      kMag = sqrt(mqc_float(2)*photoelectronEnergyHartree)
!
      return
      end subroutine padComputePhotoelectronEnergyAndK


!PROCEDURE runPADCalculation
      subroutine runPADCalculation(faf,options,results)
!
!     This routine carries out the PAD calculation for a selected Dyson orbital
!     from an already-loaded Gaussian FAF object. The molecule is kept fixed and
!     the lab-frame polarization vectors are varied. For each polarization
!     direction, the scan plane can also be rotated uniformly about epsilon
!     through a periodic chi quadrature.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      type(mqc_gaussian_unformatted_matrix_file),intent(inout)::faf
      type(pad_options),intent(in)::options
      type(pad_results),intent(out)::results
!
      integer(kind=int64)::i,j,iChi,nIntPlanes,indexTheta90
      real(kind=real64)::tStart1,tEnd1,stepSizeIntM,stepSizeTheta,  &
        thetaStart,MSquared0,MSquared90,tStart2,tEnd2,  &
        timeSingleAngles,timeThetaList,thetaMatchTol
      real(kind=real64),dimension(3)::cartStart,cartEnd,epsilonVec,  &
        epsilonUnit,uBasis,vBasis,uChi
      real(kind=real64),dimension(:),allocatable::lWeights0,lWeights90,  &
        lWeights0Tmp,lWeights90Tmp
      real(kind=real64),dimension(:),allocatable::quadWeightsM,  &
        MSquaredList
      real(kind=real64),dimension(:,:),allocatable::quadGridM,moCoeffs
      logical::found,foundTheta90
      type(mqc_basisset)::basisSet
      type(MQC_Variable)::tmp,quadTmp
!
 1100 format(/,1x,'Lab frame ',i4,': epsilon = (',f7.3,',',f7.3,',',f7.3,')',  &
        '  reference k-plane vector = (',f7.3,',',f7.3,',',f7.3,')')
 1200 format(1x,'Photoelectron model flag = ',i3)
 1210 format(1x,'Photon energy          = ',f12.6,' eV')
 1220 format(1x,'Electron binding energy= ',f12.6,' eV')
 1230 format(1x,'Electron kinetic energy= ',f12.6,' eV = ',es14.6,' Eh')
 1240 format(1x,'Photoelectron k        = ',f12.6,' a.u.')
 1250 format(1x,'Lab-frame model flag   = ',i3)
 1260 format(1x,'Lab-frame orientations = ',i8,3x,'weight sum = ',es14.6)
 1270 format(1x,'Chi quadrature points  = ',i8,3x,'weight sum = ',es14.6)
  2000 format(1x,'Data for the ',A,' grid:',/,  &
        3x,'nPoints = ',i15,' step size = ',f20.5)
 2010 format(1x,'Data for the ',A,' grid:',/,3x,'nPoints = ',i15)
 2500 format(/,1x,'Dyson orbital self-overlap = ',f12.8)
 3000 format(/,1x,55('-'),/,  &
        12x,'Intensity as a Function of Theta',/,  &
        17x,'(k = ',f10.6,' a.u.)',/,  &
        9x,'theta (rad)',15x,'I(theta)',/,  &
        1x,55('='))
 3010 format(9x,f7.3,3x,es25.12)
 3020 format(1x,55('='))
 3100 format(1x,'I(0) = ',es15.8,3x,'I(90) = ',es15.8,3x,'beta = ',f12.6)
 3200 format(1x,'l weights at theta = ',A,':',*(1x,f9.5))
 3500 format(/,1x,87('-'),/,  &
        31x,'Summary of PAD Calculation',/,  &
        34x,'(k = ',f10.6,' a.u.)',/,  &
        1x,87('='))
  3510 format(1x,'Plane: ',A,'  |  epsilon: (',f5.2,',',f5.2,',',f5.2,  &
        ')',/,18x,'theta-integrated intensity       = ',es14.6,/,  &
        18x,'solid-angle integrated intensity = ',es14.6,/,  &
        18x,'beta(ratio) = ',f10.6,'  |  beta(fit) = ',f10.6,  &
        ' (R**2 = ',f8.5,')')
  3520 format(1x,87('='))
  3530 format(1x,'Weighted mean beta(ratio) = ',f10.6,/,  &
        1x,'Weighted mean beta(fit)   = ',f10.6)
 3540 format(1x,'Averaged theta-integrated intensity = ',es14.6,/,  &
        1x,'Averaged solid-angle integrated intensity = ',es14.6,/,  &
        1x,'Beta from averaged PAD(ratio) = ',f10.6,/,  &
        1x,'Beta from averaged PAD(fit)   = ',f10.6,  &
        ' (R**2 = ',f8.5,')',/)
 8998 format(1x,'Wall time for ',A,': ',f15.3,' s')
!
!     Validate user/model options and compute the photoelectron kinematics.
!
      if(options%dysonMOIndex.lt.1)  &
        call mqc_error('PAD: requested Dyson orbital index is out of range.')
      if(options%nGridPointsTheta.lt.2)  &
        call mqc_error('PAD: nGridPointsTheta must be at least 2.')
      if(options%nGridPointsM.lt.2)  &
        call mqc_error('PAD: nGridPointsM must be at least 2.')
      if(options%nChi.lt.1)  &
        call mqc_error('PAD: nChi must be positive.')
      if(options%iPEType.lt.0.or.options%iPEType.gt.2)  &
        call mqc_error('PAD: invalid photoelectron model flag.')
      if(options%lMax.lt.0)  &
        call mqc_error('PAD: lMax must be nonnegative.')
      if(options%nLabFrameTheta.lt.2)  &
        call mqc_error('PAD: nLabFrameTheta must be at least 2.')
      if(options%nLabFramePhi.lt.1)  &
        call mqc_error('PAD: nLabFramePhi must be positive.')
      if(options%labFrameType.eq.PAD_LAB_FRAMES_AXISYMMETRIC) then
        if(options%labFrameAlignment.lt.-mqc_float(1).or.  &
          options%labFrameAlignment.gt.mqc_float(2))  &
          call mqc_error('PAD: labFrameAlignment must be between -1 and 2.')
      endIf
      call padComputePhotoelectronEnergyAndK(options,  &
        results%photoelectronEnergyEV,results%photoelectronEnergyHartree,  &
        results%kMag)
      if(options%printResults) then
        write(iOut,1200) options%iPEType
        write(iOut,1210) options%photonEnergyEV
        write(iOut,1220) options%bindingEnergyEV
        write(iOut,1230) results%photoelectronEnergyEV,  &
          results%photoelectronEnergyHartree
        write(iOut,1240) results%kMag
      endIf
!
!     Set up the lab-frame polarizations. These are rotations of the lab frame
!     relative to a fixed molecular frame.
!
      call buildPADLabFrames(options,results%epsilonVector,  &
        results%kPlaneVector,results%intPlaneLabels,  &
        results%labFrameWeights)
      nIntPlanes = Size(results%epsilonVector,2)
      results%labFrameWeightSum = SUM(results%labFrameWeights)
      if(results%labFrameWeightSum.le.mqc_small)  &
        call mqc_error('PAD: lab-frame weights sum to zero.')
      if(options%printResults) then
        write(iOut,1250) options%labFrameType
        write(iOut,1260) nIntPlanes,results%labFrameWeightSum
      endIf
      results%nLabFrames = nIntPlanes
      results%nTheta = options%nGridPointsTheta
      results%nChi = options%nChi
      call buildChiQuadrature(results%chi,results%chiWeights,results%nChi)
      results%chiWeightSum = SUM(results%chiWeights)
      if(results%chiWeightSum.le.mqc_small)  &
        call mqc_error('PAD: chi quadrature weights sum to zero.')
      if(options%printResults) write(iOut,1270) results%nChi,  &
        results%chiWeightSum
      Allocate(results%thetaIntegratedIntensity(nIntPlanes),  &
        results%solidAngleIntegratedIntensity(nIntPlanes),  &
        results%intensity0(nIntPlanes),results%intensity90(nIntPlanes),  &
        results%betaValsParaPerp(nIntPlanes),  &
        results%betaValsFit(nIntPlanes),results%rSquared(nIntPlanes),  &
        results%intensityTheta(options%nGridPointsTheta,nIntPlanes),  &
        results%averageIntensityTheta(options%nGridPointsTheta),  &
        results%theta(options%nGridPointsTheta),  &
        results%thetaWeights(options%nGridPointsTheta),  &
        results%thetaSolidAngleWeights(options%nGridPointsTheta))
      results%thetaIntegratedIntensity = mqc_float(0)
      results%solidAngleIntegratedIntensity = mqc_float(0)
      results%intensity0 = mqc_float(0)
      results%intensity90 = mqc_float(0)
      results%betaValsParaPerp = mqc_float(0)
      results%betaValsFit = mqc_float(0)
      results%rSquared = mqc_float(0)
      results%intensityTheta = mqc_float(0)
      results%averageIntensityTheta = mqc_float(0)
!
!     Read the basis set and MO coefficients from the FAF.
!
      call loadGaussianBasisSet(faf,basisSet)
      call faf%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmp)
      moCoeffs = tmp
      if(options%dysonMOIndex.gt.Size(moCoeffs,2))  &
        call mqc_error('PAD: requested Dyson orbital index is out of range.')
!
!     Use the Gaussian XC quadrature grid from the FAF when available. Otherwise
!     fall back to a simple Cartesian trapezoid grid.
!
      tStart1 = omp_get_wtime()
      call faf%getArray('3D Quadrature Grid',mqcVarOut=quadTmp,foundOut=found)
      if(found) then
        Allocate(quadWeightsM(SIZE(quadTmp,2)),quadGridM(3,SIZE(quadTmp,2)))
        do i = 1,SIZE(quadTmp,2)
          quadWeightsM(i) = quadTmp%getVal([ 1,i ])
          quadGridM(1,i) = quadTmp%getVal([ 2,i ])
          quadGridM(2,i) = quadTmp%getVal([ 3,i ])
          quadGridM(3,i) = quadTmp%getVal([ 4,i ])
        endDo
        if(options%printResults) write(iOut,2010)  &
          'FAF M quadrature',SIZE(quadWeightsM)
      else
        cartStart = [ -mqc_float(6),-mqc_float(6),-mqc_float(6) ]
        cartEnd = [ mqc_float(6),mqc_float(6),mqc_float(6) ]
        stepSizeIntM = (cartEnd(1)-cartStart(1))/  &
          mqc_float(options%nGridPointsM-1)
        Allocate(quadGridM(3,options%nGridPointsM**3),  &
          quadWeightsM(options%nGridPointsM**3))
        call setup_quadrature_trapezoid3d(options%nGridPointsM,  &
          stepSizeIntM,cartStart,quadGridM,quadWeightsM)
        if(options%printResults) write(iOut,2000) 'fallback M quadrature',  &
          options%nGridPointsM**3,stepSizeIntM
      endIf
      tEnd1 = omp_get_wtime()
      if(options%printResults) write(iOut,8998)  &
        'M quadrature setup',tEnd1-tStart1
      if(options%printResults) flush(iOut)
!
!     Fill the theta grid for the PAD scan.
!
      thetaStart = mqc_float(0)
      stepSizeTheta = Pi/mqc_float(options%nGridPointsTheta-1)
      call setup_quadrature_trapezoid1d(options%nGridPointsTheta,  &
        stepSizeTheta,thetaStart,results%theta,results%thetaWeights)
      results%thetaSolidAngleWeights = sin(results%theta)*results%thetaWeights
      if(options%printResults) write(iOut,2000)  &
        'theta',options%nGridPointsTheta,stepSizeTheta
      if(options%printResults) flush(iOut)
      foundTheta90 = .false.
      indexTheta90 = 0
      thetaMatchTol = sqrt(mqc_small)
      do j = 1,options%nGridPointsTheta
        if(abs(results%theta(j)-Pi/mqc_float(2)).lt.thetaMatchTol) then
          foundTheta90 = .true.
          indexTheta90 = j
        endIf
      endDo
!
!     Test that the chosen MO is normalized on the quadrature grid.
!
      tStart1 = omp_get_wtime()
      results%dysonSelfOverlap = moInnerProductNumericalIntegration(  &
        moCoeffs(:,options%dysonMOIndex),quadGridM,quadWeightsM,basisSet)
      if(options%printResults) write(iOut,2500) results%dysonSelfOverlap
      tEnd1 = omp_get_wtime()
      if(options%printResults) write(iOut,8998) 'norm test',tEnd1-tStart1
      if(options%printResults) flush(iOut)
!
!     Evaluate the PAD for each lab-frame orientation.
!
      Allocate(lWeights0(0:options%lMax),lWeights90(0:options%lMax),  &
        lWeights0Tmp(0:options%lMax),lWeights90Tmp(0:options%lMax),  &
        MSquaredList(options%nGridPointsTheta))
      do i = 1,nIntPlanes
        call buildTransverseBasis(results%epsilonVector(:,i),  &
          results%kPlaneVector(:,i),epsilonUnit,uBasis,vBasis)
        results%epsilonVector(:,i) = epsilonUnit
        if(options%printResults) write(iOut,1100) i,  &
          results%epsilonVector(:,i),results%kPlaneVector(:,i)
        results%intensityTheta(:,i) = mqc_float(0)
        results%intensity0(i) = mqc_float(0)
        results%intensity90(i) = mqc_float(0)
        lWeights0 = mqc_float(0)
        lWeights90 = mqc_float(0)
        timeSingleAngles = mqc_float(0)
        timeThetaList = mqc_float(0)
        tStart1 = omp_get_wtime()
!
        do iChi = 1,results%nChi
          epsilonVec = results%epsilonVector(:,i)
          uChi = cos(results%chi(iChi))*uBasis+sin(results%chi(iChi))*vBasis
          call mqc_normalizeVector(uChi)
          if(options%iPEType.eq.0) then
            tStart2 = omp_get_wtime()
            MSquaredList = dysonPlaneWaveMatrixElementSquaredThetaList(  &
              results%theta,results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM)
            tEnd2 = omp_get_wtime()
            timeThetaList = timeThetaList+tEnd2-tStart2
            MSquared0 = MSquaredList(1)
            if(foundTheta90) then
              MSquared90 = MSquaredList(indexTheta90)
            else
              tStart2 = omp_get_wtime()
              MSquared90 = dysonPlaneWaveMatrixElementSquared(  &
                Pi/mqc_float(2),results%kMag,epsilonVec,uChi,  &
                moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
                quadWeightsM)
              tEnd2 = omp_get_wtime()
              timeSingleAngles = timeSingleAngles+tEnd2-tStart2
            endIf
          else
            tStart2 = omp_get_wtime()
            call dysonMatrixElement1Angle(options%iPEType,options%lMax,  &
              mqc_float(0),results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM,MSquared0,lWeights0Tmp)
            call dysonMatrixElement1Angle(options%iPEType,options%lMax,  &
              Pi/mqc_float(2),results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM,MSquared90,lWeights90Tmp)
            tEnd2 = omp_get_wtime()
            timeSingleAngles = timeSingleAngles+tEnd2-tStart2
            tStart2 = omp_get_wtime()
            call dysonMatrixElementThetaList(options%iPEType,options%lMax,  &
              results%theta,results%kMag,epsilonVec,uChi,  &
              moCoeffs(:,options%dysonMOIndex),basisSet,quadGridM,  &
              quadWeightsM,MSquaredList)
            tEnd2 = omp_get_wtime()
            timeThetaList = timeThetaList+tEnd2-tStart2
            if(options%iPEType.eq.2) then
              lWeights0 = lWeights0+results%chiWeights(iChi)*lWeights0Tmp
              lWeights90 = lWeights90+results%chiWeights(iChi)*lWeights90Tmp
            endIf
          endIf
          results%intensityTheta(:,i) = results%intensityTheta(:,i)+  &
            results%chiWeights(iChi)*MSquaredList
          results%intensity0(i) = results%intensity0(i)+  &
            results%chiWeights(iChi)*MSquared0
          results%intensity90(i) = results%intensity90(i)+  &
            results%chiWeights(iChi)*MSquared90
        endDo
        results%intensityTheta(:,i) = results%intensityTheta(:,i)/  &
          results%chiWeightSum
        results%intensity0(i) = results%intensity0(i)/results%chiWeightSum
        results%intensity90(i) = results%intensity90(i)/results%chiWeightSum
        if(options%iPEType.eq.2) then
          lWeights0 = lWeights0/results%chiWeightSum
          lWeights90 = lWeights90/results%chiWeightSum
        endIf
!
        tEnd1 = omp_get_wtime()
        if(options%printResults) then
          write(iOut,8998) 'single-angle MSquared calls',timeSingleAngles
          write(iOut,8998) 'theta-list MSquared calls',timeThetaList
          if(results%nChi.gt.1) then
            write(iOut,8998) 'chi-averaged MSquared(theta) list',tEnd1-tStart1
          else
            write(iOut,8998) 'MSquared(theta) list',tEnd1-tStart1
          endIf
        endIf
        if(options%iPEType.eq.2.and.options%printResults) then
          write(iOut,3200) '0 ',lWeights0
          write(iOut,3200) '90',lWeights90
        endIf
!
!       Compute beta using both the parallel/perpendicular ratio and the
!       least-squares fit to the PAD shape.
!
        call betaLeastSquares(results%intensityTheta(:,i),results%theta,  &
          results%betaValsFit(i),results%rSquared(i))
        results%thetaIntegratedIntensity(i) =  &
          dot_product(results%intensityTheta(:,i),results%thetaWeights)
        results%solidAngleIntegratedIntensity(i) = results%chiWeightSum*  &
          dot_product(results%intensityTheta(:,i),  &
          results%thetaSolidAngleWeights)
        results%betaValsParaPerp(i) = betaParaPerp(results%intensity0(i),  &
          results%intensity90(i))
!
        if(options%printResults.and.options%printThetaTable) then
          write(iOut,3000) results%kMag
          do j = 1,options%nGridPointsTheta
            write(iOut,3010) results%theta(j),results%intensityTheta(j,i)
          endDo
          write(iOut,3020)
          write(iOut,3100) results%intensity0(i),results%intensity90(i),  &
            results%betaValsParaPerp(i)
          flush(iOut)
        endIf
      endDo
!
!     Build the weighted orientation-averaged PAD and associated beta values.
!
      do i = 1,nIntPlanes
        results%averageIntensityTheta = results%averageIntensityTheta +  &
          results%labFrameWeights(i)*results%intensityTheta(:,i)
      endDo
      results%averageIntensityTheta = results%averageIntensityTheta/  &
        results%labFrameWeightSum
      results%averageIntensity0 = dot_product(results%labFrameWeights,  &
        results%intensity0)/results%labFrameWeightSum
      results%averageIntensity90 = dot_product(results%labFrameWeights,  &
        results%intensity90)/results%labFrameWeightSum
      results%averageThetaIntegratedIntensity =  &
        dot_product(results%averageIntensityTheta,results%thetaWeights)
      results%averageSolidAngleIntegratedIntensity = results%chiWeightSum*  &
        dot_product(results%averageIntensityTheta,  &
        results%thetaSolidAngleWeights)
      results%averageBetaParaPerp = betaParaPerp(results%averageIntensity0,  &
        results%averageIntensity90)
      call betaLeastSquares(results%averageIntensityTheta,results%theta,  &
        results%averageBetaFit,results%averageRSquared)
      results%meanBetaParaPerp = dot_product(results%labFrameWeights,  &
        results%betaValsParaPerp)/results%labFrameWeightSum
      results%meanBetaFit = dot_product(results%labFrameWeights,  &
        results%betaValsFit)/results%labFrameWeightSum
!
!     Print the summary.
!
      if(options%printResults) then
        write(iOut,3500) results%kMag
        do i = 1,nIntPlanes
          write(iOut,3510) TRIM(results%intPlaneLabels(i)),  &
            results%epsilonVector(:,i),results%thetaIntegratedIntensity(i),  &
            results%solidAngleIntegratedIntensity(i),  &
            results%betaValsParaPerp(i),results%betaValsFit(i),  &
            results%rSquared(i)
        endDo
        write(iOut,3520)
        write(iOut,3530) results%meanBetaParaPerp,results%meanBetaFit
        write(iOut,3540) results%averageThetaIntegratedIntensity,  &
          results%averageSolidAngleIntegratedIntensity,  &
          results%averageBetaParaPerp,results%averageBetaFit,  &
          results%averageRSquared
      endIf
!
      return
      end subroutine runPADCalculation


!PROCEDURE buildPADLabFrames
      subroutine buildPADLabFrames(options,epsilonVectors,kPlaneVectors,  &
        planeLabels,labFrameWeights)
!
!     This routine builds the lab-frame polarization and reference k-plane
!     vector arrays requested in the PAD options object. The calculation driver
!     only loops over these arrays and does not need to know how they were
!     generated.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      type(pad_options),intent(in)::options
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        epsilonVectors,kPlaneVectors
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        labFrameWeights
      character(len=8),dimension(:),allocatable,intent(out)::planeLabels
!
      integer(kind=int64)::i,nLabFrames
!
!     Build the requested lab-frame vector set.
!
      select case(options%labFrameType)
      case(PAD_LAB_FRAMES_CARTESIAN)
        call buildCartesianLabFrames(epsilonVectors,kPlaneVectors,  &
          planeLabels,labFrameWeights)
!
      case(PAD_LAB_FRAMES_SPHERE)
        call buildSphereLabFrames(epsilonVectors,kPlaneVectors,  &
          planeLabels,labFrameWeights,options%nLabFrameTheta,  &
          options%nLabFramePhi)
!
      case(PAD_LAB_FRAMES_AXISYMMETRIC)
        call buildAxisymmetricLabFrames(epsilonVectors,kPlaneVectors,  &
          planeLabels,labFrameWeights,options%nLabFrameTheta,  &
          options%nLabFramePhi,options%labFrameAlignment)
!
      case(PAD_LAB_FRAMES_CUSTOM)
        if(.not.Allocated(options%labEpsilonVector))  &
          call mqc_error('buildPADLabFrames: custom epsilon vectors missing.')
        if(.not.Allocated(options%labKPlaneVector))  &
          call mqc_error('buildPADLabFrames: custom k-plane vectors missing.')
        if(Size(options%labEpsilonVector,1).ne.3)  &
          call mqc_error('buildPADLabFrames: epsilon vectors must be 3 x N.')
        if(Size(options%labKPlaneVector,1).ne.3)  &
          call mqc_error('buildPADLabFrames: k-plane vectors must be 3 x N.')
        if(Size(options%labEpsilonVector,2).ne.  &
          Size(options%labKPlaneVector,2))  &
          call mqc_error('buildPADLabFrames: custom vector counts differ.')
        nLabFrames = Size(options%labEpsilonVector,2)
        if(nLabFrames.lt.1)  &
          call mqc_error('buildPADLabFrames: no custom lab frames supplied.')
        epsilonVectors = options%labEpsilonVector
        kPlaneVectors = options%labKPlaneVector
        if(Allocated(options%labFrameLabels)) then
          if(Size(options%labFrameLabels).ne.nLabFrames)  &
            call mqc_error('buildPADLabFrames: custom label count differs.')
          planeLabels = options%labFrameLabels
        else
          Allocate(planeLabels(nLabFrames))
          do i = 1,nLabFrames
            write(planeLabels(i),'("c",i7.7)') i
          endDo
        endIf
        if(Allocated(options%labFrameWeights)) then
          if(Size(options%labFrameWeights).ne.nLabFrames)  &
            call mqc_error('buildPADLabFrames: custom weight count differs.')
          if(MINVAL(options%labFrameWeights).lt.mqc_float(0))  &
            call mqc_error('buildPADLabFrames: negative custom weight.')
          labFrameWeights = options%labFrameWeights
        else
          Allocate(labFrameWeights(nLabFrames))
          labFrameWeights = mqc_float(1)
        endIf
!
      case default
        call mqc_error('buildPADLabFrames: invalid lab-frame type.')
      end select
!
      return
      end subroutine buildPADLabFrames


!PROCEDURE buildCartesianLabFrames
      subroutine buildCartesianLabFrames(epsilonVectors,kPlaneVectors,  &
        planeLabels,labFrameWeights)
!
!     Build the default fixed-orientation lab-frame set. The molecule stays in
!     the input molecular frame; each column gives a lab electric-field
!     polarization vector and one perpendicular reference vector used to define
!     the chi-rotated k-vector scan plane.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        epsilonVectors,kPlaneVectors
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        labFrameWeights
      character(len=8),dimension(:),allocatable,intent(out)::planeLabels
!
      Allocate(epsilonVectors(3,3),kPlaneVectors(3,3),planeLabels(3),  &
        labFrameWeights(3))
      labFrameWeights = mqc_float(1)
!
      planeLabels(1) = 'xy'
      epsilonVectors(:,1) = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
      kPlaneVectors(:,1)  = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
!
      planeLabels(2) = 'yz'
      epsilonVectors(:,2) = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
      kPlaneVectors(:,2)  = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
!
      planeLabels(3) = 'zx'
      epsilonVectors(:,3) = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
      kPlaneVectors(:,3)  = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
!
      return
      end subroutine buildCartesianLabFrames


!PROCEDURE buildSphereLabFrames
      subroutine buildSphereLabFrames(epsilonVectors,kPlaneVectors,  &
        planeLabels,labFrameWeights,nLabFrameTheta,nLabFramePhi)
!
!     This routine builds a lab-frame vector set from points on the unit sphere.
!     The sphere points are used as electric-field polarization vectors, and the
!     associated orthogonal vectors define reference directions for the
!     chi-rotated k-vector scan planes.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        epsilonVectors,kPlaneVectors
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        labFrameWeights
      character(len=8),dimension(:),allocatable,intent(out)::planeLabels
      integer(kind=int64),intent(in)::nLabFrameTheta,nLabFramePhi
!
      integer(kind=int64)::i,nLabFrames
!
!     Build the vectors and labels.
!
      if(nLabFrameTheta.lt.3)  &
        call mqc_error('buildSphereLabFrames: theta grid must be at least 3.')
      call buildSphereGrid(epsilonVectors,kPlaneVectors,nLabFrameTheta,  &
        nLabFramePhi,weights=labFrameWeights)
      nLabFrames = Size(epsilonVectors,2)
      Allocate(planeLabels(nLabFrames))
      do i = 1,nLabFrames
        write(planeLabels(i),'("s",i7.7)') i
      endDo
!
      return
      end subroutine buildSphereLabFrames


!PROCEDURE buildAxisymmetricLabFrames
      subroutine buildAxisymmetricLabFrames(epsilonVectors,kPlaneVectors,  &
        planeLabels,labFrameWeights,nLabFrameTheta,nLabFramePhi,  &
        alignment)
!
!     This routine builds an axisymmetric lab-frame distribution around the
!     lab z axis. It starts from the sphere-grid orientations and applies a
!     nonnegative weight factor 1 + alignment*P2(cos(theta)).
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(:,:),allocatable,intent(out)::  &
        epsilonVectors,kPlaneVectors
      real(kind=real64),dimension(:),allocatable,intent(out)::  &
        labFrameWeights
      character(len=8),dimension(:),allocatable,intent(out)::planeLabels
      integer(kind=int64),intent(in)::nLabFrameTheta,nLabFramePhi
      real(kind=real64),intent(in)::alignment
!
      integer(kind=int64)::i,nLabFrames
      real(kind=real64)::p2Val,weightFactor,weightSum
!
!     Build the base sphere-grid vectors and reweight them by the axisymmetric
!     distribution. The accepted alignment range keeps all weights
!     nonnegative for -1/2 <= P2(cos(theta)) <= 1.
!
      if(alignment.lt.-mqc_float(1).or.alignment.gt.mqc_float(2))  &
        call mqc_error('buildAxisymmetricLabFrames: alignment out of range.')
      call buildSphereLabFrames(epsilonVectors,kPlaneVectors,planeLabels,  &
        labFrameWeights,nLabFrameTheta,nLabFramePhi)
      nLabFrames = Size(epsilonVectors,2)
      do i = 1,nLabFrames
        p2Val = (mqc_float(3)*epsilonVectors(3,i)*epsilonVectors(3,i)-  &
          mqc_float(1))/mqc_float(2)
        weightFactor = mqc_float(1)+alignment*p2Val
        if(weightFactor.lt.-mqc_small)  &
          call mqc_error('buildAxisymmetricLabFrames: negative weight factor.')
        labFrameWeights(i) = labFrameWeights(i)*max(mqc_float(0),weightFactor)
        write(planeLabels(i),'("a",i7.7)') i
      endDo
      weightSum = SUM(labFrameWeights)
      if(weightSum.le.mqc_small)  &
        call mqc_error('buildAxisymmetricLabFrames: weights sum to zero.')
      labFrameWeights = labFrameWeights*(mqc_float(4)*Pi/weightSum)
!
      return
      end subroutine buildAxisymmetricLabFrames


!PROCEDURE buildChiQuadrature
      subroutine buildChiQuadrature(chiVals,chiWeights,nChi)
!
!     This routine builds a periodic uniform quadrature over chi from 0 to
!     2*pi. The point at 2*pi is omitted because it is equivalent to the point
!     at 0 for periodic sampling.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      integer(kind=int64),intent(in)::nChi
      real(kind=real64),dimension(:),allocatable,intent(out)::chiVals,  &
        chiWeights
!
      integer(kind=int64)::i
      real(kind=real64)::chiStep
!
      if(nChi.lt.1) call mqc_error('buildChiQuadrature: nChi must be positive.')
      chiStep = mqc_float(2)*Pi/mqc_float(nChi)
      Allocate(chiVals(nChi),chiWeights(nChi))
      do i = 1,nChi
        chiVals(i) = mqc_float(i-1)*chiStep
        chiWeights(i) = chiStep
      endDo
!
      return
      end subroutine buildChiQuadrature


!PROCEDURE buildTransverseBasis
      subroutine buildTransverseBasis(epsilonVector,kPlaneReference,  &
        epsilonUnit,uBasis,vBasis)
!
!     This routine normalizes the input polarization vector and constructs an
!     orthonormal transverse basis around it. The input reference k-plane vector
!     is projected into the transverse space when possible and is otherwise
!     replaced by a simple Cartesian fallback direction.
!
!
!     H. P. Hratchian, 2026.
!
      implicit none
      real(kind=real64),dimension(3),intent(in)::epsilonVector,  &
        kPlaneReference
      real(kind=real64),dimension(3),intent(out)::epsilonUnit,uBasis,vBasis
!
      real(kind=real64),dimension(3)::fallbackVector
!
      epsilonUnit = epsilonVector
      if(dot_product(epsilonUnit,epsilonUnit).le.mqc_small)  &
        call mqc_error('buildTransverseBasis: epsilon vector has zero norm.')
      call mqc_normalizeVector(epsilonUnit)
!
!     Start from the supplied reference vector and project out any component
!     parallel to epsilon.
!
      uBasis = kPlaneReference-dot_product(kPlaneReference,epsilonUnit)*  &
        epsilonUnit
      if(dot_product(uBasis,uBasis).le.mqc_small) then
        if(abs(epsilonUnit(1)).le.abs(epsilonUnit(2)).and.  &
          abs(epsilonUnit(1)).le.abs(epsilonUnit(3))) then
          fallbackVector = [ mqc_float(1),mqc_float(0),mqc_float(0) ]
        elseIf(abs(epsilonUnit(2)).le.abs(epsilonUnit(3))) then
          fallbackVector = [ mqc_float(0),mqc_float(1),mqc_float(0) ]
        else
          fallbackVector = [ mqc_float(0),mqc_float(0),mqc_float(1) ]
        endIf
        uBasis = fallbackVector-dot_product(fallbackVector,epsilonUnit)*  &
          epsilonUnit
      endIf
      if(dot_product(uBasis,uBasis).le.mqc_small)  &
        call mqc_error('buildTransverseBasis: could not form transverse basis.')
      call mqc_normalizeVector(uBasis)
!
!     Complete the right-handed transverse basis and re-orthogonalize uBasis.
!
      vBasis = mqc_crossProduct3D_real(epsilonUnit,uBasis)
      if(dot_product(vBasis,vBasis).le.mqc_small)  &
        call mqc_error('buildTransverseBasis: degenerate transverse basis.')
      call mqc_normalizeVector(vBasis)
      uBasis = mqc_crossProduct3D_real(vBasis,epsilonUnit)
      call mqc_normalizeVector(uBasis)
!
      return
      end subroutine buildTransverseBasis



!PROCEDURE betaParaPerp
      function betaParaPerp(IPara,IPerp) result(beta)
!
!     This function evaluates the photoelectron anisotropy parameter, beta,
!     using the analytic form based on the angular intensities when the photon
!     electric field is parallel and perpendicular to the detected detachment
!     direction.
!
!
!     H. P. Hratchian, 2025, 2026.
!
      implicit none
      real(kind=real64),intent(in)::IPara,IPerp
      real(kind=real64)::beta,denominator
!
!     Do the work...
!
      denominator = (IPara/mqc_float(2))+IPerp
      if(abs(denominator).le.tiny(denominator)) then
        beta = mqc_float(0)
      else
        beta = (IPara-IPerp)/denominator
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
!     H. P. Hratchian, 2025, 2026.
!
!
      implicit none
      real(kind=real64),dimension(:),intent(in)::MSquaredList,thetaList
      real(kind=real64)::beta,rSquared
!
      integer(kind=int64)::n
      real(kind=real64)::slope,intercept,signalMax,signalRange
      real(kind=real64),dimension(:),allocatable::x,y
!
!     Before doing any work, check to see if all intensities are truly zero. If
!     they are, then beta is zero. The tests are relative to the PAD signal
!     scale so near-threshold physical intensities are not discarded.
!
      beta = mqc_float(0)
      rSquared = mqc_float(0)
      signalMax = maxval(abs(MSquaredList))
      if(signalMax.le.tiny(signalMax)) return
      signalRange = maxval(MSquaredList)-minval(MSquaredList)
      if(signalRange.le.epsilon(signalRange)*signalMax) then
        rSquared = mqc_float(1)
        return
      endIf
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
      if(abs(intercept).le.tiny(intercept)) then
        beta = mqc_float(0)
      else
        beta = slope/intercept
      endIf
!
      DeAllocate(x)
      return
      end subroutine betaLeastSquares
!
!PROCEDURE buildSphereGrid
      subroutine buildSphereGrid(xyz,orthogVector,nThetaInput,nPhiInput,  &
        weights)
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
      integer(kind=int64),intent(in),optional::nThetaInput,nPhiInput
      real(kind=real64),dimension(:),allocatable,intent(out),optional::weights
      integer(kind=int64)::i,j,nThetaPoints,nPhiPoints,ixyz
      real(kind=real64)::thetaStep,phiStep,theta,phi,x,y,z,  &
        thetaLower,thetaUpper,pointWeight
!
 1000 format(1x,i4,': theta,phi (pi)=',f5.2,',',f5.2,' | x,y,z=',f8.4,',',f8.4,',',f8.4)
!
!     Start by setting nThetaPoints and  nPhiPoints, and figuring out thetaStep
!     and phiStep. Then, allocate xyz.
!
      nThetaPoints = 5
      nPhiPoints = 8
      if(PRESENT(nThetaInput)) nThetaPoints = nThetaInput
      if(PRESENT(nPhiInput)) nPhiPoints = nPhiInput
      if(nThetaPoints.lt.2)  &
        call mqc_error('buildSphereGrid: nThetaPoints must be at least 2.')
      if(nPhiPoints.lt.1)  &
        call mqc_error('buildSphereGrid: nPhiPoints must be positive.')
      thetaStep = Pi/mqc_float(nThetaPoints-1)
      phiStep = mqc_float(2)*Pi/mqc_float(nPhiPoints)
!
!     Loop over theta and phi to build the array xyz.
!
      Allocate(xyz(3,(nThetaPoints-2)*nPhiPoints+2),  &
        orthogVector(3,(nThetaPoints-2)*nPhiPoints+2))
      if(PRESENT(weights))  &
        Allocate(weights((nThetaPoints-2)*nPhiPoints+2))
      theta = 0
      ixyz = 1
      do i = 1,nThetaPoints
        theta = mqc_float(i-1)*thetaStep
        thetaLower = max(mqc_float(0),theta-(thetaStep/mqc_float(2)))
        thetaUpper = min(Pi,theta+(thetaStep/mqc_float(2)))
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
          if(PRESENT(weights)) then
            pointWeight = mqc_float(2)*Pi*(cos(thetaLower)-cos(thetaUpper))
            weights(ixyz) = pointWeight
          endIf
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
          if(PRESENT(weights)) then
            pointWeight = mqc_float(2)*Pi*(cos(thetaLower)-cos(thetaUpper))
            weights(ixyz) = pointWeight
          endIf
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
            if(PRESENT(weights)) then
              pointWeight = phiStep*(cos(thetaLower)-cos(thetaUpper))
              weights(ixyz) = pointWeight
            endIf
            ixyz = ixyz+1
          endDo
        endIf
      endDo
!
      return
      end subroutine buildSphereGrid

      end module pad_mod
