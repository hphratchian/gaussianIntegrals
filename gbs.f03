      program gbs
!
!     This program tests the MQC_Integrals module by using a local version named
!     mqc_integrals1.F03. To run the program, an FAF is given as the only
!     command line argument. If the FAF is named myFile.faf, then this program
!     is run using:
!           ./gbs.exe myFile.faf
!
!     The program runs three sets of tests.
!
!     First, the program loads values for a couple contracted gaussian-type
!     function (CGTF) into mqc_cgtf objects.
!
!     Second, the program reads the basis set information from a Gaussian FAF.
!     Info about the basis set is written to the output and then we test the
!     MQC_Integrals overlap code to form the overlap matrix. 
!
!     The third section carries out numerical integration of the density of the
!     first five MOs.
!
!
!     Hrant P. Hratchian, 2024.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      USE OMP_LIB
      USE gbs_mod
      implicit none
      integer,parameter::lVal=0
      integer::i,j,k,l
      logical::fail=.false.,atEnd
      real(kind=real64)::moVal1,moVal2
      real(kind=real64),dimension(5,5)::tmpSij=float(0)
      character(len=256)::fafName
      logical,dimension(5,5)::haveSij=.false.
      type(mqc_cgtf)::bf1,bf2,bf3,bf4
      type(mqc_gaussian_unformatted_matrix_file)::faf
      type(mqc_basisset)::basisSet
!
      integer(kind=int64)::ixyz,jxyz
      real(kind=real64)::mu,p
      real(kind=real64),dimension(3)::xAB,xPA,xPB
!
      integer(kind=int64)::nCGTF=0,iCGTF=0,nBasis,nGridPoints=401,nOMP
      real(kind=real64)::tStart,tEnd
      real(kind=real64),dimension(:),allocatable::basisValues
      real(kind=real64),dimension(:,:),allocatable::basisIntegrals,  &
        overlapMatrix
      real(kind=real64),dimension(:),allocatable::quadWeights,quadValues
      real(kind=real64),dimension(:,:),allocatable::quadGrid,moCoeffs
      type(MQC_Variable)::tmp
!
!     Format statements.
!
 1000 format(1x,'Program gbs.')
!
!
!     Begin the program.
!
      call CPU_TIME(tStart)
      nOMP=48
      call omp_set_num_threads(nOMP)
      fail = .false.
      write(iOut,1000)
!$OMP PARALLEL
PRINT *, "Thread number:", OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
!
!     Read the FAF name from the command line and load the object faf.
!
      call get_command_argument(1,fafName)
      call faf%load(fafName)
!
!     Fill bf1-bf4 with info and then print them.
!
      call MQC_CGTF_init(bf1,lVal,[ 0.0,0.0,0.370424/angPBohr ],  &
        [ 1.0 ],[ 1.309756377 ])
      call MQC_CGTF_init(bf2,0,[ 0.0,0.0,0.370424/angPBohr ],  &
        [ 1.0 ],[ 0.2331359749 ])
      call MQC_CGTF_init(bf3,0,[ 0.0,0.0,-0.370424/angPBohr ],  &
        [ 1.0 ],[ 1.309756377 ])
      call MQC_CGTF_init(bf4,0,[ 0.0,0.0,-0.370424/angPBohr ],  &
        [ 1.0 ],[ 0.2331359749 ])
      call MQC_CGTF_print(bf1,iOut)
      call MQC_CGTF_print(bf2,iOut)
      call MQC_CGTF_print(bf3,iOut)
      call MQC_CGTF_print(bf4,iOut)
      call mqc_print(MQC_CGTF_extendLArray(bf1),iOut,header='bf1 lArray:',  &
        blank_at_top=.true.,blank_at_bottom=.true.)
      call mqc_print(MQC_CGTF_extendLArray(bf2),iOut,header='bf2 lArray:',  &
        blank_at_top=.true.,blank_at_bottom=.true.)
!
!     Try out the shell2nBasis function.
!
      write(iOut,*)
      write(IOut,*)' For basis 1, nBasis = ',bf1%shell2nBasis()
      write(IOut,*)' For basis 2, nBasis = ',bf2%shell2nBasis()
      write(IOut,*)' For basis 3, nBasis = ',bf3%shell2nBasis()
      write(IOut,*)' For basis 4, nBasis = ',bf4%shell2nBasis()
!
!     Try the self-overlap function.
!
      write(iOut,*)
      write(iOut,*)' For basis 1, self-S  = ',bf1%primitiveSelfOverlap(1,[lVal,0,0])
      write(iOut,*)' For basis 2, self-S  = ',bf2%primitiveSelfOverlap(1,[0,0,0])
      write(iOut,*)' For basis 3, self-S  = ',bf3%primitiveSelfOverlap(1,[0,0,0])
      write(iOut,*)' For basis 4, self-S  = ',bf4%primitiveSelfOverlap(1,[0,0,0])
!
!     Test the primitive overlap integral routine.
!
!     (1|1) ...
      ixyz = lVal
      jxyz = lVal
      mu = bf1%alpha(1)*bf1%alpha(1)/(bf1%alpha(1)+bf1%alpha(1)) 
      p = bf1%alpha(1)+bf1%alpha(1)
      xAB = bf1%center(1)-bf1%center(1)
      xPA = -xAB*bf1%alpha(1)/p
      xPB =  xAB*bf1%alpha(1)/p
      write(iOut,*)
      write(iOut,*)' Intermediates...1:'
      write(iOut,*)'    mu     = ',mu
      write(iOut,*)'    p      = ',p
      write(iOut,*)'    xAB(1) = ',xAB(1)
      write(iOut,*)'    xPA(1) = ',xPA(1)
      write(iOut,*)'    xPB(1) = ',xPB(1)
      call MQC_Overlap_Distribution_Primitive_XYZ_Constants(bf1,bf1, &
        1,1,mu,p,xAB,xPA,xPB)
      write(iOut,*)
      write(iOut,*)' Intermediates...2:'
      write(iOut,*)'    mu     = ',mu
      write(iOut,*)'    p      = ',p
      write(iOut,*)'    xAB(1) = ',xAB(1)
      write(iOut,*)'    xPA(1) = ',xPA(1)
      write(iOut,*)'    xPB(1) = ',xPB(1)
      call MQC_Overlap_Primitive_XYZ_OS(ixyz,jxyz,mu,p,xAB(1),xPA(1),  &
        xPB(1),tmpSij,haveSij)
      call mqc_print(haveSij,iOut,header='The HaveSij Matrix:')
      call mqc_print(tmpSij,iOut,header='The Sij Matrix:')
!
!     Test routine MQC_Overlap_CGFT.
!
      write(IOut,*)
      write(iOut,*)' ----------------------------------------------'
      write(iOut,*)' Running overlap of 1 with 1.'
      call MQC_Overlap_CGFT(bf1,bf1,basisIntegrals)
      write(IOut,*)
      write(iOut,*)' ----------------------------------------------'
      write(iOut,*)' Running overlap of 1 with 3.'
      call MQC_Overlap_CGFT(bf1,bf3,basisIntegrals)
!
!     Try reading a basis set from faf.
!
      call loadGaussianBasisSet(faf,basisSet)
      call basisSet%print(iOut)
      overlapMatrix = basisSetOverlapMatrix(basisSet)
      call mqc_print(overlapMatrix,iOut,header='overlap matrix')
      write(iOut,*)
      write(iOut,*)
      call MQC_Value_CGFT(basisSet%shells(2),[ 0.1,0.1,0.7 ],basisValues)
      call mqc_print(basisValues,iOut,header='basisValues on shell 2')
      write(iOut,*)
      write(iOut,*)
      basisValues = basisSetValuesList(basisSet,[ 0.1,0.1,0.7 ])
      call mqc_print(basisValues,iOut,header='basisValues over all basis functions')
      write(iOut,*)
      write(iOut,*)
      Allocate(quadGrid(3,nGridPoints**3),quadWeights(nGridPoints**3),  &
        quadValues(nGridPoints**3))
      call setup_quadrature_trapezoid3d(nGridPoints,  &
        0.025,  &
        [ -5.0,-5.0,-5.0 ],  &
        quadGrid,quadWeights)
      write(iOut,*)' max grid point: ',maxval(quadGrid)
!
!     The next block of code computes the integrals of the first 5 (alpha) MOs
!     with themselves. If the numerical integration is working well, the answer
!     should be very close to 1.
!
!     For a slightly different test, the value of moVal2 below can be changed to
!
!          moVal2 = dot_product(moCoeffs(:,1),basisValues)
!
!     In such a case, this should show that MO 1 is orthogonal to MOs 2, 3, 4,
!     and 5. It should also should, again, that MO1 is normalized.
!
      call faf%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmp)
      moCoeffs = tmp
      do j = 1,5
      !$omp parallel do default(none) &
      !$omp shared(j, moCoeffs, quadWeights, quadGrid, basisSet, quadValues) &
      !$omp private(i, basisValues, moVal1, moVal2)
        do i = 1, size(quadWeights)
          basisValues = basisSetValuesList(basisSet,quadGrid(:,i))
          moVal1 = dot_product(moCoeffs(:,j),basisValues)
          moVal2 = moVal1
!          moVal2 = dot_product(moCoeffs(:,1),basisValues)
          quadValues(i) = moVal1*moVal2
        end do
      !$omp end parallel do
        write(iOut,*) ' MO, Integral = ',j,dot_product(quadWeights,quadValues)
      end do
      write(iOut,*)
!
!     The end of the program.
!
  999 Continue
      call CPU_TIME(tEnd)
      write(iOut,*)' TIME = ',tEnd-tStart
      write(iOut,*)
      end program gbs
