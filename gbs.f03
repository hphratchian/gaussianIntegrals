Include "gbs_mod.f03"
      program gbs
!
!     This program tests the MQC_Integrals module.
!
!     Hrant P. Hratchian, 2024.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      USE gbs_mod
      implicit none
      integer,parameter::lVal=0
      integer::i,j,k,l
      logical::fail=.false.,atEnd
      real(kind=real64),dimension(lVal+1,lVal+1)::tmpSij=float(0)
      character(len=256)::fafName
      logical,dimension(lVal+1,lVal+1)::haveSij=.false.
      type(mqc_cgtf)::bf1,bf2,bf3,bf4
      type(mqc_gaussian_unformatted_matrix_file)::faf
      type(MQC_basisSet)::basisSet
!
      integer(kind=int64)::ixyz,jxyz
      real(kind=real64)::mu,p
      real(kind=real64),dimension(3)::xAB,xPA,xPB
!
      integer(kind=int64)::nCGTF=0,iCGTF=0,nBasis
      integer(kind=int64),dimension(:),allocatable::CGTF2IBasis,lArray
      real(kind=real64),dimension(:),allocatable::normConstants
      real(kind=real64),dimension(:,:),allocatable::basisIntegrals
!
!     Format statements.
!
 1000 format(1x,'Program gbs.')
!
!
!     Begin the program.
!
      fail = .false.
      write(iOut,1000)
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
!      write(IOut,*)
!      write(iOut,*)' Running overlap of 2 with 2.'
!      call MQC_Overlap_CGFT(bf2,bf2,basisIntegrals)
!      write(IOut,*)
!      write(iOut,*)' Running overlap of 1 with 2.'
!      call MQC_Overlap_CGFT(bf1,bf2,basisIntegrals)
      write(IOut,*)
      write(iOut,*)' ----------------------------------------------'
      write(iOut,*)' Running overlap of 1 with 3.'
      call MQC_Overlap_CGFT(bf1,bf3,basisIntegrals)
!      write(IOut,*)
!      write(iOut,*)' Running overlap of 1 with 4.'
!      call MQC_Overlap_CGFT(bf1,bf4,basisIntegrals)


!
!     Try reading a basis set from faf.
!
      call loadGaussianBasisSet(faf,basisSet)

      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Using Basis Set Object...'
      write(IOut,*)
      write(iOut,*)' ----------------------------------------------'
      write(iOut,*)' Running overlap of 1 with 1.'
      call MQC_Overlap_CGFT(basisSet%shells(1),basisSet%shells(1),basisIntegrals)
      write(IOut,*)
      write(iOut,*)' ----------------------------------------------'
      write(iOut,*)' Running overlap of 1 with 2.'
      call MQC_Overlap_CGFT(basisSet%shells(1),basisSet%shells(2),basisIntegrals)
      write(IOut,*)
      write(iOut,*)' ----------------------------------------------'
      write(iOut,*)' Running overlap of 2 with 1.'
      call MQC_Overlap_CGFT(basisSet%shells(2),basisSet%shells(1),basisIntegrals)
      write(IOut,*)
      write(iOut,*)' ----------------------------------------------'
      write(iOut,*)' Running overlap of 2 with 2.'
      call MQC_Overlap_CGFT(basisSet%shells(2),basisSet%shells(2),basisIntegrals)
      write(IOut,*)
      write(iOut,*)' ----------------------------------------------'

!
!     The end of the program.
!
  999 Continue
      end program gbs
