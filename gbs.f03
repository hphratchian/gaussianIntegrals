      program gbs
!
!     This program tests the MQC_Integrals module.
!
!     Hrant P. Hratchian, 2024.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      USE iso_fortran_env
      USE MQC_Integrals1
      implicit none
      integer,parameter::lVal=1
      integer::i,j,k,l
      logical::fail=.false.,atEnd
      real(kind=real64),dimension(lVal+1,lVal+1)::tmpSij=float(0)
      logical,dimension(lVal+1,lVal+1)::haveSij=.false.
      type(mqc_cgtf)::bf1,bf2,bf3,bf4
!
      integer(kind=int64)::ixyz,jxyz
      real(kind=real64)::mu,p,xAB,xPA,xPB
!
      integer(kind=int64)::nCGTF=0,iCGTF=0,nBasis
      integer(kind=int64),dimension(:),allocatable::CGTF2IBasis,lArray
      real(kind=real64),dimension(:),allocatable::normConstants
      type(mqc_cgtf)::bfTmp
      type(mqc_cgtf),dimension(:),allocatable::basisSetList
      type(mqc_linkedList),pointer::basisSet,basisSetCurrentNode
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
      write(iOut,*)'    mu   = ',mu
      write(iOut,*)'    p    = ',p
      write(iOut,*)'    xAB  = ',xAB
      write(iOut,*)'    xPA  = ',xPA
      write(iOut,*)'    xPB  = ',xPB
      call MQC_Overlap_Distribution_Primitive_XYZ_Constants(bf1,bf1, &
        1,1,mu,p,xAB,xPA,xPB)
      write(iOut,*)
      write(iOut,*)' Intermediates...2:'
      write(iOut,*)'    mu   = ',mu
      write(iOut,*)'    p    = ',p
      write(iOut,*)'    xAB  = ',xAB
      write(iOut,*)'    xPA  = ',xPA
      write(iOut,*)'    xPB  = ',xPB

      call MQC_Overlap_Primitive_XYZ_OS(ixyz,jxyz,mu,p,xAB,xPA,xPB,tmpSij,haveSij)
      call mqc_print(haveSij,iOut,header='The HaveSij Matrix:')
      call mqc_print(tmpSij,iOut,header='The Sij Matrix:')

!
!     Try to using the basisSet linked list.
!
      call LinkedList_Push(basisSet,bf1)
      call LinkedList_Push(basisSet,bf2)
      call LinkedList_Push(basisSet,bf3)
      call LinkedList_Push(basisSet,bf4)

      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Hrant - Trying out the linked list approach for building up the basis set...'
      basisSetCurrentNode => basisSet
      call LinkedList_Return_CGTF_Value(basisSetCurrentNode,bfTmp)
      call MQC_CGTF_print(bfTmp,iOut)
      nCGTF = nCGTF+1
      call LinkedList_GetNext(basisSetCurrentNode,atEnd,last_looks_ahead=.false.)
      do while(.not.atEnd)
        call LinkedList_Return_CGTF_Value(basisSetCurrentNode,bfTmp)
        call MQC_CGTF_print(bfTmp,iOut)
        nCGTF = nCGTF+1
        call LinkedList_GetNext(basisSetCurrentNode,atEnd,last_looks_ahead=.false.)
      endDo
      write(iOut,*)
      write(iOut,*)' Hrant - nCGTF = ',nCGTF
!
      Allocate(basisSetList(nCGTF),CGTF2IBasis(nCGTF))
      iCGTF = 0
      basisSetCurrentNode => basisSet
      call LinkedList_Return_CGTF_Value(basisSetCurrentNode,bfTmp)
      iCGTF = iCGTF+1
      basisSetList(iCGTF) = bfTmp
      call LinkedList_GetNext(basisSetCurrentNode,atEnd,last_looks_ahead=.false.)
      do while(.not.atEnd)
        call LinkedList_Return_CGTF_Value(basisSetCurrentNode,bfTmp)
        iCGTF = iCGTF+1
        basisSetList(iCGTF) = bfTmp
        call LinkedList_GetNext(basisSetCurrentNode,atEnd,last_looks_ahead=.false.)
      endDo
      nBasis = 0
      do i = 1,nCGTF
        CGTF2IBasis(i) = nBasis+1
        nBasis = nBasis+basisSetList(i)%shell2nBasis()
      endDo
      write(iOut,*)' Hrant - nBasis = ',nBasis
      Allocate(normConstants(nBasis))
      if(Allocated(lArray)) DeAllocate(lArray)


      goto 999
!hph+
      goto 999
!hph-


!
!     Work on the logic for building the set of lVectors given a total angular
!     momentum.
!
      write(iOut,*)
      l = 2
      do i = l,0,-1
        do j = (l-i),0,-1
          k = l-i-j
          write(iOut,*)i,j,k
        endDo
      endDo

!
!     The end of the program.
!
  999 Continue
      end program gbs
