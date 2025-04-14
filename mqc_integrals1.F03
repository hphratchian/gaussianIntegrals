      module MQC_Integrals1
!
!     ***********************************************************************
!     ***********************************************************************
!     **                                                                   **
!     **               The Merced Quantum Chemistry Package                **
!     **                            (MQCPack)                              **
!     **                                                                   **
!     ** Written By:                                                       **
!     **    Lee M. Thompson, Xianghai Sheng, Andrew D. Mahler, Dave        **
!     **    Mullally, and Hrant P. Hratchian                               **
!     **                                                                   **
!     **                          Version 25.2.1                           **
!     **                         February 9, 2025                          **
!     **                                                                   **
!     **                                                                   **
!     ** Modules beloning to MQCPack:                                      **
!     **    1. MQC_General                                                 **
!     **    2. MQC_Binary                                                  **
!     **    3. MQC_DataStructures                                          **
!     **    4. MQC_Integrals                                               **
!     **    5. MQC_Algebra                                                 **
!     **    6. MQC_Files                                                   **
!     **    7. MQC_Molecule                                                **
!     **    8. MQC_EST                                                     **
!     **    9. MQC_Gaussian                                                **
!     **                                                                   **
!     ***********************************************************************
!     ***********************************************************************
!
!
!     This module contains objects and procedures for working with and
!     evaluating certain types of integrals over gaussian-type basis functions.
!
!
      use MQC_General
      use MQC_DataStructures
      use iso_fortran_env, only: int64, real64
      implicit none
      integer(kind=int64),parameter::iOut=6
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!     MQC_BasisSet
!
!     The MQC_BasisSet class is an object of a list of MQC_CGTF members to build
!     up a full basis set.
      type MQC_BasisSet
        integer(kind=int64)::nBasis=0
      end type MQC_BasisSet



!
!     MQC_CGTF
!
!     The MQC_CGTF class contains a contracted Gaussian-type function.
      type MQC_CGTF
        logical,private::initialized=.false.
        integer(kind=int64),private::angularMomentum
        integer(kind=int64),private::nPrimMax,nPrim,nBasis=0
        integer(kind=int64),dimension(:,:),allocatable,private::lArrays
        real(kind=real64),dimension(3)::center
        real(kind=real64),dimension(:),allocatable::coefficient
        real(kind=real64),dimension(:),allocatable::alpha
        real(kind=real64),dimension(:,:),allocatable::primitiveNorm
      Contains
        procedure,public::init                  => MQC_CGTF_init
        procedure,public::print                 => MQC_CGTF_print
        procedure,public::shell2nBasis          => MQC_CGTF_nBasisFunctions
        procedure,public::primitiveSelfOverlap  => MQC_CGTF_primitiveSelfOverlap
      end type MQC_CGTF
!
!
      CONTAINS
!
!--------------------------------------------------------------------------
!                                                                         |
!     OBJECT SET-UP, PRINTING, AND GENERAL MANIPULATION PROCEDURES        |
!                                                                         |
!--------------------------------------------------------------------------
!
!
!PROCEDURE MQC_CGTG_init
      subroutine MQC_CGTF_init(myCGTF,angularMomentum,center,coefficients,  &
        alphas)
!
!     This routine initializes a contracted Gaussian-type basis function object.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64),intent(in)::angularMomentum
      real(kind=real64),dimension(3),intent(in)::center
      real(kind=real64),dimension(:)::coefficients,alphas
      integer(kind=int64)::i,j
!
!     Begin by making sure the object's allocatable arrays aren't already
!     allocated.
!
      if(Allocated(myCGTF%coefficient)) DeAllocate(myCGTF%coefficient)
      if(Allocated(myCGTF%alpha)) DeAllocate(myCGTF%alpha)
      if(Allocated(myCGTF%primitiveNorm)) DeAllocate(myCGTF%primitiveNorm)
!
!     Fill scalars and center array. Also, allocate the allocatable arrays in
!     the object. Also, set the initialization flage to TRUE.
!
      myCGTF%angularMomentum = angularMomentum
      myCGTF%nPrim           = Size(coefficients)
      myCGTF%center          = center
      Allocate(myCGTF%coefficient(myCGTF%nPrim),myCGTF%alpha(myCGTF%nPrim))
!
!     Load the primitive coefficients and alphas.
!
      if(Size(coefficients).ne.Size(alphas))  &
        call mqc_error('MQC_CGTF_init: Coefficients and exponents must be same size.')
      myCGTF%coefficient = coefficients
      myCGTF%alpha = alphas
!
!     Allocate primitiveNorm and lArrays arrays. Then fill in both arrays.
!
      myCGTF%nBasis = ((angularMomentum+1)*(angularMomentum+2))/2
      myCGTF%lArrays = MQC_CGTF_extendLArray(myCGTF)
      Allocate(myCGTF%primitiveNorm(myCGTF%nPrim,myCGTF%nBasis))
      do i = 1,myCGTF%nBasis
        do j = 1,myCGTF%nPrim
          myCGTF%primitiveNorm(j,i) = myCGTF%primitiveSelfOverlap(j,myCGTF%lArrays(:,i))
        endDo
      endDo
      call mqc_print(myCGTF%primitiveNorm,iOut,header='Hrant - In init, primitiveNorm:')
!
!     Set the initialization flag to TRUE.
!
      myCGTF%initialized = .true.
!
      return
      end subroutine MQC_CGTF_init

!
!PROCEDURE MQC_CGTG_print
      subroutine MQC_CGTF_print(myCGTF,iOut)
!
!     This routine prints the data in a MQC_CGTF object.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64),intent(in)::iOut
      integer(kind=int64)::i
!
 1000 format(1x,'l=',i4,' shell at (',f8.4,',',f8.4,',',f8.4,'). NPrim=',i4)
 2000 format(3x,'Primitive Exponents and Coefficients.')
 2100 format(6x,f15.8,4x,f15.8)
!
      write(iOut,1000) myCGTF%angularMomentum,myCGTF%center,myCGTF%nPrim
      write(iOut,2000)
      do i = 1,myCGTF%nPrim
        write(iOut,2100) myCGTF%alpha(i),myCGTF%coefficient(i)
      endDo
!
      return
      end subroutine MQC_CGTF_print

!
!--------------------------------------------------------------------------
!                                                                         |
!     HIGH LEVEL PROCEDURES                                               |
!                                                                         |
!--------------------------------------------------------------------------
!
!PROCEDURE LinkedList_Return_CGTF_Value
      subroutine LinkedList_Return_CGTF_Value(LL_Cur,Return_Value)
!
!     This routine returns the value stored in a linked list node
!     when the value is of type MQC_CGTF.
!
!     H. P. Hratchian, 2025.
!
      implicit none
      type(MQC_LinkedList),pointer::LL_Cur
      type(MQC_CGTF),intent(out)::Return_Value
!
!     Check that the current node is associated and contains a CGTF object.
!
      if(Associated(LL_Cur)) then
        if(Allocated(LL_Cur%Val)) then
          select type (Val => LL_Cur%Val)
          type Is (MQC_CGTF)
            Return_Value = Val
          class default
            call mqc_error('Data types not compatible in LinkedList_Return_CGTF_Value')
          end Select
        endIf
      endIf
!
      return
      end subroutine LinkedList_Return_CGTF_Value

!
!PROCEDURE MQC_CGTG_nBasisFunctions
      function MQC_CGTF_nBasisFunctions(myCGTF)  result(nBasis)
!
!     This function returns the number of basis functions associated with the
!     shell described by object myCGTF.
!
      implicit none
      class(MQC_CGTF)::myCGTF
      integer(kind=int64)::nBasis
!
      nBasis = myCGTF%nBasis
!
      return
      end function MQC_CGTF_nBasisFunctions

!
!PROCEDURE MQC_CGTF_extendLArray
      function MQC_CGTF_extendLArray(myCGTF,ordering)  result(lArrays)
!
!     This function builds the list of angular momentum arrays, which give
!     (lx,ly,lz) triples, for the contracted Gaussian type function myCGTF. The
!     returned array, lArrays, is allocated (3,nBasisFunctions).
!
!     There are multiple ways to do the ordering. The OPTIONAL dummy argument
!     <ordering> is used to chose an ordering algorithm. The available options
!     are:
!           lexicographic     This ordering is the standard used in libint,
!                             PySCF, and MolPro. As an example, Cartesian d
!                             functions are ordered as (0,0,2), (0,1,1),
!                             (0,2,0), (1,0,1), (1,1,0), and (2,0,0).
!           gaussian          This ordering is the standard used in Gaussian. As
!                             an example, Cartesian d functions are ordered as
!                             (2,0,0), (0,2,0), (0,0,2), (1,1,0), (1,0,1), and
!                             (0,1,1). This is the default.
!
!
!     - H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      character(len=64),optional::ordering
      integer(kind=int64),dimension(:,:),allocatable::lArrays
      integer(kind=int64)::nBasisFunctions,i,lx,ly,lz
      character(len=64)::myOrdering
!
!     Start by setting myOrdering.
!
      if(PRESENT(ordering)) then
        call String_Change_Case(ordering,'l',myOrdering)
      else
        myOrdering = 'gaussian'
      endIf
!
!     Allocate lArrays and then fill it.
!
      nBasisFunctions = myCGTF%shell2nBasis()
      Allocate(lArrays(3,nBasisFunctions))
      lArrays = -1
      select case(myOrdering)
      case('lexicographic')
        i = 1
        do lx = 0,myCGTF%angularMomentum
          do ly = 0,myCGTF%angularMomentum-lx
            lz = myCGTF%angularMomentum-lx-ly
            lArrays(1,i) = lx
            lArrays(2,i) = ly
            lArrays(3,i) = lz
            i = i+1
          endDo
        endDo
      case('gaussian')
        if(myCGTF%angularMomentum.eq.0) then
          lArrays(1,1) = 0
          lArrays(2,1) = 0
          lArrays(3,1) = 0
        elseIf(myCGTF%angularMomentum.eq.1) then
          lArrays(1,1) = 1
          lArrays(2,1) = 0
          lArrays(3,1) = 0
          lArrays(1,2) = 0
          lArrays(2,2) = 1
          lArrays(3,2) = 0
          lArrays(1,3) = 0
          lArrays(2,3) = 0
          lArrays(3,3) = 1
        elseIf(myCGTF%angularMomentum.eq.2) then
          lArrays(1,1) = 2
          lArrays(2,1) = 0
          lArrays(3,1) = 0
          lArrays(1,2) = 0
          lArrays(2,2) = 2
          lArrays(3,2) = 0
          lArrays(1,3) = 0
          lArrays(2,3) = 0
          lArrays(3,3) = 2
          lArrays(1,4) = 1
          lArrays(2,4) = 1
          lArrays(3,4) = 0
          lArrays(1,5) = 1
          lArrays(2,5) = 0
          lArrays(3,5) = 1
          lArrays(1,6) = 0
          lArrays(2,6) = 1
          lArrays(3,6) = 1
        elseIf(myCGTF%angularMomentum.eq.3) then
          lArrays(1,1) = 3
          lArrays(2,1) = 0
          lArrays(3,1) = 0
          lArrays(1,2) = 0
          lArrays(2,2) = 3
          lArrays(3,2) = 0
          lArrays(1,3) = 0
          lArrays(2,3) = 0
          lArrays(3,3) = 3
          lArrays(1,4) = 1
          lArrays(2,4) = 2
          lArrays(3,4) = 0
          lArrays(1,5) = 2
          lArrays(2,5) = 1
          lArrays(3,5) = 0
          lArrays(1,6) = 2
          lArrays(2,6) = 0
          lArrays(3,6) = 1
          lArrays(1,7) = 1
          lArrays(2,7) = 0
          lArrays(3,7) = 2
          lArrays(1,8) = 0
          lArrays(2,8) = 1
          lArrays(3,8) = 2
          lArrays(1,9) = 0
          lArrays(2,9) = 2
          lArrays(3,9) = 1
          lArrays(1,10) = 1
          lArrays(2,10) = 1
          lArrays(3,10) = 1
        elseIf(myCGTF%angularMomentum.gt.3) then
          call mqc_error('MQC_CGTF_extendLArray: Gaussian ordering NYI for L>3.')
        endIf
      case default
        call mqc_error('MQC_CGTF_extendLArray: Invalid ordering sent.')
      endSelect
!
      return
      end function MQC_CGTF_extendLArray

!
!PROCEDURE MQC_CGTG_primitiveSelfOverlap
      function MQC_CGTF_primitiveSelfOverlap(myCGTF,i,lVector)  result(selfOverlap)
!
!     This routine evaluates and returns the self-overlap of the i-th primitive
!     basis function in the shell of myCGTF with angular momentum vector
!     lVector.
!
!     -H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::myCGTF
      integer(kind=int64),intent(in)::i
      integer(kind=int64),dimension(3),intent(in)::lVector
      real(kind=real64)::selfOverlap
      real(kind=real64)::fourAlpha,rootPiFourA3
!
!     Check that lVector is valid.
!
      if(Sum(lVector).ne.myCGTF%angularMomentum) call mqc_error('MQC_CGTF_selfOverlap: Invalid lVector sent.')
!
!     Evaluate the self-overlap.
!
      fourAlpha = 4*myCGTF%alpha(i)
      rootPiFourA3 = (SQRT(Pi/(2*myCGTF%alpha(i))))**3
      if(lVector(1).eq.0) then
        selfOverlap = float(1)
      else
        selfOverlap = mqc_doubleFactorial((2*lVector(1)-1))/(fourAlpha**lVector(1))
      endIf
      if(lVector(2).gt.0) selfOverlap = selfOverlap*(mqc_doubleFactorial((2*lVector(2)-1))/(fourAlpha**lVector(2)))
      if(lVector(3).gt.0) selfOverlap = selfOverlap*(mqc_doubleFactorial((2*lVector(3)-1))/(fourAlpha**lVector(3)))
      selfOverlap = selfOverlap*rootPiFourA3
!
      return
      end function MQC_CGTF_primitiveSelfOverlap

!
!--------------------------------------------------------------------------
!                                                                         |
!     INTERMEDIATE AND LOW LEVEL ROUTINES                                 |
!                                                                         |
!--------------------------------------------------------------------------
!
!PROCEDURE MQC_Overlap_Distribution_Constants
      subroutine MQC_Overlap_Distribution_Primitive_XYZ_Constants(Ga,Gb, &
        iPrimA,iPrimB,mu,p,xAB,xPA,xPB)
!
!     This routine is sent two primitive contracted Gaussian basis functions and
!     returns the Gaussian Overlap Distribution intermediate constants (mu, p,
!     xAB, xPA, and xPB) for the product of primitive number iPrimA from Ga and
!     iPrimB from Gb.
!
!
!     - H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::Ga,Gb
      integer(kind=int64),intent(in)::iPrimA,iPrimB
      real(kind=real64),intent(out)::mu,p,xAB,xPA,xPB
!
!     Do the work...
!
      mu = ga%alpha(iPrimA)*gb%alpha(iPrimB)/  &
        (ga%alpha(iPrimA)+gb%alpha(iPrimB))
      p = ga%alpha(iPrimA)+gb%alpha(iPrimB)
      xAB = ga%center(iPrimA)-gb%center(iPrimB)
      xPA = -xAB*gb%alpha(iPrimB)/p
      xPB =  xAB*ga%alpha(iPrimA)/p
!
      return
      end subroutine MQC_Overlap_Distribution_Primitive_XYZ_Constants

!
!PROCEDURE MQC_Overlap_CGTF
      subroutine MQC_Overlap_CGFT(bra,ket,basisIntegrals)
!
!     This subroutine calculates an overlap integral between two contracted
!     Gaussian type functions. These are sent as MQC_CGTF objects. The output is
!     basisIntegrals, which is a matrix of the overlap integrals over the
!     expansion of basis functions that come from the contracted Gaussians in
!     bra and ket. At present, all overlap integrals are solved using the
!     Obara-Saika recursion relation.
!
!     - H. P. Hratchian, 2025.
!
      implicit none
      class(MQC_CGTF),intent(in)::bra,ket
      real(kind=real64),dimension(:,:),allocatable,intent(out)::basisIntegrals
      integer(kind=int64)::iBasisBra,iBasisKet,nBasisBra,nBasisKet
      real(kind=real64),dimension(:),allocatable::primInts,primCoeffs
      real(kind=real64),dimension(:,:),allocatable::basisIntegralsX,  &
        basisIntegralsY,basisIntegralsZ
!
!     Figure out the dimension of basisIntegrals and allocate the matrix.
!
      nBasisBra = bra%nBasis
      nBasisKet = ket%nBasis
      Allocate(basisIntegrals(nBasisBra,nBasisKet),  &
        basisIntegralsX(nBasisBra,nBasisKet),  &
        basisIntegralsY(nBasisBra,nBasisKet),  &
        basisIntegralsZ(nBasisBra,nBasisKet))
!
!     Form the overlap distribution terms and call the overlap primitives OS
!     routine for x, y, and z factorized terms.
!
      do iBasisBra = 1,nBasisBra
        do iBasisKet = 1,nBasisKet

!HRANT

      call MQC_Overlap_Distribution_Primitive_XYZ_Constants(bf1,bf1, &
        1,1,mu,p,xAB,xPA,xPB)
      call MQC_Overlap_Primitive_XYZ_OS(ixyz,jxyz,mu,p,xAB,xPA,xPB,tmpSij,haveSij)



!
!PROCEDURE MQC_Overlap_Primitive_XYZ_OS
      recursive subroutine MQC_Overlap_Primitive_XYZ_OS(ixyz,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
!
!     This recursive subroutine solves a 1-D Cartesian overlap between two
!     Cartesian Gaussian-type primitive functions using the Obara-Saika
!     recurrence relation.
!
!
!     -H. P. Hratchian, 2025.
!
      implicit none
      integer(kind=int64)::ixyz,jxyz
      real(kind=real64)::mu,p,XAB,XPA,XPB
      real(kind=real64),dimension(:,:)::Sij
      logical,dimension(:,:)::haveSij
!
!     Look at ixyz and jxyz to determine what we need to do here. There are 9
!     unique cases.
!
!     0,0
      if(ixyz.eq.0.and.jxyz.eq.0) then
        if(haveSij(1,1)) then
          write(*,*)' Hrant - asked for a value I already have (0,0).'
          return
        endIf
        Sij(1,1) = SQRT(pi/p)*exp(-mu*XAB**2)
        haveSij(1,1) = .true.
!
!     0,1
      elseIf(ixyz.eq.0.and.jxyz.eq.1) then
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(1,2) = XPB*Sij(1,1)
        haveSij(1,2) = .true.
!
!     0,>=2
      elseIf(ixyz.eq.0.and.jxyz.ge.2) then
        if(.not.haveSij(1,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,jxyz-1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-2,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(1,jxyz+1) = XPB*Sij(1,jxyz)+Sij(1,jxyz-1)*(jxyz-1)/(float(2)*p)
        haveSij(1,jxyz+1) = .true.
!
!     1,0
      elseIf(ixyz.eq.1.and.jxyz.eq.0) then
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,1) = XPA*Sij(1,1)
        haveSij(2,1) = .true.
!
!     >=2,0
      elseIf(ixyz.ge.2.and.jxyz.eq.0) then
        if(.not.haveSij(ixyz,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,1) = XPA*Sij(ixyz,1)+Sij(ixyz-1,1)*(ixyz-1)/(float(2)*p)
        haveSij(ixyz+1,1) = .true.
!
!     1,1  ==>  0,1
      elseIf(ixyz.eq.1.and.jxyz.eq.1) then
        if(.not.haveSij(1,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,2) = XPA*Sij(1,2)+Sij(1,1)/(float(2)*p)
        haveSij(2,2) = .true.
!
!     1,>=2
      elseIf(ixyz.eq.1.and.jxyz.ge.2) then
        if(.not.haveSij(1,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(1,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(0,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(2,jxyz+1) = XPA*Sij(1,jxyz+1)+jxyz*Sij(1,jxyz)/(float(2)*p)
        haveSij(2,jxyz+1) = .true.
!
!     >=2,1
      elseIf(ixyz.ge.2.and.jxyz.eq.1) then
        if(.not.haveSij(ixyz,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,2))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz,1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,0,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,2) = XPA*Sij(ixyz,2)+(ixyz*Sij(ixyz-1,2)+Sij(ixyz,1))/(float(2)*p)
        haveSij(ixyz+1,2) = .true.
!
!     >=2,>=2
      elseIf(ixyz.ge.2.and.jxyz.ge.2) then
        if(.not.haveSij(ixyz,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz-1,jxyz+1))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-2,jxyz,mu,p,XAB,XPA,XPB,Sij,haveSij)
        if(.not.haveSij(ixyz,jxyz))  &
          call MQC_Overlap_Primitive_XYZ_OS(ixyz-1,jxyz-1,mu,p,XAB,XPA,XPB,Sij,haveSij)
        Sij(ixyz+1,jxyz+1) = XPA*Sij(ixyz,jxyz+1)+((ixyz)*Sij(ixyz-1,jxyz+1)+(jxyz+1)*Sij(ixyz,jxyz))/(float(2)*p)
        haveSij(ixyz+1,jxyz+1) = .true.
!
!     Uknown combination...fail out.
      else
        call mqc_error('MQC_Overlap_Primitive_XYZ_OS: Invalid i and j.')
      endIf
!
      return
      end subroutine MQC_Overlap_Primitive_XYZ_OS

!
!
      end module MQC_Integrals1
