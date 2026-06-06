      program unitTest1
!
!     This program tests the spherical integration and helper routines used in
!     PAD. In particular, it checks the spherical quadrature normalization used
!     with spherical harmonics and the limiting values of the spherical Bessel
!     functions.
!
!     Command line arguments:
!           No command line arguments are required.
!
!     Usage:
!           ./unitTest1.exe
!
!
!     H. P. Hratchian, 2025.
!
      use iso_fortran_env
      use omp_lib
      use mqc_general
      use dyson_matrix_elements_mod
      implicit none
      integer(kind=int64),parameter::iOut=6_int64
      integer(kind=int64),parameter::nOMP=1,nTheta=181,nPhi=361,l=2,m=0
      integer(kind=int64)::i,j
      real(kind=real64)::wTotal,totalNorm,x,tmpVal
      real(kind=real64),dimension(:),allocatable::thetaVals,phiVals
      real(kind=real64),dimension(:,:),allocatable::weights
      complex(kind=real64)::totalOrthog
!
!     Format statements.
!
 1000 format(1x,'Starting unitTest1')
 1500 format(1x,i5,3x,i5,3x,f20.10)
 2000 format(/,1x,'Sum of weights = ',f20.10,3x,'w/pi = ',f20.10)
 2010 format(1x,'Integeral Value = ',f20.10)
 2500 format(1x,'l=',i3,3x,'x=',f15.10,5x,'j_l(x)=',f15.10)
 2510 format(1x,'Exact in limit=',f15.10,5x,'j_l(x)=',f15.10)
!
!     Start the unit test program...
!
      write(iOut,1000)
      call omp_set_num_threads(nOMP)
      call mqc_version_print(iOut)
      write(*,*)' Pi = ',Pi
!
!     Allocate memory.
!
      Allocate(thetaVals(nTheta),phiVals(nPhi))
      Allocate(weights(nTheta,nPhi))
!
!     Load the integration grid.
!
      call generate_sph_grid(nTheta,nPhi,thetaVals,phiVals,weights)
!
!     Test the orthonormality of the spherical harmonics function.
!
      wTotal = mqc_float(0)
      totalNorm= mqc_float(0)
      totalOrthog = mqc_float(0)
      write(*,*)' totalOrthog init = ',totalOrthog
      do i = 1,nTheta
        do j = 1,nPhi
          wTotal = wTotal+weights(i,j)
          totalNorm = totalNorm+weights(i,j)*  &
            real(Ylm_complex(l,m,thetaVals(i),phiVals(j))*  &
            conjg(Ylm_complex(l,m,thetaVals(i),phiVals(j))),kind=real64)
          totalOrthog = totalOrthog+weights(i,j)*  &
            Ylm_complex(l,m,thetaVals(i),phiVals(j))*  &
            conjg(Ylm_complex(0,0,thetaVals(i),phiVals(j)))
        endDo
      endDo
      write(iOut,2000) wTotal,wTotal/Pi
      write(iOut,2010) totalNorm
      write(*,*)' totalOrthog = ',totalOrthog
!
!     Test the limiting values of the spherical Bessel functions.
!
      write(iOut,*)
      write(iOut,2500) 0,mqc_float(0),sph_bessel_j(0,mqc_float(0))
      write(iOut,2500) 1,mqc_float(0),sph_bessel_j(1,mqc_float(0))
      write(iOut,2500) 2,mqc_float(0),sph_bessel_j(2,mqc_float(0))
      x = mqc_float(50)
      tmpVal = Sin(x-mqc_float(l)*pi/mqc_float(2))/x
      write(iOut,2510) tmpVal,sph_bessel_j(l,x)
!
      end program unitTest1
