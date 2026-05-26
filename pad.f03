      program pad
!
!     Evaluate photoelectron angular distribution intensities, I(theta), and
!     anisotropy parameters for a selected Dyson orbital read from a Gaussian
!     FAF. The production path currently uses a linearly propagating plane-wave
!     photoelectron in the length-gauge dipole approximation.
!
!     Command line arguments:
!           1. FAF filename
!           2. alpha MO number to use as the Dyson orbital
!           3. photon energy in eV
!           4. electron binding energy in eV
!           5. number of theta angles from 0 --> pi (optional; default=5)
!           6. fallback Cartesian grid points per axis (optional; default=101)
!           7. photoelectron model flag (optional; default=0)
!              0 = direct plane wave path
!              1 = plane wave path through dysonMatrixElement routines
!              2 = experimental free partial-wave path
!           8. lab-frame model flag (optional; default=0)
!              0 = 3 Cartesian lab-frame orientations
!              1 = sphere-grid lab-frame orientations
!           9. number of lab-frame theta points (optional; default=5)
!          10. number of lab-frame phi points (optional; default=8)
!          11. number of chi points from 0 --> 2*pi (optional; default=36)
!
!
!     Hrant P. Hratchian, 2025, 2026.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      use pad_mod
      implicit none
      real(kind=real64)::tStart,tEnd
      character(len=256)::fafName
      type(pad_options)::options
      type(pad_results)::results
      type(mqc_gaussian_unformatted_matrix_file)::faf
!
!     Format statements.
!
 1000 format(1x,'Program PAD.')
 8999 format(/,1x,'Job Time: ',f15.1,' s',/,1x,'PAD Complete.')
!
!
!     Begin the program.
!
      call CPU_TIME(tStart)
      call padCommandLine(options,fafName)
      call omp_set_num_threads(options%nOMP)
      write(iOut,1000)
      call mqc_version_print(iOut)
!
      if(MEMChecks) call print_memory_usage(iOut,'At top of PAD.')
!
!     Load the FAF and dispatch the reusable PAD driver.
!
      call faf%load(fafName)
      call runPADCalculation(faf,options,results)
!
      call CPU_TIME(tEnd)
      write(iOut,8999) tEnd-tStart
      if(MEMChecks) call print_memory_usage(iOut,'End of PAD.')
      end program pad
