      program pad
!
!     Evaluate photoelectron angular distribution intensities, I(theta), and
!     anisotropy parameters for a selected Dyson orbital read from a Gaussian
!     FAF. The production path currently uses a linearly propagating plane-wave
!     photoelectron in the length-gauge dipole approximation.
!
!     Command line arguments are supplied as -option value pairs:
!       -faf FILE -dyson-mo N -photon-ev EV -binding-ev EV
!       [-n-theta N] [-n-grid N] [-pe-type N]
!       [-lab-frame cartesian|sphere|axisymmetric]
!       [-lab-theta N] [-lab-phi N] [-lab-alignment A] [-n-chi N]
!     Programmatic callers may also use PAD_LAB_FRAMES_CUSTOM=-1 with
!     user-supplied lab-frame vector arrays.
!
!     Legacy positional arguments are still accepted for existing scripts.
!
!
!     Hrant P. Hratchian, 2025, 2026.
!     hhratchian@umcerced.edu
!     University of California, Merced
!
!
      use pad_mod
      use omp_lib
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
 8999 format(/,1x,'Wall Time: ',f15.3,' s',/,1x,'PAD Complete.')
!
!
!     Begin the program.
!
      tStart = omp_get_wtime()
      call padCommandLine(options,fafName)
      call omp_set_num_threads(options%nOMP)
      write(iOut,1000)
      call mqc_version_print(iOut)
      call padPrintReproducibleCommand(fafName,options)
      call padPrintOpenMPSettings(options)
!
      if(MEMChecks) call print_memory_usage(iOut,'At top of PAD.')
!
!     Load the FAF and dispatch the reusable PAD driver.
!
      call faf%load(fafName)
      call runPADCalculation(faf,options,results)
!
      tEnd = omp_get_wtime()
      write(iOut,8999) tEnd-tStart
      if(MEMChecks) call print_memory_usage(iOut,'End of PAD.')
      end program pad
