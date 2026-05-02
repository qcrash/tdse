cmake_minimum_required(VERSION 3.1)

project(tdse VERSION 0.1
        DESCRIPTION "Numerical Solution to the 1D Time Dependent Schroedinger Equation"
        LANGUAGES Fortran)
enable_language(Fortran)

# Currently setting the Fortran compiler to use -std=gnu, change this if you
# want a specific standard
set(FVERSION "-std=gnu")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION} -fdefault-integer-8")
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG "TRUE")
find_package(${CMAKE_MATH_LIBRARY} REQUIRED)
#find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)

# Source code
set(MY_MODS   dump_psi.f90 propagate.f90 scalar.f90 ham_psi.f90
segal_bargmann.f90 propagate_green.f90 propagate_ab.f90
	propagate_trap.f90 variance.f90 propagate_central.f90 pval.f90 wigner.f90
	propagate_convert.f90 pvar.f90 xval.f90 getcli.f90 check_arg_int.f90 check_arg_double.f90 dump_repr.f90 check_arg_char.f90  #propagate_fft.f90
#to do: adapt FFT to FFTW
CACHE INTERNAL  "")
       	      add_library(mylib "${MY_MODS}")
add_executable(tdse tdse.f90)
target_link_libraries(tdse PRIVATE mylib PRIVATE ${CMAKE_MATH_LIBRARY}::${CMAKE_MATH_LIBRARY_NAMESPACE} ${CMAKE_DL_LIBS}  m)

#target_link_option(tdse PRIVATE "-fdefault-integer-8")

install(TARGETS tdse)
