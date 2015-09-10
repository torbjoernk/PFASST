# Toolchain File for Jureca @ Juelich Supercomputing Center
# Intel system
#
# Make sure you have loaded only the following modules:
#  - intel-para/2015.07
#  - FFTW/3.3.4
#  - Boost/1.58.0-Python-2.7.10
#

set(MPI_ROOT /usr/local/software/jureca/Stage3/software/Toolchain/iccifort/2015.3.187-GCC-bare-4.9.3/psmpi/5.1.4-1/)
set(MPI_C_COMPILER ${MPI_ROOT}/bin/mpicc)
set(MPI_CXX_COMPILER ${MPI_ROOT}/bin/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_ROOT}/bin/mpifort)

set(FFTW3_ROOT /usr/local/software/jureca/Stage3/software/Toolchain/intel-para/2015.07/FFTW/3.3.4)
set(FFTW3_INCLUDE ${FFTW3_ROOT}/include)
set(FFTW3_LIB ${FFTW3_ROOT}/lib)
