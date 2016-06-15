
Install the dependencies:

This ought to keep you busy...

(1) Compilers:
      (a) C and an MPI C compiler
      (b) Fortran90 and an MPI Fortran90 compiler
      (c) C++ and an MPI c++ compiler
    I recommend OpenMPI (https://www.open-mpi.org/)
(2) zlib - this is a hdf5 dependency and can be substituted for szip:
      http://www.zlib.net/
(3) hdf5 (>= 1.18.16):
      https://www.hdfgroup.org/downloads/index.html
(4) metis (>= 5.1) - graph partitioner:
      http://glaros.dtc.umn.edu/gkhome/metis/metis/download
(5) parmetis - parallel graph partitioner (>= 4.0):
      http://glaros.dtc.umn.edu/gkhome/metis/metis/download
(6) scotch and ptscotch - more graph partitioners (>= 6)
      https://www.labri.fr/perso/pelegrin/scotch/
(7) BLAS and LAPACK:
      http://www.netlib.org/lapack/
(8) BLACS and ScaLapack (>= 2.0.2):
      http://www.netlib.org/scalapack/
(9) HSS matrix factorization software
      http://portal.nersc.gov/project/sparse/strumpack/
(10) ini file parser library
      https://github.com/ndevilla/iniparser

Select/modify the appropriate Makefile.inc template in make_inc.  
Copy it to the root diretory and rename it Makefile.inc

Make the application.


