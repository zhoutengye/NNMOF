#ifdef SINGLE_PRECISION
#  define setsp           4
#  define MPI_REAL_SP  MPI_REAL
#  define HDF5_REAL_SP H5T_NATIVE_REAL
#  define ToFloat      float
#else
#  define setsp           8
#  define MPI_REAL_SP  MPI_DOUBLE_PRECISION
#  define HDF5_REAL_SP H5T_NATIVE_DOUBLE
#  define ToFloat      dble
#endif

# define NOTZERO  1.0e-12
# define Pi  3.1415926535897932d0
