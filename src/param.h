#ifdef SINGLE_PRECISION
#  define sp           4
#  define SP           4
#  define MPI_REAL_SP  MPI_REAL
#  define HDF5_REAL_SP H5T_NATIVE_REAL
#  define ToFloat      float
#else
#  define sp           8
#  define SP           8
#  define MPI_REAL_SP  MPI_DOUBLE_PRECISION
#  define HDF5_REAL_SP H5T_NATIVE_DOUBLE
#  define ToFloat      dble
#endif

