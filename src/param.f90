#ifdef SINGLE_PRECISION
#  define sp          4
#  define MPI_REAL_RP MPI_REAL
#else
#  define sp          8
#  define MPI_REAL_RP MPI_DOUBLE_PRECISION
#endif

# ifdef THREE_DIMENSION
# define dim3 2
# else
# define dim3 3
# endif
