#pragma once

// Optional build flags for embedded targets:
//
//   -DRIO_EIGEN_NO_MALLOC
//   -DRIO_EIGEN_DONT_VECTORIZE
//   -DRIO_EIGEN_DONT_ALIGN_STATICALLY
//

#if defined(RIO_EIGEN_NO_MALLOC)
  #ifndef EIGEN_NO_MALLOC
    #define EIGEN_NO_MALLOC
  #endif
#endif

#if defined(RIO_EIGEN_DONT_ALIGN_STATICALLY)
  #ifndef EIGEN_DONT_ALIGN_STATICALLY
    #define EIGEN_DONT_ALIGN_STATICALLY
  #endif
#endif

#if defined(RIO_EIGEN_DONT_VECTORIZE)
  #ifndef EIGEN_DONT_VECTORIZE
    #define EIGEN_DONT_VECTORIZE
  #endif
#endif