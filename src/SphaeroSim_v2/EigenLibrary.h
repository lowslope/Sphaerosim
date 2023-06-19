#ifndef _SPHAEROSIM_EIGEN_LIBRARY_H_
#define _SPHAEROSIM_EIGEN_LIBRARY_H_

#ifndef _INCLUDED_STDINT_H_

#include <stdint.h>

#define _INCLUDED_STDINT_H_
#endif

#ifndef _INCLUDED_DENSE_H_

#include <Dense>

#define _INCLUDED_DENSE_H_
#endif

namespace Eigen {
    typedef Eigen::Matrix<std::size_t, 3, 1, 0, 3, 1> Vector3s;
    typedef Eigen::Matrix<int8_t, 3, 1, 0, 3, 1> Vector3i8;
    typedef Eigen::Matrix<int64_t, 3, 1, 0, 3, 1> Vector3i64;
}

#endif