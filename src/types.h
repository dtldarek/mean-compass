/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef __MEAN_COMPASS_TYPES_H__
#define __MEAN_COMPASS_TYPES_H__

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <Eigen/OrderingMethods>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>


namespace mean_compass {

typedef boost::multiprecision::cpp_int             Integer;
typedef boost::multiprecision::cpp_rational        Rational;
typedef boost::multiprecision::mpfr_float          Real;
typedef Eigen::Triplet<Real>                       Triplet;
typedef int                                        Index;

typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>  DenseMatrix;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1>               DenseVector;
typedef Eigen::SparseMatrix<Real, Eigen::ColMajor, Index>    SparseMatrix;
// The parameter 0 denotes that no flag is used.
typedef Eigen::SparseVector<Real, 0, Index>                  SparseVector;

typedef Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<Index>> SparseLU;

}  // namespace mean_compass

#endif  // __MEAN_COMPASS_TYPES_H__
