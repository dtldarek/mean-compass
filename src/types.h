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

namespace default_types {

// I find the `using` syntax way more natural than `typedef`.
using Index          = int;
using Integer        = boost::multiprecision::mpz_int;
using Rational       = boost::multiprecision::mpq_rational;
using Real           = boost::multiprecision::mpfr_float;
using Triplet        = Eigen::Triplet<Real>;

using DenseMatrix    = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
using DenseVector    = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using SparseMatrix   = Eigen::SparseMatrix<Real, Eigen::ColMajor, Index>;
using DiagonalMatrix = Eigen::DiagonalMatrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
// The parameter 0 denotes that no flag is used.
using SparseVector   = Eigen::SparseVector<Real, 0, Index>;

using SparseLU       = Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<Index>>;

}  // namespace default_types

template<
    bool option_use_colors = true,
    typename Index_     = default_types::Index,
    typename Integer_   = default_types::Integer,
    typename Rational_  = default_types::Rational,
    typename Real_      = default_types::Real,
    // If you use non-default Matrix_, please specify LU_ as well.
    typename Matrix_    = Eigen::SparseMatrix<Real_, Eigen::ColMajor, Index_>,
	typename Diagonal_  = Eigen::DiagonalMatrix<Real_, Eigen::Dynamic, Eigen::Dynamic>,
    typename Vector_    = Eigen::Matrix<Real_, Eigen::Dynamic, 1>,
    // Beware, the default assumes that Matrix_ works with SparseLU.
    typename LU_        = Eigen::SparseLU<Matrix_, Eigen::COLAMDOrdering<Index_>>,
    typename Weight_    = Real_>
class Config {
 public:
  using Index      = Index_;
  using Integer    = Integer_;
  using Rational   = Rational_;
  using Real       = Real_;
  using Matrix     = Matrix_;
  using Diagonal   = Diagonal_;
  using Triplet    = Eigen::Triplet<Real>;
  using Vector     = Vector_;
  using LU         = LU_;
  using Weight     = Weight_;
  static constexpr bool use_colors = option_use_colors;
};

}  // namespace mean_compass

#endif  // __MEAN_COMPASS_TYPES_H__
