/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef __MEAN_COMPASS_TYPES_H__
#define __MEAN_COMPASS_TYPES_H__

#include <algorithm>
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

// The dynamic part of the config.
struct DynamicConfig {
 public:
  DynamicConfig() = default;
 protected:
  DynamicConfig(const DynamicConfig&) = default;
  DynamicConfig& operator=(const DynamicConfig&) = default;
 public:
  // We do not delete move constructor and operator because it
  // won't copy the expensive data structures like string or vector.

# define __MEAN_COMPASS__OPTION(type, name) \
 private: \
  type name ## _; \
 public: \
  DynamicConfig& name(type param) {  /* Temporary. */ \
    using std::swap; \
    swap(name ## _, param); \
    return *this; \
  } \
  const type& name() const { return name ## _; } \
  type& name() { return name ## _; } \
  /* End of __MEAN_COMPASS__OPTION macro. */

  __MEAN_COMPASS__OPTION(std::string, config_file_name)
  __MEAN_COMPASS__OPTION(std::vector<std::string>, input_files)
  __MEAN_COMPASS__OPTION(bool,        parity)
  __MEAN_COMPASS__OPTION(int,         default_precision)
  __MEAN_COMPASS__OPTION(int,         display_precision)
  __MEAN_COMPASS__OPTION(std::string, barrier_multiplier_str)
  __MEAN_COMPASS__OPTION(bool,        barrier_adjustment)

# undef __MEAN_COMPASS__OPTION
};

template<
    bool option_verbose = false,
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
struct Config : public DynamicConfig {
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

  constexpr static const bool verbose = option_verbose;
  constexpr static const bool use_colors = option_use_colors;

  Config() = default;
  Config(const Config&) = delete;
  Config& operator=(const Config&) = delete;

  Config(const DynamicConfig& param) : DynamicConfig(param) { }
  Config(DynamicConfig&& param) :
      DynamicConfig(std::forward<DynamicConfig>(param)) {
    // Nothing here.
  }
  // We do not delete move constructor and operator because it
  // won't copy the expensive data structures like string or vector.

  // Forward everything else to DynamicConfig,
  // we explictly state the first argument, so that it does not interfere
  // with the copy constructor.
  template<typename ...Args>
  Config(std::string&& config_file_name, Args&&... args) :
    DynamicConfig(std::forward<std::string>(config_file_name),
                  std::forward<Args>(args)...) {
    // Nothing here.
  }
};

}  // namespace mean_compass

#endif  // __MEAN_COMPASS_TYPES_H__

// vim: et sw=2 ts=2
