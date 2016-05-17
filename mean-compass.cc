/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <Eigen/OrderingMethods>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

namespace {

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

// For any integer-like type T returns a random number of that type
// that is uniformly distributed in [0, size).
template<typename I = Integer> I unirand(I size) {
  static boost::random_device random_device;
  static boost::random::mt19937 generator(random_device);
  assert(size >= 1);
  boost::random::uniform_int_distribution<I> uniform_integer(0, size - 1);
  return uniform_integer(generator);
}

// Returns a random Real number uniformly distributed in [0,1).
Real unirand() {
  Integer denominator = Integer(1) << Real::default_precision();
  Integer nominator = unirand(denominator);
  Rational result(nominator, denominator);
  return static_cast<Real>(result);
}

// Check if the value is between min and max and throw an exception if not.
template<typename T, const T min, const T max>
void check_range(const T& value) {
  if (value < min || value > max) {
    std::stringstream description;
    description
      << "The value " << value << " is outside of "
      << "[" << min << ", " << max << "] range";
    throw std::out_of_range(description.str());
  }
}

}  // anonymous namespace


int main(int argc, char** argv) {
  // Parse cmdline flags. {{{
  int default_precision = 256;
  std::string config_file_name;

  try {
    namespace po = boost::program_options;
    po::options_description generic_options("Generic options");
    generic_options.add_options()
      ("help,h", "Print help and usage message.")
      ("verbose,v", "Be more verbose.")
      ("version", "Print version of the program.")
      ("config-file",
           po::value<std::string>(&config_file_name),
           "Parse options from configuration file given.");
    po::options_description config_options("Configuration");
    config_options.add_options()
      // We use int, because program_options parses "-5" with unsigned types.
      ("default-precision,p",
           po::value<int>(&default_precision)->notifier(
             check_range<int, 2, std::numeric_limits<int>::max()>),
           "Set the default precision of MPFR.")
      ("input-file",
           po::value<std::vector<std::string>>()->composing(),
           "The input files to process.");
    po::options_description all_options;
    all_options.add(generic_options).add(config_options);
    po::positional_options_description positional_description;
    positional_description.add("input_file", -1);

    po::variables_map variables_map;
    po::store(po::command_line_parser(argc, argv)
            .options(all_options)
            .positional(positional_description).run(),
        variables_map);
    po::notify(variables_map);
    if (variables_map.count("config-file")) {
      std::ifstream config_file(config_file_name);
      po::store(po::parse_config_file(config_file, all_options),
          variables_map);
      variables_map.notify();
    }

    if (variables_map.count("help")) {
      std::cout << all_options << "\n";
      return 0;
    }
  } catch (std::exception& e) {
    std::cout << e.what() << "\n" << std::flush;
    return 1;
  }
  // End of parsing cmdline flags }}}

  // Do some tests. {{{
  const size_t size = 100;

  Real::default_precision(default_precision);
  std::cout << std::setprecision(Real::default_precision());

  std::vector<Triplet> triplets;
  triplets.reserve(10*size);
  for (size_t ii = 0; ii < 10*size; ++ii) {
    triplets.push_back(Triplet(
          unirand<Index>(size),
          unirand<Index>(size),
          unirand()));
  }
  SparseMatrix A(size, size);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  DenseVector b = DenseVector::Random(size);

  SparseLU solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  DenseVector x = solver.solve(b);
  std::cout << "relative error: " << (A*x - b).norm() / b.norm() << std::endl;
  // End doing some tests. }}}

  return 0;
}

// vim: foldmethod=marker
