#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/random.hpp>
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
typedef Eigen::SparseVector<Real, 0, Index>                  SparseVector; // 0 means no flags.

typedef Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<Index>> SparseLU;

Real unirand() {
  static boost::random::mt19937 generator;
  boost::random::uniform_int_distribution<Integer> uniform_integer(0, Integer(1) << Real::default_precision());
  Integer nominator = uniform_integer(generator);
  Integer denominator = Integer(1) << Real::default_precision();
  Rational result(nominator, denominator);
  return static_cast<Real>(result);
}

} // anonymous namespace


int main() {
  const size_t size = 1000;

  Real::default_precision(256);
  //std::cout << std::setprecision(Real::default_precision());

  std::vector<Triplet> triplets;
  triplets.reserve(10*size);
  for (size_t ii = 0; ii < 10*size; ++ii) {
    triplets.push_back(Triplet(rand()%size, rand()%size, unirand()));
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

  return 0;
}
