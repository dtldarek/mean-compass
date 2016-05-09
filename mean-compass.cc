#include <iostream>
#include <boost/multiprecision/mpfr.hpp>
#include <Eigen/LU>

typedef boost::multiprecision::mpfr_float                     real_t;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1>              VectorX;

int main() {
  real_t::default_precision(256);
  MatrixX A = MatrixX::Random(100, 100);
  VectorX b = VectorX::Random(100);
  // Solve Ax=b using LU.
  VectorX x = A.lu().solve(b);
  std::cout << "relative error: " << (A*x - b).norm() / b.norm() << std::endl;

  return 0;
}
