/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef __MEAN_COMPASS_NEWTON_H__
#define __MEAN_COMPASS_NEWTON_H__

#include <boost/multiprecision/number.hpp>
#include "types.h"

namespace mean_compass {


/* SimpleNewton provides a simple **minimization** routine for
 * convex functions with equality constraints based on the Newton step
 * with infeasible start and elimination.
 *
 * Suppose we like to minimize f(x) subject to Ax = b, where f is convex.
 * The Newton step with infeasible start is defined as:
 *
 *   [ D^2 f(x) | A^T ]  [ step ]   [ D f(x) ]
 *   [----------+-----]  [------] = [--------]                (eq:step)
 *   [    A     |  0  ]  [ dual ]   [ Ax - b ]
 *
 * and then the new x is defined as
 *
 *   x <- x + t*step
 *
 * where parameter t is obtained usually via backtracking line search.
 * Instead of solving equation (eq:step) directly we use substitution.
 * For simplicity we rewrite the above as two equations where
 * H = D^2 f(x), h = D f(x) and g = Ax - b:
 *
 *   H step + A^T dual = -g,   A step = -h
 *
 * and solve for step:
 *
 *   step = - H^(-1) (g + A^T dual)
 *   A H^(-1) (g + A^T dual) = h
 *
 * thus letting S = - A H^(-1) A^T we get
 *
 *   S dual = A H^(-1) g - h
 *   H step = -A^T dual - g
 *
 * This is particularly useful, because in our setting (interior point method)
 * matrix H is diagonal, thus H^(-1) is easy to calculate and S is sparse.
 */
template<typename Config, typename Problem>
class SimpleNewton {
 public:
  using Index = typename Config::Index;
  using Real = typename Config::Real;
  using Vector = typename Config::Vector;
  using Matrix = typename Config::Matrix;

  SimpleNewton(Problem* problem) :
      problem_(problem),
      position_(problem_->init_position()),
      equality_matrix_(problem_->equality_matrix()),
      equality_vector_(problem_->equality_vector()) {
  }

  const Problem& problem() const { return *problem_; }
  Problem& problem() { return *problem_; }
  const Vector& position() const { return position_; }
  Vector& position() { return position_; }

  void init() {
    position_ = problem_->init_position();
    equality_matrix_ = problem_->equality_matrix();
    equality_vector_ = problem_->equality_vector();
  }

  //
  Real step(Vector* dual) {
    typename Config::Diagonal hessian = problem_->hessian(position_);
    typename Config::Diagonal hessian_inv = hessian.inverse();
    Vector gradient = problem_->gradient(position_);
    Matrix mid_mat = equality_matrix_ * hessian_inv;
    Matrix schur_complement = -mid_mat * equality_matrix_.transpose();
    Vector infeasibility = equality_matrix_ * position_ - equality_vector_;
    typename Config::LU solver;
    solver.compute(schur_complement);
    Vector new_dual = solver.solve(mid_mat * gradient - infeasibility);
    Vector step = hessian_inv * (equality_matrix_.transpose() * new_dual + gradient) * (-1);
    Vector dual_step = new_dual - *dual;
    Real alt = 0.01;
    while (true) {
      Vector new_position = position_ + step;
      new_dual = *dual + dual_step;
      Vector blabla = equality_matrix_ * new_position - equality_vector_;
      if (new_position.minCoeff() > 0.0 &&
        residual(new_position, new_dual) <= (Real(1) - alt) * residual(position_, *dual)) {
        position_ = new_position;
        *dual = new_dual;
        return residual(position_, *dual);
      } else {
        // Normally one prefers *= over /=, but /=2 (note 2 instead of 2.0)
        // has much better performance than *= 0.5 or *= Real(0.5).
        step /= 2;
        dual_step /= 2;
        alt /= 2;
      }
    }
  }

  Real residual2(const Vector& primal, const Vector& dual) const {
    Vector r1 = problem_->gradient(primal) + equality_matrix_.transpose() * dual;
    Vector r2 = equality_matrix_ * position_ - equality_vector_;
    Real result = 0;
    for (Index ii = 0; ii < r1.size(); ++ii) {
      const Real& value = r1(ii);
      result += value * value;
    }
    for (Index ii = 0; ii < r2.size(); ++ii) {
      const Real& value = r2(ii);
      result += value * value;
    }
    return result;
  }
  Real residual(const Vector& primal, const Vector& dual) const {
    return boost::multiprecision::sqrt(residual2(primal, dual));
  }

  void loop() {
    Vector dual = Vector::Constant(equality_matrix_.rows(), 1);
    Real yet = 0;
    do {
      yet = step(&dual);
    } while (yet > problem_->epsilon());
  }

 protected:
  Problem* problem_;
  Vector   position_;
  Matrix   equality_matrix_;
  Vector   equality_vector_;
};

}  // namespace mean_compass

#endif  // __MEAN_COMPASS_NEWTON_H__

// vim: et sw=2 ts=2
