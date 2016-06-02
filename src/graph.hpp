/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <stdexcept>
#include "graph.h"


namespace mean_compass {

template<typename Config> Graph<Config>::Graph(const Config& config, UTF8Input* input) {
  // Parse the number of vertices and edges.
  input->get(&n_min_);
  input->get(&n_max_);
  input->get(&m_);

  n_ = n_min_ + n_max_;

  int log2 = 0;
  Weight pow2 = 1;
  while (pow2 <= n_) {
    pow2 *= 2;
    log2 += 1;
  }

  // We resize the containers for better memory efficiency.
  label_.resize(n_);       // We use `resize` instead of `reserve`,
  inedges_.resize(n_);     // because we want the basic objects constructed,
  outedges_.resize(n_);    // so that we can pass their addresses to get(...).
  weight_.resize(n_);
  outdegrees_.resize(n_);
  cumulative_outdegrees_.reserve(n_ + 1);  // We start with 0 and use `reserve`
  cumulative_outdegrees_.push_back(0);     // for the convenience of `back`.

  Real min_abs_weight = 0;
  Real max_abs_weight = 0;

  // Read label and vertex-weights.
  // It would be convenient to use 1-based indexes for vertices,
  // but we would like the matrix to be compatible with all the rest,
  // so we will use 0-based indexes.
  for (Index ii = 0; ii < n_; ++ii) {
    input->get(&label_[ii]);
    index_[label_[ii]] = ii;
    if (config.parity()) {  // The weight is actually a parity game priority.
      int priority;
      input->get(&priority);
      // TODO: Depends on the type of Weight.
      weight_[ii] = boost::multiprecision::ldexp(Weight(1.0), priority * log2);
      if (priority & 0x01) {  // The priority is odd.
        weight_[ii] *= -1;  // Odd is min.
      }
    } else {
      input->get(&weight_[ii]);
    }
    const Real& abs_weight = boost::multiprecision::abs(weight_[ii]);
    if (abs_weight > 0 && max_abs_weight == 0) {
      min_abs_weight = abs_weight;
      max_abs_weight = abs_weight;
    } else if (abs_weight > 0 && min_abs_weight > abs_weight) {
      min_abs_weight = abs_weight;
    } else if (abs_weight > 0 && max_abs_weight < abs_weight) {
      max_abs_weight = abs_weight;
    }
  }

  // TODO: is this enough?
  epsilon_ = min_abs_weight / max_abs_weight / n_;

  // Read edges. Currently there are no edge-weights.
  for (Index ii = 0; ii < m_; ++ii) {
    const Index v_1 = index(input->get_string());
    const Index v_2 = index(input->get_string());
    if (v_1 == v_2) {
      std::stringstream description;
      description << "Encountered a loop edge: "
                  << v_1 << " -> " << v_2 << ", "
                  << "currently such edges are not supported";
      throw std::runtime_error(description.str());
    }
    outedges_[v_1].push_back(v_2);
    inedges_[v_2].push_back(v_1);
  }
  // Postprocessing for the min player (edge sorting, outdegrees).
  m_min_ = 0;
  for (Index ii = 0; ii < n_min_; ++ii) {
    assert(outedges_[ii].size() > 0);
    std::sort(inedges_[ii].begin(), inedges_[ii].end());
    std::sort(outedges_[ii].begin(), outedges_[ii].end());
    m_min_ += outedges_[ii].size();
    outdegrees_[ii] = outedges_[ii].size();
    const Index next_sum = cumulative_outdegrees_.back() + outedges_[ii].size();
    cumulative_outdegrees_.push_back(next_sum);
  }
  // Postprocessing for the max player (edge sorting, outdegrees).
  m_max_ = 0;
  for (Index ii = n_min_; ii < n_; ++ii) {
    assert(outedges_[ii].size() > 0);
    std::sort(inedges_[ii].begin(), inedges_[ii].end());
    std::sort(outedges_[ii].begin(), outedges_[ii].end());
    m_max_ += outedges_[ii].size();
    outdegrees_[ii] = outedges_[ii].size();
    const Index next_sum = cumulative_outdegrees_.back() + outedges_[ii].size();
    cumulative_outdegrees_.push_back(next_sum);
  }
}
template<typename Config> void Graph<Config>::init_state(UTF8Input* input) {
  // Initialize flow.
  flow_.resize(n_, n_);
  flow_.reserve(outdegrees_);
  for (Index col = 0; col < n_; ++col) {
    for (Index row : outedges_[col]) {
      flow_.insert(row, col) = input->get_real();
    }
  }
  // Initialize position.
  position_.resize(n_);
  for (Index ii = 0; ii < n_; ++ii) {
    position_(ii) = input->get_real();
  }
  // Non-input-based initializations.
  init_targets();
}
template<typename Config> void Graph<Config>::init_state(
    const Real& barrier_coef, const Real& mixing_coef) {
  init_flow();
  init_position(barrier_coef, mixing_coef);
  init_targets();
}
template<typename Config> void Graph<Config>::init_flow() {
  flow_.resize(n_, n_);
  flow_.reserve(outdegrees_);
  min_flow_log_sum_ = 0;
  max_flow_log_sum_ = 0;
  // We intialize the flow matrix with uniform split for each vertex.
  for (Index col = 0; col < n_; ++col) {
    Real equal_split = Real(1.0) / outdegrees_[col];
    for (Index row : outedges_[col]) {
      flow_.insert(row, col) = equal_split;
    }
    // Update {min|max}_flow_log_sum.
    (col < n_min_? min_flow_log_sum_ : max_flow_log_sum_) +=
      boost::multiprecision::log(equal_split) * outdegrees_[col];
  }
}
template<typename Config> void Graph<Config>::init_position(
    const Real& barrier_coef, const Real& mixing_coef) {
  static_cast<void>(barrier_coef);  // We don't need it now.

  // Mixing coefficient denotes the probability of making a random move.
  // When, `non_mixing_coef` = 1.0 then we strictly follow matrix `flow_`.
  const Real non_mixing_coef = Real(1.0) - mixing_coef;

  using Triplet = typename Config::Triplet;
  std::vector<Triplet> triplets;
  for (Index k = 0; k < flow_.outerSize(); ++k) {
    for (typename Matrix::InnerIterator it(flow_, k); it; ++it) {
      // The last row is the sum of all the variables (see the next loop),
      // so we don't want some other values to interfere.
      if (it.row() != n_ - 1) {
        triplets.push_back(Triplet(it.row(), it.col(), it.value() * non_mixing_coef));
      }
    }
  }
  for (Index ii = 0; ii < n_; ++ii) {
    if (ii != n_ - 1) {
      triplets.push_back(Triplet(ii, ii, -1.0));
    }
    triplets.push_back(Triplet(n_ - 1, ii, 1.0));
  }
  Matrix mat(n_, n_);
  mat.setFromTriplets(triplets.begin(), triplets.end());
  mat.makeCompressed();
  Vector vec = Vector::Constant(n_, mixing_coef / n_ * -non_mixing_coef);
  vec(n_ - 1) = non_mixing_coef;

  typename Config::LU solver;
  solver.compute(mat);
  position_ = solver.solve(vec);
}
template<typename Config> void Graph<Config>::init_targets() {
  min_target_.resize(n_max_ + m_min_);
  for (Index ii = 0; ii < n_max_; ++ii) {
    const Index v_max = ii + n_min_;  // The indexes of max vertices start at n_min_;
    min_target_(ii) = weight_[v_max];
  }
  for (Index ii = 0; ii < n_min_; ++ii) {
    const Index v_min = ii;  // The indexes of min vertices start at 0.
    for (Index jj = 0; jj < outdegrees_[v_min]; ++jj) {
      const Index row = n_max_ + cumulative_outdegrees_[v_min] - 0 + jj;
      min_target_(row) = weight_[v_min];
    }
  }

  // We do not negate max_target here.
  max_target_.resize(n_min_ + m_max_);
  for (Index ii = 0; ii < n_min_; ++ii) {
    const Index v_min = ii;  // The indexes of min vertices start at 0.
    max_target_(ii) = weight_[v_min];
  }
  for (Index ii = 0; ii < n_max_; ++ii) {
    const Index v_max = ii + n_min_;  // The indexes of max vertices start at n_min_;
    for (Index jj = 0; jj < outdegrees_[v_max]; ++jj) {
      const Index row = n_min_ + cumulative_outdegrees_[v_max] - m_min_ + jj;
      max_target_(row) = weight_[v_max];
    }
  }
}

template<typename Config> Graph<Config>::MinProblem::MinProblem(
    const Real& barrier_coef, const Real& mixing_coef, Graph* graph) :
    barrier_coef_(barrier_coef),
    mixing_coef_(mixing_coef),
    non_mixing_coef_(Real(1.0) - mixing_coef_),
    graph_(graph) { }

template<typename Config>
typename Graph<Config>::Vector Graph<Config>::MinProblem::init_position() const {
  typename Config::Vector min_position(graph_->n_max_ + graph_->m_min_);
  for (Index ii = 0; ii < graph_->n_max_; ++ii) {
    const Index v_max = ii + graph_->n_min_;  // The indexes of max vertices start at n_min_;
    min_position(ii) = graph_->position_(v_max);
  }
  for (Index ii = 0; ii < graph_->n_min_; ++ii) {
    const Index v_min = ii;  // The indexes of min vertices start at 0.
    const std::vector<Index>& current_outedges = graph_->outedges_[v_min];
    for (Index jj = 0; jj < graph_->outdegrees_[v_min]; ++jj) {
      const Index v_head = current_outedges[jj];
      const Index row = graph_->n_max_ + graph_->cumulative_outdegrees_[v_min] - 0 + jj;
      min_position(row) = graph_->position_(v_min) * graph_->flow_.coeff(v_head, v_min);
    }
  }
  return min_position;
}
template<typename Config>
typename Graph<Config>::Real Graph<Config>::MinProblem::value(
    const Vector& min_position) const {
  Real result = 0;
  Real barrier = graph_->min_flow_log_sum_;
  for (Index ii = 0; ii < graph_->n_max_; ++ii) {
    const Index v_max = ii + graph_->n_min_;
    result += min_position(ii) * graph_->min_target_(ii);
    barrier += boost::multiprecision::log(min_position(ii)) * graph_->outdegrees_[v_max];
  }
  for (Index ii = graph_->n_max_; ii < graph_->n_max_ + graph_->m_min_; ++ii) {
    result += min_position(ii) * graph_->min_target_(ii);
    barrier += boost::multiprecision::log(min_position(ii));
  }
  result -= barrier * barrier_coef_;
  return result;
}
template<typename Config>
typename Graph<Config>::Vector Graph<Config>::MinProblem::gradient(
    const Vector& min_position) const {
  Vector result = graph_->min_target_;
  for (Index ii = 0; ii < graph_->n_max_; ++ii) {
    const Index v_max = ii + graph_->n_min_;
    result(ii) += -Real(graph_->outdegrees_[v_max]) / min_position(ii) * barrier_coef_;
  }
  for (Index ii = graph_->n_max_; ii < graph_->n_max_ + graph_->m_min_; ++ii) {
    result(ii) += -Real(1) / min_position(ii) * barrier_coef_;
  }
  return result;
}
template<typename Config>
typename Graph<Config>::Diagonal Graph<Config>::MinProblem::hessian(
    const Vector& min_position) const {
  Vector result = Vector::Zero(graph_->n_max_ + graph_->m_min_);
  for (Index ii = 0; ii < graph_->n_max_; ++ii) {
    const Index v_max = ii + graph_->n_min_;
    result(ii) += Real(graph_->outdegrees_[v_max]) /
      min_position(ii) / min_position(ii) * barrier_coef_;
  }
  for (Index ii = graph_->n_max_; ii < graph_->n_max_ + graph_->m_min_; ++ii) {
    result(ii) += Real(1) / min_position(ii) / min_position(ii) * barrier_coef_;
  }
  return result.asDiagonal();
}
template<typename Config>
typename Graph<Config>::Matrix Graph<Config>::MinProblem::equality_matrix(
    ) const {
  /* The basic equality matrix A we are encoding looks equals
   *
   *          <-- v_max -->     <-- e_min -->
   *
   *     / [in these columns | -1 -1 -1             <--- if deg(v_min_0) = 3
   *    /  [ each cell +=    |          -1 -1       <--- if deg(v_min_1) = 2
   * v_min [ flow from v_max |                ...
   *    \  [ to v_row, i.e., |                    -1
   *     \ [ flow(row,col)   |                       -1 -1 -1 -1
   *       [-----------------+----------------------------------
   *     / [ -1              |
   *    /  [    -1           |  in these columns each cell += 1.0
   * v_max [       ...       |  if it represents an edge e_min
   *    \  [           -1    |  that ends at v_row
   *     \ [              -1 |
   *
   * However, we modify it using the mixing coeficient, that is, we use
   *
   *   non_mixing_coef * A - mixing_coef * J
   *
   * Where J is the kind-of identity matrix which ensures the -1's above
   * stay equal to -1.0 (that is, J would be identity if we were to sort it
   * and process vertices only, instead of max-vertices and min-edges).
   *
   * Finally, we substitute the last row, with [1.0, 1.0, 1.0, ..., 1.0], so
   * that it represents the sum of all the variables.
   */

  using Triplet = typename Config::Triplet;
  std::vector<Triplet> triplets;
  // WARNING: implementation dependent optimization.
  // The outer loop iterates column-related _vertices_(not exact columns).
  Index col_min = graph_->n_max_;  // No const.
  for (Index ii = 0; ii < graph_->n_; ++ii) {
    if (ii < graph_->n_max_) {
      const Index v_max = ii + graph_->n_min_;
      const Index col = ii;
      // The inner loop iterates over some rows.
      for (typename Matrix::InnerIterator it(graph_->flow_, v_max); it; ++it) {
        if (it.row() != graph_->n_ - 1) {
          // Out-flow from v_max to v_row.
          triplets.push_back(Triplet(it.row(), col, it.value() * non_mixing_coef_));
        }
      }
      // In-flow of v_max.
      if (v_max != graph_->n_ - 1) {
        triplets.push_back(Triplet(v_max, col, -1.0));
      }
      // Add the 1.0 in the last row.
      triplets.push_back(Triplet(graph_->n_ - 1, col, 1.0));
    } else {
      const Index v_min = ii - graph_->n_max_;
      // Here we iterate over out-edges, which means
      // we know both the row and the colum.
      // For each column we have at most 3 values:
      //   * the out-flow of current edge,
      //   * the in-flow of the tail vertex,
      //   * the 1.0 in the last row.
      for (Index v_row : graph_->outedges_[v_min]) {
        // Here `col_min` is actually the column.
        if (v_row != graph_->n_ - 1) {
          // Out-flow from v_min via e_min to v_row.
          triplets.push_back(Triplet(v_row, col_min, non_mixing_coef_));
        }
        // In-flow of v_min, i.e., the sum of all related e_min.
        triplets.push_back(Triplet(v_min, col_min, -1.0));
        // The last row value.
        triplets.push_back(Triplet(graph_->n_ - 1, col_min, 1.0));
        col_min++;
      }
    }
  }
  Matrix result(graph_->n_, graph_->n_max_ + graph_->m_min_);
  result.setFromTriplets(triplets.begin(), triplets.end());
  result.makeCompressed();
  return result;
}
template<typename Config>
typename Graph<Config>::Vector Graph<Config>::MinProblem::equality_vector(
    ) const {
  /* The basic equality vector would be equal to [0, 0, 0, ..., 0, 0, 1],
   * but because of the mixing_coeficient, we get that:
   *   * the non-last-row value is equal to non_mixing_coef * mixing_coef / n,
   *   * the last-row value is equal to non_mixing_coef.
   */
  Vector result = Vector::Constant(graph_->n_, mixing_coef_ / graph_->n_ * -non_mixing_coef_);
  result(graph_->n_ - 1) = non_mixing_coef_;
  return result;
}
template<typename Config> void Graph<Config>::MinProblem::update(
    const Vector& min_position) {
  for (Index ii = 0; ii < graph_->n_max_; ++ii) {
    const Index v_max = ii + graph_->n_min_;
    graph_->position_(v_max) = min_position(ii);
  }
  for (Index ii = 0; ii < graph_->n_min_; ++ii) {
    const Index v_min = ii;  // The indexes of min vertices start at 0.
    const std::vector<Index>& current_outedges = graph_->outedges_[v_min];
    graph_->position_(v_min) = 0;
    for (Index jj = 0; jj < graph_->outdegrees_[v_min]; ++jj) {
      const Index row = graph_->n_max_ + graph_->cumulative_outdegrees_[v_min] - 0 + jj;
      graph_->position_(v_min) += min_position(row);
    }
    for (Index jj = 0; jj < graph_->outdegrees_[v_min]; ++jj) {
      const Index v_head = current_outedges[jj];
      const Index row = graph_->n_max_ + graph_->cumulative_outdegrees_[v_min] - 0 + jj;
      graph_->flow_.coeffRef(v_head, v_min) = min_position(row) / graph_->position_(v_min);
    }
  }
}

// There is a lot of code share between MinProblem and MaxProblem, but
// there are tons of small changes, which makes it (at least now) infeasible.
template<typename Config> Graph<Config>::MaxProblem::MaxProblem(
    const Real& barrier_coef, const Real& mixing_coef, Graph* graph) :
    barrier_coef_(barrier_coef),
    mixing_coef_(mixing_coef),
    non_mixing_coef_(Real(1.0) - mixing_coef_),
    graph_(graph) { }

template<typename Config>
typename Graph<Config>::Vector Graph<Config>::MaxProblem::init_position() const {
  typename Config::Vector max_position(graph_->n_min_ + graph_->m_max_);
  for (Index ii = 0; ii < graph_->n_min_; ++ii) {
    const Index v_min = ii;  // The indexes of min vertices start at 0.
    max_position(ii) = graph_->position_(v_min);
  }
  for (Index ii = 0; ii < graph_->n_max_; ++ii) {
    const Index v_max = ii + graph_->n_min_;  // The indexes of max vertices start at n_min_.
    const std::vector<Index>& current_outedges = graph_->outedges_[v_max];
    for (Index jj = 0; jj < graph_->outdegrees_[v_max]; ++jj) {
      const Index v_head = current_outedges[jj];
      const Index row = graph_->n_min_ + graph_->cumulative_outdegrees_[v_max]
                      - graph_->cumulative_outdegrees_[graph_->n_min_] + jj;
      max_position(row) = graph_->position_(v_max) * graph_->flow_.coeff(v_head, v_max);
    }
  }
  return max_position;
}
template<typename Config>
typename Graph<Config>::Real Graph<Config>::MaxProblem::value(
    const Vector& max_position) const {
  Real result = 0;
  Real barrier = graph_->max_flow_log_sum_;
  for (Index ii = 0; ii < graph_->n_min_; ++ii) {
    const Index v_min = ii;
    result -= max_position(ii) * graph_->max_target_(ii);
    barrier += boost::multiprecision::log(max_position(ii)) * graph_->outdegrees_[v_min];
  }
  for (Index ii = graph_->n_min_; ii < graph_->n_min_ + graph_->m_max_; ++ii) {
    result -= max_position(ii) * graph_->max_target_(ii);
    barrier += boost::multiprecision::log(max_position(ii));
  }
  result -= barrier * barrier_coef_;
  return result;
}
template<typename Config>
typename Graph<Config>::Vector Graph<Config>::MaxProblem::gradient(
    const Vector& max_position) const {
  Vector result = -graph_->max_target_;
  for (Index ii = 0; ii < graph_->n_min_; ++ii) {
    const Index v_min = ii + graph_->n_max_;
    result(ii) += -Real(graph_->outdegrees_[v_min]) / max_position(ii) * barrier_coef_;
  }
  for (Index ii = graph_->n_min_; ii < graph_->n_min_ + graph_->m_max_; ++ii) {
    result(ii) += -Real(1) / max_position(ii) * barrier_coef_;
  }
  return result;
}
template<typename Config>
typename Graph<Config>::Diagonal Graph<Config>::MaxProblem::hessian(
    const Vector& max_position) const {
  Vector result = Vector::Zero(graph_->n_min_ + graph_->m_max_);
  for (Index ii = 0; ii < graph_->n_min_; ++ii) {
    const Index v_min = ii + graph_->n_max_;
    result(ii) += Real(graph_->outdegrees_[v_min]) /
      max_position(ii) / max_position(ii) * barrier_coef_;
  }
  for (Index ii = graph_->n_min_; ii < graph_->n_min_ + graph_->m_max_; ++ii) {
    result(ii) += Real(1) / max_position(ii) / max_position(ii) * barrier_coef_;
  }
  return result.asDiagonal();
}
template<typename Config>
typename Graph<Config>::Matrix Graph<Config>::MaxProblem::equality_matrix(
    ) const {
  /* The basic equality matrix A looks similar to the one in
   * MinProblem::equality_matrix (notice the reordering of columns
   * to match max_position):
   *
   *          <-- v_min -->     <-- e_max -->
   *
   *     / [ -1              |
   *    /  [    -1           |  in these columns each cell += 1.0
   * v_min [       ...       |  if it represents an edge e_max
   *    \  [           -1    |  that ends at v_row
   *     \ [              -1 |
   *       [-----------------+----------------------------------
   *     / [in these columns | -1 -1 -1             <--- if deg(v_max_0) = 3
   *    /  [ each cell +=    |          -1 -1       <--- if deg(v_max_1) = 2
   * v_max [ flow from v_min |                ...
   *    \  [ to v_row, i.e., |                    -1
   *     \ [ flow(row,col)   |                       -1 -1 -1 -1
   *
   */

  using Triplet = typename Config::Triplet;
  std::vector<Triplet> triplets;
  // WARNING: implementation dependent optimization.
  // The outer loop iterates column-related _vertices_(not exact columns).
  Index col_max = graph_->n_min_;  // No const.
  for (Index ii = 0; ii < graph_->n_; ++ii) {
    if (ii < graph_->n_min_) {
      const Index v_min = ii;
      const Index col = ii;
      // The inner loop iterates over some rows.
      for (typename Matrix::InnerIterator it(graph_->flow_, v_min); it; ++it) {
        if (it.row() != graph_->n_ - 1) {
          // Out-flow from v_min to v_row.
          triplets.push_back(Triplet(it.row(), col, it.value() * non_mixing_coef_));
        }
      }
      // In-flow of v_min.
      if (v_min != graph_->n_ - 1) {
        triplets.push_back(Triplet(v_min, col, -1.0));
      }
      // Add the 1.0 in the last row.
      triplets.push_back(Triplet(graph_->n_ - 1, col, 1.0));
    } else {
      const Index v_max = ii;
      // Here we iterate over out-edges, which means
      // we know both the row and the colum.
      // For each column we have at most 3 values:
      //   * the out-flow of current edge,
      //   * the in-flow of the tail vertex,
      //   * the 1.0 in the last row.
      for (Index v_row : graph_->outedges_[v_max]) {
        // Here `col_max` is actually the column.
        if (v_row != graph_->n_ - 1) {
          // Out-flow from v_max via e_max to v_row.
          triplets.push_back(Triplet(v_row, col_max, non_mixing_coef_));
        }
        if (v_max != graph_->n_ - 1) {
          // In-flow of v_max, i.e., the sum of all related e_max.
          triplets.push_back(Triplet(v_max, col_max, -1.0));
        }
        // The last row value.
        triplets.push_back(Triplet(graph_->n_ - 1, col_max, 1.0));
        col_max++;
      }
    }
  }
  Matrix result(graph_->n_, graph_->n_min_ + graph_->m_max_);
  result.setFromTriplets(triplets.begin(), triplets.end());
  result.makeCompressed();
  return result;
}
template<typename Config>
typename Graph<Config>::Vector Graph<Config>::MaxProblem::equality_vector(
    ) const {
  // See the comment at MinProblem::equality_vector.
  Vector result = Vector::Constant(graph_->n_, mixing_coef_ / graph_->n_ * -non_mixing_coef_);
  result(graph_->n_ - 1) = non_mixing_coef_;
  return result;
}
template<typename Config> void Graph<Config>::MaxProblem::update(
    const Vector& max_position) {
  for (Index ii = 0; ii < graph_->n_min_; ++ii) {
    const Index v_min = ii;
    graph_->position_(v_min) = max_position(ii);
  }
  for (Index ii = 0; ii < graph_->n_max_; ++ii) {
    const Index v_max = ii + graph_->n_min_;
    const std::vector<Index>& current_outedges = graph_->outedges_[v_max];
    graph_->position_(v_max) = 0;
    for (Index jj = 0; jj < graph_->outdegrees_[v_max]; ++jj) {
      const Index row = graph_->n_min_ + graph_->cumulative_outdegrees_[v_max] - graph_->m_min_ + jj;
      graph_->position_(v_max) += max_position(row);
    }
    for (Index jj = 0; jj < graph_->outdegrees_[v_max]; ++jj) {
      const Index v_head = current_outedges[jj];
      const Index row = graph_->n_min_ + graph_->cumulative_outdegrees_[v_max] - graph_->m_min_ + jj;
      graph_->flow_.coeffRef(v_head, v_max) = max_position(row) / graph_->position_(v_max);
    }
  }
}


}  // namespace mean_compass

// vim: et sw=2 ts=2 foldmethod=syntax
