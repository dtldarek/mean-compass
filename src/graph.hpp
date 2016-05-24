/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <stdexcept>

namespace mean_compass {

template<typename Config> Graph<Config>::Graph(UTF8Input* input) {
  // Parse the number of vertices and edges.
  input->get(&n_min_);
  input->get(&n_max_);
  input->get(&m_);

  n_ = n_min_ + n_max_;

  // We resize the containers for better memory efficiency.
  label_.resize(n_);       // We use `resize` instead of `reserve`,
  inedges_.resize(n_);     // because we want the basic objects constructed,
  outedges_.resize(n_);    // so that we can pass their addresses to get(...).
  weight_.resize(n_);
  outdegrees_.resize(n_);
  cumulative_outdegrees_.reserve(n_ + 1);  // We start with 0 and use `reserve`
  cumulative_outdegrees_.push_back(0);     // for the convenience of `back`.

  // Read label and vertex-weights.
  // It would be convenient to use 1-based indexes for vertices,
  // but we would like the matrix to be compatible with all the rest,
  // so we will use 0-based indexes.
  for (Index ii = 0; ii < n_; ++ii) {
    input->get(&label_[ii]);
    index_[label_[ii]] = ii;
    input->get(&weight_[ii]);
  }
  // Read edges. Currently there are no edge-weights.
  for (Index ii = 0; ii < m_; ++ii) {
    Index v_1 = index(input->get_string());
    Index v_2 = index(input->get_string());
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
    std::sort(inedges_[ii].begin(), inedges_[ii].end());
    std::sort(outedges_[ii].begin(), outedges_[ii].end());
    m_min_ += outedges_[ii].size();
    Index next_sum = cumulative_outdegrees_.back() + outedges_[ii].size();
    cumulative_outdegrees_.push_back(next_sum);
  }
  // Postprocessing for the max player (edge sorting, outdegrees).
  m_max_ = 0;
  for (Index ii = n_min_; ii < n_; ++ii) {
    std::sort(inedges_[ii].begin(), inedges_[ii].end());
    std::sort(outedges_[ii].begin(), outedges_[ii].end());
    m_max_ += outedges_[ii].size();
    outdegrees_[ii] = outedges_[ii].size();
    Index next_sum = cumulative_outdegrees_.back() + outedges_[ii].size();
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
  update_min_target();
  update_max_target();
}
template<typename Config> void Graph<Config>::init_state(
    const Real& barrier_coef, const Real& mixing_coef) {
  init_flow();
  init_position(barrier_coef, mixing_coef);
  init_targets();
  update_min_target();
  update_max_target();
}
template<typename Config> void Graph<Config>::init_flow() {
  flow_.resize(n_, n_);
  flow_.reserve(outdegrees_);
  for (Index col = 0; col < n_; ++col) {
    Real equal_split(1.0);
    equal_split /= outedges_[col].size();
    for (Index row : outedges_[col]) {
      flow_.insert(row, col) = equal_split;
    }
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
  solver.analyzePattern(mat);
  solver.factorize(mat);
  position_ = solver.solve(vec);
}
template<typename Config> void Graph<Config>::init_targets() {
  std::cout << n_max_ << " " << n_min_ << std::endl;
  std::cout << m_max_ << " " << m_min_ << std::endl;
  min_target_.resize(n_max_ + m_min_);
  max_target_.resize(n_min_ + m_max_);
  update_min_target();
  update_max_target();
}
template<typename Config> void Graph<Config>::update_min_target() {
  // TODO
}
template<typename Config> void Graph<Config>::update_max_target() {
  // TODO
}

template<typename Config> Graph<Config>::MinProblem::MinProblem(
    const Real& barrier_coef, const Real& mixing_coef, Graph* graph) :
    barrier_coef_(barrier_coef), mixing_coef_(mixing_coef), graph_(graph) { }

template<typename Config> void Graph<Config>::MinProblem::init_position(
    Vector* result) const {
  static_cast<void>(result);
}
template<typename Config> void Graph<Config>::MinProblem::value(
    const Vector& position, Real* result) const {
  static_cast<void>(position);
  static_cast<void>(result);
}
template<typename Config> void Graph<Config>::MinProblem::gradient(
    const Vector& position, Vector* result) const {
  static_cast<void>(position);
  static_cast<void>(result);
}
template<typename Config> void Graph<Config>::MinProblem::hessian(
    const Vector& position, Matrix* result) const {
  static_cast<void>(position);
  static_cast<void>(result);
}
template<typename Config> void Graph<Config>::MinProblem::update(
    const Vector& position) {
  static_cast<void>(position);
}

}  // namespace mean_compass

