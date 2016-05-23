/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

namespace mean_compass {

template<typename Config> Graph<Config>::Graph(UTF8Input* input) {
  input->get(&n_min_);
  input->get(&n_max_);
  input->get(&m_);

  n_ = n_min_ + n_max_;

  label_.resize(n_);
  inedges_.resize(n_);
  outedges_.resize(n_);
  weight_.resize(n_);
  outdegrees_.resize(n_);
  cumulative_outdegrees_.reserve(n_ + 1);  // We start with 0.
  cumulative_outdegrees_.push_back(0);

  // It would be convinient to use 1-based indexes for vertices,
  // but I would like the matrix to be compatible with all the rest,
  // so we will use 0-based indexes.
  for (Index ii = 0; ii < n_; ++ii) {
    input->get(&label_[ii]);
    index_[label_[ii]] = ii;
    input->get(&weight_[ii]);
  }
  for (Index ii = 0; ii < m_; ++ii) {
    Index v_1 = index(input->get_string());
    Index v_2 = index(input->get_string());
    outedges_[v_1].push_back(v_2);
    inedges_[v_2].push_back(v_1);
  }
  for (Index ii = 0, m_min = 0; ii < n_min_; ++ii) {
    std::sort(inedges_[ii].begin(), inedges_[ii].end());
    std::sort(outedges_[ii].begin(), outedges_[ii].end());
    m_min += outedges_[ii].size();
    Index next_sum = cumulative_outdegrees_.back() + outedges_[ii].size();
    cumulative_outdegrees_.push_back(next_sum);
  }
  for (Index ii = n_min_, m_max = 0; ii < n_; ++ii) {
    std::sort(inedges_[ii].begin(), inedges_[ii].end());
    std::sort(outedges_[ii].begin(), outedges_[ii].end());
    m_max += outedges_[ii].size();
    outdegrees_[ii] = outedges_[ii].size();
    Index next_sum = cumulative_outdegrees_.back() + outedges_[ii].size();
    cumulative_outdegrees_.push_back(next_sum);
  }

  init_flow();
  init_position();
  update_min_target();
  update_max_target();
}
template<typename Config> void Graph<Config>::init_flow() {
  flow_.resize(n_, n_);
  std::vector<Index> memory_requirements(n_, 0);
  for (Index ii = 0; ii < n_; ++ii) {
    memory_requirements[ii] = outedges_[ii].size();
  }
  flow_.reserve(memory_requirements);
  for (Index col = 0; col < n_; ++col) {
    Real equal_split(1.0);
    equal_split /= outedges_[col].size();
    for (Index row : outedges_[col]) {
      flow_.insert(row, col) = equal_split;
    }
  }
}
template<typename Config> void Graph<Config>::init_position() {
  typename Config::LU solver;
  Matrix mat(n_, n_);
  Vector vec(n_);
  vec(n_ - 1) = 1.0;
  //solver.analyzePattern(mat);
  //solver.factorize(mat);
  //position_ = solver.solve(b, &x);
}
template<typename Config> void Graph<Config>::init_targets() {
  min_target_.resize(n_max_ + m_min_);
  max_target_.resize(n_min_ + m_max_);
  update_min_target();
  update_max_target();
}
template<typename Config> void Graph<Config>::update_min_target() {
}
template<typename Config> void Graph<Config>::update_max_target() {
}

template<typename Config> Graph<Config>::MinProblem::MinProblem(
    const Real& barrier_coef, const Real& randdst_coef, Graph* graph) :
    barrier_coef_(barrier_coef), randdst_coef_(randdst_coef), graph_(graph) { }

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

