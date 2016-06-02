/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef __MEAN_COMPASS_GRAPH_H__
#define __MEAN_COMPASS_GRAPH_H__

#include <unordered_map>
#include <vector>
#include "types.h"
#include "utf8_io.h"

namespace mean_compass {

template<typename Config>
class Graph {
 public:
  using Index    = typename Config::Index;
  using Real     = typename Config::Real;
  using Matrix   = typename Config::Matrix;
  using Vector   = typename Config::Vector;
  using Weight   = typename Config::Weight;
  using Diagonal = typename Config::Diagonal;

  class MinProblem {
   public:
    Vector init_position() const;
    Real value(const Vector& min_position) const;
    Vector gradient(const Vector& min_position) const;
    Diagonal hessian(const Vector& min_position) const;
    Vector equality_vector() const;
    Matrix equality_matrix() const;
    void update(const Vector& min_position);

    Real epsilon() const { return graph_->epsilon_; }

   protected:
    MinProblem(const Real& barrier_coef,
               const Real& mixing_coef,
               Graph* graph);
    Real barrier_coef_;
    Real mixing_coef_;
    Real non_mixing_coef_;
    Graph* graph_;
    friend Graph;
  };
  // A lot of code between MinProblem and MaxProblem could have been shared,
  // but there are many small changes and I couldn't find a good way to do it.
  class MaxProblem {
   public:
    Vector init_position() const;
    Real value(const Vector& max_position) const;
    Vector gradient(const Vector& max_position) const;
    Diagonal hessian(const Vector& max_position) const;
    Vector equality_vector() const;
    Matrix equality_matrix() const;
    void update(const Vector& max_position);

    Real epsilon() const { return graph_->epsilon_; }

   protected:
    MaxProblem(const Real& barrier_coef,
               const Real& mixing_coef,
               Graph* graph);
    Real barrier_coef_;
    Real mixing_coef_;
    Real non_mixing_coef_;
    Graph* graph_;
    friend Graph;
  };

  /* Parse a graph from the following input format:
   *
   *     # This is a comment, anything after '#' is ignored.
   *     # Empty lines are ingored too.
   *     n1 n2 m             # |min_vertices| |max_vertices| |edges|
   *     label_u1 weight_u1  # label for the first min-vertex and its weight
   *     label_u2 weight_u2
   *     ...                 # continue for (n1 - 2) more lines
   *     label_v1 weight_v1  # label for the first max-vertex and its weight
   *     ...                 # continue for (n2 - 1) more lines
   *     label1 label2       # an edge between any two different vertices
   *     label3 label4
   *     ...                 # continue for (m - 2) more lines
   *
   *
   * To avoid long live of the UTF8Input object, one can call it like this:
   *
   *     Graph my_graph(UTF8Input(filename));
   *
   * (note: temporary objects are destroyed at the end of the full expression).
   */
  explicit Graph(UTF8Input&& input) : Graph(&input) { }
  explicit Graph(UTF8Input* input);

  // Initialize the flow matrix, position and targets.
  void init_state(UTF8Input* input);
  void init_state(const Real& barrier_coef, const Real& mixing_coef);

  Index n()       const { return n_; }
  Index n_min()   const { return n_min_; }
  Index n_max()   const { return n_max_; }
  Index m()       const { return m_; }
  Index m_min()   const { return m_min_; }
  Index m_max()   const { return m_max_; }
  Real  epsilon() const { return epsilon_; }
  Index index(const std::string& label) const { return index_.at(label); }

  const std::string& label(Index index) const { return label_[index]; }
  const Weight& weight(Index index) const { return weight_[index]; }
  const std::vector<Index>& inedges(Index index) const { return inedges_[index]; }
  const std::vector<Index>& outedges(Index index) const { return outedges_[index]; }

  const Matrix& flow() const { return flow_; }
  Matrix& flow() { return flow_; }
  const Vector& position() const { return position_; }
  Vector& position() { return position_; }

  MinProblem get_min_problem(const Real& barrier_coef, const Real& mixing_coef) {
    return MinProblem(barrier_coef, mixing_coef, this);
  }
  MaxProblem get_max_problem(const Real& barrier_coef, const Real& mixing_coef) {
    return MaxProblem(barrier_coef, mixing_coef, this);
  }

 protected:
  // Number of vertices(n) and edges(m), controlled by player min and max.
  Index n_, n_min_, n_max_;
  Index m_, m_min_, m_max_;
  Real epsilon_;

  std::unordered_map<std::string, Index>  index_;
  std::vector<std::string>                label_;
  std::vector<std::vector<Index>>         inedges_;
  std::vector<std::vector<Index>>         outedges_;
  std::vector<Index>                      outdegrees_;
  std::vector<Index>                      cumulative_outdegrees_;
  std::vector<Weight>                     weight_;

  // The matrix representing the split of flow.
  // The cell flow[row,col] corresponds to v_col -> v_row, (note the reverse),
  // in particular the matrix is row-stochastic.
  //
  // Non-zero value means an edge in the graph (this value may be _very_ small,
  // but has to be positive). If the defaults were not overriden, then the
  // Matrix type stores value in column-major order.
  Matrix flow_;

  // Roughly min_flow_log_sum = sum_{v_min} sum_{neigh} log(v_min -> neigh).
  // More precisely, suppose that a + b + c = 1 describe a flow split at some
  // vertex, then
  // log(flow * a) + log(flow * b) + log(flow * c) =
  //   3 * log(flow) + log(a) + log(b) + log(c).
  // The component log(a)+... does not change, so we
  // recalculate it don't have to recalculate it every time.
  Real min_flow_log_sum_;
  Real max_flow_log_sum_;

  // The i-th cell describes the probability of being at vertex v_i.
  Vector position_;

  // The components of min_position and min_target are as follows:
  //  * n_max cells corresponding to max vertices,
  //  * m_min cells corresponding to min edges.
  // For max_position and max_target we start with n_min and then m_max.
  Vector min_target_;
  Vector max_target_;

  void init_flow();
  void init_position(const Real& barrier_coef, const Real& mixing_coef);
  void init_targets();

  friend MinProblem;
  friend MaxProblem;
};

}  // namespace mean_compass

// Include definitions.
#include "graph.hpp"

#endif  // __MEAN_COMPASS_GRAPH_H__

// vim: et sw=2 ts=2
