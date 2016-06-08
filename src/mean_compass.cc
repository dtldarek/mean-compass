/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "mean_compass.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include "graph.h"
#include "newton.h"
#include "types.h"
#include "utils.h"
#include "utf8_io.h"

namespace {

// void check_for_sigint(config, graph, barrier_coef, mixing_coef) {{{
template<typename Config> inline void check_for_sigint(
    const Config& config,
    const mean_compass::Graph<Config>& graph,
    const typename Config::Real& barrier_coef,
    const typename Config::Real& mixing_coef) {
  using namespace mean_compass;
  using Index = typename Config::Index;
  using Matrix = typename Config::Matrix;
  if (utils::sigint_caught) {
    if (!config.dump_file().empty()) {
      std::ofstream file;
      if (config.overwrite_dump_file()) {
        file.open(config.dump_file(),
                  std::ofstream::out | std::ofstream::trunc);
      } else {
        file.open(config.dump_file(),
                  std::ofstream::out | std::ofstream::app);
      }
      file << static_cast<std::string>(barrier_coef) << ' '
           << static_cast<std::string>(mixing_coef) << '\n';
      for (Index k = 0; k < graph.flow().outerSize(); ++k) {
        for (typename Matrix::InnerIterator it(graph.flow(), k); it; ++it) {
          file << it.row() << ' '
               << it.col() << ' '
               << static_cast<std::string>(it.value()) << '\n';
        }
      }
      file << static_cast<std::string>(graph.position()(0));
      for (Index ii = 1; ii < graph.n(); ++ii) {
        file << ' ' << static_cast<std::string>(graph.position()(ii));
      }
      file << "\n\n";
    }
    throw utils::SIGINTException();
  }
}  // }}} End of check_for_sigint.


template<typename Config> inline void handle_graph(
    const Config& config,
    mean_compass::Graph<Config>&& graph_) {
  static_cast<void>(config);
  using namespace mean_compass;
  using Index   = typename Config::Index;
  using Real    = typename Config::Real;
  using Vector  = typename Config::Vector;

  // Setup {{{
  Real mixing_coef = 0.1;
  Real barrier_coef = 1.0;
  Real barrier_multiplier = 0.5;
  if (!config.barrier_multiplier_str().empty()) {
    barrier_multiplier = static_cast<Real>(config.barrier_multiplier_str());
  }

  Graph<Config> graph(graph_);
  if (config.restore_file().empty()) {
    graph.init_state(barrier_coef, mixing_coef);
  } else {
    UTF8Input state_data(config.restore_file());
    barrier_coef = state_data.get_real();
    mixing_coef = state_data.get_real();
    if (state_data) graph.init_state(&state_data);
  }
  // }}} End of setup.

  // Solve the graph. {{{
  Vector old_position = Vector::Constant(graph.n(), 0);
  Vector mid_position = Vector::Constant(graph.n(), 0);
  Vector min_dual = Vector::Constant(graph.n(), 1);
  Vector max_dual = Vector::Constant(graph.n(), 1);
  for (; barrier_coef >= graph.epsilon();
         barrier_coef *= barrier_multiplier,
         mixing_coef *= barrier_multiplier) {
    std::cout << "barrier: " << barrier_coef << (Config::verbose ? '\n' : ' ');
    int small_steps = 0;
    do {
      check_for_sigint(config, graph, barrier_coef, mixing_coef);
      old_position = graph.position();
      typename Graph<Config>::MinProblem min_problem =
          graph.get_min_problem(barrier_coef, mixing_coef);
      SimpleNewton<Config, typename Graph<Config>::MinProblem>
          min_newton(&min_problem);
      if (small_steps < 3) {
        min_newton.step_with_backtracking_line_search(&min_dual);
        min_newton.step_with_backtracking_line_search(&min_dual);
        min_newton.step_with_backtracking_line_search(&min_dual);
      } else {
        min_newton.step_with_exact_line_search(&min_dual);
      }
      min_problem.update(min_newton.position());
      if (Config::verbose) {  // Logging {{{
        //std::cout << "min: "
        //          << std::fixed << graph.position().transpose()
        //          << std::scientific << '\n';
        Vector diff = graph.position() - old_position;
        std::cout << "min: ";
        for (Index ii = 0; ii < graph.n(); ++ii) {
          char c = diff(ii) > graph.epsilon() / graph.n() ? '+' : diff(ii) < -graph.epsilon() / graph.n() ? '-' : ' ';
          std::cout << c;
        }
        std::cout << '\n';
      }  // }}} End of logging.

      if (Config::verbose) {
        mid_position = graph.position();
      }
      typename Graph<Config>::MaxProblem max_problem =
          graph.get_max_problem(barrier_coef, mixing_coef);
      SimpleNewton<Config, typename Graph<Config>::MaxProblem>
          max_newton(&max_problem);
      if (small_steps < 3) {
        min_newton.step_with_backtracking_line_search(&min_dual);
        min_newton.step_with_backtracking_line_search(&min_dual);
        min_newton.step_with_backtracking_line_search(&min_dual);
      } else {
        min_newton.step_with_exact_line_search(&min_dual);
      }
      max_problem.update(max_newton.position());
      if (Config::verbose) {  // Logging {{{
        //std::cout << "max: "
        //          << std::fixed << graph.position().transpose()
        //          << std::scientific << '\n';
        Vector diff = graph.position() - mid_position;
        std::cout << "max: ";
        for (Index ii = 0; ii < graph.n(); ++ii) {
          char c = diff(ii) > graph.epsilon() / graph.n() ? '+' : diff(ii) < -graph.epsilon() / graph.n() ? '-' : ' ';
          std::cout << c;
        }
        std::cout << '\n';
      } else {
        std::cout << '.';
      }  // End of logging. }}}
      small_steps++;
      // TODO: Move barrier adjustment to another function or class.
      if (config.barrier_adjustment()) {
        if (small_steps > 10) {
          barrier_coef /= barrier_multiplier;
          mixing_coef /= barrier_multiplier;
          Real barrier_step = Real(1) / (Real(1) - barrier_multiplier) + Real(1);
          barrier_multiplier = Real(1) - Real(1) / barrier_step;
          std::cout << " new step: " << barrier_step << ' ';
          break;
        }
      }
    } while ((graph.position() - old_position).norm() > graph.epsilon());
    if (config.barrier_adjustment() && small_steps <= 2 && barrier_multiplier > 0.5) {
      Real barrier_step = Real(1) / (Real(1) - barrier_multiplier) - Real(1);
      if (barrier_step < 2) {
        barrier_step = 2;
      }
      barrier_multiplier = Real(1) - Real(1) / barrier_step;
      std::cout << " new step: " << barrier_step << ' ';
    }
    std::cout << '\n';
  }
  // }}} The graph-solving has finished.

  // Discretize strategy and calculate the sign of infinite games. {{{
  std::vector<Index> v_positive;
  std::vector<Index> v_negative;
  std::vector<Index> strategy(graph.n(), graph.n()+1);
  std::vector<typename Graph<Config>::Weight>  score(graph.n(), 0);
  std::vector<Index> visited(graph.n(), 0);
  std::vector<Index> path;
  path.reserve(graph.n()+2);
  for (Index ii = 0; ii < graph.n(); ++ii) {
    strategy[ii] = graph.outedges(ii)[0];
    for (Index outedge : graph.outedges(ii)) {
      if (graph.flow().coeff(outedge, ii) > graph.flow().coeff(strategy[ii], ii)) {
        strategy[ii] = outedge;
      }
    }
  }
  for (Index ii = 0; ii < graph.n(); ++ii) {
    if (visited[ii] != 0) continue;
    path.push_back(ii);
    while (visited[path.back()] == 0) {
      visited[path.back()] = 1;
      path.push_back(strategy[path.back()]);
    }
    if (visited[path.back()] == 1) {  // The score wasn't calculated yet.
      score[path.back()] = graph.weight(path.back());
      for (Index jj = path.size()-2; path[jj] != path.back(); --jj) {
        score[path.back()] += graph.weight(path[jj]);
      }
      assert(score[path.back()] != 0);
      visited[path.back()] = 2;  // This vertex will be pushed back next loop.
    }
    if (visited[path.back()] == 2) {  // Has a score.
      for (Index jj : path) {
        visited[jj] = 2;
        score[jj] = score[path.back()];
      }
      if (score[path.back()] > 0) {
        path.pop_back();
        v_positive.insert(v_positive.end(), path.begin(), path.end());
      } else {
        path.pop_back();
        v_negative.insert(v_negative.end(), path.begin(), path.end());
      }
    }
    path.clear();
  }  // }}} End of discretizing the strategy.

  // Print the results. {{{
  std::sort(v_positive.begin(), v_positive.end());
  std::sort(v_negative.begin(), v_negative.end());
  if (v_negative.size() > 0) {
    std::cout << "Player " << utils::AnsiColors<Config>::YELLOW
              << "min" << utils::AnsiColors<Config>::ENDC
              << " wins from nodes:\n"
              << "  {" << graph.label(v_negative[0]);
    for (Index ii = 1; ii < static_cast<Index>(v_negative.size()); ++ii) {
      std::cout << ", " << graph.label(v_negative[ii]);
    }
    std::cout << "}\n";
  } else {
    std::cout << "Player " << utils::AnsiColors<Config>::YELLOW
              << "min" << utils::AnsiColors<Config>::ENDC
              << " wins from nodes:\n"
              << "  {}\n";
  }
  std::cout << "with strategy\n"
            << "  [" << graph.label(0)
            << "->" << graph.label(strategy[0]);
  for (Index ii = 1; ii < graph.n_min(); ++ii) {
    std::cout << "," << graph.label(ii)
      << "->" << graph.label(strategy[ii]);
  }
  std::cout << "]\n";
  if (v_positive.size() > 0) {
    std::cout << "Player " << utils::AnsiColors<Config>::YELLOW
              << "max" << utils::AnsiColors<Config>::ENDC
              << " wins from nodes:\n"
              << "  {" << graph.label(v_positive[0]);
    for (Index ii = 1; ii < static_cast<Index>(v_positive.size()); ++ii) {
      std::cout << ", " << graph.label(v_positive[ii]);
    }
    std::cout << "}\n";
  } else {
    std::cout << "Player " << utils::AnsiColors<Config>::YELLOW
              << "max" << utils::AnsiColors<Config>::ENDC
              << " wins from nodes:\n"
              << "  {}\n";

  }
  std::cout << "with strategy\n"
            << "  [" << graph.label(graph.n_min())
            << "->" << graph.label(strategy[graph.n_min()]);
  for (Index ii = graph.n_min() + 1; ii < graph.n(); ++ii) {
    std::cout << "," << graph.label(ii)
              << "->" << graph.label(strategy[ii]);
  }
  std::cout << "]\n";

  std::cout << '\n' << std::flush;
  // }}} End of printing results.
}

// Main with config. {{{
template<typename Config> inline int main_with_config(const Config& config) {
  using namespace mean_compass;
  using Real    = typename Config::Real;

  std::cout << utils::AnsiColors<Config>::GREEN
            << "Starting..."
            << utils::AnsiColors<Config>::ENDC << '\n'
            << std::flush;

  Real::default_precision(config.default_precision());

  if (config.input_files().size() == 0) {
      std::cout << utils::AnsiColors<Config>::GREEN
        << "Processing stdin"
        << utils::AnsiColors<Config>::ENDC << '\n' << std::flush;
    try {
      handle_graph(config, Graph<Config>(config, UTF8Input()));
    } catch (utils::SIGINTException& e) {
      std::cout << utils::AnsiColors<Config>::YELLOW
                << "\nSIGINT caught while processing stdin"
                << utils::AnsiColors<Config>::ENDC << '\n'
                << std::flush;
    }
  } else {
    for (const std::string& input_file : config.input_files()) {
      if (utils::sigint_caught) break;
      std::cout << utils::AnsiColors<Config>::GREEN
        << "Processing file: " << input_file
        << utils::AnsiColors<Config>::ENDC << '\n' << std::flush;
      try {
        handle_graph(config, Graph<Config>(config, UTF8Input(input_file)));
      } catch (utils::SIGINTException& e) {
        std::cout << utils::AnsiColors<Config>::YELLOW
                  << "\nSIGINT caught while processing file: " << input_file
                  << utils::AnsiColors<Config>::ENDC << '\n'
                  << std::flush;
        break;
      }
    }
  }

  std::cout << utils::AnsiColors<Config>::GREEN
            << "Exiting ;-)"
            << utils::AnsiColors<Config>::ENDC << '\n'
            << std::flush;

  return 0;
}  // }}} End of main with config.

}  // anonymous namespace

constexpr const char* version_string = "Mean Compass version 0.1.1";

int main(int argc, char** argv) {  // {{{
  // FIXME: Only the most basic things in main().
  using namespace mean_compass;

  utils::setup_sigint_handler();

  // Parse cmdline flags. {{{

  // Static options:
  bool option_verbose = false;
  bool option_use_colors = false;

  // Dynamic options:
  DynamicConfig config;
  config.parity(false)
        .default_precision(256)
        .display_precision(4);
  std::cout << std::setprecision(config.display_precision());

  try {
    namespace po = boost::program_options;
    po::options_description generic_options("Generic options");
    generic_options.add_options()
      ("help,h", "Print help and usage message.")
      ("verbose,v",
            po::bool_switch(&option_verbose),
            "Be more verbose.")
      ("color,c",
            po::bool_switch(&option_use_colors),
            "Use ANSI terminal colors.")
      ("version", "Print version of the program.")
      ("config-file",
           po::value<std::string>(&config.config_file_name()),
           "Parse options from configuration file given.")
      ("input-file,i",
           po::value<std::vector<std::string>>(&config.input_files())->composing(),
           "The input files to process, can be used multiple times.")
      ("dump-file,D",
           po::value<std::string>(&config.dump_file()),
           "The dump file, that will store the state of the last computation, "
           "if interrupted by SIGINT. By default, the dump is appended to the file.")
      ("overwrite-dump-file",
           po::bool_switch(&config.overwrite_dump_file()),
           "If this flag is set, the dump file will be overwritten. "
           "Use with caution, may cause loss of data.")
      ("restore-file,R",
           po::value<std::string>(&config.restore_file()),
           "Restore the computation state from the given dump file.");
    po::options_description config_options("Configuration");
    config_options.add_options()
      // We use int, because program_options parses "-5" with unsigned types.
      ("parity-game,g",
           po::bool_switch(&config.parity()),
           "The input graph specifies a parity game, even=max")
      ("default-precision,p",
           po::value<int>(&config.default_precision())->notifier(
             utils::check_range<int, 2, std::numeric_limits<int>::max()>),
           "Set the default precision of MPFR.")
      ("display-precision,d",
           po::value<int>(&config.display_precision())->notifier(
             utils::check_range<int, 2, std::numeric_limits<int>::max()>),
           "Set the precision of the output.")
      ("barrier-multiplier,b",
           po::value<std::string>(&config.barrier_multiplier_str()),
           "Set the value the barrier coefficient is multiplied by each turn.")
      ("barrier-adjustment,a",
           po::bool_switch(&config.barrier_adjustment()),
           "Make the barrier step variable.");
    po::options_description all_options;
    all_options.add(generic_options).add(config_options);
    po::positional_options_description positional_description;
    positional_description.add("input-file", -1);

    po::variables_map variables_map;
    po::store(po::command_line_parser(argc, argv)
            .options(all_options)
            .positional(positional_description).run(),
        variables_map);
    po::notify(variables_map);
    if (variables_map.count("config-file")) {
      std::ifstream config_file(config.config_file_name());
      po::store(po::parse_config_file(config_file, all_options),
          variables_map);
      variables_map.notify();
    }
    if (variables_map.count("version")) {
      std::cout << version_string << '\n';
      return 0;
    }
    if (variables_map.count("help")) {
      std::cout << version_string << '\n';
      std::cout << all_options << '\n';
      return 0;
    }
    if (variables_map.count("default-precision")) {
      std::cout << std::setprecision(config.default_precision());
    }
    if (variables_map.count("display-precision")) {
      std::cout << std::setprecision(config.display_precision());
    }
  } catch (std::exception& e) {
    std::cout << e.what() << '\n' << std::flush;
    return 1;
  }
  // End of parsing cmdline flags }}}

  // Call main_with_config. {{{
  if (!option_verbose && !option_use_colors) {
    return main_with_config(static_cast<Config<false, false>>(config));
  } else if (!option_verbose && option_use_colors) {
    return main_with_config(static_cast<Config<false, true>>(config));
  } else if (option_verbose && !option_use_colors) {
    return main_with_config(static_cast<Config<true, false>>(config));
  } else if (option_verbose && option_use_colors) {
    return main_with_config(static_cast<Config<true, true>>(config));
  } else {
    assert(false);
  }
  // }}} End of main_with_config calls.

}  // }}} End of main.

// vim: et sw=2 ts=2 foldmethod=marker
