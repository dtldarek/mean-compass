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

// We pass all the parameters and options directly.
// Should there be too many of them, we can put them into Config class.
template<typename Config> inline int main_with_config(
    const std::vector<std::string>& input_files,
    const int option_default_precision) {
  using namespace mean_compass;
  using Index   = typename Config::Index;
  using Real    = typename Config::Real;
  using Matrix  = typename Config::Matrix;
  using Vector  = typename Config::Vector;
  using Triplet = typename Config::Triplet;

  std::cout << utils::AnsiColors<Config>::GREEN
            << "Starting..."
            << utils::AnsiColors<Config>::ENDC << '\n'
            << std::flush;

  // Do some tests. {{{
  const size_t size = 100;

  Real::default_precision(option_default_precision);
  std::cout << std::setprecision(Real::default_precision());

  std::vector<Triplet> triplets;
  triplets.reserve(10*size);
  for (size_t ii = 0; ii < 10*size; ++ii) {
    triplets.push_back(Triplet(
          utils::unirand<Index>(size),
          utils::unirand<Index>(size),
          utils::unirand()));
  }
  Matrix A(size, size);
  A.setFromTriplets(triplets.begin(), triplets.end());
  A.makeCompressed();
  Vector b = Vector::Random(size);

  typename Config::LU solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  Vector x = solver.solve(b);
  std::cout << "relative error: " << (A*x - b).norm() / b.norm() << std::endl;
  // End doing some tests. }}}

  for (const std::string& input_file : input_files) {
    std::cout << utils::AnsiColors<Config>::GREEN
              << "Processing " << input_file
              << utils::AnsiColors<Config>::ENDC << '\n' << std::flush;
    Graph<Config> graph((UTF8Input(input_file)));
    std::cout << graph.n() << '\n';
    graph.init_state(0.001, 0.001);
    for (Index ii = 0; ii < graph.n(); ++ii) {
      std::cout << ' ' << graph.position()(ii);
    }
    std::cout << '\n';

    typename Graph<Config>::MinProblem min_problem = graph.get_min_problem(0.5, 0.5);
    SimpleNewton<Config, typename Graph<Config>::MinProblem> min_newton(&min_problem);
    min_newton.loop();
    for (Index ii = 0; ii < min_newton.position().size(); ++ii) {
      std::cout << ' ' << min_newton.position()(ii);
    }
    std::cout << '\n';
  }

  std::cout << utils::AnsiColors<Config>::GREEN
            << "Exiting ;-)"
            << utils::AnsiColors<Config>::ENDC << '\n'
            << std::flush;

  return 0;
}

}  // anonymous namespace

constexpr const char* version_string = "Mean Compass version 0.1.0";

int main(int argc, char** argv) {
  // FIXME: Only the most basic things in main().
  using namespace mean_compass;

  // Parse cmdline flags. {{{
  bool option_verbose = false;
  bool option_use_colors = false;
  std::string option_config_file_name;
  std::vector<std::string> input_files;
  int option_default_precision = 256;

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
           po::value<std::string>(&option_config_file_name),
           "Parse options from configuration file given.");
    po::options_description config_options("Configuration");
    config_options.add_options()
      // We use int, because program_options parses "-5" with unsigned types.
      ("default-precision,p",
           po::value<int>(&option_default_precision)->notifier(
             utils::check_range<int, 2, std::numeric_limits<int>::max()>),
           "Set the default precision of MPFR.")
      ("input-file",
           po::value<std::vector<std::string>>(&input_files)->composing(),
           "The input files to process, can be used multiple times.");
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
      std::ifstream config_file(option_config_file_name);
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
  } catch (std::exception& e) {
    std::cout << e.what() << '\n' << std::flush;
    return 1;
  }
  // End of parsing cmdline flags }}}

  if (option_use_colors) {
    return main_with_config<Config<true>>(input_files, option_default_precision);
  } else {
    return main_with_config<Config<false>>(input_files, option_default_precision);
  }

}

// vim: foldmethod=marker
