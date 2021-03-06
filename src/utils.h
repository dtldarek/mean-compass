/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef __MEAN_COMPASS_UTILS_H__
#define __MEAN_COMPASS_UTILS_H__

#include <csignal>
#include <sstream>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include "types.h"


namespace mean_compass {
namespace utils {

extern sig_atomic_t sigint_caught;
void setup_sigint_handler();  // Should be called only once.
class SIGINTException : public std::exception {
 public:
  SIGINTException() : std::exception() { }
};


template<typename Config> struct AnsiColors {
  static constexpr bool use_colors = Config::use_colors;
  static constexpr const char* const RED       = use_colors ? "\033[91m" : "";
  static constexpr const char* const GREEN     = use_colors ? "\033[92m" : "";
  static constexpr const char* const YELLOW    = use_colors ? "\033[93m" : "";
  static constexpr const char* const BLUE      = use_colors ? "\033[94m" : "";
  static constexpr const char* const MAGENTA   = use_colors ? "\033[95m" : "";
  static constexpr const char* const CYAN      = use_colors ? "\033[96m" : "";
  static constexpr const char* const WHITE     = use_colors ? "\033[97m" : "";
  static constexpr const char* const GRAY      = use_colors ? "\033[90m" : "";
  static constexpr const char* const ENDC      = use_colors ? "\033[0m" : "";
  static constexpr const char* const BOLD      = use_colors ? "\033[1m" : "";
  static constexpr const char* const UNDERLINE = use_colors ? "\033[4m" : "";
};

template<typename Config> constexpr const char* const AnsiColors<Config>::RED;
template<typename Config> constexpr const char* const AnsiColors<Config>::GREEN;
template<typename Config> constexpr const char* const AnsiColors<Config>::YELLOW;
template<typename Config> constexpr const char* const AnsiColors<Config>::BLUE;
template<typename Config> constexpr const char* const AnsiColors<Config>::MAGENTA;
template<typename Config> constexpr const char* const AnsiColors<Config>::CYAN;
template<typename Config> constexpr const char* const AnsiColors<Config>::WHITE;
template<typename Config> constexpr const char* const AnsiColors<Config>::GRAY;
template<typename Config> constexpr const char* const AnsiColors<Config>::ENDC;
template<typename Config> constexpr const char* const AnsiColors<Config>::BOLD;
template<typename Config> constexpr const char* const AnsiColors<Config>::UNDERLINE;

// For any integer-like type T returns a random number of that type
// that is uniformly distributed in [0, size).
template<typename I = default_types::Integer> inline I unirand(I size) {
  static boost::random_device random_device;
  static boost::random::mt19937 generator(random_device);
  assert(size >= 1);
  boost::random::uniform_int_distribution<I> uniform_integer(0, size - 1);
  return uniform_integer(generator);
}

// Returns a random Real number uniformly distributed in [0,1).
inline default_types::Real unirand() {
  default_types::Integer denominator = default_types::Integer(1) << default_types::Real::default_precision();
  default_types::Integer nominator = unirand(denominator);
  default_types::Rational result(nominator, denominator);
  return static_cast<default_types::Real>(result);
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

template<typename Real> inline auto asqrt(Real value)
    -> decltype(boost::multiprecision::scalbn(Real(0.0), int(0))) {
  int exp = -boost::multiprecision::ilogb(value) / 2;
  return boost::multiprecision::scalbn(value, exp);
}

}  // namespace utils
}  // namespace mean_compass

#endif  // __MEAN_COMPASS_UTILS_H__

// vim: et sw=2 ts=2
