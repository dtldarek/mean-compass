/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef __MEAN_COMPASS_UTILS_H__
#define __MEAN_COMPASS_UTILS_H__

#include <sstream>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include "types.h"


namespace mean_compass {
namespace utils {

struct AnsiColors {
  static constexpr const char* const RED = "\033[91m";
  static constexpr const char* const GREEN = "\033[92m";
  static constexpr const char* const YELLOW = "\033[93m";
  static constexpr const char* const BLUE = "\033[94m";
  static constexpr const char* const MAGENTA = "\033[95m";
  static constexpr const char* const CYAN = "\033[96m";
  static constexpr const char* const WHITE = "\033[97m";
  static constexpr const char* const GRAY = "\033[90m";
  static constexpr const char* const ENDC = "\033[0m";
  static constexpr const char* const BOLD = "\033[1m";
  static constexpr const char* const UNDERLINE = "\033[4m";
};

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


}  // namespace utils
}  // namespace mean_compass

#endif  // __MEAN_COMPASS_UTILS_H__
