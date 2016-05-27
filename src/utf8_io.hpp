/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "utf8_io.h"
#include <vector>

namespace mean_compass {

template<typename F, typename C>
UTF8Input& UTF8Input::get_cps(F validate, C continuation) {
  char32_t c = strip_nonanwp();
  if (c == std::char_traits<char32_t>::eof()) {
    throw EndOfFileE("Encountered an unexpected end-of-file marker");
  }
  size_t size = 0;
  std::vector<std::string> buffer = { std::string(1, c) };
  buffer[0].reserve(24);  // 2**64 has 20 digits.
  while (true) {
    const size_t capacity = buffer.back().capacity();
    for (size_t ii = buffer.back().size(); ii < capacity; ++ii) {
      if (validate(c = get_char())) {
        buffer.back() += c;
      } else {
        size += buffer.back().size();
        return continuation(std::move(buffer), size);
      }
    }
    size += buffer.back().size();
    buffer.push_back(std::string());
    buffer.back().reserve(size);  // Double the current size.
  }
}

// UTF8Input::get_nostrip (i.e. the unsigned version) {{{
template<typename T>
UTF8Input& UTF8Input::get_nostrip(T* result, char32_t c) {
  if (is_digit(c)) {
    *result = d2i(c);
    while (is_digit(c = get_char()))
      *result = *result * 10 + d2i(c);
    return *this;
  } else if (c == std::char_traits<char32_t>::eof()) {
    throw EndOfFileE("Encountered an unexpected end-of-file marker");
  } else {
    std::stringstream description;
    description << "Encountered an unexpected non-digit character: 0x"
      << std::uppercase << std::setfill(ZERO) << std::setw(2)
      << std::hex << static_cast<uint_least32_t>(c);
    throw std::runtime_error(description.str());
  }
  return *this;
}  // }}} End get_nostrip_unsigned

// UTF8Input::get_nostrip_signed {{{
template<typename T>
UTF8Input& UTF8Input::get_nostrip_signed(T* result, char32_t c) {
  if (c != '-') {
    // The is_digit check is done in the %_unsigned call.
    return get_nostrip(result, c);
  } else {
    get_nostrip(result, get_char());
    *result = -*result;
    return *this;
  }
}  // }}} End get_nostrip_signed

}  // namespace mean_compass

// vim: foldmethod=marker
