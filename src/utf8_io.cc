/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "utf8_io.h"
#include <memory>
#include <vector>

namespace mean_compass {

namespace {

class CouldNotLockFileE : public std::runtime_error {
 public:
  CouldNotLockFileE(const char* msg) : std::runtime_error(msg) { }
  CouldNotLockFileE(const std::string& msg) : std::runtime_error(msg) { }
};

FILE* open_or_throw(const std::string& filename) {
  FILE* file = fopen(filename.c_str(), "rb");
  if (!file) {
    std::stringstream description;
    description << "Could not open file: " << filename;
    throw std::runtime_error(description.str());
  }
  return file;
}

}  // anonymous namespace

// UTF8Input::UTF8Input and UTF8Input::finalize {{{
UTF8Input::UTF8Input()
try : UTF8Input(stdin) {
  // Nothing to do here.
} catch (CouldNotLockFileE& e) {
  throw CouldNotLockFileE("Could not lock stdin");
}

UTF8Input::UTF8Input(const std::string& filename)
try : UTF8Input(open_or_throw(filename)) {
} catch (CouldNotLockFileE& e) {
  std::stringstream description;
  description << "Could not lock file: " << filename;
  throw CouldNotLockFileE(description.str());
}

// Take ownership of file_, assumes the file is not locked.
UTF8Input::UTF8Input(FILE* file_) : file(nullptr), finalized(false) {
  assert(file_ != nullptr);
  if (ftrylockfile(file_) != 0) {
    fclose(file_);
    throw CouldNotLockFileE("Could not lock the given file");
  }
  file = file_;
}

void UTF8Input::finalize() noexcept {
  if (!finalized && file) {
    finalized = true;
    funlockfile(file);
    fclose(file);
    file = nullptr;
  }
}
// }}} End constructors

// UTF8Input::get_char {{{
char32_t UTF8Input::get_char() {
  // We are not using the standard library, because it's quite slow.
  // Although fread could be actually even faster, I would like to avoid
  // writtng my own buffering here.
  char c = getc_unlocked(file);
  if (c >= 0) {
    return c;
  } else {
    unsigned char u = static_cast<unsigned char>(c);
    uint_least32_t u_result;
    int size;
    if (c == std::char_traits<char>::eof()) {
      return std::char_traits<char32_t>::eof();
    } else if ((u & 0xE0) == 0xC0) {
      u_result = u & 0x1F;
      size = 1;
    } else if ((u & 0xF0) == 0xE0) {
      u_result = u & 0x0F;
      size = 2;
    } else if ((u & 0xF8) == 0xF0) {
      u_result = u & 0x07;
      size = 3;
    } else if ((u & 0xFC) == 0xF8) {
      u_result = u & 0x03;
      size = 4;
    } else if ((u & 0xFE) == 0xFC) {
      u_result = u & 0x01;
      size = 5;
    } else {
      std::stringstream description;
      description << "Encountered invalid initial UTF-8 character: 0x"
        << ((u & 0xF0) >> 4) + 'A'
        << ((u & 0x0F) >> 0) + 'A';
      throw std::runtime_error(description.str());
    }
    for (int ii = 0; ii < size; ++ii) {
      c = getc_unlocked(file);
      u = static_cast<unsigned char>(c);
      if (c == std::char_traits<char>::eof()) {
        throw std::runtime_error("Encountered EOF while reading an UTF-8 symbol");
      } else if ((u & 0xC0) == 0x80) {
        u_result = (u_result << 6) + (u & 0x3f);
      } else {
        std::stringstream description;
        description << "Encountered invalid UTF-8 character: 0x"
          << ((u & 0xF0) >> 4) + 'A'
          << ((u & 0x0F) >> 0) + 'A';
        throw std::runtime_error(description.str());
      }
    }
    return static_cast<char32_t>(u_result);
  }
}  // }}} End get_char

UTF8Input& UTF8Input::get(Integer* result) {
  std::string s;
  get(&s);
  result->assign(s);
  return *this;
}
UTF8Input& UTF8Input::get(Real* result) {
  std::string s;
  get(&s);
  if (result->precision() < s.size() + 1) {
    result->precision(s.size()+1);
  }
  result->assign(s);
  return *this;
}
UTF8Input& UTF8Input::get(std::string* result) {
  auto validate = [this](char32_t c) { return this->is_anwp(c); };
  auto continuation = [result, this](std::vector<std::string>&& buffer_,
                                     size_t size) -> UTF8Input& {
    std::vector<std::string> buffer(buffer_);
    result->clear();
    result->reserve(size);
    for (const std::string& s : buffer) {
      result->append(s);
    }
    return *this;
  };
  return get_cps(validate, continuation);
}
UTF8Input& UTF8Input::get_base64(Integer* result) {
  auto validate = [this](char32_t c) { return this->is_anwp(c); };
  auto continuation = [this, result](std::vector<std::string>&& buffer_,
                                     size_t size) -> UTF8Input& {
    std::vector<std::string> buffer(buffer_);
    uint_least32_t mid_buffer = 0;
    uint_least32_t bytes_in_mid_buffer = 0;
    std::vector<unsigned char> bytes;
    size_t ii = 0;
    size_t negated = 0;
    size_t padding = 0;
    if (buffer.front().front() == MINUS) {
      negated = 1;  // Size of the '-' sign.
      ii++;
    }
    if ((size - ii) % 4 != 0) {
      std::string description;
      description.reserve(32 + size);
      description += "Invalid base64 encoding: ";
      for (const std::string& s : buffer)
        description += s;
      throw std::runtime_error(description);
    }
    if (buffer.back().back() == EQUAL) {
      buffer.back().pop_back();
      padding += 1;
      // There are at most two '=' padding symbols.
      if (buffer.back().back() == EQUAL) {
        buffer.back().pop_back();
        padding += 1;
      }
    }
    bytes.reserve((size-ii) / 4 * 3 - padding);
    for (const std::string& s : buffer) {
      for (const char& c : s) {
        char data = base64tab[static_cast<unsigned char>(c)];
        if (data >= 0) {
          mid_buffer <<= 6;
          mid_buffer += base64tab[static_cast<unsigned char>(c)];
          bytes_in_mid_buffer += 6;
        } else {
          std::stringstream description;
          description << "Invalid base64 symbol: 0x"
                      << std::hex << std::setfill(ZERO)
                      << std::setw(2) << data;
          throw std::runtime_error(description.str());
        }
        if (bytes_in_mid_buffer >= 8) {
          bytes_in_mid_buffer -= 8;
          bytes.push_back((mid_buffer >> bytes_in_mid_buffer) & 0xFF);
        }
      }
    }
    if (bytes.size() != (size-ii) / 4 * 3 - padding) {
      std::stringstream description;
      description << "Read " << bytes.size()
                  << " of expected " << bytes.capacity() << " bytes";
      throw std::runtime_error(description.str());
    }
    *result = 0;
    mpz_import (result->backend().data(), bytes.size(), 1, sizeof(unsigned char), 0, 0, bytes.data());
    if (negated != 0) result->backend().negate();
    return *this;
  };
  return get_cps(validate, continuation);
}

constexpr const char UTF8Input::base64tab[256];

}  // namespace mean_compass

// vim: foldmethod=marker
