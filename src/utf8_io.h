/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef __MEAN_COMPASS_UTF8_IO_H__
#define __MEAN_COMPASS_UTF8_IO_H__

#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>
// By intention we do not use "types.h" and class Config.

namespace mean_compass {

class UTF8Input {
 public:
  // By intention we do not use the Config class.
  // We redefine these two types, because (unfortunatelly)
  // our implementation depends on the actual type.
  using Integer  = boost::multiprecision::mpz_int;
  using Real     = boost::multiprecision::mpfr_float;

  static constexpr const char COMMENT = 0x23;  // ASCII code for '#';
  static constexpr const char PLUS    = 0x2B;  // ASCII code for '+';
  static constexpr const char MINUS   = 0x2D;  // ASCII code for '-';
  static constexpr const char SLASH   = 0x2F;  // ASCII code for '/';
  static constexpr const char ZERO    = 0x30;  // ASCII code for '0';
  static constexpr const char EQUAL   = 0x3D;  // ASCII code for '=';
  static constexpr const char UPPERCASE_A = 0x41;  // ASCII code for 'A';
  static constexpr const char UPPERCASE_F = 0x46;  // ASCII code for 'F';
  static constexpr const char UPPERCASE_Z = 0x5A;  // ASCII code for 'Z';
  static constexpr const char LOWERCASE_A = 0x61;  // ASCII code for 'a';
  static constexpr const char LOWERCASE_F = 0x66;  // ASCII code for 'f';
  static constexpr const char LOWERCASE_Z = 0x7A;  // ASCII code for 'z';

  class EndOfFileE : public std::runtime_error {
   public:
    EndOfFileE(const char* msg) : std::runtime_error(msg) { }
    EndOfFileE(const std::string& msg) : std::runtime_error(msg) { }
  };

  UTF8Input();  // Use stdin.
  UTF8Input(const std::string& filename);

  UTF8Input(const UTF8Input&) = delete;
  UTF8Input(UTF8Input&&) = delete;
  UTF8Input& operator=(const UTF8Input&) = delete;
  UTF8Input& operator=(UTF8Input&&) = delete;

  ~UTF8Input() { finalize(); }

  bool eof() const { return feof(file); }
  bool ok() const { return !ferror(file); }
  explicit operator bool() const { return ok() && !eof(); }

  // Read a single UTF-8 symbol from FILE* using getc_unlocked.
  // Throws std::runtime_exception if an error is encountered.
  char32_t get_char();

  inline unsigned int get_uint() { return get_simple<unsigned int>(); }
  inline int get_int() { return get_simple<int>(); }
  inline unsigned long get_ulong() { return get_simple<unsigned long>(); }
  inline long get_long() { return get_simple<long>(); }
  inline unsigned long long get_ulonglong() { return get_simple<unsigned long long>(); }
  inline long long get_longlong() { return get_simple<long long>(); }

  inline Integer get_integer() { return get_simple<Integer>(); }
  inline Integer get_integer_base64() { Integer i; get_base64(&i); return i; }
  inline Real get_real() { return get_simple<Real>(); }

  // Get first anwp chunk of letters.
  std::string get_string() { return get_simple<std::string>(); }

  // We try to avoid mixing templates and ordinary overloading.
  inline UTF8Input& get(unsigned int* result) { return get_nostrip(result, strip_nonanwp()); }
  inline UTF8Input& get(int* result) { return get_nostrip_signed(result, strip_nonanwp()); }
  inline UTF8Input& get(unsigned long* result) { return get_nostrip(result, strip_nonanwp()); }
  inline UTF8Input& get(long* result) { return get_nostrip_signed(result, strip_nonanwp()); }
  inline UTF8Input& get(unsigned long long* result) { return get_nostrip(result, strip_nonanwp()); }
  inline UTF8Input& get(long long* result) { return get_nostrip_signed(result, strip_nonanwp()); }
  UTF8Input& get(Integer* result);
  UTF8Input& get(Real* result);
  UTF8Input& get(std::string* result);
  UTF8Input& get_base64(Integer* result);

  template<typename T>
    UTF8Input& operator>>(T& result) { return get(&result); }

 protected:
  // NWP stands for ASCII non-white printable.
  static inline bool is_anwp(char32_t c) { return 0x21 <= c && c <= 0x7E; }
  static inline bool is_digit(char32_t c) { return ZERO <= c && c <= (ZERO + 9); }
  static inline bool is_hex(char32_t c) { return is_digit(c) || is_hexletter(c); }
  static inline bool is_hexletter(char32_t c) {
    uint_least32_t lowercase = static_cast<uint_least32_t>(c) | 0x20;
    return LOWERCASE_A <= lowercase && lowercase <= LOWERCASE_F;
  }
  static constexpr const char base64tab[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 62, -1, -2, -1, 63,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, -1, -1, -1, -3, -1, -1,
    -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, -1, -1, -1, -1, -1,
    -1, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

  static inline bool is_base64(char32_t c) {
    uint_least32_t u = static_cast<uint_least32_t>(c);
    return (u >> 8) == 0 && base64tab[u & 0xFF] >= 0;
  }
  // We assume that the input is valid, i.e. it is a digit or hex.
  // We use int, because for types of the same rank we have
  // (signed) + (unsigned) = (unsigned).
  static inline int d2i(char32_t c) { return c - ZERO; }
  static inline int h2i(char32_t c) {
    return c <= (ZERO + 9) ? c - ZERO : (c <= UPPERCASE_F ? c - UPPERCASE_A + 10 : c - LOWERCASE_A + 10);
  }

  // Take ownership of file_.
  // Assumes file_ is not locked.
  UTF8Input(FILE* file_);

  FILE* file;
  bool finalized;
  void finalize() noexcept;

  inline bool strip_comments(char32_t c) {
    if (c == COMMENT) {  // COMMENT == 0x23 is the ASCII code for '#'
      while (ok() && !eof() && get_char() != '\n') { }
      return true;
    }
    return false;
  }
  inline char32_t strip_nonanwp() {
    char32_t c;
    while (!is_anwp(c = get_char()) || strip_comments(c)) { }
    return c;
  }

  template<typename F, typename C> UTF8Input& get_cps(F validate, C continuation);
  template<typename T> UTF8Input& get_nostrip(T* result, char32_t c);
  template<typename T> UTF8Input& get_nostrip_signed(T* result, char32_t c);
  template<typename T> inline T get_simple() { T r; get(&r); return r; }
};

}  // namespace mean_compass

#include "utf8_io.hpp"

#endif  // __MEAN_COMPASS_UTF8_IO_H__
