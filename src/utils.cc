/* Copyright (c) 2016 dtldarek
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "utils.h"

namespace mean_compass {

namespace utils {

// SIGINT handler {{{
sig_atomic_t sigint_caught = 0;

namespace {

void (*prev_sigint_handler)(int) = SIG_IGN;
void sigint_handler(int signum) {
  assert(signum == SIGINT);
  sigint_caught = 1;
  if (prev_sigint_handler != SIG_IGN && prev_sigint_handler != SIG_DFL) {
    (*prev_sigint_handler)(signum);
  }
}

}  // anonymous namespace

void setup_sigint_handler() {
  prev_sigint_handler = signal(SIGINT, sigint_handler);
}
// }}} End of SIGINT handler.


}  // namespace utils
}  // namespace mean_compass

// vim: et sw=2 ts=2 foldmethod=marker
