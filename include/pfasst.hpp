#ifndef _PFASST_HPP_
#define _PFASST_HPP_

#include <functional>

#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  /**
   * Global entry point for _PFASST_.
   *
   * This function must be called at runtime before any other methods of the _PFASST_ library.
   * In case of a MPI build, `MPI_Init` must be called prior to this function`.
   *
   * @param[in] argc  number of command line arguments as given by `main()`
   * @param[in] argv  list of command line arguments as given by `main()`
   * @param[in] opts  optional std::function initializing further additional parameters via
   *                  config::options::add_option()
   * @param[in] logs  optional std::function initializing further additional loggers via
   *                  log::add_custom_logger()
   *
   * @ingroup Utilities
   */
  inline static void init(int argc, char** argv,
                          std::function<void()> opts = nullptr,
                          std::function<void()> logs = nullptr)
  {
    if (opts) {
      opts();
    }
    config::init();
    config::read_commandline(argc, argv);
    log::start_log(argc, argv);
    auto to_log = config::check_unrecognized_args();
    if (!to_log.empty()) {
      for (const auto& l : to_log) {
        LOG(WARNING) << l;
      }
    }
    if (logs) {
      logs();
    }
  }
} // ::pfasst


/**
 * @defgroup Controllers Controllers
 *   Controllers represent the main entry point of PFASST++ as they ensemble the central algorithmic
 *   logic of PFASST and related algorithms.
 *
 * @defgroup Assistance Assistance
 *   The assistance module groups the user data wrappers, transfer operators, quadrature and
 *   communicators.
 *
 * @defgroup Internals Internals
 *   Entities listed in this module are meant to be for internal use only and should not be used
 *   outside of PFASST++ itself.
 *
 * @defgroup Traits Traits
 *   Traits offer compile-time introspection of template types.
 * @ingroup Internals
 *
 * @defgroup Tags Tags
 *   Tags are utility types to shortcut and generalize certain template expressions
 *   @ingroup Traits
 *
 * @defgroup Utilities Utilities
 *   General utility functions not directly related to PFASST++ but also handy in user code.
 *
 * @defgroup Contributed Contributed
 *   Miscellaneous utilities not necessarily required by _PFASST_, but a few examples or extensions.
 *
 * @defgroup Examples Examples
 *   A few different examples demonstrating the use and functionality of PFASST++.
 */

/**
 * PFASST++ version: closest release tag and a brief git hash.
 *
 * @var pfasst::VERSION
 * @ingroup Utilities
 */

#endif
