#ifndef _PFASST__TESTS__TEST_HELPERS_HPP_
#define _PFASST__TESTS__TEST_HELPERS_HPP_

#include <tuple>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

MATCHER(DoubleNear, "")
{
  return abs(get<0>(arg) - get<1>(arg)) < 1e-15;
}

MATCHER(MutuallyEqual, "")
{
  auto size = arg.size();
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = i + 1; j < size; ++j) {
      if (arg[i] != arg[j]) {
        return false;
      }
    }
  }
  return true;
}
#include <leathers/pop>


#include <pfasst/config.hpp>

class ConfigModder
{
  private:
    boost::program_options::variables_map backup;

  public:
    explicit ConfigModder()
      : backup(pfasst::config::options::get_instance().get_variables_map())
    {}
    ConfigModder(const ConfigModder& o) = delete;
    ConfigModder(ConfigModder&& o) = delete;
    ~ConfigModder()
    {
      pfasst::config::options::get_instance().get_variables_map() = move(this->backup);
    }

    ConfigModder& operator=(const ConfigModder& o) = delete;
    ConfigModder&& operator=(ConfigModder&& o) = delete;

    template<typename T>
    void update(const string& name, const T& value)
    {
      pfasst::config::options::update_value<T>(name, value);
    }
};


#include <pfasst/logging.hpp>

#ifdef WITH_MPI
#define TEST_MAIN() \
  int main(int argc, char** argv) { \
    MPI_Init(&argc, &argv); \
    pfasst::log::start_log(argc, argv); \
    InitGoogleTest(&argc, argv); \
    const int out = RUN_ALL_TESTS(); \
    MPI_Finalize(); \
    return out; \
  }
#else
#define TEST_MAIN() \
  int main(int argc, char** argv) { \
    pfasst::log::start_log(argc, argv); \
    InitGoogleTest(&argc, argv); \
    return RUN_ALL_TESTS(); \
  }
#endif


#include "fixtures/concepts.hpp"

#endif  // _PFASST__TESTS__TEST_HELPERS_HPP_
