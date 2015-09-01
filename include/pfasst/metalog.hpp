#ifndef METALOG_HPP_
#define METALOG_HPP_

#include <utility>

/*
 * The idea of making logging calls types was borrowed from Marc Eaddy's talk
 * "Pimp my Log()" at CppCon 2014 where he used this technique to optimize
 * logging performance.
 *
 * It can be found on YouTube:
 * https://www.youtube.com/watch?v=TS_waQZcZVc
 */

namespace metalog {

  struct None {};

  template<typename First, typename Second>
  struct Pair {
      First  first;
      Second second;
  };

  template<typename List>
  struct Data {
      List list;
  };

  template<typename Head, typename Element>
  constexpr Data<Pair<Head, const Element&>> operator<<(Data<Head>&& head, Element&& element) {
      return {{ std::forward<Head>(head.list), std::forward<Element>(element) }};
  }

  template<typename Head, std::size_t N>
  constexpr Data<Pair<Head, const char*>> operator<<(Data<Head>&& head, const char (&element)[N]) {
      return {{ std::forward<Head>(head.list), element }};
  }

  template<typename Stream>
  inline void print(Stream&, None) {}

#ifndef METALOG_NOLOG
  template<typename Stream, typename Head, typename Tail>
  void print(Stream &strm, const Pair<Head, Tail>&& data) {
      print(strm, std::move(data.first));
      strm << data.second;
  }
#else
  template<typename Stream, typename... Args>
  void print(Stream& strm, Args... args) {
  }
#endif

  template<typename Stream, typename List>
  void log(Stream& strm, Data<List>&& data) {
      print(strm, std::move(data.list));
  }

}

#ifndef ML_NOLOG
  #define ML_LOG(level, x) \
      (metalog::log((LOG(level)), metalog::Data<metalog::None>() << x))
  #define ML_CLOG(level, logger_id, x) \
      (metalog::log((CLOG(level, logger_id)), metalog::Data<metalog::None>() << x))
  #define ML_CLOG_IF(condition, level, logger_id, x) \
      CLOG_IF(condition, level, logger_id) << x
  #define ML_CVLOG(verbose_level, logger_id, x) \
      CVLOG(verbose_level, logger_id) << x
  #define ML_CVLOG_IF(condition, verbose_level, logger_id, x) \
      CVLOG_IF(condition, verbose_level, logger_id) << x
#else
  #define ML_LOG(level, x)
  #define ML_CLOG(level, logger_id, x)
  #define ML_CLOG_IF(condition, level, logger_id, x)
  #define ML_CVLOG(verbose_level, logger_id, x)
  #define ML_CVLOG_IF(condition, verbose_level, logger_id, x)
#endif

#endif
