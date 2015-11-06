#ifndef _PFASST__JSON_WRITER__HPP_
#define _PFASST__JSON_WRITER__HPP_

#include <string>
using std::string;

#include <leathers/push>
#include <leathers/all>
#include <json/json.h>
#include <leathers/pop>


namespace pfasst
{
  class JsonWriter
  {
    private:
      Json::Value _json;

      JsonWriter() = default;
      virtual ~JsonWriter() = default;

    public:
      JsonWriter(const JsonWriter& other) = delete;
      JsonWriter(JsonWriter&& other) = delete;
      JsonWriter& operator=(const JsonWriter& other) = delete;
      JsonWriter& operator=(JsonWriter&& other) = delete;

      static JsonWriter& instance();
      Json::Value& json();
      static void to_file(const string& filename);
  };
}  // ::pfasst

#include "pfasst/json_writer_impl.hpp"

#endif  // _PFASST__JSON_WRITER__HPP_
