#include "pfasst/json_writer.hpp"

#include <fstream>
#include <memory>

#include "pfasst/logging.hpp"


namespace pfasst
{
  JsonWriter&
  JsonWriter::instance()
  {
    static JsonWriter instance;
    return instance;
  }

  Json::Value&
  JsonWriter::json()
  {
    return this->_json;
  }

  void
  JsonWriter::to_file(const string& filename)
  {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "  ";
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());

    std::ofstream fh(filename, std::ios_base::out);
    if (fh.is_open()) {
      writer->write(JsonWriter::instance().json(), &fh);
    } else {
      ML_CLOG(ERROR, "OUTPUT", "Cannot open file: " << filename);
    }
    // file handle is closed via dtor as it goes out of scope
  }
}  // ::pfasst
