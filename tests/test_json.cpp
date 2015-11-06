#include "fixtures/test_helpers.hpp"
using ::testing::Eq;
using ::testing::StrEq;

#include <fstream>
#include <string>

#include <leathers/push>
#include <leathers/all>
#include <json/json.h>
#include <leathers/pop>

#include <pfasst/json_writer.hpp>


TEST(JsonWriter, is_singleton)
{
  auto& writer = pfasst::JsonWriter::instance();
  auto& other = pfasst::JsonWriter::instance();
  EXPECT_EQ(&writer, &other);
}


class ToFileTest
  : public ::testing::Test
{
  public:
    Json::Value example;
    std::string expected_string;
    
  virtual void SetUp()
  {
    this->example["string value"] = "a string";
    this->example["nested map"]["integer"] = 42;
    this->example["array"] = Json::Value(Json::arrayValue);
    this->example["array"][0] = Json::Value(1);
    this->example["array"][1] = Json::Value(2);
    this->example["array"][2] = Json::Value(3);
    
    this->expected_string = "{\n  \"array\" : \n  [\n    1,\n    2,\n    3\n  ],\n  \"nested map\" : \n  {\n    \"integer\" : 42\n  },\n  \"string value\" : \"a string\"\n}\n";
  }
};

TEST_F(ToFileTest, write)
{
  pfasst::JsonWriter::instance().json() = this->example;
  pfasst::JsonWriter::to_file("testfile.json");

  std::ifstream fh("testfile.json", std::ios_base::in);
  ASSERT_TRUE(fh.is_open());
  std::string content = "", line = "";
  while (getline(fh, line)) {
    content += line + "\n";
  };

  EXPECT_THAT(content, StrEq(this->expected_string));

  // fh is closed via dtor as it goes out of scope
}

TEST_F(ToFileTest, fails_silently_on_failed_file_creation) {
  pfasst::JsonWriter::to_file("/should/not/be/able/to/write.here");
}


TEST_MAIN()
