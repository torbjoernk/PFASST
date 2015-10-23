#include "fixtures/test_helpers.hpp"
using ::testing::Eq;

using std::make_tuple;
using std::get;
using std::tuple;

#include <pfasst/util.hpp>


class IndexManipulation2D
  : public ::testing::TestWithParam<tuple<size_t, size_t>>
{
  public:
    size_t index;
    size_t ndofs;

    void SetUp()
    {
      this->index = get<0>(GetParam());
      this->ndofs = get<1>(GetParam());
    }
};

TEST_P(IndexManipulation2D, split_and_linearize)
{
  auto tuple = pfasst::split_index<2>(index, ndofs);
  size_t i = pfasst::linearized_index(tuple, ndofs);

  EXPECT_THAT(i, Eq(index));
}

INSTANTIATE_TEST_CASE_P(IndexManipulation2DTests,
                        IndexManipulation2D,
                        ::testing::Values(
                          make_tuple<size_t, size_t>(4, 3),
                          make_tuple<size_t, size_t>(3, 3),
                          make_tuple<size_t, size_t>(1, 3),
                          make_tuple<size_t, size_t>(5, 3)
                        )
                       );


class IndexManipulation3D
  : public ::testing::TestWithParam<tuple<size_t, size_t>>
{
  public:
    size_t index;
    size_t ndofs;

    void SetUp()
    {
      this->index = get<0>(GetParam());
      this->ndofs = get<1>(GetParam());
    }
};

TEST_P(IndexManipulation3D, split_and_linearize)
{
  auto tuple = pfasst::split_index<3>(index, ndofs);
  size_t i = pfasst::linearized_index(tuple, ndofs);

  EXPECT_THAT(i, Eq(index));
}

INSTANTIATE_TEST_CASE_P(IndexManipulation3DTests,
                        IndexManipulation3D,
                        ::testing::Values(
                          make_tuple<size_t, size_t>(4, 3),
                          make_tuple<size_t, size_t>(3, 3),
                          make_tuple<size_t, size_t>(1, 3),
                          make_tuple<size_t, size_t>(5, 3),
                          make_tuple<size_t, size_t>(8, 3)
                        )
                       );


class ConfigModderTest
  : public ::testing::Test
{};

TEST_F(ConfigModderTest, makes_backup_of_initial_config)
{
  ASSERT_FALSE(pfasst::config::has_value("num_iters"));
  {
    ConfigModder config_backup;
    pfasst::config::options::update_value<size_t>("num_iters", 3);
    EXPECT_THAT(pfasst::config::get_value<size_t>("num_iters"), Eq(3));
  }
  EXPECT_FALSE(pfasst::config::has_value("num_iters"));
}

TEST_F(ConfigModderTest, provides_temporary_config_updated)
{
  ASSERT_FALSE(pfasst::config::has_value("num_iters"));
  {
    ConfigModder config_backup;
    config_backup.update<size_t>("num_iters", 3);
    EXPECT_THAT(pfasst::config::get_value<size_t>("num_iters"), Eq(3));
  }
  EXPECT_FALSE(pfasst::config::has_value("num_iters"));
}


TEST_MAIN()
