#include "fixtures/test_helpers.hpp"

#include <algorithm>
#include <memory>
#include <type_traits>
#include <vector>
using namespace std;

#include <pfasst/globals.hpp>
#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/eigen3_vector.hpp>
typedef pfasst::eigen3_encap_traits<double, double, 1>      EigenVectorEncapTrait;
typedef pfasst::encap::Encapsulation<EigenVectorEncapTrait> EigenVectorEncapsulation;

#include "comm/mocks.hpp"
typedef CommMock                                       CommunicatorType;


typedef ::testing::Types<EigenVectorEncapsulation> EncapTypes;
INSTANTIATE_TYPED_TEST_CASE_P(EigenVectorEncap, Concepts, EncapTypes);


TEST(Construction, empty_constructible)
{
  EigenVectorEncapsulation vec;
  EXPECT_THAT(vec.get_data().cols(), Eq(1));
  EXPECT_THAT(vec.get_data().rows(), Eq(0));
  EXPECT_THAT(vec.get_data(), Eq(EigenVector<double>{}));
}

TEST(Construction, data_constructible)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  EigenVectorEncapsulation vec(data);

  EXPECT_THAT(vec.get_data().rows(), Eq(3));
  EXPECT_THAT(vec.get_data().cols(), Eq(1));
  EXPECT_THAT(vec.get_data(), Eq(data));
}

TEST(DataAccession, assignable)
{
  EigenVectorEncapsulation vec;
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;

  vec = data;
  EXPECT_THAT(vec.get_data().rows(), Eq(3));
  EXPECT_THAT(vec.get_data().cols(), Eq(1));
  EXPECT_THAT(vec.get_data(), Eq(data));
}


TEST(Operation, zeroing_out)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  EigenVectorEncapsulation x(data);
  EXPECT_TRUE((x.get_data().array() == data.array()).all());

  x.zero();
  EXPECT_TRUE((x.get_data().array() == EigenVector<double>::Zero(3).array()).all());
}

TEST(Operation, in_place_axpy)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  EigenVectorEncapsulation vec_x(data);
  shared_ptr<EigenVectorEncapsulation> vec_y = \
    make_shared<EigenVectorEncapsulation>(EigenVector<double>::Ones(3));

  vec_x.scaled_add(0.5, vec_y);
  EigenVector<double> expected(3);
  expected << 1.5, 2.5, 3.5;
  EXPECT_TRUE((vec_x.get_data().array() == expected.array()).all());
}

TEST(Operation, global_axpy)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  shared_ptr<EigenVectorEncapsulation> vec_x = make_shared<EigenVectorEncapsulation>(data);
  shared_ptr<EigenVectorEncapsulation> vec_y = 
    make_shared<EigenVectorEncapsulation>(EigenVector<double>::Ones(3));

  auto result = pfasst::encap::axpy(0.5, vec_x, vec_y);

  EigenVector<double> expected(3);
  expected << 1.5, 2.0, 2.5;
  EXPECT_TRUE((result->get_data().array() == expected.array()).all());
}

TEST(Operation, norm0_as_member)
{
  EigenVector<double> data(3);
  data << 1.0, -4.0, 3.0;
  EigenVectorEncapsulation vec_x(data);
  EXPECT_THAT(vec_x.norm0(), Eq(4.0));
}

TEST(Operation, global_norm0)
{
  EigenVector<double> data(3);
  data << 1.0, -4.0, 3.0;
  shared_ptr<EigenVectorEncapsulation> vec_x = \
    make_shared<EigenVectorEncapsulation>(data);
  EXPECT_THAT(norm0(vec_x), Eq(4.0));
}


TEST(MatrixApplication, identity)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  vector<shared_ptr<EigenVectorEncapsulation>> vec(3);
  generate(vec.begin(), vec.end(),
           [data]() { return make_shared<EigenVectorEncapsulation>(data); });
  for_each(vec.cbegin(), vec.cend(), [data](const shared_ptr<EigenVectorEncapsulation>& xi) {
    EXPECT_TRUE((xi->get_data().array() == data.array()).all());
  });
  Matrix<double> mat = Matrix<double>::Identity(3, 3);

  auto result_mat_mul_vec = pfasst::encap::mat_mul_vec(1.0, mat, vec);
  for_each(result_mat_mul_vec.cbegin(), result_mat_mul_vec.cend(),
           [data](const shared_ptr<EigenVectorEncapsulation>& xi) {
             EXPECT_TRUE((xi->get_data().array() == data.array()).all());
           });

  vector<shared_ptr<EigenVectorEncapsulation>> result_mat_apply(result_mat_mul_vec);
  pfasst::encap::mat_apply(result_mat_apply, 1.0, mat, vec, true);
  for_each(result_mat_apply.cbegin(), result_mat_apply.cend(),
           [data](const shared_ptr<EigenVectorEncapsulation>& xi) {
             EXPECT_TRUE((xi->get_data().array() == data.array()).all());
           });
}

TEST(MatrixApplication, zero_matrix)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  vector<shared_ptr<EigenVectorEncapsulation>> vec(3);
  generate(vec.begin(), vec.end(),
           [data]() { return make_shared<EigenVectorEncapsulation>(data); });
  for_each(vec.cbegin(), vec.cend(), [data](const shared_ptr<EigenVectorEncapsulation>& xi) {
    EXPECT_TRUE((xi->get_data().array() == data.array()).all());
  });
  Matrix<double> mat = Matrix<double>::Zero(3, 3);

  auto result_mat_mul_vec = pfasst::encap::mat_mul_vec(1.0, mat, vec);
  for_each(result_mat_mul_vec.cbegin(), result_mat_mul_vec.cend(),
           [&](const shared_ptr<EigenVectorEncapsulation>& xi) {
             EXPECT_TRUE((xi->get_data().array() == EigenVector<double>::Zero(3).array()).all());
           });

  vector<shared_ptr<EigenVectorEncapsulation>> result_mat_apply(result_mat_mul_vec);
  pfasst::encap::mat_apply(result_mat_apply, 1.0, mat, vec, true);
  for_each(result_mat_apply.cbegin(), result_mat_apply.cend(),
           [&](const shared_ptr<EigenVectorEncapsulation>& xi) {
             EXPECT_TRUE((xi->get_data().array() == EigenVector<double>::Zero(3).array()).all());
           });
}

TEST(MatrixApplication, all_ones)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  vector<shared_ptr<EigenVectorEncapsulation>> vec(3);
  generate(vec.begin(), vec.end(),
           [data]() { return make_shared<EigenVectorEncapsulation>(data); });
  for_each(vec.cbegin(), vec.cend(), [data](const shared_ptr<EigenVectorEncapsulation>& xi) {
    EXPECT_TRUE((xi->get_data().array() == data.array()).all());
  });
  Matrix<double> mat = Matrix<double>::Ones(3, 3);

  auto result_mat_mul_vec = pfasst::encap::mat_mul_vec(1.0, mat, vec);

  EigenVector<double> expected(3);
  expected << 3.0, 6.0, 9.0;
  for_each(result_mat_mul_vec.cbegin(), result_mat_mul_vec.cend(),
           [&](const shared_ptr<EigenVectorEncapsulation>& xi) {
             EXPECT_TRUE((xi->get_data().array() == expected.array()).all());
           });

  vector<shared_ptr<EigenVectorEncapsulation>> result_mat_apply(result_mat_mul_vec);
  pfasst::encap::mat_apply(result_mat_apply, 1.0, mat, vec, true);
  for_each(result_mat_apply.cbegin(), result_mat_apply.cend(),
           [&](const shared_ptr<EigenVectorEncapsulation>& xi) {
             EXPECT_TRUE((xi->get_data().array() == expected.array()).all());
           });
}


TEST(Communication, sending)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  shared_ptr<CommunicatorType> comm = make_shared<CommunicatorType>();
  shared_ptr<EigenVectorEncapsulation> vec = \
    make_shared<EigenVectorEncapsulation>(data);
  vec->send(comm, 1, 0, true);

  vec->send(comm, 1, 0, false);
}

TEST(Communication, receiving)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  shared_ptr<CommunicatorType> comm = make_shared<CommunicatorType>();
  shared_ptr<EigenVectorEncapsulation> vec = \
    make_shared<EigenVectorEncapsulation>(data);
  vec->recv(comm, 1, 0, true);

  vec->recv(comm, 1, 0, false);
}

TEST(Communication, broadcasting)
{
  EigenVector<double> data(3);
  data << 1.0, 2.0, 3.0;
  shared_ptr<CommunicatorType> comm = make_shared<CommunicatorType>();
  shared_ptr<EigenVectorEncapsulation> vec = \
    make_shared<EigenVectorEncapsulation>(data);
  vec->bcast(comm, 0);
}


TEST(Factory, predefine_size)
{
  pfasst::encap::EncapsulationFactory<EigenVectorEncapTrait> null_factory;
  EXPECT_THAT(null_factory.size(), Eq(0));

  pfasst::encap::EncapsulationFactory<EigenVectorEncapTrait> sized_factory(3);
  EXPECT_THAT(sized_factory.size(), Eq(3));

  null_factory.set_size(3);
  EXPECT_THAT(null_factory.size(), Eq(3));
}

TEST(Factory, create_vector_encap)
{
  pfasst::encap::EncapsulationFactory<EigenVectorEncapTrait> factory(3);
  auto encap = factory.create();
  EXPECT_THAT(encap->get_data().size(), Eq(3));

  factory.set_size(5);
  EXPECT_THAT(factory.size(), Eq(5));
  auto encap_5 = factory.create();
  EXPECT_THAT(encap_5->get_data().cols(), Eq(1));
  EXPECT_THAT(encap_5->get_data().rows(), Eq(5));
}

TEST_MAIN()
