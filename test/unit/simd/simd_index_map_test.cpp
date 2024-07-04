// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <array>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <type_traits>

#include <pairwise_aligner/simd/simd_score_type.hpp>
#include <pairwise_aligner/simd/simd_index_map.hpp>


namespace pa = seqan::pairwise_aligner;

template <typename test_param_t>
struct simd_index_map_test : public testing::Test
{
    using scalar_value_t = std::tuple_element_t<0, test_param_t>;
    using scalar_key_t = std::tuple_element_t<1, test_param_t>;

    void SetUp() override
    {
    }
};

using test_types = ::testing::Types<
    std::pair<int8_t, uint8_t>//,
    // std::pair<int16_t, uint16_t>,
    // std::pair<int32_t, uint32_t>
>;

TYPED_TEST_SUITE(simd_index_map_test, test_types);

// ----------------------------------------------------------------------------
// Test Cases
// ----------------------------------------------------------------------------

TYPED_TEST(simd_index_map_test, select_size_16)
{
    using map_t = pa::simd_index_map<typename TestFixture::scalar_value_t, typename TestFixture::scalar_key_t, 5>;

    std::array data{0, 1, 2, 3, 4};
    map_t map{data};

    using key_t = typename map_t::key_type;
    key_t key{};
    auto values = map[key];

    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], 0);

    // Set all to one
    key += 1;
    values = map[key];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], 1);

    key += 1;
    values = map[key];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], 2);

    key += 1;
    values = map[key];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], 3);

    key += 1;
    values = map[key];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], 4);

    for (size_t i = 0; i < key_t::size; ++i)
        key[i] = i % 5;

    values = map[key];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], i % 5);
}

TYPED_TEST(simd_index_map_test, select_simd_size)
{
    constexpr size_t data_size = pa::detail::max_simd_size;
    std::array<typename TestFixture::scalar_value_t, data_size> data{};
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i;

    using map_t = pa::simd_index_map<typename TestFixture::scalar_value_t, typename TestFixture::scalar_key_t, data_size>;
    map_t map{data};

    using key_t = typename map_t::key_type;
    key_t key{};

    for (size_t i = 0; i < key_t::size; ++i)
        key[i] = i % data_size;

    auto values = map[key];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], i % data_size);

}

TYPED_TEST(simd_index_map_test, select_size_128)
{
    constexpr size_t data_size = 128;
    std::array<typename TestFixture::scalar_value_t, data_size> data{};
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i;

    using map_t = pa::simd_index_map<typename TestFixture::scalar_value_t, typename TestFixture::scalar_key_t, data_size>;
    map_t map{data};

    using key_t = typename map_t::key_type;
    key_t key{};

    for (size_t i = 0; i < key_t::size; ++i)
        key[i] = i % data_size;

    auto values = map[key];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], i % data_size);
}

TYPED_TEST(simd_index_map_test, select_size_256)
{
    constexpr size_t data_size = 256;
    std::array<typename TestFixture::scalar_value_t, data_size> data{};
    for (size_t i = 0; i < data.size(); ++i)
        data[i] = i;

    using map_t = pa::simd_index_map<typename TestFixture::scalar_value_t, typename TestFixture::scalar_key_t, data_size>;
    map_t map{data};

    using key_t = typename map_t::key_type;
    key_t key1{};

    for (size_t i = 0; i < key_t::size; ++i)
        key1[i] = i % data_size;

    key_t key2 = key1 + pa::detail::max_simd_size;

    auto values = map[key1];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], i % data_size);

    values = map[key2];
    for (size_t i = 0; i < key_t::size; ++i)
        EXPECT_EQ(values[i], i % data_size + pa::detail::max_simd_size);
}
