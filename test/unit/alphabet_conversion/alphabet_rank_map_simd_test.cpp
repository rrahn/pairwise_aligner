// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string_view>

#include <pairwise_aligner/alphabet_conversion/alphabet_rank_map_simd.hpp>
#include <pairwise_aligner/simd/simd_score_type.hpp>

#include "../fixture/fixture_base.hpp"

namespace alphabet_conversion::test {

// ----------------------------------------------------------------------------
// Values for test
// ----------------------------------------------------------------------------

struct values
{
    std::string_view symbol_list;
};

// ----------------------------------------------------------------------------
// Define test
// ----------------------------------------------------------------------------

template <typename fixture_t>
struct fixture : public fixture_t
{
    using simd_t = seqan::pairwise_aligner::simd_score<int8_t>;

    static constexpr size_t simd_size = simd_t::size;

    using alphabet_rank_map_t = seqan::pairwise_aligner::alphabet_rank_map_simd<simd_t>;
    using key_t = typename alphabet_rank_map_t::key_type;
    using value_t = typename alphabet_rank_map_t::value_type;

};

} // namespec alphabet_conversion::test

template <typename fixture_t>
using test_suite = alphabet_conversion::test::fixture<fixture_t>;

TYPED_TEST_SUITE_P(test_suite);

// ----------------------------------------------------------------------------
// Define test cases
// ----------------------------------------------------------------------------

// Need multi value!
// It is a value consisting of multiple elements.
// We choose by max of symbol list size and simd
TYPED_TEST_P(test_suite, smallest_rank)
{
    try {
        typename TestFixture::alphabet_rank_map_t rank_map{this->GetParam().symbol_list};

        typename TestFixture::key_t key{this->GetParam().symbol_list.front()};
        typename TestFixture::value_t rank = rank_map[key];

        for (size_t i = 0; i < TestFixture::simd_size; ++i)
            EXPECT_EQ(rank[i], 0);
    } catch(...) {
        EXPECT_THROW(std::rethrow_exception(std::current_exception()), std::invalid_argument);
    }
}

TYPED_TEST_P(test_suite, largest_rank)
{
    try {
        typename TestFixture::alphabet_rank_map_t rank_map{this->GetParam().symbol_list};

        typename TestFixture::key_t key{this->GetParam().symbol_list.back()};
        typename TestFixture::value_t rank = rank_map[key];

        for (size_t i = 0; i < TestFixture::simd_size; ++i)
            EXPECT_EQ(rank[i], this->GetParam().symbol_list.size() - 1);
    } catch(...) {
        EXPECT_THROW(std::rethrow_exception(std::current_exception()), std::invalid_argument);
    }
}

TYPED_TEST_P(test_suite, rank_in_order)
{
    try {
        auto symbol_list = this->GetParam().symbol_list;
        typename TestFixture::alphabet_rank_map_t rank_map{symbol_list};

        for (size_t j = 0; j < symbol_list.size(); j += TestFixture::simd_size) {
            auto slice = symbol_list.substr(j, j + TestFixture::simd_size);
            size_t const slice_size = slice.size();

            typename TestFixture::key_t key{};
            for (size_t i = 0; i < TestFixture::simd_size; ++i)
                key[i] = slice[i % slice_size];

            typename TestFixture::value_t rank = rank_map[key];

            for (size_t i = 0; i < TestFixture::simd_size; ++i)
                EXPECT_EQ(rank[i], j + (i % slice_size));
        }

    } catch(...) {
        EXPECT_THROW(std::rethrow_exception(std::current_exception()), std::invalid_argument);
    }
}

TYPED_TEST_P(test_suite, rank_in_reversed_order)
{
    try {
        auto symbol_list = this->GetParam().symbol_list;
        typename TestFixture::alphabet_rank_map_t rank_map{symbol_list};

        for (size_t j = 0; j < symbol_list.size(); j += TestFixture::simd_size) {
            auto slice = symbol_list.substr(j, j + TestFixture::simd_size);
            size_t const slice_size = slice.size();

            typename TestFixture::key_t key{};
            for (size_t i = 0; i < TestFixture::simd_size; ++i)
                key[TestFixture::simd_size - 1 - i] = slice[i % slice_size];

            typename TestFixture::value_t rank = rank_map[key];

            for (size_t i = 0; i < TestFixture::simd_size; ++i)
                EXPECT_EQ(rank[i], j + ((TestFixture::simd_size - 1 - i) % slice_size));
        }
    } catch (...) {
        EXPECT_THROW(std::rethrow_exception(std::current_exception()), std::invalid_argument);
    }
}

// ----------------------------------------------------------------------------
// Register test cases
// ----------------------------------------------------------------------------

REGISTER_TYPED_TEST_SUITE_P(test_suite, smallest_rank, largest_rank, rank_in_order, rank_in_reversed_order);

// ----------------------------------------------------------------------------
// Helper macro to define test values
// ----------------------------------------------------------------------------

#define DEFINE_TEST_VALUES(name, ...)  \
static auto name = [](){ return alphabet_conversion::test::values{__VA_ARGS__}; }();

// ----------------------------------------------------------------------------
// Define test values
// ----------------------------------------------------------------------------

namespace aligner = seqan::pairwise_aligner;

DEFINE_TEST_VALUES(empty_symbol_list,
    .symbol_list{""}
)

DEFINE_TEST_VALUES(symbol_list_with_single_element,
    .symbol_list{"X"}
)

DEFINE_TEST_VALUES(symbol_list_dna5,
    .symbol_list{"ACGTN"}
)

DEFINE_TEST_VALUES(symbol_list_aa20,
    .symbol_list{"ACDEFGHIKLMNPQRSTVWY"}
)

DEFINE_TEST_VALUES(symbol_list_aa27,
    .symbol_list{"ABCDEFGHIJKLMNOPQRSTUVWXYZ*"}
)

DEFINE_TEST_VALUES(symbol_list_printable_char,
    .symbol_list{" !\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"}
)

using test_types =
    ::testing::Types<
        pairwise_aligner::test::fixture<&empty_symbol_list>,
        pairwise_aligner::test::fixture<&symbol_list_with_single_element>,
        pairwise_aligner::test::fixture<&symbol_list_dna5>,
        pairwise_aligner::test::fixture<&symbol_list_aa20>,
        pairwise_aligner::test::fixture<&symbol_list_aa27>,
        pairwise_aligner::test::fixture<&symbol_list_printable_char>
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_conversion_simd_test, test_suite, test_types,);
