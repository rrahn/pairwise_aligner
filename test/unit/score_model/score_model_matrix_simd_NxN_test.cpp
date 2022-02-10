// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <array>
#include <seqan3/std/algorithm>
#include <numeric>
#include <tuple>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/score_model/score_model_matrix_simd_NxN.hpp>
#include <pairwise_aligner/simd/simd_score_type.hpp>

#include <pairwise_aligner/score_model/substitution_matrix.hpp>

namespace pa = seqan::pairwise_aligner;

template <typename test_param_t>
struct score_model_matrix_simd_test : public testing::Test
{
    using scalar_score_t = std::tuple_element_t<0, test_param_t>;
    static constexpr size_t strip_width_v = std::tuple_element_t<1, test_param_t>::value;

    using matrix_t = std::remove_cvref_t<decltype(pa::blosum62_standard<scalar_score_t>)>;

    static constexpr size_t dimension = std::tuple_size_v<matrix_t>;
    static constexpr size_t matrix_size = (dimension * (dimension + 1)) / 2;
    using score_model_t = pa::score_model_matrix_simd_NxN<pa::simd_score<scalar_score_t, pa::detail::max_simd_size>,
                                                          dimension>;
    using simd_score_t = typename score_model_t::score_type;
    using simd_rank_t = typename score_model_t::index_type;
    using simd_offset_t = std::invoke_result_t<pa::offset_transform, simd_rank_t>;

    static_assert(simd_score_t::size == simd_rank_t::size);

    static constexpr int8_t simd_size = simd_score_t::size;

    score_model_t matrix;
    pa::offset_transform _offset_fn{dimension, matrix_size};

    void SetUp() override
    {
        std::array<std::array<scalar_score_t, dimension>, dimension> tmp{};

        for (size_t i = 0; i < dimension; ++i)
            std::ranges::copy(pa::blosum62_standard<scalar_score_t>[i].second, tmp[i].data());

        matrix = score_model_t{tmp}; // construct the scoring matrix scheme.
    }

    template <typename simd_rank_t>
    simd_offset_t to_offsets(simd_rank_t const & rank) const noexcept
    {
        return _offset_fn(rank);
    }
};

using scalar_score_types = ::testing::Types<
    std::pair<int8_t, std::integral_constant<size_t, 1>>,
    std::pair<int8_t, std::integral_constant<size_t, 4>>,
    std::pair<int8_t, std::integral_constant<size_t, 8>>//,
    // std::pair<int16_t, std::integral_constant<size_t, 1>>,
    // std::pair<int16_t, std::integral_constant<size_t, 4>>,
    // std::pair<int16_t, std::integral_constant<size_t, 8>>,
    // std::pair<int32_t, std::integral_constant<size_t, 1>>,
    // std::pair<int32_t, std::integral_constant<size_t, 4>>,
    // std::pair<int32_t, std::integral_constant<size_t, 8>>,
    // std::pair<int64_t, std::integral_constant<size_t, 1>>,
    // std::pair<int64_t, std::integral_constant<size_t, 4>>,
    // std::pair<int64_t, std::integral_constant<size_t, 8>>
>;

TYPED_TEST_SUITE(score_model_matrix_simd_test, scalar_score_types);

// ----------------------------------------------------------------------------
// Test Cases
// ----------------------------------------------------------------------------

TYPED_TEST(score_model_matrix_simd_test, access_same)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using simd_rank_t = typename TestFixture::simd_rank_t;
    using scalar_rank_t = typename simd_rank_t::value_type;
    using scalar_score_t = typename TestFixture::scalar_score_t;
    using simd_offset_t = typename TestFixture::simd_offset_t;

    // Move over ranks!
    for (size_t c = 0; c < TestFixture::dimension; ++c)
    {
        simd_rank_t column_rank{static_cast<scalar_rank_t>(c)};

        // We have to get the profile first, for a sequence char:
        std::array<simd_offset_t, TestFixture::strip_width_v> symbols{};
        symbols.fill(this->to_offsets(column_rank));
        auto profile = this->matrix.initialise_profile(symbols, pa::strip_width<TestFixture::strip_width_v>);

        // Receive score for every element
        for (size_t r = 0; r < TestFixture::dimension; ++r)
        {
            simd_rank_t row_rank{static_cast<scalar_rank_t>(r)};
            auto score_array = profile.scores_for(this->to_offsets(row_rank));

            std::ranges::for_each(score_array, [&] (simd_score_t const & scores) {
                for (int8_t k = 0; k < TestFixture::simd_size; ++k) {
                    EXPECT_EQ((int8_t)scores[k], pa::blosum62_standard<scalar_score_t>[r].second[c])
                              << "k = " << (int32_t) k << " "
                              << "rank_c = " << (int32_t) c << " "
                              << "rank_r = " << (int32_t) r << "\n";
                }
            });
        }
    }
}

TYPED_TEST(score_model_matrix_simd_test, access_rotate)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using simd_rank_t = typename TestFixture::simd_rank_t;
    using scalar_score_t = typename TestFixture::scalar_score_t;
    using simd_offset_t = typename TestFixture::simd_offset_t;

    std::vector<int8_t> symbols_col{};
    symbols_col.resize(TestFixture::dimension);

    std::vector<int8_t> symbols_row{};
    symbols_row.resize(TestFixture::dimension);

    // Rotate the symbols_col in every step.
    for (size_t cycle_col = 0; cycle_col < symbols_col.size(); ++cycle_col)
    // size_t cycle_col = 0;
    {
        std::iota(symbols_col.begin(), symbols_col.end(), 0);
        std::ranges::rotate(symbols_col, symbols_col.end() - cycle_col);

        simd_rank_t column_rank{};
        for (int8_t k = 0; k < TestFixture::simd_size; ++k)
            column_rank[k] = symbols_col[k % symbols_col.size()]; // cycle around symbols

        // We have to get the profile first, for a sequence char:
        std::array<simd_offset_t, TestFixture::strip_width_v> symbols{};
        symbols.fill(this->to_offsets(column_rank));
        auto profile = this->matrix.initialise_profile(symbols, pa::strip_width<TestFixture::strip_width_v>);

        for (size_t cycle_row = 0; cycle_row < symbols_row.size(); ++cycle_row)
        {
            std::iota(symbols_row.begin(), symbols_row.end(), 0);
            std::ranges::rotate(symbols_row, symbols_row.end() - cycle_row);
            simd_rank_t row_rank{};
            for (int8_t k = 0; k < TestFixture::simd_size; ++k)
                row_rank[k] = symbols_row[k % symbols_row.size()]; // cycle around symbols

            // Which index are we now comapring for every element?
            auto score_array = profile.scores_for(this->to_offsets(row_rank));

            std::ranges::for_each(score_array, [&] (simd_score_t const & scores) {
                for (int8_t k = 0; k < TestFixture::simd_size; ++k) {
                    EXPECT_EQ(scores[k],
                            pa::blosum62_standard<scalar_score_t>[row_rank[k]].second[column_rank[k]])
                            << "k = " << (int32_t) k << " "
                            << "rank_c = " << (int32_t) column_rank[k] << " "
                            << "rank_r = " << (int32_t) row_rank[k] << " "
                            << "cycle_col = " << (int32_t) cycle_col << " "
                            << "cycle_row = " << (int32_t) cycle_row << "\n";
                }
            });
        }
    }
}
