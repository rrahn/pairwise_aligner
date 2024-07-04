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

#include <pairwise_aligner/score_model/score_model_matrix_simd_NxN.hpp>
#include <pairwise_aligner/simd/simd_score_type.hpp>

#include <pairwise_aligner/score_model/substitution_matrix.hpp>

namespace pa = seqan::pairwise_aligner;

template <typename scalar_score_t>
struct score_model_matrix_simd_test : public testing::Test
{
    using matrix_t = std::remove_cvref_t<decltype(pa::blosum62_standard<scalar_score_t>)>;

    static constexpr size_t dimension = std::tuple_size_v<matrix_t>;
    static constexpr size_t matrix_size = (dimension * (dimension + 1)) / 2;

    using simd_score_t = pa::simd_score<scalar_score_t>;
    using simd_index_t = pa::detail::make_unsigned_t<simd_score_t>;
    using score_model_t = pa::score_model_matrix_simd_NxN<simd_score_t, simd_index_t, dimension>;

    static_assert(simd_score_t::size == simd_index_t::size);

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

    template <typename simd_index_t>
    auto to_offsets(simd_index_t const & rank) const noexcept
    {
        return _offset_fn(rank);
    }
};

using scalar_score_types = ::testing::Types<int8_t
                                            // int16_t,
                                            // int32_t,
                                            // int64_t
>;

TYPED_TEST_SUITE(score_model_matrix_simd_test, scalar_score_types);

// ----------------------------------------------------------------------------
// Test Cases
// ----------------------------------------------------------------------------

TYPED_TEST(score_model_matrix_simd_test, access_same)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_score_t = typename simd_score_t::value_type;

    using simd_index_t = typename TestFixture::simd_index_t;
    using scalar_rank_t = typename simd_index_t::value_type;

    // Move over ranks!
    // for (scalar_rank_t c = 0; c < static_cast<scalar_rank_t>(TestFixture::dimension); ++c)
    scalar_rank_t c = 0;
    {
        simd_index_t column_rank{c};

        // Receive score for every element
        // for (scalar_rank_t r = 0; r < static_cast<scalar_rank_t>(TestFixture::dimension); ++r)
        scalar_rank_t r = 0;
        {
            simd_index_t row_rank{r};
            simd_score_t scores = this->matrix.score(simd_score_t{},
                                                     this->to_offsets(column_rank),
                                                     this->to_offsets(row_rank));

            for (size_t k = 0; k < simd_score_t::size; ++k) {
                EXPECT_EQ(scores[k], pa::blosum62_standard<scalar_score_t>[r].second[c])
                            << "k = " << (int32_t) k << " "
                            << "rank_c = " << (int32_t) c << " "
                            << "rank_r = " << (int32_t) r << "\n";
            }
        }
    }
}

// TYPED_TEST(score_model_matrix_simd_test, access_rotate)
// {
//     using simd_score_t = typename TestFixture::simd_score_t;
//     using scalar_score_t = typename simd_score_t::value_type;

//     using simd_index_t = typename TestFixture::simd_index_t;
//     using scalar_rank_t = typename simd_index_t::value_type;

//     std::vector<scalar_rank_t> symbols_col{};
//     symbols_col.resize(TestFixture::dimension);

//     std::vector<scalar_rank_t> symbols_row{};
//     symbols_row.resize(TestFixture::dimension);

//     // Rotate the symbols_col in every step.
//     for (size_t cycle_col = 0; cycle_col < symbols_col.size(); ++cycle_col)
//     {
//         std::iota(symbols_col.begin(), symbols_col.end(), 0);
//         std::ranges::rotate(symbols_col, symbols_col.end() - cycle_col);

//         simd_index_t column_rank{};
//         for (size_t k = 0; k < simd_score_t::size; ++k)
//             column_rank[k] = symbols_col[k % symbols_col.size()]; // cycle around symbols

//         for (size_t cycle_row = 0; cycle_row < symbols_row.size(); ++cycle_row)
//         {
//             std::iota(symbols_row.begin(), symbols_row.end(), 0);
//             std::ranges::rotate(symbols_row, symbols_row.end() - cycle_row);
//             simd_index_t row_rank{};
//             for (size_t k = 0; k < simd_score_t::size; ++k)
//                 row_rank[k] = symbols_row[k % symbols_row.size()]; // cycle around symbols

//             auto scores = this->matrix.score(simd_score_t{}, this->to_offsets(column_rank), this->to_offsets(row_rank));

//             for (size_t k = 0; k < simd_score_t::size; ++k) {
//                 EXPECT_EQ(scores[k], pa::blosum62_standard<scalar_score_t>[row_rank[k]].second[column_rank[k]])
//                         << "k = " << (int32_t) k << " "
//                         << "rank_c = " << (int32_t) column_rank[k] << " "
//                         << "rank_r = " << (int32_t) row_rank[k] << " "
//                         << "cycle_col = " << (int32_t) cycle_col << " "
//                         << "cycle_row = " << (int32_t) cycle_row << "\n";
//             }
//         }
//     }
// }
