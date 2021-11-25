// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the public pairwise aligner interface.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <iostream>

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <tuple>
#include <vector>

#include <pairwise_aligner/simd_score_type.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

using scalar_score_t = int16_t;
using simd_score_t = simd_score<scalar_score_t>;

template <typename dp_cell_t, typename dp_vector_t = std::vector<dp_cell_t>>
class intermediate_dp_vector
{
private:
    dp_vector_t _dp_vector;

public:

    using value_type = std::ranges::range_value_t<dp_vector_t>;
    using reference = std::ranges::range_reference_t<dp_vector_t>;
    using const_reference = std::ranges::range_reference_t<dp_vector_t const>;

    reference operator[](size_t const pos) noexcept(noexcept(_dp_vector[pos]))
    {
        return _dp_vector[pos];
    }

    const_reference operator[](size_t const pos) const noexcept(noexcept(_dp_vector[pos]))
    {
        return _dp_vector[pos];
    }

    template <std::ranges::viewable_range sequence_t,
              typename initialisation_strategy_t>
        requires std::ranges::forward_range<sequence_t>
    sequence_t initialise(sequence_t && sequence, initialisation_strategy_t && init_strategy)
    {
        size_t const sequence_size = std::ranges::distance(sequence);
        _dp_vector.resize(sequence_size + 1);

        std::ranges::for_each(_dp_vector, init_strategy);

        return sequence;
    }
};

// ----------------------------------------------------------------------------
// simd_intermediate_dp_vector
// ----------------------------------------------------------------------------

template <typename dp_cell_t, typename dp_vector_t = std::vector<dp_cell_t>>
class simd_intermediate_dp_vector
{
private:
    intermediate_dp_vector<dp_cell_t, dp_vector_t> _underlying_dp_vector;

public:

    using value_type = std::ranges::range_value_t<dp_vector_t>;
    using reference = std::ranges::range_reference_t<dp_vector_t>;
    using const_reference = std::ranges::range_reference_t<dp_vector_t const>;

    reference operator[](size_t const pos) noexcept(noexcept(_underlying_dp_vector[pos]))
    {
        return _underlying_dp_vector[pos];
    }

    const_reference operator[](size_t const pos) const noexcept(noexcept(_underlying_dp_vector[pos]))
    {
        return _underlying_dp_vector[pos];
    }

    template <std::ranges::forward_range sequence_t,
              typename initialisation_strategy_t>
        requires (std::ranges::forward_range<std::ranges::range_reference_t<sequence_t>> &&
                  std::ranges::viewable_range<std::ranges::range_reference_t<sequence_t>>)
    auto initialise(sequence_t && sequence, initialisation_strategy_t && init_strategy)
    {
        size_t const sequence_count = std::ranges::distance(sequence);
        size_t max_sequence_size = 0;

        std::ranges::for_each(sequence, [&] (auto && inner_sequence)
        {
            max_sequence_size = std::max<size_t>(max_sequence_size, std::ranges::distance(inner_sequence));
        });

        std::vector<simd_score_t> simd_sequence{};
        simd_sequence.resize(max_sequence_size);

        _underlying_dp_vector.initialise(simd_sequence, init_strategy);

        for (size_t j = 0; j < max_sequence_size; ++j)
            for (size_t i = 0; i < sequence_count; ++i)
                simd_sequence[j][i] = sequence[i][j];

        return simd_sequence;
    }
};

template <typename cell_t>
using dp_column_vector = intermediate_dp_vector<cell_t>;

template <typename cell_t>
using dp_row_vector = intermediate_dp_vector<cell_t>;

template <typename derived_t>
class pairwise_aligner
{
private:

    friend derived_t;

    pairwise_aligner() = default;

public:

    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
        requires (std::ranges::viewable_range<sequence1_t> &&
                  std::ranges::viewable_range<sequence2_t>)
    auto compute(sequence1_t && sequence1, sequence2_t && sequence2)
        -> int32_t
    {
        return compute(std::forward<sequence1_t>(sequence1),
                       std::forward<sequence2_t>(sequence2),
                       dp_column_vector<affine_dp_cell_t>{},
                       dp_row_vector<affine_dp_cell_t>{});
    }

    template <std::ranges::forward_range sequence1_t,
              std::ranges::forward_range sequence2_t,
              typename dp_column_vector_t,
              typename dp_row_vector_t>
        requires (std::ranges::viewable_range<sequence1_t> &&
                  std::ranges::viewable_range<sequence2_t>)
    auto compute(sequence1_t && sequence1,
                 sequence2_t && sequence2,
                 dp_column_vector_t first_dp_column,
                 dp_row_vector_t first_dp_row)
        -> int32_t
    {
        auto && [last_dp_column, last_dp_row] = run(std::forward<sequence1_t>(sequence1),
                                                    std::forward<sequence2_t>(sequence2),
                                                    std::move(first_dp_column),
                                                    std::move(first_dp_row));
        return get<0>(last_dp_column[std::ranges::size(sequence1)]);
    }

private:

    template <std::ranges::forward_range sequence1_t,
              std::ranges::forward_range sequence2_t,
              typename dp_column_vector_t,
              typename dp_row_vector_t>
        requires (std::ranges::viewable_range<sequence1_t> &&
                  std::ranges::viewable_range<sequence2_t>)
    auto run(sequence1_t && sequence1,
             sequence2_t && sequence2,
             dp_column_vector_t dp_column_vector,
             dp_row_vector_t dp_row_vector)
    {
        // ----------------------------------------------------------------------------
        // Initialisation
        // ----------------------------------------------------------------------------

        auto transformed_seq1 = as_derived().initialise_column_vector(std::forward<sequence1_t>(sequence1),
                                                                      dp_column_vector);
        auto transformed_seq2 = as_derived().initialise_row_vector(std::forward<sequence2_t>(sequence2),
                                                                   dp_row_vector);

        // auto print_col = [&] (auto col)
        // {
        //     for (size_t j = 0; j < std::ranges::size(transformed_seq1) + 1; ++j)
        //     {
        //         std::cout << get<0>(col[j]) << "\t";
        //     }
        //     std::cout << "\n";
        // };

        // std::cout << "\t\t\t";
        // for (char x : transformed_seq2)
        //     std::cout << x << "\t";
        // std::cout << "\n";
        // std::cout << "col 0:\t\t"; print_col(dp_column_vector);
        // int si = 0;

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------

        for (size_t j = 0; j < std::ranges::size(transformed_seq2); ++j)
        {
            auto cache = as_derived().initialise_cache(dp_column_vector[0], dp_row_vector[j+1]);
            dp_column_vector[0] = dp_row_vector[j+1];
            for (size_t i = 0; i < std::ranges::size(transformed_seq1); ++i)
            {
                as_derived().compute_cell(cache,
                                          dp_column_vector[i+1],
                                          transformed_seq1[i],
                                          transformed_seq2[j]);
            }
            // std::cout << "col " << j + 1 <<": " << transformed_seq1[si++] << "\t"; print_col(dp_column_vector);
        }

        return std::pair{std::move(dp_column_vector), std::move(dp_row_vector)};
    }

    derived_t const & as_derived() const noexcept
    {
        return static_cast<derived_t const &>(*this);
    }

    derived_t & as_derived() noexcept
    {
        return static_cast<derived_t &>(*this);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
