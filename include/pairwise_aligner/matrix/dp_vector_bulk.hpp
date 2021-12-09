// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_bulk.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

#include <pairwise_aligner/matrix/dp_vector_single.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
template <typename dp_cell_t, typename dp_vector_t = std::vector<dp_cell_t>>
class dp_vector_bulk
{
private:
    dp_vector_single<dp_cell_t, dp_vector_t> _underlying_dp_vector;

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

    constexpr size_t size() const noexcept
    {
        return _underlying_dp_vector.size();
    }

    template <std::ranges::forward_range sequence_t,
              typename initialisation_strategy_t>
    auto initialise(sequence_t && sequence, initialisation_strategy_t && init_strategy)
    {
        size_t const sequence_count = std::ranges::distance(sequence);
        size_t max_sequence_size = 0;

        std::ranges::for_each(sequence, [&] (auto && optional_sequence)
        {
            if (optional_sequence.has_value())
                max_sequence_size = std::max<size_t>(max_sequence_size, std::ranges::distance(*optional_sequence));
        });

        using simd_score_t = typename value_type::score_type;
        using score_t = typename simd_score_t::value_type;

        constexpr size_t bit_count = sizeof(score_t) * 8;
        constexpr score_t padding_mask = static_cast<score_t>(1 << (bit_count - 1));

        std::vector<simd_score_t> simd_sequence{};
        simd_sequence.resize(max_sequence_size, simd_score_t{padding_mask});

        _underlying_dp_vector.initialise(simd_sequence, init_strategy);

        for (size_t i = 0; i < sequence_count; ++i)
        {
            if (!sequence[i].has_value())
                continue;

            size_t j = 0;
            for (auto it = std::ranges::begin(*sequence[i]); it != std::ranges::end(*sequence[i]); ++it, ++j)
                simd_sequence[j][i] = static_cast<score_t>(*it);
        }

        return simd_sequence;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
