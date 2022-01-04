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
#include <seqan3/utility/container/aligned_allocator.hpp>
#include <seqan3/utility/simd/views/to_simd.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
template <typename dp_vector_t, typename simd_t>
class dp_vector_bulk
{
private:
    dp_vector_t _dp_vector{};

public:

    using range_type = typename dp_vector_t::range_type;
    using value_type = typename dp_vector_t::value_type;
    using reference = typename dp_vector_t::reference;
    using const_reference = typename dp_vector_t::const_reference;

    reference operator[](size_t const pos) noexcept(noexcept(_dp_vector[pos]))
    {
        return _dp_vector[pos];
    }

    const_reference operator[](size_t const pos) const noexcept(noexcept(_dp_vector[pos]))
    {
        return _dp_vector[pos];
    }

    constexpr size_t size() const noexcept
    {
        return _dp_vector.size();
    }

    decltype(auto) range() noexcept
    {
        return _dp_vector.range();
    }

    decltype(auto) range() const noexcept
    {
        return _dp_vector.range();
    }

    template <std::ranges::forward_range sequence_t, typename initialisation_strategy_t>
    auto initialise(sequence_t && sequence, initialisation_strategy_t && init_strategy)
    {
        size_t const sequence_count = std::ranges::distance(sequence);
        size_t max_sequence_size = 0;

        using optional_t = std::ranges::range_value_t<sequence_t>;
        using seq_view_t = typename optional_t::value_type;

        std::vector<seq_view_t, seqan3::aligned_allocator<seq_view_t, alignof(simd_t)>> sequences{};
        sequences.reserve(sequence_count);

        std::ranges::for_each(sequence, [&] (auto && optional_sequence)
        {
            if (optional_sequence.has_value()) {
                max_sequence_size = std::max<size_t>(max_sequence_size, std::ranges::distance(*optional_sequence));
                sequences.push_back(*optional_sequence);
            }
        });

        using score_t = typename simd_t::value_type;
        using simd_score_t = typename simd_t::simd_type::value_type;

        constexpr size_t bit_count = sizeof(score_t) * 8;
        constexpr score_t padding_mask = static_cast<score_t>(1 << (bit_count - 1));

        std::vector<simd_t, seqan3::aligned_allocator<simd_t, alignof(simd_t)>> simd_sequence{};
        simd_sequence.reserve(max_sequence_size);

        auto simd_view = sequences | seqan3::views::to_simd<simd_score_t>(padding_mask);

        for (auto && simd_vector_chunk : simd_view) {
            for (auto && simd_vector : simd_vector_chunk) {
                simd_sequence.emplace_back(std::move(simd_vector));
            }
        }
        return _dp_vector.initialise(std::move(simd_sequence), std::forward<initialisation_strategy_t>(init_strategy));
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
