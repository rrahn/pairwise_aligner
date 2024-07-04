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

#include <algorithm>
#include <ranges>
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
    using scalar_t = typename simd_t::value_type;
    using native_simd_t = typename simd_t::simd_type::value_type;

    dp_vector_t _dp_vector{};
    scalar_t _padding_symbol{};

public:

    dp_vector_bulk() = default;
    explicit dp_vector_bulk(dp_vector_t dp_vector, simd_t padding_vector) :
        _dp_vector{std::move(dp_vector)},
        _padding_symbol{padding_vector[0]}
    {}

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

    dp_vector_t & base() noexcept
    {
        return _dp_vector;
    }

    dp_vector_t const & base() const noexcept
    {
        return _dp_vector;
    }

    decltype(auto) range() noexcept
    {
        return _dp_vector.range();
    }

    decltype(auto) range() const noexcept
    {
        return _dp_vector.range();
    }

    template <std::ranges::forward_range sequence_collection_t, typename initialisation_strategy_t>
    auto initialise(sequence_collection_t && sequence_collection, initialisation_strategy_t && init_strategy)
    {
        size_t max_sequence_size = 0;

        std::ranges::for_each(sequence_collection, [&] (auto && sequence) {
            max_sequence_size = std::max<size_t>(max_sequence_size, std::ranges::distance(sequence));
        });

        std::vector<simd_t, seqan3::aligned_allocator<simd_t, alignof(simd_t)>> simd_sequence{};
        simd_sequence.reserve(max_sequence_size);

        auto simd_view = sequence_collection | seqan3::views::to_simd<native_simd_t>(_padding_symbol);

        for (auto && simd_vector_chunk : simd_view) {
            for (auto && simd_vector : simd_vector_chunk) {
                simd_sequence.emplace_back(std::move(simd_vector));
            }
        }
        return _dp_vector.initialise(std::move(simd_sequence), std::forward<initialisation_strategy_t>(init_strategy));
    }
};

namespace detail
{

struct dp_vector_bulk_factory_fn
{
    template <typename dp_vector_t, typename simd_t>
    auto operator()(dp_vector_t && dp_vector, simd_t padding_vector) const noexcept
    {
        return dp_vector_bulk<std::remove_cvref_t<dp_vector_t>, simd_t>{
            std::forward<dp_vector_t>(dp_vector),
            std::move(padding_vector)
        };
    }
};

} // namespace detail

inline constexpr detail::dp_vector_bulk_factory_fn dp_vector_bulk_factory{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
