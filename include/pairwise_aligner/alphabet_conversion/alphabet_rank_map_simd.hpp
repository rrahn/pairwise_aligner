// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::alphabet_rank_map_simd.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <ranges>
#include <vector>

#include <seqan3/utility/container/aligned_allocator.hpp>

#include <pairwise_aligner/simd/simd_rank_selector.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

// Now the rank type is a simd type.
template <typename simd_from_t, typename simd_to_t>
struct _alphabet_rank_map_simd
{
    class type;
};

template <typename simd_from_t, typename simd_to_t = simd_from_t>
using  alphabet_rank_map_simd = typename _alphabet_rank_map_simd<simd_from_t, simd_to_t>::type;

template <typename simd_from_t, typename simd_to_t>
class _alphabet_rank_map_simd<simd_from_t, simd_to_t>::type : protected detail::simd_rank_selector_t<simd_to_t>
{
private:
    using simd_rank_selector_t = detail::simd_rank_selector_t<simd_to_t>; // conversion to value_type?
    using scalar_from_t = typename simd_from_t::value_type;
    using rank_map_t = typename simd_rank_selector_t::rank_map_t;

    std::shared_ptr<rank_map_t> _rank_map_ptr{};
    simd_to_t _min_rank_offset{};

public:

    using key_type = simd_from_t;
    using value_type = simd_to_t;

    type() = default;

    template <std::ranges::random_access_range symbol_list_t>
        requires (std::convertible_to<std::ranges::range_reference_t<symbol_list_t const>, scalar_from_t>)
    explicit constexpr type(symbol_list_t const & valid_symbols)
    {
        if (valid_symbols.empty()) {
           throw std::invalid_argument{"The given set of valid symbols is empty!"};
        }

        // Assume for now the list of symbols is unique.
        // First scan list of symbols to find the minimal and maximal decimal rank value in list of symbols.
        auto [min_rank, max_rank] = std::ranges::minmax(valid_symbols,
                                                        std::ranges::less{},
                                                        [] (auto const & symbol) -> scalar_from_t {
            return static_cast<scalar_from_t>(symbol);
        });

        assert(min_rank >= 0);
        assert(max_rank >= min_rank);
        _min_rank_offset = simd_from_t{min_rank};
        // For example DNA5 = {A, C, G, N, T} as char only requires range of size (84(T) - 65(A)) + 1 = 20.
        scalar_from_t const max_rank_range = max_rank - min_rank + 1;

        assert(static_cast<int32_t>(std::numeric_limits<scalar_from_t>::max()) >
               static_cast<int32_t>(max_rank_range));

        // Resize the map to hold all ranks.
        size_t const rank_width = (simd_from_t::size_v + max_rank_range - 1) / simd_from_t::size_v;
        std::vector<simd_from_t, seqan3::aligned_allocator<simd_from_t, alignof(simd_from_t)>> ranks{};
        ranks.resize(rank_width, simd_from_t{static_cast<scalar_from_t>(-1)});

        // Fill rank map.
        for (scalar_from_t rank = 0; rank < std::ranges::ssize(valid_symbols); ++rank) {
            scalar_from_t key = static_cast<scalar_from_t>(valid_symbols[rank]);
            auto [index, position] = key_offset(key);
            ranks[index][position] = rank;
        }

        _rank_map_ptr = std::make_shared<rank_map_t>(simd_rank_selector_t::initialise_rank_map(std::move(ranks)));
    }

    // Now it looks like a map.
    constexpr value_type operator[](key_type const & keys) const noexcept
    {
        return simd_rank_selector_t::select_rank_for(*_rank_map_ptr, normalise_keys(keys));
    }

private:
    constexpr key_type normalise_keys(key_type const & keys) const noexcept
    {
        return keys - _min_rank_offset;
    }

    constexpr std::pair<scalar_from_t, scalar_from_t> key_offset(scalar_from_t key) const noexcept
    {
        scalar_from_t normalised_key = key - _min_rank_offset[0];
        return {normalised_key / simd_from_t::size_v, normalised_key % simd_from_t::size_v};
    }
};

} // inline namespace v1
} // namespace seqan::pairwise_aligner
