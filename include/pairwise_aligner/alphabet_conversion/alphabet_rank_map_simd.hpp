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

#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/utility/container/aligned_allocator.hpp>

#include <pairwise_aligner/simd/simd_rank_selector.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

// Now the rank type is a simd type.
template <typename key_t>
struct _alphabet_rank_map_simd
{
    class type;
};

template <typename key_t>
using  alphabet_rank_map_simd = typename _alphabet_rank_map_simd<key_t>::type;

template <typename key_t>
class _alphabet_rank_map_simd<key_t>::type : protected detail::simd_rank_selector_t<key_t>
{
private:
    using simd_rank_selector_t = detail::simd_rank_selector_t<key_t>;
    using rank_t = key_t;
    using scalar_key_t = typename rank_t::value_type;
    using rank_map_t = typename simd_rank_selector_t::rank_map_t;

    std::shared_ptr<rank_map_t> _rank_map_ptr{};
    rank_t _min_rank_offset{};

public:

    using key_type = key_t;
    using value_type = key_t;

    type() = default;

    template <std::ranges::random_access_range symbol_list_t>
        requires (std::convertible_to<std::ranges::range_reference_t<symbol_list_t const>, scalar_key_t>)
    explicit constexpr type(symbol_list_t const & valid_symbols)
    {
        if (valid_symbols.empty()) {
           throw std::invalid_argument{"The given set of valid symbols is empty!"};
        }

        // Assume for now the list of symbols is unique.
        // First scan list of symbols to find the minimal and maximal decimal rank value in list of symbols.
        auto [min_rank, max_rank] = std::ranges::minmax(valid_symbols,
                                                        std::ranges::less{},
                                                        [] (auto const & symbol) -> scalar_key_t {
            return static_cast<scalar_key_t>(symbol);
        });

        assert(min_rank >= 0);
        assert(max_rank >= min_rank);
        _min_rank_offset = rank_t{min_rank};
        // For example DNA5 = {A, C, G, N, T} as char only requires range of size (84(T) - 65(A)) + 1 = 20.
        scalar_key_t const max_rank_range = max_rank - min_rank + 1;

        assert(static_cast<int32_t>(std::numeric_limits<scalar_key_t>::max()) >
               static_cast<int32_t>(max_rank_range));

        // Resize the map to hold all ranks.
        size_t const rank_width = (rank_t::size + max_rank_range - 1) / rank_t::size;
        std::vector<rank_t, seqan3::aligned_allocator<rank_t, alignof(rank_t)>> ranks{};
        ranks.resize(rank_width, rank_t{static_cast<scalar_key_t>(-1)});

        // Fill rank map.
        for (scalar_key_t rank = 0; rank < std::ranges::ssize(valid_symbols); ++rank) {
            scalar_key_t key = static_cast<scalar_key_t>(valid_symbols[rank]);
            auto [index, position] = key_offset(key);
            ranks[index][position] = rank;
        }

        _rank_map_ptr = std::make_shared<rank_map_t>(simd_rank_selector_t::initialise_rank_map(std::move(ranks)));
    }

    // Now it looks like a map.
    constexpr value_type operator[](key_type const & key) const noexcept
    {
        return simd_rank_selector_t::select_rank_for(*_rank_map_ptr, normalise_key(key));
    }

private:
    constexpr key_type normalise_key(key_type const & key) const noexcept
    {
        return key - _min_rank_offset;
    }

    constexpr std::pair<scalar_key_t, scalar_key_t> key_offset(scalar_key_t key) const noexcept
    {
        scalar_key_t normalised_key = key - _min_rank_offset[0];
        return {normalised_key / rank_t::size, normalised_key % rank_t::size};
    }
};

} // inline namespace v1
} // namespace seqan::pairwise_aligner
