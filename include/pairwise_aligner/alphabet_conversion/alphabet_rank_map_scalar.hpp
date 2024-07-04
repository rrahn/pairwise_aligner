// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::alphabet_rank_map_scalar.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <ranges>
#include <vector>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

class alphabet_rank_map_scalar
{
public:
    using rank_type = int32_t;
    using value_type = rank_type;

private:
    std::vector<rank_type> _map{};
    rank_type _scale{1};

public:
    alphabet_rank_map_scalar() = default;

    template <std::ranges::forward_range symbol_list_t>
        requires (std::convertible_to<std::ranges::range_reference_t<symbol_list_t const>, rank_type>)
    explicit constexpr alphabet_rank_map_scalar(symbol_list_t const & valid_symbols)
    {
        _map.resize(256, -1); // fill with -1.
        rank_type rank{};
        std::ranges::for_each(valid_symbols, [&] (auto const & symbol) {
            assert(static_cast<rank_type>(symbol) >= 0);
            assert(static_cast<rank_type>(symbol) < std::ranges::ssize(_map));

            _map[static_cast<rank_type>(symbol)] = rank++;
        });
    }

    constexpr void set_scale(int32_t scale) noexcept
    {
        _scale = scale;
    }

    template <std::convertible_to<rank_type> alphabet_t>
    constexpr rank_type operator[](alphabet_t && symbol) const noexcept
    {
        assert(_map[static_cast<rank_type>(symbol)] > -1);
        return _map[static_cast<rank_type>(symbol)] * _scale;
    }
};

} // inline namespace v1
} // namespace seqan::pairwise_aligner
