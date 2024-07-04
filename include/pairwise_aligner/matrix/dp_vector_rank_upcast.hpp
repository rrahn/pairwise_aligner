// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_rank_upcast.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <seqan3/utility/container/aligned_allocator.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
template <typename dp_vector_t, typename target_simd_t>
class dp_vector_rank_upcast
{
private:
    dp_vector_t _dp_vector{};

public:

    dp_vector_rank_upcast() = default;
    explicit dp_vector_rank_upcast(dp_vector_t dp_vector) : _dp_vector{std::move(dp_vector)}
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

    template <std::ranges::forward_range sequence_t, typename initialisation_strategy_t>
    auto initialise(sequence_t && sequence, initialisation_strategy_t && init_strategy)
    {
        std::vector<target_simd_t, seqan3::aligned_allocator<target_simd_t, alignof(target_simd_t)>> tmp_sequence{};
        tmp_sequence.resize(std::ranges::distance(sequence));
        std::ranges::copy(sequence | std::views::transform([&] (auto const & symbol) -> target_simd_t {
            return static_cast<target_simd_t>(symbol);
        }), tmp_sequence.begin());

        return _dp_vector.initialise(std::move(tmp_sequence), std::forward<initialisation_strategy_t>(init_strategy));
    }
};

namespace detail
{

template <typename simd_t>
struct dp_vector_rank_upcast_factory_fn
{
    template <typename dp_vector_t>
    auto operator()(dp_vector_t && dp_vector) const noexcept
    {
        return dp_vector_rank_upcast<std::remove_cvref_t<dp_vector_t>, simd_t>{
            std::forward<dp_vector_t>(dp_vector)
        };
    }
};

} // namespace detail

template <typename simd_t>
inline constexpr detail::dp_vector_rank_upcast_factory_fn dp_vector_rank_upcast_factory<simd_t>{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
