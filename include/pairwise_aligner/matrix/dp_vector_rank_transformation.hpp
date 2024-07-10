// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_rank_transformation.
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
template <typename dp_vector_t, typename rank_map_t>
class dp_vector_rank_transformation
{
private:
    dp_vector_t _dp_vector{};
    rank_map_t _rank_map{};

    using rank_t = typename rank_map_t::value_type;

    static constexpr bool is_simd_rank_v = !std::integral<rank_t>;

public:

    dp_vector_rank_transformation() = default;
    explicit dp_vector_rank_transformation(dp_vector_t dp_vector, rank_map_t rank_map) :
        _dp_vector{std::move(dp_vector)},
        _rank_map{std::move(rank_map)}
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
        std::vector<rank_t, seqan3::aligned_allocator<rank_t, alignof(rank_t)>> rank_sequence{};
        rank_sequence.resize(std::ranges::distance(sequence));
        std::ranges::copy(sequence | std::views::transform([&] (auto const & symbol) -> rank_t {
            return _rank_map[symbol];
        }), rank_sequence.begin());

        return _dp_vector.initialise(std::move(rank_sequence), std::forward<initialisation_strategy_t>(init_strategy));
    }

    template <std::ranges::forward_range sequence_t, typename initialisation_strategy_t>
        requires (is_simd_rank_v && std::integral<std::ranges::range_value_t<sequence_t>>)
    auto initialise(sequence_t && sequence, initialisation_strategy_t && init_strategy)
    {
        using scalar_rank_t = typename rank_t::value_type;
        // Load the sequence into a single vector of simd values.
        std::ptrdiff_t sequence_size = std::ranges::distance(sequence);
        std::vector<scalar_rank_t> rank_sequence{};
        std::ptrdiff_t max_size = (sequence_size - 1 + rank_t::size_v) / rank_t::size_v;
        rank_sequence.reserve(max_size * rank_t::size_v);
        rank_sequence.resize(sequence_size);

        if constexpr (std::ranges::contiguous_range<sequence_t>) {
            std::ranges::for_each(std::views::iota(0, max_size), [&] (std::ptrdiff_t i) {
                std::ptrdiff_t memory_offset = i * rank_t::size_v;
                rank_t tmp{};
                tmp.load(reinterpret_cast<scalar_rank_t const *>(sequence.data()) + memory_offset);
                tmp = _rank_map[tmp]; // convert the rank.
                tmp.store(rank_sequence.data() + memory_offset);
            });
        } else {
            std::ranges::for_each(std::views::iota(0, max_size), [&] (std::ptrdiff_t i) {
                std::ptrdiff_t offset = i * rank_t::size_v;
                std::ptrdiff_t offset_end = std::min<std::ptrdiff_t>(rank_t::size_v, sequence_size);
                rank_t tmp{};
                // load from non-contiguous memory:
                for (std::ptrdiff_t k = 0; k + offset < offset_end; ++k) {
                    tmp[k] = sequence[k + offset];
                }
                tmp = _rank_map[tmp]; // convert the rank.
                tmp.store(rank_sequence.data() + offset);
            });
        }

        return _dp_vector.initialise(std::move(rank_sequence), std::forward<initialisation_strategy_t>(init_strategy));
    }
};

namespace detail
{

struct dp_vector_rank_transformation_factory_fn
{
    template <typename dp_vector_t, typename rank_map_t>
    auto operator()(dp_vector_t && dp_vector, rank_map_t rank_map) const noexcept
    {
        return dp_vector_rank_transformation<std::remove_cvref_t<dp_vector_t>, rank_map_t>{
            std::forward<dp_vector_t>(dp_vector),
            std::move(rank_map)
        };
    }
};

} // namespace detail

inline constexpr detail::dp_vector_rank_transformation_factory_fn dp_vector_rank_transformation_factory{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
