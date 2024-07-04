// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_offset_transformation.
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
template <typename dp_vector_t, typename offset_fn_t>
class dp_vector_offset_transformation
{
private:
    dp_vector_t _dp_vector{};
    offset_fn_t _offset_fn{};

public:

    dp_vector_offset_transformation() = default;
    explicit dp_vector_offset_transformation(dp_vector_t dp_vector, offset_fn_t offset_fn) :
        _dp_vector{std::move(dp_vector)},
        _offset_fn{std::move(offset_fn)}
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
        // expect simd range!
        using offset_t = std::invoke_result_t<offset_fn_t, std::ranges::range_rvalue_reference_t<sequence_t>>;

        std::vector<offset_t, seqan3::aligned_allocator<offset_t, alignof(offset_t)>> offset_sequence{};
        offset_sequence.resize(std::ranges::distance(sequence));
        std::ranges::move(sequence | std::views::transform([&] (auto && symbol) -> offset_t {
            return _offset_fn(std::move(symbol));
        }), offset_sequence.begin());

        return _dp_vector.initialise(std::move(offset_sequence),
                                     std::forward<initialisation_strategy_t>(init_strategy));
    }
};

namespace detail
{

struct dp_vector_offset_transformation_factory_fn
{
    template <typename dp_vector_t, typename offset_fn_t>
    auto operator()(dp_vector_t && dp_vector, offset_fn_t && offset_fn) const noexcept
    {
        return dp_vector_offset_transformation<std::remove_cvref_t<dp_vector_t>, offset_fn_t>{
            std::forward<dp_vector_t>(dp_vector),
            std::forward<offset_fn_t>(offset_fn)
        };
    }
};

} // namespace detail

inline constexpr detail::dp_vector_offset_transformation_factory_fn dp_vector_offset_transformation_factory{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
