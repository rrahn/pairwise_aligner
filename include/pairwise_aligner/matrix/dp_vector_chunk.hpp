// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_chunk.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_vector_t>
class dp_vector_chunk
{
private:

    dp_vector_t _dp_vector{};
    size_t _index{};

public:

    using range_type = typename dp_vector_t::range_type;
    using value_type = std::ranges::range_value_t<range_type>;
    using reference = std::ranges::range_reference_t<range_type>;
    using const_reference = std::ranges::range_reference_t<range_type const>;

    dp_vector_chunk() = default;
    explicit dp_vector_chunk(size_t const initial_index) noexcept : _index{initial_index}
    {}

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

    auto & range() noexcept
    {
        return _dp_vector.range();
    }

    auto const & range() const noexcept
    {
        return _dp_vector.range();
    }

    // initialisation interface
    template <typename predecessor_t>
    struct _factory
    {
        predecessor_t _predecessor;
        size_t _index;

        template <typename score_t, typename op_t>
        struct _op
        {
            op_t _op;
            size_t _index;

            constexpr auto operator()() noexcept
            {
                return _op(_index++);
            }
        };

        template <typename score_t>
        constexpr auto create() const noexcept
        {
            using op_t = std::remove_reference_t<decltype(std::declval<predecessor_t>().template create<score_t>())>;
            return _op<score_t, op_t>{_predecessor.template create<score_t>(), _index};
        }
    };

    template <typename sequence_t, typename initialisation_strategy_t>
    auto initialise(sequence_t && sequence, initialisation_strategy_t && init_factory)
    {
        using factory_t = _factory<initialisation_strategy_t>;
        return _dp_vector.initialise(std::forward<sequence_t>(sequence),
                                     factory_t{std::forward<initialisation_strategy_t>(init_factory), _index});
    }
};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
