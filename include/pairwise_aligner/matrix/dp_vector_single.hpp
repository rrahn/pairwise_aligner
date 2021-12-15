// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_single.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

#include <vector>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_cell_t, typename dp_vector_t = std::vector<dp_cell_t>>
class dp_vector_single
{
private:
    dp_vector_t _dp_vector;

public:

    using range_type = dp_vector_t;
    using value_type = std::ranges::range_value_t<dp_vector_t>;
    using reference = std::ranges::range_reference_t<dp_vector_t>;
    using const_reference = std::ranges::range_reference_t<dp_vector_t const>;

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

    dp_vector_t & range() noexcept
    {
        return _dp_vector;
    }

    dp_vector_t const & range() const noexcept
    {
        return _dp_vector;
    }

    // range interface

    // initialisation interface

    template <std::ranges::forward_range sequence_t, typename initialisation_strategy_t>
    sequence_t initialise(sequence_t && sequence, initialisation_strategy_t && init_factory)
    {
        size_t const sequence_size = std::ranges::distance(sequence);
        _dp_vector.resize(sequence_size + 1);

        using score_t = typename dp_cell_t::score_type;

        std::ranges::generate(_dp_vector, init_factory.template create<score_t>());

        return sequence;
    }
};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
