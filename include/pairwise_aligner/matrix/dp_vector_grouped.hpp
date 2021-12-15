// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_grouped.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/to.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_vector_t>
class dp_vector_grouped
{
private:

    std::vector<dp_vector_t> _dp_vector_group{};
    size_t _element_size{25};

public:

    using range_type = std::vector<dp_vector_t>;
    using value_type = std::ranges::range_value_t<range_type>;
    using reference = std::ranges::range_reference_t<range_type>;
    using const_reference = std::ranges::range_reference_t<range_type const>;

    dp_vector_grouped() = default;
    explicit dp_vector_grouped(size_t const element_size) noexcept : _element_size{element_size}
    {}

    reference operator[](size_t const pos) noexcept
    {
        return _dp_vector_group[pos];
    }

    const_reference operator[](size_t const pos) const noexcept
    {
        return _dp_vector_group[pos];
    }

    constexpr size_t size() const noexcept
    {
        return _dp_vector_group.size();
    }

    range_type & range() noexcept
    {
        return _dp_vector_group;
    }

    range_type const & range() const noexcept
    {
        return _dp_vector_group;
    }

    // initialisation interface
    template <typename predecessor_t>
    struct _factory
    {
        predecessor_t _predecessor;
        size_t _index;

        template <typename op_t>
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
            return _op<op_t>{_predecessor.template create<score_t>(), _index};
        }
    };

    template <std::ranges::forward_range sequence_t, typename factory_t>
    auto initialise(sequence_t && sequence, factory_t && init_factory)
    {
        using pure_factory_t = std::remove_cvref_t<factory_t>;
        // using value_t = std::ranges::range_value_t<sequence_t>;
        // using chunk_seq_t = std::vector<value_t>;
        size_t const sequence_size = std::ranges::distance(sequence);
        size_t const element_count = (sequence_size + _element_size - 1) / _element_size;

        _dp_vector_group.resize(element_count);

        // std::cout << "sequence_size = " << sequence_size << "\n";
        // std::cout << "_element_size = " << _element_size << "\n";
        // std::cout << "element_count = " << element_count << "\n";

        // std::vector<chunk_seq_t> sequence_group{};
        // sequence_group.resize(element_count);

        for (size_t i = 0; i < _dp_vector_group.size(); ++i)
        {
            size_t const first = i * _element_size;
            size_t const last = (i + 1) * _element_size;
            std::span tmp{std::ranges::next(std::ranges::begin(sequence), first),
                          std::ranges::next(std::ranges::begin(sequence), last, std::ranges::end(sequence))};

            // std::cout << "begin = " << begin << "\n";
            // std::cout << "end = " << end << "\n";
            _dp_vector_group[i].initialise(std::move(tmp), _factory<pure_factory_t>{init_factory, first});
            // std::cout << "size(sequence_group[i]) = " << std::ranges::size(sequence_group[i]) << "\n";
        }

        return std::forward<sequence_t>(sequence);
    }
};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
