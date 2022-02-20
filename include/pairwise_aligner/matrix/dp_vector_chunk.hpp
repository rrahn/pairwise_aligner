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
#include <seqan3/std/span>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_vector_t>
class dp_vector_chunk
{
private:

    std::vector<dp_vector_t> _dp_vector_chunks{};
    size_t _chunk_size{};

public:

    using range_type = std::vector<dp_vector_t>;
    using value_type = std::ranges::range_value_t<range_type>;
    using reference = std::ranges::range_reference_t<range_type>;
    using const_reference = std::ranges::range_reference_t<range_type const>;

    explicit dp_vector_chunk(dp_vector_t && dp_vector, size_t const chunk_size) noexcept :
        _dp_vector_chunks{1, std::move(dp_vector)},
        _chunk_size{chunk_size}
    {}

    reference operator[](size_t const pos) noexcept
    {
        return _dp_vector_chunks[pos];
    }

    const_reference operator[](size_t const pos) const noexcept
    {
        return _dp_vector_chunks[pos];
    }

    constexpr size_t size() const noexcept
    {
        return _dp_vector_chunks.size();
    }

    constexpr size_t chunk_size() const noexcept
    {
        return _chunk_size;
    }

    range_type & range() noexcept
    {
        return _dp_vector_chunks;
    }

    range_type const & range() const noexcept
    {
        return _dp_vector_chunks;
    }

    // initialisation interface
    template <typename predecessor_t>
    struct _factory
    {
        predecessor_t _predecessor;
        size_t _offset;

        template <typename op_t>
        struct _op
        {
            op_t _op;
            size_t _offset;

            constexpr auto operator()(size_t const index) noexcept
            {
                return _op(index + _offset);
            }
        };

        template <typename score_t>
        constexpr auto create() const noexcept
        {
            using op_t = std::remove_reference_t<decltype(std::declval<predecessor_t>().template create<score_t>())>;
            return _op<op_t>{_predecessor.template create<score_t>(), _offset};
        }
    };

    template <std::ranges::forward_range sequence_t, typename factory_t>
    auto initialise(sequence_t && sequence, factory_t && init_factory)
    {
        using pure_factory_t = std::remove_cvref_t<factory_t>;

        size_t const sequence_size = std::ranges::distance(sequence);
        size_t const chunk_size = std::min(sequence_size, _chunk_size);
        size_t const element_count = (chunk_size > 0) ? (sequence_size + chunk_size - 1) / chunk_size : 1;

        _dp_vector_chunks.resize(element_count, _dp_vector_chunks.front());

        for (size_t i = 0; i < _dp_vector_chunks.size(); ++i)
        {
            size_t const first = i * chunk_size;
            size_t const last = (i + 1) * chunk_size;
            std::span tmp{std::ranges::next(std::ranges::begin(sequence), first),
                          std::ranges::next(std::ranges::begin(sequence), last, std::ranges::end(sequence))};

            _dp_vector_chunks[i].initialise(std::move(tmp), _factory<pure_factory_t>{init_factory, first});
        }

        return std::forward<sequence_t>(sequence);
    }
};

namespace detail
{

struct dp_vector_chunk_factory_fn
{

    template <typename dp_vector_t>
    auto operator()(dp_vector_t && dp_vector, size_t const block_size = -1) const noexcept
    {
        return dp_vector_chunk<std::remove_cvref_t<dp_vector_t>>{std::forward<dp_vector_t>(dp_vector), block_size};
    }
};

} // namespace detail

inline constexpr detail::dp_vector_chunk_factory_fn dp_vector_chunk_factory{};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
