// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::aligner_result.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <limits>

#include <pairwise_aligner/configuration/end_gap_policy.hpp>
#include <pairwise_aligner/result/aligner_result.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace _chunk_factory
{

namespace detail
{
template <typename dp_vector_grouped_t>
class dp_vector_join
{
private:

    dp_vector_grouped_t _grouped_dp_vector{};
    size_t _chunk_size{};
    size_t _size{};

public:

    using range_type = typename dp_vector_grouped_t::range_type;
    using dp_vector_t = typename range_type::value_type;

    using value_type = typename dp_vector_t::value_type;
    using reference = typename dp_vector_t::reference;
    using const_reference = typename dp_vector_t::const_reference;

    explicit dp_vector_join(dp_vector_grouped_t && dp_vector) : _grouped_dp_vector{std::move(dp_vector)}
    {
        _chunk_size = _grouped_dp_vector[0].size() - 1;
        size_t const full_chunks = _grouped_dp_vector.size() - 1;
        _size = (full_chunks * _chunk_size) + _grouped_dp_vector[full_chunks].size();
        // std::cout << "_chunk_size = " << _chunk_size << "\n";
        // std::cout << "full_chunks = " << full_chunks << "\n";
        // std::cout << "_size = " << _size << "\n";
        // std::cout << "_grouped_dp_vector[full_chunks].size() = " << _grouped_dp_vector[full_chunks].size() << "\n";
    }

    reference operator[](size_t const pos) noexcept
    {
        auto [chunk_id, chunk_offset] = to_local_position(pos);
        // std::cout << "[chunk_id, chunk_offset] = [" << chunk_id << ", " << chunk_offset << "]\n";
        return _grouped_dp_vector[chunk_id][chunk_offset];
    }

    const_reference operator[](size_t const pos) const noexcept
    {
        auto [chunk_id, chunk_offset] = to_local_position(pos);
        // std::cout << "[chunk_id, chunk_offset] = [" << chunk_id << ", " << chunk_offset << "]\n";
        return _grouped_dp_vector[chunk_id][chunk_offset];
    }

    constexpr size_t size() const noexcept
    {
        // std::cout << "size = " << _size << "\n";
        return _size;
    }

    constexpr dp_vector_grouped_t & base() noexcept
    {
        return _grouped_dp_vector;
    }

    constexpr dp_vector_grouped_t const & base() const noexcept
    {
        return _grouped_dp_vector;
    }

private:

    constexpr std::pair<size_t, size_t> to_local_position(size_t const pos) const noexcept
    {
        // 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
        // 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
                                      //  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10
        size_t idx = pos / _chunk_size;
        size_t offset = pos % _chunk_size;
        size_t factor = (idx == _grouped_dp_vector.size());
        // std::cout << "factor = " << factor << "\n";
        return {idx - factor, offset + (_chunk_size * factor)};
    }

};
} // namespace detail

template <typename base_value_t>
struct _value
{
    struct type;
};

template <typename base_value_t>
using value = typename _value<base_value_t>::type;

template <typename base_value_t>
struct _value<base_value_t>::type : public base_value_t
{
    type() = delete;
    type(base_value_t && _value) : base_value_t{std::move(_value)}
    {}

    auto const & dp_column() const & noexcept
    {
        return this->dp_column().base();
    }

    auto const & dp_row() const & noexcept
    {
        return this->dp_row().base();
    }
};

} // namespace _chunk_factory

template <typename other_factory_t>
struct _result_factory_chunk
{
    struct type;
};

template <typename other_factory_t>
using result_factory_chunk = typename _result_factory_chunk<other_factory_t>::type;

template <typename other_factory_t>
struct _result_factory_chunk<other_factory_t>::type
{
    other_factory_t _factory;

    template <typename sequence1_t,
              typename sequence2_t,
              typename dp_column_t,
              typename dp_row_t>
    auto operator()(sequence1_t && sequence1,
                    sequence2_t && sequence2,
                    dp_column_t dp_column,
                    dp_row_t dp_row,
                    cfg::end_gap _column_trailing_gaps = cfg::end_gap::penalised,
                    cfg::end_gap _row_trailing_gaps = cfg::end_gap::penalised) const noexcept
    {

        auto result = _factory(std::forward<sequence1_t>(sequence1),
                               std::forward<sequence2_t>(sequence2),
                               _chunk_factory::detail::dp_vector_join<dp_column_t>{std::move(dp_column)},
                               _chunk_factory::detail::dp_vector_join<dp_row_t>{std::move(dp_row)},
                               _column_trailing_gaps,
                               _row_trailing_gaps);

        using result_t = decltype(result);
        return _chunk_factory::value<result_t>{std::move(result)};
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
