// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_column_base.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cassert>

#include <pairwise_aligner/matrix/dp_matrix_state_handle.hpp>
namespace seqan::pairwise_aligner {
inline namespace v1 {
namespace dp_matrix {
namespace detail {

template <typename block_closure_t, typename data_handle_t>
struct _column_base
{
    class type;
};

template <typename block_closure_t, typename data_handle_t>
using column_base_t = typename _column_base<block_closure_t, data_handle_t>::type;

template <typename block_closure_t, typename data_handle_t>
class _column_base<block_closure_t, data_handle_t>::type : public data_handle_t
{
protected:

    using base_t = data_handle_t;

    block_closure_t _block_closure{};

public:
    type() = delete;
    template <typename ...dp_data_t>
    constexpr explicit type(block_closure_t block_closure, dp_data_t && ...dp_data) noexcept :
        base_t{std::forward<dp_data_t>(dp_data)...},
        _block_closure{std::move(block_closure)}
    {
        rotate_row_scores_right(base_t::row());
    }

    ~type() noexcept
    {
        rotate_row_scores_left(base_t::row());
    }

    // So different implementations can have different columns
    constexpr auto operator[](size_t const index) noexcept
    {
        assert(index < size());
        return make_matrix_block(base_t::column()[index],
                                 base_t::row(),
                                 base_t::substitution_model(),
                                 base_t::tracker(),
                                 base_t::row_sequence(),
                                 base_t::lane_width());
    }

    constexpr size_t size() const noexcept // override
    {
        return base_t::column().size();
    }

protected:
    template <typename ...args_t>
    constexpr auto make_matrix_block(args_t && ...args) const noexcept
        -> std::invoke_result_t<block_closure_t, args_t...>
    {
        return std::invoke(_block_closure, std::forward<args_t>(args)...);
    }

private:
    constexpr void rotate_row_scores_right(typename base_t::row_type & dp_row) const noexcept
    {
        size_t const dp_row_size = dp_row.size() - 1;
        // cache score of last cell.
        auto tmp = std::move(dp_row[dp_row_size].score());

        // rotate scores right.
        for (size_t j = dp_row_size; j > 0; --j)
            dp_row[j].score() = dp_row[j - 1].score();

        // store last value in first cell.
        dp_row[0].score() = std::move(tmp);
    }

    constexpr void rotate_row_scores_left(typename base_t::row_type & dp_row) const noexcept
    {
        size_t const dp_row_size = dp_row.size() - 1;
        // cache score of first cell.
        auto tmp = std::move(dp_row[0].score());

        // rotate scores left.
        for (size_t j = 0; j < dp_row_size; ++j)
            dp_row[j].score() = dp_row[j + 1].score();

        // store cached score in last cell.
        dp_row[dp_row_size].score() = std::move(tmp);
    }
};

template <typename block_fn_t, typename ...dp_state_t>
class column_base : public state_handle<dp_state_t...>
{
protected:

    using base_t = state_handle<dp_state_t...>;

    block_fn_t _block_fn{};

public:
    column_base() = delete;

    constexpr explicit column_base(block_fn_t block_fn, dp_state_t ...dp_state) noexcept :
        base_t{std::forward<dp_state_t>(dp_state)...},
        _block_fn{std::move(block_fn)}
    {
        rotate_row_scores_right(base_t::dp_row());
    }

    ~column_base() noexcept
    {
        rotate_row_scores_left(base_t::dp_row());
    }

protected:

    constexpr auto column_slice_at(std::ptrdiff_t const index) const noexcept
    {
        std::ptrdiff_t const column_chunk_size = base_t::dp_column()[0].size() - 1;
        std::ptrdiff_t const column_offset = column_chunk_size * index;
        return seqan3::views::slice(base_t::column_sequence(), column_offset, column_offset + column_chunk_size);

    }

    template <typename ...args_t>
    constexpr auto make_matrix_block(args_t && ...args) const noexcept -> std::invoke_result_t<block_fn_t, args_t...>
    {
        return std::invoke(_block_fn, std::forward<args_t>(args)...);
    }

private:
    constexpr void rotate_row_scores_right(typename base_t::dp_row_type & dp_row) const noexcept
    {
        size_t const dp_row_size = dp_row.size() - 1;
        // cache score of last cell.
        auto tmp = std::move(dp_row[dp_row_size].score());

        // rotate scores right.
        for (size_t j = dp_row_size; j > 0; --j)
            dp_row[j].score() = dp_row[j - 1].score();

        // store last value in first cell.
        dp_row[0].score() = std::move(tmp);
    }

    constexpr void rotate_row_scores_left(typename base_t::dp_row_type & dp_row) const noexcept
    {
        size_t const dp_row_size = dp_row.size() - 1;
        // cache score of first cell.
        auto tmp = std::move(dp_row[0].score());

        // rotate scores left.
        for (size_t j = 0; j < dp_row_size; ++j)
            dp_row[j].score() = dp_row[j + 1].score();

        // store cached score in last cell.
        dp_row[dp_row_size].score() = std::move(tmp);
    }
};

} // namespace detail
} // namespace dp_matrix
} // inline namespace v1
} // namespace seqan::pairwise_aligner
