// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_column.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/matrix/dp_matrix_block.hpp>
#include <pairwise_aligner/matrix/dp_matrix_data_handle.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

template <typename block_closure_t, typename ...dp_data_t>
struct _column
{
    class type;
};

template <typename block_closure_t, typename ...dp_data_t>
using column_t = typename _column<block_closure_t, dp_data_t...>::type;

template <typename block_closure_t, typename ...dp_data_t>
class _column<block_closure_t, dp_data_t...>::type : public detail::dp_matrix_data_handle<dp_data_t...>
{
protected:

    using base_t = detail::dp_matrix_data_handle<dp_data_t...>;

    block_closure_t _block_closure{};

public:
    type() = delete;
    constexpr explicit type(block_closure_t block_closure, dp_data_t ...dp_data) noexcept :
        base_t{std::forward<dp_data_t>(dp_data)...},
        _block_closure{std::move(block_closure)}
    {
        rotate_row_scores_right(base_t::row());
    }

    ~type() noexcept
    {
        rotate_row_scores_left(base_t::row());
    }

    constexpr size_t size() const noexcept // override
    {
        return base_t::column().size();
    }

    // So different implementations can have different columns
    constexpr auto operator[](size_t const index) noexcept
        // -> decltype(dp_matrix_block(_dp_column[index], index))
    {
        assert(index < size());
        return make_matrix_block(base_t::column()[index],
                                 base_t::row(),
                                 base_t::substitution_model(),
                                 base_t::tracker(),
                                 base_t::row_sequence());
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

namespace cpo {
template <typename block_closure_t = dp_matrix::cpo::_block_closure<>>
struct _column_closure
{
    block_closure_t block_closure{};

    template <typename ...dp_data_t>
    constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
        using dp_column_t = dp_matrix::column_t<block_closure_t, dp_data_t...>;

        return dp_column_t{block_closure, std::forward<dp_data_t>(dp_data)...};
    }
};

} // namespace cpo
} // namespace dp_matrix

inline constexpr dp_matrix::cpo::_column_closure<> dp_matrix_column{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
