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

template <typename ...dp_data_t>
struct _column
{
    class type;
};

template <typename ...dp_data_t>
using column_t = typename _column<dp_data_t...>::type;

template <typename ...dp_data_t>
class _column<dp_data_t...>::type : public detail::dp_matrix_data_handle<dp_data_t...>
{
protected:

    using base_t = detail::dp_matrix_data_handle<dp_data_t...>;

public:
    type() = delete;
    constexpr explicit type(dp_data_t ...dp_data) noexcept : base_t{std::forward<dp_data_t>(dp_data)...}
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
        return dp_matrix_block(base_t::column()[index],
                               base_t::row(),
                               base_t::substitution_model(),
                               base_t::tracker(),
                               base_t::row_sequence());
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
struct _column_closure
{
    // template <typename dp_column_vector_t, typename dp_row_vector_t, typename substitution_model_t, typename tracker_t>
    // constexpr auto operator()(dp_column_vector_t && dp_column,
    //                           dp_row_vector_t && dp_row,
    //                           substitution_model_t && substitution_model,
    //                           tracker_t && tracker) const noexcept {
    //     using dp_column_t = dp_matrix::column_t<dp_column_vector_t, dp_row_vector_t, substitution_model_t, tracker_t>;

    //     return dp_column_t{std::forward<dp_column_vector_t>(dp_column),
    //                        std::forward<dp_row_vector_t>(dp_row),
    //                        std::forward<substitution_model_t>(substitution_model),
    //                        std::forward<tracker_t>(tracker)};
    // }

    template <typename ...dp_data_t>
    constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
        using dp_column_t = dp_matrix::column_t<dp_data_t...>;

        return dp_column_t{std::forward<dp_data_t>(dp_data)...};
    }
};

} // namespace cpo
} // namespace dp_matrix

inline constexpr dp_matrix::cpo::_column_closure dp_matrix_column{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
