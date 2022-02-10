// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_block.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/matrix/dp_matrix_lane.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

template <typename dp_column_t, typename dp_row_t, typename substitution_model_t, typename tracker_t>
struct _block
{
    class type;
};

template <typename dp_column_t, typename dp_row_t, typename substitution_model_t, typename tracker_t>
using block_t = typename _block<dp_column_t, dp_row_t, substitution_model_t, tracker_t>::type;

template <typename dp_column_t, typename dp_row_t, typename substitution_model_t, typename tracker_t>
class _block<dp_column_t, dp_row_t, substitution_model_t, tracker_t>::type
{
    static constexpr size_t _lane_width = 8; //TODO: make copnfigurable?

    dp_column_t _dp_column;
    dp_row_t _dp_row;
    substitution_model_t _substitution_model;
    tracker_t _tracker;

public:

    using dp_column_type = std::remove_reference_t<dp_column_t>;
    using dp_row_type = std::remove_reference_t<dp_row_t>;

    type() = delete;
    constexpr explicit type(dp_column_t dp_column,
                            dp_row_t dp_row,
                            substitution_model_t substitution_model,
                            tracker_t tracker) noexcept :
        _dp_column{std::forward<dp_column_t>(dp_column)},
        _dp_row{std::forward<dp_row_t>(dp_row)},
        _substitution_model{std::forward<substitution_model_t>(substitution_model)},
        _tracker{std::forward<tracker_t>(tracker)}
    {
        // Note the first column/row is not computed again, as they were already initialised.
        _dp_column[0].score() = _dp_row[0].score();
        // std::cout << "First score = " << _dp_handle.dp_row[0].score() << "\n";
    }

    ~type() noexcept
    {
        // Store score of last column in first cell of row.
        _dp_row[0].score() = _dp_column[_dp_column.size() - 1].score();
        // std::cout << "Last score = " << _dp_handle.dp_row[0].score() << "\n";
    }

    constexpr auto operator[](size_t const index) noexcept
    {
        // return dp_matrix_block(column()[index], row());
        // what can we return now?
        return dp_matrix_lane(*this, index * lane_width());
    }

    constexpr auto last_lane() noexcept
    {
        // return dp_matrix_block(column()[index], row());
        // what can we return now?
        return dp_matrix_last_lane(*this, ((row().size() - 1) / lane_width()) * lane_width());
    }

    static constexpr size_t lane_width() noexcept
    {
        return _lane_width;
    }

    constexpr size_t size() const noexcept
    {
        return (_dp_row.size() - 1 + lane_width()) / lane_width();
    }

    constexpr dp_column_type & column() noexcept
    {
        return _dp_column;
    }

    constexpr dp_row_type & row() noexcept
    {
        return _dp_row;
    }

    constexpr std::remove_reference_t<substitution_model_t> & substitution_model() noexcept
    {
        return _substitution_model;
    }

    constexpr std::remove_reference_t<tracker_t> & tracker() noexcept
    {
        return _tracker;
    }

    // constexpr size_t size() const noexcept
    // {
    //     assert(_lane_size > 0);
    //     return (_lane_size - 1 + std::ranges::distance(_dp_handle.sequence2)) / _lane_size;
    // }

    // // So different implementations can have different columns
    // constexpr auto lane_at(size_t const lane_index) noexcept
    // {
    //     assert(lane_index < lane_count());

    //     size_t const first_position = lane_index * _lane_size;
    //     size_t const last_position = first_position + _lane_size;

    //     // slice the sequence position again!
    //     return _dp_lane_closure(
    //                 dp_matrix_handle{_dp_handle.sequence1,
    //                                  seqan3::views::slice(_dp_handle.sequence2, first_position, last_position),
    //                                  _dp_handle.dp_column, // put in index wrapper
    //                                  _dp_handle.dp_row}, // put in index wrapper
    //                                  first_position);
    // }
};

namespace cpo {

struct _block_closure
{
    template <typename dp_column_vector_t, typename dp_row_vector_t, typename substitution_model_t, typename tracker_t>
    constexpr auto operator()(dp_column_vector_t && dp_column,
                              dp_row_vector_t && dp_row,
                              substitution_model_t && substitution_model,
                              tracker_t && tracker) const noexcept {
        using dp_block_t = dp_matrix::block_t<dp_column_vector_t, dp_row_vector_t, substitution_model_t, tracker_t>;

        return dp_block_t{std::forward<dp_column_vector_t>(dp_column),
                          std::forward<dp_row_vector_t>(dp_row),
                          std::forward<substitution_model_t>(substitution_model),
                          std::forward<tracker_t>(tracker)};
    }

    // template <typename dp_matrix_column_t>
    // constexpr auto operator()(dp_matrix_column_t && dp_matrix_column, size_t const index) const noexcept {
    //     assert(index < dp_matrix_column.size());

    //     using dp_block_t = dp_matrix::block_t<dp_matrix_column_t>;

    //     return dp_block_t{std::forward<dp_matrix_column_t>(dp_matrix_column), index};
    // }
};

} // namespace cpo
} // namespace dp_matrxix

inline constexpr dp_matrix::cpo::_block_closure dp_matrix_block{};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
