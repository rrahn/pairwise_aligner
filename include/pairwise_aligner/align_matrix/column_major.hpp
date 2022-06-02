// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides just alignment prompt.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <algorithm>
#include <ranges>
#include <type_traits>

// #include <pairwise_aligner/align_matrix/linear_score_matrix.hpp>

namespace align
{
inline namespace v1
{

template <typename base_dp_matrix_t>
struct linear_score_matrix : protected base_dp_matrix_t {
private:

    using typename base_dp_matrix_t::score_t;

    using vector_entry_t = std::pair<score_t, score_t>;

    using row_t = std::vector<vector_entry_t>;
    using column_t = std::vector<vector_entry_t>;

    using row_reference_t = std::ranges::range_reference_t<row_t>;
    using column_reference_t = std::ranges::range_reference_t<column_t>;

    struct entry {
        row_reference_t _row_entry;
        column_reference_t _column_entry;

private:
        constexpr friend auto tag_invoke(tag_t<align::diagonal_score>, entry & me) noexcept -> score_t & {
            return me._row_entry.first;
        }

        constexpr friend auto tag_invoke(tag_t<align::up_score>, entry & me) noexcept -> score_t & {
            return me._row_entry.second;
        }

        constexpr friend auto tag_invoke(tag_t<align::optimal_score>, entry & me) noexcept -> score_t & {
            return me._column_entry.first;
        }

        constexpr friend auto tag_invoke(tag_t<align::left_score>, entry & me) noexcept -> score_t & {
            return me._column_entry.second;
        }
    };

    row_t _row;
    column_t _column;

public:

    linear_score_matrix(base_dp_matrix_t base_dp_matrix) : base_dp_matrix_t{std::move(base_dp_matrix)}
    {}

    constexpr size_t row_dimension() const noexcept {
        return _row.size();
    }

    constexpr size_t column_dimension() const noexcept {
        return _column.size();
    }

private:

    constexpr friend auto tag_invoke(tag_t<align::entry_at>,
                                     linear_score_matrix & me,
                                     size_t const row_position,
                                     size_t const column_position) noexcept {
        return base_dp_matrix_t::entry(entry{me._row[row_position], me._column[column_position]});
    }

    template <typename row_init_fn_t, typename column_init_fn_t>
    constexpr friend void tag_invoke(tag_t<align::initialise_row>,
                                     linear_score_matrix & me,
                                     size_t const row_dimension,
                                     size_t start_index) noexcept {
        me._row.resize(row_dimension);
        std::for_each(me._row, [&] (vector_entry_t & row_entry) {
            row_entry.first = score(me.row_initialiser, start_index);
            row_entry.second = score(me.row_initialiser, ++start_index);
        });
    }

    template <typename row_init_fn_t, typename column_init_fn_t>
    constexpr friend void tag_invoke(tag_t<align::initialise_column>,
                                     linear_score_matrix & me,
                                     size_t const column_dimension,
                                     size_t start_index) noexcept {
        me._column.resize(column_dimension);
        std::for_each(me._column, [&] (vector_entry_t & column_entry) {
            column_entry.second = column_entry.first = score(me.column_initialiser, start_index++);
        });
    }
};

} // inline namespace v1
} // namespace align
