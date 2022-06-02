// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides base dp matrix.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <pairwise_aligner/align_matrix/dp_entry_concept.hpp>
#include <pairwise_aligner/align_matrix/dp_matrix_concept.hpp>

namespace align {

namespace _linear_row_column {

template <typename row_value_t, typename column_value_t>
class _matrix {
private:

    using row_t = std::vector<row_value_t>;
    using column_t = std::vector<column_value_t>;

    struct entry;

    row_t _row{};
    column_t _column{};

public:

    using row_value_type = row_value_t;
    using column_value_type = column_value_t;
    using row_type = row_t;
    using column_type = column_t;
    using entry_type = entry;

private:

    template <typename matrix_position_t>
    constexpr friend entry tag_invoke(tag_t<align::entry_at>, _matrix & me,
                                      matrix_position_t const row_position,
                                      matrix_position_t const column_position) noexcept {
        assert(row_position < static_cast<matrix_position_t>(me._row.size()));
        assert(column_position < static_cast<matrix_position_t>(me._column.size()));

        return entry{me._row[row_position], me._column[column_position]};
    }

    template <typename gap_cost_model_t>
    constexpr friend void tag_invoke(tag_t<align::initialise_row>, _matrix & me,
                                     size_t const dimension,
                                     gap_cost_model_t const & gap_cost_model) noexcept {
        me._row.resize(dimension);
        size_t i = 0;
        std::ranges::for_each(me._row, [&] (auto & row_value) {
            align::diagonal_score(row_value) = align::score(gap_cost_model, i++);
            align::up_score(row_value) = align::score(gap_cost_model, i);
        });
    }

    template <typename gap_cost_model_t>
    constexpr friend void tag_invoke(tag_t<align::initialise_column>, _matrix & me,
                                     size_t const dimension,
                                     gap_cost_model_t const & gap_cost_model) noexcept {
        me._column.resize(dimension);
        size_t i = 0;
        std::ranges::for_each(me._column, [&] (auto & column_value) {
            align::current_score(column_value) = align::score(gap_cost_model, i++);
            align::left_score(column_value) = align::current_score(column_value);
        });
    }

};

template <typename row_value_t, typename column_value_t>
struct _matrix<row_value_t, column_value_t>::entry {
private:
    friend _matrix;

    using row_reference_t = std::ranges::range_reference_t<row_t>;
    using column_reference_t = std::ranges::range_reference_t<column_t>;

    row_reference_t _row_reference;
    column_reference_t _column_reference;

public:

    entry() = delete;
    entry(row_reference_t row_reference, column_reference_t column_reference) :
        _row_reference{row_reference},
        _column_reference{column_reference}
    {}

private:

    template <typename entry_t>
        requires std::same_as<std::remove_cvref_t<entry_t>, entry>
    constexpr friend auto tag_invoke(tag_t<align::diagonal_score>, entry_t && me)
        noexcept(noexcept(align::diagonal_score(me._row_reference)))
        -> std::invoke_result_t<tag_t<align::diagonal_score>, row_reference_t> {
        return align::diagonal_score(me._row_reference);
    }

    template <typename entry_t>
        requires std::same_as<std::remove_cvref_t<entry_t>, entry>
    constexpr friend auto tag_invoke(tag_t<align::up_score>, entry_t && me)
        noexcept(noexcept(align::up_score(me._row_reference)))
        -> std::invoke_result_t<tag_t<align::up_score>, row_reference_t> {
        return align::up_score(me._row_reference);
    }

    template <typename entry_t>
        requires std::same_as<std::remove_cvref_t<entry_t>, entry>
    constexpr friend auto tag_invoke(tag_t<align::current_score>, entry_t && me)
        noexcept(noexcept(align::current_score(me._column_reference)))
        -> std::invoke_result_t<tag_t<align::current_score>, column_reference_t> {
        return align::current_score(me._column_reference);
    }

    template <typename entry_t>
        requires std::same_as<std::remove_cvref_t<entry_t>, entry>
    constexpr friend auto tag_invoke(tag_t<align::left_score>, entry_t && me)
        noexcept(noexcept(align::left_score(me._column_reference)))
        -> std::invoke_result_t<tag_t<align::left_score>, column_reference_t> {
        return align::left_score(me._column_reference);
    }
};

struct _factory {
    template <typename row_value_t, typename column_value_t>
    using invoke = _matrix<row_value_t, column_value_t>;
};

namespace _cpo {
struct _fn {
    constexpr _factory operator()() const noexcept {
        return _factory{};
    }
};
} // namespace _cpo
} // namespace _linear_row_column

inline constexpr _linear_row_column::_cpo::_fn linear_row_column_matrix{};

} // namespace align
