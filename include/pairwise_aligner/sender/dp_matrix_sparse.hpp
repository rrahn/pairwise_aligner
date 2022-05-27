// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides .
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>
#include <ranges>

#include <pairwise_aligner/sender/dp_matrix_concept.hpp>
#include <pairwise_aligner/sender/dp_matrix_dimension.hpp>
#include <pairwise_aligner/sender/dp_matrix_position.hpp>
#include <pairwise_aligner/sender/dp_attribute.hpp>
#include <pairwise_aligner/sender/tag_invoke.hpp>

namespace seqan::align
{
inline namespace v1
{

// What if void?
template <typename first_t, typename second_t>
struct dp_matrix_sparse_value : public std::pair<first_t, second_t> {
    using base_t = std::pair<first_t, second_t>;
};

template <typename row_value_t, typename column_value_t>
struct dp_matrix_sparse_entry {
    row_value_t _row_value;
    column_value_t _column_value;

private:

    template <typename entry_t>
        requires std::same_as<std::remove_cvref_t<entry_t>, dp_matrix_sparse_entry>
    constexpr friend auto tag_invoke(tag_t<align::entry_diagonal>, entry_t && me) noexcept
        -> std::tuple_element_t<0, typename std::remove_cvref_t<entry_t>::row_value_t::base_t> {
        return std::get<0>(me._row_value);
    }

    template <typename entry_t>
        requires std::same_as<std::remove_cvref_t<entry_t>, dp_matrix_sparse_entry>
    constexpr friend auto tag_invoke(tag_t<align::entry_up>, entry_t && me) noexcept
        -> std::tuple_element_t<1, typename std::remove_cvref_t<entry_t>::row_value_t::base_t> {
        return std::get<1>(me._row_value);
    }

    template <typename entry_t>
        requires std::same_as<std::remove_cvref_t<entry_t>, dp_matrix_sparse_entry>
    constexpr friend auto tag_invoke(tag_t<align::entry_score>, entry_t && me) noexcept
        -> std::tuple_element_t<0, typename std::remove_cvref_t<entry_t>::column_value_t::base_t> {
        return std::get<0>(me._column_value);
    }

    template <typename entry_t>
        requires std::same_as<std::remove_cvref_t<entry_t>, dp_matrix_sparse_entry>
    constexpr friend auto tag_invoke(tag_t<align::entry_left>, entry_t && me) noexcept
        -> std::tuple_element_t<1, typename std::remove_cvref_t<entry_t>::column_value_t::base_t> {
        return std::get<1>(me._column_value);
    }
};

template <template <typename...> typename gap_cost_model_t>
class dp_matrix_sparse {

    // initialised matrix returned after initialisation
    template <typename dp_configuration_t>
    class activated_matrix {

        using score_t = typename dp_traits_t<std::remove_cvref_t<dp_configuration_t>>::score_type;
        using gap_model_t = gap_cost_model_t<score_t>;

        using diagonal_value_t = typename gap_model_t::diagonal_value_t;
        using up_value_t = typename gap_model_t::up_value_t;
        using left_value_t = typename gap_model_t::left_value_t;

        // TODO: we can check certain if the up and left value type are void or not!
        //

        // type members
        using row_value_t = dp_matrix_sparse_value_t<diagonal_value_t, up_value_t>;
        using column_value_t = dp_matrix_sparse_value_t<diagonal_value_t, left_value_t>;

        using row_t = std::vector<row_value_t>;
        using column_t = std::vector<column_value_t>;

        using row_reference_t = std::ranges::range_reference_t<row_t>;
        using column_reference_t = std::ranges::range_reference_t<column_t>;

        // maybe some configuration type that we can query for the gap model:
        dp_configuration_t _dp_configuration{};

        row_t _dp_row;
        column_t _dp_column;

    public:
        // How can we now adapt the dp_matrix entry type?
        // So this would be the base type
        // Now we want to wrap it through a series of wrapped dp_matrices.
        using value_type = dp_matrix_sparse_entry<row_value_t, column_value_t>;
        using reference = dp_matrix_sparse_entry<row_reference_t, column_reference_t>;

        activated_matrix(dp_configuration_t dp_configuration, dp_matrix_dimension const & dimension) :
            _dp_configuration{std::move(dp_configuration)}
        {
            // initialise the matrix
            _dp_row.resize(dimension.sequence_row_size + 1);
            _dp_column.resize(dimension.sequence_column_size + 1);

            size_t index = 0;
            std::ranges::for_each(_dp_row, [&] (row_reference_t row_value) {
                align::initialise_row_at(query(dp_attribute<gap_cost>, _dp_configuration), row_value, index);
                // align::initialise_row_at(cfg::query(cfg::attribute<cfg::gap_model>, _dp_config), row_value, index);
            });
            index = 0;
            std::ranges::for_each(_dp_column, [&] (row_reference_t column_value) {
                align::initialise_column_at(query(dp_attribute<gap_cost>, _dp_configuration), column_value, index);
                // align::initialise_column_at(cfg::query(cfg::attribute<cfg::gap_model>, _dp_config), column_value, index);
            });
        }

    private:

        constexpr friend auto tag_invoke(tag_t<align::entry_at>, activated_matrix & me,
                                         dp_matrix_position const & position) noexcept -> reference {
            return reference{me._dp_row[position.row_index], me._dp_column[position.column_index]};
        }

        constexpr friend auto tag_invoke(tag_t<align::result>, activated_matrix const & me) noexcept -> reference {
            return align::entry_score(align::entry_at(me,
                                                      dp_matrix_position{me._dp_row.size() - 1,
                                                                         me._dp_column.size() - 1}));
        }
    };

public:

    dp_matrix_sparse() // create from what types?
    {}

private:

    template <typename config_t>
    constexpr friend auto tag_invoke(tag_t<align::initialise>, dp_matrix_sparse & dp_matrix,
                                     config_t && config,
                                     dp_matrix_dimension const & dimension) noexcept
            -> activated_matrix<std::tuple_element_t<0, std::remove_cvref_t<config_t>>> {
        using score_t = std::tuple_element_t<0, std::remove_cvref_t<config_t>>;
        // on the way through the initialisation we can change configuration
        return activated_matrix<config_t>{std::forward<config_t>(cfg), dimension};

        // How can we now dynamically compose an algorithm?
        // dimension could also be part of matrix allocation attribute
    }
    // or forward the data to next sender!
    // implement initialise -> return initialised matrix
};

//

} // inline namespace v1
} // namespace seqan::align
