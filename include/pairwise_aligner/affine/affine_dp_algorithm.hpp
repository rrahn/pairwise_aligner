// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::affine_dp_algorithm.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <functional>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/affine/affine_initialisation_strategy.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <template <typename> typename dp_template,
          typename score_model_t,
          typename result_factory_t,
          typename gap_model_t>
class affine_dp_algorithm : protected dp_template<affine_dp_algorithm<dp_template,
                                                                      score_model_t,
                                                                      result_factory_t,
                                                                      gap_model_t>>
{
private:

    using base_t = dp_template<affine_dp_algorithm<dp_template, score_model_t, result_factory_t, gap_model_t>>;

    friend base_t;

    score_model_t _score_model{};
    result_factory_t _result_factory{};
    gap_model_t _gap_model{};
    affine_initialisation_strategy<gap_model_t> _init{};

public:
    affine_dp_algorithm() = default;

    explicit affine_dp_algorithm(score_model_t score_model,
                                 result_factory_t result_factory,
                                 gap_model_t gap_model,
                                 initialisation_rule init_param) noexcept :
        _score_model{std::move(score_model)},
        _result_factory{std::move(result_factory)},
        _gap_model{std::move(gap_model)},
        _init{_gap_model, init_param}
    {}

protected:

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_row_vector(sequence_t && sequence, dp_vector_t & dp_vector) const
    {
        return dp_vector.initialise(std::forward<sequence_t>(sequence), _init);
    }

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_column_vector(sequence_t && sequence, dp_vector_t & dp_vector) const
    {
        return dp_vector.initialise(std::forward<sequence_t>(sequence), _init);
    }

    template <typename row_cell_t, typename column_cell_t>
    auto initialise_column(row_cell_t & current_row_cell, column_cell_t & first_column_cell) const noexcept
    {
        using std::max;
        using score_t = typename row_cell_t::score_type;

        std::pair cache{get<0>(first_column_cell), get<1>(current_row_cell)};
        get<0>(first_column_cell) = get<0>(current_row_cell);
        get<1>(first_column_cell) = max(static_cast<score_t>(cache.first + _gap_model.gap_open_score),
                                        static_cast<score_t>(get<1>(first_column_cell) +
                                                             _gap_model.gap_extension_score));
        return cache;
    }

    template <typename row_cell_t, typename column_cell_t, typename cache_t>
    void finalise_column(row_cell_t & current_row_cell,
                         column_cell_t const & last_column_cell,
                         cache_t & cache) const noexcept
    {
        get<0>(current_row_cell) = get<0>(last_column_cell);
        get<1>(current_row_cell) = std::move(cache.second);
    }

    template <typename ...args_t>
    constexpr auto make_result(args_t && ...args) const noexcept
    {
        return _result_factory(std::forward<args_t>(args)...);
    }

    template <typename cache_t, typename seq1_val_t, typename seq2_val_t, typename dp_cell_t>
    auto compute_cell(cache_t & cache,
                      dp_cell_t & column_cell,
                      seq1_val_t const & seq1_val,
                      seq2_val_t const & seq2_val) const noexcept
    {
        using std::max;
        using score_t = typename dp_cell_t::score_type;

        auto [next_diagonal, horizontal_score] = column_cell;
        cache.first += _score_model.score(seq1_val, seq2_val);
        cache.first = max(max(cache.first, cache.second), horizontal_score);
        get<0>(column_cell) = cache.first;
        cache.first += (_gap_model.gap_open_score + _gap_model.gap_extension_score);
        cache.second = max(static_cast<score_t>(cache.second + _gap_model.gap_extension_score), cache.first);
        get<1>(column_cell) = max(static_cast<score_t>(horizontal_score + _gap_model.gap_extension_score), cache.first);
        cache.first = next_diagonal; // cache score
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
