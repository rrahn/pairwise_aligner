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
#include <ranges>
#include <type_traits>

#include <pairwise_aligner/affine/affine_gap_model.hpp>
#include <pairwise_aligner/affine/affine_initialisation_strategy.hpp>
#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_attorney.hpp>
#include <pairwise_aligner/utility/math.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <template <typename> typename dp_template,
          typename ...policies_t>
class affine_dp_algorithm : public dp_template<affine_dp_algorithm<dp_template, policies_t...>>, // crtp-base
                            protected policies_t... // policies
{
private:

    using base_t = dp_template<affine_dp_algorithm<dp_template, policies_t...>>;

    friend class dp_algorithm_attorney<affine_dp_algorithm<dp_template, policies_t...>>;

public:
    affine_dp_algorithm() = default;

    explicit affine_dp_algorithm(policies_t const & ...policies) : base_t{}, policies_t{policies}...
    {}

protected:

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_row_vector(sequence_t && sequence, dp_vector_t & dp_vector) const
    {
        using gap_score_t = decltype(this->gap_open_score);
        using gap_model_t = affine_gap_model<gap_score_t>;
        using init_t = affine_initialisation_strategy<dp_vector_order::row, gap_model_t>;

        return dp_vector.initialise(std::forward<sequence_t>(sequence),
                                    init_t{gap_model_t{this->gap_open_score, this->gap_extension_score},
                                           this->first_row});
    }

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_column_vector(sequence_t && sequence, dp_vector_t & dp_vector) const
    {
        using gap_score_t = decltype(this->gap_open_score);
        using gap_model_t = affine_gap_model<gap_score_t>;
        using init_t = affine_initialisation_strategy<dp_vector_order::column, gap_model_t>;

        return dp_vector.initialise(std::forward<sequence_t>(sequence),
                                    init_t{gap_model_t{this->gap_open_score, this->gap_extension_score},
                                           this->first_column});
    }

    template <typename row_cell_t, typename column_cell_t>
    void compute_first_column(row_cell_t & first_row_cell, column_cell_t const & current_column_cell) const noexcept
    {
        get<0>(first_row_cell) = get<0>(current_column_cell);
    }

    template <typename row_cell_t, typename column_cell_t>
    auto initialise_column(row_cell_t & current_row_cell, column_cell_t & first_column_cell) const noexcept
    {
        using std::max;
        using score_t = typename row_cell_t::score_type;

        std::pair cache{get<0>(first_column_cell), get<1>(current_row_cell)};
        get<0>(first_column_cell) = get<0>(current_row_cell);
        get<1>(first_column_cell) = max(static_cast<score_t>(get<0>(current_row_cell) + (this->gap_open_score +
                                                                                         this->gap_extension_score)),
                                        static_cast<score_t>(get<1>(first_column_cell) +
                                                             this->gap_extension_score));
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
    constexpr auto initialise_substitution_scheme(args_t && ...args) const noexcept
    {
        return this->make_substitution_scheme(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr auto initialise_tracker(args_t && ...args) const noexcept
    {
        return this->make_tracker(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr auto initialise_policies(args_t && ...args) const noexcept
    {
        return this->make_policies(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr auto make_result(args_t && ...args) const noexcept
    {
        return this->operator()(std::forward<args_t>(args)..., this->last_column, this->last_row);
    }

    constexpr auto lane_width() const noexcept
    {
        return this->make_lane_width();
    }

    template <typename cache_t,
              typename dp_cell_t,
              typename scorer_t,
              typename tracker_t,
              typename seq1_val_t,
              typename seq2_val_t>
    constexpr auto compute_cell(cache_t & cache,
                                dp_cell_t & column_cell,
                                scorer_t & scorer,
                                tracker_t & tracker,
                                [[maybe_unused]] seq1_val_t const & seq1_val,
                                [[maybe_unused]] seq2_val_t const & seq2_val) const noexcept
    {
        using std::max;
        using score_t = typename dp_cell_t::score_type;

        score_t best = scorer.score(cache.first, seq1_val, seq2_val);
        best = max(max(best, cache.second), get<1>(column_cell));
        cache.first = get<0>(column_cell); // cache next diagonal score!
        get<0>(column_cell) = tracker.track(best);
        best = add(best, (this->gap_open_score + this->gap_extension_score));
        cache.second = max(static_cast<score_t>(add(cache.second, this->gap_extension_score)), best);
        get<1>(column_cell) = max(static_cast<score_t>(add(get<1>(column_cell), this->gap_extension_score)), best);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
