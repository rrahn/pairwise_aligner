// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_algorithm_template_base.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_attorney.hpp>
#include <pairwise_aligner/result/aligner_result.hpp>
#include <pairwise_aligner/matrix/dp_matrix_cpo.hpp>
#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename algorithm_impl_t>
struct _dp_algorithm_template_base
{
    class type;
};

template <typename algorithm_impl_t>
using dp_algorithm_template_base = typename _dp_algorithm_template_base<algorithm_impl_t>::type;

template <typename algorithm_impl_t>
class _dp_algorithm_template_base<algorithm_impl_t>::type
{
private:
    using algorithm_attorney_t = dp_algorithm_attorney<algorithm_impl_t>;

protected:

    auto initialise_substitution_scheme() const noexcept
    {
        return algorithm_attorney_t::initialise_substitution_scheme(as_algorithm());
    }

    auto initialise_tracker() const noexcept
    {
        return algorithm_attorney_t::initialise_tracker(as_algorithm());
    }

    template <typename sequence1_t, typename dp_column_t>
    auto initialise_column(sequence1_t && sequence1, dp_column_t && dp_column) const noexcept
    {
        return algorithm_attorney_t::initialise_column_vector(as_algorithm(), sequence1, dp_column);
    }

    template <typename sequence2_t, typename dp_row_t>
    auto initialise_row(sequence2_t && sequence2, dp_row_t && dp_row) const noexcept
    {
        return algorithm_attorney_t::initialise_row_vector(as_algorithm(), sequence2, dp_row);
    }

    template <typename ...args_t>
    auto initialise_dp_matrix(args_t && ...args) const noexcept
    {
        return algorithm_attorney_t::initialise_policies(as_algorithm(),
                                                         std::forward<args_t>(args)...,
                                                         initialise_substitution_scheme(),
                                                         initialise_tracker());
    }

    template <typename ...args_t>
    constexpr auto lane_width(args_t && ...args) const noexcept
    {
        return algorithm_attorney_t::lane_width(as_algorithm(), std::forward<args_t>(args)...);
    }

    template <typename dp_block_t>
    void compute_block(dp_block_t && dp_block) const noexcept
    {
        // Initialise bulk_cache array.
        // constexpr std::ptrdiff_t lane_width = std::remove_cvref_t<dp_block_t>::lane_width_v;
        // std::ptrdiff_t const sequence1_size = std::ranges::distance(sequence1);

        constexpr auto index_sequence = std::make_index_sequence<std::remove_reference_t<dp_block_t>::lane_width>();
        auto && tracker = dp_matrix::tracker(dp_block);
        auto && scorer = dp_matrix::substitution_model(dp_block);

        // We are moving over the sequences here.
        for (std::ptrdiff_t lane_index = 0; lane_index < dp_matrix::column_count(dp_block) - 1; ++lane_index)
        {
            auto dp_lane = dp_matrix::column_at(dp_block, lane_index);
            auto && seq2_slice = dp_matrix::row_sequence(dp_lane);

            // compute cache many cells in one row for one horizontal value.
            for (std::ptrdiff_t i = 0; i < dp_matrix::row_count(dp_lane); ++i) {
                auto cacheH = dp_matrix::dp_column(dp_lane)[i+1];
                unroll_loop(dp_matrix::dp_row(dp_lane),
                            cacheH,
                            scorer,
                            tracker,
                            dp_matrix::column_sequence(dp_lane)[i],
                            seq2_slice,
                            index_sequence);
                dp_matrix::dp_column(dp_lane)[i+1] = cacheH;
            }
        }

        // Compute remaining cells requesting explicitly last lane.
        auto final_dp_lane = dp_block.final_lane(); // Not a CPO -> last_column/row
        auto && seq2_slice = dp_matrix::row_sequence(final_dp_lane);
        // assert(seq2_slice.size() <= lane_width);

        // compute cache many cells in one row for one horizontal value.
        for (std::ptrdiff_t i = 0; i < dp_matrix::row_count(final_dp_lane); ++i) {
            auto cacheH = dp_matrix::dp_column(final_dp_lane)[i+1];
            unroll_loop(dp_matrix::dp_row(final_dp_lane),
                        cacheH,
                        scorer,
                        tracker,
                        dp_matrix::column_sequence(final_dp_lane)[i],
                        seq2_slice);
            dp_matrix::dp_column(final_dp_lane)[i+1] = cacheH;
        }
    }

    template <typename tracker_t, typename ...args_t>
    auto make_result(tracker_t const & tracker, args_t && ...args) const noexcept
    {
        auto max_score = tracker.max_score(args...);
        return aligner_result(std::forward<args_t>(args)..., std::move(max_score));
    }

private:
    template <typename row_cells_t,
              typename col_cell_t,
              typename seq1_value_t,
              typename seq2_values_t,
              typename scorer_t,
              typename tracker_t,
              size_t ...idx>
    constexpr void unroll_loop(row_cells_t & row_cells,
                               col_cell_t & col_cell,
                               scorer_t const & scorer,
                               tracker_t & tracker,
                               seq1_value_t const & seq1_value,
                               seq2_values_t const & seq2_values,
                               [[maybe_unused]] std::index_sequence<idx...> const & indices) const noexcept
    {
        (algorithm_attorney_t::compute_cell(as_algorithm(),
                                            row_cells[idx],
                                            col_cell,
                                            scorer,
                                            tracker,
                                            seq1_value,
                                            seq2_values[idx]), ...);
    }

    template <typename row_cells_t,
              typename col_cell_t,
              typename scorer_t,
              typename tracker_t,
              typename seq1_value_t,
              typename seq2_values_t>
    constexpr void unroll_loop(row_cells_t & row_cells,
                               col_cell_t & col_cell,
                               scorer_t const & scorer,
                               tracker_t & tracker,
                               seq1_value_t const & seq1_value,
                               seq2_values_t const & seq2_values) const noexcept
    {
        for (size_t idx = 0; idx < seq2_values.size(); ++idx) {
            algorithm_attorney_t::compute_cell(as_algorithm(),
                                               row_cells[idx],
                                               col_cell,
                                               scorer,
                                               tracker,
                                               seq1_value,
                                               seq2_values[idx]);
        }
    }

    constexpr algorithm_impl_t const & as_algorithm() const noexcept
    {
        return static_cast<algorithm_impl_t const &>(*this);
    }

    constexpr algorithm_impl_t & as_algorithm() noexcept
    {
        return static_cast<algorithm_impl_t &>(*this);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
