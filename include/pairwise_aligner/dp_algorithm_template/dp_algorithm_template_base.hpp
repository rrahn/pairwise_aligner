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

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_attorney.hpp>
#include <pairwise_aligner/result/aligner_result.hpp>
#include <pairwise_aligner/score_model/strip_width.hpp>
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
    auto initialise_policies(args_t && ...args) const noexcept
    {
        return algorithm_attorney_t::initialise_policies(as_algorithm(), std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr auto lane_width(args_t && ...args) const noexcept
    {
        return algorithm_attorney_t::lane_width(as_algorithm(), std::forward<args_t>(args)...);
    }

    template <typename sequence1_t, typename dp_block_t>
    void compute_block(sequence1_t && sequence1, dp_block_t && dp_block) const noexcept
    {
        // Initialise bulk_cache array.
        constexpr std::ptrdiff_t lane_width = std::remove_cvref_t<dp_block_t>::lane_width_v;
        constexpr auto index_sequence = std::make_index_sequence<lane_width>();
        std::ptrdiff_t const sequence1_size = std::ranges::distance(sequence1);

        auto && tracker = dp_block.tracker();
        auto && scorer = dp_block.substitution_model();

        for (size_t lane_index = 0; lane_index < dp_block.size() - 1; ++lane_index)
        {
            auto dp_lane = dp_block[lane_index];
            auto && seq2_slice = dp_lane.row_sequence();

            // compute cache many cells in one row for one horizontal value.
            for (std::ptrdiff_t i = 0; i < sequence1_size; ++i) {
                auto cacheH = dp_lane.column()[i+1];
                unroll_loop(dp_lane.row(), cacheH, scorer, tracker, sequence1[i], seq2_slice, index_sequence);
                dp_lane.column()[i+1] = cacheH;
            }
        }

        // Compute remaining cells requesting explicitly last lane.
        auto last_dp_lane = dp_block.last_lane();
        auto && seq2_slice = last_dp_lane.row_sequence();
        assert(seq2_slice.size() <= lane_width);

        // compute cache many cells in one row for one horizontal value.
        for (std::ptrdiff_t i = 0; i < sequence1_size; ++i) {
            auto cacheH = last_dp_lane.column()[i+1];
            unroll_loop(last_dp_lane.row(), cacheH, scorer, tracker, sequence1[i], seq2_slice);
            last_dp_lane.column()[i+1] = cacheH;
        }
    }

    template <typename sequence1_t,
              typename sequence2_t,
              typename dp_column_t,
              typename dp_row_t,
              typename scorer_t,
              typename tracker_t>
        // requires requires (scorer_t const & sc, sequence2_t && seq)
        // {
        //     { sc.initialise_profile(seq, strip_width<1>) };
        // }
    void compute_block(sequence1_t && sequence1,
                       sequence2_t && sequence2,
                       dp_column_t && dp_column,
                       dp_row_t && dp_row,
                       scorer_t && scorer,
                       tracker_t && tracker)
        const noexcept
    {
        std::ptrdiff_t const sequence1_size = std::ranges::distance(sequence1);
        std::ptrdiff_t const sequence2_size = std::ranges::distance(sequence2);

        // Precompute the offsets for the row of this block
        // using s1_value_t = typename std::ranges::range_value_t<sequence1_t>;
        // using uvalue_t = detail::make_unsigned_t<s1_value_t>;

        // auto profile_init = scorer.initialise_profile(std::span{sequence2.begin(), 1}, strip_width<1>);
        // std::vector<uvalue_t> tmp{};
        // tmp.resize(sequence1_size);
        // for (size_t i = 0; i < tmp.size(); ++i) {
        //     tmp[i] = profile_init.to_offset(sequence1[i]);
        // }

        using value_t = typename std::remove_cvref_t<dp_row_t>::value_type;

        // Store score of first row cell in first column cell.
        // Note the first column/row is not computed again, as they were already initialised.
        dp_column[0].score() = dp_row[0].score();

        // Initialise bulk_cache array.
        constexpr std::ptrdiff_t cache_size = (detail::max_simd_size == 64) ? 8 : 4;
        std::array<value_t, cache_size> bulk_cache{};

        std::ptrdiff_t j = 0;
        for (; j < sequence2_size - (cache_size - 1); j += cache_size)
        {
            // copy values into cache.
            std::span seq2_slice{sequence2.begin() + j, sequence2.begin() + j + cache_size};

            unroll_load(bulk_cache, dp_row, j + 1, std::make_index_sequence<cache_size>());
            auto profile = scorer.initialise_profile(seq2_slice, strip_width<cache_size>); // initialise_profile<cache_size>(scorer, seq2_slice);

            // compute cache many cells in one row for one horizontal value.
            for (std::ptrdiff_t i = 0; i < sequence1_size; ++i) {
                auto cacheH = dp_column[i+1];
                auto profile_scores = profile.scores_for(sequence1[i]);

                unroll_loop(bulk_cache,
                            cacheH,
                            scorer,
                            tracker,
                            sequence1[i],
                            profile_scores,
                            std::make_index_sequence<cache_size>());

                dp_column[i+1] = cacheH;
            }

            // store results of cache back in vector.
            unroll_store(dp_row, bulk_cache, j + 1, std::make_index_sequence<cache_size>());
        }

        // Compute remaining cells.
        // Run last block with less than 4/8 cached values.
        std::span seq2_slice{sequence2.begin() + j, sequence2.end()};
        assert(seq2_slice.size() <= cache_size);

        // unroll_load
        for (std::ptrdiff_t k = 0 ; k < std::ranges::ssize(seq2_slice); ++k) {
            bulk_cache[k] = dp_row[j + k + 1];
        }

        // auto profile = scorer.initialise_profile(seq2_slice, strip_width<cache_size>);

        // compute cache many cells in one row for one horizontal value.
        for (std::ptrdiff_t i = 0; i < sequence1_size; ++i) {
            auto cacheH = dp_column[i+1];
            // auto profile_scores = profile.scores_for(sequence1[i]);

            for (std::ptrdiff_t k = 0; k < std::ranges::ssize(seq2_slice); ++k) {
                algorithm_attorney_t::compute_cell(as_algorithm(),
                                                bulk_cache[k],
                                                cacheH,
                                                scorer,
                                                tracker,
                                                sequence1[i],
                                                seq2_slice[k]);
            }

            dp_column[i+1] = cacheH;
        }

        // unroll_store
        for (std::ptrdiff_t k = 0 ; k < std::ranges::ssize(seq2_slice); ++k) {
            dp_row[j + k + 1] = bulk_cache[k];
        }
        // Store score of last column in first cell of row.
        dp_row[0].score() = dp_column[dp_column.size() - 1].score();
    }

    template <typename dp_row_t>
    constexpr void rotate_row_scores_right(dp_row_t && dp_row) const noexcept
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

    template <typename dp_row_t>
    constexpr void rotate_row_scores_left(dp_row_t && dp_row) const noexcept
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

    template <typename cache_t, typename row_vector_t, size_t ...idx>
    constexpr void unroll_load(cache_t & bulk_cache,
                               row_vector_t const & row_vector,
                               size_t const offset,
                               [[maybe_unused]] std::index_sequence<idx...> const & indices) const noexcept
    {
        ((bulk_cache[idx] = row_vector[offset + idx]), ...);
    }

    template <typename row_vector_t, typename cache_t, size_t ...idx>
    constexpr void unroll_store(row_vector_t & row_vector,
                                cache_t const & bulk_cache,
                                size_t const offset,
                                [[maybe_unused]] std::index_sequence<idx...> const & indices) const noexcept
    {
        ((row_vector[offset + idx] = bulk_cache[idx]), ...);
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
