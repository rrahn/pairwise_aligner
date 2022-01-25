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

#include <seqan3/std/ranges>

#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_attorney.hpp>
#include <pairwise_aligner/result/aligner_result.hpp>

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

    template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t, typename tracker_t>
    void compute_block(sequence1_t && sequence1,
                       sequence2_t && sequence2,
                       dp_column_t && dp_column,
                       dp_row_t && dp_row,
                       tracker_t && tracker)
        const noexcept
    {
        std::ptrdiff_t const sequence1_size = std::ranges::distance(sequence1);
        std::ptrdiff_t const sequence2_size = std::ranges::distance(sequence2);

        using value_t = typename std::remove_cvref_t<dp_row_t>::value_type;

        // Initialise bulk_cache array.
        constexpr std::ptrdiff_t cache_size = 8;
        std::array<value_t, cache_size> bulk_cache{};

        std::ptrdiff_t j = 0;
        for (; j < sequence2_size - (cache_size - 1); j += cache_size)
        {
            // copy values into cache.
            std::span seq2_slice{sequence2.begin() + j, sequence2.begin() + j + cache_size};
            unroll_load(bulk_cache, dp_row, j + 1, std::make_index_sequence<cache_size>());

            // compute cache many cells in one row for one horizontal value.
            for (std::ptrdiff_t i = 0; i < sequence1_size; ++i) {
                auto cacheH = dp_column[i+1];

                unroll_loop(bulk_cache,
                            cacheH,
                            tracker,
                            sequence1[i],
                            seq2_slice,
                            std::make_index_sequence<cache_size>());

                dp_column[i+1] = cacheH;
            }

            // store results of cache back in vector.
            unroll_store(dp_row, bulk_cache, j + 1, std::make_index_sequence<cache_size>());
        }

        // Compute remaining cells.
        for (; j < sequence2_size; ++j) {
            auto cache = dp_row[j + 1];

            for (std::ptrdiff_t i = 0; i < sequence1_size; ++i) {
                algorithm_attorney_t::compute_cell(as_algorithm(),
                                                   cache,
                                                   dp_column[i+1],
                                                   tracker,
                                                   sequence1[i],
                                                   sequence2[j]);
            }

            dp_row[j + 1] = cache;
        }
    }

    template <typename dp_column_t, typename dp_row_t>
    constexpr void initialise_block(dp_column_t && dp_column, dp_row_t && dp_row) const noexcept
    {
        // H'[0].score() = V'[m].score();
        // V'[0].score() = H'[n].score();
        size_t const dp_row_size = dp_row.size() - 1;
        dp_column[0].score() = dp_row[dp_row_size].score();

        for (size_t j = dp_row_size; j > 0; --j)
            dp_row[j].score() = dp_row[j - 1].score();
    }

    template <typename dp_column_t, typename dp_row_t>
    constexpr void postprocess_block(dp_column_t const & dp_column, dp_row_t && dp_row) const noexcept
    {
        size_t const dp_row_size = dp_row.size() - 1;
        for (size_t j = 0; j < dp_row_size; ++j)
            dp_row[j].score() = dp_row[j + 1].score();

        dp_row[dp_row_size].score() = dp_column[dp_column.size() - 1].score();
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
              typename tracker_t,
              size_t ...idx>
    constexpr void unroll_loop(row_cells_t & row_cells,
                               col_cell_t & col_cell,
                               tracker_t & tracker,
                               seq1_value_t const & seq1_value,
                               seq2_values_t const & seq2_values,
                               [[maybe_unused]] std::index_sequence<idx...> const & indices) const noexcept
    {
        (algorithm_attorney_t::compute_cell(as_algorithm(),
                                            row_cells[idx],
                                            col_cell,
                                            tracker,
                                            seq1_value,
                                            seq2_values[idx]), ...);
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
