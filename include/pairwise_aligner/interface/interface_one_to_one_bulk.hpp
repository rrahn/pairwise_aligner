// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::interface_one_to_one_bulk.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <array>
#include <memory>
#include <optional>
#include <seqan3/std/ranges>

#include <pairwise_aligner/result/aligner_result_bulk.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_algorithm_t, size_t max_bulk_size>
struct _interface_one_to_one_bulk
{
    struct type;
};

template <typename dp_algorithm_t, size_t max_bulk_size>
using interface_one_to_one_bulk = typename _interface_one_to_one_bulk<dp_algorithm_t, max_bulk_size>::type;

template <typename dp_algorithm_t, size_t max_bulk_size>
struct _interface_one_to_one_bulk<dp_algorithm_t, max_bulk_size>::type : protected dp_algorithm_t
{

    type() = default;
    explicit type(dp_algorithm_t algorithm) noexcept : dp_algorithm_t{std::move(algorithm)}
    {}
    using dp_algorithm_t::dp_algorithm_t;

    using dp_algorithm_t::column_vector;
    using dp_algorithm_t::row_vector;

    template <std::ranges::forward_range sequence_bulk1_t,
              std::ranges::forward_range sequence_bulk2_t>
        requires (std::ranges::forward_range<std::ranges::range_reference_t<sequence_bulk1_t>> &&
                  std::ranges::forward_range<std::ranges::range_reference_t<sequence_bulk2_t>>) &&
                 (std::ranges::viewable_range<std::ranges::range_reference_t<sequence_bulk1_t>> &&
                  std::ranges::viewable_range<std::ranges::range_reference_t<sequence_bulk2_t>>)
    auto compute(sequence_bulk1_t && sequence_bulk1, sequence_bulk2_t && sequence_bulk2)
    {
        assert(static_cast<size_t>(std::ranges::distance(sequence_bulk1)) <= max_bulk_size);
        assert(static_cast<size_t>(std::ranges::distance(sequence_bulk2)) <= max_bulk_size);
        assert(std::ranges::distance(sequence_bulk1) == std::ranges::distance(sequence_bulk2));

        return compute(std::forward<sequence_bulk1_t>(sequence_bulk1),
                       std::forward<sequence_bulk2_t>(sequence_bulk2),
                       column_vector(),
                       row_vector());
    }

    template <std::ranges::forward_range sequence_bulk1_t,
              std::ranges::forward_range sequence_bulk2_t,
              typename dp_column_t,
              typename dp_row_t>
        requires (std::ranges::forward_range<std::ranges::range_reference_t<sequence_bulk1_t>> &&
                  std::ranges::forward_range<std::ranges::range_reference_t<sequence_bulk2_t>>) &&
                 (std::ranges::viewable_range<std::ranges::range_reference_t<sequence_bulk1_t>> &&
                  std::ranges::viewable_range<std::ranges::range_reference_t<sequence_bulk2_t>>)
    auto compute(sequence_bulk1_t && sequence_bulk1,
                 sequence_bulk2_t && sequence_bulk2,
                 dp_column_t first_dp_column,
                 dp_row_t first_dp_row)
    {
        assert(static_cast<size_t>(std::ranges::distance(sequence_bulk1)) <= max_bulk_size);
        assert(static_cast<size_t>(std::ranges::distance(sequence_bulk2)) <= max_bulk_size);
        assert(std::ranges::distance(sequence_bulk1) == std::ranges::distance(sequence_bulk2));

        // how can we handle wrong types?
        using sequence1_ref_t = std::ranges::range_reference_t<sequence_bulk1_t>;
        using sequence2_ref_t = std::ranges::range_reference_t<sequence_bulk2_t>;
        using bulk1_t = std::array<std::optional<std::views::all_t<sequence1_ref_t>>, max_bulk_size>;
        using bulk2_t = std::array<std::optional<std::views::all_t<sequence2_ref_t>>, max_bulk_size>;

        bulk1_t bulk1{};
        bulk2_t bulk2{};

        auto sequence1_it = std::ranges::begin(sequence_bulk1);
        auto sequence2_it = std::ranges::begin(sequence_bulk2);
        size_t sequence_idx{};
        for (; sequence1_it != std::ranges::end(sequence_bulk1); ++sequence_idx, ++sequence1_it, ++sequence2_it)
        {
            bulk1[sequence_idx] = *sequence1_it | std::views::all;
            bulk2[sequence_idx] = *sequence2_it | std::views::all;
        }

        auto result = dp_algorithm_t::run(std::move(bulk1),
                                          std::move(bulk2),
                                          std::move(first_dp_column),
                                          std::move(first_dp_row));

        // Transfrom bulk result to single results.
        using result_t = decltype(result);
        std::vector<aligner_result_bulk<result_t>> results{};
        results.reserve(sequence_idx);

        auto shared_result = std::make_shared<result_t>(std::move(result));

        for (size_t result_idx = 0; result_idx < sequence_idx; ++result_idx)
            results.emplace_back(shared_result, result_idx);

        return results;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
