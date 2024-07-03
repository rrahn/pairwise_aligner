// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::aligner_result_batch.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <memory>
#include <seqan3/std/ranges>
#include <tuple>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename aligner_result_t>
class aligner_result_bulk
{
    static constexpr bool is_bulk_sequence_v =
        std::conditional_t<std::ranges::range<std::ranges::range_reference_t<decltype(std::declval<aligner_result_t>().sequence1())>>,
                           std::true_type,
                           std::false_type>::value;

    std::shared_ptr<aligner_result_t> _result{};
    size_t _index{};

public:

    aligner_result_bulk() = default;
    explicit aligner_result_bulk(std::shared_ptr<aligner_result_t> result, size_t const index) noexcept :
        _result{std::move(result)},
        _index{index}
    {}

    auto const & dp_column() const & noexcept
    {
        return _result->dp_column();
    }

    auto const & dp_row() const & noexcept
    {
        return _result->dp_row();
    }

    auto const & sequence1() const noexcept
    {
        if constexpr (is_bulk_sequence_v)
            return _result->sequence1()[_index];
        else
            return _result->sequence1();
    }

    auto const & sequence2() const noexcept
    {
        return _result->sequence2()[_index];
    }

    auto score() const noexcept
    {
        return _result->score()[_index];
    }
};

// namespace cpo
// {
// struct fn
// {
//     template <typename aligner_result_t>
//     auto operator()(std::shared_ptr<aligner_result_t> aligner,
//                     index const index) const noexcept
//     {
//         return value<aligner_result>{std::move(aligner), index};
//     }
// };

// } // namespace cpo
// } // namespace _result.

// inline constexpr _aligner_result_bulk::cpo::fn make_result_bulk{};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
