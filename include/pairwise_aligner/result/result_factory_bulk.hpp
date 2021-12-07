// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::aligner_result.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/result/aligner_result.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace _bulk_factory
{
template <typename base_value_t, typename score_t>
struct _value
{
    struct type;
};

template <typename base_value_t, typename bulk_score_t>
using value = typename _value<base_value_t, bulk_score_t>::type;

template <typename base_value_t, typename bulk_score_t>
struct _value<base_value_t, bulk_score_t>::type : base_value_t
{
    using score_type = bulk_score_t;

    bulk_score_t _padding_score;
};

} // namespace _bulk_factory

template <typename score_t>
struct _result_factory_bulk
{
    struct type;
};

template <typename score_t>
using result_factory_bulk = typename _result_factory_bulk<score_t>::type;

template <typename score_t>
struct _result_factory_bulk<score_t>::type : result_factory_single
{
    score_t _padding_score;

    template <typename ...args_t>
    auto operator()(args_t && ...args) const noexcept
    {
        auto base = result_factory_single::operator()(std::forward<args_t>(args)...);
        return _bulk_factory::value<decltype(base), score_t>{{std::move(base)}, _padding_score};
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
