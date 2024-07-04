// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::score_model_unitary.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <pairwise_aligner/utility/math.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t>
struct _score_model_unitary
{
    class type;
};

template <typename score_t>
    requires (std::is_arithmetic_v<score_t>)
using  score_model_unitary = typename _score_model_unitary<score_t>::type;

template <typename score_t>
class _score_model_unitary<score_t>::type
{
private:
    score_t _match_score{};
    score_t _mismatch_score{};

public:

    using score_type = score_t;

    type() = default;
    explicit type(score_t match_score, score_t mismatch_score) :
        _match_score{std::move(match_score)},
        _mismatch_score{std::move(mismatch_score)}
    {}

    template <typename value1_t, typename value2_t>
        requires (std::equality_comparable_with<value1_t, value2_t>)
    score_t score(score_t const & last_diagonal,
                  value1_t const & value1,
                  value2_t const & value2) const noexcept
    {
        return add(last_diagonal, ((value1 == value2) ? _match_score : _mismatch_score));
    }

    // TODO: Refactor into separate factory CPO.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
