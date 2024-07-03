// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::score_model_unitary_simd.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/concepts>

#include <pairwise_aligner/simd/concept.hpp>
#include <pairwise_aligner/utility/math.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t>
struct _score_model_unitary_simd
{
    class type;
};

template <simd::simd_type score_t>
using  score_model_unitary_simd = typename _score_model_unitary_simd<score_t>::type;

template <typename score_t>
class _score_model_unitary_simd<score_t>::type
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

    template <simd::simd_type value_t>
        requires (std::same_as<typename score_type::mask_type, typename value_t::mask_type>)
    score_type score(score_type const & last_diagonal, value_t const & value1, value_t const & value2) const noexcept
    {
        return add(blend(compare(value1, value2, [] (score_type const & lhs, score_type const & rhs) {
                        return (lhs ^ rhs).le(score_type{});
                     }), _match_score, _mismatch_score),  last_diagonal);
    }

    // TODO: Refactor into separate factory CPO.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
