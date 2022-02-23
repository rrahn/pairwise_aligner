// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::score_model_unitary_local.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t>
struct _score_model_unitary_local
{
    class type;
};

template <typename score_t>
using  score_model_unitary_local = typename _score_model_unitary_local<score_t>::type;

template <typename score_t>
class _score_model_unitary_local<score_t>::type
{
private:
    score_t _match_score{};
    score_t _mismatch_score{};
    score_t _zero{};
public:

    using score_type = score_t;

    type() = default;
    explicit type(score_t match_score, score_t mismatch_score, score_t zero = {}) :
        _match_score{std::move(match_score)},
        _mismatch_score{std::move(mismatch_score)},
        _zero{std::move(zero)}
    {}

    constexpr score_t const & zero() const noexcept
    {
        return _zero;
    }

    template <typename value1_t, typename value2_t>
        requires (std::equality_comparable_with<value1_t, value2_t>)
    score_t score(score_t const & last_diagonal,
                  value1_t const & value1,
                  value2_t const & value2) const noexcept
    {
        return last_diagonal + ((value1 == value2) ? _match_score : _mismatch_score);
    }

    template <std::integral scalar_score_t, size_t bulk_size>
    simd_score<scalar_score_t, bulk_size> score(simd_score<scalar_score_t, bulk_size> const & last_diagonal,
                                                simd_score<scalar_score_t, bulk_size> const & value1,
                                                simd_score<scalar_score_t, bulk_size> const & value2) const noexcept
    {
        static_assert(std::same_as<score_t, simd_score<scalar_score_t, bulk_size>>,
                      "The simd score type does not match the score type of this score class.");
        return blend(value1.eq(value2), _match_score, _mismatch_score) + last_diagonal;
    }

    // TODO: Refactor into separate factory CPO.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
