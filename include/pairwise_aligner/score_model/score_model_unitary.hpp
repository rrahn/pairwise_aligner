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

#include <seqan3/std/concepts>

#include <pairwise_aligner/simd/simd_score_type.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t, bool is_local>
struct _score_model_unitary
{
    class type;
};

template <typename score_t, bool is_local>
using  score_model_unitary = typename _score_model_unitary<score_t, is_local>::type;

template <typename score_t, bool is_local>
class _score_model_unitary<score_t, is_local>::type
{
private:
    score_t _match_score{};
    score_t _mismatch_score{};
    score_t _absolute_mismatch_score{};

public:

    using score_type = score_t;

    type() = default;
    explicit type(score_t match_score, score_t mismatch_score) :
        _match_score{std::move(match_score)},
        _mismatch_score{std::move(mismatch_score)},
        _absolute_mismatch_score{_mismatch_score * -1}
    {}

    template <typename value1_t, typename value2_t>
        requires (std::equality_comparable_with<value1_t, value2_t>)
    score_t score(score_t const & last_diagonal,
                  value1_t const & value1,
                  value2_t const & value2) const noexcept
    {
        if constexpr (is_local)
            return (value1 == value2) ? _match_score + last_diagonal :
                        (last_diagonal < -_mismatch_score ? 0 : _mismatch_score + last_diagonal);
        else
            return last_diagonal + score_impl(value1, value2);
    }

    template <std::integral scalar_score_t, size_t bulk_size>
    simd_score<scalar_score_t, bulk_size> score(simd_score<scalar_score_t, bulk_size> const & last_diagonal,
                                                simd_score<scalar_score_t, bulk_size> const & value1,
                                                simd_score<scalar_score_t, bulk_size> const & value2) const noexcept
    {
        static_assert(std::same_as<score_t, simd_score<scalar_score_t, bulk_size>>,
                      "The simd score type does not match the score type of this score class.");

        using simd_t = simd_score<scalar_score_t, bulk_size>;

        if constexpr (is_local) {
            // auto mismatch_mask =
            //     compare(value1, value2, [] (simd_t const & lhs, simd_t const & rhs) {
            //         return (lhs ^ rhs).ne(simd_t{});
            //     });

            // auto zero_mask = mask_lt(mismatch_mask, last_diagonal, _absolute_mismatch_score);
            // return blend(zero_mask, simd_t{}, blend(mismatch_mask, _mismatch_score, _match_score) + last_diagonal);

            auto k1 = value1.eq(value2);
            return mask_add(simd_t{},
                            k1 | _absolute_mismatch_score.lt(last_diagonal),
                            blend(k1, _match_score, _mismatch_score),
                            last_diagonal);

            // auto tmp_score = blend(compare(value1, value2, [] (simd_t const & lhs, simd_t const & rhs) {
            //                 return (lhs ^ rhs).eq(simd_t{});
            //             }), _match_score, _mismatch_score) + last_diagonal;
            // return blend(tmp_score.lt(simd_t{}), simd_t{}, tmp_score);

        } else {
            return blend(compare(value1, value2, [] (simd_t const & lhs, simd_t const & rhs) {
                            return (lhs ^ rhs).le(simd_t{});
                        }), _match_score, _mismatch_score) + last_diagonal;
        }
    }

    constexpr score_t padding_score() const noexcept
    {
        return _match_score;
    }

    // creates itself.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }

private:

    template <typename value1_t, typename value2_t>
        requires (std::equality_comparable_with<value1_t, value2_t>)
    score_t score_impl(value1_t const & value1, value2_t const & value2) const noexcept
    {
        return (value1 == value2) ? _match_score : _mismatch_score;
    }

    // template <std::integral scalar_score_t, size_t bulk_size>
    // simd_score<scalar_score_t, bulk_size> score_impl(simd_score<scalar_score_t, bulk_size> const & last_diagonal,
    //                                                  simd_score<scalar_score_t, bulk_size> const & value1,
    //                                                  simd_score<scalar_score_t, bulk_size> const & value2) const noexcept
    // {
    //     static_assert(std::same_as<score_t, simd_score<scalar_score_t, bulk_size>>,
    //                   "The simd score type does not match the score type of this score class.");

    //     using simd_t = simd_score<scalar_score_t, bulk_size>;

    //     return blend(compare(value1, value2, [] (simd_t const & lhs, simd_t const & rhs) {
    //         return (lhs ^ rhs).le(simd_t{0});
    //     }), _match_score, _mismatch_score);
    // }


};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
