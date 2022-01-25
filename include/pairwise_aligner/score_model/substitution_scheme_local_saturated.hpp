// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::substitution_scheme_local_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/concepts>

#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t>
struct _substitution_scheme_local_saturated
{
    class type;
};

template <typename score_t>
using  substitution_scheme_local_saturated = typename _substitution_scheme_local_saturated<score_t>::type;

template <typename score_t>
class _substitution_scheme_local_saturated<score_t>::type
{
private:
    score_t _match_score{};
    score_t _mismatch_score{};

    struct _local_scheme
    {
        score_t _match_score{};
        score_t _mismatch_score{};
        score_t _threshold{};

        template <std::integral scalar_score_t, size_t bulk_size>
        simd_score<scalar_score_t, bulk_size> score(simd_score<scalar_score_t, bulk_size> const & last_diagonal,
                                                    simd_score<scalar_score_t, bulk_size> const & value1,
                                                    simd_score<scalar_score_t, bulk_size> const & value2) const noexcept
        {
            static_assert(std::same_as<score_t, simd_score<scalar_score_t, bulk_size>>,
                        "The simd score type does not match the score type of this score class.");

            using simd_t = simd_score<scalar_score_t, bulk_size>;

            // initialise once per block:
            // block_mask -> 1: zero_out 0: do not zero_out
            // _threshold = blend(block_mask, _absolute_mismatch_score, lowest_viable_score (never reachable by global block))
            // what changes? the absolute mismatch score

            // keep instruction but on initialisation set the _threshold
            auto k1 = value1.eq(value2);
            return mask_add(simd_t{},
                            k1 | _threshold.lt(last_diagonal),
                            blend(k1, _match_score, _mismatch_score),
                            last_diagonal);

                // auto mismatch_mask =
                //     compare(value1, value2, [] (simd_t const & lhs, simd_t const & rhs) {
                //         return (lhs ^ rhs).ne(simd_t{});
                //     });

                // auto zero_mask = mask_lt(mismatch_mask, last_diagonal, _absolute_mismatch_score);
                // return blend(zero_mask, simd_t{}, blend(mismatch_mask, _mismatch_score, _match_score) + last_diagonal);


                // auto tmp_score = blend(compare(value1, value2, [] (simd_t const & lhs, simd_t const & rhs) {
                //                 return (lhs ^ rhs).eq(simd_t{});
                //             }), _match_score, _mismatch_score) + last_diagonal;
                // return blend(tmp_score.lt(simd_t{}), simd_t{}, tmp_score);
        }
    };

    struct _global_scheme
    {
        score_t _match_score{};
        score_t _mismatch_score{};

        template <std::integral scalar_score_t, size_t bulk_size>
        simd_score<scalar_score_t, bulk_size> score(simd_score<scalar_score_t, bulk_size> const & last_diagonal,
                                                    simd_score<scalar_score_t, bulk_size> const & value1,
                                                    simd_score<scalar_score_t, bulk_size> const & value2) const noexcept
        {
            static_assert(std::same_as<score_t, simd_score<scalar_score_t, bulk_size>>,
                        "The simd score type does not match the score type of this score class.");

            using simd_t = simd_score<scalar_score_t, bulk_size>;

            return blend(compare(value1, value2, [] (simd_t const & lhs, simd_t const & rhs) {
                            return (lhs ^ rhs).le(simd_t{});
                        }), _match_score, _mismatch_score) + last_diagonal;
        }
    };
public:

    using score_type = score_t;

    type() = default;
    explicit type(score_t match_score, score_t mismatch_score) :
        _match_score{std::move(match_score)},
        _mismatch_score{std::move(mismatch_score)}
        // _absolute_mismatch_score{_mismatch_score * -1}
    {}

    constexpr _local_scheme local_scheme() const noexcept
    {
        return _local_scheme{_match_score, _mismatch_score, _mismatch_score * -1};
    }

    constexpr _global_scheme global_scheme() const noexcept
    {
        return _global_scheme{_match_score, _mismatch_score};
    }

    // creates itself.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
