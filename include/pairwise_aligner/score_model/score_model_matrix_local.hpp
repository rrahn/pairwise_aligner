// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::score_model_matrix_local.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <concepts>

#include <pairwise_aligner/simd/simd_score_type.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename substitution_matrix_t>
struct _score_model_matrix_local
{
    class type;
};

template <typename substitution_matrix_t>
using  score_model_matrix_local = typename _score_model_matrix_local<substitution_matrix_t>::type;

template <typename substitution_matrix_t>
class _score_model_matrix_local<substitution_matrix_t>::type
{
private:

    substitution_matrix_t _substitution_matrix{};

    // TODO: Generic concept for substitution matrix.
    using score_t = typename substitution_matrix_t::value_type;

public:

    using score_type = score_t;

    type() = default;
    explicit type(substitution_matrix_t substitution_matrix) : _substitution_matrix{std::move(substitution_matrix)}
    {}

    template <typename value1_t, typename value2_t>
        requires (std::equality_comparable_with<value1_t, value2_t>)
    score_t score(score_t const & last_diagonal,
                  value1_t const & value1,
                  value2_t const & value2) const noexcept
    {
        using std::max;
        return max<score_t>(last_diagonal + _substitution_matrix[value1 + value2], 0);
    }

    // TODO: Refactor into separate factory CPO.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
