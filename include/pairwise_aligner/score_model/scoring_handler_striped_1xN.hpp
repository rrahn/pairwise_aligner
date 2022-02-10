// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::scoring_handler_striped_1xN.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename substitution_matrix_t, size_t strip_width_v>
struct _scoring_handler_striped_1xN
{
    class type;
};

template <typename substitution_matrix_t, size_t strip_width_v>
using scoring_handler_striped_1xN = typename _scoring_handler_striped_1xN<substitution_matrix_t, strip_width_v>::type;

template <typename substitution_matrix_t, size_t strip_width_v>
class _scoring_handler_striped_1xN<substitution_matrix_t, strip_width_v>::type
{
    static constexpr size_t alphabet_size_v = substitution_matrix_t::dimension_v;

    using simd_score_t = typename substitution_matrix_t::score_type;
    using index_t = typename substitution_matrix_t::index_type;
    using scalar_index_t = typename index_t::value_type;

    using interleaved_scores_t = std::array<index_t, strip_width_v>;
    using profile_t = std::array<interleaved_scores_t, alphabet_size_v>;

    profile_t _interleaved_profile;

public:

    type() = delete;
    template <typename sequence_slice_t>
    explicit type(substitution_matrix_t const & matrix, sequence_slice_t && sequence) noexcept
    {
        using offset_t = typename substitution_matrix_t::offset_type;

        assert(std::ranges::size(sequence) <= strip_width_v); // can't be larger, but smaller.
        // Initialise profile: - go over all symbols in range [0..sigma)
        for_each_symbol([&] (scalar_index_t symbol) {
            interleaved_scores_t & profile = _interleaved_profile[symbol]; // fill profile for current symbol!
            for (size_t index = 0; index < std::ranges::size(sequence); ++index) { // for every symbol in sequence
                profile[index] = profile[index] | matrix[offset_t{symbol, sequence[index]}]; // scores for ranks at
            }
        }, std::make_index_sequence<alphabet_size_v>());
        // for (size_t symbol = 0; symbol < alphabet_size_v; ++symbol) {

        // }
    }

    constexpr auto scores_for(scalar_index_t const & rank) const noexcept
    {
        return _interleaved_profile[rank];
    }

private:
    template <typename fn_t, size_t ...alphabet_rank>
    void for_each_symbol(fn_t && fn, std::index_sequence<alphabet_rank...> const &) const noexcept
    {
        (fn(static_cast<scalar_index_t>(alphabet_rank)), ...);
    }
};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
