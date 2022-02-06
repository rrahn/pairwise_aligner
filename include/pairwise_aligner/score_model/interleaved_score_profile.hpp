// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::interleaved_score_profile.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename substitution_matrix_t, size_t stride_v>
struct _interleaved_score_profile
{
    class type;
};

template <typename substitution_matrix_t, size_t stride_v>
using interleaved_score_profile = typename _interleaved_score_profile<substitution_matrix_t, stride_v>::type;

template <typename substitution_matrix_t, size_t stride_v>
class _interleaved_score_profile<substitution_matrix_t, stride_v>::type
{
private:

    static constexpr size_t alphabet_size_v = substitution_matrix_t::dimension_v;

    using simd_score_t = typename substitution_matrix_t::score_type;
    using index_t = typename substitution_matrix_t::index_type;
    using offset_t = typename substitution_matrix_t::offset_type;
    using interleaved_scores_t = std::array<index_t, stride_v>;
    using profile_t = std::array<interleaved_scores_t, alphabet_size_v>;

    profile_t _profile;
    substitution_matrix_t const & _matrix;

public:

    type() = delete;
    template <typename sequence_slice_t>
    explicit type(substitution_matrix_t const & matrix, sequence_slice_t && sequence) noexcept : _matrix{matrix}
    {
        assert(std::ranges::size(sequence) <= stride_v); // can't be larger.
        // Initialise profile: - go over all symbols in range [0..sigma)
        for (size_t symbol = 0; symbol < alphabet_size_v; ++symbol) {
            interleaved_scores_t & profile = _profile[symbol]; // fill profile for current symbol!
            for (size_t index = 0; index < std::ranges::size(sequence); ++index) { // for every symbol in sequence
                profile[index] = profile[index] | _matrix[offset_t{symbol, sequence[index]}]; // scores for ranks at
            }
        }
    }

    constexpr interleaved_scores_t scores_for(index_t const & ranks) const noexcept
    {
        return for_each_symbol(ranks, std::make_index_sequence<alphabet_size_v>());
    }

private:

    template <size_t ...idx>
    constexpr auto for_each_symbol(index_t const & ranks, std::index_sequence<idx...> const &) const noexcept
    {
        constexpr auto profile_sequence = std::make_index_sequence<stride_v>();
        interleaved_scores_t scores{};
        (for_each_profile(scores, ranks.eq(index_t{idx}), _profile[idx], profile_sequence), ...);
        return scores;
    }

    template <typename simd_mask_t, size_t ...idx>
    constexpr void for_each_profile(interleaved_scores_t & scores,
                                    simd_mask_t const & mask,
                                    interleaved_scores_t const & profile,
                                    std::index_sequence<idx...> const &) const noexcept
    {
        ((scores[idx] = blend(mask, profile[idx], scores[idx])), ...);
    }

};
} // inline namespace v1
} // namespace seqan::pairwise_aligner
