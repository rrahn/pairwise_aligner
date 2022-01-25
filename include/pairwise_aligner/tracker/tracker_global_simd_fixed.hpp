// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::tracker::global_simd_fixed.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/configuration/end_gap_policy.hpp>
#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace tracker::global_simd_fixed {

template <typename score_t>
struct _tracker
{
    class type;
};

template <typename score_t>
using tracker = typename _tracker<score_t>::type;

template <typename score_t>
class _tracker<score_t>::type
{
public:

    score_t _padding_score;
    cfg::trailing_end_gap _end_gap;

    constexpr score_t const & track(score_t const & score) const noexcept {
        return score; // no-op.
    }

    // receive the optimal score.
    template <typename sequences1_t, typename sequences2_t, typename dp_column_t, typename dp_row_t>
    constexpr score_t max_score(sequences1_t && sequences1,
                                sequences2_t && sequences2,
                                dp_column_t const & dp_column,
                                dp_row_t const & dp_row) const noexcept {

        if (_end_gap.last_column == cfg::end_gap::penalised && _end_gap.last_row == cfg::end_gap::penalised) {
            score_t best_score{};
            for (std::ptrdiff_t idx = 0; idx < std::ranges::distance(sequences1); ++idx) {
                best_score[idx] = select_max_score(idx, sequences1[idx], sequences2[idx], dp_column, dp_row);
            }
            return best_score;
        }

        score_t best_score{std::numeric_limits<typename score_t::value_type>::lowest()};
        if (_end_gap.last_column == cfg::end_gap::free) {
            best_score = find_max_score(sequences1, dp_column, sequences2, dp_row);
        }

        if (_end_gap.last_row == cfg::end_gap::free) {
            best_score = max(best_score, find_max_score(sequences2, dp_row, sequences1, dp_column));
        }

        return best_score;
    }

    // TODO: optimal_coordinate()

private:


    template <typename sequence_t, typename dp_vector_t>
    constexpr auto get_offsets(sequence_t && sequence, dp_vector_t && dp_vector) const noexcept
    {
        size_t sequence_size = std::ranges::distance(sequence);
        size_t dp_vector_size = dp_vector.size();

        assert(dp_vector_size > sequence_size);
        size_t offset = dp_vector_size - 1 - sequence_size;

        return std::tuple{sequence_size, dp_vector_size, offset};
    }

    template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
    constexpr auto select_max_score(size_t const simd_idx,
                                    sequence1_t && sequence1,
                                    sequence2_t && sequence2,
                                    dp_column_t const & dp_column,
                                    dp_row_t const & dp_row) const noexcept
        -> typename score_t::value_type
    {
        auto [column_sequence_size, column_vector_size, row_offset] = get_offsets(sequence1, dp_column);
        auto [row_sequence_size, row_vector_size, column_offset] = get_offsets(sequence2, dp_row);

        using scalar_t = typename score_t::value_type;
        size_t scale = std::min(column_offset, row_offset);
        scalar_t best_score{};

        if (scale == column_offset) {
            best_score = score_at(dp_column[column_sequence_size + column_offset], simd_idx);
        } else {
            best_score = score_at(dp_row[row_sequence_size + row_offset], simd_idx);
        }

        return best_score - (_padding_score[simd_idx] * scale);
    }

    template <typename first_sequence_t, typename first_vector_t, typename second_sequence_t, typename second_vector_t>
    constexpr auto find_max_score(first_sequence_t && first_sequence,
                                  first_vector_t && first_vector,
                                  second_sequence_t && second_sequence,
                                  second_vector_t && second_vector) const noexcept
    {
        using scalar_t = typename score_t::value_type;
        using unsigned_scalar_t = std::make_unsigned_t<scalar_t>;
        using offset_simd_t = simd_score<unsigned_scalar_t, score_t::size>;

        offset_simd_t start_offset_first{};
        offset_simd_t end_offset_first{};
        score_t scale_first{};

        offset_simd_t start_offset_second{};
        offset_simd_t end_offset_second{};
        score_t scale_second{};

        // Prepare the offsets.
        for (std::ptrdiff_t idx = 0; idx < std::ranges::distance(first_sequence); ++idx) {
            auto [first_sequence_size, first_vector_size, second_offset] = get_offsets(first_sequence[idx], first_vector);
            auto [second_sequence_size, second_vector_size, first_offset] = get_offsets(second_sequence[idx], second_vector);

            start_offset_first[idx] = first_offset;
            end_offset_first[idx] = std::min<size_t>(first_offset + first_sequence_size + 1, first_vector_size);
            scale_first[idx] = first_offset;

            start_offset_second[idx] = second_sequence_size + second_offset;
            end_offset_second[idx] = second_vector_size;
            scale_second[idx] = second_offset;
        }

        // Iterate over the corresponding slice in the projected dp vector (first_vector) and subtract the scaled
        // padding score from the retrieved values of the corresponding simd index.

        score_t best_score{std::numeric_limits<scalar_t>::lowest()};
        score_t scale = _padding_score * scale_first;
        for (unsigned_scalar_t idx = 0; idx < first_vector.size(); ++idx) {
            offset_simd_t simd_idx{idx};

            auto mask = (start_offset_first.le(simd_idx) && simd_idx.lt(end_offset_first));
            best_score = mask_max(best_score, mask, best_score, first_vector[idx].score() - scale);
        }

        // Note if second_offset is less than first_offset, the last (first_offset - second_offset) elements of the
        // second vector must be considered as well; these cells contain the values of the projected first vector, which
        // breaks around the cell (n, m) of the extended simd matrix.
        // This slice starts at the projected second vector cell.
        scale = _padding_score * scale_second;
        for (unsigned_scalar_t idx = 0; idx < second_vector.size(); ++idx) {
            offset_simd_t simd_idx{idx};

            auto mask = (start_offset_second.le(simd_idx) && simd_idx.lt(end_offset_second));
            best_score = mask_max(best_score, mask, best_score, second_vector[idx].score() - scale);
            scale = mask_add(scale, mask, scale, _padding_score);
        }

        return best_score;
    }

    template <typename cell_t>
        requires requires (cell_t const & cell, size_t const idx){ { cell.score_at(idx) } -> std::integral; }
    constexpr auto score_at(cell_t const & cell, size_t const idx) const noexcept
    {
        return cell.score_at(idx);
    }

    template <typename cell_t>
    constexpr auto score_at(cell_t const & cell, size_t const idx) const noexcept
    {
        return cell.score()[idx];
    }
};

// this is the capture object which creates a new instance of a tracker every time we call it.
// the client can provide everything here.
// I need the factory here but some value to create it.
template <typename score_t>
struct _factory
{
    struct type;
};

template <typename score_t>
using factory = typename _factory<score_t>::type;

template <typename score_t>
struct _factory<score_t>::type
{
    // params for free end-gaps.
    score_t _padding_score{};
    cfg::trailing_end_gap _end_gap{};

    constexpr auto make_tracker() const noexcept {
        return tracker<score_t>{_padding_score, _end_gap};
    }
};

} // namespace tracker::global_simd_fixed

} // inline namespace v1
} // namespace seqan::pairwise_aligner
