// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::bulk_factory_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <limits>

#include <pairwise_aligner/configuration/end_gap_policy.hpp>
#include <pairwise_aligner/result/aligner_result.hpp>
#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace _bulk_factory_saturated
{
template <typename base_value_t, typename score_t>
struct _value
{
    struct type;
};

template <typename base_value_t, typename score_t>
using value = typename _value<base_value_t, score_t>::type;

template <typename base_value_t, typename score_t>
struct _value<base_value_t, score_t>::type : public base_value_t
{
    score_t _padding_score;

    constexpr score_t const & score() const noexcept
    {
        return _padding_score;
    }
};

} // namespace _bulk_factory_saturated

template <typename score_t>
struct _result_factory_bulk_saturated
{
    struct type;
};

template <typename score_t>
using result_factory_bulk_saturated = typename _result_factory_bulk_saturated<score_t>::type;

template <typename score_t>
struct _result_factory_bulk_saturated<score_t>::type
{
    score_t _padding_score;

    template <typename sequence_bulk1_t,
              typename sequence_bulk2_t,
              typename dp_column_t,
              typename dp_row_t>
    auto operator()(sequence_bulk1_t && sequence_bulk1,
                    sequence_bulk2_t && sequence_bulk2,
                    dp_column_t dp_column,
                    dp_row_t dp_row,
                    cfg::end_gap _column_trailing_gaps = cfg::end_gap::penalised,
                    cfg::end_gap _row_trailing_gaps = cfg::end_gap::penalised) const noexcept
    {
        score_t result_score = score_bulk(sequence_bulk1,
                                          sequence_bulk2,
                                          dp_column,
                                          dp_row,
                                          _column_trailing_gaps,
                                          _row_trailing_gaps);

        using aligner_result_t = _aligner_result::value<sequence_bulk1_t, sequence_bulk2_t, dp_column_t, dp_row_t>;
        aligner_result_t base{std::forward<sequence_bulk1_t>(sequence_bulk1),
                              std::forward<sequence_bulk2_t>(sequence_bulk2),
                              std::move(dp_column),
                              std::move(dp_row),
                              _column_trailing_gaps,
                              _row_trailing_gaps};

        return _bulk_factory_saturated::value<aligner_result_t, score_t>{std::move(base), result_score};
    }

private:

    template <typename sequence_bulk1_t,
              typename sequence_bulk2_t,
              typename dp_column_t,
              typename dp_row_t>
    constexpr auto score_bulk(sequence_bulk1_t && sequence_bulk1,
                              sequence_bulk2_t && sequence_bulk2,
                              dp_column_t && dp_column,
                              dp_row_t && dp_row,
                              cfg::end_gap last_column,
                              cfg::end_gap last_row) const noexcept
    {
        if (last_column == cfg::end_gap::penalised && last_row == cfg::end_gap::penalised) {
            score_t best_score{};
            for (std::ptrdiff_t idx = 0; idx < std::ranges::distance(sequence_bulk1); ++idx) {
                best_score[idx] = select_score(idx, sequence_bulk1[idx], sequence_bulk2[idx], dp_column, dp_row);
            }
            return best_score;
        }

        score_t best_score{std::numeric_limits<typename score_t::value_type>::lowest()};
        if (last_column == cfg::end_gap::free) {
            best_score = find_max_score_bulk(sequence_bulk1, dp_column, sequence_bulk2, dp_row);
        }

        if (last_row == cfg::end_gap::free) {
            best_score = max(best_score, find_max_score_bulk(sequence_bulk2, dp_row, sequence_bulk1, dp_column));
        }
        return best_score;
    }

    template <typename sequence_t, typename dp_vector_t>
    constexpr auto get_offsets(sequence_t && sequence, dp_vector_t && dp_vector) const noexcept
    {
        size_t sequence_size = std::ranges::distance(sequence);

        size_t chunk_size = dp_vector[0].size() - 1;
        size_t const full_chunks = dp_vector.size() - 1;
        size_t dp_vector_size = (full_chunks * chunk_size) + dp_vector[full_chunks].size();

        assert(dp_vector_size > sequence_size);
        size_t offset = dp_vector_size - 1 - sequence_size;

        return std::tuple{sequence_size, dp_vector_size, offset};
    }

    template <typename sequence1_t,
              typename sequence2_t,
              typename dp_column_t,
              typename dp_row_t>
    constexpr auto select_score(size_t const simd_idx,
                                sequence1_t && sequence1,
                                sequence2_t && sequence2,
                                dp_column_t && dp_column,
                                dp_row_t && dp_row) const noexcept
        -> typename score_t::value_type
    {
        auto [column_sequence_size, column_vector_size, row_offset] = get_offsets(sequence1, dp_column);
        auto [row_sequence_size, row_vector_size, column_offset] = get_offsets(sequence2, dp_row);

        using scalar_t = typename score_t::value_type;
        size_t scale = std::min(column_offset, row_offset);
        scalar_t best_score{};

        size_t const chunk_size = dp_column[0].size() - 1;

        if (scale == column_offset) {
            auto [chunk_id, chunk_offset] =
                to_local_position(column_sequence_size + column_offset, chunk_size, dp_column.size());
            best_score = score_at(dp_column[chunk_id][chunk_offset], simd_idx);
        } else {
            auto [chunk_id, chunk_offset] =
                to_local_position(row_sequence_size + row_offset, chunk_size, dp_row.size());
            best_score = score_at(dp_row[chunk_id][chunk_offset], simd_idx);
        }

        return best_score - (_padding_score[simd_idx] * scale);
    }

    template <typename first_sequence_t, typename first_vector_t, typename second_sequence_t, typename second_vector_t>
    constexpr auto find_max_score_bulk(first_sequence_t && first_sequence,
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
        size_t first_vector_size{};
        size_t second_vector_size{};

        for (std::ptrdiff_t idx = 0; idx < std::ranges::distance(first_sequence); ++idx) {
            auto [first_sequence_size, first_size, second_offset] = get_offsets(first_sequence[idx], first_vector);
            auto [second_sequence_size, second_size, first_offset] = get_offsets(second_sequence[idx], second_vector);

            first_vector_size = first_size;
            second_vector_size = second_size;

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
        offset_simd_t simd_idx{0};
        for (size_t chunk_idx = 0; chunk_idx < first_vector.size(); ++chunk_idx) {
            size_t chunk_end = first_vector[chunk_idx].size() - (chunk_idx < first_vector.size() - 1);
            for (size_t local_idx = 0; local_idx < chunk_end; ++local_idx, ++simd_idx) {
                auto mask = (start_offset_first.le(simd_idx) && simd_idx.lt(end_offset_first));
                best_score = mask_max(best_score, mask, best_score, first_vector[chunk_idx][local_idx].score() - scale);
            }
        }

        // Note if second_offset is less than first_offset, the last (first_offset - second_offset) elements of the
        // second vector must be considered as well; these cells contain the values of the projected first vector, which
        // breaks around the cell (n, m) of the extended simd matrix.
        // This slice starts at the projected second vector cell.
        scale = _padding_score * scale_second;
        simd_idx = offset_simd_t{0};
        for (size_t chunk_idx = 0; chunk_idx < second_vector.size(); ++chunk_idx) {
            size_t chunk_end = second_vector[chunk_idx].size() - (chunk_idx < second_vector.size() - 1);
            for (size_t local_idx = 0; local_idx < chunk_end; ++local_idx, ++simd_idx) {
                auto mask = (start_offset_second.le(simd_idx) && simd_idx.lt(end_offset_second));
                best_score = mask_max(best_score, mask, best_score, second_vector[chunk_idx][local_idx].score() -
                                        scale);
                scale = mask_add(scale, mask, scale, _padding_score);
            }
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

    constexpr std::pair<size_t, size_t> to_local_position(size_t const position, size_t const chunk_size, size_t const chunk_count)
        const noexcept
    {
        size_t idx = position / chunk_size;
        size_t offset = position % chunk_size;
        // The last block has one more cell to account for the initialisation cell.
        // If the last given position of the currently processed vector is a multiple of the chunk size, then
        // it is the stored inside the last chunk. The respective index and offset are corrected accordingly.
        size_t is_last_position = (idx == chunk_count);
        return {idx - is_last_position, offset + (chunk_size * is_last_position)};
    }

};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
