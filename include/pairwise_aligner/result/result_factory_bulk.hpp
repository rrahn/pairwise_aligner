// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::aligner_result.
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
namespace _bulk_factory
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

} // namespace _bulk_factory

template <typename score_t>
struct _result_factory_bulk
{
    struct type;
};

template <typename score_t>
using result_factory_bulk = typename _result_factory_bulk<score_t>::type;

template <typename score_t>
struct _result_factory_bulk<score_t>::type
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

        return _bulk_factory::value<aligner_result_t, score_t>{std::move(base), result_score};
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
        size_t dp_vector_size = dp_vector.size();

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

        if (scale == column_offset) {
            best_score = score_at(dp_column[column_sequence_size + column_offset], simd_idx);
        } else {
            best_score = score_at(dp_row[row_sequence_size + row_offset], simd_idx);
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

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
