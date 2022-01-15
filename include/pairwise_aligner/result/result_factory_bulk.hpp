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

    constexpr auto score_at(size_t const idx) const noexcept
    {
        if (base_value_t::_row_trailing_gaps == cfg::end_gap::penalised &&
            base_value_t::_column_trailing_gaps == cfg::end_gap::penalised) {
            return select_score(idx);
        }

        // Maximal score is either in last row or in last column or in both.
        using scalar_t = typename score_t::value_type;
        scalar_t best_score =  std::numeric_limits<scalar_t>::lowest();
        if (base_value_t::_column_trailing_gaps == cfg::end_gap::free) {
            best_score = find_max_score(this->sequence1()[idx], this->dp_column(),
                                        this->sequence2()[idx], this->dp_row(),
                                        idx);
        }

        if (base_value_t::_row_trailing_gaps == cfg::end_gap::free) {
            best_score = std::max(best_score, find_max_score(this->sequence2()[idx], this->dp_row(),
                                                             this->sequence1()[idx], this->dp_column(),
                                                              idx));
        }

        return best_score;
    }

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

    constexpr auto select_score(size_t const simd_idx) const noexcept
        -> typename score_t::value_type
    {
        auto [column_sequence_size, column_vector_size, row_offset] =
            get_offsets(this->sequence1()[simd_idx], this->dp_column());
        auto [row_sequence_size, row_vector_size, column_offset] =
            get_offsets(this->sequence2()[simd_idx], this->dp_row());

        using scalar_t = typename score_t::value_type;
        size_t scale = std::min(column_offset, row_offset);
        scalar_t best_score{};

        if (scale == column_offset) {
            best_score = score_at(this->dp_column()[column_sequence_size + scale], simd_idx);
        } else {
            best_score = score_at(this->dp_row()[row_sequence_size + row_offset], simd_idx);
        }

        return best_score - (_padding_score[simd_idx] * scale);
    }

    template <typename first_sequence_t, typename first_vector_t, typename second_sequence_t, typename second_vector_t>
    constexpr auto find_max_score(first_sequence_t && first_sequence,
                                  first_vector_t && first_vector,
                                  second_sequence_t && second_sequence,
                                  second_vector_t && second_vector,
                                  size_t const simd_idx) const noexcept
        -> typename score_t::value_type
    {
        using scalar_t = typename score_t::value_type;

        auto [first_sequence_size, first_vector_size, second_offset] = get_offsets(first_sequence, first_vector);
        auto [second_sequence_size, second_vector_size, first_offset] = get_offsets(second_sequence, second_vector);

        // Iterate over the corresponding slice in the projected dp vector (first_vector) and subtract the scaled
        // padding score from the retrieved values of the corresponding simd index.
        scalar_t best_score = std::numeric_limits<scalar_t>::lowest();
        size_t first_slice_end = std::min<size_t>(first_offset + first_sequence_size + 1, first_vector_size);
        size_t scale = first_offset;
        for (size_t idx = first_offset; idx < first_slice_end; ++idx) {
            best_score = std::max<scalar_t>(score_at(first_vector[idx], simd_idx) -
                                                static_cast<scalar_t>(_padding_score[simd_idx] * scale),
                                            best_score);
        }

        // Note if second_offset is less than first_offset, the last (first_offset - second_offset) elements of the
        // second vector must be considered as well; these cells contain the values of the projected first vector, which
        // breaks around the cell (n, m) of the extended simd matrix.
        // This slice starts at the projected second vector cell.
        scale = second_offset; // Reset scaling factor to second_offset.
        for (size_t idx = second_sequence_size + second_offset; idx < second_vector_size; ++idx, ++scale) {
            best_score = std::max<scalar_t>(score_at(second_vector[idx], simd_idx) -
                                                static_cast<scalar_t>(_padding_score[simd_idx] * scale),
                                            best_score);
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
        using aligner_result_t = _aligner_result::value<sequence_bulk1_t, sequence_bulk2_t, dp_column_t, dp_row_t>;
        aligner_result_t base{std::forward<sequence_bulk1_t>(sequence_bulk1),
                              std::forward<sequence_bulk2_t>(sequence_bulk2),
                              std::move(dp_column),
                              std::move(dp_row),
                              _column_trailing_gaps,
                              _row_trailing_gaps};

        return _bulk_factory::value<aligner_result_t, score_t>{std::move(base), _padding_score};
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
