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
        using scalar_t = typename score_t::value_type;

        scalar_t best_score = std::numeric_limits<scalar_t>::lowest();
        if (base_value_t::_row_trailing_gaps == cfg::end_gap::penalised &&
            base_value_t::_column_trailing_gaps == cfg::end_gap::penalised)
        {
            auto && [row_idx, col_idx, offset] = projected_coordinate(idx);

            // std::cout << "row_idx = " << row_idx << "\n";
            // std::cout << "col_idx = " << col_idx << "\n";
            // std::cout << "offset = " << offset << "\n";

            assert(row_idx == this->dp_column().size() - 1 || col_idx == this->dp_row().size() - 1);

            if (row_idx == this->dp_column().size() - 1) {
                best_score = score_at(this->dp_row()[col_idx], idx);
                // std::cout << "this->dp_row()[col_idx].score()[idx] = " << this->dp_row()[col_idx].score()[idx] << "\n";
            } else {
                best_score = score_at(this->dp_column()[row_idx], idx);
                // std::cout << "this->dp_column()[row_idx].score()[idx] = " << this->dp_column()[row_idx].score()[idx] << "\n";
            }

            // std::cout << "best_score = " << best_score << "\n";
            return static_cast<scalar_t>(best_score - static_cast<scalar_t>(_padding_score[idx] * offset));
        }

        if (base_value_t::_column_trailing_gaps == cfg::end_gap::free) {
            best_score = find_max_score(this->sequence1()[idx], this->dp_column(),
                                        this->sequence2()[idx], this->dp_row(),
                                        idx);
        }

        if (base_value_t::_row_trailing_gaps == cfg::end_gap::free) {
            best_score = find_max_score(this->sequence2()[idx], this->dp_row(),
                                        this->sequence1()[idx], this->dp_column(),
                                        idx);
        }

        return best_score;
    }

    constexpr std::tuple<size_t, size_t, size_t> projected_coordinate(size_t const idx) const noexcept
    {
        size_t const original_row_dim = std::ranges::distance(this->sequence1()[idx]);
        size_t const original_column_dim = std::ranges::distance(this->sequence2()[idx]);

        // std::cout << "original_row_dim = " << original_row_dim << "\n";
        // std::cout << "original_column_dim = " << original_column_dim << "\n";

        size_t offset = std::min<size_t>(this->dp_column().size() - 1 - original_row_dim,
                                         this->dp_row().size() -1 - original_column_dim);

        // std::cout << "this->dp_column().size() = " << this->dp_column().size() << "\n";
        // std::cout << "this->dp_row().size() = " << this->dp_row().size() << "\n";
    }

    template <typename first_sequence_t, typename first_vector_t, typename second_sequence_t, typename second_vector_t>
    constexpr auto find_max_score(first_sequence_t && first_sequence,
                                  first_vector_t && first_vector,
                                  second_sequence_t && second_sequence,
                                  second_vector_t && second_vector,
                                  size_t const simd_idx) const noexcept
    {
        using scalar_t = typename score_t::value_type;

        size_t const first_sequence_size = std::ranges::distance(first_sequence);
        size_t const second_sequence_size = std::ranges::distance(second_sequence);
        size_t const first_vector_size = first_vector.size();
        size_t const second_vector_size = second_vector.size();

        size_t const first_offset = second_vector_size - 1 - second_sequence_size;
        size_t const second_offset = first_vector_size - 1 - first_sequence_size;

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

private:

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
