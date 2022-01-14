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

        // std::cout << "index = " << idx << "\n";

        scalar_t best_score = std::numeric_limits<scalar_t>::min();
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

        if (base_value_t::_row_trailing_gaps == cfg::end_gap::free)
        {
            for (size_t cell_idx = 0; cell_idx < this->dp_row().size(); ++cell_idx)
                best_score = std::max<scalar_t>(score_at(this->dp_row()[cell_idx], idx), best_score);
        }

        if (base_value_t::_column_trailing_gaps == cfg::end_gap::free)
        {
            for (size_t cell_idx = 0; cell_idx < this->dp_column().size(); ++cell_idx)
                best_score = std::max<scalar_t>(score_at(this->dp_column()[cell_idx], idx), best_score);
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
        return std::tuple{original_row_dim + offset, original_column_dim + offset, offset};
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
