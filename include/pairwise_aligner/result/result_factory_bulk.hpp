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

#include <pairwise_aligner/dp_trailing_gaps.hpp>
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
    dp_trailing_gaps _column_trailing_gaps{};
    dp_trailing_gaps _row_trailing_gaps{};

    constexpr auto score_at(size_t const idx) const noexcept
    {
        using scalar_t = typename score_t::value_type;

        assert(this->sequence1()[idx].has_value());
        assert(this->sequence2()[idx].has_value());

        scalar_t best_score = std::numeric_limits<scalar_t>::min();
        if (_row_trailing_gaps == dp_trailing_gaps::regular && _column_trailing_gaps == dp_trailing_gaps::regular)
        {
            auto && [row_idx, col_idx, offset] = projected_coordinate(idx);
            assert(row_idx == this->dp_column().size() - 1 || col_idx == this->dp_row().size() - 1);

            if (row_idx == this->dp_column().size() - 1)
                best_score = this->dp_row()[col_idx].score()[idx];
            else
                best_score = this->dp_column()[row_idx].score()[idx];
            return static_cast<scalar_t>(best_score - static_cast<scalar_t>(_padding_score[idx] * offset));
        }

        if (_row_trailing_gaps == dp_trailing_gaps::free)
        {
            for (size_t cell_idx = 0; cell_idx < this->dp_row().size(); ++cell_idx)
                best_score = std::max<scalar_t>(this->dp_row()[cell_idx].score()[idx], best_score);
        }

        if (_column_trailing_gaps == dp_trailing_gaps::free)
        {
            for (size_t cell_idx = 0; cell_idx < this->dp_column().size(); ++cell_idx)
                best_score = std::max<scalar_t>(this->dp_column()[cell_idx].score()[idx], best_score);
        }

        return best_score;
    }

    constexpr std::tuple<size_t, size_t, size_t> projected_coordinate(size_t const idx) const noexcept
    {
        size_t const original_row_dim = std::ranges::distance(*(this->sequence1()[idx]));
        size_t const original_column_dim = std::ranges::distance(*(this->sequence2()[idx]));

        size_t offset = std::min<size_t>(this->dp_column().size() - 1 - original_row_dim,
                                         this->dp_row().size() -1 - original_column_dim);
        return std::tuple{original_row_dim + offset, original_column_dim + offset, offset};
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
    dp_trailing_gaps _column_trailing_gaps;
    dp_trailing_gaps _row_trailing_gaps;

    template <typename sequence1_t, size_t bulk1_size,
              typename sequence2_t, size_t bulk2_size,
              typename dp_column_t,
              typename dp_row_t>
    auto operator()(std::array<sequence1_t, bulk1_size> sequence_bulk1,
                    std::array<sequence2_t, bulk2_size> sequence_bulk2,
                    dp_column_t dp_column,
                    dp_row_t dp_row,
                    score_t score) const noexcept
    {
        static_assert(bulk1_size == bulk2_size, "The sequence bulks must have the same length.");

        using aligner_result_t = _aligner_result::value<std::array<sequence1_t, bulk1_size>,
                                                        std::array<sequence2_t, bulk2_size>,
                                                        dp_column_t,
                                                        dp_row_t,
                                                        score_t>;
        aligner_result_t base{std::move(sequence_bulk1),
                              std::move(sequence_bulk2),
                              std::move(dp_column),
                              std::move(dp_row),
                              std::move(score)};

        return _bulk_factory::value<aligner_result_t, score_t>{std::move(base),
                                                               _padding_score,
                                                               _column_trailing_gaps,
                                                               _row_trailing_gaps};
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
