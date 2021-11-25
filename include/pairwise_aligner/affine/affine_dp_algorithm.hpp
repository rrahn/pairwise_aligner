// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::affine_dp_algorithm.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <functional>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/affine/affine_cell.hpp>
#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_standard.hpp>
#include <pairwise_aligner/simd_score_type.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <template <typename> typename dp_template, typename gap_model_t, typename init_strategy_t>
class affine_dp_algorithm : protected dp_template<affine_dp_algorithm<dp_template, gap_model_t, init_strategy_t>>
{
private:

    using base_t = dp_template<affine_dp_algorithm<dp_template, gap_model_t, init_strategy_t>>;

    friend base_t;

    using score_t = int16_t;
    using simd_score_t = simd_score<score_t>;

    gap_model_t _gap_model{};
    init_strategy_t _initialisaion_strategy{};

    alignas(32) score_t match_score = 4;
    alignas(32) score_t mismatch_score = -5;

    simd_score_t gap_extension_score_simd{-1};
    simd_score_t gap_open_score_simd{-10};

    simd_score_t match_score_simd{match_score};
    simd_score_t mismatch_score_simd{mismatch_score};

public:
    affine_dp_algorithm() = default;
    explicit affine_dp_algorithm(gap_model_t gap_model, init_strategy_t initialisaion_strategy) noexcept :
        _gap_model{std::move(gap_model)},
        _initialisaion_strategy{std::move(initialisaion_strategy)}
    {}

protected:

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_row_vector(sequence_t && sequence, dp_vector_t & dp_vector)
    {
        return dp_vector.initialise(std::forward<sequence_t>(sequence), _initialisaion_strategy);
    }

    template <std::ranges::viewable_range sequence_t, typename dp_vector_t>
        requires std::ranges::forward_range<sequence_t>
    auto initialise_column_vector(sequence_t && sequence, dp_vector_t & dp_vector) const
    {
        return dp_vector.initialise(std::forward<sequence_t>(sequence), _initialisaion_strategy);
    }

    template <typename row_cell_t, typename column_cell_t>
    auto initialise_column(row_cell_t & current_row_cell, column_cell_t & first_column_cell) const noexcept
    {
        using std::max;
        using score_t = typename row_cell_t::score_type;

        std::pair cache{get<0>(first_column_cell), get<1>(current_row_cell)};
        get<0>(first_column_cell) = get<0>(current_row_cell);

        if constexpr (!std::integral<score_t>) {
            get<1>(first_column_cell) = max(static_cast<score_t>(cache.first + gap_open_score_simd),
                                            static_cast<score_t>(get<1>(first_column_cell) + gap_extension_score_simd));
        } else {
            get<1>(first_column_cell) = max(static_cast<score_t>(cache.first + _gap_model.gap_open_score),
                                            static_cast<score_t>(get<1>(first_column_cell) + _gap_model.gap_extension_score));
        }
        return cache;
    }

    template <typename row_cell_t, typename column_cell_t, typename cache_t>
    void finalise_column(row_cell_t & current_row_cell,
                         column_cell_t const & last_column_cell,
                         cache_t & cache) const noexcept
    {
        get<0>(current_row_cell) = get<0>(last_column_cell);
        get<1>(current_row_cell) = std::move(cache.second);
    }

    template <typename cache_t, typename seq1_val_t, typename seq2_val_t, typename dp_cell_t>
    auto compute_cell(cache_t & cache,
                      dp_cell_t & column_cell,
                      seq1_val_t const & seq1_val,
                      seq2_val_t const & seq2_val) const noexcept
    {
        using std::max;
        using score_t = typename dp_cell_t::score_type;

        auto [next_diagonal, horizontal_score] = column_cell;
        // TODO: Should be a score model -> the code depends on the score model.
        if constexpr (std::integral<score_t>) {
            cache.first += (seq1_val == seq2_val) ? match_score : mismatch_score;
        } else {
            cache.first += compare_and_blend(seq1_val, seq2_val, match_score_simd, mismatch_score_simd);
        }
        cache.first = max(max(cache.first, cache.second), horizontal_score);
        get<0>(column_cell) = cache.first;
        cache.first += (_gap_model.gap_open_score + _gap_model.gap_extension_score);
        cache.second = max(static_cast<score_t>(cache.second + _gap_model.gap_extension_score), cache.first);
        get<1>(column_cell) = max(static_cast<score_t>(horizontal_score + _gap_model.gap_extension_score), cache.first);
        cache.first = next_diagonal; // cache score
    }
};

// ----------------------------------------------------------------------------
// Define the concrete parwise aligner instances.
// ----------------------------------------------------------------------------

template <typename ...policies_t>
class pairwise_aligner_affine : public affine_dp_algorithm<dp_algorithm_template_standard, std::remove_reference_t<policies_t>...>
{
private:
    using base_t = affine_dp_algorithm<dp_algorithm_template_standard, std::remove_reference_t<policies_t>...>;
public:

    explicit pairwise_aligner_affine(policies_t && ...policies) noexcept :
        base_t{std::forward<policies_t>(policies)...}
    {}
};

template <typename ...policies_t>
pairwise_aligner_affine(policies_t && ...) -> pairwise_aligner_affine<policies_t...>;

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
