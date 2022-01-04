// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_algorithm_template_standard.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <chrono>
#include <seqan3/std/ranges>

#include <pairwise_aligner/result/aligner_result.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename derived_t>
struct _dp_algorithm_template_standard
{
    class type;
};

template <typename derived_t>
using dp_algorithm_template_standard = typename _dp_algorithm_template_standard<derived_t>::type;

template <typename derived_t>
class _dp_algorithm_template_standard<derived_t>::type
{
protected:

    friend derived_t;

    template <typename sequence1_t,
              typename sequence2_t,
              typename dp_column_vector_t,
              typename dp_row_vector_t>
    auto run(sequence1_t && sequence1,
             sequence2_t && sequence2,
             dp_column_vector_t dp_column_vector,
             dp_row_vector_t dp_row_vector) const
    {
        // ----------------------------------------------------------------------------
        // Initialisation
        // ----------------------------------------------------------------------------

        // auto start = std::chrono::high_resolution_clock::now();
        auto transformed_seq1 = as_derived().initialise_column_vector(sequence1, dp_column_vector);
        auto transformed_seq2 = as_derived().initialise_row_vector(sequence2, dp_row_vector);

        // auto end = std::chrono::high_resolution_clock::now();
        // std::cout << "\tInitialise: " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count() << " [ns]" << std::endl;

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------
        // using base_t = std::remove_cvref_t<decltype(dp_column_vector.base())>;
        // using cell_t = typename base_t::value_type;
        // using score_t = typename cell_t::score_type;

        // using cell2_t = typename dp_column_vector_t::value_type;
        // using score2_t = typename cell2_t::score_type;

        // std::cout << "dp_row_vector.size() = " << dp_row_vector.size() << "\n";
        // std::cout << "dp_column_vector.size() = " << dp_column_vector.size() << "\n";
        // std::cout << "dp_row_vector.base().offset()[0] = " << dp_row_vector.base().offset()[0] << "\n";
        // std::cout << "dp_column_vector.base().offset()[0] = " << dp_column_vector.base().offset()[0] << "\n";
        // std::cout << "dp_row_vector.offset()[0] = " << (int) dp_row_vector.offset()[0] << "\n";
        // std::cout << "dp_column_vector.offset()[0] = " << (int) dp_column_vector.offset()[0] << "\n";

        // std::vector<std::vector<score_t>> debug_matrix{};
        // debug_matrix.resize(dp_row_vector.size());
        // for (auto & col : debug_matrix)
        //     col.resize(dp_column_vector.size());

        // std::vector<std::vector<score2_t>> debug_matrix2{};
        // debug_matrix2.resize(dp_row_vector.size());
        // for (auto & col : debug_matrix2)
        //     col.resize(dp_column_vector.size());

        // size_t k = 0;

        // k = 0;
        // for (auto & init_col : debug_matrix) {
        //    init_col[0] = dp_row_vector.base()[k++].score();
        // }

        // k = 0;
        // for (auto & init_row : debug_matrix[0]) {
        //    init_row = dp_column_vector.base()[k++].score();
        // }

        // k = 0;
        // for (auto & init_row : debug_matrix2[0]) {
        //    init_row = dp_column_vector[k++].score();
        // }

        // k = 0;
        // for (auto & init_col : debug_matrix2) {
        //    init_col[0] = dp_row_vector[k++].score();
        // }
        // what is the first value?
        // we don't want to do something in the first row.
        // it has been initialised already.
        // can we compute it?
        // using time_t = std::remove_reference_t<decltype(std::chrono::high_resolution_clock::now())>;
        // start = std::chrono::high_resolution_clock::now();
        // Store the best score of the last cell of the column vector in the first cell of the row vector.
        as_derived().compute_first_column(dp_row_vector[0], dp_column_vector[dp_column_vector.size() - 1]);

        for (size_t j = 0; j < std::ranges::size(transformed_seq2); ++j)
        {
            auto cache = as_derived().initialise_column(dp_row_vector[j+1], dp_column_vector[0]);
            size_t i = 0;
            for (; i < std::ranges::size(transformed_seq1); ++i) {
                as_derived().compute_cell(cache, dp_column_vector[i+1], transformed_seq1[i], transformed_seq2[j]);
                // debug_matrix[j+1][i+1] = dp_column_vector.base()[i+1].score();
                // debug_matrix2[j+1][i+1] = dp_column_vector[i+1].score();
            }
            // total += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start).count();

            as_derived().finalise_column(dp_row_vector[j+1], dp_column_vector[i], cache);
        }
        // end = std::chrono::high_resolution_clock::now();
        // std::cout << "\tCompute: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " [ns]" << std::endl;
        // std::cout << "Debug matrix:\n";
        // for (size_t i = 0; i < std::ranges::size(debug_matrix[0]); ++i) {
        //     for (size_t j = 0; j < std::ranges::size(debug_matrix); ++j) {
        //         std::cout << (int) debug_matrix[j][i][0] << "\t";
        //     }
        //     std::cout << "\n";
        // }
        // std::cout << "\n";

        // std::cout << "Debug matrix2:\n";
        // for (size_t i = 0; i < std::ranges::size(debug_matrix2[0]); ++i) {
        //     for (size_t j = 0; j < std::ranges::size(debug_matrix2); ++j) {
        //         std::cout << (int) debug_matrix2[j][i][0] << "\t";
        //     }
        //     std::cout << "\n";
        // }
        // std::cout << "\n";

        // ----------------------------------------------------------------------------
        // Create result
        // ----------------------------------------------------------------------------

        auto best_score = get<0>(dp_column_vector[std::ranges::size(transformed_seq1)]);
        return as_derived().make_result(std::forward<sequence1_t>(sequence1),
                                        std::forward<sequence2_t>(sequence2),
                                        std::move(dp_column_vector),
                                        std::move(dp_row_vector),
                                        std::move(best_score));
    }

    derived_t const & as_derived() const noexcept
    {
        return static_cast<derived_t const &>(*this);
    }

    derived_t & as_derived() noexcept
    {
        return static_cast<derived_t &>(*this);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
