// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_algorithm_template_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <iostream>
#include <seqan3/std/span>

#include <seqan3/utility/views/slice.hpp>

#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_standard.hpp>
#include <pairwise_aligner/matrix/dp_vector_saturated.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace detail
{
template <typename dp_vector_t> // saturated vector
class saturated_wrapper
{
public:

    using range_type = std::remove_cvref_t<decltype(std::declval<dp_vector_t>().base())>;
    using value_type = typename range_type::value_type;
    using reference = typename range_type::reference;
    using const_reference = typename range_type::const_reference;

private:
    using score_t = typename value_type::score_type; //int8_t

    dp_vector_t & _dp_vector; // int8_t
    // score_t _offset{}; // int8_t

public:

    saturated_wrapper() = delete;
    explicit saturated_wrapper(dp_vector_t & dp_vector) : _dp_vector{dp_vector}
    {}

    reference operator[](size_t const pos) noexcept
    {
        return range()[pos];
    }

    const_reference operator[](size_t const pos) const noexcept
    {
        return range()[pos];
    }

    constexpr size_t size() const noexcept
    {
        return _dp_vector.size();
    }

    auto & range() noexcept
    {
        return _dp_vector.base();
    }

    auto const & range() const noexcept
    {
        return _dp_vector.base();
    }

    dp_vector_t & base() noexcept
    {
        return _dp_vector;
    }

    dp_vector_t const & base() const noexcept
    {
        return _dp_vector;
    }

    // constexpr score_t offset() const noexcept
    // {
    //     return _offset;
    // }

    // template <typename offset_score_t>
    constexpr void offset(score_t new_offset) noexcept
    {
        initialise(new_offset);
        // upcast the local offset to the large offset value
        // add the converted offset to the global offset of the current vector type.
        using large_offset_t = std::remove_cvref_t<decltype(_dp_vector.offset())>;
        _dp_vector.offset(_dp_vector.offset() + large_offset_t{new_offset}); // int32_t
    }

    // template <typename offset_score_t>
    void initialise(score_t const & new_offset) noexcept
    {
        for (size_t i = 0; i < size(); ++i)
            std::apply([&] (auto & ...values) { ((values -= new_offset), ...); }, range()[i]);
    }
};

}

template <typename derived_t>
struct _dp_algorithm_template_saturated
{
    class type;
};

template <typename derived_t>
using dp_algorithm_template_saturated = typename _dp_algorithm_template_saturated<derived_t>::type;

template <typename derived_t>
class _dp_algorithm_template_saturated<derived_t>::type
    // : dp_algorithm_template_standard<_dp_algorithm_template_saturated<derived_t>::type>
{
    // using base_t = dp_algorithm_template_standard<_dp_algorithm_template_saturated<derived_t>::type>;
protected:

    // friend base_t;
    friend derived_t;

    template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
    auto run(sequence1_t && sequence1, sequence2_t && sequence2, dp_column_t dp_column, dp_row_t dp_row) const
    {
        // ----------------------------------------------------------------------------
        // Initialisation
        // ----------------------------------------------------------------------------

        // auto start = std::chrono::high_resolution_clock::now();

        // step 1) initialise original vectors
        auto simd_seq1 = as_derived().initialise_column_vector(sequence1, dp_column);
        auto simd_seq2 = as_derived().initialise_row_vector(sequence2, dp_row);

        // auto end = std::chrono::high_resolution_clock::now();
        // std::cout << "Initialise: " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - start).count() << " [ns]" << std::endl;

        // 1) transformed sequence is a vector of chunks.
        // 2) the range is a vector of ranges (but should return proxies)
        // 3) we don't want to work with proxy but only with the pure types.
        // 4) we want to modify the offset again.

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------

        // auto & dp_column_group = dp_column.range();
        // auto & dp_row_group = dp_row.range();

        size_t const col_count = dp_column.size();
        size_t const row_count = dp_row.size();
        size_t const col_chunk_size = dp_column[0].size() - 1;
        size_t const row_chunk_size = dp_row[0].size() - 1;

        using dp_column_chunk_t = std::ranges::range_value_t<decltype(dp_column.range())>;
        using dp_row_chunk_t = std::ranges::range_value_t<decltype(dp_row.range())>;

        using value_t = std::ranges::range_value_t<decltype(simd_seq1)>;
        std::vector<std::span<value_t>> seq1_chunked{};
        seq1_chunked.reserve(col_count);

        using saturated_col_t = detail::saturated_wrapper<dp_column_chunk_t>;
        std::vector<saturated_col_t> dp_column_chunks{};
        dp_column_chunks.reserve(col_count);

        for (size_t i = 0; i < col_count; ++i)
        {
            seq1_chunked.emplace_back(std::ranges::next(std::ranges::begin(simd_seq1), (i * col_chunk_size)),
                                      std::ranges::next(std::ranges::begin(simd_seq1), ((i + 1) * col_chunk_size), std::ranges::end(simd_seq1)));
            dp_column_chunks.emplace_back(saturated_col_t{dp_column.range()[i]});
        }

        // double compute{};
        // double offset{};
        // start = std::chrono::high_resolution_clock::now();
        for (size_t j = 0; j < col_count; ++j) {
            detail::saturated_wrapper<dp_row_chunk_t> current_row_vector{dp_row.range()[j]};
            std::span transformed_seq2{std::ranges::next(std::ranges::begin(simd_seq2), (j * row_chunk_size)),
                                       std::ranges::next(std::ranges::begin(simd_seq2), ((j + 1) * row_chunk_size), std::ranges::end(simd_seq2))};
            for (size_t i = 0; i < row_count; ++i) {
                // start = std::chrono::high_resolution_clock::now();
                // std::span transformed_seq1{std::ranges::next(std::ranges::begin(simd_seq1), (i * col_chunk_size)),
                //                            std::ranges::next(std::ranges::begin(simd_seq1), ((i + 1) * col_chunk_size), std::ranges::end(simd_seq1))};
                detail::saturated_wrapper<dp_column_chunk_t> current_column_vector{dp_column.range()[i]};
                current_row_vector.offset(current_row_vector[0].score());
                dp_column_chunks[i].offset(dp_column_chunks[i][0].score());
                // auto & col = dp_column_group[i].base();
                // auto & row = dp_row_group[j].base();
                // end = std::chrono::high_resolution_clock::now();
                // offset += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
                // start = std::chrono::high_resolution_clock::now();
                _run(seq1_chunked[i], transformed_seq2, dp_column_chunks[i], current_row_vector);
                // end = std::chrono::high_resolution_clock::now();
                // compute += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

                // std::cout << "col:";
                // for (size_t k = 0; k < current_column_vector.size(); ++k)
                //     std::cout << " " << (int32_t) current_column_vector[k].score()[0];
                // std::cout << "\n";

                // std::cout << "row:";
                // for (size_t k = 0; k < current_row_vector.size(); ++k)
                //     std::cout << " " << (int32_t) current_row_vector[k].score()[0];
                // std::cout << "\n";

            }
        }
        // end = std::chrono::high_resolution_clock::now();
        // std::cout << "Compute: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " [ns]" << std::endl;
        // std::cout << "Offset: " << offset << " [ns]" << std::endl;
        // std::cout << "Compute: " << compute << " [ns]" << std::endl;

        // ----------------------------------------------------------------------------
        // Create result
        // ----------------------------------------------------------------------------

        // auto best_score = dp_column.range().back().range().back().score();
        using score_t = std::remove_cvref_t<decltype(dp_column.range().back().range().back().score())>;
        return as_derived().make_result(std::forward<sequence1_t>(sequence1),
                                        std::forward<sequence2_t>(sequence2),
                                        std::move(dp_column),
                                        std::move(dp_row),
                                        score_t{});
    }

    template <typename sequence1_t,
              typename sequence2_t,
              typename dp_column_vector_t,
              typename dp_row_vector_t>
    void _run(sequence1_t && sequence1,
             sequence2_t && sequence2,
             dp_column_vector_t && dp_column_vector,
             dp_row_vector_t && dp_row_vector) const
    {

        for (size_t i = 1; i < dp_column_vector.size(); ++i) {
            as_derived().compute_first_column(dp_row_vector[0], dp_column_vector[i]);
        }

        for (size_t j = 0; j < std::ranges::size(sequence2); ++j)
        {
            auto cache = as_derived().initialise_column(dp_row_vector[j+1], dp_column_vector[0]);
            size_t i = 0;
            for (; i < std::ranges::size(sequence1); ++i) {
                as_derived().compute_cell(cache, dp_column_vector[i+1], sequence1[i], sequence2[j]);
            }
            as_derived().finalise_column(dp_row_vector[j+1], dp_column_vector[i], cache);
        }
    }

    // template <typename sequence_t, typename dp_vector_t>
    // sequence_t initialise_row_vector(sequence_t && sequence, [[maybe_unused]] dp_vector_t & dp_vector) const
    // {
    //     // TODO: do what?
    //     // dp_vector.initialise();
    //     return std::forward<sequence_t>(sequence);
    // }

    // template <typename sequence_t, typename dp_vector_t>
    // sequence_t initialise_column_vector(sequence_t && sequence, [[maybe_unused]] dp_vector_t & dp_vector) const
    // {
    //     // dp_vector.initialise();
    //     return std::forward<sequence_t>(sequence);
    // }

    // template <typename ...args_t>
    // auto compute_first_column(args_t && ...args) const noexcept
    // {
    //     return as_derived().compute_first_column(std::forward<args_t>(args)...);
    // }

    // template <typename ...args_t>
    // auto initialise_column(args_t && ...args) const noexcept
    // {
    //     return as_derived().initialise_column(std::forward<args_t>(args)...);
    // }

    // template <typename ...args_t>
    // void finalise_column(args_t && ...args) const noexcept
    // {
    //     as_derived().finalise_column(std::forward<args_t>(args)...);
    // }

    // template <typename ...args_t>
    // constexpr void make_result(args_t && .../*args*/) const noexcept
    // {
    //     // noop
    // }

    // template <typename ...args_t>
    // void compute_cell(args_t && ...args) const noexcept
    // {
    //     as_derived().compute_cell(std::forward<args_t>(args)...);
    // }

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
