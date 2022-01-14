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

#include <array>
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
        _dp_vector.offset() += large_offset_t{new_offset}; // int32_t
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

        // first dp row is initialised correctly
        // first dp col is initialised correctly

    /*
        TODO: make_result without score!
        return make_result(...);
    */

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------

        size_t const col_count = dp_column.size();
        size_t const row_count = dp_row.size();
        size_t const col_chunk_size = dp_column[0].size() - 1;
        size_t const row_chunk_size = dp_row[0].size() - 1;

        using dp_column_chunk_t = typename dp_column_t::value_type; // std::ranges::range_value_t<decltype(dp_column.range())>;
        using dp_row_chunk_t = typename dp_row_t::value_type; // std::ranges::range_value_t<decltype(dp_row.range())>;

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
            dp_column_chunks.emplace_back(dp_column[i]);
        }

        for (size_t j = 0; j < row_count; ++j) {
            detail::saturated_wrapper<dp_row_chunk_t> current_row_vector{dp_row[j]};
            std::span transformed_seq2{std::ranges::next(std::ranges::begin(simd_seq2), (j * row_chunk_size)),
                                       std::ranges::next(std::ranges::begin(simd_seq2), ((j + 1) * row_chunk_size), std::ranges::end(simd_seq2))};

            // Initialise first block of current column.
            current_row_vector.offset(current_row_vector[0].score());
            dp_column_chunks[0].offset(dp_column_chunks[0][0].score());

            initialise_block(dp_column_chunks[0], current_row_vector);

            // Iterate over blocks in current column.
            for (size_t i = 0; i < col_count; ++i) {
                if (i > 0) {
                    current_row_vector[0].score() = dp_column_chunks[i - 1][dp_column_chunks[i - 1].size() - 1].score();
                    current_row_vector.offset(current_row_vector[1].score());
                    dp_column_chunks[i].offset(dp_column_chunks[i][0].score());
                    dp_column_chunks[i][0].score() = current_row_vector[0].score();
                }

                _run(seq1_chunked[i], transformed_seq2, dp_column_chunks[i], current_row_vector);

                // std::cout << "col:";
                // for (size_t k = 0; k < current_column_vector.size(); ++k)
                //     std::cout << " " << (int32_t) current_column_vector[k].score()[0];
                // std::cout << "\n";

                // std::cout << "row:";
                // for (size_t k = 0; k < current_row_vector.size(); ++k)
                //     std::cout << " " << (int32_t) current_row_vector[k].score()[0];
                // std::cout << "\n";

            }

            // Write back optimal score to row vector.
            postprocess_block(dp_column_chunks.back(), current_row_vector);
        }

        // ----------------------------------------------------------------------------
        // Create result
        // ----------------------------------------------------------------------------

        return as_derived().make_result(std::forward<sequence1_t>(sequence1),
                                        std::forward<sequence2_t>(sequence2),
                                        std::move(dp_column),
                                        std::move(dp_row));
    }

    template <typename sequence1_t,
              typename sequence2_t,
              typename dp_column_vector_t,
              typename dp_row_vector_t>
    void _run(sequence1_t sequence1,
              sequence2_t sequence2,
              dp_column_vector_t && dp_column_vector,
              dp_row_vector_t && dp_row_vector) const
    {
        std::ptrdiff_t const seq1_size = std::ranges::ssize(sequence1);
        std::ptrdiff_t const seq2_size = std::ranges::ssize(sequence2);

        using value_t = typename std::remove_cvref_t<dp_row_vector_t>::value_type;

        // Initialise bulk_cache array.
        constexpr std::ptrdiff_t cache_size = 8;
        std::array<value_t, cache_size> bulk_cache{};

        std::ptrdiff_t j = 0;
        for (; j < seq2_size - (cache_size - 1); j += cache_size)
        {
            // copy values into cache.
            std::span seq2_slice{sequence2.begin() + j, sequence2.begin() + j + cache_size};
            unroll_load(bulk_cache, dp_row_vector, j + 1, std::make_index_sequence<cache_size>());

            // compute cache many cells in one row for one horziontal value.
            for (std::ptrdiff_t i = 0; i < seq1_size; ++i) {
                auto cacheH = dp_column_vector[i+1];

                unroll_loop(bulk_cache, cacheH, sequence1[i], seq2_slice, std::make_index_sequence<cache_size>());

                dp_column_vector[i+1] = cacheH;
            }

            // store results of cache back in vector.
            unroll_store(dp_row_vector, bulk_cache, j + 1, std::make_index_sequence<cache_size>());
        }

        // Compute remaining cells.
        for (; j < seq2_size; ++j) {
            auto cache = dp_row_vector[j + 1];

            for (std::ptrdiff_t i = 0; i < seq1_size; ++i) {
                as_derived().compute_cell(cache, dp_column_vector[i+1], sequence1[i], sequence2[j]);
            }

            dp_row_vector[j + 1] = cache;
        }
    }

    template <typename row_cells_t,
              typename col_cell_t,
              typename seq1_value_t,
              typename seq2_values_t,
              size_t ...idx>
    constexpr void unroll_loop(row_cells_t & row_cells,
                               col_cell_t & col_cell,
                               seq1_value_t const & seq1_value,
                               seq2_values_t const & seq2_values,
                               [[maybe_unused]] std::index_sequence<idx...> const & indices) const noexcept
    {
        (as_derived().compute_cell(row_cells[idx], col_cell, seq1_value, seq2_values[idx]), ...);
    }

    template <typename cache_t, typename row_vector_t, size_t ...idx>
    constexpr void unroll_load(cache_t & bulk_cache,
                               row_vector_t const & row_vector,
                               size_t const offset,
                               [[maybe_unused]] std::index_sequence<idx...> const & indices) const noexcept
    {
        ((bulk_cache[idx] = row_vector[offset + idx]), ...);
    }

    template <typename row_vector_t, typename cache_t, size_t ...idx>
    constexpr void unroll_store(row_vector_t & row_vector,
                                cache_t const & bulk_cache,
                                size_t const offset,
                                [[maybe_unused]] std::index_sequence<idx...> const & indices) const noexcept
    {
        ((row_vector[offset + idx] = bulk_cache[idx]), ...);
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

    template <typename column_vector_t, typename row_vector_t>
    constexpr void initialise_block(column_vector_t && column_vector, row_vector_t && row_vector) const noexcept
    {
        // H'[0].score() = V'[m].score();
        // V'[0].score() = H'[n].score();
        size_t const row_vector_size = row_vector.size() - 1;
        column_vector[0].score() = row_vector[row_vector_size].score();

        for (size_t j = row_vector_size; j > 0; --j)
            row_vector[j].score() = row_vector[j - 1].score();
    }

    template <typename column_vector_t, typename row_vector_t>
    constexpr void postprocess_block(column_vector_t const & column_vector, row_vector_t && row_vector) const noexcept
    {
        size_t const row_vector_size = row_vector.size() - 1;
        for (size_t j = 0; j < row_vector_size; ++j)
            row_vector[j].score() = row_vector[j + 1].score();

        row_vector[row_vector_size].score() = column_vector[column_vector.size() - 1].score();
    }

    constexpr derived_t const & as_derived() const noexcept
    {
        return static_cast<derived_t const &>(*this);
    }

    constexpr derived_t & as_derived() noexcept
    {
        return static_cast<derived_t &>(*this);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
