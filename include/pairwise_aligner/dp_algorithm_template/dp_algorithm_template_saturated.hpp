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

#include <seqan3/std/span>

#include <seqan3/utility/views/slice.hpp>

#include <pairwise_aligner/dp_algorithm_template/dp_algorithm_template_base.hpp>

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

} // namespace detail

template <typename algorithm_impl_t>
struct _dp_algorithm_template_saturated
{
    class type;
};

template <typename algorithm_impl_t>
using dp_algorithm_template_saturated = typename _dp_algorithm_template_saturated<algorithm_impl_t>::type;

template <typename algorithm_impl_t>
class _dp_algorithm_template_saturated<algorithm_impl_t>::type : public dp_algorithm_template_base<algorithm_impl_t>
{
private:
    using base_t = dp_algorithm_template_base<algorithm_impl_t>;

protected:

    template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
    auto run(sequence1_t && sequence1, sequence2_t && sequence2, dp_column_t dp_column, dp_row_t dp_row) const
    {
        // ----------------------------------------------------------------------------
        // Initialisation
        // ----------------------------------------------------------------------------

        auto simd_seq1 = base_t::initialise_column(sequence1, dp_column);
        auto simd_seq2 = base_t::initialise_row(sequence2, dp_row);
        auto tracker = base_t::initialise_tracker();

        // ----------------------------------------------------------------------------
        // Recursion
        // ----------------------------------------------------------------------------

        size_t const column_block_count = dp_column.size();
        size_t const row_block_count = dp_row.size();
        size_t const column_block_size = dp_column[0].size() - 1;
        size_t const row_block_size = dp_row[0].size() - 1;

        using dp_column_block_t = typename dp_column_t::value_type;
        using dp_row_block_t = typename dp_row_t::value_type;

        using value_t = std::ranges::range_value_t<decltype(simd_seq1)>;
        std::vector<std::span<value_t>> seq1_blocked{};
        seq1_blocked.reserve(column_block_count);

        using saturated_col_t = detail::saturated_wrapper<dp_column_block_t>;
        std::vector<saturated_col_t> dp_column_blocks{};
        dp_column_blocks.reserve(column_block_count);

        for (size_t i = 0; i < column_block_count; ++i)
        {
            seq1_blocked.emplace_back(std::ranges::next(std::ranges::begin(simd_seq1), (i * column_block_size)),
                                      std::ranges::next(std::ranges::begin(simd_seq1),
                                                        ((i + 1) * column_block_size),
                                                        std::ranges::end(simd_seq1)));
            dp_column_blocks.emplace_back(dp_column[i]);
        }

        for (size_t j = 0; j < row_block_count; ++j) {
            detail::saturated_wrapper<dp_row_block_t> current_row_vector{dp_row[j]};
            std::span transformed_seq2{std::ranges::next(std::ranges::begin(simd_seq2), (j * row_block_size)),
                                       std::ranges::next(std::ranges::begin(simd_seq2),
                                                         ((j + 1) * row_block_size),
                                                         std::ranges::end(simd_seq2))};

            // Initialise first block of current column.
            current_row_vector.offset(current_row_vector[0].score());
            dp_column_blocks[0].offset(dp_column_blocks[0][0].score());

            base_t::initialise_block(dp_column_blocks[0], current_row_vector);

            // Iterate over blocks in current column.
            for (size_t i = 0; i < column_block_count; ++i) {
                if (i > 0) {
                    current_row_vector[0].score() = dp_column_blocks[i - 1][dp_column_blocks[i - 1].size() - 1].score();
                    current_row_vector.offset(current_row_vector[1].score());
                    dp_column_blocks[i].offset(dp_column_blocks[i][0].score());
                    dp_column_blocks[i][0].score() = current_row_vector[0].score();
                }

                base_t::compute_block(seq1_blocked[i],
                                      transformed_seq2,
                                      dp_column_blocks[i],
                                      current_row_vector,
                                      tracker);
            }

            // Write back optimal score to row vector.
            base_t::postprocess_block(dp_column_blocks.back(), current_row_vector);
        }

        // ----------------------------------------------------------------------------
        // Create result
        // ----------------------------------------------------------------------------

        return base_t::make_result(std::move(tracker),
                                   std::forward<sequence1_t>(sequence1),
                                   std::forward<sequence2_t>(sequence2),
                                   std::move(dp_column),
                                   std::move(dp_row));
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
