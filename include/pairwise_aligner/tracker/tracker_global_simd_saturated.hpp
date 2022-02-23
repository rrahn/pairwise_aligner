// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::tracker::global_simd_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/configuration/end_gap_policy.hpp>
#include <pairwise_aligner/simd/simd_score_type.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace tracker::global_simd_saturated {

namespace detail
{

template <typename score_t, typename dp_column_t, typename dp_row_t>
struct _saturated_max_score_finder
{
    class type;
};

template <typename score_t, typename dp_column_t, typename dp_row_t>
using saturated_max_score_finder = typename _saturated_max_score_finder<score_t, dp_column_t, dp_row_t>::type;

template <typename score_t, typename dp_column_t, typename dp_row_t>
class _saturated_max_score_finder<score_t, dp_column_t, dp_row_t>::type
{
private:

    using vec_int8_t = simd_score<int8_t>;
    using vec_uint8_t = simd_score<uint8_t>;
    using mask_uint8_t = typename vec_uint8_t::mask_type;
    using scalar_t = typename score_t::value_type;

    inline static constexpr size_t level1_max_size_v = (1ull << 8) - 1; // maximal size of sequence usable with level1
    inline static constexpr size_t level2_max_size_v = (1ull << 16) - 1; // maximal size of sequence usable with level2
    inline static constexpr size_t level3_max_size_v = (1ull << 24) - 1; // maximal size of sequence usable with level3

    struct level_state {
        score_t begin_position{};
        score_t end_position{};
        mask_uint8_t reached_first{~mask_uint8_t{}};
        mask_uint8_t reached_last{~mask_uint8_t{}};

        size_t chunk_idx{};
        size_t chunk_position{};
        size_t chunk_end{};
    };

    // Extracted parameters
    score_t _sequence1_sizes{};
    score_t _sequence2_sizes{};
    dp_column_t _dp_column;
    dp_row_t _dp_row;
    size_t _chunk_size{};
    size_t _dp_column_size{};
    size_t _dp_row_size{};

public:
    type() = delete;

    template <typename sequence1_bulk_t, typename sequence2_bulk_t>
    constexpr explicit type(dp_column_t dp_column,
                            dp_row_t dp_row,
                            sequence1_bulk_t && sequence1_bulk,
                            sequence2_bulk_t && sequence2_bulk,
                            size_t const chink_size) :
        _dp_column{std::forward<dp_column_t>(dp_column)},
        _dp_row{std::forward<dp_row_t>(dp_row)},
        _chunk_size{chink_size}
    {
        // Determine the maximal column and row size and store the sequence sizes in a simd vector.
        auto get_size = [] (auto && sequence) { return std::ranges::distance(sequence); };
        auto sequence1_sizes = sequence1_bulk | std::views::transform(get_size);
        auto sequence2_sizes = sequence2_bulk | std::views::transform(get_size);

        for (std::ptrdiff_t idx = 0; idx < std::ranges::distance(sequence1_sizes) ; ++idx) {
            _sequence1_sizes[idx] = sequence1_sizes[idx];
            _sequence2_sizes[idx] = sequence2_sizes[idx];
            _dp_column_size = std::max<size_t>(_dp_column_size, sequence1_sizes[idx] + 1);
            _dp_row_size = std::max<size_t>(_dp_row_size, sequence2_sizes[idx] + 1);
        }

        if (_dp_column_size > level3_max_size_v || _dp_row_size > level3_max_size_v)
            throw std::runtime_error{"The given dynamic programming matrix exceeds the maximal allowed column and/or "
                                     " row size of 2^24 - 1."};
    }

    constexpr auto in_column(score_t const & padding_score) const noexcept
    {
        auto [column_offsets, row_offsets] = get_offsets();
        return max(run_first(_dp_column, _dp_column_size, _sequence1_sizes, column_offsets, padding_score),
                   run_second(_dp_row, _dp_row_size, _sequence2_sizes, row_offsets, padding_score));
    }

    constexpr auto in_row(score_t const & padding_score) const noexcept
    {
        auto [column_offsets, row_offsets] = get_offsets();
        return max(run_first(_dp_row, _dp_row_size, _sequence2_sizes, row_offsets, padding_score),
                   run_second(_dp_column, _dp_column_size, _sequence1_sizes, column_offsets, padding_score));
    }

private:

    constexpr auto get_offsets() const noexcept
    {
        return std::pair{score_t{static_cast<scalar_t>(_dp_row_size - 1)} - _sequence2_sizes,
                         score_t{static_cast<scalar_t>(_dp_column_size - 1)} - _sequence1_sizes};
    }

    template <typename dp_vector_t>
    auto run_first(dp_vector_t const & dp_vector,
                   size_t const dp_vector_size,
                   score_t const & sequence_sizes,
                   score_t const & vector_offsets,
                   score_t const & padding_score) const noexcept
    {
        // Global states
        score_t best_score{std::numeric_limits<scalar_t>::lowest()};
        score_t const score_correction = padding_score * vector_offsets;
        score_t const mask_infinity{std::numeric_limits<int8_t>::lowest()};
        size_t const chunk_count = dp_vector.size();

        vec_int8_t local_max_score{std::numeric_limits<int8_t>::lowest()};
        mask_uint8_t reached_first{};
        mask_uint8_t reached_last{};

        auto find = [&] (level_state & state) -> bool {

            if (state.chunk_position == state.chunk_end) {
                // Get max best score.
                auto in_range = mask_infinity.lt(score_t{local_max_score});
                best_score = mask_max(best_score,
                                      in_range,
                                      best_score,
                                      (score_t{local_max_score} + dp_vector[state.chunk_idx].offset()) -
                                        score_correction);

                if (++state.chunk_idx == chunk_count) {
                    return false;
                }

                local_max_score = vec_int8_t{std::numeric_limits<int8_t>::lowest()};
                state.chunk_end = dp_vector[state.chunk_idx].size() - (state.chunk_idx < chunk_count - 1);
                state.chunk_position = 0;
            }

            reached_first |= state.reached_first;
            reached_last |= state.reached_last;
            mask_uint8_t mask{reached_first & ~reached_last};

            auto const & base_chunk = dp_vector[state.chunk_idx].base(); // TODO: No guarantee that base is what it is supposed to be!
            local_max_score = mask_max(local_max_score, mask, local_max_score, base_chunk[state.chunk_position].score());
            return true;
        };

        level_state state{
            .begin_position = vector_offsets,
            .end_position = min(vector_offsets + sequence_sizes + score_t{1},
                                score_t{static_cast<scalar_t>(dp_vector_size)}),
            .chunk_end = dp_vector[0].size() - (0 < chunk_count - 1)
        };

        if (dp_vector_size > level2_max_size_v) {
            run_level3(state, find);
        } else if (dp_vector_size > level1_max_size_v) {
            run_level2(state, find);
        } else {
            run_level1(state, find);
        }

        return best_score;
    }

    template <typename dp_vector_t>
    auto run_second(dp_vector_t const & dp_vector,
                    size_t const dp_vector_size,
                    score_t const & sequence_sizes,
                    score_t const & vector_offsets,
                    score_t const & padding_score) const noexcept
    {
        // Global state!
        score_t best_score{std::numeric_limits<scalar_t>::lowest()};
        score_t const mask_infinity{std::numeric_limits<int8_t>::lowest()};
        size_t const chunk_count = dp_vector.size();

        vec_int8_t local_max_score{std::numeric_limits<int8_t>::lowest()};
        vec_int8_t const local_padding_score{padding_score};
        vec_int8_t local_score_correction{0};
        mask_uint8_t reached_first{};

        auto find = [&] (level_state & state) -> bool {
            if (state.chunk_position == state.chunk_end) {
                auto in_range = mask_infinity.lt(score_t{local_max_score});
                score_t score_correction =
                    max(score_t{static_cast<scalar_t>(state.chunk_idx * _chunk_size)} - sequence_sizes, vector_offsets) *
                        padding_score;
                best_score = mask_max(best_score,
                                      in_range,
                                      best_score,
                                      (score_t{local_max_score} + dp_vector[state.chunk_idx].offset() -
                                        score_correction));

                if (++state.chunk_idx == chunk_count)
                    return false;

                // Reinitialise loop variables.
                state.chunk_end = dp_vector[state.chunk_idx].size() - 1;
                local_max_score = vec_int8_t{std::numeric_limits<int8_t>::lowest()};
                local_score_correction = vec_int8_t{0};
                state.chunk_position = 0;
            }

            reached_first |= state.reached_first;
            auto const & base_chunk = dp_vector[state.chunk_idx].base();
            local_max_score = mask_max(local_max_score,
                                       reached_first,
                                       local_max_score,
                                       base_chunk[state.chunk_position].score() - local_score_correction);
            local_score_correction = mask_add(local_score_correction,
                                              reached_first,
                                              local_score_correction,
                                              local_padding_score);
            return true;
        };

        level_state state{
            .begin_position = min(sequence_sizes + vector_offsets, score_t{static_cast<scalar_t>(dp_vector_size)}),
            .chunk_end = dp_vector[0].size() - 1
        };

        if (dp_vector_size > level2_max_size_v) {
            run_level3(state, find);
        } else if (dp_vector_size > level1_max_size_v) {
            run_level2(state, find);
        } else {
            run_level1(state, find);
        }

        return best_score;
    }

    template <typename fn_t>
    void run_level3(level_state & state, fn_t && fn) const noexcept
    {
        vec_uint8_t begin_position = vec_uint8_t{state.begin_position >> 16};
        vec_uint8_t end_position = vec_uint8_t{state.end_position >> 16};
        vec_uint8_t vec_level3_idx{0};
        for (vec_uint8_t vec_level3_idx{0};; ++vec_level3_idx) {
            state.reached_first = begin_position.le(vec_level3_idx);
            state.reached_last = end_position.le(vec_level3_idx);
            if (!run_level2(state, std::forward<fn_t>(fn)))
                break; // terminate.
        }
    }

    template <typename fn_t>
    auto run_level2(level_state & state, fn_t && fn) const noexcept
    {
        mask_uint8_t parent_reached_first = state.reached_first;
        mask_uint8_t parent_reached_last = state.reached_last;
        vec_uint8_t begin_position = vec_uint8_t{(state.begin_position >> 8) & score_t{static_cast<scalar_t>(255)}};
        vec_uint8_t end_position = vec_uint8_t{(state.end_position >> 8) & score_t{static_cast<scalar_t>(255)}};
        vec_uint8_t vec_level2_idx{0};
        for (size_t idx = 0; idx < 256; ++idx, ++vec_level2_idx) {
            state.reached_first = parent_reached_first & begin_position.le(vec_level2_idx);
            state.reached_last = parent_reached_last & end_position.le(vec_level2_idx);

            if (!run_level1(state, std::forward<fn_t>(fn)))
                return false; // terminate.
        }
        return true;
    }

    template <typename fn_t>
    auto run_level1(level_state & state, fn_t && fn) const noexcept
    {
        mask_uint8_t parent_reached_first = state.reached_first;
        mask_uint8_t parent_reached_last = state.reached_last;
        vec_uint8_t begin_position = vec_uint8_t{state.begin_position & score_t{static_cast<scalar_t>(255)}};
        vec_uint8_t end_position = vec_uint8_t{state.end_position & score_t{static_cast<scalar_t>(255)}};
        for (size_t level1_idx = 0; level1_idx < 256; ++level1_idx, ++state.chunk_position) {
            vec_uint8_t vec_level1_idx{static_cast<uint8_t>(level1_idx)};

            // it is a mask and we later add the mask to it.
            state.reached_first = parent_reached_first & begin_position.le(vec_level1_idx);
            state.reached_last = parent_reached_last & end_position.le(vec_level1_idx);

            if (!fn(state))
                return false; // terminate.
        }
        return true;
    }

};

} // namespace detail

template <typename score_t>
struct _tracker
{
    class type;
};

template <typename score_t>
using tracker = typename _tracker<score_t>::type;

template <typename score_t>
class _tracker<score_t>::type
{
public:
    score_t _padding_score;
    cfg::trailing_end_gap _end_gap;
    size_t _chunk_size{};

    template <typename saturated_score_t>
    constexpr saturated_score_t const & track(saturated_score_t const & score) const noexcept {
        return score; // no-op.
    }

    template <typename sequence1_t, typename sequences2_t, typename dp_column_t, typename dp_row_t>
    constexpr score_t max_score(sequence1_t && sequence1,
                                sequences2_t && sequences2,
                                dp_column_t const & dp_column,
                                dp_row_t const & dp_row) const noexcept
    {
        std::vector<std::views::all_t<sequence1_t>> sequence1_bulk{};
        sequence1_bulk.resize(std::ranges::distance(sequences2), sequence1 | std::views::all);

        return max_score(std::move(sequence1_bulk), std::forward<sequences2_t>(sequences2), dp_column, dp_row);
    }

    template <typename sequences1_t, typename sequences2_t, typename dp_column_t, typename dp_row_t>
        requires std::ranges::range<std::ranges::range_reference_t<sequences1_t>>
    constexpr score_t max_score(sequences1_t && sequences1,
                                sequences2_t && sequences2,
                                dp_column_t const & dp_column,
                                dp_row_t const & dp_row) const noexcept {
        if (_end_gap.last_column == cfg::end_gap::penalised && _end_gap.last_row == cfg::end_gap::penalised) {
            score_t best_score{};
            for (std::ptrdiff_t idx = 0; idx < std::ranges::distance(sequences1); ++idx) {
                best_score[idx] = select_max_score(idx, sequences1[idx], sequences2[idx], dp_column, dp_row);
            }
            return best_score;
        }

        using max_score_finder_t = detail::saturated_max_score_finder<score_t, dp_column_t const &, dp_row_t const &>;
        max_score_finder_t max_score_finder{dp_column, dp_row, sequences1, sequences2, _chunk_size};

        score_t best_score{std::numeric_limits<typename score_t::value_type>::lowest()};
        if (_end_gap.last_column == cfg::end_gap::free) {
            best_score = max_score_finder.in_column(_padding_score);
        }

        if (_end_gap.last_row == cfg::end_gap::free) {
            best_score = max(best_score, max_score_finder.in_row(_padding_score));
        }
        return best_score;
    }

    // TODO: optimal_coordinate()

private:

    template <typename sequence_t, typename dp_vector_t>
    constexpr auto get_offsets(sequence_t && sequence, dp_vector_t && dp_vector) const noexcept
    {
        assert(dp_vector.size() != 0);
        size_t sequence_size = std::ranges::distance(sequence);

        size_t chunk_size = dp_vector[0].size() - 1;
        size_t const full_chunks = dp_vector.size() - 1;
        size_t dp_vector_size = (full_chunks * chunk_size) + dp_vector[full_chunks].size();

        assert(dp_vector_size > sequence_size);
        size_t offset = dp_vector_size - 1 - sequence_size;

        return std::tuple{sequence_size, dp_vector_size, offset};
    }


    template <typename sequence1_t, typename sequence2_t, typename dp_column_t, typename dp_row_t>
    constexpr auto select_max_score(size_t const simd_idx,
                                    sequence1_t && sequence1,
                                    sequence2_t && sequence2,
                                    dp_column_t const & dp_column,
                                    dp_row_t const & dp_row) const noexcept
        -> typename score_t::value_type
    {
        auto [column_sequence_size, column_vector_size, row_offset] = get_offsets(sequence1, dp_column);
        auto [row_sequence_size, row_vector_size, column_offset] = get_offsets(sequence2, dp_row);

        using scalar_t = typename score_t::value_type;
        size_t scale = std::min(column_offset, row_offset);
        scalar_t best_score{};

        size_t const chunk_size = dp_column[0].size() - 1;

        if (scale == column_offset) {
            auto [chunk_id, chunk_position] =
                to_local_position(column_sequence_size + column_offset, chunk_size, dp_column.size());
            best_score = score_at(dp_column[chunk_id][chunk_position], simd_idx);
        } else {
            auto [chunk_id, chunk_position] =
                to_local_position(row_sequence_size + row_offset, chunk_size, dp_row.size());
            best_score = score_at(dp_row[chunk_id][chunk_position], simd_idx);
        }

        return best_score - (_padding_score[simd_idx] * scale);
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

    constexpr std::pair<size_t, size_t> to_local_position(size_t const position, size_t const chunk_size, size_t const chunk_count)
        const noexcept
    {
        size_t idx = position / chunk_size;
        size_t offset = position % chunk_size;
        // The last block has one more cell to account for the initialisation cell.
        // If the last given position of the currently processed vector is a multiple of the chunk size, then
        // it is the stored inside the last chunk. The respective index and offset are corrected accordingly.
        size_t is_last_position = (idx == chunk_count);
        return {idx - is_last_position, offset + (chunk_size * is_last_position)};
    }
};

template <typename score_t>
struct _factory
{
    struct type;
};

template <typename score_t>
using factory = typename _factory<score_t>::type;

template <typename score_t>
struct _factory<score_t>::type
{
    // params for free end-gaps.
    score_t _padding_score{};
    cfg::trailing_end_gap _end_gap{};
    size_t _chunk_size{};

    constexpr auto make_tracker() const noexcept {
        return tracker<score_t>{_padding_score, _end_gap, _chunk_size};
    }
};

} // namespace tracker::global_simd_saturated

} // inline namespace v1
} // namespace seqan::pairwise_aligner
