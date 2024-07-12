// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_column.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cassert>
#include <iostream>

#include <pairwise_aligner/matrix/dp_matrix_column_base.hpp>
#include <pairwise_aligner/matrix/dp_matrix_state_handle.hpp>
#include <pairwise_aligner/utility/type_list.hpp>
#include <pairwise_aligner/type_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {
namespace _column_saturated {

template <typename dp_vector_t> // saturated vector
class _wrapper
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

    _wrapper() = delete;
    explicit _wrapper(dp_vector_t & dp_vector) : _dp_vector{dp_vector}
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

    constexpr void update_offset() noexcept
    {
        score_t new_offset = (*this)[is_row_cell_v<value_type>].score();
        assert(check_saturated_arithmetic(new_offset));
        update_offset_impl(new_offset);
    }

    constexpr decltype(auto) offset() const noexcept
    {
        return base().offset();
    }

protected:

    constexpr void update_offset_impl(score_t const & new_offset) noexcept
    {
        reset(new_offset);
        _dp_vector.update_offset(new_offset);
    }

    void reset(score_t const & new_offset) noexcept
    {
        for (size_t i = 0; i < size(); ++i)
            std::apply([&] (auto & ...values) {
                ((values -= new_offset), ...);
                ((values += _dp_vector.saturated_zero_offset()), ...);
            }, range()[i]);
    }

    constexpr bool check_saturated_arithmetic(score_t const & new_offset) const noexcept
    {
        bool test = true;
        using scalar_score_t = typename score_t::value_type;
        using wide_scalar_score_t = std::conditional_t<std::is_signed_v<scalar_score_t>, int32_t, uint32_t>;
        try {
            for (size_t i = 0; i < size(); ++i) {
                using large_score_t = simd_score<wide_scalar_score_t, score_t::size_v>;
                large_score_t expected_score = large_score_t{get<0>(range()[i])} - large_score_t{new_offset};
                expected_score += large_score_t{_dp_vector.saturated_zero_offset()};

                auto real_score = get<0>(range()[i]) - new_offset;
                real_score += _dp_vector.saturated_zero_offset();

                auto throw_error = [&] (size_t k) {
                    throw std::runtime_error{" i: " + std::to_string(i) +
                                             ", k: " + std::to_string(k) +
                                             ", real_score: " + std::to_string(real_score[k]) +
                                             ", expected_score: " + std::to_string(expected_score[k]) +
                                             ", cell: <" + std::to_string(get<0>(range()[i])[k]) + ", " +
                                                         std::to_string(get<1>(range()[i])[k]) + ">" +
                                             ", offset: " + std::to_string(new_offset[k]) +
                                             ", zero_offset: " + std::to_string(_dp_vector.saturated_zero_offset()[k])};
                };

                for (size_t k = 0; k < score_t::size_v; ++k) {
                    if (expected_score[k] != real_score[k])
                        throw_error(k);
                }

                // TODO: Make generic for different alignment cell types.
                if (i > 0) { // Check also the gap costs for all i > 0.
                    expected_score = large_score_t{get<1>(range()[i])} - large_score_t{new_offset};
                    expected_score += large_score_t{_dp_vector.saturated_zero_offset()};

                    real_score = get<1>(range()[i]) - new_offset;
                    real_score += _dp_vector.saturated_zero_offset();

                    for (size_t k = 0; k < score_t::size_v; ++k) {
                        if (expected_score[k] != real_score[k])
                            throw_error(k);
                    }
                }
            }
        } catch (std::exception const & ex) {
            std::cerr << "[ERROR] Updating the offset caused an arithmetic over- or underflow! " << ex.what() << "\n";
            test = false;
        }

        return test;
    }
};

template <typename lazy_wrapper_t, typename block_fn_t, typename dp_state_t>
class _type : public dp_matrix::detail::column_base<block_fn_t, dp_state_t>
{
    using base_t = dp_matrix::detail::column_base<block_fn_t, dp_state_t>;
    using dp_inner_column_t = std::remove_reference_t<decltype(std::declval<_type &>().dp_column()[0])>;
    using dp_wrapper_t = instantiate_t<lazy_wrapper_t, dp_inner_column_t>;

public:
    using base_t::base_t;

    constexpr auto row_at(std::ptrdiff_t const index) noexcept
        -> decltype(std::declval<_type &>().make_matrix_block(dp_wrapper_t{base_t::dp_column()[index]},
                                                              base_t::dp_row(),
                                                              base_t::column_slice_at(index),
                                                              base_t::row_sequence(),
                                                              base_t::substitution_model(),
                                                              base_t::tracker().in_block_tracker(base_t::dp_row().offset())))
    {
        assert(index < base_t::row_count());

        dp_wrapper_t saturated_column{base_t::dp_column()[index]};
        saturated_column.update_offset();
        base_t::dp_row().update_offset();
        return base_t::make_matrix_block(std::move(saturated_column),
                                         base_t::dp_row(),
                                         base_t::column_slice_at(index),
                                         base_t::row_sequence(),
                                         base_t::substitution_model(),
                                         base_t::tracker().in_block_tracker(base_t::dp_row().offset()));
    }
};

template <typename lazy_wrapper_t>
struct _fn
{
    template <typename dp_block_fn_t>
    constexpr auto operator()(dp_block_fn_t && dp_block_fn) const noexcept
    {
        std::tuple<dp_block_fn_t> tmp{std::forward<dp_block_fn_t>(dp_block_fn)};
        return [fwd_capture = std::move(tmp)] <typename dp_state_t> (dp_state_t && dp_state) {
            using fwd_dp_block_fn_t = std::tuple_element_t<0, decltype(fwd_capture)>;
            using dp_row_t = typename std::remove_reference_t<dp_state_t>::dp_row_type;
            using wrapped_dp_row_t = instantiate_t<lazy_wrapper_t, dp_row_t>;

            wrapped_dp_row_t wrapped_row{dp_state.dp_row()};
            auto modified_state = dp_matrix::detail::make_dp_state(
                                    std::move(dp_state).dp_column(),
                                    std::move(wrapped_row),
                                    std::move(dp_state).column_sequence(),
                                    std::move(dp_state).row_sequence(),
                                    std::move(dp_state).substitution_model(),
                                    std::move(dp_state).tracker());

            using column_t = _column_saturated::_type<lazy_wrapper_t, fwd_dp_block_fn_t, decltype(modified_state)>;
            return column_t{std::forward<fwd_dp_block_fn_t>(std::get<0>(fwd_capture)), std::move(modified_state)};
        };
    }
};
} // namespace _column_saturated

inline namespace _cpo {
inline constexpr dp_matrix::_column_saturated::_fn<lazy_type<_column_saturated::_wrapper>> column_saturated{};

} // inline namespace _cpo
} // namespace dp_matrix
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
