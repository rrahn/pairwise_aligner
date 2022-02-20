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

#include <pairwise_aligner/matrix/dp_matrix_block.hpp>
#include <pairwise_aligner/matrix/dp_matrix_column.hpp>
#include <pairwise_aligner/type_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

namespace detail {

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

    constexpr void update_offset() noexcept
    {
        score_t new_offset = (*this)[is_row_cell_v<value_type>].score();
        assert(check_saturated_arithmetic(new_offset));
        update_offset_impl(new_offset);
    }

protected:

    constexpr void update_offset_impl(score_t const & new_offset) noexcept
    {
        reset(new_offset);
        _dp_vector.update_offset(new_offset);
    }

    void reset(score_t new_offset) noexcept
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
                using large_score_t = simd_score<wide_scalar_score_t, score_t::size>;
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

                for (size_t k = 0; k < score_t::size; ++k) {
                    if (expected_score[k] != real_score[k])
                        throw_error(k);
                }

                // TODO: Make generic for different alignment cell types.
                if (i > 0) { // Check also the gap costs for all i > 0.
                    expected_score = large_score_t{get<1>(range()[i])} - large_score_t{new_offset};
                    expected_score += large_score_t{_dp_vector.saturated_zero_offset()};

                    real_score = get<1>(range()[i]) - new_offset;
                    real_score += _dp_vector.saturated_zero_offset();

                    for (size_t k = 0; k < score_t::size; ++k) {
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
} // namespace detail

template <typename ...dp_data_t>
struct _column_saturated
{
    class type;
};

template <typename ...dp_data_t>
using column_saturated_t = typename _column_saturated<dp_data_t...>::type;

template <typename ...dp_data_t>
class _column_saturated<dp_data_t...>::type : public dp_matrix::column_t<dp_data_t...>
{
    using base_t = dp_matrix::column_t<dp_data_t...>;

public:

    type() = delete;
    constexpr explicit type(dp_data_t ...dp_data) noexcept : base_t{std::forward<dp_data_t>(dp_data)...}
    {}
    ~type() = default;

    constexpr auto operator[](size_t const index) noexcept
    {
        assert(index < base_t::size());

        detail::saturated_wrapper saturated_column{base_t::column()[index]};
        saturated_column.update_offset();
        base_t::row().update_offset();
        return dp_matrix_block(std::move(saturated_column),
                               base_t::row(),
                               base_t::substitution_model(),
                               base_t::tracker(),
                               base_t::row_sequence());
    }
};

namespace cpo {
struct _column_saturated_closure
{
    template <typename dp_column_t, typename dp_row_t, typename ...remaining_dp_data_t>
    constexpr auto operator()(dp_column_t && dp_column,
                              dp_row_t & dp_row,
                              remaining_dp_data_t && ...remaining_dp_data) const noexcept {
        using dp_saturated_column_t =
            dp_matrix::column_saturated_t<dp_column_t, detail::saturated_wrapper<dp_row_t>, remaining_dp_data_t...>;

        return dp_saturated_column_t{std::forward<dp_column_t>(dp_column),
                                     detail::saturated_wrapper{dp_row},
                                     std::forward<remaining_dp_data_t>(remaining_dp_data)...};
    }
};

} // namespace cpo
} // namespace dp_matrix

inline constexpr dp_matrix::cpo::_column_saturated_closure dp_matrix_column_saturated{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
