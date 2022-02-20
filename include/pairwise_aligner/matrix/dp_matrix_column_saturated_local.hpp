// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_column_saturated_local.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/matrix/dp_matrix_column_saturated.hpp>
#include <pairwise_aligner/type_traits.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

namespace detail {

template <typename dp_vector_t> // saturated vector
class saturated_wrapper_local : public saturated_wrapper<dp_vector_t>
{
private:
    using base_t = saturated_wrapper<dp_vector_t>;

    using score_t = typename base_t::value_type::score_type;
    using saturated_mask_t = typename score_t::mask_type;

    saturated_mask_t _is_local{true}; // Begin always in local block.

public:

    using typename base_t::range_type;
    using typename base_t::value_type;
    using typename base_t::reference;
    using typename base_t::const_reference;

    saturated_wrapper_local() = delete;
    explicit saturated_wrapper_local(dp_vector_t & dp_vector) : base_t{dp_vector}
    {}

    constexpr void update_offset() noexcept
    {
        score_t new_offset = (*this)[is_row_cell_v<value_type>].score();
        assert(base_t::check_saturated_arithmetic(new_offset));
        update_offset_impl(new_offset);
    }

    constexpr saturated_mask_t const & is_local() const noexcept
    {
        return _is_local;
    }

    constexpr auto to_regular_score(score_t const & score) const noexcept
    {
        return base_t::base().to_regular_score(score);
    }
protected:

    constexpr void update_offset_impl(score_t const & new_offset) noexcept
    {
        // 1. get absolute values
        using regular_score_t = std::remove_cvref_t<decltype(base_t::base().local_zero_offset())>;

        regular_score_t new_offset_regular{new_offset};
        auto absolute_values = base_t::base().offset() + new_offset_regular - base_t::base().local_zero_offset();
        // 2. compare with threshold
        auto is_local = absolute_values.lt(base_t::base().threshold());
        _is_local = saturated_mask_t{is_local};

        // 3. set global_offsets and local offsets
        base_t::reset(blend(_is_local, base_t::base().saturated_zero_offset(), new_offset));

        // 4. update global offset
        base_t::base().update_offset(new_offset_regular, is_local);
    }
};
} // namespace detail


template <typename ...dp_data_t>
struct _column_saturated_local
{
    class type;
};

template <typename ...dp_data_t>
using column_saturated_local_t =
    typename _column_saturated_local<dp_data_t...>::type;

template <typename ...dp_data_t>
class _column_saturated_local<dp_data_t...>::type : public dp_matrix::column_t<dp_data_t...>
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

        detail::saturated_wrapper_local saturated_column{base_t::column()[index]};
        saturated_column.update_offset();
        base_t::row().update_offset();
        // we have regular simd local tracker and we wrap it into a local wrapper.
        return dp_matrix_block(std::move(saturated_column),
                               base_t::row(),
                               base_t::substitution_model().block_scheme(base_t::row().is_local()),
                               base_t::tracker().in_block_tracker([&](auto const & in_block_score) {
                                   return base_t::row().to_regular_score(in_block_score);
                               }));
    }
};

namespace cpo {
struct _column_saturated_local_closure
{
    template <typename dp_column_t, typename dp_row_t, typename ...remaining_dp_data_t>
    constexpr auto operator()(dp_column_t && dp_column,
                              dp_row_t & dp_row,
                              remaining_dp_data_t && ...remaining_dp_data) const noexcept {
        using column_saturated_local_t =
            dp_matrix::column_saturated_local_t<
                dp_column_t, detail::saturated_wrapper_local<dp_row_t>, remaining_dp_data_t...>;

        return column_saturated_local_t{std::forward<dp_column_t>(dp_column),
                                        detail::saturated_wrapper_local{dp_row},
                                        std::forward<remaining_dp_data_t>(remaining_dp_data)...};
    }
};

} // namespace cpo
} // namespace dp_matrix

inline constexpr dp_matrix::cpo::_column_saturated_local_closure dp_matrix_column_saturated_local{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
