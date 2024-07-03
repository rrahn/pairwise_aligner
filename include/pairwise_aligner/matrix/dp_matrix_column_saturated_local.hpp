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

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {
namespace _column_saturated_local {

template <typename dp_vector_t> // saturated vector
class _wrapper : public dp_matrix::_column_saturated::_wrapper<dp_vector_t>
{
private:
    using base_t = dp_matrix::_column_saturated::_wrapper<dp_vector_t>;

    using score_t = typename base_t::value_type::score_type;
    using saturated_mask_t = typename score_t::mask_type;

    saturated_mask_t _is_local{true}; // Begin always in local block.

public:

    using typename base_t::range_type;
    using typename base_t::value_type;
    using typename base_t::reference;
    using typename base_t::const_reference;

    _wrapper() = delete;
    explicit _wrapper(dp_vector_t & dp_vector) : base_t{dp_vector}
    {}

    constexpr void update_offset() noexcept
    {
        score_t new_offset = (*this)[is_row_cell_v<value_type>].score();
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
        using regular_cell_t = typename dp_vector_t::value_type;
        using regular_score_t = typename regular_cell_t::score_type;

        regular_score_t new_offset_regular{new_offset};
        regular_score_t absolute_values = base_t::base().offset() + new_offset_regular - base_t::base().local_zero_offset();

        // 2. compare with threshold
        auto is_local = absolute_values.lt(base_t::base().threshold());
        _is_local = saturated_mask_t{is_local};

        // 3. set global_offsets and local offsets
        assert(base_t::check_saturated_arithmetic(blend(_is_local, base_t::base().saturated_zero_offset(), new_offset)));
        base_t::reset(blend(_is_local, base_t::base().saturated_zero_offset(), new_offset));

        // 4. update global offset
        base_t::base().update_offset(new_offset_regular, is_local);
    }
};

} // namespace _column_saturated_local

inline namespace _cpo {
inline constexpr dp_matrix::_column_saturated::_fn<lazy_type<_column_saturated_local::_wrapper>>
    column_saturated_local{};

} // inline namespace _cpo
} // namespace dp_matrix
} // inline namespace v1
} // namespace seqan::pairwise_aligner
