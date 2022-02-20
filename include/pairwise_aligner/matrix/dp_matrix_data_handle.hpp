// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_data_handle.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <tuple>

#include <pairwise_aligner/matrix/dp_matrix_block.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

namespace detail {

template <typename ...dp_data_t>
struct _dp_matrix_data_handle
{
    class type;
};

template <typename ...dp_data_t>
using dp_matrix_data_handle = typename _dp_matrix_data_handle<dp_data_t...>::type;

template <typename ...dp_data_t>
class _dp_matrix_data_handle<dp_data_t...>::type
{
private:

    enum accessor_id
    {
        _column = 0,
        _row = 1,
        _substitution_model = 2,
        _tracker = 3,
    };

    using data_as_tuple_t = std::tuple<dp_data_t...>;
    data_as_tuple_t _dp_data_as_tuple;

protected:

    using fwd_dp_column_t = std::tuple_element_t<accessor_id::_column, data_as_tuple_t>;
    using fwd_dp_row_t = std::tuple_element_t<accessor_id::_row, data_as_tuple_t>;
    using fwd_substitution_model_t = std::tuple_element_t<accessor_id::_substitution_model, data_as_tuple_t>;
    using fwd_tracker_t = std::tuple_element_t<accessor_id::_tracker, data_as_tuple_t>;

public:
    using column_type = std::remove_reference_t<fwd_dp_column_t>;
    using row_type = std::remove_reference_t<fwd_dp_row_t>;
    using substitution_model_type = std::remove_reference_t<fwd_substitution_model_t>;
    using tracker_type = std::remove_reference_t<fwd_tracker_t>;

    type() = default;
    explicit constexpr type(dp_data_t ...dp_data) noexcept :
        _dp_data_as_tuple{std::forward<dp_data_t>(dp_data)...}
    {}

    constexpr row_type & row() noexcept
    {
        return get<accessor_id::_row>(_dp_data_as_tuple);
    }

    constexpr row_type const & row() const noexcept
    {
        return get<accessor_id::_row>(_dp_data_as_tuple);
    }

    constexpr column_type & column() noexcept
    {
        return get<accessor_id::_column>(_dp_data_as_tuple);
    }

    constexpr column_type const & column() const noexcept
    {
        return get<accessor_id::_column>(_dp_data_as_tuple);
    }

    constexpr substitution_model_type & substitution_model() noexcept
    {
        return get<accessor_id::_substitution_model>(_dp_data_as_tuple);
    }

    constexpr substitution_model_type const & substitution_model() const noexcept
    {
        return get<accessor_id::_substitution_model>(_dp_data_as_tuple);
    }

    constexpr tracker_type & tracker() noexcept
    {
        return get<accessor_id::_tracker>(_dp_data_as_tuple);
    }

    constexpr tracker_type const & tracker() const noexcept
    {
        return get<accessor_id::_tracker>(_dp_data_as_tuple);
    }
};

} // namespace detail
} // namespace dp_matrix
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
