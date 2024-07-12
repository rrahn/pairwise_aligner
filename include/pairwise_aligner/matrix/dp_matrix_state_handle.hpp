// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix::detail::state_handle.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <tuple>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {
namespace detail {

enum struct dp_state_accessor
{
    id_dp_column = 0,
    id_dp_row = 1,
    id_column_sequence = 2,
    id_row_sequence = 3,
    id_substitution_model = 4,
    id_tracker = 5,
};

template <dp_state_accessor accessor_id>
inline constexpr size_t dp_state_accessor_id_v = static_cast<size_t>(accessor_id);

template <typename ...dp_state_t>
class state_handle
{
private:

    static_assert(sizeof...(dp_state_t) == 6, "ERROR: Too few/many arguments for the dp state!");

    using data_as_tuple_t = std::tuple<dp_state_t...>;
    data_as_tuple_t _dp_state_as_tuple;

protected:

    using fwd_dp_column_t =
        std::tuple_element_t<dp_state_accessor_id_v<dp_state_accessor::id_dp_column>, data_as_tuple_t>;
    using fwd_dp_row_t =
        std::tuple_element_t<dp_state_accessor_id_v<dp_state_accessor::id_dp_row>, data_as_tuple_t>;
    using fwd_column_sequence_t =
        std::tuple_element_t<dp_state_accessor_id_v<dp_state_accessor::id_column_sequence>, data_as_tuple_t>;
    using fwd_row_sequence_t =
        std::tuple_element_t<dp_state_accessor_id_v<dp_state_accessor::id_row_sequence>, data_as_tuple_t>;
    using fwd_substitution_model_t =
        std::tuple_element_t<dp_state_accessor_id_v<dp_state_accessor::id_substitution_model>, data_as_tuple_t>;
    using fwd_tracker_t =
        std::tuple_element_t<dp_state_accessor_id_v<dp_state_accessor::id_tracker>, data_as_tuple_t>;

public:

    using dp_column_type = std::remove_reference_t<fwd_dp_column_t>;
    using dp_row_type = std::remove_reference_t<fwd_dp_row_t>;
    using column_sequence_type = std::remove_reference_t<fwd_column_sequence_t>;
    using row_sequence_type = std::remove_reference_t<fwd_row_sequence_t>;
    using substitution_model_type = std::remove_reference_t<fwd_substitution_model_t>;
    using tracker_type = std::remove_reference_t<fwd_tracker_t>;

    state_handle() = default;
    explicit constexpr state_handle(dp_state_t ...dp_data) noexcept :
        _dp_state_as_tuple{std::forward<dp_state_t>(dp_data)...}
    {}

    constexpr dp_row_type & dp_row() noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_dp_row>>(_dp_state_as_tuple);
    }

    constexpr dp_row_type const & dp_row() const noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_dp_row>>(_dp_state_as_tuple);
    }

    constexpr dp_column_type & dp_column() noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_dp_column>>(_dp_state_as_tuple);
    }

    constexpr dp_column_type const & dp_column() const noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_dp_column>>(_dp_state_as_tuple);
    }

    constexpr column_sequence_type & column_sequence() noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_column_sequence>>(_dp_state_as_tuple);
    }

    constexpr column_sequence_type const & column_sequence() const noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_column_sequence>>(_dp_state_as_tuple);
    }

    constexpr row_sequence_type & row_sequence() noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_row_sequence>>(_dp_state_as_tuple);
    }

    constexpr row_sequence_type const & row_sequence() const noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_row_sequence>>(_dp_state_as_tuple);
    }

    constexpr substitution_model_type & substitution_model() noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_substitution_model>>(_dp_state_as_tuple);
    }

    constexpr substitution_model_type const & substitution_model() const noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_substitution_model>>(_dp_state_as_tuple);
    }

    constexpr tracker_type & tracker() noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_tracker>>(_dp_state_as_tuple);
    }

    constexpr tracker_type const & tracker() const noexcept
    {
        return get<dp_state_accessor_id_v<dp_state_accessor::id_tracker>>(_dp_state_as_tuple);
    }

    constexpr std::ptrdiff_t column_count() const noexcept
    {
        return dp_row().size();
    }

    constexpr std::ptrdiff_t row_count() const noexcept
    {
        return dp_column().size();
    }
};

template <typename ...args_t>
inline constexpr auto make_dp_state(args_t && ...args) -> state_handle<args_t...>
{
    return state_handle<args_t...>{std::forward<args_t>(args)...};
}

} // namespace detail
} // namespace dp_matrix
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
