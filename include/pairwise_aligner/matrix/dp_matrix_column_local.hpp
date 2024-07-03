// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_matrix_column_local.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/matrix/dp_matrix_block.hpp>
#include <pairwise_aligner/matrix/dp_matrix_column_base.hpp>
#include <pairwise_aligner/matrix/dp_matrix_data_handle_local.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{
namespace dp_matrix {

template <typename block_closure_t, typename ...dp_data_t>
struct _column_local
{
    class type;
};

template <typename block_closure_t, typename ...dp_data_t>
using column_local_t = typename _column_local<block_closure_t, dp_data_t...>::type;

template <typename block_closure_t, typename ...dp_data_t>
class _column_local<block_closure_t, dp_data_t...>::type :
    public detail::column_base_t<block_closure_t, detail::dp_matrix_data_handle_local<dp_data_t...>>
{
protected:
    using base_t = detail::column_base_t<block_closure_t, detail::dp_matrix_data_handle_local<dp_data_t...>>;
public:
    using base_t::base_t;
};

namespace cpo {

template <typename block_closure_t = dp_matrix::cpo::_block_closure<>>
struct _column_local_closure
{
    block_closure_t block_closure{};

    template <typename ...dp_data_t>
    constexpr auto operator()(dp_data_t && ...dp_data) const noexcept {
        using dp_column_t = dp_matrix::column_local_t<block_closure_t, dp_data_t...>;

        return dp_column_t{block_closure, std::forward<dp_data_t>(dp_data)...};
    }
};

} // namespace cpo
} // namespace dp_matrix

inline constexpr dp_matrix::cpo::_column_local_closure<> dp_matrix_column_local{};

} // namespace v1
} // namespace seqan::pairwise_aligner
