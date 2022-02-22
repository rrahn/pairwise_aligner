// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_algorithm_attorney.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace seqan::pairwise_aligner
{
inline namespace v1
{

// ----------------------------------------------------------------------------
// Forward declarations of classes granted access to algorithm via attorney.
// ----------------------------------------------------------------------------

template <typename dp_algorithm_impl_t>
struct _dp_algorithm_template_base;

// ----------------------------------------------------------------------------
// Definition of the algorithm attorney managing access to the client.
// ----------------------------------------------------------------------------

template <typename algorithm_client_t>
class dp_algorithm_attorney
{
private:
    // Classes that have been granted access (grantees) to the algorithm implementation (client/grantor).
    friend _dp_algorithm_template_base<algorithm_client_t>;

    // Member functions the grantees can access.
    template <typename ...args_t>
    constexpr static auto compute_cell(algorithm_client_t const & client, args_t && ...args)
        noexcept(noexcept(client.compute_cell(std::forward<args_t>(args)...)))
    {
        return client.compute_cell(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr static auto make_result(algorithm_client_t const & client, args_t && ...args)
        noexcept(noexcept(client.make_result(std::forward<args_t>(args)...)))
    {
        return client.make_result(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr static auto initialise_substitution_scheme(algorithm_client_t const & client, args_t && ...args)
        noexcept(noexcept(client.initialise_substitution_scheme(std::forward<args_t>(args)...)))
    {
        return client.initialise_substitution_scheme(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr static auto initialise_tracker(algorithm_client_t const & client, args_t && ...args)
        noexcept(noexcept(client.initialise_tracker(std::forward<args_t>(args)...)))
    {
        return client.initialise_tracker(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr static auto initialise_policies(algorithm_client_t const & client, args_t && ...args)
        noexcept(noexcept(client.initialise_policies(std::forward<args_t>(args)...)))
    {
        return client.initialise_policies(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr static auto initialise_column_vector(algorithm_client_t const & client, args_t && ...args)
        noexcept(noexcept(client.initialise_column_vector(std::forward<args_t>(args)...)))
    {
        return client.initialise_column_vector(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr static auto initialise_row_vector(algorithm_client_t const & client, args_t && ...args)
        noexcept(noexcept(client.initialise_row_vector(std::forward<args_t>(args)...)))
    {
        return client.initialise_row_vector(std::forward<args_t>(args)...);
    }

    template <typename ...args_t>
    constexpr static auto lane_width(algorithm_client_t const & client, args_t && ...args)
        noexcept(noexcept(client.lane_width(std::forward<args_t>(args)...)))
    {
        return client.lane_width(std::forward<args_t>(args)...);
    }
};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
