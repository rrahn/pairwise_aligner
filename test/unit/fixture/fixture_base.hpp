// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/type_traits>

namespace pairwise_aligner::test {
// ----------------------------------------------------------------------------
// Base fixture: extracts the values object it is parameterised with
// ----------------------------------------------------------------------------

template <auto _fixture>
struct fixture : public ::testing::Test
{
    using test_values_type = std::remove_cvref_t<decltype(*_fixture)>;

    // Method in same naming scheme as used by Google Test
    auto GetParam() -> test_values_type const &
    {
        return *_fixture;
    }
};

} // namespace pairwise_aligner::test
