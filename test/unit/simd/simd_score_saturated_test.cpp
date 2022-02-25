// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <pairwise_aligner/simd/simd_score_type.hpp>
#include <pairwise_aligner/utility/math.hpp>

template <typename _scalar_t>
struct saturated_simd_test : public ::testing::Test
{
    using scalar_t = _scalar_t;
    using simd_score_t = seqan::pairwise_aligner::simd_score_saturated<scalar_t>;
};

using testing_types = ::testing::Types<int8_t, uint8_t, int16_t, uint16_t>;

TYPED_TEST_SUITE(saturated_simd_test, testing_types);

TYPED_TEST(saturated_simd_test, add_above)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_t = typename TestFixture::scalar_t;

    simd_score_t a{std::numeric_limits<scalar_t>::max()};
    simd_score_t b{static_cast<scalar_t>(1)};

    simd_score_t c = seqan::pairwise_aligner::add(a, b);

    for (size_t i = 0; i < simd_score_t::size; ++i)
        EXPECT_EQ(c[i], std::numeric_limits<scalar_t>::max());
}

TYPED_TEST(saturated_simd_test, add_above_scalar)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_t = typename TestFixture::scalar_t;

    simd_score_t a{std::numeric_limits<scalar_t>::max()};

    simd_score_t c = seqan::pairwise_aligner::add(a, 1);
    for (size_t i = 0; i < simd_score_t::size; ++i)
        EXPECT_EQ(c[i], std::numeric_limits<scalar_t>::max());
}

TYPED_TEST(saturated_simd_test, add_below)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_t = typename TestFixture::scalar_t;

    simd_score_t a{std::numeric_limits<scalar_t>::lowest() + 1};
    simd_score_t b{static_cast<scalar_t>(-2)};

    simd_score_t c = seqan::pairwise_aligner::add(a, b);

    auto expected = [] () {
        if constexpr (std::is_signed_v<scalar_t>)
            return std::numeric_limits<scalar_t>::lowest();
        else
            return std::numeric_limits<scalar_t>::max();
    };

    for (size_t i = 0; i < simd_score_t::size; ++i)
        EXPECT_EQ(c[i], expected());
}

TYPED_TEST(saturated_simd_test, add_below_scalar)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_t = typename TestFixture::scalar_t;

    simd_score_t a{std::numeric_limits<scalar_t>::lowest() + 1};

    simd_score_t c = seqan::pairwise_aligner::add(a, -2);

    auto expected = [] () {
        if constexpr (std::is_signed_v<scalar_t>)
            return std::numeric_limits<scalar_t>::lowest();
        else
            return std::numeric_limits<scalar_t>::max();
    };

    for (size_t i = 0; i < simd_score_t::size; ++i)
        EXPECT_EQ(c[i], expected());
}

TYPED_TEST(saturated_simd_test, subtract_above)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_t = typename TestFixture::scalar_t;

    simd_score_t a{std::numeric_limits<scalar_t>::max()};
    simd_score_t b{static_cast<scalar_t>(-1)};

    simd_score_t c = seqan::pairwise_aligner::subtract(a, b);

    auto expected = [] () {
        if constexpr (std::is_signed_v<scalar_t>)
            return std::numeric_limits<scalar_t>::max();
        else
            return std::numeric_limits<scalar_t>::lowest();
    };

    for (size_t i = 0; i < simd_score_t::size; ++i)
        EXPECT_EQ(c[i], expected());
}

TYPED_TEST(saturated_simd_test, subtract_above_scalar)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_t = typename TestFixture::scalar_t;

    simd_score_t a{std::numeric_limits<scalar_t>::max()};
    simd_score_t c = seqan::pairwise_aligner::subtract(a, -1);

    auto expected = [] () {
        if constexpr (std::is_signed_v<scalar_t>)
            return std::numeric_limits<scalar_t>::max();
        else
            return std::numeric_limits<scalar_t>::lowest();
    };

    for (size_t i = 0; i < simd_score_t::size; ++i)
        EXPECT_EQ(c[i], expected());
}

TYPED_TEST(saturated_simd_test, subtract_below)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_t = typename TestFixture::scalar_t;

    simd_score_t a{std::numeric_limits<scalar_t>::lowest()};
    simd_score_t b{static_cast<scalar_t>(1)};

    simd_score_t c = seqan::pairwise_aligner::subtract(a, b);

    for (size_t i = 0; i < simd_score_t::size; ++i)
        EXPECT_EQ(c[i], std::numeric_limits<scalar_t>::lowest());
}

TYPED_TEST(saturated_simd_test, subtract_below_scalar)
{
    using simd_score_t = typename TestFixture::simd_score_t;
    using scalar_t = typename TestFixture::scalar_t;

    simd_score_t a{std::numeric_limits<scalar_t>::lowest()};

    simd_score_t c = seqan::pairwise_aligner::subtract(a, 1);

    for (size_t i = 0; i < simd_score_t::size; ++i)
        EXPECT_EQ(c[i], std::numeric_limits<scalar_t>::lowest());
}
