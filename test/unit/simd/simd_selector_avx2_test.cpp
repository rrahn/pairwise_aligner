// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <array>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <type_traits>

#include <pairwise_aligner/simd/simd_score_type.hpp>
#include <pairwise_aligner/simd/simd_selector.hpp>


namespace pa = seqan::pairwise_aligner;

template <typename test_param_t>
struct simd_selector_avx2_test : public testing::Test
{
    using scalar_value_t = std::tuple_element_t<0, test_param_t>;
    using scalar_key_t = std::tuple_element_t<1, test_param_t>;

    using simd_value_t = pa::simd_score<scalar_value_t>;
    using simd_key_t = pa::simd_score<scalar_key_t>;

    template <size_t size_v>
    using selector_t = pa::simd_selector<simd_value_t, simd_key_t, size_v>;

    void SetUp() override
    {
    }
};

using test_types = ::testing::Types<
    std::pair<int8_t, uint8_t>
>;

TYPED_TEST_SUITE(simd_selector_avx2_test, test_types);


TYPED_TEST(simd_selector_avx2_test, select_from_9_elements)
{
    if constexpr (pa::detail::max_simd_size == 32) {
        constexpr std::ptrdiff_t element_count_v = 9;

        using selector_t = typename TestFixture::template selector_t<element_count_v>;
        using simd_value_t = TestFixture::simd_value_t;

        EXPECT_EQ(selector_t::elements_per_select, 16ul);

        using address_t = typename selector_t::address_t;

        // Initialize the container with the data from outside.
        std::array<typename TestFixture::scalar_value_t, element_count_v> data{1, 2, 3, 4, 5, 6, 7, 8, 9};

        address_t converted_data = selector_t::load(data);

        // Select the elements from the following indices.
        typename TestFixture::simd_key_t select_keys{
            8, 7, 6, 5, 4, 3, 2, 1,
            0, 1, 6, 4, 8, 7, 0, 1,
            1, 1, 0, 7, 8, 8, 5, 2,
            3, 4, 5, 6, 7, 8, 0, 1
        };

        auto selector = selector_t::selector_for(select_keys);
        __m256i actual = selector(converted_data);

        simd_value_t expected_values{
            9, 8, 7, 6, 5, 4, 3, 2,
            1, 2, 7, 5, 9, 8, 1, 2,
            2, 2, 1, 8, 9, 9, 6, 3,
            4, 5, 6, 7, 8, 9, 1, 2
        };

        using native_simd_t = typename simd_value_t::native_simd_type;

        for (size_t i = 0; i < expected_values.size(); ++i) {
            EXPECT_EQ(static_cast<int>(expected_values[i]),
                      static_cast<int>(reinterpret_cast<native_simd_t const &>(actual)[i])) << "at index" << i;
        }
    } else {
        GTEST_SKIP() << "Test only available for AVX2.";
    }
}

TYPED_TEST(simd_selector_avx2_test, select_from_32_elements_reverse)
{
    if constexpr (pa::detail::max_simd_size == 32) {
        constexpr std::ptrdiff_t element_count_v = 32;

        using selector_t = typename TestFixture::template selector_t<element_count_v>;
        using simd_value_t = TestFixture::simd_value_t;

        EXPECT_EQ(selector_t::elements_per_select, 32ul);

        using address_t = typename selector_t::address_t;

        // Initialize the container with the data from outside.
        std::array<typename TestFixture::scalar_value_t, element_count_v> data{
             1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
            17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

        address_t converted_data = selector_t::load(data);

        // Scenario 1: Select the elements in reverse order.
        typename TestFixture::simd_key_t select_keys{
            31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
            15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1,  0
        };

        auto selector = selector_t::selector_for(select_keys);
        __m256i actual = selector(converted_data);

        simd_value_t expected_values{
            32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,
            16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1
        };

        using native_simd_t = typename simd_value_t::native_simd_type;

        for (size_t i = 0; i < expected_values.size(); ++i) {
            EXPECT_EQ(static_cast<int>(expected_values[i]),
                      static_cast<int>(reinterpret_cast<native_simd_t const &>(actual)[i])) << "at index " << i;
        }

    } else {
        GTEST_SKIP() << "Test only available for AVX2.";
    }
}

TYPED_TEST(simd_selector_avx2_test, select_from_32_elements_interleaved) {
        if constexpr (pa::detail::max_simd_size == 32) {
        constexpr std::ptrdiff_t element_count_v = 32;

        using selector_t = typename TestFixture::template selector_t<element_count_v>;
        using simd_value_t = TestFixture::simd_value_t;

        EXPECT_EQ(selector_t::elements_per_select, 32ul);

        using address_t = typename selector_t::address_t;

        // Initialize the container with the data from outside.
        std::array<typename TestFixture::scalar_value_t, element_count_v> data{
             1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
            17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

        address_t converted_data = selector_t::load(data);

        // Scenario 2: Select the elements from interleaved indices.
        typename TestFixture::simd_key_t select_keys{
            0, 16, 1, 17,  2, 18,  3, 19,  4, 20,  5, 21,  6, 22,  7, 23,
            8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31
        };

        auto selector = selector_t::selector_for(select_keys);
        __m256i actual = selector(converted_data);

        simd_value_t expected_values{
            1, 17,  2, 18,  3, 19,  4, 20,  5, 21,  6, 22,  7, 23,  8, 24,
            9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31, 16, 32
        };

        using native_simd_t = typename simd_value_t::native_simd_type;

        for (size_t i = 0; i < expected_values.size(); ++i) {
            EXPECT_EQ(static_cast<int>(expected_values[i]),
                      static_cast<int>(reinterpret_cast<native_simd_t const &>(actual)[i])) << "at index " << i;
        }
    } else {
        GTEST_SKIP() << "Test only available for AVX2.";
    }
}
