// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/std/ranges>
#include <string_view>
#include <seqan3/utility/simd/all.hpp>

#include <pairwise_aligner/simd_score_type.hpp>

template <typename score_t>
void simd_compare_and_blend(benchmark::State& state) {
    namespace pa = seqan::pairwise_aligner;

    using simd_type = pa::simd_score<score_t>;
    simd_type a{};
    simd_type b{};

    simd_type t{};
    simd_type f{};

    for (size_t i = 0; i < simd_type::size; ++i)
    {
        a[i] = (std::rand() % (sizeof(score_t) << 3));
        b[i] = (std::rand() % (sizeof(score_t) << 3));
        t[i] = (std::rand() % (sizeof(score_t) << 3));
        f[i] = (std::rand() % (sizeof(score_t) << 3));
    }

    simd_type c{};
    for (auto _ : state) {
        c = blend(compare(a, b, [] (auto const & x, auto const & y) { return (x ^ y).le(simd_type{0}); }), t, f);
    }

    int32_t score{};
    for (size_t i = 0; i < simd_type::size; ++i)
        score += c[i];

    state.counters["score"] = score;
}

// C++11 or newer, you can use the BENCHMARK macro with template parameters:
BENCHMARK_TEMPLATE(simd_compare_and_blend, int32_t);
BENCHMARK_TEMPLATE(simd_compare_and_blend, int16_t);
BENCHMARK_TEMPLATE(simd_compare_and_blend, int8_t);
