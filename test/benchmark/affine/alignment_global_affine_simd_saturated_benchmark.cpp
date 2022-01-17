// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/core/configuration/configuration.hpp>

#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd_saturated.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>

#include "alignment_benchmark_fixture.hpp"

namespace aligner::benchmark::saturated_simd {
namespace pa = seqan::pairwise_aligner;

using score_t = int16_t;

inline constexpr auto base_configurator =
    pa::cfg::gap_model_affine(pa::cfg::score_model_unitary_simd_saturated(static_cast<score_t>(4),
                                                                          static_cast<score_t>(-5)),
                              -10, -1);

DEFINE_BENCHMARK_VALUES(standard_unitary_same_size,
    .configurator = base_configurator,
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .sequence_size_mean = 150,
    .sequence_size_variance = 0,
    .sequence_count = pa::detail::max_simd_size
)

DEFINE_BENCHMARK_VALUES(semi_first_unitary_same_size,
    .configurator = pa::cfg::method_global(base_configurator,
                        pa::cfg::leading_end_gap{.first_column = pa::cfg::end_gap::free },
                        pa::cfg::trailing_end_gap{.last_column = pa::cfg::end_gap::free }),
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .sequence_size_mean = 150,
    .sequence_size_variance = 0,
    .sequence_count = seqan::pairwise_aligner::detail::max_simd_size
)

DEFINE_BENCHMARK_VALUES(semi_second_unitary_same_size,
    .configurator = pa::cfg::method_global(base_configurator,
                        pa::cfg::leading_end_gap{.first_row = pa::cfg::end_gap::free },
                        pa::cfg::trailing_end_gap{.last_row = pa::cfg::end_gap::free }),
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .sequence_size_mean = 150,
    .sequence_size_variance = 0,
    .sequence_count = seqan::pairwise_aligner::detail::max_simd_size
)

DEFINE_BENCHMARK_VALUES(overlap_unitary_same_size,
    .configurator = pa::cfg::method_global(base_configurator,
                        pa::cfg::leading_end_gap{pa::cfg::end_gap::free, pa::cfg::end_gap::free},
                        pa::cfg::trailing_end_gap{pa::cfg::end_gap::free, pa::cfg::end_gap::free}),
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .sequence_size_mean = 150,
    .sequence_size_variance = 0,
    .sequence_count = seqan::pairwise_aligner::detail::max_simd_size
)

ALIGNER_BENCHMARK(saturated_simd, standard_unitary_same_size)
ALIGNER_BENCHMARK(saturated_simd, semi_first_unitary_same_size)
ALIGNER_BENCHMARK(saturated_simd, semi_second_unitary_same_size)
ALIGNER_BENCHMARK(saturated_simd, overlap_unitary_same_size)

} // namespace aligner::benchmark::saturated_simd


BENCHMARK_MAIN();
