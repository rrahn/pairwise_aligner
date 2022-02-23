// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------


#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/core/configuration/configuration.hpp>

#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_matrix_simd_NxN.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>

#include <pairwise_aligner/score_model/substitution_matrix.hpp>

#include "alignment_benchmark_fixture.hpp"

namespace aligner::benchmark::fixed_simd {
namespace pa = seqan::pairwise_aligner;

using score_t = int8_t;

inline constexpr size_t max_sequence_count = pa::simd_score<int8_t>::size;
inline constexpr auto base_configurator =
    pa::cfg::gap_model_affine(pa::cfg::score_model_matrix_simd_NxN(pa::blosum62_standard<score_t>), -10, -1);

DEFINE_BENCHMARK_VALUES(standard_matrix_same_size,
    .configurator = base_configurator,
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .alphabet = seqan3::aa20{},
    .sequence_size_mean = aligner::benchmark::sequence_size,
    .sequence_size_variance = 0,
    .sequence_count = max_sequence_count
)

DEFINE_BENCHMARK_VALUES(semi_first_matrix_same_size,
    .configurator = pa::cfg::method_global(base_configurator,
                        pa::cfg::leading_end_gap{.first_column = pa::cfg::end_gap::free },
                        pa::cfg::trailing_end_gap{.last_column = pa::cfg::end_gap::free }),
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .alphabet = seqan3::aa20{},
    .sequence_size_mean = aligner::benchmark::sequence_size,
    .sequence_size_variance = 0,
    .sequence_count = max_sequence_count
)

DEFINE_BENCHMARK_VALUES(semi_second_matrix_same_size,
    .configurator = pa::cfg::method_global(base_configurator,
                        pa::cfg::leading_end_gap{.first_row = pa::cfg::end_gap::free },
                        pa::cfg::trailing_end_gap{.last_row = pa::cfg::end_gap::free }),
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .alphabet = seqan3::aa20{},
    .sequence_size_mean = aligner::benchmark::sequence_size,
    .sequence_size_variance = 0,
    .sequence_count = max_sequence_count
)

DEFINE_BENCHMARK_VALUES(overlap_matrix_same_size,
    .configurator = pa::cfg::method_global(base_configurator,
                        pa::cfg::leading_end_gap{pa::cfg::end_gap::free, pa::cfg::end_gap::free},
                        pa::cfg::trailing_end_gap{pa::cfg::end_gap::free, pa::cfg::end_gap::free}),
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .alphabet = seqan3::aa20{},
    .sequence_size_mean = aligner::benchmark::sequence_size,
    .sequence_size_variance = 0,
    .sequence_count = max_sequence_count
)

ALIGNER_BENCHMARK(fixed_simd, standard_matrix_same_size)
// ALIGNER_BENCHMARK(fixed_simd, semi_first_matrix_same_size)
// ALIGNER_BENCHMARK(fixed_simd, semi_second_matrix_same_size)
// ALIGNER_BENCHMARK(fixed_simd, overlap_matrix_same_size)

} // namespace aligner::benchmark::fixed_simd

BENCHMARK_MAIN();
