// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/core/configuration/configuration.hpp>

#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>

#include "alignment_benchmark_fixture.hpp"

namespace aligner::benchmark::fixed_simd {
namespace pa = seqan::pairwise_aligner;

DEFINE_BENCHMARK_VALUES(standard_unitary_same_size,
    .configurator = pa::cfg::gap_model_affine(
                       pa::cfg::score_model_unitary_simd(static_cast<int8_t>(4), static_cast<int8_t>(-5)),
                    -10, -1),
    .seqan_configurator = seqan3::configuration{} | seqan3::align_cfg::method_global{},
    .sequence_size_mean = 150,
    .sequence_size_variance = 0,
    .sequence_count = pa::detail::max_simd_size
)

ALIGNER_BENCHMARK(fixed_simd, standard_unitary_same_size)

} // namespace aligner::benchmark::fixed_simd

BENCHMARK_MAIN();
