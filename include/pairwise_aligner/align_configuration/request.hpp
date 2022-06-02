// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides request.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/align_configuration/configure_cpo.hpp>
#include <pairwise_aligner/sender/tag_invoke.hpp>

#include <pairwise_aligner/align_configuration/request.hpp>
#include <pairwise_aligner/score_model/score_model_unitary.hpp>
#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/gap_model_affine.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

namespace align {

namespace _get_first_sequence {

inline const struct _fn {
    template<typename object_t>
        requires tag_invocable<_fn, object_t>
    constexpr auto operator()(object_t && object) const
        noexcept(is_nothrow_tag_invocable_v<_fn, object_t>)
        -> tag_invoke_result_t<_fn, object_t> {
        return align::tag_invoke(_fn{}, std::forward<object_t>(object));
    }

} get_first_sequence{};

} // namespace _get_first_sequence

using _get_first_sequence::get_first_sequence;

namespace _get_second_sequence {

inline const struct _fn {
    template<typename object_t>
        requires tag_invocable<_fn, object_t>
    constexpr auto operator()(object_t && object) const
        noexcept(is_nothrow_tag_invocable_v<_fn, object_t>)
        -> tag_invoke_result_t<_fn, object_t> {
        return align::tag_invoke(_fn{}, std::forward<object_t>(object));
    }
} get_second_sequence{};

} // namespace _get_second_sequence

using _get_second_sequence::get_second_sequence;

template <typename sequence1_t, typename sequence2_t, typename align_config_t>
struct request {
    sequence1_t _seq1;
    sequence2_t _seq2;
    align_config_t _config;

    template <typename aligner_t>
    class aligner {
    private:
        aligner_t _aligner;

    public:

        aligner() = delete;
        constexpr aligner(aligner_t aligner) : _aligner{std::move(aligner)}
        {}

        template <typename ...args_t>
        constexpr auto operator()(args_t && ...args) -> decltype(_aligner.compute(std::forward<args_t>(args)...)) {
            return _aligner.compute(std::forward<args_t>(args)...);
        }
    };

private:

    constexpr friend auto tag_invoke(tag_t<align::configure>, request const &) {
        // TODO Make configurable by user!
        auto aligner_config =
            seqan::pairwise_aligner::cfg::method_global(
                seqan::pairwise_aligner::cfg::gap_model_affine(
                    seqan::pairwise_aligner::cfg::score_model_unitary(4, -5),
                    -10, -1
                ),
                seqan::pairwise_aligner::cfg::leading_end_gap{
                        .first_column = seqan::pairwise_aligner::cfg::end_gap::penalised,
                        .first_row = seqan::pairwise_aligner::cfg::end_gap::free},
                seqan::pairwise_aligner::cfg::trailing_end_gap{
                    .last_column = seqan::pairwise_aligner::cfg::end_gap::penalised,
                    .last_row = seqan::pairwise_aligner::cfg::end_gap::free}
            );

        return aligner{seqan::pairwise_aligner::cfg::configure_aligner(aligner_config)};
    }

    constexpr friend sequence1_t & tag_invoke(tag_t<align::get_first_sequence>, request & me) noexcept {
        return me._seq1;
    }

    constexpr friend sequence2_t & tag_invoke(tag_t<align::get_second_sequence>, request & me) noexcept {
        return me._seq2;
    }
};

} // namespace align
