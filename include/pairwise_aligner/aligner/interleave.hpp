// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides regular alignment prompt.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <type_traits>

// #include <pairwise_aligner/align_matrix/interleave_matrix.hpp>
#include <pairwise_aligner/utility/type_list.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <template <typename ...> typename prev_matrix_t, typename interleaved_sequence1_t, typename interleaved_sequence2_t>
struct interleave_matrix_prompt {
    interleaved_sequence1_t _seq1;
    interleaved_sequence2_t _seq2;

    // will check if works here.
    template <typename base_dp_matrix_t>
    using matrix_type = interleave_matrix<prev_matrix_t<wrapped_matrix_t<base_dp_matrix_t>>>;

private:
    constexpr friend interleaved_sequence1_t & tag_invoke(tag_t<align::first_sequence>,
                                                          interleave_matrix_prompt & me) noexcept {
        return me._seq1; // ref view?
    }

    constexpr friend interleaved_sequence2_t & tag_invoke(tag_t<align::second_sequence>,
                                                          interleave_matrix_prompt & me) noexcept {
        return me._seq2; // ref view
    }
};

template <typename aligner_t>
struct _interleave_aligner {

    aligner_t _aligner;

    template <typename prompt_t>
        // requires some specific input required here.
    constexpr void set_value(prompt_t && prompt) {
        // access the value type
        // transform sequences here with simd specification.
        align::set_value(std::forward<aligner_t>(_aligner),
                         interleave_matrix_prompt{std::forward<sequence1_t>(seq1),
                                                  std::forward<sequence2_t>(seq2)});
    }

    template <typename ...errors_t>
        // requires some specific input required here.
    constexpr void set_error(errors_t && ...errors) {
        align::set_error(std::forward<aligner_t>(_aligner), std::forward<errors_t>(errors)...);
    }

    constexpr void set_done() {
        align::set_done(std::forward<aligner_t>(_aligner));
    }
};

template <typename prompt_t, typename interleave_count>
struct _interleave_prompt {

    prompt_t _prompt;

    // we get the sender types from sender here!
    // we send a prompt
    // we need to get the types from the sender_t
    using _prompt_t = typename std::remove_cvref_t<sender_t>::template
                        value_types<concat_type_lists_t, type_list>::apply<interleave_matrix_prompt>;

    // getting prompt type from previous prompt type etc.

    template <template <typename ...> typename Variant,  template <typename ...> typename Tuple>
    using value_types = Variant<Tuple<_prompt_t>>;

    template <template <typename ...> typename Tuple>
    using error_types = typename std::remove_cvref_t<sender_t>::template error_types<Tuple>;

    inline constexpr bool sends_done = true;

    template <typename aligner_t>
        // requires aligner_for<prompt_t>
    constexpr auto connect_to(aligner_t && aligner) {
        return align::connect_to(std::forward<prompt_t>(_prompt),
                                 _interleave_aligner<aligner_t>(std::forward<aligner_t>(aligner)));
    }
};

template <typename prompt_t, typename interleave_count>
inline constexpr auto interleave(prompt_t && prompt, interleave_count) {
    return _interleave_prompt<prompt_t, interleave_count>{std::forward<prompt_t>(prompt)};
}
} // inline namespace v1
} // namespace seqan::pairwise_aligner
