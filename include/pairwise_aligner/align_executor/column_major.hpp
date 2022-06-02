// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides just alignment prompt.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

namespace align
{
inline namespace v1
{
// why do I wrap this into a separate operation state.
// because I need to be able to handle this as a separate operation.
template<typename align_prompt_t, typename receiver_t>
struct column_major_sender_operation {

    align_prompt_t _prompt; // local state generated here.
    receiver_t _receiver; // the one we need to call once we obtained the alignment.

    constexpr void start() noexcept {
        try {
            using base_dp_t = void;
            using dp_matrix_t = typename align_prompt_t::template matrix_type<base_dp_t>;

            // if we get the prompt -> we want to execute it here.
            // the prompt is prepared but not initialised yet.
            auto && seq1 = align::first_sequence(_prompt);
            auto && seq2 = align::second_sequence(_prompt); // transformed already.

            substituion_model substitution_costs; // initialised?
            dp_matrix_t matrix{};
            auto scout = align::initialise(matrix, seq1, seq2); // returns some scout for the final tracking phase.

            for (std::ptrdiff_t i = 0; i < std::ranges::distance(seq1); ++i) {
                // auto && column = column_at(matrix, i); // cached the column and possible initialisations.
                for (std::ptrdiff_t j = 0; j < std::ranges::distance(seq2); ++j) {
                    align::compute(align::entry_at(matrix, i, j),
                                   align::score(substitution_costs, seq1[i], seq2[j]),
                                   scout);
                }
            }

            // scout is based on policies about finding the optimal score.
            // like from regular matrix we need to get a scout that gets the proper value from the list
            // we need to overload it for interleaved matrix for example.
            // we need to overload it again for the tiled dp matrices, etc..
            auto result = align::find_best(scout);
            unifex::set_value(std::move(_receiver), std::move(result));
        } catch(...){
            unifex::set_error(std::move(_receiver), std::current_exception());
        }
    }
};

// convert this to a sender type
template <typename align_prompt_t>
struct column_major_sender {

    align_prompt_t _prompt;

    typename alignment_t = typename align_traits_t<std::remove_cvref_t<align_prompt_t>>::alignment_type;

    // typed sender traits
    template <template <typename ...> Variant, template <typename ...> Tuple>
    using value_types = Variant<Tuple<alignment_t>>;

    template <template <typename ...> Tuple>
    using error_types = Tuple<std::exception_ptr>;

    inline constexpr bool sends_done = false;

    template <typename receiver_t>
    auto connect_to(receiver_t && receiver) {
        using operation_t = column_major_sender_operation<align_prompt_t, receiver_t>;
        return operation_t{std::forward<align_prompt_t>(_prompt),
                           std::forward<receiver_t>(receiver)};
    }
};

template <typename align_prompt_t>
inline constexpr auto column_major(align_prompt_t && prompt) {
    return column_major_sender<align_prompt_t>{std::forward<align_prompt_t>(prompt)};
}

} // inline namespace v1
} // namespace align
