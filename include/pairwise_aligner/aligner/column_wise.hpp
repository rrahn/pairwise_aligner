// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides column wise aligner.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <iostream>

#include <unifex/bind_back.hpp>
#include <unifex/sender_concepts.hpp>

#include <pairwise_aligner/sender/tag_invoke.hpp>
#include <pairwise_aligner/sender/type_traits.hpp>
#include <pairwise_aligner/align_configuration/alignment_begin.hpp>
#include <pairwise_aligner/align_configuration/request.hpp>
#include <pairwise_aligner/align_scoring_scheme/zero_cost.hpp>

#include <pairwise_aligner/align_matrix/linear_row_column_matrix.hpp>
#include <pairwise_aligner/align_matrix/score_matrix.hpp>

namespace align {
namespace _column_wise {

// wrap the score type
// make new score


// TODO: create aligner factory!
template <typename operation_t>
struct _receiver {

    using next_receiver_t = typename operation_t::next_receiver_type;

    operation_t & _operation;

    template <typename request_t>
    void set_value(request_t && request) && noexcept {
        try {
            std::cout << "3. Receive their request!\n";

            // request the score type from the gap cost model
            // TODO: proper metafunction style for this
            // these are basically traits of the request!
            auto const & gap_cost = align::get_gap_cost_model(request._config);
            auto const & substitution_cost = align::get_substitution_cost_model(request._config);
            using gap_cost_model_t = std::remove_cvref_t<decltype(gap_cost)>;
            using score_t = typename gap_cost_model_t::score_type;

            // we have used the bind back option to modify the matrix structure
            // now we need to specify how to compute it
            auto && first_sequence = align::get_first_sequence(request);
            auto && second_sequence = align::get_second_sequence(request);

            size_t const seq1_size = std::ranges::size(first_sequence);
            size_t const seq2_size = std::ranges::size(second_sequence);

            auto matrix = align::score_matrix<score_t>(align::linear_row_column_matrix());

            if (align::get_alignment_begin(request._config) == align::alignment_begin_position::first_row ||
                align::get_alignment_begin(request._config) == align::alignment_begin_position::first_row_or_column ||
                align::get_alignment_begin(request._config) == align::alignment_begin_position::any) {

                align::initialise_row(matrix, seq1_size + 1, align::zero_cost<decltype(gap_cost)>{gap_cost});
            } else {
                align::initialise_row(matrix, seq1_size + 1, gap_cost);
            }

            if (align::get_alignment_begin(request._config) == align::alignment_begin_position::first_column ||
                align::get_alignment_begin(request._config) == align::alignment_begin_position::first_row_or_column ||
                align::get_alignment_begin(request._config) == align::alignment_begin_position::any) {

                align::initialise_column(matrix, seq2_size + 1, align::zero_cost<decltype(gap_cost)>{gap_cost});
            } else {
                align::initialise_column(matrix, seq2_size + 1, gap_cost);
            }

            for (size_t i = 0; i < seq1_size; ++i) {
                for (size_t j = 0; j < seq2_size; ++j) {

                    align::compute(gap_cost,
                                   align::entry_at(matrix, i, j + 1),
                                   align::score(substitution_cost, first_sequence[i], second_sequence[j]));
                }
            }

            std::cout << "the final result = " << align::current_score(align::entry_at(matrix, seq1_size, seq2_size)) << "\n";
            // // why this indirection, just to wrap again how the entry must be computed!
            // for (auto && column : matrix) {
            //     for (auto && entry : column) {
            //         align::compute(entry, gap_cost, substitution_cost); // shall we pass this through?
            //     }
            // }

            // 4. we process the request
            std::cout << "4. Process their request!\n";
            // We get this from the request!
            // out source into aligner factory!
            // allocate memory and
            auto aligner = align::configure(request);
            auto alignment = aligner(align::get_first_sequence(request), align::get_second_sequence(request));
            // std::cout << "seq1 " << request._seq1 << "\n";
            // std::cout << "seq2 " << request._seq2 << "\n";
            // 5. we send the next our result
            std::cout << "5. Send our request to next!\n";
            unifex::set_value(std::move(next_receiver()), std::move(alignment));
        } catch (...) {
            unifex::set_error(std::move(next_receiver()), std::current_exception());
        }
    }

    // pass through in case the error was set
    template <typename error_t>
    void set_error(error_t && error) && noexcept {

        unifex::set_error(std::move(next_receiver()), std::forward<error_t>(error));
    }

    void set_done() && noexcept {
        unifex::set_done(std::move(next_receiver()));
    }

private:

    next_receiver_t & next_receiver() const noexcept {
        return _operation._next_receiver;
    }
};

template <typename predecessor_t, typename next_receiver_t>
struct _operation {

    using next_receiver_type = next_receiver_t;
    using current_receiver_t = _receiver<_operation>;
    // when invoked predecessor sends request to current_receiver
    using previous_operation_t = unifex::connect_result_t<predecessor_t, current_receiver_t>;

    friend current_receiver_t;

private:

    next_receiver_t _next_receiver;
    previous_operation_t _previous_operation; // local to this operation

public:
    // now we have connected the previous operation to this operations.
    _operation(predecessor_t predecessor, next_receiver_t next_receiver) :
        _next_receiver{std::forward<next_receiver_t>(next_receiver)},
        _previous_operation{unifex::connect(std::forward<predecessor_t>(predecessor), current_receiver_t{*this})}
    {}
    // when this operation is started we need to get the values and send it to the next operation
    // so we store the connected result between the predecessor and our own receiver
    // once this operation is started we start the operation that we get from the left sender
    // and then we use it to define the score matrix and then send off to the next problem

    // when we
    // 1. our operation starts
    void start() noexcept {
        std::cout << "1. Our operation starts!\n";
        // We keep the state of this active in our local storage!
        std::cout << "2. Start their operation!\n";
        unifex::start(_previous_operation);
    }
};

template <typename predecessor_t>
struct _aligner {

    predecessor_t _predecessor;

    // TODO: simplify!
    template <typename request_t>
    using alignment_t =
        seqan::pairwise_aligner::type_list<
            std::invoke_result_t<
                tag_invoke_result_t<tag_t<align::configure>, request_t>,
                tag_invoke_result_t<tag_t<align::get_first_sequence>, request_t &>,
                tag_invoke_result_t<tag_t<align::get_second_sequence>, request_t &>
            >
        >;

    template <template <typename ...> typename Variant, template <typename ...> typename Tuple>
    using value_types = seqan::pairwise_aligner::type_list_nested_apply_t<
                            unifex::sender_value_types_t<predecessor_t,
                                                         seqan::pairwise_aligner::type_list,
                                                         alignment_t>,
                            Variant,
                            Tuple>;

    template <template <typename ...> typename Variant>
    using error_types = unifex::sender_error_types_t<predecessor_t, Variant>;

    static constexpr bool sends_done = unifex::sender_traits<predecessor_t>::sends_done;

    // Shortcut definition for the operation.
    template <typename receiver_t>
    using operation_t = _operation<predecessor_t, receiver_t>;

private:
    template <typename aligner_t, typename next_receiver_t>
        requires std::same_as<std::remove_cvref_t<aligner_t>, _aligner> &&
                 unifex::receiver<std::remove_cvref_t<next_receiver_t>>
    friend auto tag_invoke(tag_t<unifex::connect>, aligner_t && me, next_receiver_t && next_receiver)
        ->  operation_t<next_receiver_t> {
            // now we could basically check the conditions of the configuration and check if the aligner can be used.
        return operation_t<next_receiver_t>{std::forward<aligner_t>(me)._predecessor,
                                            std::forward<next_receiver_t>(next_receiver)};
    }
};

namespace _cpo {
struct _fn {
    template <typename aligner_t>
    auto operator()(aligner_t && aligner) const
        noexcept(std::is_nothrow_constructible_v<_column_wise::_aligner<aligner_t>, aligner_t>)
        -> _column_wise::_aligner<aligner_t> {
        return _column_wise::_aligner<aligner_t>{(aligner_t &&) aligner};
    }

    constexpr auto operator()() const
        noexcept(unifex::is_nothrow_callable_v<tag_t<unifex::bind_back>, _fn>)
        -> unifex::bind_back_result_t<_fn> {
        return unifex::bind_back(*this);
    }
};
} // namespace _cpo
} // namespace _column_wise

inline constexpr _column_wise::_cpo::_fn column_wise{};

} // namespace align
