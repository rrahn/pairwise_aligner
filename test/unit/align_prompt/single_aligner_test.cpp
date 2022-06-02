// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>

#include <unifex/just.hpp>
#include <unifex/sync_wait.hpp>

#include <pairwise_aligner/align_configuration/alignment_kind.hpp>
#include <pairwise_aligner/align_configuration/alignment_begin.hpp>
#include <pairwise_aligner/align_configuration/alignment_end.hpp>
#include <pairwise_aligner/align_scoring_scheme/affine_gap_cost.hpp>
#include <pairwise_aligner/align_scoring_scheme/unitary_cost.hpp>
#include <pairwise_aligner/align_configuration/configuration.hpp>
#include <pairwise_aligner/aligner/column_wise.hpp>
#include <pairwise_aligner/aligner/single.hpp>

#include <pairwise_aligner/align_matrix/linear_row_column_matrix.hpp>
#include <pairwise_aligner/align_matrix/score_matrix.hpp>

TEST(basic_sender_test, align_pair)
{
    std::string_view seq1{"ACGTGACTGACACTACGACT"};
    std::string_view seq2{"ACGTGACTGACACTACGACT"};

    align::affine_cost_model<int> cost_model{-10, -1};
    align::unitary_cost_model<int> unitary_cost{4, -5};
    align::configuration cfg{align::alignment_kind{align::alignment_kind_option::score},
                             std::move(cost_model),
                             std::move(unitary_cost),
                             align::alignment_begin_position::first_row_and_column,
                             align::alignment_end_position::last_row_and_column};

    auto gap_model = align::get_gap_cost_model(cfg);
    EXPECT_EQ(align::score(gap_model, 1), -11);
    auto kind = align::get_alignment_kind(cfg);
    EXPECT_TRUE(align::is_option_set(kind, align::alignment_kind_option::score));
    auto substitution_model = align::get_substitution_cost_model(cfg);
    EXPECT_EQ(align::score(substitution_model, 'a', 'a'), 4);
    EXPECT_EQ(align::get_alignment_begin(cfg), align::alignment_begin_position::first_row_and_column);
    EXPECT_EQ(align::get_alignment_end(cfg), align::alignment_end_position::last_row_and_column);


    auto s1 = unifex::just(seq1, seq2) | align::single(cfg) | align::column_wise();
    auto res = unifex::sync_wait(std::move(s1));
    EXPECT_EQ(res->score(), 80);
}

TEST(basic_sender_test, score_matrix)
{
    align::affine_cost_model<int> cost_model{-10, -1};
    auto matrix = align::score_matrix<int>(align::linear_row_column_matrix());
    align::initialise_row(matrix, 10, cost_model);
    align::initialise_column(matrix, 10, cost_model);

    auto e = align::entry_at(matrix, 0, 0);
    EXPECT_EQ(align::current_score(e), 0);

    align::compute(cost_model,
                   align::entry_at(matrix, 0, 0),
                   4);
}
    // Alignment attributes
    // gap_cost
    // align::cfg::alignment_begin{row | column | origin | any}
    // align::cfg::alignment_end{row | column | sink | any}
    // align::cfg::gap_cost{}
    // align::cfg::substitution_cost
    // align::cfg::display{score | position | transcript} -> policy to define how to display the traceback

    // at the moment we just assume everything is given.
    // two states!
    // auto align_request = unifex::just(seq1, seq2) | align::linear_score

    // recursion:
        // depends on gap cost model
        // score() // depends on gap cost model
        // redefine score for row initialisation
        // redefine score for column initialisation

    // how can we get this here?

    // dp_matrix: CPOs
        // initialise_column(, policy)
        // initialise_row(, policy)
        // entry_at(i, j) -> entry_t

    // entry_t:
        // dp_matrix_entry:
        // how to initialise?

    // compute:
        // scout depends on matrix type
        // what if I want to know the coordinates
        // what if I want to know the traceback
        // what if I want to know the score only?


    // initialisation -> row/column => can be adapted how?
    // recursion -> compute
    // linear score matrix + quadratic trace matrix => basic matrix?
    // after that we might want to refine it further

    // instead of saying I want this or that.
    // I am adding these information on top of that.
    // So I want a
        // what is the interface for these?
        // How do I initialise the interleaved matrix?
        // How do I get interleaved trace matrix type? -> can be obtained generically?
        // storage + layout
        // scalar_score_matrix<linear, quadratic>
        // scalar_trace_matrix<quadratic> // linear might or might not exist!
        // scalar_trace_matrix<sampled> // linear might or might not exist!
        // interleaved_score_matrix<linear> // -- needs to implement hirschberg algorithm
        // interleaved_trace_matrix<quadratic_matrix>
        // interleaved_trace_matrix<sampled_matrix>

        // compute matrix, track pointers in middle and continue
        // get the score from the end
        // to get the alignment, the aligner recomputes by splitting the matrix in separate types.
        // but this means the aligner depends on the trace policy
        // so can we get this any different?
        // what if we require a trace?
        // how do we manage this?
        // this is not dependent on the sequences but part of the aligner itself
        // so we can implement an hirschberg_aligner
        // but can we do this for inter_simd as well?
            // interleaved score is foremost a vector of strings
            // hirschberg_simd_aligner etc. ...
            // so what is the trick that makes them combinable?
        // it doesn't matter. If you want this to work, you need some special aligner
        // this might or might not work.
        // we could try to organise this via the scout.
        // then the aligner is responsible for setting the correct scout
        // that does mean the scout is actually tracking the  score etc.
        // can we then have an interleaved hirschberg?
            // if the scout supports it, why not.

        // so that gives us:
            // the score matrix which is the data model that can be wrapped and refined
            // this is req

        // what about finding the optimum in case of local aligner?
        //

        // how do we implement a linear trace matrix

        // linear_score_matrix | quadratic_trace_matrix | interleaved_matrix |

    // so what I want to show is the flexibility of this design.

    // what is the actual aligner and what is the more refined specification?
    // // // aligner is just a handler to some executor.

    // std::cout << res->_seq1 << "\n";

    // using Receiver = unifex::_sync_wait::_receiver<
    //                     align::_single::_request<
    //                         std::basic_string_view<char, std::char_traits<char> >,
    //                         std::basic_string_view<char, std::char_traits<char> >,
    //                         align::configuration<align::alignment_kind>&
    //                     >
    //                 >::type;
    // using Sender = align::_single::_aligner<
    //                     unifex::_just::_sender<
    //                         std::basic_string_view<char, std::char_traits<char> >,
    //                         std::basic_string_view<char, std::char_traits<char> >
    //                     >::type,
    //                     align::configuration<align::alignment_kind>&
    //                 >;

    // static_assert(unifex::tag_invocable<unifex::_connect::_cpo::_fn, Sender, Receiver>);


    // EXPECT_EQ(alignment->_seq1, seq1);



    // auto s2 = unifex::just(std::vector<std::string_view>{16, seq1}, std::vector<std::string_view>{16, seq2})
    //         | align::interleave(align::max_score_range<int16_t>)
    //         | align::direct_column_wise(align_config);  // setting the substitution model

    // auto s3 = unifex::just(std::vector<std::string_view>{32, seq1}, std::vector<std::string_view>{32, seq2})
    //         | align::interleave(align::max_score_range<int8_t>) | align::saturate(align::max_score_range<int16_t>)
    //         | align::split_column_wise | align::join_column_wise(align_config); // does not call set_next anymore



    // we pair the split and the join operation:
        // what about align config?
        // gap cost model
        // scout depends on the configuration method
        // how do I know that the scout is comparable?
        // we have to specify the correct one.
    // auto s2 = unifex::just(seq1[], seq2[]) | align::interleaved | align::striped(unitary_cost_local_simd{});
    // auto s3 = unifex::just(seq1[], seq2[]) | align::interleave
    //                                        | align::saturate
    //                                        | align::column_major | align::striped_join(unitary_cost_global_simd{});

    // column_major sends matrix type to interleave
    // wraps matrix type does something with the score function
    //


    // now how can we model prompts to get from simple direct prompt to pipeline of prompts with
    // some possibilities to modify it?


    // Number one we get a single source!
    // auto aligner = align::one(seq1, seq2)
    //              | align::request(align::cfg::gap_cost_linear{.extension_cost}) // => forms an alignment prompt
    //              | align::prompt()
    //              | align::prompt2()
    //              | align::prompt3()
    //              | align::aligner;

    // prompt connects to aligner ->
        // chain this by setting the prompt of the next aligner
    //


    // designing the pipeline:
