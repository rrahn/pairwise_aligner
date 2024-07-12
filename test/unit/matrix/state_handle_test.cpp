// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <pairwise_aligner/matrix/dp_matrix_state_handle.hpp>

namespace pa = seqan::pairwise_aligner;

struct object {
    int value{42};
};

TEST(state_handle_test, lvalue_state_with_lvalue_ref_members) {

    object obj{};
    auto state = pa::dp_matrix::detail::make_dp_state(obj, obj, obj, obj, obj, obj);
    EXPECT_TRUE((std::same_as<decltype(state),
                 pa::dp_matrix::detail::state_handle<object &, object &, object &, object &, object &, object &>>));

    EXPECT_TRUE((std::same_as<decltype(state.dp_column()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.dp_row()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.column_sequence()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.row_sequence()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.substitution_model()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.tracker()), object &>));

}

TEST(state_handle_test, const_lvalue_state_with_lvalue_ref_members) {
    object obj{};
    auto state = pa::dp_matrix::detail::make_dp_state(obj, obj, obj, obj, obj, obj);
    EXPECT_TRUE((std::same_as<decltype(state),
                 pa::dp_matrix::detail::state_handle<object &, object &, object &, object &, object &, object &>>));

    auto const & c_state = std::as_const(state);
    EXPECT_TRUE((std::same_as<decltype(c_state.dp_column()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.dp_row()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.column_sequence()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.row_sequence()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.substitution_model()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.tracker()), object const &>));
}

TEST(state_handle_test, rvalue_state_with_lvalue_ref_members) {
    object obj{};
    auto state = pa::dp_matrix::detail::make_dp_state(obj, obj, obj, obj, obj, obj);
    EXPECT_TRUE((std::same_as<decltype(state),
                 pa::dp_matrix::detail::state_handle<object &, object &, object &, object &, object &, object &>>));

    EXPECT_TRUE((std::same_as<decltype(std::move(state).dp_column()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).dp_row()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).column_sequence()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).row_sequence()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).substitution_model()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).tracker()), object &>));
}

TEST(state_handle_test, rvalue_const_state_with_lvalue_ref_members) {
    object obj{};
    auto state = pa::dp_matrix::detail::make_dp_state(obj, obj, obj, obj, obj, obj);
    EXPECT_TRUE((std::same_as<decltype(state),
                 pa::dp_matrix::detail::state_handle<object &, object &, object &, object &, object &, object &>>));

    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).dp_column()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).dp_row()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).column_sequence()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).row_sequence()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).substitution_model()), object &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).tracker()), object &>));
}

TEST(state_handle_test, lvalue_state_with_moved_members) {

    auto state = pa::dp_matrix::detail::make_dp_state(object{}, object{}, object{}, object{}, object{}, object{});
    EXPECT_TRUE((std::same_as<decltype(state),
                 pa::dp_matrix::detail::state_handle<object, object, object, object, object, object>>));

    EXPECT_TRUE((std::same_as<decltype(state.dp_column()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.dp_row()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.column_sequence()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.row_sequence()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.substitution_model()), object &>));
    EXPECT_TRUE((std::same_as<decltype(state.tracker()), object &>));

}

TEST(state_handle_test, const_lvalue_state_with_moved_members) {

    auto state = pa::dp_matrix::detail::make_dp_state(object{}, object{}, object{}, object{}, object{}, object{});
    EXPECT_TRUE((std::same_as<decltype(state),
                 pa::dp_matrix::detail::state_handle<object, object, object, object, object, object>>));

    auto const & c_state = std::as_const(state);
    EXPECT_TRUE((std::same_as<decltype(c_state.dp_column()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.dp_row()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.column_sequence()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.row_sequence()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.substitution_model()), object const &>));
    EXPECT_TRUE((std::same_as<decltype(c_state.tracker()), object const &>));
}

TEST(state_handle_test, rvalue_state_with_moved_members) {

    auto state = pa::dp_matrix::detail::make_dp_state(object{}, object{}, object{}, object{}, object{}, object{});
    EXPECT_TRUE((std::same_as<decltype(state),
                 pa::dp_matrix::detail::state_handle<object, object, object, object, object, object>>));

    EXPECT_TRUE((std::same_as<decltype(std::move(state).dp_column()), object &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).dp_row()), object &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).column_sequence()), object &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).row_sequence()), object &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).substitution_model()), object &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(state).tracker()), object &&>));
}

TEST(state_handle_test, rvalue_const_state_with_moved_members) {

    auto state = pa::dp_matrix::detail::make_dp_state(object{}, object{}, object{}, object{}, object{}, object{});
    EXPECT_TRUE((std::same_as<decltype(state),
                 pa::dp_matrix::detail::state_handle<object, object, object, object, object, object>>));

    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).dp_column()), object const &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).dp_row()), object const &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).column_sequence()), object const &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).row_sequence()), object const &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).substitution_model()), object const &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(state)).tracker()), object const &&>));
}
