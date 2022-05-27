# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

macro (configure_libunifex)
    enable_testing ()

    # set (gtest_git_tag "main")
    message (STATUS "Fetch libunifex")

    include (FetchContent)
    FetchContent_Declare (
        libunifex_fetch_content
        GIT_REPOSITORY "https://github.com/facebookexperimental/libunifex.git"
        GIT_TAG d7d191e4dc61b67f21cf18220751973025ae520e # last main commit
    )
    FetchContent_MakeAvailable(libunifex_fetch_content)

    # if (NOT TARGET gtest_build)
    #     add_custom_target (gtest_build DEPENDS gtest_main gtest)
    # endif ()
endmacro ()
