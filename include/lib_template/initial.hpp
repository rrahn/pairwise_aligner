// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/lib_template/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides initial file.
 * \author name <name [AT] domain>
 */

#pragma once

#include <iostream>

namespace seqan::lib_template
{
inline namespace v1
{

//!\brief Prints hello!
inline void hello()
{
    std::cout << "Hello!\n";
}
}  // inline namespace v1
}  // namespace seqan::lib_template
