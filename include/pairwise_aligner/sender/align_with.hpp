// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides aligner::align_with.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <unifex/receiver_concepts.hpp>
#include <unifex/sender_concepts.hpp>

#include <pairwise_aligner/configuration/configure_aligner.hpp>

namespace aligner {
namespace _align_with {
namespace detail {

template <typename aligner_config_t>
using aligner_t = decltype(pa::cfg::configure_aligner(std::declval<aligner_config_t>()));

} // namespace detail

template <typename receiver_t, typename aligner_config_t>
struct _receiver<receiver_t, aligner_config_t>{ class type; };

template <typename receiver_t, typename aligner_config_t>
using receiver = typename _receiver<receiver_t, aligner_config_t>::type;

template <typename receiver_t, typename aligner_config_t>
class _receiver<receiver_t, aligner_config_t>::type
{
    receiver_t _receiver;
    aligner_config_t _aligner_config;

public:

    type(receiver_t receiver, aligner_config_t aligner_config) :
        _receiver(std::move(receiver)),
        _aligner_config(std::forward<aligner_config_t>(aligner_config))
    {}

    template <typename sequence1_t, typename sequence2_t>
    void set_value(sequence1_t && seq1, sequence2_t && seq2) && noexcept {
        try {
            // Try to configure aligner and set value together with sequences.
            unifex::set_value(std::move(_receiver),
                            std::forward<sequence1_t>(seq1),
                            std::forward<sequence2_t>(seq2),
                            pairwise_aligner::configure_aligner(std::forward<aligner_config_t>(aligner_config_)));
        } catch (...) {
        unifex::set_error(std::move(_receiver), std::current_exception());
        }
    }

    template <typename error_t>
    void set_error(error_t && error) && noexcept {
        unifex::set_error(std::move(_receiver), std::forward<error_t>(error));
    }

    void set_done() && noexcept {
        unifex::set_done(std::move(_receiver));
    }
};

template <typename inner_sender_t, typename aligner_config_t>
struct _sender
{
    class type;
};

template <typename inner_sender_t, typename aligner_config_t>
using sender = typename _sender<inner_sender_t, aligner_config_t>::type;

// Sender implementation
template <typename inner_sender_t, typename aligner_config_t>
class _sender<inner_sender_t, aligner_config_t>::type
{
    // We need to decide what the configured alignment algorithm would be!

    // What are the types send by inner_sender?

    // Typed sender traits:
    // What the inner sender passes through and the aligner
    template <template <typename...> typename Variant,
              template <typename...> typename Tuple>;
    using value_types = Variant<Tuple<detail::aligner_t<aligner_config_t>>>;

    // set either exception pointer or whatever the other sender passes through?
    template <template <typename...> typename Variant>
    using error_types = Variant<std::exception_ptr>;

    static constexpr bool sends_done = true;

    // Closure parameters
    inner_sender_t _inner_sender;
    aligner_config_t _aligner_config;

public:
    // Constructor
    type(inner_sender_t inner_sender, aligner_config_t aligner_config) :
        _inner_sender{std::forward<inner_sender_t>(inner_sender)},
        _aligner_config{std::forward<inner_sender_t>(aligner_config)}
    {}

private:
    // Connect
    template <typename this_t, typename receiver_t>
        requires sender_to<inner_sender_t, receiver>
    friend auto unifex::tag_invoke(unifex::tag_t<connect>, this_t && me, receiver_t && receiver)
        // noexcept(?)
        // use standard operation?
        // get inner sender type from this_t!
        -> operation<this_t, >
    {
        return operation(me._inner_sender,
                         receiver<receiver_t, aligner_config_t>(std::move(receiver),
                                                                std::forward<aligner_config_t>(_aligner_config)));
    }
};

// Sender algorithm
struct _fn
{
    template <typename inner_sender_t, typename aligner_config_t>
        // requires ?
    auto operator()(inner_sender_t && inner_sender, aligner_config_t && aligner_config)
        // noexcept ?
    {
        return sender<inner_sender_t, aligner_config_t>{std::forward<inner_sender_t>(inner_sender),
                                                        std::forward<aligner_config_t>(aligner_config)};
    }
};

} // namespace _align_with

inline constexpr _align_with::_fn align_with;
} // namespace aligner
