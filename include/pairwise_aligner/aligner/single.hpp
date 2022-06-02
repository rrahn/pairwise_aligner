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

#include <unifex/bind_back.hpp>
#include <unifex/sender_concepts.hpp>

#include <pairwise_aligner/sender/tag_invoke.hpp>
#include <pairwise_aligner/utility/type_list.hpp>

namespace align {

namespace _single {

template <typename next_receiver_t, typename aligner_config_t>
struct _receiver {

    next_receiver_t _next_receiver;
    aligner_config_t _config;

    template <typename sequence1_t, typename sequence2_t>
    void set_value(sequence1_t && seq1, sequence2_t && seq2) && noexcept {
        unifex::set_value(std::move(_next_receiver),
                          request<sequence1_t, sequence2_t, aligner_config_t>{
                                std::forward<sequence1_t>(seq1),
                                std::forward<sequence2_t>(seq2),
                                std::forward<aligner_config_t>(_config)});
    }

    template <typename error_t>
    void set_error(error_t && error) && noexcept {
        unifex::set_error(std::move(_next_receiver), std::forward<error_t>(error));
    }

    void set_done() && noexcept {
        unifex::set_done(std::move(_next_receiver));
    }
};

template <typename sender_t, typename config_t>
struct _aligner {

    sender_t _sender;
    config_t _config;

    template <typename ...args_t>
    using request_t = seqan::pairwise_aligner::type_list<request<args_t..., config_t>>;

    template <template <typename ...> typename Variant, template <typename ...> typename Tuple>
    using value_types = seqan::pairwise_aligner::type_list_nested_apply_t<
                            unifex::sender_value_types_t<sender_t, seqan::pairwise_aligner::type_list, request_t>,
                            Variant,
                            Tuple>;

    template <template <typename ...> typename Variant>
    using error_types = unifex::sender_error_types_t<sender_t, Variant>;

    static constexpr bool sends_done = true;

private:
    template <typename aligner_t, typename receiver_t>
        requires std::same_as<std::remove_cvref_t<aligner_t>, _aligner> &&
                 unifex::receiver<std::remove_cvref_t<receiver_t>> &&
                 unifex::sender_to<sender_t, _receiver<std::remove_cvref_t<receiver_t>, config_t>>
    friend auto tag_invoke(tag_t<unifex::connect>, aligner_t && me, receiver_t && receiver)
        -> unifex::connect_result_t<unifex::member_t<aligner_t, sender_t>,
                                    _receiver<std::remove_cvref_t<receiver_t>, config_t>> {
        return unifex::connect(std::forward<aligner_t>(me)._sender,
                               _receiver<std::remove_cvref_t<receiver_t>, config_t>{
                                    (receiver_t &&) receiver,
                                    std::forward<aligner_t>(me)._config});
    }
};

template <typename source_sender_t, typename config_t>
inline constexpr auto single(source_sender_t && sequence_sender, config_t && config) {
    return _single::_aligner<std::remove_cvref_t<source_sender_t>, config_t>{
                std::forward<source_sender_t>(sequence_sender),
                std::forward<config_t>(config)};
}

namespace _cpo {
  struct _fn {
    template <typename source_sender_t, typename config_t>
    auto operator()(source_sender_t && source, config_t && config) const
        noexcept(std::is_nothrow_constructible_v<
                    _single::_aligner<source_sender_t, config_t>, source_sender_t, config_t>)
        -> _single::_aligner<source_sender_t, config_t> {
      return _single::_aligner<source_sender_t, config_t> {
          (source_sender_t &&) source,
          (config_t &&) config};
    }

    template <typename config_t>
    constexpr auto operator()(config_t && config) const
        noexcept(unifex::is_nothrow_callable_v<tag_t<unifex::bind_back>, _fn, config_t>)
        -> unifex::bind_back_result_t<_fn, config_t> {
      return unifex::bind_back(*this, (config_t &&) config);
    }
  };
} // namespace _cpo
} // namespace _single

inline constexpr _single::_cpo::_fn single {};

} // namespace align
