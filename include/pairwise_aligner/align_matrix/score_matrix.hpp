// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides base dp matrix.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <type_traits>
#include <utility>

#include <pairwise_aligner/align_matrix/dp_entry_concept.hpp>

namespace align {

template <typename t>
struct underlying_type;

template <typename t>
using underlying_type_t = typename underlying_type<t>::type;

namespace _score_matrix {

template <typename score_t>
struct _entry : public std::pair<score_t, score_t>
{
    using base_t = std::pair<score_t, score_t>;
    using base_t::base_t;
};

template <typename score_t>
struct row_value : public _entry<score_t> {

    using _entry<score_t>::_entry;
    using score_type = score_t;
private:

    constexpr friend score_type & tag_invoke(tag_t<align::diagonal_score>, row_value & me) noexcept {
        return me.first;
    }

    constexpr friend score_type const & tag_invoke(tag_t<align::diagonal_score>, row_value const & me) noexcept {
        return me.first;
    }

    constexpr friend score_type & tag_invoke(tag_t<align::up_score>, row_value & me) noexcept {
        return me.second;
    }

    constexpr friend score_type const & tag_invoke(tag_t<align::up_score>, row_value const & me) noexcept {
        return me.second;
    }
};

template <typename score_t>
struct column_value : public _entry<score_t> {

    using _entry<score_t>::_entry;
    using score_type = score_t;

private:
    constexpr friend score_type & tag_invoke(tag_t<align::current_score>, column_value & me) noexcept {
        return me.first;
    }

    constexpr friend score_type const & tag_invoke(tag_t<align::current_score>, column_value const & me) noexcept {
        return me.first;
    }

    constexpr friend score_type & tag_invoke(tag_t<align::left_score>, column_value & me) noexcept {
        return me.second;
    }

    constexpr friend score_type const & tag_invoke(tag_t<align::left_score>, column_value const & me) noexcept {
        return me.second;
    }
};

template <typename score_t, typename buffer_factory_t>
class _matrix {
private:
    using row_value_t = row_value<score_t>;
    using column_value_t = column_value<score_t>;

    using buffer_t = typename buffer_factory_t::template invoke<row_value_t, column_value_t>;

    buffer_t _buffer{};
public:
    _matrix() = default;

private:

    template <typename cpo_t, typename matrix_t, typename ...args_t>
        requires std::same_as<std::remove_cvref_t<matrix_t>, _matrix>
    constexpr friend auto tag_invoke(cpo_t cpo, matrix_t && me, args_t && ...args)
        noexcept(noexcept(cpo(std::forward<matrix_t>(me)._buffer, std::forward<args_t>(args)...)))
        -> std::invoke_result_t<cpo_t, unifex::member_t<matrix_t, buffer_t>, args_t...> {
        return cpo(std::forward<matrix_t>(me)._buffer, std::forward<args_t>(args)...);
    }
};

namespace _cpo {

template <typename score_t>
struct _fn {
    template <typename buffer_factory_t>
    _matrix<score_t, buffer_factory_t> operator()(buffer_factory_t const &) const noexcept {
        return _matrix<score_t, buffer_factory_t>{};
    }
};

} // namespace _cpo
} // namespace _score_matrix

template <typename score_t>
inline constexpr _score_matrix::_cpo::_fn<score_t> score_matrix{};

template <typename score_t>
struct underlying_type<_score_matrix::_entry<score_t>> :
    public std::type_identity<typename _score_matrix::_entry<score_t>::base_t>
{};

} // namespace align

namespace std {

template <typename score_t>
struct tuple_size<align::_score_matrix::_entry<score_t>> :
    public tuple_size<align::underlying_type_t<align::_score_matrix::_entry<score_t>>>
{};

template <size_t idx, typename score_t>
struct tuple_element<idx, align::_score_matrix::_entry<score_t>> :
    public tuple_element_t<idx, align::underlying_type_t<align::_score_matrix::_entry<score_t>>>
{};
} // namespace std
