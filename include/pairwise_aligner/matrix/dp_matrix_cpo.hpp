// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides basic CPOs for modelling dp matrix.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <pairwise_aligner/utility/priority_tag.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

namespace dp_matrix {

// ============================================================================
// Concepts for navigating through dp matrix
// ============================================================================

// ----------------------------------------------------------------------------
// CPO: column_count
// ----------------------------------------------------------------------------

namespace _column_count {

constexpr auto column_count(...) noexcept = delete;

struct _fn
{
    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<0> const &) const &
        noexcept(noexcept(column_count(std::declval<object_t &&>())))
        -> decltype(column_count(std::declval<object_t>()))
    {
        return column_count(std::forward<object_t>(object));
    }

    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().column_count()))
        -> decltype(std::declval<object_t>().column_count())
    {
        return std::forward<object_t>(object).column_count();
    }
};
} // namespace _column_count

inline namespace _cpo {

inline constexpr auto column_count = [] (auto && object)
    noexcept(noexcept(std::declval<_column_count::_fn const &>()(std::declval<decltype(object)>(), detail::priority_tag<1>{})))
    -> std::invoke_result_t<_column_count::_fn, decltype(object), detail::priority_tag<1>>
{
    return std::invoke(_column_count::_fn{},
                       std::forward<decltype(object)>(object),
                       detail::priority_tag<1>{});
};
} // namespace _cpo

// ----------------------------------------------------------------------------
// CPO: column_at
// ----------------------------------------------------------------------------

namespace _column_at {

constexpr auto column_at(...) noexcept = delete;

struct _fn
{
    template <typename object_t, typename index_t>
    constexpr auto operator()(object_t && object, index_t && index, detail::priority_tag<0> const &) const &
        noexcept(noexcept(column_at(std::declval<object_t &&>(std::forward<index_t>(index)))))
        -> decltype(column_at(std::declval<object_t>(std::forward<index_t>(index))))
    {
        return column_at(std::forward<object_t>(object), std::forward<index_t>(index));
    }

    template <typename object_t, typename index_t>
    constexpr auto operator()(object_t && object, index_t && index, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().column_at(std::forward<index_t>(index))))
        -> decltype(std::declval<object_t>().column_at(std::forward<index_t>(index)))
    {
        return std::forward<object_t>(object).column_at(std::forward<index_t>(index));
    }
};
} // namespace _column_at

inline namespace _cpo {

inline constexpr auto column_at = [] (auto && object, auto && index)
    noexcept(noexcept(std::declval<_column_at::_fn const &>()(std::declval<decltype(object)>(), index, detail::priority_tag<1>{})))
    -> std::invoke_result_t<_column_at::_fn, decltype(object), decltype(index), detail::priority_tag<1>>
{
    return std::invoke(_column_at::_fn{},
                       std::forward<decltype(object)>(object),
                       std::forward<decltype(index)>(index),
                       detail::priority_tag<1>{});
};
} // namespace _cpo

// ----------------------------------------------------------------------------
// CPO: row_count
// ----------------------------------------------------------------------------

namespace _row_count {

constexpr auto row_count(...) noexcept = delete;

struct _fn
{
    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<0> const &) const &
        noexcept(noexcept(row_count(std::declval<object_t &&>())))
        -> decltype(row_count(std::declval<object_t>()))
    {
        return row_count(std::forward<object_t>(object));
    }

    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().row_count()))
        -> decltype(std::declval<object_t>().row_count())
    {
        return std::forward<object_t>(object).row_count();
    }
};
} // namespace _row_count

inline namespace _cpo {

inline constexpr auto row_count = [] (auto && object)
    noexcept(noexcept(std::declval<_row_count::_fn const &>()(std::declval<decltype(object)>(), detail::priority_tag<1>{})))
    -> std::invoke_result_t<_row_count::_fn, decltype(object), detail::priority_tag<1>>
{
    return std::invoke(_row_count::_fn{}, std::forward<decltype(object)>(object), detail::priority_tag<1>{});
};
} // namespace _cpo


// ----------------------------------------------------------------------------
// CPO: row_at
// ----------------------------------------------------------------------------

namespace _row_at {

constexpr auto row_at(...) noexcept = delete;

struct _fn
{
    template <typename object_t, typename index_t>
    constexpr auto operator()(object_t && object, index_t && index, detail::priority_tag<0> const &) const &
        noexcept(noexcept(row_at(std::declval<object_t &&>(std::forward<index_t>(index)))))
        -> decltype(row_at(std::declval<object_t>(std::forward<index_t>(index))))
    {
        return row_at(std::forward<object_t>(object), std::forward<index_t>(index));
    }

    template <typename object_t, typename index_t>
    constexpr auto operator()(object_t && object, index_t && index, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().row_at(std::forward<index_t>(index))))
        -> decltype(std::declval<object_t>().row_at(std::forward<index_t>(index)))
    {
        return std::forward<object_t>(object).row_at(std::forward<index_t>(index));
    }
};
} // namespace _row_at

inline namespace _cpo {

inline constexpr auto row_at = [] (auto && object, auto && index)
    noexcept(noexcept(std::declval<_row_at::_fn const &>()(std::declval<decltype(object)>(), index, detail::priority_tag<1>{})))
    -> std::invoke_result_t<_row_at::_fn, decltype(object), decltype(index), detail::priority_tag<1>>
{
    return std::invoke(_row_at::_fn{},
                       std::forward<decltype(object)>(object),
                       std::forward<decltype(index)>(index),
                       detail::priority_tag<1>{});
};
} // namespace _cpo

// ============================================================================
// Concepts for accessing dp matrix state
// ============================================================================

// ----------------------------------------------------------------------------
// CPO: dp_column
// ----------------------------------------------------------------------------

namespace _dp_column {

constexpr auto dp_column(...) noexcept = delete;

struct _fn
{
    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<0> const &) const &
        noexcept(noexcept(dp_column(std::declval<object_t &&>())))
        -> decltype(dp_column(std::declval<object_t>()))
    {
        return dp_column(std::forward<object_t>(object));
    }

    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().dp_column()))
        -> decltype(std::declval<object_t>().dp_column())
    {
        return std::forward<object_t>(object).dp_column();
    }
};
} // namespace _dp_column

inline namespace _cpo {

inline constexpr auto dp_column = [] (auto && object)
    noexcept(noexcept(std::declval<_dp_column::_fn const &>()(std::declval<decltype(object)>(), detail::priority_tag<1>{})))
    -> std::invoke_result_t<_dp_column::_fn, decltype(object), detail::priority_tag<1>>
{
    return std::invoke(_dp_column::_fn{},
                       std::forward<decltype(object)>(object),
                       detail::priority_tag<1>{});
};
} // namespace _cpo

// ----------------------------------------------------------------------------
// CPO: dp_row
// ----------------------------------------------------------------------------

namespace _dp_row {

constexpr auto dp_row(...) noexcept = delete;

struct _fn
{
    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<0> const &) const &
        noexcept(noexcept(dp_row(std::declval<object_t &&>())))
        -> decltype(dp_row(std::declval<object_t>()))
    {
        return dp_row(std::forward<object_t>(object));
    }

    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().dp_row()))
        -> decltype(std::declval<object_t>().dp_row())
    {
        return std::forward<object_t>(object).dp_row();
    }
};
} // namespace _dp_row

inline namespace _cpo {

inline constexpr auto dp_row = [] (auto && object)
    noexcept(noexcept(std::declval<_dp_row::_fn const &>()(std::declval<decltype(object)>(), detail::priority_tag<1>{})))
    -> std::invoke_result_t<_dp_row::_fn, decltype(object), detail::priority_tag<1>>
{
    return std::invoke(_dp_row::_fn{},
                       std::forward<decltype(object)>(object),
                       detail::priority_tag<1>{});
};
} // namespace _cpo

// ----------------------------------------------------------------------------
// CPO: column_sequence
// ----------------------------------------------------------------------------

namespace _column_sequence {

constexpr auto column_sequence(...) noexcept = delete;

struct _fn
{
    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<0> const &) const &
        noexcept(noexcept(column_sequence(std::declval<object_t &&>())))
        -> decltype(column_sequence(std::declval<object_t>()))
    {
        return column_sequence(std::forward<object_t>(object));
    }

    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().column_sequence()))
        -> decltype(std::declval<object_t>().column_sequence())
    {
        return std::forward<object_t>(object).column_sequence();
    }
};
} // namespace _column_sequence

inline namespace _cpo {

inline constexpr auto column_sequence = [] (auto && object)
    noexcept(noexcept(std::declval<_column_sequence::_fn const &>()(std::declval<decltype(object)>(), detail::priority_tag<1>{})))
    -> std::invoke_result_t<_column_sequence::_fn, decltype(object), detail::priority_tag<1>>
{
    return std::invoke(_column_sequence::_fn{},
                       std::forward<decltype(object)>(object),
                       detail::priority_tag<1>{});
};
} // namespace _cpo

// ----------------------------------------------------------------------------
// CPO: row_sequence
// ----------------------------------------------------------------------------

namespace _row_sequence {

constexpr auto row_sequence(...) noexcept = delete;

struct _fn
{
    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<0> const &) const &
        noexcept(noexcept(row_sequence(std::declval<object_t &&>())))
        -> decltype(row_sequence(std::declval<object_t>()))
    {
        return row_sequence(std::forward<object_t>(object));
    }

    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().row_sequence()))
        -> decltype(std::declval<object_t>().row_sequence())
    {
        return std::forward<object_t>(object).row_sequence();
    }
};
} // namespace _row_sequence

inline namespace _cpo {

inline constexpr auto row_sequence = [] (auto && object)
    noexcept(noexcept(std::declval<_row_sequence::_fn const &>()(std::declval<decltype(object)>(), detail::priority_tag<1>{})))
    -> std::invoke_result_t<_row_sequence::_fn, decltype(object), detail::priority_tag<1>>
{
    return std::invoke(_row_sequence::_fn{},
                       std::forward<decltype(object)>(object),
                       detail::priority_tag<1>{});
};
} // namespace _cpo

// ----------------------------------------------------------------------------
// CPO: substitution_model
// ----------------------------------------------------------------------------

namespace _substitution_model {

constexpr auto substitution_model(...) noexcept = delete;

struct _fn
{
    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<0> const &) const &
        noexcept(noexcept(substitution_model(std::declval<object_t &&>())))
        -> decltype(substitution_model(std::declval<object_t>()))
    {
        return substitution_model(std::forward<object_t>(object));
    }

    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().substitution_model()))
        -> decltype(std::declval<object_t>().substitution_model())
    {
        return std::forward<object_t>(object).substitution_model();
    }
};
} // namespace _substitution_model

inline namespace _cpo {

inline constexpr auto substitution_model = [] (auto && object)
    noexcept(noexcept(std::declval<_substitution_model::_fn const &>()(std::declval<decltype(object)>(), detail::priority_tag<1>{})))
    -> std::invoke_result_t<_substitution_model::_fn, decltype(object), detail::priority_tag<1>>
{
    return std::invoke(_substitution_model::_fn{},
                       std::forward<decltype(object)>(object),
                       detail::priority_tag<1>{});
};
} // namespace _cpo

// ----------------------------------------------------------------------------
// CPO: tracker
// ----------------------------------------------------------------------------

namespace _tracker {

constexpr auto tracker(...) noexcept = delete;

struct _fn
{
    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<0> const &) const &
        noexcept(noexcept(tracker(std::declval<object_t &&>())))
        -> decltype(tracker(std::declval<object_t>()))
    {
        return tracker(std::forward<object_t>(object));
    }

    template <typename object_t>
    constexpr auto operator()(object_t && object, detail::priority_tag<1> const &) const &
        noexcept(noexcept(std::declval<object_t &&>().tracker()))
        -> decltype(std::declval<object_t>().tracker())
    {
        return std::forward<object_t>(object).tracker();
    }
};
} // namespace _tracker

inline namespace _cpo {

inline constexpr auto tracker = [] (auto && object)
    noexcept(noexcept(std::declval<_tracker::_fn const &>()(std::declval<decltype(object)>(), detail::priority_tag<1>{})))
    -> std::invoke_result_t<_tracker::_fn, decltype(object), detail::priority_tag<1>>
{
    return std::invoke(_tracker::_fn{},
                       std::forward<decltype(object)>(object),
                       detail::priority_tag<1>{});
};
} // namespace _cpo
} // namespace dp_matrix
} // inline namespace v1
} // namespace seqan::pairwise_aligner
