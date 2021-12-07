// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::optional_range and utility type traits.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <optional>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <std::ranges::view underlying_range_t>
class optional_range : public std::ranges::view_base
{
private:
    using optional_range_t = std::optional<underlying_range_t>;

    optional_range_t _underlying_range{};

    using iterator_t = std::ranges::iterator_t<underlying_range_t>;
    using sentinel_t = std::ranges::sentinel_t<underlying_range_t>;
    using const_iterator_t = std::ranges::iterator_t<underlying_range_t const>;
    using const_sentinel_t = std::ranges::sentinel_t<underlying_range_t const>;

    static_assert(std::semiregular<iterator_t>);
    static_assert(std::semiregular<sentinel_t>);

public:

    template <bool is_const>
    struct iterator;

    optional_range() = default;
    optional_range(underlying_range_t underlying_range) noexcept : _underlying_range{std::move(underlying_range)}
    {}

    std::iter_reference_t<iterator<false>> operator[](size_t const idx) noexcept
        requires std::ranges::random_access_range<underlying_range_t>
    {
        return *std::ranges::next(begin(), idx);
    }

    std::iter_reference_t<iterator<true>> operator[](size_t const idx) const noexcept
        requires std::ranges::random_access_range<underlying_range_t>
    {
        return *std::ranges::next(begin(), idx);
    }

    constexpr iterator<false> begin() noexcept
    {
        return iterator<false>{_underlying_range, false};
    }

    constexpr iterator<true> begin() const noexcept
    {
        return iterator<true>{_underlying_range, false};
    }

    constexpr iterator<false> end() noexcept
    {
        return iterator<false>{_underlying_range, true};
    }

    constexpr iterator<true> end() const noexcept
    {
        return iterator<true>{_underlying_range, true};
    }
};

template <std::ranges::view underlying_range_t>
template <bool is_const>
class optional_range<underlying_range_t>::iterator
{
    using maybe_const_optional_t = std::conditional_t<is_const, optional_range_t const, optional_range_t>;
    using maybe_const_range_t = std::conditional_t<is_const, underlying_range_t const, underlying_range_t>;
    using underlying_iter_t = std::ranges::iterator_t<maybe_const_range_t>;

    underlying_iter_t _iter{};
    bool _is_empty{false};

public:

    using value_type = std::iter_value_t<underlying_iter_t>;
    using reference = std::iter_reference_t<underlying_iter_t>;
    using pointer = typename std::iterator_traits<underlying_iter_t>::pointer;
    using difference_type = std::iter_difference_t<underlying_iter_t>;
    using iterator_category = typename std::iterator_traits<underlying_iter_t>::iterator_category;

    iterator() = default;
    explicit iterator(maybe_const_optional_t & optional_range, bool at_end)
    {
        _is_empty = !optional_range.has_value();
        if (!_is_empty)
        {
            if (at_end)
                _iter = std::ranges::end(optional_range.value());
            else
                _iter = std::ranges::begin(optional_range.value());
        }
    }

    constexpr reference operator*() const noexcept
    {
        return *_iter;
    }

    constexpr reference operator[](difference_type const offset) const noexcept
        requires std::random_access_iterator<underlying_iter_t>
    {
        return *(_iter + offset);
    }

    constexpr iterator & operator++() noexcept
    {
        ++_iter;
        return *this;
    }

    constexpr auto operator++(int) noexcept
    {
        if constexpr (std::forward_iterator<underlying_iter_t>)
        {
            iterator tmp{*this};
            ++_iter;
            return tmp;
        }
        else
        {
            ++_iter;
        }
    }

    constexpr iterator & operator+=(difference_type const offset) noexcept
        requires std::random_access_iterator<underlying_iter_t>
    {
        _iter += offset;
        return *this;
    }

    constexpr iterator operator+(difference_type const offset) const noexcept
        requires std::random_access_iterator<underlying_iter_t>
    {
        iterator tmp{*this};
        tmp += offset;
        return tmp;
    }

    constexpr friend iterator operator+(difference_type const offset, iterator const & iter) noexcept
        requires std::random_access_iterator<underlying_iter_t>
    {
        return iter + offset;
    }

    constexpr iterator & operator--() noexcept
        requires std::bidirectional_iterator<underlying_iter_t>
    {
        --_iter;
        return *this;
    }

    constexpr auto operator--(int) noexcept
        requires std::bidirectional_iterator<underlying_iter_t>
    {
        iterator tmp{*this};
        --_iter;
        return tmp;
    }

    constexpr iterator & operator-=(difference_type const offset) noexcept
        requires std::random_access_iterator<underlying_iter_t>
    {
        _iter -= offset;
        return *this;
    }

    constexpr iterator operator-(difference_type const offset) const noexcept
        requires std::random_access_iterator<underlying_iter_t>
    {
        iterator tmp{*this};
        tmp -= offset;
        return tmp;
    }

    constexpr iterator operator-(iterator const & rhs) const noexcept
        requires std::sized_sentinel_for<underlying_iter_t, underlying_iter_t>
    {
        return _iter - rhs._iter;
    }

    constexpr bool operator==(iterator const & rhs) const noexcept
    {
        return _is_empty || (_iter == rhs._iter);
    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
