// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::score_model_matrix_simd_1xN.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <cassert>
#include <concepts>

#include <seqan3/utility/container/aligned_allocator.hpp>

#include <pairwise_aligner/simd/simd_base.hpp>
#include <pairwise_aligner/simd/simd_rank_selector.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename score_t, size_t dimension>
struct _score_model_matrix_simd_1xN
{
    class type;
};

template <typename score_t, size_t dimension>
using score_model_matrix_simd_1xN = typename _score_model_matrix_simd_1xN<score_t, dimension>::type;

template <typename score_t, size_t dimension>
class _score_model_matrix_simd_1xN<score_t, dimension>::type :
    protected detail::simd_rank_selector_t<simd_score<int8_t, score_t::size_v>>
{
private:

    using index_t = simd_score<int8_t, score_t::size_v>;
    using simd_rank_selector_t = detail::simd_rank_selector_t<index_t>;
    using typename simd_rank_selector_t::rank_map_t;

    struct _interleaved_substitution_profile;

    std::array<rank_map_t, dimension> _matrix{};
    score_t _zero{};

public:

    static constexpr size_t dimension_v = dimension;
    using score_type = score_t;
    using index_type = index_t;
    using offset_type = std::pair<int32_t, index_t>;

    using profile_type = _interleaved_substitution_profile;

    type() = default;

    template <typename substitution_matrix_t> // TODO: does this remain scalar?
    constexpr explicit type(substitution_matrix_t const & matrix, score_t zero = score_t{}) : _zero{zero}
    {
        constexpr size_t chunk_size = (dimension_v - 1 + index_type::size_v) / index_type::size_v;
        for (size_t symbol_rank = 0; symbol_rank < dimension_v; ++symbol_rank) { // we move over the substitution_matrix
            std::array<index_type, chunk_size> tmp;
            for (size_t i = 0; i < dimension_v; ++i) {
                auto [index, offset] = std::pair{i / index_type::size_v, i % index_type::size_v};
                tmp[index][offset] = matrix[symbol_rank][i];
            }
            _matrix[symbol_rank] = simd_rank_selector_t::initialise_rank_map(tmp);
        }
    }

    constexpr score_t const & zero() const noexcept
    {
        return _zero;
    }

    template <typename strip_t>
        requires (std::same_as<std::ranges::range_value_t<strip_t>, index_type>)
    constexpr auto initialise_profile(strip_t && sequence_strip) const noexcept
    {
        return profile_type{_matrix, std::forward<strip_t>(sequence_strip)};
    }

    template <typename value1_t, typename interleaved_profile_t>
    score_type score(score_type const & last_diagonal,
                     value1_t const & value1,
                     interleaved_profile_t const & profile) const noexcept
    {
        // Upcasting the index scores to the score type.
        return last_diagonal + profile[value1];
    }

    // TODO: Refactor into separate factory CPO.
    constexpr type make_substitution_scheme() const noexcept
    {
        return *this;
    }
};

template <typename score_t, size_t dimension>
class _score_model_matrix_simd_1xN<score_t, dimension>::type::_interleaved_substitution_profile
{
    static constexpr size_t alphabet_size_v = dimension;

    using simd_score_t = score_type;
    using scalar_index_t = typename index_t::value_type;

    static constexpr size_t min_align_v = std::max<size_t>(alignof(index_t), 16);

    using interleaved_scores_t = std::vector<index_t, seqan3::aligned_allocator<index_t, min_align_v>>;
    using profile_t = std::vector<interleaved_scores_t,
                                  seqan3::aligned_allocator<interleaved_scores_t, min_align_v>>;

    struct proxy_reference
    {
        profile_t const & _interleaved_profile;
        size_t const _profile_index;

        constexpr auto operator[](size_t const rank) const noexcept
        {
            return _interleaved_profile[rank][_profile_index];
        }
    };

    class iterator;

    profile_t _interleaved_profile;
    size_t _size;

public:

    _interleaved_substitution_profile() = default;
    template <typename matrix_t, typename sequence_slice_t>
    explicit _interleaved_substitution_profile(matrix_t const & matrix, sequence_slice_t && sequence) noexcept
    {
        _size = std::ranges::distance(sequence);

        // Initialise profile: - go over all symbols in range [0..sigma)
        _interleaved_profile.resize(dimension);
        for_each_symbol([&] (scalar_index_t const rank) {
            interleaved_scores_t & profile = _interleaved_profile[rank]; // fill profile for current symbol!
            profile.resize(_size);
            for (size_t index = 0; index < _size; ++index) { // for every symbol in sequence
                profile[index] |= simd_rank_selector_t::select_rank_for(matrix[rank], sequence[index]);
            }
        }, std::make_index_sequence<dimension>());
    }

    iterator begin() const noexcept {
        return iterator{*this, 0};
    }

    iterator end() const noexcept {
        return iterator{*this, size()};
    }

    constexpr auto operator[](size_t const offset) const noexcept
    {
        return proxy_reference{_interleaved_profile, offset};
    }

    constexpr size_t size() const noexcept
    {
        return _size;
    }

    constexpr auto scores_for(scalar_index_t const & rank) const noexcept
    {
        return _interleaved_profile[rank];
    }

private:
    template <typename fn_t, size_t ...alphabet_rank>
    void for_each_symbol(fn_t && fn, std::index_sequence<alphabet_rank...> const &) const noexcept
    {
        (fn(static_cast<scalar_index_t>(alphabet_rank)), ...);
    }
};

template <typename score_t, size_t dimension>
class _score_model_matrix_simd_1xN<score_t, dimension>::type::_interleaved_substitution_profile::iterator
{
    _interleaved_substitution_profile const * _profile;
    size_t _index{0};

public:

    using value_type = proxy_reference;
    using reference = proxy_reference;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using iterator_category = std::random_access_iterator_tag;

    iterator() = default;
    iterator(_interleaved_substitution_profile const & profile, size_t const index) :
        _profile{std::addressof(profile)},
        _index{index}
    {}

    reference operator*() const noexcept {
        assert(_profile != nullptr);
        return _profile->operator[](_index);
    }

    reference operator[](size_t const offset) const noexcept {
        assert(_profile != nullptr);
        return _profile->operator[](_index + offset);
    }

    iterator & operator++() noexcept {
        ++_index;
        return *this;
    }

    iterator operator++(int) noexcept {
        iterator tmp{*this};
        ++_index;
        return tmp;
    }

    iterator & operator--() noexcept {
        --_index;
        return *this;
    }

    iterator operator--(int) noexcept {
        iterator tmp{*this};
        --_index;
        return tmp;
    }

    iterator & operator+=(difference_type const offset) noexcept {
        _index += offset;
        return *this;
    }

    iterator & operator-=(difference_type const offset) noexcept {
        _index -= offset;
        return *this;
    }

    friend constexpr iterator operator+(iterator lhs, difference_type const offset) noexcept {
        lhs += offset;
        return lhs;
    }

    friend constexpr iterator operator+(difference_type const offset, iterator const & rhs) noexcept {
        return rhs + offset;
    }

    friend constexpr iterator operator-(iterator lhs, difference_type const offset) noexcept {
        lhs -= offset;
        return lhs;
    }

    friend constexpr difference_type operator-(iterator const & lhs, iterator const & rhs) noexcept {
        return lhs._index - rhs._index;
    }

    friend constexpr bool operator==(iterator const & lhs, iterator const & rhs) noexcept {
        return lhs._index == rhs._index;
    }

    friend constexpr auto operator<=>(iterator const & lhs, iterator const & rhs) noexcept {
        return lhs._index <=> rhs._index;

    }
};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
