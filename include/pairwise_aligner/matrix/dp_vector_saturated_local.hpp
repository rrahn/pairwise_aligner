// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_saturated_local.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/std/type_traits>

#include <pairwise_aligner/simd/simd_base.hpp>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_vector_t, typename regular_cell_t> // int32_t
class dp_vector_saturated_local
{
private:
    using regular_score_t = typename regular_cell_t::score_type;
    using saturated_score_t = simd_score<int8_t, regular_score_t::size>;

    template <bool is_const>
    struct _proxy
    {
        using value_type = typename dp_vector_t::value_type;
        using reference = std::conditional_t<is_const,
                                             typename dp_vector_t::const_reference,
                                             typename dp_vector_t::reference>;
        using score_type = regular_score_t;

    private:

        reference _saturated_value;
        regular_score_t const & _regular_offset;

    public:

        _proxy() = delete;
        explicit _proxy(reference saturated_value, regular_score_t const & offset) :
            _saturated_value{saturated_value},
            _regular_offset{offset}
        {}

        // assignable from actual type.
        _proxy & operator=(_proxy other) noexcept
        {
            using std::swap;
            swap(_saturated_value, other._saturated_value);
            return *this;
        }

        _proxy & operator=(value_type cell) noexcept
        {
            using std::swap;
            swap(_saturated_value, cell);
            return *this;
        }

        // TODO: cast into original cell type
        constexpr operator regular_cell_t() const noexcept
        {
            regular_cell_t cell{_saturated_value};
            std::apply([this] (auto & ...values) { ((values += _regular_offset), ...); }, cell);
            return cell;
        }

        constexpr score_type score() const noexcept
        {
            return score_type{_saturated_value.score()} + _regular_offset;
        }

        constexpr typename score_type::value_type score_at(size_t const position) const noexcept
        {
            return _saturated_value.score()[position] + _regular_offset[position];
        }
    };

    dp_vector_t _dp_vector{}; // int8_t
    regular_score_t _local_zero_offset{};
    regular_score_t _threshold{};
    regular_score_t _regular_offset{}; // int32_t
    regular_score_t _regular_zero_offset{}; // the zero offset in regular score.
    saturated_score_t _saturated_zero_offset{}; // the zero offset in saturated score.

public:

    using range_type = typename dp_vector_t::range_type;
    using value_type = regular_cell_t;
    using reference = _proxy<false>;
    using const_reference = _proxy<true>;

    dp_vector_saturated_local() = default;
    explicit dp_vector_saturated_local(dp_vector_t dp_vector,
                                 int8_t local_zero_offset,
                                 int8_t global_zero_offset,
                                 int8_t threshold) noexcept :
        _dp_vector{std::move(dp_vector)},
        _local_zero_offset{static_cast<typename regular_score_t::value_type>(local_zero_offset * -1)},
        _threshold{static_cast<typename regular_score_t::value_type>(threshold)},
        _regular_offset{static_cast<typename regular_score_t::value_type>(local_zero_offset * -1)},
        _regular_zero_offset{static_cast<typename regular_score_t::value_type>(global_zero_offset)},
        _saturated_zero_offset{global_zero_offset}
    {}

    // return a proxy!
    reference operator[](size_t const pos) noexcept
    {
        return reference{_dp_vector[pos], _regular_offset};
    }

    const_reference operator[](size_t const pos) const noexcept
    {
        return const_reference{_dp_vector[pos], _regular_offset};
    }

    constexpr size_t size() const noexcept
    {
        return _dp_vector.size();
    }

    dp_vector_t & base() noexcept
    {
        return _dp_vector;
    }

    dp_vector_t const & base() const noexcept
    {
        return _dp_vector;
    }

    decltype(auto) range() noexcept
    {
        return _dp_vector.range() | std::views::transform([this] (auto & cell) {
            return reference{cell, _regular_offset};
        });
    }

    decltype(auto) range() const noexcept
    {
        return _dp_vector.range() | std::views::transform([this] (auto const & cell) {
            return const_reference{cell, _regular_offset};
        });
    }

    constexpr saturated_score_t const & saturated_zero_offset() const noexcept
    {
        return _saturated_zero_offset;
    }

    constexpr regular_score_t const & regular_zero_offset() const noexcept
    {
        return _regular_zero_offset;
    }

    constexpr regular_score_t const & local_zero_offset() const noexcept
    {
        return _local_zero_offset;
    }

    constexpr regular_score_t const & threshold() const noexcept
    {
        return _threshold;
    }

    constexpr regular_score_t to_regular_score(saturated_score_t const & score) const noexcept
    {
        // TODO: Update me!
        return _regular_offset + regular_score_t{score};
    }

    constexpr regular_score_t const & offset() const noexcept
    {
        return _regular_offset;
    }

    template <typename mask_t>
    constexpr void update_offset(regular_score_t const & offset, mask_t const & is_local) noexcept
    {
        _regular_offset = blend(is_local,
                                local_zero_offset(),
                                _regular_offset + (offset - _regular_zero_offset));
    }

    // initialisation interface
    template <typename predecessor_t>
    struct _factory
    {
        predecessor_t _predecessor;
        regular_score_t & _regular_offset;

        template <typename op_t>
        struct _op
        {
            using small_cell_t = typename dp_vector_t::value_type;
            op_t _op;
            regular_score_t & _regular_offset;
            bool _first_call{true};

            constexpr small_cell_t operator()(size_t const index) noexcept
            {
                auto scalar_cell = _op(index);
                if (_first_call)
                {
                    _regular_offset += regular_score_t{scalar_cell.score()}; // TODO: synchronise with global alignment!
                    _first_call = false;
                }

                std::apply([this] (auto & ...values) { ((values -= _regular_offset[0]), ...); }, scalar_cell);
                return small_cell_t{scalar_cell}; // construct simd type with relative scores.
            }
        };

        template <typename score_t>
        constexpr auto create() const noexcept
        {
            using score_value_t = typename regular_score_t::value_type;
            using op_t = std::remove_reference_t<decltype(std::declval<predecessor_t>().template
                create<score_value_t>())>;
            return _op<op_t>{_predecessor.template create<score_value_t>(), _regular_offset};
        }
    };

    template <typename sequence_t, typename factory_t>
    auto initialise(sequence_t && sequence, factory_t && init_factory)
    {
        return _dp_vector.initialise(std::forward<sequence_t>(sequence),
                                     _factory<factory_t>{std::forward<factory_t>(init_factory), _regular_offset});
    }
};

namespace detail
{

template <typename regular_cell_t>
struct dp_vector_saturated_local_factory_fn
{

    template <typename dp_vector_t, typename offset_t>
    auto operator()(dp_vector_t && dp_vector,
                    offset_t local_zero_offset,
                    offset_t zero_offset,
                    offset_t threshold) const noexcept
    {
        using pure_dp_vector_t = std::remove_cvref_t<dp_vector_t>;
        return dp_vector_saturated_local<pure_dp_vector_t, regular_cell_t>{
                std::forward<dp_vector_t>(dp_vector),
                std::move(local_zero_offset),
                std::move(zero_offset),
                std::move(threshold)
        };
    }
};

} // namespace detail

template <typename regular_cell_t>
inline constexpr detail::dp_vector_saturated_local_factory_fn<regular_cell_t> dp_vector_saturated_local_factory{};

} // inline namespace v1
}  // namespace seqan::pairwise_aligner
