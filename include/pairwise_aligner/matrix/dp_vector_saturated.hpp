// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/rrahn/pairwise_aligner/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::pairwise_aligner::dp_vector_saturated.
 * \author Rene Rahn <rahn AT molgen.mpg.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/std/type_traits>

namespace seqan::pairwise_aligner
{
inline namespace v1
{

template <typename dp_vector_t, typename original_cell_t> // int32_t
class dp_vector_saturated
{
private:
    using offset_score_t = typename original_cell_t::score_type;

    dp_vector_t _dp_vector{}; // int8_t
    offset_score_t _offset{}; // int32_t

    template <bool is_const>
    struct _proxy
    {
        using value_type = typename dp_vector_t::value_type;
        using reference = std::conditional_t<is_const,
                                             typename dp_vector_t::const_reference,
                                             typename dp_vector_t::reference>;
        using score_type = offset_score_t;

    private:

        reference _saturated_value;
        offset_score_t const & _offset;

    public:

        _proxy() = delete;
        explicit _proxy(reference saturated_value, offset_score_t const & offset) :
            _saturated_value{saturated_value},
            _offset{offset}
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
        constexpr operator original_cell_t() const noexcept
        {
            original_cell_t cell{_saturated_value};
            std::apply([this] (auto & ...values) { ((values += _offset), ...); }, cell);
            return cell;
        }

        constexpr score_type score() const noexcept
        {
            // std::cout << "_saturated_value[0] = " << (int)_saturated_value.score()[0] << "\n";
            // std::cout << "_offset[0] = " << _offset[0] << "\n";
            // std::cout << "original_cell_t{*this}.score()[0] = " << original_cell_t{*this}.score()[0] << "\n";
            return original_cell_t{*this}.score();
        }

        constexpr typename score_type::value_type score_at(size_t const position) const noexcept
        {
            return _saturated_value.score()[position] + _offset[position];
        }
    };

public:

    using range_type = typename dp_vector_t::range_type;
    using value_type = original_cell_t;
    using reference = _proxy<false>;
    using const_reference = _proxy<true>;

    dp_vector_saturated() = default;

    // return a proxy!
    reference operator[](size_t const pos) noexcept
    {
        return reference{_dp_vector[pos], _offset};
    }

    const_reference operator[](size_t const pos) const noexcept
    {
        return const_reference{_dp_vector[pos], _offset};
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
            return reference{cell, _offset};
        });
    }

    decltype(auto) range() const noexcept
    {
        return _dp_vector.range() | std::views::transform([this] (auto const & cell) {
            return const_reference{cell, _offset};
        });
    }

    constexpr offset_score_t offset() const noexcept
    {
        return _offset;
    }

    constexpr offset_score_t & offset(offset_score_t new_offset) noexcept
    {
        return _offset = std::move(new_offset);
    }

    // initialisation interface
    template <typename predecessor_t>
    struct _factory
    {
        predecessor_t _predecessor;
        offset_score_t & _offset;

        template <typename op_t>
        struct _op
        {
            using small_cell_t = typename dp_vector_t::value_type;
            op_t _op;
            offset_score_t & _offset;
            bool _first_call{true};

            template <typename ...args_t>
            constexpr small_cell_t operator()(args_t && ...args) noexcept
            {
                auto scalar_cell = _op(std::forward<args_t>(args)...);
                if (_first_call)
                {
                    _offset = offset_score_t{scalar_cell.score()};
                    _first_call = false;
                }

                std::apply([this] (auto & ...values) { ((values -= _offset[0]), ...); }, scalar_cell);
                return small_cell_t{scalar_cell}; // construct simd type with relative scores.
            }
        };

        template <typename score_t>
        constexpr auto create() const noexcept
        {
            using score_value_t = typename offset_score_t::value_type;
            using op_t = std::remove_reference_t<decltype(std::declval<predecessor_t>().template
                create<score_value_t>())>;
            return _op<op_t>{_predecessor.template create<score_value_t>(), _offset};
        }
    };

    template <typename sequence_t, typename factory_t>
    auto initialise(sequence_t && sequence, factory_t && init_factory)
    {
        return _dp_vector.initialise(std::forward<sequence_t>(sequence),
                                     _factory<factory_t>{std::forward<factory_t>(init_factory), _offset});
    }
};
} // inline namespace v1
}  // namespace seqan::pairwise_aligner
