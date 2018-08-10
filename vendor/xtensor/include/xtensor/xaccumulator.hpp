/***************************************************************************
* Copyright (c) 2016, Johan Mabille, Sylvain Corlay and Wolf Vollprecht    *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#ifndef XTENSOR_ACCUMULATOR_HPP
#define XTENSOR_ACCUMULATOR_HPP

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <type_traits>

#include "xexpression.hpp"
#include "xstrides.hpp"
#include "xtensor_forward.hpp"

namespace xt
{

#define DEFAULT_STRATEGY_ACCUMULATORS evaluation_strategy::immediate

    /**************
     * accumulate *
     **************/

    template <class ACCUMULATE_FUNC, class INIT_FUNC = xtl::identity>
    struct xaccumulator_functor
        : public std::tuple<ACCUMULATE_FUNC, INIT_FUNC>
    {
        using self_type = xaccumulator_functor<ACCUMULATE_FUNC, INIT_FUNC>;
        using base_type = std::tuple<ACCUMULATE_FUNC, INIT_FUNC>;
        using accumulate_functor_type = ACCUMULATE_FUNC;
        using init_functor_type = INIT_FUNC;

        xaccumulator_functor()
            : base_type()
        {
        }

        template <class RF>
        xaccumulator_functor(RF&& accumulate_func)
            : base_type(std::forward<RF>(accumulate_func), INIT_FUNC())
        {
        }

        template <class RF, class IF>
        xaccumulator_functor(RF&& accumulate_func, IF&& init_func)
            : base_type(std::forward<RF>(accumulate_func), std::forward<IF>(init_func))
        {
        }
    };

    template <class RF>
    auto make_xaccumulator_functor(RF&& accumulate_func)
    {
        using accumulator_type = xaccumulator_functor<std::remove_reference_t<RF>>;
        return accumulator_type(std::forward<RF>(accumulate_func));
    }

    template <class RF, class IF>
    auto make_xaccumulator_functor(RF&& accumulate_func, IF&& init_func)
    {
        using accumulator_type = xaccumulator_functor<std::remove_reference_t<RF>, std::remove_reference_t<IF>>;
        return accumulator_type(std::forward<RF>(accumulate_func), std::forward<IF>(init_func));
    }

    namespace detail
    {
        template <class F, class E, class EVS>
        xarray<typename std::decay_t<E>::value_type> accumulator_impl(F&&, E&&, std::size_t, EVS)
        {
            static_assert(!std::is_same<evaluation_strategy::lazy, EVS>::value, "Lazy accumulators not yet implemented.");
        }

        template <class F, class E, class EVS>
        xarray<typename std::decay_t<E>::value_type> accumulator_impl(F&&, E&&, EVS)
        {
            static_assert(!std::is_same<evaluation_strategy::lazy, EVS>::value, "Lazy accumulators not yet implemented.");
        }

        template <class T, class R>
        struct xaccumulator_return_type
        {
            using type = xarray<R>;
        };

        template <class T, std::size_t N, class R>
        struct xaccumulator_return_type<xtensor<T, N>, R>
        {
            using type = xtensor<R, N>;
        };

        template <class T, class R>
        using xaccumulator_return_type_t = typename xaccumulator_return_type<T, R>::type;

        template <class F, class E>
        inline auto accumulator_init_with_f(F&& f, E& e, std::size_t axis)
        {
            // this function is the equivalent (but hopefully faster) to (if axis == 1)
            // e[:, 0, :, :, ...] = f(e[:, 0, :, :, ...])
            // so that all "first" values are initialized in a first pass

            std::size_t outer_loop_size, inner_loop_size, outer_stride, inner_stride, pos = 0;

            auto set_loop_sizes = [&outer_loop_size, &inner_loop_size](auto first, auto last, std::ptrdiff_t ax) {
                outer_loop_size = std::accumulate(first, first + ax,
                                                  std::size_t(1), std::multiplies<std::size_t>());
                inner_loop_size = std::accumulate(first + ax + 1, last,
                                                  std::size_t(1), std::multiplies<std::size_t>());
            };

            auto set_loop_strides = [&outer_stride, &inner_stride](auto first, auto last, std::ptrdiff_t ax) {
                outer_stride = ax == 0 ? 1 : *std::min_element(first, first + ax);
                inner_stride = (ax == std::distance(first, last) - 1) ? 1 : *std::min_element(first + ax + 1, last);
            };

            set_loop_sizes(e.shape().begin(), e.shape().end(), static_cast<std::ptrdiff_t>(axis));
            set_loop_strides(e.strides().begin(), e.strides().end(), static_cast<std::ptrdiff_t>(axis));

            if (e.layout() == layout_type::column_major)
            {
                // swap for better memory locality (smaller stride in the inner loop)
                std::swap(outer_loop_size, inner_loop_size);
                std::swap(outer_stride, inner_stride);
            }

            for (std::size_t i = 0; i < outer_loop_size; ++i)
            {
                pos = i * outer_stride;
                for (std::size_t j = 0; j < inner_loop_size; ++j)
                {
                    e.storage()[pos] = f(e.storage()[pos]);
                    pos += inner_stride;
                }
            }
        }

        template <class F, class E>
        inline auto accumulator_impl(F&& f, E&& e, std::size_t axis, evaluation_strategy::immediate)
        {
            using accumulate_functor = std::decay_t<decltype(std::get<0>(f))>;
            using function_return_type = typename accumulate_functor::result_type;
            using result_type = xaccumulator_return_type_t<std::decay_t<E>, function_return_type>;

            if (axis >= e.dimension())
            {
                throw std::runtime_error("Axis larger than expression dimension in accumulator.");
            }

            result_type result = e;  // assign + make a copy, we need it anyways

            std::size_t inner_stride = result.strides()[axis];
            std::size_t outer_stride = 1;  // this is either going row- or column-wise (strides.back / strides.front)
            std::size_t outer_loop_size = 0;
            std::size_t inner_loop_size = 0;

            auto set_loop_sizes = [&outer_loop_size, &inner_loop_size](auto first, auto last, std::ptrdiff_t ax) {
                outer_loop_size = std::accumulate(first,
                                                  first + ax,
                                                  std::size_t(1), std::multiplies<std::size_t>());

                inner_loop_size = std::accumulate(first + ax,
                                                  last,
                                                  std::size_t(1), std::multiplies<std::size_t>());
            };

            if (result_type::static_layout == layout_type::row_major)
            {
                set_loop_sizes(result.shape().cbegin(), result.shape().cend(), static_cast<std::ptrdiff_t>(axis));
            }
            else
            {
                set_loop_sizes(result.shape().cbegin(), result.shape().cend(), static_cast<std::ptrdiff_t>(axis + 1));
                std::swap(inner_loop_size, outer_loop_size);
            }

            std::size_t pos = 0;

            inner_loop_size = inner_loop_size - inner_stride;

            // activate the init loop if we have an init function other than identity
            if (!std::is_same<decltype(std::get<1>(f)), xtl::identity>::value)
            {
                accumulator_init_with_f(std::get<1>(f), result, axis);
            }

            pos = 0;
            for (std::size_t i = 0; i < outer_loop_size; ++i)
            {
                for (std::size_t j = 0; j < inner_loop_size; ++j)
                {
                    result.storage()[pos + inner_stride] = std::get<0>(f)(result.storage()[pos],
                                                                       result.storage()[pos + inner_stride]);
                    pos += outer_stride;
                }
                pos += inner_stride;
            }
            return result;
        }

        template <class F, class E>
        inline auto accumulator_impl(F&& f, E&& e, evaluation_strategy::immediate)
        {
            using accumulate_functor = std::decay_t<decltype(std::get<0>(f))>;
            using T = typename accumulate_functor::result_type;

            using result_type = xtensor<T, 1>;
            std::size_t sz = e.size();
            auto result = result_type::from_shape({sz});

            auto it = e.template begin<XTENSOR_DEFAULT_LAYOUT>();

            result.storage()[0] = std::get<1>(f)(*it);
            ++it;

            for (std::size_t idx = 0; it != e.template end<XTENSOR_DEFAULT_LAYOUT>(); ++it)
            {
                result.storage()[idx + 1] = std::get<0>(f)(result.storage()[idx], *it);
                ++idx;
            }
            return result;
        }
    }

    /**
     * Accumulate and flatten array
     * **NOTE** This function is not lazy!
     *
     * @param f functor to use for accumulation
     * @param e xexpression to be accumulated
     * @param evaluation_strategy evaluation strategy of the accumulation
     *
     * @return returns xarray<T> filled with accumulated values
     */
    template <class F, class E, class EVS = DEFAULT_STRATEGY_ACCUMULATORS,
              typename std::enable_if_t<!std::is_integral<EVS>::value, int> = 0>
    inline auto accumulate(F&& f, E&& e, EVS evaluation_strategy = EVS())
    {
        // Note we need to check is_integral above in order to prohibit EVS = int, and not taking the std::size_t
        // overload below!
        return detail::accumulator_impl(std::forward<F>(f), std::forward<E>(e), evaluation_strategy);
    }

    /**
     * Accumulate over axis
     * **NOTE** This function is not lazy!
     *
     * @param f Functor to use for accumulation
     * @param e xexpression to accumulate
     * @param axis Axis to perform accumulation over
     * @param evaluation_strategy evaluation strategy of the accumulation
     *
     * @return returns xarray<T> filled with accumulated values
     */
    template <class F, class E, class EVS = DEFAULT_STRATEGY_ACCUMULATORS>
    inline auto accumulate(F&& f, E&& e, std::size_t axis, EVS evaluation_strategy = EVS())
    {
        return detail::accumulator_impl(std::forward<F>(f), std::forward<E>(e), axis, evaluation_strategy);
    }
}

#endif
