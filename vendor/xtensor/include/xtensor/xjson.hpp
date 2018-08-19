/***************************************************************************
* Copyright (c) 2016, Johan Mabille, Sylvain Corlay and Wolf Vollprecht    *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#ifndef XTENSOR_JSON_HPP
#define XTENSOR_JSON_HPP

#include <cstddef>
#include <stdexcept>
#include <utility>

#include <nlohmann/json.hpp>

#include "xstrided_view.hpp"

namespace xt
{
    /*************************************
     * to_json and from_json declaration *
     *************************************/

    template <class E>
    enable_xexpression<E> to_json(nlohmann::json&, const E&);

    template <class E>
    enable_xcontainer_semantics<E> from_json(const nlohmann::json&, E&);

    /// @cond DOXYGEN_INCLUDE_SFINAE
    template <class E>
    enable_xview_semantics<E> from_json(const nlohmann::json&, E&);
    /// @endcond

    /****************************************
     * to_json and from_json implementation *
     ****************************************/

    namespace detail
    {
        template <class D>
        void to_json_impl(nlohmann::json& j, const xexpression<D>& e, xstrided_slice_vector& slices)
        {
            const auto view = strided_view(e.derived_cast(), slices);
            if (view.dimension() == 0)
            {
                j = view();
            }
            else
            {
                j = nlohmann::json::array();
                using size_type = typename D::size_type;
                size_type nrows = view.shape()[0];
                for (size_type i = 0; i != nrows; ++i)
                {
                    slices.push_back(i);
                    nlohmann::json k;
                    to_json_impl(k, e, slices);
                    j.push_back(std::move(k));
                    slices.pop_back();
                }
            }
        }

        template <class D>
        inline void from_json_impl(const nlohmann::json& j, xexpression<D>& e, xstrided_slice_vector& slices)
        {
            auto view = strided_view(e.derived_cast(), slices);

            if (view.dimension() == 0)
            {
                view() = j;
            }
            else
            {
                using size_type = typename D::size_type;
                size_type nrows = view.shape()[0];
                for (size_type i = 0; i != nrows; ++i)
                {
                    slices.push_back(i);
                    const nlohmann::json& k = j[i];
                    from_json_impl(k, e, slices);
                    slices.pop_back();
                }
            }
        }

        inline unsigned int json_dimension(const nlohmann::json& j)
        {
            if (j.is_array() && j.size())
            {
                return 1 + json_dimension(j[0]);
            }
            else
            {
                return 0;
            }
        }

        template <class S>
        inline void json_shape(const nlohmann::json& j, S& s, std::size_t pos = 0)
        {
            if (j.is_array())
            {
                auto size = j.size();
                s[pos] = size;
                if (size)
                {
                    json_shape(j[0], s, pos + 1);
                }
            }
        }
    }

    /**
     * @brief JSON serialization of an xtensor expression.
     *
     * The to_json method is used by the nlohmann_json package for automatic
     * serialization of user-defined types. The method is picked up by
     * argument-dependent lookup.
     *
     * @param j a JSON object
     * @param e a const \ref xexpression
     */
    template <class E>
    inline enable_xexpression<E> to_json(nlohmann::json& j, const E& e)
    {
        auto sv = xstrided_slice_vector();
        detail::to_json_impl(j, e, sv);
    }

    /**
     * @brief JSON deserialization of a xtensor expression with a container or
     * a view semantics.
     *
     * The from_json method is used by the nlohmann_json library for automatic
     * serialization of user-defined types. The method is picked up by
     * argument-dependent lookup.
     *
     * Note: for converting a JSON object to a value, nlohmann_json requiress
     * the value type to be default constructible, which is typically not the
     * case for expressions with a view semantics. In this case, from_json can
     * be called directly.
     *
     * @param j a const JSON object
     * @param e an \ref xexpression
     */
    template <class E>
    inline enable_xcontainer_semantics<E> from_json(const nlohmann::json& j, E& e)
    {
        auto dimension = detail::json_dimension(j);
        auto s = xtl::make_sequence<typename E::shape_type>(dimension);
        detail::json_shape(j, s);

        // In the case of a container, we resize the container.
        e.resize(s);

        auto sv = xstrided_slice_vector();
        detail::from_json_impl(j, e, sv);
    }

    /// @cond DOXYGEN_INCLUDE_SFINAE
    template <class E>
    inline enable_xview_semantics<E> from_json(const nlohmann::json& j, E& e)
    {
        typename E::shape_type s;
        detail::json_shape(j, s);

        // In the case of a view, we check the size of the container.
        if (!std::equal(s.cbegin(), s.cend(), e.shape().cbegin()))
        {
            throw std::runtime_error("Shape mismatch when deserializing JSON to view");
        }

        auto sv = xstrided_slice_vector();
        detail::from_json_impl(j, e, sv);
    }
    /// @endcond
}

#endif
