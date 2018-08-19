/***************************************************************************
* Copyright (c) 2016, Johan Mabille, Sylvain Corlay and Wolf Vollprecht    *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#ifndef XTENSOR_SEMANTIC_HPP
#define XTENSOR_SEMANTIC_HPP

#include <functional>
#include <utility>

#include "xassign.hpp"
#include "xexpression.hpp"

namespace xt
{
    /**
     * @class xsemantic_base
     * @brief Base interface for assignable xexpressions.
     *
     * The xsemantic_base class defines the interface for assignable
     * xexpressions.
     *
     * @tparam D The derived type, i.e. the inheriting class for which xsemantic_base
     *           provides the interface.
     */
    template <class D>
    class xsemantic_base : public xexpression<D>
    {
    public:

        using base_type = xexpression<D>;
        using derived_type = typename base_type::derived_type;

        using temporary_type = typename xcontainer_inner_types<D>::temporary_type;

        template <class E>
        disable_xexpression<E, derived_type&> operator+=(const E&);

        template <class E>
        disable_xexpression<E, derived_type&> operator-=(const E&);

        template <class E>
        disable_xexpression<E, derived_type&> operator*=(const E&);

        template <class E>
        disable_xexpression<E, derived_type&> operator/=(const E&);

        template <class E>
        disable_xexpression<E, derived_type&> operator%=(const E&);

        template <class E>
        disable_xexpression<E, derived_type&> operator&=(const E&);

        template <class E>
        disable_xexpression<E, derived_type&> operator|=(const E&);

        template <class E>
        disable_xexpression<E, derived_type&> operator^=(const E&);

        template <class E>
        derived_type& operator+=(const xexpression<E>&);

        template <class E>
        derived_type& operator-=(const xexpression<E>&);

        template <class E>
        derived_type& operator*=(const xexpression<E>&);

        template <class E>
        derived_type& operator/=(const xexpression<E>&);

        template <class E>
        derived_type& operator%=(const xexpression<E>&);

        template <class E>
        derived_type& operator&=(const xexpression<E>&);

        template <class E>
        derived_type& operator|=(const xexpression<E>&);

        template <class E>
        derived_type& operator^=(const xexpression<E>&);

        template <class E>
        derived_type& assign(const xexpression<E>&);

        template <class E>
        derived_type& plus_assign(const xexpression<E>&);

        template <class E>
        derived_type& minus_assign(const xexpression<E>&);

        template <class E>
        derived_type& multiplies_assign(const xexpression<E>&);

        template <class E>
        derived_type& divides_assign(const xexpression<E>&);

        template <class E>
        derived_type& modulus_assign(const xexpression<E>&);

    protected:

        xsemantic_base() = default;
        ~xsemantic_base() = default;

        xsemantic_base(const xsemantic_base&) = default;
        xsemantic_base& operator=(const xsemantic_base&) = default;

        xsemantic_base(xsemantic_base&&) = default;
        xsemantic_base& operator=(xsemantic_base&&) = default;

        template <class E>
        derived_type& operator=(const xexpression<E>&);
    };

    /**
     * @class xcontainer_semantic
     * @brief Implementation of the xsemantic_base interface
     * for dense multidimensional containers.
     *
     * The xcontainer_semantic class is an implementation of the
     * xsemantic_base interface for dense multidimensional
     * containers.
     *
     * @tparam D the derived type
     */
    template <class D>
    class xcontainer_semantic : public xsemantic_base<D>
    {
    public:

        using base_type = xsemantic_base<D>;
        using derived_type = D;
        using temporary_type = typename base_type::temporary_type;

        derived_type& assign_temporary(temporary_type&&);

        template <class E>
        derived_type& assign_xexpression(const xexpression<E>& e);

        template <class E>
        derived_type& computed_assign(const xexpression<E>& e);

        template <class E, class F>
        derived_type& scalar_computed_assign(const E& e, F&& f);

    protected:

        xcontainer_semantic() = default;
        ~xcontainer_semantic() = default;

        xcontainer_semantic(const xcontainer_semantic&) = default;
        xcontainer_semantic& operator=(const xcontainer_semantic&) = default;

        xcontainer_semantic(xcontainer_semantic&&) = default;
        xcontainer_semantic& operator=(xcontainer_semantic&&) = default;

        template <class E>
        derived_type& operator=(const xexpression<E>&);
    };

    namespace detail
    {
        template <class E>
        struct has_container_semantics_impl : std::is_base_of<xcontainer_semantic<std::decay_t<E>>, std::decay_t<E>>
        {
        };

        template <class E>
        struct has_container_semantics_impl<xcontainer_semantic<E>> : std::true_type
        {
        };
    }

    template <class E>
    using has_container_semantics = detail::has_container_semantics_impl<E>;

    template <class E, class R = void>
    using enable_xcontainer_semantics = typename std::enable_if<has_container_semantics<E>::value, R>::type;

    template <class E, class R = void>
    using disable_xcontainer_semantics = typename std::enable_if<!has_container_semantics<E>::value, R>::type;

    /**
     * @class xview_semantic
     * @brief Implementation of the xsemantic_base interface for
     * multidimensional views
     *
     * The xview_semantic is an implementation of the xsemantic_base
     * interface for multidimensional views.
     *
     * @tparam D the derived type
     */
    template <class D>
    class xview_semantic : public xsemantic_base<D>
    {
    public:

        using base_type = xsemantic_base<D>;
        using derived_type = D;
        using temporary_type = typename base_type::temporary_type;

        derived_type& assign_temporary(temporary_type&&);

        template <class E>
        derived_type& assign_xexpression(const xexpression<E>& e);

        template <class E>
        derived_type& computed_assign(const xexpression<E>& e);

        template <class E, class F>
        derived_type& scalar_computed_assign(const E& e, F&& f);

    protected:

        xview_semantic() = default;
        ~xview_semantic() = default;

        xview_semantic(const xview_semantic&) = default;
        xview_semantic& operator=(const xview_semantic&) = default;

        xview_semantic(xview_semantic&&) = default;
        xview_semantic& operator=(xview_semantic&&) = default;

        template <class E>
        derived_type& operator=(const xexpression<E>&);
    };

    namespace detail
    {
        template <class E>
        struct has_view_semantics_impl : std::is_base_of<xview_semantic<std::decay_t<E>>, std::decay_t<E>>
        {
        };

        template <class E>
        struct has_view_semantics_impl<xview_semantic<E>> : std::true_type
        {
        };
    }

    template <class E>
    using has_view_semantics = detail::has_view_semantics_impl<E>;

    template <class E, class R = void>
    using enable_xview_semantics = typename std::enable_if<has_view_semantics<E>::value, R>::type;

    template <class E, class R = void>
    using disable_xview_semantics = typename std::enable_if<!has_view_semantics<E>::value, R>::type;

    /*********************************
     * xsemantic_base implementation *
     *********************************/

    /**
     * @name Computed assignement
     */
    //@{
    /**
     * Adds the scalar \c e to \c *this.
     * @param e the scalar to add.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator+=(const E& e) -> disable_xexpression<E, derived_type&>
    {
        return this->derived_cast().scalar_computed_assign(e, std::plus<>());
    }

    /**
     * Subtracts the scalar \c e from \c *this.
     * @param e the scalar to subtract.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator-=(const E& e) -> disable_xexpression<E, derived_type&>
    {
        return this->derived_cast().scalar_computed_assign(e, std::minus<>());
    }

    /**
     * Multiplies \c *this with the scalar \c e.
     * @param e the scalar involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator*=(const E& e) -> disable_xexpression<E, derived_type&>
    {
        return this->derived_cast().scalar_computed_assign(e, std::multiplies<>());
    }

    /**
     * Divides \c *this by the scalar \c e.
     * @param e the scalar involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator/=(const E& e) -> disable_xexpression<E, derived_type&>
    {
        return this->derived_cast().scalar_computed_assign(e, std::divides<>());
    }

    /**
     * Computes the remainder of \c *this after division by the scalar \c e.
     * @param e the scalar involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator%=(const E& e) -> disable_xexpression<E, derived_type&>
    {
        return this->derived_cast().scalar_computed_assign(e, std::modulus<>());
    }

    /**
     * Computes the bitwise and of \c *this and the scalar \c e and assigns it to \c *this.
     * @param e the scalar involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator&=(const E& e) -> disable_xexpression<E, derived_type&>
    {
        return this->derived_cast().scalar_computed_assign(e, std::bit_and<>());
    }

    /**
     * Computes the bitwise or of \c *this and the scalar \c e and assigns it to \c *this.
     * @param e the scalar involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator|=(const E& e) -> disable_xexpression<E, derived_type&>
    {
        return this->derived_cast().scalar_computed_assign(e, std::bit_or<>());
    }

    /**
     * Computes the bitwise xor of \c *this and the scalar \c e and assigns it to \c *this.
     * @param e the scalar involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator^=(const E& e) -> disable_xexpression<E, derived_type&>
    {
        return this->derived_cast().scalar_computed_assign(e, std::bit_xor<>());
    }

    /**
     * Adds the xexpression \c e to \c *this.
     * @param e the xexpression to add.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator+=(const xexpression<E>& e) -> derived_type&
    {
        return operator=(this->derived_cast() + e.derived_cast());
    }

    /**
     * Subtracts the xexpression \c e from \c *this.
     * @param e the xexpression to subtract.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator-=(const xexpression<E>& e) -> derived_type&
    {
        return operator=(this->derived_cast() - e.derived_cast());
    }

    /**
     * Multiplies \c *this with the xexpression \c e.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator*=(const xexpression<E>& e) -> derived_type&
    {
        return operator=(this->derived_cast() * e.derived_cast());
    }

    /**
     * Divides \c *this by the xexpression \c e.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator/=(const xexpression<E>& e) -> derived_type&
    {
        return operator=(this->derived_cast() / e.derived_cast());
    }

    /**
     * Computes the remainder of \c *this after division by the xexpression \c e.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator%=(const xexpression<E>& e) -> derived_type&
    {
        return operator=(this->derived_cast() % e.derived_cast());
    }

    /**
     * Computes the bitwise and of \c *this and the xexpression \c e and assigns it to \c *this.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator&=(const xexpression<E>& e) -> derived_type&
    {
        return operator=(this->derived_cast() & e.derived_cast());
    }

    /**
     * Computes the bitwise or of \c *this and the xexpression \c e and assigns it to \c *this.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator|=(const xexpression<E>& e) -> derived_type&
    {
        return operator=(this->derived_cast() | e.derived_cast());
    }

    /**
     * Computes the bitwise xor of \c *this and the xexpression \c e and assigns it to \c *this.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator^=(const xexpression<E>& e) -> derived_type&
    {
        return operator=(this->derived_cast() ^ e.derived_cast());
    }
    //@}

    /**
     * @name Assign functions
     */
    /**
     * Assigns the xexpression \c e to \c *this. Ensures no temporary
     * will be used to perform the assignment.
     * @param e the xexpression to assign.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::assign(const xexpression<E>& e) -> derived_type&
    {
        return this->derived_cast().assign_xexpression(e);
    }

    /**
     * Adds the xexpression \c e to \c *this. Ensures no temporary
     * will be used to perform the assignment.
     * @param e the xexpression to add.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::plus_assign(const xexpression<E>& e) -> derived_type&
    {
        return this->derived_cast().computed_assign(this->derived_cast() + e.derived_cast());
    }

    /**
     * Subtracts the xexpression \c e to \c *this. Ensures no temporary
     * will be used to perform the assignment.
     * @param e the xexpression to subtract.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::minus_assign(const xexpression<E>& e) -> derived_type&
    {
        return this->derived_cast().computed_assign(this->derived_cast() - e.derived_cast());
    }

    /**
     * Multiplies \c *this with the xexpression \c e. Ensures no temporary
     * will be used to perform the assignment.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::multiplies_assign(const xexpression<E>& e) -> derived_type&
    {
        return this->derived_cast().computed_assign(this->derived_cast() * e.derived_cast());
    }

    /**
     * Divides \c *this by the xexpression \c e. Ensures no temporary
     * will be used to perform the assignment.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::divides_assign(const xexpression<E>& e) -> derived_type&
    {
        return this->derived_cast().computed_assign(this->derived_cast() / e.derived_cast());
    }

    /**
     * Computes the remainder of \c *this after division by the xexpression \c e.
     * Ensures no temporary will be used to perform the assignment.
     * @param e the xexpression involved in the operation.
     * @return a reference to \c *this.
     */
    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::modulus_assign(const xexpression<E>& e) -> derived_type&
    {
        return this->derived_cast().computed_assign(this->derived_cast() % e.derived_cast());
    }

    template <class D>
    template <class E>
    inline auto xsemantic_base<D>::operator=(const xexpression<E>& e) -> derived_type&
    {
        temporary_type tmp(e);
        return this->derived_cast().assign_temporary(std::move(tmp));
    }

    /**************************************
     * xcontainer_semantic implementation *
     **************************************/

    /**
     * Assigns the temporary \c tmp to \c *this.
     * @param tmp the temporary to assign.
     * @return a reference to \c *this.
     */
    template <class D>
    inline auto xcontainer_semantic<D>::assign_temporary(temporary_type&& tmp) -> derived_type&
    {
        return (this->derived_cast() = std::move(tmp));
    }

    template <class D>
    template <class E>
    inline auto xcontainer_semantic<D>::assign_xexpression(const xexpression<E>& e) -> derived_type&
    {
        xt::assign_xexpression(*this, e);
        return this->derived_cast();
    }

    template <class D>
    template <class E>
    inline auto xcontainer_semantic<D>::computed_assign(const xexpression<E>& e) -> derived_type&
    {
        xt::computed_assign(*this, e);
        return this->derived_cast();
    }

    template <class D>
    template <class E, class F>
    inline auto xcontainer_semantic<D>::scalar_computed_assign(const E& e, F&& f) -> derived_type&
    {
        xt::scalar_computed_assign(*this, e, std::forward<F>(f));
        return this->derived_cast();
    }

    template <class D>
    template <class E>
    inline auto xcontainer_semantic<D>::operator=(const xexpression<E>& e) -> derived_type&
    {
        return base_type::operator=(e);
    }

    /*********************************
     * xview_semantic implementation *
     *********************************/

    /**
     * Assigns the temporary \c tmp to \c *this.
     * @param tmp the temporary to assign.
     * @return a reference to \c *this.
     */
    template <class D>
    inline auto xview_semantic<D>::assign_temporary(temporary_type&& tmp) -> derived_type&
    {
        this->derived_cast().assign_temporary_impl(std::move(tmp));
        return this->derived_cast();
    }

    template <class D>
    template <class E>
    inline auto xview_semantic<D>::assign_xexpression(const xexpression<E>& e) -> derived_type&
    {
        xt::assert_compatible_shape(*this, e);
        xt::assign_data(*this, e, false);
        return this->derived_cast();
    }

    template <class D>
    template <class E>
    inline auto xview_semantic<D>::computed_assign(const xexpression<E>& e) -> derived_type&
    {
        xt::assert_compatible_shape(*this, e);
        xt::assign_data(*this, e, false);
        return this->derived_cast();
    }

    template <class D>
    template <class E, class F>
    inline auto xview_semantic<D>::scalar_computed_assign(const E& e, F&& f) -> derived_type&
    {
        D& d = this->derived_cast();
        std::transform(d.begin(), d.end(), d.begin(),
                       [e, &f](const auto& v) { return f(v, e); });
        return this->derived_cast();
    }

    template <class D>
    template <class E>
    inline auto xview_semantic<D>::operator=(const xexpression<E>& e) -> derived_type&
    {
        bool cond = (e.derived_cast().shape().size() == this->derived_cast().dimension()) &&
            std::equal(this->derived_cast().shape().begin(),
                       this->derived_cast().shape().end(),
                       e.derived_cast().shape().begin());
        if (!cond)
        {
            base_type::operator=(broadcast(e.derived_cast(), this->derived_cast().shape()));
        }
        else
        {
            base_type::operator=(e);
        }
        return this->derived_cast();
    }
}

#endif
