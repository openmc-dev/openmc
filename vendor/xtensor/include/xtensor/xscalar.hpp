/***************************************************************************
* Copyright (c) 2016, Johan Mabille, Sylvain Corlay and Wolf Vollprecht    *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#ifndef XTENSOR_SCALAR_HPP
#define XTENSOR_SCALAR_HPP

#include <array>
#include <cstddef>
#include <utility>

#include <xtl/xtype_traits.hpp>

#include "xexpression.hpp"
#include "xiterable.hpp"
#include "xlayout.hpp"
#include "xtensor_simd.hpp"

namespace xt
{

    /***********
     * xscalar *
     ***********/

    // xscalar is a cheap wrapper for a scalar value as an xexpression.

    template <bool is_const, class CT>
    class xscalar_stepper;

    template <bool is_const, class CT>
    class xdummy_iterator;

    template <class CT>
    class xscalar;

    template <class CT>
    struct xiterable_inner_types<xscalar<CT>>
    {
        using value_type = std::decay_t<CT>;
        using inner_shape_type = std::array<std::size_t, 0>;
        using const_stepper = xscalar_stepper<true, CT>;
        using stepper = xscalar_stepper<false, CT>;
    };

#define DL XTENSOR_DEFAULT_LAYOUT
    template <class CT>
    class xscalar : public xexpression<xscalar<CT>>,
                    private xiterable<xscalar<CT>>
    {
    public:

        using self_type = xscalar<CT>;

        using value_type = std::decay_t<CT>;
        using reference = value_type&;
        using const_reference = const value_type&;
        using pointer = value_type*;
        using const_pointer = const value_type*;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using simd_value_type = xsimd::simd_type<value_type>;

        using iterable_base = xiterable<self_type>;
        using inner_shape_type = typename iterable_base::inner_shape_type;
        using shape_type = inner_shape_type;
        using expression_tag = xscalar_expression_tag;

        using stepper = typename iterable_base::stepper;
        using const_stepper = typename iterable_base::const_stepper;

        template <layout_type L>
        using layout_iterator = typename iterable_base::template layout_iterator<L>;
        template <layout_type L>
        using const_layout_iterator = typename iterable_base::template const_layout_iterator<L>;

        template <layout_type L>
        using reverse_layout_iterator = typename iterable_base::template reverse_layout_iterator<L>;
        template <layout_type L>
        using const_reverse_layout_iterator = typename iterable_base::template const_reverse_layout_iterator<L>;

        template <class S, layout_type L>
        using broadcast_iterator = typename iterable_base::template broadcast_iterator<S, L>;
        template <class S, layout_type L>
        using const_broadcast_iterator = typename iterable_base::template const_broadcast_iterator<S, L>;

        template <class S, layout_type L>
        using reverse_broadcast_iterator = typename iterable_base::template reverse_broadcast_iterator<S, L>;
        template <class S, layout_type L>
        using const_reverse_broadcast_iterator = typename iterable_base::template const_reverse_broadcast_iterator<S, L>;

        using iterator = value_type*;
        using const_iterator = const value_type*;
        using reverse_iterator = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

        using dummy_iterator = xdummy_iterator<false, CT>;
        using const_dummy_iterator = xdummy_iterator<true, CT>;

        static constexpr layout_type static_layout = layout_type::any;
        static constexpr bool contiguous_layout = true;

        xscalar() noexcept;
        xscalar(CT value) noexcept;

        operator value_type&() noexcept;
        operator const value_type&() const noexcept;

        size_type size() const noexcept;
        size_type dimension() const noexcept;
        const shape_type& shape() const noexcept;
        layout_type layout() const noexcept;

        template <class... Args>
        reference operator()(Args...) noexcept;
        template <class... Args>
        reference at(Args...);
        template <class... Args>
        reference unchecked(Args...) noexcept;
        template <class S>
        disable_integral_t<S, reference> operator[](const S&) noexcept;
        template <class I>
        reference operator[](std::initializer_list<I>) noexcept;
        reference operator[](size_type) noexcept;

        template <class... Args>
        const_reference operator()(Args...) const noexcept;
        template <class... Args>
        const_reference at(Args...) const;
        template <class... Args>
        const_reference unchecked(Args...) const noexcept;
        template <class S>
        disable_integral_t<S, const_reference> operator[](const S&) const noexcept;
        template <class I>
        const_reference operator[](std::initializer_list<I>) const noexcept;
        const_reference operator[](size_type) const noexcept;

        template <class It>
        reference element(It, It) noexcept;

        template <class It>
        const_reference element(It, It) const noexcept;

        template <class S>
        bool broadcast_shape(S& shape, bool reuse_cache = false) const noexcept;

        template <class S>
        bool is_trivial_broadcast(const S& strides) const noexcept;

        template <layout_type L = DL>
        iterator begin() noexcept;
        template <layout_type L = DL>
        iterator end() noexcept;

        template <layout_type L = DL>
        const_iterator begin() const noexcept;
        template <layout_type L = DL>
        const_iterator end() const noexcept;
        template <layout_type L = DL>
        const_iterator cbegin() const noexcept;
        template <layout_type L = DL>
        const_iterator cend() const noexcept;

        template <layout_type L = DL>
        reverse_iterator rbegin() noexcept;
        template <layout_type L = DL>
        reverse_iterator rend() noexcept;

        template <layout_type L = DL>
        const_reverse_iterator rbegin() const noexcept;
        template <layout_type L = DL>
        const_reverse_iterator rend() const noexcept;
        template <layout_type L = DL>
        const_reverse_iterator crbegin() const noexcept;
        template <layout_type L = DL>
        const_reverse_iterator crend() const noexcept;

        template <class S, layout_type L = DL>
        broadcast_iterator<S, L> begin(const S& shape) noexcept;
        template <class S, layout_type L = DL>
        broadcast_iterator<S, L> end(const S& shape) noexcept;

        template <class S, layout_type L = DL>
        const_broadcast_iterator<S, L> begin(const S& shape) const noexcept;
        template <class S, layout_type L = DL>
        const_broadcast_iterator<S, L> end(const S& shape) const noexcept;
        template <class S, layout_type L = DL>
        const_broadcast_iterator<S, L> cbegin(const S& shape) const noexcept;
        template <class S, layout_type L = DL>
        const_broadcast_iterator<S, L> cend(const S& shape) const noexcept;


        template <class S, layout_type L = DL>
        reverse_broadcast_iterator<S, L> rbegin(const S& shape) noexcept;
        template <class S, layout_type L = DL>
        reverse_broadcast_iterator<S, L> rend(const S& shape) noexcept;

        template <class S, layout_type L = DL>
        const_reverse_broadcast_iterator<S, L> rbegin(const S& shape) const noexcept;
        template <class S, layout_type L = DL>
        const_reverse_broadcast_iterator<S, L> rend(const S& shape) const noexcept;
        template <class S, layout_type L = DL>
        const_reverse_broadcast_iterator<S, L> crbegin(const S& shape) const noexcept;
        template <class S, layout_type L = DL>
        const_reverse_broadcast_iterator<S, L> crend(const S& shape) const noexcept;

        template <layout_type L = DL>
        iterator storage_begin() noexcept;
        template <layout_type L = DL>
        iterator storage_end() noexcept;

        template <layout_type L = DL>
        const_iterator storage_begin() const noexcept;
        template <layout_type L = DL>
        const_iterator storage_end() const noexcept;
        template <layout_type L = DL>
        const_iterator storage_cbegin() const noexcept;
        template <layout_type L = DL>
        const_iterator storage_cend() const noexcept;

        template <layout_type L = DL>
        reverse_iterator storage_rbegin() noexcept;
        template <layout_type L = DL>
        reverse_iterator storage_rend() noexcept;

        template <layout_type L = DL>
        const_reverse_iterator storage_rbegin() const noexcept;
        template <layout_type L = DL>
        const_reverse_iterator storage_rend() const noexcept;
        template <layout_type L = DL>
        const_reverse_iterator storage_crbegin() const noexcept;
        template <layout_type L = DL>
        const_reverse_iterator storage_crend() const noexcept;

        template <class S>
        stepper stepper_begin(const S& shape) noexcept;
        template <class S>
        stepper stepper_end(const S& shape, layout_type l) noexcept;

        template <class S>
        const_stepper stepper_begin(const S& shape) const noexcept;
        template <class S>
        const_stepper stepper_end(const S& shape, layout_type l) const noexcept;

        dummy_iterator dummy_begin() noexcept;
        dummy_iterator dummy_end() noexcept;

        const_dummy_iterator dummy_begin() const noexcept;
        const_dummy_iterator dummy_end() const noexcept;

        reference data_element(size_type i) noexcept;
        const_reference data_element(size_type i) const noexcept;

        template <class align, class simd = simd_value_type>
        void store_simd(size_type i, const simd& e);
        template <class align, class simd = simd_value_type>
        simd load_simd(size_type i) const;

    private:

        CT m_value;

        friend class xconst_iterable<self_type>;
        friend class xiterable<self_type>;
    };
#undef DL

    namespace detail
    {
        template <class E>
        struct is_xscalar_impl : std::false_type
        {
        };

        template <class E>
        struct is_xscalar_impl<xscalar<E>> : std::true_type
        {
        };
    }

    template <class E>
    using is_xscalar = detail::is_xscalar_impl<E>;

    namespace detail
    {
        template <class... E>
        struct all_xscalar
        {
            static constexpr bool value = xtl::conjunction<is_xscalar<std::decay_t<E>>...>::value;
        };
    }

    // Note: MSVC bug workaround. Cannot just define
    // template <class... E>
    // using all_xscalar = xtl::conjunction<is_xscalar<std::decay_t<E>>...>;

    template <class... E>
    using all_xscalar = detail::all_xscalar<E...>;

    /******************
     * xref and xcref *
     ******************/

    template <class T>
    xscalar<T&> xref(T& t);

    template <class T>
    xscalar<const T&> xcref(T& t);

    /*******************
     * xscalar_stepper *
     *******************/

    template <bool is_const, class CT>
    class xscalar_stepper
    {
    public:

        using self_type = xscalar_stepper<is_const, CT>;
        using storage_type = std::conditional_t<is_const,
                                                const xscalar<CT>,
                                                xscalar<CT>>;

        using value_type = typename storage_type::value_type;
        using reference = std::conditional_t<is_const,
                                             typename storage_type::const_reference,
                                             typename storage_type::reference>;
        using pointer = std::conditional_t<is_const,
                                           typename storage_type::const_pointer,
                                           typename storage_type::pointer>;
        using size_type = typename storage_type::size_type;
        using difference_type = typename storage_type::difference_type;

        xscalar_stepper(storage_type* c) noexcept;

        reference operator*() const noexcept;

        void step(size_type dim, size_type n = 1) noexcept;
        void step_back(size_type dim, size_type n = 1) noexcept;
        void reset(size_type dim) noexcept;
        void reset_back(size_type dim) noexcept;

        void to_begin() noexcept;
        void to_end(layout_type l) noexcept;

        template <class R>
        R step_simd();

        value_type step_leading();

    private:

        storage_type* p_c;
    };

    /*******************
     * xdummy_iterator *
     *******************/

    namespace detail
    {
        template <bool is_const, class CT>
        using dummy_reference_t = std::conditional_t<is_const,
                                                     typename xscalar<CT>::const_reference,
                                                     typename xscalar<CT>::reference>;

        template <bool is_const, class CT>
        using dummy_pointer_t = std::conditional_t<is_const,
                                                   typename xscalar<CT>::const_pointer,
                                                   typename xscalar<CT>::pointer>;
    }

    template <bool is_const, class CT>
    class xdummy_iterator
        : public xtl::xrandom_access_iterator_base<xdummy_iterator<is_const, CT>,
                                                   typename xscalar<CT>::value_type,
                                                   typename xscalar<CT>::difference_type,
                                                   detail::dummy_pointer_t<is_const, CT>,
                                                   detail::dummy_reference_t<is_const, CT>>
    {
    public:

        using self_type = xdummy_iterator<is_const, CT>;
        using storage_type = std::conditional_t<is_const,
                                                const xscalar<CT>,
                                                xscalar<CT>>;

        using value_type = typename storage_type::value_type;
        using reference = detail::dummy_reference_t<is_const, CT>;
        using pointer = detail::dummy_pointer_t<is_const, CT>;
        using difference_type = typename storage_type::difference_type;
        using iterator_category = std::random_access_iterator_tag;

        explicit xdummy_iterator(storage_type* c) noexcept;

        self_type& operator++() noexcept;
        self_type& operator--() noexcept;

        self_type& operator+=(difference_type n) noexcept;
        self_type& operator-=(difference_type n) noexcept;

        difference_type operator-(const self_type& rhs) const noexcept;

        reference operator*() const noexcept;

        bool equal(const self_type& rhs) const noexcept;
        bool less_than(const self_type& rhs) const noexcept;

    private:

        storage_type* p_c;
    };

    template <bool is_const, class CT>
    bool operator==(const xdummy_iterator<is_const, CT>& lhs,
                    const xdummy_iterator<is_const, CT>& rhs) noexcept;

    template <bool is_const, class CT>
    bool operator<(const xdummy_iterator<is_const, CT>& lhs,
                   const xdummy_iterator<is_const, CT>& rhs) noexcept;

    /*******************************
     * trivial_begin / trivial_end *
     *******************************/

    namespace detail
    {
        template <class CT>
        constexpr auto trivial_begin(xscalar<CT>& c) noexcept -> decltype(c.dummy_begin())
        {
            return c.dummy_begin();
        }

        template <class CT>
        constexpr auto trivial_end(xscalar<CT>& c) noexcept -> decltype(c.dummy_end())
        {
            return c.dummy_end();
        }

        template <class CT>
        constexpr auto trivial_begin(const xscalar<CT>& c) noexcept -> decltype(c.dummy_begin())
        {
            return c.dummy_begin();
        }

        template <class CT>
        constexpr auto trivial_end(const xscalar<CT>& c) noexcept -> decltype(c.dummy_end())
        {
            return c.dummy_end();
        }
    }

    /**************************
     * xscalar implementation *
     **************************/

    // This constructor will not compile when CT is a reference type.
    template <class CT>
    inline xscalar<CT>::xscalar() noexcept
        : m_value()
    {
    }

    template <class CT>
    inline xscalar<CT>::xscalar(CT value) noexcept
        : m_value(value)
    {
    }

    template <class CT>
    inline xscalar<CT>::operator value_type&() noexcept
    {
        return m_value;
    }

    template <class CT>
    inline xscalar<CT>::operator const value_type&() const noexcept
    {
        return m_value;
    }

    template <class CT>
    inline auto xscalar<CT>::size() const noexcept -> size_type
    {
        return 1;
    }

    template <class CT>
    inline auto xscalar<CT>::dimension() const noexcept -> size_type
    {
        return 0;
    }

    template <class CT>
    inline auto xscalar<CT>::shape() const noexcept -> const shape_type&
    {
        static std::array<size_type, 0> zero_shape;
        return zero_shape;
    }

    template <class CT>
    inline layout_type xscalar<CT>::layout() const noexcept
    {
        return static_layout;
    }

    template <class CT>
    template <class... Args>
    inline auto xscalar<CT>::operator()(Args...) noexcept -> reference
    {
        XTENSOR_CHECK_DIMENSION((std::array<int, 0>()), Args()...);
        return m_value;
    }

    template <class CT>
    template <class... Args>
    inline auto xscalar<CT>::at(Args... args) -> reference
    {
        check_dimension(shape(), args...);
        return this->operator()(args...);
    }

    template <class CT>
    template <class... Args>
    inline auto xscalar<CT>::unchecked(Args...) noexcept -> reference
    {
        return m_value;
    }

    template <class CT>
    template <class S>
    inline auto xscalar<CT>::operator[](const S&) noexcept
        -> disable_integral_t<S, reference>
    {
        return m_value;
    }

    template <class CT>
    template <class I>
    inline auto xscalar<CT>::operator[](std::initializer_list<I>) noexcept
        -> reference
    {
        return m_value;
    }

    template <class CT>
    inline auto xscalar<CT>::operator[](size_type) noexcept -> reference
    {
        return m_value;
    }

    template <class CT>
    template <class... Args>
    inline auto xscalar<CT>::operator()(Args...) const noexcept -> const_reference
    {
        XTENSOR_CHECK_DIMENSION((std::array<int, 0>()), Args()...);
        return m_value;
    }

    template <class CT>
    template <class... Args>
    inline auto xscalar<CT>::at(Args... args) const -> const_reference
    {
        check_dimension(shape(), args...);
        return this->operator()(args...);
    }

    template <class CT>
    template <class... Args>
    inline auto xscalar<CT>::unchecked(Args...) const noexcept -> const_reference
    {
        return m_value;
    }

    template <class CT>
    template <class S>
    inline auto xscalar<CT>::operator[](const S&) const noexcept
        -> disable_integral_t<S, const_reference>
    {
        return m_value;
    }

    template <class CT>
    template <class I>
    inline auto xscalar<CT>::operator[](std::initializer_list<I>) const noexcept
        -> const_reference
    {
        return m_value;
    }

    template <class CT>
    inline auto xscalar<CT>::operator[](size_type) const noexcept -> const_reference
    {
        return m_value;
    }

    template <class CT>
    template <class It>
    inline auto xscalar<CT>::element(It, It) noexcept -> reference
    {
        return m_value;
    }

    template <class CT>
    template <class It>
    inline auto xscalar<CT>::element(It, It) const noexcept -> const_reference
    {
        return m_value;
    }

    template <class CT>
    template <class S>
    inline bool xscalar<CT>::broadcast_shape(S&, bool) const noexcept
    {
        return true;
    }

    template <class CT>
    template <class S>
    inline bool xscalar<CT>::is_trivial_broadcast(const S&) const noexcept
    {
        return true;
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::begin() noexcept -> iterator
    {
        return &m_value;
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::end() noexcept -> iterator
    {
        return &m_value + 1;
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::begin() const noexcept -> const_iterator
    {
        return &m_value;
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::end() const noexcept -> const_iterator
    {
        return &m_value + 1;
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::cbegin() const noexcept -> const_iterator
    {
        return &m_value;
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::cend() const noexcept -> const_iterator
    {
        return &m_value + 1;
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::rbegin() noexcept -> reverse_iterator
    {
        return reverse_storage_iterator(end());
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::rend() noexcept -> reverse_iterator
    {
        return reverse_storage_iterator(begin());
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::rbegin() const noexcept -> const_reverse_iterator
    {
        return crbegin();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::rend() const noexcept -> const_reverse_iterator
    {
        return crend();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::crbegin() const noexcept -> const_reverse_iterator
    {
        return const_reverse_storage_iterator(cend());
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::crend() const noexcept -> const_reverse_iterator
    {
        return const_reverse_storage_iterator(cbegin());
    }

    /*****************************
    * Broadcasting iterator api *
    *****************************/

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::begin(const S& shape) noexcept -> broadcast_iterator<S, L>
    {
        return iterable_base::template begin<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::end(const S& shape) noexcept -> broadcast_iterator<S, L>
    {
        return iterable_base::template end<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::begin(const S& shape) const noexcept -> const_broadcast_iterator<S, L>
    {
        return iterable_base::template begin<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::end(const S& shape) const noexcept -> const_broadcast_iterator<S, L>
    {
        return iterable_base::template end<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::cbegin(const S& shape) const noexcept -> const_broadcast_iterator<S, L>
    {
        return iterable_base::template cbegin<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::cend(const S& shape) const noexcept -> const_broadcast_iterator<S, L>
    {
        return iterable_base::template cend<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::rbegin(const S& shape) noexcept -> reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template rbegin<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::rend(const S& shape) noexcept -> reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template rend<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::rbegin(const S& shape) const noexcept -> const_reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template rbegin<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::rend(const S& shape) const noexcept -> const_reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template rend<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::crbegin(const S& shape) const noexcept -> const_reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template crbegin<S, L>(shape);
    }

    template <class CT>
    template <class S, layout_type L>
    inline auto xscalar<CT>::crend(const S& shape) const noexcept -> const_reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template crend<S, L>(shape);
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_begin() noexcept -> iterator
    {
        return this->template begin<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_end() noexcept -> iterator
    {
        return this->template end<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_begin() const noexcept -> const_iterator
    {
        return this->template begin<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_end() const noexcept -> const_iterator
    {
        return this->template end<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_cbegin() const noexcept -> const_iterator
    {
        return this->template cbegin<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_cend() const noexcept -> const_iterator
    {
        return this->template cend<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_rbegin() noexcept -> reverse_iterator
    {
        return this->template rbegin<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_rend() noexcept -> reverse_iterator
    {
        return this->template rend<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_rbegin() const noexcept -> const_reverse_iterator
    {
        return this->template rbegin<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_rend() const noexcept -> const_reverse_iterator
    {
        return this->template rend<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_crbegin() const noexcept -> const_reverse_iterator
    {
        return this->template crbegin<L>();
    }

    template <class CT>
    template <layout_type L>
    inline auto xscalar<CT>::storage_crend() const noexcept -> const_reverse_iterator
    {
        return this->template crend<L>();
    }

    template <class CT>
    template <class S>
    inline auto xscalar<CT>::stepper_begin(const S&) noexcept -> stepper
    {
        return stepper(this, false);
    }

    template <class CT>
    template <class S>
    inline auto xscalar<CT>::stepper_end(const S&, layout_type) noexcept -> stepper
    {
        return stepper(this);
    }

    template <class CT>
    template <class S>
    inline auto xscalar<CT>::stepper_begin(const S&) const noexcept -> const_stepper
    {
        return const_stepper(this);
    }

    template <class CT>
    template <class S>
    inline auto xscalar<CT>::stepper_end(const S&, layout_type) const noexcept -> const_stepper
    {
        return const_stepper(this);
    }

    template <class CT>
    inline auto xscalar<CT>::dummy_begin() noexcept -> dummy_iterator
    {
        return dummy_iterator(this);
    }

    template <class CT>
    inline auto xscalar<CT>::dummy_end() noexcept -> dummy_iterator
    {
        return dummy_iterator(this);
    }

    template <class CT>
    inline auto xscalar<CT>::dummy_begin() const noexcept -> const_dummy_iterator
    {
        return const_dummy_iterator(this);
    }

    template <class CT>
    inline auto xscalar<CT>::dummy_end() const noexcept -> const_dummy_iterator
    {
        return const_dummy_iterator(this);
    }

    template <class CT>
    inline auto xscalar<CT>::data_element(size_type) noexcept -> reference
    {
        return m_value;
    }

    template <class CT>
    inline auto xscalar<CT>::data_element(size_type) const noexcept -> const_reference
    {
        return m_value;
    }

    template <class CT>
    template <class align, class simd>
    inline void xscalar<CT>::store_simd(size_type, const simd& e)
    {
        m_value = static_cast<value_type>(e[0]);
    }

    template <class CT>
    template <class align, class simd>
    inline auto xscalar<CT>::load_simd(size_type) const -> simd
    {
        return xsimd::set_simd<value_type, typename simd::value_type>(m_value);
    }

    template <class T>
    inline xscalar<T&> xref(T& t)
    {
        return xscalar<T&>(t);
    }

    template <class T>
    inline xscalar<const T&> xcref(T& t)
    {
        return xscalar<const T&>(t);
    }

    /**********************************
     * xscalar_stepper implementation *
     **********************************/

    template <bool is_const, class CT>
    inline xscalar_stepper<is_const, CT>::xscalar_stepper(storage_type* c) noexcept
        : p_c(c)
    {
    }

    template <bool is_const, class CT>
    inline auto xscalar_stepper<is_const, CT>::operator*() const noexcept -> reference
    {
        return p_c->operator()();
    }

    template <bool is_const, class CT>
    inline void xscalar_stepper<is_const, CT>::step(size_type /*dim*/, size_type /*n*/) noexcept
    {
    }

    template <bool is_const, class CT>
    inline void xscalar_stepper<is_const, CT>::step_back(size_type /*dim*/, size_type /*n*/) noexcept
    {
    }

    template <bool is_const, class CT>
    inline void xscalar_stepper<is_const, CT>::reset(size_type /*dim*/) noexcept
    {
    }

    template <bool is_const, class CT>
    inline void xscalar_stepper<is_const, CT>::reset_back(size_type /*dim*/) noexcept
    {
    }

    template <bool is_const, class CT>
    inline void xscalar_stepper<is_const, CT>::to_begin() noexcept
    {
    }

    template <bool is_const, class CT>
    inline void xscalar_stepper<is_const, CT>::to_end(layout_type /*l*/) noexcept
    {
    }

    template <bool is_const, class CT>
    template <class R>
    inline R xscalar_stepper<is_const, CT>::step_simd()
    {
        return R(p_c->operator()());
    }

    template <bool is_const, class CT>
    inline auto xscalar_stepper<is_const, CT>::step_leading()
        -> value_type
    {
        return p_c->operator()();
    }

    /**********************************
     * xdummy_iterator implementation *
     **********************************/

    template <bool is_const, class CT>
    inline xdummy_iterator<is_const, CT>::xdummy_iterator(storage_type* c) noexcept
        : p_c(c)
    {
    }

    template <bool is_const, class CT>
    inline auto xdummy_iterator<is_const, CT>::operator++() noexcept -> self_type&
    {
        return *this;
    }

    template <bool is_const, class CT>
    inline auto xdummy_iterator<is_const, CT>::operator--() noexcept -> self_type&
    {
        return *this;
    }

    template <bool is_const, class CT>
    inline auto xdummy_iterator<is_const, CT>::operator+=(difference_type) noexcept -> self_type&
    {
        return *this;
    }

    template <bool is_const, class CT>
    inline auto xdummy_iterator<is_const, CT>::operator-=(difference_type) noexcept -> self_type&
    {
        return *this;
    }

    template <bool is_const, class CT>
    inline auto xdummy_iterator<is_const, CT>::operator-(const self_type&) const noexcept -> difference_type
    {
        return 0;
    }

    template <bool is_const, class CT>
    inline auto xdummy_iterator<is_const, CT>::operator*() const noexcept -> reference
    {
        return p_c->operator()();
    }

    template <bool is_const, class CT>
    inline bool xdummy_iterator<is_const, CT>::equal(const self_type& rhs) const noexcept
    {
        return p_c == rhs.p_c;
    }

    template <bool is_const, class CT>
    inline bool xdummy_iterator<is_const, CT>::less_than(const self_type& rhs) const noexcept
    {
        return p_c < rhs.p_c;
    }

    template <bool is_const, class CT>
    inline bool operator==(const xdummy_iterator<is_const, CT>& lhs,
                           const xdummy_iterator<is_const, CT>& rhs) noexcept
    {
        return lhs.equal(rhs);
    }

    template <bool is_const, class CT>
    inline bool operator<(const xdummy_iterator<is_const, CT>& lhs,
                          const xdummy_iterator<is_const, CT>& rhs) noexcept
    {
        return lhs.less_than(rhs);
    }
}

#endif
