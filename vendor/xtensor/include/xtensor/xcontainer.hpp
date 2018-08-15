/***************************************************************************
* Copyright (c) 2016, Johan Mabille, Sylvain Corlay and Wolf Vollprecht    *
*                                                                          *
* Distributed under the terms of the BSD 3-Clause License.                 *
*                                                                          *
* The full license is in the file LICENSE, distributed with this software. *
****************************************************************************/

#ifndef XTENSOR_CONTAINER_HPP
#define XTENSOR_CONTAINER_HPP

#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include <xtl/xmeta_utils.hpp>
#include <xtl/xsequence.hpp>

#include "xiterable.hpp"
#include "xiterator.hpp"
#include "xmath.hpp"
#include "xoperation.hpp"
#include "xstrides.hpp"
#include "xtensor_forward.hpp"

namespace xt
{
    template <class D>
    struct xcontainer_iterable_types
    {
        using inner_shape_type = typename xcontainer_inner_types<D>::inner_shape_type;
        using storage_type = typename xcontainer_inner_types<D>::storage_type;
        using stepper = xstepper<D>;
        using const_stepper = xstepper<const D>;
    };

#define DL XTENSOR_DEFAULT_LAYOUT

    namespace detail
    {
        template <class T>
        struct allocator_type_impl
        {
            using type = typename T::allocator_type;
        };

        template <class T, std::size_t N>
        struct allocator_type_impl<std::array<T, N>>
        {
            using type = std::allocator<T>;  // fake allocator for testing
        };
    }

    template <class T>
    using allocator_type_t = typename detail::allocator_type_impl<T>::type;


    /**
     * @class xcontainer
     * @brief Base class for dense multidimensional containers.
     *
     * The xcontainer class defines the interface for dense multidimensional
     * container classes. It does not embed any data container, this responsibility
     * is delegated to the inheriting classes.
     *
     * @tparam D The derived type, i.e. the inheriting class for which xcontainer
     *           provides the interface.
     */
    template <class D>
    class xcontainer : private xiterable<D>
    {
    public:

        using derived_type = D;

        using inner_types = xcontainer_inner_types<D>;
        using storage_type = typename inner_types::storage_type;
        using allocator_type = allocator_type_t<std::decay_t<storage_type>>;
        using value_type = typename storage_type::value_type;
        using reference = std::conditional_t<std::is_const<storage_type>::value,
                                             typename storage_type::const_reference,
                                             typename storage_type::reference>;
        using const_reference = typename storage_type::const_reference;
        using pointer = typename storage_type::pointer;
        using const_pointer = typename storage_type::const_pointer;
        using size_type = typename storage_type::size_type;
        using difference_type = typename storage_type::difference_type;
        using simd_value_type = xsimd::simd_type<value_type>;

        using shape_type = typename inner_types::shape_type;
        using strides_type = typename inner_types::strides_type;
        using backstrides_type = typename inner_types::backstrides_type;

        using inner_shape_type = typename inner_types::inner_shape_type;
        using inner_strides_type = typename inner_types::inner_strides_type;
        using inner_backstrides_type = typename inner_types::inner_backstrides_type;

        using iterable_base = xiterable<D>;
        using stepper = typename iterable_base::stepper;
        using const_stepper = typename iterable_base::const_stepper;

        static constexpr layout_type static_layout = inner_types::layout;
        static constexpr bool contiguous_layout = static_layout != layout_type::dynamic;
        using data_alignment = xsimd::container_alignment_t<storage_type>;
        using simd_type = xsimd::simd_type<value_type>;

        static_assert(static_layout != layout_type::any, "Container layout can never be layout_type::any!");

        size_type size() const noexcept;

        constexpr size_type dimension() const noexcept;

        constexpr const inner_shape_type& shape() const noexcept;
        constexpr const inner_strides_type& strides() const noexcept;
        constexpr const inner_backstrides_type& backstrides() const noexcept;

        template <class T>
        void fill(const T& value);

        template <class... Args>
        reference operator()(Args... args);

        template <class... Args>
        const_reference operator()(Args... args) const;

        template <class... Args>
        reference at(Args... args);

        template <class... Args>
        const_reference at(Args... args) const;

        template <class... Args>
        reference unchecked(Args... args);

        template <class... Args>
        const_reference unchecked(Args... args) const;

        template <class S>
        disable_integral_t<S, reference> operator[](const S& index);
        template <class I>
        reference operator[](std::initializer_list<I> index);
        reference operator[](size_type i);

        template <class S>
        disable_integral_t<S, const_reference> operator[](const S& index) const;
        template <class I>
        const_reference operator[](std::initializer_list<I> index) const;
        const_reference operator[](size_type i) const;

        template <class It>
        reference element(It first, It last);
        template <class It>
        const_reference element(It first, It last) const;

        storage_type& storage() noexcept;
        const storage_type& storage() const noexcept;

        value_type* data() noexcept;
        const value_type* data() const noexcept;
        const size_type data_offset() const noexcept;

        template <class S>
        bool broadcast_shape(S& shape, bool reuse_cache = false) const;

        template <class S>
        bool is_trivial_broadcast(const S& strides) const noexcept;

        template <class S>
        stepper stepper_begin(const S& shape) noexcept;
        template <class S>
        stepper stepper_end(const S& shape, layout_type l) noexcept;

        template <class S>
        const_stepper stepper_begin(const S& shape) const noexcept;
        template <class S>
        const_stepper stepper_end(const S& shape, layout_type l) const noexcept;

        reference data_element(size_type i);
        const_reference data_element(size_type i) const;

        template <class simd>
        using simd_return_type = xsimd::simd_return_type<value_type, typename simd::value_type>;

        template <class align, class simd = simd_value_type>
        void store_simd(size_type i, const simd& e);
        template <class align, class simd = simd_value_type>
        simd_return_type<simd> load_simd(size_type i) const;

#if defined(_MSC_VER) && _MSC_VER >= 1910
        // Workaround for compiler bug in Visual Studio 2017 with respect to alias templates with non-type parameters.
        template <layout_type L>
        using layout_iterator = xiterator<typename iterable_base::stepper, typename iterable_base::inner_shape_type*, L>;
        template <layout_type L>
        using const_layout_iterator = xiterator<typename iterable_base::const_stepper, typename iterable_base::inner_shape_type*, L>;
        template <layout_type L>
        using reverse_layout_iterator = std::reverse_iterator<layout_iterator<L>>;
        template <layout_type L>
        using const_reverse_layout_iterator = std::reverse_iterator<const_layout_iterator<L>>;
#else
        template <layout_type L>
        using layout_iterator = typename iterable_base::template layout_iterator<L>;
        template <layout_type L>
        using const_layout_iterator = typename iterable_base::template const_layout_iterator<L>;
        template <layout_type L>
        using reverse_layout_iterator = typename iterable_base::template reverse_layout_iterator<L>;
        template <layout_type L>
        using const_reverse_layout_iterator = typename iterable_base::template const_reverse_layout_iterator<L>;
#endif

        template <class S, layout_type L>
        using broadcast_iterator = typename iterable_base::template broadcast_iterator<S, L>;
        template <class S, layout_type L>
        using const_broadcast_iterator = typename iterable_base::template const_broadcast_iterator<S, L>;
        template <class S, layout_type L>
        using reverse_broadcast_iterator = typename iterable_base::template reverse_broadcast_iterator<S, L>;
        template <class S, layout_type L>
        using const_reverse_broadcast_iterator = typename iterable_base::template const_reverse_broadcast_iterator<S, L>;

        using storage_iterator = typename storage_type::iterator;
        using const_storage_iterator = typename storage_type::const_iterator;
        using reverse_storage_iterator = typename storage_type::reverse_iterator;
        using const_reverse_storage_iterator = typename storage_type::const_reverse_iterator;

        template <layout_type L, class It1, class It2>
        using select_iterator_impl = std::conditional_t<L == static_layout, It1, It2>;

        template <layout_type L>
        using select_iterator = select_iterator_impl<L, storage_iterator, layout_iterator<L>>;
        template <layout_type L>
        using select_const_iterator = select_iterator_impl<L, const_storage_iterator, const_layout_iterator<L>>;
        template <layout_type L>
        using select_reverse_iterator = select_iterator_impl<L, reverse_storage_iterator, reverse_layout_iterator<L>>;
        template <layout_type L>
        using select_const_reverse_iterator = select_iterator_impl<L, const_reverse_storage_iterator, const_reverse_layout_iterator<L>>;

        using iterator = select_iterator<DL>;
        using const_iterator = select_const_iterator<DL>;
        using reverse_iterator = select_reverse_iterator<DL>;
        using const_reverse_iterator = select_const_reverse_iterator<DL>;

        template <layout_type L = DL>
        select_iterator<L> begin() noexcept;
        template <layout_type L = DL>
        select_iterator<L> end() noexcept;

        template <layout_type L = DL>
        select_const_iterator<L> begin() const noexcept;
        template <layout_type L = DL>
        select_const_iterator<L> end() const noexcept;
        template <layout_type L = DL>
        select_const_iterator<L> cbegin() const noexcept;
        template <layout_type L = DL>
        select_const_iterator<L> cend() const noexcept;

        template <layout_type L = DL>
        select_reverse_iterator<L> rbegin() noexcept;
        template <layout_type L = DL>
        select_reverse_iterator<L> rend() noexcept;

        template <layout_type L = DL>
        select_const_reverse_iterator<L> rbegin() const noexcept;
        template <layout_type L = DL>
        select_const_reverse_iterator<L> rend() const noexcept;
        template <layout_type L = DL>
        select_const_reverse_iterator<L> crbegin() const noexcept;
        template <layout_type L = DL>
        select_const_reverse_iterator<L> crend() const noexcept;

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
        storage_iterator storage_begin() noexcept;
        template <layout_type L = DL>
        storage_iterator storage_end() noexcept;

        template <layout_type L = DL>
        const_storage_iterator storage_begin() const noexcept;
        template <layout_type L = DL>
        const_storage_iterator storage_end() const noexcept;
        template <layout_type L = DL>
        const_storage_iterator storage_cbegin() const noexcept;
        template <layout_type L = DL>
        const_storage_iterator storage_cend() const noexcept;

        template <layout_type L = DL>
        reverse_storage_iterator storage_rbegin() noexcept;
        template <layout_type L = DL>
        reverse_storage_iterator storage_rend() noexcept;

        template <layout_type L = DL>
        const_reverse_storage_iterator storage_rbegin() const noexcept;
        template <layout_type L = DL>
        const_reverse_storage_iterator storage_rend() const noexcept;
        template <layout_type L = DL>
        const_reverse_storage_iterator storage_crbegin() const noexcept;
        template <layout_type L = DL>
        const_reverse_storage_iterator storage_crend() const noexcept;

        using container_iterator = storage_iterator;
        using const_container_iterator = const_storage_iterator;

    protected:

        xcontainer() = default;
        ~xcontainer() = default;

        xcontainer(const xcontainer&) = default;
        xcontainer& operator=(const xcontainer&) = default;

        xcontainer(xcontainer&&) = default;
        xcontainer& operator=(xcontainer&&) = default;

        container_iterator data_xbegin() noexcept;
        const_container_iterator data_xbegin() const noexcept;
        container_iterator data_xend(layout_type l) noexcept;
        const_container_iterator data_xend(layout_type l) const noexcept;

    protected:

        derived_type& derived_cast() & noexcept;
        const derived_type& derived_cast() const & noexcept;
        derived_type derived_cast() && noexcept;

    private:

        friend class xiterable<D>;
        friend class xconst_iterable<D>;

        template <class C>
        friend class xstepper;

        template <class It>
        It data_xend_impl(It end, layout_type l) const noexcept;

        inner_shape_type& mutable_shape();
        inner_strides_type& mutable_strides();
        inner_backstrides_type& mutable_backstrides();
    };

#undef DL

    /**
     * @class xstrided_container
     * @brief Partial implementation of xcontainer that embeds the strides and the shape
     *
     * The xstrided_container class is a partial implementation of the xcontainer interface
     * that embed the strides and the shape of the multidimensional container. It does
     * not embed the data container, this responsibility is delegated to the inheriting
     * classes.
     *
     * @tparam D The derived type, i.e. the inheriting class for which xstrided_container
     *           provides the partial imlpementation of xcontainer.
     */
    template <class D>
    class xstrided_container : public xcontainer<D>
    {
    public:

        using base_type = xcontainer<D>;
        using storage_type = typename base_type::storage_type;
        using value_type = typename base_type::value_type;
        using reference = typename base_type::reference;
        using const_reference = typename base_type::const_reference;
        using pointer = typename base_type::pointer;
        using const_pointer = typename base_type::const_pointer;
        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;
        using strides_type = typename base_type::strides_type;
        using inner_shape_type = typename base_type::inner_shape_type;
        using inner_strides_type = typename base_type::inner_strides_type;
        using inner_backstrides_type = typename base_type::inner_backstrides_type;

        template <class S = shape_type>
        void resize(S&& shape, bool force = false);
        template <class S = shape_type>
        void resize(S&& shape, layout_type l);
        template <class S = shape_type>
        void resize(S&& shape, const strides_type& strides);

        template <class S = shape_type>
        void reshape(S&& shape, layout_type layout = base_type::static_layout);

        layout_type layout() const noexcept;

    protected:

        xstrided_container() noexcept;
        ~xstrided_container() = default;

        xstrided_container(const xstrided_container&) = default;
        xstrided_container& operator=(const xstrided_container&) = default;

        xstrided_container(xstrided_container&&) = default;
        xstrided_container& operator=(xstrided_container&&) = default;

        explicit xstrided_container(inner_shape_type&&, inner_strides_type&&) noexcept;

        inner_shape_type& shape_impl() noexcept;
        const inner_shape_type& shape_impl() const noexcept;

        inner_strides_type& strides_impl() noexcept;
        const inner_strides_type& strides_impl() const noexcept;

        inner_backstrides_type& backstrides_impl() noexcept;
        const inner_backstrides_type& backstrides_impl() const noexcept;

    private:

        inner_shape_type m_shape;
        inner_strides_type m_strides;
        inner_backstrides_type m_backstrides;
        layout_type m_layout = base_type::static_layout;
    };

    /******************************
     * xcontainer implementation *
     ******************************/

    template <class D>
    template <class It>
    inline It xcontainer<D>::data_xend_impl(It end, layout_type l) const noexcept
    {
        return strided_data_end(*this, end, l);
    }

    template <class D>
    inline auto xcontainer<D>::mutable_shape() -> inner_shape_type&
    {
        return derived_cast().shape_impl();
    }

    template <class D>
    inline auto xcontainer<D>::mutable_strides() -> inner_strides_type&
    {
        return derived_cast().strides_impl();
    }

    template <class D>
    inline auto xcontainer<D>::mutable_backstrides() -> inner_backstrides_type&
    {
        return derived_cast().backstrides_impl();
    }

    /**
     * @name Size and shape
     */
    //@{
    /**
     * Returns the number of element in the container.
     */
    template <class D>
    inline auto xcontainer<D>::size() const noexcept -> size_type
    {
        return contiguous_layout ? storage().size() : compute_size(shape());
    }

    /**
     * Returns the number of dimensions of the container.
     */
    template <class D>
    inline constexpr auto xcontainer<D>::dimension() const noexcept -> size_type
    {
        return shape().size();
    }

    /**
     * Returns the shape of the container.
     */
    template <class D>
    constexpr inline auto xcontainer<D>::shape() const noexcept -> const inner_shape_type&
    {
        return derived_cast().shape_impl();
    }

    /**
     * Returns the strides of the container.
     */
    template <class D>
    constexpr inline auto xcontainer<D>::strides() const noexcept -> const inner_strides_type&
    {
        return derived_cast().strides_impl();
    }

    /**
     * Returns the backstrides of the container.
     */
    template <class D>
    constexpr inline auto xcontainer<D>::backstrides() const noexcept -> const inner_backstrides_type&
    {
        return derived_cast().backstrides_impl();
    }
    //@}

    /**
     * @name Data
     */
    //@{

    /**
     * Fills the container with the given value.
     * @param value the value to fill the container with.
     */
    template <class D>
    template <class T>
    inline void xcontainer<D>::fill(const T& value)
    {
        std::fill(storage_begin(), storage_end(), value);
    }

    /**
     * Returns a reference to the element at the specified position in the container.
     * @param args a list of indices specifying the position in the container. Indices
     * must be unsigned integers, the number of indices should be equal or greater than
     * the number of dimensions of the container.
     */
    template <class D>
    template <class... Args>
    inline auto xcontainer<D>::operator()(Args... args) -> reference
    {
        XTENSOR_TRY(check_index(shape(), args...));
        XTENSOR_CHECK_DIMENSION(shape(), args...);
        size_type index = xt::data_offset<size_type>(strides(), static_cast<size_type>(args)...);
        return storage()[index];
    }

    /**
     * Returns a constant reference to the element at the specified position in the container.
     * @param args a list of indices specifying the position in the container. Indices
     * must be unsigned integers, the number of indices should be equal or greater than
     * the number of dimensions of the container.
     */
    template <class D>
    template <class... Args>
    inline auto xcontainer<D>::operator()(Args... args) const -> const_reference
    {
        XTENSOR_TRY(check_index(shape(), args...));
        XTENSOR_CHECK_DIMENSION(shape(), args...);
        size_type index = xt::data_offset<size_type>(strides(), static_cast<size_type>(args)...);
        return storage()[index];
    }

    /**
     * Returns a reference to the element at the specified position in the container,
     * after dimension and bounds checking.
     * @param args a list of indices specifying the position in the container. Indices
     * must be unsigned integers, the number of indices should be equal to the number of dimensions
     * of the container.
     * @exception std::out_of_range if the number of argument is greater than the number of dimensions
     * or if indices are out of bounds.
     */
    template <class D>
    template <class... Args>
    inline auto xcontainer<D>::at(Args... args) -> reference
    {
        check_access(shape(), static_cast<size_type>(args)...);
        return this->operator()(args...);
    }

    /**
     * Returns a constant reference to the element at the specified position in the container,
     * after dimension and bounds checking.
     * @param args a list of indices specifying the position in the container. Indices
     * must be unsigned integers, the number of indices should be equal to the number of dimensions
     * of the container.
     * @exception std::out_of_range if the number of argument is greater than the number of dimensions
     * or if indices are out of bounds.
     */
    template <class D>
    template <class... Args>
    inline auto xcontainer<D>::at(Args... args) const -> const_reference
    {
        check_access(shape(), static_cast<size_type>(args)...);
        return this->operator()(args...);
    }

    /**
     * Returns a reference to the element at the specified position in the container.
     * @param args a list of indices specifying the position in the container. Indices
     * must be unsigned integers, the number of indices must be equal to the number of
     * dimensions of the container, else the behavior is undefined.
     *
     * @warning This method is meant for performance, for expressions with a dynamic
     * number of dimensions (i.e. not known at compile time). Since it may have
     * undefined behavior (see parameters), operator() should be prefered whenever
     * it is possible.
     * @warning This method is NOT compatible with broadcasting, meaning the following
     * code has undefined behavior:
     * \code{.cpp}
     * xt::xarray<double> a = {{0, 1}, {2, 3}};
     * xt::xarray<double> b = {0, 1};
     * auto fd = a + b;
     * double res = fd.uncheked(0, 1);
     * \endcode
     */
    template <class D>
    template <class... Args>
    inline auto xcontainer<D>::unchecked(Args... args) -> reference
    {
        size_type index = xt::unchecked_data_offset<size_type>(strides(), static_cast<size_type>(args)...);
        return storage()[index];
    }

    /**
     * Returns a constant reference to the element at the specified position in the container.
     * @param args a list of indices specifying the position in the container. Indices
     * must be unsigned integers, the number of indices must be equal to the number of
     * dimensions of the container, else the behavior is undefined.
     *
     * @warning This method is meant for performance, for expressions with a dynamic
     * number of dimensions (i.e. not known at compile time). Since it may have
     * undefined behavior (see parameters), operator() should be prefered whenever
     * it is possible.
     * @warning This method is NOT compatible with broadcasting, meaning the following
     * code has undefined behavior:
     * \code{.cpp}
     * xt::xarray<double> a = {{0, 1}, {2, 3}};
     * xt::xarray<double> b = {0, 1};
     * auto fd = a + b;
     * double res = fd.uncheked(0, 1);
     * \endcode
     */
    template <class D>
    template <class... Args>
    inline auto xcontainer<D>::unchecked(Args... args) const -> const_reference
    {
        size_type index = xt::unchecked_data_offset<size_type>(strides(), static_cast<size_type>(args)...);
        return storage()[index];
    }

    /**
     * Returns a reference to the element at the specified position in the container.
     * @param index a sequence of indices specifying the position in the container. Indices
     * must be unsigned integers, the number of indices in the list should be equal or greater
     * than the number of dimensions of the container.
     */
    template <class D>
    template <class S>
    inline auto xcontainer<D>::operator[](const S& index)
        -> disable_integral_t<S, reference>
    {
        return element(index.cbegin(), index.cend());
    }

    template <class D>
    template <class I>
    inline auto xcontainer<D>::operator[](std::initializer_list<I> index) -> reference
    {
        return element(index.begin(), index.end());
    }

    template <class D>
    inline auto xcontainer<D>::operator[](size_type i) -> reference
    {
        return operator()(i);
    }

    /**
     * Returns a constant reference to the element at the specified position in the container.
     * @param index a sequence of indices specifying the position in the container. Indices
     * must be unsigned integers, the number of indices in the list should be equal or greater
     * than the number of dimensions of the container.
     */
    template <class D>
    template <class S>
    inline auto xcontainer<D>::operator[](const S& index) const
        -> disable_integral_t<S, const_reference>
    {
        return element(index.cbegin(), index.cend());
    }

    template <class D>
    template <class I>
    inline auto xcontainer<D>::operator[](std::initializer_list<I> index) const -> const_reference
    {
        return element(index.begin(), index.end());
    }

    template <class D>
    inline auto xcontainer<D>::operator[](size_type i) const -> const_reference
    {
        return operator()(i);
    }

    /**
     * Returns a reference to the element at the specified position in the container.
     * @param first iterator starting the sequence of indices
     * @param last iterator ending the sequence of indices
     * The number of indices in the sequence should be equal to or greater
     * than the number of dimensions of the container.
     */
    template <class D>
    template <class It>
    inline auto xcontainer<D>::element(It first, It last) -> reference
    {
        XTENSOR_TRY(check_element_index(shape(), first, last));
        return storage()[element_offset<size_type>(strides(), first, last)];
    }

    /**
     * Returns a reference to the element at the specified position in the container.
     * @param first iterator starting the sequence of indices
     * @param last iterator ending the sequence of indices
     * The number of indices in the sequence should be equal to or greater
     * than the number of dimensions of the container.
     */
    template <class D>
    template <class It>
    inline auto xcontainer<D>::element(It first, It last) const -> const_reference
    {
        XTENSOR_TRY(check_element_index(shape(), first, last));
        return storage()[element_offset<size_type>(strides(), first, last)];
    }

    /**
     * Returns a reference to the buffer containing the elements of the container.
     */
    template <class D>
    inline auto xcontainer<D>::storage() noexcept -> storage_type&
    {
        return derived_cast().storage_impl();
    }

    /**
     * Returns a constant reference to the buffer containing the elements of the
     * container.
     */
    template <class D>
    inline auto xcontainer<D>::storage() const noexcept -> const storage_type&
    {
        return derived_cast().storage_impl();
    }

    /**
     * Returns a pointer to the underlying array serving as element storage. The pointer
     * is such that range [data(); data() + size()] is always a valid range, even if the
     * container is empty (data() is not is not dereferenceable in that case)
     */
    template <class D>
    inline auto xcontainer<D>::data() noexcept -> value_type*
    {
        return storage().data();
    }

    /**
    * Returns a constant pointer to the underlying array serving as element storage. The pointer
    * is such that range [data(); data() + size()] is always a valid range, even if the
    * container is empty (data() is not is not dereferenceable in that case)
    */
    template <class D>
    inline auto xcontainer<D>::data() const noexcept -> const value_type*
    {
        return storage().data();
    }

    /**
     * Returns the offset to the first element in the container.
     */
    template <class D>
    inline auto xcontainer<D>::data_offset() const noexcept -> const size_type
    {
        return size_type(0);
    }
    //@}

    /**
     * @name Broadcasting
     */
    //@{
    /**
     * Broadcast the shape of the container to the specified parameter.
     * @param shape the result shape
     * @param reuse_cache parameter for internal optimization
     * @return a boolean indicating whether the broadcasting is trivial
     */
    template <class D>
    template <class S>
    inline bool xcontainer<D>::broadcast_shape(S& shape, bool) const
    {
        return xt::broadcast_shape(this->shape(), shape);
    }

    /**
     * Compares the specified strides with those of the container to see whether
     * the broadcasting is trivial.
     * @return a boolean indicating whether the broadcasting is trivial
     */
    template <class D>
    template <class S>
    inline bool xcontainer<D>::is_trivial_broadcast(const S& str) const noexcept
    {
        return str.size() == strides().size() &&
            std::equal(str.cbegin(), str.cend(), strides().begin());
    }
    //@}

    /****************
     * Iterator api *
     ****************/

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::begin() noexcept -> select_iterator<L>
    {
        return xtl::mpl::static_if<L == static_layout>([&](auto self)
        {
            return self(*this).template storage_begin<L>();
        }, /*else*/ [&](auto self)
        {
            return self(*this).iterable_base::template begin<L>();
        });
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::end() noexcept -> select_iterator<L>
    {
        return xtl::mpl::static_if<L == static_layout>([&](auto self)
        {
            return self(*this).template storage_end<L>();
        }, /*else*/ [&](auto self)
        {
            return self(*this).iterable_base::template end<L>();
        });
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::begin() const noexcept -> select_const_iterator<L>
    {
        return this->template cbegin<L>();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::end() const noexcept -> select_const_iterator<L>
    {
        return this->template cend<L>();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::cbegin() const noexcept -> select_const_iterator<L>
    {
        return xtl::mpl::static_if<L == static_layout>([&](auto self)
        {
            return self(*this).template storage_cbegin<L>();
        }, /*else*/ [&](auto self)
        {
            return self(*this).iterable_base::template cbegin<L>();
        });
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::cend() const noexcept -> select_const_iterator<L>
    {
        return xtl::mpl::static_if<L == static_layout>([&](auto self)
        {
            return self(*this).template storage_cend<L>();
        }, /*else*/ [&](auto self)
        {
            return self(*this).iterable_base::template cend<L>();
        });
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::rbegin() noexcept -> select_reverse_iterator<L>
    {
        return xtl::mpl::static_if<L == static_layout>([&](auto self)
        {
            return self(*this).template storage_rbegin<L>();
        }, /*else*/ [&](auto self)
        {
            return self(*this).iterable_base::template rbegin<L>();
        });
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::rend() noexcept -> select_reverse_iterator<L>
    {
        return xtl::mpl::static_if<L == static_layout>([&](auto self)
        {
            return self(*this).template storage_rend<L>();
        }, /*else*/ [&](auto self)
        {
            return self(*this).iterable_base::template rend<L>();
        });
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::rbegin() const noexcept -> select_const_reverse_iterator<L>
    {
        return this->template crbegin<L>();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::rend() const noexcept -> select_const_reverse_iterator<L>
    {
        return this->template crend<L>();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::crbegin() const noexcept -> select_const_reverse_iterator<L>
    {
        return xtl::mpl::static_if<L == static_layout>([&](auto self)
        {
            return self(*this).template storage_crbegin<L>();
        }, /*else*/ [&](auto self)
        {
            return self(*this).iterable_base::template crbegin<L>();
        });
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::crend() const noexcept -> select_const_reverse_iterator<L>
    {
        return xtl::mpl::static_if<L == static_layout>([&](auto self)
        {
            return self(*this).template storage_crend<L>();
        }, /*else*/ [&](auto self)
        {
            return self(*this).iterable_base::template crend<L>();
        });
    }

    /*****************************
     * Broadcasting iterator api *
     *****************************/

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::begin(const S& shape) noexcept -> broadcast_iterator<S, L>
    {
        return iterable_base::template begin<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::end(const S& shape) noexcept -> broadcast_iterator<S, L>
    {
        return iterable_base::template end<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::begin(const S& shape) const noexcept -> const_broadcast_iterator<S, L>
    {
        return iterable_base::template begin<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::end(const S& shape) const noexcept -> const_broadcast_iterator<S, L>
    {
        return iterable_base::template end<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::cbegin(const S& shape) const noexcept -> const_broadcast_iterator<S, L>
    {
        return iterable_base::template cbegin<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::cend(const S& shape) const noexcept -> const_broadcast_iterator<S, L>
    {
        return iterable_base::template cend<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::rbegin(const S& shape) noexcept -> reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template rbegin<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::rend(const S& shape) noexcept -> reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template rend<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::rbegin(const S& shape) const noexcept -> const_reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template rbegin<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::rend(const S& shape) const noexcept -> const_reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template rend<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::crbegin(const S& shape) const noexcept -> const_reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template crbegin<S, L>(shape);
    }

    template <class D>
    template <class S, layout_type L>
    inline auto xcontainer<D>::crend(const S& shape) const noexcept -> const_reverse_broadcast_iterator<S, L>
    {
        return iterable_base::template crend<S, L>(shape);
    }

    /***********************
     * Linear iterator api *
     ***********************/

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_begin() noexcept -> storage_iterator
    {
        return storage().begin();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_end() noexcept -> storage_iterator
    {
        return storage().end();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_begin() const noexcept -> const_storage_iterator
    {
        return storage().begin();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_end() const noexcept -> const_storage_iterator
    {
        return storage().end();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_cbegin() const noexcept -> const_storage_iterator
    {
        return storage().cbegin();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_cend() const noexcept -> const_storage_iterator
    {
        return storage().cend();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_rbegin() noexcept -> reverse_storage_iterator
    {
        return storage().rbegin();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_rend() noexcept -> reverse_storage_iterator
    {
        return storage().rend();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_rbegin() const noexcept -> const_reverse_storage_iterator
    {
        return storage().rbegin();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_rend() const noexcept -> const_reverse_storage_iterator
    {
        return storage().rend();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_crbegin() const noexcept -> const_reverse_storage_iterator
    {
        return storage().crbegin();
    }

    template <class D>
    template <layout_type L>
    inline auto xcontainer<D>::storage_crend() const noexcept -> const_reverse_storage_iterator
    {
        return storage().crend();
    }

    /***************
     * stepper api *
     ***************/

    template <class D>
    template <class S>
    inline auto xcontainer<D>::stepper_begin(const S& shape) noexcept -> stepper
    {
        size_type offset = shape.size() - dimension();
        return stepper(static_cast<derived_type*>(this), data_xbegin(), offset);
    }

    template <class D>
    template <class S>
    inline auto xcontainer<D>::stepper_end(const S& shape, layout_type l) noexcept -> stepper
    {
        size_type offset = shape.size() - dimension();
        return stepper(static_cast<derived_type*>(this), data_xend(l), offset);
    }

    template <class D>
    template <class S>
    inline auto xcontainer<D>::stepper_begin(const S& shape) const noexcept -> const_stepper
    {
        size_type offset = shape.size() - dimension();
        return const_stepper(static_cast<const derived_type*>(this), data_xbegin(), offset);
    }

    template <class D>
    template <class S>
    inline auto xcontainer<D>::stepper_end(const S& shape, layout_type l) const noexcept -> const_stepper
    {
        size_type offset = shape.size() - dimension();
        return const_stepper(static_cast<const derived_type*>(this), data_xend(l), offset);
    }

    template <class D>
    inline auto xcontainer<D>::data_xbegin() noexcept -> container_iterator
    {
        return storage().begin();
    }

    template <class D>
    inline auto xcontainer<D>::data_xbegin() const noexcept -> const_container_iterator
    {
        return storage().begin();
    }

    template <class D>
    inline auto xcontainer<D>::data_xend(layout_type l) noexcept -> container_iterator
    {
        return data_xend_impl(storage().end(), l);
    }

    template <class D>
    inline auto xcontainer<D>::data_xend(layout_type l) const noexcept -> const_container_iterator
    {
        return data_xend_impl(storage().end(), l);
    }

    template <class D>
    inline auto xcontainer<D>::derived_cast() & noexcept -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    template <class D>
    inline auto xcontainer<D>::derived_cast() const & noexcept -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

    template <class D>
    inline auto xcontainer<D>::derived_cast() && noexcept -> derived_type
    {
        return *static_cast<derived_type*>(this);
    }

    template <class D>
    inline auto xcontainer<D>::data_element(size_type i) -> reference
    {
        return storage()[i];
    }

    template <class D>
    inline auto xcontainer<D>::data_element(size_type i) const -> const_reference
    {
        return storage()[i];
    }

    template <class D>
    template <class alignment, class simd>
    inline void xcontainer<D>::store_simd(size_type i, const simd& e)
    {
        using align_mode = driven_align_mode_t<alignment, data_alignment>;
        xsimd::store_simd<value_type, typename simd::value_type>(&(storage()[i]), e, align_mode());
    }

    template <class D>
    template <class alignment, class simd>
    inline auto xcontainer<D>::load_simd(size_type i) const -> simd_return_type<simd>
    {
        using align_mode = driven_align_mode_t<alignment, data_alignment>;
        return xsimd::load_simd<value_type, typename simd::value_type>(&(storage()[i]), align_mode());
    }

    /*************************************
     * xstrided_container implementation *
     *************************************/

    template <class D>
    inline xstrided_container<D>::xstrided_container() noexcept
        : base_type()
    {
        m_shape = xtl::make_sequence<inner_shape_type>(base_type::dimension(), 1);
        m_strides = xtl::make_sequence<inner_strides_type>(base_type::dimension(), 0);
        m_backstrides = xtl::make_sequence<inner_backstrides_type>(base_type::dimension(), 0);
    }

    template <class D>
    inline xstrided_container<D>::xstrided_container(inner_shape_type&& shape, inner_strides_type&& strides) noexcept
        : base_type(), m_shape(std::move(shape)), m_strides(std::move(strides))
    {
        m_backstrides = xtl::make_sequence<inner_backstrides_type>(m_shape.size(), 0);
        adapt_strides(m_shape, m_strides, m_backstrides);
    }

    template <class D>
    inline auto xstrided_container<D>::shape_impl() noexcept -> inner_shape_type&
    {
        return m_shape;
    }

    template <class D>
    inline auto xstrided_container<D>::shape_impl() const noexcept -> const inner_shape_type&
    {
        return m_shape;
    }

    template <class D>
    inline auto xstrided_container<D>::strides_impl() noexcept -> inner_strides_type&
    {
        return m_strides;
    }

    template <class D>
    inline auto xstrided_container<D>::strides_impl() const noexcept -> const inner_strides_type&
    {
        return m_strides;
    }

    template <class D>
    inline auto xstrided_container<D>::backstrides_impl() noexcept -> inner_backstrides_type&
    {
        return m_backstrides;
    }

    template <class D>
    inline auto xstrided_container<D>::backstrides_impl() const noexcept -> const inner_backstrides_type&
    {
        return m_backstrides;
    }

    /**
     * Return the layout_type of the container
     * @return layout_type of the container
     */
    template <class D>
    layout_type xstrided_container<D>::layout() const noexcept
    {
        return m_layout;
    }

    namespace detail
    {
        template <class C, class S>
        inline void resize_data_container(C& c, S size)
        {
            xt::resize_container(c, size);
        }

        template <class C, class S>
        inline void resize_data_container(const C& c, S size)
        {
            (void) c;  // remove unused parameter warning
            (void) size;
            XTENSOR_ASSERT_MSG(c.size() == size, "Trying to resize const data container with wrong size.");
        }
    }

    /**
     * resizes the container.
     * @param shape the new shape
     * @param force force reshaping, even if the shape stays the same (default: false)
     */
    template <class D>
    template <class S>
    inline void xstrided_container<D>::resize(S&& shape, bool force)
    {
        if (m_shape.size() != shape.size() || !std::equal(std::begin(shape), std::end(shape), std::begin(m_shape)) || force)
        {
            if (m_layout == layout_type::dynamic || m_layout == layout_type::any)
            {
                m_layout = XTENSOR_DEFAULT_LAYOUT;  // fall back to default layout
            }
            m_shape = xtl::forward_sequence<shape_type>(shape);
            resize_container(m_strides, m_shape.size());
            resize_container(m_backstrides, m_shape.size());
            size_type data_size = compute_strides(m_shape, m_layout, m_strides, m_backstrides);
            detail::resize_data_container(this->storage(), data_size);
        }
    }

    /**
     * resizes the container.
     * @param shape the new shape
     * @param l the new layout_type
     */
    template <class D>
    template <class S>
    inline void xstrided_container<D>::resize(S&& shape, layout_type l)
    {
        if (base_type::static_layout != layout_type::dynamic && l != base_type::static_layout)
        {
            throw std::runtime_error("Cannot change layout_type if template parameter not layout_type::dynamic.");
        }
        m_layout = l;
        resize(std::forward<S>(shape), true);
    }

    /**
     * Resizes the container.
     * @param shape the new shape
     * @param strides the new strides
     */
    template <class D>
    template <class S>
    inline void xstrided_container<D>::resize(S&& shape, const strides_type& strides)
    {
        if (base_type::static_layout != layout_type::dynamic)
        {
            throw std::runtime_error("Cannot resize with custom strides when layout() is != layout_type::dynamic.");
        }
        m_shape = xtl::forward_sequence<shape_type>(shape);
        m_strides = strides;
        resize_container(m_backstrides, m_strides.size());
        adapt_strides(m_shape, m_strides, m_backstrides);
        m_layout = layout_type::dynamic;
        detail::resize_data_container(this->storage(), compute_size(m_shape));
    }

    /**
     * Reshapes the container and keeps old elements
     * @param shape the new shape (has to have same number of elements as the original container)
     * @param layout the layout to compute the strides (defaults to static layout of the container,
     *               or for a container with dynamic layout to XTENSOR_DEFAULT_LAYOUT)
     */
    template <class D>
    template <class S>
    inline void xstrided_container<D>::reshape(S&& shape, layout_type layout)
    {
        if (compute_size(shape) != this->size())
        {
            throw std::runtime_error("Cannot reshape with incorrect number of elements. Do you mean to resize?");
        }
        if (layout == layout_type::dynamic || layout == layout_type::any)
        {
            layout = XTENSOR_DEFAULT_LAYOUT;  // fall back to default layout
        }
        if (layout != base_type::static_layout && base_type::static_layout != layout_type::dynamic)
        {
            throw std::runtime_error("Cannot reshape with different layout if static layout != dynamic.");
        }
        m_layout = layout;
        m_shape = xtl::forward_sequence<shape_type>(shape);
        resize_container(m_strides, m_shape.size());
        resize_container(m_backstrides, m_shape.size());
        compute_strides(m_shape, m_layout, m_strides, m_backstrides);
    }
}

#endif
