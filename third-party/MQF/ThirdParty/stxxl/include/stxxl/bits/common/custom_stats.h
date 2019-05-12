/***************************************************************************
 *  include/stxxl/bits/common/custom_stats.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014 Thomas Keh <thomas.keh@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_CUSTOM_STATS_HEADER
#define STXXL_COMMON_CUSTOM_STATS_HEADER

#include <string>
#include <utility>
#include <algorithm>
#include <stxxl/bits/common/timer.h>
#include <stxxl/bits/common/utils.h>

STXXL_BEGIN_NAMESPACE

/*!
 * This class provides a statistical counter that can easily be deactivated
 * using a typedef to dummy_custom_stats_counter.  It's basically a wrapper for
 * a unsigned long long value.
 *
 * \see dummy_custom_stats_counter
 */
template <typename ValueType>
class custom_stats_counter
{
public:
    //! The counter's value type
    typedef ValueType value_type;

protected:
    //! The counter's value
    value_type m_value;

public:
    //! The constructor. Initializes the counter to 0.
    custom_stats_counter()
        : m_value(0)
    { }

    //! Increases the counter by right.
    //! \param right The corresponding integer value
    custom_stats_counter& operator += (const value_type& right)
    {
        m_value += right;
        return *this;
    }
    //! Increases the counter by 1 (prefix).
    custom_stats_counter& operator ++ ()
    {
        ++m_value;
        return *this;
    }
    //! Increases the counter by 1 (postfix).
    custom_stats_counter operator ++ (int)
    {
        custom_stats_counter copy = *this;
        ++m_value;
        return copy;
    }
    //! Assignment operator
    //! \param other The corresponding integer value
    custom_stats_counter& operator = (const value_type& other)
    {
        m_value = other;
        return *this;
    }
    /*!
     * Set the counter to other if other is larger than the current counter
     * value.
     *
     * \param other The corresponding integer value
     */
    void set_max(const value_type& other)
    {
        m_value = std::max(m_value, other);
    }
    /*!
     * Return the counter value interpreted as a memory amount in IEC units as
     * a string.  For that purpose the counter value is multiplied with the
     * byte_per_element argument.
     *
     * \param byte_per_element The memory amount per "counter element".
     */
    std::string as_memory_amount(const value_type& byte_per_element) const
    {
        return format_IEC_size(m_value * byte_per_element) + "B";
    }
    /*!
     * Cast to counter_type: Returns the counter's value as a regular integer
     * value.  This can be used as a getter as well as for printing with
     * std::out.
     */
    operator value_type () const
    {
        return m_value;
    }
};

/*!
 * Dummy class for custom_stats_counter. The methods do nothing.  The compiler
 * should optimize out the code.
 *
 * \see custom_stats_counter
 */
template <typename ValueType>
class dummy_custom_stats_counter
{
public:
    typedef ValueType value_type;

public:
    dummy_custom_stats_counter() { }
    dummy_custom_stats_counter& operator += (const value_type&)
    {
        return *this;
    }
    dummy_custom_stats_counter& operator ++ ()
    {
        return *this;
    }
    dummy_custom_stats_counter& operator ++ (int)
    {
        return *this;
    }
    dummy_custom_stats_counter& operator = (const value_type&)
    {
        return *this;
    }
    void set_max(value_type) { }
    std::string as_memory_amount(const value_type&) const
    {
        return "";
    }
    operator value_type () const
    {
        return value_type();
    }
};

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_CUSTOM_STATS_HEADER
