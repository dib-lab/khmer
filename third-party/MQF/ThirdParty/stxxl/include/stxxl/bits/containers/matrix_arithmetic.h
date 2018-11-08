/***************************************************************************
 *  include/stxxl/bits/containers/matrix_arithmetic.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2010-2011 Raoul Steffen <R-Steffen@gmx.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_MATRIX_ARITHMETIC_HEADER
#define STXXL_CONTAINERS_MATRIX_ARITHMETIC_HEADER

#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/containers/matrix_low_level.h>

STXXL_BEGIN_NAMESPACE

#ifndef STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_MAX_NUM_LEVELS
#define STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_MAX_NUM_LEVELS 3
#endif
#ifndef STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_BASE_CASE
#define STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_BASE_CASE 2
#endif

template <typename ValueType>
class column_vector;

template <typename ValueType>
class row_vector;

template <typename ValueType, unsigned BlockSideLength>
class swappable_block_matrix;

//! \addtogroup matrix
//! \{

struct matrix_operation_statistic_dataset
{
    int_type block_multiplication_calls,
        block_multiplications_saved_through_zero,
        block_addition_calls,
        block_additions_saved_through_zero;

    matrix_operation_statistic_dataset()
        : block_multiplication_calls(0),
          block_multiplications_saved_through_zero(0),
          block_addition_calls(0),
          block_additions_saved_through_zero(0) { }

    matrix_operation_statistic_dataset operator + (const matrix_operation_statistic_dataset& stat)
    {
        matrix_operation_statistic_dataset res(*this);
        res.block_multiplication_calls += stat.block_multiplication_calls;
        res.block_multiplications_saved_through_zero += stat.block_multiplications_saved_through_zero;
        res.block_addition_calls += stat.block_addition_calls;
        res.block_additions_saved_through_zero += stat.block_additions_saved_through_zero;
        return res;
    }

    matrix_operation_statistic_dataset operator - (const matrix_operation_statistic_dataset& stat)
    {
        matrix_operation_statistic_dataset res(*this);
        res.block_multiplication_calls -= stat.block_multiplication_calls;
        res.block_multiplications_saved_through_zero -= stat.block_multiplications_saved_through_zero;
        res.block_addition_calls -= stat.block_addition_calls;
        res.block_additions_saved_through_zero -= stat.block_additions_saved_through_zero;
        return res;
    }
};

struct matrix_operation_statistic
    : public singleton<matrix_operation_statistic>, public matrix_operation_statistic_dataset
{ };

struct matrix_operation_statistic_data : public matrix_operation_statistic_dataset
{
    matrix_operation_statistic_data(const matrix_operation_statistic& stat = * matrix_operation_statistic::get_instance())
        : matrix_operation_statistic_dataset(stat) { }

    matrix_operation_statistic_data(const matrix_operation_statistic_dataset& stat)
        : matrix_operation_statistic_dataset(stat) { }

    matrix_operation_statistic_data& operator = (const matrix_operation_statistic& stat)
    {
        return *this = matrix_operation_statistic_data(stat);
    }

    void set()
    { operator = (*matrix_operation_statistic::get_instance()); }

    matrix_operation_statistic_data operator + (const matrix_operation_statistic_data& stat)
    { return matrix_operation_statistic_data(matrix_operation_statistic_dataset(*this) + matrix_operation_statistic_dataset(stat)); }

    matrix_operation_statistic_data operator - (const matrix_operation_statistic_data& stat)
    { return matrix_operation_statistic_data(matrix_operation_statistic_dataset(*this) - matrix_operation_statistic_dataset(stat)); }
};

std::ostream& operator << (std::ostream& o, const matrix_operation_statistic_data& statsd)
{
    o << "matrix operation statistics" << std::endl;
    o << "block multiplication calls                     : "
      << statsd.block_multiplication_calls << std::endl;
    o << "block multiplications saved through zero blocks: "
      << statsd.block_multiplications_saved_through_zero << std::endl;
    o << "block multiplications performed                : "
      << statsd.block_multiplication_calls - statsd.block_multiplications_saved_through_zero << std::endl;
    o << "block addition calls                           : "
      << statsd.block_addition_calls << std::endl;
    o << "block additions saved through zero blocks      : "
      << statsd.block_additions_saved_through_zero << std::endl;
    o << "block additions performed                      : "
      << statsd.block_addition_calls - statsd.block_additions_saved_through_zero << std::endl;
    return o;
}

//! \}

//! matrix low-level operations and tools
namespace matrix_local {

//! A static_quadtree holds 4^Level elements arranged in a quad tree.
//!
//! Static quad trees are useful for recursive algorithms with fixed depth
//!   that partition the in- and output and perform pre- and postcalculations on the partitions.
//! The four children of one node are denoted as ul (up left), ur (up right), dl (down left), and dr (down right).
template <typename ValueType, unsigned Level>
struct static_quadtree
{
    typedef static_quadtree<ValueType, Level - 1> smaller_static_quadtree;

    smaller_static_quadtree ul, ur, dl, dr;

    static_quadtree(smaller_static_quadtree ul, smaller_static_quadtree ur,
                    smaller_static_quadtree dl, smaller_static_quadtree dr)
        : ul(ul), ur(ur), dl(dl), dr(dr) { }

    static_quadtree() { }

    static_quadtree& operator &= (const static_quadtree& right)
    {
        ul &= right.ul, ur &= right.ur;
        dl &= right.dl, dr &= right.dr;
        return *this;
    }

    static_quadtree& operator += (const static_quadtree& right)
    {
        ul += right.ul, ur += right.ur;
        dl += right.dl, dr += right.dr;
        return *this;
    }

    static_quadtree& operator -= (const static_quadtree& right)
    {
        ul -= right.ul, ur -= right.ur;
        dl -= right.dl, dr -= right.dr;
        return *this;
    }

    static_quadtree operator & (const static_quadtree& right) const
    { return static_quadtree(ul & right.ul, ur & right.ur, dl & right.dl, dr & right.dr); }

    static_quadtree operator + (const static_quadtree& right) const
    { return static_quadtree(ul + right.ul, ur + right.ur, dl + right.dl, dr + right.dr); }

    static_quadtree operator - (const static_quadtree& right) const
    { return static_quadtree(ul - right.ul, ur - right.ur, dl - right.dl, dr - right.dr); }
};

template <typename ValueType>
struct static_quadtree<ValueType, 0>
{
    ValueType val;

    static_quadtree(const ValueType& v)
        : val(v) { }

    static_quadtree() { }

    operator const ValueType& () const
    { return val; }

    operator ValueType& ()
    { return val; }

    static_quadtree& operator &= (const static_quadtree& right)
    {
        val &= right.val;
        return *this;
    }

    static_quadtree& operator += (const static_quadtree& right)
    {
        val += right.val;
        return *this;
    }

    static_quadtree& operator -= (const static_quadtree& right)
    {
        val -= right.val;
        return *this;
    }

    static_quadtree operator ! () const
    { return static_quadtree(! val); }

    static_quadtree operator & (const static_quadtree& right) const
    { return val & right.val; }

    static_quadtree operator + (const static_quadtree& right) const
    { return val + right.val; }

    static_quadtree operator - (const static_quadtree& right) const
    { return val - right.val; }
};

template <typename ValueType, unsigned BlockSideLength, unsigned Level, bool AExists, bool BExists>
struct feedable_strassen_winograd
{
    typedef static_quadtree<bool, Level> zbt;     // true <=> is a zero-block
    typedef static_quadtree<ValueType, Level> vt;

    typedef feedable_strassen_winograd<ValueType, BlockSideLength, Level - 1, AExists, BExists> smaller_feedable_strassen_winograd_ab;
    typedef feedable_strassen_winograd<ValueType, BlockSideLength, Level - 1, AExists, false> smaller_feedable_strassen_winograd_a;
    typedef feedable_strassen_winograd<ValueType, BlockSideLength, Level - 1, false, BExists> smaller_feedable_strassen_winograd_b;
    typedef feedable_strassen_winograd<ValueType, BlockSideLength, Level - 1, false, false> smaller_feedable_strassen_winograd_n;

    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename swappable_block_matrix_type::size_type size_type;

    const size_type n, m, l;
    smaller_feedable_strassen_winograd_ab p1, p2;
    smaller_feedable_strassen_winograd_n p3, p4, p5;
    smaller_feedable_strassen_winograd_b p6;
    smaller_feedable_strassen_winograd_a p7;

    feedable_strassen_winograd(
        const swappable_block_matrix_type& existing_a, const size_type a_from_row, const size_type a_from_col,
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l,
        const swappable_block_matrix_type& existing_b, const size_type b_from_row, const size_type b_from_col)
        : n(n), m(m), l(l),
          p1(existing_a, a_from_row,       a_from_col,       bs_c, n/2, m/2, l/2, existing_b, b_from_row,       b_from_col),
          p2(existing_a, a_from_row,       a_from_col + l/2, bs_c, n/2, m/2, l/2, existing_b, b_from_row + l/2, b_from_col),
          p3(                                                bs_c, n/2, m/2, l/2),
          p4(                                                bs_c, n/2, m/2, l/2),
          p5(                                                bs_c, n/2, m/2, l/2),
          p6(                                                bs_c, n/2, m/2, l/2, existing_b, b_from_row + l/2, b_from_col + m/2),
          p7(existing_a, a_from_row + n/2, a_from_col + l/2, bs_c, n/2, m/2, l/2) {}

    feedable_strassen_winograd(
        const swappable_block_matrix_type& existing_a, const size_type a_from_row, const size_type a_from_col,
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l)
        : n(n), m(m), l(l),
          p1(existing_a, a_from_row,       a_from_col,       bs_c, n/2, m/2, l/2),
          p2(existing_a, a_from_row,       a_from_col + l/2, bs_c, n/2, m/2, l/2),
          p3(                                                bs_c, n/2, m/2, l/2),
          p4(                                                bs_c, n/2, m/2, l/2),
          p5(                                                bs_c, n/2, m/2, l/2),
          p6(                                                bs_c, n/2, m/2, l/2),
          p7(existing_a, a_from_row + n/2, a_from_col + l/2, bs_c, n/2, m/2, l/2) {}

    feedable_strassen_winograd(
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l,
        const swappable_block_matrix_type& existing_b, const size_type b_from_row, const size_type b_from_col)
        : n(n), m(m), l(l),
          p1(bs_c, n/2, m/2, l/2, existing_b, b_from_row,       b_from_col),
          p2(bs_c, n/2, m/2, l/2, existing_b, b_from_row + l/2, b_from_col),
          p3(bs_c, n/2, m/2, l/2),
          p4(bs_c, n/2, m/2, l/2),
          p5(bs_c, n/2, m/2, l/2),
          p6(bs_c, n/2, m/2, l/2, existing_b, b_from_row + l/2, b_from_col + m/2),
          p7(bs_c, n/2, m/2, l/2) {}

    feedable_strassen_winograd(
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l)
        : n(n), m(m), l(l),
          p1(bs_c, n / 2, m / 2, l / 2),
          p2(bs_c, n / 2, m / 2, l / 2),
          p3(bs_c, n / 2, m / 2, l / 2),
          p4(bs_c, n / 2, m / 2, l / 2),
          p5(bs_c, n / 2, m / 2, l / 2),
          p6(bs_c, n / 2, m / 2, l / 2),
          p7(bs_c, n / 2, m / 2, l / 2) { }

    void begin_feeding_a_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        typename zbt::smaller_static_quadtree
        s1 = zb.dl & zb.dr,
        s2 = s1 & zb.ul,
        s3 = zb.ul & zb.dl,
        s4 = zb.ur & s2;
        p1.begin_feeding_a_block(block_row, block_col, zb.ul);
        p2.begin_feeding_a_block(block_row, block_col, zb.ur);
        p3.begin_feeding_a_block(block_row, block_col, s1);
        p4.begin_feeding_a_block(block_row, block_col, s2);
        p5.begin_feeding_a_block(block_row, block_col, s3);
        p6.begin_feeding_a_block(block_row, block_col, s4);
        p7.begin_feeding_a_block(block_row, block_col, zb.dr);
    }

    void feed_a_element(const int_type element_num, const vt v)
    {
        typename vt::smaller_static_quadtree
        s1 = v.dl + v.dr,
        s2 = s1 - v.ul,
        s3 = v.ul - v.dl,
        s4 = v.ur - s2;
        p1.feed_a_element(element_num, v.ul);
        p2.feed_a_element(element_num, v.ur);
        p3.feed_a_element(element_num, s1);
        p4.feed_a_element(element_num, s2);
        p5.feed_a_element(element_num, s3);
        p6.feed_a_element(element_num, s4);
        p7.feed_a_element(element_num, v.dr);
    }

    void end_feeding_a_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        typename zbt::smaller_static_quadtree
        s1 = zb.dl & zb.dr,
        s2 = s1 & zb.ul,
        s3 = zb.ul & zb.dl,
        s4 = zb.ur & s2;
        p1.end_feeding_a_block(block_row, block_col, zb.ul);
        p2.end_feeding_a_block(block_row, block_col, zb.ur);
        p3.end_feeding_a_block(block_row, block_col, s1);
        p4.end_feeding_a_block(block_row, block_col, s2);
        p5.end_feeding_a_block(block_row, block_col, s3);
        p6.end_feeding_a_block(block_row, block_col, s4);
        p7.end_feeding_a_block(block_row, block_col, zb.dr);
    }

    void begin_feeding_b_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        typename zbt::smaller_static_quadtree
        t1 = zb.ur & zb.ul,
        t2 = zb.dr & t1,
        t3 = zb.dr & zb.ur,
        t4 = zb.dl & t2;
        p1.begin_feeding_b_block(block_row, block_col, zb.ul);
        p2.begin_feeding_b_block(block_row, block_col, zb.dl);
        p3.begin_feeding_b_block(block_row, block_col, t1);
        p4.begin_feeding_b_block(block_row, block_col, t2);
        p5.begin_feeding_b_block(block_row, block_col, t3);
        p6.begin_feeding_b_block(block_row, block_col, zb.dr);
        p7.begin_feeding_b_block(block_row, block_col, t4);
    }

    void feed_b_element(const int_type element_num, const vt v)
    {
        typename vt::smaller_static_quadtree
        t1 = v.ur - v.ul,
        t2 = v.dr - t1,
        t3 = v.dr - v.ur,
        t4 = v.dl - t2;
        p1.feed_b_element(element_num, v.ul);
        p2.feed_b_element(element_num, v.dl);
        p3.feed_b_element(element_num, t1);
        p4.feed_b_element(element_num, t2);
        p5.feed_b_element(element_num, t3);
        p6.feed_b_element(element_num, v.dr);
        p7.feed_b_element(element_num, t4);
    }

    void end_feeding_b_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        typename zbt::smaller_static_quadtree
        t1 = zb.ur & zb.ul,
        t2 = zb.dr & t1,
        t3 = zb.dr & zb.ur,
        t4 = zb.dl & t2;
        p1.end_feeding_b_block(block_row, block_col, zb.ul);
        p2.end_feeding_b_block(block_row, block_col, zb.dl);
        p3.end_feeding_b_block(block_row, block_col, t1);
        p4.end_feeding_b_block(block_row, block_col, t2);
        p5.end_feeding_b_block(block_row, block_col, t3);
        p6.end_feeding_b_block(block_row, block_col, zb.dr);
        p7.end_feeding_b_block(block_row, block_col, t4);
    }

    void multiply()
    {
        p1.multiply();
        p2.multiply();
        p3.multiply();
        p4.multiply();
        p5.multiply();
        p6.multiply();
        p7.multiply();
    }

    zbt begin_reading_block(const size_type& block_row, const size_type& block_col)
    {
        zbt r;
        r.ur = r.ul = p1.begin_reading_block(block_row, block_col);
        r.ul &= p2.begin_reading_block(block_row, block_col);
        r.ur &= p4.begin_reading_block(block_row, block_col);
        r.dr = r.dl = p5.begin_reading_block(block_row, block_col);
        r.dl &= r.ur;
        r.dl &= p7.begin_reading_block(block_row, block_col);
        r.ur &= p3.begin_reading_block(block_row, block_col);
        r.dr &= r.ur;
        r.ur &= p6.begin_reading_block(block_row, block_col);
        return r;
    }

    vt read_element(int_type element_num)
    {
        vt r;
        r.ur = r.ul = p1.read_element(element_num);
        r.ul += p2.read_element(element_num);
        r.ur += p4.read_element(element_num);
        r.dr = r.dl = p5.read_element(element_num);
        r.dl += r.ur;
        r.dl += p7.read_element(element_num);
        r.ur += p3.read_element(element_num);
        r.dr += r.ur;
        r.ur += p6.read_element(element_num);
        return r;
    }

    zbt end_reading_block(const size_type& block_row, const size_type& block_col)
    {
        zbt r;
        r.ur = r.ul = p1.end_reading_block(block_row, block_col);
        r.ul &= p2.end_reading_block(block_row, block_col);
        r.ur &= p4.end_reading_block(block_row, block_col);
        r.dr = r.dl = p5.end_reading_block(block_row, block_col);
        r.dl &= r.ur;
        r.dl &= p7.end_reading_block(block_row, block_col);
        r.ur &= p3.end_reading_block(block_row, block_col);
        r.dr &= r.ur;
        r.ur &= p6.end_reading_block(block_row, block_col);
        return r;
    }
};

template <typename ValueType, unsigned BlockSideLength, bool AExists, bool BExists>
struct feedable_strassen_winograd<ValueType, BlockSideLength, 0, AExists, BExists>
{
    typedef static_quadtree<bool, 0> zbt;     // true <=> is a zero-block
    typedef static_quadtree<ValueType, 0> vt;

    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename swappable_block_matrix_type::size_type size_type;

    swappable_block_matrix_type a, b, c;
    const size_type n, m, l;
    internal_block_type* iblock;

    feedable_strassen_winograd(
        const swappable_block_matrix_type& existing_a, const size_type a_from_row, const size_type a_from_col,
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l,
        const swappable_block_matrix_type& existing_b, const size_type b_from_row, const size_type b_from_col)
        : a(existing_a, n, l, a_from_row, a_from_col),
          b(existing_b, n, l, b_from_row, b_from_col),
          c(bs_c, n, m),
          n(n), m(m), l(l),
          iblock(0) { }

    feedable_strassen_winograd(
        const swappable_block_matrix_type& existing_a, const size_type a_from_row, const size_type a_from_col,
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l)
        : a(existing_a, n, l, a_from_row, a_from_col),
          b(bs_c, n, l),
          c(bs_c, n, m),
          n(n), m(m), l(l),
          iblock(0) { }

    feedable_strassen_winograd(
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l,
        const swappable_block_matrix_type& existing_b, const size_type b_from_row, const size_type b_from_col)
        : a(bs_c, n, l),
          b(existing_b, n, l, b_from_row, b_from_col),
          c(bs_c, n, m),
          n(n), m(m), l(l),
          iblock(0) { }

    feedable_strassen_winograd(
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l)
        : a(bs_c, n, l),
          b(bs_c, n, l),
          c(bs_c, n, m),
          n(n), m(m), l(l),
          iblock(0) { }

    void begin_feeding_a_block(const size_type& block_row, const size_type& block_col, const zbt)
    {
        if (! AExists)
            iblock = &a.bs.acquire(a(block_row, block_col), true);
    }

    void feed_a_element(const int_type element_num, const vt v)
    {
        if (! AExists)
            (*iblock)[element_num] = v;
    }

    void end_feeding_a_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        if (! AExists)
        {
            a.bs.release(a(block_row, block_col), ! zb);
            iblock = 0;
        }
    }

    void begin_feeding_b_block(const size_type& block_row, const size_type& block_col, const zbt)
    {
        if (! BExists)
            iblock = &b.bs.acquire(b(block_row, block_col), true);
    }

    void feed_b_element(const int_type element_num, const vt v)
    {
        if (! BExists)
            (*iblock)[element_num] = v;
    }

    void end_feeding_b_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        if (! BExists)
        {
            b.bs.release(b(block_row, block_col), ! zb);
            iblock = 0;
        }
    }

    void multiply()
    { matrix_operations<ValueType, BlockSideLength>::choose_level_for_feedable_sw(a, b, c); }

    zbt begin_reading_block(const size_type& block_row, const size_type& block_col)
    {
        bool zb = ! c.bs.is_initialized(c(block_row, block_col));
        iblock = &c.bs.acquire(c(block_row, block_col));
        return zb;
    }

    vt read_element(const int_type element_num)
    { return (*iblock)[element_num]; }

    zbt end_reading_block(const size_type& block_row, const size_type& block_col)
    {
        c.bs.release(c(block_row, block_col), false);
        iblock = 0;
        return ! c.bs.is_initialized(c(block_row, block_col));
    }
};

template <typename ValueType, unsigned BlockSideLength, unsigned Level>
struct matrix_to_quadtree
{
    typedef static_quadtree<bool, Level> zbt;     // true <=> is a zero-block
    typedef static_quadtree<ValueType, Level> vt;

    typedef matrix_to_quadtree<ValueType, BlockSideLength, Level - 1> smaller_matrix_to_quadtree;

    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename swappable_block_matrix_type::size_type size_type;

    smaller_matrix_to_quadtree ul, ur, dl, dr;

    matrix_to_quadtree(const swappable_block_matrix_type & matrix)
        : ul(matrix, matrix.get_height()/2, matrix.get_width()/2,                     0,                    0),
          ur(matrix, matrix.get_height()/2, matrix.get_width()/2,                     0, matrix.get_width()/2),
          dl(matrix, matrix.get_height()/2, matrix.get_width()/2, matrix.get_height()/2,                    0),
          dr(matrix, matrix.get_height()/2, matrix.get_width()/2, matrix.get_height()/2, matrix.get_width()/2)
    { assert(! (matrix.get_height() % 2 | matrix.get_width() % 2)); }

    matrix_to_quadtree(const swappable_block_matrix_type & matrix,
            const size_type height, const size_type width, const size_type from_row, const size_type from_col)
        : ul(matrix, height/2, width/2, from_row,            from_col),
          ur(matrix, height/2, width/2, from_row,            from_col + width/2),
          dl(matrix, height/2, width/2, from_row + height/2, from_col),
          dr(matrix, height/2, width/2, from_row + height/2, from_col + width/2)
    { assert(! (height % 2 | width % 2)); }

    void begin_feeding_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        ul.begin_feeding_block(block_row, block_col, zb.ul);
        ur.begin_feeding_block(block_row, block_col, zb.ur);
        dl.begin_feeding_block(block_row, block_col, zb.dl);
        dr.begin_feeding_block(block_row, block_col, zb.dr);
    }

    void feed_element(const int_type element_num, const vt v)
    {
        ul.feed_element(element_num, v.ul);
        ur.feed_element(element_num, v.ur);
        dl.feed_element(element_num, v.dl);
        dr.feed_element(element_num, v.dr);
    }

    void feed_and_add_element(const int_type element_num, const vt v)
    {
        ul.feed_and_add_element(element_num, v.ul);
        ur.feed_and_add_element(element_num, v.ur);
        dl.feed_and_add_element(element_num, v.dl);
        dr.feed_and_add_element(element_num, v.dr);
    }

    void end_feeding_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        ul.end_feeding_block(block_row, block_col, zb.ul);
        ur.end_feeding_block(block_row, block_col, zb.ur);
        dl.end_feeding_block(block_row, block_col, zb.dl);
        dr.end_feeding_block(block_row, block_col, zb.dr);
    }

    zbt begin_reading_block(const size_type& block_row, const size_type& block_col)
    {
        zbt zb;
        zb.ul = ul.begin_reading_block(block_row, block_col);
        zb.ur = ur.begin_reading_block(block_row, block_col);
        zb.dl = dl.begin_reading_block(block_row, block_col);
        zb.dr = dr.begin_reading_block(block_row, block_col);
        return zb;
    }

    vt read_element(const int_type element_num)
    {
        vt v;
        v.ul = ul.read_element(element_num);
        v.ur = ur.read_element(element_num);
        v.dl = dl.read_element(element_num);
        v.dr = dr.read_element(element_num);
        return v;
    }

    zbt end_reading_block(const size_type& block_row, const size_type& block_col)
    {
        zbt zb;
        zb.ul = ul.end_reading_block(block_row, block_col);
        zb.ur = ur.end_reading_block(block_row, block_col);
        zb.dl = dl.end_reading_block(block_row, block_col);
        zb.dr = dr.end_reading_block(block_row, block_col);
        return zb;
    }

    const size_type & get_height_in_blocks()
    { return ul.get_height_in_blocks(); }

    const size_type & get_width_in_blocks()
    { return ul.get_width_in_blocks(); }
};

template <typename ValueType, unsigned BlockSideLength>
struct matrix_to_quadtree<ValueType, BlockSideLength, 0>
{
    typedef static_quadtree<bool, 0> zbt;     // true <=> is a zero-block
    typedef static_quadtree<ValueType, 0> vt;

    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename swappable_block_matrix_type::size_type size_type;

    swappable_block_matrix_type m;
    internal_block_type* iblock;

    matrix_to_quadtree(const swappable_block_matrix_type& matrix)
        : m(matrix, matrix.get_height(), matrix.get_width(), 0, 0),
          iblock(0) { }

    matrix_to_quadtree(const swappable_block_matrix_type& matrix,
                       const size_type height, const size_type width, const size_type from_row, const size_type from_col)
        : m(matrix, height, width, from_row, from_col),
          iblock(0) { }

    void begin_feeding_block(const size_type& block_row, const size_type& block_col, const zbt)
    { iblock = &m.bs.acquire(m(block_row, block_col)); }

    void feed_element(const int_type element_num, const vt v)
    { (*iblock)[element_num] = v; }

    void feed_and_add_element(const int_type element_num, const vt v)
    { (*iblock)[element_num] += v; }

    void end_feeding_block(const size_type& block_row, const size_type& block_col, const zbt zb)
    {
        m.bs.release(m(block_row, block_col), ! zb);
        iblock = 0;
    }

    zbt begin_reading_block(const size_type& block_row, const size_type& block_col)
    {
        zbt zb = ! m.bs.is_initialized(m(block_row, block_col));
        iblock = &m.bs.acquire(m(block_row, block_col));
        return zb;
    }

    vt read_element(const int_type element_num)
    { return (*iblock)[element_num]; }

    zbt end_reading_block(const size_type& block_row, const size_type& block_col)
    {
        m.bs.release(m(block_row, block_col), false);
        iblock = 0;
        return ! m.bs.is_initialized(m(block_row, block_col));
    }

    const size_type & get_height_in_blocks()
    { return m.get_height(); }

    const size_type & get_width_in_blocks()
    { return m.get_width(); }
};

template <typename ValueType, unsigned BlockSideLength, unsigned Level, bool AExists, bool BExists>
struct feedable_strassen_winograd_block_grained
{
    typedef static_quadtree<bool, Level> zbt;     // true <=> is a zero-block
    typedef static_quadtree<ValueType, Level> vt;

    typedef feedable_strassen_winograd_block_grained<ValueType, BlockSideLength, Level - 1, AExists, BExists> smaller_feedable_strassen_winograd_ab;
    typedef feedable_strassen_winograd_block_grained<ValueType, BlockSideLength, Level - 1, AExists, false> smaller_feedable_strassen_winograd_a;
    typedef feedable_strassen_winograd_block_grained<ValueType, BlockSideLength, Level - 1, false, BExists> smaller_feedable_strassen_winograd_b;
    typedef feedable_strassen_winograd_block_grained<ValueType, BlockSideLength, Level - 1, false, false> smaller_feedable_strassen_winograd_n;

    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename swappable_block_matrix_type::size_type size_type;
    typedef matrix_operations<ValueType, BlockSideLength> Ops;

    const size_type n, m, l;
    smaller_feedable_strassen_winograd_ab p1, p2;
    smaller_feedable_strassen_winograd_n p3, p4, p5;
    smaller_feedable_strassen_winograd_b p6;
    smaller_feedable_strassen_winograd_a p7;

    inline feedable_strassen_winograd_block_grained(
        const swappable_block_matrix_type& existing_a, const size_type a_from_row, const size_type a_from_col,
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l,
        const swappable_block_matrix_type& existing_b, const size_type b_from_row, const size_type b_from_col)
        : n(n), m(m), l(l),
          p1(existing_a, a_from_row,       a_from_col,       bs_c, n/2, m/2, l/2, existing_b, b_from_row,       b_from_col),
          p2(existing_a, a_from_row,       a_from_col + l/2, bs_c, n/2, m/2, l/2, existing_b, b_from_row + l/2, b_from_col),
          p3(                                                bs_c, n/2, m/2, l/2),
          p4(                                                bs_c, n/2, m/2, l/2),
          p5(                                                bs_c, n/2, m/2, l/2),
          p6(                                                bs_c, n/2, m/2, l/2, existing_b, b_from_row + l/2, b_from_col + m/2),
          p7(existing_a, a_from_row + n/2, a_from_col + l/2, bs_c, n/2, m/2, l/2) {}

    inline feedable_strassen_winograd_block_grained(
        const swappable_block_matrix_type& existing_a, const size_type a_from_row, const size_type a_from_col,
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l)
        : n(n), m(m), l(l),
          p1(existing_a, a_from_row,       a_from_col,       bs_c, n/2, m/2, l/2),
          p2(existing_a, a_from_row,       a_from_col + l/2, bs_c, n/2, m/2, l/2),
          p3(                                                bs_c, n/2, m/2, l/2),
          p4(                                                bs_c, n/2, m/2, l/2),
          p5(                                                bs_c, n/2, m/2, l/2),
          p6(                                                bs_c, n/2, m/2, l/2),
          p7(existing_a, a_from_row + n/2, a_from_col + l/2, bs_c, n/2, m/2, l/2) {}

    inline feedable_strassen_winograd_block_grained(
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l,
        const swappable_block_matrix_type& existing_b, const size_type b_from_row, const size_type b_from_col)
        : n(n), m(m), l(l),
          p1(bs_c, n/2, m/2, l/2, existing_b, b_from_row,       b_from_col),
          p2(bs_c, n/2, m/2, l/2, existing_b, b_from_row + l/2, b_from_col),
          p3(bs_c, n/2, m/2, l/2),
          p4(bs_c, n/2, m/2, l/2),
          p5(bs_c, n/2, m/2, l/2),
          p6(bs_c, n/2, m/2, l/2, existing_b, b_from_row + l/2, b_from_col + m/2),
          p7(bs_c, n/2, m/2, l/2) {}

    inline feedable_strassen_winograd_block_grained(
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l)
        : n(n), m(m), l(l),
          p1(bs_c, n / 2, m / 2, l / 2),
          p2(bs_c, n / 2, m / 2, l / 2),
          p3(bs_c, n / 2, m / 2, l / 2),
          p4(bs_c, n / 2, m / 2, l / 2),
          p5(bs_c, n / 2, m / 2, l / 2),
          p6(bs_c, n / 2, m / 2, l / 2),
          p7(bs_c, n / 2, m / 2, l / 2) { }

    inline void feed_a(const size_type& row, const size_type& col, const swappable_block_matrix_type& bl)
    {
        // partition bl
        typename Ops::swappable_block_matrix_quarterer qbl(bl);
        // preadditions
        swappable_block_matrix_type
            s1(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed()),
        s2(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed()),
        s3(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed()),
        s4(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed());
        Ops::strassen_winograd_preaddition_a(qbl.ul, qbl.ur, qbl.dl, qbl.dr, s1, s2, s3, s4);
        // feed recursive
        p1.feed_a(row, col, qbl.ul);
        p2.feed_a(row, col, qbl.ur);
        p3.feed_a(row, col, s1);
        p4.feed_a(row, col, s2);
        p5.feed_a(row, col, s3);
        p6.feed_a(row, col, s4);
        p7.feed_a(row, col, qbl.dr);
    }

    inline void feed_b(const size_type& row, const size_type& col, const swappable_block_matrix_type& bl)
    {
        // partition bl
        typename Ops::swappable_block_matrix_quarterer qbl(bl);
        // preadditions
        swappable_block_matrix_type
            t1(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed()),
        t2(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed()),
        t3(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed()),
        t4(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed());
        Ops::strassen_winograd_preaddition_b(qbl.ul, qbl.ur, qbl.dl, qbl.dr, t1, t2, t3, t4);
        // feed recursive
        p1.feed_b(row, col, qbl.ul);
        p2.feed_b(row, col, qbl.dl);
        p3.feed_b(row, col, t1);
        p4.feed_b(row, col, t2);
        p5.feed_b(row, col, t3);
        p6.feed_b(row, col, qbl.dr);
        p7.feed_b(row, col, t4);
    }

    inline void multiply()
    {
        p1.multiply();
        p2.multiply();
        p3.multiply();
        p4.multiply();
        p5.multiply();
        p6.multiply();
        p7.multiply();
    }

    inline void read_and_add(const size_type& row, const size_type& col, const swappable_block_matrix_type& bl)
    {
        // partition bl
        typename Ops::swappable_block_matrix_quarterer qbl(bl);
        // postadditions
        swappable_block_matrix_type px(bl.bs, qbl.ul.get_height(), qbl.ul.get_width(), qbl.ul.is_transposed());
        p2.read_and_add(row, col, qbl.ul);
        p1.read_and_add(row, col, px);
        Ops::element_op(qbl.ul, px, typename Ops::addition());
        p4.read_and_add(row, col, px);
        Ops::element_op(qbl.ur, px, typename Ops::addition());
        p5.read_and_add(row, col, px);
        Ops::element_op_twice_nontransposed(qbl.dl, qbl.dr, px, typename Ops::addition());
        px.set_zero();
        p7.read_and_add(row, col, qbl.dl);
        p3.read_and_add(row, col, px);
        Ops::element_op_twice_nontransposed(qbl.dr, qbl.ur, px, typename Ops::addition());
        p6.read_and_add(row, col, qbl.ur);
    }

    inline static unsigned_type get_num_temp_grains()
    { return smaller_feedable_strassen_winograd_ab::get_num_temp_grains() + (4 ^ Level) * 2; }
};

template <typename ValueType, unsigned BlockSideLength, bool AExists, bool BExists>
struct feedable_strassen_winograd_block_grained<ValueType, BlockSideLength, 0, AExists, BExists>
{
    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename swappable_block_matrix_type::swappable_block_identifier_type swappable_block_identifier_type;
    typedef typename swappable_block_matrix_type::size_type size_type;
    typedef matrix_operations<ValueType, BlockSideLength> Ops;

    typedef static_quadtree<swappable_block_identifier_type, 0> bt;

    swappable_block_matrix_type a, b, c;

    inline feedable_strassen_winograd_block_grained(
        const swappable_block_matrix_type& existing_a, const size_type a_from_row, const size_type a_from_col,
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l,
        const swappable_block_matrix_type& existing_b, const size_type b_from_row, const size_type b_from_col)
        : a(existing_a, n, l, a_from_row, a_from_col),
          b(existing_b, n, l, b_from_row, b_from_col),
          c(bs_c, n, m) { }

    inline feedable_strassen_winograd_block_grained(
        const swappable_block_matrix_type& existing_a, const size_type a_from_row, const size_type a_from_col,
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l)
        : a(existing_a, n, l, a_from_row, a_from_col),
          b(bs_c, n, l),
          c(bs_c, n, m) { }

    inline feedable_strassen_winograd_block_grained(
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l,
        const swappable_block_matrix_type& existing_b, const size_type b_from_row, const size_type b_from_col)
        : a(bs_c, n, l),
          b(existing_b, n, l, b_from_row, b_from_col),
          c(bs_c, n, m) { }

    inline feedable_strassen_winograd_block_grained(
        block_scheduler_type& bs_c, const size_type n, const size_type m, const size_type l)
        : a(bs_c, n, l),
          b(bs_c, n, l),
          c(bs_c, n, m) { }

    inline void feed_a(const size_type& row, const size_type& col, const swappable_block_matrix_type& bl)
    {
        if (! AExists)
        {
            // copy bl to a from (row, col) (assuming a from (row, col) == 0)
            swappable_block_matrix_type at(a, bl.get_height(), bl.get_width(), row, col);
            Ops::element_op(at, bl, typename Ops::addition());
        }
    }

    inline void feed_b(const size_type& row, const size_type& col, const swappable_block_matrix_type& bl)
    {
        if (! BExists)
        {
            // copy bl(0,0) to b(row, col) (assuming b from (row, col) == 0)
            swappable_block_matrix_type bt(b, bl.get_height(), bl.get_width(), row, col);
            Ops::element_op(bt, bl, typename Ops::addition());
        }
    }

    inline void multiply()
    {
        matrix_operations<ValueType, BlockSideLength>::
        multi_level_strassen_winograd_multiply_and_add_block_grained(a, b, c);
        if (! AExists)
            a.set_zero();
        if (! BExists)
            b.set_zero();
    }


    inline void read_and_add(const size_type& row, const size_type& col, swappable_block_matrix_type& bl)
    {
        // add c from (row, col) to bl
        swappable_block_matrix_type ct(c, bl.get_height(), bl.get_width(), row, col);
        Ops::element_op(bl, ct, typename Ops::addition());
        ct.set_zero();
    }

    inline static unsigned_type get_num_temp_grains()
    { return 0; }
};

template <typename ValueType, unsigned BlockSideLength, unsigned Level, unsigned Granularity>
struct matrix_to_quadtree_block_grained
{
    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::size_type size_type;

    typedef matrix_to_quadtree_block_grained<ValueType, BlockSideLength, Level - 1, Granularity> smaller_matrix_to_quadtree_block_grained;

    smaller_matrix_to_quadtree_block_grained ul, ur, dl, dr;

    inline matrix_to_quadtree_block_grained(const swappable_block_matrix_type & matrix)
        : ul(matrix, matrix.get_height()/2, matrix.get_width()/2,                     0,                    0),
          ur(matrix, matrix.get_height()/2, matrix.get_width()/2,                     0, matrix.get_width()/2),
          dl(matrix, matrix.get_height()/2, matrix.get_width()/2, matrix.get_height()/2,                    0),
          dr(matrix, matrix.get_height()/2, matrix.get_width()/2, matrix.get_height()/2, matrix.get_width()/2)
    { assert(! (matrix.get_height() % 2 | matrix.get_width() % 2)); }

    inline matrix_to_quadtree_block_grained(const swappable_block_matrix_type & matrix,
            const size_type height, const size_type width, const size_type from_row, const size_type from_col)
        : ul(matrix, height/2, width/2, from_row,            from_col),
          ur(matrix, height/2, width/2, from_row,            from_col + width/2),
          dl(matrix, height/2, width/2, from_row + height/2, from_col),
          dr(matrix, height/2, width/2, from_row + height/2, from_col + width/2)
    { assert(! (height % 2 | width % 2)); }

    inline swappable_block_matrix_type operator () (const size_type& row, const size_type& col)
    {
        return swappable_block_matrix_type(ul(row, col), ur(row, col), dl(row, col), dr(row, col));
    }

    inline const size_type get_height()
    { return ul.get_height(); }

    inline const size_type get_width()
    { return ul.get_width(); }
};

template <typename ValueType, unsigned BlockSideLength, unsigned Granularity>
struct matrix_to_quadtree_block_grained<ValueType, BlockSideLength, 0, Granularity>
{
    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::size_type size_type;

    swappable_block_matrix_type m;

    inline matrix_to_quadtree_block_grained(const swappable_block_matrix_type& matrix)
        : m(matrix, matrix.get_height(), matrix.get_width(), 0, 0)
    { assert(! (matrix.get_height() % Granularity | matrix.get_width() % Granularity)); }

    inline matrix_to_quadtree_block_grained(const swappable_block_matrix_type& matrix,
                                            const size_type height, const size_type width, const size_type from_row, const size_type from_col)
        : m(matrix, height, width, from_row, from_col)
    { assert(! (matrix.get_height() % Granularity | matrix.get_width() % Granularity)); }

    inline swappable_block_matrix_type operator () (const size_type& row, const size_type& col)
    {
        return swappable_block_matrix_type(m, Granularity, Granularity, row * Granularity, col * Granularity);
    }

    inline const size_type get_height()
    { return m.get_height() / Granularity; }

    inline const size_type get_width()
    { return m.get_width() / Granularity; }
};

template <typename ValueType, unsigned BlockSideLength>
struct matrix_operations
{
    // tuning-parameter: Only matrices larger than this (in blocks) are processed by Strassen-Winograd.
    // you have to adapt choose_level_for_feedable_sw, too
    static const int_type strassen_winograd_base_case_size;

    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename swappable_block_matrix_type::swappable_block_identifier_type swappable_block_identifier_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename swappable_block_matrix_type::size_type size_type;
    typedef column_vector<ValueType> column_vector_type;
    typedef row_vector<ValueType> row_vector_type;
    typedef typename column_vector_type::size_type vector_size_type;

    // +-+-+-+ addition +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    struct addition
    {
        /* op(c,a,b) means c = a <op> b  e.g. assign sum
         * op(c,a)   means c <op>= a     e.g. add up
         * op(a)     means <op>a         e.g. sign
         *
         * it should hold:
         * op(c,0,0) equivalent c = 0
         * op(c=0,a) equivalent c = op(a)
         * op(c,0)   equivalent {}
         */

        inline ValueType& operator () (ValueType& c, const ValueType& a, const ValueType& b) { return c = a + b; }
        inline ValueType& operator () (ValueType& c, const ValueType& a) { return c += a; }
        inline ValueType operator () (const ValueType& a) { return +a; }
    };

    struct subtraction
    {
        inline ValueType& operator () (ValueType& c, const ValueType& a, const ValueType& b) { return c = a - b; }
        inline ValueType& operator () (ValueType& c, const ValueType& a) { return c -= a; }
        inline ValueType operator () (const ValueType& a) { return -a; }
    };

    struct scalar_multiplication
    {
        inline scalar_multiplication(const ValueType scalar = 1) : s(scalar) { }
        inline ValueType& operator () (ValueType& c, const ValueType& a) { return c = a * s; }
        inline ValueType operator () (const ValueType& a) { return a * s; }
        inline operator const ValueType& () { return s; }
        const ValueType s;
    };

    // element_op<Op>(C,A,B) calculates C = A <Op> B
    template <class Op>
    static swappable_block_matrix_type&
    element_op(swappable_block_matrix_type& C,
               const swappable_block_matrix_type& A,
               const swappable_block_matrix_type& B, Op op = Op())
    {
        for (size_type row = 0; row < C.get_height(); ++row)
            for (size_type col = 0; col < C.get_width(); ++col)
                element_op_swappable_block(
                    C(row, col), C.is_transposed(), C.bs,
                    A(row, col), A.is_transposed(), A.bs,
                    B(row, col), B.is_transposed(), B.bs, op);
        return C;
    }

    // element_op<Op>(C,A) calculates C <Op>= A
    template <class Op>
    static swappable_block_matrix_type&
    element_op(swappable_block_matrix_type& C,
               const swappable_block_matrix_type& A, Op op = Op())
    {
        for (size_type row = 0; row < C.get_height(); ++row)
            for (size_type col = 0; col < C.get_width(); ++col)
                element_op_swappable_block(
                    C(row, col), C.is_transposed(), C.bs,
                    A(row, col), A.is_transposed(), A.bs, op);
        return C;
    }

    // element_op<Op>(C) calculates C = <Op>C
    template <class Op>
    static swappable_block_matrix_type&
    element_op(swappable_block_matrix_type& C, Op op = Op())
    {
        for (size_type row = 0; row < C.get_height(); ++row)
            for (size_type col = 0; col < C.get_width(); ++col)
                element_op_swappable_block(
                    C(row, col), C.bs, op);
        return C;
    }

    // calculates c = a <Op> b
    template <class Op>
    static void
    element_op_swappable_block(
        const swappable_block_identifier_type c, const bool c_is_transposed, block_scheduler_type& bs_c,
        const swappable_block_identifier_type a, bool a_is_transposed, block_scheduler_type& bs_a,
        const swappable_block_identifier_type b, bool b_is_transposed, block_scheduler_type& bs_b, Op op = Op())
    {
        if (! bs_c.is_simulating())
            ++matrix_operation_statistic::get_instance()->block_addition_calls;
        // check if zero-block (== ! initialized)
        if (! bs_a.is_initialized(a) && ! bs_b.is_initialized(b))
        {
            // => a and b are zero -> set c zero
            bs_c.deinitialize(c);
            if (! bs_c.is_simulating())
                ++matrix_operation_statistic::get_instance()->block_additions_saved_through_zero;
            return;
        }
        a_is_transposed = a_is_transposed != c_is_transposed;
        b_is_transposed = b_is_transposed != c_is_transposed;
        if (! bs_a.is_initialized(a))
        {
            // a is zero -> copy b
            internal_block_type& ic = bs_c.acquire(c, true),
            & ib = bs_b.acquire(b);
            if (! bs_c.is_simulating())
            {
                if (b_is_transposed)
                    low_level_matrix_binary_ass_op<ValueType, BlockSideLength, false, true, Op>(&ic[0], 0, &ib[0], op);
                else
                    low_level_matrix_binary_ass_op<ValueType, BlockSideLength, false, false, Op>(&ic[0], 0, &ib[0], op);
            }
            bs_b.release(b, false);
            bs_c.release(c, true);
        }
        else if (! bs_b.is_initialized(b))
        {
            // b is zero -> copy a
            internal_block_type& ic = bs_c.acquire(c, true),
            & ia = bs_a.acquire(a);
            if (! bs_c.is_simulating())
            {
                if (a_is_transposed)
                    low_level_matrix_binary_ass_op<ValueType, BlockSideLength, true, false, Op>(&ic[0], &ia[0], 0, op);
                else
                    low_level_matrix_binary_ass_op<ValueType, BlockSideLength, false, false, Op>(&ic[0], &ia[0], 0, op);
            }
            bs_a.release(a, false);
            bs_c.release(c, true);
        }
        else
        {
            internal_block_type& ic = bs_c.acquire(c, true),
            & ia = bs_a.acquire(a),
            & ib = bs_b.acquire(b);
            if (! bs_c.is_simulating())
            {
                if (a_is_transposed)
                {
                    if (b_is_transposed)
                        low_level_matrix_binary_ass_op<ValueType, BlockSideLength, true, true, Op>(&ic[0], &ia[0], &ib[0], op);
                    else
                        low_level_matrix_binary_ass_op<ValueType, BlockSideLength, true, false, Op>(&ic[0], &ia[0], &ib[0], op);
                }
                else
                {
                    if (b_is_transposed)
                        low_level_matrix_binary_ass_op<ValueType, BlockSideLength, false, true, Op>(&ic[0], &ia[0], &ib[0], op);
                    else
                        low_level_matrix_binary_ass_op<ValueType, BlockSideLength, false, false, Op>(&ic[0], &ia[0], &ib[0], op);
                }
            }
            bs_a.release(a, false);
            bs_b.release(b, false);
            bs_c.release(c, true);
        }
    }

    // calculates c <op>= a
    template <class Op>
    static void
    element_op_swappable_block(
        const swappable_block_identifier_type c, const bool c_is_transposed, block_scheduler_type& bs_c,
        const swappable_block_identifier_type a, const bool a_is_transposed, block_scheduler_type& bs_a, Op op = Op())
    {
        if (! bs_c.is_simulating())
            ++matrix_operation_statistic::get_instance()->block_addition_calls;
        // check if zero-block (== ! initialized)
        if (! bs_a.is_initialized(a))
        {
            // => b is zero => nothing to do
            if (! bs_c.is_simulating())
                ++matrix_operation_statistic::get_instance()->block_additions_saved_through_zero;
            return;
        }
        const bool c_is_zero = ! bs_c.is_initialized(c);
        // acquire
        internal_block_type& ic = bs_c.acquire(c, c_is_zero),
        & ia = bs_a.acquire(a);
        // add
        if (! bs_c.is_simulating())
        {
            if (c_is_zero) {
                if (c_is_transposed == a_is_transposed)
                    low_level_matrix_unary_op<ValueType, BlockSideLength, false, Op>(&ic[0], &ia[0], op);
                else
                    low_level_matrix_unary_op<ValueType, BlockSideLength, true, Op>(&ic[0], &ia[0], op);
            }
            else {
                if (c_is_transposed == a_is_transposed)
                    low_level_matrix_unary_ass_op<ValueType, BlockSideLength, false, Op>(&ic[0], &ia[0], op);
                else
                    low_level_matrix_unary_ass_op<ValueType, BlockSideLength, true, Op>(&ic[0], &ia[0], op);
            }
        }
        // release
        bs_c.release(c, true);
        bs_a.release(a, false);
    }

    // calculates c = <op>c
    template <class Op>
    static void
    element_op_swappable_block(
        const swappable_block_identifier_type c, block_scheduler_type& bs_c, Op op = Op())
    {
        if (! bs_c.is_simulating())
            ++matrix_operation_statistic::get_instance()->block_addition_calls;
        // check if zero-block (== ! initialized)
        if (! bs_c.is_initialized(c))
        {
            // => c is zero => nothing to do
            if (! bs_c.is_simulating())
                ++matrix_operation_statistic::get_instance()->block_additions_saved_through_zero;
            return;
        }
        // acquire
        internal_block_type& ic = bs_c.acquire(c);
        // add
        if (! bs_c.is_simulating())
            low_level_matrix_unary_op<ValueType, BlockSideLength, false, Op>(&ic[0], &ic[0], op);
        // release
        bs_c.release(c, true);
    }

    // additions for strassen-winograd

    inline static void
    strassen_winograd_preaddition_a(swappable_block_matrix_type& a11,
                                    swappable_block_matrix_type& a12,
                                    swappable_block_matrix_type& a21,
                                    swappable_block_matrix_type& a22,
                                    swappable_block_matrix_type& s1,
                                    swappable_block_matrix_type& s2,
                                    swappable_block_matrix_type& s3,
                                    swappable_block_matrix_type& s4)
    {
        for (size_type row = 0; row < a11.get_height(); ++row)
            for (size_type col = 0; col < a11.get_width(); ++col)
            {
                op_swappable_block_nontransposed(s3, a11, subtraction(), a21, row, col);
                op_swappable_block_nontransposed(s1, a21,    addition(), a22, row, col);
                op_swappable_block_nontransposed(s2,  s1, subtraction(), a11, row, col);
                op_swappable_block_nontransposed(s4, a12, subtraction(),  s2, row, col);
            }
    }

    inline static void
    strassen_winograd_preaddition_b(swappable_block_matrix_type& b11,
                                    swappable_block_matrix_type& b12,
                                    swappable_block_matrix_type& b21,
                                    swappable_block_matrix_type& b22,
                                    swappable_block_matrix_type& t1,
                                    swappable_block_matrix_type& t2,
                                    swappable_block_matrix_type& t3,
                                    swappable_block_matrix_type& t4)
    {
        for (size_type row = 0; row < b11.get_height(); ++row)
            for (size_type col = 0; col < b11.get_width(); ++col)
            {
                op_swappable_block_nontransposed(t3, b22, subtraction(), b12, row, col);
                op_swappable_block_nontransposed(t1, b12, subtraction(), b11, row, col);
                op_swappable_block_nontransposed(t2, b22, subtraction(),  t1, row, col);
                op_swappable_block_nontransposed(t4, b21, subtraction(),  t2, row, col);
            }
    }

    inline static void
    strassen_winograd_postaddition(swappable_block_matrix_type& c11,  // = p2
                                   swappable_block_matrix_type& c12,  // = p6
                                   swappable_block_matrix_type& c21,  // = p7
                                   swappable_block_matrix_type& c22,  // = p4
                                   swappable_block_matrix_type& p1,
                                   swappable_block_matrix_type& p3,
                                   swappable_block_matrix_type& p5)
    {
        for (size_type row = 0; row < c11.get_height(); ++row)
            for (size_type col = 0; col < c11.get_width(); ++col)
            {
                op_swappable_block_nontransposed(c11,     addition(),  p1, row, col); // (u1)
                op_swappable_block_nontransposed( p1,     addition(), c22, row, col); // (u2)
                op_swappable_block_nontransposed( p5,     addition(),  p1, row, col); // (u3)
                op_swappable_block_nontransposed(c21,     addition(),  p5, row, col); // (u4)
                op_swappable_block_nontransposed(c22, p5, addition(),  p3, row, col); // (u5)
                op_swappable_block_nontransposed( p1,     addition(),  p3, row, col); // (u6)
                op_swappable_block_nontransposed(c12,     addition(),  p1, row, col); // (u7)
            }
    }

    // calculates c1 += a; c2 += a
    template <class Op>
    inline static void
    element_op_twice_nontransposed(swappable_block_matrix_type& c1,
                                   swappable_block_matrix_type& c2,
                                   const swappable_block_matrix_type& a, Op op = Op())
    {
        for (size_type row = 0; row < a.get_height(); ++row)
            for (size_type col = 0; col < a.get_width(); ++col)
            {
                element_op_swappable_block(
                    c1(row, col), false, c1.bs,
                    a(row, col), false, a.bs, op);
                element_op_swappable_block(
                    c2(row, col), false, c2.bs,
                    a(row, col), false, a.bs, op);
            }
    }

    template <class Op>
    inline static void
    op_swappable_block_nontransposed(swappable_block_matrix_type& c,
                                     swappable_block_matrix_type& a, Op op, swappable_block_matrix_type& b,
                                     size_type& row, size_type& col)
    {
        element_op_swappable_block(
            c(row, col), false, c.bs,
            a(row, col), false, a.bs,
            b(row, col), false, b.bs, op);
    }

    template <class Op>
    inline static void
    op_swappable_block_nontransposed(swappable_block_matrix_type& c, Op op, swappable_block_matrix_type& a,
                                     size_type& row, size_type& col)
    {
        element_op_swappable_block(
            c(row, col), false, c.bs,
            a(row, col), false, a.bs, op);
    }

    // +-+ end addition +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    // +-+-+-+ matrix multiplication +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    /*  n, m and l denote the three dimensions of a matrix multiplication, according to the following ascii-art diagram:
     *
     *                 +--m--+
     *  +----l-----+   |     |   +--m--+
     *  |          |   |     |   |     |
     *  n    A     | • l  B  | = n  C  |
     *  |          |   |     |   |     |
     *  +----------+   |     |   +-----+
     *                 +-----+
     *
     * The index-variables are called i, j, k for dimension
     *                                n, m, l .
     */

    // requires height and width divisible by 2
    struct swappable_block_matrix_quarterer
    {
        swappable_block_matrix_type  upleft,   upright,
                                    downleft, downright,
                & ul, & ur, & dl, & dr;

        swappable_block_matrix_quarterer(const swappable_block_matrix_type & whole)
            : upleft   (whole, whole.get_height()/2, whole.get_width()/2,                    0,                   0),
              upright  (whole, whole.get_height()/2, whole.get_width()/2,                    0, whole.get_width()/2),
              downleft (whole, whole.get_height()/2, whole.get_width()/2, whole.get_height()/2,                   0),
              downright(whole, whole.get_height()/2, whole.get_width()/2, whole.get_height()/2, whole.get_width()/2),
              ul(upleft), ur(upright), dl(downleft), dr(downright)
        { assert(! (whole.get_height() % 2 | whole.get_width() % 2)); }
    };

    struct swappable_block_matrix_padding_quarterer
    {
        swappable_block_matrix_type  upleft,   upright,
                                    downleft, downright,
                & ul, & ur, & dl, & dr;

        swappable_block_matrix_padding_quarterer(const swappable_block_matrix_type & whole)
            : upleft   (whole, div_ceil(whole.get_height(),2), div_ceil(whole.get_width(),2),                              0,                             0),
              upright  (whole, div_ceil(whole.get_height(),2), div_ceil(whole.get_width(),2),                              0, div_ceil(whole.get_width(),2)),
              downleft (whole, div_ceil(whole.get_height(),2), div_ceil(whole.get_width(),2), div_ceil(whole.get_height(),2),                             0),
              downright(whole, div_ceil(whole.get_height(),2), div_ceil(whole.get_width(),2), div_ceil(whole.get_height(),2), div_ceil(whole.get_width(),2)),
              ul(upleft), ur(upright), dl(downleft), dr(downright) {}
    };

    struct swappable_block_matrix_approximative_quarterer
    {
        swappable_block_matrix_type upleft, upright,
            downleft, downright,
        & ul, & ur, & dl, & dr;

        swappable_block_matrix_approximative_quarterer(const swappable_block_matrix_type & whole)
            : upleft   (whole,                      whole.get_height()/2,                     whole.get_width()/2,                    0,                   0),
              upright  (whole,                      whole.get_height()/2, whole.get_width() - whole.get_width()/2,                    0, whole.get_width()/2),
              downleft (whole, whole.get_height() - whole.get_height()/2,                     whole.get_width()/2, whole.get_height()/2,                   0),
              downright(whole, whole.get_height() - whole.get_height()/2, whole.get_width() - whole.get_width()/2, whole.get_height()/2, whole.get_width()/2),
              ul(upleft), ur(upright), dl(downleft), dr(downright) {}
    };

    //! calculates C = A * B + C
    // requires fitting dimensions
    static swappable_block_matrix_type&
    multi_level_strassen_winograd_multiply_and_add_block_grained(const swappable_block_matrix_type& A,
                                                                 const swappable_block_matrix_type& B,
                                                                 swappable_block_matrix_type& C)
    {
        int_type num_levels = ilog2_ceil(std::min(A.get_width(), std::min(C.get_width(), C.get_height())));
        if (num_levels > STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_BASE_CASE)
        {
            if (num_levels > STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_MAX_NUM_LEVELS)
                num_levels = STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_MAX_NUM_LEVELS;
            swappable_block_matrix_type padded_a(A, round_up_to_power_of_two(A.get_height(), num_levels),
                                                 round_up_to_power_of_two(A.get_width(), num_levels), 0, 0),
            padded_b(B, round_up_to_power_of_two(B.get_height(), num_levels),
                     round_up_to_power_of_two(B.get_width(), num_levels), 0, 0),
            padded_c(C, round_up_to_power_of_two(C.get_height(), num_levels),
                     round_up_to_power_of_two(C.get_width(), num_levels), 0, 0);
            switch (num_levels)
            {
            #if (STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_MAX_NUM_LEVELS >= 5 && 5 > STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_BASE_CASE)
            case 5:
                use_feedable_sw_block_grained<5>(padded_a, padded_a, padded_c);
                break;
            #endif
            #if (STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_MAX_NUM_LEVELS >= 4 && 4 > STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_BASE_CASE)
            case 4:
                use_feedable_sw_block_grained<4>(padded_a, padded_a, padded_c);
                break;
            #endif
            #if (STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_MAX_NUM_LEVELS >= 3 && 3 > STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_BASE_CASE)
            case 3:
                use_feedable_sw_block_grained<3>(padded_a, padded_a, padded_c);
                break;
            #endif
            #if (STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_MAX_NUM_LEVELS >= 2 && 2 > STXXL_MATRIX_MULTI_LEVEL_STRASSEN_WINOGRAD_BASE_CASE)
            case 2:
                use_feedable_sw_block_grained<2>(padded_a, padded_a, padded_c);
                break;
            #endif
            default:  // only here in case of wrong bounds
                strassen_winograd_multiply_and_add_interleaved(A, B, C);
                break;
            }
        }
        else
            // base case
            strassen_winograd_multiply_and_add_interleaved(A, B, C);
        return C;
    }

    // input matrices have to be padded
    template <unsigned Level>
    static void use_feedable_sw_block_grained(const swappable_block_matrix_type& A,
                                              const swappable_block_matrix_type& B,
                                              swappable_block_matrix_type& C)
    {
        const unsigned granularity = 1;

        feedable_strassen_winograd_block_grained<ValueType, BlockSideLength, Level, true, true>
        fsw(A, 0, 0, C.bs, C.get_height(), C.get_width(), A.get_width(), B, 0, 0);
        // preadditions for A
        {
            matrix_to_quadtree_block_grained<ValueType, BlockSideLength, Level, granularity>
            mtq_a(A);
            for (size_type row = 0; row < mtq_a.get_height(); ++row)
                for (size_type col = 0; col < mtq_a.get_width(); ++col)
                    fsw.feed_a(row, col, mtq_a(row, col));
        }
        // preadditions for B
        {
            matrix_to_quadtree_block_grained<ValueType, BlockSideLength, Level, granularity>
            mtq_b(B);
            for (size_type row = 0; row < mtq_b.get_height(); ++row)
                for (size_type col = 0; col < mtq_b.get_width(); ++col)
                    fsw.feed_b(row, col, mtq_b(row, col));
        }
        // recursive multiplications
        fsw.multiply();
        // postadditions
        {
            matrix_to_quadtree_block_grained<ValueType, BlockSideLength, Level, granularity>
            mtq_c(C);
            for (size_type row = 0; row < mtq_c.get_height(); ++row)
                for (size_type col = 0; col < mtq_c.get_width(); ++col)
                    fsw.read_and_add(row, col, mtq_c(row, col));
        }
    }

    //! calculates C = A * B + C
    // requires fitting dimensions
    static swappable_block_matrix_type&
    multi_level_strassen_winograd_multiply_and_add(const swappable_block_matrix_type& A,
                                                   const swappable_block_matrix_type& B,
                                                   swappable_block_matrix_type& C)
    {
        int_type p = ilog2_ceil(std::min(A.get_width(), std::min(C.get_width(), C.get_height())));

        swappable_block_matrix_type padded_a(A, round_up_to_power_of_two(A.get_height(), p),
                                             round_up_to_power_of_two(A.get_width(), p), 0, 0),
        padded_b(B, round_up_to_power_of_two(B.get_height(), p),
                 round_up_to_power_of_two(B.get_width(), p), 0, 0),
        padded_c(C, round_up_to_power_of_two(C.get_height(), p),
                 round_up_to_power_of_two(C.get_width(), p), 0, 0);
        choose_level_for_feedable_sw(padded_a, padded_b, padded_c);
        return C;
    }

    // input matrices have to be padded
    static void choose_level_for_feedable_sw(const swappable_block_matrix_type& A,
                                             const swappable_block_matrix_type& B,
                                             swappable_block_matrix_type& C)
    {
        switch (ilog2_ceil(std::min(A.get_width(), std::min(C.get_width(), C.get_height()))))
        {
        default:
            /*
            use_feedable_sw<4>(A, B, C);
            break;
        case 3:
            use_feedable_sw<3>(A, B, C);
            break;
        case 2:*/
            use_feedable_sw<2>(A, B, C);
            break;
        case 1:
        /*use_feedable_sw<1>(A, B, C);
        break;*/
        case 0:
            // base case
            recursive_multiply_and_add(A, B, C);
            break;
        }
    }

    // input matrices have to be padded
    template <unsigned Level>
    static void use_feedable_sw(const swappable_block_matrix_type& A,
                                const swappable_block_matrix_type& B,
                                swappable_block_matrix_type& C)
    {
        feedable_strassen_winograd<ValueType, BlockSideLength, Level, true, true>
        fsw(A, 0, 0, C.bs, C.get_height(), C.get_width(), A.get_width(), B, 0, 0);
        // preadditions for A
        matrix_to_quadtree<ValueType, BlockSideLength, Level>
        mtq_a(A);
        for (size_type block_row = 0; block_row < mtq_a.get_height_in_blocks(); ++block_row)
            for (size_type block_col = 0; block_col < mtq_a.get_width_in_blocks(); ++block_col)
            {
                fsw.begin_feeding_a_block(block_row, block_col,
                                          mtq_a.begin_reading_block(block_row, block_col));
                #if STXXL_PARALLEL
                #pragma omp parallel for
                #endif
                for (int_type element_row_in_block = 0; element_row_in_block < int_type(BlockSideLength); ++element_row_in_block)
                    for (int_type element_col_in_block = 0; element_col_in_block < int_type(BlockSideLength); ++element_col_in_block)
                        fsw.feed_a_element(element_row_in_block * BlockSideLength + element_col_in_block,
                                           mtq_a.read_element(element_row_in_block * BlockSideLength + element_col_in_block));
                fsw.end_feeding_a_block(block_row, block_col,
                                        mtq_a.end_reading_block(block_row, block_col));
            }
        // preadditions for B
        matrix_to_quadtree<ValueType, BlockSideLength, Level>
        mtq_b(B);
        for (size_type block_row = 0; block_row < mtq_b.get_height_in_blocks(); ++block_row)
            for (size_type block_col = 0; block_col < mtq_b.get_width_in_blocks(); ++block_col)
            {
                fsw.begin_feeding_b_block(block_row, block_col,
                                          mtq_b.begin_reading_block(block_row, block_col));
                #if STXXL_PARALLEL
                #pragma omp parallel for
                #endif
                for (int_type element_row_in_block = 0; element_row_in_block < int_type(BlockSideLength); ++element_row_in_block)
                    for (int_type element_col_in_block = 0; element_col_in_block < int_type(BlockSideLength); ++element_col_in_block)
                        fsw.feed_b_element(element_row_in_block * BlockSideLength + element_col_in_block,
                                           mtq_b.read_element(element_row_in_block * BlockSideLength + element_col_in_block));
                fsw.end_feeding_b_block(block_row, block_col,
                                        mtq_b.end_reading_block(block_row, block_col));
            }
        // recursive multiplications
        fsw.multiply();
        // postadditions
        matrix_to_quadtree<ValueType, BlockSideLength, Level>
        mtq_c(C);
        for (size_type block_row = 0; block_row < mtq_c.get_height_in_blocks(); ++block_row)
            for (size_type block_col = 0; block_col < mtq_c.get_width_in_blocks(); ++block_col)
            {
                mtq_c.begin_feeding_block(block_row, block_col,
                                          fsw.begin_reading_block(block_row, block_col));
                #if STXXL_PARALLEL
                #pragma omp parallel for
                #endif
                for (int_type element_row_in_block = 0; element_row_in_block < int_type(BlockSideLength); ++element_row_in_block)
                    for (int_type element_col_in_block = 0; element_col_in_block < int_type(BlockSideLength); ++element_col_in_block)
                        mtq_c.feed_and_add_element(element_row_in_block * BlockSideLength + element_col_in_block,
                                                   fsw.read_element(element_row_in_block * BlockSideLength + element_col_in_block));
                mtq_c.end_feeding_block(block_row, block_col,
                                        fsw.end_reading_block(block_row, block_col));
            }
    }

    //! calculates C = A * B
    // assumes fitting dimensions
    static swappable_block_matrix_type&
    strassen_winograd_multiply(const swappable_block_matrix_type& A,
                               const swappable_block_matrix_type& B,
                               swappable_block_matrix_type& C)
    {
        // base case
        if (C.get_height() <= strassen_winograd_base_case_size
            || C.get_width() <= strassen_winograd_base_case_size
            || A.get_width() <= strassen_winograd_base_case_size)
        {
            C.set_zero();
            return recursive_multiply_and_add(A, B, C);
        }

        // partition matrix
        swappable_block_matrix_padding_quarterer qa(A), qb(B), qc(C);
        // preadditions
        swappable_block_matrix_type s1(C.bs, qa.ul.get_height(), qa.ul.get_width(), qa.ul.is_transposed()),
                                    s2(C.bs, qa.ul.get_height(), qa.ul.get_width(), qa.ul.is_transposed()),
                                    s3(C.bs, qa.ul.get_height(), qa.ul.get_width(), qa.ul.is_transposed()),
                                    s4(C.bs, qa.ul.get_height(), qa.ul.get_width(), qa.ul.is_transposed()),
                                    t1(C.bs, qb.ul.get_height(), qb.ul.get_width(), qb.ul.is_transposed()),
                                    t2(C.bs, qb.ul.get_height(), qb.ul.get_width(), qb.ul.is_transposed()),
                                    t3(C.bs, qb.ul.get_height(), qb.ul.get_width(), qb.ul.is_transposed()),
                                    t4(C.bs, qb.ul.get_height(), qb.ul.get_width(), qb.ul.is_transposed());
        strassen_winograd_preaddition_a(qa.ul, qa.ur, qa.dl, qa.dr, s1, s2, s3, s4);
        strassen_winograd_preaddition_b(qb.ul, qb.ur, qb.dl, qb.dr, t1, t2, t3, t4);
        // recursive multiplications
        swappable_block_matrix_type p1(C.bs, qc.ul.get_height(), qc.ul.get_width(), qc.ul.is_transposed()),
                                 // p2 stored in qc.ul
                                    p3(C.bs, qc.ul.get_height(), qc.ul.get_width(), qc.ul.is_transposed()),
                                 // p4 stored in qc.dr
                                    p5(C.bs, qc.ul.get_height(), qc.ul.get_width(), qc.ul.is_transposed());
                                 // p6 stored in qc.ur
                                 // p7 stored in qc.dl
        strassen_winograd_multiply(qa.ul, qb.ul,    p1);
        strassen_winograd_multiply(qa.ur, qb.dl, qc.ul);
        strassen_winograd_multiply(   s1,    t1,    p3);
        strassen_winograd_multiply(   s2,    t2, qc.dr);
        strassen_winograd_multiply(   s3,    t3,    p5);
        strassen_winograd_multiply(   s4, qb.dr, qc.ur);
        strassen_winograd_multiply(qa.dr,    t4, qc.dl);
        // postadditions
        strassen_winograd_postaddition(qc.ul, qc.ur, qc.dl, qc.dr, p1, p3, p5);
        return C;
    }

    //! calculates C = A * B + C
    // assumes fitting dimensions
    static swappable_block_matrix_type&
    strassen_winograd_multiply_and_add_interleaved(const swappable_block_matrix_type& A,
                                                   const swappable_block_matrix_type& B,
                                                   swappable_block_matrix_type& C)
    {
        // base case
        if (C.get_height() <= strassen_winograd_base_case_size
            || C.get_width() <= strassen_winograd_base_case_size
            || A.get_width() <= strassen_winograd_base_case_size)
            return recursive_multiply_and_add(A, B, C);

        // partition matrix
        swappable_block_matrix_padding_quarterer qa(A), qb(B), qc(C);
        // preadditions
        swappable_block_matrix_type s1(C.bs, qa.ul.get_height(), qa.ul.get_width(), qa.ul.is_transposed()),
                                    s2(C.bs, qa.ul.get_height(), qa.ul.get_width(), qa.ul.is_transposed()),
                                    s3(C.bs, qa.ul.get_height(), qa.ul.get_width(), qa.ul.is_transposed()),
                                    s4(C.bs, qa.ul.get_height(), qa.ul.get_width(), qa.ul.is_transposed()),
                                    t1(C.bs, qb.ul.get_height(), qb.ul.get_width(), qb.ul.is_transposed()),
                                    t2(C.bs, qb.ul.get_height(), qb.ul.get_width(), qb.ul.is_transposed()),
                                    t3(C.bs, qb.ul.get_height(), qb.ul.get_width(), qb.ul.is_transposed()),
                                    t4(C.bs, qb.ul.get_height(), qb.ul.get_width(), qb.ul.is_transposed());
        strassen_winograd_preaddition_a(qa.ul, qa.ur, qa.dl, qa.dr, s1, s2, s3, s4);
        strassen_winograd_preaddition_b(qb.ul, qb.ur, qb.dl, qb.dr, t1, t2, t3, t4);
        // recursive multiplications and postadditions
        swappable_block_matrix_type px(C.bs, qc.ul.get_height(), qc.ul.get_width(), qc.ul.is_transposed());
        strassen_winograd_multiply_and_add_interleaved(qa.ur, qb.dl, qc.ul); // p2
        strassen_winograd_multiply_and_add_interleaved(qa.ul, qb.ul, px); // p1
        element_op<addition>(qc.ul, px);
        strassen_winograd_multiply_and_add_interleaved(s2, t2, px); // p4
        s2.set_zero();
        t2.set_zero();
        element_op<addition>(qc.ur, px);
        strassen_winograd_multiply_and_add_interleaved(s3, t3, px); // p5
        s3.set_zero();
        t3.set_zero();
        element_op_twice_nontransposed<addition>(qc.dl, qc.dr, px);
        px.set_zero();
        strassen_winograd_multiply_and_add_interleaved(qa.dr, t4, qc.dl); // p7
        t4.set_zero();
        strassen_winograd_multiply_and_add_interleaved(s1, t1, px); // p3
        s1.set_zero();
        t1.set_zero();
        element_op_twice_nontransposed<addition>(qc.dr, qc.ur, px);
        px.set_zero();
        strassen_winograd_multiply_and_add_interleaved(s4, qb.dr, qc.ur); // p6
        return C;
    }

    //! calculates C = A * B + C
    // assumes fitting dimensions
    static swappable_block_matrix_type&
    strassen_winograd_multiply_and_add(const swappable_block_matrix_type& A,
                                       const swappable_block_matrix_type& B,
                                       swappable_block_matrix_type& C)
    {
        // base case
        if (C.get_height() <= strassen_winograd_base_case_size
            || C.get_width() <= strassen_winograd_base_case_size
            || A.get_width() <= strassen_winograd_base_case_size)
            return recursive_multiply_and_add(A, B, C);

        // partition matrix
        swappable_block_matrix_padding_quarterer qa(A), qb(B), qc(C);
        // preadditions
        swappable_block_matrix_type s1(C.bs, qa.ul.get_height(), qa.ul.get_width()),
        s2(C.bs, qa.ul.get_height(), qa.ul.get_width()),
        s3(C.bs, qa.ul.get_height(), qa.ul.get_width()),
        s4(C.bs, qa.ul.get_height(), qa.ul.get_width()),
        t1(C.bs, qb.ul.get_height(), qb.ul.get_width()),
        t2(C.bs, qb.ul.get_height(), qb.ul.get_width()),
        t3(C.bs, qb.ul.get_height(), qb.ul.get_width()),
        t4(C.bs, qb.ul.get_height(), qb.ul.get_width());
        element_op<subtraction>(s3, qa.ul, qa.dl);
        element_op<addition>(s1, qa.dl, qa.dr);
        element_op<subtraction>(s2, s1, qa.ul);
        element_op<subtraction>(s4, qa.ur, s2);
        element_op<subtraction>(t3, qb.dr, qb.ur);
        element_op<subtraction>(t1, qb.ur, qb.ul);
        element_op<subtraction>(t2, qb.dr, t1);
        element_op<subtraction>(t4, qb.dl, t2);
        // recursive multiplications and postadditions
        swappable_block_matrix_type px(C.bs, qc.ul.get_height(), qc.ul.get_width());
        strassen_winograd_multiply_and_add(qa.ur, qb.dl, qc.ul); // p2
        strassen_winograd_multiply_and_add(qa.ul, qb.ul, px); // p1
        element_op<addition>(qc.ul, px);
        strassen_winograd_multiply_and_add(s2, t2, px); // p4
        element_op<addition>(qc.ur, px);
        strassen_winograd_multiply_and_add(s3, t3, px); // p5
        element_op<addition>(qc.dl, px);
        element_op<addition>(qc.dr, px);
        px.set_zero();
        strassen_winograd_multiply_and_add(qa.dr, t4, qc.dl); // p7
        strassen_winograd_multiply_and_add(s1, t1, px); // p3
        element_op<addition>(qc.dr, px);
        element_op<addition>(qc.ur, px);
        strassen_winograd_multiply_and_add(s4, qb.dr, qc.ur); // p6
        return C;
    }

    //! calculates C = A * B + C
    // assumes fitting dimensions
    static swappable_block_matrix_type&
    recursive_multiply_and_add(const swappable_block_matrix_type& A,
                               const swappable_block_matrix_type& B,
                               swappable_block_matrix_type& C)
    {
        // catch empty intervals
        if (C.get_height() * C.get_width() * A.get_width() == 0)
            return C;
        // base case
        if ((C.get_height() == 1) + (C.get_width() == 1) + (A.get_width() == 1) >= 2)
            return naive_multiply_and_add(A, B, C);

        // partition matrix
        swappable_block_matrix_approximative_quarterer qa(A), qb(B), qc(C);
        // recursive multiplication
        // The order of recursive calls is optimized to enhance locality. C has priority because it has to be read and written.
        recursive_multiply_and_add(qa.ul, qb.ul, qc.ul);
        recursive_multiply_and_add(qa.ur, qb.dl, qc.ul);
        recursive_multiply_and_add(qa.ur, qb.dr, qc.ur);
        recursive_multiply_and_add(qa.ul, qb.ur, qc.ur);
        recursive_multiply_and_add(qa.dl, qb.ur, qc.dr);
        recursive_multiply_and_add(qa.dr, qb.dr, qc.dr);
        recursive_multiply_and_add(qa.dr, qb.dl, qc.dl);
        recursive_multiply_and_add(qa.dl, qb.ul, qc.dl);

        return C;
    }

    //! calculates C = A * B + C
    // requires fitting dimensions
    static swappable_block_matrix_type&
    naive_multiply_and_add(const swappable_block_matrix_type& A,
                           const swappable_block_matrix_type& B,
                           swappable_block_matrix_type& C)
    {
        const size_type& n = C.get_height(),
        & m = C.get_width(),
        & l = A.get_width();
        for (size_type i = 0; i < n; ++i)
            for (size_type j = 0; j < m; ++j)
                for (size_type k = 0; k < l; ++k)
                    multiply_and_add_swappable_block(A(i, k), A.is_transposed(), A.bs,
                                                     B(k, j), B.is_transposed(), B.bs,
                                                     C(i, j), C.is_transposed(), C.bs);
        return C;
    }

    static void multiply_and_add_swappable_block(
        const swappable_block_identifier_type a, const bool a_is_transposed, block_scheduler_type& bs_a,
        const swappable_block_identifier_type b, const bool b_is_transposed, block_scheduler_type& bs_b,
        const swappable_block_identifier_type c, const bool c_is_transposed, block_scheduler_type& bs_c)
    {
        if (! bs_c.is_simulating())
            ++matrix_operation_statistic::get_instance()->block_multiplication_calls;
        // check if zero-block (== ! initialized)
        if (! bs_a.is_initialized(a) || ! bs_b.is_initialized(b))
        {
            // => one factor is zero => product is zero
            if (! bs_c.is_simulating())
                ++matrix_operation_statistic::get_instance()->block_multiplications_saved_through_zero;
            return;
        }
        // acquire
        ValueType* ap = bs_a.acquire(a).begin(),
        * bp = bs_b.acquire(b).begin(),
        * cp = bs_c.acquire(c).begin();
        // multiply
        if (! bs_c.is_simulating())
            low_level_matrix_multiply_and_add<ValueType, BlockSideLength>
                (ap, a_is_transposed, bp, b_is_transposed, cp, c_is_transposed);
        // release
        bs_a.release(a, false);
        bs_b.release(b, false);
        bs_c.release(c, true);
    }

    // +-+ end matrix multiplication +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    // +-+-+-+ matrix-vector multiplication +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    //! calculates z = A * x
    static column_vector_type&
    recursive_matrix_col_vector_multiply_and_add(const swappable_block_matrix_type& A,
                                                 const column_vector_type& x, column_vector_type& z,
                                                 const vector_size_type offset_x = 0, const vector_size_type offset_z = 0)
    {
        // catch empty intervals
        if (A.get_height() * A.get_width() == 0)
            return z;
        // base case
        if (A.get_height() == 1 || A.get_width() == 1)
            return naive_matrix_col_vector_multiply_and_add(A, x, z, offset_x, offset_z);

        // partition matrix
        swappable_block_matrix_approximative_quarterer qa(A);
        // recursive multiplication
        // The order of recursive calls is optimized to enhance locality.
        recursive_matrix_col_vector_multiply_and_add(qa.ul, x, z, offset_x,                     offset_z                     );
        recursive_matrix_col_vector_multiply_and_add(qa.ur, x, z, offset_x + qa.ul.get_width(), offset_z                     );
        recursive_matrix_col_vector_multiply_and_add(qa.dr, x, z, offset_x + qa.ul.get_width(), offset_z + qa.ul.get_height());
        recursive_matrix_col_vector_multiply_and_add(qa.dl, x, z, offset_x,                     offset_z + qa.ul.get_height());

        return z;
    }

    static column_vector_type&
    naive_matrix_col_vector_multiply_and_add(const swappable_block_matrix_type& A,
                                             const column_vector_type& x, column_vector_type& z,
                                             const vector_size_type offset_x = 0, const vector_size_type offset_z = 0)
    {
        for (size_type row = 0; row < A.get_height(); ++row)
            for (size_type col = 0; col < A.get_width(); ++col)
                matrix_col_vector_multiply_and_add_swappable_block(A(row, col), A.is_transposed(), A.bs,
                                                                   x, z, (offset_x + col) * BlockSideLength, (offset_z + row) * BlockSideLength);
        return z;
    }

    static void matrix_col_vector_multiply_and_add_swappable_block(
        const swappable_block_identifier_type a, const bool a_is_transposed, block_scheduler_type& bs_a,
        const column_vector_type& x, column_vector_type& z,
        const vector_size_type offset_x = 0, const vector_size_type offset_z = 0)
    {
        // check if zero-block (== ! initialized)
        if (! bs_a.is_initialized(a))
        {
            // => matrix is zero => product is zero
            return;
        }
        // acquire
        internal_block_type& ia = bs_a.acquire(a);
        // multiply
        if (! bs_a.is_simulating())
        {
            int_type row_limit = std::min(BlockSideLength, unsigned(z.size() - offset_z)),
                col_limit = std::min(BlockSideLength, unsigned(x.size() - offset_x));
            if (a_is_transposed)
                for (int_type col = 0; col < col_limit; ++col)
                    for (int_type row = 0; row < row_limit; ++row)
                        z[offset_z + row] += x[offset_x + col] * ia[row + col * BlockSideLength];
            else
                for (int_type row = 0; row < row_limit; ++row)
                    for (int_type col = 0; col < col_limit; ++col)
                        z[offset_z + row] += x[offset_x + col] * ia[row * BlockSideLength + col];
        }
        // release
        bs_a.release(a, false);
    }

    //! calculates z = y * A
    static row_vector_type&
    recursive_matrix_row_vector_multiply_and_add(const row_vector_type& y,
                                                 const swappable_block_matrix_type& A, row_vector_type& z,
                                                 const vector_size_type offset_y = 0, const vector_size_type offset_z = 0)
    {
        // catch empty intervals
        if (A.get_height() * A.get_width() == 0)
            return z;
        // base case
        if (A.get_height() == 1 || A.get_width() == 1)
            return naive_matrix_row_vector_multiply_and_add(y, A, z, offset_y, offset_z);

        // partition matrix
        swappable_block_matrix_approximative_quarterer qa(A);
        // recursive multiplication
        // The order of recursive calls is optimized to enhance locality.
        recursive_matrix_row_vector_multiply_and_add(y, qa.ul, z, offset_y,                      offset_z                    );
        recursive_matrix_row_vector_multiply_and_add(y, qa.dl, z, offset_y + qa.ul.get_height(), offset_z                    );
        recursive_matrix_row_vector_multiply_and_add(y, qa.dr, z, offset_y + qa.ul.get_height(), offset_z + qa.ul.get_width());
        recursive_matrix_row_vector_multiply_and_add(y, qa.ur, z, offset_y,                      offset_z + qa.ul.get_width());

        return z;
    }

    static row_vector_type&
    naive_matrix_row_vector_multiply_and_add(const row_vector_type& y, const swappable_block_matrix_type& A,
                                             row_vector_type& z,
                                             const vector_size_type offset_y = 0, const vector_size_type offset_z = 0)
    {
        for (size_type row = 0; row < A.get_height(); ++row)
            for (size_type col = 0; col < A.get_width(); ++col)
                matrix_row_vector_multiply_and_add_swappable_block(y, A(row, col), A.is_transposed(), A.bs,
                                                                   z, (offset_y + row) * BlockSideLength, (offset_z + col) * BlockSideLength);
        return z;
    }

    static void matrix_row_vector_multiply_and_add_swappable_block(const row_vector_type& y,
                                                                   const swappable_block_identifier_type a, const bool a_is_transposed, block_scheduler_type& bs_a,
                                                                   row_vector_type& z,
                                                                   const vector_size_type offset_y = 0, const vector_size_type offset_z = 0)
    {
        // check if zero-block (== ! initialized)
        if (! bs_a.is_initialized(a))
        {
            // => matrix is zero => product is zero
            return;
        }
        // acquire
        internal_block_type& ia = bs_a.acquire(a);
        // multiply
        if (! bs_a.is_simulating())
        {
            int_type row_limit = std::min(BlockSideLength, unsigned(y.size() - offset_y)),
                col_limit = std::min(BlockSideLength, unsigned(z.size() - offset_z));
            if (a_is_transposed)
                for (int_type col = 0; col < col_limit; ++col)
                    for (int_type row = 0; row < row_limit; ++row)
                        z[offset_z + col] += ia[row + col * BlockSideLength] * y[offset_y + row];
            else
                for (int_type row = 0; row < row_limit; ++row)
                    for (int_type col = 0; col < col_limit; ++col)
                        z[offset_z + col] += ia[row * BlockSideLength + col] * y[offset_y + row];
        }
        // release
        bs_a.release(a, false);
    }

    // +-+ end matrix-vector multiplication +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    // +-+-+-+ vector-vector multiplication +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    static void recursive_matrix_from_vectors(swappable_block_matrix_type A, const column_vector_type& l,
                                              const row_vector_type& r, vector_size_type offset_l = 0, vector_size_type offset_r = 0)
    {
        // catch empty intervals
        if (A.get_height() * A.get_width() == 0)
            return;
        // base case
        if (A.get_height() == 1 || A.get_width() == 1)
        {
            naive_matrix_from_vectors(A, l, r, offset_l, offset_r);
            return;
        }

        // partition matrix
        swappable_block_matrix_approximative_quarterer qa(A);
        // recursive creation
        // The order of recursive calls is optimized to enhance locality.
        recursive_matrix_from_vectors(qa.ul, l, r, offset_l,                      offset_r                    );
        recursive_matrix_from_vectors(qa.ur, l, r, offset_l,                      offset_r + qa.ul.get_width());
        recursive_matrix_from_vectors(qa.dr, l, r, offset_l + qa.ul.get_height(), offset_r + qa.ul.get_width());
        recursive_matrix_from_vectors(qa.dl, l, r, offset_l + qa.ul.get_height(), offset_r                    );
    }

    static void naive_matrix_from_vectors(swappable_block_matrix_type A, const column_vector_type& l,
                                          const row_vector_type& r, vector_size_type offset_l = 0, vector_size_type offset_r = 0)
    {
        for (size_type row = 0; row < A.get_height(); ++row)
            for (size_type col = 0; col < A.get_width(); ++col)
                matrix_from_vectors_swappable_block(A(row, col), A.is_transposed(), A.bs,
                                                    l, r, (offset_l + row) * BlockSideLength, (offset_r + col) * BlockSideLength);
    }

    static void matrix_from_vectors_swappable_block(swappable_block_identifier_type a,
                                                    const bool a_is_transposed, block_scheduler_type& bs_a,
                                                    const column_vector_type& l, const row_vector_type& r,
                                                    vector_size_type offset_l, vector_size_type offset_r)
    {
        // acquire
        internal_block_type& ia = bs_a.acquire(a, true);
        // multiply
        if (! bs_a.is_simulating())
        {
            int_type row_limit = std::min(BlockSideLength, unsigned(l.size() - offset_l)),
                col_limit = std::min(BlockSideLength, unsigned(r.size() - offset_r));
            if (a_is_transposed)
                for (int_type col = 0; col < col_limit; ++col)
                    for (int_type row = 0; row < row_limit; ++row)
                        ia[row + col * BlockSideLength] = l[row + offset_l] * r[col + offset_r];
            else
                for (int_type row = 0; row < row_limit; ++row)
                    for (int_type col = 0; col < col_limit; ++col)
                        ia[row * BlockSideLength + col] = l[row + offset_l] * r[col + offset_r];
        }
        // release
        bs_a.release(a, true);
    }

    // +-+ end vector-vector multiplication +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
};

// Adjust choose_level_for_feedable_sw, too!
template <typename ValueType, unsigned BlockSideLength>
const int_type matrix_operations<ValueType, BlockSideLength>::strassen_winograd_base_case_size = 3;

} // namespace matrix_local

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_MATRIX_ARITHMETIC_HEADER
