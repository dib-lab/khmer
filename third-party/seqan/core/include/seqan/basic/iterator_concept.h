// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef CORE_INCLUDE_SEQAN_BASIC_ITERATOR_CONCEPT_H_
#define CORE_INCLUDE_SEQAN_BASIC_ITERATOR_CONCEPT_H_

namespace seqan {

/**
.Metafunction.Pointer
..cat:Basic
..summary:Returns pointer to an object, required for @Function.operator->@, for example.
..signature:Pointer<T>::Type
..param.T:The type to query.
..returns:Pointer type.
..include:seqan/basic.h
 */

// Forward Declaration / Prototype.
template <typename T> struct Pointer;

/*!
 * @concept IteratorAssociatedTypesConcept
 * @headerfile <seqan/basic.h>
 * @brief Requires metafunctions for the associated types used in the iterator concepts.
 * 
 * @signature IteratorAssociatedTypesConcept<T>
 * 
 * @section Remarks
 * 
 * The SeqAn iterators mirror the definitions from <a href="http://generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator">ConceptC++</a>.
 */

/*!
 * @mfn IteratorAssociatedTypesConcept#Value
 * @brief The value type of the iterator (same as the value type of the underlying container).
 *
 * @signature Value<TIter>::Type
 *
 * @tparam TIter The <tt>TIter</tt> class to query for its value type.
 *
 * @return Type The value type of <tt>TIter</tt>
 */

/*!
 * @mfn IteratorAssociatedTypesConcept#GetValue
 * @brief The get-value type of the iterator (same as the get-value type of the underlying container).
 *
 * @signature GetValue<TIter>::Type
 *
 * @tparam TIter The <tt>TIter</tt> class to query for its get-value type.
 *
 * @return Type The get-value type of <tt>TIter</tt>
 */

/*!
 * @mfn IteratorAssociatedTypesConcept#Reference
 * @brief The reference type of the iterator (same as the reference type of the underlying container).
 *
 * @signature Reference<TIter>::Type
 *
 * @tparam TIter The <tt>TIter</tt> class to query for its reference type.
 *
 * @return Type The reference type of <tt>TIter</tt>
 */

/*!
 * @mfn IteratorAssociatedTypesConcept#Difference
 * @brief Returns the type for distances between two iterators.
 *
 * @signature Size<TContainer>::Type
 *
 * @tparam TContainer The Container type to query.
 * @return Type       The type to use for storing iterator distances sizes.
 *
 * This must be the same type as the distance type of the containers iterators.
 */

/*!
 * @mfn IteratorAssociatedTypesConcept#Pointer
 * @brief Returns pointer to an object, required for <tt>operator-></tt>, for example.
 * 
 * @signature Pointer<TIter>::Type
 * 
 * @tparam TIter The type to query.
 * 
 * @return Type Pointer type.
 */

/*!
 * @fn IteratorAssociatedTypesConcept#operator*
 * @brief Returns reference to the pointed-to value.
 * 
 * @signature TReference operator*(it);
 * 
 * @param it The iterator to dereference.
 * 
 * @return TReference The reference type.
 */

/*!
 * @fn IteratorAssociatedTypesConcept#value
 * @brief Returns reference to the pointed-to value.
 * @deprecated Use <tt>operator*()</tt> instead.
 * 
 * @signature TReference value(it);
 * 
 * @param it The iterator to dereference.
 * 
 * @return TReference The reference type.
 */

/*!
 * @fn IteratorAssociatedTypesConcept#getValue
 * @brief Returns get-value of pointed-to character.
 * @deprecated Use <tt>operator*()</tt> instead.
 * 
 * @signature TGetValue getValue(it);
 * 
 * @param it The iterator to get get-value from.
 * 
 * @return TGetValue The get-value that is pointed to.
 */

/**
.Concept.IteratorAssociatedTypesConcept
..cat:Iterators
..summary:Requires metafunctions for the associated types used in the iterator concepts.
..signature:IteratorAssociatedTypesConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h

.Metafunction.Value.concept:Concept.IteratorAssociatedTypesConcept
.Metafunction.GetValue.concept:Concept.IteratorAssociatedTypesConcept
.Metafunction.Difference.concept:Concept.IteratorAssociatedTypesConcept
.Metafunction.Reference.concept:Concept.IteratorAssociatedTypesConcept
.Metafunction.Pointer.concept:Concept.IteratorAssociatedTypesConcept
 */

SEQAN_CONCEPT(IteratorAssociatedTypesConcept, (T))
{
    typedef typename Value<T>::Type      TValue;
    typedef typename GetValue<T>::Type   TGetValue;
    typedef typename Difference<T>::Type TDifference;
    typedef typename Reference<T>::Type  TReference;
    typedef typename Pointer<T>::Type    TPointer;
    
    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TDifference>));

    SEQAN_CONCEPT_USAGE(IteratorAssociatedTypesConcept)
    {
    }
};

/*!
 * @concept InputIteratorConcept
 * @extends IteratorAssociatedTypesConcept
 * @extends CopyConstructibleConcept
 * @extends EqualityComparableConcept
 * @headerfile seqan/basic.h
 * @brief Iterator that allows dereferenced reading.
 * 
 * @signature InputIteratorConcept<T>
 * 
 * @section Remarks
 * 
 * The SeqAn iterators mirror the definitions from <a href="http://generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator">ConceptC++</a>.
 * 
 * @see BasicOutputIteratorConcept
 */

/*!
 * @fn InputIteratorConcept#operator++(suffix)
 * @brief C++ built-in suffix increment operator.
 *
 * @signature TIterator operator++(it, i)
 *
 * @param it The iterator to increment.
 * @param i  An integer, used to mark suffix decrement.
 *
 * @return TIterator A copy of the original iterator.
 */

/*!
 * @fn InputIteratorConcept#operator++(prefix)
 * @brief C++ built-in prefix increment operator.
 *
 * @signature TIterator operator++(it)
 *
 * @param it The iterator to increment.
 *
 * @return TIterator A copy of the original iterator.
 */

/*!
 * @fn InputIteratorConcept#operator->
 * @brief C++ built-in structure dereference operator.
 *
 * @signature TResult operator->(it)
 *
 * @param it The iterator to structure-dereference.
 *
 * @return TResult Either a pointer or another type.  If it is another type, the <tt>operator-></tt> is called
 *                 recursively.
 */

/*!
 * @fn InputIteratorConcept#goNext
 * @brief Iterates to next position.
 *
 * @signature void goNext(it);
 * 
 * @param it The iterator to increment.
 * 
 * @section Remarks
 * 
 * This function is equivalent to <tt>++iterator</tt>.
 */

/**
.Concept.InputIteratorConcept
..cat:Iterators
..summary:Iterator that allows dereferenced reading.
..baseconcept:Concept.IteratorAssociatedTypesConcept
..baseconcept:Concept.CopyConstructibleConcept
..baseconcept:Concept.EqualityComparableConcept
..signature:InputIteratorConcept<T>
..see:Concept.BasicOutputIteratorConcept
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h

.Function.operator->.concept:Concept.InputIteratorConcept
.Function.operator++ (prefix).concept:Concept.InputIteratorConcept
.Function.operator++ (suffix).concept:Concept.InputIteratorConcept
.Function.goNext.concept:Concept.InputIteratorConcept
.Function.operator*.concept:Concept.InputIteratorConcept
.Function.operator!=.concept.Concept.InputIteratorConcept
 */

SEQAN_CONCEPT_REFINE(InputIteratorConcept, (T), (IteratorAssociatedTypesConcept)(CopyConstructible)(EqualityComparable))
{
    typedef typename Value<T>::Type      TValue;
    typedef typename GetValue<T>::Type   TGetValue;
    typedef typename Difference<T>::Type TDifference;
    typedef typename Reference<T>::Type  TReference;
    typedef typename Pointer<T>::Type    TPointer;

    TValue v;
    T      x, y;

    SEQAN_CONCEPT_USAGE(InputIteratorConcept)
    {
        TReference & rv = v;
        T & rx =          x;

        SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TDifference>));
        // TODO(holtgrew): requires Convertible<reference, value_type>;
        // TODO(holtgrew): requires Convertible<pointer, cont value_type*>;

        // TODO(holtgrew): requires Dereferenceable<postincrement_result>;
        // TODO(holtgrew): requires requires Dereferenceable<postincrement_result>;

        // operator->: Cannot check this, need to know member for this.
        
        sameType(++x, rx);
        sameType(x++, y);
        goNext(x);

        sameType(rv, *x);

        x != x;
    }
};

/*!
 * @concept BasicOutputIteratorConcept
 * @extends IteratorAssociatedTypesConcept
 * @extends CopyConstructibleConcept
 * @headerfile <seqan/basic.h>
 * 
 * @brief Iterator that allows dereferenced writing.
 * 
 * @signature OutputIteratorConcept<T>
 * 
 * @section Remarks
 * 
 * The SeqAn iterators mirror the definitions from <a href="http://generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator">ConceptC++</a>. * 
 *
 * @section Examples
 * 
 * In the following, <tt>x</tt> is an iterator to type <tt>X</tt>, <tt>t</tt> is
 * a valid rvalue of type <tt>X</tt>.
 * 
 * The following expressions must be valid.
 * 
 * @code{.cpp}
 * *x = t     // Dereference assignment.
 * ++x        // Preincrement.
 * (void)x++  // Postincrement.
 * *x++ = t   // Postincrement and assign.
 *  
 * assignValue(x, t);
 * @endcode
 *
 * @see InputIteratorConcept
 */

/*!
 * @fn BasicOutputIteratorConcept#assignValue
 * @brief Assigns value to iterator.
 * @deprecated Use dereferencement and assignment instead.
 * 
 * @signature void assignValue(it, value);
 * @signature assignValue(container, pos, value)
 * 
 * @param it The iterator to assign value to.
 * @param value A value that is assigned to the item <tt>it</tt> points to.
 */

/*!
 * @fn BasicOutputIteratorConcept#operator++(suffix)
 * @brief C++ built-in suffix increment operator.
 * 
 * @signature TIterator operator++(it, i)
 * 
 * @param it The iterator to increment.
 * @param i An integer, used to mark suffix decrement.
 * 
 * @return TIterator A copy of the original iterator.
 */

/*!
 * @fn BasicOutputIteratorConcept#operator++(prefix)
 * @brief C++ built-in prefix increment operator.
 * 
 * @signature TIterator operator++(it)
 * 
 * @param it The iterator to increment.
 * 
 * @return TIterator A copy of the original iterator.
 */

/*!
 * @fn BasicOutputIteratorConcept#goNext
 * @brief Iterates to next position.
 * 
 * @signature void goNext(it);
 * 
 * @param it The iterator to increment.
 * 
 * @section Remarks
 * 
 * This function is equivalent to <tt>++iterator</tt>.
 */

/**
.Concept.BasicOutputIteratorConcept
..cat:Iterators
..summary:Iterator that allows dereferenced writing.
..baseconcept:Concept.IteratorAssociatedTypesConcept
..baseconcept:Concept.CopyConstructibleConcept
..signature:OutputIteratorConcept<T>
..see:Concept.InputIteratorConcept
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
..example.text:In the following, $x$ is an iterator to type $X$, $t$ is a valid rvalue of type $X$.
..example.text:The following expressions must be valid.
..example.code:
*x = t     // Dereference assignment.
++x        // Preincrement.
(void)x++  // Postincrement.
*x++ = t   // Postincrement and assign.

assignValue(x, t);

.Function.assignValue.concept:Concept.BasicOutputIteratorConcept
.Function.operator++ (prefix).concept:Concept.BasicOutputIteratorConcept
.Function.operator++ (suffix).concept:Concept.BasicOutputIteratorConcept
.Function.goNext.concept:Concept.BasicOutputIteratorConcept
 */

SEQAN_CONCEPT_REFINE(BasicOutputIteratorConcept, (T), (CopyConstructible))
{
    typedef typename Value<T>::Type TValue;

    SEQAN_CONCEPT_ASSERT((Is<Assignable<TValue> >));
    
    T      x;
    TValue v;

    SEQAN_CONCEPT_USAGE(BasicOutputIteratorConcept)
    {
        *x = v;
        assignValue(x, v);
        value(x) = v;

        ++x;
        ignoreUnusedVariableWarning(x++);
        goNext(x);
        *x++ = v;

        ignoreUnusedVariableWarning(x);
    }
};

/*!
 * @concept ForwardIteratorConcept
 * @extends InputIteratorConcept
 * @extends DefaultConstructibleConcept
 * @headerfile <seqan/basic.h>
 * @brief Iterator that allows passing over a linear sequence multiple times.
 * 
 * @signature ForwardIteratorConcept<T>
 * 
 * @section Remarks
 * 
 * The SeqAn iterators mirror the definitions from <a href="http://generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator">ConceptC++</a>.
 * 
 * @section Examples
 * 
 * In the following, <tt>x</tt> is an iterator to type <tt>X</tt>.
 * 
 * The following expressions must be valid.
 * 
 * @code{.cpp}
 * ++x  // Preincrement.
 * x++  // Postincrement.
 * @endcode
 */

/**
.Concept.ForwardIteratorConcept
..cat:Iterators
..summary:Iterator that allows passing over a linear sequence multiple times.
..baseconcept:Concept.InputIteratorConcept
..baseconcept:Concept.DefaultConstructibleConcept
..signature:ForwardIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
..example.text:In the following, $x$ is an iterator to type $X$.
..example.text:The following expressions must be valid.
..example.code:
++x  // Preincrement.
x++  // Postincrement.
 */

SEQAN_CONCEPT_REFINE(ForwardIteratorConcept, (T), (InputIteratorConcept)(DefaultConstructible))
{
    typedef typename Value<T>::Type TValue;

    T x;

    SEQAN_CONCEPT_USAGE(ForwardIteratorConcept)
    {
        ++x;
        ignoreUnusedVariableWarning(*x++);
    }
};

/*!
 * @concept MutableForwardIteratorConcept
 * @extends ForwardIteratorConcept
 * @extends BasicOutputIteratorConcept
 * @headerfile <seqan/basic.h>
 * 
 * @brief A @link ForwardIteratorConcept Forward Iterator @endlink that allows dereferenced assignment.
 * 
 * @section Remarks
 * 
 * The SeqAn iterators mirror the definitions from <a href="http://generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator">ConceptC++</a>.
 */

/**
.Concept.MutableForwardIteratorConcept
..cat:Iterators
..summary:A @Concept.ForwardIteratorConcept|Forward Iterator@ that allows dereferenced assignment.
..baseconcept:Concept.ForwardIteratorConcept
..baseconcept:Concept.BasicOutputIteratorConcept
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableForwardIteratorConcept, (T), (ForwardIteratorConcept)(BasicOutputIteratorConcept))
{
    typedef typename Value<T>::Type      TValue;

    T      x;
    TValue v;

    SEQAN_CONCEPT_USAGE(MutableForwardIteratorConcept)
    {
        *x = v;
    }
};

/*!
 * @concept BidirectionalIteratorConcept
 * 
 * @headerfile seqan/basic.h
 * 
 * @extends ForwardIteratorConcept
 * 
 * @brief Iterator that can be both incremented and decremented.
 * 
 * @signature BidirectionalIteratorConcept<T>
 * 
 * @section Remarks
 * 
 * The SeqAn iterators mirror the definitions from <a href="http://generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator">ConceptC++</a>.
 * 
 * @section Examples
 * 
 * In the following, <tt>x</tt> is an iterator to type <tt>X</tt>.
 * 
 * The following expressions must be valid.
 * 
 * @code{.cpp}
 * --x  // Predecrement.
 * x--  // Postdecrement.
 * @endcode
 */

/*!
 * @fn BidirectionalIteratorConcept#operator--(prefix)
 * @brief C++ built-in prefix decrement operator.
 *
 * @signature TIterator operator--(it);
 *
 * @param it The iterator to increment.
 *
 * @return TIterator Reference to the incremented iterator.
 */

/*!
 * @fn BidirectionalIteratorConcept#operator--(suffix)
 * @brief C++ built-in suffix decrement operator.
 *
 * @signature TIterator operator--(it, i);
 *
 * @param it The iterator to increment.
 * @param i  An int value to mark the operator as suffix decrement.
 *
 * @return TIterator Reference to the incremented iterator.
 */

/*!
 * @fn BidirectionalIteratorConcept#goPrevious
 * @brief Iterates to pevious position.
 * 
 * @signature void goPrevious(it);
 * 
 * @param it Iterator to move to previous position.
 * 
 * @section Remarks
 * 
 * This function is equivalent to <tt>--iterator</tt>.
 */

/**
.Concept.BidirectionalIteratorConcept
..cat:Iterators
..summary:Iterator that can be both incremented and decremented.
..baseconcept:Concept.ForwardIteratorConcept
..signature:BidirectionalIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
..example.text:In the following, $x$ is an iterator to type $X$.
..example.text:The following expressions must be valid.
..example.code:
--x  // Predecrement.
x--  // Postdecrement.

.Function.operator-- (prefix).concept:Concept.BidirectionalIteratorConcept
.Function.operator-- (suffix).concept:Concept.BidirectionalIteratorConcept
.Function.goPrevious.concept:Concept.BidirectionalIteratorConcept
 */

SEQAN_CONCEPT_REFINE(BidirectionalIteratorConcept, (T), (ForwardIteratorConcept))
{
    T x;

    SEQAN_CONCEPT_USAGE(BidirectionalIteratorConcept)
    {
        --x;
        x--;
        goPrevious(x);
    }
};

/*!
 * @concept MutableBidirectionalIteratorConcept
 * @extends BidirectionalIteratorConcept
 * @brief Bidirectional iterator that also allows writing of dereferenced values.
 * @headerfile <seqan/basic.h>
 *
 * @signature MutableBidirectionalIteratorConcept<T>
 */

/**
.Concept.MutableBidirectionalIteratorConcept
..cat:Iterators
..summary:A @Concept.BidirectionalIteratorConcept|Bidirectional Iterator@ that allows dereferenced assignment
..baseconcept:Concept.ForwardIteratorConcept
..signature:MutableBidirectionalIteratorConcept<T>
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableBidirectionalIteratorConcept, (T), (BidirectionalIteratorConcept)(MutableForwardIteratorConcept))
{
    typedef typename Value<T>::Type TValue;

    T      x;
    TValue v;

    SEQAN_CONCEPT_USAGE(MutableBidirectionalIteratorConcept)
    {
        *x = v;
    }
};

/*!
 * @concept RandomAccessIteratorConcept
 * @extends BidirectionalIteratorConcept
 * @extends LessThanComparableConcept
 * @brief An iterator allowing random access.
 * @headerfile <seqan/basic.h>
 *
 * @signature RandomAccessIteratorConcept<T>
 *
 * @section Remarks
 * 
 * This function is equivalent to <tt>--iterator</tt>.
 *
 * @section Examples
 * 
 * In the following, <tt>x</tt> is an iterator to type <tt>X</tt>, <tt>t</tt> is
 * a valid rvalue of type <tt>X</tt>, <tt>n</tt> is a distance type.
 * 
 * The following expressions must be valid.
 * 
 * @code{.cpp}
 * x += n    // Iterator addition assignment.
 * x + n     // Iterator addition.
 * n + i     // Iterator addition.
 * x -= n    // Iterator subtraction assignment.
 * x - n     // Iterator subtraction.
 * x - a     // Difference.
 * x[n]      // Element operator.
 * @endcode
 */

/*!
 * @mfn RandomAccessIteratorConcept#Difference
 * @brief Type of an object that stores the difference between two iterators.
 * 
 * @signature Difference<T>::Type
 * 
 * @tparam T Type for which the difference type is determined.
 * 
 * @return Type The difference type.
 * 
 * @section Remarks
 * 
 * In most cases this type is <tt>ptrdiff_t</tt>.
 */

/*!
 * @fn RandomAccessIteratorConcept#difference
 * @brief The difference between two iterators.
 * 
 * @signature TDifference difference(begin, end);
 * 
 * @param begin Iterator to the first position of a range.
 * @param end  Iterator behind the last position of a range.
 * 
 * @return TDifference Length of the range between <tt>begin</tt> and <tt>end</tt>, type from 
 * 
 * @section Remarks
 * 
 * This function is equivalent to <tt>end - begin</tt>.
 * 
 * Usually, <tt>begin</tt> and <tt>end</tt> have the same type.
 */

/*!
 * @fn RandomAccessIteratorConcept#operator+=
 * @brief C++ built-in addition assignment operator.
 *
 * @signature TIterator operator+=(it, diff);
 */

/*!
 * @fn RandomAccessIteratorConcept#operator+
 * @brief C++ built-in addition operator.
 *
 * @signature TIterator operator+(lhs, rhs);
 */

/*!
 * @fn RandomAccessIteratorConcept#operator-=
 * @brief C++ built-in subtraction assignment operator.
 *
 * @signature TIterator operator-=(it, diff);
 */

/*!
 * @fn RandomAccessIteratorConcept#goNil
 * @brief Moves iterator to nil position.
 * 
 * @signature void goNil(it);
 * 
 * @param it The iterator that will be moved.
 *
 * 
 * @section Remarks
 * 
 * <tt>it</tt> is set to an invalid position, e.g. <tt>NULL</tt> for pointer types.
 */

/*!
 * @fn RandomAccessIteratorConcept#goFurther
 * @brief Iterates some steps further.
 * 
 * @signature void goFurther(iterator, steps);
 * 
 * @param it    The iterator to move.
 * @param steps Number of steps <tt>iterator</tt> should be moved further.
 * 
 * @section Remarks
 * 
 * This function is equivalent to <tt>iterator += steps</tt>.
 */

/*!
 * @fn RandomAccessIteratorConcept::operator[]
 * @brief C++ built-in array subscript operator.
 *
 * @signature TReference TIterator::operator[](pos);
 *
 * @tparam TReference The referenced element.
 * @tparam TIterator  The iterator to subscript.
 *
 * @param pos Position to get element at.
 */

/*!
 * @fn RandomAccessIteratorConcept#atNil
 * @brief Tests whether iterator is at nil position.
 * 
 * @signature bool atNil(it);
 * 
 * @param it The iterator to query.
 *
 * @return bool Whether or not the iterator is at the nil positions (NULL for pointers).
 */

/**
.Concept.RandomAccessIteratorConcept
..cat:Iterators
..summary:An iterator allowing random access.
..baseconcept:Concept.BidirectionalIteratorConcept
..baseconcept:Concept.LessThanComparableConcept
..signature:RandomAccessIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
..example.text:In the following, $x$ is an iterator to type $X$, $t$ is a valid rvalue of type $X$, $n$ is a distance type.
..example.text:The following expressions must be valid.
..example.code:
x += n    // Iterator addition assignment.
x + n     // Iterator addition.
n + i     // Iterator addition.
x -= n    // Iterator subtraction assignment.
x - n     // Iterator subtraction.
x - a     // Difference.
x[n]      // Element operator.

.Metafunction.Difference.concept:Concept.RandomAccessIteratorConcept
.Function.operator+=.concept:Concept.RandomAccessIteratorConcept
.Function.operator+.concept:Concept.RandomAccessIteratorConcept
.Function.operator-=.concept:Concept.RandomAccessIteratorConcept
.Function.operator-.concept:Concept.RandomAccessIteratorConcept
.Function.difference.concept:Concept.RandomAccessIteratorConcept
.Function.operator[].concept:Concept.RandomAccessIteratorConcept
.Function.goFurther.concept:Concept.RandomAccessIteratorConcept

.Function.operator>=.concept:Concept.RandomAccessIteratorConcept
.Function.operator>.concept:Concept.RandomAccessIteratorConcept
.Function.operator<=.concept:Concept.RandomAccessIteratorConcept
 */

SEQAN_CONCEPT_REFINE(RandomAccessIteratorConcept, (T), (BidirectionalIteratorConcept)(LessThanComparable))
{
    typedef typename Difference<T>::Type TDifference;
    typedef typename Value<T>::Type      TValue;

    T           x, y, a;
    TValue      t;
    TDifference n;

    SEQAN_CONCEPT_USAGE(RandomAccessIteratorConcept)
    {

        x += n;
        goFurther(x, n);        
        ignoreUnusedVariableWarning(x + n);
        ignoreUnusedVariableWarning(n + x);
        x -= n;
        y = x - n;
        n = x - a;
        n = difference(x, a);
        ignoreUnusedVariableWarning(x[n]);

        ignoreUnusedVariableWarning(x);
        ignoreUnusedVariableWarning(y);
    }
};

/*!
 * @concept MutableRandomAccessIteratorConcept
 * @extends RandomAccessIteratorConcept
 * @headerfile <seqan/basic.h>
 * @brief A random access iterator whose dereferenced values can be assigned.
 *
 * @signature MutableRandomAccessIteratorConcept<T>
 */

/**
.Concept.MutableRandomAccessIteratorConcept
..cat:Iterators
..summary:A @Concept.RandomAccessIteratorConcept@ that allows assignable derefentiation.
..baseconcept:Concept.BidirectionalIteratorConcept
..baseconcept:Concept.LessThanComparableConcept
..signature:RandomAccessIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..example.text:The following expressions should be valid.
..example.code:
value(x, n) = t
x[n] = t
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableRandomAccessIteratorConcept, (T), (RandomAccessIteratorConcept)(MutableBidirectionalIteratorConcept))
{
    typedef typename Difference<T>::Type TDifference;
    typedef typename Value<T>::Type      TValue;

    T           x;
    TValue      t;
    TDifference n;

    SEQAN_CONCEPT_USAGE(MutableRandomAccessIteratorConcept)
    {
        value(x, n) = t;  // TODO(holtgrew): Not supported?
        x[n] = t;
    }
};

/*!
 * @concept RootedIteratorConcept
 * @extends ForwardIteratorConcept
 * @brief Iterator that knows its container.
 *
 * @signature RootedIteratorConcept<T>
 */

/*!
 * @mfn RootedIteratorConcept#Container
 * @brief Metafunction that returns the container of an iterator.
 *
 * @signature Container<TIterator>::Type
 *
 * @tparam TIterator The type of the iterator to query for its container.
 *
 * @return Type The type of the container for <tt>TIterator</tt>
 */

/*!
 * @fn RootedIteratorConcept#container
 * @brief Returns the container.
 *
 * @signature TContainer container(it);
 *
 * @param it The iterator to get the container of.
 *
 * @return TContainer The container of the iterat.r
 *
 * TODO(holtgrew): Need to document Reference_ and Parameter_.
 */

/*!
 * @fn RootedIteratorConcept#atBegin
 * @brief Queries whether the rooted iterator is at the beginning of the container or not.
 *
 * @signature bool atBegin(it);
 *
 * @param it The rooted iterator to query.
 *
 * @return bool Whether or not the iterator is at the beginning.
 */

/*!
 * @fn RootedIteratorConcept#atEnd
 * @brief Queries whether the rooted iterator is at the end of the container or not.
 *
 * @signature bool atEnd(it);
 *
 * @param it The rooted iterator to query.
 *
 * @return bool Whether or not the iterator is at the end.
 */

/**
.Concept.RootedIteratorConcept
..cat:Iterators
..summary:Iterator that knows its container.
..baseconcept:Concept.ForwardIteratorConcept
..signature:RootedIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h

.Metafunction.Container.concept:Concept.RootedIteratorConcept
.Function.container.concept:Concept.RootedIteratorConcept
.Function.atBegin.concept:Concept.RootedIteratorConcept
.Function.atEnd.concept:Concept.RootedIteratorConcept
 */

SEQAN_CONCEPT_REFINE(RootedIteratorConcept, (T), (IteratorAssociatedTypesConcept))
{
    typedef typename Container<T>::Type TContainer;

    T x;

    SEQAN_CONCEPT_USAGE(RootedIteratorConcept)
    {
        T xs;
        ignoreUnusedVariableWarning(xs);

        TContainer & c = container(x);
        atBegin(x);
        atEnd(x);
        ignoreUnusedVariableWarning(c);
    }
};


/*!
 * @concept MutableRootedIteratorConcept
 * @extends RootedIteratorConcept
 * @extends MutableForwardIteratorConcept
 * @brief Rooted iterator that allows mutation after dereferencing.
 *
 * @signature MutableRootedIteratorConcept<T>
 */

/**
.Concept.MutableRootedIteratorConcept
..cat:Iterators
..summary:A @Concept.RootedIteratorConcept|Rooted Iterator@ that allows dereferenced assignment.
..baseconcept:Concept.ForwardIteratorConcept
..signature:RootedIteratorConcept<T>
..remarks:The SeqAn iterators mirror the definitions from @http://www.generic-programming.org/languages/conceptcpp/concept_web.php?header=iterator|ConceptC++@.
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableRootedIteratorConcept, (T), (RootedIteratorConcept)(MutableForwardIteratorConcept))
{
    SEQAN_CONCEPT_USAGE(MutableRootedIteratorConcept)
    {
    }
};

/*!
 * @concept RootedRandomAccessIteratorConcept
 * @extends RootedIteratorConcept
 * @extends RandomAccessIteratorConcept
 * @brief Rooted iterator with random access.
 *
 * @signature RootedRandomAccessIteratorConcept<T>
 */

/*!
 * @mfn RootedRandomAccessIteratorConcept#Position
 * @brief Metafunction to get Position type of a rooted random access iterator.
 *
 * @signature Position<TIter>::Type
 *
 * @tparam TIter Iterator to query for its position type.
 *
 * @return Type The position type of the iterator.
 */

/*!
 * @fn RootedRandomAccessIteratorConcept#position
 * @brief Function to get the position of a rooted random access iterator.
 *
 * @signature TPosition position(it);
 *
 * @param it The iterator to query for its position.
 *
 * @return TPosition The position of <tt>it</tt>
 */

/*!
 * @fn RootedRandomAccessIteratorConcept#setPosition
 * @brief Set position of a rooted random access iterator.
 *
 * @signature void setPosition(it, pos);
 *
 * @param it  The iterator to set the position of.
 * @param pos The position to set <tt>it</tt> to.
 */

/*!
 * @fn RootedRandomAccessIteratorConcept#goBegin
 * @brief Set position of rooted random access iterator to the beginning of the container.
 *
 * @signature void goBegin(it);
 *
 * @param it  The iterator to set the position of.
 */

/*!
 * @fn RootedRandomAccessIteratorConcept#goEnd
 * @brief Set position of rooted random access iterator to the end of the container.
 *
 * @signature void goEnd(it, pos);
 *
 * @param it  The iterator to set the position of.
 */

/**
.Concept.RootedRandomAccessIteratorConcept
..cat:Iterators
..summary:An iterator that is both rooted and random access, allowing to implement @Function.position@.
..signature:RootedRandomAccessIteratorConcept<T>
..baseconcept:Concept.RandomAccessIteratorConcept
..baseconcept:Concept.RootedIteratorConcept

.Metafunction.Position.concept:Concept.RootedRandomAccessIteratorConcept
.Function.position.concept:Concept.RootedRandomAccessIteratorConcept
.Function.setPosition.concept:Concept.RootedRandomAccessIteratorConcept
.Function.goBegin.concept:Concept.RootedRandomAccessIteratorConcept
.Function.goEnd.concept:Concept.RootedRandomAccessIteratorConcept
 */

SEQAN_CONCEPT_REFINE(RootedRandomAccessIteratorConcept, (T), (RootedIteratorConcept)(RandomAccessIteratorConcept))
{
    typedef typename Position<T>::Type TPosition;

    SEQAN_CONCEPT_USAGE(RootedRandomAccessIteratorConcept)
    {
        T x;

        TPosition p = position(x);
        setPosition(x, p);
        goBegin(x);
        goEnd(x);
    }
};

/*!
 * @concept MutableRootedRandomAccessIteratorConcept
 * @extends RootedRandomAccessIteratorConcept
 * @extends MutableBidirectionalIteratorConcept
 * @brief Rooted iterator with random access that allows the mutation of dereferenced value.
 *
 * @signature RootedRandomAccessIteratorConcept<T>
 */

/**
.Concept.MutableRootedRandomAccessIteratorConcept
..cat:Iterators
..baseconcept:Concept.RootedRandomAccessIteratorConcept
..baseconcept:Concept.MutableBidirectionalIteratorConcept
..summary:A @Concept.RootedIteratorConcept|Rooted Iterator@ that allows dereferenced assignment.
..signature:MutableRootedRandomAccessIteratorConcept<T>
..include:seqan/basic.h
 */

SEQAN_CONCEPT_REFINE(MutableRootedRandomAccessIteratorConcept, (T), (RootedRandomAccessIteratorConcept)(MutableBidirectionalIteratorConcept))
{
    SEQAN_CONCEPT_USAGE(MutableRootedRandomAccessIteratorConcept)
    {
    }
};

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BASIC_ITERATOR_CONCEPT_H_
