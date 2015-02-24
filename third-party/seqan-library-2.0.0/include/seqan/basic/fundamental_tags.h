// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Global (future: generic) tag definitions.
// ==========================================================================

// TODO(holtgrew): This should probably be minimalized, tags should be moved to single modules whenever possible.

#ifndef SEQAN_BASIC_BASIC_TAGS_H_
#define SEQAN_BASIC_BASIC_TAGS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Tag<T>
// ----------------------------------------------------------------------------

/*!
 * @class Tag
 * @headerfile <seqan/basic.h>
 * @brief Template for tag definition.
 *
 * @signature template <typename T>
 *            struct Tag;
 *
 * @tparam T Any parameterless types.
 *
 * This <tt>struct</tt> is defined such that parameter less tags are easier recognizeable.  This is best explained with
 * the example below.
 *
 * @section Examples
 *
 * Usually, tags are defined in the following way.
 *
 * @code{.cpp}
 * struct SomeTag_;
 * typedef Tag<SomeTag_> SomeTag;
 * @endcode
 *
 * They are then used as follows.
 *
 * @code{.cpp}
 * template <typename T>
 * void f(T const & x, SomeTag const & tag)
 * {
 *     // ...
 * }
 *
 * // Somewhere else:
 * f(3, SomeTag());
 * @endcode
 *
 * This has the advantages that (1) the type of tag parameters is printed as <tt>Tag&lt;SomeTag_&gt;</tt> in compiler error
 * traces.  Furthermore, (2) parameter less tags can be defined redundantly in multiple headers and we can still
 * instantiate them anywhere where they are declared.  The latter (2) cannot be achieved with only forward declaration
 * (<tt>struct SomeTag;</tt>) or full declarations (<tt>struct SomeTag {};</tt>) everywhere.
 */

template <typename T>
struct Tag {};

// ----------------------------------------------------------------------------
// Tag Default
// ----------------------------------------------------------------------------

/*!
 * @tag Default
 * @headerfile <seqan/basic.h>
 * @brief Tag that specifies default behaviour.
 *
 * @signature typedef Tag<Default_> Default;
 */

struct Default_;
typedef Tag<Default_> Default;

// ----------------------------------------------------------------------------
// Tag Nothing
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should use Tag<>.

/*!
 * @tag Nothing
 * @headerfile <seqan/basic.h>
 * @brief Tag tha trepresents an absent parameter or absent type.
 *
 * @signature struct Nothing {};
 */

struct Nothing_;
typedef Tag<Nothing_> Nothing;

// --------------------------------------------------------------------------
// Tag Raw
// --------------------------------------------------------------------------

struct Raw_;
typedef Tag<Raw_> Raw;

// ----------------------------------------------------------------------------
// Tag Move
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should probably go into basic_transport.h.

/*!
 * @tag Move
 * @headerfile <seqan/align.h>
 * @brief Switch to force move.
 *
 * @signature typedef Tag<Move_> Move;
 *
 * The difference between move constructor and copy constructor is that the source object is not copied but moved into
 * the target object.  The source object can lose its content and will be empty after this operation in this case.  A
 * move constructor can sigificantly faster than a copy constructor.
 *
 * @section Examples
 *
 * @code{.cpp}
 * String source("hello");
 * String target(source, Move()); // source is moved to target
 * std::cout << source; //nothing printed since source lost content
 * std::cout << target; //"hello"
 * @endcode
 *
 * Move constructors are like copy-constructors.  However, their argument is not const.
 *
 * @code{.cpp}
 * class Klass
 * {
 * public:
 *     seqan::String m;
 *
 *     // Copy constructor, other is left untouched.
 *     Klass(Klass const & other) { ... }
 *
 *     // Move constructor, leaves other and its members in an "empty" state.
 *     Klass(Klass & other, seqan::Move const &) { ... }
 * };
 * @endcode
 *
 * @see AssignableConcept#move
 */

struct Move_;
typedef Tag<Move_> Move;

// ----------------------------------------------------------------------------
// Tag MinimalCtor
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should go into initialization part of alphabet header set.

//construct without initializing
struct MinimalCtor_;
typedef Tag<MinimalCtor_> MinimalCtor;

// ----------------------------------------------------------------------------
// Tag MinimalCtor
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should go into initialization part of alphabet header set.

//construct with initializing
struct NonMinimalCtor_;
typedef Tag<NonMinimalCtor_> NonMinimalCtor;

// ----------------------------------------------------------------------------
// Tag MinimalCtor
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should go into the iterators header set?

//Pass to c'tor of iterator to move it to the end
struct GoEnd_;
typedef Tag<GoEnd_> GoEnd;

// ----------------------------------------------------------------------------
// Tag Serial
// ----------------------------------------------------------------------------

struct Serial_;
typedef Tag<Serial_> Serial;

// ----------------------------------------------------------------------------
// Tag TagList<TTag, TNext>
// ----------------------------------------------------------------------------

/*!
 * @class TagList
 * @headerfile <seqan/basic.h>
 * @brief A structure to represent a list of tags.
 *
 * @signature template <[typename TTag[, typename TSubList]]>
 *            struct TagList;
 *
 * @tparam TTag     The tag of the front for the list.
 * @tparam TSubList Nested list.
 */

template <typename TTag = void, typename TSubList = void>
struct TagList
{
    typedef TTag Type;
};

// ----------------------------------------------------------------------------
// Class TagSelector
// ----------------------------------------------------------------------------

/*!
 * @class TagSelector
 * @headerfile <seqan/basic.h>
 * @brief A structure to select a tag from a @link TagList @endlink.
 *
 * @signature template <typename TTagList>
 *            struct TagSelector;
 *
 * @tparam TTagList A tag list.
 */

/*!
 * @var T TagSelector::tagId
 * @headerfile <seqan/basic.h>
 * @brief Stores the index of a Tag in the tag list.
 */

template <typename TTagList = void>
struct TagSelector
{
    int tagId;

    TagSelector() :
        tagId(-1) {}    // -1 is an important initialization to signal a not yet selected tag (used for file format auto-detection)

    inline bool
    operator==(TagSelector const & other) const
    {
        return other.tagId == tagId;
    }
};

template <typename TTag, typename TSubList>
struct TagSelector< TagList<TTag, TSubList> >
        : TagSelector<TSubList>
{
    typedef TTag                    Type;
    typedef TagSelector<TSubList>   Base;
};


template <typename TTagList>
struct Value<TagSelector<TTagList> >
{
    typedef int Type;
};

template <typename TTagList>
inline typename Reference<TagSelector<TTagList> >::Type
value(TagSelector<TTagList> &selector)
{
    return selector.tagId;
}

template <typename TTagList>
inline typename Reference<TagSelector<TTagList> const>::Type
value(TagSelector<TTagList> const &selector)
{
    return selector.tagId;
}

// ----------------------------------------------------------------------------
// Tag DotDrawing
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/*!
 * @tag DotDrawing
 * @headerfile <seqan/basic.h>
 * @brief Switch to trigger drawing in dot format.
 *
 * @signature typedef Tag<DotDrawing_> DotDrawing;
 */

struct DotDrawing_;
typedef Tag<DotDrawing_> DotDrawing;

// TODO(holtgrew): Should probably not be defined here.
// TODO(holtgrew): Are these used at all?

/*!
 * @tag HammingDistance
 * @headerfile <seqan/basic.h>
 * @brief Hamming distance.
 *
 * @signature typedef Tag<HammingDistance_> HammingDistance;
 */

// TODO(holtgrew): Why ambiguous here? Edit distance is the more common name.

/*!
 * @tag LevenshteinDistance
 * @headerfile <seqan/basic.h>
 * @brief Levenshtein distance.
 *
 * @signature typedef Tag<LevenshteinDistance_> LevenshteinDistance;
 */

/*!
 * @tag EditDistance
 * @headerfile <seqan/basic.h>
 * @brief Edit distance.
 *
 * @signature typedef Tag<LevenshteinDistance_> EditDistance;
 */

struct HammingDistance_;
struct LevenshteinDistance_;

typedef Tag<HammingDistance_>       HammingDistance;
typedef Tag<LevenshteinDistance_>   LevenshteinDistance;
typedef Tag<LevenshteinDistance_>   EditDistance;


// ----------------------------------------------------------------------------
// Tag Blat
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

struct Blat_;
typedef Tag<Blat_> Blat;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction LENGTH for TagLists
// ----------------------------------------------------------------------------

/*!
 * @mfn TagList#LENGTH
 * @brief Return the length of a tag list.
 *
 * @signature LENGTH<TTagList>::VALUE;
 *
 * @tparam TTagList The TagList to query for its length.
 *
 * @return VALUE The length of the TagList.
 */

template <>
struct LENGTH<void>
{
    enum { VALUE = 0 };
};

template <typename TTag>
struct LENGTH<TagList<TTag, void> >
{
    enum { VALUE = 1 };
};

template <typename TTag, typename TSubList>
struct LENGTH<TagList<TTag, TSubList> >
{
    enum { VALUE = LENGTH<TSubList>::VALUE + 1 };
};

// ----------------------------------------------------------------------------
// Metafunction TagListValue
// ----------------------------------------------------------------------------

/*!
 * @mfn TagList#TagListValue
 * @brief A metafunction to retrieve a tag from a TagList.
 *
 * @signature TagListValue<TagList, TAG_ID>::Type;
 *
 * @tparam TagList The TagList to query.
 * @tparam TAG_ID  An index of a tag in the tag list (<tt>int</tt>.  This value must be in
 *                 <tt>0..LENGTH&lt;TTagList&gt;::VALUE -1</tt>.
 */

template <typename TList, int I>
struct TagListValue
{
    typedef void Type;
};

template <typename TTag, typename TSubList, int I>
struct TagListValue<TagList<TTag, TSubList>, I>:
    public If<Eval<I == LENGTH<TSubList>::VALUE>,
              TTag,
              typename TagListValue<TSubList, I>::Type> {};

//template <typename TTag, typename TSubList, int I>
//struct TagListValue<TagList<TTag, TSubList>, I> :
//    public typename TagListValue<TSubList, I - 1> {};
//
//template <typename TTagList, int I>
//struct TagListValue<TagSelector<TTagList>, I> :
//    public typename TagListValue<TTagList, I> {};

// ----------------------------------------------------------------------------
// Metafunction Find
// ----------------------------------------------------------------------------

/*!
 * @mfn Find
 * @headerfile <seqan/basic.h>
 * @brief A metafunction to retrieve the index of a tag in the TagList.
 *
 * @signature Find<TTagList, TSearchTag>::VALUE;
 *
 * @tparam TSearchTag A tag to retrieve the index of.
 * @tparam TTagList   A tag list.
 *
 * @return VALUE This meta-function can be used to test whether the value of a TagSelector equals a specific tag.
 *
 * @section Examples
 *
 * @code{.cpp}
 * AutoSeqFormat format;
 * if (format.tagId == Find<AutoSeqFormat, Fasta>::VALUE)
 * {
 *     // do something specific to Fasta format
 * }
 *
 * // or even shorter:
 *
 * if (isEqual(format.tagId, Fasta()))
 * {
 *     // do something specific to Fasta format
 * }
 * @endcode
 */

template <typename TList, typename TSearchTag>
struct Find;

template <typename TTag, typename TSearchTag>
struct Find<TagList<TTag, void>, TSearchTag>
{
    enum { VALUE = -1 };    // not found
};

template <typename TSearchTag>
struct Find<TagList<TSearchTag, void>, TSearchTag>
{
    enum { VALUE = 0 };
};

template <typename TTag, typename TSubList, typename TSearchTag>
struct Find<TagList<TTag, TSubList>, TSearchTag>
{
    enum { VALUE = Find<TSubList, TSearchTag>::VALUE };
};

template <typename TSubList, typename TSearchTag>
struct Find<TagList<TSearchTag, TSubList>, TSearchTag>
{
    enum { VALUE = LENGTH<TSubList>::VALUE };
};

template <typename TTagList, typename TSearchTag>
struct Find<TagSelector<TTagList>, TSearchTag>:
    public Find<TTagList, TSearchTag> {};

template <typename TTagList, typename TSearchTag>
inline int find(TagSelector<TTagList> const &, TSearchTag const &)
{
    return Find<TTagList, TSearchTag>::VALUE;
}

// ============================================================================
// Functions
// ============================================================================

// isEqual()
template <typename TTagList, typename TTag>
inline bool isEqual(TagSelector<TTagList> const &selector, TTag const &)
{
    return selector.tagId == Find<TTagList, TTag>::VALUE;
}

// assign()
template <typename TTagList, typename TTag>
inline void assign(TagSelector<TTagList> &selector, TTag &)
{
    SEQAN_ASSERT_NEQ(int(Find<TTagList, TTag>::VALUE), -1);
    selector.tagId = Find<TTagList, TTag>::VALUE;
}

template <typename TTagList, typename TTag>
inline void assign(TagSelector<TTagList> &selector, TTag const &)
{
    SEQAN_ASSERT_NEQ(int(Find<TTagList, TTag>::VALUE), -1);
    selector.tagId = Find<TTagList, TTag>::VALUE;
}

template <typename TTagList>
inline void assign(TagSelector<TTagList> &selector, TagSelector<TTagList> &other)
{
    selector.tagId = other.tagId;
}

template <typename TTagList>
inline void assign(TagSelector<TTagList> &selector, TagSelector<TTagList> const &other)
{
    selector.tagId = other.tagId;
}

// tagApply()
template <typename TFunctor, typename TTag>
inline bool
tagApply(TFunctor &func, TagList<TTag>)
{
    return func(TTag());
}

template <typename TFunctor, typename TTag, typename TSubList>
inline bool
tagApply(TFunctor &func, TagList<TTag, TSubList>)
{
    if (func(TTag()))
        return true;
    return tagApply(func, TSubList());
}


}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_TAGS_H_
