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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_PIPE_EDIT_ENVIRONMENT_H
#define SEQAN_HEADER_PIPE_EDIT_ENVIRONMENT_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    template < typename TDistanceSpec, unsigned STEP_SIZE = 1 >
    struct EditEnvironment;

/*!
 * @class EditEnvironment
 * @extends Pipe
 * @headerfile <seqan/pipe.h>
 * @brief Outputs tuples of the <tt>tupleLen</tt> consecutive elements of the input stream.
 *
 * @signature template <typename TInput, unsigned TUPLE_LEN, bool OMIT_LAST>
 *            class Pipe<TInput, Tupler<TUPLE_LEN, OMIT_LAST> >;
 *
 * @tparam OMIT_LAST Omit half filled tuples.  If <tt>true</tt>, the output stream is <tt>tupleLen-1</tt> elements
 *                   shorter than the input stream.  If <tt>false</tt>, the lengths are identical and the last
 *                   tuples are filled with blanks (default constructed elements) for undefined entries.
 * @tparam TInput    The type of the pipeline module this module reads from.
 * @tparam TUPLE_LEN The tuple length.The tuples contain elements <tt>in[i]in[i+1]...in[i+(tupleLen-1)]</tt>.
 *
 * The output type is a @link Tuple @endlink of input elements and length <tt>tupleLen</tt> (i.e.
 * <tt>Tuple&lt;Value&lt;TInput&gt;::Type, tupleLen&gt;</tt>).
 *
 * The tuples are sequences of the form <tt>in[i]in[i-1]in[i-2]..in[i-tupleLen+1]</tt>. For <tt>omitLast=false</tt>
 * <tt>i</tt> begins with 0 and for <tt>omitLast=true</tt> <tt>i</tt> begins with <tt>tupleLen-1</tt>.
 */

    //////////////////////////////////////////////////////////////////////////////
    // pipe to enumerate the hamming 1-environment
    template < typename TInput, unsigned STEP_SIZE >
    struct Pipe< TInput, EditEnvironment< Tag<HammingDistance_>, STEP_SIZE > >
    {
        typedef typename Value< typename Value<TInput>::Type, 2 >::Type    TTuple;
        typedef typename Value<TTuple>::Type                            TValue;

        TInput                      &in;
        typename Value<Pipe>::Type    tmp, orig;
        unsigned                    errorPos;        // position of substitution
        unsigned                    character;        // replacement character
        unsigned                    skipChar;        // skip the original character

        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
            do {
                if (++character < ValueSize<TValue>::VALUE)
                    // next replacement value
                    assignValue(tmp.i2, errorPos, (TValue) character);
                else {
                    // next error position
                    assignValue(tmp.i2, errorPos, orig.i2[errorPos]);
                    character = 0;
                    if (++errorPos < length(tmp.i2)) {
                        skipChar = (unsigned) orig.i2[errorPos];
                        assignValue(tmp.i2, errorPos, (TValue) 0);
                    } else {
                        // next tuple
                        errorPos = 0;
                        ++in;
                        for(unsigned i = 1; i < STEP_SIZE && !eof(in); ++i)
                            ++in;
                        if (!eof(in)) {
                            tmp = orig = *in;
                            assignValue(tmp.i2, 0, (TValue) 0);
                        }
                    }
                }
            // output the original tuple only once
            } while ((errorPos > 0) && (character == skipChar));

            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // pipe to enumerate the levenshtein 1-environment
    template < typename TInput, unsigned STEP_SIZE >
    struct Pipe< TInput, EditEnvironment< Tag<LevenshteinDistance_>, STEP_SIZE > >
    {
        typedef typename Value< typename Value<TInput>::Type, 2 >::Type    TTuple;
        typedef typename Value<TTuple>::Type                            TValue;

        enum TState { SUBST_, DELETE_, INSERT_, INSERT_LAST_, Eof_, INSERT_EOS_ };

        TInput                      &in;
        typename Value<Pipe>::Type    tmp, orig, prev;
        unsigned                    errorPos;        // position of substitution
        unsigned                    character;        // replacement character
        unsigned                    skipChar;        // skip the original character
        TState                        state;

        Pipe(TInput& _in):
            in(_in),
            state(Eof_) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
            switch (state) {
            case SUBST_:
                // before SUBST_ (tmp[1..] == orig[1..] and tmp[0] == 0) holds
                do {
                    if (++character < ValueSize<TValue>::VALUE) {
                        // next replacement value
                        assignValue(tmp.i2, errorPos, (TValue) character);
                    } else {
                        // next substitution position
                        assignValue(tmp.i2, errorPos, orig.i2[errorPos]);
                        character = 0;
                        if (++errorPos < length(tmp.i2)) {
                            skipChar = (unsigned) orig.i2[errorPos];
                            assignValue(tmp.i2, errorPos, (TValue) 0);
                        } else {
                            // NEXT TUPLE
                            // now (tmp == orig) holds
                            ++in;
                            for(unsigned i = 1; i < STEP_SIZE && !eof(in); ++i)
                                ++in;
                            if (!eos(in)) {
                                prev = orig;
                                orig = *in;
                                tmp.i2 = orig.i2;
                                assignValue(tmp.i2, 0, prev.i2[0]);
                                assignValue(tmp.i2, 1, prev.i2[1]);
                                if (length(tmp.i2) >= 4) {
                                    errorPos = 2;
                                    state = DELETE_;
                                    //std::cerr << std::endl << "_DELETIONS____" << std::endl;
                                    return *this;
                                }
                            } else {
                                // LAST TUPLE
                                shiftLeft(orig.i2);
                                assignValue(tmp.i2, 0, orig.i2[0]);
                                assignValue(tmp.i2, 1, (TValue) 0);
                                character = 0;
                                errorPos = 1;
                                state = INSERT_LAST_;
                                //std::cerr << std::endl << "_INSERTS______" << std::endl;
                                return *this;
                            }
                        }
                    }
                // output the original tuple only once
                } while ((errorPos > 0) && (character == skipChar));
                break;
            case DELETE_:
                // before DELETE_ (prev=orig, ++in; tmp=orig=*in) holds
                assignValue(tmp.i2, errorPos, prev.i2[errorPos]);
                if (++errorPos >= length(tmp.i2) - 1) {
                    assignValue(tmp.i2, length(tmp.i2)-1, prev.i2[length(tmp.i2)-1]);
                    assignValue(tmp.i2, 0, orig.i2[0]);
                    assignValue(tmp.i2, 1, (TValue) 0);
                    character = 0;
                    errorPos = 1;
                    state = INSERT_;
                    //std::cerr << std::endl << "_INSERTS______" << std::endl;
                }
                break;

            case INSERT_EOS_:
                state = INSERT_;
            case INSERT_:
            case INSERT_LAST_:
                // before INSERT_ (prev=orig, ++in; tmp=prev) holds
                if (++character < ValueSize<TValue>::VALUE)
                    // next replacement value
                    assignValue(tmp.i2, errorPos, (TValue) character);
                else {
                    // next insert position
                    assignValue(tmp.i2, errorPos, orig.i2[errorPos]);
                    character = 0;
                    if (++errorPos >= length(tmp.i2) - 1 && state == INSERT_) {
                        tmp = orig;
                        state = SUBST_;
                        //std::cerr << std::endl << "_REPLACEMENTS_" << std::endl;
                        errorPos = 0;
                        assignValue(tmp.i2, 0, (TValue) 0);
                        break;
                    }
                    if (errorPos >= length(tmp.i2)) {
                        if (eof(in))
                            state = Eof_;
                        else {
                            tmp = orig = *in;

                            // begin to insert the first char at position 0
                            shiftRight(tmp.i2);
                            assignValue(tmp.i2, 1, (TValue) 0);
                            errorPos = 0;
                            state = INSERT_EOS_;
                            //std::cerr << std::endl << "_INSERTS______" << std::endl;
                        }
                        break;
                    }
                    assignValue(tmp.i2, errorPos, (TValue) 0);
                }
            default:;
            }
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned STEP_SIZE >
    inline bool
    control(
        Pipe< TInput, EditEnvironment< Tag<HammingDistance_>, STEP_SIZE > > &me,
        ControlBeginRead const &command)
    {
        if (!control(me.in, command)) return false;

        me.tmp = me.orig = *me.in;
        me.errorPos = 0;
        me.character = 0;

        return true;
    }

    template < typename TInput, unsigned STEP_SIZE >
    inline bool
    control(
        Pipe< TInput, EditEnvironment< Tag<LevenshteinDistance_>, STEP_SIZE > > &me,
        ControlBeginRead const &command)
    {
        if (!control(me.in, command)) return false;

        if (eof(me.in)) {
            me.state = me.Eof_;
            return true;
        }

        typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
        typedef typename Value<TTuple>::Type                            TValue;

        me.tmp = me.orig = *me.in;

        // begin to insert the first char at position 0
        shiftRight(me.tmp.i2);
        assignValue(me.tmp.i2, 0, (TValue) 0);
        me.character = 0;
        me.errorPos = 0;
        me.state = me.INSERT_;
        //std::cerr << std::endl << "_INSERTS______" << std::endl;

        return true;
    }

    template < typename TInput, unsigned STEP_SIZE >
    inline bool
    control(
        Pipe< TInput, EditEnvironment< Tag<LevenshteinDistance_>, STEP_SIZE > > &me,
        ControlEof const &)
    {
        return me.state == me.Eof_;
    }

    template < typename TInput, unsigned STEP_SIZE >
    inline bool
    control(
        Pipe< TInput, EditEnvironment< Tag<LevenshteinDistance_>, STEP_SIZE > > &me,
        ControlEos const &)
    {
        return me.state == me.Eof_ || me.state == me.INSERT_EOS;
    }

    template < typename TInput, unsigned STEP_SIZE >
    inline typename Size< Pipe< TInput, Pipe< TInput, EditEnvironment< Tag<HammingDistance_>, STEP_SIZE > > > >::Type
    length(Pipe< TInput, EditEnvironment< Tag<HammingDistance_>, STEP_SIZE > > const &me) {
        typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
        typedef typename Size< Pipe< TInput, Pipe< TInput, EditEnvironment< Tag<HammingDistance_>, STEP_SIZE > > > >::Type TSize;

        unsigned alphabetSize = ValueSize< typename Value<TTuple>::Type >::VALUE;
        /*unsigned seqs = countSequences(me.in);*/
        TSize sum = 0;
/*        for(unsigned i = 0; i < seqs; ++i)
            sum += (length((*me.in.in.in.in.set)[i]) / STEP_SIZE) * (1 + length(me.tmp.i2) * (alphabetSize - 1));
        return sum;
*/        return (length(me.in) / STEP_SIZE) * (1 + length(me.tmp.i2) * (alphabetSize - 1));
    }

    template < typename TInput, unsigned STEP_SIZE >
    inline typename Size< Pipe< TInput, Pipe< TInput, EditEnvironment< Tag<LevenshteinDistance_>, STEP_SIZE > > > >::Type
    length(Pipe< TInput, EditEnvironment< Tag<LevenshteinDistance_>, STEP_SIZE > > const &me) {
        typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
        unsigned alphabetSize = ValueSize< typename Value<TTuple>::Type >::VALUE;
        unsigned seqs = countSequences(me.in);
        // TODO: We run into problems when one sequence contains 1 or less tuples
        // length should be ommitted in future, but Pools or the skew algorithm needs to know the stream length
        if (length(me.in) >= seqs)
            return
                  (length(me.in) / STEP_SIZE)     * (1 + length(me.tmp.i2) * (alphabetSize - 1)) +            // substitutions and original
                 ((length(me.in) / STEP_SIZE) - seqs) * (length(me.tmp.i2) - 3) +                            // deletions
                (((length(me.in) / STEP_SIZE) + seqs) * (length(me.tmp.i2) - 2) + 2 * seqs) * alphabetSize;    // insertions
        else
            return 0;
    }
//}

}

#endif
