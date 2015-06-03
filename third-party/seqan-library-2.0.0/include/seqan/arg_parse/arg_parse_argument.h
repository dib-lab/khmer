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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_
#define SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_

#include <seqan/arg_parse/arg_parse_exceptions.h>
#include <seqan/arg_parse/arg_parse_type_support.h>

#include <string>
#include <vector>

#include <sstream>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// Required for _checkStringRestrictions().
inline std::string getFileExtension(ArgParseArgument const & me, unsigned pos);

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// ----------------------------------------------------------------------------
// Class ArgParseArgument
// ----------------------------------------------------------------------------

/*!
 * @class ArgParseArgument
 * @implements AssignableConcept
 * @headerfile <seqan/arg_parse.h>
 * @brief Information for a specific command line argument.
 *
 * @signature class ArgParseArgument
 */

/*!
 * @enum ArgParseArgument::ArgumentType
 * @headerfile <seqan/arg_parse.h>
 * @brief Define the type of an @link ArgParseArgument @endlink.
 *
 * @signature enum ArgParseArgument::ArgumentType;
 *
 * @section Examples
 *
 * In the following example, the types <tt>INPUT_FILE</tt>, <tt>OUTPUT_FILE</tt>, and <tt>DOUBLE</tt> are used.
 *
 * @include demos/arg_parse/argument_parser.cpp
 */

/*!
 * @val ArgParseArgument::ArgumentType STRING
 * @brief Argument is a string.
 *
 * @val ArgParseArgument::ArgumentType ArgParseArgument::INTEGER;
 * @brief Argument is a signed 32 bit integer.
 *
 * @val ArgParseArgument::ArgumentType ArgParseArgument::INT64;
 * @brief Argument is a signed 64 bit integer.
 *
 * @val ArgParseArgument::ArgumentType ArgParseArgument::DOUBLE;
 * @brief Argument is a floating point number stored as double.
 *
 * @val ArgParseArgument::ArgumentType ArgParseArgument::INPUT_FILE;
 * @brief Argument is an input file.
 *
 * @val ArgParseArgument::ArgumentType ArgParseArgument::OUTPUT_FILE;
 * @brief Argument is an output file.
 */

/*!
 * @fn ArgParseArgument::ArgParseArgument
 * @brief Constructor
 *
 * @signature ArgParseArgument::ArgParseArgument(argumentType[, argumentLabel[, isListArgument[, numberOfArgument]]]);
 *
 * @param[in] argumentType      Type of the argument (<tt>ArgParseArgument::ArgumentType</tt>).
 * @param[in] argumentLabel     Label for the argument (<tt>char const *</tt>).
 * @param[in] isListArgument    Whether or not this argument can be given multiple times (<tt>bool</tt>).
 * @param[in] numberOfArguments Number of times the argument must be given.  E.g. set to 2 for the parser to always
 *                              expect two values (<tt>int</tt>, default is 1).
 */

class ArgParseArgument
{
public:
    enum ArgumentType
    {
        // argument is
        STRING,      // .. a string
        INTEGER,     // .. an integer
        INT64,       // .. a 64 bit integer
        DOUBLE,      // .. a float
        INPUT_FILE,   // .. an inputfile (implicitly also a string)
        OUTPUT_FILE,  // .. an outputfile (implicitly also a string)
        INPUTPREFIX, // .. an inputprefix (implicitly also a string)
        OUTPUT_PREFIX // .. an outoutprefix (implicitly also a string)
    };


    // ----------------------------------------------------------------------------
    // Members to store type information
    // ----------------------------------------------------------------------------
    ArgumentType _argumentType;
    unsigned     _numberOfValues;
    std::string  _argumentLabel;
    bool         _isListArgument;

    // ----------------------------------------------------------------------------
    // Members to store the values
    // ----------------------------------------------------------------------------
    std::vector<std::string>  defaultValue;
    std::vector<std::string>  value;
    // Override for the file extension;  only used for input/output file arguments.
    std::vector<std::string> _fileExtensions;

    // ----------------------------------------------------------------------------
    // Members for restrictions
    // ----------------------------------------------------------------------------
    std::string           minValue;
    std::string           maxValue;
    std::vector<std::string> validValues;

    // ----------------------------------------------------------------------------
    // Tags
    // ----------------------------------------------------------------------------
    // Tags can be used to attach hints to the arguments (and options).  Currently,
    // this is used for tagging the "-file-ext" arguments as "file-ext-override"
    // and "gkn-ignore" for ignoring in GKN.
    std::vector<std::string> tags;

    // ----------------------------------------------------------------------------
    // Members to help text
    // ----------------------------------------------------------------------------
    std::string         _helpText;    // The help text shown on the command line

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    ArgParseArgument(ArgumentType argumentType,
                     std::string const & argumentLabel = "",
                     bool isListArgument = false,
                     unsigned numberOfValues = 1) :
        _argumentType(argumentType),
        _numberOfValues(numberOfValues),
        _argumentLabel(argumentLabel),
        _isListArgument(isListArgument),
        minValue(""),
        maxValue(""),
        _helpText("")
    {}
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Helper Function _typeToString()
// ----------------------------------------------------------------------------

inline std::string _typeToString(ArgParseArgument const & me)
{
    std::string typeName = "";

    switch (me._argumentType)
    {
    case ArgParseArgument::DOUBLE:
        typeName = "double";
        break;

    case ArgParseArgument::INTEGER:
        typeName = "integer";
        break;

    case ArgParseArgument::INT64:
        typeName = "int64";
        break;

    case ArgParseArgument::STRING:
        typeName = "string";
        break;

    case ArgParseArgument::INPUT_FILE:
        typeName = "inputfile";
        break;

    case ArgParseArgument::OUTPUT_FILE:
        typeName = "outputfile";
        break;

    case ArgParseArgument::INPUTPREFIX:
        typeName = "inputprefix";
        break;

    case ArgParseArgument::OUTPUT_PREFIX:
        typeName = "outputprefix";
        break;

    default:
        typeName = "unknown";
        break;
    }

    return typeName;
}

// ----------------------------------------------------------------------------
// Function isListArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isListArgument
 * @headerfile <seqan/arg_parse.h>
 *
 * @brief Returns whether the argument can be given more than one time.
 *
 * @signature bool isListArgument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it can be given multiple times, <tt>false</tt> otherwise.
 */

inline bool isListArgument(ArgParseArgument const & me)
{
    return me._isListArgument;
}

// ----------------------------------------------------------------------------
// Function isStringArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isStringArgument
 * @headerfile <seqan/arg_parse.h>
 *
 * @brief Returns whether the argument is a string.
 *
 * @signature bool isStringArgument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a string, <tt>false</tt> otherwise.
 */

inline bool isStringArgument(ArgParseArgument const & me)
{
    return (me._argumentType == ArgParseArgument::STRING) ||
           (me._argumentType == ArgParseArgument::INPUT_FILE) ||
           (me._argumentType == ArgParseArgument::OUTPUT_FILE) ||
           (me._argumentType == ArgParseArgument::INPUTPREFIX) ||
           (me._argumentType == ArgParseArgument::OUTPUT_PREFIX) ;
}

// ----------------------------------------------------------------------------
// Function isIntegerArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isIntegerArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is an integer.
 *
 * @signature bool isIntegerArgument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is an integer, <tt>false</tt> otherwise.
 */

inline bool isIntegerArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::INTEGER;
}

// ----------------------------------------------------------------------------
// Function isInt64Argument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isInt64Argument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is a 64 bit integer.
 *
 * @signature bool isInt64Argument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a 64 bit integer, <tt>false</tt> otherwise.
 */

inline bool isInt64Argument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::INT64;
}

// ----------------------------------------------------------------------------
// Function isDoubleArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isDoubleArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is a double integer.
 *
 * @signature bool isDoubleArgument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a double argument, <tt>false</tt> otherwise.
 */

inline bool isDoubleArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::DOUBLE;
}

// ----------------------------------------------------------------------------
// Function isInputFileArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isInputFileArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is a input file.
 *
 * @signature bool isInputFileArgument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a input file argument, <tt>false</tt> otherwise.
 */

inline bool isInputFileArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::INPUT_FILE;
}

// ----------------------------------------------------------------------------
// Function isOutputFileArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isOutputFileArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is a output file.
 *
 * @signature bool isOutputFileArgument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a output file argument, <tt>false</tt> otherwise.
 */

inline bool isOutputFileArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::OUTPUT_FILE;
}

// ----------------------------------------------------------------------------
// Function isOutputPrefixArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isOutputPrefixArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is an output prefix.
 *
 * @signature bool isOutputPrefixArgument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is an output prefix argument, <tt>false</tt> otherwise.
 */

inline bool isOutputPrefixArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::OUTPUT_PREFIX;
}

// ----------------------------------------------------------------------------
// Function isOutputFileArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isInputPrefixArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is an input prefix argument.
 *
 * @signature bool isInputPrefixArgument(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is an input prefix argument, <tt>false</tt> otherwise.
 */

inline bool isInputPrefixArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::INPUTPREFIX;
}

// ----------------------------------------------------------------------------
// Function getArgumentLabel()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#getArgumentLabel
 * @headerfile <seqan/arg_parse.h>
 * @brief Return argument label.
 *
 * @signature std::string getArgumentLabel(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return std::string The argument label as a STL string.
 */

inline std::string const getArgumentLabel(ArgParseArgument const & me)
{
    if (me._argumentLabel != "")
    {
        return me._argumentLabel;
    }
    else
    {
        // infer from argument type
        std::string baseLabel = "";
        if (isInputFileArgument(me) || isOutputFileArgument(me))
            baseLabel = "FILE";
        else if (isInputPrefixArgument(me) || isOutputPrefixArgument(me))
            baseLabel = "PREFIX";
        else if (isStringArgument(me))
            baseLabel = "STR";
        else if (isIntegerArgument(me) || isDoubleArgument(me))
            baseLabel = "NUM";

        std::string finalLabel;

        if (me._numberOfValues != 1)
        {
            for (unsigned i = 0; i < me._numberOfValues; ++i)
            {
                if (i != 0)
                    append(finalLabel, " ");
                append(finalLabel, baseLabel);
            }
        }
        else if (isListArgument(me))
            finalLabel = baseLabel;                         // maybe we want to customize list labels
        else
            finalLabel = baseLabel;

        return finalLabel;
    }
}

// ----------------------------------------------------------------------------
// Helper Function _intervalAssert()
// ----------------------------------------------------------------------------

// this methods ensures that the given arguments define a non emtpy value interval
// otherwise it will trigger a SEQAN_CHECK failure
template <typename TIntervalBorder>
inline void _intervalAssert(const std::string minValueAsString, const std::string maxValueAsString)
{
    if (minValueAsString != "" && maxValueAsString != "")
        SEQAN_CHECK(_cast<TIntervalBorder>(minValueAsString) < _cast<TIntervalBorder>(maxValueAsString),
                    "The interval [%s:%s] is empty. Please specify a valid, non-empty interval.",
                    minValueAsString.c_str(),
                    maxValueAsString.c_str());
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#setMinValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Set smallest allowed value for the argument.
 *
 * @signature void setMinValue(arg, minValue);
 *
 * @param[in,out] arg      The ArgParseArgument to set the smallest value of.
 * @param[in]     minValue The smallest value to set (<tt>std::string</tt>).
 */

inline void setMinValue(ArgParseArgument & me, const std::string minValue)
{
    if (isDoubleArgument(me))
    {
        SEQAN_CHECK(_isCastable<double>(minValue), "The maximal value for a double argument must be double.");
        _intervalAssert<double>(minValue, me.maxValue);
        me.minValue = minValue;
    }
    else if (isIntegerArgument(me))
    {
        SEQAN_CHECK(_isCastable<int>(minValue), "The maximal value for an integer argument must be an integer");
        _intervalAssert<int>(minValue, me.maxValue);
        me.minValue = minValue;
    }
    else if (isInt64Argument(me))
    {
        SEQAN_CHECK(_isCastable<__int64>(minValue), "The maximal value for a 64 integer argument must be a 64 bit integer");
        _intervalAssert<__int64>(minValue, me.maxValue);
        me.minValue = minValue;
    }
    else
        SEQAN_FAIL("min/max values are not applicable to non numeric arguments");
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#setMaxValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Set smallest allowed value for the argument.
 *
 * @signature void setMaxValue(arg, maxValue);
 *
 * @param[in,out] arg      The ArgParseArgument to set the smallest value of.
 * @param[in]     maxValue The largest value to set (<tt>std::string</tt>).
 */

inline void setMaxValue(ArgParseArgument & me, const std::string maxValue)
{
    if (isDoubleArgument(me))
    {
        SEQAN_CHECK(_isCastable<double>(maxValue), "The maximal value for a double argument must be double.");
        _intervalAssert<double>(me.minValue, maxValue);
        me.maxValue = maxValue;
    }
    else if (isIntegerArgument(me))
    {
        SEQAN_CHECK(_isCastable<int>(maxValue), "The maximal value for an integer argument must be an integer");
        _intervalAssert<int>(me.minValue, maxValue);
        me.maxValue = maxValue;
    }
    else if (isInt64Argument(me))
    {
        SEQAN_CHECK(_isCastable<int>(maxValue), "The maximal value for a 64 bit integer argument must be an 64 bit integer");
        _intervalAssert<int>(me.minValue, maxValue);
        me.maxValue = maxValue;
    }
    else
        SEQAN_FAIL("min/max values are not applicable to non numeric arguments");
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#setValidValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Set list of valid values.
 *
 * @signature void setValidValues(arg, values);
 *
 * @param[in,out] arg    The ArgParseArgument to set the valid values for.
 * @param[in]     values Either a <tt>std::string</tt> containing all valid entries, separated by spaces or a
 *                       <tt>std::vector&lt;std::string&gt;</tt> with the valid entries.
 *
 * If the argument is of type string then the list of valid values is the case-sensitive list of string values
 * allowed for this argument.  If it is an input or output file then the list of valid values is a list of
 * case-insentive file extensions identifying the allowed types.
 *
 * @section Examples
 *
 * An example of setting allowed values for a string option.
 *
 * @code{.cpp}
 * seqan::ArgParseArgument stringArg(seqan::ArgParseArgument::STRING);
 * setValidValues(stringArg, "one two three");  // one of {"one", "two", "three"}
 *
 * std::vector<std::string> values;
 * values.push_back("four");
 * values.push_back("five");
 * setValidValues(stringArg, values);  // one of {"four", "five"}
 * @endcode
 *
 * An example for an input file option.  Note that by changing <tt>INPUT_FILE</tt> to <tt>OUTPUT_FILE</tt> below,
 * the example would be the same for output files.
 *
 * @code{.cpp}
 * seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUT_FILE);
 * setValidValues(fileArg, "fq fastq");  // file must end in ".fq" or ".fastq"
 *
 * std::vector<std::string> values;
 * values.push_back("sam");
 * values.push_back("bam");
 * setValidValues(fileArg, values);  // file must end in ".sam" or ".bam"
 * @endcode
 */

inline void setValidValues(ArgParseArgument & me, std::vector<std::string> const & values)
{
    if (isDoubleArgument(me) || isIntegerArgument(me))
        SEQAN_FAIL("ArgParseArgument does not support setting valid values for numeric arguments.");

    me.validValues = values;
}

inline void setValidValues(ArgParseArgument & me, std::string const & valuesString)
{
    // convert array to String<std::string>
    std::vector<std::string> values;
    std::string current_argument;

    for (std::string::const_iterator ch  = valuesString.begin(); ch != valuesString.end(); ++ch)
    {
        if (*ch == ' ')
        {
            appendValue(values, current_argument);
            current_argument = "";
        }
        else
        {
            append(current_argument, *ch);
        }
    }
    if (current_argument != "")
        appendValue(values, current_argument);

    setValidValues(me, values);
}

// ----------------------------------------------------------------------------
// Function setHelpText()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#setHelpText
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the help text for an ArgParseArgument.
 *
 * @signature void setHelpText(arg, text);
 *
 * @param[in,out] arg  The ArgParseArgument to set the help text for.
 * @param[in]     text The text to display as the description of the argument (<tt>std::string</tt>).
 */

inline void setHelpText(ArgParseArgument & me, std::string const & text)
{
    me._helpText = text;
}

// ----------------------------------------------------------------------------
// Helper Function _isInInterval()
// ----------------------------------------------------------------------------

// check if the given value is in the provided interval
template <typename TTarget, typename TString>
inline bool _isInInterval(TString value, TString lowerIntervalBound, TString upperIntervalBound)
{
    bool isInInterval = true;

    if (lowerIntervalBound != "")
        isInInterval &= (_cast<TTarget>(lowerIntervalBound) <= _cast<TTarget>(value));
    if (upperIntervalBound != "")
        isInInterval &= (_cast<TTarget>(value) <= _cast<TTarget>(upperIntervalBound));

    return isInInterval;
}

// ----------------------------------------------------------------------------
// Helper Function _checkNumericArgument()
// ----------------------------------------------------------------------------

// test if the values can be assigned to the option and is in the given boundaries
template <typename TNumerical>
inline void _checkNumericArgument(ArgParseArgument const & me, std::string const & value)
{
    if (!_isCastable<TNumerical>(value))
    {
        std::stringstream what;
        what << "the given value '" << value << "' cannot be casted to " << _typeToString(me);
        SEQAN_THROW(ParseError(what.str()));
    }

    if (!_isInInterval<TNumerical>(value, me.minValue, me.maxValue))
    {
        std::stringstream what;
        what << "the given value '" << value << "' is not in the interval ["
             << (me.minValue != "" ? me.minValue : "-inf") << ":"
             << (me.maxValue != "" ? me.maxValue : "+inf") << "]";

        SEQAN_THROW(ParseError(what.str()));
    }
}

// ----------------------------------------------------------------------------
// Helper Function _compareExtension()
// ----------------------------------------------------------------------------

inline bool _compareExtension(std::string const & str, std::string const & ext)
{
    std::string str_ext = str.substr(str.size() - ext.size());
    for (size_t i = 0; i < str_ext.size() && i < ext.size(); ++i)
    {
        if (tolower(str_ext[i]) != tolower(ext[i]))
            return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Helper Function _checkStringRestrictions()
// ----------------------------------------------------------------------------

// The parameter i gives the index of the value in the argument.

inline void _checkStringRestrictions(ArgParseArgument const & me, std::string const &value,
                                     unsigned i)
{
    typedef std::vector<std::string>::const_iterator TVectorIterator;

    // we only check valid values for files and string arguments, but not for prefix arguments
    if (!empty(me.validValues) && !(isInputPrefixArgument(me) || isOutputPrefixArgument(me)))
    {
        // The file name "-" is reserved for stdin or stdout
        if ((isInputFileArgument(me) || isOutputFileArgument(me)) && value == "-")
            return;

        // Allow the filename to be a pipe (without checking its extension)
        if (isInputFileArgument(me) && _isPipe(value.c_str()))
            return;

        bool isContained = false;
        for (TVectorIterator validValue = me.validValues.begin();
             validValue != me.validValues.end();
             ++validValue)
        {
            // if it is an input or output file, we only check the file endings
            if (isInputFileArgument(me) || isOutputFileArgument(me))
            {
                if (length(*validValue) > length(getFileExtension(me, i)))
                    continue;
                else
                    isContained |= _compareExtension(getFileExtension(me, i), *validValue);
            }
            else
            {
                isContained |= (*validValue == value);
            }
            if (isContained)
                break;
        }
        if (!isContained)
        {
            std::stringstream what;
            if (isInputFileArgument(me) || isOutputFileArgument(me))
                what << "the given path '" << value << "' does not have one of the valid file extensions [";
            else
                what << "the given value '" << value << "' is not in the list of allowed values [";
            for (TVectorIterator validValue = me.validValues.begin();
                 validValue != me.validValues.end();
                 ++validValue)
            {
                if (validValue != me.validValues.begin())
                    what << ", ";
                what << ((isInputFileArgument(me) || isOutputFileArgument(me)) ? "*" : "") << *validValue;
            }
            what << "]";
            if (i < me._fileExtensions.size())
                what << "; the file extension was overridden to be '" << getFileExtension(me, i) << "'";
            SEQAN_THROW(ParseError(what.str()));
        }
    }
}

// ----------------------------------------------------------------------------
// Function _checkValue()
// ----------------------------------------------------------------------------

// Check value before or after assignment.
//
// The parameter i gives the index of the value to check for overriding the extension in case of file arguments.

inline void _checkValue(ArgParseArgument const & me, std::string val, unsigned i = 0)
{
    // type checks
    if (isIntegerArgument(me))
        _checkNumericArgument<int>(me, val);

    if (isInt64Argument(me))
        _checkNumericArgument<__int64>(me, val);

    if (isDoubleArgument(me))
        _checkNumericArgument<double>(me, val);

    // check valid values
    if (isStringArgument(me))
        _checkStringRestrictions(me, val, i);
}

inline void _checkValue(ArgParseArgument const & me)
{
    unsigned i = 0;
    for (std::vector<std::string>::const_iterator it = me.value.begin(); it != me.value.end(); ++it, ++i)
        _checkValue(me, *it, i);
}

// ----------------------------------------------------------------------------
// Function _assignArgumentValue()
// ----------------------------------------------------------------------------

inline void _assignArgumentValue(ArgParseArgument & me, std::string const & value)
{
    // assignment
    if (isListArgument(me)) // just append
    {
        appendValue(me.value, value, Exact());
    }
    else
    {
        // check if we already set all expected arguments
        if (length(me.value) == me._numberOfValues)
            clear(me.value);
        appendValue(me.value, value, Exact());
    }
}

// ----------------------------------------------------------------------------
// Function getArgumentValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#getArgumentValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Return the value of the argument.
 *
 * @signature std::string getArgumentValue(arg[, argNo]);
 *
 * @param[in,out] arg   The ArgParseArgument to query.
 * @param[in]     argNo In case that the ArgParseArgument allowed multiple values, give the index of the argument
 *                      that you want to retrieve (<tt>unsigned</tt>, starts at 0).
 *
 * @return std::string Const-reference to the argument value.
 */

inline std::string const & getArgumentValue(ArgParseArgument const & me, unsigned argNo)
{
    SEQAN_CHECK(argNo < me.value.size() || argNo < me.defaultValue.size(),
                "ArgParseArgument: No value set for index %d", argNo);

    if (argNo < me.value.size())
        return me.value[argNo];
    else
        return me.defaultValue[argNo];
}

inline std::string const & getArgumentValue(ArgParseArgument const & me)
{
    return getArgumentValue(me, 0);
}

// ----------------------------------------------------------------------------
// Function getArgumentValues()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#getArgumentValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Return all values of the argument.
 *
 * @signature std::vector<std::string> getArgumentValue(arg);
 *
 * @param[in] arg   The ArgParseArgument to query.
 *
 * @return std::vector<std::string> Const-reference to the argument values.
 */

inline std::vector<std::string> const & getArgumentValues(ArgParseArgument const & me)
{
    if (!me.value.empty())
        return me.value;
    else
        return me.defaultValue;
}

// ----------------------------------------------------------------------------
// Function hasArgumentValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#hasArgumentValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Return whether a value is available.
 *
 * @signature bool hasValue(arg[, pos]);
 *
 * @param[in] arg The ArgParseArgument to query.
 * @param[in] pos The position of the argument in case of being a list (<tt>unsigned</tt>, 0-based, default is 0).
 *
 * @return bool <tt>true</tt> if <tt>pos</tt> is less than the size and the argument is non-empty.
 */

inline bool hasValue(ArgParseArgument const & arg, unsigned position)
{
    return arg.value.size() > position || arg.defaultValue.size() > position;
}

inline bool hasValue(ArgParseArgument const & arg)
{
    return hasValue(arg, 0);
}

// ----------------------------------------------------------------------------
// Function isSet()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isSet
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns true if a value was assigned to the argument.
 *
 * @signature bool isSet(arg):
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if a value was assigned, <tt>false</tt> otherwise.
 */

inline bool isSet(ArgParseArgument const & me)
{
    return !me.value.empty();
}

// ----------------------------------------------------------------------------
// Function hasDefault()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#hasDefault
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument has a default value.
 *
 * @signature bool hasDefault(arg);
 *
 * @param[in] arg The argument to query.
 *
 * @return bool <tt>true</tt> if the argument has a default value and <tt>false</tt> if not.
 */

inline bool hasDefault(ArgParseArgument const & me)
{
    return !me.defaultValue.empty();
}

// ----------------------------------------------------------------------------
// Function numberOfArguments
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#numberOfAllowedValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns the number of allowed values for this ArgParseArgument.
 *
 * @signature unsigned numberOfAllowedValues(arg);
 *
 * @param[in] arg The ArgParseArgument to query.
 *
 * @return unsigned The number of allowed values.
 */

inline unsigned numberOfAllowedValues(ArgParseArgument const & me)
{
    return me._numberOfValues;
}

// ----------------------------------------------------------------------------
// Function getFileExtension()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#getFileExtension
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns the file extension for the given file argument.
 *
 * Only valid when argument is an INPUT_FILE or OUTPUT_FILE.
 *
 * Halts the program if not an input or output file argument.
 *
 * Can be overridden with special hidden options.
 * For arguments, you can pass <tt>--arg-&lt;num&gt;-file-ext</tt> for argument <tt>num</tt>.
 * For parameters, you can pass <tt>--&lt;param-name&gt;-file-ext</tt> for the option named <tt>param-name</tt>.
 *
 * @signature std::string getFileExtension(arg[, pos]);
 *
 * @param[in] arg The ArgParseArgument to query.
 * @param[in] pos The position of the value to retrieve if multiple values (<tt>unsigned</tt>).
 *
 * @return std::string The file extension, empty if no extension or not set.
 */

inline std::string getFileExtension(ArgParseArgument const & me, unsigned pos = 0)
{
    if (me._argumentType != ArgParseArgument::INPUT_FILE &&
        me._argumentType != ArgParseArgument::OUTPUT_FILE)
        SEQAN_FAIL("Cannot get file extension from non-file argument/option.");

    // Short-circuit to override file extension if set.
    if (pos < me._fileExtensions.size())
    {
        std::string result = me._fileExtensions[pos];
        if (!result.empty() && result[0] != '.')
            result.insert(0, ".");
        return result;
    }

    // Get argument value and break if empty.
    std::string value = getArgumentValue(me, pos);
    if (value.empty())
        return "";

    // If there is a list of valid values then we look for each of these in the path.
    if (!me.validValues.empty())
    {
        for (unsigned i = 0; i < length(me.validValues); ++i)
        {
            unsigned len = std::min(me.validValues[i].size(), value.size());
            std::string tmp = value.substr(value.size() - len);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            if (tmp == me.validValues[i])
                return me.validValues[i];
        }
        return "";
    }

    // Fall back to finding position of last dot.  Return suffix if found and empty string if not.
    size_t dotPos = value.find_last_of('.');
    if (dotPos == std::string::npos)
        return "";
    return value.substr(dotPos + 1);
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_
