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

#ifndef SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_
#define SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_

#include <seqan/arg_parse/tool_doc.h>
#include <seqan/arg_parse/argument_parser.h>

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function getAppName()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getAppName
 * @brief Return program name of ArgumentParser.
 *
 * @signature TCharStringRef getAppName(parser);
 *
 * @param[in] parser The ArgumentParser to get the app name for.
 *
 * @return TCharStringRef The app name, const-ref to @link CharString @endlink.
 */

inline CharString const & getAppName(ArgumentParser const & parser)
{
    return getName(parser._toolDoc);
}

// ----------------------------------------------------------------------------
// Helper Function _parseAppName()
// ----------------------------------------------------------------------------

inline void _parseAppName(ArgumentParser & parser, std::string const & candidate)
{
    //IOREV _notio_ irrelevant for io-revision
    int i = length(candidate) - 1;

    for (; i >= 0; --i)
        if (candidate[i] == '\\' || candidate[i] == '/')
            break;

    setName(parser._toolDoc, candidate.substr(i + 1));
}

// ----------------------------------------------------------------------------
// Helper Function _addLine()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addLine
 * @brief Adds a line of text to the help output of the ArgumentParser.
 *
 * The line of text will be added to the block of the options.
 *
 * @signature void addLine(parser, line);
 *
 * @param[in,out] parser The ArgumentParser to add the line to.
 * @param[in]     line   The line of text to add, @link StringConcept @endlink of <tt>char</tt>.
 */

template <typename TString>
inline void addLine(ArgumentParser & me, TString const & line)
{
    addOption(me, ArgParseOption("", "", line));
}

// ----------------------------------------------------------------------------
// Function addSection()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addSection
 * @brief Begins a new section of the option block of the ArgumentParser help output.
 *
 * @signature void addSection(parser, title);
 *
 * @param[in,out] parser The ArgumentParser to add the line to.
 * @param[in]     title  The title to add, @link StringConcept @endlink of <tt>char</tt>.
 *
 * @code{.cpp}
 * ArgumentParser parser;
 *
 * [...] // init parser
 *
 * addSection(parser, "In-/Output-Options");
 * addOption("i", ... );
 * addOption("o", ... );
 *
 * addSection(parser, "Other Options");
 * addOption("x", ... );
 * @endcode
 */

template <typename TString>
inline void addSection(ArgumentParser & me, TString const & line)
{
    addLine(me, "");
    addLine(me, line);
}

// ----------------------------------------------------------------------------
// Function addUsageLine()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addUseLine
 * @brief Adds a line of text to the usage output of the ArgumentParser.
 *
 * @signature void addUsageLine(parser, line);
 *
 * @param[in,out] parser The ArgumentParser to add the line to.
 * @param[in]     line   The line to add, a <tt>std::string</tt>.
 */

inline void addUsageLine(ArgumentParser & me, std::string const & line)
{
    me._usageText.push_back(line);
}

// ----------------------------------------------------------------------------
// Helper Function _addUsage()
// ----------------------------------------------------------------------------

inline void _addUsage(ToolDoc & toolDoc, ArgumentParser const & me)
{
    for (unsigned i = 0; i < length(me._usageText); ++i)
    {
        std::string text = "\\fB";
        append(text, getAppName(me));
        append(text, "\\fP ");
        append(text, me._usageText[i]);
        addText(toolDoc, text, false);
    }
}

// ----------------------------------------------------------------------------
// Function addDescription()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addDescription
 * @brief Appends a description paragraph to the ArgumentParser documentation.
 *
 * @signature void addDescription(parser, description);
 *
 * @param[in,out] parser      The ArgumentParser to add the line to.
 * @param[in]     description The description text, a <tt>std::string</tt>.
 */

inline void addDescription(ArgumentParser & me, std::string const & description)
{
    me._description.push_back(description);
}

// ----------------------------------------------------------------------------
// Function setAppName()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setAppName
 * @brief Sets application name of ArgumentParser.
 *
 * @signature void setAppName(parser, name);
 *
 * @param[in,out] parser The ArgumentParser to set the name of.
 * @param[in]     name   The application name, <tt>std::string</tt>.
 */

inline void setAppName(ArgumentParser & me, std::string const & name)
{
    setName(me._toolDoc, name);
}

// ----------------------------------------------------------------------------
// Function setShortDescription()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setShortDescription
 * @brief Sets shortDescription of ArgumentParser.
 *
 * @signature void setShortDescription(parser, desc);
 *
 * @param[in,out] parser The ArgumentParser to set the short description of.
 * @param[in]     desc   The short description, <tt>std::string</tt>.
 */

inline void setShortDescription(ArgumentParser & me, std::string const & description)
{
    setShortDescription(me._toolDoc, description);
}

// ----------------------------------------------------------------------------
// Function getShortDescription()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getShortDescription
 * @brief Returns the short description.
 *
 * @signature CharString getShortDescription(parser);
 *
 * @param[in,out] parser The ArgumentParser to get short description for.
 *
 * @return CharString A @link CharString @endlink with the short description.
 */

inline CharString getShortDescription(ArgumentParser const & me)
{
    return getShortDescription(me._toolDoc);
}

// ----------------------------------------------------------------------------
// Function setVersion()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setVersion
 * @brief Sets version of ArgumentParser.
 *
 * @signature void setVersion(parser, version);
 *
 * @param[in,out] parser  The ArgumentParser to set the version of.
 * @param[in]     version The version string to set, <tt>std::string</tt>.
 */

inline void setVersion(ArgumentParser & me, std::string const & versionString)
{
    setVersion(me._toolDoc, versionString);
    if (!hasOption(me, "version"))
        addOption(me, ArgParseOption("", "version", "Display version information."));
}

// --------------------------------------------------------------------------
// Function getVersion()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getVersion
 * @brief Returns the version string.
 *
 * @signature TCharStringRef getVersion(parser);
 *
 * @param[in,out] parser The ArgumentParser to get the version string from.
 *
 * @return TCharString A const-ref to a @link CharString @endlink with the version string.
 */

inline CharString const & getVersion(ArgumentParser const & me)
{
    return getVersion(me._toolDoc);
}

// --------------------------------------------------------------------------
// Function setCategory()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setCategory
 * @brief Sets category of ArgumentParser.
 *
 * @signature void setCategory(parser, category);
 *
 * @param[in,out] parser  The ArgumentParser to set the category of.
 * @param[in]     category The category to set, <tt>std::string</tt>.
 */

inline void setCategory(ArgumentParser & parser, CharString const & category)
{
    setCategory(parser._toolDoc, category);
}

// --------------------------------------------------------------------------
// Function getCategory()
// --------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#getCategory
 * @brief Returns the category.
 *
 * @signature TCharStringRef getCategory(parser);
 *
 * @param[in,out] parser The ArgumentParser to get the category from.
 *
 * @return TCharString A const-ref to a @link CharString @endlink with the category.
 */

inline CharString const & getCategory(ArgumentParser const & parser)
{
    return getCategory(parser._toolDoc);
}

// ----------------------------------------------------------------------------
// Function setDate()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#setDate
 * @brief Sets date string of ArgumentParser.
 *
 * @signature void setDate(parser, date);
 *
 * @param[in,out] parser The ArgumentParser to set the date string of.
 * @param[in]     date   The date string to set, <tt>std::string</tt>.
 */

inline void setDate(ArgumentParser & me, std::string const & date)
{
    setDate(me._toolDoc, date);
}

// ----------------------------------------------------------------------------
// Function addTextSection()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addTextSection
 * @brief Add a text section to the ArgumentParser.
 *
 * @signature void addTextSection(parser, title);
 *
 * @param[in,out] parser The ArgumentParser to add the text section title to.
 * @param[in]     title  The section title to add, <tt>std::string</tt>.
 */

inline void addTextSection(ArgumentParser & me, std::string const & title)
{
    addSection(me._toolDoc, title);
}

// ----------------------------------------------------------------------------
// Function addTextSubSection()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addTextSubSection
 * @brief Add a text sub section to the ArgumentParser.
 *
 * @signature void addTextSubSection(parser, title);
 *
 * @param[in,out] parser The ArgumentParser add the subsection title to of.
 * @param[in]     title  The sub section title to add, <tt>std::string</tt>.
 */

inline void addTextSubSection(ArgumentParser & me, std::string const & title)
{
    addSubSection(me._toolDoc, title);
}

// ----------------------------------------------------------------------------
// Function addText()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addText
 * @brief Add text to an ArgumentParser.
 *
 * @signature void addText(parser, text);
 *
 * @param[in,out] parser ArgumentParser to add text to.
 * @param[in]     text   The <tt>std::string</tt> to add to the parser.
 */

inline void addText(ArgumentParser & me, std::string const & text)
{
    addText(me._toolDoc, text);
}

// ----------------------------------------------------------------------------
// Function addListItem()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#addListItem
 * @brief Appends a list item to the ArgumentParser
 *
 * @signature void addListItem(parser, item, description);
 *
 * @param[in,out] parser      The ArgumentParser to add the list item to.
 * @param[in]     item        The item to add, <tt>std::string</tt>.
 * @param[in]     description The item to add, <tt>std::string</tt>.
 */

inline void addListItem(ArgumentParser & me, std::string const & item, std::string const & description)
{
    addListItem(me._toolDoc, item, description);
}

// ----------------------------------------------------------------------------
// Function printShortHelp()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#printShortHelp
 * @brief Prints a short help message for the parser to a stream.
 *
 * @signature void printShortHelp(parser, out);
 *
 * @param[in,out] parser The ArgumentParser to print help for.
 * @param[in,out] out    The <tt>std::ostream</tt> to print help to.
 */

inline void printShortHelp(ArgumentParser const & me, std::ostream & stream)
{
    // TODO: maybe we can get this a bit prettier
    ToolDoc shortDoc(me._toolDoc);
    clearEntries(shortDoc);

    _addUsage(shortDoc, me);

    std::stringstream shortHelp;
    shortHelp << "Try '" << getAppName(me) << " --help' for more information.\n";
    addText(shortDoc, shortHelp.str());

    print(stream, shortDoc, "txt");
}

inline void printShortHelp(ArgumentParser const & me)
{
    printShortHelp(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function printVersion()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#printVersion
 * @brief Prints the version information of the parser to a stream.
 *
 * @signature void printVersion(parser, stream);
 *
 * @param[in,out] parser The ArgumenParser to print for.
 * @param[in,out] stream The <tt>std::ostream</tt> to print to.
 */

inline void printVersion(ArgumentParser const & me, std::ostream & stream)
{
    stream << getAppName(me) << " version " << getVersion(me) << std::endl;
}

inline void printVersion(ArgumentParser const & me)
{
    printVersion(me, std::cerr);
}

// ----------------------------------------------------------------------------
// Function _addNumericalRestriction()
// ----------------------------------------------------------------------------


inline void _addNumericalRestriction(std::string & text, ArgParseOption const & opt)
{
    // expand min/max restrictions
    if (!empty(opt.minValue) || !empty(opt.maxValue))
    {
        append(text, " In range [");

        if (!empty(opt.minValue))
            append(text, opt.minValue);
        else
            append(text, "-inf");

        append(text, "..");

        if (!empty(opt.maxValue))
            append(text, opt.maxValue);
        else
            append(text, "inf");

        append(text, "].");
    }
}

// ----------------------------------------------------------------------------
// Function _expandList()
// ----------------------------------------------------------------------------

// expands the given vector as text in the form v1, v2, and v3, while respecting
// the size with respect to the used commas and "and"s
inline void _expandList(std::string & text, std::vector<std::string> const & list)
{
    for (std::vector<std::string>::size_type i = 0; i < list.size(); ++i)
    {
        if (i + 1 == list.size() && list.size() == 2u)
            append(text, " and ");
        else if (i + 1 == list.size()  && list.size() > 2u)
            append(text, ", and ");
        else if (i != 0)
            append(text, ", ");

        append(text, "\\fI");
        append(text, list[i]);
        append(text, "\\fP");

    }
}

// ----------------------------------------------------------------------------
// Function _addDefaultValues()
// ----------------------------------------------------------------------------

inline void _addDefaultValues(std::string & text, ArgParseOption const & opt)
{
    if (!empty(opt.defaultValue) && !isBooleanOption(opt))
    {
        append(text, " Default: ");
        _expandList(text, opt.defaultValue);
        append(text, ".");
    }
}

// ----------------------------------------------------------------------------
// Function _addValidValuesRestrictions()
// ----------------------------------------------------------------------------

inline void _addValidValuesRestrictions(std::string & text, ArgParseOption const & opt)
{
    if (!empty(opt.validValues) && !isBooleanOption(opt))
    {
        if (isInputFileArgument(opt) || isOutputFileArgument(opt))
        {
            append(text, " Valid filetype");

            if (opt.validValues.size() > 1)
                append(text, "s are: ");
            else
                append(text, " is: ");
        }
        else
        {
            append(text, " One of ");
        }

        _expandList(text, opt.validValues);
        append(text, ".");
    }
}

// ----------------------------------------------------------------------------
// Function printHelp()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Parameter order.

/*!
 * @fn ArgumentParser#printHelp
 * @brief Prints the complete help message for the parser.
 *
 * @signature void printHelp(parser, out, format);
 *
 * @param[in,out] parser The ArgumentParser print the help for.
 * @param[out]    out    The output stream to print to (<tt>std::ostream</tt>).
 * @param[in]     format The format to print, one of "html", "man", and "txt".
 */

inline void printHelp(ArgumentParser const & me, std::ostream & stream, CharString const & format)
{
    ToolDoc toolDoc(me._toolDoc);
    clearEntries(toolDoc);  // We will append me._toolDoc later.

    // Build synopsis section.
    addSection(toolDoc, "Synopsis");
    _addUsage(toolDoc, me);

    // Add description to tool documentation.
    addSection(toolDoc, "Description");
    for (unsigned i = 0; i < me._description.size(); ++i)
        addText(toolDoc, me._description[i]);

    // Add options to description section.
    for (unsigned i = 0; i < length(me.optionMap); ++i)
    {
        ArgParseOption const & opt = me.optionMap[i];
        if (empty(opt.shortName) && empty(opt.longName))  // this is not an option but a text line
        {
            if (empty(opt._helpText))  // TODO(holtgrew): Should go away in future.
                continue;  // Skip empty lines.

            // Is command line parser section, maps to ToolDoc subsection.
            std::string title = opt._helpText;
            append(title, ":");
            addSubSection(toolDoc, title);
        }
        else if (!isHidden(opt))
        {
            // Build list item term.
            std::string term;
            if (!empty(opt.shortName))
            {
                term = "\\fB-";
                append(term, opt.shortName);
                append(term, "\\fP");
            }
            if (!empty(opt.shortName) && !empty(opt.longName))
                append(term, ", ");
            if (!empty(opt.longName))
            {
                append(term, "\\fB--");
                append(term, opt.longName);
                append(term, "\\fP");
            }
            // Get arguments, autogenerate if necessary.
            std::string arguments = getArgumentLabel(opt);

            // Write arguments to term line -> only exception, boolean flags
            if (!empty(arguments))
            {
                // Tokenize argument names.
                std::istringstream iss(toCString(arguments));
                std::vector<std::string> tokens;
                std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                          std::back_inserter<std::vector<std::string> >(tokens));
                // Append them, formatted in italic.
                for (unsigned i = 0; i < length(tokens); ++i)
                {
                    append(term, " \\fI");
                    append(term, tokens[i]);
                    append(term, "\\fP");
                }
            }

            std::string helpText = opt._helpText;

            // expand min/max restrictions
            _addNumericalRestriction(helpText, opt);

            // expand validValues restrictions
            _addValidValuesRestrictions(helpText, opt);

            // expand defaultValue
            _addDefaultValues(helpText, opt);

            // Add list item.
            addListItem(toolDoc, term, helpText);
        }
    }

    append(toolDoc, me._toolDoc);
    print(stream, toolDoc, format);
}

inline void printHelp(ArgumentParser const & me, std::ostream & stream)
{
    printHelp(me, stream, "txt");
}

inline void printHelp(ArgumentParser const & me)
{
    printHelp(me, std::cerr, "txt");
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_DOC_H_
