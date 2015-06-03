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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_ARG_PARSE_TOOL_DOC_H_
#define SEQAN_INCLUDE_ARG_PARSE_TOOL_DOC_H_

#include <iterator>

#include <seqan/misc/terminal.h>
#include <seqan/arg_parse/xml_support.h>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

class ToolDoc;
inline void clearEntries(ToolDoc & doc);
inline void append(ToolDoc & a, ToolDoc const & b);

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Class ToolDocEntry_
// --------------------------------------------------------------------------

// Base class for entries in the ToolDoc class.

class ToolDocEntry_
{
public:
    enum EntryType
    {
        SUBSECTION, SECTION, LINE, LIST_ITEM
    };

    virtual ~ToolDocEntry_() {}

    virtual EntryType getType() const = 0;
};

// --------------------------------------------------------------------------
// Class ToolDocSection_
// --------------------------------------------------------------------------

// A section in the ToolDoc class.

class ToolDocSection_ :
    public ToolDocEntry_
{
public:
    CharString _title;

    ToolDocSection_(ToolDocSection_ const & sec) :
        ToolDocEntry_(sec), _title(sec._title)
    {}

    ToolDocSection_(CharString const & title) :
        _title(title)
    {}

    virtual EntryType getType() const
    {
        return SECTION;
    }

};

// --------------------------------------------------------------------------
// Class ToolDocSubSection_
// --------------------------------------------------------------------------

// A subsection.

class ToolDocSubSection_ :
    public ToolDocSection_
{
public:
    ToolDocSubSection_(ToolDocSubSection_ const & sec) :
        ToolDocSection_(sec)
    {}

    ToolDocSubSection_(CharString const & title) :
        ToolDocSection_(title)
    {}

    virtual EntryType getType() const
    {
        return SUBSECTION;
    }

};

// --------------------------------------------------------------------------
// Class ToolDocLine__
// --------------------------------------------------------------------------

// A Line/Paragraph in ToolDoc documents.

class ToolDocLine_ :
    public ToolDocEntry_
{
public:
    CharString _text;
    // Evaluates to true if it is a paragraph and separated by .sp, false if it is a line, separated by .br.
    bool _isPar;

    ToolDocLine_(ToolDocLine_ const & line) :
        ToolDocEntry_(line), _text(line._text), _isPar(line._isPar)
    {}

    ToolDocLine_(CharString const & text, bool isParagraph) :
        _text(text), _isPar(isParagraph)
    {}

    virtual EntryType getType() const
    {
        return LINE;
    }

    bool isParagraph() const
    {
        return _isPar;
    }

};

// --------------------------------------------------------------------------
// Class ToolDocListItem_
// --------------------------------------------------------------------------

// A list item in ToolDoc documents.

class ToolDocListItem_ :
    public ToolDocEntry_
{
public:
    CharString _term;
    CharString _description;

    ToolDocListItem_(ToolDocListItem_ const & item) :
        ToolDocEntry_(item), _term(item._term), _description(item._description)
    {}

    ToolDocListItem_(CharString const & term, CharString const & description) :
        _term(term), _description(description)
    {}

    virtual EntryType getType() const
    {
        return LIST_ITEM;
    }

};

// --------------------------------------------------------------------------
// Class ToolDocPrinter_
// --------------------------------------------------------------------------

// Abstrace base class for ToolDoc printers.

class ToolDocPrinter_
{
public:
    virtual void print(std::ostream & stream, ToolDoc const & doc) = 0;
};

// --------------------------------------------------------------------------
// Class ManToolDocPrinter_
// --------------------------------------------------------------------------

// Print ToolDoc instance in man format.

class ManToolDocPrinter_ :
    public ToolDocPrinter_
{
public:
    virtual void print(std::ostream & stream, ToolDoc const & doc);
};

// --------------------------------------------------------------------------
// Class HtmlToolDocPrinter_
// --------------------------------------------------------------------------

// Print ToolDoc instance in HTML format.

class HtmlToolDocPrinter_ :
    public ToolDocPrinter_
{
public:
    void _maybeCloseList(std::ostream & stream, bool & isListOpen)
    {
        if (!isListOpen)
            return;

        stream << "</dl>\n";
        isListOpen = false;
    }

    void _maybeCloseParagraph(std::ostream & stream, bool & isParOpen)
    {
        if (!isParOpen)
            return;

        stream << "</p>\n";
        isParOpen = false;
    }

    // Converts formatting with \fI, \fB, and \fP to HTML.
    template <typename TSequence>
    TSequence _toHtml(TSequence const & input) const
    {
        TSequence buffer = xmlEscape(input);
        TSequence result;
        String<TSequence> openTags;

        typedef typename Iterator<TSequence const, Standard>::Type TIterator;
        TIterator endIt = end(input, Standard());
        for (TIterator it = begin(input, Standard()); it != endIt; goNext(it))
        {
            if (*it == '\\')
            {
                // Handle escape sequence, we interpret only "\-", "\fI", and "\fB".
                goNext(it);
                SEQAN_ASSERT_NOT(it == endIt);
                if (*it == '-')
                {
                    appendValue(result, *it);
                }
                else if (*it == 'f')
                {
                    goNext(it);
                    SEQAN_ASSERT_NOT(it == endIt);
                    if (*it == 'I')
                    {
                        appendValue(openTags, "em");
                        append(result, "<em>");
                    }
                    else if (*it == 'B')
                    {
                        appendValue(openTags, "strong");
                        append(result, "<strong>");
                    }
                    else if (*it == 'P')
                    {
                        SEQAN_ASSERT_NOT(empty(openTags));
                        append(result, "</");
                        append(result, back(openTags));
                        append(result, ">");
                        eraseBack(openTags);
                    }
                    else
                    {
                        append(result, "\\f");
                        appendValue(result, *it);
                    }
                }
                else
                {
                    appendValue(result, '\\');
                    appendValue(result, *it);
                }
            }
            else
            {
                appendValue(result, *it);
            }
        }

        return result;
    }

    virtual void print(std::ostream & stream, ToolDoc const & doc);
};

// --------------------------------------------------------------------------
// Class TextToolDocPrinter_
// --------------------------------------------------------------------------

// Print ToolDoc objects in Text format, suitable for console output.

class TextToolDocPrinter_ :
    public ToolDocPrinter_
{
public:
    // Stores the relevant parameters of the documentation on the screen.
    struct Layout_
    {
        unsigned screenWidth;
        unsigned defaultScreenWidth;
        unsigned maximalScreenWidth;
        unsigned minimalScreenWidth;
        unsigned leftPadding;
        unsigned centerPadding;
        unsigned rightPadding;
        unsigned leftColumnWidth;
        unsigned rightColumnWidth;
        unsigned rightColumnTab;

        Layout_() :
            screenWidth(0), defaultScreenWidth(80), maximalScreenWidth(120), minimalScreenWidth(40),
            leftPadding(4), centerPadding(2), rightPadding(2), leftColumnWidth(4), rightColumnWidth(0)
        {
            // Guess terminal screen width and set into layout.
            unsigned cols = 0, rows = 0;
            bool success = getTerminalSize(cols, rows);
            screenWidth = (success && cols > 0) ? cols : defaultScreenWidth;
            screenWidth = std::max(screenWidth, minimalScreenWidth);
            screenWidth = std::min(screenWidth, maximalScreenWidth);
            screenWidth -= rightPadding;

            rightColumnWidth = screenWidth - leftPadding - leftColumnWidth - centerPadding - rightPadding;
            rightColumnTab = leftPadding + leftColumnWidth + centerPadding;

            // std::cerr << "screen width\t" << screenWidth << std::endl;
        }

    };

    Layout_ _layout;

    virtual void print(std::ostream & stream, ToolDoc const & doc);

    void _printSection(std::ostream & stream, ToolDocSection_ const & section)
    {
        std::ostream_iterator<char> out(stream);
        stream << '\n' << _toText("\\fB");
        std::transform(begin(section._title), end(section._title), out, toupper);
        stream << _toText("\\fP") << '\n';
    }

    void _printSubSection(std::ostream & stream, ToolDocSection_ const & section)
    {
        std::ostream_iterator<char> out(stream);
        stream << '\n' << _toText("\\fB");
        std::fill_n(out, _layout.leftPadding / 2, ' ');
        stream << section._title << _toText("\\fP") << '\n';
    }

    void _printLine(std::ostream & stream, ToolDocLine_ const & line)
    {
        std::ostream_iterator<char> out(stream);
        std::fill_n(out, _layout.leftPadding, ' ');

        _printText(stream, line._text, _layout.leftPadding);
    }

    CharString _toText(CharString const & str) const
    {
        CharString result;

        typedef Iterator<CharString const, Rooted>::Type TIterator;
        for (TIterator it = begin(str, Rooted()); !atEnd(it); goNext(it))
        {
            if (*it == '\\')
            {
                // Handle escape sequence, we interpret only "\-", "\fI", and "\fB".
                goNext(it);
                SEQAN_ASSERT_NOT(atEnd(it));
                if (*it == '-')
                {
                    appendValue(result, *it);
                }
                else if (*it == 'f')
                {
                    goNext(it);
                    SEQAN_ASSERT_NOT(atEnd(it));
                    if (*it == 'I')
                    {
                        if (isTerminal())
                            append(result, "\033[4m");
                    }
                    else if (*it == 'B')
                    {
                        if (isTerminal())
                            append(result, "\033[1m");
                    }
                    else if (*it == 'P')
                    {
                        if (isTerminal())
                            append(result, "\033[0m");
                    }
                    else
                    {
                        append(result, "\\f");
                        appendValue(result, *it);
                    }
                }
                else
                {
                    appendValue(result, '\\');
                    appendValue(result, *it);
                }
            }
            else
            {
                appendValue(result, *it);
            }
        }

        return result;
    }

    void _printListItem(std::ostream & stream, ToolDocListItem_ const & listItem)
    {
        std::ostream_iterator<char> out(stream);

        // Print term.
        std::fill_n(out, _layout.leftPadding, ' ');
        stream << _toText(listItem._term);
        unsigned pos = _layout.leftPadding + length(listItem._term);
        if (pos + _layout.centerPadding > _layout.rightColumnTab)
        {
            stream << '\n';
            pos = 0;
        }
        std::fill_n(out, _layout.rightColumnTab - pos, ' ');

        _printText(stream, listItem._description, _layout.rightColumnTab);
    }

    // Returns width of text if printed, i.e. "\-" has length 1, "\fI", "\fB", "\fP" have length 0.
    unsigned _textWidth(CharString const & text)
    {
        unsigned result = 0;

        for (unsigned i = 0; i < length(text); ++i)
        {
            if (text[i] != '\\')
            {
                result += 1;
                continue;
            }

            if (i + 1 == length(text))
            {
                result += 1;  // Will print "\\".
                continue;
            }

            if (text[i + 1] == '\\' || text[i + 1] == '-')
            {
                i += 1;
                result += 1;
                continue;  // Will print '\\' or '-'.
            }

            if (i + 2 == length(text))
            {
                i += 1;
                result += 2;  // Will print two chars.
                continue;
            }

            if (text[i + 1] == 'f')
            {
                if (text[i + 2] == 'B' || text[i + 2] == 'I' || text[i + 2] == 'P')
                    i += 2;  // Skip f and {B, I, P}.
                else
                    result += 1;
            }
        }

        return result;
    }

    // Print text, must be at tab already.
    void _printText(std::ostream & stream, CharString const & text, unsigned tab)
    {
        unsigned pos = tab;
        std::ostream_iterator<char> out(stream);

        // Tokenize the text.
        std::istringstream iss(toCString(text));
        std::vector<std::string> tokens;
        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                  std::back_inserter<std::vector<std::string> >(tokens));

        // Print the text.
        SEQAN_ASSERT_LEQ(pos, tab);
        std::fill_n(out, tab - pos, ' ');  // go to tab

        pos = tab;
        typedef std::vector<std::string>::const_iterator TConstIter;
        for (TConstIter it = tokens.begin(); it != tokens.end(); ++it)
        {
            if (it == tokens.begin())
            {
                stream << _toText(*it);
                pos += _textWidth(*it);
                if (pos > _layout.screenWidth)
                {
                    stream << '\n';
                    std::fill_n(out, tab, ' ');
                    pos = tab;
                }
            }
            else
            {
                if (pos + 1 + _textWidth(*it) > _layout.screenWidth)
                {
                    // Would go over screen with next, print current word on next line.
                    stream << '\n';
                    fill_n(out, tab, ' ');
                    stream << _toText(*it);
                    pos = tab + _textWidth(*it);
                }
                else
                {
                    stream << ' ';
                    stream << _toText(*it);
                    pos += _textWidth(*it) + 1;
                }
            }
        }
        if (!empty(tokens))
            stream << '\n';
    }

};

// --------------------------------------------------------------------------
// Class ToolDoc
// --------------------------------------------------------------------------

/*!
 * @class ToolDoc
 * @implements AssignableConcept
 * @headerfile <seqan/arg_parse.>
 * @brief Container for string documentation on a command line tool.
 *
 * @signature class ToolDoc;
 *
 * @section Remarks
 *
 * This class is generally not used directly by the user but through @link ArgumentParser @endlink. It allows to store
 * and represent all information related to a command line tool that would normally go into a man page. It can be
 * printed to STL streams in different formats, currently plain text, HTML and man pages are supported.
 *
 * You can also use basic formatting in text. This formatting is tailored to the usage on the command line. Use
 * <tt>\fB</tt> to start bold font, <tt>\fI</tt> to start italic font and <tt>\fP</tt> to use the previous font (of
 * course, use correct escaping of the backslash in C strings, so use <tt>"\\fB"</tt>, <tt>"\\fI"</tt>, and
 * <tt>"\\fP"</tt> in your code.
 *
 * @section Examples
 *
 * The following shows a brief example of how to use @link ToolDoc @endlink.
 *
 * @code{.cpp}
 * ToolDoc doc;
 * setName(doc, "RazerS");
 * setShortDescription(doc, "Read mapping with controllable sensitivity.");
 * setDate(doc, "04 March 2012");
 * setVersion(doc, "1.0");
 * setCategory(doc, "Read Mapping");
 * setManTitle(doc, "SeqAn Apps Reference Manual");
 *
 * addSection(doc, "Synopsis");
 * addText(doc, "\\fBrazers\\fP [\\fIOPTIONS\\fP] \\fIREFERENCE\\fP \\fIREADS\\fP", false);
 * addText(doc,
 *         "\\fBrazers\\fP [\\fIOPTIONS\\fP] \\fIREFERENCE\\fP \\fILEFT_READS\\fP "
 *         "\\fIRIGHT_READS\\fP", false);
 *
 * addSection(doc, "Description");
 * addText(doc,
 *         "RazerS is a read mapper with controllable, sensitivity.  This "
 *         "means that you can find all read matches in the reference sequence "
 *         "and optionally, you can trade lower sensitivity for better "
 *         "performance.");
 * addText(doc,
 *         "What's special about RazerS is that you can control the sensitivity.");
 *
 * addSection(doc, "Options");
 * addSubSection(doc, "Main Options");
 * addListItem(doc, "\\fB-id\\fP, \\fB--indels\\fP",
 *             "Enable mapping with indels enabled.");
 * addListItem(doc, "\\fB-i\\fP, \\fB--identity\\fP \\fIIDENTITY\\fP",
 *             "Set minimal identity of matches to find.");
 *
 * print(std::cout, doc, "text");
 * @endcode
 *
 * @see ToolDoc#addText
 * @see ToolDoc#addListItem
 * @see ArgumentParser
 */

/*!
 * @fn ToolDoc::ToolDoc
 * @brief Constructor
 *
 * @signature ToolDoc::ToolDoc()
 */

class ToolDoc
{
public:
    CharString _name;
    CharString _shortDescription;
    CharString _date;
    CharString _version;
    CharString _manTitle;
    CharString _category;
    unsigned _manSection;

    String<ToolDocEntry_ *> _entries;

    enum Format {FORMAT_HTML, FORMAT_MAN, FORMAT_TXT};

    ToolDoc() :
        _manSection(1) {}

    ToolDoc(ToolDoc const & toolDoc) :
        _name(toolDoc._name), _shortDescription(toolDoc._shortDescription),
        _date(toolDoc._date), _version(toolDoc._version), _manTitle(toolDoc._manTitle),
        _category(toolDoc._category), _manSection(1)
    {
        append(*this, toolDoc);
    }

    ~ToolDoc()
    {
        clearEntries(*this);
    }

    void print(std::ostream & stream, Format format) const
    {
        switch (format)
        {
        case FORMAT_HTML:
        {
            HtmlToolDocPrinter_ p;
            p.print(stream, *this);
        }
        break;

        case FORMAT_MAN:
        {
            ManToolDocPrinter_ p;
            p.print(stream, *this);
        }
        break;

        case FORMAT_TXT:
        {
            TextToolDocPrinter_ p;
            p.print(stream, *this);
        }
        break;
        }
    }

    void print(std::ostream & stream, CharString const & format) const
    {
        if (format == "html")
            print(stream, FORMAT_HTML);
        else if (format == "man")
            print(stream, FORMAT_MAN);
        else
            print(stream, FORMAT_TXT);
    }

};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function append()                                                  ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#append
 * @headerfile <seqan/arg_parse.h>
 * @brief Append two @link ToolDoc @endlink objects.
 *
 * @signature void append(a, b);
 *
 * @param[in,out] a This object is updated
 * @param[in]     b This object is appended to <tt>b</tt>.
 */

inline void append(ToolDoc & a, ToolDoc const & b)
{
    for (unsigned i = 0; i < length(b._entries); ++i)
    {
        // TODO(holtgrew): This is quite ugly and can be made better with shared_ptr<> in C++11.
        switch (b._entries[i]->getType())
        {
        case ToolDocEntry_::SUBSECTION:
            appendValue(a._entries, new ToolDocSubSection_(*static_cast<ToolDocSubSection_ *>(b._entries[i])));
            break;

        case ToolDocEntry_::SECTION:
            appendValue(a._entries, new ToolDocSection_(*static_cast<ToolDocSection_ *>(b._entries[i])));
            break;

        case ToolDocEntry_::LINE:
            appendValue(a._entries, new ToolDocLine_(*static_cast<ToolDocLine_ *>(b._entries[i])));
            break;

        case ToolDocEntry_::LIST_ITEM:
            appendValue(a._entries, new ToolDocListItem_(*static_cast<ToolDocListItem_ *>(b._entries[i])));
            break;
        }
    }
}

// --------------------------------------------------------------------------
// Function setName()                                                 ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#setName
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the tool name.
 *
 * @signature void setName(toolDoc, name);
 *
 * @param[in,out] toolDoc The ToolDoc object to the set the name for.
 * @param[in]     name    The name of the tool (@link CharString @endlink).
 */

inline void setName(ToolDoc & doc, CharString const & name)
{
    doc._name = name;
}

// --------------------------------------------------------------------------
// Function getName()                                                 ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#getName
 * @headerfile <seqan/arg_parse.h>
 * @brief Get the tool name.
 *
 * @signature CharString getName(toolDoc);
 *
 * @param[in] toolDoc The ToolDoc object to the get the name for.
 *
 * @return CharString Resulting name (@link CharString @endlink).
 */

inline CharString const & getName(ToolDoc const & doc)
{
    return doc._name;
}

// --------------------------------------------------------------------------
// Function setCategory()                                             ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#setCategory
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the tool name.
 *
 * @signature void setName(toolDoc, name);
 *
 * @param[in,out] toolDoc The ToolDoc object to the set the name for.
 * @param[in]     name    The name of the tool (@link CharString @endlink).
 */

inline void setCategory(ToolDoc & doc, CharString const & category)
{
    doc._category = category;
}

// --------------------------------------------------------------------------
// Function getCategory()                                             ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#getCategory
 * @headerfile <seqan/arg_parse.h>
 * @brief Get the tool category.
 *
 * @signature CharString getCategory(toolDoc);
 *
 * @param[in] toolDoc The ToolDoc object to the get the category for.
 *
 * @return CharString Resulting category (@link CharString @endlink).
 */

inline CharString const & getCategory(ToolDoc const & doc)
{
    return doc._category;
}

// --------------------------------------------------------------------------
// Function setShortDescription()                                     ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#setShortDescription
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the tool short description.
 *
 * @signature void setShortDescription(toolDoc, text);
 *
 * @param[in,out] toolDoc The ToolDoc object to the set the short description for.
 * @param[in]     text    The short description of the tool (@link CharString @endlink).
 */

inline void setShortDescription(ToolDoc & doc, CharString const & shortDescription)
{
    doc._shortDescription = shortDescription;
}

// --------------------------------------------------------------------------
// Function getShortDescription()                                     ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#getShortDescription
 * @headerfile <seqan/arg_parse.h>
 * @brief Get the tool short description.
 *
 * @signature CharString getShortDescription(toolDoc);
 *
 * @param[in] toolDoc The ToolDoc object to the get the short description for.
 *
 * @return CharString Resulting short description (@link CharString @endlink).
 */

inline CharString const & getShortDescription(ToolDoc const & doc)
{
    return doc._shortDescription;
}

// --------------------------------------------------------------------------
// Function setDate()                                                 ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#setDate
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the date string.
 *
 * @signature void setName(toolDoc, str);
 *
 * @param[in,out] toolDoc The ToolDoc object to the set the date string for.
 * @param[in]     str     The date string of the tool (@link CharString @endlink).
 */

inline void setDate(ToolDoc & doc, CharString const & date)
{
    doc._date = date;
}

// --------------------------------------------------------------------------
// Function getDate()                                                 ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#getDate
 * @headerfile <seqan/arg_parse.h>
 * @brief Get the date.
 *
 * @signature CharString getDate(toolDoc);
 *
 * @param[in] toolDoc The ToolDoc object to the get the date from.
 *
 * @return CharString Resulting date string (@link CharString @endlink).
 */

inline CharString const & getDate(ToolDoc const & doc)
{
    return doc._date;
}

// --------------------------------------------------------------------------
// Function setVersion()                                              ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#setVersion
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the tool version string.
 *
 * @signature void setName(toolDoc, str);
 *
 * @param[in,out] toolDoc The ToolDoc object to the set the version string for.
 * @param[in]     str     The version string of the tool (@link CharString @endlink).
 */

inline void setVersion(ToolDoc & doc, CharString const & version)
{
    doc._version = version;
}

// --------------------------------------------------------------------------
// Function getVersion()                                              ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#getVersion
 * @headerfile <seqan/arg_parse.h>
 * @brief Get the tool version string.
 *
 * @signature CharString getVersion(toolDoc);
 *
 * @param[in] toolDoc The ToolDoc object to the get the version string.
 *
 * @return CharString Resulting version string (@link CharString @endlink).
 */

inline CharString const & getVersion(ToolDoc const & doc)
{
    return doc._version;
}

// --------------------------------------------------------------------------
// Function setManTitle()                                             ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#setManTitle
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the man title.
 *
 * @signature void setTitle(toolDoc, title);
 *
 * @param[in,out] toolDoc The ToolDoc object to the set the title for.
 * @param[in]     title   The title of the tool (@link CharString @endlink).
 */

inline void setManTitle(ToolDoc & doc, CharString const & title)
{
    doc._manTitle = title;
}

// --------------------------------------------------------------------------
// Function getManTitle()                                             ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#getManTitle
 * @headerfile <seqan/arg_parse.h>
 * @brief Get the tool man page title of.
 *
 * @signature CharString getManTitle(toolDoc);
 *
 * @param[in] toolDoc The ToolDoc object to the get the man page title.
 *
 * @return CharString Resulting man page title (@link CharString @endlink).
 */

inline CharString const & getManTitle(ToolDoc & doc)
{
    return doc._manTitle;
}

// --------------------------------------------------------------------------
// Function addSection()                                              ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#addSection
 * @headerfile <seqan/arg_parse.h>
 * @brief Add a section with the given title.
 *
 * @signature void addSection(toolDoc, title);
 *
 * @param[in,out] toolDoc The ToolDoc object to add a section for.
 * @param[in]     title   The section title (@link CharString @endlink).
 */

inline void addSection(ToolDoc & doc, CharString const & title)
{
    appendValue(doc._entries, new ToolDocSection_(title));
}

// --------------------------------------------------------------------------
// Function addSubSection()                                           ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#addSubSection
 * @headerfile <seqan/arg_parse.h>
 * @brief Add a subsection with the given title.
 *
 * @signature void addSubSection(toolDoc, title);
 *
 * @param[in,out] toolDoc The ToolDoc object to add a subsection for.
 * @param[in]     title   The subsection title (@link CharString @endlink).
 */

inline void addSubSection(ToolDoc & doc, CharString const & title)
{
    appendValue(doc._entries, new ToolDocSubSection_(title));
}

// --------------------------------------------------------------------------
// Function addText()                                                 ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#addText
 * @headerfile <seqan/arg_parse.h>
 * @brief Add a text line/paragraph to ToolDoc.
 *
 * @signature void addText(toolDoc, text[, isParagraph]);
 *
 * @param[in,out] toolDoc     The ToolDoc to add the text to.
 * @param[in]     text        The text to add (@link CharString @endlink).
 * @param[in]     isParagraph Whether to insert as paragraph or just a line (only one line break if not a paragraph).
 */

inline void addText(ToolDoc & doc, CharString const & text, bool isParagraph)
{
    appendValue(doc._entries, new ToolDocLine_(text, isParagraph));
}

inline void addText(ToolDoc & doc, CharString const & text)
{
    addText(doc, text, true);
}

// --------------------------------------------------------------------------
// Function addListItem()                                             ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#addListItem
 * @headerfile <seqan/arg_parse.h>
 * @brief Add a list item to a ToolDoc.
 *
 * @signature void addListItem(toolDoc, key, value);
 *
 * @param[in,out] toolDoc The ToolDoc object to add the list item to.
 * @param[in]     key     The key for the list (@link CharString @endlink).
 * @param[in]     value   The value for the list (@link CharString @endlink).
 */

inline void addListItem(ToolDoc & doc, CharString const & key, CharString const & value)
{
    appendValue(doc._entries, new ToolDocListItem_(key, value));
}

// --------------------------------------------------------------------------
// Function print()                                                   ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#print
 * @headerfile <seqan/arg_parse.h>
 * @brief Print ToolDoc object in a given format.
 *
 * @signature void print(stream, toolDoc, format);
 *
 * @param[in,out] stream  The <tt>std::ostream</tt> to write to.
 * @param[in]     toolDoc The ToolDoc to print.
 * @param[in]     format  The format, one of {"html", "man", "txt"}.
 */

inline void print(std::ostream & stream, ToolDoc const & doc, CharString const & format)
{
    doc.print(stream, format);
}

// --------------------------------------------------------------------------
// Function clearEntries()                                            ToolDoc
// --------------------------------------------------------------------------

/*!
 * @fn ToolDoc#clearEntries
 * @headerfile <seqan/arg_parse.h>
 * @brief Clear entries from ToolDoc.
 *
 * @signature void clearEntries(toolDoc);
 *
 * @param[in,out] toolDoc The ToolDoc object to clear entries from.
 */

inline void clearEntries(ToolDoc & doc)
{
    typedef Iterator<String<ToolDocEntry_ *>, Rooted>::Type TIter;
    for (TIter it = begin(doc._entries, Rooted()); !atEnd(it); goNext(it))
        delete *it;
    clear(doc._entries);
}

// --------------------------------------------------------------------------
// Function HtmlDocPrinter_::print()
// --------------------------------------------------------------------------

inline
void HtmlToolDocPrinter_::print(std::ostream & stream, ToolDoc const & doc)
{
    // Print HTML boilerplate header.
    stream << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" "
           << "http://www.w3.org/TR/html4/strict.dtd\">\n"
           << "<html lang=\"en\">\n"
           << "<head>\n"
           << "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n"
           << "<title>" << xmlEscape(doc._name) << " &mdash; " << xmlEscape(doc._shortDescription) << "</title>\n"
           << "</head>\n"
           << "<body>\n";

    stream << "<h1>" << _toHtml(doc._name) << "</h1>\n"
           << "<div>" << _toHtml(doc._shortDescription) << "</div>\n";

    typedef Iterator<String<ToolDocEntry_ *> const, Rooted>::Type TIter;
    bool isDl = false;
    bool isP = false;
    for (TIter it = begin(doc._entries, Rooted()); !atEnd(it); goNext(it))
    {
        SEQAN_ASSERT_NOT_MSG(isDl && isP, "Current <dl> and <p> are mutually exclusive.");
        ToolDocEntry_ * entry = *it;

        switch (entry->getType())
        {
        case ToolDocEntry_::SECTION:
        {
            _maybeCloseList(stream, isDl);
            _maybeCloseParagraph(stream, isP);
            ToolDocSection_ const * sec = static_cast<ToolDocSection_ const *>(entry);
            stream << "<h2>" << _toHtml(sec->_title) << "</h2>\n";
        }
        break;

        case ToolDocEntry_::SUBSECTION:
        {
            _maybeCloseList(stream, isDl);
            _maybeCloseParagraph(stream, isP);
            ToolDocSection_ const * sec = static_cast<ToolDocSection_ const *>(entry);
            stream << "<h3>" << _toHtml(sec->_title) << "</h3>\n";
        }
        break;

        case ToolDocEntry_::LINE:
        {
            ToolDocLine_ const * line = static_cast<ToolDocLine_ const *>(entry);
            _maybeCloseList(stream, isDl);
            if (!isP)
            {
                stream << "<p>\n";
                isP = true;
            }
            stream << _toHtml(line->_text) << "\n";
            if (line->isParagraph())
                _maybeCloseParagraph(stream, isP);
            else
                stream << "<br />\n";
        }
        break;

        case ToolDocEntry_::LIST_ITEM:
        {
            _maybeCloseParagraph(stream, isP);
            ToolDocListItem_ const * item = static_cast<ToolDocListItem_ const *>(entry);
            if (!isDl)
            {
                stream << "<dl>\n";
                isDl = true;
            }
            stream << "<dt>" << _toHtml(item->_term) << "</dt>\n"
                   << "<dd>" << _toHtml(item->_description) << "</dd>\n";
        }
        break;
        }
    }
    _maybeCloseList(stream, isDl);

    // Print version and date.
    stream << "<h2>Version</h2>\n"
           << "<p>Last update: " << _toHtml(doc._date) << ", " << doc._name
           << " version: " << doc._version << "</p>\n";

    // Print HTML boilerplate footer.
    stream << "</body></html>";
}

// --------------------------------------------------------------------------
// Function TextDocPrinter_::print()
// --------------------------------------------------------------------------

inline
void TextToolDocPrinter_::print(std::ostream & stream, ToolDoc const & doc)
{
    std::ostream_iterator<char> out(stream);

    stream << doc._name;
    if (!empty(doc._shortDescription))
        stream << " - " << doc._shortDescription;
    stream << "\n";
    unsigned len = _textWidth(doc._name) + (empty(doc._shortDescription) ? 0 : 3) + _textWidth(doc._shortDescription);
    std::fill_n(out, len, '=');
    stream << '\n';

    typedef Iterator<String<ToolDocEntry_ *> const, Rooted>::Type TIter;
    bool prevWasParagraph = false;  // Stores whether to add a line break.
    for (TIter it = begin(doc._entries, Rooted()); !atEnd(it); goNext(it))
    {
        ToolDocEntry_ * entry = *it;

        switch (entry->getType())
        {
        case ToolDocEntry_::SECTION:
        {
            ToolDocSection_ const * sec = static_cast<ToolDocSection_ const *>(entry);
            _printSection(stream, *sec);
            prevWasParagraph = false;
        }
        break;

        case ToolDocEntry_::SUBSECTION:
        {
            ToolDocSection_ const * sec = static_cast<ToolDocSection_ const *>(entry);
            _printSubSection(stream, *sec);
            prevWasParagraph = false;
        }
        break;

        case ToolDocEntry_::LINE:
        {
            if (prevWasParagraph)
                stream << '\n';
            ToolDocLine_ const * line = static_cast<ToolDocLine_ const *>(entry);
            _printLine(stream, *line);
            prevWasParagraph = line->isParagraph();
        }
        break;

        case ToolDocEntry_::LIST_ITEM:
        {
            if (prevWasParagraph)
                stream << '\n';
            ToolDocListItem_ const * item = static_cast<ToolDocListItem_ const *>(entry);
            _printListItem(stream, *item);
            prevWasParagraph = false;
        }
        break;
        }
    }

    // Print version and date.
    stream << "\n" << _toText("\\fB") << "VERSION" << _toText("\\fP") << "\n";
    std::fill_n(out, _layout.leftPadding, ' ');
    stream << doc._name << " version: " << doc._version << "\n";
    std::fill_n(out, _layout.leftPadding, ' ');
    stream << "Last update " << doc._date << "\n";
}

inline
void ManToolDocPrinter_::print(std::ostream & stream, ToolDoc const & doc)
{
    std::ostream_iterator<char> out(stream);

    // Print .TH line.
    stream << ".TH ";
    std::transform(begin(doc._name), end(doc._name), out, toupper);
    stream << " " << doc._manSection << " \"" << doc._date << "\" \"";
    std::transform(begin(doc._name), end(doc._name), out, tolower);
    stream << " " << doc._version << "\" \"" << doc._manTitle << "\"\n";

    // Print NAME section.
    stream << ".SH NAME\n"
           << doc._name << " \\- " << doc._shortDescription << std::endl;

    // Print text.
    typedef Iterator<String<ToolDocEntry_ *> const, Rooted>::Type TIter;
    bool isFirstInSection = true;
    for (TIter it = begin(doc._entries, Rooted()); !atEnd(it); goNext(it))
    {
        ToolDocEntry_ * entry = *it;
        switch (entry->getType())
        {
        case ToolDocEntry_::SUBSECTION:
        {
            ToolDocSubSection_ const * sec = static_cast<ToolDocSubSection_ const *>(entry);
            stream << ".SS " << sec->_title << "\n";
            isFirstInSection = true;
        }
        break;

        case ToolDocEntry_::SECTION:
        {
            ToolDocSection_ const * sec = static_cast<ToolDocSection_ const *>(entry);
            stream << ".SH ";
            std::transform(begin(sec->_title), end(sec->_title), out, toupper);
            stream << "\n";
            isFirstInSection = true;
        }
        break;

        case ToolDocEntry_::LINE:
        {
            ToolDocLine_ const * line = static_cast<ToolDocLine_ const *>(entry);
            if (!isFirstInSection && line->isParagraph())
                stream << ".sp\n";
            else if (!isFirstInSection && !line->isParagraph())
                stream << ".br\n";
            stream << line->_text << "\n";
            isFirstInSection = false;
        }
        break;

        case ToolDocEntry_::LIST_ITEM:
        {
            ToolDocListItem_ const * item = static_cast<ToolDocListItem_ const *>(entry);
            stream << ".TP\n"
                   << item->_term << "\n"
                   << item->_description << "\n";
            isFirstInSection = false;
        }
        break;
        }
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_MISC_TOOL_DOC_H_
