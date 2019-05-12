/***************************************************************************
 *  lib/common/cmdline.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/cmdline.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/namespace.h>

#include <iomanip>
#include <cstring>

STXXL_BEGIN_NAMESPACE

void cmdline_parser::output_wrap(std::ostream& os, const std::string& text, size_t wraplen,
                                 size_t indent_first, size_t indent_rest,
                                 size_t current, size_t indent_newline)
{
    std::string::size_type t = 0;
    size_t indent = indent_first;
    while (t != text.size())
    {
        std::string::size_type to = t, lspace = t;

        // scan forward in text until we hit a newline or wrap point
        while (to != text.size() &&
               to + current + indent < t + wraplen &&
               text[to] != '\n')
        {
            if (text[to] == ' ') lspace = to;
            ++to;
        }

        // go back to last space
        if (to != text.size() && text[to] != '\n' &&
            lspace != t) to = lspace + 1;

        // output line
        os << std::string(indent, ' ')
           << text.substr(t, to - t) << std::endl;

        current = 0;
        indent = indent_rest;

        // skip over last newline
        if (to != text.size() && text[to] == '\n') {
            indent = indent_newline;
            ++to;
        }

        t = to;
    }
}

void cmdline_parser::print_usage(std::ostream& os)
{
    std::ios state(NULL);
    state.copyfmt(os);

    os << "Usage: " << m_progname
       << (m_optlist.size() ? " [options]" : "");

    for (arglist_type::const_iterator it = m_paramlist.begin();
         it != m_paramlist.end(); ++it)
    {
        const argument* arg = *it;

        os << (arg->m_required ? " <" : " [")
           << arg->m_longkey
           << (arg->m_repeated ? " ..." : "")
           << (arg->m_required ? '>' : ']');
    }

    os << std::endl;

    if (m_description.size())
    {
        os << std::endl;
        output_wrap(os, m_description, m_linewrap);
    }
    if (m_author.size())
    {
        os << "Author: " << m_author << std::endl;
    }

    if (m_description.size() || m_author.size())
        os << std::endl;

    if (m_paramlist.size())
    {
        os << "Parameters:" << std::endl;

        for (arglist_type::const_iterator it = m_paramlist.begin();
             it != m_paramlist.end(); ++it)
        {
            const argument* arg = *it;

            os << "  " << std::setw(m_param_maxlong) << std::left
               << arg->param_text();
            output_wrap(os, arg->m_desc, m_linewrap,
                        0, m_param_maxlong + 2, m_param_maxlong + 2, 8);
        }
    }

    if (m_optlist.size())
    {
        os << "Options:" << std::endl;

        for (arglist_type::const_iterator it = m_optlist.begin();
             it != m_optlist.end(); ++it)
        {
            const argument* arg = *it;

            os << "  " << std::setw(m_opt_maxlong) << std::left
               << arg->option_text();
            output_wrap(os, arg->m_desc, m_linewrap,
                        0, m_opt_maxlong + 2, m_opt_maxlong + 2, 8);
        }
    }

    os.copyfmt(state);
}

void cmdline_parser::print_option_error(
    int argc, const char* const* argv, const argument* arg,
    std::ostream& os)
{
    os << "Error: Argument ";
    if (argc != 0)
        os << '"' << argv[0] << '"';

    os << " for " << arg->type_name() << " option " << arg->option_text()
       << (argc == 0 ? " is missing!" : " is invalid!")
       << std::endl << std::endl;

    print_usage(os);
}

void cmdline_parser::print_param_error(int argc, const char* const* argv, const argument* arg,
                                       std::ostream& os)
{
    os << "Error: Argument ";
    if (argc != 0)
        os << '"' << argv[0] << '"';

    os << " for " << arg->type_name() << " parameter " << arg->param_text()
       << (argc == 0 ? " is missing!" : " is invalid!")
       << std::endl << std::endl;

    print_usage(os);
}

bool cmdline_parser::process(int argc, const char* const* argv, std::ostream& os)
{
    m_progname = argv[0];
    --argc, ++argv;

    // search for help string and output help
    for (int i = 0; i < argc; ++i)
    {
        if (strcmp(argv[i], "-h") == 0 ||
            strcmp(argv[i], "--help") == 0)
        {
            print_usage(os);
            return false;
        }
    }

    arglist_type::iterator argi = m_paramlist.begin(); // current argument in m_paramlist
    bool end_optlist = false;

    while (argc != 0)
    {
        const char* arg = argv[0];

        if (arg[0] == '-' && !end_optlist)
        {
            // option, advance to argument
            --argc, ++argv;
            if (arg[1] == '-')
            {
                if (arg[2] == '-') {
                    end_optlist = true;
                }
                else {
                    // long option
                    arglist_type::const_iterator oi = m_optlist.begin();
                    for ( ; oi != m_optlist.end(); ++oi)
                    {
                        if ((arg + 2) == (*oi)->m_longkey)
                        {
                            if (!(*oi)->process(argc, argv))
                            {
                                print_option_error(argc, argv, *oi, os);
                                return false;
                            }
                            else if (m_verbose_process)
                            {
                                os << "Option " << (*oi)->option_text() << " set to ";
                                (*oi)->print_value(os);
                                os << '.' << std::endl;
                            }
                            break;
                        }
                    }
                    if (oi == m_optlist.end())
                    {
                        os << "Error: Unknown option \"" << arg << "\"." << std::endl << std::endl;
                        print_usage(os);
                        return false;
                    }
                }
            }
            else
            {
                // short option
                if (arg[1] == 0) {
                    os << "Invalid option \"" << arg << '"' << std::endl;
                }
                else {
                    arglist_type::const_iterator oi = m_optlist.begin();
                    for ( ; oi != m_optlist.end(); ++oi)
                    {
                        if (arg[1] == (*oi)->m_key)
                        {
                            if (!(*oi)->process(argc, argv))
                            {
                                print_option_error(argc, argv, *oi, os);
                                return false;
                            }
                            else if (m_verbose_process)
                            {
                                os << "Option " << (*oi)->option_text() << " set to ";
                                (*oi)->print_value(os);
                                os << '.' << std::endl;
                            }
                            break;
                        }
                    }
                    if (oi == m_optlist.end())
                    {
                        os << "Error: Unknown option \"" << arg << "\"." << std::endl << std::endl;
                        print_usage(os);
                        return false;
                    }
                }
            }
        }
        else
        {
            if (argi != m_paramlist.end())
            {
                if (!(*argi)->process(argc, argv))
                {
                    print_param_error(argc, argv, *argi, os);
                    return false;
                }
                else if (m_verbose_process)
                {
                    os << "Parameter " << (*argi)->param_text() << " set to ";
                    (*argi)->print_value(os);
                    os << '.' << std::endl;
                }
                (*argi)->m_found = true;
                if (!(*argi)->m_repeated)
                    ++argi;
            }
            else
            {
                os << "Error: Unexpected extra argument \"" << argv[0] << "\"."
                   << std::endl << std::endl;
                --argc, ++argv;
                print_usage(os);
                return false;
            }
        }
    }

    bool good = true;

    for (arglist_type::const_iterator it = m_paramlist.begin();
         it != m_paramlist.end(); ++it)
    {
        if ((*it)->m_required && !(*it)->m_found) {
            os << "Error: Argument for parameter "
               << (*it)->m_longkey << " is required!" << std::endl;
            good = false;
        }
    }

    if (!good) {
        os << std::endl;
        print_usage(os);
    }

    return good;
}

void cmdline_parser::print_result(std::ostream& os)
{
    std::ios state(NULL);
    state.copyfmt(os);

    int maxlong = STXXL_MAX(m_param_maxlong, m_opt_maxlong);

    if (m_paramlist.size())
    {
        os << "Parameters:" << std::endl;

        for (arglist_type::const_iterator it = m_paramlist.begin();
             it != m_paramlist.end(); ++it)
        {
            const argument* arg = *it;

            os << "  " << std::setw(maxlong) << std::left << arg->param_text();

            std::string typestr = "(" + std::string(arg->type_name()) + ")";
            os << std::setw(m_maxtypename + 4) << typestr;

            arg->print_value(os);

            os << std::endl;
        }
    }

    if (m_optlist.size())
    {
        os << "Options:" << std::endl;

        for (arglist_type::const_iterator it = m_optlist.begin();
             it != m_optlist.end(); ++it)
        {
            const argument* arg = *it;

            os << "  " << std::setw(maxlong) << std::left << arg->option_text();

            std::string typestr = "(" + std::string(arg->type_name()) + ")";
            os << std::setw(m_maxtypename + 4) << std::left << typestr;

            arg->print_value(os);

            os << std::endl;
        }
    }

    os.copyfmt(state);
}

STXXL_END_NAMESPACE
