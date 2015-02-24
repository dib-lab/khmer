#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>      // For printing SeqAn Strings.

#include <seqan/arg_parse.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    // Initialize ArgumentParser.
    ArgumentParser parser("arg_parse_demo");
    setCategory(parser, "Demo");
    setShortDescription(parser, "Just a demo of the new ArgumentParser!");
    setVersion(parser, "0.1");
    setDate(parser, "Mar 2012");

    // Add use and description lines.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIIN\\fP \\fIOUT\\fP ");

    addDescription(
        parser,
        "This is just a little demo to show what ArgumentParser is "
        "able to do.  \\fIIN\\fP is a multi-FASTA input file.  \\fIOUT\\fP is a "
        "txt output file.");

    // Add positional arguments and set their valid file types.
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "IN"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUT_FILE, "OUT"));
    setValidValues(parser, 0, "FASTA fa");
    setValidValues(parser, 1, "txt");

    // Add a section with some options.
    addSection(parser, "Important Tool Parameters");
    addOption(parser, ArgParseOption("", "id", "Sequence identity between [0.0:1.0]",
                                     ArgParseArgument::DOUBLE, "ID"));
    setRequired(parser, "id", true);
    setMinValue(parser, "id", "0.0");
    setMaxValue(parser, "id", "1.0");

    // Adding a verbose and a hidden option.
    addSection(parser, "Miscellaneous");
    addOption(parser, ArgParseOption("v", "verbose", "Turn on verbose output."));
    addOption(parser, ArgParseOption("H", "hidden", "Super mysterious flag that will not be shown in "
                                                    "the help screen or man page."));
    hideOption(parser, "H");

    // Add a Reference section.
    addTextSection(parser, "References");
    addText(parser, "http://www.seqan.de");

    // Parse the arguments.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    // Return if there was an error or a built-in command was triggered (e.g. help).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // Extract and print the options.
    bool verbose = false;
    getOptionValue(verbose, parser, "verbose");
    std::cout << "Verbose:     " << (verbose ? "on" : "off") << std::endl;

    double identity = -1.0;
    getOptionValue(identity, parser, "id");
    std::cout << "Identity:    " << identity << std::endl;

    CharString inputFile, outputFile;
    getArgumentValue(inputFile, parser, 0);
    getArgumentValue(outputFile, parser, 1);

    std::cout << "Input-File:  " << inputFile << std::endl;
    std::cout << "Output-File: " << outputFile << std::endl;

    return 0;
}
