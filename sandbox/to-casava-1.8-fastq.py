#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#


import functools
import re
import argparse

from khmer import ReadParser


resub_read_1 = functools.partial(re.sub, r"^(.*)/1$", r"\1 1:N:0:NNNNN")
resub_read_2 = functools.partial(re.sub, r"^(.*)/2$", r"\1 2:N:0:NNNNN")


def setup_cl_parser():

    parser = \
        argparse.ArgumentParser(
            description=
            "Convert the older FASTQ format to the Casava >= 1.8 FASTQ format."
        )
    parser.add_argument("input_filename")
    parser.add_argument("output_filename")

    return parser


def main():

    cl_parser = setup_cl_parser()
    cl_args = cl_parser.parse_args()

    # Note: Only use 1 thread to ensure same ordering of reads.
    rparser = ReadParser(cl_args.input_filename, 1)

    with open(cl_args.output_filename, "w") as output_file:

        for read in rparser:

            new_name = resub_read_1(read.name)
            new_name = resub_read_2(new_name)

            output_file.write(
                "@{name}\n{sequence}\n+\n{accuracy}\n".format(
                    name=new_name,
                    sequence=read.sequence,
                    accuracy=read.accuracy,
                )
            )


if "__main__" == __name__:
    main()

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
