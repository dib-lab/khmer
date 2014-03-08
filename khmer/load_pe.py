#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#


def _is_1(name):
    return name.endswith('/1')


def _is_2(name):
    return name.endswith('/2')


def is_pe(name1, name2):
    return name1[:-1] == name2[:-1]


def load_pe(screed_handle):
    last_record = None

    screed_iter = iter(screed_handle)

    while 1:
        try:
            this_record = screed_iter.next()
        except StopIteration:
            if last_record:
                yield last_record, None

            raise StopIteration

        if _is_2(this_record.name):
            # PE!
            if last_record:
                if is_pe(last_record.name, this_record.name):
                    yield last_record, this_record
                    last_record = None
                else:
                    # both records exist but they do not match as PEs
                    yield last_record, None
                    yield this_record, None
                    last_record = None

            # first sequence (/1) is missing?
            else:
                yield this_record, None
        else:
            if last_record:
                yield last_record, None
            last_record = this_record
