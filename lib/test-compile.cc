//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

// Author:  Kevin Murray, spam@kdmurray.id.au
// This file is used to test compilation with libkhmer.a/libkhmer.so

#include  <counting.hh>

int main()
{
    khmer::CountingHash test(1,1);
    return 0;
}
