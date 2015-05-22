//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2015. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

// File to be used to check for openmp
#include <omp.h>
#include <stdio.h>

int main()
{
	#pragma omp parallel
	printf("Hello from thread %d, nthreads %d\n",
	       omp_get_thread_num(), omp_get_num_threads());
}
