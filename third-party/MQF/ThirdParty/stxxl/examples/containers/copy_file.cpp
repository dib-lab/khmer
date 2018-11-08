/***************************************************************************
 *  examples/containers/copy_file.cpp
 *
 *  This example shows three methods to copy a file to another file using
 *  stxxl:vectors.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2004-2006 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/io>
#include <stxxl/vector>
#include <stxxl/stream>

void copy_file(const char* input_path, const char* output_path, unsigned int method)
{
    using stxxl::file;

    file::unlink(output_path); // delete output file

    stxxl::timer tm(true);     // start a timer

    // input file object
    stxxl::syscall_file InputFile(input_path, file::RDONLY | file::DIRECT);
    // output file object
    stxxl::syscall_file OutputFile(output_path, file::RDWR | file::CREAT | file::DIRECT);

    typedef stxxl::vector<unsigned char> vector_type;

    std::cout << "Copying file " << input_path << " to " << output_path << std::endl;

    // InputVector is mapped to InputFile
    vector_type InputVector(&InputFile);
    vector_type OutputVector(&OutputFile);                      // OutputVector is mapped to OutputFile

    std::cout << "File " << input_path << " has size " << InputVector.size() << " bytes." << std::endl;

    if (method == 1)
    {
        // First method: copy vector elements. This is rather slow, because no
        // prefetching is used and the vector's pager does lots of internal
        // work!

        std::cout << "Using first method: copying vector elements." << std::endl;

        for (vector_type::const_iterator it = InputVector.begin(); // iterate through InputVector
             it != InputVector.end(); ++it)
        {
            OutputVector.push_back(*it);                           // add the value pointed by 'it' to OutputVector
        }
    }
    else if (method == 2)
    {
        // Second method: use a vector_iterator2stream to prefetch blocks in
        // the input vector and a vector_bufwriter for buffered writing to the
        // output.

        std::cout << "Using second method: vector_iterator2stream and vector_bufwriter." << std::endl;

        // prepare output vector's size
        OutputVector.resize(InputVector.size());

        // construct prefetching input stream
        stxxl::stream::vector_iterator2stream<vector_type::const_iterator>
        input(InputVector.begin(), InputVector.end());

        // construct buffered output writer
        vector_type::bufwriter_type writer(OutputVector.begin());

        while (!input.empty()) // iterate through InputVector
        {
            writer << *input;
            ++input;
        }
        writer.finish(); // flush buffers
    }
    else if (method == 3)
    {
        // Third method: use a vector_iterator2stream to prefetch blocks in the
        // input vector and materialize to write the stream to the output.

        std::cout << "Using third method: vector_iterator2stream and materialize." << std::endl;

        // prepare output vector's size
        OutputVector.resize(InputVector.size());

        // construct prefetching input stream
        stxxl::stream::vector_iterator2stream<vector_type::const_iterator>
        input(InputVector.begin(), InputVector.end());

        // materilize intput directly into output vector
        stxxl::stream::materialize(input, OutputVector.begin(), OutputVector.end());
    }

    std::cout << "Copied in " << tm.seconds() << " at "
              << (double)InputVector.size() / tm.seconds() / 1024 / 1024 << " MiB/s" << std::endl;
}

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " input_file output_file [method 1-3]" << std::endl;
        return -1;
    }

    int method = (argc >= 4) ? atoi(argv[3]) : 3;

    copy_file(argv[1], argv[2], method);

    return 0;
}
