..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2010-2015 Michigan State University
   Copyright (C) 2015 The Regents of the University of California.
   It is licensed under the three-clause BSD license; see LICENSE.
   Contact: khmer-project@idyll.org

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   Contact: khmer-project@idyll.org

============================
Installing and running khmer
============================


Build requirements
------------------

You'll need a 64-bit operating system, internet access, a C++11 compatible compiler (e.g. GCC 4.8 or greater), GNU Make, and Python version 3.3 or greater.

.. note::

    The khmer package is no longer compatible with Python 2!

.. note::

    If you are running khmer in a HPC environment or for other reasons do not have administrative privileges, we strongly suggest installing khmer in a virtual environment.
    See the relevant instructions below.

.. _user_install_prereqs:

Prerequisites
-------------

OS X
^^^^

#) Start by installing the Xcode command line tools if they are not already installed.

      xcode-select --install

#) If it's not already installed, install Python version 3 using `Homebrew <http://brew.sh/>`_ or `Anaconda <https://www.anaconda.com/download/>`.

Linux
^^^^^

#) Use your Linux distribution's package manager to install Python 3 and essential build tools such as ``Make`` and ``g++``.

   - On recent versions of Debian and Ubuntu this can be done with::

         sudo apt-get install python3-dev python3-venv build-essential

   - For recent versions of Red Hat, Fedora, and CentOS you can invoke::

         sudo yum install -y python3-devel gcc-c++ make


Create a virtual environment
----------------------------

Anaconda Python
^^^^^^^^^^^^^^^

If you are using the Anaconda Python distribution::

    conda create --name khmerEnv python=3.6
    source activate khmerEnv

The first command *creates* the virtual environment in a dedicated location in your home directory and only needs to be invoked once.
The second command *activates* the environment and must be invoked every time you begin a new terminal session.
The *activate* command can be invoked from any directory in your system.

System/Homebrew Python
^^^^^^^^^^^^^^^^^^^^^^

If you are using the system or Homebrew Python distribution::

    python3 -m venv khmerEnv
    source khmerEnv/bin/activate

The first command *creates* the virtual environment in your current directory and only needs to be invoked once.
The second command *activates* the environment and must be invoked every time you begin a new terminal session.
The *activate* command will only work if you provide the correct relative or absolute path of the *activate* file.


Installing khmer
----------------

Once the virtual environment is created and activated, you can use pip to install khmer and its dependencies.

    pip install khmer
