..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2016 The Regents of the University of California.
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

Intermediate-to-advanced khmer development
==========================================

.. contents::



This document is for people who have made contributions to khmer and want some
tips to improve their development process. If you are just getting started with
khmer development, please refer to the introductory guide at
":doc:`getting-started`" and return to this document only when you are
comfortable with all of the topics discussed in the introductory guide.

----


Your second contribution...
---------------------------

Here are a few pointers on getting started on your second (or third,
or fourth, or nth contribution).

So, assuming you've found an issue you'd like to work on there are a
couple things to do to make sure your local copy of the repository is
ready for a new issue--specifically, we need to make sure it's in sync
with the remote repository so you aren't working on a old copy. So::

        git checkout master
        git fetch --all
        git pull

This puts you on the latest master branch and pulls down updates from
GitHub with any changes that may have been made since your last
contribution (usually including the merge of your last
contribution). Then we merge those changes into your local copy of the
master branch.

Now, you can go back to **Claiming an issue and starting to develop** in
:doc:`getting-started`.


Resolving merge conflicts
-------------------------

It is possible that when you do a `git pull` you will get a "merge
conflict" -- This is what happens when something changed in the branch you're
pulling in in the same place you made a change in your local copy. This
frequently happens in the `ChangeLog` file.

Git will complain loudly about merges and tell you specifically in which
files they occurred. If you open the file, you'll see something vaguely
like this in the place where the merge occurred::

   <<<<<<< HEAD
   Changes made on the branch that is being merged into. In most cases,
   this is the branch that you have currently checked out
   =======
   Changes made on the branch that is being merged in, almost certainly
   master.
   >>>>>>> abcde1234

Though there are a variety of tools to assist with resolving merge
conflicts they can be quite complicated at first glance and it is usually
easy enough to manually resolve the conflict.

To resolve the conflict you simply have to manually 'meld' the changes
together and remove the merge markers.

After this you'll have to add and commit the merge just like any other set
of changes. It's also recommended that you run tests.


Virtual environments
--------------------

FIXME FIXME

Pull request cleanup (commit squashing)
---------------------------------------

Submitters are invited to reduce the numbers of commits in their pull requests
either via `git rebase -i dib/master` or this recipe::

        git pull # make sure the local is up to date
        git pull dib master # get up to date
        # fix any merge conflicts
        git status # sanity check
        git diff dib/master # does the diff look correct? (no merge markers)
        git reset --soft dib/master # un-commit the differences from dib/master
        git status # sanity check
        git commit --all # package all differences in one commit
        git status # sanity check
        git push # should fail
        git push --force # override what's in GitHub's copy of the branch/pull request
