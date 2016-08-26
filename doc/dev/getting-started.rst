..
   This file is part of khmer, https://github.com/dib-lab/khmer/, and is
   Copyright (C) 2014-2015 Michigan State University
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

Getting started with khmer development
======================================

.. contents::

This document is for people who would like to contribute to khmer.  It
walks first-time contributors through making their own copy of khmer,
building it, and submitting changes for review and merge into the master
copy of khmer.

----

Start by making your own copy of khmer and setting yourself up for
development; then, build khmer and run the tests; and finally, claim
an issue and start developing!

If you're unfamiliar with git and branching in particular, check out
the `git-scm book <http://git-scm.com/book/en/Git-Branching>`__.

We've provided a quick guide to the khmer code base here:
:doc:`codebase-guide`.

One-time Preparation
--------------------

#. Install the dependencies.

   OS X users

   a.  Install Xcode from the `Mac App Store (requires root)
       <https://developer.apple.com/xcode/>`_.
   #.  `Register as an Apple Developer
       <https://developer.apple.com/register>`__.
   #.  Install the Xcode command-line tools: Xcode -> preferences ->
       Downloads -> Command Line Tools (requires root).

   Linux users

   a.  Install the python development environment, virtualenv, pip, gcc, and
       g++.

       On recent Debian and Ubuntu this can be done with::

           sudo apt-get install python2.7-dev python-virtualenv python-pip gcc \
           g++ git astyle gcovr cppcheck enchant

       For RHEL6::

           sudo yum install -y python-devel python-pip git gcc gcc-c++ make enchant
           sudo pip install virtualenv

      For Arch Linux::

          sudo pacman -S python2 python2-pip python2-virtualenv gcc make enchant

#. Get a `GitHub <http://github.com>`__ account.

   (We use GitHub to manage khmer contributions.)

#. Fork `github.com/dib-lab/khmer <https://github.com/dib-lab/khmer>`__.

   Visit that page, and then click on the 'fork' button (upper right).

   (This makes a copy of the khmer source code in your own GitHub account.)

#. Clone your copy of khmer to your local development environment.

   Your clone URL should look something like this::

       https://github.com/empty-titus/khmer.git

   and the UNIX shell command should be::

       git clone https://github.com/empty-titus/khmer.git

   (This makes a local copy of khmer on your development machine.)

#. Add a git reference to the khmer dib-lab repository::

       cd khmer
       git remote add dib https://github.com/dib-lab/khmer.git
       cd ../

   (This makes it easy for you to pull down the latest changes in the
   main repository.)

#. Create a virtual Python environment within which to work with
   `virtualenv <https://pypi.python.org/pypi/virtualenv>`__::

       python2.7 -m virtualenv env

   This gives you a place to install packages necessary for running khmer.

   OS X users and others may need to download virtualenv first::

	curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz
	tar xzf virtualenv*
	cd virtualenv-*; python2.7 virtualenv.py ../env; cd ..

   `Mac ports <https://www.macports.org/>`__ users on the OS X platform can
   install pip by execution from the command line::

       sudo port install py27-pip

   `Homebrew <http://brew.sh/>`__ users on the OS X platform will have pip
   already installed


   `Conda <https://github.com/conda/conda>`__ users on any platform
   should instead create a separate Conda environment::

       conda create -n khmer anaconda

#. Activate the virtualenv and install a few packages::

       source env/bin/activate
       cd khmer
       make install-dependencies

   (This installs `Sphinx <http://sphinx-doc.org/>`__ and `pytest
   <https://pytest.org/>`__, packages we use for
   building the documentation and running the tests.)

   In Conda to activate the previously created environment and install
   dependencies::

       source activate khmer
       cd khmer
       make install-dependencies

#. Cppcheck installation:

   `Debian <https://www.debian.org/>`__ and
   `Ubuntu <http://www.ubuntu.com/>`__ Linux distro users can
   install cppcheck by executing from the command line::

       sudo apt-get install cppcheck

   `Mac ports <https://www.macports.org/>`__ users on the OS X platform can
   install cppcheck by executing from the command line::

       sudo port install cppcheck

   `Homebrew <http://brew.sh/>`__ users on the OS X platform can
   install cppcheck by executing from the command line::

       sudo brew install cppcheck


#. ccache installation:

   Debian and Ubuntu Linux distro users can install ``ccache`` to speed up
   their compile times::

       sudo apt-get install ccache
       echo 'export PATH="/usr/lib/ccache:$PATH" # enable ccache' >> ~/.bashrc
       export PATH="/usr/lib/ccache:$PATH"


Building khmer and running the tests
------------------------------------

#. Activate (or re-activate) the virtualenv::

      source ../env/bin/activate

   ... or for Conda users::

      source activate khmer

   You can run this many times without any ill effects.

   (This puts you in the development environment.)

#. Build khmer::

      make

   If this fails, we apologize -- please `go create a new issue
   <https://github.com/dib-lab/khmer/issues?direction=desc&sort=created&state=open>`__,
   paste in the failure message, and we'll try to help you work through it!

   (This takes the C++ source code and compiles it into something that Python
   can run.)

#. Run the tests::

      make test

   You should see lots of output, with something like::

      ====== 658 passed, 22 deselected in 40.93 seconds =======

   at the end.

   (This will run all of the Python tests in the tests/ directory.)

Congratulations! You're ready to develop!

Claiming an issue and starting to develop
------------------------------------------

#. Find an open issue and claim it.

   Go to `the list of open khmer issues
   <https://github.com/dib-lab/khmer/issues?direction=desc&sort=created&state=open>`__
   and find one you like; we suggest starting with `the low-hanging fruit issues <https://github.com/dib-lab/khmer/issues?direction=desc&labels=low-hanging-fruit&page=1&sort=created&state=open>`__).

   Once you've found an issue you like, make sure that no one has been
   assigned to it (see "assignee", bottom right near "notifications").
   Then, add a comment "I am working on this issue." You've staked
   your claim!

   (We're trying to avoid having multiple people working on the same issue.)

#. In your local copy of the source code, update your master branch
   from the main khmer master branch::

      git checkout master
      git pull dib master

   (This pulls in all of the latest changes from whatever we've been
   doing on dib-lab.)

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

#. Create a new branch and link it to your fork on GitHub::

      git checkout -b fix/brief_issue_description
      git push -u origin fix/brief_issue_description

   where you replace "brief_issue_description" with 2-3 words, separated
   by underscores, describing the issue.

   (This is the set of changes you're going to ask to be merged into khmer.)

#. Make some changes and commit them.

   Though this will largely be issue-dependent the basics of committing are
   simple. After you've made a cohesive set of changes, run the command `git
   status`. This will display a list of all the files git has noticed you
   changed. A file in the 'untracked' section are files that haven't existed
   previously in the repository but git has noticed.

   To commit changes you have to 'stage' them--this is done by issuing the
   following command::

      git add path/to/file

   If you have a large quantity of changes and you don't want to add each file
   manually you can do ``git add --patch`` which will display each set of
   changes to you before staging them for commit.

   Once you have staged your changes, it's time to make a commit::

      git commit

   Git will then open your default console text editor to write a commit
   message -- this is a short (typically 1-3 sentence) description of the
   changes you've made. Please make your commit message informative but
   concise -- these messages become part of the 'official' history of the
   project.

   Once your changes have been committed, push them up to the remote branch::

      git push

   If this is your first commit on a new branch git will error out, telling
   you the remote branch doesn't exist -- This is fine, as it will also provide
   the command to create the branch. Copy/paste/run and you should be set.

   You should also visit and read :doc:`coding-guidelines-and-review`.

#. Periodically update your branch from the main khmer master branch::

      git pull dib master

   (This pulls in all of the latest changes from whatever we've been
   doing on dib-lab - important especially during periods of fast change
   or for long-running pull requests.

#. Run the tests and/or build the docs *before* pushing to GitHub::

      make doc test pep8 diff-cover

   Make sure they all pass!

#. Push your branch to your own GitHub fork::

      git push origin

   (This pushes all of your changes to your own fork.)

#. Repeat until you're ready to merge your changes into "official" khmer.

#. Set up a Pull Request asking to merge things into the central khmer
   repository.

   In a Web browser, go to your GitHub fork of khmer, e.g.::

      https://github.com/empty-titus/khmer

   and you will see a list of "recently pushed branches" just above the
   source code listing.  On the right side of that should be a
   "Compare & pull request" green button.  Click on it!

   Now:

     * add a descriptive title ("updated tests for XXX")
     * put the issue number in the comment ("fixes issue #532")

   then click "Create pull request."

   (This creates a new issue where we can all discuss your proposed
   changes; the khmer team will be automatically notified and you will
   receive e-mail notifications as we add comments.  See `GitHub flow
   <http://scottchacon.com/2011/08/31/github-flow.html>`__ for more
   info.)

#. Paste in the committer checklist from :doc:`coding-guidelines-and-review`
   and, after its pasted in, check off as many of the boxes as you can.

#. As you add new commits to address bugs or formatting issues, you can keep
   pushing your changes to the pull request by doing::

      git push origin

#. When you are ready to have the pull request reviewed, please mention
   @luizirber, @camillescott, @mr-c, or @ctb with a comment 'Ready for review!'

#. The khmer team will now review your pull request and communicate
   with you through the pull request page.  Please feel free to add
   'ping!' and an @ in the comments if you are looking for feedback
   -- this will alert us that you are still on the line -- but we will
   automatically get notified of your pull request and any new
   comments, so use sparingly.

   If this is still your first issue, please *don't* take another issue until
   we've merged your first one - thanks!

#. If we request changes, return to the step "Make some changes and
   commit them" and go from there.  Any additional commits you make and
   push to your branch will automatically be added to the pull request
   (which is pretty dang cool.)

After your first issue is successfully merged...
------------------------------------------------

You're now an experienced GitHub user!  Go ahead and take some more
tasks; you can broaden out beyond the low hanging fruit if you like.

Here are a few suggestions:

* If you're knowledgeable in C++ and/or Python and/or documentation
  and/or biology, we'd love to attract further contributions to khmer.
  Please visit the issues list and browse about and find something
  interesting looking.

* One general thing we'd like to do is increase our test coverage.
  You can go find test coverage information `on our continuous
  integration server
  <http://ci.ged.msu.edu/job/khmer-master/label=linux/cobertura>`__ by
  clicking down to individual files; or, ask us on
  khmer-project@idyll.org for suggestions.

* Ask us! Ask khmer-project@idyll.org for suggestions on what to do next.
  We can suggest particularly ripe low-hanging fruit, or find some other
  issues that suit your interests and background.

* You can also help other people out by watching for new issues or
  looking at pull requests.  Remember to be nice and polite!

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

Now, you can go back to `Claiming an issue and starting to develop`_.

Advanced merging with git-merge-changelog
-----------------------------------------

Often one can get a merge conflict due to updates in the ChangeLog. To teach
Git how to handle these on its own you can install a special merge driver.

On Debian & Ubuntu systems you'll need the `git-merge-changelog` package::

        sudo apt-get install git-merge-changelog

Ubuntu 14.04 LTS users will need to add an external repository that contains a
backport of the package first before installing::

        sudo apt-add-repository ppa:misterc/gedlab
        sudo apt-get update
        sudo apt-get install git-merge-changelog

Everyone should then update their `~/.gitconfig` file with the following::

        [merge "merge-changelog"]
                  name = GNU-style ChangeLog merge driver
                  driver = /usr/bin/git-merge-changelog %O %A %B

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
